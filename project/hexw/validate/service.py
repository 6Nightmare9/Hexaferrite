"""Compare immutable predictions with explicitly matched observables; never refit."""
from pathlib import Path
import hashlib
import json
import numpy as np
from hexw.io import paths
from hexw.registry.freeze import verify_manifest
from hexw.registry.prior_use import apply_prior_use, load_prior_use
from hexw.registry.schemas import DataRef, Prediction, ValidationResult, Metric, Conclusion
from hexw.registry.versioning import ValidationStore
from hexw.workflow import atomic_json, digest, summary, utc


def compare(run, targets_file):
    root=paths.get_config().root; run=Path(run).resolve()
    frozen=verify_manifest(run/"freeze.json")
    active=verify_manifest(root/"registry"/"FROZEN_MODEL.json")
    if frozen["model_version"]!=active["model_version"]:
        raise ValueError("activate this frozen run before opening its experimental targets")
    with paths.open_data(targets_file,"rb",purpose="validation") as f:
        raw=f.read()
    targets=json.loads(raw)
    if targets.get("template_only"):
        raise ValueError("fill the target template and remove template_only before comparison")
    check_target_origin(targets,targets_file,root)
    req=json.loads((run/"request.json").read_text(encoding="utf-8"))
    mat=json.loads((run/"material.json").read_text(encoding="utf-8"))
    refs=[DataRef(p["record"],sha256=p.get("source_sha256"),group=p.get("measurement_group"),source_role="literature_fit") for p in mat["parameter_sources"]]
    prediction=Prediction(frozen["run_id"],frozen["model_version"],refs)
    prior=load_prior_use(root/"registry") if (root/"registry"/"prior_use.yaml").exists() else {}
    sha=hashlib.sha256(raw).hexdigest(); reports=[]
    with np.load(run/"parameters.npz",allow_pickle=False) as f:
        parameters={k:f[k] for k in f.files}
    for target in targets["targets"]:
        tid=target["id"]
        result=ValidationResult("vr-"+digest({"run":frozen["model_version"],"target":target,"sha":sha})[:20],
            frozen["model_version"],frozen["run_id"],
            [DataRef(tid,sha256=sha,group=target.get("observation_id",tid),source_role="laboratory_validation")],
            Metric("absolute_relative_error",None,"1",.1,"relative"),created_utc=utc())
        apply_prior_use(result,prior)
        try:
            predicted=_match(target,req,parameters,run)
            value=float(target["value"])
            if not np.isfinite(value) or value==0:
                raise ValueError("relative-error metric requires a finite nonzero target; specify another metric for zero")
            if predicted["median"] is None or predicted["n_valid"]<.9*predicted["n_total"]:
                raise ValueError("fewer than 90% valid prediction draws: comparison is conditional/inconclusive")
            err=abs(predicted["median"]-value)/abs(value)
            result.metric.value=float(err)
            result.conclusion=Conclusion.PASS if err<=.1 else Conclusion.FAIL
            result.notes="One matched observable only. PASS does not establish universal <=10% accuracy. Prediction not recalibrated."
            unc=target.get("uncertainty") or {}
            exp_lo,exp_hi=unc.get("lower"),unc.get("upper")
            if "std" in unc:
                sd=float(unc["std"])
                if not np.isfinite(sd) or sd<0:
                    raise ValueError("experimental std must be finite/nonnegative in the declared target unit")
                exp_lo,exp_hi=value-1.96*sd,value+1.96*sd
            interval=None if exp_lo is None or exp_hi is None else [float(exp_lo),float(exp_hi)]
            if interval is not None and (not np.isfinite(interval).all() or interval[0]>interval[1]):
                raise ValueError("invalid experimental interval")
            report={"target":target,"prediction":predicted,"signed_residual":predicted["median"]-value,
                "absolute_error":abs(predicted["median"]-value),"relative_interval_width":(predicted["hi95"]-predicted["lo95"])/abs(value),
                "experimental_interval":interval,
                "intervals_overlap":None if interval is None else max(predicted["lo95"],interval[0])<=min(predicted["hi95"],interval[1]),
                "interval_covers_target":predicted["lo95"]<=value<=predicted["hi95"],
                "compatibility_note":"interval overlap is not an accuracy proof; non-Gaussian scenario distributions are retained"}
        except (ValueError,KeyError,IndexError,TypeError) as e:
            result.conclusion=Conclusion.NOT_EVALUATED; result.metric.value=None
            result.notes=str(e); report={"target":target,"prediction":None}
        store=ValidationStore(root/"runs"/"validation"/"results.jsonl")
        if not store.verify():
            raise ValueError("validation store integrity failure")
        if result.result_id not in {r["record"]["result_id"] for r in store.all()}:
            store.append(result,prediction)
        from hexw.workflow import plain
        report["result"]=plain(result); reports.append(report)
    out=root/"runs"/"validation"/("comparison-"+digest({"run":frozen["model_version"],"sha":sha})[:20])
    out.mkdir(parents=True,exist_ok=True)
    atomic_json(out/"comparison.json",{"schema":"hexw.comparison.v1","model_version":frozen["model_version"],
        "comparisons":reports,"theory_updated":False,"scope":"observable-specific; prior-use disclosure retained"})
    return out


def check_target_origin(document,path,root):
    p=Path(path).resolve();root=Path(root).resolve()
    forbidden=(root/"vault/legacy",root/"runs/predictions",root/"runs/parameter_cache")
    if any(p.is_relative_to(x) for x in forbidden):
        raise ValueError("model predictions and legacy artifacts are not laboratory validation targets")
    if document.get("source_role")!="laboratory_validation" or document.get("derived_from_model"):
        raise ValueError("explicit laboratory_validation source_role and measured provenance required; computed predictions are not measurements")
    if any(t.get("source_role",document["source_role"])!="laboratory_validation" for t in document["targets"]):
        raise ValueError("non-laboratory target role rejected")


def _match(target,req,params,run):
    if target["composition"]!=req["composition"]:
        raise ValueError("composition differs from the frozen request")
    def index(key,val):
        matches=np.flatnonzero(np.isclose(req[key],val,rtol=0,atol=1e-8))
        if len(matches)!=1:
            raise ValueError(f"{key} not uniquely present in frozen prediction")
        return int(matches[0])
    name=target["observable"]; unit=target["unit"]
    ti=index("temperatures_K",target["T_K"])
    if name=="Tc_K":
        data=params["Tc_K"]; expected="K"
    elif name=="Ms_A_m":
        data=params[name][ti]; expected="A/m"
    elif name=="Ms_A_m2_kg":
        data=params[name][ti]; expected="A*m^2/kg"
    elif name in {"K1_J_m3","K2_J_m3","K3_J_m3","K4_J_m3"}:
        if target.get("quantity_definition")!="signed_anisotropy_coefficient":
            raise ValueError("signed K coefficient cannot be compared to an unsigned anisotropy-field magnitude")
        data=params[name][ti]; expected="J/m^3"
    elif name in {"lowest_ground_mode_Hz","peak_Hz","fwhm_Hz"}:
        for key in ("theta_deg","phi_deg","shape","g","alpha"):
            if key not in target or target[key]!=req[key]:
                raise ValueError(f"missing/different condition: {key}")
        if target.get("quantity_definition")!={"lowest_ground_mode_Hz":"undamped_lowest_ground_eigenfrequency",
            "peak_Hz":"external_chi_xx_absorption_peak","fwhm_Hz":"external_chi_xx_absorption_fwhm"}[name]:
            raise ValueError("different resonance/linewidth definition; no automatic peak/derivative/eigenmode substitution")
        hi=index("fields_T",target["mu0H_T"])
        with np.load(run/f"T{ti:03d}_H{hi:04d}.npz",allow_pickle=False) as a:
            data=a[name]
        expected="Hz"
    else:
        raise ValueError("observable has no implemented frozen forward prediction")
    if unit!=expected:
        raise ValueError(f"unit {unit!r} differs from required {expected!r}; convert explicitly")
    return summary(data)
