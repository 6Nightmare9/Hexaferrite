"""Post-freeze literature comparison. Outputs never feed the parameterization."""
import json
import hashlib
from pathlib import Path
import numpy as np
import yaml
from hexw.io import paths
from hexw.registry.freeze import verify_manifest
from hexw.registry.schemas import DataRef,Prediction,ValidationResult,Metric,Conclusion
from hexw.registry.versioning import ValidationStore
from hexw.workflow import atomic_json,digest,summary,utc

PROTOCOL="magnetometry_phasecorrected_v2"


def admissible_magnetization(name,unit):
    # Ordered neutron moments, bulk impure-powder VSM and phase-corrected
    # saturation magnetization must not be treated as the same target.
    return name=="Ms_VSM_impurity_corrected" and unit=="A*m^2/kg"


def compare_reserved(run):
    cfg=paths.get_config();root=cfg.root;run=Path(run)
    frozen=verify_manifest(run/"freeze.json")
    active=verify_manifest(cfg.freeze_file)
    if active["model_version"]!=frozen["model_version"]:
        raise ValueError("activate the matching frozen run first")
    candidates=list(cfg.literature_holdout.glob("records_h*.yaml"))
    if len(candidates)!=1:
        raise ValueError("reserved numerical corpus is ambiguous/missing")
    with paths.open_data(candidates[0],"rb",purpose="validation") as f:
        corpus=f.read()
    corpus_sha=hashlib.sha256(corpus).hexdigest()
    claims=yaml.safe_load(corpus)["claims"]
    req=json.loads((run/"request.json").read_text(encoding="utf-8"))
    material=json.loads((run/"material.json").read_text(encoding="utf-8"))
    with np.load(run/"parameters.npz",allow_pickle=False) as a:
        p={k:a[k] for k in a.files}
    formula="Sr"+"".join(f"{k}{v:g}" for k,v in sorted(req["composition"].items()))+f"Fe{18-sum(req['composition'].values()):g}O27"
    refs=[DataRef(c["record"],group=c.get("measurement_group"),source_role="literature_fit") for c in material["parameter_sources"]]
    pred=Prediction(frozen["run_id"],frozen["model_version"],refs)
    reports=[];skipped=[]
    for c in claims:
        source=c["source"];estimate=c.get("estimate")
        if not estimate or source.get("formula_in_source")!=formula:
            continue
        name=c["name"];unit=estimate["unit"];kind=None;values=None
        cond=source.get("measurement_conditions") or {};T=cond.get("T")
        T=T.get("K") if isinstance(T,dict) else T
        if admissible_magnetization(name,unit) and isinstance(T,(int,float)):
            ix=np.flatnonzero(np.isclose(req["temperatures_K"],T))
            if len(ix)==1 and "Ms_A_m2_kg" in p:
                values=p["Ms_A_m2_kg"][ix[0]];kind="Ms_A_m2_kg"
        elif ("Tc" in name or "T_C" in name) and unit=="K" and "dmax" in (source.get("quantity_definition") or "").lower():
            values=p["Tc_K"];kind="Tc_dmax_K"
        if kind is None:
            if "Ms" in name and unit=="A*m^2/kg":
                old="lit-"+digest({"model":frozen["model_version"],"rec":c["rec_id"],"corpus":corpus_sha})[:20]
                correction=ValidationResult("lit-"+digest({"supersedes":old,"protocol":PROTOCOL})[:20],
                    frozen["model_version"],frozen["run_id"],
                    [DataRef(c["rec_id"],sha256=corpus_sha,group=source.get("measurement_group"),source_role="literature_holdout")],
                    Metric("absolute_relative_error",None,"1"),created_utc=utc(),
                    notes="Supersedes "+old+": mismatched observable definition (neutron ordered moment or uncorrected impure-sample magnetization); no accuracy verdict.")
                store=ValidationStore(root/"runs/validation/literature_results.jsonl")
                if not store.verify():raise ValueError("validation chain changed")
                ids={r["record"]["result_id"] for r in store.all()}
                if old in ids and correction.result_id not in ids:
                    store.append(correction,pred)
                reason=correction.notes if old in ids else "Different observable definition; not evaluated against intrinsic phase-corrected magnetization."
                skipped.append({"record":c["rec_id"],"name":name,"reason":reason})
            continue
        s=summary(values);target=float(estimate["value"])
        err=abs(s["median"]-target)/abs(target) if target else None
        result=ValidationResult("lit-"+digest({"model":frozen["model_version"],"rec":c["rec_id"],"corpus":corpus_sha,"protocol":PROTOCOL})[:20],
            frozen["model_version"],frozen["run_id"],
            [DataRef(c["rec_id"],sha256=corpus_sha,group=source.get("measurement_group"),source_role="literature_holdout")],
            Metric("absolute_relative_error",err,"1",.1,"relative"),
            conclusion=Conclusion.PASS if err is not None and err<=.1 else Conclusion.FAIL,created_utc=utc(),
            notes="Reserved publication comparison; prior developer exposure already disclosed. Retrospective, not a fully blind proof. No fitting performed.")
        store=ValidationStore(root/"runs"/"validation"/"literature_results.jsonl")
        if not store.verify():
            raise ValueError("literature validation chain changed")
        if result.result_id not in {r["record"]["result_id"] for r in store.all()}:
            store.append(result,pred)
        from hexw.workflow import plain
        reports.append({"record":c["rec_id"],"observable":kind,"prediction":s,"observed":target,"unit":unit,
            "record_name":name,"quantity_definition":source.get("quantity_definition"),"measurement_conditions":cond,
            "source_doi":source.get("doi"),"source_locator":source.get("locator"),"source_group_kind":source.get("holdout_kind"),
            "observed_uncertainty":estimate.get("uncertainty"),"interval_covers_observation":s["lo95"]<=target<=s["hi95"],"result":plain(result)})
    out=root/"runs"/"validation"/("literature-"+frozen["run_id"]+"-"+corpus_sha[:8]+"-v2")
    atomic_json(out/"comparison.json",{"model_version":frozen["model_version"],"formula":formula,"comparisons":reports,
        "n_compared":len(reports),"protocol":PROTOCOL,"not_evaluated":skipped,"not_tested":["FMR accuracy without matched reserved FMR measurements", "all unsupported/mismatched definitions"],
        "theory_updated":False,"corpus_sha256":corpus_sha,"validation_mode":"retrospective: earlier developer exposure disclosed"})
    return out
