"""Analysis entry point separated from theoretical parameter fitting."""
from pathlib import Path
import hashlib
import json
import numpy as np
import yaml
from hexw.io import paths
from hexw.registry.freeze import verify_manifest
from hexw.workflow import atomic_json, digest, utc
from hexw.experiment import readers


def analyze(config):
    root=paths.get_config().root
    verify_manifest(root/"registry"/"FROZEN_MODEL.json")
    cfg=yaml.safe_load(Path(config).read_text(encoding="utf-8")) if not isinstance(config,dict) else dict(config)
    source=Path(cfg["file"]).resolve()
    # Builtin open is covered by the mandatory hexw audit hook and the experiment
    # frame. Read once, checksum once, then parse the same immutable byte buffer.
    with open(source,"rb") as f:
        raw=f.read()
    if cfg.get("expected_sha256") and hashlib.sha256(raw).hexdigest()!=cfg["expected_sha256"]:
        raise ValueError("measurement checksum differs from its registered source")
    code_hash={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob("*.py")}
    auxiliary={}
    for key in ("structure_cif","mask_file"):
        if cfg.get(key):
            with paths.open_data(cfg[key],"rb",purpose="experiment") as f:
                auxiliary[key]=f.read()
    auxiliary_hash={k:hashlib.sha256(v).hexdigest() for k,v in auxiliary.items()}
    ident=digest({"config":cfg,"source_sha256":hashlib.sha256(raw).hexdigest(),"code":code_hash,"auxiliary":auxiliary_hash})[:20]
    out=root/"runs"/"analysis"/ident
    out.mkdir(parents=True,exist_ok=True)
    done=out/"analysis.json"
    if done.exists():
        prior=json.loads(done.read_text(encoding="utf-8"))
        if any(hashlib.sha256((out/name).read_bytes()).hexdigest()!=sha for name,sha in prior["artifacts"].items()):
            raise ValueError("completed analysis artifact changed")
        return out
    method=cfg["method"].lower(); options=cfg.get("options",{})
    curve=None
    if method=="xrd":
        from hexw.experiment.xrd import analyze as calc
        x,y,header=readers.xrd(source,raw=raw)
        wavelength=options.get("wavelength_angstrom",header.get("wavelength_angstrom"))
        if wavelength is None:
            raise ValueError("supply XRD wavelength_angstrom")
        result=calc(x,y,wavelength_angstrom=wavelength,instrument_fwhm_deg=options.get("instrument_fwhm_deg"))
        if cfg.get("structure_cif"):
            from pymatgen.core import Structure
            from hexw.experiment.xrd import screen_structure
            structure=Structure.from_str(auxiliary["structure_cif"].decode("utf-8"),fmt="cif")
            info,curve=screen_structure(x,y,structure,wavelength)
            result["pattern_screening"]=info
        result["instrument_header"]=header
        xlabel="2 theta (deg)"
    elif method=="mossbauer":
        from hexw.experiment.mossbauer import fit_alternatives
        if source.suffix.lower()==".spc":
            channels,y,meta=readers.spc_channels(source,raw=raw)
            calibration=options.get("velocity_calibration")
            calibrated=bool(calibration and calibration.get("confirmed") and calibration.get("reference"))
            if calibrated:
                x=channels*float(calibration["mm_s_per_channel"])+float(calibration["zero_mm_s"])
                if np.any(np.diff(x)<=0):
                    x,y=x[::-1],y[::-1]
            else:
                x=channels
        else:
            x,y=readers.xy(source,raw=raw)
            calibrated=bool(options.get("calibrated_mm_s") and options.get("calibration_reference"))
            meta={}
        result,curve=fit_alternatives(x,y,calibrated=calibrated,
            n_components=options.get("n_components",[1,2,3]),starts=options.get("starts",3))
        result["instrument_header"]=meta
        xlabel="Velocity (mm/s)" if calibrated else "Uncalibrated channel"
    elif method=="sem":
        from PIL import Image
        from io import BytesIO
        from hexw.experiment.sem import analyze as calc
        image=np.asarray(Image.open(BytesIO(raw)))
        if options.get("roi_px"):
            l,t,r,b=map(int,options["roi_px"]); image=image[t:b,l:r]
        mask=None
        if cfg.get("mask_file"):
            mask=np.asarray(Image.open(BytesIO(auxiliary["mask_file"])).convert("L"))>0
        result,mask=calc(image,um_per_pixel=options["um_per_pixel"],
            foreground=options.get("foreground","bright"),mask=mask,min_area_px=options.get("min_area_px",16))
        Image.fromarray((mask*255).astype("uint8")).save(out/"segmentation.png")
        _image_plot(image,mask,out)
    elif method=="fmr":
        from hexw.theory.magnet.linewidth import fwhm
        from hexw.workflow import plain
        x,y=readers.xy(source,raw=raw)
        if options.get("axis_unit") not in {"GHz","Hz","T","mT","Oe"}:
            raise ValueError("explicit FMR axis_unit required (GHz,Hz,T,mT,Oe)")
        signal=options.get("signal_kind")
        if signal=="absorption":
            if "baseline" not in options:
                raise ValueError("an explicit absorption baseline is required")
            result=plain(fwhm(x,y,baseline=float(options["baseline"]),smooth=None))
        elif signal=="derivative":
            # This file already contains d(absorption)/dx. Do NOT differentiate it again.
            from scipy.signal import find_peaks
            z=y-float(options.get("baseline",0.))
            mx,_=find_peaks(z,prominence=.1*np.ptp(z)); mn,_=find_peaks(-z,prominence=.1*np.ptp(z))
            if len(mx)==len(mn)==1 and mx[0]<mn[0] and z[mx[0]]>0 and z[mn[0]]<0:
                result={"status":"valid","delta_pp":float(x[mn[0]]-x[mx[0]]),
                        "x_max_slope":float(x[mx[0]]),"x_min_slope":float(x[mn[0]])}
            else:
                result={"status":"multipeak_or_censored_or_wrong_sign","delta_pp":None}
        else:
            raise ValueError("signal_kind must be absorption or derivative")
        result.update(axis_unit=options["axis_unit"],signal_kind=signal,
            note="measured peak/width is not an undamped eigenfrequency; no conversion without a declared line shape and field-frequency slope")
        xlabel=options["axis_unit"]
    elif method=="eds":
        import csv,io
        from hexw.experiment.eds import analyze as calc
        content=readers.text(source,raw=raw)
        dialect=csv.Sniffer().sniff(content,delimiters=",;\t")
        rows=list(csv.DictReader(io.StringIO(content),dialect=dialect))
        result=calc(rows,normalize_element=options.get("normalize_element","Sr"))
    else:
        raise ValueError("method must be xrd, mossbauer, sem, eds or fmr")
    if method not in {"sem","eds"}:
        np.savetxt(out/"data_and_fit.csv",np.column_stack([x,y,curve if curve is not None else np.full(len(x),np.nan)]),
                   delimiter=",",header="axis,observed,fit",comments="")
        _curve_plot(x,y,curve,xlabel,out)
    frozen=json.loads((root/"registry"/"FROZEN_MODEL.json").read_text(encoding="utf-8"))
    result.update(schema="hexw.analysis.v1",method=method,created_utc=utc(),config=cfg,
        source_sha256=hashlib.sha256(raw).hexdigest(),model_version_before_access=frozen["model_version"],
        analysis_code_sha256=code_hash,auxiliary_sha256=auxiliary_hash,
        theory_parameters_updated=False)
    from hexw.registry.prior_use import load_prior_use,mode_for_observations
    prior=load_prior_use(root/"registry") if (root/"registry/prior_use.yaml").exists() else {}
    mode,status=mode_for_observations(prior,[cfg.get("observation_id","")])
    result.update(validation_mode=mode.value,prior_use_status=status)
    result["artifacts"]={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file() and p.name!="analysis.json"}
    atomic_json(out/"analysis.json",result)
    return out


def _curve_plot(x,y,fit,xlabel,out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8,4.5)); ax.plot(x,y,lw=.8,label="Observed")
    if fit is not None:
        ax.plot(x,fit,lw=.9,label="Model"); ax.plot(x,y-fit,label="Residual",lw=.6)
    ax.set(xlabel=xlabel,ylabel="Signal / counts"); ax.legend(); fig.tight_layout()
    fig.savefig(out/"analysis.png",dpi=150); plt.close(fig)


def _image_plot(image,mask,out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7,5)); ax.imshow(image,cmap="gray")
    ax.contour(mask,levels=[.5],colors="red",linewidths=.5); ax.axis("off"); fig.tight_layout()
    fig.savefig(out/"overlay.png",dpi=150); plt.close(fig)
