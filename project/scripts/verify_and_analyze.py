"""Post-freeze checks, reserved-paper comparison, native-data processing. Never refits theory."""
import os
for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[key]="2"
from pathlib import Path
import json,yaml
from hexw.io import paths
from hexw.workflow import friendly_resources,atomic_json,utc
from hexw.registry.freeze import freeze_run
from hexw.experiment.service import analyze
from hexw.validate.literature import compare_reserved

friendly_resources();root=paths.get_config().root
batch=root/"runs/material_series_2026-10-26"
series=json.loads((batch/"checkpoint.json").read_text(encoding="utf-8"))
assert series["status"]=="complete"
log={"started_utc":utc(),"literature":{},"experiment":{},"theory_refitted":False}
for name in ("cozn_co2","cozn_zn2"):
    run=Path(series["completed"][name]["run"]);freeze_run(run)
    try:
        out=compare_reserved(run);d=json.loads((out/"comparison.json").read_text(encoding="utf-8"))
        log["literature"][name]={"out":str(out),"n_compared":d["n_compared"],"comparisons":d["comparisons"]}
        print(name,"reserved comparisons",d["n_compared"],flush=True)
    except Exception as e:
        log["literature"][name]={"error":str(e),"type":type(e).__name__};print(name,str(e),flush=True)
freeze_run(Path(series["completed"]["cozn_mix"]["run"]))
obs=yaml.safe_load((root/"registry/observations.yaml").read_text(encoding="utf-8"))["observations"]
for o in obs:
    if not o["observation_id"].startswith("OBS-") or o["method"] not in ("XRD","Mossbauer"):
        continue
    f=o["files"][0];source=root/f["zone"]/f["archive_path"]
    cfg={"method":"xrd" if o["method"]=="XRD" else "mossbauer","file":str(source),
         "observation_id":o["observation_id"],"expected_sha256":f["sha256"],"options":{}}
    if o["method"]=="XRD":
        cfg["options"]["wavelength_angstrom"]=o["header_metadata"]["wavelength_angstrom"]
    try:
        out=analyze(cfg);d=json.loads((out/"analysis.json").read_text(encoding="utf-8"))
        result={"out":str(out),"status":d.get("status",d.get("method")),"n_peaks":d.get("n_peaks"),"validation_mode":d["validation_mode"]}
    except Exception as e:
        result={"error":str(e),"type":type(e).__name__}
    log["experiment"][o["observation_id"]]=result
    atomic_json(batch/"validation_and_analysis.json",log)
    print(o["observation_id"],result,flush=True)
log["finished_utc"]=utc();atomic_json(batch/"validation_and_analysis.json",log)
print("ANALYSIS COMPLETE",flush=True)
