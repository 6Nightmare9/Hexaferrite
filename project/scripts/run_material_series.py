"""Reproducible sequential laptop benchmark; can be rerun after interruption."""
import os
for name in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[name]="2"
import json
from pathlib import Path
from hexw.workflow import friendly_resources,atomic_json,utc,core_fingerprint
from hexw.io import paths
from hexw.theory.request import ForecastRequest
from hexw.theory.pipeline import run_forecast
from hexw.registry.freeze import freeze_run

friendly_resources()
root=paths.get_config().root
batch=root/"runs"/"material_series_2026-10-26"
batch.mkdir(parents=True,exist_ok=True)
state={"started_utc":utc(),"source_fingerprint":core_fingerprint(root),"completed":{},"status":"running"}
names=["cozn_co2","cozn_mix","cozn_zn2","mn1","mn2","mg1","mg2","ni1","ni2"]
for name in names:
    print("SERIES "+name,flush=True)
    req=ForecastRequest.load(root/"project"/"config"/(name+".yaml"))
    out=run_forecast(req,progress=lambda s:print(s,flush=True))
    frozen=freeze_run(out)
    state["completed"][name]={"run":str(out),"model_version":frozen["model_version"],"finished_utc":utc()}
    atomic_json(batch/"checkpoint.json",state)
    p=root/"audit"/"PROJECT_STATE.json"
    project=json.loads(p.read_text(encoding="utf-8"))
    project["active_work"]={"stage":"representative material series", "completed":list(state["completed"]),"last_run":str(out),"next":"finish sequential series, experimental smoke checks, delivery"}
    atomic_json(p,project)
state["status"]="complete"; state["finished_utc"]=utc()
atomic_json(batch/"checkpoint.json",state)
print("SERIES COMPLETE",flush=True)
