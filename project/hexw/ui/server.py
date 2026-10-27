"""Loopback-only local UI, one background worker, no external services or GUI takeover."""
from __future__ import annotations
import argparse
import json
import mimetypes
import re
import secrets
import socket
import threading
from concurrent.futures import ThreadPoolExecutor
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse,unquote
from hexw.io import paths
from hexw.workflow import friendly_resources,atomic_json,utc,core_fingerprint


def forecast_progress(message,request):
    match=re.match(r"^T(\d+)_H(\d+):",message)
    if match:
        ti,hi=map(int,match.groups())
        if ti<len(request.temperatures_K) and hi<len(request.fields_T):
            point=ti*len(request.fields_T)+hi+1
            total=len(request.temperatures_K)*len(request.fields_T)
            return (f"Рассчитываю точку {point} из {total}: "
                    f"{request.temperatures_K[ti]:g} К, поле {request.fields_T[hi]:g} Тл").replace('.',',')
    if "parameter ensemble" in message:
        return "Подготавливаю параметры материала. Пожалуйста, подождите."
    if message.startswith("complete:"):
        return "Таблицы и графики сохранены."
    return "Выполняю расчёт. Пожалуйста, подождите."


class Application:
    def __init__(self):
        self.root=paths.get_config().root
        self.current_source=core_fingerprint(self.root)
        self.token=secrets.token_urlsafe(24)
        self.jobs={}
        saved=self.root/"runs"/"ui_jobs.json"
        if saved.exists():
            try:
                self.jobs=json.loads(saved.read_text(encoding="utf-8"))["jobs"]
                for job in self.jobs.values():
                    if job["status"] in {"queued","running"}:
                        job.update(status="interrupted",progress="Работа прервана. Возобновите расчёт из списка результатов.")
            except (ValueError,KeyError):
                self.jobs={}
        self.lock=threading.Lock()
        self.worker=ThreadPoolExecutor(max_workers=1,thread_name_prefix="hexw-calculation")

    def submit(self,kind,payload):
        with self.lock:
            if sum(j["status"] in {"queued","running"} for j in self.jobs.values())>=4:
                raise ValueError("already four queued/running tasks; wait for completion")
            ident=secrets.token_hex(8)
            self.jobs[ident]={"id":ident,"kind":kind,"status":"queued","progress":"В очереди", "created_utc":utc()}
        self.worker.submit(self.execute,ident,kind,payload)
        return self.jobs[ident].copy()

    def execute(self,ident,kind,payload):
        self.update(ident,status="running",progress="Выполняется")
        try:
            if kind in {"predict","resume"}:
                from hexw.theory.request import ForecastRequest,composition_from_formula
                from hexw.theory.pipeline import run_forecast
                if kind=="resume":
                    payload=json.loads((self.run_path(payload["run"])/"request.json").read_text(encoding="utf-8"))
                payload=dict(payload)
                if "formula" in payload:
                    payload["composition"]=composition_from_formula(payload.pop("formula"))
                request=ForecastRequest(**payload)
                out=run_forecast(request,progress=lambda s:self.update(ident,progress=forecast_progress(s,request)))
            elif kind=="analyze":
                from hexw.experiment.service import analyze
                out=analyze(payload)
            elif kind=="freeze":
                from hexw.registry.freeze import freeze_run
                out=self.run_path(payload["run"])
                freeze_run(out)
            elif kind=="validate":
                from hexw.validate.service import compare
                out=compare(self.run_path(payload["run"]),payload["targets_file"])
            else:
                raise ValueError("unknown task")
            self.update(ident,status="complete",progress="Готово",result=out.relative_to(self.root/"runs").as_posix())
        except Exception as e:
            self.update(ident,status="failed",progress=str(e),error_type=type(e).__name__)

    def update(self,ident,**kwargs):
        with self.lock:
            self.jobs[ident].update(kwargs)
            snapshot={"jobs":self.jobs,"updated_utc":utc()}
            atomic_json(self.root/"runs"/"ui_jobs.json",snapshot)

    def run_path(self,rel):
        p=(self.root/"runs"/rel).resolve()
        if not p.is_relative_to((self.root/"runs").resolve()):
            raise ValueError("path escapes results folder")
        return p

    def runs(self):
        results=[]
        for p in (self.root/"runs"/"predictions").glob("*/checkpoint.json"):
            try:
                d=json.loads(p.read_text(encoding="utf-8"))
                summary=p.parent/"summary.json"
                material=json.loads(summary.read_text(encoding="utf-8")) if summary.exists() else {}
                results.append({"run":p.parent.relative_to(self.root/"runs").as_posix(),"run_id":d["run_id"],
                    "label":d["request"].get("label"),"composition":d["request"]["composition"],
                    "status":d["status"],"material_status":material.get("material_status",d["status"]),"created_utc":d.get("created_utc","")})
                results[-1]["source_matches_current"]=d.get("source_fingerprint")==self.current_source
            except (ValueError,KeyError,OSError):
                continue
        return sorted(results,key=lambda r:r["created_utc"],reverse=True)

    def analyses(self):
        results=[]
        for p in (self.root/"runs"/"analysis").glob("*/analysis.json"):
            try:
                d=json.loads(p.read_text(encoding="utf-8"))
                config=d.get("config",{})
                results.append({"run":p.parent.relative_to(self.root/"runs").as_posix(),
                    "label":Path(config["file"]).name if config.get("file") else config.get("observation_id",p.parent.name),
                    "method":d.get("method"),"status":d.get("status","complete"),
                    "_source_key":(d.get("method"),config.get("observation_id") or config.get("file") or str(p)),
                    "created_utc":d.get("created_utc","")})
            except (ValueError,KeyError,OSError,TypeError):
                continue
        results=sorted(results,key=lambda r:r["created_utc"],reverse=True)
        seen=set()
        for result in results:
            key=result.pop("_source_key")
            result["latest_for_source"]=key not in seen
            seen.add(key)
        return results


class LocalServer(ThreadingHTTPServer):
    """One owner per port, including Windows where SO_REUSEADDR can allow duplicates."""
    allow_reuse_address=False

    def server_bind(self):
        if hasattr(socket,"SO_EXCLUSIVEADDRUSE"):
            self.socket.setsockopt(socket.SOL_SOCKET,socket.SO_EXCLUSIVEADDRUSE,1)
        super().server_bind()


def handler(app):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self,*args):
            pass

        def send(self,status,data,content_type="application/json; charset=utf-8"):
            if not isinstance(data,bytes):
                data=json.dumps(data,ensure_ascii=False).encode("utf-8")
            try:
                self.send_response(status)
                self.send_header("Content-Type",content_type)
                self.send_header("Content-Length",str(len(data)))
                self.send_header("X-Content-Type-Options","nosniff")
                self.send_header("Cache-Control","no-store")
                self.end_headers(); self.wfile.write(data)
            except (BrokenPipeError,ConnectionResetError,ConnectionAbortedError):
                # A tab may close or reload while its response is being written.
                self.close_connection=True

        def do_GET(self):
            route=unquote(urlparse(self.path).path)
            if route=="/":
                data=(Path(__file__).parent/"index.html").read_text(encoding="utf-8").replace("__HEXW_TOKEN__",app.token)
                return self.send(200,data.encode("utf-8"),"text/html; charset=utf-8")
            if route=="/api/state":
                with app.lock:
                    jobs=[dict(x) for x in app.jobs.values()]
                return self.send(200,{"jobs":jobs,"runs":app.runs(),"analyses":app.analyses(),"root":str(app.root),"status":"research prototype; <=10% not established"})
            if route.startswith("/results/"):
                try:
                    p=app.run_path(route[len("/results/"):])
                    if p.is_dir():
                        return self.send(200,{"files":[f.name for f in p.iterdir() if f.is_file()]})
                    if p.suffix.lower() not in {".json",".csv",".png",".cif",".npz",".txt"}:
                        return self.send(403,{"error":"unsupported result type"})
                    return self.send(200,p.read_bytes(),mimetypes.guess_type(p.name)[0] or "application/octet-stream")
                except (ValueError,OSError):
                    return self.send(404,{"error":"result not found"})
            self.send(404,{"error":"not found"})

        def do_POST(self):
            if self.headers.get("X-Hexw-Token")!=app.token:
                return self.send(403,{"error":"reload the local interface"})
            try:
                length=int(self.headers.get("Content-Length",0))
                if not 0<length<1024*1024:
                    raise ValueError("invalid request size")
                data=json.loads(self.rfile.read(length))
                if not isinstance(data,dict):
                    raise ValueError("request must be a JSON object")
                kind=urlparse(self.path).path.removeprefix("/api/")
                if kind not in {"predict","resume","analyze","freeze","validate"}:
                    raise ValueError("unknown task")
                self.send(202,app.submit(kind,data))
            except (ValueError,KeyError) as e:
                self.send(400,{"error":str(e)})
    return Handler


def main(args=None):
    ap=argparse.ArgumentParser(prog="hexw serve")
    ap.add_argument("--port",type=int,default=8765)
    ns=ap.parse_args(args)
    friendly_resources()
    app=Application()
    try:
        server=LocalServer(("127.0.0.1",ns.port),handler(app))
    except OSError as e:
        print(f"Не удалось открыть локальный порт {ns.port}: {e}. Если интерфейс уже запущен, откройте http://127.0.0.1:{ns.port}",flush=True)
        app.worker.shutdown(wait=False)
        return 1
    print(f"HexW local interface: http://127.0.0.1:{server.server_port}",flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close(); app.worker.shutdown(wait=True)
    return 0
