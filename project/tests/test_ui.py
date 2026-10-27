"""Offline UI routing/security tests. No browser or external network access."""
from pathlib import Path
import io
import json
import pytest
from hexw.ui.server import Application
from hexw.theory.request import composition_from_formula


@pytest.mark.parametrize("f,c",[("SrCoZnFe16O27",{"Co":1,"Zn":1}),("SrMn₂Fe₁₆O₂₇",{"Mn":2}),("SrFe18O27",{})])
def test_formula_parser(f,c):
    assert composition_from_formula(f)==c


@pytest.mark.parametrize("f",["BaNi2Fe16O27","SrNi3Fe15O27","SrAlFe17O27","SrFe16O26","Sr(CoZn)Fe16O27"])
def test_formula_outside_domain_rejected(f):
    with pytest.raises(ValueError):
        composition_from_formula(f)


def test_ui_result_paths_cannot_escape_root(root):
    app=Application()
    try:
        assert app.run_path("predictions/test")==root/"runs"/"predictions"/"test"
        with pytest.raises(ValueError):
            app.run_path("../../vault/lab_validation/d.spc")
    finally:
        app.worker.shutdown()


def test_ui_analysis_failure_is_reported_without_worker_crash(root):
    app=Application()
    try:
        app.jobs["test"]={"status":"queued"}
        app.execute("test","analyze",{"file":"missing.txt","method":"xrd"})
        assert app.jobs["test"]["status"]=="failed"
        assert "FROZEN_MODEL" in app.jobs["test"]["progress"]
    finally:
        app.worker.shutdown()


def test_html_contains_four_workflows_and_local_csrf_token():
    import hexw.ui.server as s
    text=(Path(s.__file__).parent/"index.html").read_text(encoding="utf-8")
    for kind in ("predict","analyze","freeze","validate"):
        assert f"post('{kind}'" in text
    assert "__HEXW_TOKEN__" in text and "X-Hexw-Token" in text


def test_ui_lists_saved_analyses_and_ignores_incomplete_json(root):
    out=root/"runs/analysis/record"
    out.mkdir(parents=True)
    (out/"analysis.json").write_text(json.dumps({"method":"mossbauer","status":"calibration_required",
        "created_utc":"2026-10-27T00:00:00Z","config":{"observation_id":"OBS-MOS-001"}}),encoding="utf-8")
    bad=root/"runs/analysis/partial";bad.mkdir()
    (bad/"analysis.json").write_text('{',encoding="utf-8")
    app=Application()
    try:
        records=app.analyses()
        assert len(records)==1
        assert records[0]["run"]=="analysis/record"
        assert records[0]["label"]=="OBS-MOS-001"
        assert records[0]["status"]=="calibration_required"
    finally:
        app.worker.shutdown()


def test_analysis_history_retains_both_versions_and_marks_latest(root):
    for ident,date in [("earlier","2026-10-26T20:00:00Z"),("later","2026-10-27T00:00:00Z")]:
        out=root/"runs/analysis"/ident;out.mkdir(parents=True)
        (out/"analysis.json").write_text(json.dumps({"method":"xrd","created_utc":date,
            "config":{"observation_id":"OBS-XRD-001"}}),encoding="utf-8")
    app=Application()
    try:
        records=app.analyses()
        assert [r["run"] for r in records]==["analysis/later","analysis/earlier"]
        assert [r["latest_for_source"] for r in records]==[True,False]
        assert records[0]["status"]=="complete"
    finally:
        app.worker.shutdown()


def test_progress_explains_actual_point_and_hides_internal_paths():
    from types import SimpleNamespace
    from hexw.ui.server import forecast_progress
    request=SimpleNamespace(temperatures_K=[298.15,400.0],fields_T=[0.0,0.4,0.8])
    assert forecast_progress("T000_H0000: T=298.15 K, mu0 H=0 T",request)=="Рассчитываю точку 1 из 6: 298,15 К, поле 0 Тл"
    assert forecast_progress("T001_H0002: T=400 K, mu0 H=0.8 T",request)=="Рассчитываю точку 6 из 6: 400 К, поле 0,8 Тл"
    for message in ["run abc: parameter ensemble (48 draws)","complete: C:/private/result", "unfamiliar internal message"]:
        shown=forecast_progress(message,request)
        assert shown and "abc" not in shown and "C:/" not in shown and "internal" not in shown


def test_local_ui_has_one_port_owner():
    # Binding only: no requests, clients, browser access or external network.
    from http.server import BaseHTTPRequestHandler
    from hexw.ui.server import LocalServer
    first=LocalServer(("127.0.0.1",0),BaseHTTPRequestHandler)
    try:
        with pytest.raises(OSError):
            second=LocalServer(first.server_address,BaseHTTPRequestHandler)
            second.server_close()
    finally:
        first.server_close()


def _response_handler(writer):
    from hexw.ui.server import handler
    response=object.__new__(handler(None))
    response.request_version="HTTP/1.1"
    response.requestline="GET /api/state HTTP/1.1"
    response.command="GET"
    response.close_connection=False
    response.wfile=writer
    return response


def test_ui_response_keeps_json_and_headers():
    writer=io.BytesIO()
    response=_response_handler(writer)
    response.send(200,{"ready":True})
    headers,body=writer.getvalue().split(b"\r\n\r\n",1)
    assert headers.startswith(b"HTTP/1.0 200 OK\r\n")
    assert b"Content-Type: application/json; charset=utf-8" in headers
    assert f"Content-Length: {len(body)}".encode() in headers
    assert b"X-Content-Type-Options: nosniff" in headers
    assert json.loads(body)=={"ready":True}
    assert not response.close_connection


@pytest.mark.parametrize("failure_at",[1,2],ids=["headers","body"])
@pytest.mark.parametrize("error",[BrokenPipeError,ConnectionResetError,ConnectionAbortedError])
def test_ui_response_accepts_disconnected_client(failure_at,error):
    class DisconnectedWriter:
        writes=0
        def write(self,data):
            self.writes+=1
            if self.writes==failure_at:
                raise error("client disconnected")
            return len(data)
    writer=DisconnectedWriter()
    response=_response_handler(writer)
    response.send(200,{"ready":True})
    assert response.close_connection
    assert writer.writes==failure_at


def test_ui_response_does_not_hide_other_write_errors():
    class FailedWriter:
        def write(self,data):
            raise OSError("unexpected write failure")
    response=_response_handler(FailedWriter())
    with pytest.raises(OSError,match="unexpected write failure"):
        response.send(200,{"ready":True})
