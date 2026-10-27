"""Queue / checkpoint / resume logic with a fake pw.x (no real calculation): run_external is replaced."""
import json
import subprocess
import threading
import time
from datetime import datetime, timedelta, timezone

import pytest

from hexw.dft import runner as R
from _dftsynth import synthetic_pw_out


class FakeRun:
    def __init__(self, behaviour=None, delay=0.0):
        self.calls, self.active, self.max_active = [], 0, 0
        self.pid = 4242
        self.lock = threading.Lock()
        self.behaviour = behaviour or {}
        self.delay = delay

    def __call__(self, exe, args=(), *, cwd=None, timeout=3600.0, env=None, stdin_path=None, stdout_path=None, stderr_path=None,
                 pid_callback=None):
        name = cwd.name
        assert stdout_path == cwd / "pw.out" and stderr_path == cwd / "pw.err"        # streamed to files, not captured
        if pid_callback is not None:
            pid_callback(self.pid)
        def out(text):
            stdout_path.write_text(text, encoding="utf-8")
            return subprocess.CompletedProcess([], 0, None, None)
        with self.lock:
            self.calls.append(name)
            self.active += 1
            self.max_active = max(self.max_active, self.active)
        try:
            time.sleep(self.delay)
            b = self.behaviour.get(name, "ok")
            assert env["OMP_NUM_THREADS"] == "1" and str(cwd) in str(cwd)
            if b == "timeout":
                stdout_path.write_text(synthetic_pw_out(-1.0, [1.0], converged=False, job_done=False), encoding="utf-8")
                raise subprocess.TimeoutExpired("pw", timeout)
            if b == "crash":
                stderr_path.write_text("boom", encoding="utf-8")
                return subprocess.CompletedProcess([], 1, None, None)
            if b == "unconverged":
                return out(synthetic_pw_out(-1.0, [1.0], converged=False))
            return out(synthetic_pw_out(-5.0, [3.7, -3.7]))
        finally:
            with self.lock:
                self.active -= 1


@pytest.fixture
def jobs(tmp_path):
    out = []
    for k in range(6):
        d = tmp_path / f"c{k:02d}"
        d.mkdir()
        (d / "pw.in").write_text(f"&CONTROL\n/ ! job {k}\n")
        out.append(d)
    return out


@pytest.fixture
def fake(monkeypatch):
    def install(**kw):
        f = FakeRun(**kw)
        monkeypatch.setattr(R, "run_external", f)
        return f
    return install


def test_runs_all_with_bounded_parallelism(cfg, jobs, fake):
    f = fake(delay=0.15)
    out = R.run_queue(cfg, jobs, max_parallel=2, timeout_s=10)
    assert set(out.values()) == {"done"} and len(f.calls) == 6 and f.max_active == 2
    for d in jobs:
        st = R.read_status(d)
        assert st["state"] == "done" and st["attempts"] == 1 and st["config_sha256"] == cfg.sha256
        assert len(st["input_sha256"]) == 64 and len(st["output_sha256"]) == 64 and st["qe_version"] == "7.6"
        assert not (d / "tmp").exists() and (d / "pw.out").is_file()


def test_parallelism_is_capped_by_the_config(cfg, jobs, fake):
    f = fake(delay=0.1)
    R.run_queue(cfg, jobs, max_parallel=99, timeout_s=10)
    assert f.max_active <= cfg.raw["engine"]["max_parallel"]


def test_resume_does_not_touch_finished_jobs(cfg, jobs, fake):
    fake()
    R.run_queue(cfg, jobs[:3], timeout_s=10)
    before = {d.name: (d / "status.json").read_bytes() for d in jobs[:3]}
    f = fake()
    out = R.run_queue(cfg, jobs, timeout_s=10)
    assert sorted(f.calls) == sorted(d.name for d in jobs[3:])
    assert all(out[d.name] == "skip:done" for d in jobs[:3])
    assert {d.name: (d / "status.json").read_bytes() for d in jobs[:3]} == before


def test_failures_are_recorded_and_retried_only_on_request(cfg, jobs, fake):
    f = fake(behaviour={"c01": "crash", "c02": "timeout", "c03": "unconverged"})
    out = R.run_queue(cfg, jobs[:4], timeout_s=5)
    assert out["c00"] == "done" and out["c01"] == "failed" and out["c02"] == "timeout" and out["c03"] == "not_converged"
    assert "rc=1" in R.read_status(jobs[1])["message"] and "timeout" in R.read_status(jobs[2])["message"]
    assert (jobs[2] / "pw.out").is_file()                                      # partial output of the killed run is kept
    f2 = fake()
    out2 = R.run_queue(cfg, jobs[:4], timeout_s=5)
    assert not f2.calls and out2["c01"] == "skip:failed" and out2["c02"] == "skip:timeout"
    out3 = R.run_queue(cfg, jobs[:4], timeout_s=5, retry_failed=True)
    assert out3["c01"] == out3["c02"] == out3["c03"] == "done" and sorted(f2.calls) == ["c01", "c02", "c03"]
    assert R.read_status(jobs[1])["attempts"] == 2


def test_running_job_with_fresh_heartbeat_is_left_alone_stale_is_restarted(cfg, jobs, fake):
    now = datetime.now(timezone.utc)
    R.write_status(jobs[0], state="running", heartbeat_utc=now.isoformat())
    R.write_status(jobs[1], state="running", heartbeat_utc=(now - timedelta(seconds=10_000)).isoformat())
    f = fake()
    out = R.run_queue(cfg, jobs[:2], timeout_s=5)
    assert out["c00"] == "skip:running_elsewhere" and out["c01"] == "done" and f.calls == ["c01"]


def test_changed_input_of_a_finished_job_is_not_silently_rerun(cfg, jobs, fake):
    fake()
    R.run_queue(cfg, jobs[:1], timeout_s=5)
    (jobs[0] / "pw.in").write_text("&CONTROL\n/ ! edited\n")
    f = fake()
    assert R.run_queue(cfg, jobs[:1], timeout_s=5)["c00"] == "skip:input_changed" and not f.calls
    assert R.read_status(jobs[0])["state"] == "input_changed"
    assert R.run_queue(cfg, jobs[:1], timeout_s=5, force=True)["c00"] == "done"


def test_not_feasible_is_never_run(cfg, jobs, fake):
    R.write_status(jobs[0], state="not_feasible_on_laptop")
    f = fake()
    assert R.run_queue(cfg, jobs[:1], timeout_s=5, force=True, retry_failed=True)["c00"] == "skip:not_feasible_on_laptop"
    assert not f.calls


def test_missing_input_is_a_failure(cfg, tmp_path, fake):
    d = tmp_path / "c99"
    d.mkdir()
    fake()
    assert R.run_queue(cfg, [d], timeout_s=5)["c99"] == "failed"


def test_journal_lines_are_json(cfg, jobs, fake, tmp_path):
    fake()
    j = R.Journal(tmp_path / "j" / "journal.jsonl")
    R.run_queue(cfg, jobs[:2], timeout_s=5, journal=j)
    lines = [json.loads(l) for l in (tmp_path / "j" / "journal.jsonl").read_text().splitlines()]
    assert {l["event"] for l in lines} == {"queue", "keep_awake", "start", "finish"}


def test_bad_state_rejected(tmp_path):
    with pytest.raises(ValueError):
        R.write_status(tmp_path, state="weird")


def test_only_run_external_is_used_to_start_programs():
    import hexw.dft.runner as mod
    src = open(mod.__file__, encoding="utf-8").read()
    assert "subprocess" not in src
    assert "run_external" in src and "os.system" not in src and "Popen" not in src
