"""Queue, checkpoint/resume and journal for pw.x runs.

One directory per configuration (`<runs>/<campaign>/<composition>/<config_id>/`): `pw.in`, `job.json` (what the job is), `status.json`
(state machine), `pw.out`, `result.json` (after collect). Programs are started ONLY through
`hexw.util.external.run_external` (white list), up to N at a time from threads of this process (1 OpenMP thread each).

States: pending -> running -> done | failed | timeout | not_converged; `input_changed` (a finished run whose pw.in changed
afterwards, never silently re-run), `not_feasible_on_laptop` (SOC jobs; never run), `orphan_running` (a `running` job whose
heartbeat is stale but whose runner process or pw.x process is still alive: it is NOT restarted, two writers in one directory
would corrupt it). Resume: `done` is never touched; `running` with a fresh heartbeat belongs to another runner and is
skipped; with a stale heartbeat it is restarted only if neither the recorded runner PID nor the recorded pw.x PID is alive;
failed / timeout / not_converged are restarted only with `retry_failed`. A restarted run starts the SCF from scratch (a
half-finished pw.x run has no usable restart on this build).

pw.x writes `pw.out` (and `pw.err`) straight into the job directory while it runs (run_external stdout_path): no capture in
memory, the output can be followed with `tail -f`, and it survives a crash of the runner. The queue holds a Windows
keep-awake request (`hexw.util.keepawake`) for the time it runs.
"""
from __future__ import annotations

import json
import os
import shutil
import threading
import time
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Optional

from hexw.dft import parse as P
from hexw.dft import settings as ST
from hexw.util.external import run_external
from hexw.util.keepawake import keep_awake
from hexw.util.proc import pid_alive

STATES = ("pending", "running", "done", "failed", "timeout", "not_converged", "input_changed", "not_feasible_on_laptop",
          "orphan_running")
_LOCKS: dict = {}
_LOCKS_GUARD = threading.Lock()
_ACTIVE: set = set()                      # job directories that THIS process is running right now


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _lock(d: Path) -> threading.Lock:
    with _LOCKS_GUARD:
        return _LOCKS.setdefault(str(d), threading.Lock())


def read_status(d: Path) -> dict:
    p = Path(d) / "status.json"
    if not p.is_file():
        return {"state": "pending"}
    last = None
    for _ in range(40):                       # Windows: a reader can meet an os.replace of the writer; retry briefly
        try:
            with _lock(Path(d)):
                return json.loads(p.read_text(encoding="utf-8"))
        except (PermissionError, FileNotFoundError, json.JSONDecodeError) as e:
            last = e
            time.sleep(0.05)
    return {"state": "failed", "message": f"unreadable status.json: {last!r}"}


def write_status(d: Path, _only_if_state: Optional[tuple] = None, **kw) -> dict:
    """Update status.json atomically. With `_only_if_state` the update is made only if the CURRENT state is one of them (checked
    under the lock: the heartbeat must not overwrite the final state of a job that has just finished)."""
    d = Path(d)
    with _lock(d):
        p = d / "status.json"
        st = json.loads(p.read_text(encoding="utf-8")) if p.is_file() else {"state": "pending"}
        if _only_if_state is not None and st["state"] not in _only_if_state:
            return st
        st.update(kw)
        if st["state"] not in STATES:
            raise ValueError(f"bad state {st['state']!r}")
        tmp = d / "status.json.tmp"
        tmp.write_text(json.dumps(st, indent=1, sort_keys=True), encoding="utf-8")
        for attempt in range(40):             # another process may be reading the file (Windows refuses the replace)
            try:
                os.replace(tmp, p)
                break
            except PermissionError:
                if attempt == 39:
                    raise
                time.sleep(0.05)
    return st


class Journal:
    """Append-only JSON-lines journal of a composition directory."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self._lock = threading.Lock()

    def log(self, event: str, **kw) -> None:
        rec = {"t": now_utc(), "event": event, **kw}
        with self._lock:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            with open(self.path, "a", encoding="utf-8") as f:
                f.write(json.dumps(rec, sort_keys=True) + "\n")


def qe_env(cfg: ST.DftConfig) -> dict:
    env = dict(os.environ)
    env["PATH"] = str(cfg.bin_dir()) + os.pathsep + env.get("PATH", "")
    env["OMP_NUM_THREADS"] = str(int(cfg.raw["engine"]["omp_threads"]))
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["OMP_STACKSIZE"] = "256M"
    return env


def job_input_sha(d: Path) -> Optional[str]:
    p = Path(d) / "pw.in"
    return ST.sha256_file(p) if p.is_file() else None


def run_job(cfg: ST.DftConfig, d: Path, timeout_s: float, journal: Optional[Journal] = None) -> dict:
    """Run pw.x in directory d (blocking). Returns the final status."""
    d = Path(d)
    pw_in = d / "pw.in"
    if not pw_in.is_file():
        return write_status(d, state="failed", message="pw.in missing", finished_utc=now_utc())
    attempts = int(read_status(d).get("attempts", 0)) + 1
    in_sha = ST.sha256_file(pw_in)
    write_status(d, state="running", started_utc=now_utc(), heartbeat_utc=now_utc(), attempts=attempts, input_sha256=in_sha,
                 config_sha256=cfg.sha256, runner_pid=os.getpid(), pw_pid=None, message="", finished_utc=None)
    if journal:
        journal.log("start", job=d.name, attempt=attempts)
    t0 = time.perf_counter()
    rc, state, msg = None, "failed", ""
    out_p, err_p = d / "pw.out", d / "pw.err"
    _ACTIVE.add(str(d))
    try:
        cp = run_external(cfg.executable_path(), [], cwd=d, timeout=timeout_s, env=qe_env(cfg), stdin_path=pw_in,
                          stdout_path=out_p, stderr_path=err_p,
                          pid_callback=lambda pid: write_status(d, pw_pid=pid))
        rc = cp.returncode
    except Exception as e:                                              # noqa: BLE001
        if type(e).__name__ == "TimeoutExpired":
            state, msg = "timeout", f"timeout after {timeout_s:.0f} s"
        else:
            state, msg = "failed", f"{type(e).__name__}: {e}"
    finally:
        _ACTIVE.discard(str(d))
    wall = time.perf_counter() - t0
    if err_p.is_file() and err_p.stat().st_size == 0:
        err_p.unlink()
    stdout = out_p.read_text(encoding="utf-8", errors="replace") if out_p.is_file() else ""
    res = P.parse_pw_output(stdout)
    out_sha = ST.sha256_file(out_p) if out_p.is_file() else None
    if state != "timeout" and rc is not None:
        if rc == 0 and res.job_done and res.converged:
            state = "done"
        elif rc == 0 and res.job_done:
            state, msg = "not_converged", "SCF did not converge"
        else:
            state, msg = "failed", f"rc={rc}; " + ("; ".join(res.errors) if res.errors else "no JOB DONE in the output")
    if state == "done" and not cfg.raw["dft"].get("keep_scratch", False):
        shutil.rmtree(d / "tmp", ignore_errors=True)
    st = write_status(d, state=state, message=msg, rc=rc, wall_s=round(wall, 2), finished_utc=now_utc(),
                      output_sha256=out_sha, qe_version=res.qe_version, n_iterations=res.n_iterations,
                      energy_Ry=res.energy_Ry, pw_pid=None)
    if journal:
        journal.log("finish", job=d.name, state=state, wall_s=round(wall, 1), message=msg)
    return st


def alive_processes(d: Path, st: dict) -> list:
    """Descriptions of the processes recorded in the status of a `running` job that are still alive (runner, pw.x). The runner of
    this very process counts only while it is really running that directory (a stale record of a dead thread does not)."""
    out = []
    rp, wp = st.get("runner_pid"), st.get("pw_pid")
    if rp is not None:
        if int(rp) == os.getpid():
            if str(Path(d)) in _ACTIVE:
                out.append(f"runner pid {rp} (this process)")
        elif pid_alive(rp):
            out.append(f"runner pid {rp}")
    if wp is not None and int(wp) != os.getpid() and pid_alive(wp):
        out.append(f"pw.x pid {wp}")
    return out


def classify(cfg: ST.DftConfig, d: Path, retry_failed: bool = False, force: bool = False, now: Optional[float] = None) -> str:
    """'run' | 'skip:<reason>' for a job directory."""
    st = read_status(d)
    s = st["state"]
    if s == "not_feasible_on_laptop":
        return "skip:not_feasible_on_laptop"
    if s == "done":
        cur = job_input_sha(d)
        if cur is not None and st.get("input_sha256") not in (None, cur) and not force:
            write_status(d, state="input_changed", message="pw.in changed after the run finished")
            return "skip:input_changed"
        return "skip:done"
    if s == "input_changed":
        return "run" if force else "skip:input_changed"
    if s in ("running", "orphan_running"):
        hb = st.get("heartbeat_utc")
        age = (datetime.now(timezone.utc) - datetime.fromisoformat(hb)).total_seconds() if hb else float("inf")
        if s == "running" and age < float(cfg.raw["runtime"]["stale_after_s"]):
            return "skip:running_elsewhere"
        alive = alive_processes(d, st)
        if alive:                                  # a stale heartbeat is also what a sleeping laptop produces: never restart a live job
            write_status(d, state="orphan_running",
                         message="heartbeat is stale but " + ", ".join(alive) + " is still alive: not restarted "
                                 "(wait for it or kill it; for a recycled PID delete status.json)")
            return "skip:orphan_running"
        return "run"
    if s in ("failed", "timeout", "not_converged"):
        return "run" if retry_failed else f"skip:{s}"
    return "run"


def run_queue(cfg: ST.DftConfig, job_dirs: list, *, max_parallel: Optional[int] = None, timeout_s: float,
              retry_failed: bool = False, force: bool = False, journal: Optional[Journal] = None,
              progress: Optional[Callable[[str], None]] = None) -> dict:
    """Run the not yet finished jobs, at most max_parallel at a time, holding a keep-awake request. Returns {job: final state or
    'skip:...'}."""
    n = int(max_parallel or cfg.raw["engine"]["max_parallel"])
    n = max(1, min(n, int(cfg.raw["engine"]["max_parallel"])))
    say = progress or (lambda s: None)
    outcome: dict = {}
    todo = []
    for d in map(Path, job_dirs):
        c = classify(cfg, d, retry_failed, force)
        if c == "run":
            todo.append(d)
        else:
            outcome[d.name] = c
            say(f"{d.name}: {c}")
    if journal:
        journal.log("queue", n_parallel=n, to_run=[d.name for d in todo], skipped=outcome, runner_pid=os.getpid())
    if not todo:
        return outcome
    with keep_awake(bool(cfg.raw["runtime"].get("keep_awake", True))) as ka:
        if journal:
            journal.log("keep_awake", active=ka.active, supported=ka.supported, error=ka.error)
        say(f"keep-awake: {'active' if ka.active else 'NOT active'}" + (f" ({ka.error})" if ka.error else ""))
        _run_todo(cfg, todo, n, timeout_s, journal, say, outcome)
    return outcome


def _run_todo(cfg, todo, n, timeout_s, journal, say, outcome) -> None:
    hb = float(cfg.raw["runtime"]["heartbeat_s"])
    ex = ThreadPoolExecutor(max_workers=n)
    futs: dict = {}
    try:
        for d in todo:
            futs[ex.submit(run_job, cfg, d, timeout_s, journal)] = d
            # the status must show 'running' only for jobs that really started; submit order == start order
        pending = set(futs)
        while pending:
            done, pending = wait(pending, timeout=hb, return_when=FIRST_COMPLETED)
            for f in done:
                d = futs[f]
                try:
                    outcome[d.name] = f.result()["state"]
                except Exception as e:                                  # noqa: BLE001
                    write_status(d, state="failed", message=f"{type(e).__name__}: {e}", finished_utc=now_utc())
                    outcome[d.name] = "failed"
                say(f"{d.name}: {outcome[d.name]}")
            for f in pending:                                           # heartbeat of the running ones
                d = futs[f]
                write_status(d, ("running", "orphan_running"), state="running", heartbeat_utc=now_utc())
    except KeyboardInterrupt:
        ex.shutdown(wait=False, cancel_futures=True)
        if journal:
            journal.log("interrupted")
        raise
    finally:
        ex.shutdown(wait=True)
