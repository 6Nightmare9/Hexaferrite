"""Driver of the k-point test (docs/DFT.md section 12): five pw.x calculations through the CLI `python -m hexw dft ... run`, at most
3 at a time (memory: Gamma 0.92 GB, 2x2x1 2.67 GB per process, QE estimates). Every calculation is one CLI call with --only, so every one
is checkpointed by the DFT queue (status.json, pw.out streamed into the job directory, keep-awake held by the CLI process).

Started detached (Win32_Process.Create), survives the end of the launching session. Safe to start again: finished jobs are skipped by the
queue, running ones are recognised through the recorded PIDs (orphan_running).
"""
import datetime
import os
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor

ROOT = r"D:\Work\Hexx"
PROJ = ROOT + r"\project"
PY = ROOT + r"\.venv\Scripts\python.exe"
CFG = PROJ + r"\config\dft\ktest_v1.yaml"
LOGDIR = ROOT + r"\dft_runs\ktest_v1\logs"
MAX_PARALLEL = 3
C00 = "c00_ground"
C01 = "c01_flip_W3_W4_W5_W6_W7"
# longest first: 2x2x1 (~5 h) start immediately, Gamma (~1.3 h) fill the remaining slot and follow
JOBS = [
    ("k221_c00", ["--kpoints", "2,2,1", "--tag", "k221"], C00),
    ("k221_c01", ["--kpoints", "2,2,1", "--tag", "k221"], C01),
    ("g_c00", ["--kpoints", "1,1,1", "--tag", "g"], C00),
    ("g_c01", ["--kpoints", "1,1,1", "--tag", "g"], C01),
    ("g_alt_c00", ["--kpoints", "1,1,1", "--tag", "g_alt", "--fe2-start", "none"], C00),
]
_lock = threading.Lock()


def log(msg):
    line = f"{datetime.datetime.now().isoformat(timespec='seconds')}  {msg}\n"
    with _lock, open(os.path.join(LOGDIR, "driver.log"), "a", encoding="utf-8") as f:
        f.write(line)


def run_one(job):
    name, opts, cid = job
    env = dict(os.environ, HEXW_ROOT=ROOT)
    cmd = [PY, "-m", "hexw", "dft", "--config", CFG, "run", "--system", "w", *opts, "--only", cid, "--i-accept-runtime", "--parallel", "1"]
    log(f"start {name}: {' '.join(cmd)}")
    with open(os.path.join(LOGDIR, f"{name}.log"), "ab") as out:
        rc = subprocess.run(cmd, cwd=PROJ, env=env, stdout=out, stderr=subprocess.STDOUT, check=False).returncode
    log(f"finish {name}: rc={rc}")
    return name, rc


def main():
    os.makedirs(LOGDIR, exist_ok=True)
    log(f"driver pid {os.getpid()} started; {len(JOBS)} jobs, at most {MAX_PARALLEL} at a time")
    with ThreadPoolExecutor(max_workers=MAX_PARALLEL) as ex:
        res = list(ex.map(run_one, JOBS))
    log(f"driver done: {res}")


if __name__ == "__main__":
    main()
