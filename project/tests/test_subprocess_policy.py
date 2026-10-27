"""Subprocess policy for the theory-side packages (stage 4A): static check in `hexw audit` + the white-listed launcher."""
import subprocess
import sys
from pathlib import Path

import pytest

import hexw
from hexw import audit
from hexw.util import external

BAD = {
    "subprocess_run_executable": "import subprocess, sys\nsubprocess.run([sys.executable, 'a.py'])\n",
    "popen_python_string": "import subprocess\nsubprocess.Popen('python -c pass', shell=True)\n",
    "os_system": "import os\nos.system('python a.py')\n",
    "alias_module": "import subprocess as sp\nsp.check_output(['python', 'a.py'])\n",
    "from_import": "from subprocess import Popen as P\nP(['x'])\n",
    "from_os": "from os import system as s, popen\npopen('x')\n",
    "os_exec": "import os, sys\nos.execv(sys.executable, [sys.executable])\n",
    "os_spawn": "import os\nos.spawnl(os.P_WAIT, 'x', 'x')\n",
    "dunder_import": "m = __import__('subprocess')\nm.run(['x'])\n",
    "import_module": "import importlib\nimportlib.import_module('subprocess')\n",
    "getattr_system": "import os\ngetattr(os, 'system')('x')\n",
    "mp_process": "import multiprocessing as mp\nmp.Process(target=print).start()\n",
    "mp_pool": "from multiprocessing import Pool\nPool(2)\n",
    "executor": "import concurrent.futures\nconcurrent.futures.ProcessPoolExecutor(2)\n",
    "executor_from": "from concurrent.futures import ProcessPoolExecutor as P\nP()\n",
    "asyncio": "import asyncio\nasyncio.create_subprocess_exec('x')\n",
    "syntax_error": "def (:\n",
}
GOOD = {
    "worker_pool": "from hexw.util.pool import WorkerPool\nwith WorkerPool(2) as p:\n    p.map(abs, [1])\n",
    "run_external": "from hexw.util.external import run_external\nrun_external('pw.x', ['-in', 'a.in'])\n",
    "plain": "import os, sys\nprint(os.getcwd(), sys.executable)\nos.path.join('a', 'b')\n",
    "threads": "from concurrent.futures import ThreadPoolExecutor\nThreadPoolExecutor(2)\n",
    "unrelated_run": "class R:\n    def run(self): ...\nR().run()\n",
}


@pytest.mark.parametrize("name", sorted(BAD))
def test_static_check_flags_child_process_calls(name):
    assert audit.find_subprocess_violations(BAD[name]), name


@pytest.mark.parametrize("name", sorted(GOOD))
def test_static_check_allows_the_sanctioned_wrappers(name):
    assert audit.find_subprocess_violations(GOOD[name]) == [], name


def test_python_in_the_call_is_named():
    (_, msg), = audit.find_subprocess_violations("import subprocess, sys\nsubprocess.run([sys.executable, 'a.py'])\n")
    assert "sys.executable" in msg
    (_, msg), = audit.find_subprocess_violations("import subprocess\nsubprocess.run(['python3', 'a.py'])\n")
    assert "python" in msg


def test_theory_side_packages_are_clean():
    root = Path(hexw.__file__).resolve().parent
    assert audit.scan_theory_sources(root) == []


def test_scan_finds_a_violation_in_a_synthetic_package(tmp_path):
    for pkg in audit.THEORY_PACKAGES:
        (tmp_path / pkg).mkdir()
    (tmp_path / "dft" / "ok.py").write_text(GOOD["run_external"])
    (tmp_path / "params" / "sub").mkdir()
    (tmp_path / "params" / "sub" / "bad.py").write_text(BAD["subprocess_run_executable"])
    (tmp_path / "util").mkdir()
    (tmp_path / "util" / "wrapper.py").write_text(BAD["subprocess_run_executable"])      # outside the scanned packages
    hits = audit.scan_theory_sources(tmp_path)
    assert [h[0].replace("\\", "/") for h in hits] == ["params/sub/bad.py"]


def test_audit_reports_the_policy_checks():
    rep = audit.Report()
    audit.check_subprocess_policy(rep, None)
    assert len(rep.checks) == 3 and rep.ok, rep.checks


# ---------------------------------------------------------------- run_external
@pytest.mark.parametrize("exe", ["python", sys.executable, "python.exe", "cmd", "cmd.exe", "powershell", "sh", "bash",
                                 "", " pw.x", "pw.x\n", "../python"])
def test_run_external_refuses_interpreters_and_shells(exe):
    with pytest.raises(external.ExternalProgramError):
        external.run_external(exe, [])


def test_run_external_white_list_and_missing_binary():
    assert "pw.x" in external.ALLOWED_EXECUTABLES and "python" not in external.ALLOWED_EXECUTABLES
    with pytest.raises(external.ExternalProgramError, match="not found"):
        external.resolve_executable("pw.x") if external.shutil.which("pw.x") is None else pytest.skip("pw.x installed")


class FakePopen:
    """Stands in for subprocess.Popen: records the call, writes into the stdout / stderr FILES the way a child does."""
    calls: list = []

    def __init__(self, argv, **kw):
        self.argv, self.kw, self.pid, self.returncode = argv, kw, 777, 0
        FakePopen.calls.append(self)
        if kw["stdout"] not in (subprocess.PIPE, None):
            kw["stdout"].write(b"line 1\n")
            kw["stdout"].flush()                                    # visible on disk while the "child" is still running
        if kw["stderr"] not in (subprocess.PIPE, None):
            kw["stderr"].write(b"warn\n")

    def communicate(self, timeout=None):
        return (b"captured" if self.kw["stdout"] is subprocess.PIPE else None,
                b"" if self.kw["stderr"] is subprocess.PIPE else None)

    def kill(self):
        pass

    def wait(self):
        return 0


def _fake_pw(monkeypatch, tmp_path):
    fake = tmp_path / "pw.x"
    fake.write_text("")
    FakePopen.calls = []
    monkeypatch.setattr(external.subprocess, "Popen", FakePopen)
    monkeypatch.setattr(external.shutil, "which", lambda s: str(fake))
    return fake


def test_run_external_never_uses_a_shell_and_validates_arguments(monkeypatch, tmp_path):
    fake = _fake_pw(monkeypatch, tmp_path)
    r = external.run_external(str(fake), ["-in", "a.in"], cwd=tmp_path, timeout=5)
    c = FakePopen.calls[0]
    assert r.stdout == "captured" and c.argv == [str(fake), "-in", "a.in"]                  # default: captured as text
    assert c.kw["shell"] is False and c.kw["cwd"] == str(tmp_path) and c.kw["stdout"] is subprocess.PIPE
    with pytest.raises(external.ExternalProgramError):
        external.run_external(str(fake), "-in a.in")                    # one command-line string
    with pytest.raises(external.ExternalProgramError):
        external.run_external(str(fake), ["-c", "import os"])
    with pytest.raises(external.ExternalProgramError):
        external.run_external(str(fake), [], timeout=0)
    # `which` resolving to something else than a white-listed name is refused
    monkeypatch.setattr(external.shutil, "which", lambda s: str(tmp_path / "python.exe"))
    with pytest.raises(external.ExternalProgramError):
        external.run_external("pw.x", [])


def test_run_external_streams_stdout_to_a_file_without_capture(monkeypatch, tmp_path):
    """stdout_path: the child writes into the file while it runs; nothing is captured in memory; no shell; pid reported."""
    fake = _fake_pw(monkeypatch, tmp_path)
    seen = {}
    out, err, inp = tmp_path / "pw.out", tmp_path / "pw.err", tmp_path / "pw.in"
    inp.write_text("&CONTROL\n/\n")
    r = external.run_external(str(fake), [], cwd=tmp_path, timeout=5, stdin_path=inp, stdout_path=out, stderr_path=err,
                              pid_callback=lambda pid: seen.update(pid=pid, size=out.stat().st_size))
    c = FakePopen.calls[0]
    assert c.kw["shell"] is False and c.kw["stdout"] is not subprocess.PIPE and c.kw["stderr"] is not subprocess.PIPE
    assert r.stdout is None and r.stderr is None and r.returncode == 0                      # not captured
    assert out.read_bytes() == b"line 1\n" and err.read_bytes() == b"warn\n"
    assert seen["pid"] == 777 and seen["size"] > 0                                          # the file already had output at the start
    assert c.kw["stdout"].closed and c.kw["stdin"].closed                                   # handles closed afterwards


def test_run_external_source_has_no_capture_output_or_shell():
    import inspect
    src = inspect.getsource(external.run_external)
    assert "capture_output" not in src and "shell=True" not in src and "shell=False" in src
