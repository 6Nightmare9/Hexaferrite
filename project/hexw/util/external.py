"""The ONLY sanctioned way for theory-side code (hexw.theory, params, structures, dft, spin) to start an external program.

Policy (docs/ISOLATION.md, "Subprocess policy"):
  * Python child processes are started only through `hexw.util.pool.WorkerPool` (spawn, with the data guard installed
    in the child before any user code);
  * every other program only through `run_external`, which accepts the executables of a white list (the Quantum
    ESPRESSO `pw.x`) and never a shell, an interpreter (`python`, `sys.executable`, `cmd`, `powershell`, `sh`, ...) or
    a command line given as one string;
  * `hexw.audit` checks the theory-side sources statically and fails on any other use of subprocess / os.system / Popen.

The data guard of `hexw.io.paths` still applies: paths to the vault / holdout in the arguments or the cwd are refused
by the audit hook (the call comes from a theory-side frame).
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Callable, Optional, Sequence

__all__ = ["ALLOWED_EXECUTABLES", "ExternalProgramError", "resolve_executable", "run_external"]

# base names (case-insensitive, with or without .exe) of the executables that may be started; extend deliberately
ALLOWED_EXECUTABLES = frozenset({"pw.x", "pw.exe", "pw.x.exe", "pw_big.exe"})   # pw_big.exe: pw.exe with a larger stack (Windows QE build)
_FORBIDDEN_ARG_START = ("-c",)          # never allowed as an argument even to a white-listed program


class ExternalProgramError(RuntimeError):
    pass


def resolve_executable(exe: str | os.PathLike) -> str:
    """Absolute path of a white-listed executable (name looked up in PATH, or an explicit path)."""
    if not isinstance(exe, (str, os.PathLike)) or not str(exe).strip():
        raise ExternalProgramError("executable must be a non-empty string / path")
    s = str(exe)
    if any(ch in s for ch in "\r\n\t") or s != s.strip():
        raise ExternalProgramError(f"bad executable string {s!r}")
    base = Path(s).name.lower()
    if base not in ALLOWED_EXECUTABLES:
        raise ExternalProgramError(
            f"executable {base!r} is not on the white list {sorted(ALLOWED_EXECUTABLES)}. Python children go through "
            "hexw.util.pool.WorkerPool; other programs must be added to the list deliberately.")
    found = shutil.which(s)
    if found is None:
        raise ExternalProgramError(f"white-listed executable {base!r} not found (PATH / given path): {s!r}")
    if Path(found).name.lower() not in ALLOWED_EXECUTABLES:             # `which` resolved to something else
        raise ExternalProgramError(f"{found!r} is not a white-listed executable")
    return found


def run_external(exe: str | os.PathLike, args: Sequence[str] = (), *, cwd: Optional[str | os.PathLike] = None,
                 timeout: float = 3600.0, env: Optional[dict] = None, stdin_path: Optional[str | os.PathLike] = None,
                 stdout_path: Optional[str | os.PathLike] = None, stderr_path: Optional[str | os.PathLike] = None,
                 pid_callback: Optional[Callable[[int], None]] = None) -> subprocess.CompletedProcess:
    """Run a white-listed program without a shell. Raises on timeout (the child is killed) and does not raise on a non-zero
    exit code (the caller decides).

    Output: without `stdout_path` / `stderr_path` the stream is captured as text (CompletedProcess.stdout / .stderr). With
    `stdout_path` the child writes straight into that file WHILE it runs (opened "wb", i.e. truncated; nothing is held in
    memory, a long run can be followed with `tail`, and the file survives a timeout / a crash of the caller); the same for
    `stderr_path`. In that case the corresponding CompletedProcess attribute is None. `pid_callback(pid)` is called right after
    the start (the queue records the child PID for the orphan check)."""
    if isinstance(args, (str, bytes)):
        raise ExternalProgramError("args must be a list of strings, not one command-line string")
    argv = [resolve_executable(exe)]
    for a in args:
        if not isinstance(a, (str, os.PathLike)):
            raise ExternalProgramError(f"argument {a!r} is not a string / path")
        if str(a) in _FORBIDDEN_ARG_START:
            raise ExternalProgramError(f"argument {str(a)!r} is not allowed")
        argv.append(str(a))
    if timeout is None or not timeout > 0:
        raise ExternalProgramError("a positive timeout is required")
    opened = []
    try:
        stdin = open(stdin_path, "rb") if stdin_path is not None else None      # noqa: SIM115 - closed below
        if stdin is not None:
            opened.append(stdin)
        out = err = subprocess.PIPE
        if stdout_path is not None:
            out = open(stdout_path, "wb")                                       # noqa: SIM115
            opened.append(out)
        if stderr_path is not None:
            err = open(stderr_path, "wb")                                       # noqa: SIM115
            opened.append(err)
        proc = subprocess.Popen(argv, shell=False, cwd=None if cwd is None else str(cwd), env=env, stdin=stdin,
                                stdout=out, stderr=err, text=False)
        if pid_callback is not None:
            try:
                pid_callback(proc.pid)
            except Exception:                                                   # noqa: BLE001 - bookkeeping must not kill the run
                pass
        try:
            so, se = proc.communicate(timeout=timeout)
        except subprocess.TimeoutExpired as e:
            proc.kill()
            so, se = proc.communicate()
            e.stdout = so
            e.stderr = se
            raise
        except BaseException:
            proc.kill()
            proc.wait()
            raise
        dec = lambda b: None if b is None else b.decode("utf-8", "replace")     # noqa: E731
        return subprocess.CompletedProcess(argv, proc.returncode, dec(so), dec(se))
    finally:
        for f in opened:
            f.close()
