"""Is a process with a given PID alive? (read-only query; used by the DFT queue before it restarts a `running` job.)

Windows: OpenProcess(PROCESS_QUERY_LIMITED_INFORMATION) + GetExitCodeProcess == STILL_ACTIVE (kernel32 through ctypes; it
reads no file and starts no process). NOT `os.kill(pid, 0)`: on Windows that terminates the process. Elsewhere: `os.kill(pid, 0)`.
A recycled PID is indistinguishable from the original process: the caller treats "alive" as "do not restart", the safe side.
"""
from __future__ import annotations

import os
import sys

_PROCESS_QUERY_LIMITED_INFORMATION = 0x1000
_STILL_ACTIVE = 259


def pid_alive(pid) -> bool:
    try:
        pid = int(pid)
    except (TypeError, ValueError):
        return False
    if pid <= 0:
        return False
    if sys.platform == "win32":
        import ctypes
        from ctypes import wintypes
        k32 = ctypes.WinDLL("kernel32", use_last_error=True)
        k32.OpenProcess.restype = wintypes.HANDLE
        k32.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
        k32.GetExitCodeProcess.argtypes = [wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD)]
        k32.CloseHandle.argtypes = [wintypes.HANDLE]
        h = k32.OpenProcess(_PROCESS_QUERY_LIMITED_INFORMATION, False, pid)
        if not h:
            return ctypes.get_last_error() == 5                                  # 5 = access denied: the process exists
        try:
            code = wintypes.DWORD()
            ok = k32.GetExitCodeProcess(h, ctypes.byref(code))
            return bool(ok) and code.value == _STILL_ACTIVE
        finally:
            k32.CloseHandle(h)
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True
