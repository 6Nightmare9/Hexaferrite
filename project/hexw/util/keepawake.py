"""Keep Windows from going to idle sleep while a long calculation queue runs (`with keep_awake(): ...`).

`SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED)` (kernel32 through ctypes) asks the system not to enter
automatic idle sleep for as long as the calling THREAD lives (the flag is per thread: call it from the thread that stays
alive for the whole queue, which `hexw.dft.runner.run_queue` does) or until the context is left, which restores
`ES_CONTINUOUS` alone. The display may still turn off (ES_DISPLAY_REQUIRED is NOT set).

What this does NOT do: it does not change any power plan / `powercfg` setting, and it cannot stop a sleep that the user or the
hardware forces: CLOSING THE LID of a laptop (the default lid action is "sleep") or pressing the power / sleep button still
suspends the machine, and a suspended machine stops the calculation. Long runs need the lid open (or the lid action set to
"do nothing" by the user, which this code never does).

The function opens and reads no file and starts no process: it only calls one kernel32 function (docs/ISOLATION.md). On
platforms other than Windows it does nothing (`active` stays False).
"""
from __future__ import annotations

import sys
from contextlib import contextmanager

ES_CONTINUOUS = 0x80000000
ES_SYSTEM_REQUIRED = 0x00000001


class KeepAwakeState:
    """`active` is True while the request is set; `previous` is the value returned by the first call (0 on failure)."""

    def __init__(self):
        self.active = False
        self.supported = sys.platform == "win32"
        self.error = None


def _set(flags: int) -> int:
    import ctypes
    fn = ctypes.windll.kernel32.SetThreadExecutionState
    fn.restype = ctypes.c_uint32
    fn.argtypes = [ctypes.c_uint32]
    return int(fn(ctypes.c_uint32(flags)))


@contextmanager
def keep_awake(enabled: bool = True):
    """Context manager; yields a KeepAwakeState. Failure to set the request is recorded in `.error` and never raises."""
    st = KeepAwakeState()
    if enabled and st.supported:
        try:
            st.active = _set(ES_CONTINUOUS | ES_SYSTEM_REQUIRED) != 0       # 0 = failure
            if not st.active:
                st.error = "SetThreadExecutionState returned 0"
        except Exception as e:                                              # noqa: BLE001
            st.error = f"{type(e).__name__}: {e}"
    try:
        yield st
    finally:
        if st.active:
            try:
                _set(ES_CONTINUOUS)
            except Exception:                                               # noqa: BLE001
                pass
            st.active = False
