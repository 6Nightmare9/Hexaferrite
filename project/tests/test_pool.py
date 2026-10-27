"""hexw.util.pool: timeout per task and a clear error instead of endless respawning (synthetic root only)."""
import sys
import textwrap
import time

import pytest

from hexw.util import pool as P


@pytest.fixture
def mods(tmp_path, monkeypatch, root):
    (tmp_path / "poolmods").mkdir()
    d = tmp_path / "poolmods"
    (d / "pm_ok.py").write_text(textwrap.dedent('''
        import time
        def square(x): return x * x
        def slow(x):
            time.sleep(x)
            return x
        def boom(x): raise ValueError("bad item %s" % x)
        def die(x):
            import os
            os._exit(7)
    '''))
    # a module that fails at import time in CHILD processes only (the parent can import it to pickle references)
    (d / "pm_badimport.py").write_text(textwrap.dedent('''
        import multiprocessing, os
        if multiprocessing.current_process().name != "MainProcess":
            os._exit(3)            # simulates "touched forbidden data at import" -> child dies
        def f(x): return x
    '''))
    monkeypatch.syspath_prepend(str(d))
    yield
    for n in ("pm_ok", "pm_badimport"):
        sys.modules.pop(n, None)


def test_map_in_order(mods):
    import pm_ok
    with P.WorkerPool(2, timeout=60) as pool:
        assert pool.map(pm_ok.square, [3, 1, 2]) == [9, 1, 4]


def test_task_timeout_is_reported(mods):
    import pm_ok
    t0 = time.monotonic()
    with pytest.raises(P.TaskTimeoutError):
        with P.WorkerPool(1, timeout=1.0) as pool:
            pool.map(pm_ok.slow, [30])
    assert time.monotonic() - t0 < 25                    # killed, not waited for


def test_exception_of_task_propagates(mods):
    import pm_ok
    with pytest.raises(ValueError, match="bad item"):
        with P.WorkerPool(1, timeout=30) as pool:
            pool.map(pm_ok.boom, [1])


def test_worker_death_during_task(mods):
    import pm_ok
    with pytest.raises(P.WorkerCrashedError, match="import time"):
        with P.WorkerPool(1, timeout=30) as pool:
            pool.map(pm_ok.die, [1])


def test_import_time_crash_is_detected_at_start(mods):
    import pm_badimport
    # initializer lives in the failing module -> the child dies while importing it
    t0 = time.monotonic()
    with pytest.raises(P.WorkerCrashedError) as ei:
        P.WorkerPool(2, startup_timeout=60, initializer=pm_badimport.f, initargs=(1,)).start()
    assert "__main__" in str(ei.value) or P.RULE in str(ei.value)
    assert time.monotonic() - t0 < 50


def test_rule_is_documented():
    assert "import" in P.__doc__ and "__main__" in P.__doc__
