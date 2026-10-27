"""Process pool with a per-task timeout and a clear error when a worker dies (instead of endless respawning).

THE RULES (read before using multiprocessing in hexw)
  1. The 'spawn' start method re-imports the worker's module in every child. A module-level statement that touches
     data (opens a file in vault/ or holdout, reads a big CSV, starts another pool, ...) runs in EVERY child at
     import time. If it fails, the child dies before it can serve a task; stock `multiprocessing.Pool` then spawns
     a replacement that dies the same way, forever. => DO NOT TOUCH THE VAULT (or any data) AT IMPORT TIME.
  2. Create the pool only under `if __name__ == "__main__":` (or inside a function that is itself called from there).
     Without that guard the child re-executes the pool creation: `RuntimeError: An attempt has been made to start a
     new process before the current process has finished its bootstrapping phase`.
  3. Worker functions must be importable top-level functions (picklable); they get data through their arguments.
  Children of a hexw process inherit HEXW_ROOT/HEXW_GUARD and install the data-guard hook before any user code
  (see hexw/io/paths.py); the vault/holdout rules therefore apply inside workers too.

`WorkerPool` is a thin wrapper over `concurrent.futures.ProcessPoolExecutor` (spawn context):
  * a worker that crashes (also at import) -> `WorkerCrashedError` with the hint above, the pool is torn down
    (the stdlib executor does not respawn, it marks itself broken);
  * a task that runs longer than `timeout` seconds -> `TaskTimeoutError`; the whole pool is killed (a running task
    cannot be cancelled), pending tasks are dropped;
  * start-up handshake (`start()`): every worker must answer a ping within `startup_timeout`, so an import-time
    failure is reported before the first real task.
"""
from __future__ import annotations

import concurrent.futures as cf
import multiprocessing as mp
import os
import time
from concurrent.futures.process import BrokenProcessPool
from typing import Any, Callable, Iterable

__all__ = ["WorkerPool", "PoolError", "WorkerCrashedError", "TaskTimeoutError", "RULE"]

RULE = ("do not touch the vault (or any data) at import time of a worker module; create the pool under "
        "`if __name__ == '__main__':`")


class PoolError(RuntimeError):
    pass


class WorkerCrashedError(PoolError):
    """A worker process died (typically: its module failed or exited at import, or the task segfaulted/os._exit)."""


class TaskTimeoutError(PoolError, TimeoutError):
    """A task did not finish within the timeout."""


def _ping(_=None) -> int:
    return os.getpid()


def _crash_hint(extra: str = "") -> str:
    return ("a worker process died. Typical cause: the worker's module raised or called exit/touched forbidden data "
            f"at import time. Rule: {RULE}." + (" " + extra if extra else ""))


class WorkerPool:
    def __init__(self, workers: int | None = None, *, timeout: float = 120.0, startup_timeout: float = 60.0,
                 initializer: Callable | None = None, initargs: tuple = ()):
        if timeout <= 0 or startup_timeout <= 0:
            raise ValueError("timeouts must be positive")
        self.workers = max(1, int(workers or min(4, os.cpu_count() or 1)))
        self.timeout = float(timeout)
        self.startup_timeout = float(startup_timeout)
        self._init = (initializer, initargs)
        self._ex: cf.ProcessPoolExecutor | None = None

    # ---- lifecycle ----
    def start(self) -> "WorkerPool":
        if self._ex is not None:
            return self
        ini, args = self._init
        self._ex = cf.ProcessPoolExecutor(self.workers, mp_context=mp.get_context("spawn"),
                                          initializer=ini, initargs=args)
        # the executor starts its workers lazily: ping them all to learn about import-time failures NOW
        futs = [self._ex.submit(_ping) for _ in range(self.workers)]
        try:
            done, pending = cf.wait(futs, timeout=self.startup_timeout)
            if pending:
                raise TaskTimeoutError(f"workers did not start within {self.startup_timeout} s. Rule: {RULE}.")
            for f in done:
                f.result()
        except BrokenProcessPool as e:
            self._kill()
            raise WorkerCrashedError(_crash_hint("(detected at start-up)")) from e
        except BaseException:
            self._kill()
            raise
        return self

    def _kill(self):
        ex, self._ex = self._ex, None
        if ex is None:
            return
        for p in list((getattr(ex, "_processes", None) or {}).values()):
            try:
                p.kill()
            except Exception:                                   # noqa: already dead
                pass
        try:
            ex.shutdown(wait=False, cancel_futures=True)
        except Exception:                                       # noqa
            pass

    def close(self):
        ex, self._ex = self._ex, None
        if ex is not None:
            ex.shutdown(wait=True, cancel_futures=True)

    def __enter__(self):
        return self.start()

    def __exit__(self, et, ev, tb):
        if et is None:
            self.close()
        else:
            self._kill()
        return False

    # ---- work ----
    def map(self, fn: Callable[[Any], Any], items: Iterable[Any], *, timeout: float | None = None) -> list:
        """Results in input order. At most `workers` tasks are in flight, so the clock of a task starts when it
        really starts; each task gets `timeout` seconds. Exceptions raised by `fn` are re-raised unchanged."""
        if self._ex is None:
            self.start()
        tmo = self.timeout if timeout is None else float(timeout)
        items = list(items)
        results: list = [None] * len(items)
        queue = list(range(len(items)))[::-1]
        inflight: dict[cf.Future, tuple[int, float]] = {}
        try:
            while queue or inflight:
                while queue and len(inflight) < self.workers:
                    i = queue.pop()
                    inflight[self._ex.submit(fn, items[i])] = (i, time.monotonic() + tmo)
                nearest = min(d for _, d in inflight.values())
                done, _ = cf.wait(list(inflight), timeout=max(0.0, nearest - time.monotonic()),
                                  return_when=cf.FIRST_COMPLETED)
                for f in done:
                    i, _d = inflight.pop(f)
                    results[i] = f.result()
                now = time.monotonic()
                late = [i for f, (i, d) in inflight.items() if d <= now and not f.done()]
                if late:
                    raise TaskTimeoutError(f"task(s) {late} (item index) exceeded {tmo} s; the pool was killed")
        except BrokenProcessPool as e:
            self._kill()
            raise WorkerCrashedError(_crash_hint()) from e
        except BaseException:
            self._kill()
            raise
        return results
