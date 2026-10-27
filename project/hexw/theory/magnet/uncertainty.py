"""Monte-Carlo propagation of parameter uncertainty (sequential; no Pool, no I/O).

Distributions are supplied by the caller: scipy.stats frozen distributions (object with .rvs), callables
rng -> float, or plain numbers (constants). Correlated draws: pass `sampler(rng) -> dict`.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Mapping

import numpy as np


@dataclass
class QuantitySummary:
    median: float | None
    lo: float | None
    hi: float | None
    n_valid: int
    n_total: int

    @property
    def fraction_valid(self) -> float:
        return self.n_valid / self.n_total if self.n_total else float("nan")


@dataclass
class UncertaintyResult:
    summary: dict
    samples: dict                 # key -> array with NaN where invalid
    failures: dict                # reason -> count
    n: int
    interval: tuple
    seed: int | None


def _draw(spec, rng):
    if hasattr(spec, "rvs"):
        return float(spec.rvs(random_state=rng))
    if callable(spec):
        return float(spec(rng))
    return float(spec)


def monte_carlo(func: Callable[[dict], Mapping], distributions: Mapping | None = None, *, n: int,
                seed: int | None = None, sampler: Callable | None = None, fixed: Mapping | None = None,
                interval: tuple = (0.025, 0.975)) -> UncertaintyResult:
    """Run func(params) n times; params = fixed + draws. func returns {name: float | None | NaN}.
    None/NaN/non-finite and exceptions (ValueError, ArithmeticError, RuntimeError, LinAlgError) are counted as
    invalid -- reported via n_valid and `failures`, never filled in. Summary: median and the `interval`
    quantiles of the VALID samples."""
    if (distributions is None) == (sampler is None):
        raise ValueError("pass exactly one of distributions or sampler")
    rng = np.random.default_rng(seed)
    rows, fail = [], {}
    for _ in range(n):
        p = dict(fixed or {})
        p.update(sampler(rng) if sampler is not None else {k: _draw(v, rng) for k, v in distributions.items()})
        try:
            out = dict(func(p))
        except (ValueError, ArithmeticError, RuntimeError, np.linalg.LinAlgError) as e:
            fail[type(e).__name__] = fail.get(type(e).__name__, 0) + 1
            rows.append(None)
            continue
        rows.append(out)
    keys = sorted({k for r in rows if r for k in r})
    samples = {k: np.full(n, np.nan) for k in keys}
    for i, r in enumerate(rows):
        if r is None:
            continue
        for k in keys:
            v = r.get(k)
            if v is None or not isinstance(v, (int, float, np.integer, np.floating)) or not math.isfinite(v):
                fail[f"{k}:invalid"] = fail.get(f"{k}:invalid", 0) + 1
            else:
                samples[k][i] = float(v)
    summary = {}
    for k in keys:
        v = samples[k][np.isfinite(samples[k])]
        if len(v) == 0:
            summary[k] = QuantitySummary(None, None, None, 0, n)
        else:
            lo, med, hi = np.quantile(v, [interval[0], 0.5, interval[1]])
            summary[k] = QuantitySummary(float(med), float(lo), float(hi), int(len(v)), n)
    return UncertaintyResult(summary, samples, fail, n, tuple(interval), seed)
