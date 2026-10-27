"""Access to registered laboratory measurements (registry/observations.yaml), prior-analysis artefacts excluded.

Prior analyses (role prior_lab_analysis: XRD fit tables, plots, .mtd projects, reports) are derived artefacts, not
measurements. `iter_measurements()` returns ONLY original measurements unless include_prior=True is passed
explicitly. It reads the registry only (no vault bytes).
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterator

import yaml

from hexw.io import paths
from hexw.registry.prior_use import LAB_METHODS, ROLE_PRIOR_LAB_ANALYSIS


def _is_prior(o: dict) -> bool:
    return (o.get("role") == ROLE_PRIOR_LAB_ANALYSIS or o.get("method") == "prior_analysis_artifact"
            or (o.get("independent") is False and o.get("method") not in LAB_METHODS))


def iter_measurements(registry_dir=None, *, include_prior: bool = False, methods=None) -> Iterator[dict]:
    """Yield observation records (plain dicts). Default: original measurements only (XRD, Mossbauer),
    no prior_lab_analysis artefacts. `methods` optionally restricts the measurement method."""
    reg = Path(registry_dir) if registry_dir is not None else paths.get_config().registry
    doc = yaml.safe_load((reg / "observations.yaml").read_text(encoding="utf-8"))
    for o in doc["observations"]:
        if not include_prior and _is_prior(o):
            continue
        if methods is not None and o.get("method") not in methods:
            continue
        yield o
