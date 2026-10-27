"""Literature base for the theory side: ONLY the fit part (literature/fit/records_fit.yaml).

`load_fit_records()` is callable from hexw.theory/params: the fit files are outside every protected zone. The holdout
records (literature_holdout/) and the holdout publications are blocked by the data guard (hexw.io.paths, zone HOLDOUT);
this module never touches them. The excluded / methodology_only records (literature/excluded) are not loaded here:
they cannot be used in a parametrization.
"""
from __future__ import annotations

from pathlib import Path

import yaml

from hexw.io import paths
from hexw.registry import schemas as S

FIT_FILE = Path("fit") / "records_fit.yaml"


class FitRecordsError(RuntimeError):
    pass


def fit_records_path() -> Path:
    return paths.get_config().literature / FIT_FILE


def load_fit_records() -> list[S.LiteratureClaim]:
    """All fit claims as validated `LiteratureClaim` objects (role == literature_fit, split == fit)."""
    p = fit_records_path()
    if not p.is_file():
        raise FitRecordsError(f"{p} not found: run `python -m hexw.ingest.literature_split` first")
    with open(p, encoding="utf-8") as f:
        doc = yaml.safe_load(f)
    out = [S.from_dict(S.LiteratureClaim, d) for d in doc["claims"]]
    bad = [c.rec_id for c in out if c.role is not S.SourceRole.LITERATURE_FIT or c.split != "fit"]
    if bad:
        raise FitRecordsError(f"records with split != fit in the fit file: {bad[:5]}")
    return out


def load_fit_estimates() -> list[S.ParameterEstimate]:
    """Numeric fit parameters only (claims that carry a ParameterEstimate)."""
    return [c.estimate for c in load_fit_records() if c.estimate is not None]
