"""Simple chemical validators: site overflow, charge neutrality, Fe2+/Fe3+ model of the ТЗ (section 7).

Model: SrFe(18-x)A(x)O27 with Sr2+, A2+, O27 and only Fe2+/Fe3+:
    n(Fe2+) = 2 - x,  n(Fe3+) = 16.
Not applicable for x > 2 (n(Fe2+) < 0), nor to mixed valence / oxygen defects.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from hexw.registry.schemas import StructureCandidate

TOL = 1e-9
RESERVED_NAMES = frozenset({"Sr", "O", "Fe", "Fe2+", "Fe3+"})


class ChemistryError(ValueError):
    pass


class SiteOverflowError(ChemistryError):
    pass


class ChargeError(ChemistryError):
    pass


class ModelNotApplicableError(ChemistryError):
    pass


def check_site_occupancy(candidate: StructureCandidate, tol: float = TOL) -> None:
    for s in candidate.sites:
        for sp, occ in s.species.items():
            if not math.isfinite(occ):
                raise SiteOverflowError(f"site {s.site}: non-finite occupancy of {sp}: {occ!r}")
            if not occ >= -tol:
                raise SiteOverflowError(f"site {s.site}: negative occupancy of {sp}")
        total = sum(s.species.values())
        if not total <= 1 + tol:          # written so that NaN is rejected
            raise SiteOverflowError(f"site {s.site}: total occupancy {total:.6g} > 1")


def formula_charge(counts: dict[str, float], oxidation: dict[str, float]) -> float:
    missing = [el for el in counts if el not in oxidation]
    if missing:
        raise ChargeError(f"no oxidation state given for {missing}")
    for el, n in counts.items():
        if not math.isfinite(n) or not math.isfinite(oxidation[el]):
            raise ChargeError(f"non-finite count or oxidation state for {el}: {n!r}, {oxidation[el]!r}")
    return sum(n * oxidation[el] for el, n in counts.items())


def check_charge_neutral(counts: dict[str, float], oxidation: dict[str, float], tol: float = 1e-6) -> None:
    q = formula_charge(counts, oxidation)
    if not abs(q) <= tol:                 # written so that NaN is rejected
        raise ChargeError(f"formula charge {q:+.6g} != 0")


@dataclass(frozen=True)
class FeValenceModel:
    x: float
    A: str
    n_fe2: float
    n_fe3: float
    n_fe: float
    counts: dict
    oxidation: dict


def fe2_fe3_model(x: float, A: str = "A") -> FeValenceModel:
    """Conditional model; raises ModelNotApplicableError for x outside [0, 2] (or non-finite) and
    ChemistryError if A is one of the reserved species names."""
    if not isinstance(A, str) or not A or A in RESERVED_NAMES:
        raise ChemistryError(f"substituent name {A!r} is not allowed: it must be a non-empty string different from "
                             f"the reserved species {sorted(RESERVED_NAMES)} (they would silently merge in `counts`)")
    if not math.isfinite(x):
        raise ModelNotApplicableError(f"x={x!r} is not a finite number")
    if x < -TOL:
        raise ModelNotApplicableError(f"x={x}: negative substituent content")
    if x > 2 + TOL:
        raise ModelNotApplicableError(
            f"x={x} > 2: n(Fe2+) = 2 - x < 0. The Fe2+/Fe3+ model with Sr2+, A2+, O27 is not applicable; "
            "mixed valence (e.g. Fe4+, A3+, Mn3+) or oxygen/cation defects need a different model.")
    n2, n3 = 2 - x, 16.0
    counts = {"Sr": 1.0, A: x, "Fe2+": n2, "Fe3+": n3, "O": 27.0}
    oxidation = {"Sr": 2.0, A: 2.0, "Fe2+": 2.0, "Fe3+": 3.0, "O": -2.0}
    check_charge_neutral(counts, oxidation)
    return FeValenceModel(x, A, n2, n3, n2 + n3, counts, oxidation)
