"""Curie temperature model T_C(composition) and the reduced magnetisation m(T/T_C) (fast-model approximation).

Data: ONLY literature/fit (L02 Table 3, abstract numbers as an alternative, L01 Mg/Ni thermogravimetry). Everything
that is not in a fit record is an explicit ASSUMPTION (named below, with its justification in docs/ANISOTROPY_MODEL.md).

T_C model (linear dilution of the exchange network, `x_A` = A atoms per formula unit SrFe(18-x)A(x)O27):

    T_C(x) = T0 - sum_A s_A x_A                         [K]
    s_Co   = (T0 - T_Co2)/2 =: free dilution slope of a magnetic Co2+      (assumption: 0 <= s_Co <= 50 K)
    s_Zn   = (T0 - T_Zn2)/2                                                (T_Co2, T_Zn2: line through the L02 series)
    s_Mg, s_Ni = (T0 - T_A2^TG - d_def)/2                                  (L01 TG intercept; d_def = definition offset)
    s_Mn   = f_Mn s_Zn, f_Mn in U(0,1)                                     (no Mn T_C in fit: wide, identifiable = False)

For the Co-Zn series (x_Co + x_Zn = 2) T0 cancels: the model is the straight line through the five L02 points
("interpolation inside the series"). T0 (the host Fe2W T_C, absent from fit) only matters outside x_tot = 2.

T_C DEFINITIONS (never mixed silently): L02 -- maximum of |dM/dT| of an isofield M(T), 50-100 Oe ("dmax");
L01 -- intersection of tangents of a thermogravimetric curve in a permanent magnet ("TG"). The model scale is "dmax";
TG values enter through the offset d_def ~ N(0, 15 K) (ASSUMPTION; the two definitions are not measured on one sample).

Reduced magnetisation: Kuz'min-type form m = [1 - s t^(3/2) - (1-s) t^(5/2)]^(1/3), t = T/T_C, shape s in U(0.3,1.5)
(ASSUMPTION: the fit has no digitised M(T); L02-036 is figure-only/qualitative). Alternative: mean-field Brillouin S=5/2.
This is a single effective magnetisation of the whole ferrimagnet -- the multisublattice thermodynamics (TZ 9.2) is the
next stage; here it is only the fast-model approximation.

Pure computation: the fit records are read once through hexw.params.literature (fit zone only); no Pool, no file writes.
"""
from __future__ import annotations

import functools
import math
import warnings
import zlib
from dataclasses import dataclass, field

import numpy as np

from hexw.registry import schemas as S

# --------------------------------------------------------------------------------------------- composition
KNOWN_A = ("Co", "Zn", "Mn", "Ni", "Mg")
X_TOT_MAX = 2.0
_X_TOL = 1e-9

#: assumptions (named; justification in docs/ANISOTROPY_MODEL.md section 2)
SIGMA_TC_READ_K = 10.0          # reading of a derivative maximum on a 10 K/min heating run
S_CO_MAX_K = 50.0               # Co2+ (3 muB) dilutes at most half as strongly as non-magnetic Zn (100 K / x)
TG_DEF_SIGMA_K = 15.0           # TG-tangent T_C vs dM/dT-maximum T_C
R_FIT = (1.71 + 1.799) / 4.0    # mean (actual / nominal) A content of the L02 Co and Zn samples (L02-019, L02-022: 1.71, 1.799 of 2);
#                                 the T_C line and the K parameters are calibrated on samples with this content, a prediction
#                                 for a sample with content ratio r is rescaled by r / R_FIT (shared r of JointDraws)
KUZMIN_S_RANGE = (0.3, 1.5)
KUZMIN_P = 2.5
KUZMIN_BETA = 1.0 / 3.0

SAMPLE_X_ZN = {"SrCo2Fe16O27": 0.0, "SrCo1.5Zn0.5Fe16O27": 0.5, "SrCoZnFe16O27": 1.0,
               "SrCo0.5Zn1.5Fe16O27": 1.5, "SrZn2Fe16O27": 2.0}
TC_TABLE_IDS = ("L02-023", "L02-025", "L02-027", "L02-029", "L02-031")      # x_Zn = 0, .5, 1, 1.5, 2 (Table 3, dmax)
TC_ABSTRACT_IDS = {"L02-034": 0.0, "L02-035": 1.0}                           # alternative numbers (abstract)
TT_DMAX_IDS = ("L02-024", "L02-026", "L02-028", "L02-030")                   # T_T for x_Zn = 0, .5, 1, 1.5
TT_ONSET_IDS = {"L02-032": 0.0, "L02-033": 1.0}                              # text "around 520 / 360 K"
NO_TRANSITION_ID = "L02-032"                                                  # text: absent in SrZn2Fe16O27
BOUNDARY_ID = "L02-037"                                                       # x_Co = 0.65 at 300 K (assumed_by_authors)
STRUCTURE_IDS = {"L02-045": ("SrCo2Fe16O27", "planar"), "L02-046": ("SrCoZnFe16O27", "planar"),
                 "L02-047": ("SrZn2Fe16O27", "uniaxial")}
TC_TG_IDS = {"Mg": "L01-070", "Ni": "L01-078"}


def normalize_composition(comp) -> dict:
    """{'Co': x, 'Zn': y, ...} -> validated dict (A atoms per f.u. of SrFe(18-x)A(x)O27). `{}` is the unsubstituted host
    SrFe18O27 (x = 0), supported explicitly; zero entries are dropped."""
    if not isinstance(comp, dict):
        raise ValueError("composition must be a dict {A: x_A}, A in " + str(KNOWN_A) + " ({} = SrFe18O27)")
    out = {}
    for k, v in comp.items():
        if k not in KNOWN_A:
            raise ValueError(f"unknown substituent {k!r}; supported {KNOWN_A}")
        v = float(v)
        if not math.isfinite(v) or v < 0:
            raise ValueError(f"x_{k} must be finite and >= 0, got {v!r}")
        if v > 0:
            out[k] = out.get(k, 0.0) + v
    tot = sum(out.values())
    if tot > X_TOT_MAX + _X_TOL:
        raise ValueError(f"x_tot = {tot} > {X_TOT_MAX}: outside the validity range of the model")
    return out


def co_zn(x_zn: float) -> dict:
    """The L02 series SrCo(2-x)Zn(x)Fe16O27 (x_tot = 2)."""
    if not (-_X_TOL <= x_zn <= 2 + _X_TOL):
        raise ValueError("x_Zn must lie in [0, 2]")
    return {"Co": max(2.0 - x_zn, 0.0), "Zn": max(x_zn, 0.0)}


# ------------------------------------------------------------------------------------------------ results
@dataclass
class Dist:
    """Sample-based distribution of a model output. `samples[..., i]` belongs to parameter draw i; Dists coming from
    one model object share the draws, i.e. they are JOINT (K1..K4, T_C, m are mutually consistent per index).
    status: 'ok' | 'not_predicted'. `identifiable`: False when the absolute scale/value is not determined by fit data."""
    name: str
    unit: str
    samples: np.ndarray | None
    status: str = "ok"
    identifiable: bool = True
    basis: str = ""
    notes: str = ""
    extra: dict = field(default_factory=dict)

    def _q(self, q):
        if self.samples is None:
            return np.nan
        s = np.asarray(self.samples, float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)       # all-NaN slice (e.g. no transition) -> NaN
            return np.nanquantile(s, q, axis=-1)

    @property
    def median(self):
        return self._q(0.5)

    @property
    def lo(self):
        return self._q(0.025)

    @property
    def hi(self):
        return self._q(0.975)

    @property
    def interval(self):
        return (self.lo, self.hi)

    def prob(self, pred) -> float:
        if self.samples is None:
            return float("nan")
        return float(np.mean(pred(np.asarray(self.samples))))

    def __repr__(self):
        if self.samples is None:
            return f"Dist({self.name}: {self.status}; {self.basis})"
        return f"Dist({self.name}: {np.asarray(self.median)!r} [{np.asarray(self.lo)!r}, {np.asarray(self.hi)!r}] {self.unit}; {self.status})"


def not_predicted(name, unit, why, basis="") -> Dist:
    return Dist(name, unit, None, status="not_predicted", identifiable=False, basis=basis, notes=why)


# --------------------------------------------------------------------------------------------- fit inputs
@dataclass(frozen=True)
class FitInputs:
    x_zn: np.ndarray                  # 0, .5, 1, 1.5, 2
    tc_dmax: np.ndarray               # K (L02 Table 3)
    tc_sigma: np.ndarray              # K (assumption SIGMA_TC_READ_K (+) |table - abstract|)
    tc_ids: tuple
    tt_dmax: np.ndarray               # K for x_zn 0, .5, 1, 1.5
    tt_onset: dict                    # x_zn -> K
    tt_ids: tuple
    tc_tg: dict                       # 'Mg'/'Ni' -> (K, sd)
    boundary_x_co_300K: float
    structures: dict                  # sample -> 'planar' | 'uniaxial'
    rec_ids: dict


def _claims():
    from hexw.params.literature import load_fit_records
    return {c.rec_id: c for c in load_fit_records()}


@functools.lru_cache(maxsize=1)
def fit_inputs() -> FitInputs:
    cl = _claims()

    def est(rid, name_prefix):
        c = cl[rid]
        if not c.name.startswith(name_prefix):
            raise RuntimeError(f"fit record {rid} is {c.name!r}, expected {name_prefix!r}: records changed")
        if c.estimate is None:
            raise RuntimeError(f"fit record {rid} has no estimate")
        return c

    x, tc = [], []
    for rid in TC_TABLE_IDS:
        c = est(rid, "Curie_temperature_Table3")
        x.append(SAMPLE_X_ZN[c.source.sample_id_in_source])
        tc.append(c.estimate.value)
    x, tc = np.array(x), np.array(tc)
    if not np.allclose(x, [0, 0.5, 1, 1.5, 2]):
        raise RuntimeError("unexpected order of the L02 T_C records")
    alt = {}
    for rid, xz in TC_ABSTRACT_IDS.items():
        c = est(rid, "Curie_temperature_abstract")
        alt[xz] = c.estimate.value
    sig = np.array([math.hypot(SIGMA_TC_READ_K, tc[i] - alt.get(x[i], tc[i])) for i in range(5)])
    tt = []
    for rid, xz in zip(TT_DMAX_IDS, x[:4]):
        c = est(rid, "transition_temperature_Table3_derivative_max")
        if SAMPLE_X_ZN[c.source.sample_id_in_source] != xz:
            raise RuntimeError(f"{rid}: sample mismatch")
        tt.append(c.estimate.value)
    onset = {}
    for rid, xz in TT_ONSET_IDS.items():
        onset[xz] = est(rid, "transition_temperature_onset_text").estimate.value
    tg = {}
    for a, rid in TC_TG_IDS.items():
        c = est(rid, "Curie_temperature")
        tg[a] = (c.estimate.value, c.estimate.uncertainty.std)
    b = est(BOUNDARY_ID, "composition_of_planar_to_uniaxial_boundary_at_RT").estimate.value
    st = {}
    for rid, (smp, typ) in STRUCTURE_IDS.items():
        c = cl[rid]
        if c.name != "magnetic_structure_type" or c.qualitative != typ or c.source.sample_id_in_source != smp:
            raise RuntimeError(f"{rid}: structure record changed")
        st[smp] = typ
    ids = {"tc": TC_TABLE_IDS, "tc_abstract": tuple(TC_ABSTRACT_IDS), "tt": TT_DMAX_IDS, "tt_onset": tuple(TT_ONSET_IDS),
           "tg": dict(TC_TG_IDS), "boundary": BOUNDARY_ID, "structure": tuple(STRUCTURE_IDS)}
    return FitInputs(x, tc, sig, TC_TABLE_IDS, np.array(tt), onset, TT_DMAX_IDS, tg, b, st, ids)


# ------------------------------------------------------------------------------------ reduced magnetisation
def _kuz_g(t, s):
    return s * t ** 1.5 + (1.0 - s) * t ** KUZMIN_P


def reduced_magnetization(t, shape=0.8, form: str = "kuzmin"):
    """m(t), t = T/T_C (array-like, broadcast with `shape`). m(0) = 1, m(t >= 1) = 0; monotone decreasing."""
    t = np.asarray(t, float)
    if form == "kuzmin":
        s = np.asarray(shape, float)
        tc_ = np.clip(t, 0.0, 1.0)
        a = 1.0 - _kuz_g(tc_, s)
        return np.where(a > 0, np.abs(a) ** KUZMIN_BETA, 0.0)
    if form == "brillouin":
        tt, mm = _brillouin_table(2.5)
        return np.interp(np.clip(t, 0.0, 1.0), tt, mm)
    raise ValueError(f"unknown magnetisation form {form!r} ('kuzmin' | 'brillouin')")


def brillouin_m(t, spin: float):
    """Mean-field Brillouin reduced magnetisation of spin `spin` (3/2 for Co2+, 5/2 for Fe3+) at t = T/T_ordering;
    m(0) = 1, m(t >= 1) = 0."""
    tt, mm = _brillouin_table(float(spin))
    return np.interp(np.clip(np.asarray(t, float), 0.0, 1.0), tt, mm)


def inverse_reduced_magnetization(m, shape=0.8, form: str = "kuzmin"):
    """t with m(t) = m (0 < m <= 1), vectorised."""
    m = np.asarray(m, float)
    if form == "kuzmin":
        s = np.asarray(shape, float)
        target = 1.0 - np.clip(m, 0.0, 1.0) ** 3         # g(t) = 1 - m^3 (beta = 1/3)
        target, s = np.broadcast_arrays(target, s)
        lo = np.zeros_like(target)
        hi = np.ones_like(target)
        for _ in range(48):
            mid = 0.5 * (lo + hi)
            up = _kuz_g(mid, s) < target
            lo = np.where(up, mid, lo)
            hi = np.where(up, hi, mid)
        return 0.5 * (lo + hi)
    if form == "brillouin":
        tt, mm = _brillouin_table(2.5)
        return np.interp(np.clip(m, 0.0, 1.0), mm[::-1], tt[::-1])
    raise ValueError(f"unknown magnetisation form {form!r}")


@functools.lru_cache(maxsize=8)
def _brillouin_table(S_spin: float = 2.5):
    """Mean-field Brillouin magnetisation m(t) for spin S (t = T/T_C): m = B_S(3S/(S+1) m / t).

    Parametrised by y = (3S/(S+1)) m / t: m = B_S(y), t = (3S/(S+1)) m / y. Limits: y -> 0 gives m -> 0 and t -> 1
    (the Curie point, appended exactly); y -> infinity gives m -> 1 and t -> 0. Small y is not evaluated through the
    difference of two cotangents (cancellation): the grid starts at y = 0.02, where the formula is accurate to ~1e-12.
    Returns (t ascending, m descending), with m(0) = 1 and m(1) = 0 exactly."""
    y = np.geomspace(2e-2, 400.0, 8001)
    a, b = (2 * S_spin + 1) / (2 * S_spin), 1 / (2 * S_spin)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        m = a / np.tanh(a * y) - b / np.tanh(b * y)
    t = (3 * S_spin / (S_spin + 1)) * m / y
    tt = np.concatenate([[0.0], t[::-1], [1.0]])
    mm = np.concatenate([[1.0], m[::-1], [0.0]])
    keep = np.concatenate([[True], np.diff(tt) > 0])
    tt, mm = tt[keep], mm[keep]
    mm = np.minimum.accumulate(mm)
    return tt, mm


# --------------------------------------------------------------------------------------- T_C parameter draws
def _wls_line(x, y, sig):
    w = 1.0 / sig ** 2
    A = np.stack([np.ones_like(x), x], axis=-1)
    Aw = A * w[..., None]
    cov_inv = np.einsum("...ni,...nj->...ij", Aw, A)
    cov = np.linalg.inv(cov_inv)
    beta = np.einsum("...ij,...j->...i", cov, np.einsum("...ni,...n->...i", Aw, y))
    return beta, cov


@functools.lru_cache(maxsize=1)
def tc_line_fit() -> dict:
    """Weighted straight line T_C(x_Zn) = a + b x_Zn through the five L02 Table-3 points (dmax definition)."""
    f = fit_inputs()
    beta, cov = _wls_line(f.x_zn, f.tc_dmax, f.tc_sigma)
    res = f.tc_dmax - (beta[0] + beta[1] * f.x_zn)
    chi2 = float(np.sum((res / f.tc_sigma) ** 2))
    dof = len(f.x_zn) - 2
    s_res = float(math.sqrt(np.sum(res ** 2) / dof))
    return {"a": float(beta[0]), "b": float(beta[1]), "cov": cov, "chi2": chi2, "dof": dof, "residuals": res,
            "s_res": s_res, "x": f.x_zn}


def sample_tc_params(n: int, rng: np.random.Generator) -> dict:
    """Joint draws of the T_C parameters (independent of the reorientation data). Keys (arrays of length n):
    tc_co2, tc_zn2 [K, dmax definition], s_co, f_mn, d_def, tg_mg, tg_ni, tc_scatter (interpolation scatter [K])."""
    f, ln = fit_inputs(), tc_line_fit()
    infl = max(1.0, ln["chi2"] / ln["dof"])                   # inflate the line covariance if chi2/dof > 1
    L = np.linalg.cholesky(ln["cov"] * infl)
    z = rng.standard_normal((n, 2)) @ L.T
    a, b = ln["a"] + z[:, 0], ln["b"] + z[:, 1]
    return {
        "tc_co2": a, "tc_zn2": a + 2.0 * b,
        "s_co": rng.uniform(0.0, S_CO_MAX_K, n),
        "f_mn": rng.uniform(0.0, 1.0, n),
        "d_def": rng.normal(0.0, TG_DEF_SIGMA_K, n),
        "tg_mg": rng.normal(f.tc_tg["Mg"][0], f.tc_tg["Mg"][1], n),
        "tg_ni": rng.normal(f.tc_tg["Ni"][0], f.tc_tg["Ni"][1], n),
        # interpolation scatter: one INDEPENDENT N(0, s_res) per composition (tc_samples draws it from this sigma and
        # a seed keyed by the composition), not one number shared by all compositions of a draw
        "tc_scatter_sigma": np.full(n, ln["s_res"]),
        "scatter_seed": np.full(n, int(rng.integers(0, 2**31 - 1))),
    }


def dilution_slopes(p: dict) -> dict:
    """s_A [K per unit x_A] and the host T0 from a parameter-draw dict."""
    t0 = p["tc_co2"] + 2.0 * p["s_co"]
    s_zn = 0.5 * (t0 - p["tc_zn2"])
    return {"T0": t0, "Co": p["s_co"], "Zn": s_zn,
            "Mg": 0.5 * (t0 - (p["tg_mg"] - p["d_def"])),      # TG value minus offset -> dmax scale
            "Ni": 0.5 * (t0 - (p["tg_ni"] - p["d_def"])),
            "Mn": p["f_mn"] * s_zn}


def composition_key(comp: dict) -> str:
    return ";".join(f"{a}={comp[a]:.6g}" for a in sorted(comp) if comp[a] != 0)      # zero entries do not change the composition


def tc_scatter(comp: dict, p: dict) -> np.ndarray:
    """Interpolation scatter [K] of a composition: independent between compositions, the same for every call with the
    same composition and draws (reproducible)."""
    n = len(p["tc_scatter_sigma"])
    seed = (int(p["scatter_seed"][0]), zlib.crc32(composition_key(comp).encode()))
    return np.random.default_rng(seed).standard_normal(n) * p["tc_scatter_sigma"]


def tc_samples(comp: dict, p: dict) -> np.ndarray:
    """T_C [K] (dmax scale) for each draw in `p`. If `p` carries 'content_scale' (= r / R_FIT of the shared content ratio)
    the A contents are rescaled by it."""
    sl = dilution_slopes(p)
    cs = p.get("content_scale", 1.0)
    out = sl["T0"].copy()
    for a, x in comp.items():
        out = out - x * cs * sl[a]
    return out + tc_scatter(comp, p)


def tc_status(comp: dict) -> tuple[str, bool, str]:
    """(status, identifiable, basis) -- what the fit data actually say about this composition."""
    tot = sum(comp.values())
    keys = set(comp)
    if keys <= {"Co", "Zn"} and abs(tot - 2.0) < 1e-6:
        return ("interpolated_in_series", True,
                "straight line through L02 Table 3 (x_Zn = 0..2), 'interpolation inside the series'")
    if "Mn" in keys:
        return ("wide_prior", False, "no Mn T_C in fit: dilution slope f_Mn*s_Zn, f_Mn in U(0,1) (assumption)")
    if keys <= {"Mg", "Ni"} and len(keys) == 1 and abs(tot - 2.0) < 1e-6:
        return ("measured_other_definition", True, "L01 TG-intercept T_C shifted by the definition offset d_def")
    if tot == 0:
        return ("wide_prior", False, "host T0 not in fit: T0 = T_Co2 + 2 s_Co, s_Co in U(0,50 K) (assumption)")
    return ("extrapolated_from_endpoints", False,
            "x_tot != 2 or Mg/Ni mixtures: linear dilution from the host T0 whose value is NOT in fit (assumption)")


# --------------------------------------------------------------------------------------------- public API
@functools.lru_cache(maxsize=4)
def _tc_draws(n: int, seed: int) -> dict:
    return sample_tc_params(n, np.random.default_rng(seed))


DEFAULT_N = 4000
DEFAULT_SEED = 20260101


def tc(composition, *, n: int = DEFAULT_N, seed: int = DEFAULT_SEED) -> Dist:
    """T_C [K] (definition: maximum of |dM/dT| at 50-100 Oe, L02) as a sample distribution.
    Status: interpolated_in_series | measured_other_definition | extrapolated_from_endpoints | wide_prior."""
    comp = normalize_composition(composition)
    st, ident, basis = tc_status(comp)
    s = np.maximum(tc_samples(comp, _tc_draws(n, seed)), 0.0)
    return Dist("T_C", "K", s, status=st, identifiable=ident, basis=basis,
                notes="definition: dmax (L02); L01 Mg/Ni are TG-intercept values shifted by d_def")


def ms_reduced(T, composition, *, n: int = DEFAULT_N, seed: int = DEFAULT_SEED, draws=None) -> Dist:
    """Reduced magnetisation m(T) = M(T)/M(0) (fast-model approximation, single effective magnetisation).
    Uses the joint draws of hexw.params.anisotropy (shape parameter conditioned on the same fit) unless `draws` is
    given. T scalar or array [K] -> samples shape (len(T), n)."""
    from hexw.params import anisotropy as _an
    comp = normalize_composition(composition)
    d = draws if draws is not None else _an.joint_draws(n, seed)
    tcs = np.maximum(tc_samples(comp, d.tc_params), 1.0)
    Tarr = np.atleast_1d(np.asarray(T, float))
    if np.any(Tarr < 0) or not np.all(np.isfinite(Tarr)):
        raise ValueError("T must be finite and >= 0 K")
    m = reduced_magnetization(Tarr[:, None] / tcs[None, :], d.shape_s, d.form)
    st, ident, basis = tc_status(comp)
    return Dist("m_reduced", "1", m[0] if np.ndim(T) == 0 else m, status=st, identifiable=ident,
                basis=basis + "; shape from Kuz'min-type form (assumption), single effective sublattice")


# ------------------------------------------------------------------------------------ parameter estimates
ASSUMPTION_ID = "ASSUMPTION"


def _prov(ids, note, status):
    return S.Provenance(source_id=";".join(ids), locator="literature/fit/records_fit.yaml", note=note, status=status)


def _ps(name, value, unit, kind, ids, note, std=None, ident=True, uniform=None, status=None):
    """ParameterEstimate. An item whose only source is a modelling assumption gets the role `assumption` (source_id
    'ASSUMPTION', status `assumption`); sourced items are literature_fit with status `verified` (a number read from a
    record) or `derived_by_us` (computed here from record numbers). Mixed items are not built: split them."""
    ids = list(ids)
    is_assumption = ids == [ASSUMPTION_ID]
    if ASSUMPTION_ID in ids and not is_assumption:
        raise ValueError(f"{name}: mixed ASSUMPTION / sourced estimate; make two estimates")
    if is_assumption:
        role, st = S.SourceRole.ASSUMPTION, S.ProvenanceStatus.ASSUMPTION
    else:
        role, st = S.SourceRole.LITERATURE_FIT, S.ProvenanceStatus(status or "derived_by_us")
    if uniform is not None:
        u = S.Uncertainty(kind="distribution", distribution="uniform", params={"low": uniform[0], "high": uniform[1]})
    elif std is not None:
        u = S.Uncertainty(kind="std", std=float(std))
    else:
        u = S.Uncertainty(kind="unknown")
    return S.ParameterEstimate(name=name, value=float(value), unit=unit, quantity_kind=kind, source_role=role,
                               provenance=_prov(ids, note, st), uncertainty=u, identifiable=ident)


def curie_parameter_estimates(n: int = 20000, seed: int = DEFAULT_SEED) -> list:
    """T_C / m(T) parameters as ParameterEstimate. Sourced items: role literature_fit (status derived_by_us for numbers
    computed from the L02 line, verified for numbers read from a record); pure assumptions: role assumption."""
    from hexw.units import Kind
    f, ln = fit_inputs(), tc_line_fit()
    p = sample_tc_params(n, np.random.default_rng(seed))
    ids_l02 = list(f.tc_ids) + list(TC_ABSTRACT_IDS)
    A_ = ASSUMPTION_ID
    out = [
        _ps("Tc_Co2_dmax", np.median(p["tc_co2"]), "K", Kind.TEMPERATURE, ids_l02,
            "line through L02 Table 3 at x_Zn=0; abstract numbers (780/680) enter as sigma inflation", std=np.std(p["tc_co2"])),
        _ps("Tc_Zn2_dmax", np.median(p["tc_zn2"]), "K", Kind.TEMPERATURE, ids_l02,
            "line through L02 Table 3 at x_Zn=2 (dmax definition)", std=np.std(p["tc_zn2"])),
        _ps("Tc_interpolation_scatter", ln["s_res"], "K", Kind.TEMPERATURE, list(f.tc_ids),
            "rms residual of the 5-point line (3 dof), added to every predicted T_C (independent per composition)",
            std=ln["s_res"] * 0.5),
        _ps("Tc_Mg2_TG", f.tc_tg["Mg"][0], "K", Kind.TEMPERATURE, [TC_TG_IDS["Mg"]],
            "L01 thermogravimetry, tangent intersection (different definition from L02)", std=f.tc_tg["Mg"][1],
            status="verified"),
        _ps("Tc_Ni2_TG", f.tc_tg["Ni"][0], "K", Kind.TEMPERATURE, [TC_TG_IDS["Ni"]],
            "L01 thermogravimetry, tangent intersection (different definition from L02)", std=f.tc_tg["Ni"][1],
            status="verified"),
        _ps("A_content_ratio_fit_samples", R_FIT, "1", Kind.DIMENSIONLESS, ["L02-019", "L02-022"],
            "mean refined/nominal A content of the L02 SrCo2W and SrZn2W samples (1.71, 1.799 of 2): the T_C line and "
            "the K parameters are calibrated on samples of this content", std=0.02),
        _ps("A_content_ratio_prediction", 0.925, "1", Kind.DIMENSIONLESS, [A_],
            "prior r ~ U(0.85, 1) (value = its mean): content ratio of a predicted sample (shared by Ms0, T_C and K in "
            "JointDraws)", ident=False),
        _ps("Tc_def_offset_TG_minus_dmax", 0.0, "K", Kind.TEMPERATURE, [A_],
            f"N(0, {TG_DEF_SIGMA_K} K): the two T_C definitions are never measured on one sample", ident=False),
        _ps("Tc_dilution_slope_Co", S_CO_MAX_K / 2, "K", Kind.TEMPERATURE, [A_],
            "U(0,50 K) per unit x_Co; fixes the host T0 = T_Co2 + 2 s_Co (T0 is not in fit)", ident=False,
            uniform=(0.0, S_CO_MAX_K)),
        _ps("Tc_Mn_dilution_fraction", 0.5, "1", Kind.DIMENSIONLESS, [A_],
            "s_Mn = f_Mn s_Zn, f_Mn in U(0,1): no Mn T_C in fit", ident=False, uniform=(0.0, 1.0)),
        _ps("m_shape_s", 0.9, "1", Kind.DIMENSIONLESS, [A_],
            f"Kuz'min-type shape parameter prior U{KUZMIN_S_RANGE}: M(T) of L02 is figure-only (L02-036, not digitised)",
            ident=False, uniform=KUZMIN_S_RANGE),
    ]
    return out
