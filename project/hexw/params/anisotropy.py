"""Magnetocrystalline anisotropy K1..K4 (x_A, T) of SrFe(18-x)A(x)O27 from fit literature + physics. See docs/ANISOTROPY_MODEL.md.

Conventions of the magnet core (hexw.theory.magnet.core): F = K1 s^2 + K2 s^4 + K3 s^6 + K4 s^6 cos6phi, s = sin(theta),
theta from the c axis, J/m^3.  K1 > 0: easy axis, K1 < 0: easy plane.

MODEL (host: Callen-Callen l = 2 law in the reduced magnetisation m(T) of hexw.params.curie; Co2+: single-ion term):

    K1(x, T) = [K_Fe + (x_Co + x_Zn)(A - K_Fe)/2 + x_Mn (A_Mn - K_Fe)/2] m^3  -  rho A x_Co F_Co(T)  (+ Mn terms)
    K2(x, T) = q2 Kscale m^10   (l = 4)      K3 = 0 (assumed)       K4(x, T) = q4 Kscale m^21   (l = 6)
    Kscale   = |bracket| + rho A x_Co / 2    (full single-ion scale, not only the host)

  A = K1 of the (axial) SrZn2W host at T = 0, K_Fe = K1 of SrFe18O27 at T = 0, A_Mn = Mn2 end member (magnitude);
  m^3 = l = 2 single-ion law of the host (Fe3+ S-state, dipolar + weak SOC);
  F_Co(T), two STRUCTURAL ALTERNATIVES (co_form, never averaged; both are returned):
      'ion_mf'  (default)  F = m_ion(T / Theta_Co)^3, m_ion = mean-field Brillouin S = 3/2, Theta_Co = f T_C, f ~ U(0.15, 1):
                           for S = 3/2 the single-ion operators have rank l = 2 only, so K_Co ~ <O_2> ~ m_ion^3 (l(l+1)/2 = 3);
                           the Co2+ sublattice follows ITS OWN molecular field, hence the second parameter f (rho, f)
      'cc_power'           F = m^p, p ~ U(3, 21): an EFFECTIVE exponent, not a consequence of the Callen-Callen theorem
                           (p is identifiable = False; the combination of (rho, p) that the T_T data fix is)
  rho    = |Co2+ contribution| / A per Co at T = 0, log-uniform (0.5, 50), TRUNCATED by an explicit ASSUMPTION: the single-ion
           energy of Co2+ at 0 K, rho A V_fu, does not exceed E_CO_MAX = 20 meV (a general scale of spin-orbit energy of a
           3d ion in a distorted oxygen environment, not a literature number)
  q2, q4 = dimensionless prior-only nuisances (K2, K4 are not identifiable).

For the Co-Zn series (x_Co + x_Zn = 2) the bracket is A, so K1 = A [m^3 - rho x_Co F]: A is only an overall scale, the sign
and the zero crossing T_T (the T at which the planar state is left on heating) depend on (rho, F-parameter, m(T), q2) only.

IDENTIFICATION (L02): the planar -> non-planar transition on heating is the point where the planar state loses stability,
dE/d(s) at s = 1:  K1 + 2 K2 = 0 (the cone phase between K1+2K2 = 0 and K1 = 0 is not resolved by the experiment, so the
criterion K1 + 2 K2 = 0 is the first exit from the plane and is what the table T_T measures). Constraints used:
4 T_T (table, dmax; sigma = read-off (+) |table - text onset| definition spread), no transition in SrZn2W, planar at 300 K
for x_Zn = 0 and 1 (L02-045/046), axial at 300 K for x_Zn = 2 (L02-047).  The posterior is a self-normalised importance
sample from the prior (SIR); a chi^2 profile on the (rho, p | f) grid is provided for control.
Optional branch l11_constraint = 'co_zn_ratio': the ratio |K1(Co2)|/|K1(Zn2)| at 300 K of L11 (mu0 Ha Ms/2, sign-blind,
valid only for H >> Ha, which is NOT the case: hence sigma_ln = 0.7) multiplies the likelihood. Base forecast: 'none'.

VALIDITY IN T: K(x, T) at T < 250 K is `extrapolated_in_T` (the fit constrains only T_T >= 265 K and 300 K states) and the
interval is widened (log-normal factor, sigma growing linearly to 0.7 at 0 K).

Mn: the sign of the Mn contribution follows the valence hypothesis of stage 4A (shared draws): Mn2+ (S ion) -- axial with weight
0.85, magnitude as the host (L11 Mn); Mn3+ (Jahn-Teller) -- single-ion term of the SAME form as Co, sign 50/50 (all weights ASSUMPTIONS).

Everything not in a fit record is a named ASSUMPTION (constants below, justified in the document).  Pure computation:
no Pool, no file writes; fit records are read once.
"""
from __future__ import annotations

import functools
import math
from dataclasses import dataclass, field

import numpy as np

from hexw.params import curie as C
from hexw.params.curie import Dist, not_predicted, normalize_composition
from hexw.registry import schemas as S

CO_FORMS = ("ion_mf", "cc_power")
L11_CONSTRAINTS = ("none", "co_zn_ratio")
HOSTS = ("dft", "l11")
DEFAULT_CO_FORM = "ion_mf"
DEFAULT_L11 = "none"
DEFAULT_HOST = "dft"

# ------------------------------------------------------------------------------------------------ assumptions
EXP_HOST = 3.0                   # Callen-Callen l = 2 : K ~ m^(l(l+1)/2) = m^3
EXP_K2 = 10.0                    # l = 4
EXP_K4 = 21.0                    # l = 6
P_RANGE = (3.0, 21.0)            # cc_power: EFFECTIVE Co2+ exponent prior (l = 2 .. 6); not a theorem
F_RANGE = (0.15, 1.0)            # ion_mf: Theta_Co / T_C
CO_SPIN = 1.5                    # Co2+ high spin S = 3/2 (3 muB): single-ion operators of rank l = 2 only
RHO_RANGE = (0.5, 50.0)          # log-uniform prior on rho (ratio of Co2+ to host anisotropy scale at 0 K), before the energy cap
E_CO_MAX_MEV = 20.0              # ASSUMPTION: ceiling of the Co2+ single-ion energy rho*A*V_fu at 0 K (spin-orbit energy scale)
MEV_J = 1.602176634e-22          # 1 meV in J
Q2_SIGMA = 0.3                   # K2 = q2 Kscale m^10, q2 ~ N(0, 0.3)  (|K2| <~ |K1| scale; sign unknown)
Q4_SIGMA = 0.03                  # K4 = q4 Kscale m^21, q4 ~ N(0, 0.03) (basal anisotropy of hexaferrites << K1; no fit data)
SIGMA_TT_READ_K = 10.0           # read-off of a derivative maximum (same as Tc)
SIGMA_LN_DFT = 0.5               # 1-sigma log-normal width of the GGA+U+SOC Ku (3 methods agree only in trend, L05-103)
SIGMA_LN_ZN_RATIO = 0.5          # K1(Zn2)/K1(Fe18) at 0 K: no fit number, bracketed by DFT (1) and L11 (see doc)
SIGMA_LN_L11 = 0.8               # tier 3 (Ba -> Sr, sign-blind law of approach to saturation, formula not valid at H ~ Ha)
SIGMA_LN_L11_MN = 1.0            # Mn: additionally inflated (Mn2+ moment disputed in L11-016; Ms of the L11 Mn sample anomalously low)
SIGMA_LN_BA_SR = 0.4             # ASSUMPTION: width of ln(Ba->Sr factor): exponent q of K(y) ~ y^q in [0.5, 2] gives 0.20..0.41
SIGMA_LN_L11_RATIO = 0.7         # co_zn_ratio constraint: |K1(Co2)|/|K1(Zn2)| (L11, 300 K); formula inapplicable at H ~ Ha
SIGMA_DENSITY_REL = 0.05         # density used to convert sigma -> Ms for L11 (cell mass of the Ba compound / DFT Sr cell)
MN2_AXIAL_WEIGHT = 0.85          # ASSUMPTION: P(Mn2+ contribution is axial, as the host); S ion, small SOC
MN3_SIGN_WEIGHT = 0.5            # ASSUMPTION: P(Mn3+ single-ion term is axial)
LOW_T_K = 250.0                  # K(x, T) below this T is extrapolated_in_T
SIGMA_LN_LOWT = 0.7              # ASSUMPTION: extra log-normal width of K at 0 K (linear ramp from LOW_T_K): a factor ~2 (1 sigma)
N_AVOGADRO = 6.02214076e23
ATOMIC_MASS = {"Sr": 87.62, "Ba": 137.327, "Fe": 55.845, "Zn": 65.38, "Co": 58.933, "Mn": 54.938,
               "Ni": 58.6934, "Mg": 24.305, "O": 15.999}

L05_KU_ID, L05_VOL_ID, L05_MS_ID = "L05-103", "L05-007", "L05-094"
L05_KU_SRBA_ID = "L05-112"       # Sr0.5Ba0.5Fe18O27, the same DFT method: the in-fit sensitivity to Ba
L11_IDS = {"Fe": ("L11-001", "L11-003"), "Zn": ("L11-010", "L11-012"), "Mn": ("L11-013", "L11-015"),
           "Co": ("L11-004", "L11-006"), "Ni": ("L11-007", "L11-009")}  # (Ms rec, Ha rec)
SEEDS = {"posterior": 20260202}
DEFAULT_N_PRIOR = 200_000
_GRID_T = np.linspace(0.0, 0.999, 150)


# ------------------------------------------------------------------------------------------- fit-derived inputs
@dataclass(frozen=True)
class HostInputs:
    ku_dft: float              # J/m^3 (SrFe18O27, 0 K)
    ku_srba: float             # J/m^3 (Sr0.5Ba0.5Fe18O27, 0 K, same method)
    vol_dft_m3: float          # m^3 per unit cell (2 f.u.)
    ms_dft: float              # A/m
    l11: dict                  # Me -> (sigma [A m^2/kg], mu0 Ha [T])
    ids: dict


@functools.lru_cache(maxsize=1)
def host_inputs() -> HostInputs:
    from hexw.params.literature import load_fit_records
    cl = {c.rec_id: c for c in load_fit_records()}

    def val(rid, prefix, unit):
        c = cl[rid]
        if not c.name.startswith(prefix) or c.estimate is None or c.estimate.unit != unit:
            raise RuntimeError(f"fit record {rid} changed: {c.name!r} / {c.estimate}")
        return c.estimate.value

    ku = val(L05_KU_ID, "Ku_uniaxial_MAE_constant", "MJ/m^3") * 1e6
    ku_ba = val(L05_KU_SRBA_ID, "Ku_uniaxial_MAE_constant", "MJ/m^3") * 1e6
    vol = val(L05_VOL_ID, "cell_volume_DFT_relaxed", "angstrom^3") * 1e-30
    ms = val(L05_MS_ID, "mu0Ms_DFT", "T") / (4e-7 * math.pi)
    if cl[L05_KU_ID].source.sample_id_in_source != "SrW":
        raise RuntimeError("L05-103 is not the SrFe18O27 (SrW) entry")
    if cl[L05_KU_SRBA_ID].source.formula_in_source != "Sr0.5Ba0.5Fe18O27":
        raise RuntimeError("L05-112 is not the Sr0.5Ba0.5Fe18O27 entry")
    l11 = {}
    for me, (rs, rh) in L11_IDS.items():
        sig = val(rs, "Ms_RT", "A*m^2/kg")
        ha = val(rh, "Ha_anisotropy_field_RT", "Oe") * 1e-4          # the record is mu0*Ha in Oe -> tesla (1 Oe <-> 1e-4 T)
        l11[me] = (sig, ha)
        if f"Me={me}" not in cl[rs].source.sample_id_in_source:
            raise RuntimeError(f"{rs}: expected Me={me}")
    for sp,rid in (("Mg","L01-088"),("Ni","L01-089")):
        c=cl[rid]
        if c.name!="moment_angle_from_c_axis" or c.estimate.value!=0 or c.source.formula_in_source!=f"Sr{sp}2Fe16O27":
            raise RuntimeError(f"axial end-member constraint changed: {rid}")
    ids = {"ku": L05_KU_ID, "ku_srba": L05_KU_SRBA_ID, "vol": L05_VOL_ID, "ms": L05_MS_ID, "l11": dict(L11_IDS)}
    return HostInputs(ku, ku_ba, vol, ms, l11, ids)


def ba_to_sr_factor(hin: HostInputs | None = None) -> float:
    """Multiplier K(Sr compound)/K(Ba compound) from the IN-FIT sensitivity of the L05 DFT Ku to Ba (same method, same
    functional): K(y) = K_Sr + (K_Ba - K_Sr) y, y = Ba fraction, with K(0) = L05-103 (SrW), K(0.5) = L05-112 (Sr0.5Ba0.5W):
    K_Ba = 2 K(0.5) - K(0), factor = K(0)/K_Ba = 0.326. The linearity in y is an ASSUMPTION covered by SIGMA_LN_BA_SR."""
    hin = hin or host_inputs()
    kba = 2.0 * hin.ku_srba - hin.ku_dft
    if kba <= 0:
        raise RuntimeError("L05 Ba sensitivity is not monotonic: the transfer factor is undefined")
    return hin.ku_dft / kba


def ba_compound_density(me: str, x: float, vol_m3: float) -> float:
    """kg/m^3 of BaFe(18-x)Me(x)O27 using the DFT Sr cell volume (assumption: the Ba cell is ~1-2 % larger, absorbed in
    SIGMA_DENSITY_REL).  Z = 2 f.u. per cell."""
    m = ATOMIC_MASS
    mass = m["Ba"] + (18 - x) * m["Fe"] + x * m[me] + 27 * m["O"] if me != "Fe" else m["Ba"] + 18 * m["Fe"] + 27 * m["O"]
    return 2.0 * mass * 1e-3 / N_AVOGADRO / vol_m3


def l11_k_300(me: str) -> float:
    """|K1| [J/m^3] at 300 K of the BARIUM compound from L11 (K = mu0 Ha Ms / 2, uniaxial, K2 = 0): central value, density from
    the cell. The law of approach to saturation gives only |K|, and is valid for H >> Ha, which L11 does not satisfy."""
    h = host_inputs()
    sig, ha = h.l11[me]
    x = 0.0 if me == "Fe" else 2.0
    rho = ba_compound_density(me, x, h.vol_dft_m3)
    return 0.5 * ha * sig * rho


def l11_co_zn_ratio() -> float:
    """|K1(BaCo2W)| / |K1(BaZn2W)| at 300 K from L11 (tier 3, Ba; sign-blind): the target of the 'co_zn_ratio' branch."""
    return l11_k_300("Co") / l11_k_300("Zn")


def v_fu_m3() -> float:
    return host_inputs().vol_dft_m3 / 2.0


def rho_cap(a_host) -> np.ndarray:
    """Largest rho with rho |A| V_fu <= E_CO_MAX (the explicit ASSUMPTION on the Co2+ single-ion energy at 0 K)."""
    return E_CO_MAX_MEV * MEV_J / (np.abs(np.asarray(a_host, float)) * v_fu_m3())


# -------------------------------------------------------------------------------------- the Co2+ temperature function
def co_F(co_form: str, m, t, p, f):
    """F_Co(T), the temperature factor of the Co2+ single-ion term (F(0) = 1). m: host reduced magnetisation, t = T/T_C of the
    sample; p (cc_power) / f (ion_mf) broadcast against m."""
    m = np.asarray(m, float)
    if co_form == "cc_power":
        with np.errstate(divide="ignore", invalid="ignore"):
            lm = np.log(np.where(m > 0, m, 1.0))
        return np.where(m > 0, np.exp(np.asarray(p, float) * lm), 0.0)
    if co_form == "ion_mf":
        return C.brillouin_m(np.asarray(t, float) / np.asarray(f, float), CO_SPIN) ** 3
    raise ValueError(f"co_form must be one of {CO_FORMS}")


def _check_branches(host, co_form, l11):
    if host not in HOSTS:
        raise ValueError(f"host must be one of {HOSTS}")
    if co_form not in CO_FORMS:
        raise ValueError(f"co_form must be one of {CO_FORMS}")
    if l11 not in L11_CONSTRAINTS:
        raise ValueError(f"l11_constraint must be one of {L11_CONSTRAINTS}")


# ------------------------------------------------------------------------------------------------ the T_T fit
def _tt_sigma(fi: "C.FitInputs") -> np.ndarray:
    d = {0.0: abs(fi.tt_dmax[0] - fi.tt_onset[0.0]), 1.0: abs(fi.tt_dmax[2] - fi.tt_onset[1.0])}
    typical = float(np.mean(list(d.values())))
    spread = np.array([d.get(float(x), typical) for x in fi.x_zn[:4]])
    return np.hypot(SIGMA_TT_READ_K, spread)


@dataclass(frozen=True)
class TTData:
    x_zn: np.ndarray          # compositions with an observed T_T
    tc: np.ndarray            # K (table values: the measured T_C of the same sample)
    tt: np.ndarray            # K
    sigma: np.ndarray         # K
    planar_300: tuple         # x_Zn with planar structure at 300 K
    axial_300: tuple          # x_Zn uniaxial at 300 K
    tc_zn2: float             # K (Zn2 end member: used for the 'no transition' and 'axial' constraints)


def tt_data() -> TTData:
    fi = C.fit_inputs()
    planar = tuple(C.SAMPLE_X_ZN[s] for s, t in fi.structures.items() if t == "planar")
    axial = tuple(C.SAMPLE_X_ZN[s] for s, t in fi.structures.items() if t == "uniaxial")
    return TTData(fi.x_zn[:4].copy(), fi.tc_dmax[:4].copy(), fi.tt_dmax.copy(), _tt_sigma(fi), planar, axial,
                  float(fi.tc_dmax[4]))


def _k12(t, rho, p, f, s, q2, x_co, form="kuzmin", co_form=DEFAULT_CO_FORM):
    """(K1, K2) / A of the Co-Zn series at reduced temperature t; arrays broadcast (draws on axis 0)."""
    m = C.reduced_magnetization(t, s, form)
    m = np.maximum(m, 1e-12)
    lm = np.log(m)
    big_f = co_F(co_form, m, t, p, f)
    scale2 = 1.0 + 0.5 * rho * x_co                          # Kscale / A of the series
    return np.exp(EXP_HOST * lm) - rho * x_co * big_f, q2 * scale2 * np.exp(EXP_K2 * lm)


def _g(t, rho, p, f, s, q2, x_co, form="kuzmin", co_form=DEFAULT_CO_FORM):
    """(K1 + 2 K2) / A of the Co-Zn series at reduced temperature t; arrays broadcast (draws on axis 0)."""
    k1, k2 = _k12(t, rho, p, f, s, q2, x_co, form, co_form)
    return k1 + 2.0 * k2


def first_exit_t(rho, p, f, s, q2, x_co, form="kuzmin", co_form=DEFAULT_CO_FORM):
    """t = T_T/T_C at which the planar state (K1 + 2K2 < 0 at low T) is left on heating; NaN if there is no planar state
    at 0 K or the planar state survives up to T_C.  rho, p, f, s, q2: arrays (n,)."""
    rho, p, f, s, q2 = (np.asarray(a, float)[:, None] for a in (rho, p, f, s, q2))
    t = _GRID_T[None, :]
    k1, k2 = _k12(t, rho, p, f, s, q2, x_co, form, co_form)
    return _crossing(k1 + 2.0 * k2, _GRID_T, planar0=(easy_state(k1[:, 0], k2[:, 0]) == 1))


def _crossing(g, grid, planar0=None):
    """First upward zero of g[i, :] on `grid` (linear interpolation); NaN if g[i, 0] >= 0 or no crossing. `planar0`: boolean (n,)
    -- the plane must also be the GLOBAL minimum at the first grid point (a metastable plane of an axial ground state is not
    a state that is 'left on heating')."""
    planar0 = (g[:, 0] < 0) if planar0 is None else ((g[:, 0] < 0) & planar0)
    up = g >= 0
    j = np.argmax(up, axis=1)
    has = up.any(axis=1) & planar0 & (j > 0)
    jj = np.clip(j, 1, len(grid) - 1)
    g0 = g[np.arange(len(g)), jj - 1]
    g1 = g[np.arange(len(g)), jj]
    w = np.divide(-g0,g1-g0,out=np.zeros_like(g0),where=g1!=g0)
    tcross = grid[jj - 1] + w * (grid[jj] - grid[jj - 1])
    return np.where(has, tcross, np.nan)


def tt_log_likelihood(rho, p, f, s, q2, data: TTData, form="kuzmin", co_form=DEFAULT_CO_FORM, l11_constraint=DEFAULT_L11):
    """log L (Gaussian T_T + hard structure constraints [+ L11 Co/Zn ratio]) for draws; -inf where a constraint is violated."""
    n = len(rho)
    chi2 = np.zeros(n)
    ok = np.ones(n, bool)
    for xz, tcv, tt, sg in zip(data.x_zn, data.tc, data.tt, data.sigma):
        tcross = first_exit_t(rho, p, f, s, q2, 2.0 - xz, form, co_form)
        bad = np.isnan(tcross)
        ok &= ~bad
        chi2 += np.where(bad, 0.0, ((tcross * tcv - tt) / sg) ** 2)
    # structures at 300 K (x_Co = 2 - x_Zn); Zn2 end member: no planar state at any T
    # the observed structure must be the GLOBAL minimum (easy_state), not merely a local one (plane stable, axis lower)
    for xz in data.planar_300:
        tcv = float(np.interp(xz, data.x_zn, data.tc)) if xz in data.x_zn else None
        if tcv is not None:
            k1, k2 = _k12(300.0 / tcv, rho, p, f, s, q2, 2.0 - xz, form, co_form)
            ok &= (easy_state(k1, k2) == 1) & (k1 + 2.0 * k2 < 0)
    for xz in data.axial_300:
        if abs(xz - 2.0) < 1e-9:
            col = lambda a: np.asarray(a)[:, None]    # noqa: E731
            k1, k2 = _k12(_GRID_T[None, :], col(rho), col(p), col(f), col(s), col(q2), 0.0, form, co_form)
            ok &= np.all(easy_state(k1, k2) == 0, axis=1)
        elif xz in data.x_zn:
            tcv = float(data.tc[list(data.x_zn).index(xz)])
            k1, k2 = _k12(300.0 / tcv, rho, p, f, s, q2, 2.0 - xz, form, co_form)
            ok &= easy_state(k1, k2) == 0
    ll = -0.5 * chi2
    if l11_constraint == "co_zn_ratio":
        ll = ll - 0.5 * (np.log(k1_ratio_co_zn_300(rho, p, f, s, data, form, co_form) / l11_co_zn_ratio())
                         / SIGMA_LN_L11_RATIO) ** 2
    return np.where(ok, ll, -np.inf)


def k1_ratio_co_zn_300(rho, p, f, s, data: TTData, form="kuzmin", co_form=DEFAULT_CO_FORM):
    """|K1(Co2)| / |K1(Zn2)| at 300 K of the model (A cancels), measured T_C of the samples (L02 Table 3)."""
    tco, tzn = float(data.tc[0]), float(data.tc_zn2)
    mc = np.maximum(C.reduced_magnetization(300.0 / tco, s, form), 1e-12)
    mz = np.maximum(C.reduced_magnetization(300.0 / tzn, s, form), 1e-12)
    fc = co_F(co_form, mc, 300.0 / tco, p, f)
    k_co = np.abs(mc ** EXP_HOST - 2.0 * rho * fc)
    return np.maximum(k_co, 1e-9) / mz ** EXP_HOST


@dataclass
class TTPosterior:
    rho: np.ndarray
    p: np.ndarray
    f: np.ndarray
    s: np.ndarray
    q2: np.ndarray
    weights: np.ndarray        # normalised importance weights of the prior draws
    ess: float
    n_prior: int
    chi2_min: float
    best: dict
    co_form: str = DEFAULT_CO_FORM
    l11_constraint: str = DEFAULT_L11
    evidence: float = float("nan")      # mean likelihood over the prior (comparison of the structural forms)
    frac_feasible: float = float("nan")


def fit_tt(data: TTData, *, n_prior: int = DEFAULT_N_PRIOR, seed: int = SEEDS["posterior"], chunk: int = 20000,
           form: str = "kuzmin", co_form: str = DEFAULT_CO_FORM, l11_constraint: str = DEFAULT_L11) -> TTPosterior:
    """Importance sampling of (rho, p, f, s, q2) from their priors weighted by the T_T likelihood (4 numeric points +
    structure constraints [+ the L11 ratio]): 2 fit parameters (rho, p | f) + 2 nuisances kept at their priors (s, q2).
    The parameter that the chosen co_form does not use (f for cc_power, p for ion_mf) stays at its prior."""
    _check_branches("dft", co_form, l11_constraint)
    rng = np.random.default_rng(seed)
    R, Pp, Ff, Ss, Q = [], [], [], [], []
    L = []
    for i in range(0, n_prior, chunk):
        m = min(chunk, n_prior - i)
        rho = np.exp(rng.uniform(math.log(RHO_RANGE[0]), math.log(RHO_RANGE[1]), m))
        p = rng.uniform(*P_RANGE, m)
        f = rng.uniform(*F_RANGE, m)
        s = rng.uniform(*C.KUZMIN_S_RANGE, m)
        q2 = rng.normal(0.0, Q2_SIGMA, m)
        L.append(tt_log_likelihood(rho, p, f, s, q2, data, form, co_form, l11_constraint))
        R.append(rho), Pp.append(p), Ff.append(f), Ss.append(s), Q.append(q2)
    rho, p, f, s, q2, ll = map(np.concatenate, (R, Pp, Ff, Ss, Q, L))
    fin = np.isfinite(ll)
    if not fin.any():
        raise RuntimeError("no prior draw is compatible with the T_T data: the model is falsified by the fit")
    w = np.exp(ll - np.max(ll[fin]))
    w[~fin] = 0.0
    evidence = float(np.sum(np.exp(ll[fin])) / len(ll))
    w /= w.sum()
    ess = 1.0 / float(np.sum(w ** 2))
    k = int(np.argmax(ll))
    best = {"rho": float(rho[k]), "p": float(p[k]), "f": float(f[k]), "s": float(s[k]), "q2": float(q2[k])}
    return TTPosterior(rho, p, f, s, q2, w, ess, n_prior, float(-2 * ll[k]), best, co_form, l11_constraint, evidence,
                       float(fin.mean()))


@functools.lru_cache(maxsize=8)
def _posterior(n_prior: int = DEFAULT_N_PRIOR, seed: int = SEEDS["posterior"], co_form: str = DEFAULT_CO_FORM,
               l11_constraint: str = DEFAULT_L11) -> TTPosterior:
    return fit_tt(tt_data(), n_prior=n_prior, seed=seed, co_form=co_form, l11_constraint=l11_constraint)


def chi2_profile(data: TTData | None = None, s: float = 0.9, q2: float = 0.0, ng: int = 121,
                 co_form: str = "cc_power", second: float | None = None):
    """chi^2 (T_T only) on a (log rho, p) [cc_power] or (log rho, f) [ion_mf] grid at fixed nuisances; for the other
    parameter of the form `second` is not used. Returns (rho_grid, par_grid, chi2[ng, ng], feasible)."""
    data = data or tt_data()
    rg = np.exp(np.linspace(math.log(RHO_RANGE[0]), math.log(RHO_RANGE[1]), ng))
    pg = np.linspace(*(P_RANGE if co_form == "cc_power" else F_RANGE), ng)
    RR, PP = np.meshgrid(rg, pg, indexing="ij")
    one = np.full(RR.size, 1.0)
    p, f = (PP.ravel(), one * 0.5) if co_form == "cc_power" else (one * 13.0, PP.ravel())
    ll = tt_log_likelihood(RR.ravel(), p, f, np.full(RR.size, s), np.full(RR.size, q2), data, "kuzmin", co_form)
    chi2 = (-2 * ll).reshape(ng, ng)
    return rg, pg, chi2, np.isfinite(chi2)


# ------------------------------------------------------------------------------------------------ joint draws
@dataclass
class JointDraws:
    """All parameter draws, index i = one consistent parameter set (shared by T_C, m, K1..K4 and by hexw.params.moments.Ms0):
    the content ratio r, the Mn valence branch, T_C parameters, shape s, (rho, p, f, q2, q4), host scales. The branch
    (host, co_form, l11_constraint) is fixed per object; random numbers of the quantities that do not depend on the branch
    (T_C, r, Mn valence, low-T factor) are COMMON to all branches of the same (n, seed)."""
    n: int
    host: str
    co_form: str
    l11_constraint: str
    tc_params: dict            # includes 'r' and 'content_scale' = r / R_FIT
    r: np.ndarray              # content ratio actual/nominal ~ U(0.85, 1): the same r for Ms0, T_C and K
    mn_valence: np.ndarray     # 2 | 3 per draw (hypothesis of 4A, shared by every quantity)
    shape_s: np.ndarray
    form: str
    rho: np.ndarray
    p: np.ndarray
    f: np.ndarray
    q2: np.ndarray
    q4: np.ndarray
    mn: dict                   # sign2, sign3, rho, p, f (Mn3+ single-ion term / Mn2+ sign)
    hosts: dict                # {host: {'K_fe','A','A_mn'}} (J/m^3 at 0 K); only the branch host is present
    low_t_z: np.ndarray        # standard normal of the low-T widening (one number per draw)
    posterior: TTPosterior = field(repr=False, default=None)
    n_cap_resampled: int = 0   # draws resampled because of the E_CO_MAX truncation


def _lognormal(rng, sigma, n):
    return np.exp(rng.normal(0.0, sigma, n))


def _hosts(rng, n, tcp, s, hin: HostInputs) -> dict:
    sl = C.dilution_slopes(tcp)
    t0, s_zn = sl["T0"], sl["Zn"]
    # --- alternative 1: DFT L05 (SrFe18O27, 0 K); Zn2 and Mn2 ends by ratios (assumption, wide)
    k_fe = hin.ku_dft * _lognormal(rng, SIGMA_LN_DFT, n)
    a_dft = k_fe * _lognormal(rng, SIGMA_LN_ZN_RATIO, n)
    ratio_mn = (hin.l11["Mn"][0] * hin.l11["Mn"][1]) / (hin.l11["Zn"][0] * hin.l11["Zn"][1]) * \
        ba_compound_density("Mn", 2, hin.vol_dft_m3) / ba_compound_density("Zn", 2, hin.vol_dft_m3)
    a_mn_dft = a_dft * ratio_mn * _lognormal(rng, SIGMA_LN_L11_MN, n)
    # --- alternative 2: L11 (Ba, 300 K) -> 0 K by the host m^3 and Ba -> Sr by the in-fit L05 sensitivity
    dens = _lognormal(rng, SIGMA_DENSITY_REL, n)
    f_ba = ba_to_sr_factor(hin) * _lognormal(rng, SIGMA_LN_BA_SR, n)          # one transfer factor per draw

    def k0(me, tc_me, sig_ln):
        sg, ha = hin.l11[me]
        x = 0.0 if me == "Fe" else 2.0
        k300 = 0.5 * ha * sg * ba_compound_density(me, x, hin.vol_dft_m3) * dens * _lognormal(rng, sig_ln, n)
        m300 = C.reduced_magnetization(300.0 / np.maximum(tc_me, 301.0), s)
        return f_ba * k300 / np.maximum(m300, 1e-6) ** EXP_HOST

    tc_fe = t0
    tc_zn2 = tcp["tc_zn2"]
    tc_mn2 = t0 - 2.0 * sl["Mn"]
    l11 = {"K_fe": k0("Fe", tc_fe, SIGMA_LN_L11), "A": k0("Zn", tc_zn2, SIGMA_LN_L11),
           "A_mn": k0("Mn", tc_mn2, SIGMA_LN_L11_MN)}
    # Ni magnitude is available in FIT L11 (Ba analogue); the sign constraint is
    # independently available in FIT L01-089: axial SrNi2W at 300 K.
    tc_ni = C.tc_samples({"Ni": 2.0}, tcp)
    a_ni = k0("Ni", tc_ni, SIGMA_LN_L11)
    # Mg: no usable anisotropy magnitude in the current fit corpus. This is an
    # explicitly uncalibrated host-transfer SCENARIO (factor 10.5 at 95% in each
    # direction, before host uncertainty); never a learned Mg anisotropy.
    a_mg = a_dft * _lognormal(rng, 1.2, n)
    l11.update(A_ni=a_ni, A_mg=a_mg)
    return {"dft": {"K_fe": k_fe, "A": a_dft, "A_mn": a_mn_dft, "A_ni": a_ni, "A_mg": a_mg}, "l11": l11}


@functools.lru_cache(maxsize=16)
def joint_draws(n: int = C.DEFAULT_N, seed: int = C.DEFAULT_SEED, n_prior: int = DEFAULT_N_PRIOR, *,
                host: str = DEFAULT_HOST, co_form: str = DEFAULT_CO_FORM, l11_constraint: str = DEFAULT_L11) -> JointDraws:
    """Draws for T_C, m and K1..K4 (and for Ms0): (rho, p, f, s, q2) from the T_T posterior of the branch (resampled; rho
    truncated by the Co2+ energy ceiling of the host), T_C line, r, Mn branch and host/nuisance priors independent of them.
    Cached; the arrays are shared, never mutate them."""
    _check_branches(host, co_form, l11_constraint)
    from hexw.structures.ensemble import MN_VALENCE_PRIOR
    hin = host_inputs()
    post = _posterior(n_prior, SEEDS["posterior"], co_form, l11_constraint)
    tcp = C.sample_tc_params(n, np.random.default_rng(seed))
    rc = np.random.default_rng([seed, 11])                   # branch-independent quantities
    r = rc.uniform(0.85, 1.0, n)
    vals = sorted(MN_VALENCE_PRIOR)
    w = np.array([MN_VALENCE_PRIOR[v] for v in vals], float)
    mnv = rc.choice(vals, size=n, p=w / w.sum())
    sign2 = np.where(rc.random(n) < MN2_AXIAL_WEIGHT, 1.0, -1.0)
    sign3 = np.where(rc.random(n) < MN3_SIGN_WEIGHT, 1.0, -1.0)
    u_rho = rc.random(n)
    p_mn = rc.uniform(*P_RANGE, n)
    f_mn = rc.uniform(*F_RANGE, n)
    z_lowt = rc.standard_normal(n)
    q4 = rc.normal(0.0, Q4_SIGMA, n)
    tcp["r"] = r
    tcp["content_scale"] = r / C.R_FIT
    # the host scales (L11 branch) depend on the shape s, the cap on the host scale A: iterate until no draw violates the
    # ceiling (the host random numbers are the same in every iteration)
    rp = np.random.default_rng([seed, 13])
    idx = rp.choice(len(post.weights), size=n, replace=True, p=post.weights)
    n_res = 0
    for _ in range(60):
        hosts_all = _hosts(np.random.default_rng([seed, 12]), n, tcp, post.s[idx], hin)
        cap = rho_cap(hosts_all[host]["A"])
        bad = post.rho[idx] > cap
        if not bad.any():
            break
        n_res += int(bad.sum())
        idx[bad] = rp.choice(len(post.weights), size=int(bad.sum()), replace=True, p=post.weights)
    hs = hosts_all[host]
    rho = post.rho[idx].copy()
    bad = rho > cap
    if bad.mean() > 0.01:
        raise RuntimeError(f"the Co2+ energy ceiling ({E_CO_MAX_MEV} meV) excludes {bad.mean():.1%} of the posterior")
    rho[bad] = cap[bad]
    s = post.s[idx]
    # Mn3+ single-ion amplitude: log-uniform prior truncated by the same ceiling (analytic)
    lo, hi = math.log(RHO_RANGE[0]), np.log(np.minimum(RHO_RANGE[1], np.maximum(cap, RHO_RANGE[0] * 1.0000001)))
    rho_mn = np.exp(lo + u_rho * (hi - lo))
    mn = {"sign2": sign2, "sign3": sign3, "rho": rho_mn, "p": p_mn, "f": f_mn}
    return JointDraws(n, host, co_form, l11_constraint, tcp, r, mnv, s, "kuzmin", rho, post.p[idx], post.f[idx],
                      post.q2[idx], q4, mn, {host: hs}, z_lowt, post, n_res)


# --------------------------------------------------------------------------------------------------- the model
def anisotropy_status(comp: dict) -> tuple[str, bool, str]:
    """(status, scale_identifiable, basis); scale = the absolute K1 magnitude."""
    keys = set(comp)
    tot = sum(comp.values())
    if "Mg" in keys:
        return ("uncalibrated_scenario", False,
                "Mg anisotropy: uncalibrated host-transfer magnitude; axial sign supported by FIT L01-088 at 300 K; axial branch/K2 truncation away from that endpoint is an assumption; NOT a material-specific magnitude prediction")
    if "Ni" in keys:
        return ("transfer_tier3", False,
                "Ni magnitude from FIT L11-007/009 (BaNi2W, 300 K); broad Ba->Sr and temperature transfer; axial sign from FIT L01-089")
    if keys <= {"Co", "Zn"} and abs(tot - 2.0) < 1e-6 and keys:
        return ("series_interpolation", False,
                "L02 Co-Zn series: sign and T_T fitted (identifiable); absolute scale A from host priors (L05 DFT | L11), not identifiable")
    if "Mn" in keys:
        return ("transfer_tier3", False, "Mn only from L11 BaMn2W Ha (tier 3, Ba->Sr, sign-blind law of approach to saturation); "
                "sign from the valence hypothesis (assumption)")
    return ("extrapolated", False, "outside the Co-Zn series (x_tot != 2): linear host interpolation + Co2+ term are assumptions")


def _k_arrays(comp: dict, Tmat, d: JointDraws, widen: bool = True):
    """K1, K2, K3, K4, m, T_C for temperatures `Tmat` (shape (nT, 1) or (nT, n))."""
    host = d.host
    tcs = np.maximum(C.tc_samples(comp, d.tc_params), 1.0)
    t = Tmat / tcs[None, :]
    m = np.maximum(C.reduced_magnetization(t, d.shape_s[None, :], d.form), 0.0)
    h = d.hosts[host]
    cs = d.tc_params["content_scale"]
    x_co, x_zn, x_mn = (comp.get(k, 0.0) * cs for k in ("Co", "Zn", "Mn"))
    v3 = d.mn_valence == 3
    a_mn = np.where(v3, h["A"], d.mn["sign2"] * np.abs(h["A_mn"]))
    kh = h["K_fe"] + (x_co + x_zn) * (h["A"] - h["K_fe"]) / 2.0 + x_mn * (a_mn - h["K_fe"]) / 2.0
    for sp, endpoint in (("Ni", "A_ni"), ("Mg", "A_mg")):
        kh = kh + comp.get(sp, 0.0) * cs * (h[endpoint] - h["K_fe"]) / 2.0
    if set(comp)<={"Mg","Ni"} and comp:
        # The uncertain nominal/actual-content ratio cannot make an x<=2
        # divalent solid solution overshoot its endpoint and acquire a negative
        # host fraction. Preserve convex mixing for this axial transfer branch.
        fraction=np.minimum(sum(comp.values())*cs/2.,1.)
        end=sum(comp.get(sp,0.)*h[key] for sp,key in (("Mg","A_mg"),("Ni","A_ni")))/sum(comp.values())
        kh=(1-fraction)*h["K_fe"]+fraction*end
    with np.errstate(divide="ignore", invalid="ignore"):
        lm = np.log(np.where(m > 0, m, 1.0))
    mp = lambda e: np.where(m > 0, np.exp(e * lm), 0.0)       # noqa: E731
    big_f = co_F(d.co_form, m, t, d.p[None, :], d.f[None, :])
    k1 = kh[None, :] * mp(EXP_HOST) - (d.rho * h["A"] * x_co)[None, :] * big_f
    scale = np.abs(kh) + d.rho * np.abs(h["A"]) * x_co / 2.0
    if np.any(x_mn > 0):
        f_mn = co_F(d.co_form, m, t, d.mn["p"][None, :], d.mn["f"][None, :])
        k1 = k1 + np.where(v3, d.mn["sign3"] * d.mn["rho"] * h["A"] * x_mn, 0.0)[None, :] * f_mn
        scale = scale + np.where(v3, d.mn["rho"] * np.abs(h["A"]) * x_mn / 2.0, 0.0)
    k2 = (d.q2 * scale)[None, :] * mp(EXP_K2)
    k4 = (d.q4 * scale)[None, :] * mp(EXP_K4)
    if set(comp)<={"Mg","Ni"} and comp:
        # L01 Mg/Ni end-member axial magnetic order constrains the sign at 300 K.
        # Keep prior nuisance magnitudes, but enforce an axial host branch rather
        # than introducing a plane solely through the Co-series q2 nuisance.
        k2=np.maximum(k2,-.49*np.maximum(k1,0))
    if widen:
        sig = SIGMA_LN_LOWT * np.clip(1.0 - np.asarray(Tmat, float) / LOW_T_K, 0.0, 1.0)
        fac = np.exp(sig * d.low_t_z[None, :])
        k1, k2, k4 = k1 * fac, k2 * fac, k4 * fac
    return k1, k2, np.zeros_like(k1), k4, m, tcs


def easy_state(k1, k2):
    """Global minimiser of E(s) = K1 s + K2 s^2, s = sin^2(theta) in [0,1]: 0 = easy axis, 1 = easy plane, 2 = cone."""
    k1, k2 = np.broadcast_arrays(np.asarray(k1, float), np.asarray(k2, float))
    e0, e1 = np.zeros_like(k1), k1 + k2
    with np.errstate(divide="ignore", invalid="ignore"):
        sstar = np.where(k2 != 0, -k1 / (2 * k2), np.nan)
    interior = (k2 > 0) & (sstar > 0) & (sstar < 1)
    ei = np.where(interior, k1 * sstar + k2 * sstar ** 2, np.inf)
    best = np.argmin(np.stack([e0, e1, ei]), axis=0)
    return best


STATE_NAMES = {0: "axis", 1: "plane", 2: "cone"}


@dataclass
class AnisotropyResult:
    K1: Dist
    K2: Dist
    K3: Dist
    K4: Dist
    T: np.ndarray
    composition: dict
    host: str
    status: str
    p_state: dict              # state name -> probability per temperature (array)
    notes: str = ""
    co_form: str = DEFAULT_CO_FORM
    l11_constraint: str = DEFAULT_L11

    def __iter__(self):
        return iter((self.K1, self.K2, self.K3, self.K4))

    def most_likely_state(self):
        names = list(self.p_state)
        arr = np.stack([self.p_state[k] for k in names])
        return [names[i] for i in np.argmax(arr, axis=0)]


def _resolve_branch(host, co_form, l11_constraint, draws):
    if draws is not None:
        for nm, v, dv in (("host", host, draws.host), ("co_form", co_form, draws.co_form),
                          ("l11_constraint", l11_constraint, draws.l11_constraint)):
            if v is not None and v != dv:
                raise ValueError(f"{nm}={v!r} contradicts the supplied draws ({dv!r})")
        return draws.host, draws.co_form, draws.l11_constraint
    host = DEFAULT_HOST if host is None else host
    co_form = DEFAULT_CO_FORM if co_form is None else co_form
    l11_constraint = DEFAULT_L11 if l11_constraint is None else l11_constraint
    _check_branches(host, co_form, l11_constraint)
    return host, co_form, l11_constraint


def anisotropy(composition, T, *, host: str | None = None, co_form: str | None = None, l11_constraint: str | None = None,
               n: int = C.DEFAULT_N, seed: int = C.DEFAULT_SEED, draws: JointDraws | None = None) -> AnisotropyResult:
    """(K1, K2, K3, K4) [J/m^3, core convention] with sample distributions for T (K, scalar or array -> samples (nT, n)).
    Branches, ALTERNATIVES that are never averaged (use `anisotropy_branches` for all of them):
      host: 'dft' (L05 SrFe18O27 0 K prior, default) | 'l11' (L11 Ba->Sr, tier 3);
      co_form: 'ion_mf' (default) | 'cc_power' (structural form of the Co2+ term);
      l11_constraint: 'none' (default) | 'co_zn_ratio' (L11 ratio |K1(Co2)|/|K1(Zn2)| enters the likelihood).
    T < 250 K: status 'extrapolated_in_T' and a widened interval. Mg is an uncalibrated scenario; Ni is tier-3 transfer; the absolute scale is
    never identifiable (flag on each Dist; sign/T_T: see `extra`)."""
    host, co_form, l11_constraint = _resolve_branch(host, co_form, l11_constraint, draws)
    comp = normalize_composition(composition)
    Tarr = np.atleast_1d(np.asarray(T, float))
    if np.any(Tarr < 0) or not np.all(np.isfinite(Tarr)):
        raise ValueError("T must be finite and >= 0 K")
    status, _, basis = anisotropy_status(comp)
    if status == "not_predicted":
        nm = [not_predicted(f"K{i}", "J/m^3", basis) for i in (1, 2, 3, 4)]
        return AnisotropyResult(*nm, Tarr, comp, host, status, {}, notes=basis, co_form=co_form,
                                l11_constraint=l11_constraint)
    d = draws if draws is not None else joint_draws(n, seed, host=host, co_form=co_form, l11_constraint=l11_constraint)
    k1, k2, k3, k4, m, tcs = _k_arrays(comp, Tarr[:, None], d)
    k1n, k2n, _, _, _, _ = _k_arrays(comp, Tarr[:, None], d, widen=False)
    scalar = np.ndim(T) == 0
    sq = (lambda a: a[0]) if scalar else (lambda a: a)
    low_t = Tarr < LOW_T_K
    ident_sign = status == "series_interpolation"
    st_out = "extrapolated_in_T" if low_t.any() else status
    if low_t.any():
        basis = basis + f"; T < {LOW_T_K:g} K is outside the fitted range (T_T >= 265 K, 300 K states): extrapolated, interval widened"
    ex = {"scale_identifiable": False, "sign_and_T_T_identifiable": ident_sign and not low_t.all(), "host": host,
          "co_form": co_form, "l11_constraint": l11_constraint, "composition_status": status,
          "T_extrapolated": low_t.copy() if not scalar else bool(low_t[0]),
          "posterior_ess": d.posterior.ess if d.posterior is not None else None}
    mk = lambda name, a, ident, nt: Dist(name, "J/m^3", sq(a), status=st_out, identifiable=ident, basis=basis,  # noqa: E731
                                         notes=nt, extra=dict(ex))
    out = (mk("K1", k1, False, "K1 = A[m^3 - rho x_Co F_Co]; abs. scale A not identifiable, sign/T_T from T_T fit"),
           mk("K2", k2, False, "K2 = q2 Kscale m^10, q2 ~ N(0, 0.3): prior only; Mg/Ni-only branch truncated at K2>=-0.49 K1 to retain the L01 axial host hypothesis"),
           Dist("K3", "J/m^3", sq(k3), status="assumed_zero", identifiable=False, basis=basis,
                notes="no fit information on the sin^6 term: set to 0 (assumption)", extra=dict(ex)),
           mk("K4", k4, False, "basal-plane constant (literature 'K3 sin^6 cos6phi'): prior only, N(0,0.03 Kscale m^21); CRITICAL for planar FMR"))
    st = easy_state(k1n, k2n)
    p_state = {nm: np.mean(st == i, axis=1) for i, nm in STATE_NAMES.items()}
    return AnisotropyResult(*out, Tarr, comp, host, st_out, p_state, notes=basis, co_form=co_form,
                            l11_constraint=l11_constraint)


def anisotropy_branches(composition, T, *, n: int = C.DEFAULT_N, seed: int = C.DEFAULT_SEED) -> dict:
    """All 8 alternatives {(host, co_form, l11_constraint): AnisotropyResult}; nothing is averaged. The forecast of record is
    (dft, ion_mf, none)."""
    return {(h, cf, lc): anisotropy(composition, T, host=h, co_form=cf, l11_constraint=lc, n=n, seed=seed)
            for h in HOSTS for cf in CO_FORMS for lc in L11_CONSTRAINTS}


def reorientation_temperature(composition, *, host: str | None = None, co_form: str | None = None,
                              l11_constraint: str | None = None, n: int = C.DEFAULT_N, seed: int = C.DEFAULT_SEED,
                              draws: JointDraws | None = None) -> Dist:
    """T at which the planar state is left on heating (K1 + 2K2 = 0 first crossing). Samples are NaN where the
    model has no planar low-T state (extra['p_none'])."""
    host, co_form, l11_constraint = _resolve_branch(host, co_form, l11_constraint, draws)
    comp = normalize_composition(composition)
    status, _, basis = anisotropy_status(comp)
    if status == "not_predicted":
        return not_predicted("T_T", "K", basis)
    d = draws if draws is not None else joint_draws(n, seed, host=host, co_form=co_form, l11_constraint=l11_constraint)
    tcs = np.maximum(C.tc_samples(comp, d.tc_params), 1.0)
    Tmat = _GRID_T[:, None] * tcs[None, :]
    k1, k2, _, _, _, _ = _k_arrays(comp, Tmat, d, widen=False)
    tt = _crossing((k1 + 2.0 * k2).T, _GRID_T, planar0=(easy_state(k1[0], k2[0]) == 1)) * tcs
    return Dist("T_T", "K", tt, status=status, identifiable=(status == "series_interpolation"), basis=basis,
                notes="first exit from the easy-plane state on heating (K1+2K2=0)",
                extra={"p_none": float(np.mean(np.isnan(tt))), "host": host, "co_form": co_form,
                       "l11_constraint": l11_constraint})


# ------------------------------------------------------------------------------------------------ residuals
def fit_residuals(draws: JointDraws | None = None, *, co_form: str | None = None, l11_constraint: str | None = None) -> dict:
    """Posterior predictive check on the fit data: T_T (4 points), structures at 300 K, authors' boundary (not fitted)."""
    d = draws if draws is not None else joint_draws(host=DEFAULT_HOST, co_form=co_form or DEFAULT_CO_FORM,
                                                    l11_constraint=l11_constraint or DEFAULT_L11)
    data = tt_data()
    rows = []
    for xz, tcv, tt, sg in zip(data.x_zn, data.tc, data.tt, data.sigma):
        t = first_exit_t(d.rho, d.p, d.f, d.shape_s, d.q2, 2.0 - xz, d.form, d.co_form) * tcv
        lo, med, hi = np.nanquantile(t, [0.025, 0.5, 0.975])
        rows.append({"x_Zn": float(xz), "T_T_obs": float(tt), "sigma_obs": float(sg), "pred_median": float(med),
                     "pred_lo": float(lo), "pred_hi": float(hi), "pull": float((med - tt) / sg),
                     "inside_95": bool(lo <= tt <= hi)})
    # structure at 300 K from the model (compositions of L02 samples)
    st = {}
    for xz in (0.0, 1.0, 2.0):
        r = anisotropy(C.co_zn(xz), 300.0, draws=d)
        st[xz] = max(r.p_state, key=lambda k: float(r.p_state[k][0])) if r.p_state else None
        st[(xz, "p")] = {k: float(v[0]) for k, v in r.p_state.items()}
    # boundary at 300 K: x_Co at which P(plane) = 1/2
    xs = np.linspace(0.0, 2.0, 81)
    pp = []
    for xc in xs:
        r = anisotropy(C.co_zn(2.0 - xc), 300.0, draws=d)
        pp.append(float(r.p_state["plane"][0]))
    pp = np.array(pp)
    k = np.where(pp >= 0.5)[0]
    boundary = float(xs[k[0]]) if len(k) else float("nan")
    return {"T_T": rows, "structure_300K": st, "boundary_x_Co_300K": boundary,
            "boundary_authors": C.fit_inputs().boundary_x_co_300K}


def l11_ratio_tension(co_form: str = DEFAULT_CO_FORM, n: int = C.DEFAULT_N, seed: int = C.DEFAULT_SEED) -> dict:
    """The posterior |K1(Co2)|/|K1(Zn2)| at 300 K of the model WITHOUT the L11 constraint versus the L11 number: a RESULT
    (tension), not hidden. Returns the quantiles of the model ratio, the L11 value and P(model ratio <= L11 ratio)."""
    d = joint_draws(n, seed, host=DEFAULT_HOST, co_form=co_form, l11_constraint="none")
    r = anisotropy(C.co_zn(0.0), 300.0, draws=d).K1.samples
    z = anisotropy(C.co_zn(2.0), 300.0, draws=d).K1.samples
    ratio = np.abs(r) / np.abs(z)
    tgt = l11_co_zn_ratio()
    return {"quantiles": dict(zip(("p2.5", "p16", "p50", "p84", "p97.5"),
                                  (float(v) for v in np.quantile(ratio, [.025, .16, .5, .84, .975])))),
            "L11_ratio": tgt, "P_below_L11": float(np.mean(ratio <= tgt)),
            "z_ln": float((np.log(np.median(ratio)) - math.log(tgt)) / SIGMA_LN_L11_RATIO)}


# ------------------------------------------------------------------------------------ parameter estimates
def _ids(*a):
    return ";".join(a)


def _prior_samples(n: int, seed: int) -> dict:
    rng = np.random.default_rng([seed, 99])
    return {"rho": np.exp(rng.uniform(math.log(RHO_RANGE[0]), math.log(RHO_RANGE[1]), n)),
            "p": rng.uniform(*P_RANGE, n), "f": rng.uniform(*F_RANGE, n), "s": rng.uniform(*C.KUZMIN_S_RANGE, n),
            "q2": rng.normal(0.0, Q2_SIGMA, n), "q4": rng.normal(0.0, Q4_SIGMA, n)}


NARROWING_MIN = 2.0       # an estimate is identifiable only if the 95 % posterior width is < 1/2 of the prior width


def _w95(a):
    lo, hi = np.quantile(a, [0.025, 0.975])
    return float(hi - lo)


def anisotropy_parameter_estimates(n: int = 20000, seed: int = C.DEFAULT_SEED) -> list:
    """K-model parameters as ParameterEstimate. role literature_fit (status verified for a number read from a record,
    derived_by_us for a value computed here) for sourced items; role assumption (source_id 'ASSUMPTION') for pure assumptions.
    Identifiable = the 95 % posterior width is narrower than half the prior width (data must narrow the prior at least 2x);
    a non-identifiable estimate reports the PRIOR interval."""
    from hexw.units import Kind
    hin = host_inputs()
    f = C.fit_inputs()
    A_ = "ASSUMPTION"
    tt_ids = list(f.tt_ids) + list(C.TT_ONSET_IDS) + list(C.TC_TABLE_IDS) + list(C.STRUCTURE_IDS) + [C.NO_TRANSITION_ID]
    pri = _prior_samples(n, seed)
    ed = Kind.ENERGY_DENSITY

    def est(name, post, prior, unit, kind, ids, note, status="derived_by_us"):
        post = np.asarray(post, float)
        prior = np.asarray(prior, float)
        med = float(np.median(post))
        narrowing = _w95(prior) / max(_w95(post), 1e-300)
        assumption = list(ids) == [A_]
        ident = bool(narrowing >= NARROWING_MIN) and not assumption
        src = prior if not ident else post
        lo, hi = (float(x) for x in np.quantile(src, [0.025, 0.975]))
        note = note + f"; 95% width posterior/prior = {1.0 / narrowing:.2f} (identifiable requires <= {1.0 / NARROWING_MIN:.2f})"
        if not ident:
            note += "; the interval reported is the PRIOR interval"
        u = S.Uncertainty(kind="interval", lower=lo, upper=hi, declared_wide=not ident)
        role = S.SourceRole.ASSUMPTION if assumption else S.SourceRole.LITERATURE_FIT
        st = S.ProvenanceStatus.ASSUMPTION if assumption else S.ProvenanceStatus(status)
        return S.ParameterEstimate(name=name, value=med, unit=unit, quantity_kind=kind, source_role=role,
                                   provenance=S.Provenance(source_id=_ids(*ids), locator="literature/fit/records_fit.yaml",
                                                           note=note, status=st),
                                   uncertainty=u, identifiable=ident)

    def const(name, value, unit, kind, ids, note, status="derived_by_us"):
        """A single number (no sample): an assumption (role assumption) or a value computed from records."""
        assumption = list(ids) == [A_]
        role = S.SourceRole.ASSUMPTION if assumption else S.SourceRole.LITERATURE_FIT
        st = S.ProvenanceStatus.ASSUMPTION if assumption else S.ProvenanceStatus(status)
        return S.ParameterEstimate(name=name, value=float(value), unit=unit, quantity_kind=kind, source_role=role,
                                   provenance=S.Provenance(source_id=_ids(*ids), locator="literature/fit/records_fit.yaml",
                                                           note=note, status=st),
                                   uncertainty=S.Uncertainty(kind="unknown"), identifiable=not assumption)

    out = []
    # the host scales (DFT branch)
    dd = joint_draws(n, seed, host="dft")
    dl = joint_draws(n, seed, host="l11")
    out.append(S.ParameterEstimate(
        "aniso_K1_SrFe18O27_0K_L05", hin.ku_dft, "J/m^3", ed, S.SourceRole.LITERATURE_FIT,
        S.Provenance(source_id=hin.ids["ku"], locator="literature/fit/records_fit.yaml", status=S.ProvenanceStatus.VERIFIED,
                     note="L05 GGA+U+SOC Ku (0 K, effective K1+2K2), the number of the record; the log-normal width is a separate "
                          "ASSUMPTION; the excluded secondary experimental citation is NOT used"),
        S.Uncertainty(kind="unknown"), identifiable=True))
    out.append(const("aniso_sigma_ln_DFT_Ku", SIGMA_LN_DFT, "1", Kind.DIMENSIONLESS, [A_],
                   "log-normal 1-sigma of the DFT Ku (3 methods of L05 agree only in trend)"))
    out.append(est("aniso_ratio_Zn2_over_Fe18_0K", dd.hosts["dft"]["A"] / dd.hosts["dft"]["K_fe"] * 1.0,
                   dd.hosts["dft"]["A"] / dd.hosts["dft"]["K_fe"], "1", Kind.DIMENSIONLESS, [A_],
                   "K1(SrZn2W)/K1(SrFe18O27) at 0 K = LN(0, 0.5): no Zn2 DFT in fit (bracketed by DFT = 1 and L11)"))
    out.append(est("aniso_K1_SrFe18O27_0K_L11", dl.hosts["l11"]["K_fe"], dl.hosts["l11"]["K_fe"], "J/m^3", ed,
                   list(hin.ids["l11"]["Fe"]) + [hin.ids["ku"], hin.ids["ku_srba"]],
                   "L11 BaFe18 300 K: mu0 Ha Ms/2 / m^3 x Ba->Sr factor (L05-103/L05-112 sensitivity); tier 3, ALTERNATIVE host; "
                   "|K| only (sign-blind), formula valid for H >> Ha"))
    out.append(est("aniso_A_Zn2_0K_L11", dl.hosts["l11"]["A"], dl.hosts["l11"]["A"], "J/m^3", ed, list(hin.ids["l11"]["Zn"]),
                   "L11 BaZn2W 300 K -> 0 K, Ba -> Sr by the L05 sensitivity; tier 3; ALTERNATIVE"))
    out.append(est("aniso_A_Mn2_0K_L11", dl.hosts["l11"]["A_mn"], dl.hosts["l11"]["A_mn"], "J/m^3", ed, list(hin.ids["l11"]["Mn"]),
                   f"L11 BaMn2W 300 K -> 0 K; tier 3, sigma_ln {SIGMA_LN_L11_MN}; the L11 Ms of this sample is anomalously low "
                   "(47.8 A m2/kg against 73-95 of the other samples); sign from the Mn valence hypothesis (assumption)"))
    out.append(const("aniso_Ba_to_Sr_factor", ba_to_sr_factor(hin), "1",
                   Kind.DIMENSIONLESS, [hin.ids["ku"], hin.ids["ku_srba"]],
                   "K(Sr)/K(Ba) = K(0)/(2K(0.5)-K(0)): in-fit L05 DFT sensitivity to Ba, linear in the Ba fraction"))
    out.append(const("aniso_sigma_ln_Ba_to_Sr", SIGMA_LN_BA_SR, "1", Kind.DIMENSIONLESS, [A_],
                   "width of ln(factor): K(y) ~ y^q with q in [0.5, 2] changes the factor between 0.20 and 0.41"))
    out.append(const("aniso_L11_ratio_Co2_over_Zn2_300K", l11_co_zn_ratio(), "1",
                   Kind.DIMENSIONLESS, [*hin.ids["l11"]["Co"], *hin.ids["l11"]["Zn"]],
                   "|K1(BaCo2W)|/|K1(BaZn2W)| = Ha Ms rho ratio (L11, 300 K): target of the optional 'co_zn_ratio' branch, "
                   "sigma_ln 0.7 (assumption); NOT used in the forecast of record"))
    out.append(const("aniso_E_Co_max_meV", E_CO_MAX_MEV, "1", Kind.DIMENSIONLESS, [A_],
                   "ceiling of the Co2+ single-ion energy rho A V_fu at 0 K, meV (general spin-orbit energy scale; truncates "
                   "the prior of rho)"))
    out.append(const("aniso_Mn2_axial_weight", MN2_AXIAL_WEIGHT, "1", Kind.DIMENSIONLESS, [A_],
                   "P(Mn2+ contribution axial like the host); S ion"))
    out.append(const("aniso_Mn3_axial_weight", MN3_SIGN_WEIGHT, "1", Kind.DIMENSIONLESS, [A_],
                   "P(Mn3+ Jahn-Teller single-ion term axial); sign 50/50"))
    out.append(est("aniso_q2_K2_over_Kscale", dd.q2, pri["q2"], "1", Kind.DIMENSIONLESS, [A_],
                   "prior N(0,0.3) only; enters T_T through K1+2K2=0; not identifiable"))
    out.append(est("aniso_q4_K4_over_Kscale", dd.q4, pri["q4"], "1", Kind.DIMENSIONLESS, [A_],
                   "prior N(0,0.03) only (basal constant, 'K3 sin^6 cos6phi' of papers); not identifiable"))
    # the two Co forms, each against its own T_T posterior
    for cf in CO_FORMS:
        d = joint_draws(n, seed, host="dft", co_form=cf)
        out.append(est(f"aniso_rho_Co_over_A[{cf}]", d.rho, pri["rho"], "1", Kind.DIMENSIONLESS, tt_ids,
                       f"FIT (T_T of 4 samples + structure constraints), co_form={cf}: |Co2+ term| / A per Co at 0 K, truncated by "
                       "the 20 meV ceiling; correlated with the shape parameter"))
        if cf == "cc_power":
            out.append(est("aniso_p_Co_exponent[cc_power]", d.p, pri["p"], "1", Kind.DIMENSIONLESS, tt_ids,
                           "EFFECTIVE exponent of F = m^p, prior U(3,21); NOT identifiable (identifiable is the combination, see "
                           "aniso_CoTerm_over_host_300K_Co2); not a consequence of the Callen-Callen theorem"))
        else:
            out.append(est("aniso_theta_Co_over_Tc[ion_mf]", d.f, pri["f"], "1", Kind.DIMENSIONLESS, tt_ids,
                           "Theta_Co / T_C of the Co2+ sublattice (mean-field Brillouin S=3/2), prior U(0.15,1)"))
        # the identifiable combination: the Co share of the 300 K anisotropy of Co2 (T_C of the sample = L02 Table 3, 790 K)
        tco = float(C.fit_inputs().tc_dmax[0])

        def share(rho, p, f_, s, cf=cf):
            m3 = np.maximum(C.reduced_magnetization(300.0 / tco, s), 1e-12)
            return 2.0 * rho * co_F(cf, m3, 300.0 / tco, p, f_) / m3 ** EXP_HOST

        out.append(est(f"aniso_CoTerm_over_host_300K_Co2[{cf}]", share(d.rho, d.p, d.f, d.shape_s),
                       share(pri["rho"], pri["p"], pri["f"], pri["s"]), "1", Kind.DIMENSIONLESS, tt_ids,
                       f"(Co term)/(host term) of SrCo2W at 300 K, co_form={cf}: the combination the T_T data fix"))
    return out
