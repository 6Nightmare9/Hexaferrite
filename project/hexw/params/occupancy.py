"""Site-preference model of the A cations (Co, Zn, Mg, Ni, Mn) on the seven W-sites, fitted to the FIT literature only.

Model (a lattice gas with Fe as the reference species, frozen at an effective temperature T_f):

    p[s, i] = w[s, i] / (1 + sum_s' w[s', i]),    w[s, i] = exp(m_s - d[s, i]),    d[s, i] = (eps[s, i] - eps[s, W3]) / (k T_f)

p[s, i] is the fraction of site i occupied by species s (the rest is Fe). The chemical-potential-like shift m_s is NOT a
parameter: it is solved for every composition from sum_i N_i p[s, i] = n_s (N_i = sites per formula unit; Lagrange
multiplier / iterative proportional fitting), so capacity (p summed over species <= 1) and the A content are exact.
Only the dimensionless d[s, i] (6 per species, reference site W3) are fitted; eps and T_f separately are degenerate
(T_f is an explicit non-identifiable parameter and does not change any prediction).

Data (records of literature/fit only, through hexw.params.literature):
  * L02 Co (SrCo2W), Zn (SrZn2W) and Co+Zn (SrCoZnW), sites W1..W7 (Cmcm sub-sites weighted by multiplicity);
  * L01 Mg and Ni, sites W1..W6.
Co on W3 (6g): the only source is L02; there is no cross-check in the fit.

Weighted least squares in p-space; sigma_eff^2 = sigma_CIF^2 + sigma_sys^2 (sigma_sys is an ASSUMPTION: Rietveld standard
uncertainties are not systematic errors; it is a PARAMETER with the prior U(0.01, 0.03), marginalised in the bootstrap: every
replicate draws its own sigma_sys, resamples the fractions with it and refits with it). Intervals: parametric bootstrap
(resampling of the observed fractions); the same bootstrap sample is the ensemble used for propagation to the structure
candidates and Ms0. A relative energy is "robustly identified" if its site fraction is >= 3 sigma_eff over the WHOLE range
of sigma_sys (i.e. at 0.03); `n_identified_at_min` counts it at 0.01 only.

Weights of this ensemble describe the distribution FROZEN at synthesis/quench (T_f), not the one at the measuring temperature.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Optional

import numpy as np
from scipy.optimize import least_squares

from hexw.params.literature import load_fit_records
from hexw.registry import schemas as S
from hexw.structures.sites import N_FU, SITES
from hexw.units import Kind

A_MIN, A_MAX = -9.0, 9.0          # bounds of the logit a = m - d in the fit (p in [1.2e-4, 0.9999])
A_ABSENT = -40.0                  # "no such site for this species" (no data): p = 0
SIGMA_SYS_RANGE = (0.01, 0.03)    # ASSUMPTION: prior U(low, high) of the systematic floor on a refined site fraction (absolute)
SIGMA_SYS = 0.5 * (SIGMA_SYS_RANGE[0] + SIGMA_SYS_RANGE[1])      # central value (point fit, residuals, upper-limit bounds)
Z_IDENT = 3.0                     # a site fraction is "identified" if p / sigma_eff >= Z_IDENT in the reference sample
REF_SITE = 2                      # W3: the reference site (d = 0)
N_BOOT = 300
BOOT_SEED = 20261026
SPECIES_FIT = ("Co", "Zn", "Mg", "Ni")
MN_D_HALF = 3.0                   # assumption: Mn prior d[Mn, i] ~ U(-3, +3) kT_f for i in W1,W2,W4,W5,W6
T_FREEZE = (1100.0, 400.0, 1523.0)    # assumption: value, low, high (K); see docstring of freezing_temperature_estimate

# CIF sub-site multiplicities (Cmcm cell of L02): 6g -> 8d + 4b, 12k -> 16h + 8f
_SPLIT_W = {2: (8.0, 4.0), 4: (16.0, 8.0)}


class OccupancyDataError(RuntimeError):
    pass


# --------------------------------------------------------------------------------------------- data
@dataclass
class SpeciesObs:
    p: np.ndarray                       # observed fraction of the site, nan = no datum
    sigma: np.ndarray                   # standard uncertainty of the CIF (0 where not stated), nan where no datum
    rec_ids: list                       # list (len 7) of tuples of rec_ids behind each value

    def n(self) -> float:
        m = ~np.isnan(self.p)
        return float((self.p[m] * N_FU[m]).sum())

    def sigma_eff(self, sigma_sys: Optional[float] = None) -> np.ndarray:
        ss = SIGMA_SYS if sigma_sys is None else sigma_sys
        return np.sqrt(np.nan_to_num(self.sigma) ** 2 + ss ** 2)


@dataclass
class OccSample:
    name: str
    lit_id: str
    species: dict                       # species -> SpeciesObs


def _site_index(s: str) -> int:
    return int(s) - 1


def load_occupancy_samples(records=None) -> list[OccSample]:
    """Build the observation set from the fit records (L02: Co/Zn/CoZn, L01: Mg/Ni)."""
    records = load_fit_records() if records is None else records
    raw: dict = {}              # (sample, species) -> {site_idx -> {sub: (p, sigma, rec_id)}}
    lit_of = {}
    for r in records:
        e = r.estimate
        if e is None or not r.name.startswith("site_occupancy_") or r.name.endswith("_derived"):
            continue
        sample = r.source.sample_id_in_source
        m = re.match(r"^site_occupancy_(Co|Zn)(\d)(?:_(\d))?$", r.name)
        if m:
            sp, site, sub = m.group(1), _site_index(m.group(2)), int(m.group(3) or 0)
        else:
            m = re.match(r"^site_occupancy_(Mg|Ni)_W(\d)$", r.name)
            if not m:
                raise OccupancyDataError(f"unrecognised occupancy record name {r.name!r} ({r.rec_id})")
            sp, site, sub = m.group(1), _site_index(m.group(2)), 0
        sig = e.uncertainty.std if e.uncertainty.kind == "std" else 0.0
        raw.setdefault((sample, sp), {}).setdefault(site, {})[sub] = (float(e.value), float(sig or 0.0), r.rec_id)
        lit_of[sample] = r.lit_id
    samples: dict = {}
    for (sample, sp), sites in raw.items():
        p = np.full(7, np.nan)
        sg = np.full(7, np.nan)
        ids = [()] * 7
        for i, subs in sites.items():
            if len(subs) == 1 and 0 in subs:
                p[i], sg[i], rid = subs[0]
                ids[i] = (rid,)
            else:
                if i not in _SPLIT_W or set(subs) != {1, 2}:
                    raise OccupancyDataError(f"{sample} {sp} W{i + 1}: unexpected sub-sites {sorted(subs)}")
                w1, w2 = _SPLIT_W[i]
                (p1, s1, r1), (p2, s2, r2) = subs[1], subs[2]
                tot = w1 + w2
                p[i] = (w1 * p1 + w2 * p2) / tot
                sg[i] = math.hypot(w1 * s1, w2 * s2) / tot
                ids[i] = (r1, r2)
        samples.setdefault(sample, OccSample(sample, lit_of[sample], {}))
        samples[sample].species[sp] = SpeciesObs(p, sg, ids)
    out = [samples[k] for k in sorted(samples)]
    if not out:
        raise OccupancyDataError("no site-occupancy records found in the fit file")
    return out


# --------------------------------------------------------------------------------------------- site-fraction solver
def solve_site_fractions(a, n, tol: float = 1e-11, maxit: int = 200) -> np.ndarray:
    """Site fractions p[..., s, i] of species s on site i given the base logits a[..., s, i] and the contents n[..., s]
    (atoms per formula unit). Fe fills the rest of every site. Solves sum_i N_i p[s, i] = n_s for the shifts m_s
    (Newton on log n with the exact Jacobian; this is the Lagrange-multiplier step of the constrained Boltzmann
    distribution). Raises ValueError if a content cannot be accommodated."""
    a = np.asarray(a, float)
    n = np.asarray(n, float)
    ns = a.shape[-2]
    if np.any(n <= 0):
        raise ValueError("contents must be positive (drop absent species before the call)")
    shift = np.zeros(a.shape[:-1])

    def model(sh):
        w = np.exp(np.clip(a + sh[..., None], -60, 60))
        p = w / (1.0 + w.sum(-2, keepdims=True))
        return p, (p * N_FU).sum(-1)

    p, nm = model(shift)
    for _ in range(maxit):
        g = np.log(np.maximum(nm, 1e-300)) - np.log(n)
        if np.max(np.abs(nm - n)) < tol * max(1.0, float(np.max(n))):
            break
        # Jacobian dn_k/dm_l = sum_i N_i (delta_kl p_k - p_k p_l)
        J = -np.einsum("...ki,...li,i->...kl", p, p, N_FU)
        J = J + nm[..., :, None] * np.eye(ns)
        J = J / np.maximum(nm, 1e-300)[..., :, None]
        J = J + 1e-12 * np.eye(ns)
        step = -np.linalg.solve(J, g[..., None])[..., 0]
        step = np.clip(step, -4.0, 4.0)
        shift = shift + step
        p, nm = model(shift)
    else:
        if np.max(np.abs(nm - n)) > 1e-7 * max(1.0, float(np.max(n))):
            raise ValueError("site-fraction solver did not converge: content exceeds the capacity of the sites "
                             f"(max |n_model - n| = {np.max(np.abs(nm - n)):.3g})")
    return p


# --------------------------------------------------------------------------------------------- the fit
@dataclass
class OccupancyFit:
    species: tuple                      # fitted species, order of the arrays below
    a_hat: np.ndarray                   # (S, 7) best-fit logits in the REFERENCE sample of each species
    boot_a: np.ndarray                  # (B, S, 7) bootstrap logits (the ensemble)
    d_hat: np.ndarray                   # (S, 7)  d = a[W3] - a[i]
    boot_d: np.ndarray                  # (B, S, 7)
    has_site: np.ndarray                # (S, 7) bool: the data say something about (species, site)
    identified: np.ndarray              # (S, 7) bool
    d_lower: np.ndarray                 # (S, 7) 95%-style lower bound of d for the non-identified sites (nan otherwise)
    samples: list
    chi2: float
    dof: int
    n_fitted: int                       # number of fitted logits (free parameters of the least-squares problem)
    n_identified: int                   # number of ROBUSTLY identified relative energies d (>= 3 sigma_eff over the whole sigma_sys range)
    n_identified_at_min: int = 0        # the same at the smallest sigma_sys of the prior only
    identified_at_min: Optional[np.ndarray] = None
    boot_sigma_sys: Optional[np.ndarray] = None     # (B,) sigma_sys of every bootstrap replicate
    estimates: list = field(default_factory=list)
    residuals: list = field(default_factory=list)       # (sample, species, site, p_obs, p_fit, sigma_eff)
    ref_sample: dict = field(default_factory=dict)

    def index(self, sp: str) -> int:
        return self.species.index(sp)


def _pick_reference(samples: list[OccSample]) -> dict:
    ref = {}
    for sp in {sp for s in samples for sp in s.species}:
        cands = [(len(s.species) == 1, s.species[sp].n(), s.name) for s in samples if sp in s.species]
        ref[sp] = max(cands)[2]
    return ref


def _fit_once(samples: list[OccSample], species: tuple, ref: dict, p_over: Optional[dict] = None,
              start: Optional[np.ndarray] = None, sig_sys: Optional[float] = None):
    """One weighted LSQ fit. `p_over[(sample, species)]` replaces the observed fractions (bootstrap); `sig_sys` is the
    systematic floor of this replicate (None = the central SIGMA_SYS)."""
    S_ = len(species)
    mask = np.zeros((S_, 7), bool)
    for s in samples:
        for sp, o in s.species.items():
            mask[species.index(sp)] |= ~np.isnan(o.p)
    free = np.argwhere(mask)

    def p_obs(s, sp):
        return (p_over or {}).get((s.name, sp), s.species[sp].p)

    # initial logits from the reference samples
    a0 = np.full((S_, 7), A_ABSENT)
    for k, i in free:
        sp = species[k]
        s = next(x for x in samples if x.name == ref[sp])
        v = p_obs(s, sp)[i]
        v = 0.5 if np.isnan(v) else float(np.clip(v, 1e-4, 1 - 1e-4))
        a0[k, i] = float(np.clip(math.log(v / (1 - v)), A_MIN, A_MAX))
    x0 = np.array([a0[k, i] for k, i in free]) if start is None else np.array([start[k, i] for k, i in free])

    def unpack(x):
        a = np.full((S_, 7), A_ABSENT)
        a[free[:, 0], free[:, 1]] = x
        return a

    def resid(x):
        a = unpack(x)
        out = []
        for s in samples:
            sps = list(s.species)
            idx = [species.index(sp) for sp in sps]
            if len(sps) == 1 and ref[sps[0]] == s.name:
                p = 1.0 / (1.0 + np.exp(-a[idx[0]]))
                pp = [p]
            else:
                n = np.array([max(_n_of(p_obs(s, sp)), 1e-6) for sp in sps])
                pf = solve_site_fractions(a[idx], n)
                pp = [pf[j] for j in range(len(sps))]
            for j, sp in enumerate(sps):
                o = s.species[sp]
                po = p_obs(s, sp)
                m = ~np.isnan(po)
                out.append(((pp[j] - po)[m]) / o.sigma_eff(sig_sys)[m])
        return np.concatenate(out)

    lo = np.full(x0.shape, A_MIN)
    hi = np.full(x0.shape, A_MAX)
    sol = least_squares(resid, np.clip(x0, A_MIN + 1e-6, A_MAX - 1e-6), bounds=(lo, hi), xtol=1e-9, ftol=1e-9,
                        gtol=1e-9, method="trf")
    return unpack(sol.x), sol, mask


def _n_of(p):
    m = ~np.isnan(p)
    return float((p[m] * N_FU[m]).sum())


def fit_occupancy(samples: Optional[list] = None, n_boot: int = N_BOOT, seed: int = BOOT_SEED) -> OccupancyFit:
    """Fit the logits of every species/site, bootstrap them, classify identifiability, build the ParameterEstimates."""
    samples = load_occupancy_samples() if samples is None else samples
    species = tuple(sp for sp in SPECIES_FIT if any(sp in s.species for s in samples))
    species = species + tuple(sorted({sp for s in samples for sp in s.species} - set(species)))
    ref = _pick_reference(samples)
    a_hat, sol, mask = _fit_once(samples, species, ref)
    chi2 = float(np.sum(sol.fun ** 2))
    n_data = sol.fun.size
    n_fit = int(mask.sum())
    dof = n_data - n_fit

    rng = np.random.default_rng(seed)
    boot = np.empty((n_boot, len(species), 7))
    boot_ss = rng.uniform(*SIGMA_SYS_RANGE, n_boot)                  # sigma_sys is marginalised, not fixed
    for b in range(n_boot):
        over = {}
        for s in samples:
            for sp, o in s.species.items():
                over[(s.name, sp)] = np.clip(o.p + rng.normal(size=7) * o.sigma_eff(boot_ss[b]), 0.0, 1.0)
        boot[b], _, _ = _fit_once(samples, species, ref, over, start=a_hat, sig_sys=boot_ss[b])

    def to_d(a):
        return a[..., REF_SITE:REF_SITE + 1] - a
    d_hat = to_d(a_hat)
    boot_d = to_d(boot)

    identified = np.zeros_like(mask)
    identified_min = np.zeros_like(mask)
    d_lower = np.full(mask.shape, np.nan)
    for k, sp in enumerate(species):
        rs = next(x for x in samples if x.name == ref[sp])
        o = rs.species[sp]
        se = o.sigma_eff()
        for i in range(7):
            if not mask[k, i]:
                continue
            po = o.p[i]
            if np.isnan(po):                       # datum only in a non-reference sample
                identified[k, i] = False
                continue
            se_hi = o.sigma_eff(SIGMA_SYS_RANGE[1])[i]
            se_lo = o.sigma_eff(SIGMA_SYS_RANGE[0])[i]
            identified_min[k, i] = po / se_lo >= Z_IDENT
            if po / se_hi >= Z_IDENT:                       # robust: identified over the whole sigma_sys range
                identified[k, i] = True
            else:
                pul = min(max(po, 0.0) + 2.0 * se[i], 0.999)
                a_ul = math.log(pul / (1 - pul))
                d_lower[k, i] = a_hat[k, REF_SITE] - a_ul
    fit = OccupancyFit(species=species, a_hat=a_hat, boot_a=boot, d_hat=d_hat, boot_d=boot_d, has_site=mask,
                       identified=identified, d_lower=d_lower, samples=samples, chi2=chi2, dof=dof, n_fitted=n_fit,
                       n_identified=int((identified & (np.arange(7) != REF_SITE)).sum()), ref_sample=ref,
                       n_identified_at_min=int((identified_min & (np.arange(7) != REF_SITE)).sum()),
                       identified_at_min=identified_min, boot_sigma_sys=boot_ss)
    fit.estimates = _build_estimates(fit)
    fit.residuals = _residual_table(fit)
    return fit


def _residual_table(fit: OccupancyFit) -> list:
    rows = []
    for s in fit.samples:
        sps = list(s.species)
        idx = [fit.index(sp) for sp in sps]
        if len(sps) == 1 and fit.ref_sample[sps[0]] == s.name:
            pf = [1.0 / (1.0 + np.exp(-fit.a_hat[idx[0]]))]
        else:
            n = np.array([max(s.species[sp].n(), 1e-6) for sp in sps])
            pm = solve_site_fractions(fit.a_hat[idx], n)
            pf = [pm[j] for j in range(len(sps))]
        for j, sp in enumerate(sps):
            o = s.species[sp]
            for i in range(7):
                if not np.isnan(o.p[i]):
                    rows.append((s.name, sp, SITES[i], float(o.p[i]), float(pf[j][i]), float(o.sigma_eff()[i])))
    return rows


def _rec_ids(fit: OccupancyFit, sp: str, i: int) -> list[str]:
    out = []
    for s in fit.samples:
        if sp in s.species:
            out.extend(s.species[sp].rec_ids[i])
    return sorted(set(out))


def _build_estimates(fit: OccupancyFit) -> list:
    est = []
    for k, sp in enumerate(fit.species):
        for i, site in enumerate(SITES):
            if i == REF_SITE or not fit.has_site[k, i]:
                continue
            ids = _rec_ids(fit, sp, i)
            prov = S.Provenance(
                source_id=";".join(ids), locator=f"reference sample {fit.ref_sample[sp]}", status=S.ProvenanceStatus.DERIVED_BY_US,
                note=("relative site-preference energy d = (eps[%s,%s] - eps[%s,W3])/(k T_f) of a frozen lattice gas "
                      "fitted by weighted LSQ to the fit-literature site fractions" % (sp, site, sp)))
            name = f"occ_dE_{sp}_{site}"
            if fit.identified[k, i]:
                u = S.Uncertainty(kind="std", std=float(max(np.std(fit.boot_d[:, k, i], ddof=1), 1e-6)),
                                  convention="sd")
                est.append(S.ParameterEstimate(name, float(fit.d_hat[k, i]), "1", Kind.DIMENSIONLESS,
                                               S.SourceRole.LITERATURE_FIT, prov, u, identifiable=True))
            else:
                lo = float(fit.d_lower[k, i]) if not np.isnan(fit.d_lower[k, i]) else float(fit.d_hat[k, i])
                width = 1.01 * max(8.0, abs(lo))
                prov.note += ("; NOT robustly identified: the site fraction is < %g sigma_eff somewhere in the sigma_sys range "
                              "(%s identified at sigma_sys = %g only); only a lower bound of d (upper limit of p) is "
                              "known where it is consistent with 0" % (Z_IDENT, "is" if fit.identified_at_min[k, i] else "not",
                                                                       SIGMA_SYS_RANGE[0]))
                u = S.Uncertainty(kind="interval", lower=lo, upper=lo + width, declared_wide=True, inequality="ge")
                est.append(S.ParameterEstimate(name, lo, "1", Kind.DIMENSIONLESS, S.SourceRole.LITERATURE_FIT, prov, u,
                                               identifiable=False))
    est.append(freezing_temperature_estimate())
    est.append(sigma_sys_estimate())
    return est


def sigma_sys_estimate() -> S.ParameterEstimate:
    """sigma_sys: ASSUMPTION, a parameter with the prior U(0.01, 0.03), marginalised in the bootstrap."""
    lo, hi = SIGMA_SYS_RANGE
    prov = S.Provenance(source_id="ASSUMPTION:sigma_sys", status=S.ProvenanceStatus.ASSUMPTION,
                        locator="Rietveld standard uncertainties of the CIFs are not systematic errors",
                        note="systematic floor of a refined site fraction, absolute; prior U(0.01, 0.03) (value = its mean), marginalised in the bootstrap")
    return S.ParameterEstimate("occ_sigma_sys", 0.5 * (lo + hi), "1", Kind.DIMENSIONLESS, S.SourceRole.ASSUMPTION, prov,
                               S.Uncertainty(kind="unknown"), identifiable=False)


def freezing_temperature_estimate() -> S.ParameterEstimate:
    """T_f is degenerate with the energy scale (only d = eps/kT_f enters any prediction): not identifiable, wide."""
    v, lo, hi = T_FREEZE
    prov = S.Provenance(
        source_id="ASSUMPTION:T_f", status=S.ProvenanceStatus.ASSUMPTION,
        locator="upper end = highest calcination temperature in L02-003/-005/-007 (1523 K); L03-005/-007 (1423-1473 K)",
        note=("effective cation-freezing temperature during cooling (below the synthesis temperature); not measurable "
              "from the fit data; degenerate with eps. It converts d to eV (eps = d k T_f) and changes no prediction"))
    return S.ParameterEstimate("occ_T_freeze", v, "K", Kind.TEMPERATURE, S.SourceRole.ASSUMPTION, prov,
                               S.Uncertainty(kind="interval", lower=lo, upper=hi, declared_wide=True),
                               identifiable=False)


def energy_eV(d: float, T_f: float = T_FREEZE[0]) -> float:
    """eps(species, site) - eps(species, W3) in eV for an assumed freezing temperature."""
    return d * 8.617333262e-5 * T_f


# --------------------------------------------------------------------------------------------- Mn prior
@dataclass(frozen=True)
class MnPrior:
    """No Mn site fractions exist in the fit corpus (only L11 BaMn2W, tier 3, no site occupancies). The prior is an
    ASSUMPTION built from general considerations only: Mn has no strong, literature-supported site preference on W1..W6,
    so d[Mn, i] ~ U(-MN_D_HALF, +MN_D_HALF) (kT_f); on W7 no fitted species (Co, Zn) has ever been found -> the
    smallest lower bound of d among them is used as the lower end."""
    half: float
    w7_lower: float
    w7_width: float


def mn_prior(fit: OccupancyFit) -> MnPrior:
    lows = []
    for sp in ("Co", "Zn"):
        if sp in fit.species:
            k = fit.index(sp)
            if fit.has_site[k, 6] and not fit.identified[k, 6]:
                lows.append(float(fit.d_lower[k, 6]))
    lo = min(lows) if lows else 4.0
    return MnPrior(MN_D_HALF, lo, 1.01 * max(8.0, lo))


def sample_mn_a(fit: OccupancyFit, size: int, rng: np.random.Generator) -> np.ndarray:
    """(size, 7) base logits a = -d of Mn (reference d[W3] = 0)."""
    pr = mn_prior(fit)
    d = rng.uniform(-pr.half, pr.half, size=(size, 7))
    d[:, REF_SITE] = 0.0
    d[:, 6] = rng.uniform(pr.w7_lower, pr.w7_lower + pr.w7_width, size=size)
    return -d


def mn_estimates(fit: OccupancyFit) -> list:
    pr = mn_prior(fit)
    prov = S.Provenance(
        source_id="ASSUMPTION:Mn_site_preference", status=S.ProvenanceStatus.ASSUMPTION,
        locator="no site occupancies of Mn in the fit corpus (only L11 BaMn2W, tier 3: Ba->Sr transfer unjustified)",
        note="uninformative prior from general considerations; weight/prior not determined by the Sr-W literature")
    out = []
    for i, site in enumerate(SITES):
        if i == REF_SITE:
            continue
        if i == 6:
            u = S.Uncertainty(kind="interval", lower=pr.w7_lower, upper=pr.w7_lower + pr.w7_width, declared_wide=True,
                              inequality="ge")
            v = pr.w7_lower
        else:
            u = S.Uncertainty(kind="interval", lower=-pr.half, upper=pr.half, declared_wide=True)
            v = 0.0
        out.append(S.ParameterEstimate(f"occ_dE_Mn_{site}", v, "1", Kind.DIMENSIONLESS, S.SourceRole.ASSUMPTION,
                                       prov, u, identifiable=False))
    return out


# --------------------------------------------------------------------------------------------- cached default fit
@lru_cache(maxsize=1)
def default_fit() -> OccupancyFit:
    """The fit of the full fit-literature set (computed on first use, never at import)."""
    return fit_occupancy()


def parameter_table(fit: Optional[OccupancyFit] = None) -> list:
    fit = default_fit() if fit is None else fit
    return list(fit.estimates) + mn_estimates(fit)


# --------------------------------------------------------------------------------------------- ensemble draws
SUPPORTED_A = ("Co", "Zn", "Mg", "Ni", "Mn")


def clean_composition(composition: dict) -> dict:
    """Validated {species: x per formula unit}; zero entries dropped."""
    if not isinstance(composition, dict):
        raise ValueError("composition must be a dict {A: x} (empty = SrFe18O27)")
    out = {}
    for sp, x in composition.items():
        if sp not in SUPPORTED_A:
            raise ValueError(f"unsupported substituent {sp!r}; supported: {SUPPORTED_A}")
        if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) or x < 0:
            raise ValueError(f"bad content of {sp}: {x!r}")
        if x > 0:
            out[sp] = float(x)
    return out


def site_fraction_draws(composition: dict, fit: Optional[OccupancyFit] = None, n: Optional[int] = None,
                        seed: int = 0, scale: Optional[np.ndarray] = None):
    """Draws of the frozen distribution of the A cations: (species tuple, p[n, S, 7]) with p the fraction of each
    site occupied by the species (Fe = the rest). Co/Zn/Mg/Ni: one draw per bootstrap fit (cycled if n > B);
    Mn: the (assumption) prior. Capacity and A content are exact for every draw. `scale` (n,) multiplies the A contents
    draw by draw (actual / nominal content)."""
    fit = default_fit() if fit is None else fit
    comp = clean_composition(composition)
    species = tuple(sp for sp in SUPPORTED_A if sp in comp)
    nb = fit.boot_a.shape[0]
    n = nb if n is None else int(n)
    rng = np.random.default_rng(seed)
    if not species:
        return species, np.zeros((n, 0, 7))
    a = np.empty((n, len(species), 7))
    for j, sp in enumerate(species):
        if sp == "Mn":
            a[:, j] = sample_mn_a(fit, n, rng)
        else:
            if sp not in fit.species:
                raise ValueError(f"no fit data for {sp}")
            a[:, j] = fit.boot_a[np.arange(n) % nb, fit.index(sp)]
    tgt = np.broadcast_to(np.array([comp[sp] for sp in species]), (n, len(species)))
    if scale is not None:
        scale = np.asarray(scale, float)
        if scale.shape != (n,) or np.any(scale <= 0):
            raise ValueError("scale must be positive with shape (n,)")
        tgt = tgt * scale[:, None]
    return species, solve_site_fractions(a, tgt)
