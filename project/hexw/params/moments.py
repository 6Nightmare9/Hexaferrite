"""Saturation magnetisation at 0 K, Ms0, of SrFe(18-x)A(x)O27 in the collinear Gorter ferrimagnet model.

    M0 [muB / f.u.] = sum_i sigma_i * sum_species N_i * p[species, i] * m(species)

sigma_i = +1 (spin up: W3, W4, W5, W7) / -1 (down: W1, W2, W6) -- the ferrimagnetic scheme refined in L01 (Fig. 4) and
L02 (mCIF) and used in the L05 DFT; p = frozen cation distribution (hexw.params.occupancy), Fe on the rest of each site.

Ionic moments (ASSUMPTIONS, identifiable=False, never presented as literature numbers):
    spin-only g S muB: Fe3+ 5, Fe2+ 4, Co2+ 3 (L02 also fixes Fe3+:Co2+ = 5:3), Ni2+ 2, Mn2+ 5, Mn3+ 4, Zn2+/Mg2+ 0;
    orbital contribution (added parallel to the spin for more-than-half-filled shells): ranges Co2+ U(0, 1), Fe2+ U(0, 0.5),
    Ni2+ U(0, 0.3) muB; Mn3+ U(-0.3, 0). No fit-corpus value of an orbital moment of Co exists (L05 gives Fe only, ~0.02 muB).
The authors' 'assumed' dopant moments of L11 (Co 4.0, Ni 3.0, Mn 2.7) are NOT used (tier 3, flagged in CONFLICTS F); the Co range
simply contains 4.0.

The fit moments (L02 neutron, L05 DFT, L01 NPD) are used for CHECKS only (`sanity_checks`): they have different moment
definitions (Fe ion / site average) and temperatures (300 K vs 0 K) and are never mixed in a calibration.

Ms(T) = eta * Ms0 * m(T / T_C)   (`ms_T`). eta is an EFFECTIVE moment-reduction multiplier, estimated from the three fit samples
that have Ms at a known temperature (SrFe18O27: C004-002/-004/-006, 300 K; SrMg2W: L01-069, SrNi2W: L01-077, 300 K, VSM,
impurity-corrected; T_C of the L01 samples from L01-070/-078): eta_j = Ms_obs / (Ms0 m(300 K / T_C)). The inter-sample scatter is
carried into the prediction (3 samples -> 1 parameter + between-sample scatter, Student-t predictive). eta absorbs covalent
reduction of the moments, non-collinearity and the under-saturation of the VSM -- it is NOT a physical constant of the Fe ion.

Joint draws: `Ms0(..., draws=JointDraws)` takes the content ratio r and the Mn valence branch from hexw.params.anisotropy.JointDraws, the
same that T_C and K use.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Optional

import numpy as np
from scipy import constants as C

from hexw.params import occupancy as O
from hexw.params.literature import load_fit_records
from hexw.registry import schemas as S
from hexw.structures.ensemble import FE2_SCENARIOS, MN_VALENCE_PRIOR, Hypothesis, fe2_scenario_weights, valence_counts
from hexw.structures.sites import N_FU, SITES, SPIN_FE2, SPIN_FE3, SPIN_SIGN
from hexw.units import Kind

MU_B = C.physical_constants["Bohr magneton"][0]            # A m^2
N_A = C.N_A
ATOMIC_MASS = {"Sr": 87.62, "Fe": 55.845, "Co": 58.933, "Zn": 65.38, "Mg": 24.305, "Ni": 58.693, "Mn": 54.938,
               "O": 15.999}                                  # g/mol, IUPAC standard atomic weights

SPIN_MOMENT = {"Fe3+": SPIN_FE3, "Fe2+": SPIN_FE2, "Co": 3.0, "Ni": 2.0, "Mn2": 5.0, "Mn3": 4.0, "Zn": 0.0, "Mg": 0.0}
ORBITAL_RANGE = {"Fe3+": (0.0, 0.0), "Fe2+": (0.0, 0.5), "Co": (0.0, 1.0), "Ni": (0.0, 0.3), "Mn2": (0.0, 0.0),
                 "Mn3": (-0.3, 0.0), "Zn": (0.0, 0.0), "Mg": (0.0, 0.0)}
V_REL_SYS = 0.004             # assumption: x-dependence of the cell volume (relative, 1 sigma)
# Actual / nominal A content. The refined contents in the fit samples are Co 1.71, Zn 1.80 (L02), Mg 1.69 (L01) per f.u. for
# nominal 2 (ratio 0.84-0.90; L03: Me < 2 in air). ASSUMPTION: a sample of nominal x contains r*x with r ~ U(0.85, 1).
DEFAULT_CONTENT_RATIO = (0.85, 1.0)


@dataclass
class Ms0Result:
    composition: dict
    hypothesis: Hypothesis
    mu_B_per_fu: np.ndarray
    spin_only_mu_B_per_fu: np.ndarray
    A_m2_per_kg: np.ndarray
    A_per_m: np.ndarray
    mass_g_per_mol: float
    volume_fu_A3: np.ndarray
    mn_valence: np.ndarray                 # per draw (0 if no Mn)
    fe2_scenario: np.ndarray               # per draw, index into scenario names
    scenario_names: tuple
    notes: list = field(default_factory=list)
    content_ratio: Optional[np.ndarray] = None      # r per draw (actual / nominal A content): shared with T_C and K via JointDraws

    def summary(self) -> dict:
        def st(v):
            q = np.percentile(v, [2.5, 16, 50, 84, 97.5])
            return {"mean": float(v.mean()), "std": float(v.std(ddof=1)), "p2.5": float(q[0]), "p16": float(q[1]),
                    "median": float(q[2]), "p84": float(q[3]), "p97.5": float(q[4])}
        return {"muB_per_fu": st(self.mu_B_per_fu), "A*m^2/kg": st(self.A_m2_per_kg), "A/m": st(self.A_per_m),
                "muB_per_fu_spin_only": st(self.spin_only_mu_B_per_fu)}

    def by_hypothesis(self) -> dict:
        out = {}
        for v in sorted(set(self.mn_valence.tolist())):
            for k, name in enumerate(self.scenario_names):
                m = (self.mn_valence == v) & (self.fe2_scenario == k)
                if m.any():
                    key = (f"Mn{int(v)}+" if v else "A2+", name)
                    out[key] = {"n": int(m.sum()), "median_muB_per_fu": float(np.median(self.mu_B_per_fu[m])),
                                "p16": float(np.percentile(self.mu_B_per_fu[m], 16)),
                                "p84": float(np.percentile(self.mu_B_per_fu[m], 84))}
        return out


def formula_mass_g_per_mol(composition: dict) -> float:
    comp = O.clean_composition(composition)
    x = sum(comp.values())
    return (ATOMIC_MASS["Sr"] + sum(v * ATOMIC_MASS[sp] for sp, v in comp.items()) + (18.0 - x) * ATOMIC_MASS["Fe"]
            + 27.0 * ATOMIC_MASS["O"])


# ---- lattice (geometry) from the fit records ------------------------------------------------------------
@lru_cache(maxsize=1)
def _default_lattice_points() -> dict:
    return _lattice_points(load_fit_records())


def _lattice_points(records=None) -> dict:
    """{sample_id: (a, c, sigma_a, sigma_c)} from the L01/L02 lattice records (cached for the default record set)."""
    if records is None:
        return _default_lattice_points()
    d: dict = {}
    for r in records:
        e = r.estimate
        if e is None or r.name not in ("lattice_a", "lattice_c"):
            continue
        sid = r.source.sample_id_in_source
        ent = d.setdefault(sid, {})
        ent[r.name] = (float(e.value), float(e.uncertainty.std or 0.0))
    return {k: (v["lattice_a"][0], v["lattice_c"][0], v["lattice_a"][1], v["lattice_c"][1])
            for k, v in d.items() if "lattice_a" in v and "lattice_c" in v}


def cell_volume_A3(a: float, c: float) -> float:
    return math.sqrt(3.0) / 2.0 * a * a * c


_LATTICE_OF = {"Co": "SrCo2Fe16O27", "Zn": "SrZn2Fe16O27", "Mg": "SrMg2Fe16O27", "Ni": "SrNi2Fe16O27"}


def volume_fu_model(composition: dict, records=None):
    """(V_fu mean, sigma) in A^3. Rule of mixtures of the measured x~2 end-member volumes (L02 Co/Zn, L01 Mg/Ni) by the
    share of each A; Mn (no lattice in the fit corpus): mean of all measured volumes. sigma: half the spread of the
    measured volumes in quadrature with V_REL_SYS (the x dependence is not measured)."""
    pts = _lattice_points(records)
    vol = {k: cell_volume_A3(v[0], v[1]) / 2.0 for k, v in pts.items()}            # Z = 2 per hexagonal cell
    comp = O.clean_composition(composition)
    x = sum(comp.values())
    allv = np.array(list(vol.values()))
    if x == 0:
        v0 = float(allv.mean())
    else:
        v0 = 0.0
        for sp, v in comp.items():
            sid = _LATTICE_OF.get(sp)
            v0 += (v / x) * (vol[sid] if sid in vol else float(allv.mean()))
    sig = math.hypot(0.5 * float(allv.max() - allv.min()), V_REL_SYS * v0)
    return float(v0), float(sig)


# ---- Ms0 -------------------------------------------------------------------------------------------------
def Ms0(composition: dict, hypothesis: Optional[Hypothesis] = None, n_draws: int = 1000, seed: int = 0,
        fit: Optional[O.OccupancyFit] = None, records=None, content_ratio: tuple = DEFAULT_CONTENT_RATIO,
        draws=None) -> Ms0Result:
    """Distribution of Ms0 for a composition {A: x} (`{}` = SrFe18O27) and a valence/Fe2+-placement hypothesis (None = all, with
    their weights). `content_ratio` = (lo, hi): the actual A content of the phase is r*x with r ~ U(lo, hi) (default: the
    refined contents of the fit samples are 0.84-0.90 of nominal); (1, 1) = exactly the nominal composition.
    `draws` (hexw.params.anisotropy.JointDraws): the SHARED draws -- r (the same as in T_C and K) and the Mn valence branch are
    taken from them, n_draws = draws.n, `content_ratio` is ignored."""
    hyp = hypothesis or Hypothesis()
    comp = O.clean_composition(composition)
    lo_r, hi_r = content_ratio
    if not 0 < lo_r <= hi_r <= 1.0 + 1e-12:
        raise ValueError("content_ratio must satisfy 0 < lo <= hi <= 1")
    if draws is not None:
        n = int(draws.n)
        r = np.asarray(draws.r, float)
        lo_r, hi_r = float(r.min()), float(r.max())
    else:
        n = int(n_draws)
        r = np.random.default_rng(seed + 2).uniform(lo_r, hi_r, size=n) if hi_r > lo_r else np.full(n, float(lo_r))
    valence_counts(comp, 2 if "Mn" in comp else None)                       # applicability (x <= 2) at the nominal content
    species, p = O.site_fraction_draws(comp, fit, n, seed, scale=r if comp else None)      # (n, S, 7)
    rng = np.random.default_rng(seed + 1)
    has_mn = "Mn" in comp
    if has_mn:
        vals = [hyp.mn_valence] if hyp.mn_valence else [2, 3]
        w = np.array([MN_VALENCE_PRIOR[v] for v in vals], float)
        if draws is not None and not hyp.mn_valence:
            mnv = np.asarray(draws.mn_valence, int)             # the shared Mn branch
        else:
            mnv = rng.choice(vals, size=n, p=w / w.sum())
    else:
        if hyp.mn_valence is not None:
            raise ValueError("mn_valence given, but no Mn in the composition")
        mnv = np.zeros(n, int)
    names = tuple(FE2_SCENARIOS)
    if hyp.fe2_scenario:
        if hyp.fe2_scenario not in FE2_SCENARIOS:
            raise ValueError(f"unknown Fe2+ scenario {hyp.fe2_scenario!r}")
        sc = np.full(n, names.index(hyp.fe2_scenario))
    else:
        fw = fe2_scenario_weights()
        sc = rng.choice(len(names), size=n, p=np.array([fw[k] for k in names]))

    x2 = np.full(n, sum(v for sp, v in comp.items() if sp != "Mn"), float)       # A2+ content (nominal)
    if has_mn:
        x2 = x2 + np.where(mnv == 2, comp["Mn"], 0.0)
    n2 = 2.0 - r * x2                                                        # charge neutrality at the actual content
    if np.any(n2 < -1e-9):
        raise ValueError("n(Fe2+) < 0: composition outside the Fe2+/Fe3+ model")
    n2 = np.maximum(n2, 0.0)
    sel = np.zeros((len(names), 7), bool)
    for k, nm in enumerate(names):
        sel[k, list(FE2_SCENARIOS[nm][0])] = True

    fe_site = N_FU * (1.0 - p.sum(1))                                       # (n, 7) Fe atoms / f.u. on the site
    cap = (fe_site * sel[sc]).sum(1)
    if np.any(n2 > cap + 1e-9):
        raise ValueError("Fe2+ scenario cannot host the required Fe2+ for some draws")
    fe2 = np.zeros_like(fe_site)
    nz = n2 > 0
    fe2[nz] = n2[nz, None] * fe_site[nz] * sel[sc][nz] / cap[nz, None]
    fe3 = fe_site - fe2

    def orb(key, size):
        lo, hi = ORBITAL_RANGE[key]
        return rng.uniform(lo, hi, size=size) if hi > lo else np.zeros(size)

    def moment_per_species(orbital: bool):
        m = np.zeros((n, len(species)))
        for j, sp in enumerate(species):
            if sp == "Mn":
                s = np.where(mnv == 3, SPIN_MOMENT["Mn3"], SPIN_MOMENT["Mn2"])
                o = np.where(mnv == 3, orb("Mn3", n), 0.0) if orbital else 0.0
            else:
                s = SPIN_MOMENT[sp]
                o = orb(sp, n) if orbital and sp in ("Co", "Ni") else 0.0
            m[:, j] = s + o
        return m

    def total(orbital: bool):
        ms = moment_per_species(orbital)
        a_part = np.einsum("nsi,ns,i->ni", p * N_FU, ms, np.ones(7)) if species else np.zeros((n, 7))
        mfe2 = SPIN_MOMENT["Fe2+"] + (orb("Fe2+", n) if orbital else np.zeros(n))
        fe_part = fe3 * SPIN_MOMENT["Fe3+"] + fe2 * mfe2[:, None]
        return ((a_part + fe_part) * SPIN_SIGN).sum(1)

    M = total(True)
    M_spin = total(False)
    mass_n = (np.array([formula_mass_g_per_mol({sp: v * ri for sp, v in comp.items()}) for ri in r]) if hi_r > lo_r
              else np.full(n, formula_mass_g_per_mol({sp: v * lo_r for sp, v in comp.items()})))
    mass = float(mass_n.mean())
    v0, vs = volume_fu_model(comp, records)
    V = np.maximum(rng.normal(v0, vs, size=n), 1e-3)
    sigma = M * MU_B * N_A / (mass_n * 1e-3)                                # A m^2 / kg
    Mv = M * MU_B / (V * 1e-30)                                             # A / m
    notes = [
        "T = 0 K, collinear Gorter ferrimagnet; ionic moments and orbital ranges are assumptions (identifiable=False)",
        "orbital moments: Co U(0,1), Fe2+ U(0,0.5), Ni U(0,0.3), Mn3+ U(-0.3,0) muB; spin_only has them all zero",
        f"formula-unit mass {mass:.3f} g/mol from stoichiometry; V_fu = {v0:.2f} +- {vs:.2f} A^3 (fit lattices L01/L02)",
        "the A distribution is the one frozen at synthesis; actual A content = r*x, r ~ U(%.2f, %.2f) (assumption from the "
        "refined contents of the fit samples)" % (lo_r, hi_r),
    ]
    return Ms0Result(comp, hyp, M, M_spin, sigma, Mv, mass, V, mnv, sc, names, notes, content_ratio=np.asarray(r, float))


def moment_parameter_estimates() -> list:
    """The moment assumptions as ParameterEstimates (identifiable=False, wide, provenance = assumption)."""
    est = []
    for key, val in SPIN_MOMENT.items():
        prov = S.Provenance(source_id=f"ASSUMPTION:spin_only_{key}", status=S.ProvenanceStatus.ASSUMPTION,
                            note="spin-only g*S*muB (g=2); not from a literature record" +
                                 ("; L02 fixes the ratio Fe3+:Co2+ = 5:3" if key == "Co" else ""))
        est.append((key, val, prov))
    out = []
    for key, val, prov in est:
        lo, hi = ORBITAL_RANGE[key]
        if hi > lo:
            prov2 = S.Provenance(source_id=f"ASSUMPTION:orbital_{key}", status=S.ProvenanceStatus.ASSUMPTION,
                                 note="orbital contribution added to the spin-only moment: a range, not a fit value")
            mid = 0.5 * (lo + hi)
            out.append(S.ParameterEstimate(f"moment_orbital_{key}", mid, "muB", Kind.MAGNETIC_MOMENT,
                                           S.SourceRole.ASSUMPTION, prov2,
                                           S.Uncertainty(kind="interval", lower=lo, upper=hi, declared_wide=True),
                                           identifiable=False))
    return out


# ---- SiteMoment helpers (definitions are never mixed) -------------------------------------------------------
def site_average_moment(site_index: int, A_fraction: dict, fe2_fraction: float = 0.0) -> S.SiteMoment:
    """Gorter spin-only moment of a site averaged over all occupants (diamagnetic ones count as 0), muB, signed."""
    fe = 1.0 - sum(A_fraction.values())
    m = fe * ((1 - fe2_fraction) * SPIN_MOMENT["Fe3+"] + fe2_fraction * SPIN_MOMENT["Fe2+"])
    for sp, f in A_fraction.items():
        m += f * SPIN_MOMENT[sp if sp in SPIN_MOMENT else "Mn2"]
    return S.SiteMoment(float(SPIN_SIGN[site_index] * m), S.MomentDefinition.SITE_AVERAGE)


def fe_ion_moment(site_index: int, fe2_share: float = 0.0) -> S.SiteMoment:
    m = (1 - fe2_share) * SPIN_MOMENT["Fe3+"] + fe2_share * SPIN_MOMENT["Fe2+"]
    return S.SiteMoment(float(SPIN_SIGN[site_index] * m), S.MomentDefinition.FE_ION)


def attach_site_moments(candidate: S.StructureCandidate) -> S.StructureCandidate:
    """Adds the spin-only Gorter moments (site average and Fe ion, separately defined) to the sites of a candidate."""
    for i, s in enumerate(candidate.sites):
        a = {k: v for k, v in s.species.items() if k not in ("Fe2+", "Fe3+")}
        fe2, fe3 = s.species.get("Fe2+", 0.0), s.species.get("Fe3+", 0.0)
        share = fe2 / (fe2 + fe3) if fe2 + fe3 > 0 else 0.0
        s.moments = [site_average_moment(i, {k: v for k, v in a.items()}, share), fe_ion_moment(i, share)]
    return candidate


# ---- sanity checks against the FIT corpus (never holdout) ------------------------------------------------------
def sanity_checks(records=None, n_draws: int = 1000, seed: int = 0) -> dict:
    """Comparison of the model with fit-literature moments. Returns {section: [rows]}. These are checks, not calibration:
    L01/L02 moments are at 300 K (the model is 0 K), L05 is DFT at 0 K."""
    records = load_fit_records() if records is None else records
    byid = {r.rec_id: r for r in records}
    rep: dict = {"L05_total": [], "L01_site_average": [], "L01_Ms": [], "L02_ZnW_Fe_ion": [], "C004_Ms_RT": []}

    # L05: SrFe18O27 (x = 0) total spin moment per cell (2 f.u.), 0 K
    r = Ms0({}, n_draws=n_draws, seed=seed, records=records, content_ratio=(1.0, 1.0))
    dft = byid["L05-093"].estimate.value / 2.0
    lo, hi = r.mu_B_per_fu.min(), r.mu_B_per_fu.max()
    so = {nm: float(np.median(r.spin_only_mu_B_per_fu[r.fe2_scenario == k])) for k, nm in enumerate(r.scenario_names)}
    rep["L05_total"].append({"rec": "L05-093", "dft_muB_per_fu": dft, "model_spin_only_by_scenario": so,
                             "model_min": float(lo), "model_max": float(hi), "model_median": float(np.median(r.mu_B_per_fu)),
                             "dft_below_model_min_by": float(lo - dft)})

    # L01 Mg / Ni: site-average moments (same definition as the Fig. 4 numbers), Gorter 0 K vs 300 K
    for sp, sid, c_site in (("Mg", "SrMg2Fe16O27", None), ("Ni", "SrNi2Fe16O27", None)):
        nA = {"Mg": 1.689, "Ni": 2.113}[sp]
        x = min(nA, 2.0)
        ens_comp = {sp: x}
        # data fractions of the sample itself (not the model): moment from the observed occupancies
        for i, site in enumerate(SITES):
            rid = next((q for q in byid if byid[q].name == f"site_magnetic_moment_{site}" and
                        byid[q].source.sample_id_in_source == sid), None)
            if rid is None:
                continue
            obs = byid[rid].estimate
            om = S.SiteMoment(obs.value, S.MomentDefinition.SITE_AVERAGE, "muB", None, obs.uncertainty.std)
            fr = next((q for q in byid.values() if q.name == f"site_occupancy_{sp}_{site}"), None)
            pA = fr.estimate.value if fr is not None else 0.0
            mm = site_average_moment(i, {sp: pA})
            rep["L01_site_average"].append({"sample": sid, "site": site, "obs_300K": obs.value, "sigma": obs.uncertainty.std,
                                            "model_0K_from_observed_occupancy": mm.value, "diff_obs_minus_model": om - mm,
                                            "z": (om - mm) / obs.uncertainty.std, "rec": rid})
        # Ms: observed totals
        ms = {r2.name: r2 for r2 in records if r2.source.sample_id_in_source == sid and r2.name.startswith("Ms_")}
        res = Ms0(ens_comp, hypothesis=Hypothesis(fe2_scenario="uniform"), n_draws=n_draws, seed=seed, records=records,
                  content_ratio=(1.0, 1.0))
        row = {"sample": sid, "model_Ms0_A_m2_kg": res.summary()["A*m^2/kg"], "x_used": x}
        for nm, rr in ms.items():
            row[nm] = (rr.rec_id, rr.estimate.value)
        rep["L01_Ms"].append(row)

    # L02 SrZn2W: Fe3+ ion moments (300 K) must not exceed the ionic limit 5 muB
    for i, site in enumerate(SITES):
        rid = next((q for q in byid if byid[q].name.startswith(f"Fe_moment_Mz_{site}_") and
                    byid[q].source.sample_id_in_source == "SrZn2Fe16O27"), None)
        if rid is None:
            continue
        obs = byid[rid].estimate
        om = S.SiteMoment(abs(obs.value), S.MomentDefinition.FE_ION, "muB", "Fe", obs.uncertainty.std)
        mm = fe_ion_moment(i)
        mm = S.SiteMoment(abs(mm.value), S.MomentDefinition.FE_ION)
        rep["L02_ZnW_Fe_ion"].append({"site": site, "obs_abs_300K": om.value, "sigma": obs.uncertainty.std,
                                      "model_Fe3+_0K": mm.value, "ratio": om.value / mm.value, "rec": rid,
                                      "within_ionic_limit": bool(om.value - 3 * (om.uncertainty or 0) <= mm.value)})

    # C004: SrFe18O27 Ms at RT (81 A m2/kg, 300 K) vs Ms0
    for rid in ("C004-002", "C004-004", "C004-006"):
        if rid in byid:
            v = byid[rid].estimate.value
            q = r.summary()["A*m^2/kg"]
            rep["C004_Ms_RT"].append({"rec": rid, "obs_300K": v, "model_Ms0_median": q["median"],
                                      "model_p16": q["p16"], "model_p84": q["p84"], "ratio_obs_over_median": v / q["median"]})
    return rep


# ---- Ms(T) = eta Ms0 m(T): the moment-reduction multiplier eta from the fit samples with Ms at a known T ----------
ETA_T_K = 300.0                              # all three Ms values were measured at room temperature (taken as 300 K)
ETA_SAMPLES = {
    # name: (composition used for Ms0 and T_C, rec_ids of Ms [A m^2/kg], rec_id of T_C or None = host T_C of the model)
    "C004_SrFe18O27": ({}, ("C004-002", "C004-004", "C004-006"), None),
    "L01_SrMg2W": ({"Mg": 1.689}, ("L01-069",), "L01-070"),          # refined content 1.69 of 2 (L01-033); impurity-corrected VSM
    "L01_SrNi2W": ({"Ni": 2.0}, ("L01-077",), "L01-078"),            # refined 2.11 > 2: the model is defined for x <= 2
}
ETA_RECORD_NAMES = {"C004-002": "Ms_RT", "C004-004": "Ms_RT", "C004-006": "Ms_RT", "L01-069": "Ms_VSM_impurity_corrected",
                    "L01-077": "Ms_VSM_impurity_corrected", "L01-070": "Curie_temperature", "L01-078": "Curie_temperature"}


@dataclass
class EtaResult:
    eta: np.ndarray                          # (n,) predictive draws of eta for a NEW sample
    per_sample: dict                         # name -> (n,) eta_j draws of the fit samples
    table: list                              # rows: name, Ms_obs, Ms0 model (median), m(300), eta_j (median, p16, p84)
    mean_ln: float
    scatter_ln: float                        # between-sample sd of ln eta (median over draws)
    rec_ids: list


def eta_draws(draws=None, *, n: int = 4000, seed: int = 20260101) -> EtaResult:
    """Effective multiplier eta of Ms(T) = eta Ms0 m(T), estimated from the fit samples with Ms at a known T (see the module
    docstring). Per draw i of the shared draws: ln eta_j = ln(Ms_obs_j) - ln(Ms0_j m_j(300 K)) for the 3 samples; the predictive
    ln eta for a new sample is mean + s sqrt(1 + 1/3) t_2 (Student-t, 3 samples, 2 dof), s = between-sample sd."""
    from hexw.params import anisotropy as AN
    from hexw.params import curie as C
    d = draws if draws is not None else AN.joint_draws(n, seed)
    recs = {c.rec_id: c for c in load_fit_records()}
    lnj = {}
    table = []
    ids = []
    for name, (comp, ms_ids, tc_id) in ETA_SAMPLES.items():
        for rid in ms_ids:
            if recs[rid].name != ETA_RECORD_NAMES[rid] or recs[rid].estimate.unit != "A*m^2/kg":
                raise RuntimeError(f"fit record {rid} changed: {recs[rid].name!r}")
        ms_obs = float(np.mean([recs[rid].estimate.value for rid in ms_ids]))
        ids.extend(ms_ids)
        res = Ms0(comp, n_draws=d.n, seed=7, content_ratio=(1.0, 1.0))
        if tc_id is None:
            tcs = np.maximum(C.tc_samples({}, d.tc_params), 1.0)             # the host T0 of the model (not measured in fit)
        else:
            if recs[tc_id].name != ETA_RECORD_NAMES[tc_id]:
                raise RuntimeError(f"fit record {tc_id} changed")
            tcs = np.full(d.n, recs[tc_id].estimate.value)
            ids.append(tc_id)
        m = C.reduced_magnetization(ETA_T_K / tcs, d.shape_s, d.form)
        e = ms_obs / (res.A_m2_per_kg * m)
        lnj[name] = np.log(e)
        q = np.percentile(e, [16, 50, 84])
        table.append({"sample": name, "Ms_obs_A_m2_kg": ms_obs, "Ms0_model_median": float(np.median(res.A_m2_per_kg)),
                      "m_300K_median": float(np.median(m)), "eta_p16": float(q[0]), "eta_median": float(q[1]),
                      "eta_p84": float(q[2]), "rec_ids": list(ms_ids) + ([tc_id] if tc_id else [])})
    L = np.stack(list(lnj.values()))                                           # (3, n)
    k = L.shape[0]
    mean = L.mean(0)
    sd = L.std(0, ddof=1)
    rng = np.random.default_rng([seed, 41])
    t2 = rng.standard_t(k - 1, size=d.n)
    ln_new = mean + sd * np.sqrt(1.0 + 1.0 / k) * t2
    return EtaResult(np.exp(ln_new), {kk: np.exp(v) for kk, v in lnj.items()}, table, float(np.median(mean)),
                     float(np.median(sd)), sorted(set(ids)))


@dataclass
class MsTResult:
    A_per_m: object                          # Dist [A/m]
    A_m2_per_kg: object                      # Dist [A m^2/kg]
    ms0: Ms0Result
    eta: np.ndarray
    m: np.ndarray                            # (nT, n) reduced magnetisation


def ms_T(T, composition, *, n: int = 4000, seed: int = 20260101, draws=None, eta: Optional[EtaResult] = None) -> MsTResult:
    """Ms(T) = eta * Ms0 * m(T / T_C) with the SHARED draws (r, Mn valence, T_C, shape s) of hexw.params.anisotropy.JointDraws.
    T scalar or array [K] -> samples (nT, n) (scalar: (n,)). Valid for x <= 2; Ms0 is the collinear Gorter value at 0 K."""
    from hexw.params import anisotropy as AN
    from hexw.params import curie as C
    d = draws if draws is not None else AN.joint_draws(n, seed)
    comp = O.clean_composition(composition)
    Tarr = np.atleast_1d(np.asarray(T, float))
    if np.any(Tarr < 0) or not np.all(np.isfinite(Tarr)):
        raise ValueError("T must be finite and >= 0 K")
    res = Ms0(comp, n_draws=d.n, seed=int(seed) % (2**31), draws=d)
    ex = eta if eta is not None else eta_draws(d)
    tcs = np.maximum(C.tc_samples(comp, d.tc_params), 1.0)
    m = C.reduced_magnetization(Tarr[:, None] / tcs[None, :], d.shape_s[None, :], d.form)
    ms = ex.eta[None, :] * res.A_per_m[None, :] * m
    sg = ex.eta[None, :] * res.A_m2_per_kg[None, :] * m
    sq = (lambda a: a[0]) if np.ndim(T) == 0 else (lambda a: a)
    notes = ("Ms(T) = eta Ms0 m(T/T_C): eta is an EFFECTIVE multiplier (covalent moment reduction, non-collinearity, VSM "
             "under-saturation) estimated from 3 fit samples; Ms0 collinear Gorter, 0 K")
    return MsTResult(C.Dist("Ms", "A/m", sq(ms), basis=notes, extra={"eta_median": float(np.median(ex.eta))}),
                     C.Dist("sigma_s", "A m^2/kg", sq(sg), basis=notes), res, ex.eta, m)


def eta_estimates(res: Optional[EtaResult] = None) -> list:
    """eta (moment-reduction multiplier) and its between-sample scatter as ParameterEstimates: role literature_fit, status
    derived_by_us (computed from the Ms values of C004 and L01 Mg/Ni at 300 K). It is an effective, composition-independent
    parameter; identifiable = True only as a fit to these three samples."""
    res = res or eta_draws()
    lo, hi = np.percentile(res.eta, [2.5, 97.5])
    prov = S.Provenance(source_id=";".join(res.rec_ids), locator="literature/fit/records_fit.yaml",
                        status=S.ProvenanceStatus.DERIVED_BY_US,
                        note="eta = Ms_obs(300 K) / (Ms0 m(300 K/T_C)) over SrFe18O27 (C004), SrMg2W, SrNi2W (L01 Mg/Ni, fit); "
                             "absorbs covalent reduction, non-collinearity and VSM under-saturation: an effective parameter")
    out = [S.ParameterEstimate("Ms_eta_moment_reduction", float(np.median(res.eta)), "1", Kind.DIMENSIONLESS,
                               S.SourceRole.LITERATURE_FIT, prov,
                               S.Uncertainty(kind="interval", lower=float(lo), upper=float(hi)), identifiable=True)]
    return out
