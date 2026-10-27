"""Ensemble of cation distributions for SrFe(18-x)A(x)O27 (TZ section 7): candidates, weights, ranges, ambiguity.

A candidate is the MEAN partially occupied description of one hypothesis (valence hypothesis x Fe2+ placement scenario):
it is flagged partial_occupancy_cif_not_dft_config=True; ordered DFT configurations come from `hexw.structures.configs`.

Constraints (checked by hexw.chem.validators): capacity of every site <= 1; charge neutrality of Sr2+, A(q+), Fe2+, Fe3+, O27:
    n(Fe2+) = 2 - sum_{A in 2+} x_A,     n(Fe3+) = 18 - x - n(Fe2+).
For A2+ only this is the TZ model n(Fe2+)=2-x, n(Fe3+)=16 (not applicable for x > 2). Mn3+ needs no Fe2+ compensation
FOR ITSELF (the extra Fe2+ relative to the Mn2+ model); it is a separate candidate. Both are defined for x <= 2 only (the
fit and the TZ model do not cover more).

Weights: priors of the HYPOTHESES (assumptions, not literature numbers) times the frozen A distribution (fit ensemble).
The A distribution reflects the cation distribution frozen at synthesis/quench (T_f), not equilibrium at the measuring T.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from hexw.chem import validators as V
from hexw.params import occupancy as O
from hexw.registry import schemas as S
from hexw.structures.sites import ENV, N_CELL, N_FU, SITES, SPIN_FE2, SPIN_FE3, SPIN_SIGN, WYCKOFF

FE2_SCENARIOS = {
    # name: (site indices that can host Fe2+, description)
    "uniform": ((0, 1, 2, 3, 4, 5, 6), "Fe2+ statistically on all Fe positions (proportional to the Fe content of the site)"),
    "oct_S_up": ((2, 3), "Fe2+ on W3 (6g) and W4 (4f_VI, S block, spin up): the placement mentioned by Nakai (L03/C004: "
                         "6g and/or '4f3' = W4 if the Fig. 1 labelling inferred in site_notation.md is right)"),
    "oct_R": ((4, 5), "Fe2+ on the octahedra of the R block / R-S boundary, W5 (12k) and W6 (4f_VI*)"),
}
FE2_PRIOR_NOTE = ("weights from the compatibility of the Gorter moment of SrFe18O27 with the DFT moment of L05 (L05-093, "
                  "sigma 1 muB/f.u., Gaussian likelihood; an estimate of the literature DFT, not a measurement); the scenarios "
                  "themselves are ASSUMPTIONS: the fit corpus does not localise Fe2+ (site_notation.md warning 6; L02 refinements "
                  "keep Fe3+ everywhere)")
L05_TOTAL_MOMENT_ID = "L05-093"            # total spin moment per cell of SrFe18O27 (DFT, 0 K), 54.991 muB / 2 f.u.
SIGMA_DFT_MOMENT = 1.0                      # ASSUMPTION: muB / f.u., covers the spheres-vs-interstitial definition of the DFT total
MN_VALENCE_PRIOR = {2: 0.5, 3: 0.5}
MN_PRIOR_NOTE = "weight not determined by the Sr-W literature (assumption: equal prior); only L11 BaMn2W (tier 3) exists"
SPACE_GROUP = "P6_3/mmc"
OX_A = {"Co": 2.0, "Zn": 2.0, "Mg": 2.0, "Ni": 2.0, "Mn": 2.0}


class EnsembleError(ValueError):
    pass


def gorter_moment_x0(scenario: str) -> float:
    """Spin-only Gorter moment (muB / f.u.) of SrFe18O27 (n(Fe2+) = 2) with Fe2+ placed by `scenario`: the number the DFT
    moment of L05 is compared with (analytic: no A cations, nothing to sample)."""
    sel = np.zeros(7, bool)
    sel[list(FE2_SCENARIOS[scenario][0])] = True
    fe2 = np.zeros(7)
    fe2[sel] = 2.0 * N_FU[sel] / N_FU[sel].sum()
    fe3 = N_FU - fe2
    return float(((fe3 * SPIN_FE3 + fe2 * SPIN_FE2) * SPIN_SIGN).sum())


def fe2_scenario_weights(sigma: float = SIGMA_DFT_MOMENT) -> dict:
    """Weights of the Fe2+ placement scenarios: Gaussian likelihood of the L05 DFT total moment (per f.u.) given the
    Gorter moment of each scenario (sigma ~ 1 muB/f.u.), normalised. Provenance L05 (an estimate of the DFT literature)."""
    from hexw.params.literature import load_fit_records
    rec = next(c for c in load_fit_records() if c.rec_id == L05_TOTAL_MOMENT_ID)
    if rec.name != "total_spin_moment_per_cell" or rec.estimate is None or rec.estimate.unit != "muB":
        raise EnsembleError(f"{L05_TOTAL_MOMENT_ID} changed: {rec.name!r}")
    dft = rec.estimate.value / 2.0
    ll = np.array([-0.5 * ((gorter_moment_x0(k) - dft) / sigma) ** 2 for k in FE2_SCENARIOS])
    w = np.exp(ll - ll.max())
    w /= w.sum()
    return dict(zip(FE2_SCENARIOS, (float(x) for x in w)))


@dataclass(frozen=True)
class Hypothesis:
    """Selection of hypotheses; None = all of them (with their prior weights)."""
    mn_valence: Optional[int] = None          # 2 | 3 (only meaningful if Mn is present)
    fe2_scenario: Optional[str] = None        # key of FE2_SCENARIOS


@dataclass
class EnsembleMember:
    candidate: S.StructureCandidate
    weight: float
    weight_basis: str
    mn_valence: Optional[int]
    fe2_scenario: str
    counts_per_fu: dict                       # {species: atoms per f.u.}
    site_counts: dict                         # {site: {species: atoms per f.u. on the site}}
    ranges: dict                              # {(site, species): (lo, hi)} 5-95 % of the A fractions over the fit ensemble
    ambiguity: list = field(default_factory=list)


@dataclass
class StructureEnsemble:
    composition: dict
    x: float
    members: list
    ambiguity: list
    A_mean_fraction: dict                     # {species: array(7)} mean site fractions of the A cations
    A_draws: np.ndarray                       # (n, S, 7) the fit-ensemble draws
    A_species: tuple
    notes: list = field(default_factory=list)

    def weights(self) -> np.ndarray:
        return np.array([m.weight for m in self.members])

    def best(self) -> EnsembleMember:
        return max(self.members, key=lambda m: m.weight)


def valence_counts(composition: dict, mn_valence: Optional[int] = None):
    """(n_Fe2, n_Fe3, oxidation dict, counts dict) of one valence hypothesis; the model is defined for x <= 2 only."""
    comp = O.clean_composition(composition)
    x = sum(comp.values())
    if x > 2 + V.TOL:
        raise V.ModelNotApplicableError(
            f"x={x} > 2: the Fe2+/Fe3+ model with Sr2+, A2+, O27 (n(Fe2+)=2-x<0) is not applicable; this module is "
            "defined for x <= 2 (fit corpus and TZ section 7)")
    q = {sp: OX_A[sp] for sp in comp}
    if "Mn" in comp:
        if mn_valence not in (2, 3):
            raise EnsembleError("Mn present: mn_valence must be 2 or 3")
        q["Mn"] = float(mn_valence)
    n2 = 2.0 - sum(v for sp, v in comp.items() if q[sp] == 2.0)
    if n2 < -V.TOL:
        raise V.ModelNotApplicableError(f"n(Fe2+) = {n2:.3f} < 0 for composition {comp}")
    n2 = max(n2, 0.0)
    n3 = 18.0 - x - n2
    counts = {"Sr": 1.0, **comp, "Fe2+": n2, "Fe3+": n3, "O": 27.0}
    ox = {"Sr": 2.0, **q, "Fe2+": 2.0, "Fe3+": 3.0, "O": -2.0}
    V.check_charge_neutral(counts, ox)
    return n2, n3, ox, counts


def _formula(comp: dict) -> str:
    x = sum(comp.values())
    parts = "".join(f"{sp}{v:g}" for sp, v in sorted(comp.items()))
    return f"Sr{parts}Fe{18 - x:g}O27"


def candidate_counts(c: S.StructureCandidate, z: float = 2.0) -> dict:
    """Atoms per formula unit by species, from the site occupancies (multiplicities refer to a cell with z formula units;
    z = 2 for the hexagonal cell)."""
    out = {"Sr": 1.0, "O": 27.0}
    for s in c.sites:
        for sp, f in s.species.items():
            out[sp] = out.get(sp, 0.0) + f * s.multiplicity / z
    return out


def check_candidate(c: S.StructureCandidate, z: float = 2.0) -> None:
    """Capacity (every site <= 1) and charge neutrality (Sr2+, A, Fe2+, Fe3+, O27) of a candidate."""
    V.check_site_occupancy(c)
    counts = candidate_counts(c, z)
    V.check_charge_neutral(counts, {"Sr": 2.0, "O": -2.0, **c.oxidation_states})


def _fe2_allocation(fe_site: np.ndarray, n2: float, sites: tuple) -> np.ndarray:
    """Fe2+ atoms per f.u. on each site: proportional to the Fe content inside the scenario's site set."""
    sel = np.zeros(7, bool)
    sel[list(sites)] = True
    cap = fe_site[sel].sum()
    if n2 > cap + 1e-12:
        raise EnsembleError(f"scenario cannot host {n2:.3f} Fe2+ (Fe available {cap:.3f})")
    out = np.zeros(7)
    if n2 > 0:
        out[sel] = n2 * fe_site[sel] / cap
    return out


def build_ensemble(composition: dict, hypothesis: Optional[Hypothesis] = None, fit: Optional[O.OccupancyFit] = None,
                   n_draws: Optional[int] = None, seed: int = 0) -> StructureEnsemble:
    """Candidates for a composition: valence hypotheses x Fe2+ placement scenarios, each with the mean frozen A
    distribution, 5-95 % ranges and reasons of ambiguity."""
    hyp = hypothesis or Hypothesis()
    comp = O.clean_composition(composition)
    x = sum(comp.values())
    fit = O.default_fit() if fit is None else fit
    species, draws = O.site_fraction_draws(comp, fit, n_draws, seed)        # (n, S, 7), fractions of the site
    pmean = draws.mean(0)
    lo, hi = np.percentile(draws, [5, 95], axis=0) if species else (np.zeros((0, 7)), np.zeros((0, 7)))

    has_mn = "Mn" in comp
    vals = [hyp.mn_valence] if hyp.mn_valence is not None else ([2, 3] if has_mn else [None])
    if has_mn and any(v not in (2, 3) for v in vals):
        raise EnsembleError("mn_valence must be 2 or 3")
    if not has_mn and hyp.mn_valence is not None:
        raise EnsembleError("mn_valence given, but the composition has no Mn")
    scen = [hyp.fe2_scenario] if hyp.fe2_scenario else list(FE2_SCENARIOS)
    for s in scen:
        if s not in FE2_SCENARIOS:
            raise EnsembleError(f"unknown Fe2+ scenario {s!r}; allowed {sorted(FE2_SCENARIOS)}")

    amb_common = _ambiguity_common(comp, fit, species)
    fe2w_all = fe2_scenario_weights()
    tot_w = sum(fe2w_all[s_] for s_ in scen)
    fe2w = {s_: fe2w_all[s_] / tot_w for s_ in scen}
    members = []
    for v in vals:
        n2, n3, ox, counts = valence_counts(comp, v)
        for sc in scen:
            fe_site = N_FU * (1.0 - pmean.sum(0))                    # Fe atoms per f.u. on each site
            try:
                fe2 = _fe2_allocation(fe_site, n2, FE2_SCENARIOS[sc][0])
            except EnsembleError:
                continue                                              # scenario physically impossible here: dropped
            fe3 = fe_site - fe2
            sites, site_counts = [], {}
            for i, sid in enumerate(SITES):
                sp = {}
                for j, a in enumerate(species):
                    if pmean[j, i] > 1e-12:
                        sp[a] = float(pmean[j, i])
                if fe3[i] > 1e-12:
                    sp["Fe3+"] = float(fe3[i] / N_FU[i])
                if fe2[i] > 1e-12:
                    sp["Fe2+"] = float(fe2[i] / N_FU[i])
                sites.append(S.SiteOccupancy(site=sid, species=sp, multiplicity=int(N_CELL[i]), wyckoff=WYCKOFF[i],
                                             site_id=sid, space_group=SPACE_GROUP))
                site_counts[sid] = {k: f * N_FU[i] for k, f in sp.items()}
            ox_c = {k: val for k, val in ox.items() if k not in ("Sr", "O")}
            w_val = MN_VALENCE_PRIOR[v] / sum(MN_VALENCE_PRIOR[u] for u in vals) if has_mn else 1.0
            w_sc = fe2w[sc]
            basis = f"weight(Fe2+ scenario)={w_sc:.3g} [{FE2_PRIOR_NOTE}]"
            if has_mn:
                basis = f"prior(Mn{v}+)={w_val:.3g} [{MN_PRIOR_NOTE}]; " + basis
            cand = S.StructureCandidate(
                candidate_id=f"{_formula(comp)}|Mn{v}+|{sc}" if has_mn else f"{_formula(comp)}|{sc}",
                formula=_formula(comp), sites=sites, oxidation_states=ox_c,
                magnetic_hypothesis="collinear Gorter ferrimagnet (W1,W2,W6 down; W3,W4,W5,W7 up)",
                source="hexw.structures.ensemble: site-preference ensemble fitted to fit literature (L02 Co/Zn, L01 Mg/Ni)",
                status="candidate", partial_occupancy_cif_not_dft_config=True)
            check_candidate(cand)
            rng_ = {(SITES[i], a): (float(lo[j, i]), float(hi[j, i])) for j, a in enumerate(species) for i in range(7)
                    if hi[j, i] > 1e-6}
            amb = list(amb_common)
            amb.append(f"Fe2+ placement ({sc}): {FE2_SCENARIOS[sc][1]}; {FE2_PRIOR_NOTE}")
            if has_mn:
                amb.append(f"Mn valence hypothesis Mn{v}+; {MN_PRIOR_NOTE}")
            members.append(EnsembleMember(cand, w_val * w_sc, basis, v, sc,
                                          {k: float(vv) for k, vv in counts.items()}, site_counts, rng_, amb))
    if not members:
        raise EnsembleError("no admissible hypothesis for this composition")
    tot = sum(m.weight for m in members)
    for m in members:
        m.weight /= tot
    return StructureEnsemble(comp, x, members, amb_common, {a: pmean[j].copy() for j, a in enumerate(species)},
                             draws, species, notes=[
                                 "weights describe the distribution FROZEN at synthesis/quench (T_f), not at the "
                                 "measuring temperature",
                                 "fractions refer to the site; W7 = the 2d bipyramid (the CIF has a half-occupied split 4f)"])


def _ambiguity_common(comp: dict, fit: O.OccupancyFit, species: tuple) -> list:
    out = [("A distribution comes from a near-saturated fit: %d fitted logits, %d identified relative energies, "
            "chi2/dof = %.1f/%d (bootstrap ensemble of %d fits)" % (fit.n_fitted, fit.n_identified, fit.chi2, fit.dof,
                                                                       fit.boot_a.shape[0]))]
    for sp in species:
        if sp == "Mn":
            out.append("Mn: NO site-occupancy data in the fit corpus; uninformative prior (assumption), identifiable=False")
            continue
        k = fit.index(sp)
        bad = [SITES[i] for i in range(7) if fit.has_site[k, i] and not fit.identified[k, i] and i != O.REF_SITE]
        if bad:
            out.append(f"{sp}: sites {', '.join(bad)} not identified (fraction consistent with 0): only upper limits")
        if not fit.has_site[k, 6]:
            out.append(f"{sp}: no datum for W7, taken as unoccupied (absent in all fit samples of other species)")
    if "Co" in comp:
        out.append("Co on W3 (6g): the only source is L02; there is no cross-check in the fit")
    out.append("refined A content in the fit samples is below nominal (Co 1.71, Zn 1.80, Mg 1.69 per f.u.); the "
               "prediction uses the nominal content x")
    return out
