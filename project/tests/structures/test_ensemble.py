"""Ensemble of cation distributions: candidates, constraints (capacity, charge), Mn hypotheses, weights, ambiguity."""
import numpy as np
import pytest

from hexw.chem import validators as V
from hexw.params import occupancy as O
from hexw.registry import schemas as S
from hexw.structures import ensemble as E
from hexw.structures.sites import N_FU

CO_ZN_SERIES = [{"Co": 2.0 - x, "Zn": x} for x in (0.0, 0.5, 1.0, 1.5, 2.0)]
SERIES = [{a: x} for a in ("Co", "Zn", "Mn") for x in (0.5, 1.0, 1.5, 2.0)]


@pytest.fixture(scope="module", autouse=True)
def _warm():
    O.default_fit()


@pytest.mark.parametrize("comp", [c for c in CO_ZN_SERIES if sum(c.values()) > 0] + SERIES,
                         ids=lambda c: "_".join(f"{k}{v:g}" for k, v in c.items()))
def test_candidates_returned_and_validated(comp):
    ens = E.build_ensemble(comp, n_draws=60)
    assert ens.members and abs(ens.weights().sum() - 1.0) < 1e-12
    has_mn = "Mn" in comp
    assert len(ens.members) == len(E.FE2_SCENARIOS) * (2 if has_mn else 1)
    x = sum(comp.values())
    for m in ens.members:
        c = m.candidate
        assert c.partial_occupancy_cif_not_dft_config is True             # the mean candidate is NOT a DFT configuration
        with pytest.raises(S.PartialOccupancyError):
            c.require_dft_configuration()
        V.check_site_occupancy(c)
        counts = E.candidate_counts(c)
        V.check_charge_neutral(counts, {"Sr": 2.0, "O": -2.0, **c.oxidation_states})
        fe2, fe3 = counts.get("Fe2+", 0.0), counts.get("Fe3+", 0.0)
        assert abs(fe2 + fe3 + x - 18.0) < 1e-9                           # all 18 TM sites are full
        for sp, n in comp.items():
            if n > 0:
                assert abs(counts[sp] - n) < 1e-9
        if m.mn_valence == 3:
            assert abs(fe2 - 2.0) < 1e-9                                  # Mn3+: the extra Fe2+ of the hypothesis
        else:
            assert abs(fe2 - (2.0 - x)) < 1e-9 and abs(fe3 - 16.0) < 1e-9     # TZ section 7
        assert all(s.site_id in S.SITE_IDS and s.multiplicity for s in c.sites)
        assert m.ambiguity and any("Fe2+ placement" in a for a in m.ambiguity)
        assert m.ranges and all(lo <= hi for lo, hi in m.ranges.values())


def test_capacity_and_content_hold_for_every_draw():
    for comp in ({"Co": 2.0}, {"Zn": 2.0}, {"Co": 1.0, "Zn": 1.0}, {"Co": 0.5, "Mn": 1.5}):
        ens = E.build_ensemble(comp, n_draws=80)
        assert np.all(ens.A_draws.sum(1) <= 1 + 1e-12)
        assert np.allclose((ens.A_draws * N_FU).sum(-1), [comp[s] for s in ens.A_species], atol=1e-7)


def test_zn_goes_to_tetrahedra_and_co_to_6g():
    zn = E.build_ensemble({"Zn": 2.0}, n_draws=60).A_mean_fraction["Zn"]
    co = E.build_ensemble({"Co": 2.0}, n_draws=60).A_mean_fraction["Co"]
    assert zn[0] > 0.4 and zn[0] > zn[1] > zn[2] > zn[3] and zn[6] < 0.02    # W1 first, W7 empty
    assert np.argmax(co) == 2 and co[2] > 0.4                                # Co on 6g (W3), L02 only


def test_mn_hypotheses_are_separate_candidates_with_prior_weights():
    ens = E.build_ensemble({"Mn": 1.0}, n_draws=60)
    by_val = {}
    for m in ens.members:
        by_val.setdefault(m.mn_valence, []).append(m)
    assert set(by_val) == {2, 3}
    w2, w3 = sum(m.weight for m in by_val[2]), sum(m.weight for m in by_val[3])
    assert abs(w2 - 0.5) < 1e-12 and abs(w3 - 0.5) < 1e-12
    for m in ens.members:
        assert "not determined by the Sr-W literature" in m.weight_basis
        assert any("Mn" in a and "NO site-occupancy data" in a for a in m.ambiguity)
        assert m.candidate.oxidation_states["Mn"] == float(m.mn_valence)
    ids = [m.candidate.candidate_id for m in ens.members]
    assert len(set(ids)) == len(ids)
    # a hypothesis can be selected
    only3 = E.build_ensemble({"Mn": 1.0}, E.Hypothesis(mn_valence=3, fe2_scenario="oct_R"), n_draws=40)
    assert len(only3.members) == 1 and only3.members[0].weight == 1.0
    with pytest.raises(E.EnsembleError):
        E.build_ensemble({"Co": 1.0}, E.Hypothesis(mn_valence=3))
    with pytest.raises(E.EnsembleError):
        E.valence_counts({"Mn": 1.0}, None)


def test_mn3_adds_fe2_relative_to_mn2():
    n2_a, n3_a, *_ = E.valence_counts({"Mn": 1.0}, 2)
    n2_b, n3_b, *_ = E.valence_counts({"Mn": 1.0}, 3)
    assert (n2_a, n3_a) == (1.0, 16.0) and (n2_b, n3_b) == (2.0, 15.0)
    # Co + Mn3+: only the Co2+ takes part in n(Fe2+) = 2 - x_2+
    n2_c, n3_c, ox, counts = E.valence_counts({"Co": 0.5, "Mn": 1.0}, 3)
    assert abs(n2_c - 1.5) < 1e-12 and abs(n3_c - 15.0) < 1e-12 and ox["Mn"] == 3.0


def test_model_not_applicable_beyond_x2_and_unknown_species():
    with pytest.raises(V.ModelNotApplicableError):
        E.build_ensemble({"Co": 2.5})
    with pytest.raises(V.ModelNotApplicableError):
        E.valence_counts({"Mn": 2.5}, 3)
    with pytest.raises(ValueError):
        E.build_ensemble({"Cu": 1.0})


def test_end_member_without_a_cation():
    ens = E.build_ensemble({}, n_draws=10)
    assert len(ens.members) == 3
    for m in ens.members:
        c = E.candidate_counts(m.candidate)
        assert abs(c["Fe2+"] - 2.0) < 1e-9 and abs(c["Fe3+"] - 16.0) < 1e-9


def test_fe2_scenarios_place_fe2_on_their_sites():
    ens = E.build_ensemble({"Co": 1.0}, n_draws=40)
    for m in ens.members:
        sites_with = {s.site for s in m.candidate.sites if s.species.get("Fe2+", 0) > 0}
        allowed = {f"W{i + 1}" for i in E.FE2_SCENARIOS[m.fe2_scenario][0]}
        assert sites_with and sites_with <= allowed


def test_fit_check_mg_ni_mean_matches_observed_fractions():
    """Check on the fit corpus: the ensemble at the refined content reproduces the L01 Mg / Ni site fractions."""
    for sp, n in (("Mg", 1.689), ("Ni", 2.0)):
        ens = E.build_ensemble({sp: n}, n_draws=100)
        obs = next(s for s in O.load_occupancy_samples() if sp in s.species).species[sp]
        # the observed totals differ slightly from the requested content (Ni refined 2.113 > 2), so compare shapes
        mean = ens.A_mean_fraction[sp]
        for i in range(6):
            tol = 3 * obs.sigma_eff()[i] + (0.02 if sp == "Ni" else 0.0)
            assert abs(mean[i] - obs.p[i]) < tol, (sp, i, mean[i], obs.p[i])


def test_schema_roundtrip_and_ranges():
    ens = E.build_ensemble({"Co": 1.0, "Zn": 1.0}, n_draws=40)
    m = ens.best()
    plain = S.dump_json(m.candidate)
    assert "partial_occupancy_cif_not_dft_config" in plain
    # 5-95 % range brackets the mean for the dominant occupancies
    lo, hi = m.ranges[("W1", "Zn")]
    assert lo <= ens.A_mean_fraction["Zn"][0] <= hi
