"""Site-preference model of the A cations: solver invariants, synthetic recovery, fit parameters, provenance. Fit data are the
fit-literature records only (L02 Co/Zn/CoZn, L01 Mg/Ni)."""
import subprocess
import sys

import numpy as np
import pytest

from hexw.params import occupancy as O
from hexw.params.literature import load_fit_records
from hexw.registry import schemas as S
from hexw.structures.sites import N_FU


# ---------------------------------------------------------------- solver: capacity and content are exact
def test_solver_conserves_content_and_capacity():
    rng = np.random.default_rng(1)
    a = rng.uniform(-6, 3, size=(50, 3, 7))
    n = rng.uniform(0.1, 2.0, size=(50, 3))
    n = n / np.maximum(1.0, n.sum(1, keepdims=True) / 4.0)           # at most 4 A atoms per f.u. in total
    p = O.solve_site_fractions(a, n)
    assert np.allclose((p * N_FU).sum(-1), n, atol=1e-8)             # content: the Lagrange constraint is met
    assert np.all(p >= 0) and np.all(p.sum(-2) <= 1 + 1e-12)         # capacity of every site


def test_solver_rejects_overfull_and_nonpositive():
    a = np.full((1, 7), -3.0)
    with pytest.raises(ValueError):
        O.solve_site_fractions(a, np.array([19.0]))                  # more atoms than the 18 sites
    with pytest.raises(ValueError):
        O.solve_site_fractions(a, np.array([0.0]))


def test_solver_single_site_limit():
    # one site type reachable (others absent): p = n / N_i
    a = np.full((1, 7), O.A_ABSENT)
    a[0, 4] = 0.0
    p = O.solve_site_fractions(a, np.array([3.0]))
    assert abs(p[0, 4] - 3.0 / 6.0) < 1e-8 and np.allclose(np.delete(p[0], 4), 0, atol=1e-12)


# ---------------------------------------------------------------- synthetic recovery
def _synthetic_samples(d_co, d_zn, rng, sig=0.004):
    def obs(p, sg=sig):
        p = np.clip(p + rng.normal(size=7) * sg, 0, 1)
        return O.SpeciesObs(p, np.full(7, sg), [()] * 7)
    a_co, a_zn = -np.asarray(d_co, float), -np.asarray(d_zn, float)
    p_co = O.solve_site_fractions(a_co[None], np.array([1.7]))[0]
    p_zn = O.solve_site_fractions(a_zn[None], np.array([1.8]))[0]
    mix = O.solve_site_fractions(np.stack([a_co, a_zn]), np.array([0.85, 0.8]))
    return [O.OccSample("pureCo", "SYN", {"Co": obs(p_co)}), O.OccSample("pureZn", "SYN", {"Zn": obs(p_zn)}),
            O.OccSample("mix", "SYN", {"Co": obs(mix[0]), "Zn": obs(mix[1])})]


def test_synthetic_parameter_recovery():
    d_co = np.array([2.5, 4.0, 0.0, 3.5, 5.0, 8.0, 8.0])
    d_zn = np.array([-2.5, -0.6, 0.0, 6.0, 6.0, 1.2, 6.0])
    rng = np.random.default_rng(7)
    fit = O.fit_occupancy(_synthetic_samples(d_co, d_zn, rng), n_boot=40, seed=3)
    assert fit.dof > 0
    for sp, truth in (("Co", d_co), ("Zn", d_zn)):
        k = fit.index(sp)
        for i in range(7):
            if i == O.REF_SITE:
                continue
            if fit.identified[k, i]:
                sd = fit.boot_d[:, k, i].std(ddof=1)
                assert abs(fit.d_hat[k, i] - truth[i]) < 3 * sd + 0.2, (sp, i)
            else:
                # not identified: only a lower bound, and the true d is above it (the data allow it)
                assert fit.d_lower[k, i] <= truth[i] + 1e-6, (sp, i)
    # sites with a clear signal are identified at the smallest sigma_sys, sites with p ~ 0 are not; robust => at the minimum too
    assert fit.identified_at_min[fit.index("Co"), 0] and fit.identified_at_min[fit.index("Zn"), 0]
    assert not fit.identified[fit.index("Co"), 5] and not fit.identified[fit.index("Zn"), 3]
    assert not np.any(fit.identified & ~fit.identified_at_min)
    assert fit.n_identified <= fit.n_identified_at_min


def test_synthetic_fit_is_exact_without_noise():
    d_co = np.array([2.0, 3.0, 0.0, 3.0, 4.0, 5.0, 6.0])
    d_zn = np.array([-2.0, -0.5, 0.0, 5.0, 5.0, 1.5, 5.0])
    rng = np.random.default_rng(0)
    s = _synthetic_samples(d_co, d_zn, rng, sig=0.0)
    fit = O.fit_occupancy(s, n_boot=2)
    k = fit.index("Co")
    assert np.allclose(fit.d_hat[k, [0, 1, 3]], d_co[[0, 1, 3]], atol=5e-3)
    assert fit.chi2 < 1e-3


# ---------------------------------------------------------------- the real fit-literature data
@pytest.fixture(scope="module")
def fit():
    return O.default_fit()


def test_samples_come_only_from_the_fit_literature():
    samples = O.load_occupancy_samples()
    names = {s.name for s in samples}
    assert names == {"SrCo2Fe16O27", "SrCoZnFe16O27", "SrZn2Fe16O27", "SrMg2Fe16O27", "SrNi2Fe16O27"}
    fit_ids = {r.rec_id for r in load_fit_records()}
    used = {rid for s in samples for o in s.species.values() for ids in o.rec_ids for rid in ids}
    assert used and used <= fit_ids
    # L01 contributes Mg and Ni only
    assert {sp for s in samples if s.lit_id == "L01" for sp in s.species} == {"Mg", "Ni"}


def test_weighted_cmcm_sites_reproduce_the_derived_records():
    byid = {r.rec_id: r for r in load_fit_records()}
    samples = {s.name: s for s in O.load_occupancy_samples()}
    for rid, sample, sp, site in (("L02-098", "SrCo2Fe16O27", "Co", 2), ("L02-099", "SrCo2Fe16O27", "Co", 4),
                                  ("L02-100", "SrCoZnFe16O27", "Co", 2)):
        assert abs(samples[sample].species[sp].p[site] - byid[rid].estimate.value) < 1.5e-3, rid


def test_fit_quality_and_parameter_count(fit):
    n_data = len(fit.residuals)
    assert fit.dof > 0 and fit.n_fitted < n_data                      # not more free logits than data
    assert fit.n_identified < fit.n_fitted
    assert fit.chi2 / fit.dof < 3.0
    assert fit.boot_a.shape[0] == O.N_BOOT


def test_parameter_estimates_are_valid_with_provenance_and_units(fit):
    ids = {r.rec_id for r in load_fit_records()}
    est = fit.estimates
    assert est and all(isinstance(e, S.ParameterEstimate) for e in est)
    for e in est:
        assert e.unit in ("1", "K") and e.provenance.source_id
        # an assumption is an assumption (T_f, sigma_sys); the fitted energies are literature_fit, derived by us
        assert ("ASSUMPTION" in e.provenance.source_id) == (e.source_role is S.SourceRole.ASSUMPTION)
        if e.name.startswith("occ_dE_"):
            assert e.source_role is S.SourceRole.LITERATURE_FIT and e.provenance.status is S.ProvenanceStatus.DERIVED_BY_US
        if e.name.startswith("occ_dE_"):
            cited = e.provenance.source_id.split(";")
            assert set(cited) <= ids, e.name
        if e.identifiable:
            assert e.uncertainty.kind == "std" and e.uncertainty.std > 0
        else:
            hw = e.uncertainty.half_width()
            assert e.uncertainty.kind in ("unknown", "interval")
            assert e.uncertainty.kind == "unknown" or e.uncertainty.declared_wide or hw / abs(e.value) >= 0.5


def test_non_identified_parameters_have_wide_intervals_and_are_flagged(fit):
    non = [e for e in fit.estimates if not e.identifiable and e.name.startswith("occ_dE_")]
    assert non
    for e in non:
        assert e.uncertainty.inequality == "ge"                       # a bound, not a value
        assert e.uncertainty.upper - e.uncertainty.lower >= 8.0 - 1e-9
        assert "NOT robustly identified" in e.provenance.note
    by = {e.name: e for e in fit.estimates}
    assert not by["occ_dE_Co_W7"].identifiable and not by["occ_dE_Zn_W4"].identifiable
    assert by["occ_dE_Zn_W1"].identifiable
    # Co on W1 is identified at sigma_sys = 0.01 but NOT over the whole range up to 0.03: the honest (robust) status is False
    assert fit.identified_at_min[fit.index("Co"), 0] and not by["occ_dE_Co_W1"].identifiable


def test_physical_sign_pattern(fit):
    # Zn prefers the tetrahedron W1 over W3; Co prefers W3 over W1
    z, c = fit.index("Zn"), fit.index("Co")
    assert fit.d_hat[z, 0] < 0 and fit.d_hat[z, 1] < fit.d_hat[c, 1]
    assert fit.d_hat[c, 0] > 1.0


def test_freezing_temperature_is_a_non_identifiable_assumption():
    t = O.freezing_temperature_estimate()
    assert not t.identifiable and t.provenance.status == "assumption" and t.unit == "K"
    assert t.uncertainty.upper == 1523.0 and t.uncertainty.lower < 600
    assert abs(O.energy_eV(1.0, 1000.0) - 8.617333262e-2) < 1e-9


def test_mn_prior_is_an_unidentifiable_assumption(fit):
    est = O.mn_estimates(fit)
    assert len(est) == 6
    for e in est:
        assert not e.identifiable and e.provenance.status == "assumption"
        assert "Sr-W literature" in e.provenance.note
        assert e.uncertainty.kind == "interval" and e.uncertainty.declared_wide
    a = O.sample_mn_a(fit, 500, np.random.default_rng(0))
    assert a.shape == (500, 7) and np.all(a[:, 2] == 0)               # reference site W3
    assert np.all(-a[:, 6] >= O.mn_prior(fit).w7_lower - 1e-9)        # W7: never found for other cations


def test_site_fraction_draws_content_and_capacity(fit):
    for comp in ({"Co": 1.5}, {"Zn": 2.0}, {"Co": 1.0, "Zn": 1.0}, {"Mg": 1.0}, {"Co": 0.5, "Mn": 1.0}):
        sp, p = O.site_fraction_draws(comp, fit, 60)
        assert sp == tuple(s for s in O.SUPPORTED_A if s in comp)
        assert np.allclose((p * N_FU).sum(-1), [comp[s] for s in sp], atol=1e-7)
        assert np.all(p.sum(1) <= 1 + 1e-12)


def test_clean_composition_rejects_garbage():
    for bad in ({"Cu": 1.0}, {"Co": -1.0}, {"Co": float("nan")}, {"Co": True}, [("Co", 1)]):
        with pytest.raises(ValueError):
            O.clean_composition(bad)
    assert O.clean_composition({"Co": 0.0, "Zn": 1.0}) == {"Zn": 1.0}


def test_cozn_is_predicted_from_the_pure_samples():
    """Out-of-sample check inside the fit corpus: fit Co and Zn WITHOUT the CoZn sample, predict its site fractions."""
    ss = O.load_occupancy_samples()
    pure = [s for s in ss if s.name != "SrCoZnFe16O27"]
    f = O.fit_occupancy(pure, n_boot=2)
    cz = next(s for s in ss if s.name == "SrCoZnFe16O27")
    sp = ["Co", "Zn"]
    p = O.solve_site_fractions(f.a_hat[[f.index(x) for x in sp]], np.array([cz.species[x].n() for x in sp]))
    z = np.array([(p[j] - cz.species[x].p) / cz.species[x].sigma_eff() for j, x in enumerate(sp)])
    assert np.all(np.abs(z) < 4.5)                                    # worst: Zn on W2 (predicted 0.05, observed 0.10)
    assert abs(z[0, 2]) < 1.5 and abs(z[1, 0]) < 1.5                  # Co on W3 and Zn on W1: the dominant occupancies


def test_import_does_not_start_pools_or_fit():
    code = ("import sys\nimport hexw.params.occupancy as o, hexw.params.moments, hexw.structures.ensemble, "
            "hexw.structures.configs\n"
            "assert 'hexw.util.pool' not in sys.modules\n"
            "assert o.default_fit.cache_info().currsize == 0\n"
            "assert 'pymatgen' not in sys.modules\n")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
