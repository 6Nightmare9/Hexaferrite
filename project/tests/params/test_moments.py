"""Ms0 in the Gorter model: scheme vs the fit literature, analytic limits, units, intervals, hypotheses, sanity checks."""
import numpy as np
import pytest

from hexw.params import moments as M
from hexw.params import occupancy as O
from hexw.params.literature import load_fit_records
from hexw.registry import schemas as S
from hexw.structures import ensemble as E
from hexw.structures.sites import N_FU, SITES, SPIN_SIGN

N = 200


@pytest.fixture(scope="module", autouse=True)
def _warm():
    O.default_fit()


# ---------------------------------------------------------------- the ferrimagnetic scheme is the one in the fit data
def _sign(x):
    return 1.0 if x > 0 else -1.0


def test_gorter_scheme_matches_the_fit_moments():
    recs = {r.rec_id: r for r in load_fit_records() if r.estimate is not None}
    by = {}
    for r in recs.values():
        by.setdefault((r.name, r.source.sample_id_in_source), r.estimate.value)
    for sample in ("SrMg2Fe16O27", "SrNi2Fe16O27"):                    # L01 Fig. 4: + parallel to the net moment
        for i, s in enumerate(SITES):
            assert _sign(by[(f"site_magnetic_moment_{s}", sample)]) == SPIN_SIGN[i], (sample, s)
    for i, s in enumerate(SITES):                                      # L05 DFT, collinear Gorter
        if s == "W2":
            continue                                                   # (W2 DFT row not in the table of the same name)
        v = next((r.estimate.value for r in recs.values() if r.name.startswith(f"Fe_spin_moment_{s}_")
                  and r.source.sample_id_in_source == "SrW"), None)
        if v is not None:
            assert _sign(v) == SPIN_SIGN[i], s
    # L02 SrZn2W (hexagonal axes, collinear): W4's printed sign is the documented extraction error (CONFLICTS B.3)
    for i, s in enumerate(SITES):
        if s == "W4":
            continue
        v = next(r.estimate.value for r in recs.values() if r.name.startswith(f"Fe_moment_Mz_{s}_")
                 and r.source.sample_id_in_source == "SrZn2Fe16O27")
        assert _sign(v) == SPIN_SIGN[i], s
    assert (N_FU * (SPIN_SIGN > 0)).sum() == 12 and (N_FU * (SPIN_SIGN < 0)).sum() == 6      # 12 up, 6 down per f.u.


# ---------------------------------------------------------------- analytic limits
def test_srfe18o27_spin_only_gorter_values():
    for sc, want in (("oct_S_up", 28.0), ("uniform", 30.0 - 2 * (12 - 6) / 18), ("oct_R", 29.0)):
        r = M.Ms0({}, E.Hypothesis(fe2_scenario=sc), n_draws=N)
        # 12 sites up x 5 - 6 down x 5 = 30; each Fe2+ (4 muB) on an up site: -1, on a down site: +1
        assert np.allclose(r.spin_only_mu_B_per_fu, want, atol=1e-9), sc
        assert np.all(r.mu_B_per_fu >= r.spin_only_mu_B_per_fu - 1e-12)            # orbital parts are >= 0 here


def test_unit_conversions():
    r = M.Ms0({}, E.Hypothesis(fe2_scenario="oct_S_up"), n_draws=N)
    m_fu = 87.62 + 18 * 55.845 + 27 * 15.999
    assert abs(r.mass_g_per_mol - m_fu) < 1e-9
    one_muB = 9.2740100783e-24 * 6.02214076e23 / (m_fu * 1e-3)           # A m^2/kg per muB/f.u.
    assert abs(one_muB - 3.659) < 0.01
    assert np.allclose(r.A_m2_per_kg / r.mu_B_per_fu, one_muB, rtol=1e-9)
    # A/m from the lattice of L02/L01: V_fu ~ 496 A^3
    v, s = M.volume_fu_model({})
    assert 490 < v < 502 and 0 < s < 0.02 * v
    assert np.allclose(r.A_per_m * r.volume_fu_A3 * 1e-30 / 9.2740100783e-24, r.mu_B_per_fu, rtol=1e-6)
    assert 4e5 < np.median(r.A_per_m) < 7e5


def test_volume_uses_fit_lattices_only():
    pts = M._lattice_points()
    assert set(pts) == {"SrMg2Fe16O27", "SrNi2Fe16O27", "SrCo2Fe16O27", "SrCoZnFe16O27", "SrZn2Fe16O27"}
    a, c = pts["SrZn2Fe16O27"][:2]
    assert abs(M.volume_fu_model({"Zn": 2.0})[0] - M.cell_volume_A3(a, c) / 2) < 1e-9


# ---------------------------------------------------------------- Ms0 for the series
SERIES = [{"Co": 2.0 - x, "Zn": x} for x in (0.0, 0.5, 1.0, 1.5, 2.0) if 2.0 - x > 0 or x > 0] + \
         [{a: x} for a in ("Co", "Zn", "Mn") for x in (0.5, 1.0, 1.5, 2.0)]


@pytest.mark.parametrize("comp", SERIES, ids=lambda c: "_".join(f"{k}{v:g}" for k, v in c.items()))
def test_ms0_has_ordered_finite_intervals(comp):
    comp = {k: v for k, v in comp.items() if v > 0}
    s = M.Ms0(comp, n_draws=N).summary()
    for key in ("muB_per_fu", "A*m^2/kg", "A/m"):
        q = s[key]
        assert q["p2.5"] <= q["p16"] <= q["median"] <= q["p84"] <= q["p97.5"]
        assert q["std"] > 0 and np.isfinite(q["mean"]) and q["median"] > 0
    assert 15 < s["muB_per_fu"]["median"] < 45


def test_physical_trends():
    med = lambda c: M.Ms0(c, n_draws=N).summary()["muB_per_fu"]["median"]      # noqa: E731
    zn = [med({"Zn": x}) for x in (0.5, 1.0, 1.5, 2.0)]
    assert zn == sorted(zn) and zn[-1] > med({}) + 5                            # Zn on the down tetrahedra raises Ms0
    co = [med({"Co": x}) for x in (0.5, 1.0, 2.0)]
    assert co == sorted(co, reverse=True) and co[-1] < med({})                  # Co on up 6g (3 muB instead of 5) lowers it
    # the Co+Zn series lies between the Zn and the Co end member
    assert med({"Co": 2.0}) < med({"Co": 1.0, "Zn": 1.0}) < med({"Zn": 2.0})


def test_co_orbital_range_widens_the_interval():
    r = M.Ms0({"Co": 2.0}, E.Hypothesis(fe2_scenario="oct_S_up"), n_draws=400)
    shift = r.mu_B_per_fu - r.spin_only_mu_B_per_fu
    assert shift.min() >= 0 and 1.0 < shift.max() < 4.0 and shift.std() > 0.1   # Co L in U(0,1) x ~2.9 Co on up sites


def test_mn_hypotheses_and_wide_interval():
    r = M.Ms0({"Mn": 2.0}, n_draws=600)
    by = r.by_hypothesis()
    assert {k[0] for k in by} == {"Mn2+", "Mn3+"}
    med = {v: np.median(r.mu_B_per_fu[r.mn_valence == v]) for v in (2, 3)}
    assert med[2] != med[3]
    # the pure Mn3+ hypothesis has 4 muB ions and 2 Fe2+: distinct from Mn2+
    only3 = M.Ms0({"Mn": 2.0}, E.Hypothesis(mn_valence=3), n_draws=N)
    assert np.all(only3.mn_valence == 3)
    # Mn2+ (5 muB) is magnetically identical to Fe3+: at x = 2 (no Fe2+) the Mn2+ moment does not depend on the
    # (unknown) Mn distribution at all, so the uncertainty of Mn comes from the valence hypothesis (Mn3+: 4 muB + 2 Fe2+)
    m2 = M.Ms0({"Mn": 2.0}, E.Hypothesis(mn_valence=2), n_draws=N, content_ratio=(1.0, 1.0))
    assert np.allclose(m2.spin_only_mu_B_per_fu, 30.0)
    w = lambda res: np.subtract(*np.percentile(res.mu_B_per_fu, [97.5, 2.5]))     # noqa: E731
    assert w(r) > w(M.Ms0({"Zn": 2.0}, n_draws=N)) and abs(med[2] - med[3]) > 1.0
    with pytest.raises(ValueError):
        M.Ms0({"Co": 1.0}, E.Hypothesis(mn_valence=2))


def test_content_ratio_models_the_deficit_of_the_a_cation():
    hyp = E.Hypothesis(fe2_scenario="oct_S_up")                      # one scenario: the weighted mixture of three is bimodal (noisy std)
    exact = M.Ms0({"Zn": 2.0}, hyp, n_draws=N, content_ratio=(1.0, 1.0)).summary()["muB_per_fu"]
    deficit = M.Ms0({"Zn": 2.0}, hyp, n_draws=N).summary()["muB_per_fu"]         # default r ~ U(0.85, 1)
    assert deficit["median"] < exact["median"] - 0.2                              # less Zn on the down sites: lower Ms0
    # r is a source of variation (the bootstrap scatter of the moment itself scales with the Zn content, so the std alone is not
    # a measure of it): the correlation of Ms0 with r is positive (more Zn: higher Ms0)
    rr = M.Ms0({"Zn": 2.0}, hyp, n_draws=N)
    assert np.corrcoef(rr.content_ratio, rr.mu_B_per_fu)[0, 1] > 0.0          # positive, if weak (bootstrap noise dominates)
    r = M.Ms0({"Zn": 2.0}, n_draws=50, content_ratio=(0.9, 0.9))
    assert r.mass_g_per_mol == pytest.approx(M.formula_mass_g_per_mol({"Zn": 1.8}))
    for bad in ((0.0, 1.0), (1.0, 0.8), (0.9, 1.2)):
        with pytest.raises(ValueError):
            M.Ms0({"Zn": 1.0}, n_draws=10, content_ratio=bad)


def test_x_above_2_is_refused():
    from hexw.chem.validators import ModelNotApplicableError
    with pytest.raises(ModelNotApplicableError):
        M.Ms0({"Zn": 2.5}, n_draws=10)


def test_moment_assumptions_are_flagged_unidentifiable():
    est = M.moment_parameter_estimates()
    assert {e.name for e in est} == {"moment_orbital_Fe2+", "moment_orbital_Co", "moment_orbital_Ni", "moment_orbital_Mn3"}
    for e in est:
        assert not e.identifiable and e.provenance.status == "assumption" and e.unit == "muB"
        assert e.uncertainty.declared_wide


# ---------------------------------------------------------------- moment definitions are never mixed
def test_site_moment_definitions_are_kept_apart():
    cand = E.build_ensemble({"Zn": 2.0}, n_draws=40).best().candidate
    M.attach_site_moments(cand)
    for s in cand.sites:
        avg = s.moment(S.MomentDefinition.SITE_AVERAGE)
        fe = s.moment(S.MomentDefinition.FE_ION)
        assert avg is not None and fe is not None
        with pytest.raises(S.MomentDefinitionError):
            avg - fe
        assert abs(avg.value) <= abs(fe.value) + 1e-12                          # diamagnetic Zn only lowers the average
    w1 = cand.sites[0]
    assert w1.moment(S.MomentDefinition.SITE_AVERAGE).value < -1.0 and w1.moment(S.MomentDefinition.FE_ION).value < 0


# ---------------------------------------------------------------- checks against the FIT literature (not holdout)
@pytest.fixture(scope="module")
def sanity():
    return M.sanity_checks(n_draws=N)


def test_sanity_l05_total_moment_close_to_the_ensemble(sanity):
    row = sanity["L05_total"][0]
    assert row["rec"] == "L05-093" and abs(row["dft_muB_per_fu"] - 27.4955) < 1e-6
    # the DFT total (0 K) lies within 1 muB of the Gorter spin-only range 28.0 ... 29.33
    assert -1.0 < row["dft_below_model_min_by"] < 1.0
    assert abs(row["model_spin_only_by_scenario"]["oct_S_up"] - 28.0) < 1e-9


def test_sanity_l01_site_moments_use_the_site_average_definition(sanity):
    rows = sanity["L01_site_average"]
    assert len(rows) == 14                                                       # Mg and Ni, W1..W7
    for r in rows:
        # 300 K moments are reduced w.r.t. the 0 K Gorter value, never larger beyond the uncertainty
        assert abs(r["obs_300K"]) <= abs(r["model_0K_from_observed_occupancy"]) + 3 * r["sigma"]


def test_sanity_l02_fe_moments_do_not_exceed_the_ionic_value(sanity):
    rows = sanity["L02_ZnW_Fe_ion"]
    assert len(rows) == 7 and all(r["within_ionic_limit"] and 0.55 < r["ratio"] < 1.0 for r in rows)


def test_sanity_c004_room_temperature_ms_is_below_ms0(sanity):
    for r in sanity["C004_Ms_RT"]:
        assert 0.6 < r["ratio_obs_over_median"] < 0.95        # 300 K / 0 K reduction of the SrFe18O27 magnetisation


def test_sanity_l01_ms_is_reported_not_fitted(sanity):
    for row in sanity["L01_Ms"]:
        q = row["model_Ms0_A_m2_kg"]
        assert q["p16"] < q["median"] < q["p84"]
        assert row["Ms_VSM_raw"][1] < q["median"]                                # 300 K VSM below the 0 K model
