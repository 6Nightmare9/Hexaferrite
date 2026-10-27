"""K1..K4(x, T): fit-literature model, identification from T_T, synthetic recovery, branches (co_form, l11_constraint, host),
Mn by valence, cap of the Co2+ energy, validity in T, shared draws. Real fit file; skipped if absent."""
import math

import numpy as np
import pytest

N_PRIOR = 60_000


def _need_fit():
    try:
        from hexw.params.literature import fit_records_path
        if not fit_records_path().is_file():
            pytest.skip("literature split not built")
    except Exception:
        pytest.skip("no project root")


@pytest.fixture(scope="module")
def A():
    _need_fit()
    from hexw.params import anisotropy
    return anisotropy


@pytest.fixture(scope="module")
def C():
    _need_fit()
    from hexw.params import curie
    return curie


@pytest.fixture(scope="module")
def d(A):
    """The forecast of record: host dft, co_form ion_mf, l11_constraint none."""
    return A.joint_draws(2000, 7, N_PRIOR)


@pytest.fixture(scope="module")
def d_cc(A):
    return A.joint_draws(2000, 7, N_PRIOR, co_form="cc_power")


def test_defaults_are_ion_mf_dft_none(A, d):
    assert (A.DEFAULT_CO_FORM, A.DEFAULT_L11, A.DEFAULT_HOST) == ("ion_mf", "none", "dft")
    assert (d.co_form, d.l11_constraint, d.host) == ("ion_mf", "none", "dft")


@pytest.mark.parametrize("which", ["d", "d_cc"])
def test_fit_reproduces_T_T_for_both_forms(A, which, request):
    dd = request.getfixturevalue(which)
    res = A.fit_residuals(dd)
    assert dd.posterior.ess > 300
    assert len(res["T_T"]) == 4 and all(r["inside_95"] for r in res["T_T"]), res["T_T"]
    assert max(abs(r["pull"]) for r in res["T_T"]) < 1.5
    st = res["structure_300K"]
    assert st[0.0] == "plane" and st[1.0] == "plane" and st[2.0] == "axis"
    assert 0.3 < res["boundary_x_Co_300K"] < 1.0


def test_no_transition_in_zn2_and_ordering(A, C, d):
    z = A.reorientation_temperature(C.co_zn(2.0), draws=d)
    assert np.all(np.isnan(z.samples)) and z.extra["p_none"] == 1.0
    meds = [A.reorientation_temperature(C.co_zn(x), draws=d).median for x in (0.0, 0.5, 1.0, 1.5)]
    assert all(np.diff(meds) < 0)


def test_k_convention_and_identifiability_flags(A, C, d):
    r = A.anisotropy(C.co_zn(1.0), [300.0, 600.0, 900.0], draws=d)
    K1, K2, K3, K4 = r
    assert K1.samples.shape == (3, d.n)
    assert np.all(K3.samples == 0) and K3.status == "assumed_zero"
    assert not K1.identifiable and K1.extra["sign_and_T_T_identifiable"] and K1.status == "series_interpolation"
    assert not K2.identifiable and not K4.identifiable
    assert np.all(K1.samples[2] == 0) and np.all(K4.samples[2] == 0)          # above T_C of nearly all draws
    assert r.p_state["plane"][0] > 0.9
    # sign: SrZn2W axial, SrCo2W planar at 300 K
    assert A.anisotropy(C.co_zn(2.0), 300.0, draws=d).K1.prob(lambda k: k > 0) == 1.0
    # K1 > 0 with a strongly negative K2 (K1 + K2 < 0) is still a planar state: P(K1 < 0) < P(plane)
    assert A.anisotropy(C.co_zn(0.0), 300.0, draws=d).K1.prob(lambda k: k < 0) > 0.8


# -------------------------------------------------------------------------------------------- the two Co forms
def test_both_co_forms_are_returned_and_not_averaged(A, C):
    res = A.anisotropy_branches(C.co_zn(0.0), 300.0, n=1500, seed=5)
    assert len(res) == 8 and set(res) == {(h, f, l) for h in A.HOSTS for f in A.CO_FORMS for l in A.L11_CONSTRAINTS}
    a, b = res[("dft", "ion_mf", "none")], res[("dft", "cc_power", "none")]
    assert a.co_form == "ion_mf" and b.co_form == "cc_power"
    assert a.K1.extra["co_form"] == "ion_mf"
    assert not np.array_equal(a.K1.samples, b.K1.samples)
    # the default of record is ion_mf
    dflt = A.anisotropy(C.co_zn(0.0), 300.0, n=1500, seed=5)
    assert np.array_equal(dflt.K1.samples, a.K1.samples)
    with pytest.raises(ValueError):
        A.anisotropy(C.co_zn(0.0), 300.0, co_form="mean")
    with pytest.raises(ValueError):
        A.anisotropy(C.co_zn(0.0), 300.0, l11_constraint="maybe")


def _series_k1(A, C, dd, x_zn, T):
    """K1 of the Co-Zn series recomputed from the draws by hand (no use of _k_arrays)."""
    comp = C.co_zn(x_zn)
    h = dd.hosts[dd.host]
    tcs = np.maximum(C.tc_samples(comp, dd.tc_params), 1.0)
    t = T / tcs
    m = C.reduced_magnetization(t, dd.shape_s, dd.form)
    cs = dd.tc_params["content_scale"]
    x_co = comp["Co"] * cs
    if dd.co_form == "ion_mf":
        f_co = C.brillouin_m(t / dd.f, 1.5) ** 3
    else:
        f_co = m ** dd.p
    kh = h["K_fe"] + cs * (h["A"] - h["K_fe"])           # (x_Co + x_Zn) cs / 2 = cs for the series
    return kh * m ** 3 - dd.rho * h["A"] * x_co * f_co, kh, m


@pytest.mark.parametrize("which", ["d", "d_cc"])
def test_k1_formula_of_each_co_form(A, C, which, request):
    dd = request.getfixturevalue(which)
    for x_zn in (0.0, 1.0):
        k1, _, _ = _series_k1(A, C, dd, x_zn, 300.0)
        got = A.anisotropy(C.co_zn(x_zn), 300.0, draws=dd).K1.samples
        assert np.allclose(got, k1, rtol=1e-9, atol=1e-3)


def test_ion_mf_uses_s32_brillouin_cubed_and_theta_fraction(A, C, d):
    """F = m_ion(T / (f T_C))^3 with the S = 3/2 mean-field table; F(0) = 1, F = 0 above Theta_Co."""
    assert 0.15 <= d.f.min() and d.f.max() <= 1.0
    f_hi = np.where(d.f > 0, 1.0, 0.0)
    assert np.allclose(A.co_F("ion_mf", 1.0, 0.0, 1.0, d.f), 1.0)
    assert np.all(A.co_F("ion_mf", 0.5, 1.2 * d.f, 1.0, d.f) == 0.0) and f_hi.all()
    assert float(A.co_F("ion_mf", 1.0, 0.5, 0.0, 1.0)) == pytest.approx(float(C.brillouin_m(0.5, 1.5)) ** 3)


def test_p_is_effective_and_not_identifiable(A):
    ests = {e.name: e for e in A.anisotropy_parameter_estimates(n=1500)}
    assert not ests["aniso_p_Co_exponent[cc_power]"].identifiable
    assert "NOT identifiable" in ests["aniso_p_Co_exponent[cc_power]"].provenance.note
    comb = ests["aniso_CoTerm_over_host_300K_Co2[cc_power]"]
    assert comb.identifiable                                  # the combination is what the data fix


# ------------------------------------------------------------------------------------ K2/K4 scale, T validity, cap
def test_k2_k4_scale_with_the_full_single_ion_contribution(A, C, d):
    comp = C.co_zn(0.0)
    r = A.anisotropy(comp, 300.0, draws=d)
    h = d.hosts[d.host]
    cs = d.tc_params["content_scale"]
    kh = h["K_fe"] + 2.0 * cs * (h["A"] - h["K_fe"]) / 2.0
    scale = np.abs(kh) + d.rho * np.abs(h["A"]) * 2.0 * cs / 2.0
    tcs = np.maximum(C.tc_samples(comp, d.tc_params), 1.0)
    m = C.reduced_magnetization(300.0 / tcs, d.shape_s, d.form)
    assert np.allclose(r.K2.samples, d.q2 * scale * m ** 10, rtol=1e-9, atol=1e-6)
    assert np.allclose(r.K4.samples, d.q4 * scale * m ** 21, rtol=1e-9, atol=1e-6)
    # the old scaling with A only would be much smaller for a Co-rich composition
    assert np.median(scale / np.abs(h["A"])) > 2.0


def test_low_T_is_extrapolated_with_a_wider_interval(A, C, d):
    r = A.anisotropy(C.co_zn(1.0), [0.0, 100.0, 300.0], draws=d)
    assert r.K1.status == "extrapolated_in_T" and r.K1.extra["composition_status"] == "series_interpolation"
    assert list(r.K1.extra["T_extrapolated"]) == [True, True, False]
    ok = A.anisotropy(C.co_zn(1.0), [300.0, 400.0], draws=d)
    assert ok.K1.status == "series_interpolation" and not ok.K1.extra["T_extrapolated"].any()
    assert A.anisotropy(C.co_zn(1.0), 250.0, draws=d).K1.status == "series_interpolation"
    assert A.anisotropy(C.co_zn(1.0), 249.0, draws=d).K1.status == "extrapolated_in_T"
    # the interval at 0 K is wider than without the extrapolation factor (same draws)
    k1w = A._k_arrays(C.co_zn(1.0), np.array([[0.0]]), d, widen=True)[0][0]
    k1n = A._k_arrays(C.co_zn(1.0), np.array([[0.0]]), d, widen=False)[0][0]
    w = lambda k: np.subtract(*np.quantile(k, [0.975, 0.025]))     # noqa: E731
    assert w(k1w) > 1.15 * w(k1n)
    assert np.array_equal(np.sign(k1w), np.sign(k1n))                  # the factor never flips a sign
    # no restriction from 300 K up
    assert np.allclose(A._k_arrays(C.co_zn(1.0), np.array([[300.0]]), d, widen=True)[0],
                       A._k_arrays(C.co_zn(1.0), np.array([[300.0]]), d, widen=False)[0])


@pytest.mark.parametrize("host", ["dft", "l11"])
def test_co_single_ion_energy_is_capped(A, host):
    dd = A.joint_draws(1500, 11, N_PRIOR, host=host)
    e_mev = dd.rho * np.abs(dd.hosts[host]["A"]) * A.v_fu_m3() / A.MEV_J
    assert e_mev.max() <= A.E_CO_MAX_MEV * (1 + 1e-9)
    e_mn = dd.mn["rho"] * np.abs(dd.hosts[host]["A"]) * A.v_fu_m3() / A.MEV_J
    assert e_mn.max() <= A.E_CO_MAX_MEV * (1 + 1e-9)
    assert dd.rho.min() >= A.RHO_RANGE[0] * 0.999


def test_cap_cuts_the_prior_when_it_binds(A):
    huge = np.array([1e7])                       # 10 MJ/m3 host: only rho < 0.65 would fit
    assert A.rho_cap(huge)[0] < A.RHO_RANGE[0] * 2
    assert A.rho_cap(np.array([1e5]))[0] == pytest.approx(A.E_CO_MAX_MEV * A.MEV_J / (1e5 * A.v_fu_m3()))


# ----------------------------------------------------------------------------------------------- L11 branches
def test_l11_ratio_branch_changes_the_posterior_toward_the_l11_ratio(A, C):
    base = A.joint_draws(2000, 9, N_PRIOR)
    cons = A.joint_draws(2000, 9, N_PRIOR, l11_constraint="co_zn_ratio")
    assert cons.l11_constraint == "co_zn_ratio" and cons.posterior.ess > 150

    def ratio(dd):
        r = A.anisotropy(C.co_zn(0.0), 300.0, draws=dd).K1.samples
        z = A.anisotropy(C.co_zn(2.0), 300.0, draws=dd).K1.samples
        return np.abs(r) / np.abs(z)

    tgt = A.l11_co_zn_ratio()
    assert 0.5 < tgt < 1.5                                   # the L11 number (~0.9) from the fit records
    dev = lambda dd: abs(np.median(np.log(ratio(dd) / tgt)))        # noqa: E731
    assert dev(cons) <= dev(base) + 0.05
    assert A.SIGMA_LN_L11_RATIO == 0.7
    ten = A.l11_ratio_tension("ion_mf", n=1500, seed=9)
    assert set(ten) >= {"quantiles", "L11_ratio", "P_below_L11", "z_ln"} and ten["quantiles"]["p50"] > 0


def test_l11_formula_note_and_sigma(A):
    import inspect
    assert A.SIGMA_LN_L11 >= 0.8 and A.SIGMA_LN_L11_MN >= A.SIGMA_LN_L11
    src = inspect.getsource(A)
    assert "H >> Ha" in src and "sign-blind" in src
    assert "tier 3" in src


def test_ba_to_sr_factor_from_the_in_fit_l05_sensitivity(A):
    hin = A.host_inputs()
    assert hin.ids["ku"] == "L05-103" and hin.ids["ku_srba"] == "L05-112"
    k_ba = 2 * hin.ku_srba - hin.ku_dft
    assert A.ba_to_sr_factor(hin) == pytest.approx(hin.ku_dft / k_ba)
    assert 0.2 < A.ba_to_sr_factor(hin) < 0.45
    dl = A.joint_draws(1500, 3, N_PRIOR, host="l11")
    dd = A.joint_draws(1500, 3, N_PRIOR, host="dft")
    # with the transfer the two hosts agree within a factor ~2 in the median (no longer the 2.4x of the Ba value)
    r = np.median(dl.hosts["l11"]["K_fe"]) / np.median(dd.hosts["dft"]["K_fe"])
    assert 0.5 < r < 2.0


def test_alternatives_not_averaged(A, C):
    a = A.anisotropy(C.co_zn(2.0), 300.0, host="dft", n=1500, seed=3).K1
    b = A.anisotropy(C.co_zn(2.0), 300.0, host="l11", n=1500, seed=3).K1
    assert a.extra["host"] == "dft" and b.extra["host"] == "l11" and not np.array_equal(a.samples, b.samples)
    with pytest.raises(ValueError):
        A.anisotropy(C.co_zn(2.0), 0.0, host="mean")


# --------------------------------------------------------------------------------------------------- Mn
def test_mn_sign_follows_the_valence_branch(A, C, d):
    h = d.hosts[d.host]
    v2, v3 = d.mn_valence == 2, d.mn_valence == 3
    assert 0.4 < v3.mean() < 0.6
    assert abs(np.mean(d.mn["sign2"] > 0) - A.MN2_AXIAL_WEIGHT) < 0.04
    assert abs(np.mean(d.mn["sign3"] > 0) - 0.5) < 0.05
    T = np.array([[300.0]])
    k1 = A._k_arrays({"Mn": 1.0}, T, d, widen=False)[0][0]
    cs = d.tc_params["content_scale"]
    tcs = np.maximum(C.tc_samples({"Mn": 1.0}, d.tc_params), 1.0)
    m = C.reduced_magnetization(300.0 / tcs, d.shape_s, d.form)
    # Mn2+: the host-like axial/planar term with the sign of the 0.85/0.15 draw, no single-ion term
    a_mn = d.mn["sign2"] * np.abs(h["A_mn"])
    k2p = (h["K_fe"] + cs * (a_mn - h["K_fe"]) / 2.0) * m ** 3
    assert np.allclose(k1[v2], k2p[v2], rtol=1e-9)
    # Mn3+: Zn-like dilution + a single-ion term of the Co form with a 50/50 sign
    t = 300.0 / tcs
    f3 = C.brillouin_m(t / d.mn["f"], 1.5) ** 3
    k3p = (h["K_fe"] + cs * (h["A"] - h["K_fe"]) / 2.0) * m ** 3 + d.mn["sign3"] * d.mn["rho"] * h["A"] * cs * f3
    assert np.allclose(k1[v3], k3p[v3], rtol=1e-9)
    # the branches really differ
    assert abs(np.median(k1[v3]) - np.median(k1[v2])) > 0 and np.mean(k1[v2] < 0) > 0.05


def test_mn_valence_is_shared_with_ms0(A, d):
    from hexw.params.moments import Ms0
    res = Ms0({"Mn": 1.0}, draws=d, fit=None)
    assert np.array_equal(res.mn_valence, d.mn_valence)
    assert np.array_equal(res.content_ratio, d.r)


def test_ni_transfer_and_mg_uncalibrated_scenario(A, d):
    r = A.anisotropy({"Ni": 2.0}, 300.0, draws=d)
    assert r.status == "transfer_tier3" and all(np.isfinite(k.samples).all() for k in r)
    assert not r.K1.identifiable and "L01-089" in r.notes
    mg = A.anisotropy({"Mg": 2.0}, 300., draws=d)
    assert mg.status == "uncalibrated_scenario" and not mg.K1.identifiable
    assert np.all(mg.K1.samples>0) and mg.p_state["axis"]==pytest.approx([1.])
    mn = A.anisotropy({"Mn": 1.0}, 300.0, draws=d)
    assert mn.status == "transfer_tier3" and not mn.K1.identifiable
    zn = A.anisotropy({"Zn": 1.0}, 300.0, draws=d)
    assert (mn.K1.hi - mn.K1.lo) >= 0.9 * (zn.K1.hi - zn.K1.lo)


# ----------------------------------------------------------------------------------------------- x = 0, joint
def test_x0_is_supported_without_a_hack(A, C, d):
    a = A.anisotropy({}, [300.0, 500.0], draws=d)
    b = A.anisotropy({"Co": 0.0}, [300.0, 500.0], draws=d)
    assert a.status == "extrapolated" and np.array_equal(a.K1.samples, b.K1.samples)
    assert np.all(np.isfinite(a.K1.samples))
    h = d.hosts[d.host]
    tcs = np.maximum(C.tc_samples({}, d.tc_params), 1.0)
    m = C.reduced_magnetization(300.0 / tcs, d.shape_s, d.form)
    assert np.allclose(a.K1.samples[0], h["K_fe"] * m ** 3, rtol=1e-9)
    assert A.reorientation_temperature({}, draws=d).extra["p_none"] > 0.9          # axial host: the plane is the global minimum only if q2 < -1


def test_shared_content_ratio_in_tc_and_k(A, C, d):
    assert np.array_equal(d.tc_params["r"], d.r) and np.allclose(d.tc_params["content_scale"], d.r / C.R_FIT)
    assert 0.85 <= d.r.min() and d.r.max() <= 1.0
    # a larger r dilutes T_C more (Zn): same draws, so the correlation is through r
    t = C.tc_samples({"Zn": 1.0}, d.tc_params)
    assert np.corrcoef(t, d.r)[0, 1] < -0.05


def test_branch_independent_quantities_are_common(A):
    a = A.joint_draws(1000, 4, N_PRIOR)
    b = A.joint_draws(1000, 4, N_PRIOR, co_form="cc_power")
    c = A.joint_draws(1000, 4, N_PRIOR, host="l11")
    for o in (b, c):
        assert np.array_equal(a.r, o.r) and np.array_equal(a.mn_valence, o.mn_valence)
        assert np.array_equal(a.tc_params["tc_co2"], o.tc_params["tc_co2"])
        assert np.array_equal(a.low_t_z, o.low_t_z)


def test_acceptance_grid_finite_or_not_predicted(A, C, d):
    Ts = [0.0, 100.0, 300.0, 500.0, 800.0]
    for xz in np.arange(0.0, 2.01, 0.25):
        r = A.anisotropy(C.co_zn(xz), Ts, draws=d)
        for k in r:
            assert np.all(np.isfinite(k.samples)), (xz, k.name)
    for a in ("Co", "Zn", "Mn"):
        for x in (0.5, 1.0, 1.5, 2.0):
            r = A.anisotropy({a: x}, Ts, draws=d)
            assert r.status != "not_predicted" and np.all(np.isfinite(r.K1.samples))
            assert np.all(np.isfinite(C.tc({a: x}).samples))
    assert A.anisotropy({"Mg": 1.0}, Ts, draws=d).status == "extrapolated_in_T"


def test_easy_state(A):
    assert A.easy_state(1.0, 0.0) == 0 and A.easy_state(-1.0, 0.0) == 1
    assert A.easy_state(-1.0, 1.0) == 2            # K1 < 0, K1 + 2 K2 > 0 : cone
    assert A.easy_state(-1.0, 0.2) == 1 and A.easy_state(0.5, -2.0) == 1


@pytest.mark.parametrize("co_form,truth", [("cc_power", dict(rho=5.0, p=12.0, f=0.5, s=0.8, q2=0.0)),
                                           ("ion_mf", dict(rho=3.5, p=12.0, f=0.67, s=0.8, q2=0.0))])
def test_synthetic_recovery(A, co_form, truth):
    """The truth must itself satisfy the structure constraints kept from the real data (planar at 300 K for x_Zn = 0, 1: T_T > 300 K).
    Generate T_T with known parameters from the model, add noise, refit: truth inside the 95 % intervals of the observable
    (held-out T_T at x_Zn = 1.25) and of the fitted parameter combination."""
    base = A.tt_data()

    def one(x, tc):
        a = {k: np.array([v]) for k, v in truth.items()}
        return A.first_exit_t(a["rho"], a["p"], a["f"], a["s"], a["q2"], 2.0 - x, "kuzmin", co_form)[0] * tc

    tt_true = np.array([one(x, tc) for x, tc in zip(base.x_zn, base.tc)])
    rng = np.random.default_rng(3)
    sig = np.full(4, 15.0)
    syn = A.TTData(base.x_zn, base.tc, tt_true + rng.normal(0, sig), sig, base.planar_300, base.axial_300, base.tc_zn2)
    post = A.fit_tt(syn, n_prior=N_PRIOR, seed=11, co_form=co_form)
    assert post.ess > 100
    idx = np.random.default_rng(5).choice(len(post.weights), 3000, p=post.weights)
    lo, hi = np.quantile(post.rho[idx], [0.025, 0.975])
    assert lo <= truth["rho"] <= hi, ("rho", lo, hi)
    tcv = float(np.interp(1.25, base.x_zn, base.tc))
    held = A.first_exit_t(post.rho[idx], post.p[idx], post.f[idx], post.s[idx], post.q2[idx], 2.0 - 1.25, "kuzmin", co_form) * tcv
    lo, hi = np.nanquantile(held, [0.025, 0.975])
    assert lo <= one(1.25, tcv) <= hi
    assert (hi - lo) < 150.0                       # the observable is much better determined than the individual parameters


def test_chi2_profile_runs_for_both_forms(A):
    for cf in A.CO_FORMS:
        rg, pg, chi2, feas = A.chi2_profile(ng=31, co_form=cf)
        assert chi2.shape == (31, 31) and feas.any() and np.nanmin(chi2[feas]) < 10.0


# ------------------------------------------------------------------------------------------- estimates / schema
def test_parameter_estimates_roles_and_identifiability(A):
    from hexw.registry import schemas as S
    ests = A.anisotropy_parameter_estimates(n=2000)
    nm = {e.name: e for e in ests}
    for e in ests:
        has = "ASSUMPTION" in e.provenance.source_id.upper()
        assert has == (e.source_role is S.SourceRole.ASSUMPTION), e.name
        assert not (has and ";" in e.provenance.source_id), f"mixed record {e.name}"
        if has:
            assert e.provenance.status is S.ProvenanceStatus.ASSUMPTION and not e.identifiable
        else:
            assert e.source_role is S.SourceRole.LITERATURE_FIT
            assert e.provenance.status in (S.ProvenanceStatus.VERIFIED, S.ProvenanceStatus.DERIVED_BY_US)
    assert nm["aniso_K1_SrFe18O27_0K_L05"].provenance.source_id == "L05-103"
    assert nm["aniso_K1_SrFe18O27_0K_L05"].provenance.status is S.ProvenanceStatus.VERIFIED
    assert not nm["aniso_q2_K2_over_Kscale"].identifiable and not nm["aniso_q4_K4_over_Kscale"].identifiable
    assert "L02-024" in nm["aniso_rho_Co_over_A[ion_mf]"].provenance.source_id
    assert nm["aniso_A_Zn2_0K_L11"].provenance.source_id == "L11-010;L11-012"
    assert nm["aniso_E_Co_max_meV"].value == 20.0 and nm["aniso_E_Co_max_meV"].source_role is S.SourceRole.ASSUMPTION
    assert nm["aniso_Ba_to_Sr_factor"].provenance.source_id == "L05-103;L05-112"
    # identifiable = the data narrow the prior by >= 2x (95 % widths); a non-identifiable estimate reports the prior interval
    for e in ests:
        if "posterior/prior" in (e.provenance.note or ""):
            ratio = float(e.provenance.note.split("posterior/prior = ")[1].split(" ")[0])
            if e.source_role is S.SourceRole.LITERATURE_FIT:
                assert e.identifiable == (ratio <= 0.5 + 1e-9), (e.name, ratio)
