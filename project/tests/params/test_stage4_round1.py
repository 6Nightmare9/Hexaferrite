"""Stage 4 round 1 (review 4): eta from the fit samples, Fe2+ scenario weights from the L05 DFT moment, sigma_sys marginalisation,
shared draws of Ms0 / T_C / K, Ms(T). Real fit file; skipped if absent."""
import numpy as np
import pytest


def _need_fit():
    try:
        from hexw.params.literature import fit_records_path
        if not fit_records_path().is_file():
            pytest.skip("literature split not built")
    except Exception:
        pytest.skip("no project root")


@pytest.fixture(scope="module")
def mods():
    _need_fit()
    from hexw.params import anisotropy, curie, moments, occupancy
    from hexw.structures import ensemble
    return anisotropy, curie, moments, occupancy, ensemble


@pytest.fixture(scope="module")
def draws(mods):
    return mods[0].joint_draws(1500, 5, 60_000)


# ------------------------------------------------------------------------------------------- Fe2+ scenario weights
def test_fe2_weights_follow_the_l05_dft_moment(mods):
    E = mods[4]
    w = E.fe2_scenario_weights()
    assert set(w) == set(E.FE2_SCENARIOS) and abs(sum(w.values()) - 1.0) < 1e-12
    # Gorter moments of the scenarios (analytic): 29.33 / 28.0 / 29.0 muB per f.u.; the DFT moment 54.991 / 2 = 27.50
    g = {k: E.gorter_moment_x0(k) for k in E.FE2_SCENARIOS}
    assert g["uniform"] == pytest.approx(29.3333, abs=1e-3) and g["oct_S_up"] == pytest.approx(28.0) and g["oct_R"] == pytest.approx(29.0)
    assert w["oct_S_up"] > w["oct_R"] > w["uniform"]
    ll = {k: np.exp(-0.5 * ((g[k] - 54.991 / 2) / 1.0) ** 2) for k in g}
    z = sum(ll.values())
    for k in g:
        assert w[k] == pytest.approx(ll[k] / z, rel=1e-9)
    assert max(w.values()) < 0.7 and min(w.values()) > 0.05        # a Gaussian with sigma 1 muB: informative, not decisive
    # the scenarios stay, only the weights change; the members carry them
    ens = E.build_ensemble({"Co": 1.0, "Zn": 1.0})
    assert len(ens.members) == 3
    for m in ens.members:
        assert m.weight == pytest.approx(w[m.fe2_scenario]) and "L05-093" in m.weight_basis


def test_ms0_samples_scenarios_with_the_weights(mods):
    _, _, M, _, E = mods
    r = M.Ms0({"Zn": 1.0}, n_draws=6000, seed=2, content_ratio=(1.0, 1.0))
    w = E.fe2_scenario_weights()
    frac = {nm: np.mean(r.fe2_scenario == k) for k, nm in enumerate(r.scenario_names)}
    for k in w:
        assert abs(frac[k] - w[k]) < 0.02


# ------------------------------------------------------------------------------------------------- sigma_sys, bootstrap
def test_sigma_sys_is_a_parameter_with_a_prior_and_the_bootstrap_marginalises_it(mods):
    O = mods[3]
    assert O.SIGMA_SYS_RANGE == (0.01, 0.03) and O.N_BOOT == 300
    fit = O.default_fit()
    ss = fit.boot_sigma_sys
    assert ss.shape == (300,) and 0.0095 < ss.min() and ss.max() < 0.0305 and abs(ss.mean() - 0.02) < 0.003
    est = {e.name: e for e in fit.estimates}
    e = est["occ_sigma_sys"]
    from hexw.registry import schemas as S
    assert e.source_role is S.SourceRole.ASSUMPTION and e.provenance.status is S.ProvenanceStatus.ASSUMPTION and not e.identifiable
    assert "0.01" in e.provenance.note and "0.03" in e.provenance.note


def test_robust_identifiability_is_over_the_whole_sigma_sys_range(mods):
    O = mods[3]
    fit = O.default_fit()
    assert fit.n_identified <= fit.n_identified_at_min
    assert fit.n_identified < fit.n_identified_at_min                   # sigma_sys = 0.03 loses some sites that 0.01 identifies
    assert not np.any(fit.identified & ~fit.identified_at_min)
    # the robust status is what the estimates carry
    by = {e.name: e for e in fit.estimates}
    for k, sp in enumerate(fit.species):
        for i, site in enumerate(O.SITES):
            nm = f"occ_dE_{sp}_{site}"
            if nm in by:
                assert by[nm].identifiable == bool(fit.identified[k, i])
    # recompute at the two ends of the prior independently from the reference samples
    n_hi = n_lo = 0
    for k, sp in enumerate(fit.species):
        o = next(s for s in fit.samples if s.name == fit.ref_sample[sp]).species[sp]
        for i in range(7):
            if i == O.REF_SITE or not fit.has_site[k, i] or np.isnan(o.p[i]):
                continue
            n_hi += o.p[i] / o.sigma_eff(0.03)[i] >= O.Z_IDENT
            n_lo += o.p[i] / o.sigma_eff(0.01)[i] >= O.Z_IDENT
    assert n_hi == fit.n_identified and n_lo == fit.n_identified_at_min


# -------------------------------------------------------------------------------------------------------- eta, Ms(T)
def test_eta_is_estimated_from_the_three_fit_samples(mods, draws):
    _, _, M, _, _ = mods
    r = M.eta_draws(draws)
    assert set(r.per_sample) == {"C004_SrFe18O27", "L01_SrMg2W", "L01_SrNi2W"}
    assert r.eta.shape == (draws.n,) and np.all(r.eta > 0)
    ids = set(r.rec_ids)
    assert {"C004-002", "C004-004", "C004-006", "L01-069", "L01-077", "L01-070", "L01-078"} <= ids
    med = {k: float(np.median(v)) for k, v in r.per_sample.items()}
    assert all(0.6 < v < 0.95 for v in med.values())                     # ~0.74 .. 0.83
    assert med["L01_SrNi2W"] < med["L01_SrMg2W"]
    # the predictive eta carries the between-sample scatter: wider than the spread of the per-sample medians alone
    lo, hi = np.percentile(r.eta, [2.5, 97.5])
    assert lo < min(med.values()) and hi > max(med.values())
    assert 0.02 < r.scatter_ln < 0.15
    # independent arithmetic for the C004 row: eta = Ms_obs / (Ms0 m(300 K/T0))
    row = next(t for t in r.table if t["sample"] == "C004_SrFe18O27")
    assert row["Ms_obs_A_m2_kg"] == pytest.approx(np.mean([81.1, 81.1, 81.0]))
    assert row["eta_median"] == pytest.approx(row["Ms_obs_A_m2_kg"] / (row["Ms0_model_median"] * row["m_300K_median"]), rel=0.06)


def test_ms_of_t_is_eta_ms0_m(mods, draws):
    A, C, M, _, _ = mods
    res = M.ms_T([0.0, 300.0, 1100.0], {"Zn": 2.0}, draws=draws)
    ms = res.A_per_m.samples
    assert ms.shape == (3, draws.n)
    ms0 = M.Ms0({"Zn": 2.0}, seed=20260101, draws=draws)          # the seed of ms_T
    assert np.allclose(ms[0], res.eta * ms0.A_per_m)                      # m(0) = 1
    assert np.all(ms[2] == 0.0)                                           # above every T_C
    tcs = np.maximum(C.tc_samples({"Zn": 2.0}, draws.tc_params), 1.0)
    m300 = C.reduced_magnetization(300.0 / tcs, draws.shape_s, draws.form)
    assert np.allclose(ms[1], res.eta * ms0.A_per_m * m300)
    assert np.allclose(res.A_m2_per_kg.samples[1], res.eta * ms0.A_m2_per_kg * m300)
    assert np.array_equal(ms0.mu_B_per_fu, res.ms0.mu_B_per_fu)
    assert M.ms_T(300.0, {}, draws=draws).A_per_m.samples.ndim == 1       # x = 0, no hack
    # the 300 K / 0 K ratio is NOT m(T) alone: eta (~0.8) multiplies it (the old statement of a consistent ratio was wrong)
    ratio = np.median(ms[1] / ms[0])
    assert ratio == pytest.approx(np.median(m300), rel=1e-6)              # eta cancels in the ratio at fixed draw
    sr = np.median(ms[1]) / np.median(ms0.A_per_m)
    assert sr < 0.9 * np.median(m300)


def test_eta_estimate_is_a_derived_fit_parameter(mods, draws):
    from hexw.registry import schemas as S
    M = mods[2]
    (e,) = M.eta_estimates(M.eta_draws(draws))
    assert e.source_role is S.SourceRole.LITERATURE_FIT and e.provenance.status is S.ProvenanceStatus.DERIVED_BY_US
    assert "C004-002" in e.provenance.source_id and "L01-069" in e.provenance.source_id
    assert "effective" in e.provenance.note and 0.5 < e.value < 1.0


# --------------------------------------------------------------------------------------------------- shared draws
def test_ms0_shares_r_and_mn_branch_with_tc_and_k(mods, draws):
    A, C, M, _, _ = mods
    r = M.Ms0({"Mn": 1.0}, draws=draws)
    assert np.array_equal(r.content_ratio, draws.r) and np.array_equal(r.mn_valence, draws.mn_valence)
    assert r.mu_B_per_fu.shape == (draws.n,)
    # a different fixed hypothesis overrides the shared branch
    r2 = M.Ms0({"Mn": 1.0}, hypothesis=M.Hypothesis(mn_valence=3), draws=draws)
    assert np.all(r2.mn_valence == 3)
    # without draws: independent, but the same prior
    r3 = M.Ms0({"Zn": 1.0}, n_draws=2000)
    assert 0.85 <= r3.content_ratio.min() and r3.content_ratio.max() <= 1.0


def test_x0_tc_is_the_empty_composition(mods, draws):
    A, C, *_ = mods
    t0 = C.tc_samples({}, draws.tc_params)
    assert np.all(np.isfinite(t0))
    assert np.array_equal(C.tc({}).samples, C.tc({"Zn": 0.0}).samples)
