"""T_C(composition) and m(T) from the fit literature (L02 Table 3, L01 Mg/Ni TG). Real fit file; skipped if absent."""
import numpy as np
import pytest


def _need_fit():
    try:
        from hexw.params.literature import fit_records_path
        if not fit_records_path().is_file():
            pytest.skip("literature split not built")
    except Exception:           # no configured root
        pytest.skip("no project root")


@pytest.fixture(scope="module")
def C():
    _need_fit()
    from hexw.params import curie
    return curie


def test_series_reproduces_table(C):
    f = C.fit_inputs()
    for xz, obs, sg in zip(f.x_zn, f.tc_dmax, f.tc_sigma):
        d = C.tc(C.co_zn(xz))
        assert d.status == "interpolated_in_series" and d.identifiable
        assert d.lo <= obs <= d.hi, (xz, obs, d.interval)
        assert abs(d.median - obs) < 2.5 * max(sg, 10.0)


def test_table_vs_abstract_goes_into_uncertainty(C):
    f = C.fit_inputs()
    # abstract numbers 780 / 680 differ from the table (790 / 700): sigma > read-off sigma at x_Zn = 0, 1 only
    assert f.tc_sigma[0] > C.SIGMA_TC_READ_K and f.tc_sigma[2] > f.tc_sigma[0]
    assert f.tc_sigma[1] == pytest.approx(C.SIGMA_TC_READ_K)
    assert f.tc_dmax[0] == 790.0 and f.tc_dmax[2] == 700.0       # main value: the table


def test_outside_domain_is_wide_and_flagged(C):
    mn = C.tc({"Mn": 2.0})
    assert mn.status == "wide_prior" and not mn.identifiable
    zn = C.tc(C.co_zn(2.0))
    assert (mn.hi - mn.lo) > 2 * (zn.hi - zn.lo)
    host = C.tc({"Co": 1.0})
    assert not host.identifiable and host.status in ("extrapolated_from_endpoints", "wide_prior")


def test_x0_is_supported_explicitly(C):
    """x = 0 (SrFe18O27) is the empty composition; {} and {"Co": 0} are the same thing, no hack needed."""
    a, b = C.tc({}), C.tc({"Co": 0.0})
    assert np.array_equal(a.samples, b.samples)
    assert a.status == "wide_prior" and not a.identifiable
    assert np.all(C.ms_reduced(300.0, {}).samples > 0)
    with pytest.raises(ValueError):
        C.normalize_composition("Co")           # not a dict


def test_tc_scatter_independent_between_compositions(C):
    """The interpolation scatter is one independent draw per composition (not one number shared by all of them)."""
    from hexw.params import anisotropy as A
    d = A.joint_draws(1000, 3, 30_000)
    s1 = C.tc_scatter(C.co_zn(0.0), d.tc_params)
    s2 = C.tc_scatter(C.co_zn(1.0), d.tc_params)
    assert abs(np.corrcoef(s1, s2)[0, 1]) < 0.15
    assert np.array_equal(s1, C.tc_scatter(C.co_zn(0.0), d.tc_params))        # reproducible
    assert 5.0 < np.std(s1) < 30.0


def test_mg_ni_other_definition(C):
    for a in ("Mg", "Ni"):
        d = C.tc({a: 2.0})
        assert d.status == "measured_other_definition" and "dmax" in d.notes


def test_bad_compositions(C):
    for bad in ({"Xx": 1.0}, {"Co": -1.0}, {"Co": float("nan")}, {"Co": 1.5, "Zn": 1.5}):
        with pytest.raises(ValueError):
            C.tc(bad)


@pytest.mark.parametrize("form", ["kuzmin", "brillouin"])
def test_m_reduced_shape(C, form):
    t = np.linspace(0, 1.2, 121)
    m = C.reduced_magnetization(t, 0.8, form)
    assert m[0] == 1.0 and np.all(m[t >= 1.0] == 0.0) and np.all(np.diff(m) <= 1e-12)
    sel = (m > 0.05) & (m < 1)
    back = C.inverse_reduced_magnetization(m[sel], 0.8, form)
    assert np.allclose(C.reduced_magnetization(back, 0.8, form), m[sel], atol=2e-3)


def _mf_solution(S, t):
    """Independent solution of the mean-field self-consistency m = B_S(3S/(S+1) m / t) (scipy fsolve)."""
    from scipy.optimize import fsolve

    def B(y):
        a, b = (2 * S + 1) / (2 * S), 1 / (2 * S)
        return a / np.tanh(a * y) - b / np.tanh(b * y)

    return float(fsolve(lambda m: m - B(3 * S / (S + 1) * m / t), 0.9, xtol=1e-13)[0])


@pytest.mark.parametrize("spin", [1.5, 2.5])
def test_brillouin_table_matches_the_self_consistency_equation(C, spin):
    """Regression (review): the y -> 0 end of the table used to collapse it (m[0] = 1e-12 was paired with t ~ 0). The
    limits are m(0) = 1, m(1) = 0 (y -> 0 is the Curie point t = 1, m = 0); values from fsolve."""
    tt, mm = C._brillouin_table(spin)
    assert tt[0] == 0.0 and mm[0] == 1.0 and tt[-1] == 1.0 and mm[-1] == 0.0
    assert np.all(np.diff(tt) > 0) and np.all(np.diff(mm) <= 0)
    for t in (0.2, 0.5, 0.8, 0.95, 0.99):
        assert float(C.brillouin_m(t, spin)) == pytest.approx(_mf_solution(spin, t), abs=2e-4)
    if spin == 2.5:
        assert 0.87 < float(C.brillouin_m(0.5, 2.5)) < 0.90           # ~0.889
        assert np.allclose(C.reduced_magnetization(np.array([0.0, 0.5, 1.0]), 0.8, "brillouin"),
                           [1.0, _mf_solution(2.5, 0.5), 0.0], atol=2e-4)
    assert float(C.brillouin_m(0.0, spin)) == 1.0 and float(C.brillouin_m(1.0, spin)) == 0.0
    assert float(C.brillouin_m(1.7, spin)) == 0.0


def test_ms_reduced_api(C):
    d = C.ms_reduced([0.0, 300.0, 600.0, 1000.0], C.co_zn(1.0))
    assert d.samples.shape[0] == 4
    assert np.allclose(d.samples[0], 1.0) and np.all(d.samples[3] == 0.0)
    assert np.all(np.diff(d.median) <= 0)
    assert C.ms_reduced(300.0, C.co_zn(1.0)).samples.ndim == 1
    with pytest.raises(ValueError):
        C.ms_reduced(-1.0, C.co_zn(1.0))


def test_parameter_estimates(C):
    from hexw.registry import schemas as S
    ests = C.curie_parameter_estimates(n=4000)
    names = {e.name: e for e in ests}
    assert names["Tc_Co2_dmax"].provenance.source_id.startswith("L02-")
    assert names["Tc_Co2_dmax"].source_role is S.SourceRole.LITERATURE_FIT
    assert names["Tc_Co2_dmax"].provenance.status is S.ProvenanceStatus.DERIVED_BY_US        # a line fit, not a record value
    assert names["Tc_Mg2_TG"].provenance.source_id == "L01-070"
    assert names["Tc_Mg2_TG"].provenance.status is S.ProvenanceStatus.VERIFIED
    assert not names["Tc_Mn_dilution_fraction"].identifiable
    # assumptions are assumptions: role, id and status; no assumption claims the status 'verified'; no mixed ids
    for e in ests:
        assert ("ASSUMPTION" in e.provenance.source_id) == (e.source_role is S.SourceRole.ASSUMPTION)
        if e.source_role is S.SourceRole.ASSUMPTION:
            assert e.provenance.source_id == "ASSUMPTION" and e.provenance.status is S.ProvenanceStatus.ASSUMPTION
        else:
            assert e.provenance.status is not S.ProvenanceStatus.ASSUMPTION
    for nm in ("Tc_dilution_slope_Co", "Tc_def_offset_TG_minus_dmax", "m_shape_s", "A_content_ratio_prediction"):
        assert names[nm].source_role is S.SourceRole.ASSUMPTION
