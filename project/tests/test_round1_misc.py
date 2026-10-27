import math

import pytest

from hexw import units
from hexw.chem import validators as V
from hexw.registry import schemas as S
from hexw.registry.versioning import parameter_set_hash


def P(v=10.0, **kw):
    return S.ParameterEstimate("Ms", v, "A/m", "magnetization", S.SourceRole.DFT, S.Provenance("src"), **kw)


# ---- 6: uncertainty
@pytest.mark.parametrize("u", [
    S.Uncertainty("std", std=0.01),
    S.Uncertainty("interval", lower=9.9, upper=10.1),
    S.Uncertainty("distribution", distribution="normal", params={"sigma": 0.1}),
    S.Uncertainty("distribution", distribution="uniform", params={"low": 9.9, "high": 10.1}),
])
def test_nonidentifiable_rejects_narrow_for_all_kinds(u):
    with pytest.raises(ValueError):
        P(identifiable=False, uncertainty=u)
    P(identifiable=True, uncertainty=u)


def test_nonidentifiable_wide_unknown_and_flag():
    P(identifiable=False, uncertainty=S.Uncertainty("std", std=6.0))
    P(identifiable=False, uncertainty=S.Uncertainty("interval", lower=1.0, upper=19.0))
    P(identifiable=False)                                              # unknown
    with pytest.raises(ValueError):                                    # unrecognised distribution width
        P(identifiable=False, uncertainty=S.Uncertainty("distribution", distribution="lognormal"))
    P(identifiable=False, uncertainty=S.Uncertainty("distribution", distribution="lognormal", declared_wide=True))
    with pytest.raises(ValueError):                                    # flag never overrides a narrow width
        P(identifiable=False, uncertainty=S.Uncertainty("std", std=0.01, declared_wide=True))
    with pytest.raises(ValueError):                                    # value 0: needs the flag
        P(0.0, identifiable=False, uncertainty=S.Uncertainty("std", std=1e-9))


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -float("inf")])
def test_uncertainty_rejects_nan_inf(bad):
    for kw in ({"kind": "std", "std": bad}, {"kind": "interval", "lower": bad, "upper": 1.0},
               {"kind": "interval", "lower": 0.0, "upper": bad}, {"kind": "unknown", "std": bad},
               {"kind": "distribution", "distribution": "normal", "params": {"sigma": bad}}):
        with pytest.raises(ValueError):
            S.Uncertainty(**kw)


# ---- 7: validators with NaN/inf
def test_validators_reject_nan_inf():
    nan, inf = float("nan"), float("inf")
    for occ in (nan, inf, -inf):
        c = S.StructureCandidate("c", "x", [S.SiteOccupancy("12k", {"Fe": occ})])
        with pytest.raises(V.SiteOverflowError):
            V.check_site_occupancy(c)
    with pytest.raises(V.SiteOverflowError):
        V.check_site_occupancy(S.StructureCandidate("c", "x", [S.SiteOccupancy("a", {"Fe": 0.5, "Co": nan})]))
    for n in (nan, inf):
        with pytest.raises(V.ChargeError):
            V.check_charge_neutral({"Sr": n, "O": 1}, {"Sr": 2, "O": -2})
    with pytest.raises(V.ChargeError):
        V.check_charge_neutral({"Sr": 1}, {"Sr": nan})
    for x in (nan, inf, -inf):
        with pytest.raises(V.ModelNotApplicableError):
            V.fe2_fe3_model(x, "Co")


# ---- 8: hash canonical numbers
def test_parameter_set_hash_int_float_equal():
    a = [S.ParameterEstimate("n", 1, "1", "dimensionless", "dft", S.Provenance("s"))]
    b = [S.ParameterEstimate("n", 1.0, "1", "dimensionless", "dft", S.Provenance("s"))]
    assert parameter_set_hash(a) == parameter_set_hash(b)
    c = [S.ParameterEstimate("n", 1.0000001, "1", "dimensionless", "dft", S.Provenance("s"))]
    assert parameter_set_hash(a) != parameter_set_hash(c)
    u1 = [S.ParameterEstimate("n", 1, "1", "dimensionless", "dft", S.Provenance("s"), S.Uncertainty("std", std=1))]
    u2 = [S.ParameterEstimate("n", 1, "1", "dimensionless", "dft", S.Provenance("s"), S.Uncertainty("std", std=1.0))]
    assert parameter_set_hash(u1) == parameter_set_hash(u2)


# ---- 9: reserved names
@pytest.mark.parametrize("A", ["Sr", "O", "Fe", "Fe2+", "Fe3+", "", None])
def test_fe2_fe3_model_reserved_names(A):
    with pytest.raises(V.ChemistryError, match="not allowed"):
        V.fe2_fe3_model(1.0, A)
    assert V.fe2_fe3_model(1.0, "Co").n_fe2 == 1.0


# ---- 10: units
def test_gamma_e_over_2pi_ghz_per_tesla():
    g = units.gamma_from_g(2.00231930436)
    f = units.gamma_over_2pi_from_gamma(g).value / 1e9
    assert f == pytest.approx(28.0249514, rel=1e-8)


def test_cgs_units_and_conventions():
    q = units.gamma_over_2pi(28.0249514e9)
    assert units.Quantity(1.0, "MHz/Oe", units.Kind.GAMMA_OVER_2PI).canonical().value == pytest.approx(1e10)
    assert units.Quantity(1.0, "GHz/kOe", units.Kind.GAMMA_OVER_2PI).canonical().value == pytest.approx(1e10)
    # 2.8 MHz/Oe == 28 GHz/T
    assert units.Quantity(2.8, "MHz/Oe", units.Kind.GAMMA_OVER_2PI).canonical().value == pytest.approx(28e9)
    assert q.unit == "Hz/T"
    w = units.angular_frequency_from_frequency(units.Quantity(1.0, "GHz", units.Kind.FREQUENCY))
    assert w.unit == "rad/s" and w.value == pytest.approx(2 * math.pi * 1e9)
    m = units.Quantity(1.0, "emu/cm3", units.Kind.MAGNETIZATION)
    assert m.canonical().value == pytest.approx(1000.0)
    f4 = units.fourpiM_from_M(m)
    assert f4.unit == "G" and f4.value == pytest.approx(4 * math.pi)
    assert units.M_from_fourpiM(f4).value == pytest.approx(1000.0)
    assert units.mu0H_from_H(units.H(1e5)).value * 1e4 == pytest.approx(1e5 * units.MU0 * 1e4)
    # B_hf: 1 kOe <-> 0.1 T
    assert units.B_hf(30.0, "kOe").canonical().value == pytest.approx(3.0)


def test_constants_come_from_scipy_codata():
    import scipy.constants as c
    assert units.MU0 == c.mu_0 and units.H_PLANCK == c.h and units.HBAR == c.hbar and units.E_CHARGE == c.e
    assert units.CODATA_VERSION == "CODATA 2022"
