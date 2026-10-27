import pytest

from hexw.chem.validators import (ChargeError, ModelNotApplicableError, SiteOverflowError, check_charge_neutral,
                                  check_site_occupancy, fe2_fe3_model)
from hexw.registry.schemas import SiteOccupancy, StructureCandidate


def cand(*sites):
    return StructureCandidate("c1", "SrFe16Co2O27", list(sites))


def test_site_ok_and_partial():
    check_site_occupancy(cand(SiteOccupancy("12k", {"Fe": 0.75, "Co": 0.25}), SiteOccupancy("4e", {"Fe": 0.5})))


def test_site_overflow():
    with pytest.raises(SiteOverflowError, match="12k"):
        check_site_occupancy(cand(SiteOccupancy("12k", {"Fe": 0.8, "Co": 0.4})))


def test_negative_occupancy():
    with pytest.raises(SiteOverflowError):
        check_site_occupancy(cand(SiteOccupancy("12k", {"Fe": 1.2, "Co": -0.2})))


def test_charge_violation_and_ok():
    ox = {"Sr": 2, "Fe": 3, "O": -2}
    check_charge_neutral({"Sr": 1, "Fe": 52 / 3, "O": 27}, ox)
    with pytest.raises(ChargeError):
        check_charge_neutral({"Sr": 1, "Fe": 18, "O": 27}, ox)
    with pytest.raises(ChargeError, match="no oxidation"):
        check_charge_neutral({"Sr": 1, "Xx": 1}, ox)


@pytest.mark.parametrize("x", [0.0, 0.5, 1.1, 1.5, 2.0])
def test_fe_model_valid(x):
    m = fe2_fe3_model(x, "Co")
    assert m.n_fe2 == pytest.approx(2 - x) and m.n_fe3 == 16 and m.n_fe == pytest.approx(18 - x)


@pytest.mark.parametrize("x", [2.0001, 2.5, 3.0])
def test_x_above_two_rejected(x):
    with pytest.raises(ModelNotApplicableError) as e:
        fe2_fe3_model(x, "Mn")
    msg = str(e.value)
    assert "not applicable" in msg and "mixed valence" in msg and "defects" in msg


def test_negative_x_rejected():
    with pytest.raises(ModelNotApplicableError):
        fe2_fe3_model(-0.1)
