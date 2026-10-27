import json
import math

import pytest
import yaml

from hexw import units
from hexw.registry import schemas as S


def test_field_like_kinds_are_distinct_and_not_addable():
    h, m = units.H(1e5), units.mu0H(0.1)
    assert h.kind != m.kind
    with pytest.raises(units.UnitError):
        _ = h + m
    with pytest.raises(units.UnitError):
        units.B_hf(30, "A/m")
    with pytest.raises(units.UnitError):
        units.mu0H(1.0, "A/m")


def test_mu0H_conversion_explicit():
    m = units.mu0H_from_H(units.H(1e5))
    assert m.kind == units.Kind.MU0H and m.unit == "T" and m.value == pytest.approx(0.1256637, rel=1e-6)
    assert units.H_from_mu0H(m).value == pytest.approx(1e5)


def test_gamma_vs_gamma_over_2pi():
    g = units.gamma(1.76085963e11)
    g2 = units.gamma_over_2pi_from_gamma(g)
    assert g2.unit == "Hz/T" and g2.value == pytest.approx(1.76085963e11 / (2 * math.pi))
    assert units.gamma_from_gamma_over_2pi(g2).value == pytest.approx(g.value)
    with pytest.raises(units.UnitError):
        units.gamma(1.0, "Hz/T")


def test_unit_is_mandatory():
    with pytest.raises(units.UnitError):
        S.ParameterEstimate("Ms", 1.0, "", "magnetization", "dft", S.Provenance("x"))
    with pytest.raises(units.UnitError):
        S.ParameterEstimate("B", 1.0, "Hz/T", "B_hf", "dft", S.Provenance("x"))
    with pytest.raises(ValueError):
        S.ParameterEstimate("Ms", 1.0, "A/m", "magnetization", "made_up_role", S.Provenance("x"))


def test_source_roles_exact():
    assert {r.value for r in S.SourceRole} == {"literature_fit", "literature_holdout", "dft",
                                               "laboratory_validation", "legacy_prediction", "explicit_calibration",
                                               # stage 3: never usable in a parametrization (tests/test_stage3_schemas)
                                               "excluded_secondary", "excluded_low_quality", "methodology_only",
                                               "assumption"}


def test_yaml_json_round_trip():
    p = S.ParameterEstimate("gamma_over_2pi", 28.0e9, "Hz/T", "gamma_over_2pi", "literature_fit",
                            S.Provenance("L01", "Table 2"), S.Uncertainty("interval", lower=27e9, upper=29e9))
    d = yaml.safe_load(S.dump_yaml(p))
    assert S.from_dict(S.ParameterEstimate, d) == p
    assert S.from_dict(S.ParameterEstimate, json.loads(S.dump_json(p))) == p
    pr = S.Prediction("p", "mv", [S.DataRef("a", source_role="dft")],
                      [S.OutputValue("f", 1e10, "Hz", "frequency")])
    assert S.from_dict(S.Prediction, json.loads(S.dump_json(pr))) == pr


def test_other_entities_construct():
    S.RunManifest("r", "audit", "2026-10-25T00:00:00Z")
    S.LiteratureRecord("L01", "Morch 2019", claims=[S.Claim("Ms", 1.0, "A/m", "magnetization")])
    with pytest.raises(ValueError):
        S.Claim("number without unit", 1.0)
    S.SampleRecord("SMP-X", "BAT-X", "x")
    S.Observation("OBS-X", None, "XRD", [])
    S.StructureCandidate("c", "SrFe18O27", [S.SiteOccupancy("12k", {"Fe": 1.0})])
    S.ValidationResult("r", "mv", "p", [S.DataRef("t")], S.Metric("e"))
    with pytest.raises(ValueError):
        S.Metric("e", 0.1)   # value without unit
