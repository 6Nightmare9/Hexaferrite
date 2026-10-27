import pytest

from hexw.registry.schemas import (DataRef, Metric, ParameterEstimate, Prediction, Provenance, SourceRole,
                                   Uncertainty, ValidationResult)
from hexw.registry.versioning import (AppendOnlyError, LeakageError, ValidationStore, check_inputs_targets_disjoint,
                                      check_prediction_inputs_roles, make_model_version, parameter_set_hash)


def P(name="Ms", v=1.0, role=SourceRole.DFT, **kw):
    return ParameterEstimate(name, v, "A/m", "magnetization", role, Provenance("src"), **kw)


def pred(inputs, mv="mv-1"):
    return Prediction("p1", mv, inputs)


def vr(targets, rid="r1", mv="mv-1"):
    return ValidationResult(rid, mv, "p1", targets, Metric("rel_err", 0.1, "1"))


def test_disjoint_ok():
    check_inputs_targets_disjoint(pred([DataRef("passport:x"), DataRef("L01")]), vr([DataRef("OBS-XRD-001")]))


@pytest.mark.parametrize("inp,tgt", [
    (DataRef("OBS-1"), DataRef("OBS-1")),                                   # same id
    (DataRef("a.mdi", sha256="aa"), DataRef("b.txt", sha256="aa")),         # same bytes
    (DataRef("a.mdi", group="OBS-1"), DataRef("b.txt", group="OBS-1")),     # .mdi vs .txt of one observation
    (DataRef("OBS-1"), DataRef("b.txt", group="OBS-1")),                    # id vs group
])
def test_overlap_raises(inp, tgt):
    with pytest.raises(LeakageError, match="overlap"):
        check_inputs_targets_disjoint(pred([inp]), vr([tgt]))


def test_different_model_versions_rejected():
    with pytest.raises(LeakageError):
        check_inputs_targets_disjoint(pred([DataRef("a")], "mv-1"), vr([DataRef("b")], mv="mv-2"))


def test_laboratory_inputs_rejected():
    with pytest.raises(LeakageError):
        check_prediction_inputs_roles(pred([DataRef("x", source_role="laboratory_validation")]))
    with pytest.raises(LeakageError):
        check_prediction_inputs_roles(pred([DataRef("x", source_role="legacy_prediction")]))


def test_version_changes_with_parameter():
    v1 = make_model_version([P("Ms", 1.0), P("K1", 2.0)])
    assert make_model_version([P("K1", 2.0), P("Ms", 1.0)]).version_id == v1.version_id   # order-independent
    assert make_model_version([P("Ms", 1.0000001), P("K1", 2.0)]).version_id != v1.version_id
    unit_changed = ParameterEstimate("Ms", 1.0, "kA/m", "magnetization", SourceRole.DFT, Provenance("src"))
    assert make_model_version([unit_changed, P("K1", 2.0)]).version_id != v1.version_id
    prov_changed = ParameterEstimate("Ms", 1.0, "A/m", "magnetization", SourceRole.DFT, Provenance("other"))
    assert make_model_version([prov_changed, P("K1", 2.0)]).version_id != v1.version_id
    assert parameter_set_hash([P("Ms", 1.0), P("K1", 2.0)]) == v1.param_hash


def test_lab_and_legacy_params_cannot_define_blind_version():
    with pytest.raises(LeakageError):
        make_model_version([P(role=SourceRole.LABORATORY_VALIDATION)])
    with pytest.raises(LeakageError):
        make_model_version([P(role=SourceRole.LEGACY_PREDICTION)])
    with pytest.raises(LeakageError):
        make_model_version([P(role=SourceRole.EXPLICIT_CALIBRATION)])
    assert make_model_version([P(role=SourceRole.EXPLICIT_CALIBRATION)], allow_explicit_calibration=True).allows_calibration


def test_nonidentifiable_parameter_cannot_get_narrow_interval():
    with pytest.raises(ValueError):
        P(v=10.0, identifiable=False, uncertainty=Uncertainty("std", std=0.01))


def test_validation_store_append_only(tmp_path):
    st = ValidationStore(tmp_path / "v.jsonl")
    v1 = make_model_version([P(v=1.0)])
    v2 = make_model_version([P(v=1.1)])
    assert v1.version_id != v2.version_id
    r1 = vr([DataRef("OBS-9")], "r1", v1.version_id)
    st.append(r1, pred([DataRef("fit-input")], v1.version_id))
    snapshot = (tmp_path / "v.jsonl").read_text(encoding="utf-8")
    with pytest.raises(AppendOnlyError):
        st.append(r1, pred([DataRef("fit-input")], v1.version_id))
    with pytest.raises(AppendOnlyError):                       # same id, different content: still refused
        st.append(vr([DataRef("OBS-10")], "r1", v1.version_id), pred([DataRef("fit-input")], v1.version_id))
    st.append(vr([DataRef("OBS-9")], "r2", v2.version_id), pred([DataRef("fit-input")], v2.version_id))
    assert (tmp_path / "v.jsonl").read_text(encoding="utf-8").startswith(snapshot)
    assert [r["result_id"] for r in st.for_model(v1.version_id)] == ["r1"]
    assert st.verify()
    f = tmp_path / "v.jsonl"                                   # tampering is detected
    f.write_text(f.read_text(encoding="utf-8").replace("OBS-9", "OBS-8", 1), encoding="utf-8")
    assert not st.verify()


def test_store_refuses_overlapping_prediction(tmp_path):
    st = ValidationStore(tmp_path / "v.jsonl")
    with pytest.raises(LeakageError):
        st.append(vr([DataRef("OBS-1")]), prediction=pred([DataRef("OBS-1")]))
    assert not (tmp_path / "v.jsonl").exists()
