"""Stage 4 round 1: SourceRole.ASSUMPTION and Provenance.status as an enumeration (hexw.registry.schemas)."""
import pytest

from hexw.registry import schemas as S
from hexw.registry.versioning import make_model_version


def _pe(role, source_id, status="unverified", **kw):
    return S.ParameterEstimate("x", 1.0, "1", "dimensionless", role, S.Provenance(source_id, status=status), **kw)


def test_assumption_role_is_not_leaky_and_not_non_parametric():
    assert S.SourceRole.ASSUMPTION.value == "assumption"
    assert S.SourceRole.ASSUMPTION not in S.LEAKY_ROLES and S.SourceRole.ASSUMPTION not in S.NON_PARAMETRIC_ROLES
    p = _pe("assumption", "ASSUMPTION", "assumption", identifiable=False)
    assert make_model_version([p]).n_params == 1                 # usable in a blind model version


def test_source_id_naming_an_assumption_requires_the_role():
    for sid in ("ASSUMPTION", "ASSUMPTION:T_f", "assumption:spin_only_Fe3+"):
        with pytest.raises(ValueError, match="role must be 'assumption'"):
            _pe("literature_fit", sid, "assumption")
        _pe("assumption", sid, "assumption", identifiable=False)


def test_mixed_records_are_refused_split_them():
    for role in ("literature_fit", "assumption"):
        with pytest.raises(ValueError):
            _pe(role, "ASSUMPTION;L02-036", "assumption" if role == "assumption" else "derived_by_us")
    _pe("literature_fit", "L02-036", "derived_by_us")            # the sourced half
    _pe("assumption", "ASSUMPTION", "assumption", identifiable=False)       # the assumption half


def test_role_assumption_needs_the_id_and_the_status():
    with pytest.raises(ValueError, match="ASSUMPTION"):
        _pe("assumption", "L02-036", "assumption")
    with pytest.raises(ValueError, match="status"):
        _pe("assumption", "ASSUMPTION", "verified")
    with pytest.raises(ValueError, match="requires role"):
        _pe("literature_fit", "L02-036", "assumption")           # an assumption cannot be a literature_fit


def test_provenance_status_is_an_enumeration():
    assert [m.value for m in S.ProvenanceStatus] == ["verified", "unverified", "derived_by_us", "assumption"]
    assert S.Provenance("a").status is S.ProvenanceStatus.UNVERIFIED
    assert S.Provenance("a", status="derived_by_us").status is S.ProvenanceStatus.DERIVED_BY_US
    for bad in ("preliminary_from_filename", "Verified!", ""):
        with pytest.raises(ValueError):
            S.Provenance("a", status=bad)
    p = _pe("literature_fit", "L02-036", S.ProvenanceStatus.VERIFIED)
    d = S.to_plain(p)
    assert d["provenance"]["status"] == "verified"
    q = S.from_dict(S.ParameterEstimate, d)
    assert q.provenance.status is S.ProvenanceStatus.VERIFIED and q == p


def test_modelled_assumptions_carry_the_role():
    """Every parameter table of the stage-4 models: ASSUMPTION <=> role, no assumption claims 'verified'."""
    try:
        from hexw.params import curie, moments, occupancy
        from hexw.params.literature import fit_records_path
        if not fit_records_path().is_file():
            pytest.skip("literature split not built")
    except Exception:
        pytest.skip("no project root")
    ests = list(curie.curie_parameter_estimates(n=1000)) + list(moments.moment_parameter_estimates())
    ests += [occupancy.freezing_temperature_estimate(), occupancy.sigma_sys_estimate()]
    assert ests
    for e in ests:
        is_a = "ASSUMPTION" in e.provenance.source_id.upper()
        assert is_a == (e.source_role is S.SourceRole.ASSUMPTION), e.name
        if is_a:
            assert e.provenance.status is S.ProvenanceStatus.ASSUMPTION
