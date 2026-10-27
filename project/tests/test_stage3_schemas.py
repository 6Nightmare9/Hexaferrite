"""Stage 3: units, roles, literature trace, uncertainty conventions, site IDs and moment definitions."""
import math

import pytest

from hexw import units as U
from hexw.registry import schemas as S

P = S.Provenance


def test_new_units():
    assert U.Quantity(10.0, "emu/g", U.Kind.SPECIFIC_MAGNETIZATION).canonical().value == 10.0
    assert U.Quantity(1.0, "erg/cm^3", U.Kind.ENERGY_DENSITY).canonical().value == pytest.approx(0.1)
    assert U.Quantity(1.0, "MJ/m^3", U.Kind.ENERGY_DENSITY).canonical().value == 1e6
    assert U.Quantity(1.0, "atm", U.Kind.PRESSURE).canonical().value == 101325.0
    assert U.Quantity(2.0, "MA/m", U.Kind.MAGNETIZATION).canonical().value == 2e6
    assert U.Quantity(3.0, "kA/m", U.Kind.H_FIELD).canonical().value == 3e3
    for k, u in ((U.Kind.MAGNETIC_MOMENT, "muB"), (U.Kind.VOLUME, "angstrom^3"), (U.Kind.DENSITY, "g/cm^3"),
                 (U.Kind.ENERGY_PER_ATOM, "eV/atom"), (U.Kind.ENERGY_PER_ATOM, "meV/atom"),
                 (U.Kind.ENERGY_PER_FORMULA_UNIT, "eV/f.u."), (U.Kind.ENERGY_PER_FORMULA_UNIT, "meV/f.u."),
                 (U.Kind.MU0M, "T"), (U.Kind.REFLECTION_LOSS, "dB")):
        U.check_unit(k, u)


def test_kinds_do_not_mix():
    e_atom = U.Quantity(1.0, "eV/atom", U.Kind.ENERGY_PER_ATOM)
    e_fu = U.Quantity(46.0, "eV/f.u.", U.Kind.ENERGY_PER_FORMULA_UNIT)
    with pytest.raises(U.UnitError):
        U.Quantity(1.0, "eV/f.u.", U.Kind.ENERGY_PER_ATOM)
    assert U.energy_per_fu_from_per_atom(e_atom, 46).canonical().value == 46.0
    assert U.energy_per_atom_from_per_fu(e_fu, 46).canonical().value == pytest.approx(1.0)
    with pytest.raises(U.UnitError):
        U.energy_per_fu_from_per_atom(e_fu, 46)
    with pytest.raises(U.UnitError):                       # mass fraction is not a spectral-area fraction
        U.Quantity(5.0, "area%", U.Kind.FRACTION_MASS)


def test_specific_magnetization_roundtrip():
    s = U.specific_magnetization(30.0, "emu/g")
    rho = U.Quantity(5.1, "g/cm^3", U.Kind.DENSITY)
    M = U.M_from_specific(s, rho)
    assert M.value == pytest.approx(30.0 * 5100.0)
    assert U.specific_from_M(M, rho).value == pytest.approx(30.0)


def test_legacy_constants_kept():
    for n in ("MU0", "MU_B", "HBAR", "gamma_from_g"):
        assert hasattr(U, n)


@pytest.mark.parametrize("role", ["excluded_secondary", "methodology_only", "excluded_low_quality"])
def test_non_parametric_roles_cannot_parametrize(role):
    assert S.SourceRole(role) in S.NON_PARAMETRIC_ROLES
    with pytest.raises(ValueError, match="cannot be used in a parametrization"):
        S.ParameterEstimate("Ms", 1.0, "A/m", "magnetization", role, P("x"))


def _src(**kw):
    d = dict(locator=S.Locator(table="Table 1"), data_status="measured", tier=1, doi="10.1/abc")
    d.update(kw)
    return S.LiteratureSource(**d)


def _est(role="literature_fit"):
    return S.ParameterEstimate("a", 5.9, "angstrom", "length", role, P("lit"))


def test_claim_role_split_consistency_and_excluded_no_numbers():
    c = S.LiteratureClaim("R1", "L02", "a", "literature_fit", _src(), estimate=_est())
    assert c.split == "fit" and c.estimate.provenance.literature is c.source
    with pytest.raises(ValueError):
        S.LiteratureClaim("R1", "L02", "a", "literature_fit", _src(), estimate=_est(), split="holdout")
    with pytest.raises(ValueError):
        S.LiteratureClaim("R2", "C1", "a", "excluded_secondary", _src(), qualitative="x",
                          estimate=S.ParameterEstimate("a", 1.0, "angstrom", "length", "literature_fit", P("l")))


def test_source_trace_rules():
    with pytest.raises(ValueError):
        _src(tier=5)
    with pytest.raises(ValueError):
        _src(figure_only=True)                               # axes required
    with pytest.raises(ValueError):
        _src(data_status="digitized")                        # digitized -> figure_only
    with pytest.raises(ValueError):
        _src(transfer="Ba_to_Sr", tier=2)                    # transfer needs tier >= 3
    assert _src(transfer="Ba_to_Sr", tier=3).transfer == "Ba_to_Sr"
    with pytest.raises(ValueError):
        S.Locator()
    with pytest.raises(ValueError):
        _src(holdout_kind="nonsense")


def test_uncertainty_conventions():
    v, u = S.parse_reported("5.8955(1)")
    assert v == 5.8955 and u.convention == "parenthetical_se" and u.std == pytest.approx(1e-4)
    v, u = S.parse_reported("60-65")
    assert v == 62.5 and (u.lower, u.upper, u.convention) == (60.0, 65.0, "range")
    v, u = S.parse_reported("≥ 5")
    assert v == 5.0 and u.inequality == "ge"
    v, u = S.parse_reported("<= 5")
    assert u.inequality == "le"
    with pytest.raises(ValueError):
        S.Uncertainty(kind="std", std=0.1, convention="bogus")
    with pytest.raises(ValueError):
        S.Uncertainty(inequality="gt")


def test_old_objects_serialize_unchanged():
    d = S.to_plain(S.Uncertainty())
    assert "convention" not in d and "inequality" not in d
    d = S.to_plain(S.SiteOccupancy("12k", {"Fe": 1.0}))
    assert "moments" not in d and "site_id" not in d


def test_site_ids_and_split():
    assert S.site_id("4f3", "islam2025") == "W6" and S.site_id("4f3", "nakai") == "W4"   # schemes differ!
    with pytest.raises(KeyError):
        S.site_id("4f", "islam2025")                          # ambiguous
    with pytest.raises(ValueError):
        S.SiteOccupancy("12k", {"Fe": 1.0}, site_id="W9")
    with pytest.raises(ValueError):
        S.SiteOccupancy("12k", {"Fe": 1.0}, split_label="8d")   # split needs W-site and space group
    s = S.SiteOccupancy("12k", {"Fe": 1.0}, site_id="W5", space_group="Cmcm", split_label="8d")
    assert s.split_label == "8d"


def test_moment_definitions_not_mixable():
    fe = S.SiteMoment(4.1, "fe_ion")
    avg = S.SiteMoment(3.5, "site_average_incl_diamagnetic")
    co = S.SiteMoment(2.0, "species_ion", species="Co")
    assert fe - S.SiteMoment(3.9, "fe_ion") == pytest.approx(0.2)
    for a, b in ((fe, avg), (fe, co), (avg, co)):
        with pytest.raises(S.MomentDefinitionError):
            a - b
        with pytest.raises(S.MomentDefinitionError):
            a < b
    assert co.to_site_average(0.5).value == pytest.approx(1.0)
    site = S.SiteOccupancy("12k", {"Fe": 0.9, "Co": 0.1}, site_id="W5", moments=[fe, avg])
    assert site.moment("fe_ion") is fe
    with pytest.raises(ValueError):
        S.SiteOccupancy("12k", {"Fe": 1.0}, moments=[fe, S.SiteMoment(3.0, "fe_ion")])


def test_partial_occupancy_flag():
    ordered = S.StructureCandidate("c1", "SrFe18O27", [S.SiteOccupancy("12k", {"Fe": 1.0})])
    ordered.require_dft_configuration()
    flagged = S.StructureCandidate("c2", "x", [S.SiteOccupancy("12k", {"Fe": 1.0})],
                                   partial_occupancy_cif_not_dft_config=True)
    with pytest.raises(S.PartialOccupancyError):
        flagged.require_dft_configuration()
    frac = S.StructureCandidate("c3", "x", [S.SiteOccupancy("12k", {"Fe": 0.8, "Co": 0.2})])
    with pytest.raises(S.PartialOccupancyError):
        frac.require_dft_configuration()


def test_yaml_round_trip_claim():
    c = S.LiteratureClaim("R1", "L11", "a", "literature_fit", _src(transfer="Ba_to_Sr", tier=3),
                          estimate=_est(), value_original="5.9(1)")
    d = S.to_plain(c)
    d["estimate"]["provenance"].pop("literature", None)
    c2 = S.from_dict(S.LiteratureClaim, d)
    assert c2.source.transfer == "Ba_to_Sr" and c2.estimate.provenance.literature is c2.source
