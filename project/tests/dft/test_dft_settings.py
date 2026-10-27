import copy

import pytest

from hexw.dft import settings as ST


def test_default_config_loads_and_matches_literature(cfg):
    assert cfg.campaign == "campaign_v2" and len(cfg.sha256) == 64
    lit = cfg.literature["U_eff_Fe_eV"]
    assert lit["rec_id"] == "L05-001" and lit["value"] == 3.0 and lit["lit_id"] == "L05"
    assert cfg.hubbard_U("Fe") == 3.0 and cfg.hubbard_U("Zn") is None
    # Co, Mn are assumptions, never labelled as literature
    st = cfg.raw["dft"]["hubbard"]["U_eff_status"]
    assert st["Fe"] == "literature" and st["Co"] == st["Mn"] == "assumption" and st["Zn"] == "not_applied"
    assert cfg.raw["dft"]["occupations"]["tot_magnetization"] is None
    assert cfg.executable() == "pw_big.exe" and cfg.raw["engine"]["omp_threads"] == 1


def test_hash_is_stable_and_sensitive(cfg):
    again = ST.load_config()
    assert again.sha256 == cfg.sha256
    other = cfg.with_overrides(runs_root="elsewhere")
    assert other.sha256 != cfg.sha256


def test_u_value_must_match_the_literature_record(tmp_path, cfg):
    raw = copy.deepcopy(cfg.raw)
    raw["dft"]["hubbard"]["U_eff_eV"]["Fe"] = 4.0
    import yaml
    p = tmp_path / "c.yaml"
    p.write_text(yaml.safe_dump(raw), encoding="utf-8")
    with pytest.raises(ST.DftConfigError, match="L05-001"):
        ST.load_config(p)
    raw = copy.deepcopy(cfg.raw)
    raw["literature_basis"]["U_eff_Fe_eV"]["value"] = 5.0
    p.write_text(yaml.safe_dump(raw), encoding="utf-8")
    with pytest.raises(ST.DftConfigError, match="L05-001"):
        ST.load_config(p)


@pytest.mark.parametrize("path,value,msg", [
    (("engine", "omp_threads"), 4, "omp_threads"),
    (("dft", "nspin"), 1, "nspin"),
    (("dft", "hubbard", "projector"), "wrong", "projector"),
    (("dft", "hubbard", "U_eff_status", "Co"), "literature", "assumption"),
    (("dft", "hubbard", "U_eff_eV", "Mn"), 40.0, "out of range"),
    (("exchange", "n_holdout"), 0, "hold-out"),
    (("exchange", "group_by"), "bogus", "group_by"),
    (("soc", "allow_run"), "yes", "boolean"),
])
def test_bad_values_are_refused(cfg, path, value, msg):
    raw = copy.deepcopy(cfg.raw)
    d = raw
    for k in path[:-1]:
        d = d[k]
    d[path[-1]] = value
    with pytest.raises(ST.DftConfigError, match=msg):
        ST._validate(raw)


def test_pseudo_profiles(cfg):
    pil = ST.load_profile(cfg, "gbrv_pilot")
    assert (pil.ecutwfc, pil.ecutrho) == (40.0, 240.0) and pil.cutoff_status == "pilot"
    ss = ST.load_profile(cfg, "sssp_efficiency", ("Sr", "Fe", "O"))
    assert (ss.ecutwfc, ss.ecutrho) == (90.0, 1080.0)                 # Fe PAW dominates (SSSP 1.3.0 efficiency metadata)
    assert "not_feasible" in ss.cutoff_status
    mn = ST.load_profile(cfg, "sssp_efficiency", ("Sr", "Fe", "O", "Mn", "Co", "Zn"))
    assert mn.ecutwfc == 90.0 and "Mn 65/780" in mn.cutoff_basis
    for p in (pil, ss):
        for el, f in p.files.items():
            assert (p.directory / f).is_file(), (p.name, el, f)
    with pytest.raises(ST.DftConfigError):
        ST.load_profile(cfg, "nope")


def test_upf_valence_reads_both_formats(cfg):
    pil = ST.load_profile(cfg, "gbrv_pilot")
    assert ST.upf_valence(pil.directory / pil.files["Fe"]) == pytest.approx(16.0)          # GBRV: 3s 3p 3d 4s
    assert ST.upf_valence(pil.directory / pil.files["O"]) == pytest.approx(6.0)
    ss = ST.load_profile(cfg, "sssp_efficiency")
    assert ST.upf_valence(ss.directory / ss.files["Fe"]) == pytest.approx(16.0)             # UPF v2 attribute


def test_system_overrides(cfg):
    assert cfg.exchange("W92")["n_holdout"] == 2
    assert cfg.exchange("hematite")["group_by"] == "shell" and cfg.exchange("hematite")["n_holdout"] == 1
    with pytest.raises(ST.DftConfigError):
        cfg.system("nothing")
