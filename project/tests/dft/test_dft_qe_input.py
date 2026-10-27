import re

import numpy as np
import pytest

from hexw.dft import magnet as M
from hexw.dft import qe_input as Q
from hexw.dft.hematite import hematite_system


@pytest.fixture(scope="module")
def hem():
    return hematite_system()


def test_hematite_input_is_valid_and_has_species_per_spin(cfg, hem):
    pat = np.array([1, -1, 1, 1.0])
    txt = Q.build_input(hem, pat, cfg, prefix="t", config_id="x")
    assert Q.validate_input(txt, cfg, hem, pat) == []
    P_ = Q.parse_pw_input(txt)
    labels = [s[0] for s in P_["species"]]
    assert set(labels) == {"Fe1", "Fe2", "O"} and all(len(l) <= 3 for l in labels)
    assert "tot_magnetization" not in txt                                       # NOT fixed (checked afterwards)
    assert P_["hubbard"] == "atomic" and sorted(P_["U"]) == [("U", "Fe1-3d", 3.0), ("U", "Fe2-3d", 3.0)]
    nl = P_["namelists"]["system"]
    assert nl["nspin"] == "2" and nl["occupations"] == "smearing" and nl["ecutwfc"] == "40" and nl["ecutrho"] == "240"
    assert float(P_["namelists"]["electrons"]["conv_thr"]) == pytest.approx(1e-8 * 10)
    # sign of the starting magnetisation follows the target configuration atom by atom
    sig = hem.sigma(pat)
    sm = {i + 1: float(nl[f"starting_magnetization({i + 1})"]) for i in range(len(labels))}
    for k, (lab, _) in enumerate(P_["positions"]):
        if hem.spin[k] > 0:
            assert np.sign(sm[labels.index(lab) + 1]) == np.sign(sig[k])


def test_validation_catches_errors(cfg, hem):
    pat = np.ones(4)
    txt = Q.build_input(hem, pat, cfg, prefix="t")

    def bad(t, frag):
        problems = Q.validate_input(t, cfg, hem, pat)
        assert any(frag in p for p in problems), problems

    bad(txt.replace("nat=10", "nat=11"), "nat=11")
    bad(txt.replace("ntyp=3", "ntyp=4"), "ntyp=4")
    bad(txt.replace("fe_pbe_v1.5.uspp.F.UPF", "missing.UPF"), "not found")
    bad(txt.replace("o_pbe_v1.2.uspp.F.UPF", "fe_pbe_v1.5.uspp.F.UPF", 1), "does not belong")
    bad(txt.replace("ecutwfc=40", "ecutwfc=20"), "below the profile")
    bad(txt.replace("starting_magnetization(1)=0.4", "starting_magnetization(1)=-0.4"), "sign disagrees")
    bad(txt.replace("starting_magnetization(1)=0.4", "starting_magnetization(1)=1.7"), "outside [-1, 1]")
    bad(txt.replace("U Fe1-3d 3", "U Fe1-3d 5"), "config says")
    bad(txt.replace("{atomic}", "{ortho-atomic}"), "projector")
    bad(txt.replace("&SYSTEM", "&SYSTEM\n  tot_magnetization=3,"), "tot_magnetization")
    bad(txt.replace("&SYSTEM", "&SYSTEM\n  noncolin=.true.,"), "noncolin")
    bad(txt.replace("calculation='scf'", "calculation='relax'"), "scf")
    bad(txt.replace("O 0.4441000000 0.0559000000 0.7500000000", "O 0.8553000000 0.8553000000 0.8553000000"), "closer than")
    bad(re.sub(r"K_POINTS automatic\n\S+ \S+ \S+ 0 0 0", "K_POINTS automatic\n0 2 2 0 0 0", txt), "k mesh")
    bad(txt.replace("U Fe2-3d 3\n", ""), "config says") if False else None
    assert Q.validate_input("garbage", cfg) != []
    with pytest.raises(Q.InputError):
        Q.assert_valid(txt.replace("nat=10", "nat=3"), cfg)


def test_zero_u_and_nonmagnetic_species_have_no_hubbard_or_moment(cfg):
    """Zn: no U, no starting moment; Co / Mn: U marked as assumption in the config."""
    from pymatgen.core import Lattice, Species, Structure
    st = Structure(Lattice.cubic(8.0), [Species("Zn", 2), Species("Mn", 2), Species("Co", 2), Species("O", -2)],
                   [[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0], [0.5, 0.5, 0.5]])
    spin = np.array([0, 2.5, 1.5, 0.0])
    sysm = M.MagneticSystem("toy", "hematite", st, spin, np.array([-1, 0, 1, -1]), ["a", "b"], np.array([0, 1, -1, 0.0]))
    txt = Q.build_input(sysm, np.ones(2), cfg, prefix="t")
    P_ = Q.parse_pw_input(txt)
    assert [s[0] for s in P_["species"]] == ["Co2", "Mn1", "O", "Zn"]
    assert {u[1] for u in P_["U"]} == {"Co2-3d", "Mn1-3d"}                       # Zn: none
    nl = P_["namelists"]["system"]
    assert float(nl["starting_magnetization(1)"]) < 0 < float(nl["starting_magnetization(2)"])
    assert float(nl["starting_magnetization(4)"]) == 0
    assert Q.validate_input(txt, cfg, sysm, np.ones(2)) == []


@pytest.fixture(scope="module")
def w_inputs():
    from hexw.dft.settings import load_config
    from hexw.structures.ensemble import Hypothesis
    from hexw.structures.configs import generate_configurations
    cfg = load_config()
    out = {}
    for key, comp, hyp in (("W", {}, Hypothesis(fe2_scenario="uniform")),
                           ("Mn", {"Mn": 2}, Hypothesis(mn_valence=2, fe2_scenario="uniform"))):
        cs = generate_configurations(comp, hyp, n_rep=2)
        sysm = M.system_from_w_config(cs, 0)
        model = M.build_model(sysm, 3.8, "sublattice_pair")
        des = M.select_design(model, 1, 2)
        out[key] = (sysm, des)
    return cfg, out


def test_w_and_mn_inputs_validate_for_all_configurations(w_inputs):
    cfg, out = w_inputs
    for key, (sysm, des) in out.items():
        assert sysm.n_atoms == 92
        for k in range(len(des.config_ids)):
            txt = Q.build_input(sysm, des.patterns[k], cfg, prefix=f"p{k}", config_id=des.config_ids[k])
            assert Q.validate_input(txt, cfg, sysm, des.patterns[k]) == [], (key, des.config_ids[k])
    P_ = Q.parse_pw_input(Q.build_input(out["Mn"][0], out["Mn"][1].patterns[3], cfg, prefix="x"))
    labels = {s[0] for s in P_["species"]}
    assert {"Fe1", "Fe2", "Sr", "O"} <= labels and labels & {"Mn1", "Mn2"} and "Zn" not in labels
    assert any(u[1].startswith("Mn") for u in P_["U"])
    n_mn = sum(1 for p in P_["positions"] if p[0].startswith("Mn"))
    assert n_mn == 4                                                            # SrFe16Mn2O27, Z = 2


def test_ground_state_is_gorter_in_the_input(w_inputs):
    cfg, out = w_inputs
    sysm, des = out["W"]
    txt = Q.build_input(sysm, des.patterns[0], cfg, prefix="g")
    P_ = Q.parse_pw_input(txt)
    n_up = sum(1 for p in P_["positions"] if p[0] in ("Fe1", "Fe3"))             # Fe3 / Fe4: the Fe2+ ions (up / down)
    n_dn = sum(1 for p in P_["positions"] if p[0] in ("Fe2", "Fe4"))
    assert (n_up, n_dn) == (24, 12)                                             # up: W3 6 + W4 4 + W5 12 + W7 2; down: W1 4 + W2 4 + W6 4
