"""Ordered DFT configurations from the ensemble: integer counts, charge, symmetry-inequivalence, geometry source."""
import numpy as np
import pytest

from hexw.chem import validators as V
from hexw.params import occupancy as O
from hexw.structures import configs as C
from hexw.structures import ensemble as E


@pytest.fixture(scope="module", autouse=True)
def _warm():
    O.default_fit()


def test_geometry_comes_from_the_l02_cif_only():
    assert "L02_Morch2023" in str(C.CIF_REL) and "L01" not in str(C.CIF_REL)
    t = C.load_template()
    assert t.source == str(C.CIF_REL)
    assert len(t.tm_frac) == 36 and np.bincount(t.tm_site).tolist() == [4, 4, 6, 4, 12, 4, 2]
    assert abs(t.lattice.a - 5.908361) < 1e-5 and abs(t.lattice.c - 32.83268) < 1e-5
    # W7: the split 4f pair was merged onto the mirror planes
    assert sorted(np.round(t.tm_frac[t.tm_site == 6][:, 2], 3)) == [0.25, 0.75]


def test_symmetry_group_is_p63mmc():
    perms = C.symmetry_permutations(C.load_template())
    assert len(perms) == 24
    # a permutation of the TM positions that keeps the site type
    t = C.load_template()
    assert all((t.tm_site[p] == t.tm_site).all() for p in perms)


def test_orbit_reduction_burnside():
    """2 Co on the six 6g positions: C(6,2) = 15 placements fall into 3 orbits (sizes 6, 6, 3)."""
    cs = C.generate_configurations({"Co": 1.0}, n_rep=5, n_draws=60)
    assert cs.n_orbits_total == 3 and abs(cs.coverage - 1.0) < 1e-12
    assert sorted(c.orbit_size for c in cs.configs) == [3, 6, 6]
    assert sum(c.orbit_size for c in cs.configs) == 15


def test_configs_are_ordered_charge_neutral_and_have_92_atoms():
    for comp in ({"Zn": 2.0}, {"Co": 1.0, "Zn": 1.0}, {"Co": 2.0}, {"Mn": 1.0}):
        cs = C.generate_configurations(comp, n_draws=60)
        assert C.N_REP_MIN <= len(cs.configs) <= C.N_REP_MAX
        assert abs(sum(c.weight for c in cs.configs) - 1.0) < 1e-12
        assert 0 < cs.coverage <= 1.0 + 1e-12
        x = sum(comp.values())
        for c in cs.configs:
            assert c.n_atoms == 92 and len(c.structure) == 92
            c.candidate.require_dft_configuration()                      # integer occupancies
            assert c.candidate.partial_occupancy_cif_not_dft_config is False
            V.check_site_occupancy(c.candidate)
            counts = E.candidate_counts(c.candidate, z=2.0)
            V.check_charge_neutral(counts, {"Sr": 2.0, "O": -2.0, **c.candidate.oxidation_states})
            for sp, n in comp.items():
                assert abs(counts[sp] - n) < 1e-9                        # integer atoms: n per f.u. x 2
            assert abs(counts.get("Fe2+", 0.0) - (2.0 - x)) < 1e-9
            # capacity: every TM position is occupied by exactly one atom
            assert sum(1 for s in c.candidate.sites) == 36
            assert len(c.structure.site_properties["magmom"]) == 92


def test_representatives_are_symmetry_inequivalent():
    cs = C.generate_configurations({"Zn": 2.0}, n_rep=5, n_draws=60)
    t = C.load_template()
    perms = C.symmetry_permutations(t)
    # rebuild the Zn pattern of every representative and check that no two are related by a symmetry operation
    patterns = []
    for c in cs.configs:
        lab = np.zeros(36, np.int8)
        for j, s in enumerate(c.candidate.sites):
            lab[j] = 1 if "Zn" in s.species else 0
        patterns.append(lab)
    canon = []
    for lab in patterns:
        images = [tuple(np.asarray(lab)[np.argsort(p)]) for p in perms]
        canon.append(min(images))
    assert len(set(canon)) == len(canon)


def test_count_vectors_follow_the_ensemble():
    cs = C.generate_configurations({"Zn": 2.0}, n_draws=100)
    # 4 Zn in the cell: mostly on W1 (4e) and W2, as in the L02 refinement
    cv = cs.configs[0].count_vector["Zn"]
    assert sum(cv) == 4 and cv[0] >= 2 and cv[6] == 0


def test_non_commensurate_content_needs_a_supercell():
    with pytest.raises(C.ConfigError, match="supercell"):
        C.generate_configurations({"Co": 0.25})


def test_n_rep_is_clipped_to_2_5():
    cs = C.generate_configurations({"Zn": 2.0}, n_rep=50, n_draws=40)
    assert len(cs.configs) == 5
    cs = C.generate_configurations({"Zn": 2.0}, n_rep=1, n_draws=40)
    assert len(cs.configs) == 2


def test_fe2_placement_follows_the_scenario():
    cs = C.generate_configurations({"Co": 1.0}, E.Hypothesis(fe2_scenario="oct_R"), n_draws=40)
    for c in cs.configs:
        fe2_sites = {s.site_id for s in c.candidate.sites if "Fe2+" in s.species}
        assert fe2_sites <= {"W5", "W6"} and fe2_sites
        assert sum(1 for s in c.candidate.sites if "Fe2+" in s.species) == 2          # n(Fe2+) = 1 per f.u. x 2
    with pytest.raises(C.ConfigError):
        C.generate_configurations({"Co": 1.0}, E.Hypothesis(fe2_scenario="nope"))


def test_mn_configurations_carry_the_valence_hypothesis():
    for v, n2 in ((2, 1.0), (3, 2.0)):
        cs = C.generate_configurations({"Mn": 1.0}, E.Hypothesis(mn_valence=v), n_draws=40)
        c = cs.configs[0]
        assert c.candidate.oxidation_states["Mn"] == float(v)
        assert abs(E.candidate_counts(c.candidate)["Fe2+"] - n2) < 1e-9
