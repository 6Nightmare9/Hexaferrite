import numpy as np
import pytest

from hexw.dft import hamiltonian as H
from hexw.dft import magnet as M
from hexw.dft.hematite import hematite_system


@pytest.fixture(scope="module")
def hem():
    return hematite_system()


def test_hematite_system(hem):
    assert hem.n_atoms == 10 and (hem.spin > 0).sum() == 4 and set(hem.spin[hem.spin > 0]) == {2.5}
    assert list(hem.base_sign[hem.spin > 0]) == [1, -1, -1, 1] or sorted(hem.base_sign[hem.spin > 0]) == [-1, -1, 1, 1]
    assert hem.sigma(np.ones(4))[hem.spin > 0].sum() == 0                       # AFM ground state: zero net moment
    assert hem.meta["geometry"]["a_A"] == 5.038


def test_hematite_bonds_and_shells(hem):
    b = M.find_bonds(hem, 3.8)
    d = sorted({round(x.distance, 2) for x in b})
    assert d == [2.90, 2.97, 3.36, 3.71]
    by = {r: [x for x in b if round(x.distance, 2) == r] for r in d}
    assert {x.n_shared_O for x in by[2.90]} == {3} and {x.n_shared_O for x in by[2.97]} == {2}      # face, edge sharing
    assert {x.n_shared_O for x in by[3.36]} == {1} and {x.n_shared_O for x in by[3.71]} == {1}      # corner sharing
    assert 100 < np.mean([x.angle_deg for x in by[3.71]]) < 180
    assert len(b) == 2 + 6 + 6 + 12 and all(x.i < x.j for x in b)


def test_hematite_model_merges_proportional_columns(hem):
    m = M.build_model(hem, 3.8, "shell")
    labels = m.class_labels
    assert len(labels) == 3 and any("+" in l for l in labels)                   # 2.90 and 3.36 A join the same Fe pair
    merged = next(c for c in m.class_info if "combination" in c)
    assert "3*J[" in merged["combination"]
    assert any("proportional" in a["reason"] for a in m.absorbed)
    assert H.conditioning(m.pool_X)["rank"] == 4


def test_design_needs_distinct_rows_for_holdout(hem):
    m = M.build_model(hem, 3.8, "shell")
    with pytest.raises(M.DesignError, match="distinct design rows"):
        M.select_design(m, 1, 1)                                                # 3 J + E0 + 1 + hold-out > 5 distinct rows
    d3 = M.select_design(m, 0, 1)                                               # exactly determined + one new row for the hold-out
    assert d3.roles.count("fit") == 4 and d3.roles.count("holdout") == 1 and d3.cond_fit["rank"] == 4
    m2 = M.build_model(hem, 3.5, "shell")
    d = M.select_design(m2, 1, 1)
    assert d.roles.count("fit") == 4 and d.roles.count("holdout") == 1 and d.config_ids[0] == "c00_ground"
    fit_rows = d.X[d.fit_rows]
    hold = d.X[d.holdout_rows][0]
    assert not any(np.allclose(hold, r) for r in fit_rows)                      # hold-out is a new design row
    assert d.cond_fit["cond_scaled"] < 5 and d.cond_fit["rank"] == 3
    assert len(set(d.config_ids)) == len(d.config_ids)


def test_max_pairs_keeps_shortest(hem):
    m = M.build_model(hem, 3.8, "shell", max_pairs=2)
    assert m.n_classes == 2
    assert any("neglected" in a["reason"] for a in m.absorbed)


# ---- the W phase (92 atoms) -------------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def w_system():
    from hexw.structures.configs import generate_configurations
    cs = generate_configurations({}, n_rep=2)
    return M.system_from_w_config(cs, 0)


@pytest.fixture(scope="module")
def w_design(w_system):
    m = M.build_model(w_system, 3.8, "sublattice_pair")
    return M.select_design(m, 1, 2)


def test_w_system_is_gorter(w_system):
    s = w_system
    assert s.n_atoms == 92 and s.group_names == ["W1", "W2", "W3", "W4", "W5", "W6", "W7"]
    mag = s.spin > 0
    assert mag.sum() == 36
    sig = s.sigma(np.ones(7))
    up = {int(g) for g, sg in zip(s.group[mag], sig[mag]) if sg > 0}
    dn = {int(g) for g, sg in zip(s.group[mag], sig[mag]) if sg < 0}
    assert up == {2, 3, 4, 6} and dn == {0, 1, 5}                              # W3 W4 W5 W7 up; W1 W2 W6 down
    assert set(s.spin[mag]) == {2.5, 2.0}                                       # Fe3+ and the Fe2+ placeholder


def test_w_exchange_classes(w_design):
    m = w_design.model
    assert m.n_classes == 13                                                    # inter-sublattice pairs within 3.8 A
    labels = set(m.class_labels)
    assert {"W3-W4", "W4-W5", "W1-W3", "W5-W6", "W6-W7", "W2-W6"} <= labels
    assert not any(l in labels for l in ("W1-W1", "W3-W3", "W5-W5", "W6-W6"))
    assert {a["label"] for a in m.absorbed} == {"W1-W1", "W3-W3", "W5-W5", "W6-W6"}
    info = {c["label"]: c for c in m.class_info}
    assert "edge (2 O)" in info["W3-W4"]["paths"] and "corner (1 O)" in info["W1-W3"]["paths"]
    assert 100 < info["W1-W3"]["mean_Fe_O_Fe_angle_corner_deg"] < 150


def test_w_design_is_well_conditioned_with_holdout(w_design):
    d = w_design
    assert d.roles.count("fit") == 15 and d.roles.count("holdout") == 2          # n_J + 1 + 1 extra, 2 hold-out
    assert d.cond_fit["rank"] == 14 and d.cond_fit["cond_scaled"] < 10
    assert d.config_ids[0] == "c00_ground" and len(set(d.config_ids)) == 17
    assert d.pool_cond["cond_scaled"] == pytest.approx(1.0)                     # distinct characters are orthogonal
    sv = d.cond_fit["singular_values_scaled"]
    assert len(sv) == 14 and min(sv) > 0
    fit_rows = d.X[d.fit_rows]
    for h in d.X[d.holdout_rows]:
        assert not any(np.allclose(h, r) for r in fit_rows)
    # every fit pattern is a flip pattern relative to the Gorter ground state
    assert set(np.unique(d.patterns)) == {-1.0, 1.0} and np.all(d.patterns[0] == 1)


def test_w_design_reproduces_energies_of_a_known_model(w_design):
    """End-to-end identifiability: synthetic Heisenberg energies -> J recovered from the chosen configurations only."""
    d = w_design
    rng = np.random.default_rng(0)
    J = rng.normal(0, 3.0, d.model.n_classes)
    E = 500.0 + d.X @ J
    fit = H.fit_exchange(d.X[d.fit_rows], E[d.fit_rows], d.model.class_labels, noise_floor_meV=0.1,
                         X_holdout=d.X[d.holdout_rows], E_holdout_meV=E[d.holdout_rows])
    assert np.allclose(fit.J, J, atol=1e-8) and fit.holdout["max_abs_meV"] < 1e-8


@pytest.fixture(scope="module")
def w_far(w_system):
    m = M.build_model(w_system, 3.8, "sublattice_pair", None, 8.0)
    return m, M.select_design(m, 1, 2)


def test_w_pairs_without_a_short_bond_become_nuisance_columns(w_far):
    m, d = w_far
    assert m.far_labels == ["far:W1-W6", "far:W1-W7", "far:W2-W7", "far:W3-W5", "far:W3-W6", "far:W4-W6", "far:W4-W7"]
    assert m.pool_X_far.shape == (64, 7) and d.X_far.shape == (17, 7)
    assert all(i["n_bonds_far"] > 0 and i["d_min"] > 3.8 for i in m.far_info)
    assert m.n_classes == 13                                                    # the fitted classes are unchanged


def test_effective_j_of_a_class_absorbs_all_pairs_of_the_two_sublattices(w_system, w_far):
    """Random J for EVERY bond below 8 A: the energy of any sublattice-flip pattern is E0 + X J_eff + X_far J_far, with
    J_eff = sum over all pairs of the two sublattices of J_ij S_i S_j / sum of S_i S_j of the bonds < 3.8 A (docs/DFT.md section 3)."""
    m, _ = w_far
    sysm = w_system
    bonds = M.find_bonds(sysm, 8.0)
    rng = np.random.default_rng(5)
    g = lambda b: tuple(sorted((int(sysm.group[b.i]), int(sysm.group[b.j]))))
    types = sorted({g(b) for b in bonds if g(b)[0] != g(b)[1]})
    near_types = {g(b) for b in bonds if b.distance <= 3.8 + 1e-6}
    Jt = {t: rng.normal(0, 3.0) for t in types}
    Jb = [(b, Jt[g(b)] + rng.normal(0, 0.5)) if g(b) in near_types else (b, Jt[g(b)]) for b in bonds if g(b)[0] != g(b)[1]]
    pats = m.pool_patterns
    sig = np.array([sysm.sigma(p) for p in pats])
    E = np.array([-sum(j * b.S_i * b.S_j * s[b.i] * s[b.j] for b, j in Jb) for s in sig])
    eff = {}
    for t in near_types:
        if t[0] == t[1]:
            continue
        num = sum(j * b.S_i * b.S_j for b, j in Jb if g(b) == t)
        den = sum(b.S_i * b.S_j for b, j in Jb if g(b) == t and b.distance <= 3.8 + 1e-6)
        eff[t] = num / den
    Jeff = np.array([eff[tuple(int(x[1:]) - 1 for x in lab.split("-"))] for lab in m.class_labels])
    Jfar = np.array([Jt[tuple(int(x[1:]) - 1 for x in lab[4:].split("-"))] for lab in m.far_labels])
    pred = m.pool_X @ Jeff + m.pool_X_far @ Jfar
    assert np.allclose(E - E[0], pred - pred[0], atol=1e-8)
