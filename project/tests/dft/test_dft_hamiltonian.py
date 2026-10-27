import math

import numpy as np
import pytest

from hexw.dft import hamiltonian as H


# ---- conventions ------------------------------------------------------------------------------------------------------
def test_convention_is_documented_and_fixed():
    assert H.OURS == H.Convention(-1, "once", "S")
    assert "H = -sum_{i<j}" in H.CONVENTION_TEXT and "J>0 ferromagnetic" in H.CONVENTION_TEXT
    assert H.RY_TO_MEV == pytest.approx(13605.693, abs=1e-3)


def test_spin_lengths():
    assert H.spin_length("Fe", 3) == 2.5 and H.spin_length("Fe", 2) == 2.0
    assert H.spin_length("Mn", 2) == 2.5 and H.spin_length("Co", 2) == 1.5 and H.spin_length("Zn", 2) == 0.0
    with pytest.raises(ValueError):
        H.spin_length("Xx", 2)


def test_convert_unit_vectors_and_double_counting():
    J, Si, Sj = -4.0, 2.5, 2.5
    unit = H.Convention(-1, "once", "unit")
    assert H.convert_J(J, Si, Sj, H.OURS, unit) == pytest.approx(J * Si * Sj)         # J_unit = J S_i S_j
    twice = H.Convention(-1, "twice", "S")
    assert H.convert_J(J, Si, Sj, H.OURS, twice) == pytest.approx(J / 2)              # sum over i != j counts a pair twice
    plus = H.Convention(+1, "once", "S")
    assert H.convert_J(J, Si, Sj, H.OURS, plus) == pytest.approx(-J)                  # H = +sum J: J>0 antiferro
    # the pair energy is the same in every convention
    for c in (unit, twice, plus, H.Convention(1, "twice", "unit")):
        Jc = H.convert_J(J, Si, Sj, H.OURS, c)
        mult = 2.0 if c.counting == "twice" else 1.0
        norm = Si * Sj if c.spin == "S" else 1.0
        assert c.sign * Jc * mult * norm == pytest.approx(-J * Si * Sj)


@pytest.mark.parametrize("src", [H.OURS, H.Convention(1, "twice", "unit"), H.Convention(-1, "twice", "S")])
@pytest.mark.parametrize("dst", [H.OURS, H.Convention(1, "once", "unit"), H.Convention(-1, "once", "unit")])
def test_convert_round_trip(src, dst):
    for J in (-7.3, 0.4, 12.0):
        back = H.convert_J(H.convert_J(J, 2.5, 2.0, src, dst), 2.5, 2.0, dst, src)
        assert back == pytest.approx(J)


def test_convert_rejects_zero_spin_and_bad_convention():
    with pytest.raises(ValueError):
        H.convert_J(1.0, 0.0, 2.5, H.OURS, H.Convention(-1, "once", "S"))
    with pytest.raises(ValueError):
        H.Convention(0, "once", "S")
    with pytest.raises(ValueError):
        H.Convention(-1, "thrice", "S")


def test_sign_ferro_is_lower_energy():
    # J > 0: parallel pair lower than antiparallel in our convention
    assert H.energy_per_pair_parallel(+3.0, 2.5, 2.5) < 0 < H.energy_per_pair_parallel(-3.0, 2.5, 2.5)
    assert H.J_to_kelvin(8.617333262e-2) == pytest.approx(1.0)


# ---- design matrix against a brute-force Heisenberg energy ----------------------------------------------------------------
def _toy(n_atoms=6, n_bonds=14, n_classes=3, seed=1):
    rng = np.random.default_rng(seed)
    S = rng.choice([2.5, 2.0, 1.5], n_atoms)
    bonds = []
    for _ in range(n_bonds):
        i, j = sorted(rng.choice(n_atoms, 2, replace=False))
        bonds.append(H.Bond(int(i), int(j), (0, 0, 0), 3.0, int(rng.integers(n_classes)), float(S[i]), float(S[j])))
    return bonds, n_atoms, n_classes


def test_design_matrix_equals_brute_force_energy():
    bonds, n, K = _toy()
    J = np.array([4.0, -2.0, 0.7])
    rng = np.random.default_rng(5)
    sig = rng.choice([-1.0, 1.0], (9, n))
    X = H.design_matrix(bonds, sig, K)
    for c in range(len(sig)):
        e = sum(-J[b.cls] * b.S_i * b.S_j * sig[c, b.i] * sig[c, b.j] for b in bonds)     # H = -sum J S_i S_j sigma_i sigma_j
        assert X[c] @ J == pytest.approx(e)
        assert np.allclose(H.design_row(bonds, sig[c], K), X[c])


def test_absorbed_bonds_do_not_enter():
    b = H.Bond(0, 1, (0, 0, 0), 3.0, -1, 2.5, 2.5)
    assert np.all(H.design_matrix([b], np.array([[1.0, -1.0]]), 1) == 0)


# ---- fit ------------------------------------------------------------------------------------------------------------------------
def _problem(n_cfg=14, K=4, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n_cfg, K)) * 30
    J = np.array([5.0, -3.0, 1.5, -8.0])[:K]
    E0 = -1234.5
    return X, J, E0


def test_fit_recovers_noise_free_parameters():
    X, J, E0 = _problem()
    E = E0 + X @ J
    fr = H.fit_exchange(X, E, list("abcd"), noise_floor_meV=0.1)
    assert fr.E0 == pytest.approx(E0) and np.allclose(fr.J, J)
    assert fr.dof == 14 - 5 and fr.rms_meV < 1e-9
    assert fr.sigma_basis == "noise_floor" and np.all(fr.sigma_J > 0)
    assert np.nanmax(np.abs(fr.loo_residuals)) < 1e-8


def test_fit_uncertainties_cover_truth_statistically():
    X, J, E0 = _problem(n_cfg=20)
    rng = np.random.default_rng(3)
    z = []
    for _ in range(300):
        E = E0 + X @ J + rng.normal(0, 1.0, len(X))
        fr = H.fit_exchange(X, E, list("abcd"))
        z.append((fr.J - J) / fr.sigma_J)
    z = np.array(z)
    assert abs(z.std() - 1.0) < 0.12 and np.abs(z.mean()) < 0.1          # residual-based sigma is calibrated


def test_fit_holdout_prediction_and_flagging():
    X, J, E0 = _problem(n_cfg=14)
    Xh = np.random.default_rng(9).normal(size=(3, 4)) * 30
    E = E0 + X @ J
    Eh = E0 + Xh @ J
    fr = H.fit_exchange(X, E, list("abcd"), noise_floor_meV=0.2, X_holdout=Xh, E_holdout_meV=Eh)
    assert fr.holdout["max_abs_meV"] < 1e-8
    # a hold-out generated by a different model (a missing term) is badly predicted
    fr2 = H.fit_exchange(X, E, list("abcd"), noise_floor_meV=0.2, X_holdout=Xh, E_holdout_meV=Eh + 25.0)
    assert fr2.holdout["rms_meV"] == pytest.approx(25.0, rel=1e-6)
    assert np.all(np.abs(fr2.holdout["z"]) > 3)
    # a disagreeing hold-out inflates the J uncertainties, a consistent one does not
    assert fr.holdout_consistent and np.allclose(fr.sigma_J_total, fr.sigma_J) and fr.sigma_model_meV == 0
    assert fr2.holdout_consistent is False and fr2.sigma_model_meV == pytest.approx(25.0, rel=1e-6)
    assert np.all(fr2.sigma_J_total > 10 * fr2.sigma_J) and any("hold-out disagrees" in n for n in fr2.notes)


def test_fit_rank_deficient_and_underdetermined_raise():
    X, J, E0 = _problem()
    X2 = np.hstack([X, X[:, :1] * 2.0])                               # proportional column
    with pytest.raises(ValueError, match="rank deficient"):
        H.fit_exchange(X2, E0 + X @ J, list("abcde"))
    with pytest.raises(ValueError):
        H.fit_exchange(X[:3], np.zeros(3), list("abcd"))


def test_exactly_determined_fit_needs_noise_floor():
    X, J, E0 = _problem(n_cfg=5)
    E = E0 + X @ J
    with pytest.raises(ValueError, match="no error scale"):
        H.fit_exchange(X, E, list("abcd"), noise_floor_meV=0.0)
    fr = H.fit_exchange(X, E, list("abcd"), noise_floor_meV=0.5)
    assert fr.dof == 0 and fr.sigma_basis == "noise_floor" and any("dof = 0" in n for n in fr.notes)


def test_conditioning_reports_svd():
    X, _, _ = _problem()
    c = H.conditioning(X)
    assert c["rank"] == 5 and len(c["singular_values_scaled"]) == 5
    assert c["cond_scaled"] == pytest.approx(c["singular_values_scaled"][0] / c["singular_values_scaled"][-1])
    Xd = np.hstack([X, X[:, :1]])
    assert H.conditioning(Xd)["rank"] == 5 and math.isinf(H.conditioning(Xd)["cond_scaled"])
