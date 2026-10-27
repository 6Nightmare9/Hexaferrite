"""Susceptibility: analytic closed form, loss sign, passivity, Kramers-Kronig, linewidth link. SYNTHETIC inputs."""
import math

import numpy as np
import pytest

from hexw.theory.magnet import (MagnetModel, MU0, find_equilibria, field_vector, fwhm, kramers_kronig_real, mu_r,
                                passivity_margin, sphere, susceptibility, thin_disk_perp_c)
from _synth import GAMMA_2PI, GY, SYNTH_K1, SYNTH_K2, SYNTH_K4, SYNTH_MS

MS = SYNTH_MS
G = GY.gamma_rad_per_s_T


def sphere_chi_closed_form(B, f, alpha):
    """Independent derivation: isotropic sphere, field B along z (m0 = z), circular components of the
    exp(-i w t) Gilbert LLG:  chi_xx = (mu0 Ms / 2) gamma [ 1/(gamma B - w(1+i a)) + 1/(gamma B + w(1-i a)) ]."""
    w = 2 * math.pi * f
    return 0.5 * MU0 * MS * G * (1 / (G * B - w * (1 + 1j * alpha)) + 1 / (G * B + w * (1 - 1j * alpha)))


def test_chi_matches_closed_form_sphere():
    m = MagnetModel(Ms=MS, N=sphere())
    B = 0.4
    fr = np.linspace(0.1e9, 30e9, 60)
    chi = susceptibility(m, [0, 0, B], [0, 0, 1], fr, alpha=0.03, gyro=GY)
    ref = np.array([sphere_chi_closed_form(B, f, 0.03) for f in fr])
    assert np.allclose(chi[:, 0, 0], ref, rtol=1e-10, atol=0)
    assert np.allclose(chi[:, 1, 1], ref, rtol=1e-10, atol=0)
    assert np.allclose(chi[:, 2, :], 0) and np.allclose(chi[:, :, 2], 0)       # no response along m0
    # antisymmetric (gyrotropic) part exists and is purely off-diagonal
    assert np.abs(chi[:, 0, 1]).max() > 0.1 * np.abs(chi[:, 0, 0]).max()


def test_static_limit_and_mu_r():
    m = MagnetModel(Ms=MS, N=sphere())
    chi0 = susceptibility(m, [0, 0, 0.4], [0, 0, 1], [0.0], alpha=0.02, gyro=GY)
    assert abs(chi0[0, 0, 0] - MU0 * MS / 0.4) < 1e-10
    mu = mu_r(chi0)
    assert np.allclose(mu[0] - chi0[0], np.eye(3))


def test_loss_sign_convention_exp_minus_i_omega_t():
    m = MagnetModel(Ms=MS, K1=SYNTH_K1, N=sphere())
    f0 = GAMMA_2PI * 2 * SYNTH_K1 / MS
    fr = np.linspace(0.2 * f0, 3 * f0, 400)
    chi = susceptibility(m, np.zeros(3), [0, 0, 1], fr, alpha=0.05, gyro=GY)
    assert np.all(chi[:, 0, 0].imag > 0) and np.all(chi[:, 1, 1].imag > 0)       # chi'' > 0 for omega > 0
    # the same medium with the opposite time convention would have chi'' < 0: complex conjugate flips the sign
    assert np.all(np.conj(chi[:, 0, 0]).imag < 0)


@pytest.mark.parametrize("kind", ["external", "internal"])
@pytest.mark.parametrize("seed", range(4))
def test_passivity_antihermitian_part_psd(seed, kind):
    rng = np.random.default_rng(seed)
    N = np.diag(rng.dirichlet([2, 2, 2]))
    m = MagnetModel(Ms=MS, K1=rng.uniform(-1, 1) * 1e5, K2=SYNTH_K2, K4=SYNTH_K4, N=N)
    B = field_vector(rng.uniform(0.1, 0.8), rng.uniform(0.2, 2.9), rng.uniform(0, 6))
    fr = np.linspace(0.05e9, 80e9, 700)
    for e in find_equilibria(m, B).minima:                      # metastable branches are passive as well
        chi = susceptibility(m, B, e.m, fr, alpha=0.02, gyro=GY, kind=kind)
        assert passivity_margin(chi) > -1e-9


def test_passivity_detects_negative_damping():
    m = MagnetModel(Ms=MS, N=sphere())
    with pytest.raises(ValueError):
        susceptibility(m, [0, 0, 0.4], [0, 0, 1], [1e9], alpha=-0.01, gyro=GY)
    # forcing a gain medium by conjugating the time convention violates the check
    chi = susceptibility(m, [0, 0, 0.4], [0, 0, 1], np.linspace(1e9, 30e9, 50), alpha=0.02, gyro=GY)
    assert passivity_margin(chi) > -1e-9
    assert passivity_margin(np.conj(chi)) < -1e-3


def test_non_equilibrium_rejected():
    m = MagnetModel(Ms=MS, K1=SYNTH_K1, N=sphere())
    with pytest.raises(ValueError):
        susceptibility(m, [0, 0, 0.0], [1, 0, 1], [1e9], alpha=0.02, gyro=GY)


def test_external_internal_relation_and_thin_film():
    N = thin_disk_perp_c()
    m = MagnetModel(Ms=MS, N=N)
    B = 0.9       # perpendicular saturation of a film: B > mu0 Ms
    fr = np.linspace(1e9, 40e9, 30)
    ce = susceptibility(m, [0, 0, B], [0, 0, 1], fr, alpha=0.02, gyro=GY, kind="external")
    ci = susceptibility(m, [0, 0, B], [0, 0, 1], fr, alpha=0.02, gyro=GY, kind="internal")
    for a, b in zip(ce, ci):
        assert np.allclose(a, b @ np.linalg.inv(np.eye(3) + N @ b), rtol=1e-9, atol=1e-12)
    # for the perpendicular film the dynamic demag in the plane is zero (N_xx = 0): ext == int
    assert np.allclose(ce, ci, rtol=1e-9, atol=1e-12)
    # in-plane magnetised film: dynamic demag matters, ext != int
    ce2 = susceptibility(m, [0.3, 0, 0], [1, 0, 0], fr, alpha=0.02, gyro=GY, kind="external")
    ci2 = susceptibility(m, [0.3, 0, 0], [1, 0, 0], fr, alpha=0.02, gyro=GY, kind="internal")
    assert not np.allclose(ce2, ci2, rtol=1e-3)


def test_fwhm_of_Im_chi_matches_2_alpha_f0():
    alpha = 0.005
    m = MagnetModel(Ms=MS, N=sphere())
    B = 0.4
    f0 = GAMMA_2PI * B
    fr = np.linspace(0.5 * f0, 1.5 * f0, 20001)
    chi = susceptibility(m, [0, 0, B], [0, 0, 1], fr, alpha=alpha, gyro=GY)
    r = fwhm(fr, chi[:, 0, 0].imag, mode_x=f0)
    assert r.status == "valid"
    assert abs(r.width - 2 * alpha * f0) < 0.02 * 2 * alpha * f0
    assert abs(r.peak_x - f0) < 3 * alpha * f0 and r.mode_x == f0


def test_kramers_kronig_numerical():
    alpha = 0.05
    m = MagnetModel(Ms=MS, K1=SYNTH_K1, K2=SYNTH_K2, N=np.diag([0.2, 0.3, 0.5]))
    B = field_vector(0.1, 0.6, 0.0)
    e = find_equilibria(m, B).ground
    f0 = 20e9
    fr = np.linspace(0.0, 300e9, 6001)            # very wide: far beyond the resonance
    chi = susceptibility(m, B, e.m, fr, alpha=alpha, gyro=GY)
    w = 2 * math.pi * fr
    for (i, j) in ((0, 0), (1, 1)):
        im = chi[:, i, j].imag
        kk = kramers_kronig_real(w, im)
        re = chi[:, i, j].real
        sel = (fr > 1e9) & (fr < 120e9)
        err = np.abs(kk[sel] - re[sel]).max() / np.abs(re[sel]).max()
        assert err < 2e-2, (i, j, err)
    # sum rule at w=0 : chi'(0) = (2/pi) int chi''/w dw
    im = chi[1:, 0, 0].imag
    s0 = (2 / math.pi) * np.trapezoid(im / w[1:], w[1:])
    assert abs(s0 - chi[0, 0, 0].real) < 2e-2 * abs(chi[0, 0, 0].real)
