"""Regressions from the independent review of the magnet core (round 1). All parameters are SYNTHETIC.

Sections: anisotropy fields / notation, equilibrium search, derivative kernels, H_res, neutral modes, LLG in time,
polycrystal (azimuth, demag modes), apparent permeability, line widths, uncertainty."""
import math
import time
import warnings

import numpy as np
import pytest
from scipy import optimize
from scipy.integrate import solve_ivp
from scipy.spatial import cKDTree

from hexw.theory.magnet import (Gyro, H_res, MagnetModel, MU0, Texture, apparent_mu_ellipsoid, apparent_mu_shape,
                                derivative_peak_to_peak, effective_medium_mixing, ellipsoid, energy, f_res,
                                field_vector, find_equilibria, fwhm, llg_omega, map_literature_constants, monte_carlo,
                                mu_r, natural_fmr, no_shape, polycrystal_chi, sample_demag_correction,
                                schloemann_porosity_linewidth, sphere, sphere_derivs, susceptibility,
                                sweep_branches, thin_disk_perp_c, width_freq_to_field)
from hexw.theory.magnet._fast import ScalarKernel
from hexw.theory.magnet.core import (_aniso_terms, eig2_batch, grad_hess, grad_hess_batch, sphere_derivs_batch,
                                     tangent_basis)
from hexw.theory.magnet.equilibrium import fibonacci_sphere, relax_batch
from hexw.theory.magnet.fmr import modes_at
from hexw.theory.magnet.linewidth import noise_sigma
from _synth import GAMMA_2PI, GY, SYNTH_K1, SYNTH_K2, SYNTH_K3, SYNTH_K4, SYNTH_MS

MS, G2 = SYNTH_MS, GAMMA_2PI
MU0MS = MU0 * MS


# =============================================================================== 1. H_A, notation

@pytest.mark.parametrize("K4", [-8e3, +8e3])
def test_basal_state_fields_and_minima_for_both_signs_of_K4(K4):
    K1, K2, K3 = -SYNTH_K1, SYNTH_K2, SYNTH_K3
    m = MagnetModel(Ms=MS, K1=K1, K2=K2, K3=K3, K4=K4, N=sphere())
    eqs = find_equilibria(m, np.zeros(3))
    assert len(eqs.minima) == 6 and all(e.kind == "planar" for e in eqs.minima)
    off = 0.0 if K4 < 0 else 30.0                                   # cos6phi0 = -sgn(K4)
    for e in eqs.minima:
        assert abs(((math.degrees(e.phi) - off + 30.0) % 60.0) - 30.0) < 1e-3
    Hth, Hph = m.mu0_H_theta_basal, m.mu0_H_phi_basal
    assert abs(Hth - (-2 * K1 - 4 * K2 - 6 * K3 + 6 * abs(K4)) / MS) < 1e-15
    assert abs(Hph - 36 * abs(K4) / MS) < 1e-15
    f_ref = G2 * math.sqrt(Hth * Hph)         # sphere: the shape adds the same constant to both
    for x in natural_fmr(m, gyro=GY).modes:
        assert abs(x.f_Hz - f_ref) < 1e-9 * f_ref                     # the natural frequency depends on |K4| only
    # ground-state energy at theta = pi/2: K1 + K2 + K3 + K4 cos6phi0 = K1 + K2 + K3 - |K4| (+ sphere shape constant)
    assert abs(eqs.ground.energy - (K1 + K2 + K3 - abs(K4) + MU0 * MS ** 2 / 6)) < 1e-6 * abs(K1)


def test_K4_positive_equals_K4_negative_in_natural_frequency():
    kw = dict(Ms=MS, K1=-SYNTH_K1, K2=SYNTH_K2, N=sphere())
    fm = natural_fmr(MagnetModel(K4=-SYNTH_K4 * 0.8, **kw), gyro=GY).modes[0].f_Hz
    fp = natural_fmr(MagnetModel(K4=+SYNTH_K4 * 0.8, **kw), gyro=GY).modes[0].f_Hz
    assert abs(fm - fp) < 1e-9 * fm


def test_c_state_anisotropy_field_is_2K1_over_Ms_independent_of_higher_constants():
    m = MagnetModel(Ms=MS, K1=SYNTH_K1, K2=SYNTH_K2, K3=SYNTH_K3, K4=SYNTH_K4, N=sphere())
    modes = natural_fmr(m, gyro=GY).modes
    assert len(modes) == 2 and all(x.state.kind == "axial" for x in modes)
    for x in modes:
        assert abs(x.f_Hz - G2 * m.mu0_HA) < 1e-9 * x.f_Hz and abs(m.mu0_HA - 2 * SYNTH_K1 / MS) < 1e-15


@pytest.mark.parametrize("K4,B", [(-8e3, 1.2), (+8e3, 1.2), (0.0, 0.6), (-8e3, 0.3)])
def test_planar_kittel_by_hand_field_in_plane(K4, B):
    """Independent hand Kittel formula (reviewer probe p1): film (disk perp c), field along phi = 0, m || B."""
    K1, K2, K3 = -1.2e5, 1e4, 2e3
    m = MagnetModel(Ms=MS, K1=K1, K2=K2, K3=K3, K4=K4, N=thin_disk_perp_c())
    Ft = B * MS + MU0 * MS ** 2 + (-2 * K1 - 4 * K2 - 6 * K3 - 6 * K4)       # cos(6 phi) = 1 at phi = 0
    Fp = B * MS + (-36 * K4)
    f_hand = G2 * math.sqrt(Ft * Fp) / MS
    ground = [x for x in f_res(m, B, math.pi / 2, 0.0, gyro=GY).modes if x.is_ground]
    assert ground and all(x.state.kind == "planar" for x in ground)
    assert abs(ground[0].f_Hz - f_hand) < 1e-9 * f_hand


def test_literature_notation_map_is_explicit():
    # a paper's "K3 sin^6 cos6phi" is our K4
    assert map_literature_constants({"K1": -1e5, "K2": 3e4, "K3": 8e3}, {"K1": "K1", "K2": "K2", "K3": "K4"}) == \
        {"K1": -1e5, "K2": 3e4, "K4": 8e3}
    with pytest.raises(KeyError):
        map_literature_constants({"K3": 1.0}, {})
    with pytest.raises(ValueError):
        map_literature_constants({"a": 1.0, "b": 2.0}, {"a": "K4", "b": "K4"})
    with pytest.raises(ValueError):
        map_literature_constants({"K5": 1.0}, {"K5": "K5"})


# =============================================================================== 2. equilibrium search

def test_twelve_cone_minima_at_n_grid_200():
    m = MagnetModel(Ms=MS, K1=-1e5, K2=3e5, K4=-1e3, N=sphere())
    for ng in (100, 200, 600):
        e = find_equilibria(m, np.zeros(3), n_grid=ng)
        assert len(e.minima) == 12 and e.status == "ok" and not e.possibly_incomplete
    th = sorted(set(round(math.degrees(x.theta), 2) for x in e.minima))
    assert th == [24.09, 155.91] or all(abs(t - a) < 0.05 for t, a in zip(th, [24.09, 155.91]))


def test_coarse_grid_flags_possibly_incomplete_and_returns_denser_result():
    m = MagnetModel(Ms=MS, K1=-1e5, K2=3e5, K4=-1e3, N=sphere())
    e = find_equilibria(m, np.zeros(3), n_grid=3)
    assert e.status == "possibly_incomplete" and e.possibly_incomplete
    assert any("possibly_incomplete" in n for n in e.notes)
    assert len(e.minima) >= 1


def _reference_minima(m, B, n=20000):
    """Independent reference: grid minima of a 20000-point grid polished by scipy BFGS in a tangent chart."""
    from hexw.theory.magnet.core import energy_many
    pts = fibonacci_sphere(n)
    sc = m.energy_scale(float(np.linalg.norm(B)))
    E = energy_many(m, B, pts)
    _, nb = cKDTree(pts).query(pts, k=9)
    ismin = np.all(E[:, None] <= E[nb[:, 1:]] + 1e-14 * sc, axis=1)
    out = []
    for i in np.nonzero(ismin)[0]:
        x0 = pts[i]
        Eb = tangent_basis(x0)
        res = optimize.minimize(lambda s: energy(m, B, x0 + Eb @ s) / sc, np.zeros(2), method="BFGS",
                                options={"gtol": 1e-10})
        x = (x0 + Eb @ res.x)
        x = x / np.linalg.norm(x)
        _, Hs, _ = sphere_derivs(m, B, x)
        if np.linalg.eigvalsh(Hs)[0] < -1e-7 * sc:
            continue
        if not any(np.linalg.norm(x - o) < 1e-2 for o in out):
            out.append(x)
    return out


def test_sixty_random_models_no_missed_minima_against_20000_point_reference():
    rng = np.random.default_rng(7)
    missed = 0
    n_ref = 0
    for _ in range(60):
        K1 = rng.uniform(-1.5, 1.5) * 1e5
        K2 = rng.uniform(-1, 3) * 1e5 * rng.choice([0.1, 1])
        K3 = rng.uniform(-1, 1) * 3e4
        K4 = rng.uniform(-1, 1) * 3e4
        N = np.diag(rng.dirichlet([1, 1, 1]))
        m = MagnetModel(Ms=MS, K1=K1, K2=K2, K3=K3, K4=K4, N=N)
        B = field_vector(rng.uniform(0, 0.3) * rng.choice([0, 1, 1]), rng.uniform(0, math.pi), rng.uniform(0, 6.28))
        ref = _reference_minima(m, B)
        got = [x.m for x in find_equilibria(m, B, n_grid=200).minima]
        n_ref += len(ref)
        missed += sum(1 for r in ref if not any(np.linalg.norm(r - g) < 2e-2 for g in got))
    assert n_ref > 100
    assert missed == 0


# =============================================================================== derivative kernels

def _ext_energy(m, B, r):
    """Cartesian extension of the energy at an arbitrary (non-unit) point r, as used by grad_hess."""
    return _aniso_terms(m, r)[0] - m.Ms * float(np.dot(B, r)) + 0.5 * MU0 * m.Ms ** 2 * float(r @ m.N @ r)


@pytest.mark.parametrize("seed", range(5))
def test_grad_hess_matches_finite_differences(seed):
    rng = np.random.default_rng(seed)
    m = MagnetModel(Ms=MS, K1=rng.uniform(-1, 1) * 1e5, K2=rng.uniform(-1, 1) * 5e4, K3=rng.uniform(-1, 1) * 3e4,
                    K4=rng.uniform(-1, 1) * 3e4, N=np.diag(rng.dirichlet([1, 1, 1])))
    B = rng.normal(size=3) * 0.2
    for _ in range(4):
        r = rng.normal(size=3)
        r /= np.linalg.norm(r)
        F, g, H = grad_hess(m, B, r)
        assert abs(F - _ext_energy(m, B, r)) < 1e-9 * m.energy_scale(0.2)
        h = 1e-5
        gfd = np.zeros(3)
        Hfd = np.zeros((3, 3))
        for i in range(3):
            ei = np.eye(3)[i] * h
            gfd[i] = (_ext_energy(m, B, r + ei) - _ext_energy(m, B, r - ei)) / (2 * h)
            for j in range(3):
                ej = np.eye(3)[j] * h
                Hfd[i, j] = (_ext_energy(m, B, r + ei + ej) - _ext_energy(m, B, r + ei - ej)
                             - _ext_energy(m, B, r - ei + ej) + _ext_energy(m, B, r - ei - ej)) / (4 * h * h)
        sc = m.energy_scale(0.2)
        assert np.abs(g - gfd).max() < 1e-6 * sc
        assert np.abs(H - Hfd).max() < 1e-4 * sc


@pytest.mark.parametrize("seed", range(4))
def test_sphere_derivs_match_finite_differences_of_retraction(seed):
    """Riemannian gradient/Hessian = derivatives of F(normalise(m + E s)) at s = 0 (a second-order retraction)."""
    rng = np.random.default_rng(10 + seed)
    m = MagnetModel(Ms=MS, K1=rng.uniform(-1, 1) * 1e5, K2=rng.uniform(-1, 1) * 5e4, K4=rng.uniform(-1, 1) * 3e4,
                    N=np.diag(rng.dirichlet([1, 1, 1])))
    B = rng.normal(size=3) * 0.2
    r = rng.normal(size=3)
    r /= np.linalg.norm(r)
    gt, Hs, Eb = sphere_derivs(m, B, r)
    f = lambda s: energy(m, B, r + Eb @ s)
    h = 1e-4
    gfd = np.array([(f(np.eye(2)[i] * h) - f(-np.eye(2)[i] * h)) / (2 * h) for i in range(2)])
    Hfd = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            ei, ej = np.eye(2)[i] * h, np.eye(2)[j] * h
            Hfd[i, j] = (f(ei + ej) - f(ei - ej) - f(-ei + ej) + f(-ei - ej)) / (4 * h * h)
    sc = m.energy_scale(0.2)
    assert np.abs(gt - gfd).max() < 1e-6 * sc and np.abs(Hs - Hfd).max() < 2e-4 * sc


def test_batch_and_scalar_kernels_agree():
    rng = np.random.default_rng(5)
    m = MagnetModel(Ms=MS, K1=-1.2e5, K2=2e4, K3=5e3, K4=-9e3, N=np.diag([0.2, 0.3, 0.5]))
    B = field_vector(0.35, 1.1, 0.7)
    M = rng.normal(size=(6, 3))
    M /= np.linalg.norm(M, axis=1)[:, None]
    gt, Hs, Eb = sphere_derivs_batch(m, B, M)
    lmin, lmax, vmin = eig2_batch(Hs)
    kern = ScalarKernel(m, B)
    for i in range(len(M)):
        g0, g1, a, b, d, e1, e2 = kern.derivs(*M[i])
        gts, Hss, _ = sphere_derivs(m, B, M[i])
        sc = m.energy_scale(0.35)
        # numpy bases of different kernels may differ only by a rotation; compare invariants
        assert abs(np.hypot(g0, g1) - np.linalg.norm(gts)) < 1e-10 * sc
        assert abs(np.linalg.norm(gt[i]) - np.linalg.norm(gts)) < 1e-10 * sc
        ls = np.linalg.eigvalsh(Hss)
        assert abs(kern.eig(a, b, d)[0] - ls[0]) < 1e-9 * sc and abs(kern.eig(a, b, d)[1] - ls[1]) < 1e-9 * sc
        assert abs(lmin[i] - ls[0]) < 1e-9 * sc and abs(lmax[i] - ls[1]) < 1e-9 * sc
        assert abs(kern.energy(*M[i]) - energy(m, B, M[i])) < 1e-9 * sc
        # eigenvector of lmin in the batch basis
        assert np.linalg.norm(Hs[i] @ vmin[i] - lmin[i] * vmin[i]) < 1e-9 * sc
    # relaxations agree
    Mb, ok = relax_batch(m, B, M)
    assert ok.all()
    for i in range(len(M)):
        ms_, oks, lm, _ = kern.relax(tuple(M[i]))
        assert oks and np.linalg.norm(np.array(ms_) - Mb[i]) < 1e-6 or abs(
            energy(m, B, ms_) - energy(m, B, Mb[i])) < 1e-9 * m.energy_scale(0.35)


def test_trust_region_newton_never_raises_energy_from_far_starts():
    m = MagnetModel(Ms=MS, K1=-1e5, K2=3e5, K4=-1e3, N=sphere())
    from hexw.theory.magnet.core import energy_many
    P = fibonacci_sphere(300)
    M, ok = relax_batch(m, np.zeros(3), P)
    assert (energy_many(m, np.zeros(3), M) <= energy_many(m, np.zeros(3), P) + 1e-9).all() and ok.mean() > 0.99


# =============================================================================== 3. H_res regressions

def _unique(xs, tol=1e-4):
    out = []
    for x in sorted(xs):
        if not out or abs(x - out[-1]) > tol:
            out.append(x)
    return out


@pytest.mark.parametrize("ft_GHz", [0.3, 1.0, 3.0])
def test_H_res_narrow_dip_uniaxial_sphere_field_perpendicular_to_c(ft_GHz):
    """Probe p3: both analytic roots must be found although the dip f(B) -> 0 at B = H_A lies between samples."""
    K1, Ms = 1e5, 3e5
    HA = 2 * K1 / Ms
    g = 28e9
    gy = Gyro.from_gamma_over_2pi(g)
    m = MagnetModel(Ms=Ms, K1=K1, N=sphere())
    ft = ft_GHz * 1e9
    r = H_res(m, ft, math.pi / 2, 0.0, B_range_T=(0.0, 1.5), gyro=gy, n_points=241)
    assert r.status == "ok"
    ana = [(HA + math.sqrt(HA ** 2 + 4 * (ft / g) ** 2)) / 2]
    if ft < g * HA:
        ana.append(math.sqrt(HA ** 2 - (ft / g) ** 2))
    got = _unique([x.B_res_T for x in r.roots], tol=1e-6)
    assert len(got) == len(ana)
    assert np.allclose(sorted(got), sorted(ana), atol=2e-6)
    for x in r.roots:                                # every root is a true resonance of its state
        assert x.f_Hz == ft


def test_H_res_metastable_branch_ending_at_spinodal():
    K1, Ms = 1e5, 3e5
    HA = 2 * K1 / Ms
    g = 28e9
    gy = Gyro.from_gamma_over_2pi(g)
    m = MagnetModel(Ms=Ms, K1=K1, N=sphere())
    r = H_res(m, 0.2e9, 0.0, 0.0, B_range_T=(-1.0, 0.0), gyro=gy, n_points=121)
    assert r.status == "ok" and len(r.roots) == 1
    assert abs(r.roots[0].B_res_T - (-HA + 0.2e9 / g)) < 1e-6 and not r.roots[0].is_ground
    # sweeping upward from -1 T the +z state is BORN at the spinodal field -H_A (found by backward tracking), f(-H_A) = 0
    br = [b for b in r.branches if b.B_birth is not None]
    assert len(br) == 1 and abs(br[0].B_birth + HA) < 1e-5
    assert abs(r.roots[0].df_dB_Hz_per_T - g) < 1e-3 * g          # linear branch f = g (B + HA)


def _random_hres_config(rng):
    K1 = rng.uniform(-1.5, 1.5) * 1e5
    K2 = rng.uniform(-0.5, 1) * 1e5 * rng.choice([0.1, 1])
    K3 = rng.uniform(-1, 1) * 2e4
    K4 = rng.uniform(-1, 1) * 2e4 * rng.choice([0, 1])
    N = np.diag(rng.dirichlet([2, 2, 2])) if rng.random() < 0.6 else sphere()
    m = MagnetModel(Ms=MS, K1=K1, K2=K2, K3=K3, K4=K4, N=N)
    th = rng.choice([0.0, math.pi / 2, rng.uniform(0, math.pi)])
    ph = rng.uniform(0, 6.28)
    Bhi = rng.uniform(0.5, 1.5)
    Blo = rng.choice([0.0, -Bhi * rng.uniform(0.3, 1)])
    ft = rng.uniform(3, 40) * 1e9
    return m, th, ph, Blo, Bhi, ft


def _dense_root_count(m, th, ph, Blo, Bhi, ft, gy, n=500):
    d = field_vector(1.0, th, ph)
    cnt = [sum(1 for x in modes_at(m, b * d, gy, n_grid=100).modes if x.f_Hz > ft) for b in np.linspace(Blo, Bhi, n)]
    return int(np.abs(np.diff(cnt)).sum())


def test_H_res_flags_possibly_incomplete_on_a_coarse_field_grid():
    """13 field points: the curvature of F = f^2 between nodes does not allow to exclude hidden roots -> the status says so
    (the roots found are still returned, and agree with a dense scan of the mode set)."""
    gy = Gyro.from_gamma_over_2pi(28e9)
    m, th, ph, Blo, Bhi, ft = _random_hres_config(np.random.default_rng(11))
    r = H_res(m, ft, th, ph, B_range_T=(Blo, Bhi), gyro=gy, n_points=13)
    assert r.status == "possibly_incomplete" and r.possibly_incomplete and len(r.roots) >= 1
    assert any("possibly_incomplete" in n for n in r.notes)
    assert len(r.roots) == _dense_root_count(m, th, ph, Blo, Bhi, ft, gy, n=300)


def test_H_res_no_resonance_when_branch_is_smooth_and_far():
    m = MagnetModel(Ms=MS, N=sphere())
    r = H_res(m, 80e9, 0.3, 0.0, B_range_T=(0.0, 1.0), gyro=GY, n_points=41)
    assert r.status == "no_resonance_in_range" and not r.possibly_incomplete


@pytest.mark.slow
def test_H_res_root_count_matches_dense_scan_of_mode_set():
    """Number of roots = number of changes of #(modes above f_target) along a dense field scan (branches born/dying at
    folds have f -> 0 and do not change the count)."""
    rng = np.random.default_rng(3)
    gy = Gyro.from_gamma_over_2pi(28e9)
    for _ in range(4):
        m, th, ph, Blo, Bhi, ft = _random_hres_config(rng)
        r = H_res(m, ft, th, ph, B_range_T=(Blo, Bhi), gyro=gy, n_points=161)
        assert len(r.roots) == _dense_root_count(m, th, ph, Blo, Bhi, ft, gy), (th, Blo, Bhi, ft)


def test_H_res_speed_regression():
    """Planar film, 241 points: was ~3-5 s (global search at every step); now a few tenths of a second."""
    gy = Gyro.from_gamma_over_2pi(28e9)
    m = MagnetModel(Ms=MS, K1=-1.2e5, K2=1e4, K4=-8e3, N=sphere())
    t0 = time.perf_counter()
    r = H_res(m, 30e9, math.pi / 2, 0.0, B_range_T=(0.0, 1.5), gyro=gy)
    dt = time.perf_counter() - t0
    assert r.status == "ok" and len(r.roots) == 1 and r.roots[0].is_ground
    assert dt < 2.5, dt


def test_pitchfork_branch_merges_and_is_not_a_spurious_root():
    """Uniaxial sphere, B perp c: the two tilted states fuse into the saturated state at B = H_A."""
    gy = Gyro.from_gamma_over_2pi(28e9)
    m = MagnetModel(Ms=3e5, K1=1e5, N=sphere())
    r = H_res(m, 1e9, math.pi / 2, 0.0, B_range_T=(0.0, 1.5), gyro=gy, n_points=241)
    merged = [b for b in r.branches if b.end_reason == "merged"]
    assert len(merged) == 1 and abs(merged[0].B_death - 2e5 / 3e5) < 1e-5
    assert not any("discontinuity" in n for n in r.notes)


# =============================================================================== 4. neutral modes / NaN

def test_cone_with_field_along_c_has_no_nan_and_zero_frequencies():
    gy = Gyro.from_gamma_over_2pi(28e9)
    m = MagnetModel(Ms=3e5, K1=-1e5, K2=3e5, N=sphere())
    brs = sweep_branches(m, [0, 0, 1.0], np.linspace(0.0, 0.3, 13), gy)
    assert brs
    for b in brs:
        f = np.array(b.f_Hz)
        assert np.isfinite(f).all()
        assert all(s in ("ok", "neutral_zero_mode") for s in b.mode_status)
    cone0 = [b for b in brs if b.kind[0] == "cone"]
    assert cone0 and all(s == "neutral_zero_mode" for s in cone0[0].mode_status) and (np.array(cone0[0].f_Hz) == 0).all()
    # the same state through llg_omega and modes_at: identical zero mode, no NaN
    e = find_equilibria(m, np.zeros(3)).ground
    assert llg_omega(m, np.zeros(3), e.m, gy) == 0.0
    assert all(x.f_Hz == 0.0 and x.status == "neutral_zero_mode" for x in modes_at(m, np.zeros(3), gy).modes)


def test_llg_omega_nan_only_for_non_minimum():
    m = MagnetModel(Ms=MS, K1=SYNTH_K1, N=sphere())
    assert math.isnan(llg_omega(m, np.zeros(3), np.array([1.0, 0, 0]), GY))        # hard direction: a saddle
    assert llg_omega(m, np.zeros(3), np.array([0, 0, 1.0]), GY) > 0


# =============================================================================== LLG in time (reviewer probe p2)

@pytest.mark.slow
def test_chi_against_time_integration_of_gilbert_llg():
    K1, K2, K3, K4 = -1.2e5, 1e4, 2e3, -8e3
    gam, alpha = 2 * math.pi * 28e9, 0.05
    N = np.diag([0.2, 0.3, 0.5])
    model = MagnetModel(Ms=MS, K1=K1, K2=K2, K3=K3, K4=K4, N=N)
    gy = Gyro(gam)
    Bs = field_vector(0.4, 1.2, 0.5)

    def Fm(mv):                              # independent energy in (theta, phi)
        th = math.acos(mv[2])
        ph = math.atan2(mv[1], mv[0])
        s = math.sin(th)
        return (K1 * s ** 2 + K2 * s ** 4 + K3 * s ** 6 + K4 * s ** 6 * math.cos(6 * ph) - MS * Bs @ mv
                + 0.5 * MU0 * MS ** 2 * mv @ N @ mv)

    def gradF(mv, h=1e-6):
        g = np.zeros(3)
        for i in range(3):
            d = np.zeros(3)
            d[i] = h
            mp, mm = mv + d, mv - d
            g[i] = (Fm(mp / np.linalg.norm(mp)) - Fm(mm / np.linalg.norm(mm))) / (2 * h)
        return g

    m0 = find_equilibria(model, Bs).ground.m
    f0 = llg_omega(model, Bs, m0, gy) / (2 * math.pi)
    H0 = 100.0
    hdir = np.array([0.3, 0.8, 0.5])
    hdir /= np.linalg.norm(hdir)

    def run(f):
        w = 2 * math.pi * f

        def rhs(t, mv):
            mm = mv / np.linalg.norm(mv)
            Be = -gradF(mm) / MS + MU0 * H0 * hdir * math.cos(w * t)
            a = -gam * np.cross(mm, Be)
            return (a + alpha * np.cross(mm, a)) / (1 + alpha ** 2)
        T = 2 * math.pi / w
        tf = T * (int(2.0 / alpha) + 10)
        ts = np.linspace(tf - 4 * T, tf, 200, endpoint=False)
        sol = solve_ivp(rhs, [0, tf], m0, t_eval=ts, rtol=1e-10, atol=1e-12, method="DOP853")
        M1 = MS * (sol.y.T - m0)
        A = np.stack([np.cos(w * ts), np.sin(w * ts)], 1)
        return np.array([complex(*np.linalg.lstsq(A, M1[:, i], rcond=None)[0]) for i in range(3)]) / H0

    for f in (0.98 * f0, 1.5 * f0):
        chi = susceptibility(model, Bs, m0, [f], alpha=alpha, gyro=gy)[0]
        ref = run(f)
        assert np.abs(chi @ hdir - ref).max() / np.abs(ref).max() < 2e-3


# =============================================================================== 5. polycrystal

def _trace(res):
    return np.trace(res.chi_lab, axis1=1, axis2=2) / 3


@pytest.mark.slow
def test_azimuth_average_p5_converges_monotonically_and_matches_reference():
    """Probe p5: K4 != 0, sphere, B = 0.3 T along z, alpha = 0.02 -- the azimuth of the crystal must be averaged."""
    gy = Gyro.from_gamma_over_2pi(28e9)
    fs = np.linspace(6e9, 60e9, 55)
    m = MagnetModel(Ms=MS, K1=-1.2e5, K2=1e4, K4=-8e3, N=sphere())
    B = np.array([0.0, 0.0, 0.3])
    kw = dict(alpha=0.02, gyro=gy, state="field_aligned")
    ref = polycrystal_chi(m, B, fs, n_c=100, n_beta=24, **kw)
    assert ref.converged and ref.quadrature.startswith("gauss_legendre") and ref.n_distinct_azimuths_mod_pi3 == 24
    rt = _trace(ref)
    errs = []
    for nc in (12, 25, 50):
        errs.append(np.linalg.norm(_trace(polycrystal_chi(m, B, fs, n_c=nc, **kw)) - rt) / np.linalg.norm(rt))
    assert errs[0] > errs[1] > errs[2] and errs[1] < 0.03
    # independent estimate: random orientations with quasi-random azimuths (>= 24 distinct) agree within 3%
    q = polycrystal_chi(m, B, fs, n_c=1600, quadrature="quasi", **kw)
    assert q.n_distinct_azimuths_mod_pi3 >= 24 and q.quadrature == "quasi_random"
    assert np.linalg.norm(_trace(q) - rt) / np.linalg.norm(rt) < 0.03


def test_zero_azimuth_sampling_is_wrong_and_warns():
    """The old sampling (every grain at azimuth 0) misses the K4 average by tens of percent; it must warn."""
    gy = Gyro.from_gamma_over_2pi(28e9)
    fs = np.linspace(6e9, 60e9, 28)
    m = MagnetModel(Ms=MS, K1=-1.2e5, K2=1e4, K4=-8e3, N=sphere())
    B = np.array([0.0, 0.0, 0.3])
    kw = dict(alpha=0.02, gyro=gy, state="field_aligned")
    ref = _trace(polycrystal_chi(m, B, fs, n_c=40, **kw))
    with pytest.warns(RuntimeWarning, match="distinct crystal azimuths"):
        bad = polycrystal_chi(m, B, fs, n_c=200, quadrature="quasi", azimuth="zero", **kw)
    assert bad.n_distinct_azimuths_mod_pi3 == 1
    assert np.linalg.norm(_trace(bad) - ref) / np.linalg.norm(ref) > 0.1


def test_gauss_legendre_equals_random_orientation_average_for_K4_zero():
    gy = GY
    fs = np.linspace(8e9, 50e9, 20)
    m = MagnetModel(Ms=MS, K1=SYNTH_K1, K2=SYNTH_K2, N=sphere())
    B = np.array([0.0, 0.0, 0.5])
    kw = dict(alpha=0.2, gyro=gy, state="field_aligned")
    gl = polycrystal_chi(m, B, fs, n_c=40, **kw)
    assert gl.quadrature == "gauss_legendre_cos_psi"
    q = polycrystal_chi(m, B, fs, n_c=2000, quadrature="quasi", **kw)
    assert np.abs(gl.chi_lab - q.chi_lab).max() / np.abs(gl.chi_lab).max() < 5e-3
    # the lab-axisymmetric average is a tensor diag(a, a, b) without gyrotropic cross terms along the axis
    c = gl.chi_lab
    assert np.abs(c[:, 0, 0] - c[:, 1, 1]).max() < 1e-12 * np.abs(c).max() and np.abs(c[:, 2, :2]).max() < 1e-12


def test_non_axisymmetric_lab_falls_back_to_quasi_random():
    gy = GY
    m = MagnetModel(Ms=MS, K1=SYNTH_K1, N=sphere())
    r = polycrystal_chi(m, [0.1, 0.0, 0.5], np.linspace(10e9, 40e9, 4), alpha=0.2, gyro=gy, n_c=30,
                        N_sample=np.diag([0.2, 0.3, 0.5]), state="field_aligned")
    assert r.quadrature == "quasi_random" and r.converged


def test_macro_meanfield_sphere_pole_shifts_between_internal_and_external():
    gy = Gyro.from_gamma_over_2pi(28e9)
    m = MagnetModel(Ms=MS, N=sphere())                    # isotropic grains
    Bapp = 0.5
    fs = np.linspace(8e9, 20e9, 1201)
    r = polycrystal_chi(m, [0, 0, Bapp], fs, alpha=0.01, gyro=gy, n_c=4, state="field_aligned")
    assert r.demag == "macro_meanfield" and not r.physically_invalid_for_bulk and r.converged
    assert r.kind == "internal" and r.chi_lab is r.chi_int_lab
    Bint = Bapp - MU0MS / 3.0
    assert np.allclose(r.B_int_lab, [0, 0, Bint], atol=1e-6)
    f_int = fs[np.argmax(r.chi_int_lab[:, 0, 0].imag)]
    f_ext = fs[np.argmax(r.chi_ext_lab[:, 0, 0].imag)]
    assert abs(f_int - 28e9 * Bint) < 0.02e9                # response to the internal field: pole at gamma B_int
    assert abs(f_ext - 28e9 * Bapp) < 0.02e9                # sphere Kittel: gamma B_app after the dynamic demag field
    assert r.chi_ext_lab.shape == r.chi_int_lab.shape
    ext = polycrystal_chi(m, [0, 0, Bapp], fs[::100], alpha=0.01, gyro=gy, n_c=4, state="field_aligned", kind="external")
    assert ext.kind == "external" and np.allclose(ext.chi_lab, sample_demag_correction(ext.chi_int_lab, sphere()))


def test_macro_meanfield_self_consistency_thin_disk():
    gy = GY
    m = MagnetModel(Ms=MS, K1=SYNTH_K1, N=sphere())
    N_s = thin_disk_perp_c()
    Bapp = np.array([0.0, 0.0, 0.9])
    r = polycrystal_chi(m, Bapp, np.linspace(10e9, 40e9, 4), alpha=0.1, gyro=gy, n_c=16, N_sample=N_s,
                        state="field_aligned", macro_tol_T=1e-10)
    assert r.converged
    Mz = r.M_avg_A_per_m[2]
    assert abs(r.B_int_lab[2] - (Bapp[2] - MU0 * Mz)) < 1e-8               # B_int = B_app - mu0 N M
    assert 0 < Mz <= MS and np.allclose(r.N_sample, N_s)


def test_per_grain_mode_is_flagged_and_warns():
    gy = GY
    m = MagnetModel(Ms=MS, K1=SYNTH_K1, N=sphere())
    with pytest.warns(UserWarning, match="physically invalid for a bulk sample"):
        r = polycrystal_chi(m, [0, 0, 1.2], np.linspace(10e9, 60e9, 4), alpha=0.2, gyro=gy, n_c=10,
                            demag="per_grain", state="field_aligned")
    assert r.physically_invalid_for_bulk and r.demag == "per_grain" and r.kind == "external"
    assert any("per_grain" in n for n in r.notes)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ok = polycrystal_chi(m, [0, 0, 1.2], np.linspace(10e9, 60e9, 4), alpha=0.2, gyro=gy, n_c=10,
                             state="field_aligned")
    assert not ok.physically_invalid_for_bulk and ok.demag == "macro_meanfield"
    with pytest.raises(ValueError):
        polycrystal_chi(m, [0, 0, 1.2], [10e9], alpha=0.2, gyro=gy, n_c=4, demag="nope")


def test_grain_without_own_shape_is_allowed_and_n_must_otherwise_have_trace_one():
    m = MagnetModel(Ms=MS, K1=SYNTH_K1, N=no_shape())
    assert not np.any(m.N)
    with pytest.raises(ValueError):
        MagnetModel(Ms=MS, K1=SYNTH_K1, N=np.diag([0.2, 0.2, 0.2]))
    e = find_equilibria(m, field_vector(0.3, 0.4))
    assert e.ground is not None


def test_effective_medium_and_schloemann_are_explicit_stubs():
    for fn in (effective_medium_mixing, schloemann_porosity_linewidth):
        r = fn(1, 2, x=3)
        assert r.status == "not_implemented" and r.value is NotImplemented and r.reason


# =============================================================================== 6. apparent permeability

def test_mu_r_is_material_and_refuses_external():
    chi = np.zeros((2, 3, 3), complex) + 0.5 * np.eye(3)
    assert np.allclose(mu_r(chi), 1.5 * np.eye(3))
    with pytest.raises(ValueError, match="apparent_mu"):
        mu_r(chi, kind="external")


def test_apparent_mu_shape_requires_a_label():
    chi = np.ones((1, 3, 3), complex)
    a = apparent_mu_shape(chi, "sphere")
    assert a.shape_label == "sphere" and np.allclose(a.mu, np.eye(3) + chi)
    for bad in ("", "  ", None):
        with pytest.raises(ValueError):
            apparent_mu_shape(chi, bad)


def test_apparent_mu_ellipsoid_is_I_plus_I_minus_N_chi_ext():
    N = ellipsoid(1.0, 2.0, 3.0)
    m = MagnetModel(Ms=MS, K1=SYNTH_K1, N=N)
    B = field_vector(0.9, 0.0)
    e = find_equilibria(m, B).ground
    freqs = np.array([5e9, 20e9, 45e9])
    chi_ext = susceptibility(m, B, e.m, freqs, alpha=0.1, gyro=GY, kind="external")
    chi_int = susceptibility(m, B, e.m, freqs, alpha=0.1, gyro=GY, kind="internal")
    mu = apparent_mu_ellipsoid(chi_ext, N, "ellipsoid 1:2:3").mu
    I = np.eye(3)
    for k in range(len(freqs)):
        assert np.allclose(mu[k], I + (I - N) @ chi_ext[k])
        # <B>/(mu0 H_app) = (I + chi_int)(I - N chi_ext): flux density of the internal field H_int = (I - N chi_ext) H_app
        assert np.allclose(mu[k], (I + chi_int[k]) @ (I - N @ chi_ext[k]), atol=1e-9 * np.abs(chi_int).max())
    sph = apparent_mu_ellipsoid(chi_ext, sphere()).mu
    assert np.allclose(sph, I + (2.0 / 3.0) * chi_ext)             # sphere: 1 + (2/3) chi, not 1 + chi


# =============================================================================== 7. linewidth

def _lor(x, x0, w, a=1.0):
    return a / (1.0 + ((x - x0) / (w / 2.0)) ** 2)


def test_noisy_lorentzian_is_valid_and_within_five_percent():
    x = np.linspace(0, 10, 2001)
    worst = 0.0
    for seed in range(25):
        y = _lor(x, 5, 0.6) + 0.03 * np.random.default_rng(seed).normal(size=x.size)
        r = fwhm(x, y)
        assert r.status == "valid", (seed, r.status)
        assert r.smoothed_window > 0 and abs(r.noise_sigma - 0.03) < 0.005
        worst = max(worst, abs(r.width - 0.6) / 0.6)
    assert worst < 0.05
    # without smoothing the maximum is biased upward: the unsmoothed estimate is clearly worse in the mean
    ws = [fwhm(x, _lor(x, 5, 0.6) + 0.03 * np.random.default_rng(s).normal(size=x.size), smooth=None)
          for s in range(25)]
    assert np.mean([abs(w.width - 0.6) / 0.6 for w in ws if w.width]) > 0.02 or any(w.status != "valid" for w in ws)


def test_noise_sigma_estimate_ignores_the_line():
    x = np.linspace(-5, 5, 4001)
    y = _lor(x, 0, 1.0) + 0.02 * np.random.default_rng(0).normal(size=x.size)
    assert abs(noise_sigma(y) - 0.02) < 0.002
    assert noise_sigma(_lor(x, 0, 1.0)) < 1e-6


def test_derivative_pp_on_noise_smoothed_or_declared_unreliable():
    x = np.linspace(-20, 20, 8001)
    exact = 2.0 / math.sqrt(3.0)
    # light noise: valid and close; heavy noise: unreliable, no value
    light = [derivative_peak_to_peak(x, _lor(x, 0, 2.0) + 0.003 * np.random.default_rng(s).normal(size=x.size))
             for s in range(8)]
    ok = [d for d in light if d.status == "valid"]
    assert len(ok) >= 5 and all(abs(d.delta_pp - exact) < 0.2 * exact and d.smoothed_window > 0 for d in ok)
    heavy = [derivative_peak_to_peak(x, _lor(x, 0, 2.0) + 0.3 * np.random.default_rng(s).normal(size=x.size))
             for s in range(8)]
    assert all(d.status != "valid" and d.delta_pp is None for d in heavy)
    assert sum(d.status == "noisy_unreliable" for d in heavy) >= 6
    assert all(len(d.estimates) == 2 for d in heavy if d.status == "noisy_unreliable")
    # a clean line keeps the exact old algorithm
    d = derivative_peak_to_peak(x, _lor(x, 0, 2.0))
    assert d.status == "valid" and d.smoothed_window == 0 and abs(d.delta_pp - exact) < 1e-3


def test_slope_too_small_even_with_zero_tolerance():
    c = width_freq_to_field(1e9, 1e-3, slope_tol=0.0)
    assert c.status == "slope_too_small" and c.value is None
    assert width_freq_to_field(1e9, 0.0).status == "slope_zero"
    assert width_freq_to_field(1e9, float("nan")).status == "slope_zero"
    assert width_freq_to_field(1e9, 28e9).status in ("ok", "ok_with_warnings")
    assert width_freq_to_field(1e9, 5e5, slope_floor=1e3).value == pytest.approx(1e9 / 5e5)


def test_unresolved_asymmetric_doublet_is_flagged():
    x = np.linspace(0, 10, 4001)
    y = _lor(x, 5.0, 1.0) + 0.45 * _lor(x, 5.9, 1.0)           # strong overlap, no second maximum
    r = fwhm(x, y)
    assert r.status == "valid" and r.asymmetric_possible_doublet
    assert any("asymmetric_possible_doublet" in n for n in r.notes)
    single = fwhm(x, _lor(x, 5.0, 1.0))
    assert single.status == "valid" and not single.asymmetric_possible_doublet


# =============================================================================== 8. uncertainty

def test_monte_carlo_accepts_numpy_integers_and_floats():
    res = monte_carlo(lambda p: {"i": np.int64(3), "j": np.int32(4), "f": np.float32(1.5), "g": np.float64(2.0)},
                      {"x": 1.0}, n=5, seed=1)
    for k in ("i", "j", "f", "g"):
        assert res.summary[k].n_valid == 5, (k, res.failures)
    assert res.summary["i"].median == 3.0 and not res.failures


# =============================================================================== import hygiene (strict)

def test_import_does_not_load_process_pools():
    import subprocess
    import sys
    code = ("import sys\nimport hexw.theory.magnet\n"
            "assert 'multiprocessing.pool' not in sys.modules, 'multiprocessing.pool'\n"
            "assert 'concurrent.futures.process' not in sys.modules, 'concurrent.futures.process'\n")
    r = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr
