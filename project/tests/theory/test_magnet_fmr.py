"""FMR: analytic limits, Smit-Beljers vs LLG, rotation invariance, gamma conventions, H_res. SYNTHETIC inputs."""
import math

import numpy as np
import pytest

from hexw import units
from hexw.theory.magnet import (Gyro, H_res, MagnetModel, MU0, ellipsoid, f_res, field_vector, find_equilibria,
                                llg_omega, natural_fmr, resolve_gamma, smit_beljers_omega, sphere,
                                thin_disk_perp_c, thin_disk_parallel_c)
from hexw.theory.magnet.fmr import modes_at
from hexw.theory.magnet.polycrystal import to_crystal_frame
from _synth import GAMMA_2PI, GY, SYNTH_K1, SYNTH_K2, SYNTH_K3, SYNTH_K4, SYNTH_MS, rand_rotation

MS, G2 = SYNTH_MS, GAMMA_2PI
MU0MS = MU0 * MS


def f_of(model, B, th=0.0, ph=0.0, which="ground"):
    ms = f_res(model, B, th, ph, gyro=GY).modes
    return [m.f_Hz for m in ms if m.is_ground][0]


def test_kittel_sphere_no_anisotropy():
    for B in (0.05, 0.3, 1.0):
        f = f_of(MagnetModel(Ms=MS, N=sphere()), B, 0.4, 0.2)
        assert abs(f - G2 * B) < 1e-9 * G2 * B


def test_axial_natural_frequency_sphere():
    m = MagnetModel(Ms=MS, K1=SYNTH_K1, N=sphere())
    modes = natural_fmr(m, gyro=GY).modes
    HA_T = 2 * SYNTH_K1 / MS
    assert all(abs(x.f_Hz - G2 * HA_T) < 1e-9 * G2 * HA_T for x in modes) and len(modes) == 2


def test_axial_natural_with_K2_K3_and_disk_shape():
    K1, K2, K3 = SYNTH_K1, SYNTH_K2, SYNTH_K3
    m = MagnetModel(Ms=MS, K1=K1, K2=K2, K3=K3, N=thin_disk_perp_c())
    f = natural_fmr(m, gyro=GY).modes[0].f_Hz
    # c-axis state: H_eff = 2(K1+...)/Ms (K2,K3 do not contribute at theta=0) minus mu0 Ms for a disk normal to c
    assert abs(f - G2 * (2 * K1 / MS - MU0MS)) < 1e-9 * G2 * MU0MS


def test_film_in_plane_and_out_of_plane_kittel():
    B = 0.6
    inplane = MagnetModel(Ms=MS, N=thin_disk_perp_c())
    f = f_of(inplane, B, math.pi / 2, 0.0)
    assert abs(f - G2 * math.sqrt(B * (B + MU0MS))) < 1e-9 * f
    f = f_of(inplane, B + MU0MS, 0.0, 0.0)              # perpendicular: f = gamma (B - mu0 Ms), saturated
    assert abs(f - G2 * B) < 1e-9 * f


def test_general_ellipsoid_kittel_field_along_x():
    N = ellipsoid(1.0, 2.0, 3.0)
    m = MagnetModel(Ms=MS, N=N)
    B = 0.5
    f = f_of(m, B, math.pi / 2, 0.0)
    Nx, Ny, Nz = np.diag(N)
    ref = G2 * math.sqrt((B + (Ny - Nx) * MU0MS) * (B + (Nz - Nx) * MU0MS))
    assert abs(f - ref) < 1e-9 * ref


def test_planar_state_with_basal_anisotropy_natural():
    K1, K4 = -SYNTH_K1, SYNTH_K4              # K1<0 easy plane, K4<0 -> m at phi=0
    m = MagnetModel(Ms=MS, K1=K1, K2=SYNTH_K2, K3=SYNTH_K3, K4=K4, N=sphere())
    # independent closed form: H_theta, H_phi from theta-/phi-second derivatives at theta=pi/2, phi=0 (K4<0)
    H_th = (-2 * K1 - 4 * SYNTH_K2 - 6 * SYNTH_K3 - 6 * K4) / MS
    H_ph = (-36 * K4) / MS
    ref = G2 * math.sqrt(H_th * H_ph)
    modes = natural_fmr(m, gyro=GY).modes
    assert len(modes) == 6 and all(x.state.kind == "planar" for x in modes)
    for x in modes:
        assert abs(x.f_Hz - ref) < 1e-9 * ref


def test_sphere_shape_gives_no_shape_contribution_disk_does():
    kw = dict(Ms=MS, K1=SYNTH_K1)
    fs = natural_fmr(MagnetModel(N=sphere(), **kw), gyro=GY).modes[0].f_Hz
    fd = natural_fmr(MagnetModel(N=thin_disk_perp_c(), **kw), gyro=GY).modes[0].f_Hz
    assert abs(fs - G2 * 2 * SYNTH_K1 / MS) < 1e-6 * fs and fd < fs
    assert abs((fs - fd) - G2 * MU0MS) < 1e-6 * fs


def test_cone_goldstone_mode_flagged():
    m = MagnetModel(Ms=MS, K1=-1e5, K2=3e5, N=sphere())
    modes = natural_fmr(m, gyro=GY).modes
    assert all(x.status == "neutral_zero_mode" and x.f_Hz == 0.0 for x in modes)


@pytest.mark.parametrize("seed", range(6))
def test_smit_beljers_equals_llg(seed):
    rng = np.random.default_rng(seed)
    N = np.diag(rng.dirichlet([2, 2, 2]))
    m = MagnetModel(Ms=MS, K1=rng.uniform(-1, 1) * 1e5, K2=rng.uniform(0, 3e4), K3=rng.uniform(-1, 1) * 5e3,
                    K4=rng.uniform(-1, 1) * 1e4, N=N)
    B = field_vector(rng.uniform(0.05, 0.8), rng.uniform(0.1, 3.0), rng.uniform(0, 6.28))
    eqs = find_equilibria(m, B)
    assert eqs.minima
    for e in eqs.minima:
        wl = llg_omega(m, B, e.m, GY)
        ws, info = smit_beljers_omega(m, B, e.m, GY, pole="auto")
        assert abs(ws - wl) < 1e-9 * wl, info
        ws2, _ = smit_beljers_omega(m, B, e.m, GY, pole="local")
        assert abs(ws2 - wl) < 1e-9 * wl


def test_smit_beljers_singularity_and_continuity_near_pole():
    m = MagnetModel(Ms=MS, K1=SYNTH_K1, K4=1e3, N=np.diag([0.3, 0.3, 0.4]))
    B0 = field_vector(0.4, 0.0)
    e0 = find_equilibria(m, B0).ground
    w_ref = llg_omega(m, B0, e0.m, GY)
    ws, info = smit_beljers_omega(m, B0, np.array([0, 0, 1.0]), GY, pole="crystal")
    assert math.isnan(ws) and info["status"] == "singular_pole"      # explicit, not silently regularised
    ws, info = smit_beljers_omega(m, B0, np.array([0, 0, 1.0]), GY, pole="auto")
    assert info["pole_used"] == "local" and abs(ws - w_ref) < 1e-9 * w_ref
    # continuity as the field tilts away from c
    prev = None
    for tilt in (1e-9, 1e-7, 1e-5, 1e-3, 1e-2, 3e-2):
        B = field_vector(0.4, tilt, 0.3)
        e = find_equilibria(m, B).ground
        wl = llg_omega(m, B, e.m, GY)
        ws, info = smit_beljers_omega(m, B, e.m, GY, pole="auto")
        assert abs(ws - wl) < 1e-8 * wl
        assert abs(wl - w_ref) / w_ref < 1e-2 * (1 + tilt)          # no jump as theta -> 0
        prev = wl


@pytest.mark.parametrize("seed", range(4))
def test_frequency_invariant_under_lab_rotation(seed):
    rng = np.random.default_rng(100 + seed)
    m = MagnetModel(Ms=MS, K1=SYNTH_K1, K2=SYNTH_K2, K4=SYNTH_K4)
    N_lab = np.diag(rng.dirichlet([2, 2, 2]))
    B_lab = rng.normal(size=3) * 0.2
    R = rand_rotation(rng)                     # crystal -> lab
    S = rand_rotation(rng)                     # change of the laboratory frame
    f = []
    for (Rr, Br, Nr) in ((R, B_lab, N_lab), (S @ R, S @ B_lab, S @ N_lab @ S.T)):
        Bc, Nc = to_crystal_frame(Rr, Br, Nr)
        mm = m.with_N(Nc)
        f.append(sorted(x.f_Hz for x in modes_at(mm, Bc, GY).modes))
    assert len(f[0]) == len(f[1])
    assert np.allclose(f[0], f[1], rtol=1e-8)


def test_gamma_vs_gamma_over_2pi():
    g_rad = resolve_gamma(g=2.0)
    assert abs(g_rad - 2.0 * units.MU_B / units.HBAR) < 1e-6 * g_rad
    assert abs(resolve_gamma(gamma_over_2pi_Hz_per_T=g_rad / (2 * math.pi)) - g_rad) < 1e-6 * g_rad
    assert abs(g_rad / (2 * math.pi) / 1e9 - 27.99) < 0.02            # CODATA electron value (sanity of units only)
    m = MagnetModel(Ms=MS, N=sphere())
    f_a = f_res(m, 0.5, 0, 0, gamma_over_2pi_Hz_per_T=28e9).modes[0].f_Hz
    f_b = f_res(m, 0.5, 0, 0, gamma_rad_per_s_T=28e9 * 2 * math.pi).modes[0].f_Hz
    assert abs(f_a - 14e9) < 1 and abs(f_b - 14e9) < 1
    f_wrong = f_res(m, 0.5, 0, 0, gamma_rad_per_s_T=28e9).modes[0].f_Hz     # treating gamma/2pi as gamma
    assert abs(f_wrong * 2 * math.pi - f_a) < 1e-3 * f_a


def test_gamma_is_mandatory_and_unique():
    m = MagnetModel(Ms=MS, N=sphere())
    with pytest.raises(ValueError):
        f_res(m, 0.5, 0, 0)
    with pytest.raises(ValueError):
        f_res(m, 0.5, 0, 0, g=2.0, gamma_over_2pi_Hz_per_T=28e9)


def test_dimensions_of_outputs():
    # B = 1 T, free spin: ~ 28 GHz (Hz, not rad/s)
    f = f_res(MagnetModel(Ms=MS, N=sphere()), 1.0, 0, 0, g=2.0).modes[0].f_Hz
    assert 27.9e9 < f < 28.1e9
    # anisotropy field in tesla from J/m^3 and A/m
    assert abs(MagnetModel(Ms=MS, K1=SYNTH_K1).mu0_HA - 2 * SYNTH_K1 / MS) < 1e-15
    assert abs(MagnetModel(Ms=MS).mu0Ms - MU0 * MS) < 1e-15


def test_H_res_kittel_sphere_exact_root():
    m = MagnetModel(Ms=MS, N=sphere())
    r = H_res(m, 10e9, 0.3, 0.0, B_range_T=(0.0, 1.0), gyro=GY, n_points=41)
    assert r.status == "ok" and len(r.roots) == 1
    assert abs(r.roots[0].B_res_T - 10e9 / G2) < 1e-9
    assert abs(r.roots[0].df_dB_Hz_per_T - G2) < 1e-4 * G2
    assert abs(r.roots[0].H_res_A_per_m * MU0 - r.roots[0].B_res_T) < 1e-15


def test_H_res_no_root_returns_status_not_edge():
    m = MagnetModel(Ms=MS, N=sphere())
    r = H_res(m, 80e9, 0.3, 0.0, B_range_T=(0.0, 1.0), gyro=GY, n_points=41)      # needs B ~ 2.9 T
    assert r.status == "no_resonance_in_range" and r.roots == []
    r = H_res(m, 1e9, 0.3, 0.0, B_range_T=(0.5, 1.0), gyro=GY, n_points=21)       # below the range
    assert r.status == "no_resonance_in_range" and r.roots == []


def test_H_res_multiple_branches_flagged():
    HA = 2 * SYNTH_K1 / MS
    m = MagnetModel(Ms=MS, K1=SYNTH_K1, N=sphere())
    f_t = G2 * 0.5 * HA
    r = H_res(m, f_t, 0.0, 0.0, B_range_T=(-0.9 * HA, 0.9 * HA), gyro=GY, n_points=61)
    assert r.status == "ok" and len(r.roots) == 2
    Bs = sorted(x.B_res_T for x in r.roots)
    assert abs(Bs[0] + 0.5 * HA) < 1e-8 and abs(Bs[1] - 0.5 * HA) < 1e-8
    assert len({x.branch_id for x in r.roots}) == 2
    assert all(not x.is_ground for x in r.roots)              # both are metastable (antiparallel) states
    f_t2 = G2 * 1.5 * HA
    r2 = H_res(m, f_t2, 0.0, 0.0, B_range_T=(-0.9 * HA, 0.9 * HA), gyro=GY, n_points=61)
    assert r2.status == "ok" and len(r2.roots) == 2 and all(x.is_ground for x in r2.roots)
