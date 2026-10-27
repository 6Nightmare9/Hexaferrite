"""Equilibria, stability, classification, cone, branch following. All parameters are SYNTHETIC."""
import math

import numpy as np
import pytest

from hexw.theory.magnet import (MagnetModel, cone_angle_analytic, ellipsoid, field_vector, find_equilibria,
                                hysteresis_loop, sphere, thin_disk_parallel_c, thin_disk_perp_c, MU0)
from _synth import SYNTH_K1, SYNTH_K2, SYNTH_K4, SYNTH_MS


def test_demag_tensors_trace_one():
    for N in (sphere(), thin_disk_perp_c(), thin_disk_parallel_c(), ellipsoid(1, 2, 3), ellipsoid(1, 1, 50)):
        assert abs(np.trace(N) - 1) < 1e-9
    assert np.allclose(ellipsoid(1, 1, 1), np.eye(3) / 3, atol=1e-9)


def test_ellipsoid_prolate_analytic():
    m = 4.0                       # c/a
    e = math.sqrt(1 - 1 / m ** 2)
    Nz = (1 - e ** 2) / e ** 3 * (math.atanh(e) - e)
    N = ellipsoid(1, 1, m)
    assert abs(N[2, 2] - Nz) < 1e-9
    assert abs(N[0, 0] - (1 - Nz) / 2) < 1e-9


def test_ellipsoid_flat_limit_approaches_disk():
    N = ellipsoid(100, 100, 1)
    assert abs(N[2, 2] - 1) < 0.02 and N[0, 0] < 0.01


def test_cone_analytic():
    K1, K2 = -1.0e5, 3.0e5       # synthetic: K1<0, K2>0, sin^2 = -K1/(2K2) = 1/6
    m = MagnetModel(Ms=SYNTH_MS, K1=K1, K2=K2, N=sphere())
    eqs = find_equilibria(m, np.zeros(3))
    th0 = cone_angle_analytic(K1, K2)
    assert th0 is not None
    assert eqs.ground.kind == "cone"
    assert abs(math.sin(eqs.ground.theta) ** 2 - (-K1 / (2 * K2))) < 1e-10
    assert abs(eqs.ground.theta - th0) < 1e-10 or abs(math.pi - eqs.ground.theta - th0) < 1e-10
    assert eqs.ground.azimuthal_degenerate and eqs.ground.neutral        # ring: Goldstone direction
    # the axial and planar states are not minima here
    assert all(e.kind == "cone" for e in eqs.minima)


def test_cone_analytic_none_outside_region():
    assert cone_angle_analytic(1e5, 1e5) is None
    assert cone_angle_analytic(-1e5, 4e4) is None       # -K1/(2K2) >= 1


@pytest.mark.parametrize("K1,kind", [(+SYNTH_K1, "axial"), (-SYNTH_K1, "planar")])
def test_axial_planar_switch_with_sign_of_K1(K1, kind):
    m = MagnetModel(Ms=SYNTH_MS, K1=K1, N=sphere())
    eqs = find_equilibria(m, np.zeros(3))
    assert eqs.ground.kind == kind
    assert eqs.status == "ok"


def test_shape_overrules_weak_axial_anisotropy_for_disk():
    # mu0 Ms^2/2 = 0.5*MU0*9e10 ~ 5.7e4 < K1 -> axial survives; with K1 smaller the disk is planar
    weak = MagnetModel(Ms=SYNTH_MS, K1=0.2 * 0.5 * MU0 * SYNTH_MS ** 2, N=thin_disk_perp_c())
    assert find_equilibria(weak, np.zeros(3)).ground.kind == "planar"
    strong = MagnetModel(Ms=SYNTH_MS, K1=3 * 0.5 * MU0 * SYNTH_MS ** 2, N=thin_disk_perp_c())
    assert find_equilibria(strong, np.zeros(3)).ground.kind == "axial"


def test_metastable_state_and_energies_axial_antiparallel():
    m = MagnetModel(Ms=SYNTH_MS, K1=SYNTH_K1, N=sphere())
    HA = 2 * SYNTH_K1 / SYNTH_MS
    B = field_vector(-0.5 * HA, 0.0)              # field antiparallel to +c, below the switching field
    eqs = find_equilibria(m, B)
    assert len(eqs.minima) == 2 and len(eqs.metastable) == 1
    gs, (ms_, dE) = eqs.ground, eqs.metastable[0]
    assert gs.m[2] < -0.99 and ms_.m[2] > 0.99
    assert abs(dE - 2 * SYNTH_MS * 0.5 * HA) < 1e-6 * SYNTH_K1      # E(+z) - E(-z) = 2 Ms |B|
    assert all(e.stable for e in eqs.minima)
    # above the switching field the metastable state disappears (stability check of the Hessian)
    eqs2 = find_equilibria(m, field_vector(-1.2 * HA, 0.0))
    assert len(eqs2.minima) == 1 and not eqs2.metastable


def test_tilted_by_field_classification():
    m = MagnetModel(Ms=SYNTH_MS, K1=SYNTH_K1, N=sphere())
    eqs = find_equilibria(m, field_vector(0.1, math.pi / 3))
    assert eqs.ground.kind == "tilted_by_field"
    # planar easy-plane with in-plane field stays planar
    mp = MagnetModel(Ms=SYNTH_MS, K1=-SYNTH_K1, N=sphere())
    assert find_equilibria(mp, field_vector(0.1, math.pi / 2)).ground.kind == "planar"


def test_hexagonal_basal_minima_count_and_energy():
    # K1<0 (easy plane) and K4<0: minima at phi = k*pi/3 (cos6phi = 1 is the minimum for K4<0) -> 6 minima
    m = MagnetModel(Ms=SYNTH_MS, K1=-SYNTH_K1, K4=SYNTH_K4, N=sphere())
    eqs = find_equilibria(m, np.zeros(3))
    assert len(eqs.minima) == 6 and all(e.kind == "planar" for e in eqs.minima)
    assert abs(eqs.ground.energy - (-SYNTH_K1 + SYNTH_K4 + MU0 * SYNTH_MS ** 2 / 6)) < 1e-6 * SYNTH_K1
    phis = sorted(round(e.phi / (math.pi / 3), 6) % 6 for e in eqs.minima)
    assert np.allclose(phis, [0, 1, 2, 3, 4, 5], atol=1e-5)


def test_isotropic_degenerate_status():
    m = MagnetModel(Ms=SYNTH_MS, N=sphere())
    eqs = find_equilibria(m, np.zeros(3))
    assert eqs.status == "degenerate_isotropic" and eqs.ground.neutral


def test_hysteresis_switching_field_axial():
    HA = 2 * SYNTH_K1 / SYNTH_MS
    m = MagnetModel(Ms=SYNTH_MS, K1=SYNTH_K1, N=sphere())
    Bs, pts, proj = hysteresis_loop(m, [0, 0, 1], 1.5 * HA, 301)
    jumps = [(Bs[i], p.angle_step) for i, p in enumerate(pts) if p.jumped]
    assert len(jumps) == 2                                   # one switching on each sweep
    assert abs(jumps[0][0] + HA) < 0.02 * HA                 # down sweep switches near -HA
    assert abs(jumps[1][0] - HA) < 0.02 * HA                 # up sweep near +HA
    # the branch is remembered: at B=0 on the down sweep m is still +c, on the up sweep -c
    i0 = int(np.argmin(np.abs(Bs[:301])))
    assert proj[i0] > 0.99
    j0 = 301 + int(np.argmin(np.abs(Bs[301:])))
    assert proj[j0] < -0.99
