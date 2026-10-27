"""Grain ensemble, texture, stubs, Monte-Carlo, import hygiene. SYNTHETIC inputs."""
import math
import subprocess
import sys

import numpy as np
import pytest
from scipy import stats

from hexw.theory.magnet import (LEVEL_ENSEMBLE, LEVEL_IDEAL, MagnetModel, Texture, crystal_chi, f_res, find_equilibria,
                                intergranular_interaction, monte_carlo, passivity_margin, polycrystal_chi,
                                porosity_correction, rotation_from_c_axis, sample_demag_correction, sphere,
                                susceptibility)
from hexw.theory.magnet.polycrystal import orientation_grid, rotation_matrix
from _synth import GAMMA_2PI, GY, SYNTH_K1, SYNTH_MS

MS = SYNTH_MS


def _setup():
    # field above the anisotropy field (2K1/Ms = 0.667 T) -> a single minimum in every grain
    return MagnetModel(Ms=MS, K1=SYNTH_K1, N=sphere()), np.array([0.0, 0.0, 1.2]), np.linspace(10e9, 60e9, 7)


def test_levels_labelled():
    m, B, fr = _setup()
    r1 = crystal_chi(m, rotation_from_c_axis([0, 0, 1]), B, fr, alpha=0.2, gyro=GY)
    assert r1.level == LEVEL_IDEAL == "ideal_crystal"
    r2 = polycrystal_chi(m, B, fr, alpha=0.2, gyro=GY, n_c=30, state="field_aligned")
    assert r2.level == LEVEL_ENSEMBLE == "grain_ensemble_noninteracting"
    assert r2.n_failed == 0


def test_isotropic_average_converges_and_is_grid_independent():
    m, B, fr = _setup()
    ref = polycrystal_chi(m, B, fr, alpha=0.2, gyro=GY, n_c=1000, state="field_aligned").chi_lab
    scale = np.abs(ref).max()
    errs = []
    for n in (40, 160, 640):
        c = polycrystal_chi(m, B, fr, alpha=0.2, gyro=GY, n_c=n, state="field_aligned").chi_lab
        errs.append(np.abs(c - ref).max() / scale)
    assert errs[0] > errs[2] and errs[2] < 5e-3 and errs[1] < 2e-2
    # a differently rotated grid and a different azimuth sampling give the same average
    Rg = rotation_matrix([1, 2, 3], 0.7)
    c2 = polycrystal_chi(m, B, fr, alpha=0.2, gyro=GY, n_c=640, n_phi=1, state="field_aligned",
                         global_rotation=Rg).chi_lab
    c3 = polycrystal_chi(m, B, fr, alpha=0.2, gyro=GY, n_c=200, n_phi=3, state="field_aligned").chi_lab
    assert np.abs(c2 - ref).max() / scale < 5e-3
    assert np.abs(c3 - ref).max() / scale < 1e-2
    # isotropic average of a uniaxial grain with a lab field along z is isotropic in xx/yy
    assert abs(ref[3, 0, 0] - ref[3, 1, 1]) < 5e-3 * scale
    assert passivity_margin(ref) > -1e-9


def test_global_state_ensemble_runs_and_is_passive():
    m = MagnetModel(Ms=MS, K1=SYNTH_K1, N=sphere())
    r = polycrystal_chi(m, [0.0, 0.0, 0.3], np.linspace(5e9, 60e9, 5), alpha=0.2, gyro=GY, n_c=24)
    assert r.state_selection == "global" and passivity_margin(r.chi_lab) > -1e-9


def test_sample_demag_correction_inverse_of_internal_definition():
    m, B, fr = _setup()
    N = np.diag([0.2, 0.3, 0.5])
    R = np.eye(3)
    ci = crystal_chi(m.with_N(np.zeros((3, 3)) + np.eye(3) / 3), R, B, fr, alpha=0.1, gyro=GY, kind="internal")
    ce = sample_demag_correction(ci.chi_lab, np.eye(3) / 3)
    ref = crystal_chi(m, R, B, fr, alpha=0.1, gyro=GY, kind="external")
    assert np.allclose(ce, ref.chi_lab, rtol=1e-9, atol=1e-12)


def test_march_dollase_weights():
    t = Texture(r=0.5, axis=(0, 0, 1))
    c = np.array([[0, 0, 1.0], [1.0, 0, 0]])
    w = t.weights(c)
    assert abs(w.sum() - 1) < 1e-15
    assert abs(w[0] / w[1] - 0.5 ** -4.5) < 1e-9          # (r^2)^-1.5 / (1/r)^-1.5
    assert np.allclose(Texture(r=1.0).weights(c), 0.5)
    with pytest.raises(ValueError):
        Texture(r=0.0).weights(c)
    # strong texture moves the ensemble towards the single crystal along the axis
    m, B, fr = _setup()
    B = np.array([0.0, 0.0, 1.2])
    single = crystal_chi(m, rotation_from_c_axis([0, 0, 1]), B, fr, alpha=0.2, gyro=GY).chi_lab
    iso = polycrystal_chi(m, B, fr, alpha=0.2, gyro=GY, n_c=200, state="field_aligned").chi_lab
    tex = polycrystal_chi(m, B, fr, alpha=0.2, gyro=GY, n_c=200, state="field_aligned",
                          texture=Texture(r=0.1)).chi_lab
    assert np.abs(tex - single).max() < np.abs(iso - single).max()


def test_orientation_grid_proper_rotations():
    Rs, c = orientation_grid(20, 3)
    assert len(Rs) == 60
    for R in Rs[:10]:
        assert np.allclose(R @ R.T, np.eye(3), atol=1e-12) and abs(np.linalg.det(R) - 1) < 1e-12


def test_stubs_return_not_implemented_with_status():
    for fn in (intergranular_interaction, porosity_correction):
        r = fn(model=None)
        assert r.status == "not_implemented" and r.value is NotImplemented and r.reason


def test_monte_carlo_axial_frequency():
    HA_K1 = lambda K1: GAMMA_2PI * 2 * K1 / MS

    def func(p):
        m = MagnetModel(Ms=p["Ms"], K1=p["K1"], N=sphere())
        mode = f_res(m, 0.0, 0.0, 0.0, gamma_over_2pi_Hz_per_T=p["g2pi"], n_grid=100).modes[0]
        return {"f_Hz": mode.f_Hz}

    dists = {"K1": stats.norm(SYNTH_K1, 0.05 * SYNTH_K1), "Ms": MS, "g2pi": GAMMA_2PI}
    res = monte_carlo(func, dists, n=150, seed=7)
    s = res.summary["f_Hz"]
    assert s.n_valid == 150 and s.lo < s.median < s.hi
    assert abs(s.median - HA_K1(SYNTH_K1)) < 0.01 * s.median
    assert abs((s.hi - s.lo) / s.median - 2 * 1.96 * 0.05) < 0.03
    res2 = monte_carlo(func, dists, n=150, seed=7)
    assert res2.summary["f_Hz"].median == s.median                       # reproducible


def test_monte_carlo_counts_invalid_without_filling():
    def func(p):
        if p["x"] > 0.5:
            return {"y": None}
        if p["x"] < 0.1:
            raise ValueError("bad draw")
        return {"y": p["x"]}

    res = monte_carlo(func, {"x": stats.uniform(0, 1)}, n=400, seed=1)
    s = res.summary["y"]
    assert s.n_valid < 400 and s.fraction_valid < 0.6
    assert res.failures.get("ValueError", 0) > 0 and res.failures.get("y:invalid", 0) > 0
    assert np.isnan(res.samples["y"]).sum() == 400 - s.n_valid
    with pytest.raises(ValueError):
        monte_carlo(func, None, n=3)


def test_import_has_no_filesystem_writes_or_pools():
    code = (
        "import sys\n"
        "ev=[]\n"
        "def hook(e,a):\n"
        "    if e=='open':\n"
        "        p,mode=str(a[0]),a[1] or 'r'\n"
        "        if any(c in str(mode) for c in 'wax+') and not p.endswith(('.pyc','.pyo')): ev.append(p)\n"
        "    if e in ('os.mkdir','os.remove','shutil.rmtree','subprocess.Popen','multiprocessing.pool'): ev.append(e)\n"
        "sys.addaudithook(hook)\n"
        "import hexw.theory.magnet\n"
        "assert 'multiprocessing.pool' not in sys.modules, 'multiprocessing.pool imported'\n"
        "assert 'concurrent.futures.process' not in sys.modules, 'process pool imported'\n"
        "assert not ev, ev\n"
    )
    # -I ignores PYTHONDONTWRITEBYTECODE; disable bytecode explicitly so the check covers model I/O.
    r = subprocess.run([sys.executable, "-I", "-B", "-c", code], capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr
