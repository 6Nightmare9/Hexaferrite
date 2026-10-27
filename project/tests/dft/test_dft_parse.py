import numpy as np
import pytest

from hexw.dft import parse as P
from _dftsynth import synthetic_pw_out


def test_parse_basic_fields():
    r = P.parse_pw_output(synthetic_pw_out(-1194.89823655, [3.77, -3.77, 3.77, -3.77, 0.0, 0.0], wall="1m 10.40s"))
    assert r.converged and r.job_done and r.n_iterations == 11 and r.qe_version == "7.6"
    assert r.energy_Ry == pytest.approx(-1194.89823655) and r.energy_meV == pytest.approx(-1194.89823655 * 13605.693122994)
    assert r.site_moments == [3.77, -3.77, 3.77, -3.77, 0.0, 0.0]
    assert r.wall_s == pytest.approx(70.4) and r.ram_estimate_MB == pytest.approx(203.39)
    assert r.n_kpoints == 8 and r.n_electrons == 100.0 and r.nat == 6


def test_unconverged_has_no_energy():
    r = P.parse_pw_output(synthetic_pw_out(-1.0, [1.0], converged=False))
    assert not r.converged and r.energy_Ry is None and r.energy_meV is None        # an unconverged energy is not a result


def test_hours_in_wall_time_and_errors():
    r = P.parse_pw_output(synthetic_pw_out(-1.0, [1.0], wall="2h 3m 1.50s"))
    assert r.wall_s == pytest.approx(7381.5)
    txt = "     Error in routine read_cards (1):\n     bad card\n"
    assert P.parse_pw_output(txt).errors == ["read_cards: bad card"]


def test_real_output_snippet_hubbard_moments():
    txt = """     Tr[ns(  1)] (up, down, total) =   4.94572  0.83783  5.78355
     Tr[ns(  2)] (up, down, total) =   0.83783  4.94572  5.78355
"""
    r = P.parse_pw_output(txt)
    assert r.hubbard_d_moments[1] == pytest.approx(4.10789) and r.hubbard_d_moments[2] == pytest.approx(-4.10789)


def _check(moments, target, spin=None):
    spin = np.array([2.5, 2.5, 2.5, 2.5, 0.0]) if spin is None else spin
    r = P.parse_pw_output(synthetic_pw_out(-1.0, moments))
    return P.check_configuration(r, np.array(target, float), spin, min_fraction=0.5, start_moment=2 * spin)


def test_configuration_preserved():
    c = _check([3.7, -3.7, -3.7, 3.7, 0.0], [1, -1, -1, 1, 0])
    assert c["ok"] and c["reason"] == "preserved" and c["min_abs_moment"] == pytest.approx(3.7)


def test_flipped_moment_is_detected():
    c = _check([3.7, 3.6, -3.7, 3.7, 0.0], [1, -1, -1, 1, 0])
    assert not c["ok"] and c["flipped"] == [2]


def test_collapsed_moment_is_detected():
    c = _check([3.7, -3.7, -1.0, 3.7, 0.0], [1, -1, -1, 1, 0])                    # 1.0 < 0.5 * 5
    assert not c["ok"] and c["collapsed"] == [3] and not c["flipped"]


def test_missing_moments_cannot_be_verified():
    r = P.parse_pw_output("!    total energy              =  -1.0 Ry\n convergence has been achieved in 3 iterations\n JOB DONE.\n")
    c = P.check_configuration(r, np.array([1.0]), np.array([2.5]), min_fraction=0.5, start_moment=np.array([5.0]))
    assert not c["ok"] and c["reason"] == "no_site_moments"
