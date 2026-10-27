"""Peaks and widths on constructed lines (Lorentzians). SYNTHETIC."""
import math

import numpy as np

from hexw.theory.magnet import (derivative_peak_to_peak, find_spectrum_peaks, fwhm, lorentzian_fwhm_from_dpp,
                                width_field_to_freq, width_freq_to_field)
from hexw.theory.magnet.linewidth import lorentzian_dpp_from_fwhm


def lor(x, x0, a, h=1.0):
    return h / (1.0 + ((x - x0) / a) ** 2)


def test_valid_lorentzian_fwhm():
    x = np.linspace(-10, 10, 4001)
    r = fwhm(x, lor(x, 0.7, 0.5), mode_x=0.7)
    assert r.status == "valid" and abs(r.width - 1.0) < 1e-4
    assert abs(r.peak_x - 0.7) < 1e-4 and abs(r.center_half_max - 0.7) < 1e-4
    assert r.n_points_in_fwhm > 100 and r.censored_sides == ()


def test_three_positions_are_distinct_for_asymmetric_line():
    x = np.linspace(-10, 10, 8001)
    y = lor(x, 0.0, 0.5) * (1 + 0.5 * np.tanh(x / 0.5))        # skewed line
    r = fwhm(x, y, mode_x=0.0)
    assert r.status == "valid"
    assert r.mode_x == 0.0
    assert abs(r.peak_x - r.center_half_max) > 1e-3 and abs(r.peak_x - r.mode_x) > 1e-3


def test_left_censored_no_edge_substitution():
    x = np.linspace(-1.2, 10, 3001)
    r = fwhm(x, lor(x, 0.0, 1.0))                    # left half-maximum point at -1 is... inside; shift the line
    assert r.status == "valid"
    x = np.linspace(-0.5, 10, 3001)                  # left half-max point (x=-1) outside the grid
    r = fwhm(x, lor(x, 0.0, 1.0))
    assert r.status == "left_censored" and r.width is None and r.x_left is None
    assert abs(r.x_right - 1.0) < 1e-3


def test_right_censored():
    x = np.linspace(-10, 0.5, 3001)
    r = fwhm(x, lor(x, 0.0, 1.0))
    assert r.status == "right_censored" and r.width is None and r.x_right is None
    assert abs(r.x_left + 1.0) < 1e-3


def test_both_censored_and_peak_on_edge():
    x = np.linspace(-0.2, 0.2, 401)
    assert fwhm(x, lor(x, 0.0, 1.0)).status == "both_censored"
    x = np.linspace(0.0, 10, 2001)
    r = fwhm(x, lor(x, -0.3, 1.0))                   # maximum at the grid edge
    assert r.status == "left_censored" and r.width is None


def test_multipeak():
    x = np.linspace(-20, 20, 8001)
    y = lor(x, -3.0, 0.5) + lor(x, 3.0, 0.5, 0.9)
    r = fwhm(x, y)
    assert r.status == "multipeak" and r.width is None and len(r.peaks) == 2
    # a weak distant satellite below half level does not spoil the main FWHM
    y = lor(x, 0.0, 0.5) + lor(x, 8.0, 0.5, 0.2)
    r = fwhm(x, y)
    assert r.status == "valid" and abs(r.width - 1.0) < 2e-3 and len(r.peaks) == 2


def test_insufficient_resolution_coarse_grid():
    x = np.linspace(-10, 10, 21)                    # step 1 with HWHM 0.5
    r = fwhm(x, lor(x, 0.0, 0.5))
    assert r.status == "insufficient_resolution" and r.width is None
    assert r.width_estimate is not None and r.n_points_in_fwhm < 5


def test_adaptive_expansion_callback():
    calls = []

    def spectrum(lo, hi):
        calls.append((lo, hi))
        x = np.linspace(lo, hi, 4001)
        return x, lor(x, 0.0, 3.0)

    x0, y0 = spectrum(-0.5, 0.5)                     # far too narrow a window
    calls.clear()
    r = fwhm(x0, y0, expand=spectrum, max_expand=8)
    assert r.status == "valid" and abs(r.width - 6.0) < 2e-2 and r.n_expansions >= 1 and calls
    # limited expansion stays censored (no fake width)
    r = fwhm(x0, y0, expand=spectrum, x_limits=(-1.0, 1.0))
    assert r.status == "both_censored" and r.width is None
    r = fwhm(x0, y0, expand=spectrum, max_expand=1)
    assert r.status in ("both_censored", "left_censored", "right_censored") and r.width is None


def test_find_peaks_ignores_noise_level():
    x = np.linspace(-10, 10, 4001)
    y = lor(x, 0, 0.5) + 0.004 * np.random.default_rng(3).normal(size=x.size)
    assert len(find_spectrum_peaks(x, y)) == 1


def test_dpp_and_lorentzian_relation():
    x = np.linspace(-20, 20, 40001)
    a = 0.8
    y = lor(x, 0.0, a)
    d = derivative_peak_to_peak(x, y)
    w = fwhm(x, y).width
    assert d.status == "valid"
    assert abs(d.delta_pp - w / math.sqrt(3)) < 1e-3 * w
    assert abs(lorentzian_fwhm_from_dpp(d.delta_pp) - w) < 1e-3 * w
    assert abs(lorentzian_dpp_from_fwhm(2 * a) - d.delta_pp) < 1e-3 * w
    # Gaussian is NOT sqrt(3): the Lorentzian relation must not be applied blindly
    s = 1.0
    g = np.exp(-x ** 2 / (2 * s ** 2))
    dg = derivative_peak_to_peak(x, g).delta_pp
    wg = fwhm(x, g).width
    assert abs(wg / dg - math.sqrt(3)) > 0.3


def test_dpp_censored():
    x = np.linspace(-0.1, 10, 3001)
    assert derivative_peak_to_peak(x, lor(x, 0.0, 0.5)).status == "left_censored"
    x = np.linspace(-10, 0.1, 3001)
    assert derivative_peak_to_peak(x, lor(x, 0.0, 0.5)).status == "right_censored"


def test_width_conversion_and_warnings():
    c = width_freq_to_field(1.4e9, 28e9)             # Hz, Hz/T
    assert abs(c.value - 0.05) < 1e-12 and c.status in ("ok", "ok_with_warnings") and c.warnings
    back = width_field_to_freq(c.value, 28e9)
    assert abs(back.value - 1.4e9) < 1e-3
    assert width_freq_to_field(1e9, 0.0).status == "slope_zero" and width_freq_to_field(1e9, 0.0).value is None
    # strongly curved branch: warning
    c = width_freq_to_field(5e9, 10e9, d2fdB2=500e9)
    assert c.status == "ok_with_warnings" and any("slope changes" in w for w in c.warnings)
    # negative slope gives a positive width
    assert width_freq_to_field(1e9, -20e9).value > 0
