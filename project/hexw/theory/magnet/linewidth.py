"""Peaks, FWHM with censoring statuses, derivative peak-to-peak width, frequency <-> field conversion.

Three DIFFERENT positions are kept apart in LinewidthResult:
    mode_x         eigenfrequency / resonance field of the mode (from the model; given by the caller)
    peak_x         maximum of the (possibly broadened) spectrum (parabolic interpolation)
    center_half_max  midpoint of the half-maximum interval  (x_left + x_right)/2
For an asymmetric or overlapping line they differ.

NOISE. The white-noise level is estimated from the MAD of the second difference (`noise_sigma`). A "peak" must have a
prominence >= max(prominence_rel * max, NOISE_K * sigma) and be >= `min_separation` away from a higher one (default
3 samples), so noise ripples on top of a line are not peaks. With smooth='auto' a cubic Savitzky-Golay filter (window
~ FWHM/4 samples) is applied when the noise is noticeable (sigma/max > 2e-3): the maximum of a noisy line is otherwise
biased UPWARD and the half level with it. The window actually used is returned.
`asymmetric_possible_doublet` is set when |center_half_max - peak_x| > 10% of the width: a skewed line or an
unresolved doublet -- the width is then not the width of one symmetric line (a symmetric unresolved doublet is not
detectable by this test).

Statuses of fwhm(): valid | left_censored | right_censored | both_censored | multipeak |
insufficient_resolution | no_peak.  (both_censored and no_peak are additions to the five required ones.)
Precedence: censored > multipeak > insufficient_resolution > valid. A grid edge is NEVER substituted for a
half-maximum crossing; width is None unless status == 'valid'. For insufficient_resolution the interpolated
interval is given as `width_estimate` only.

Delta B_pp (distance between the extrema of d(absorption)/dx) is a different quantity from FWHM:
for a LORENTZIAN absorption line only, FWHM = sqrt(3) * Delta_pp. For other shapes the factor differs
(Gaussian: FWHM = 1.177 Delta_pp); no universal relation exists.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from scipy.signal import find_peaks, savgol_filter

NOISE_K = 5.0                          # prominence must exceed NOISE_K x sigma
NOISE_REL_FOR_SMOOTHING = 2e-3         # sigma / max above which smooth='auto' smooths
DOUBLET_ASYMMETRY_FRACTION = 0.1
SLOPE_FLOOR_HZ_PER_T = 1.0e6           # |df/dB| below this (1 MHz/T, ~4e-5 of the free-spin 28 GHz/T) is "too small"


@dataclass
class LinewidthResult:
    status: str
    peak_x: float | None = None
    peak_y: float | None = None
    half_level: float | None = None
    x_left: float | None = None
    x_right: float | None = None
    width: float | None = None
    width_estimate: float | None = None
    center_half_max: float | None = None
    mode_x: float | None = None
    n_points_in_fwhm: int = 0
    censored_sides: tuple = ()
    peaks: list = field(default_factory=list)
    n_expansions: int = 0
    x_range: tuple | None = None
    notes: list = field(default_factory=list)
    noise_sigma: float | None = None
    smoothed_window: int = 0                # 0 = no smoothing
    asymmetric_possible_doublet: bool = False


def noise_sigma(y) -> float:
    """White-noise standard deviation from the MAD of the second difference (var(d2) = 6 sigma^2); robust against
    the smooth line itself."""
    y = np.asarray(y, float)
    if len(y) < 5:
        return 0.0
    d2 = y[2:] - 2.0 * y[1:-1] + y[:-2]
    return float(1.4826 * np.median(np.abs(d2 - np.median(d2))) / math.sqrt(6.0))


def _odd_window(w: int, n: int) -> int:
    """Odd window >= 5 and <= n (SG needs window <= len)."""
    w = max(5, int(w))
    if w % 2 == 0:
        w += 1
    top = n if n % 2 == 1 else n - 1
    return max(5, min(w, top))


def _sg(y, window: int):
    w = _odd_window(window, len(y))
    return savgol_filter(y, w, 3), w


def _parabola(x, y, i):
    """Vertex of the parabola through (i-1, i, i+1); falls back to the sample when degenerate."""
    if i <= 0 or i >= len(x) - 1:
        return float(x[i]), float(y[i])
    x0, x1, x2 = x[i - 1], x[i], x[i + 1]
    y0, y1, y2 = y[i - 1], y[i], y[i + 1]
    d = (x0 - x1) * (x0 - x2) * (x1 - x2)
    if d == 0:
        return float(x1), float(y1)
    A = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / d
    Bq = (x2 ** 2 * (y0 - y1) + x1 ** 2 * (y2 - y0) + x0 ** 2 * (y1 - y2)) / d
    if A >= 0:
        return float(x1), float(y1)
    xv = -Bq / (2 * A)
    if not (x0 <= xv <= x2):
        return float(x1), float(y1)
    yv = A * xv ** 2 + Bq * xv + (y1 - A * x1 ** 2 - Bq * x1)
    return float(xv), float(yv)


def find_spectrum_peaks(x, y, baseline: float = 0.0, prominence_rel: float = 0.1, sigma: float | None = None,
                        min_separation: float | None = None) -> list:
    """Significant maxima [(x_peak, y_peak_above_baseline)] with prominence >= max(prominence_rel * global maximum,
    NOISE_K * sigma) (sigma: noise level, estimated by `noise_sigma` if None) and mutual distance >= min_separation
    (x units; default 3 samples). Grid-edge samples are not reported as peaks unless they are local maxima of an
    interior shoulder."""
    x = np.asarray(x, float)
    yb = np.asarray(y, float) - baseline
    ymax = yb.max()
    if ymax <= 0:
        return []
    sig = noise_sigma(yb) if sigma is None else sigma
    dist = 3 if min_separation is None else max(1, int(round(min_separation / float(np.median(np.diff(x))))))
    idx, _ = find_peaks(yb, prominence=max(prominence_rel * ymax, NOISE_K * sig), distance=dist)
    out = []
    for i in idx:
        xv, yv = _parabola(x, yb, int(i))
        out.append((xv, yv))
    return out


def _check_grid(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if x.shape != y.shape or x.ndim != 1 or len(x) < 3:
        raise ValueError("x, y must be 1-D arrays of equal length >= 3")
    if np.any(np.diff(x) <= 0):
        raise ValueError("x must be strictly increasing")
    return x, y


def _evaluate(x, y, baseline, min_points, prominence_rel, mode_x, sigma=None) -> LinewidthResult:
    yb = y - baseline
    res = LinewidthResult(status="no_peak", mode_x=mode_x, x_range=(float(x[0]), float(x[-1])))
    ip = int(np.argmax(yb))
    if yb[ip] <= 0:
        return res
    px, py = _parabola(x, yb, ip)
    half = 0.5 * py
    res.peak_x, res.peak_y, res.half_level = px, py, half
    # crossings (no substitution of edges)
    left = None
    i = ip
    while i > 0:
        if yb[i - 1] < half <= yb[i]:
            left = x[i - 1] + (half - yb[i - 1]) * (x[i] - x[i - 1]) / (yb[i] - yb[i - 1])
            break
        i -= 1
    right = None
    i = ip
    while i < len(x) - 1:
        if yb[i + 1] < half <= yb[i]:
            right = x[i] + (half - yb[i]) * (x[i + 1] - x[i]) / (yb[i + 1] - yb[i])
            break
        i += 1
    cens = tuple(s for s, v in (("left", left), ("right", right)) if v is None)
    res.censored_sides = cens
    pk = find_spectrum_peaks(x, y, baseline, prominence_rel, sigma=sigma)
    res.peaks = [p[0] for p in pk]
    dx_main = float(np.max(np.diff(x)))
    others = [p for p in pk if abs(p[0] - px) > 2 * dx_main and p[1] >= half]
    # a second significant maximum reaching the half level means the half-maximum interval is contaminated
    if len(cens) == 2:
        res.status = "both_censored"
        return res
    if len(cens) == 1:
        res.status = f"{cens[0]}_censored"
        if left is not None:
            res.x_left = float(left)
        if right is not None:
            res.x_right = float(right)
        return res
    res.x_left, res.x_right = float(left), float(right)
    res.center_half_max = 0.5 * (res.x_left + res.x_right)
    if others:
        res.status = "multipeak"
        res.notes.append(f"{len(others)} other significant maxima at/above half level")
        return res
    n_in = int(np.count_nonzero((x >= left) & (x <= right)))
    res.n_points_in_fwhm = n_in
    if n_in < min_points:
        res.status = "insufficient_resolution"
        res.width_estimate = float(right - left)
        res.notes.append(f"only {n_in} grid points inside the FWHM interval (< {min_points})")
        return res
    res.status = "valid"
    res.width = float(right - left)
    if abs(res.center_half_max - px) > DOUBLET_ASYMMETRY_FRACTION * res.width:
        res.asymmetric_possible_doublet = True
        res.notes.append("asymmetric_possible_doublet: |center_half_max - peak_x| > 10% of the width (skewed line or "
                         "unresolved doublet): the width is not the width of one symmetric line")
    return res


def _fwhm_samples(x, y, baseline) -> int:
    """Approximate number of samples inside the FWHM from a heavily pre-smoothed pass (for choosing windows)."""
    n = len(y)
    pre, _ = _sg(y, max(5, n // 20))
    r0 = _evaluate(x, pre, baseline, 1, 0.1, None, sigma=None)
    if r0.x_left is not None and r0.x_right is not None:
        return max(5, int(np.count_nonzero((x >= r0.x_left) & (x <= r0.x_right))))
    return max(5, n // 20)


def _eval_smoothed(x, y, baseline, min_points, prominence_rel, mode_x, smooth):
    sigma = noise_sigma(y)
    ymax = float(np.max(y - baseline))
    win, ys = 0, y
    if smooth == "auto":
        if ymax > 0 and sigma / ymax > NOISE_REL_FOR_SMOOTHING and len(y) >= 25:
            ys, win = _sg(y, max(5, _fwhm_samples(x, y, baseline) // 4))
    elif smooth is not None:
        ys, win = _sg(y, int(smooth))
    return _evaluate(x, ys, baseline, min_points, prominence_rel, mode_x, sigma=sigma), sigma, win


def fwhm(x, y, *, baseline: float = 0.0, mode_x: float | None = None, min_points: int = 5,
         prominence_rel: float = 0.1,
         expand: Callable[[float, float], tuple] | None = None, max_expand: int = 6,
         expand_factor: float = 2.0, x_limits: tuple | None = None, smooth="auto") -> LinewidthResult:
    """FWHM of the main maximum of y(x) above `baseline` (default 0, i.e. an absorption-type spectrum).

    expand(x_min, x_max) -> (x, y): optional spectrum callback. If a half-maximum crossing is missing on one
    side, the range is extended on THAT side by (expand_factor - 1) x current span (clipped to x_limits) and
    the spectrum is re-requested, up to max_expand times. If the limit is reached the status stays censored.

    smooth: 'auto' (Savitzky-Golay only if the noise is noticeable, sigma/max > 2e-3), None (never), or a window
    length in samples. The noise estimate and the window used are returned in the result."""
    x, y = _check_grid(x, y)
    n_exp = 0
    notes = []
    res, sigma, win = _eval_smoothed(x, y, baseline, min_points, prominence_rel, mode_x, smooth)
    while expand is not None and res.status in ("left_censored", "right_censored", "both_censored") \
            and n_exp < max_expand:
        lo, hi = float(x[0]), float(x[-1])
        span = hi - lo
        nlo, nhi = lo, hi
        if "left" in res.censored_sides:
            nlo = lo - (expand_factor - 1.0) * span
        if "right" in res.censored_sides:
            nhi = hi + (expand_factor - 1.0) * span
        if x_limits is not None:
            nlo, nhi = max(nlo, x_limits[0]), min(nhi, x_limits[1])
        if nlo == lo and nhi == hi:
            notes.append("x_limits reached; range cannot be extended")
            break
        x, y = _check_grid(*expand(nlo, nhi))
        n_exp += 1
        res, sigma, win = _eval_smoothed(x, y, baseline, min_points, prominence_rel, mode_x, smooth)
    res.n_expansions = n_exp
    res.noise_sigma, res.smoothed_window = sigma, win
    if win:
        notes.append(f"Savitzky-Golay smoothing (window {win} samples) applied; noise sigma = {sigma:.3g}")
    res.notes.extend(notes)
    return res


# --------------------------------------------------------------------------------- derivative peak-to-peak

@dataclass
class DerivativeWidth:
    status: str                   # valid | left_censored | right_censored | insufficient_resolution | no_peak |
    #                               noisy_unreliable
    x_max_slope: float | None = None
    x_min_slope: float | None = None
    delta_pp: float | None = None
    noise_sigma: float | None = None
    smoothed_window: int = 0
    estimates: tuple = ()         # delta_pp obtained with the two smoothing widths (noisy data)
    notes: list = field(default_factory=list)


def _dpp_core(x, y, baseline) -> DerivativeWidth:
    yb = y - baseline
    if yb.max() <= 0:
        return DerivativeWidth("no_peak")
    ip = int(np.argmax(yb))
    d = np.gradient(yb, x)
    il = int(np.argmax(d[: ip + 1])) if ip > 0 else 0
    ir = ip + int(np.argmin(d[ip:]))
    if il == 0 and d[0] >= d[1]:
        return DerivativeWidth("left_censored")
    if ir == len(x) - 1 and d[-1] <= d[-2]:
        return DerivativeWidth("right_censored")
    xl, _ = _parabola(x, d, il)
    xr, _ = _parabola(x, -d, ir)
    out = DerivativeWidth("valid", xl, xr, float(xr - xl))
    if ir - il < 3:
        out.status = "insufficient_resolution"
        out.delta_pp = None
    return out


def derivative_peak_to_peak(x, y, baseline: float = 0.0, smooth="auto", stability_tol: float = 0.1) -> DerivativeWidth:
    """Distance between the maximum (left of the line peak) and minimum (right) of dy/dx for an
    absorption-type line y(x): Delta_pp (e.g. Delta H_pp). Extrema at the grid edge -> censored.

    A derivative amplifies noise. If the noise is noticeable (sigma/max > 1e-3; smooth='auto') the line is smoothed
    (cubic Savitzky-Golay) with TWO window widths (w and 2w+1, w ~ FWHM/4 samples); if the two Delta_pp differ by more
    than `stability_tol` (relative) the status is 'noisy_unreliable' and no value is given (the two estimates are kept in
    `estimates`). Clean data are differentiated directly."""
    x, y = _check_grid(x, y)
    yb = y - baseline
    if yb.max() <= 0:
        return DerivativeWidth("no_peak")
    sigma = noise_sigma(y)
    rel = sigma / float(yb.max())
    if smooth is None or (smooth == "auto" and rel <= 0.5 * NOISE_REL_FOR_SMOOTHING) or len(y) < 25:
        out = _dpp_core(x, y, baseline)
        out.noise_sigma = sigma
        return out
    w1 = int(smooth) if smooth != "auto" else max(9, _fwhm_samples(x, y, baseline) // 4)
    outs = []
    for w in (w1, 2 * w1 + 1):
        ys, wu = _sg(y, w)
        o = _dpp_core(x, ys, baseline)
        o.smoothed_window = wu
        outs.append(o)
    res = outs[0]
    res.noise_sigma = sigma
    for o in outs:
        if o.status in ("left_censored", "right_censored", "no_peak"):
            o.noise_sigma = sigma
            return o
    est = tuple(o.delta_pp for o in outs)
    res.estimates = est
    if any(e is None for e in est) or abs(est[0] - est[1]) > stability_tol * max(abs(est[1]), 1e-300):
        res.status = "noisy_unreliable"
        res.delta_pp = None
        res.notes.append(f"Delta_pp unstable against the smoothing width ({est}); noise sigma = {sigma:.3g}")
    return res


def lorentzian_fwhm_from_dpp(dpp: float) -> float:
    """FWHM = sqrt(3) * Delta_pp. ONLY for a Lorentzian absorption line."""
    return math.sqrt(3.0) * dpp


def lorentzian_dpp_from_fwhm(w: float) -> float:
    return w / math.sqrt(3.0)


# --------------------------------------------------------------------------------- frequency <-> field

@dataclass
class WidthConversion:
    value: float | None
    status: str                   # ok | ok_with_warnings | slope_zero | slope_too_small
    warnings: list = field(default_factory=list)


_APPLICABILITY = ("conversion assumes the line is a rigid shift along one branch f(B) with a locally constant "
                  "slope; it does not convert intrinsic (Gilbert) broadening between swept variables for "
                  "strongly non-linear branches or lines superposed from several grains/modes")


def _convert(width_in, slope, d2, direction, curv_tol, slope_tol, slope_floor):
    if not (slope is not None and math.isfinite(slope)) or slope == 0.0:
        return WidthConversion(None, "slope_zero", ["df/dB = 0 (turning point of the branch): conversion undefined"])
    floor = max(slope_tol, slope_floor if direction == "f->B" else 0.0)
    if abs(slope) <= floor:
        return WidthConversion(None, "slope_too_small",
                               [f"|df/dB| = {abs(slope):.3g} Hz/T <= {floor:.3g} Hz/T: the conversion would blow up "
                                "(close to a turning point of the branch)"])
    wr = abs(slope)
    out = width_in / wr if direction == "f->B" else width_in * wr
    warns = [_APPLICABILITY]
    wB = width_in / wr if direction == "f->B" else width_in
    if d2 is not None and abs(d2) * wB / wr > curv_tol:
        warns.append(f"slope changes by {abs(d2) * wB / wr:.2f} (> {curv_tol}) across the line: linear conversion unreliable")
    return WidthConversion(float(out), "ok_with_warnings" if len(warns) > 1 else "ok", warns)


def width_freq_to_field(width_Hz: float, dfdB_Hz_per_T: float, d2fdB2: float | None = None,
                        curvature_tol: float = 0.1, slope_tol: float = 0.0,
                        slope_floor: float = SLOPE_FLOOR_HZ_PER_T) -> WidthConversion:
    """Delta(mu0 H) [T] = Delta f / |df/dB|. Warns if the branch curvature d2f/dB2 changes the slope by more
    than curvature_tol (relative) across the line. df/dB = 0 -> status slope_zero; |df/dB| <= max(slope_tol,
    slope_floor) (default floor 1 MHz/T, applied even when slope_tol = 0) -> slope_too_small. No value in either case."""
    return _convert(width_Hz, dfdB_Hz_per_T, d2fdB2, "f->B", curvature_tol, slope_tol, slope_floor)


def width_field_to_freq(width_T: float, dfdB_Hz_per_T: float, d2fdB2: float | None = None,
                        curvature_tol: float = 0.1, slope_tol: float = 0.0) -> WidthConversion:
    """Delta f [Hz] = |df/dB| Delta(mu0 H) (slope 0 -> slope_zero; |slope| <= slope_tol -> slope_too_small)."""
    return _convert(width_T, dfdB_Hz_per_T, d2fdB2, "B->f", curvature_tol, slope_tol, 0.0)
