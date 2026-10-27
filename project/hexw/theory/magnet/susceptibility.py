"""High-frequency susceptibility tensor chi(f, H) from the linearised Gilbert LLG.

CONVENTIONS (fixed, do not mix):
  * time dependence exp(-i omega t), omega = 2 pi f > 0.  chi = chi' + i chi''.  A passive medium has
    chi'' > 0 on the diagonal for omega > 0; in general the anti-Hermitian part
    (chi - chi^dagger)/(2i) is positive semidefinite (absorbed power  P = (omega mu0 / 2) h^dagger [...] h).
  * Linear response  M1 = chi H1  where H1 is the small AC field APPLIED from outside (vacuum field, uniform
    sample, shape through the static N). kind='external' (default): demag fields of the dynamic magnetisation
    are included (N enters the Hessian) -> chi_ext (a property of material AND shape).
    kind='internal': M1 = chi_int (H1 - N M1), the response to the field inside the sample ignoring dynamic
    demag:  chi_int = chi_ext (I - N chi_ext)^-1,  and  chi_ext = chi_int (I + N chi_int)^-1.
  * Tensors are returned in the CRYSTAL frame (z = c); chi is rank 2 (no response along m0 in linear order).

Equation in the local frame (e1 x e2 = m0, W = Hs/Ms [T], R = [[0,-1],[1,0]]):
    [ i omega (I - alpha R) + gamma R W ] delta = gamma R h,   M1 = Ms delta,  h = mu0 H1 (tesla),
    chi = mu0 Ms * E A E^T,  A = gamma [ ... ]^-1 R.
alpha is the Gilbert damping, a mandatory argument (no default).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .core import MU0, Gyro, MagnetModel, sphere_derivs, unit

R2 = np.array([[0.0, -1.0], [1.0, 0.0]])


def susceptibility(model: MagnetModel, B, m0, freqs_Hz, *, alpha: float, gyro: Gyro,
                   kind: str = "external") -> np.ndarray:
    """chi(f) at the equilibrium m0 for field vector B (tesla, crystal frame); freqs in Hz. Shape (nf, 3, 3).
    m0 must be a stable equilibrium (checked loosely: gradient and Hessian)."""
    if alpha is None or alpha < 0:
        raise ValueError("alpha (Gilbert damping, >= 0) is mandatory")
    freqs = np.atleast_1d(np.asarray(freqs_Hz, float))
    m0 = unit(m0)
    gt, Hs, E = sphere_derivs(model, B, m0)
    scale = model.energy_scale(float(np.linalg.norm(B)))
    if np.linalg.norm(gt) > 1e-6 * scale:
        raise ValueError("m0 is not an equilibrium (tangential gradient not zero)")
    if np.linalg.eigvalsh(Hs)[0] < -1e-9 * scale:
        raise ValueError("m0 is not a stable minimum")
    W = Hs / model.Ms
    g = gyro.gamma_rad_per_s_T
    om = 2.0 * math.pi * freqs
    I2 = np.eye(2)
    Mm = 1j * om[:, None, None] * (I2 - alpha * R2)[None] + (g * (R2 @ W))[None]
    A = g * np.linalg.solve(Mm, np.broadcast_to(R2, Mm.shape).astype(complex))
    chi = MU0 * model.Ms * np.einsum("ia,fab,jb->fij", E, A, E)
    if kind == "external":
        return chi
    if kind == "internal":
        I3 = np.eye(3)
        return np.array([c @ np.linalg.inv(I3 - model.N @ c) for c in chi])
    raise ValueError("kind must be 'external' or 'internal'")


def mu_r(chi_int: np.ndarray, kind: str = "internal") -> np.ndarray:
    """Material relative permeability  mu_r = I + chi_int  (chi_int: response to the field INSIDE the body, no
    sample shape). kind='external' is refused: I + chi_ext mixes material and shape; use `apparent_mu_shape`
    (with an explicit shape label) or `apparent_mu_ellipsoid`."""
    if kind != "internal":
        raise ValueError("mu_r is defined for chi_int (material). I + chi_ext depends on the shape: use "
                         "apparent_mu_shape(chi_ext, shape_label) or apparent_mu_ellipsoid(chi_ext, N)")
    return np.eye(3)[None] + chi_int


@dataclass(frozen=True)
class ApparentMu:
    """A shape-dependent 'apparent permeability' together with the shape it is valid for."""
    mu: np.ndarray
    shape_label: str
    definition: str


def apparent_mu_shape(chi_ext: np.ndarray, shape_label: str) -> ApparentMu:
    """I + chi_ext: ratio of the dynamic magnetisation-plus-field response to the APPLIED field for the particular
    sample shape whose demagnetising tensor is already inside chi_ext (kind='external'). This is NOT a material
    property; the mandatory `shape_label` (e.g. 'sphere', 'thin_disk_perp_c', 'ellipsoid 1:2:3') is carried with
    the result. For the flux density averaged over an ellipsoid use `apparent_mu_ellipsoid`."""
    if not isinstance(shape_label, str) or not shape_label.strip():
        raise ValueError("shape_label (non-empty string naming the sample shape) is mandatory")
    return ApparentMu(np.eye(3)[None] + chi_ext, shape_label.strip(), "I + chi_ext (shape-dependent, not material)")


def apparent_mu_ellipsoid(chi_ext: np.ndarray, N: np.ndarray, shape_label: str = "ellipsoid") -> ApparentMu:
    """Apparent permeability of a uniformly magnetised ELLIPSOID: <B> / (mu0 H_app) = I + (I - N) chi_ext.

    Derivation: M = chi_ext H_app, H_int = H_app - N M, B = mu0 (H_int + M) = mu0 (H_app + (I - N) M). For a sphere
    (N = I/3) this is I + (2/3) chi_ext, not I + chi_ext. chi_ext must be the external susceptibility of the same
    ellipsoid (N in the model equal to this N)."""
    N = np.asarray(N, float)
    mu = np.eye(3)[None] + np.einsum("ij,fjk->fik", np.eye(3) - N, chi_ext)
    return ApparentMu(mu, shape_label, "<B>/(mu0 H_app) = I + (I - N) chi_ext (ellipsoid)")


def antihermitian_part(chi: np.ndarray) -> np.ndarray:
    return (chi - np.conj(np.swapaxes(chi, -1, -2))) / 2j


def passivity_margin(chi: np.ndarray) -> float:
    """min over the grid of the smallest eigenvalue of (chi - chi^dagger)/(2i), divided by max |chi|.
    Passive (omega > 0, exp(-i omega t)):  margin >= -tol (rounding)."""
    ah = antihermitian_part(chi)
    ev = np.linalg.eigvalsh(ah)
    return float(ev.min() / max(np.abs(chi).max(), 1e-300))


def kramers_kronig_real(omega: np.ndarray, chi_imag: np.ndarray) -> np.ndarray:
    """Real part from the imaginary part of a causal response with chi(-w) = conj chi(w):
        chi'(w) = (2/pi) P int_0^inf w' chi''(w') / (w'^2 - w^2) dw'
    evaluated by singularity subtraction on the (non-uniform) positive grid `omega` (trapezoid rule).
    O(n^2) cost. The integral is truncated at omega[-1] (tail neglected: caller must provide a wide grid; the
    analytic piece of the principal value on [0, Omega] is kept exactly)."""
    w = np.asarray(omega, float)
    y = np.asarray(chi_imag, float)
    out = np.empty_like(y)
    Om = w[-1]
    for i, wi in enumerate(w):
        num = w * y
        den = w ** 2 - wi ** 2
        sing = den == 0.0
        integrand = np.where(sing, 0.0, (num - wi * y[i]) / np.where(sing, 1.0, den))
        if sing.any():              # derivative limit at w' = w (smooth part), central estimate
            j = int(np.nonzero(sing)[0][0])
            jm, jp = max(j - 1, 0), min(j + 1, len(w) - 1)
            integrand[j] = 0.5 * (integrand[jm] + integrand[jp])
        smooth = np.trapezoid(integrand, w)
        if wi >= Om:
            out[i] = np.nan          # no data above the last point to close the principal value
            continue
        pv = 0.5 * y[i] * math.log(abs((Om - wi) / (Om + wi))) if wi > 0 else 0.0
        out[i] = (2.0 / math.pi) * (smooth + pv)
    return out
