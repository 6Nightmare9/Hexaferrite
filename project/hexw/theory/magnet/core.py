"""Macrospin free energy, gradient and sphere-Hessian (crystal frame, z = c axis).

SI throughout. Free energy density F [J/m^3] of the unit magnetisation vector m:

    F(m) = K1 s^2 + K2 s^4 + K3 s^6 + K4 s^6 cos(6 phi)            (magnetocrystalline; s = sin(theta))
           - Ms (B . m) + (mu0 Ms^2 / 2) m.N.m                       (Zeeman + shape)

B denotes mu0*H_ext in TESLA (applied field, vacuum-equivalent); it is NOT the flux density inside the body.
K_i are magnetocrystalline constants ONLY; the shape contribution is carried separately by N.
Anisotropy fields (as mu0*H, tesla), shape NOT included (it enters through N):
  c-axis state (theta = 0):   mu0 H_A = 2 K1 / Ms   (K2, K3, K4 are O(theta^4) and do not contribute)
  basal state (theta = pi/2): effective fields of the two transverse modes
        mu0 H_theta = (-2 K1 - 4 K2 - 6 K3 + 6 |K4|) / Ms,      mu0 H_phi = 36 |K4| / Ms,
        f = (gamma/2pi) sqrt(H_theta H_phi) (+ Zeeman/shape terms). The basal minimum of K4 sin^6(theta) cos(6 phi)
        has cos(6 phi0) = -sgn(K4): phi0 = k*pi/3 for K4 < 0 and phi0 = pi/6 + k*pi/3 for K4 > 0.
Do NOT lump the shape term into K1 and also pass N (double counting).

NOTATION. In this code the basal-plane constant is K4 (coefficient of sin^6(theta) cos(6 phi)). Papers often call
the same term "K3 sin^6(theta) cos(6 phi)" (with K3 meaning the 6th-order basal constant, and then usually no
separate sin^6 theta term). Parameters taken from a paper MUST go through `map_literature_constants` with an
explicit mapping; there is no implicit import.

Cartesian extension: s^2 = 1 - mz^2 and s^6 cos6phi = Re (mx + i my)^6, valid on the unit sphere. The
Riemannian Hessian  P (Hess F) P - (m . grad F) P  does not depend on the extension.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from hexw import units as _u
from .demag import check_demag, sphere

MU0 = _u.MU0


@dataclass(frozen=True)
class MagnetModel:
    """Macrospin parameters supplied by the caller (no material defaults).
    Ms [A/m] mandatory; K1..K4 [J/m^3] (0 = term absent); N demag tensor in crystal frame (tr = 1, or the zero tensor = 'no own shape')."""
    Ms: float
    K1: float = 0.0
    K2: float = 0.0
    K3: float = 0.0
    K4: float = 0.0
    N: np.ndarray = field(default_factory=sphere)
    label: str = ""

    def __post_init__(self):
        if not (self.Ms > 0 and math.isfinite(self.Ms)):
            raise ValueError("Ms must be positive [A/m]")
        object.__setattr__(self, "N", check_demag(np.array(self.N, dtype=float)))

    @property
    def mu0Ms(self) -> float:
        """mu0*Ms [T]."""
        return MU0 * self.Ms

    @property
    def mu0_HA(self) -> float:
        """mu0*H_A = 2 K1 / Ms [T]: anisotropy field of the c-axis state (K2, K3, K4 do not contribute)."""
        return 2.0 * self.K1 / self.Ms

    @property
    def mu0_H_theta_basal(self) -> float:
        """mu0*H_theta of the basal-plane state [T] (shape not included):
        (-2 K1 - 4 K2 - 6 K3 + 6 |K4|)/Ms, valid at the basal minimum cos(6 phi0) = -sgn(K4)."""
        return (-2.0 * self.K1 - 4.0 * self.K2 - 6.0 * self.K3 + 6.0 * abs(self.K4)) / self.Ms

    @property
    def mu0_H_phi_basal(self) -> float:
        """mu0*H_phi of the basal-plane state [T]: 36 |K4| / Ms (at cos(6 phi0) = -sgn(K4))."""
        return 36.0 * abs(self.K4) / self.Ms

    def energy_scale(self, B_norm: float = 0.0) -> float:
        return max(abs(self.K1), abs(self.K2), abs(self.K3), abs(self.K4),
                   0.5 * MU0 * self.Ms ** 2 * float(np.abs(self.N).max()), self.Ms * B_norm, 1e-30)

    def with_N(self, N) -> "MagnetModel":
        return MagnetModel(self.Ms, self.K1, self.K2, self.K3, self.K4, np.array(N, dtype=float), self.label)


_ALLOWED_LIT_TARGETS = ("K1", "K2", "K3", "K4")


def map_literature_constants(constants: dict, mapping: dict) -> dict:
    """Translate anisotropy constants quoted in a paper into this code's K1..K4 through an EXPLICIT map.

    constants: {name used in the paper: value [J/m^3]};  mapping: {name used in the paper: 'K1'|'K2'|'K3'|'K4'}.
    Example: a paper writing  K1 sin^2 + K2 sin^4 + K3 sin^6 cos(6 phi)  maps {'K3': 'K4'} (so our K3 stays 0).
    Every constant must be mapped, targets must be unique, unknown targets are an error (no guessing); NaN/inf values are rejected."""
    out: dict = {}
    for name, value in constants.items():
        if name not in mapping:
            raise KeyError(f"constant {name!r} has no entry in the notation map (explicit mapping required)")
        tgt = mapping[name]
        if tgt not in _ALLOWED_LIT_TARGETS:
            raise ValueError(f"mapping target {tgt!r} must be one of {_ALLOWED_LIT_TARGETS}")
        if tgt in out:
            raise ValueError(f"two literature constants map to {tgt}")
        v = float(value)
        if not math.isfinite(v):
            raise ValueError(f"constant {name!r} is not finite ({value!r}); NaN/inf are rejected")
        out[tgt] = v
    return out


def unit(v):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if n == 0:
        raise ValueError("zero vector")
    return v / n


def field_vector(B_T: float, theta_H: float, phi_H: float = 0.0) -> np.ndarray:
    """B = mu0 H [T] from magnitude and crystal-frame angles (theta from c axis, phi azimuth), radians."""
    st = math.sin(theta_H)
    return B_T * np.array([st * math.cos(phi_H), st * math.sin(phi_H), math.cos(theta_H)])


def angles(m) -> tuple[float, float]:
    """(theta from +z in [0, pi], phi in [0, 2pi))."""
    m = unit(m)
    return math.acos(max(-1.0, min(1.0, m[2]))), math.atan2(m[1], m[0]) % (2 * math.pi)


def _aniso_terms(model: MagnetModel, m):
    """energy, gradient (3), Hessian (3x3) of the magnetocrystalline part in the Cartesian extension."""
    x, y, u = m
    s2 = 1.0 - u * u
    K1, K2, K3, K4 = model.K1, model.K2, model.K3, model.K4
    g = K1 * s2 + K2 * s2 ** 2 + K3 * s2 ** 3
    p = K1 + 2 * K2 * s2 + 3 * K3 * s2 ** 2
    dp = 2 * K2 + 6 * K3 * s2
    grad = np.array([0.0, 0.0, -2 * u * p])
    hess = np.zeros((3, 3))
    hess[2, 2] = -2 * p + 4 * u * u * dp
    E = g
    if K4 != 0.0:
        w = complex(x, y)
        w4, w5, w6 = w ** 4, w ** 5, w ** 6
        E += K4 * w6.real
        grad[0] += K4 * 6 * w5.real
        grad[1] += K4 * (-6 * w5.imag)
        hess[0, 0] += K4 * 30 * w4.real
        hess[0, 1] += K4 * (-30 * w4.imag)
        hess[1, 0] += K4 * (-30 * w4.imag)
        hess[1, 1] += K4 * (-30 * w4.real)
    return E, grad, hess


def energy(model: MagnetModel, B, m) -> float:
    m = unit(m)
    E = _aniso_terms(model, m)[0]
    E += -model.Ms * float(np.dot(B, m)) + 0.5 * MU0 * model.Ms ** 2 * float(m @ model.N @ m)
    return float(E)


def energy_many(model: MagnetModel, B, pts) -> np.ndarray:
    """Vectorised energy [J/m^3] for unit vectors pts (n,3)."""
    pts = np.asarray(pts, float)
    x, y, u = pts[:, 0], pts[:, 1], pts[:, 2]
    s2 = 1.0 - u * u
    E = model.K1 * s2 + model.K2 * s2 ** 2 + model.K3 * s2 ** 3
    if model.K4 != 0.0:
        E = E + model.K4 * ((x + 1j * y) ** 6).real
    E = E - model.Ms * (pts @ np.asarray(B, float)) + 0.5 * MU0 * model.Ms ** 2 * np.einsum("ni,ij,nj->n", pts, model.N, pts)
    return E


def grad_hess(model: MagnetModel, B, m):
    """(F, grad F, Hess F) in the 3D Cartesian extension at unit m."""
    m = unit(m)
    E, g, h = _aniso_terms(model, m)
    Ms = model.Ms
    E = E - Ms * float(np.dot(B, m)) + 0.5 * MU0 * Ms ** 2 * float(m @ model.N @ m)
    g = g - Ms * np.asarray(B, float) + MU0 * Ms ** 2 * (model.N @ m)
    h = h + MU0 * Ms ** 2 * model.N
    return float(E), g, h


def tangent_basis(m) -> np.ndarray:
    """3x2 matrix [e1 e2] with e1 x e2 = m (right-handed). Deterministic."""
    m = unit(m)
    a = np.array([1.0, 0.0, 0.0]) if abs(m[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = unit(a - m * np.dot(a, m))
    e2 = np.cross(m, e1)
    return np.stack([e1, e2], axis=1)


def sphere_derivs(model: MagnetModel, B, m):
    """Tangential gradient (2,), Riemannian Hessian (2x2) [J/m^3], basis E (3x2)."""
    m = unit(m)
    _, g, h = grad_hess(model, B, m)
    E = tangent_basis(m)
    gt = E.T @ g
    Hs = E.T @ h @ E - float(np.dot(m, g)) * np.eye(2)
    return gt, 0.5 * (Hs + Hs.T), E


# --------------------------------------------------------------------------------- vectorised (batch) kernels

def tangent_basis_batch(M: np.ndarray) -> np.ndarray:
    """(n,3,2) tangent bases [e1 e2] with e1 x e2 = m, same convention as `tangent_basis`."""
    M = np.asarray(M, float)
    use_x = np.abs(M[:, 0]) < 0.9
    a = np.zeros_like(M)
    a[use_x, 0] = 1.0
    a[~use_x, 1] = 1.0
    e1 = a - np.einsum("ni,ni->n", a, M)[:, None] * M
    e1 /= np.sqrt(np.einsum("ni,ni->n", e1, e1))[:, None]
    e2 = np.stack([M[:, 1] * e1[:, 2] - M[:, 2] * e1[:, 1],
                   M[:, 2] * e1[:, 0] - M[:, 0] * e1[:, 2],
                   M[:, 0] * e1[:, 1] - M[:, 1] * e1[:, 0]], axis=1)
    return np.stack([e1, e2], axis=2)


def grad_hess_batch(model: MagnetModel, B, M: np.ndarray):
    """Vectorised (F, grad F, Hess F) in the Cartesian extension at unit vectors M (n,3): (n,), (n,3), (n,3,3)."""
    M = np.asarray(M, float)
    n = len(M)
    x, y, u = M[:, 0], M[:, 1], M[:, 2]
    s2 = 1.0 - u * u
    K1, K2, K3, K4 = model.K1, model.K2, model.K3, model.K4
    E = K1 * s2 + K2 * s2 ** 2 + K3 * s2 ** 3
    p = K1 + 2 * K2 * s2 + 3 * K3 * s2 ** 2
    dp = 2 * K2 + 6 * K3 * s2
    g = np.zeros((n, 3))
    h = np.zeros((n, 3, 3))
    g[:, 2] = -2 * u * p
    h[:, 2, 2] = -2 * p + 4 * u * u * dp
    if K4 != 0.0:
        w = x + 1j * y
        w4 = w ** 4
        w5 = w4 * w
        E = E + K4 * (w5 * w).real
        g[:, 0] += K4 * 6 * w5.real
        g[:, 1] += K4 * (-6 * w5.imag)
        h[:, 0, 0] += K4 * 30 * w4.real
        h[:, 0, 1] += K4 * (-30 * w4.imag)
        h[:, 1, 0] += K4 * (-30 * w4.imag)
        h[:, 1, 1] += K4 * (-30 * w4.real)
    Ms = model.Ms
    Bv = np.asarray(B, float)
    NM = M @ model.N                      # N symmetric
    E = E - Ms * (M @ Bv) + 0.5 * MU0 * Ms ** 2 * np.einsum("ni,ni->n", M, NM)
    g = g - Ms * Bv[None, :] + MU0 * Ms ** 2 * NM
    h = h + MU0 * Ms ** 2 * model.N[None]
    return E, g, h


def sphere_derivs_batch(model: MagnetModel, B, M: np.ndarray, with_energy: bool = False):
    """Tangential gradient (n,2), Riemannian Hessian (n,2,2) [J/m^3], bases (n,3,2) for unit vectors M (n,3);
    with_energy=True additionally returns the energies (n,) as the last element."""
    E, g, h = grad_hess_batch(model, B, M)
    Eb = tangent_basis_batch(M)
    gt = np.einsum("nia,ni->na", Eb, g)
    Hs = np.einsum("nia,nij,njb->nab", Eb, h, Eb)
    rad = np.einsum("ni,ni->n", M, g)
    Hs[:, 0, 0] -= rad
    Hs[:, 1, 1] -= rad
    Hs = 0.5 * (Hs + np.swapaxes(Hs, 1, 2))
    if with_energy:
        return gt, Hs, Eb, E
    return gt, Hs, Eb


def eig2_batch(Hs: np.ndarray):
    """Closed-form eigen-decomposition of symmetric 2x2 matrices (n,2,2): (lam_min, lam_max, v_min (n,2))."""
    a, b, d = Hs[:, 0, 0], Hs[:, 0, 1], Hs[:, 1, 1]
    mean = 0.5 * (a + d)
    r = np.hypot(0.5 * (a - d), b)
    lmin, lmax = mean - r, mean + r
    # eigenvector of the larger eigenvalue at angle t/2 (tan t = 2b/(a-d)); the smaller one is perpendicular
    t = 0.5 * np.arctan2(2.0 * b, a - d)
    vmin = np.stack([-np.sin(t), np.cos(t)], axis=1)
    return lmin, lmax, vmin


def resolve_gamma(g: float | None = None, gamma_over_2pi_Hz_per_T: float | None = None,
                  gamma_rad_per_s_T: float | None = None) -> float:
    """Return gamma [rad/(s*T)]. EXACTLY ONE of the three must be given (no universal default).
    g -> gamma = g mu_B / hbar;  gamma/2pi [Hz/T] -> *2pi;  gamma [rad/(s T)] as is."""
    given = [v is not None for v in (g, gamma_over_2pi_Hz_per_T, gamma_rad_per_s_T)]
    if sum(given) != 1:
        raise ValueError("specify exactly one of g, gamma_over_2pi_Hz_per_T, gamma_rad_per_s_T")
    if g is not None:
        return float(_u.gamma_from_g(g).value)
    if gamma_over_2pi_Hz_per_T is not None:
        return float(gamma_over_2pi_Hz_per_T) * 2.0 * math.pi
    return float(gamma_rad_per_s_T)


@dataclass(frozen=True)
class Gyro:
    """Gyromagnetic ratio carrier: gamma [rad/(s*T)] and gamma/2pi [Hz/T] are DIFFERENT attributes."""
    gamma_rad_per_s_T: float

    @classmethod
    def from_g(cls, g): return cls(resolve_gamma(g=g))

    @classmethod
    def from_gamma_over_2pi(cls, v_Hz_per_T): return cls(resolve_gamma(gamma_over_2pi_Hz_per_T=v_Hz_per_T))

    @property
    def gamma_over_2pi_Hz_per_T(self) -> float:
        return self.gamma_rad_per_s_T / (2.0 * math.pi)


def as_gyro(g=None, gamma_over_2pi_Hz_per_T=None, gamma_rad_per_s_T=None, gyro: Gyro | None = None) -> Gyro:
    if gyro is not None:
        if any(v is not None for v in (g, gamma_over_2pi_Hz_per_T, gamma_rad_per_s_T)):
            raise ValueError("pass either gyro or g/gamma, not both")
        return gyro
    return Gyro(resolve_gamma(g, gamma_over_2pi_Hz_per_T, gamma_rad_per_s_T))
