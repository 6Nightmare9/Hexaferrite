"""Demagnetising tensors N (tr N = 1) in the CRYSTAL frame (z = hexagonal c axis).

All functions are pure; nothing is read from disk.
  * sphere()                    N = I/3
  * thin_disk_perp_c()          disk plane perpendicular to c  -> N = diag(0, 0, 1)
  * thin_disk_parallel_c()      disk plane contains c, normal along crystal x -> N = diag(1, 0, 0)
  * ellipsoid(a, b, c)          semi-axes a || x, b || y, c || z (Osborn 1945, numerical quadrature of
                                N_i = (abc/2) int_0^inf ds / ((s + a_i^2) sqrt((s+a^2)(s+b^2)(s+c^2))))
"""
from __future__ import annotations

import numpy as np
from scipy import integrate


def sphere() -> np.ndarray:
    return np.eye(3) / 3.0


def thin_disk_perp_c() -> np.ndarray:
    return np.diag([0.0, 0.0, 1.0])


def thin_disk_parallel_c() -> np.ndarray:
    return np.diag([1.0, 0.0, 0.0])


def ellipsoid(a: float, b: float, c: float) -> np.ndarray:
    """Demagnetising tensor of a general ellipsoid with semi-axes (a, b, c) along crystal (x, y, z).
    Only ratios matter. Returns diag(Nx, Ny, Nz), trace = 1 (checked)."""
    ax = np.array([a, b, c], dtype=float)
    if np.any(ax <= 0) or not np.all(np.isfinite(ax)):
        raise ValueError("semi-axes must be positive and finite")
    ax = ax / ax.max()
    a2 = ax ** 2

    def n_i(i: int) -> float:
        # substitution s = t/(1-t) maps [0, inf) onto [0, 1)
        def f(t):
            s = t / (1.0 - t)
            r = np.sqrt(np.prod(s + a2))
            return 0.5 * np.prod(ax) / ((s + a2[i]) * r) / (1.0 - t) ** 2
        val, _ = integrate.quad(f, 0.0, 1.0, epsabs=1e-13, epsrel=1e-12, limit=400)
        return val

    n = np.array([n_i(0), n_i(1), n_i(2)])
    if abs(n.sum() - 1.0) > 1e-6:
        raise ArithmeticError(f"demag trace {n.sum()} != 1")
    return np.diag(n / n.sum())


def no_shape() -> np.ndarray:
    """N = 0: the object has no shape of its own (the field acting on it is already the INTERNAL field).
    The only tensor allowed with trace != 1; used for grains in the macroscopic mean-field mode."""
    return np.zeros((3, 3))


def check_demag(N: np.ndarray, tol: float = 1e-9) -> np.ndarray:
    N = np.asarray(N, dtype=float)
    if N.shape != (3, 3):
        raise ValueError("N must be 3x3")
    if not np.any(N):
        return N                   # explicit "no own shape" (internal-field mode)
    if not np.allclose(N, N.T, atol=tol):
        raise ValueError("N must be symmetric")
    if abs(np.trace(N) - 1.0) > 1e-6:
        raise ValueError(f"tr N must be 1, got {np.trace(N)}")
    if np.linalg.eigvalsh(N).min() < -1e-9:
        raise ValueError("N must be positive semidefinite")
    return N


def rotate_tensor(N: np.ndarray, R: np.ndarray) -> np.ndarray:
    """N' = R N R^T (R maps the frame of N into the new frame)."""
    return R @ N @ R.T
