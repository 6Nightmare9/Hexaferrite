"""Pure-Python scalar kernel for the sequential steps of branch tracking (one state at a time).

For a single state numpy has a call overhead of ~10 us per operation and ~40 operations per Newton iteration;
plain `math` on floats is ~10x faster. This module duplicates the formulas of `core.grad_hess_batch` /
`equilibrium.relax_batch` (tests check agreement to rounding); it is used ONLY where one state is followed along a
field sweep. Vectorised global searches use the numpy kernels.
"""
from __future__ import annotations

import math

from .core import MU0, MagnetModel


class ScalarKernel:
    """Energy, sphere derivatives and a trust-region saddle-free Newton relaxation for ONE state at a fixed field."""

    def __init__(self, model: MagnetModel, Bv):
        self.K1, self.K2, self.K3, self.K4 = model.K1, model.K2, model.K3, model.K4
        Ms = model.Ms
        c = MU0 * Ms * Ms
        N = model.N
        self.n00, self.n01, self.n02 = c * N[0, 0], c * N[0, 1], c * N[0, 2]
        self.n11, self.n12, self.n22 = c * N[1, 1], c * N[1, 2], c * N[2, 2]
        self.zx, self.zy, self.zz = Ms * float(Bv[0]), Ms * float(Bv[1]), Ms * float(Bv[2])
        self.scale = model.energy_scale(math.sqrt(float(Bv[0]) ** 2 + float(Bv[1]) ** 2 + float(Bv[2]) ** 2))

    def energy(self, x, y, u):
        s2 = 1.0 - u * u
        E = self.K1 * s2 + self.K2 * s2 * s2 + self.K3 * s2 * s2 * s2
        if self.K4:
            w = complex(x, y)
            w2 = w * w
            E += self.K4 * (w2 * w2 * w2).real
        E -= self.zx * x + self.zy * y + self.zz * u
        E += 0.5 * (self.n00 * x * x + self.n11 * y * y + self.n22 * u * u
                    + 2.0 * (self.n01 * x * y + self.n02 * x * u + self.n12 * y * u))
        return E

    def derivs(self, x, y, u):
        """(g0, g1, a, b, d, e1, e2): tangential gradient and symmetric 2x2 Riemannian Hessian in the basis (e1, e2)."""
        K4 = self.K4
        s2 = 1.0 - u * u
        p = self.K1 + 2.0 * self.K2 * s2 + 3.0 * self.K3 * s2 * s2
        dp = 2.0 * self.K2 + 6.0 * self.K3 * s2
        gx = gy = hxx = hxy = hyy = 0.0
        gz = -2.0 * u * p
        hzz = -2.0 * p + 4.0 * u * u * dp
        if K4:
            w = complex(x, y)
            w2 = w * w
            w4 = w2 * w2
            w5 = w4 * w
            gx = 6.0 * K4 * w5.real
            gy = -6.0 * K4 * w5.imag
            hxx = 30.0 * K4 * w4.real
            hxy = -30.0 * K4 * w4.imag
            hyy = -hxx
        gx += -self.zx + self.n00 * x + self.n01 * y + self.n02 * u
        gy += -self.zy + self.n01 * x + self.n11 * y + self.n12 * u
        gz += -self.zz + self.n02 * x + self.n12 * y + self.n22 * u
        Hxx, Hxy, Hxz = hxx + self.n00, hxy + self.n01, self.n02
        Hyy, Hyz, Hzz = hyy + self.n11, self.n12, hzz + self.n22
        rad = gx * x + gy * y + gz * u
        if abs(x) < 0.9:
            e1x, e1y, e1z = 1.0 - x * x, -x * y, -x * u
        else:
            e1x, e1y, e1z = -y * x, 1.0 - y * y, -y * u
        nrm = math.sqrt(e1x * e1x + e1y * e1y + e1z * e1z)
        e1x, e1y, e1z = e1x / nrm, e1y / nrm, e1z / nrm
        e2x, e2y, e2z = y * e1z - u * e1y, u * e1x - x * e1z, x * e1y - y * e1x
        g0 = e1x * gx + e1y * gy + e1z * gz
        g1 = e2x * gx + e2y * gy + e2z * gz
        h1x = Hxx * e1x + Hxy * e1y + Hxz * e1z
        h1y = Hxy * e1x + Hyy * e1y + Hyz * e1z
        h1z = Hxz * e1x + Hyz * e1y + Hzz * e1z
        h2x = Hxx * e2x + Hxy * e2y + Hxz * e2z
        h2y = Hxy * e2x + Hyy * e2y + Hyz * e2z
        h2z = Hxz * e2x + Hyz * e2y + Hzz * e2z
        a = e1x * h1x + e1y * h1y + e1z * h1z - rad
        d = e2x * h2x + e2y * h2y + e2z * h2z - rad
        b = 0.5 * (e1x * h2x + e1y * h2y + e1z * h2z + e2x * h1x + e2y * h1y + e2z * h1z)
        return g0, g1, a, b, d, (e1x, e1y, e1z), (e2x, e2y, e2z)

    @staticmethod
    def eig(a, b, d):
        """(lmin, lmax, (v0, v1) eigenvector of lmin) for the symmetric matrix [[a, b], [b, d]]."""
        mean = 0.5 * (a + d)
        r = math.hypot(0.5 * (a - d), b)
        t = 0.5 * math.atan2(2.0 * b, a - d)
        return mean - r, mean + r, (-math.sin(t), math.cos(t))

    def relax(self, m0, maxit: int = 25, outer: int = 2, trust: float = 0.25, tol_rel: float = 1e-13,
              floor_rel: float = 1e-7, esc_amp: float = 1e-4, escape: bool = True):
        """Returns (m (3-tuple), ok, lmin, lmax). Same algorithm as `relax_batch` for one start."""
        scale = self.scale
        tol_g, tol_neg, floor = tol_rel * scale, 1e-9 * scale, floor_rel * scale
        x, y, u = m0
        nrm = math.sqrt(x * x + y * y + u * u)
        x, y, u = x / nrm, y / nrm, u / nrm
        for rnd in range(outer):
            for _ in range(maxit):
                g0, g1, a, b, d, e1, e2 = self.derivs(x, y, u)
                if math.hypot(g0, g1) < tol_g:
                    break
                lmin, lmax, (v0, v1) = self.eig(a, b, d)
                c1 = (v0 * g0 + v1 * g1) / max(abs(lmin), floor)
                c2 = (-v1 * g0 + v0 * g1) / max(abs(lmax), floor)
                s0 = -(c1 * v0 - c2 * v1)
                s1 = -(c1 * v1 + c2 * v0)
                ln = math.hypot(s0, s1)
                if ln > trust:
                    s0, s1 = s0 * trust / ln, s1 * trust / ln
                F0 = self.energy(x, y, u)
                moved = False
                for _bt in range(30):
                    L = math.hypot(s0, s1)
                    cL = math.cos(L)
                    sL = math.sin(L) / L if L > 1e-300 else 1.0
                    dx = s0 * e1[0] + s1 * e2[0]
                    dy = s0 * e1[1] + s1 * e2[1]
                    dz = s0 * e1[2] + s1 * e2[2]
                    nx, ny, nu = x * cL + dx * sL, y * cL + dy * sL, u * cL + dz * sL
                    nn = math.sqrt(nx * nx + ny * ny + nu * nu)
                    nx, ny, nu = nx / nn, ny / nn, nu / nn
                    if self.energy(nx, ny, nu) <= F0 + 1e-14 * scale:
                        x, y, u = nx, ny, nu
                        moved = True
                        break
                    s0, s1 = 0.5 * s0, 0.5 * s1
                if not moved:
                    break
            g0, g1, a, b, d, e1, e2 = self.derivs(x, y, u)
            lmin, lmax, (v0, v1) = self.eig(a, b, d)
            if not escape or rnd == outer - 1 or lmin >= -tol_neg:
                break
            best = None
            for sg in (1.0, -1.0):
                cx = x + sg * esc_amp * (v0 * e1[0] + v1 * e2[0])
                cy = y + sg * esc_amp * (v0 * e1[1] + v1 * e2[1])
                cz = u + sg * esc_amp * (v0 * e1[2] + v1 * e2[2])
                cn = math.sqrt(cx * cx + cy * cy + cz * cz)
                cx, cy, cz = cx / cn, cy / cn, cz / cn
                ce = self.energy(cx, cy, cz)
                if best is None or ce < best[0]:
                    best = (ce, cx, cy, cz)
            _, x, y, u = best
        g0, g1, a, b, d, e1, e2 = self.derivs(x, y, u)
        lmin, lmax, _ = self.eig(a, b, d)
        return (x, y, u), math.hypot(g0, g1) < 1e-8 * scale, lmin, lmax
