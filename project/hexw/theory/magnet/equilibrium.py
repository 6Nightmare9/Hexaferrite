"""Equilibria of the macrospin energy: all local minima, stability, classification, branch following.

Search: relaxation from EVERY point of a Fibonacci grid on the sphere plus symmetry seeds (+-x, +-y, +-z, +-B)
(vectorised trust-region saddle-free Riemannian Newton, steps <= 0.25 rad, saddle escape along negative curvature)
-> clustering by angle -> control run on a doubled grid (different count => status 'possibly_incomplete').
Stability from the eigenvalues of the 2x2 Riemannian Hessian [J/m^3].

Classification of a state (angles from the c axis, tolerance `ang_tol`):
  axial            m || +-c
  planar           m in the basal plane (mz ~ 0)
  cone             intermediate theta and the field has no component perpendicular to c (incl. H = 0)
  tilted_by_field  intermediate theta with a transverse field component
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from ._fast import ScalarKernel
from .core import (MagnetModel, angles, eig2_batch, energy, energy_many, sphere_derivs, sphere_derivs_batch,
                   tangent_basis, unit)


def fibonacci_sphere(n: int) -> np.ndarray:
    """n quasi-uniform (equal-area) points on the unit sphere, shape (n, 3)."""
    i = np.arange(n) + 0.5
    z = 1.0 - 2.0 * i / n
    r = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    ph = i * math.pi * (3.0 - math.sqrt(5.0))
    return np.stack([r * np.cos(ph), r * np.sin(ph), z], axis=1)


@dataclass(frozen=True)
class Equilibrium:
    m: np.ndarray
    theta: float                  # rad, from +c
    phi: float                    # rad
    energy: float                 # J/m^3
    hess_eigs: tuple              # eigenvalues of the sphere Hessian [J/m^3], ascending
    stable: bool                  # lambda_min >= -tol (local minimum, possibly neutral)
    neutral: bool                 # some |lambda| <= tol: continuous degeneracy / marginal stability
    kind: str                     # axial | planar | cone | tilted_by_field
    azimuthal_degenerate: bool = False   # representative of a continuous ring (axial symmetry)


@dataclass
class EquilibriumSet:
    status: str                   # ok | possibly_incomplete | degenerate_isotropic | no_minimum_found
    minima: list = field(default_factory=list)       # all stable states, sorted by energy
    ground: Equilibrium | None = None                # lowest-energy state
    ground_degenerate: list = field(default_factory=list)   # states with energy equal to the ground state
    metastable: list = field(default_factory=list)   # [(Equilibrium, dE [J/m^3] > 0)]
    notes: list = field(default_factory=list)
    possibly_incomplete: bool = False                # minima count changed on the doubled grid (see notes)


def has_axial_symmetry(model: MagnetModel, B, tol: float = 1e-12) -> bool:
    """True if the energy is invariant under rotation about c: K4 = 0, N axisymmetric, B_perp = 0."""
    N = model.N
    ax = (model.K4 == 0.0 and abs(N[0, 0] - N[1, 1]) <= 1e-12 and abs(N[0, 1]) <= 1e-12
          and np.abs(N[2, :2]).max() <= 1e-12 and np.hypot(B[0], B[1]) <= tol * max(1.0, np.linalg.norm(B)))
    return bool(ax)


def classify(model: MagnetModel, B, m, ang_tol: float = 1e-3) -> str:
    th, _ = angles(m)
    Bperp = math.hypot(B[0], B[1])
    if th < ang_tol or math.pi - th < ang_tol:
        return "axial"
    if abs(th - math.pi / 2) < ang_tol:
        return "planar"
    if Bperp <= 1e-9 * max(np.linalg.norm(B), 1e-30) or np.linalg.norm(B) == 0.0:
        return "cone"
    return "tilted_by_field"


def relax_batch(model: MagnetModel, B, M0, tol_rel: float = 1e-13, maxit: int = 150, trust: float = 0.25,
                escape: bool = True, floor_rel: float = 1e-7, outer: int = 4, esc_amp: float = 1e-4):
    """Vectorised trust-region saddle-free Riemannian Newton descent from every row of M0 (n,3).

    Each step is computed in the eigenbasis of the 2x2 sphere Hessian with |lambda| (floored at floor_rel*scale) in
    the denominator, capped at `trust` rad (exponential map on the sphere) and accepted only if the energy does not
    increase (backtracking). With escape=True a converged saddle (lambda_min < -1e-9*scale) is left along the
    negative eigenvector and descended again (up to `outer` rounds).
    Returns (M (n,3), ok (n,)); ok marks |tangential gradient| < 1e-8*scale."""
    B = np.asarray(B, float)
    scale = model.energy_scale(float(np.linalg.norm(B)))
    M = np.array(M0, float, copy=True).reshape(-1, 3)
    M /= np.linalg.norm(M, axis=1)[:, None]
    tol_g = tol_rel * scale
    tol_neg = 1e-9 * scale
    floor = floor_rel * scale
    todo = np.arange(len(M))
    kern = None
    switch = 12
    for rnd in range(outer):
        act = todo
        for _ in range(maxit):
            if act.size == 0:
                break
            m = M[act]
            gt, Hs, Eb, F0 = sphere_derivs_batch(model, B, m, with_energy=True)
            keep = np.sqrt(np.einsum("na,na->n", gt, gt)) >= tol_g
            if not keep.all():
                act, m, gt, Hs, Eb, F0 = act[keep], m[keep], gt[keep], Hs[keep], Eb[keep], F0[keep]
                if act.size == 0:
                    break
            lmin, lmax, vmin = eig2_batch(Hs)
            vmax = np.stack([-vmin[:, 1], vmin[:, 0]], axis=1)
            c1 = np.einsum("na,na->n", vmin, gt) / np.maximum(np.abs(lmin), floor)
            c2 = np.einsum("na,na->n", vmax, gt) / np.maximum(np.abs(lmax), floor)
            st = -(c1[:, None] * vmin + c2[:, None] * vmax)
            ln = np.sqrt(np.einsum("na,na->n", st, st))
            st *= np.minimum(1.0, trust / np.maximum(ln, 1e-300))[:, None]
            t = np.ones(len(act))
            newm = m.copy()
            pend = np.arange(len(act))
            for _bt in range(30):
                step = st[pend] * t[pend, None]
                L = np.sqrt(np.einsum("na,na->n", step, step))
                dv = np.einsum("nia,na->ni", Eb[pend], step)
                mn = m[pend] * np.cos(L)[:, None] + dv * np.sinc(L / math.pi)[:, None]
                mn /= np.sqrt(np.einsum("ni,ni->n", mn, mn))[:, None]
                ok = energy_many(model, B, mn) <= F0[pend] + 1e-14 * scale
                newm[pend[ok]] = mn[ok]
                pend = pend[~ok]
                if pend.size == 0:
                    break
                t[pend] *= 0.5
            M[act] = newm
            stalled = np.zeros(len(act), bool)
            stalled[pend] = True
            act = act[~stalled]
            if act.size <= switch:
                break
        if 0 < act.size <= switch:           # stragglers: pure-Python Newton is cheaper than a numpy iteration
            kern = kern or ScalarKernel(model, B)
            for i in act:
                M[i] = kern.relax(M[i], maxit=maxit, outer=1, trust=trust, tol_rel=tol_rel, floor_rel=floor_rel,
                                  escape=False)[0]
        if not escape or rnd == outer - 1:
            break
        gt, Hs, Eb = sphere_derivs_batch(model, B, M[todo])
        lmin, _, vmin = eig2_batch(Hs)
        sad = lmin < -tol_neg
        if not sad.any():
            break
        if sad.sum() <= switch:              # few saddles: scalar kernel does the escape and the descent
            kern = kern or ScalarKernel(model, B)
            for i in todo[sad]:
                M[i] = kern.relax(M[i], maxit=maxit, outer=outer - rnd, trust=trust, tol_rel=tol_rel,
                                  floor_rel=floor_rel, esc_amp=esc_amp, escape=True)[0]
            break
        todo = todo[sad]
        Eb, vmin = Eb[sad], vmin[sad]
        dv = esc_amp * np.einsum("nia,na->ni", Eb, vmin)
        cand = []
        for sgn in (1.0, -1.0):
            c = M[todo] + sgn * dv
            cand.append(c / np.linalg.norm(c, axis=1)[:, None])
        e_p, e_m = energy_many(model, B, cand[0]), energy_many(model, B, cand[1])
        M[todo] = np.where((e_p <= e_m)[:, None], cand[0], cand[1])
    gt, _, _ = sphere_derivs_batch(model, B, M)
    return M, np.linalg.norm(gt, axis=1) < 1e-8 * scale


def relax(model: MagnetModel, B, m0, tol_rel: float = 1e-13, maxit: int = 200, escape: bool = True):
    """Relaxation from a single start m0 (thin wrapper over `relax_batch`). Returns unit m (a critical point;
    check the Hessian for its type). With escape=True a saddle is left along the negative eigenvector."""
    M, _ = relax_batch(model, B, np.asarray(m0, float)[None, :], tol_rel=tol_rel, maxit=maxit, escape=escape)
    return M[0]


def _make_eq(model, B, m, tol, ang_tol, azdeg=False) -> Equilibrium:
    gt, Hs, _ = sphere_derivs(model, B, m)
    lam = np.linalg.eigvalsh(Hs)
    th, ph = angles(m)
    return Equilibrium(m=unit(m), theta=th, phi=ph, energy=energy(model, B, m),
                       hess_eigs=(float(lam[0]), float(lam[1])), stable=bool(lam[0] >= -tol),
                       neutral=bool(abs(lam[0]) <= tol or abs(lam[1]) <= tol),
                       kind=classify(model, B, m, ang_tol), azimuthal_degenerate=azdeg)


def _cluster(M: np.ndarray, tol: float) -> list:
    """Representatives of clusters of points closer than tol (greedy, order-preserving)."""
    left = np.ones(len(M), bool)
    reps = []
    while left.any():
        i = int(np.argmax(left))
        reps.append(M[i])
        left &= np.linalg.norm(M - M[i], axis=1) >= tol
    return reps


def _grid_search(model, B, n, tol, sym, cluster_tol):
    """Minima found by relaxing from EVERY point of an n-point Fibonacci grid plus the symmetry seeds.
    Returns (list of unit vectors, number of seeds that did not converge)."""
    Bn = float(np.linalg.norm(B))
    sd = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
    if Bn > 0:
        sd += [list(B / Bn), list(-B / Bn)]
    M, ok = relax_batch(model, B, np.vstack([fibonacci_sphere(n), np.array(sd, float)]))
    if sym:                                          # canonical azimuth: one representative of each ring
        M = np.stack([np.hypot(M[:, 0], M[:, 1]), np.zeros(len(M)), M[:, 2]], axis=1)
    _, Hs, _ = sphere_derivs_batch(model, B, M)
    lmin, _, _ = eig2_batch(Hs)
    good = ok & (lmin >= -tol)
    return _cluster(M[good], cluster_tol), int((~ok).sum())


def find_equilibria(model: MagnetModel, B, n_grid: int = 600, ang_tol: float = 1e-3,
                    cluster_tol: float = 1e-4, stab_tol_rel: float = 1e-9, verify: bool = True) -> EquilibriumSet:
    """All local minima of F on the sphere for field vector B (= mu0 H, tesla, crystal frame).

    Seeds: EVERY point of an n_grid Fibonacci grid plus the symmetry seeds. With verify=True the search is repeated
    on a grid of 2*n_grid points; if the number of distinct minima differs, the result of the denser grid is returned
    with status 'possibly_incomplete' (flag `possibly_incomplete`, note): narrow basins may still be missing."""
    B = np.asarray(B, dtype=float)
    Bn = float(np.linalg.norm(B))
    scale = model.energy_scale(Bn)
    tol = stab_tol_rel * scale
    out = EquilibriumSet(status="ok")
    Eg = energy_many(model, B, fibonacci_sphere(max(n_grid, 50)))
    if np.ptp(Eg) <= 1e-12 * scale:
        out.status = "degenerate_isotropic"
        out.notes.append("energy is constant on the sphere: every direction is an equilibrium (Goldstone)")
        eq = _make_eq(model, B, np.array([0.0, 0.0, 1.0]), tol, ang_tol, azdeg=True)
        out.minima, out.ground, out.ground_degenerate = [eq], eq, [eq]
        return out
    sym = has_axial_symmetry(model, B)
    found, nbad = _grid_search(model, B, n_grid, tol, sym, cluster_tol)
    if verify:
        found2, nbad2 = _grid_search(model, B, 2 * n_grid, tol, sym, cluster_tol)
        if len(found2) != len(found):
            out.status = "possibly_incomplete"
            out.possibly_incomplete = True
            out.notes.append(f"possibly_incomplete: {len(found)} minima on {n_grid} grid points vs {len(found2)} on "
                             f"{2 * n_grid}; the denser-grid result is returned")
        found, nbad = found2, nbad2
    if nbad:
        out.notes.append(f"{nbad} seeds did not converge and were dropped")
    eqs = sorted((_make_eq(model, B, m, tol, ang_tol, azdeg=sym and math.hypot(m[0], m[1]) > 1e-9) for m in found),
                 key=lambda e: e.energy)
    if not eqs:
        out.status = "no_minimum_found"
        return out
    out.minima = eqs
    out.ground = eqs[0]
    etol = 1e-9 * scale
    out.ground_degenerate = [e for e in eqs if e.energy - eqs[0].energy <= etol]
    out.metastable = [(e, e.energy - eqs[0].energy) for e in eqs if e.energy - eqs[0].energy > etol]
    if sym:
        out.notes.append("axial symmetry: states are representatives of azimuthal rings (phi=0 shown)")
    return out


# --------------------------------------------------------------------------------- analytic cone

def cone_angle_analytic(K1: float, K2: float) -> float | None:
    """theta0 [rad] of the zero-field cone (K3 = K4 = 0, N-free): sin^2(theta0) = -K1/(2 K2) for K2 > 0,
    K1 < 0 and -K1/(2K2) < 1. None otherwise."""
    if K2 > 0 and K1 < 0:
        s2 = -K1 / (2.0 * K2)
        if s2 < 1.0:
            return math.asin(math.sqrt(s2))
    return None


# --------------------------------------------------------------------------------- branch following

@dataclass
class BranchPoint:
    B_T: float
    state: Equilibrium
    jumped: bool
    angle_step: float


def trace_branch(model: MagnetModel, direction, B_values, m_init, jump_angle: float = 0.5,
                 ang_tol: float = 1e-3) -> list[BranchPoint]:
    """Follow the magnetisation branch for a field sweep B_values (scalar tesla along unit `direction`,
    crystal frame), starting from state m_init at B_values[0] (relaxed there). At each step the previous
    state is relaxed (no global search), so the branch is kept until it ceases to be a local minimum, at
    which point the descent leaves it: `jumped` is set when the direction changes by more than jump_angle
    (an irreversible switching event) -- the branch is then continued on the new minimum."""
    d = unit(direction)
    out: list[BranchPoint] = []
    m_prev = None
    for k, b in enumerate(B_values):
        Bv = b * d
        m_start = unit(m_init) if m_prev is None else m_prev
        # a deterministic tiny transverse nudge so that exact-symmetry saddles are not mistaken for branches
        e = tangent_basis(m_start)[:, 0]
        m_new = relax(model, Bv, unit(m_start + 1e-9 * e), escape=True)
        step = 0.0 if m_prev is None else math.acos(max(-1.0, min(1.0, float(np.dot(m_new, m_prev)))))
        scale = model.energy_scale(abs(b))
        eq = _make_eq(model, Bv, m_new, 1e-9 * scale, ang_tol)
        out.append(BranchPoint(float(b), eq, bool(step > jump_angle), step))
        m_prev = m_new
    return out


def hysteresis_loop(model: MagnetModel, direction, B_max: float, n: int, m_init=None, jump_angle: float = 0.5):
    """Field loop +B_max -> -B_max -> +B_max along `direction` from the state `m_init` (default: along +direction).
    Returns (B_values, branch points, projection of m on direction)."""
    d = unit(direction)
    Bs = np.concatenate([np.linspace(B_max, -B_max, n), np.linspace(-B_max, B_max, n)[1:]])
    pts = trace_branch(model, d, Bs, d if m_init is None else m_init, jump_angle=jump_angle)
    proj = np.array([float(np.dot(p.state.m, d)) for p in pts])
    return Bs, pts, proj
