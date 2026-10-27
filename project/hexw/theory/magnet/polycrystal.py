"""Ideal single crystal and non-interacting grain ensemble.

Result levels (always attached to the result, never dropped):
    'ideal_crystal'                 one crystal, given orientation, no averaging
    'grain_ensemble_noninteracting' orientational average of independent single-domain grains.

Frames: crystal frame (z = c) inside MagnetModel; lab frame for fields and N_lab. A grain with rotation R
(columns = crystal axes expressed in the lab frame, v_lab = R v_crystal) sees
    B_c = R^T B_lab,    N_c = R^T N_lab R,    chi_lab = R chi_c R^T.

DEMAGNETISING MODES of `polycrystal_chi` (argument `demag`):
  'macro_meanfield' (default)  grains have NO shape of their own (N_grain = 0, or given explicitly); the sample shape
        acts through the macroscopic field  B_int = B_app - mu0 Ms N_sample <m>  (mean-field, solved self-consistently
        with a convergence check). The dynamic response is averaged as the INTERNAL susceptibility chi_int (response
        to the macroscopic field inside the sample); the apparent chi_ext = chi_int (I + N_sample chi_int)^-1 of the
        sample is returned as well.
  'per_grain'  every grain is treated as a sample-shaped ellipsoid in the applied field. Physically WRONG for a bulk
        ceramic (each grain would carry the whole sample shape): the result is flagged
        `physically_invalid_for_bulk=True` and a warning is issued. Kept for isolated particles / comparisons.

Not modelled (explicit stubs returning a status, never numbers): intergranular coupling, porosity, effective-medium
mixing (Bruggeman / Maxwell-Garnett), Schloemann porosity linewidth.
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field

import numpy as np

from ._fast import ScalarKernel
from .core import MU0, Gyro, MagnetModel, unit
from .demag import no_shape
from .equilibrium import find_equilibria, fibonacci_sphere, relax
from .susceptibility import susceptibility

LEVEL_IDEAL = "ideal_crystal"
LEVEL_ENSEMBLE = "grain_ensemble_noninteracting"
DEMAG_MODES = ("macro_meanfield", "per_grain")
# Kronecker-sequence increment (in turns) for the crystal azimuth about c; irrational and unrelated to the
# Fibonacci longitude increment (1 - 1/phi = 0.381966 turns) so that (c direction, azimuth) fill SO(3) without correlation.
AZIMUTH_TURNS = 0.5698402909980532      # 1/rho^2, rho = plastic number (R2 sequence, second component)
MIN_AZIMUTHS_FOR_K4 = 24


def rotation_from_c_axis(c_dir, azimuth: float = 0.0) -> np.ndarray:
    """Proper rotation R with R[:,2] = c_dir (lab) and rotation `azimuth` about c."""
    z = unit(c_dir)
    a = np.array([1.0, 0.0, 0.0]) if abs(z[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    x = unit(a - z * np.dot(a, z))
    y = np.cross(z, x)
    ca, sa = math.cos(azimuth), math.sin(azimuth)
    return np.stack([ca * x + sa * y, -sa * x + ca * y, z], axis=1)


def rotation_matrix(axis, angle: float) -> np.ndarray:
    a = unit(axis)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + math.sin(angle) * K + (1 - math.cos(angle)) * (K @ K)


@dataclass(frozen=True)
class Texture:
    """Texture of the c axes. r = 1: isotropic. March-Dollase: w(psi) ~ (r^2 cos^2 psi + sin^2 psi / r)^(-3/2),
    psi = angle between the c axis and `axis` (lab); r < 1 platelets-like alignment of c along `axis`
    (r is the March parameter, an INPUT)."""
    r: float = 1.0
    axis: tuple = (0.0, 0.0, 1.0)

    def weights(self, c_dirs: np.ndarray) -> np.ndarray:
        if self.r <= 0:
            raise ValueError("March parameter r must be > 0")
        if self.r == 1.0:
            return np.ones(len(c_dirs)) / len(c_dirs)
        cp = c_dirs @ unit(self.axis)
        w = (self.r ** 2 * cp ** 2 + (1 - cp ** 2) / self.r) ** -1.5
        return w / w.sum()


def quasi_azimuths(n: int) -> np.ndarray:
    """Crystal azimuths about c [rad], phi_i = 2 pi frac(0.5 + i a) with an irrational increment a unrelated to the
    Fibonacci longitudes of the c directions (a quasi-random sequence independent of the polar angle)."""
    return 2.0 * math.pi * np.mod(0.5 + np.arange(n) * AZIMUTH_TURNS, 1.0)


def orientation_grid(n_c: int, n_phi: int = 1, global_rotation: np.ndarray | None = None, azimuth: str = "quasi"):
    """Quasi-uniform orientations: n_c Fibonacci c-axis directions x n_phi azimuths about c.
    azimuth='quasi' (default): every c direction gets its own quasi-random azimuth offset (`quasi_azimuths`) and the
    n_phi azimuths are offset copies spaced by 2 pi / n_phi; azimuth='zero': all offsets 0 (the old behaviour; exact
    only when the azimuth is a symmetry: K4 = 0 and axisymmetric lab field / N).
    Returns (list of R, c_dirs (n_c,3)); the azimuths used are available through `orientation_azimuths`."""
    Rs, cd, _ = _orientations(n_c, n_phi, global_rotation, azimuth)
    return Rs, cd


def _orientations(n_c, n_phi, global_rotation, azimuth):
    pts = fibonacci_sphere(n_c)
    if global_rotation is not None:
        pts = pts @ global_rotation.T
    if azimuth == "quasi":
        off = quasi_azimuths(n_c)
    elif azimuth == "zero":
        off = np.zeros(n_c)
    else:
        raise ValueError("azimuth must be 'quasi' or 'zero'")
    az = np.array([off[i] + 2.0 * math.pi * k / n_phi for i in range(n_c) for k in range(n_phi)])
    Rs = [rotation_from_c_axis(pts[i // n_phi], az[i]) for i in range(n_c * n_phi)]
    return Rs, pts, az


def distinct_azimuths_mod(az: np.ndarray, period: float = math.pi / 3.0, decimals: int = 6) -> int:
    """Number of distinct crystal azimuths modulo `period` (pi/3 for the cos 6phi term)."""
    return int(len(np.unique(np.round(np.mod(az, period) / period, decimals))))


@dataclass
class CrystalResult:
    level: str
    freqs_Hz: np.ndarray
    chi_lab: np.ndarray
    kind: str
    notes: list = field(default_factory=list)


@dataclass
class EnsembleResult:
    level: str
    freqs_Hz: np.ndarray
    chi_lab: np.ndarray           # (nf, 3, 3) weighted mean over grains, lab frame; meaning given by `kind`
    kind: str                     # 'internal' (material, macro_meanfield default) or 'external' (shape-dependent)
    n_grains: int
    n_failed: int
    state_selection: str
    texture: Texture
    notes: list = field(default_factory=list)
    demag: str = "per_grain"
    physically_invalid_for_bulk: bool = False
    converged: bool = True        # macro_meanfield: self-consistent internal field found
    n_iter: int = 0
    B_int_lab: np.ndarray | None = None          # macroscopic internal field mu0 H_int [T] (macro_meanfield)
    M_avg_A_per_m: np.ndarray | None = None      # mean magnetisation [A/m] (macro_meanfield)
    chi_int_lab: np.ndarray | None = None        # <chi_int>
    chi_ext_lab: np.ndarray | None = None        # chi_int (I + N_sample chi_int)^-1 for the sample shape
    N_sample: np.ndarray | None = None
    quadrature: str = "quasi_random"
    n_distinct_azimuths_mod_pi3: int | None = None


def to_crystal_frame(R: np.ndarray, B_lab, N_lab):
    return R.T @ np.asarray(B_lab, float), R.T @ np.asarray(N_lab, float) @ R


def _relax_fast(gm: MagnetModel, Bc, m):
    kern = ScalarKernel(gm, Bc)
    return np.array(kern.relax(tuple(float(c) for c in m), maxit=60, outer=3)[0])


def _initial_states(gm, Bc, R, B_lab, m_init_lab, state, n_grid):
    """List of equilibrium directions (crystal frame) of one grain: global ground state(s) or field-aligned."""
    if state == "global":
        eqs = find_equilibria(gm, Bc, n_grid=n_grid, verify=False)
        return [e.m for e in eqs.ground_degenerate] if eqs.ground else []
    if state == "field_aligned":
        ref = np.asarray(B_lab, float) if np.linalg.norm(B_lab) > 0 else np.asarray(m_init_lab, float)
        return [relax(gm, Bc, R.T @ unit(ref))]
    raise ValueError("state must be 'global' or 'field_aligned'")


def _grain_chi_from_states(gm, R, Bc, ms, freqs, alpha, gyro, kind):
    chis = [susceptibility(gm, Bc, m, freqs, alpha=alpha, gyro=gyro, kind=kind) for m in ms]
    chic = sum(chis) / len(chis)
    return np.einsum("ia,fab,jb->fij", R, chic, R)


def _grain_chi(model, R, B_lab, N_lab, freqs, alpha, gyro, state, m_init_lab, kind, n_grid):
    Bc, Nc = to_crystal_frame(R, B_lab, N_lab)
    gm = model.with_N(Nc)
    ms = _initial_states(gm, Bc, R, B_lab, m_init_lab, state, n_grid)
    if not ms:
        return None
    return _grain_chi_from_states(gm, R, Bc, ms, freqs, alpha, gyro, kind)


def crystal_chi(model: MagnetModel, R: np.ndarray, B_lab, freqs_Hz, *, alpha: float, gyro: Gyro,
                N_lab=None, state: str = "global", m_init_lab=(0, 0, 1), kind: str = "external",
                n_grid: int = 300) -> CrystalResult:
    """chi of ONE crystal with orientation R (level 'ideal_crystal'); the crystal has the shape N_lab (default model.N)."""
    N_lab = model.N if N_lab is None else np.asarray(N_lab, float)
    chi = _grain_chi(model, R, B_lab, N_lab, np.atleast_1d(freqs_Hz), alpha, gyro, state, m_init_lab, kind, n_grid)
    if chi is None:
        raise RuntimeError("no stable equilibrium found")
    return CrystalResult(LEVEL_IDEAL, np.atleast_1d(freqs_Hz), chi, kind)


# --------------------------------------------------------------------------------- axisymmetric quadrature

def _axis_for_quadrature(model, B_lab, N_s, N_g, texture):
    """Common symmetry axis a (unit vector) if field, N_sample, N_grain and texture are axisymmetric about a: the
    rotation of the c axis about a is then a symmetry, so only the polar angle psi (and, for K4 != 0, the crystal
    azimuth beta about c) matter. Else None."""
    Bn = np.linalg.norm(B_lab)
    if Bn > 0:
        a = np.asarray(B_lab, float) / Bn
    elif texture.r != 1.0:
        a = unit(texture.axis)
    else:
        a = np.array([0.0, 0.0, 1.0])
    if texture.r != 1.0 and abs(abs(np.dot(unit(texture.axis), a)) - 1.0) > 1e-9:
        return None
    P = np.eye(3) - np.outer(a, a)
    for N in (N_s, N_g):
        Na = N @ a
        if np.linalg.norm(Na - np.dot(a, Na) * a) > 1e-9:
            return None
        Np = P @ N @ P
        mu = np.trace(Np) / 2.0
        if np.abs(Np - mu * P).max() > 1e-9:
            return None
    return a


def _frame_with_axis(a):
    R = rotation_from_c_axis(a)          # columns (x', y', a)
    return R


def _gl_orientations(n_c, a, texture, n_beta):
    """Gauss-Legendre nodes in u = cos(psi) (psi from the axis a) x n_beta equally spaced crystal azimuths beta in one
    period pi/3 of the cos(6 phi) term (n_beta = 1, beta = 0 when K4 = 0). Weights include the texture."""
    u, w = np.polynomial.legendre.leggauss(n_c)
    if texture.r != 1.0:
        w = w * (texture.r ** 2 * u ** 2 + (1 - u ** 2) / texture.r) ** -1.5
    w = w / w.sum()
    Fa = _frame_with_axis(a)
    cd = np.stack([np.sqrt(np.maximum(0.0, 1 - u ** 2)), np.zeros(n_c), u], axis=1) @ Fa.T
    betas = (np.arange(n_beta) + 0.5) * (math.pi / 3.0) / n_beta if n_beta > 1 else np.zeros(1)
    Rs = [rotation_from_c_axis(cd[i], betas[k]) for i in range(n_c) for k in range(n_beta)]
    wts = np.repeat(w, n_beta) / n_beta
    return Rs, cd, wts, Fa, np.tile(betas, n_c)


def _average_about_axis(chi: np.ndarray, Fa: np.ndarray) -> np.ndarray:
    """Average chi (nf,3,3, lab) over rotations about the axis (third column of Fa)."""
    T = np.einsum("ai,fab,bj->fij", Fa, chi, Fa)        # axis along z'
    out = np.zeros_like(T)
    out[:, 0, 0] = out[:, 1, 1] = 0.5 * (T[:, 0, 0] + T[:, 1, 1])
    out[:, 0, 1] = 0.5 * (T[:, 0, 1] - T[:, 1, 0])
    out[:, 1, 0] = -out[:, 0, 1]
    out[:, 2, 2] = T[:, 2, 2]
    return np.einsum("ia,fab,jb->fij", Fa, out, Fa)


# --------------------------------------------------------------------------------- the ensemble

def polycrystal_chi(model: MagnetModel, B_lab, freqs_Hz, *, alpha: float, gyro: Gyro, N_lab=None,
                    N_sample=None, N_grain=None, demag: str = "macro_meanfield",
                    n_c: int = 200, n_phi: int = 1, texture: Texture = Texture(),
                    state: str = "global", m_init_lab=(0, 0, 1), kind: str | None = None,
                    global_rotation: np.ndarray | None = None, n_grid: int = 200, azimuth: str = "quasi",
                    quadrature: str = "auto", n_beta: int | None = None, macro_tol_T: float = 1e-7, macro_maxit: int = 80,
                    min_azimuths: int = MIN_AZIMUTHS_FOR_K4) -> EnsembleResult:
    """Orientation average of non-interacting grains (level 'grain_ensemble_noninteracting').

    demag='macro_meanfield' (default): grains carry N_grain (default 0 = no own shape); the sample shape N_sample
        (default: N_lab, else model.N, interpreted in the LAB frame) enters through the self-consistent macroscopic
        field B_int = B_app - mu0 Ms N_sample <m>. chi_lab is the averaged INTERNAL susceptibility (kind='internal',
        material); kind='external' returns chi_ext of the sample shape instead. Always available in the result:
        chi_int_lab, chi_ext_lab, B_int_lab, M_avg_A_per_m, converged, n_iter.
    demag='per_grain': every grain gets N_grain (default N_sample) as its own shape; default kind='external'; flagged
        physically_invalid_for_bulk=True with a warning.
    state='global': each grain sits in its global energy minimum (ties averaged); state='field_aligned': each grain
        relaxes from the field direction (or m_init_lab at B = 0) -- a history-like (saturated/remanent) state, may be
        metastable.
    Orientation sampling: quadrature='auto' uses a Gauss-Legendre rule in cos(psi) (n_c nodes) when the field, N_sample,
        N_grain and the texture are axisymmetric about a common axis (the rotation of c about that axis is a symmetry,
        averaged exactly); for K4 != 0 it is a product rule with n_beta (default 12) equally spaced crystal azimuths in
        one period pi/3. Otherwise n_c Fibonacci c directions with INDEPENDENT quasi-random crystal azimuths
        (`orientation_grid`); with K4 != 0 fewer than `min_azimuths` distinct azimuths modulo pi/3 triggers a warning
        and a note.
    Small damping (alpha < ~0.02) gives a comb of resonances for a finite number of grains/nodes: use a dense
    sampling (n_c >= several hundred) and check convergence against n_c."""
    if demag not in DEMAG_MODES:
        raise ValueError(f"demag must be one of {DEMAG_MODES}")
    freqs = np.atleast_1d(np.asarray(freqs_Hz, float))
    B_lab = np.asarray(B_lab, float)
    N_s = np.asarray(N_sample if N_sample is not None else (N_lab if N_lab is not None else model.N), float)
    notes: list[str] = []
    if demag == "per_grain":
        N_g = np.asarray(N_grain if N_grain is not None else N_s, float)
        kind = "external" if kind is None else kind
        msg = ("demag='per_grain': every grain is given the sample shape -- physically invalid for a bulk sample "
               "(use demag='macro_meanfield')")
        warnings.warn(msg, UserWarning, stacklevel=2)
        notes.append(msg)
    else:
        N_g = np.asarray(N_grain if N_grain is not None else no_shape(), float)
        kind = "internal" if kind is None else kind
        if kind not in ("internal", "external"):
            raise ValueError("kind must be 'internal' or 'external'")

    # ---- orientation set
    axis = _axis_for_quadrature(model, B_lab, N_s, N_g, texture) if quadrature in ("auto", "gauss_legendre") else None
    if quadrature == "gauss_legendre" and axis is None:
        raise ValueError("gauss_legendre requires a common symmetry axis of field, N_sample, N_grain and texture")
    if quadrature not in ("auto", "gauss_legendre", "quasi"):
        raise ValueError("quadrature must be auto|gauss_legendre|quasi")
    if axis is not None:
        nb = (12 if n_beta is None else int(n_beta)) if model.K4 != 0.0 else 1
        Rs, cdirs, w_all, Fa, az = _gl_orientations(n_c, axis, texture, nb)
        quad = "gauss_legendre_cos_psi" + ("_x_beta" if nb > 1 else "")
    else:
        Rs, cdirs, az = _orientations(n_c, n_phi, global_rotation, azimuth)
        w_c = texture.weights(cdirs)
        w_all = np.repeat(w_c, n_phi) / n_phi
        Fa = None
        quad = "quasi_random"
    n_az = None
    if az is not None:
        n_az = distinct_azimuths_mod(az)
        # deterministic quadrature in beta is exact for the periodic integrand with a few points; random sampling is not
        need = 6 if axis is not None else min_azimuths
        if model.K4 != 0.0 and n_az < need:
            msg = (f"K4 != 0 but only {n_az} distinct crystal azimuths modulo pi/3 (< {min_azimuths}): the basal "
                   f"anisotropy is under-sampled; increase n_c / n_phi")
            warnings.warn(msg, RuntimeWarning, stacklevel=2)
            notes.append(msg)

    # ---- static states (and the macroscopic field for macro_meanfield)
    grains = []
    for R in Rs:
        grains.append({"R": R})
    B_app = B_lab
    B_int = B_lab.copy()
    converged, n_iter = True, 0
    M_avg = None

    def build(Bl):
        """(gm, Bc, states) for every grain at the lab field Bl; None where no stable equilibrium."""
        out = []
        for g in grains:
            Bc, Nc = to_crystal_frame(g["R"], Bl, N_g)
            gm = model.with_N(Nc)
            out.append((gm, Bc))
        return out

    gms = build(B_int)
    states = []
    for g, (gm, Bc) in zip(grains, gms):
        states.append(_initial_states(gm, Bc, g["R"], B_int, m_init_lab, state, n_grid))
    nfail = sum(1 for s in states if not s)

    def mean_m(sts):
        acc, ws = np.zeros(3), 0.0
        for g, s, w in zip(grains, sts, w_all):
            if not s:
                continue
            mm = g["R"] @ np.mean(np.array(s), axis=0)
            acc += w * mm
            ws += w
        acc = acc / ws if ws > 0 else acc
        if axis is not None:                           # lon-average: only the component along the axis survives
            acc = np.dot(acc, axis) * axis
        return acc

    if demag == "macro_meanfield" and np.any(N_s):
        damp = 0.7
        last_err = math.inf
        for outer in range(3):
            for n_iter in range(1, macro_maxit + 1):
                M_avg = mean_m(states)
                B_new = B_app - MU0 * model.Ms * (N_s @ M_avg)
                err = float(np.linalg.norm(B_new - B_int))
                if err > last_err:
                    damp = max(0.1, 0.5 * damp)
                last_err = err
                B_int = B_int + damp * (B_new - B_int)
                if err < macro_tol_T:
                    break
                gms = build(B_int)
                for k, ((gm, Bc), s) in enumerate(zip(gms, states)):
                    if s:
                        states[k] = [_relax_fast(gm, Bc, m) for m in s]
            else:
                converged = False
            if not converged or state != "global":
                break
            # 'global': re-select the ground state(s) of every grain at the converged field; repeat if anything changed
            gms = build(B_int)
            changed = False
            for k, ((gm, Bc), g) in enumerate(zip(gms, grains)):
                new = _initial_states(gm, Bc, g["R"], B_int, m_init_lab, state, n_grid)
                old = states[k]
                if new and (not old or abs(min(ScalarKernel(gm, Bc).energy(*map(float, m)) for m in new) -
                                          min(ScalarKernel(gm, Bc).energy(*map(float, m)) for m in old))
                            > 1e-9 * gm.energy_scale(float(np.linalg.norm(Bc)))):
                    states[k] = new
                    changed = True
            if not changed:
                break
            last_err = math.inf
        if not converged:
            notes.append(f"macro_meanfield: self-consistent internal field NOT converged in {macro_maxit} iterations "
                         f"(|dB| = {last_err:.3g} T): result is unreliable")
            warnings.warn(notes[-1], RuntimeWarning, stacklevel=2)
        M_avg = mean_m(states)
    elif demag == "macro_meanfield":
        M_avg = mean_m(states)                         # N_sample = 0: nothing to iterate

    # ---- dynamics
    gms = build(B_int)
    acc = np.zeros((len(freqs), 3, 3), complex)
    wsum = 0.0
    for k, ((gm, Bc), g, w) in enumerate(zip(gms, grains, w_all)):
        if not states[k]:
            continue
        # per_grain: the requested kind; macro_meanfield: the grain's own N_grain (0 => external == internal)
        chi = _grain_chi_from_states(gm, g["R"], Bc, states[k], freqs, alpha, gyro,
                                     kind if demag == "per_grain" else "external")
        acc += w * chi
        wsum += w
    if wsum <= 0:
        raise RuntimeError("all grains failed")
    chi_avg = acc / wsum
    if Fa is not None:
        chi_avg = _average_about_axis(chi_avg, Fa)
    if nfail:
        notes.append(f"{nfail} grains without stable equilibrium were dropped; weights renormalised")
    res = EnsembleResult(LEVEL_ENSEMBLE, freqs, chi_avg, kind, len(Rs), nfail, state, texture, notes,
                         demag=demag, quadrature=quad, n_distinct_azimuths_mod_pi3=n_az)
    if demag == "per_grain":
        res.physically_invalid_for_bulk = True
        return res
    res.converged, res.n_iter = converged, n_iter
    res.B_int_lab = B_int
    res.M_avg_A_per_m = None if M_avg is None else model.Ms * M_avg
    res.N_sample = N_s
    res.chi_int_lab = chi_avg
    res.chi_ext_lab = sample_demag_correction(chi_avg, N_s)
    if kind == "external":
        res.chi_lab = res.chi_ext_lab
        notes.append("kind='external': chi_lab is the apparent response of the sample SHAPE N_sample (not a material property)")
    return res


def sample_demag_correction(chi_int: np.ndarray, N_sample: np.ndarray) -> np.ndarray:
    """chi_ext = chi_int (I + N chi_int)^-1 for the sample-shape N (apply to an averaged INTERNAL chi)."""
    I3 = np.eye(3)
    return np.array([c @ np.linalg.inv(I3 + N_sample @ c) for c in chi_int])


# --------------------------------------------------------------------------------- explicit non-implementations

@dataclass(frozen=True)
class NotImplementedResult:
    status: str
    value: object
    reason: str


def _stub(reason: str) -> NotImplementedResult:
    return NotImplementedResult("not_implemented", NotImplemented, reason)


def intergranular_interaction(*args, **kwargs) -> NotImplementedResult:
    """Placeholder: intergranular (dipolar/exchange) coupling is not modelled. Never returns numbers."""
    return _stub("intergranular interactions are outside the non-interacting grain ensemble level")


def porosity_correction(*args, **kwargs) -> NotImplementedResult:
    """Placeholder: porosity/effective-medium correction is not modelled. Never returns numbers."""
    return _stub("porosity / effective-medium mixing is not implemented at this level")


def effective_medium_mixing(*args, **kwargs) -> NotImplementedResult:
    """Placeholder: Bruggeman / Maxwell-Garnett effective-medium mixing of grains and pores. Never returns numbers."""
    return _stub("effective-medium theory (Bruggeman / Maxwell-Garnett) is not implemented")


def schloemann_porosity_linewidth(*args, **kwargs) -> NotImplementedResult:
    """Placeholder: Schloemann porosity / inhomogeneous line broadening. Never returns numbers."""
    return _stub("Schloemann porosity line broadening is not implemented")
