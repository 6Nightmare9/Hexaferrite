"""FMR frequencies of the macrospin: linearised LLG (main method) and Smit-Beljers (cross-check).

Dynamics (Gilbert, SI):  dm/dt = -gamma m x B_eff + alpha m x dm/dt,   B_eff = -(1/Ms) dF/dm  [tesla],
gamma in rad/(s*T) (NOT gamma/2pi). Undamped eigenfrequency about a minimum with sphere Hessian Hs [J/m^3]:

    omega = (gamma / Ms) sqrt(det Hs),        f = omega / (2 pi)                        (LLG, local frame)
    omega = gamma / (Ms sin(theta)) sqrt(F_tt F_pp - F_tp^2)                           (Smit-Beljers)

Smit-Beljers is singular at sin(theta) -> 0 when theta is measured from the crystal c axis; here the
polar axis is switched to a direction perpendicular to m0 (pole='auto'/'local'), where sin(theta') = 1.
The frequency is a property of the point on the sphere, so it does not depend on the choice of pole.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy import optimize

from .core import (Gyro, MagnetModel, as_gyro, eig2_batch, energy_many, field_vector, grad_hess, sphere_derivs,
                   sphere_derivs_batch, tangent_basis, tangent_basis_batch, unit)
from ._fast import ScalarKernel
from .equilibrium import (Equilibrium, _make_eq, classify, find_equilibria, has_axial_symmetry, relax,
                          relax_batch)


# --------------------------------------------------------------------------------- LLG (local frame)

DEBUG_RISK = False
NEUTRAL_TOL_REL = 1e-9        # single tolerance (x energy scale) for "neutral / unstable" everywhere


def omega_from_eigs(lmin, lmax, Ms: float, gamma: float, tol: float) -> float:
    """omega [rad/s] from the two sphere-Hessian eigenvalues with ONE tolerance rule:
    lmin < -tol -> nan (not a minimum); min(|lambda|) <= tol -> 0.0 (neutral zero mode); else gamma/Ms sqrt(lmin lmax)."""
    if lmin < -tol:
        return float("nan")
    if abs(lmin) <= tol or abs(lmax) <= tol:
        return 0.0
    return gamma / Ms * math.sqrt(lmin * lmax)


def llg_omega(model: MagnetModel, B, m, gyro: Gyro, tol_rel: float = NEUTRAL_TOL_REL) -> float:
    """Undamped angular eigenfrequency omega [rad/s] at the point m. lambda_min < -tol: not a minimum -> nan;
    |lambda| <= tol (tol = tol_rel x energy scale): neutral zero mode -> exactly 0.0 (as in `modes_at`)."""
    gt, Hs, _ = sphere_derivs(model, B, m)
    lam = np.linalg.eigvalsh(Hs)
    tol = tol_rel * model.energy_scale(float(np.linalg.norm(B)))
    return omega_from_eigs(float(lam[0]), float(lam[1]), model.Ms, gyro.gamma_rad_per_s_T, tol)


def llg_frequency_Hz(model: MagnetModel, B, m, gyro: Gyro) -> float:
    return llg_omega(model, B, m, gyro) / (2.0 * math.pi)


def smit_beljers_omega(model: MagnetModel, B, m0, gyro: Gyro, pole: str = "auto",
                       sin_switch: float = 0.05) -> tuple[float, dict]:
    """Smit-Beljers omega [rad/s] with explicit treatment of sin(theta) -> 0.

    pole='crystal'  polar axis = c (z). If sin(theta) < 1e-9 the formula is singular -> (nan, info['status']=
                    'singular_pole'); no silent regularisation.
    pole='local'    polar axis perpendicular to m0 (theta'=pi/2): regular everywhere.
    pole='auto'     crystal pole unless sin(theta) < sin_switch, then local.
    Returns (omega, info) with info = {pole_used, sin_theta, status, F_tt, F_pp, F_tp}.
    """
    m0 = unit(m0)
    use = pole
    sin_t_c = math.hypot(m0[0], m0[1])
    if pole == "auto":
        use = "local" if sin_t_c < sin_switch else "crystal"
    if use == "crystal":
        Rp = np.eye(3)
    elif use == "local":
        # right-handed frame (a, b, p): p perpendicular to m0, a = m0, b = p x a  (a x b = p)
        p = tangent_basis(m0)[:, 0]
        a = m0.copy()
        b = np.cross(p, a)
        Rp = np.stack([a, b, p], axis=1)
    else:
        raise ValueError("pole must be crystal|local|auto")
    v = Rp.T @ m0
    th = math.acos(max(-1.0, min(1.0, float(v[2]))))
    ph = math.atan2(v[1], v[0])
    st, ct, sp, cp = math.sin(th), math.cos(th), math.sin(ph), math.cos(ph)
    info = {"pole_used": use, "sin_theta": st}
    if st < 1e-9:
        info["status"] = "singular_pole"
        return float("nan"), info
    m_t = Rp @ np.array([ct * cp, ct * sp, -st])
    m_p = Rp @ np.array([-st * sp, st * cp, 0.0])
    m_tt = -unit(m0)
    m_pp = Rp @ np.array([-st * cp, -st * sp, 0.0])
    m_tp = Rp @ np.array([-ct * sp, ct * cp, 0.0])
    _, g, h = grad_hess(model, B, m0)
    Ftt = float(m_t @ h @ m_t + g @ m_tt)
    Fpp = float(m_p @ h @ m_p + g @ m_pp)
    Ftp = float(m_t @ h @ m_p + g @ m_tp)
    info.update(F_tt=Ftt, F_pp=Fpp, F_tp=Ftp)
    det = Ftt * Fpp - Ftp ** 2
    if det < 0 or Ftt < 0:
        info["status"] = "not_a_minimum"
        return float("nan"), info
    info["status"] = "ok"
    return gyro.gamma_rad_per_s_T / (model.Ms * st) * math.sqrt(det), info


# --------------------------------------------------------------------------------- modes at fixed field

@dataclass(frozen=True)
class Mode:
    f_Hz: float                   # eigenfrequency (undamped), NaN if not a stable minimum
    omega: float
    state: Equilibrium
    is_ground: bool               # lowest-energy state (within tolerance)
    dE: float                     # energy above the ground state [J/m^3]
    status: str                   # ok | neutral_zero_mode | unstable


@dataclass
class ModeSet:
    B_vec: np.ndarray
    modes: list                   # one Mode per local minimum (all branches)
    eq_status: str


def modes_at(model: MagnetModel, B, gyro: Gyro, n_grid: int = 600) -> ModeSet:
    B = np.asarray(B, float)
    eqs = find_equilibria(model, B, n_grid=n_grid)
    out = []
    E0 = eqs.ground.energy if eqs.ground else 0.0
    scale = model.energy_scale(float(np.linalg.norm(B)))
    for e in eqs.minima:
        om = llg_omega(model, B, e.m, gyro)
        st = "ok" if (e.stable and not e.neutral) else ("neutral_zero_mode" if e.stable else "unstable")
        if e.neutral and e.stable:
            om = 0.0          # a zero Hessian eigenvalue (Goldstone / marginal): exact zero mode of the macrospin
        out.append(Mode(om / (2 * math.pi), om, e, abs(e.energy - E0) <= 1e-9 * scale, e.energy - E0, st))
    return ModeSet(B, out, eqs.status)


def f_res(model: MagnetModel, B_T: float, theta_H: float, phi_H: float = 0.0, *, g=None,
          gamma_over_2pi_Hz_per_T=None, gamma_rad_per_s_T=None, gyro: Gyro | None = None,
          n_grid: int = 600) -> ModeSet:
    """Resonance frequencies of ALL local-minimum branches at field magnitude B_T = mu0 H [T], direction
    (theta_H from c, phi_H) in the crystal frame. gamma (g or gamma/2pi or gamma) is mandatory."""
    gy = as_gyro(g, gamma_over_2pi_Hz_per_T, gamma_rad_per_s_T, gyro)
    return modes_at(model, field_vector(B_T, theta_H, phi_H), gy, n_grid=n_grid)


def natural_fmr(model: MagnetModel, *, g=None, gamma_over_2pi_Hz_per_T=None, gamma_rad_per_s_T=None,
                gyro: Gyro | None = None) -> ModeSet:
    """Zero-field FMR (H = 0). A cone state with exact axial symmetry has a Goldstone mode (f = 0,
    status 'neutral_zero_mode'): the uniform macrospin has no finite natural frequency there."""
    gy = as_gyro(g, gamma_over_2pi_Hz_per_T, gamma_rad_per_s_T, gyro)
    return modes_at(model, np.zeros(3), gy)


# --------------------------------------------------------------------------------- branch tracking

@dataclass
class Branch:
    """One local-minimum branch along a field sweep. Samples are contiguous grid indices first_index..last_index.
    B_birth / B_death: spinodal fields (lambda_min = 0, f = 0) at which the branch appears / disappears between grid
    points, found by bisection; None if the branch is cut only by the end of the range or was present at the start."""
    branch_id: int
    B: list = field(default_factory=list)
    f_Hz: list = field(default_factory=list)
    m: list = field(default_factory=list)
    kind: list = field(default_factory=list)
    energy: list = field(default_factory=list)
    first_index: int = 0
    last_index: int = 0
    ended_by_jump: bool = False
    mode_status: list = field(default_factory=list)      # ok | neutral_zero_mode per sample
    B_birth: float | None = None
    B_death: float | None = None
    end_reason: str = "range_end"                         # range_end | spinodal | merged (pitchfork into another branch)


def _nudged(m):
    """Deterministic tiny transverse nudge (1e-9 rad): an exact-symmetry saddle must not be mistaken for a branch."""
    x, y, u = float(m[0]), float(m[1]), float(m[2])
    if abs(x) < 0.9:
        ex, ey, ez = 1.0 - x * x, -x * y, -x * u
    else:
        ex, ey, ez = -y * x, 1.0 - y * y, -y * u
    n = math.sqrt(ex * ex + ey * ey + ez * ez)
    return (x + 1e-9 * ex / n, y + 1e-9 * ey / n, u + 1e-9 * ez / n)


def _step_one(model, Bv, m_prev, jump_angle, kern=None):
    """Continue ONE state to field vector Bv: (m (ndarray), alive, lmin, lmax). `alive` = converged, lambda_min >=
    -tol (stable, possibly neutral) and moved by <= jump_angle."""
    kern = kern or ScalarKernel(model, Bv)
    m, ok, lmin, lmax = kern.relax(_nudged(m_prev))
    cosang = m[0] * m_prev[0] + m[1] * m_prev[1] + m[2] * m_prev[2]
    ang = math.acos(max(-1.0, min(1.0, cosang)))
    alive = ok and lmin >= -NEUTRAL_TOL_REL * kern.scale and ang <= jump_angle
    return np.array(m), bool(alive), lmin, lmax


def _advance(model, d, Ba, ma, Bb, jump_angle):
    """Follow the branch state ma (a minimum at B = Ba, along direction d) to Bb by adaptive continuation. A step is
    refused if the relaxed state is not a minimum or moved by > jump_angle; the step is then halved. Returns
    (True, Bb, m) if the branch reaches Bb, else (False, B_end, m_end) with B_end the last field at which it is
    still a minimum (the spinodal field, lambda_min -> 0, to ~1e-7 of the interval)."""
    cur, mcur = Ba, np.asarray(ma, float)
    total = Bb - Ba
    if total == 0.0:
        return True, Bb, mcur
    m, alive, _, _ = _step_one(model, Bb * d, mcur, jump_angle)
    if alive:
        return True, Bb, m
    sgn = 1.0 if total > 0 else -1.0
    step = abs(total) * 0.5
    tiny = 1e-7 * abs(total)
    for _ in range(200):
        left = abs(Bb - cur)
        if left <= tiny * 1e-3:
            return True, Bb, mcur
        trial = cur + sgn * min(step, left)
        m, alive, _, _ = _step_one(model, trial * d, mcur, jump_angle)
        if alive:
            cur, mcur = trial, m
            step = min(2.0 * step, abs(total))
        else:
            step *= 0.5
            if step < tiny:
                return False, cur, mcur
    return False, cur, mcur


def _merge_field(model, d, Ba, ma, Bb, other, Bo, jump_angle, match_tol, sym):
    """Field in [Ba, Bb] at which the branch state (minimum at Ba) merges into the state of branch `other` (pitchfork
    where a symmetric pair of minima fuses with the central state), or stops being a minimum -- by bisection."""
    lo, mlo, hi = Ba, np.asarray(ma, float), Bb
    mo = other.m[-2]                      # the other branch's state at the previous grid point (Bo)
    span = abs(Bb - Ba)
    for _ in range(40):
        if abs(hi - lo) <= 1e-7 * span:
            break
        mid = 0.5 * (lo + hi)
        ms, alive, _, _ = _step_one(model, mid * d, mlo, jump_angle)
        merged = not alive
        if alive:
            ok_o, _, m_o = _advance(model, d, Bo, mo, mid, jump_angle)
            merged = ok_o and np.linalg.norm(_canon(ms, sym) - _canon(m_o, sym)) < match_tol
        if merged:
            hi = mid
        else:
            lo, mlo = mid, ms
    return float(lo)


def _classify_py(m, Bv, ang_tol: float = 1e-3) -> str:
    """Same rule as `classify`, plain floats."""
    th = math.acos(max(-1.0, min(1.0, float(m[2]))))
    bn = math.sqrt(float(Bv[0]) ** 2 + float(Bv[1]) ** 2 + float(Bv[2]) ** 2)
    if th < ang_tol or math.pi - th < ang_tol:
        return "axial"
    if abs(th - math.pi / 2) < ang_tol:
        return "planar"
    if bn == 0.0 or math.hypot(float(Bv[0]), float(Bv[1])) <= 1e-9 * max(bn, 1e-30):
        return "cone"
    return "tilted_by_field"


def _sample_info(model, Bv, m, gyro, kern=None):
    """(f_Hz, mode_status, kind, energy) of the minimum m at field vector Bv."""
    kern = kern or ScalarKernel(model, Bv)
    g0, g1, a, b, d, _, _ = kern.derivs(float(m[0]), float(m[1]), float(m[2]))
    lmin, lmax, _ = kern.eig(a, b, d)
    om = omega_from_eigs(lmin, lmax, model.Ms, gyro.gamma_rad_per_s_T, NEUTRAL_TOL_REL * kern.scale)
    st = "neutral_zero_mode" if om == 0.0 else "ok"
    return (om / (2.0 * math.pi), st, _classify_py(m, Bv), kern.energy(float(m[0]), float(m[1]), float(m[2])))


def _append(br, b, m, info, i, front=False):
    f, st, kd, en = info
    if front:
        br.B.insert(0, float(b)); br.f_Hz.insert(0, f); br.m.insert(0, np.array(m))
        br.kind.insert(0, kd); br.energy.insert(0, en); br.mode_status.insert(0, st)
        br.first_index = i
    else:
        br.B.append(float(b)); br.f_Hz.append(f); br.m.append(np.array(m))
        br.kind.append(kd); br.energy.append(en); br.mode_status.append(st)
        br.last_index = i


def _canon(m, sym):
    return np.array([math.hypot(m[0], m[1]), 0.0, m[2]]) if sym else np.asarray(m, float)


def _sweep(model, direction, B_grid, gyro, n_grid, jump_angle, match_tol, coarse_step, lam_event_rel,
           verify_global):
    d = unit(direction)
    Bg = np.asarray(B_grid, float)
    n = len(Bg)
    branches: list[Branch] = []
    active: list[Branch] = []
    lam_prev: dict[int, float] = {}
    notes: list[str] = []
    incomplete = False
    coarse = set(range(0, n, max(1, coarse_step))) | {n - 1}
    for i in range(n):
        b = float(Bg[i])
        Bv = b * d
        kern = ScalarKernel(model, Bv)
        sym = has_axial_symmetry(model, Bv)
        taken: list[np.ndarray] = []
        owner: list[Branch] = []
        event = False
        if active:
            still = []
            for br in active:
                mj, ok, lmin, lmax = _step_one(model, Bv, br.m[-1], jump_angle, kern)
                dup = None
                if ok:
                    for t, o in zip(taken, owner):
                        if np.linalg.norm(_canon(mj, sym) - t) < match_tol:
                            dup = o
                            break
                if dup is None and not ok:
                    ok2, B_end, m_end = _advance(model, d, Bg[i - 1], br.m[-1], b, jump_angle)
                    if ok2:
                        mj, ok = m_end, True
                        lmin = kern.eig(*kern.derivs(*[float(c) for c in mj])[2:5])[0]
                        for t, o in zip(taken, owner):
                            if np.linalg.norm(_canon(mj, sym) - t) < match_tol:
                                dup = o
                                break
                    else:
                        br.B_death, br.end_reason, br.ended_by_jump = float(B_end), "spinodal", True
                        event = True
                        continue
                if dup is not None:                         # pitchfork: this branch merges into the one that took the state
                    br.B_death = _merge_field(model, d, Bg[i - 1], br.m[-1], b, dup, Bg[i - 1], jump_angle,
                                              match_tol, sym)
                    br.end_reason = "merged"
                    event = True
                    continue
                _append(br, b, mj, _sample_info(model, Bv, mj, gyro, kern), i)
                taken.append(_canon(mj, sym))
                owner.append(br)
                still.append(br)
                # a soft mode that keeps softening: a partner minimum may be about to appear
                if lmin < lam_event_rel * kern.scale and lmin < lam_prev.get(br.branch_id, math.inf):
                    event = True
                lam_prev[br.branch_id] = lmin
            active = still
        if i in coarse or event or not active:
            eqs = find_equilibria(model, Bv, n_grid=n_grid, verify=verify_global)
            if eqs.possibly_incomplete:
                incomplete = True
                notes.append(f"global search flagged possibly_incomplete at B={b:.6g} T")
            for e in eqs.minima:
                if any(np.linalg.norm(_canon(e.m, sym) - t) < match_tol for t in taken):
                    continue
                br = Branch(branch_id=len(branches), first_index=i, last_index=i)
                branches.append(br)
                _append(br, b, e.m, _sample_info(model, Bv, e.m, gyro, kern), i)
                taken.append(_canon(e.m, sym))
                owner.append(br)
                # backward tracking to the spinodal field where the branch appears (or to the first grid point)
                j, m_cur = i, e.m
                while j > 0:
                    ok2, B_end, m_end = _advance(model, d, Bg[j], m_cur, Bg[j - 1], jump_angle)
                    if not ok2:
                        br.B_birth = float(B_end)
                        break
                    j -= 1
                    m_cur = m_end
                    Bj = Bg[j] * d
                    symj = has_axial_symmetry(model, Bj)
                    if any(o is not br and o.first_index <= j <= o.last_index and
                           np.linalg.norm(_canon(o.m[j - o.first_index], symj) - _canon(m_end, symj)) < match_tol
                           for o in branches):
                        notes.append(f"branch {br.branch_id} merges with an existing branch at index {j}")
                        break
                    _append(br, Bg[j], m_end, _sample_info(model, Bj, m_end, gyro), j, front=True)
                active.append(br)
    return branches, notes, incomplete


def sweep_branches(model: MagnetModel, direction, B_grid, gyro: Gyro, n_grid: int = 200,
                   jump_angle: float = 0.5, match_tol: float = 1e-3, coarse_step: int = 8,
                   lam_event_rel: float = 0.02, verify_global: bool = False) -> list[Branch]:
    """Frequencies of every local-minimum branch along the field sweep B_grid (tesla, along `direction`).

    Branches are continued from grid point to grid point (relaxation of the previous state, path bisection if the
    direct step fails -> spinodal field B_death). The GLOBAL search for minima runs only on every `coarse_step`-th
    grid point and on events (a branch has lambda_min < lam_event_rel x scale, or a branch disappeared); a newly found
    minimum is tracked backwards to its spinodal field B_birth. Limitation: a pair of minima that exist only between
    two global searches and are not near an active branch's soft mode may be missed (use coarse_step=1 to be exhaustive).
    Samples exist only where the branch is a local minimum; neutral modes have f = 0, status 'neutral_zero_mode'."""
    return _sweep(model, direction, B_grid, gyro, n_grid, jump_angle, match_tol, coarse_step, lam_event_rel,
                  verify_global)[0]


# --------------------------------------------------------------------------------- resonance fields

@dataclass(frozen=True)
class ResonanceRoot:
    B_res_T: float                # mu0 H at resonance [T]
    H_res_A_per_m: float
    branch_id: int
    state_kind: str
    is_ground: bool               # ground state at the resonance field (False = metastable branch)
    df_dB_Hz_per_T: float         # local slope of the branch (for width conversion)
    f_Hz: float


@dataclass
class ResonanceResult:
    status: str                   # ok | possibly_incomplete | no_resonance_in_range
    roots: list
    f_target_Hz: float
    B_range_T: tuple
    branches: list
    notes: list = field(default_factory=list)
    possibly_incomplete: bool = False


class _BranchEval:
    """f(B) of one branch by continuation from the nearest stored sample (valid for B up to the spinodal field)."""

    def __init__(self, model, d, br, gyro, jump_angle):
        self.model, self.d, self.br, self.gyro, self.ja = model, d, br, gyro, jump_angle
        self.Bs = np.array(br.B)

    def state(self, B):
        k = int(np.argmin(np.abs(self.Bs - B)))
        ok, Bend, m = _advance(self.model, self.d, self.Bs[k], self.br.m[k], B, self.ja)
        return ok, m

    def covers(self, B):
        """B lies inside the field interval where the branch is known to exist (samples + spinodal ends)."""
        lo, hi = float(self.Bs.min()), float(self.Bs.max())
        if self.br.B_birth is not None:
            lo, hi = min(lo, self.br.B_birth), max(hi, self.br.B_birth)
        if self.br.B_death is not None:
            lo, hi = min(lo, self.br.B_death), max(hi, self.br.B_death)
        return lo <= B <= hi

    def energy_at(self, B):
        """Energy of the branch state at B (None where the branch does not exist)."""
        if not self.covers(B):
            return None
        ok, m = self.state(B)
        if not ok:
            return None
        return ScalarKernel(self.model, B * self.d).energy(float(m[0]), float(m[1]), float(m[2]))

    def f(self, B):
        ok, m = self.state(B)
        if not ok:
            return 0.0                                     # beyond the spinodal field the branch is absent: f -> 0
        return llg_omega(self.model, B * self.d, m, self.gyro) / (2.0 * math.pi)


def _branch_roots(model, d, br, ft, gyro, jump_angle, Bspan, n_grid, notes, extra_fields=()):
    """Roots f_branch(B) = ft on one branch: (list of B, evaluator, risk flag). extra_fields: spinodal / merging fields of
    OTHER branches; where this branch exists there it is evaluated too (a branch that merges into this one at a
    pitchfork has a soft mode, hence a dip of f(B), at exactly that field)."""
    ev = _BranchEval(model, d, br, gyro, jump_angle)
    nodes = [(b, f, False) for b, f in zip(br.B, br.f_Hz) if math.isfinite(f)]
    if len(nodes) < 1:
        return [], None, False
    if br.B_birth is not None:
        nodes.insert(0, (br.B_birth, 0.0, True))
    if br.B_death is not None:
        nodes.append((br.B_death, 0.0, True))
    sign = 1.0 if nodes[-1][0] >= nodes[0][0] else -1.0
    for Bx in extra_fields:
        if ev.covers(Bx) and min(abs(Bx - n[0]) for n in nodes) > 1e-9 * max(Bspan, 1e-12):
            nodes.append((float(Bx), ev.f(Bx), False))
    nodes.sort(key=lambda t: sign * t[0])
    # refine discrete extrema that sit on the wrong side of the target (a dip/bump can hide a pair of roots)
    out_nodes = []
    resolved = set()                       # nodes next to a refined extremum: that interval is resolved
    for k, nd in enumerate(nodes):
        out_nodes.append(nd)
        if 0 < k < len(nodes) - 1 and not nd[2]:
            f0, f1, f2 = nodes[k - 1][1], nd[1], nodes[k + 1][1]
            is_min = f1 < f0 and f1 < f2 and f1 > ft
            is_max = f1 > f0 and f1 > f2 and f1 < ft
            if is_min or is_max:
                lo, hi = sorted((nodes[k - 1][0], nodes[k + 1][0]))
                sgn = 1.0 if is_min else -1.0
                res = optimize.minimize_scalar(lambda x: sgn * ev.f(x), bounds=(lo, hi), method="bounded",
                                               options={"xatol": 1e-12 * max(Bspan, 1e-12)})
                fe = ev.f(float(res.x))
                if (is_min and fe < ft) or (is_max and fe > ft):
                    out_nodes.append((float(res.x), fe, False))
                    resolved.update((nd[0], float(res.x)))
    nodes = sorted(out_nodes, key=lambda t: sign * t[0])
    virt = {nd[0]: nd[1] for nd in nodes if nd[2]}

    def g(B):
        return (virt[B] if B in virt else ev.f(B)) - ft

    roots = []
    for (Ba, fa, _), (Bb, fb, _) in zip(nodes[:-1], nodes[1:]):
        if (fa - ft >= 0) == (fb - ft >= 0):
            continue
        ga, gb = fa - ft, fb - ft
        try:
            x = optimize.brentq(lambda B: ga if B == Ba else (gb if B == Bb else g(B)), Ba, Bb,
                                xtol=1e-13 * max(1.0, abs(Ba), abs(Bb)), rtol=1e-13, maxiter=200)
        except ValueError:
            notes.append(f"branch {br.branch_id}: bracketing failed near B={Ba:.6g} T")
            continue
        if abs(g(x)) > 1e-5 * ft:
            notes.append(f"branch {br.branch_id}: bracket near B={x:.6g} T closed on a discontinuity, root rejected")
            continue
        roots.append(x)
    # Resolution check: could a pair of roots hide inside an interval without a sign change? Work with F = f^2 =
    # (gamma/Ms)^2 det(Hs): it is smooth through the spinodal square-root cusps of f(B). The chord deviation of F over
    # an interval is bounded by h^2/8 |F''| (|F''| from the second divided differences of the neighbouring nodes).
    risk = False
    Ft = ft * ft
    for j in range(len(nodes) - 1):
        (Ba, fa, _), (Bb, fb, _) = nodes[j], nodes[j + 1]
        if (fa - ft >= 0) != (fb - ft >= 0) or Ba in resolved or Bb in resolved:
            continue
        h = abs(Bb - Ba)
        d2 = []
        for (p, q, r) in ((j - 1, j, j + 1), (j, j + 1, j + 2)):
            if p >= 0 and r < len(nodes):
                x0, x1, x2 = nodes[p][0], nodes[q][0], nodes[r][0]
                F0, F1, F2 = nodes[p][1] ** 2, nodes[q][1] ** 2, nodes[r][1] ** 2
                if x0 != x1 and x1 != x2 and x0 != x2:
                    d2.append(abs(2.0 * ((F2 - F1) / (x2 - x1) - (F1 - F0) / (x1 - x0)) / (x2 - x0)))
        if d2 and min(abs(fa * fa - Ft), abs(fb * fb - Ft)) < max(d2) * h * h / 8.0:
            risk = True
            if DEBUG_RISK:
                print('RISK', br.branch_id, j, (Ba, fa), (Bb, fb), d2, h)
    return roots, ev, risk


def H_res(model: MagnetModel, f_target_Hz: float, theta_H: float, phi_H: float = 0.0, *,
          B_range_T: tuple, g=None, gamma_over_2pi_Hz_per_T=None, gamma_rad_per_s_T=None,
          gyro: Gyro | None = None, n_points: int = 241, n_grid: int = 200, coarse_step: int = 8,
          verify_global: bool = False, jump_angle: float = 0.5) -> ResonanceResult:
    """All resonance fields mu0 H_res [T] with f_branch(B) = f_target on B_range_T, over ALL branches
    (stable and metastable, flagged). B_range_T is mu0 H, tesla.

    Per branch: sign changes of f - f_target between samples are refined by Brent; discrete extrema of f(B) on the
    wrong side of the target are located (Brent minimisation) so that a narrow dip/bump hiding a pair of roots is
    resolved; the end of a branch (spinodal field, lambda_min = 0, found by bisection) is used as a node with f = 0.
    Status: 'ok' | 'no_resonance_in_range' (no branch crosses the target and nothing suspicious) |
    'possibly_incomplete' (the resolution cannot exclude hidden roots: curvature of f(B) between samples is large
    compared with the distance to the target, or a global search was flagged). The end of the range is never
    substituted for a root. Completeness of the branch SET relies on `n_grid` (see `sweep_branches`)."""
    gy = as_gyro(g, gamma_over_2pi_Hz_per_T, gamma_rad_per_s_T, gyro)
    d = field_vector(1.0, theta_H, phi_H)
    Bg = np.linspace(B_range_T[0], B_range_T[1], n_points)
    branches, notes, incomplete = _sweep(model, d, Bg, gy, n_grid, jump_angle, 1e-3, coarse_step, 0.02, verify_global)
    Bspan = float(abs(Bg[-1] - Bg[0]))
    roots: list[ResonanceRoot] = []
    risk_any = False
    evs = {br.branch_id: _BranchEval(model, d, br, gy, jump_angle) for br in branches}
    for br in branches:
        extra = sorted({float(o.B_death) for o in branches if o is not br and o.B_death is not None} |
                       {float(o.B_birth) for o in branches if o is not br and o.B_birth is not None})
        xs, ev, risk = _branch_roots(model, d, br, f_target_Hz, gy, jump_angle, Bspan, n_grid, notes, extra)
        risk_any = risk_any or risk
        for bres in xs:
            ok, mm = ev.state(bres)
            hstep = 1e-5 * max(Bspan, 1e-12)
            fp = ev.f(bres + hstep)
            fm = ev.f(bres - hstep)
            okp, _ = ev.state(bres + hstep)
            okm, _ = ev.state(bres - hstep)
            if okp and okm:
                slope = (fp - fm) / (2 * hstep)
            elif okp:
                slope = (fp - f_target_Hz) / hstep
            else:
                slope = (f_target_Hz - fm) / hstep
            scale = model.energy_scale(abs(bres))
            eq = _make_eq(model, bres * d, mm, NEUTRAL_TOL_REL * scale, 1e-3)
            e_all = [e for e in (evs[o.branch_id].energy_at(bres) for o in branches) if e is not None]
            is_g = (not e_all) or eq.energy - min(e_all) <= 1e-9 * scale
            roots.append(ResonanceRoot(float(bres), float(bres / _mu0()), br.branch_id, eq.kind, bool(is_g),
                                       float(slope), float(f_target_Hz)))
    roots.sort(key=lambda r: (r.B_res_T, r.branch_id))
    flag = bool(risk_any or incomplete)
    if flag:
        notes.append("possibly_incomplete: resolution cannot exclude additional roots (refine n_points / n_grid)")
    status = "possibly_incomplete" if flag else ("ok" if roots else "no_resonance_in_range")
    return ResonanceResult(status, roots, float(f_target_Hz), (float(B_range_T[0]), float(B_range_T[1])),
                           branches, notes, flag)


def _mu0() -> float:
    from .core import MU0
    return MU0


def df_dB(model: MagnetModel, B_T: float, theta_H: float, phi_H: float, m0, gyro: Gyro, h: float = 1e-6) -> float:
    """Branch slope df/d(mu0 H) [Hz/T] at the state continuing from m0 (central difference)."""
    d = field_vector(1.0, theta_H, phi_H)

    def f(b):
        mm = relax(model, b * d, m0)
        return llg_omega(model, b * d, mm, gyro) / (2 * math.pi)
    return (f(B_T + h) - f(B_T - h)) / (2 * h)
