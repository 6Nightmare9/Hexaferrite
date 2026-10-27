"""Heisenberg convention of the DFT adapter and the least-squares mapping of collinear DFT energies onto it.

CONVENTION (fixed, TZ section 8):

    H = - sum_{i<j} J_ij  S_i . S_j

* every pair (i, j) is counted ONCE (i < j; a pair with a periodic image of the partner is a separate bond);
* S_i are spin VECTORS of length |S_i| = S (dimensionless spin quantum number: Fe3+ 5/2, Fe2+ 2, Mn2+ 5/2, Mn3+ 2,
  Co2+ 3/2, Ni2+ 1, Zn2+/Mg2+ 0), NOT unit vectors;
* J > 0 is ferromagnetic (parallel spins lower the energy), J < 0 antiferromagnetic;
* energies and J in meV (1 Ry = 13605.693122994 meV); the fit is for the energy PER SIMULATION CELL.

Collinear configuration c with signs sigma_i = +-1:
    E_c = E0 + sum_g J_g X_{c,g},     X_{c,g} = - sum_{bonds b in class g} S_i S_j sigma_i sigma_j.

A class g is a set of bonds that share one J (same pair of sublattices, or the same distance shell). The J of a class
is the bond-weighted mean of the true J_b inside the class when the bonds are not strictly equivalent.
Intra-class bonds that join the same flip group have sigma_i sigma_j = +1 in every configuration: they only shift E0
(unobservable by sublattice flips).

Conversion to other conventions: `convert_J`.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

RY_TO_MEV = 13605.693122994           # CODATA 2018 Rydberg energy in eV * 1000
KB_MEV_PER_K = 8.617333262e-2          # meV/K

CONVENTION_TEXT = ("H = -sum_{i<j} J_ij S_i.S_j; each pair once; |S| = spin quantum number "
                   "(Fe3+ 5/2, Fe2+ 2, Mn2+ 5/2, Co2+ 3/2); J>0 ferromagnetic; J in meV; energies per simulation cell")

# spin quantum numbers by (element, formal charge); Zn/Mg are non-magnetic
SPIN_S = {("Fe", 3): 2.5, ("Fe", 2): 2.0, ("Mn", 2): 2.5, ("Mn", 3): 2.0, ("Co", 2): 1.5, ("Ni", 2): 1.0,
          ("Zn", 2): 0.0, ("Mg", 2): 0.0}


def spin_length(element: str, charge: int) -> float:
    try:
        return SPIN_S[(element, int(charge))]
    except KeyError:
        raise ValueError(f"no spin length for {element}{charge:+d}") from None


# ---- conventions ----------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Convention:
    """H = sign * sum over pairs [counting] J * (S_i.S_j if spin == 'S' else e_i.e_j)."""
    sign: int = -1                    # -1: H = -sum J ... (J>0 ferro);  +1: H = +sum J ... (J>0 antiferro)
    counting: str = "once"            # 'once': sum_{i<j};  'twice': sum_{i!=j} (each pair twice)
    spin: str = "S"                   # 'S': vectors of length |S|;  'unit': unit vectors

    def __post_init__(self):
        if self.sign not in (-1, 1) or self.counting not in ("once", "twice") or self.spin not in ("S", "unit"):
            raise ValueError(f"bad convention {self}")


OURS = Convention(-1, "once", "S")


def convert_J(J: float, S_i: float, S_j: float, src: Convention, dst: Convention) -> float:
    """J of the pair (i, j) in convention `src` expressed in convention `dst` (same physical energy of the pair)."""
    def pair_energy_factor(c: Convention) -> float:
        # energy of one pair per unit cos(theta) = sign * J * mult * norm
        return c.sign * (2.0 if c.counting == "twice" else 1.0) * (S_i * S_j if c.spin == "S" else 1.0)
    f_src, f_dst = pair_energy_factor(src), pair_energy_factor(dst)
    if f_dst == 0:
        raise ValueError("zero spin: the pair has no exchange energy")
    return J * f_src / f_dst


def J_to_kelvin(J_meV: float) -> float:
    return J_meV / KB_MEV_PER_K


def energy_per_pair_parallel(J: float, S_i: float, S_j: float) -> float:
    """Exchange energy (meV) of a parallel pair in our convention: -J S_i S_j."""
    return -J * S_i * S_j


# ---- design matrix --------------------------------------------------------------------------------------------------
@dataclass
class Bond:
    i: int
    j: int
    image: tuple                       # lattice vector of the image of j
    distance: float
    cls: int                           # class index (column of X), -1: absorbed (constant in all configurations)
    S_i: float
    S_j: float
    n_shared_O: int = 0
    angle_deg: Optional[float] = None


def design_row(bonds: list[Bond], sigma: np.ndarray, n_classes: int) -> np.ndarray:
    """X_{c,g} = - sum_{b in g} S_i S_j sigma_i sigma_j for one collinear configuration (signs per atom)."""
    x = np.zeros(n_classes)
    for b in bonds:
        if b.cls >= 0:
            x[b.cls] -= b.S_i * b.S_j * sigma[b.i] * sigma[b.j]
    return x


def design_matrix(bonds: list[Bond], sigmas: np.ndarray, n_classes: int) -> np.ndarray:
    """(n_config, n_classes), vectorised: sigmas is (n_config, n_atoms)."""
    idx = [b for b in bonds if b.cls >= 0]
    X = np.zeros((len(sigmas), n_classes))
    if not idx:
        return X
    bi = np.array([b.i for b in idx]); bj = np.array([b.j for b in idx])
    w = np.array([b.S_i * b.S_j for b in idx]); cl = np.array([b.cls for b in idx])
    prod = sigmas[:, bi] * sigmas[:, bj] * w[None, :]
    for g in range(n_classes):
        X[:, g] = -prod[:, cl == g].sum(1)
    return X


def augmented(X: np.ndarray) -> np.ndarray:
    """[1, X]: the first column belongs to E0."""
    return np.hstack([np.ones((len(X), 1)), X])


def scaled(A: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Column-normalised matrix (unit 2-norm columns) and the scale factors."""
    s = np.linalg.norm(A, axis=0)
    s[s == 0] = 1.0
    return A / s, s


def conditioning(X: np.ndarray) -> dict:
    """Condition number and singular values of the scaled augmented matrix [1, X] (columns of unit norm) and of the
    raw one. rank is numerical (tol = max(shape) * eps * s_max)."""
    A = augmented(X)
    An, _ = scaled(A)
    s_raw = np.linalg.svd(A, compute_uv=False)
    s = np.linalg.svd(An, compute_uv=False)
    tol = max(A.shape) * np.finfo(float).eps * s[0]
    rank = int((s > tol).sum())
    cond = float(s[0] / s[-1]) if s[-1] > tol else math.inf
    cond_raw = float(s_raw[0] / s_raw[-1]) if s_raw[-1] > 0 else math.inf
    return {"cond_scaled": cond, "cond_raw": cond_raw, "singular_values_scaled": [float(v) for v in s],
            "singular_values_raw": [float(v) for v in s_raw], "rank": rank, "n_columns": A.shape[1],
            "n_rows": A.shape[0]}


# ---- least squares with uncertainties ---------------------------------------------------------------------------------
@dataclass
class FitResult:
    names: list
    E0: float
    J: np.ndarray
    sigma_J: np.ndarray               # 1 sigma (meV); with nuisance columns the aliasing of the omitted pairs is marginalised in
    sigma_E0: float
    cov: np.ndarray                   # covariance of [E0, J...]
    residuals: np.ndarray             # fit energies - model (meV), fit set
    rms_meV: float
    dof: int
    s_resid_meV: float                # sqrt(RSS / dof) (nan if dof == 0)
    sigma_used_meV: float             # max(s_resid, floor): the scale of the covariance
    sigma_basis: str                  # 'residual' | 'noise_floor'
    cond: dict
    loo_residuals: Optional[np.ndarray] = None          # leave-one-out prediction errors (meV), nan if the fit is singular
    holdout: Optional[dict] = None
    notes: list = field(default_factory=list)
    holdout_consistent: Optional[bool] = None           # hold-out rms / energy range < holdout_rel_tol (None: no hold-out)
    sigma_model_meV: float = 0.0                        # hold-out rms error added in quadrature when the hold-out disagrees
    sigma_J_total: Optional[np.ndarray] = None          # sigma_J with the model error included (== sigma_J if consistent)
    nuisance: Optional[dict] = None                     # far-pair nuisance columns: names, prior sigma, posterior mean / sigma, inflation

    def predict(self, X: np.ndarray) -> np.ndarray:
        return self.E0 + np.asarray(X) @ self.J


def _lstsq(A: np.ndarray, y: np.ndarray):
    sol, *_ = np.linalg.lstsq(A, y, rcond=None)
    return sol


def _ridge_solve(A: np.ndarray, N: Optional[np.ndarray], y: np.ndarray, sigma: float, tau: Optional[float]):
    """MAP estimate of [theta (flat prior), gamma ~ N(0, tau^2)] for y = A theta + N gamma + noise(sigma). Returns
    (theta, gamma, M^-1) where sigma^2 M^-1 is the posterior covariance of [theta, gamma]."""
    p = A.shape[1]
    if N is None or N.shape[1] == 0:
        return _lstsq(A, y), np.zeros(0), np.linalg.inv(A.T @ A)
    B = np.hstack([A, N])
    pen = np.zeros(B.shape[1])
    pen[p:] = (sigma / tau) ** 2
    M = B.T @ B + np.diag(pen)
    Minv = np.linalg.inv(M)
    sol = Minv @ (B.T @ y)
    return sol[:p], sol[p:], Minv


def fit_exchange(X: np.ndarray, E_meV: np.ndarray, names: list, *, noise_floor_meV: float = 0.0,
                 X_holdout: Optional[np.ndarray] = None, E_holdout_meV: Optional[np.ndarray] = None,
                 max_cond: float = 1e3, X_nuis: Optional[np.ndarray] = None, X_nuis_holdout: Optional[np.ndarray] = None,
                 nuis_names: Optional[list] = None, nuis_prior_sigma_meV: Optional[float] = None,
                 holdout_rel_tol: float = 0.02) -> FitResult:
    """OLS fit E = E0 + X J. Covariance = sigma^2 (A^T A)^-1, sigma = max(sqrt(RSS/dof), noise floor); with dof == 0 only the
    noise floor can be used (the fit is exactly determined and says nothing about its own quality). LOO cross-validation on
    the fit set and, when given, prediction errors on the hold-out configurations (never used in the fit).

    Nuisance columns (`X_nuis`, one per omitted sublattice pair that has no short bond): E = E0 + X J + X_nuis g, with the
    Gaussian prior g ~ N(0, nuis_prior_sigma_meV^2) (an ASSUMPTION of the caller). The covariance of [E0, J] is the marginal of the
    joint posterior, so sigma(J) includes the aliasing of these pairs onto the fitted classes. dof stays n - p (the nuisance
    parameters are shrunk by the prior and are not counted as fitted).

    Hold-out criterion (relative, not z): rms of the hold-out errors / (max - min of the fit energies) < holdout_rel_tol
    (2 %). A z-score against the noise floor is not used: the Heisenberg-model error itself (7.3 meV on hematite) exceeds any
    SCF noise."""
    X = np.asarray(X, float); E = np.asarray(E_meV, float)
    n, p = X.shape[0], X.shape[1] + 1
    A = augmented(X)
    cond = conditioning(X)
    notes = []
    if cond["rank"] < p:
        raise ValueError(f"design matrix is rank deficient (rank {cond['rank']} < {p} parameters): add configurations or merge classes")
    if n < p:
        raise ValueError(f"{n} configurations for {p} parameters")
    if cond["cond_scaled"] > max_cond:
        notes.append(f"scaled condition number {cond['cond_scaled']:.3g} exceeds {max_cond:g}")
    Nm = None
    if X_nuis is not None and np.size(X_nuis):
        Nm = np.asarray(X_nuis, float)
        if Nm.shape[0] != n:
            raise ValueError("X_nuis must have one row per fit configuration")
        if not (nuis_prior_sigma_meV and nuis_prior_sigma_meV > 0):
            raise ValueError("nuisance columns need a positive prior sigma")
    dof = n - p
    tau = nuis_prior_sigma_meV
    sig_ridge = max(noise_floor_meV, 1e-12)
    for _ in range(30):                              # sigma enters the ridge strength: fixed point (converges in a few steps)
        sol, gam, Minv = _ridge_solve(A, Nm, E, sig_ridge, tau)
        res = E - A @ sol - (Nm @ gam if Nm is not None else 0.0)
        rss = float(res @ res)
        s_res = math.sqrt(rss / dof) if dof > 0 else math.nan
        sigma = max(s_res, noise_floor_meV) if dof > 0 else noise_floor_meV
        if Nm is None or abs(sigma - sig_ridge) < 1e-9 * max(sigma, 1e-12):
            break
        sig_ridge = max(sigma, 1e-12)
    if dof > 0:
        basis = "residual" if s_res >= noise_floor_meV else "noise_floor"
    else:
        basis = "noise_floor"
        notes.append("dof = 0: exactly determined fit, uncertainties come from the assumed noise floor only")
    if sigma <= 0:
        raise ValueError("no error scale: dof = 0 and noise floor = 0")
    covf = sigma ** 2 * Minv
    cov = covf[:p, :p]
    loo = np.full(n, math.nan)
    if n > p:
        for k in range(n):
            m = np.arange(n) != k
            Ak = A[m]
            if np.linalg.matrix_rank(Ak) < p:
                continue
            solk, gk, _ = _ridge_solve(Ak, None if Nm is None else Nm[m], E[m], sigma, tau)
            loo[k] = E[k] - A[k] @ solk - (Nm[k] @ gk if Nm is not None else 0.0)
    out = FitResult(names=list(names), E0=float(sol[0]), J=sol[1:].copy(), sigma_J=np.sqrt(np.diag(cov))[1:],
                    sigma_E0=float(math.sqrt(cov[0, 0])), cov=cov, residuals=res,
                    rms_meV=float(math.sqrt(rss / n)), dof=dof, s_resid_meV=s_res, sigma_used_meV=sigma,
                    sigma_basis=basis, cond=cond, loo_residuals=loo, notes=notes)
    if Nm is not None:
        base_cov = sigma ** 2 * np.linalg.inv(A.T @ A)
        sj0 = np.sqrt(np.diag(base_cov))[1:]
        out.nuisance = {"names": list(nuis_names or [f"n{k}" for k in range(Nm.shape[1])]), "prior_sigma_meV": float(tau),
                        "posterior_mean_meV": gam.tolist(), "posterior_sigma_meV": np.sqrt(np.diag(covf)[p:]).tolist(),
                        "sigma_J_without_nuisance_meV": sj0.tolist(),
                        "sigma_inflation": (out.sigma_J / sj0).tolist()}
    span = float(np.ptp(E))
    if X_holdout is not None and E_holdout_meV is not None and len(E_holdout_meV):
        Xh = np.asarray(X_holdout, float)
        Ah = augmented(Xh)
        Bh = Ah if Nm is None else np.hstack([Ah, np.asarray(X_nuis_holdout, float)])
        theta = np.concatenate([sol, gam]) if Nm is not None else sol
        pred = Bh @ theta
        err = np.asarray(E_holdout_meV, float) - pred
        var = np.einsum("ij,jk,ik->i", Bh, covf, Bh) + sigma ** 2
        rms_h = float(math.sqrt(np.mean(err ** 2)))
        out.holdout = {"error_meV": err.tolist(), "predicted_meV": pred.tolist(), "sigma_pred_meV": np.sqrt(var).tolist(),
                       "z": (err / np.sqrt(var)).tolist(), "rms_meV": rms_h,
                       "max_abs_meV": float(np.max(np.abs(err))), "energy_range_meV": span,
                       "rel_rms": rms_h / span if span > 0 else math.inf, "rel_tol": holdout_rel_tol,
                       "criterion": f"rms of hold-out errors / energy range of the fit < {holdout_rel_tol:g}"}
        out.holdout_consistent = bool(span > 0 and rms_h / span < holdout_rel_tol)
        if not out.holdout_consistent:
            # the Heisenberg model with these classes does not predict configurations outside the fit within the tolerance: the
            # hold-out rms is a lower bound of the model error, added in quadrature to the energy error scale of the fit
            out.sigma_model_meV = rms_h
            tot = math.sqrt(sigma ** 2 + out.sigma_model_meV ** 2)
            out.sigma_J_total = tot / sigma * out.sigma_J
            notes.append(f"hold-out disagrees (rms error {rms_h:.3g} meV = {100 * out.holdout['rel_rms']:.2g} % of the energy range "
                         f"{span:.4g} meV, tolerance {100 * holdout_rel_tol:g} %): the spin model is incomplete for these classes; "
                         "model error added to the J uncertainties")
    if out.sigma_J_total is None:
        out.sigma_J_total = out.sigma_J.copy()
    return out
