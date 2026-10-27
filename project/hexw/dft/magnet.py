"""Magnetic systems, exchange bonds / classes and the choice of well-conditioned collinear configurations for the J fit.

A `MagneticSystem` is an ORDERED structure (a configuration of `hexw.structures.configs.generate_configurations`, or the
hematite smoke cell) plus the magnetic bookkeeping: spin length per atom, flip group per atom (the unit that is turned over
as a whole: a W-sublattice W1..W7, or a single Fe atom in hematite) and the sign of every atom in the ground state
(Gorter scheme for the W phase).

Bonds are Fe/Co/Mn pairs closer than r_cut. They are collected into classes (columns of the design matrix) either by the pair
of flip groups ("sublattice_pair") or by the distance shell ("shell"). Classes whose column does not change over the pool of
flip patterns are absorbed into E0 (typical: intra-sublattice bonds under sublattice flips).
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from hexw.dft import hamiltonian as H

PATH_NAMES = {0: "no shared O", 1: "corner (1 O)", 2: "edge (2 O)", 3: "face (3 O)"}


class DesignError(ValueError):
    pass


@dataclass
class MagneticSystem:
    name: str
    kind: str                           # key of config 'systems' (W92 / hematite)
    structure: object                   # pymatgen Structure, ordered, with oxidation states
    spin: np.ndarray                    # (n_atoms,) S, 0 for non-magnetic atoms
    group: np.ndarray                   # (n_atoms,) flip-group index, -1 for non-magnetic atoms
    group_names: list
    base_sign: np.ndarray               # (n_atoms,) +-1 in the ground state, 0 for non-magnetic
    meta: dict = field(default_factory=dict)

    @property
    def n_atoms(self) -> int:
        return len(self.structure)

    def sigma(self, s: np.ndarray) -> np.ndarray:
        """Per-atom signs for the group flip pattern s (+-1 per group)."""
        out = np.zeros(self.n_atoms)
        mag = self.group >= 0
        out[mag] = self.base_sign[mag] * np.asarray(s)[self.group[mag]]
        return out


# ---- builders -------------------------------------------------------------------------------------------------------
def spin_of_site(site) -> float:
    sp = site.specie
    try:
        el, ox = sp.symbol, int(round(sp.oxi_state))
    except Exception:                                            # no oxidation state
        return 0.0
    if el in ("Sr", "O"):
        return 0.0
    return H.spin_length(el, ox)


def system_from_w_config(cfgset, index: int, name: Optional[str] = None) -> MagneticSystem:
    """One ordered W-phase configuration -> MagneticSystem (flip groups = the seven W-sublattices)."""
    from hexw.structures.configs import load_template
    from hexw.structures.sites import SITES, SPIN_SIGN
    cfg = cfgset.configs[index]
    st = cfg.structure
    t = load_template(cfgset.supercell)
    n0 = len(t.other_species)
    if len(st) != n0 + len(t.tm_frac):
        raise DesignError("configuration does not match the template (atom count)")
    spin = np.zeros(len(st)); group = -np.ones(len(st), int); base = np.zeros(len(st))
    for j in range(len(t.tm_frac)):
        k = n0 + j
        d = st[k].frac_coords - t.tm_frac[j]
        if np.abs(d - np.round(d)).max() > 1e-6:
            raise DesignError(f"atom {k} does not sit on template position {j}")
        S = spin_of_site(st[k])
        if S > 0:
            spin[k] = S
            group[k] = int(t.tm_site[j])
            base[k] = SPIN_SIGN[int(t.tm_site[j])]
    return MagneticSystem(name or cfg.candidate.candidate_id.replace("|", "_"), "W92", st, spin, group, list(SITES),
                          base, meta={"candidate_id": cfg.candidate.candidate_id, "weight": cfg.weight,
                                      "count_vector": cfg.count_vector, "supercell": list(cfgset.supercell)})


# ---- bonds ----------------------------------------------------------------------------------------------------------
def _cluster_distances(d: np.ndarray, tol: float = 0.03) -> np.ndarray:
    order = np.argsort(d)
    lab = np.zeros(len(d), int)
    k = 0
    for a, b in zip(order[:-1], order[1:]):
        if d[b] - d[a] > tol:
            k += 1
        lab[b] = k
    return lab


def find_bonds(system: MagneticSystem, r_cut: float, r_O: float = 2.6) -> list[H.Bond]:
    """All magnetic pairs closer than r_cut (periodic images included; each pair once; a pair of an atom with its own
    image is skipped: it is constant in every collinear configuration)."""
    st = system.structure
    mag = [i for i in range(len(st)) if system.spin[i] > 0]
    nbrs = st.get_all_neighbors(r_cut + 1e-6, sites=[st[i] for i in mag])
    # oxygen neighbours (home cell coordinates) for the path analysis
    O_idx = [i for i, s in enumerate(st) if s.specie.symbol == "O"]
    o_nb = st.get_all_neighbors(r_O, sites=[st[i] for i in mag])
    O_home = {}
    for i, lst in zip(mag, o_nb):
        O_home[i] = np.array([n.coords for n in lst if st[n.index].specie.symbol == "O"]).reshape(-1, 3)
    lat = st.lattice.matrix
    bonds = []
    for i, lst in zip(mag, nbrs):
        for n in lst:
            j = int(n.index)
            if system.spin[j] <= 0 or j <= i:           # j < i is found from j's list; j == i is constant
                continue
            R = tuple(int(v) for v in np.round(n.image))
            # shared oxygens: O of i (home) vs O of j (home + R)
            Oj = O_home[j] + np.array(R) @ lat
            n_sh, ang = 0, []
            if len(O_home[i]) and len(Oj):
                D = np.linalg.norm(O_home[i][:, None, :] - Oj[None, :, :], axis=2)
                for a, b in zip(*np.nonzero(D < 0.1)):
                    n_sh += 1
                    po = O_home[i][a]
                    v1 = st[i].coords - po
                    v2 = n.coords - po
                    ang.append(math.degrees(math.acos(np.clip(v1 @ v2 / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1, 1))))
            bonds.append(H.Bond(i, j, R, float(n.nn_distance), -1, float(system.spin[i]), float(system.spin[j]),
                                n_sh, float(np.mean(ang)) if ang else None))
    return bonds


@dataclass
class ExchangeModel:
    system: MagneticSystem
    bonds: list
    class_keys: list                    # key per class
    class_labels: list
    class_info: list                    # dict per class
    absorbed: list                      # info on classes absorbed into E0 / neglected
    pool_patterns: np.ndarray           # (n_pool, n_groups) +-1, row 0 = ground
    pool_X: np.ndarray                  # (n_pool, n_classes)
    active_groups: list
    group_by: str
    r_cut: float
    notes: list = field(default_factory=list)
    far_labels: list = field(default_factory=list)      # nuisance columns: sublattice pairs WITHOUT a short bond (see far_columns)
    far_info: list = field(default_factory=list)
    pool_X_far: Optional[np.ndarray] = None             # (n_pool, n_far): column = -sum over far bonds of the pair S_i S_j s_i s_j

    @property
    def n_classes(self) -> int:
        return len(self.class_labels)


def far_columns(system: MagneticSystem, bonds_near: list, r_cut: float, r_far: float, sig: np.ndarray) -> tuple:
    """Nuisance columns for the sublattice pairs a-b (a != b) that have NO bond shorter than r_cut but do have bonds between r_cut
    and r_far. Their energy -sum S_i S_j s_i s_j J_ij is exactly proportional (per pair of sublattices) to a column the model does
    not contain, so it aliases onto the fitted classes. Pairs that also have short bonds need no column: all bonds of one
    pair of sublattices have the same s_i s_j in every sublattice-flip pattern, so their far part is absorbed EXACTLY into the
    effective J of the class (docs/DFT.md, section 3). Returns (labels, X (n_patterns, q), info)."""
    gname = system.group_names
    near_types = {tuple(sorted((int(system.group[b.i]), int(system.group[b.j])))) for b in bonds_near}
    far = [b for b in find_bonds(system, r_far) if b.distance > r_cut + 1e-6]
    by_type: dict = {}
    for b in far:
        t = tuple(sorted((int(system.group[b.i]), int(system.group[b.j]))))
        if t[0] != t[1] and t not in near_types:
            by_type.setdefault(t, []).append(b)
    labels, cols, info = [], [], []
    for t in sorted(by_type):
        bl = by_type[t]
        for k, b in enumerate(bl):
            b.cls = 0
        col = H.design_matrix(bl, sig, 1)[:, 0]
        if (col.max() - col.min()) < 1e-9:                     # constant over the flip pool: absorbed into E0, not identifiable either
            continue
        labels.append(f"far:{gname[t[0]]}-{gname[t[1]]}")
        cols.append(col)
        d = np.array([b.distance for b in bl])
        info.append({"label": labels[-1], "n_bonds_far": len(bl), "d_min": float(d.min()), "d_max": float(d.max()),
                     "sum_SS": float(sum(b.S_i * b.S_j for b in bl))})
    X = np.array(cols).T if cols else np.zeros((len(sig), 0))
    return labels, X, info


def build_model(system: MagneticSystem, r_cut: float, group_by: str = "sublattice_pair",
                max_pairs: Optional[int] = None, r_far: Optional[float] = None) -> ExchangeModel:
    bonds = find_bonds(system, r_cut)
    if not bonds:
        raise DesignError("no magnetic bonds found")
    gname = system.group_names
    # provisional keys
    if group_by == "sublattice_pair":
        keys = [tuple(sorted((int(system.group[b.i]), int(system.group[b.j])))) for b in bonds]
    elif group_by == "shell":
        d = np.array([b.distance for b in bonds])
        lab = _cluster_distances(d)
        keys = [("shell", int(l)) for l in lab]
    else:
        raise DesignError(f"unknown group_by {group_by!r}")
    uniq = sorted(set(keys))
    # pool of flip patterns over the groups that take part in a bond
    groups_in = sorted({int(system.group[b.i]) for b in bonds} | {int(system.group[b.j]) for b in bonds})
    G = len(system.group_names)
    pats = []
    rest = list(itertools.product([1, -1], repeat=len(groups_in) - 1))
    for r in rest:
        s = np.ones(G)
        for g, v in zip(groups_in, (1,) + tuple(r)):
            s[g] = v
        pats.append(s)
    pats = np.array(pats)                                          # row 0 = all +1 = ground state
    sig = np.array([system.sigma(s) for s in pats])
    # provisional classes
    for b, k in zip(bonds, keys):
        b.cls = uniq.index(k)
    X = H.design_matrix(bonds, sig, len(uniq))
    const = (X.max(0) - X.min(0)) < 1e-9
    absorbed, kept = [], []
    for c, k in enumerate(uniq):
        info = _class_info(system, [b for b, kk in zip(bonds, keys) if kk == k], k, group_by, gname)
        (absorbed if const[c] else kept).append((c, k, info))
    # classes with proportional columns over the whole pool cannot be told apart by any flip pattern: keep the first one
    # (master) and report the combination J_master + sum r_m J_m that the data determine
    merged_notes = []
    masters: list = []
    for t in kept:
        col = X[:, t[0]]
        for m in masters:
            colm = X[:, m[0][0]]
            r = float(col @ colm / (colm @ colm))
            if np.linalg.norm(col - r * colm) <= 1e-8 * max(np.linalg.norm(col), 1e-30):
                m.append((t, r))
                break
        else:
            masters.append([t])
    kept = []
    for m in masters:
        t0 = m[0]
        if len(m) > 1:
            comb = " + ".join([f"J[{t0[2]['label']}]"] + [f"{r:.4g}*J[{t[2]['label']}]" for t, r in m[1:]])
            t0[2]["combination"] = comb
            t0[2]["merged_with"] = [t[2]["label"] for t, _ in m[1:]]
            t0[2]["label"] = t0[2]["label"] + "+" + "+".join(t[2]["label"] for t, _ in m[1:])
            merged_notes.append(f"proportional columns merged: the fit determines {comb}")
            for t, _ in m[1:]:
                absorbed.append((t[0], t[1], dict(t[2], reason="proportional to class " + t0[2]["label"] + " (merged)")))
        kept.append(t0)
    # truncation to the leading pairs (shortest distance)
    neglected = []
    if max_pairs is not None and len(kept) > max_pairs:
        order = sorted(kept, key=lambda t: (t[2]["d_min"], -t[2]["n_bonds"]))
        neglected = order[max_pairs:]
        kept = [t for t in kept if t not in neglected]
    new_index = {c: n for n, (c, _, _) in enumerate(kept)}
    for b, k in zip(bonds, keys):
        c = uniq.index(k)
        b.cls = new_index.get(c, -1)
    X = H.design_matrix(bonds, sig, len(kept))
    notes = list(merged_notes)
    if neglected:
        notes.append("classes beyond max_pairs are neglected (their energy goes into the residual): "
                     + ", ".join(t[2]["label"] for t in neglected))
    mixed = [t[2]["label"] for t in kept if t[2].get("n_distinct_shells", 1) > 1 and group_by == "sublattice_pair"]
    if mixed:
        notes.append("J of these classes is the bond-weighted mean over several distance shells (not separable by "
                     "sublattice flips): " + ", ".join(mixed))
    far_labels, far_X, far_info = [], np.zeros((len(sig), 0)), []
    if r_far is not None and group_by == "sublattice_pair" and r_far > r_cut:
        far_labels, far_X, far_info = far_columns(system, bonds, r_cut, float(r_far), sig)
        if far_labels:
            notes.append(f"{len(far_labels)} nuisance columns for sublattice pairs without a bond < {r_cut:g} A (bonds up to "
                         f"{r_far:g} A): " + ", ".join(far_labels))
    return ExchangeModel(system, bonds, [t[1] for t in kept], [t[2]["label"] for t in kept], [t[2] for t in kept],
                         [t[2] if "reason" in t[2] else dict(t[2], reason="constant under the flip pool (absorbed into E0)")
                          for t in absorbed]
                         + [dict(t[2], reason="beyond max_pairs (neglected)") for t in neglected],
                         pats, X, groups_in, group_by, r_cut, notes, far_labels, far_info, far_X)


def _class_info(system, blist, key, group_by, gname) -> dict:
    d = np.array([b.distance for b in blist])
    lab = _cluster_distances(d) if len(d) > 1 else np.zeros(1, int)
    shells = [(float(np.mean(d[lab == l])), int((lab == l).sum())) for l in sorted(set(lab))]
    paths = {}
    for b in blist:
        paths[b.n_shared_O] = paths.get(b.n_shared_O, 0) + 1
    angs = [b.angle_deg for b in blist if b.angle_deg is not None and b.n_shared_O == 1]
    if group_by == "sublattice_pair":
        label = f"{gname[key[0]]}-{gname[key[1]]}"
    else:
        label = f"shell{key[1]}_{np.mean(d):.2f}A"
    return {"label": label, "n_bonds": len(blist), "d_min": float(d.min()), "d_max": float(d.max()), "shells": shells,
            "n_distinct_shells": len(shells), "paths": {PATH_NAMES.get(k, str(k)): v for k, v in sorted(paths.items())},
            "mean_Fe_O_Fe_angle_corner_deg": float(np.mean(angs)) if angs else None,
            "spin_products": sorted({round(b.S_i * b.S_j, 3) for b in blist})}


# ---- choice of configurations ---------------------------------------------------------------------------------------
def _scaled_cond(X: np.ndarray) -> float:
    c = H.conditioning(X)
    return c["cond_scaled"]


def _logdet(A: np.ndarray) -> float:
    An, _ = H.scaled(A)
    M = An.T @ An
    M = M + 1e-9 * np.eye(len(M))
    sign, ld = np.linalg.slogdet(M)
    return ld if sign > 0 else -math.inf


@dataclass
class Design:
    model: ExchangeModel
    pool_index: list                    # indices into the pool, in run order (row 0 = ground)
    roles: list                         # 'fit' | 'holdout'
    config_ids: list
    flips: list                         # list of group-name lists
    patterns: np.ndarray                # (n, n_groups) +-1
    X: np.ndarray                       # (n, n_classes)
    cond_fit: dict
    pool_cond: dict
    leverage_holdout: list
    notes: list = field(default_factory=list)
    X_far: Optional[np.ndarray] = None  # (n, q) nuisance columns of the chosen configurations (far sublattice pairs)

    @property
    def fit_rows(self):
        return [k for k, r in enumerate(self.roles) if r == "fit"]

    @property
    def holdout_rows(self):
        return [k for k, r in enumerate(self.roles) if r == "holdout"]


def _hamming(a, b):
    d = int((a != b).sum())
    return min(d, len(a) - d)               # s and -s are the same state


def select_design(model: ExchangeModel, n_extra_fit: int = 1, n_holdout: int = 2) -> Design:
    X = model.pool_X
    n_pool, p_cls = X.shape
    p = p_cls + 1
    pool_cond = H.conditioning(X)
    if pool_cond["rank"] < p:
        raise DesignError(f"the full pool of {n_pool} flip patterns has rank {pool_cond['rank']} < {p} parameters "
                          "(classes are not separable by these flips): merge classes / reduce max_pairs")
    n_fit = p + n_extra_fit
    # patterns with the same X row are equivalent for the model (and, by the crystal symmetry, usually for the DFT energy):
    # they carry no new information, so only one of them can be run; a hold-out must have a row the fit set does not have
    uniq: list = []
    for c in range(n_pool):
        if not any(np.allclose(X[c], X[u], rtol=1e-9, atol=1e-9 * max(1.0, np.abs(X).max())) for u in uniq):
            uniq.append(c)
    if n_fit + n_holdout > len(uniq):
        raise DesignError(f"only {len(uniq)} distinct design rows (of {n_pool} patterns) for {n_fit} fit + {n_holdout} hold-out "
                          f"configurations ({p} parameters): reduce the number of J classes (r_cut / max_pairs) or n_holdout")
    A_all = H.augmented(X)
    chosen = [0]
    while len(chosen) < n_fit:                       # greedy D-optimal
        best, best_v = None, -math.inf
        for c in uniq:
            if c in chosen:
                continue
            v = _logdet(A_all[chosen + [c]])
            if v > best_v:
                best, best_v = c, v
        chosen.append(best)
    # exchange refinement on the scaled condition number (ground state stays)
    cur = _scaled_cond(X[chosen])
    for _ in range(30):
        improved = False
        for pos in range(1, len(chosen)):
            for c in uniq:
                if c in chosen:
                    continue
                trial = chosen.copy(); trial[pos] = c
                v = _scaled_cond(X[trial])
                if v < cur * (1 - 1e-9):
                    chosen, cur, improved = trial, v, True
        if not improved:
            break
    rest = [c for c in uniq if c not in chosen]
    hold = []
    An = H.augmented(X[chosen])
    Minv = np.linalg.pinv(An.T @ An)
    while len(hold) < n_holdout:
        def score(c):
            dmin = min(_hamming(model.pool_patterns[c], model.pool_patterns[f]) for f in chosen)
            lev = float(A_all[c] @ Minv @ A_all[c])
            return (dmin, -lev)
        c = max((c for c in rest if c not in hold), key=score)
        hold.append(c)
    order = chosen + hold
    gname = model.system.group_names
    ids, flips = [], []
    for k, c in enumerate(order):
        fl = [gname[g] for g in model.active_groups if model.pool_patterns[c][g] < 0]
        flips.append(fl)
        ids.append(f"c{k:02d}_" + ("ground" if not fl else "flip_" + "_".join(fl)))
    lev = [float(A_all[c] @ Minv @ A_all[c]) for c in hold]
    cond_fit = H.conditioning(X[chosen])
    notes = []
    if cond_fit["cond_scaled"] > 1e3:
        notes.append(f"scaled condition number {cond_fit['cond_scaled']:.3g} is large")
    return Design(model, order, ["fit"] * len(chosen) + ["holdout"] * len(hold), ids, flips,
                  model.pool_patterns[order], X[order], cond_fit, pool_cond, lev, notes,
                  None if model.pool_X_far is None else model.pool_X_far[order])
