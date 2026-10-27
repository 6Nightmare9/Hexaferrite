"""Ordered configurations for DFT from the A-cation ensemble (TZ section 7: a partially occupied CIF is NOT a DFT configuration).

Procedure
  1. Geometry ONLY from the L02 CIF of SrZn2Fe16O27 (hexagonal P6_3/mmc, 300 K), `raw/theory_side/06_literature/L02_Morch2023`.
     The half-occupied split 4f pair of W7 (0.33 A apart) is replaced by the
     2d point on the mirror plane (shift 0.16 A); relaxation is expected to restore the displacement.
  2. Integer cell counts of every A species on every site type are drawn from the fit ensemble (one apportionment per
     bootstrap draw; largest remainder, site capacity respected); the frequency of a count vector is its weight.
  3. For every count vector all placements on the 36 TM positions are generated and reduced to symmetry-inequivalent
     ones with the P6_3/mmc operations found by spglib (pymatgen); an orbit weighs (its size / number of placements):
     ideal mixing inside a count vector (no short-range-order energy is modelled).
  4. The best `n_rep` (2..5) orbits are returned with renormalised weights and the covered share of the probability.
  5. Fe2+ (n = 2 - sum x_2+ [-1 per Mn3+ ... see ensemble]) is placed on Fe positions of the chosen scenario by a farthest-point
     rule (Fe2+ ions as far apart as possible): a placeholder, not an enumeration - flagged in the result.
Spin-only starting moments (Gorter signs) are attached as the `magmom` site property.
"""
from __future__ import annotations

import itertools
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

from hexw.io import paths
from hexw.params import occupancy as O
from hexw.params.moments import SPIN_MOMENT
from hexw.registry import schemas as S
from hexw.structures import ensemble as E
from hexw.structures.sites import N_CELL, SITES, SPIN_SIGN

CIF_REL = Path("raw/theory_side/06_literature/L02_Morch2023/supporting/CIFs/SrZn2Fe16O27.cif")
MAX_ARRANGEMENTS = 600_000
N_REP_MIN, N_REP_MAX = 2, 5


class ConfigError(ValueError):
    pass


@dataclass
class Template:
    lattice: object                    # pymatgen Lattice
    tm_frac: np.ndarray                # (36, 3)
    tm_site: np.ndarray                # (36,) 0..6 (W1..W7)
    other_species: list                # [(symbol, frac(3))]
    source: str


@dataclass
class DftConfiguration:
    candidate: S.StructureCandidate    # ordered: integer occupancies, partial_occupancy flag False
    structure: object                  # pymatgen Structure (oxidation states, magmom)
    weight: float                      # renormalised over the returned representatives
    weight_raw: float                  # probability under the ensemble (before truncation)
    count_vector: dict                 # {species: [atoms on W1..W7 in the cell]}
    orbit_size: int
    n_atoms: int
    notes: list = field(default_factory=list)


@dataclass
class ConfigSet:
    composition: dict
    hypothesis: E.Hypothesis
    configs: list
    coverage: float                    # sum of raw weights of the returned configurations
    n_orbits_total: int
    supercell: tuple
    notes: list = field(default_factory=list)


def load_template(supercell: tuple = (1, 1, 1)) -> Template:
    """The L02 SrZn2W CIF as an ordered template (all 36 TM positions, Sr, O)."""
    import warnings

    from pymatgen.io.cif import CifParser
    p = paths.get_config().root / CIF_REL
    if "L02" not in str(p) or "L01" in str(p):
        raise ConfigError("geometry must come from the L02 CIF")
    with paths.open_data(p, "r", purpose="structure") as f:
        txt = f.read()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        st = CifParser.from_str(txt).parse_structures(primitive=False, on_error="raise")[0]
    tm, other = [], []
    w7 = []
    for s in st:
        names = {sp.symbol for sp in s.species}
        if names <= {"Sr", "O"}:
            other.append((next(iter(names)), np.array(s.frac_coords)))
            continue
        m = re.search(r"(\d)$", s.label)
        if not m:
            raise ConfigError(f"cannot read the W index from label {s.label!r}")
        idx = int(m.group(1)) - 1
        if idx == 6:
            w7.append(np.array(s.frac_coords))
        else:
            tm.append((idx, np.array(s.frac_coords)))
    used = [False] * len(w7)
    for i, a in enumerate(w7):                 # merge the split pairs on the mirror plane z = 1/4, 3/4
        if used[i]:
            continue
        for j in range(i + 1, len(w7)):
            if not used[j] and np.allclose(a[:2], w7[j][:2], atol=1e-3):
                used[i] = used[j] = True
                tm.append((6, np.array([a[0], a[1], 0.5 * (a[2] + w7[j][2])])))
                break
    if not all(used):
        raise ConfigError("unpaired W7 positions in the CIF")
    tm.sort(key=lambda t: t[0])
    tm_site = np.array([t[0] for t in tm])
    tm_frac = np.array([t[1] for t in tm])
    if len(tm) != 36 or any((tm_site == i).sum() != N_CELL[i] for i in range(7)):
        raise ConfigError(f"unexpected number of TM positions: {np.bincount(tm_site)}")
    t = Template(st.lattice, tm_frac, tm_site, other, str(CIF_REL))
    if supercell != (1, 1, 1):
        t = _supercell(t, supercell)
    return t


def _supercell(t: Template, sc: tuple) -> Template:
    from pymatgen.core import Lattice
    na, nb, nc = sc
    shifts = [np.array(s) for s in itertools.product(range(na), range(nb), range(nc))]
    scale = np.array([na, nb, nc], float)

    def rep(fr):
        return [(fr + s) / scale for s in shifts]
    tm_frac, tm_site = [], []
    for fr, si in zip(t.tm_frac, t.tm_site):
        for f2 in rep(fr):
            tm_frac.append(f2)
            tm_site.append(si)
    other = [(sym, f2) for sym, fr in t.other_species for f2 in rep(fr)]
    lat = Lattice(t.lattice.matrix * scale[:, None])
    return Template(lat, np.array(tm_frac), np.array(tm_site), other, t.source)


def symmetry_permutations(t: Template) -> np.ndarray:
    """(n_ops, n_tm): perm[k, j] = index of the TM position onto which operation k maps position j."""
    from pymatgen.core import Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    species = [e for e, _ in t.other_species] + ["Fe"] * len(t.tm_frac)
    coords = [f for _, f in t.other_species] + list(t.tm_frac)
    st = Structure(t.lattice, species, coords)
    sga = SpacegroupAnalyzer(st, symprec=1e-3)
    ops = sga.get_symmetry_operations(cartesian=False)
    perms = []
    for op in ops:
        new = op.operate_multi(t.tm_frac)
        d = new[:, None, :] - t.tm_frac[None, :, :]
        d -= np.round(d)
        match = np.linalg.norm(d @ t.lattice.matrix, axis=2) < 1e-2
        if not (match.sum(0) == 1).all() or not (match.sum(1) == 1).all():
            continue                       # operation not a symmetry of the TM sub-lattice (should not happen)
        perms.append(match.argmax(1))
    if not perms:
        raise ConfigError("no symmetry operations found")
    return np.unique(np.array(perms), axis=0)


# ---- integer counts from the ensemble ---------------------------------------------------------------------------
def _apportion(expected: np.ndarray, totals: np.ndarray, capacity: np.ndarray) -> np.ndarray:
    """expected (S,7) atoms in the cell, totals (S,) integers, capacity (7,) sites of each type in the cell."""
    k = np.floor(expected + 1e-9).astype(int)
    frac = expected - k
    for s in range(expected.shape[0]):
        need = int(totals[s] - k[s].sum())
        order = np.argsort(-frac[s])
        guard = 0
        while need > 0 and guard < 10_000:
            moved = False
            for i in order:
                if need == 0:
                    break
                if k[:, i].sum() < capacity[i]:
                    k[s, i] += 1
                    need -= 1
                    moved = True
            guard += 1
            if not moved:
                raise ConfigError("cannot accommodate the A cations (capacity)")
    return k


def count_vector_distribution(species: tuple, draws: np.ndarray, x: dict, z_cell: float, capacity: np.ndarray):
    """{count-vector tuple: probability} from the (n, S, 7) fraction draws."""
    totals = np.array([x[sp] * z_cell for sp in species])
    if np.any(np.abs(totals - np.round(totals)) > 1e-6):
        raise ConfigError(f"content {x} is not an integer number of atoms in the cell (Z = {z_cell:g}); use a supercell")
    totals = np.round(totals).astype(int)
    out: dict = {}
    for d in draws:
        kv = _apportion(d * capacity[None, :], totals, capacity)
        key = tuple(map(tuple, kv))
        out[key] = out.get(key, 0) + 1
    n = len(draws)
    return {k: v / n for k, v in out.items()}


def _site_placements(positions: np.ndarray, counts: list) -> np.ndarray:
    """All ways to put counts[s] atoms of species s+1 on `positions` (len n): rows of species labels (0 = none)."""
    n = len(positions)
    rows = [np.zeros(n, np.int8)]
    for s, c in enumerate(counts, start=1):
        new = []
        for r in rows:
            free = np.flatnonzero(r == 0)
            for comb in itertools.combinations(free, c):
                r2 = r.copy()
                r2[list(comb)] = s
                new.append(r2)
        rows = new
    return np.array(rows)


def _canonical_orbits(arr: np.ndarray, perms: np.ndarray, base: int):
    """Orbit representatives of the rows of arr under perms; returns (representatives, orbit sizes)."""
    n = arr.shape[1]
    half = n // 2
    w_hi = base ** np.arange(half - 1, -1, -1, dtype=np.int64)
    w_lo = base ** np.arange(n - half - 1, -1, -1, dtype=np.int64)
    best_hi = np.full(len(arr), np.iinfo(np.int64).max, np.int64)
    best_lo = best_hi.copy()
    arg = np.zeros(len(arr), int)
    for k, perm in enumerate(perms):
        tr = np.empty_like(arr)
        tr[:, perm] = arr                   # atom on position j moves to perm[j]
        hi = tr[:, :half].astype(np.int64) @ w_hi
        lo = tr[:, half:].astype(np.int64) @ w_lo
        better = (hi < best_hi) | ((hi == best_hi) & (lo < best_lo))
        best_hi = np.where(better, hi, best_hi)
        best_lo = np.where(better, lo, best_lo)
    keys = np.stack([best_hi, best_lo], 1)
    uniq, first, cnt = np.unique(keys, axis=0, return_index=True, return_counts=True)
    return arr[first], cnt


def _fe2_positions(t: Template, occupied_A: np.ndarray, n2: int, allowed_sites: tuple) -> list:
    """Farthest-point choice of n2 Fe positions among the allowed site types."""
    cand = [j for j in range(len(t.tm_frac)) if occupied_A[j] == 0 and t.tm_site[j] in allowed_sites]
    if n2 > len(cand):
        raise ConfigError("not enough Fe positions for Fe2+ in the chosen scenario")
    if n2 == 0:
        return []
    chosen = [cand[0]]
    f = t.tm_frac
    def dist(i, j):
        d = f[i] - f[j]
        d -= np.round(d)
        return float(np.linalg.norm(d @ t.lattice.matrix))
    while len(chosen) < n2:
        best = max((c for c in cand if c not in chosen), key=lambda c: min(dist(c, q) for q in chosen))
        chosen.append(best)
    return chosen


def _build(t: Template, species_names: tuple, mn_valence: Optional[int], labels: np.ndarray, fe2: list, z_cell: float):
    from pymatgen.core import Species, Structure
    symbols, coords, magmom = [], [], []
    for sym, fr in t.other_species:
        symbols.append(Species(sym, {"Sr": 2, "O": -2}[sym]))
        coords.append(fr)
        magmom.append(0.0)
    sites = []
    for j, (fr, si) in enumerate(zip(t.tm_frac, t.tm_site)):
        lab = int(labels[j])
        if lab == 0:
            key = "Fe2+" if j in fe2 else "Fe3+"
            sp = Species("Fe", 2 if j in fe2 else 3)
            mom = SPIN_MOMENT[key]
        else:
            name = species_names[lab - 1]
            q = float(mn_valence) if name == "Mn" else 2.0
            sp = Species(name, int(q))
            key = "Mn3" if (name == "Mn" and mn_valence == 3) else ("Mn2" if name == "Mn" else name)
            mom = SPIN_MOMENT[key]
        symbols.append(sp)
        coords.append(fr)
        magmom.append(float(SPIN_SIGN[si] * mom))
        sites.append((j, si, sp, key))
    st = Structure(t.lattice, symbols, coords, site_properties={"magmom": magmom})
    cand_sites = []
    for j, si, sp, key in sites:
        label = key if key in ("Fe2+", "Fe3+") else sp.symbol
        cand_sites.append(S.SiteOccupancy(site=f"{SITES[si]}#{j}", species={label: 1.0}, multiplicity=1,
                                          wyckoff=None, site_id=SITES[si]))
    return st, cand_sites


def generate_configurations(composition: dict, hypothesis: Optional[E.Hypothesis] = None, n_rep: int = 5,
                            supercell: tuple = (1, 1, 1), fit: Optional[O.OccupancyFit] = None, n_draws: int = 200,
                            seed: int = 0) -> ConfigSet:
    """2-5 symmetry-inequivalent ordered configurations (92 atoms for the hexagonal cell) with weights."""
    comp = O.clean_composition(composition)
    hyp = hypothesis or E.Hypothesis()
    has_mn = "Mn" in comp
    mnv = (hyp.mn_valence or 2) if has_mn else None
    scen = hyp.fe2_scenario or "uniform"
    if scen not in E.FE2_SCENARIOS:
        raise ConfigError(f"unknown Fe2+ scenario {scen!r}")
    n_rep = int(min(max(n_rep, N_REP_MIN), N_REP_MAX))
    n2, n3, ox, counts = E.valence_counts(comp, mnv)
    t = load_template(supercell)
    vol = int(np.prod(supercell))
    z_cell = 2.0 * vol
    if abs(n2 * z_cell - round(n2 * z_cell)) > 1e-6:
        raise ConfigError(f"n(Fe2+) = {n2:g} per f.u. is not an integer in the cell; use a supercell")
    species, draws = O.site_fraction_draws(comp, fit, n_draws, seed)
    capacity = np.array([(t.tm_site == i).sum() for i in range(7)])
    dist = count_vector_distribution(species, draws, comp, z_cell, capacity) if species else {(): 1.0}
    perms = symmetry_permutations(t)

    items = []                                   # (raw weight, count key, orbit rep, orbit size)
    n_orbits = 0
    for key, prob in sorted(dist.items(), key=lambda kv: -kv[1]):
        kv = np.array(key) if species else np.zeros((0, 7), int)
        blocks = []
        total = 1
        for i in range(7):
            pos = np.flatnonzero(t.tm_site == i)
            blocks.append((pos, _site_placements(pos, [int(c) for c in kv[:, i]]) if species else np.zeros((1, len(pos)), np.int8)))
            total *= len(blocks[-1][1])
            if total > MAX_ARRANGEMENTS:
                break
        if total > MAX_ARRANGEMENTS:
            continue                             # too many placements for this (rare) count vector: skipped, reported
        full = np.zeros((1, len(t.tm_frac)), np.int8)
        for pos, rows in blocks:
            na, nr = len(full), len(rows)
            full = np.repeat(full, nr, 0)
            full[:, pos] = np.tile(rows, (na, 1))
        reps, cnt = _canonical_orbits(full, perms, base=len(species) + 1)
        n_orbits += len(reps)
        for r, c in zip(reps, cnt):
            items.append((prob * c / total, key, r, int(c)))
    if not items:
        raise ConfigError("no configuration could be generated (too many placements)")
    items.sort(key=lambda it: -it[0])
    tot_raw = sum(it[0] for it in items)
    top = items[:n_rep]
    cov = sum(it[0] for it in top) / tot_raw
    wsum = sum(it[0] for it in top)
    configs = []
    for w, key, rep, osize in top:
        fe2 = _fe2_positions(t, rep, int(round(n2 * z_cell)), E.FE2_SCENARIOS[scen][0])
        st, cs = _build(t, species, mnv, rep, fe2, z_cell)
        oxs = {k: v for k, v in ox.items() if k not in ("Sr", "O")}
        cand = S.StructureCandidate(
            candidate_id=f"{E._formula(comp)}|{'Mn%d+|' % mnv if has_mn else ''}{scen}|cfg{len(configs)}",
            formula=E._formula(comp), sites=cs, oxidation_states=oxs,
            magnetic_hypothesis="collinear Gorter ferrimagnet (W1,W2,W6 down; W3,W4,W5,W7 up)",
            source="hexw.structures.configs (geometry: L02 SrZn2Fe16O27 CIF)", status="candidate",
            partial_occupancy_cif_not_dft_config=False)
        cand.require_dft_configuration()
        E.check_candidate(cand, z_cell)
        cvec = {sp: [int(c) for c in np.array(key)[j]] for j, sp in enumerate(species)} if species else {}
        configs.append(DftConfiguration(cand, st, w / wsum, w / tot_raw, cvec, osize, len(st), notes=[
            "Fe2+ positions: farthest-point placement in scenario %r (placeholder, not enumerated)" % scen,
            "magmom = signed spin-only starting moments (Gorter)",
            "W7: 2d point on the mirror plane instead of the half-occupied split 4f pair (0.16 A shift)"]))
    return ConfigSet(comp, hyp, configs, cov, n_orbits, supercell, notes=[
        "geometry only from the L02 SrZn2Fe16O27 CIF (not relaxed for the chosen composition)",
        "weights: frequency of the integer count vector in the fit ensemble x orbit size / placements (ideal mixing "
        "inside a count vector); they describe the frozen distribution, not the equilibrium at the measuring T",
        f"{n_orbits} symmetry-inequivalent placements in total; the {len(configs)} returned cover {cov:.1%}"])
