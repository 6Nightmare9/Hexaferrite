"""Parsing of a pw.x output (Quantum ESPRESSO 7.x) and the check that the run stayed in the target magnetic configuration."""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field, asdict
from typing import Optional

import numpy as np

from hexw.dft.hamiltonian import RY_TO_MEV


@dataclass
class PwResult:
    qe_version: Optional[str] = None
    job_done: bool = False
    converged: bool = False
    n_iterations: Optional[int] = None
    energy_Ry: Optional[float] = None
    energy_meV: Optional[float] = None
    smearing_contrib_Ry: Optional[float] = None
    estimated_scf_accuracy_Ry: Optional[float] = None
    total_magnetization: Optional[float] = None
    absolute_magnetization: Optional[float] = None
    fermi_eV: Optional[float] = None
    site_moments: list = field(default_factory=list)         # integrated on the atomic sphere (muB), by atom order
    site_charges: list = field(default_factory=list)
    hubbard_d_moments: dict = field(default_factory=dict)    # atom -> (up - down) of Tr[ns]
    hubbard_ns: dict = field(default_factory=dict)           # atom (1-based) -> {"up", "down", "total"}: full Tr[ns] of the LAST SCF step
    hubbard_ns_start: dict = field(default_factory=dict)     # the same for the FIRST printed block (starting occupations, iteration 0)
    nat: Optional[int] = None
    n_electrons: Optional[float] = None
    n_kpoints: Optional[int] = None
    cpu_s: Optional[float] = None
    wall_s: Optional[float] = None
    ram_estimate_MB: Optional[float] = None
    errors: list = field(default_factory=list)
    iteration_energies_Ry: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def _hms(txt: str) -> float:
    t = 0.0
    for val, unit in re.findall(r"([0-9.]+)\s*([dhms])", txt):
        t += float(val) * {"d": 86400, "h": 3600, "m": 60, "s": 1}[unit]
    return t


def parse_pw_output(text: str) -> PwResult:
    r = PwResult()
    m = re.search(r"Program PWSCF\s+v\.?\s*([0-9][0-9.\w]*)", text)
    r.qe_version = m.group(1) if m else None
    r.job_done = "JOB DONE." in text
    r.converged = "convergence has been achieved" in text
    m = re.search(r"convergence has been achieved in\s+(\d+)\s+iterations", text)
    if m:
        r.n_iterations = int(m.group(1))
    elif "convergence NOT achieved" in text:
        r.n_iterations = len(re.findall(r"iteration #\s*\d+", text))
    bang = re.findall(r"^!\s+total energy\s+=\s+(-?[0-9.]+)\s+Ry", text, re.M)
    if bang and r.converged:
        r.energy_Ry = float(bang[-1])
        r.energy_meV = r.energy_Ry * RY_TO_MEV
    r.iteration_energies_Ry = [float(x) for x in re.findall(r"^\s*total energy\s+=\s+(-?[0-9.]+)\s+Ry", text, re.M)]
    m = re.findall(r"smearing contrib\. \(-TS\)\s+=\s+(-?[0-9.]+)\s+Ry", text)
    if m:
        r.smearing_contrib_Ry = float(m[-1])
    m = re.findall(r"estimated scf accuracy\s+<\s+([0-9.Ee+-]+)\s+Ry", text)
    if m:
        r.estimated_scf_accuracy_Ry = float(m[-1])
    m = re.findall(r"total magnetization\s+=\s+(-?[0-9.]+)\s+Bohr", text)
    if m:
        r.total_magnetization = float(m[-1])
    m = re.findall(r"absolute magnetization\s+=\s+(-?[0-9.]+)\s+Bohr", text)
    if m:
        r.absolute_magnetization = float(m[-1])
    m = re.findall(r"the Fermi energy is\s+(-?[0-9.]+)\s+ev", text, re.I)
    if m:
        r.fermi_eV = float(m[-1])
    # moments per site: last block
    blocks = text.split("Magnetic moment per site")
    if len(blocks) > 1:
        last = blocks[-1]
        for mm in re.finditer(r"atom\s+(\d+)\s+\(R=[0-9.]+\)\s+charge=\s*(-?[0-9.]+)\s+magn=\s*(-?[0-9.]+)", last):
            r.site_charges.append(float(mm.group(2)))
            r.site_moments.append(float(mm.group(3)))
    ns, full, first = {}, {}, {}
    for mm in re.finditer(r"Tr\[ns\(\s*(\d+)\)\]\s+\(up, down, total\)\s+=\s+(-?[0-9.]+)\s+(-?[0-9.]+)\s+(-?[0-9.]+)", text):
        a = int(mm.group(1))
        up, dn, tot = float(mm.group(2)), float(mm.group(3)), float(mm.group(4))
        ns[a] = up - dn
        full[a] = {"up": up, "down": dn, "total": tot}                    # later blocks overwrite: the last one is the final state
        first.setdefault(a, dict(full[a]))
    r.hubbard_d_moments = ns
    r.hubbard_ns = full
    r.hubbard_ns_start = first
    m = re.search(r"number of atoms/cell\s+=\s+(\d+)", text)
    r.nat = int(m.group(1)) if m else None
    m = re.search(r"number of electrons\s+=\s+([0-9.]+)", text)
    r.n_electrons = float(m.group(1)) if m else None
    m = re.search(r"number of k points=\s*(\d+)", text)
    r.n_kpoints = int(m.group(1)) if m else None
    m = re.search(r"PWSCF\s*:\s*(.*?)\s+CPU\s+(.*?)\s+WALL", text)
    if m:
        r.cpu_s, r.wall_s = _hms(m.group(1)), _hms(m.group(2))
    m = re.search(r"Estimated max dynamical RAM per process >\s+([0-9.]+)\s+(MB|GB)", text)
    if m:
        r.ram_estimate_MB = float(m.group(1)) * (1024.0 if m.group(2) == "GB" else 1.0)
    for mm in re.finditer(r"Error in routine\s+(\S+)\s+\((\d+)\):\s*\n\s*(.+)", text):
        r.errors.append(f"{mm.group(1)}: {mm.group(3).strip()}")
    return r


def parse_pw_file(path) -> PwResult:
    with open(path, encoding="utf-8", errors="replace") as f:
        return parse_pw_output(f.read())


# ---- is the target configuration preserved? --------------------------------------------------------------------------
def check_configuration(result: PwResult, target_sign: np.ndarray, spin: np.ndarray, *, min_fraction: float,
                        start_moment: np.ndarray, require_sign: bool = True) -> dict:
    """Compare the sphere moments with the target configuration. An atom with spin S > 0 must keep the target sign and
    |m| >= min_fraction * start_moment (start_moment: spin-only ionic moment 2S, muB). Returns {ok, reason, flipped,
    collapsed, ...}. A run without per-site moments cannot be verified: ok=False ('no_site_moments')."""
    n = len(spin)
    if len(result.site_moments) != n:
        return {"ok": False, "reason": "no_site_moments", "flipped": [], "collapsed": [], "n_site_moments": len(result.site_moments)}
    m = np.array(result.site_moments)
    mag = np.flatnonzero(spin > 0)
    flipped, collapsed = [], []
    for k in mag:
        if require_sign and m[k] * target_sign[k] <= 0:
            flipped.append(int(k) + 1)
        elif abs(m[k]) < min_fraction * start_moment[k]:
            collapsed.append(int(k) + 1)
    ok = not flipped and not collapsed
    reason = "preserved" if ok else ("flipped atoms %s" % flipped[:8] if flipped else "") + (
        "; collapsed atoms %s" % collapsed[:8] if collapsed else "")
    return {"ok": ok, "reason": reason.strip("; "), "flipped": flipped, "collapsed": collapsed,
            "min_abs_moment": float(np.min(np.abs(m[mag]))) if len(mag) else None,
            "mean_abs_moment": float(np.mean(np.abs(m[mag]))) if len(mag) else None,
            "total_target_spin_moment": float(np.sum(target_sign[mag] * start_moment[mag]))}


# ---- is the valence state preserved? (Fe2+ location, total moment, smearing, Mn2+) ---------------------------------------
VALENCE_DEFAULTS = {"total_magnetization_tol_muB": 0.5, "max_abs_TS_meV": 1.0, "mn_min_fraction_of_fe3_moment": 0.9}


def _ion(a: dict) -> str:
    """'Fe2+' / 'Fe3+' / 'Mn2+' ... of an atom record of design.json (older records: spin 2 on Fe = Fe2+)."""
    if a.get("ion"):
        return str(a["ion"])
    if a["species"] == "Fe":
        return "Fe2+" if abs(float(a["spin_S"]) - 2.0) < 1e-6 else "Fe3+"
    return str(a["species"])


def check_valence(result: PwResult, atoms: list, target_sign: np.ndarray, opts: Optional[dict] = None) -> dict:
    """Did the SCF keep the assigned charge state? A collinear DFT+U run can move the extra electron of Fe2+ to another Fe
    (or spread it), which changes the magnetic energy by far more than the exchange. Checks:
      (a) the N Fe atoms with the LARGEST d occupation Tr[ns] (total) are the N assigned Fe2+ atoms; without Tr[ns] the N Fe
          atoms with the SMALLEST |m| (sphere moment) must be them;
      (b) |M_tot - M_expected| <= tol, M_expected = sum(sign * 2S) over the magnetic atoms (muB per cell);
      (c) |TS| (smearing contribution) < max_abs_TS_meV; a non-zero TS means fractional occupations at the Fermi level;
      (d) Mn assigned as Mn2+ must keep a sphere moment >= mn_min_fraction_of_fe3_moment x the mean |m| of the Fe3+ atoms
          (Mn3+ d4 would be ~0.8 of Mn2+/Fe3+ d5; the fraction is an ASSUMPTION, between 0.8 and 1).
    Returns {ok, valence_ok, smearing_ok, reasons, ...}."""
    o = dict(VALENCE_DEFAULTS, **(opts or {}))
    reasons, out = [], {}
    n = len(atoms)
    ion = [_ion(a) for a in atoms]
    spin = np.array([float(a["spin_S"]) for a in atoms])
    fe = [k for k in range(n) if atoms[k]["species"] == "Fe" and spin[k] > 0]
    fe2 = sorted(k for k in fe if ion[k] == "Fe2+")
    out["fe2_assigned"] = [k + 1 for k in fe2]
    m = np.array(result.site_moments) if len(result.site_moments) == n else None
    ns = result.hubbard_ns
    if fe2:
        basis, found = None, None
        if ns and all((k + 1) in ns for k in fe):
            nd = {k: ns[k + 1]["total"] for k in fe}
            order = sorted(fe, key=lambda k: -nd[k])
            found, basis = sorted(order[:len(fe2)]), "largest Tr[ns] (n_d)"
            rest = [nd[k] for k in fe if k not in fe2]
            out["nd_gap"] = float(min(nd[k] for k in fe2) - max(rest)) if rest else None
            out["nd_fe2_mean"] = float(np.mean([nd[k] for k in fe2]))
            out["nd_fe3_mean"] = float(np.mean(rest)) if rest else None
        elif m is not None:
            order = sorted(fe, key=lambda k: abs(m[k]))
            found, basis = sorted(order[:len(fe2)]), "smallest |m|"
        if found is None:
            reasons.append("Fe2+ location cannot be checked (no Tr[ns] and no site moments)")
        else:
            out["fe2_found"], out["fe2_basis"] = [k + 1 for k in found], basis
            if found != fe2:
                reasons.append(f"Fe2+ moved: assigned atoms {[k + 1 for k in fe2]}, {basis} gives {[k + 1 for k in found]}")
        if m is not None:
            out["abs_m_fe2_mean"] = float(np.mean([abs(m[k]) for k in fe2]))
    mag = np.flatnonzero(spin > 0)
    m_exp = float(np.sum(target_sign[mag] * 2.0 * spin[mag]))
    out["M_expected_muB"] = m_exp
    out["M_total_muB"] = result.total_magnetization
    if result.total_magnetization is None:
        reasons.append("no total magnetization in the output")
    else:
        out["dM_muB"] = float(result.total_magnetization - m_exp)
        if abs(out["dM_muB"]) > o["total_magnetization_tol_muB"]:
            reasons.append(f"total magnetization {result.total_magnetization:+.2f} muB, expected {m_exp:+.2f} "
                           f"(|diff| {abs(out['dM_muB']):.2f} > {o['total_magnetization_tol_muB']:g})")
    mn = [k for k in range(n) if atoms[k]["species"] == "Mn" and spin[k] > 0 and ion[k] == "Mn2+"]
    if mn and m is not None:
        fe3 = [abs(m[k]) for k in fe if ion[k] == "Fe3+"]
        if fe3:
            ref = float(np.mean(fe3))
            low = [k + 1 for k in mn if abs(m[k]) < o["mn_min_fraction_of_fe3_moment"] * ref]
            out["mn_abs_m"] = [float(abs(m[k])) for k in mn]
            out["fe3_mean_abs_m"] = ref
            if low:
                reasons.append(f"Mn atoms {low} have |m| < {o['mn_min_fraction_of_fe3_moment']:g} x mean Fe3+ moment "
                               f"({ref:.2f} muB): Mn2+ not preserved")
    elif mn:
        reasons.append("Mn2+ cannot be checked (no site moments)")
    valence_ok = not reasons
    ts = None if result.smearing_contrib_Ry is None else abs(result.smearing_contrib_Ry) * RY_TO_MEV
    out["TS_meV"] = ts
    smearing_ok = ts is not None and ts < o["max_abs_TS_meV"]
    if ts is None:
        out["smearing_reason"] = "no smearing contribution in the output"
    elif not smearing_ok:
        out["smearing_reason"] = f"|TS| = {ts:.3g} meV >= {o['max_abs_TS_meV']:g} meV (fractional occupations at the Fermi level)"
    out.update(ok=valence_ok and smearing_ok, valence_ok=valence_ok, smearing_ok=smearing_ok, reasons=reasons)
    return out
