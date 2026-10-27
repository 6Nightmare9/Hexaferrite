"""Non-collinear spin-orbit inputs for the angular dependence of the energy (TZ section 8.6) -- generator only.

NOT FEASIBLE ON THE LAPTOP (dft_feasibility/FEASIBILITY.md, section 3 and 4): a non-collinear SOC SCF of the 92-atom cell costs
x3-4 of the collinear one (about 7-10 h and 6-8 GB per direction on Gamma), and the magnetocrystalline energy (micro-eV..meV)
does not converge on a Gamma-only mesh. Therefore every generated job has the status `not_feasible_on_laptop`; the runner never
starts it, and `soc.allow_run` in the config is false. Also, the pseudopotentials of the campaign are scalar-relativistic
(no fully relativistic files in tools/pseudo): a real SOC run would need new pseudopotentials. The inputs are written for
completeness and for a better machine; they were never executed here, so their syntax is checked only structurally.
No number is produced by this module.
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np

from hexw.dft import qe_input as Q
from hexw.dft import runner as R
from hexw.dft import settings as ST
from hexw.dft.magnet import MagneticSystem

STATUS = "not_feasible_on_laptop"
REFERENCE = "dft_feasibility/FEASIBILITY.md (sections 3 and 4)"
REASON = ("non-collinear SOC SCF of the 92-atom cell: ~7-10 h and 6-8 GB per direction on Gamma, no k-convergence of the "
          "anisotropy energy on a Gamma-only mesh; scalar-relativistic pseudopotentials only")

# polar angle theta from the z axis, azimuth phi (degrees) of the magnetization direction
DIRECTIONS = [(0.0, 0.0), (90.0, 0.0), (90.0, 90.0), (45.0, 0.0), (45.0, 90.0), (90.0, 45.0)]


class SocNotFeasible(RuntimeError):
    pass


def soc_input(system: MagneticSystem, cfg: ST.DftConfig, theta: float, phi: float, *, prefix: str = "soc") -> str:
    """Collinear ground-state input turned into a non-collinear SOC one with the spins of the up sublattices along (theta, phi)
    and those of the down sublattices along (180 - theta, phi + 180)."""
    base = Q.build_input(system, np.ones(len(system.group_names)), cfg, prefix=prefix)
    lines = base.splitlines()
    out = []
    species = []
    for ln in lines:
        m = re.match(r"\s*starting_magnetization\((\d+)\)\s*=\s*(-?[0-9.eE+-]+),?", ln)
        if m:
            i, v = int(m.group(1)), float(m.group(2))
            species.append(i)
            out.append(f"  starting_magnetization({i})={abs(v):g},")
            out.append(f"  angle1({i})={theta if v >= 0 else 180.0 - theta:g}, angle2({i})={phi if v >= 0 else (phi + 180.0) % 360.0:g},")
            continue
        if ln.strip().startswith("nspin=2"):
            out.append("  noncolin=.true., lspinorb=.true.,")
            out.append("  occupations='smearing', smearing='gaussian', degauss=0.005,")
            out.append("  constrained_magnetization='atomic direction', lambda=1.0,")
            continue
        out.append(ln)
    text = "\n".join(out) + "\n"
    # a trailing comma before '/' is legal in a namelist; leave as is
    return "! SOC input, NEVER RUN (not_feasible_on_laptop); see hexw.dft.soc\n" + text


def generate_soc_inputs(system: MagneticSystem, cfg: ST.DftConfig, comp_dir: Path, directions=None) -> dict:
    """Write soc/<dir>/pw.in with status not_feasible_on_laptop. Returns {'state': ..., 'reason': ..., 'jobs': [...]}"""
    comp_dir = Path(comp_dir)
    jobs = []
    for theta, phi in (directions or DIRECTIONS)[: int(cfg.raw["soc"]["n_directions"])]:
        name = f"soc_t{theta:g}_p{phi:g}"
        d = comp_dir / "soc" / name
        d.mkdir(parents=True, exist_ok=True)
        (d / "pw.in").write_text(soc_input(system, cfg, theta, phi, prefix=name), encoding="utf-8")
        R.write_status(d, state=STATUS, message=REASON, reference=REFERENCE, config_sha256=cfg.sha256)
        jobs.append(str(d))
    return {"state": STATUS, "reason": REASON, "reference": REFERENCE, "jobs": jobs,
            "pseudo_fully_relativistic": False}


def soc_status(cfg: ST.DftConfig) -> dict:
    return {"state": STATUS if not cfg.raw["soc"]["allow_run"] else "allowed_by_config", "reason": REASON,
            "reference": REFERENCE}


def run_soc(*_, **__):
    """Always refuses: see the module docstring. There is no override that produces numbers on this machine."""
    raise SocNotFeasible(f"{STATUS}: {REASON}. See {REFERENCE}.")
