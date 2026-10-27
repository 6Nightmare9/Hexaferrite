"""alpha-Fe2O3 (hematite, R-3c) cell for the smoke test of the chain generate -> run -> collect -> fit.

Geometry: space group R-3c (hexagonal setting a = 5.038 A, c = 13.772 A, Fe 12c z = 0.3553, O 18e x = 0.3059), the
parameters of Blake et al., Am. Mineral. 51, 123 (1966), as used in the resource pilot `dft_feasibility/pilots/make_a.py`
(that pilot is not a derivative of a holdout file: it builds the cell from these numbers, no CIF). PROVENANCE CAVEAT: the
numbers were taken from the pilot script and from memory of the cited paper; they were not re-verified against the paper
in this run. The cell is the 10-atom primitive rhombohedral cell with four Fe atoms. The Fe atoms are ordered along the
trigonal axis and the antiferromagnetic ground state is +,-,-,+ (Fe1 Fe2 Fe3 Fe4 along c), which is the known
easy-plane/c-axis collinear order of hematite.

This is a smoke system for the CHAIN, not a result about hematite: pilot cutoffs, small k mesh, U = 3 eV.
"""
from __future__ import annotations

import numpy as np

from hexw.dft.magnet import MagneticSystem
from hexw.dft import hamiltonian as H

GEOMETRY = {"space_group": "R-3c", "a_A": 5.038, "c_A": 13.772, "Fe_z": 0.3553, "O_x": 0.3059,
            "source": "Blake et al., Am. Mineral. 51 (1966) 123, via dft_feasibility/pilots/make_a.py (not re-verified here)"}


def hematite_structure():
    from pymatgen.core import Lattice, Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    g = GEOMETRY
    s = Structure.from_spacegroup(g["space_group"], Lattice.hexagonal(g["a_A"], g["c_A"]), ["Fe", "O"],
                                  [[0, 0, g["Fe_z"]], [g["O_x"], 0, 0.25]])
    p = SpacegroupAnalyzer(s, symprec=1e-3).get_primitive_standard_structure()
    return p


def hematite_system() -> MagneticSystem:
    p = hematite_structure()
    if len(p) != 10:
        raise ValueError(f"unexpected hematite primitive cell: {len(p)} atoms")
    p.add_oxidation_state_by_element({"Fe": 3, "O": -2})
    fe = [i for i, x in enumerate(p) if x.specie.symbol == "Fe"]
    axis = p.lattice.matrix.sum(0)
    axis = axis / np.linalg.norm(axis)                      # the trigonal axis of the rhombohedral cell
    fe.sort(key=lambda i: p[i].coords @ axis)
    spin = np.zeros(len(p)); group = -np.ones(len(p), int); base = np.zeros(len(p))
    signs = [1, -1, -1, 1]
    for g, (i, sg) in enumerate(zip(fe, signs)):
        spin[i] = H.spin_length("Fe", 3)
        group[i] = g
        base[i] = sg
    return MagneticSystem("alpha-Fe2O3", "hematite", p, spin, group, [f"Fe{k + 1}" for k in range(4)], base,
                          meta={"geometry": dict(GEOMETRY), "order_axis": "trigonal axis, Fe sorted along it",
                                "ground_state": "+--+"})
