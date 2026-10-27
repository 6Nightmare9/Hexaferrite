"""Write structures/site_geometry.csv (pymatgen coordination of the unified sites) and structures/README.md (tables)."""
import warnings, csv, sys, collections
warnings.filterwarnings("ignore")
sys.dont_write_bytecode = True
import numpy as np
from pymatgen.io.cif import CifParser
LIT = r"D:/Work/Hexx/raw/theory_side/06_literature"
OUT = r"D:/Work/Hexx/literature/structures"
CASES = [("L02", "SrZn2Fe16O27 (hex P63/mmc, 300 K)", LIT + "/L02_Morch2023/supporting/CIFs/SrZn2Fe16O27.cif"),
         ("L02", "SrCo2Fe16O27 (Cmcm cell, 300 K)", LIT + "/L02_Morch2023/supporting/CIFs/SrCo2Fe16O27.cif")]
name_by_z = {}
rows = []
for src, tag, path in CASES:
    s = CifParser(path, occupancy_tolerance=1.0).parse_structures(primitive=False, on_error="ignore")[0]
    seen = {}
    for i, site in enumerate(s):
        els = {sp.symbol for sp in site.species}
        if "O" in els or "Sr" in els: continue
        z = site.frac_coords[2] % 0.5 if True else 0
        nn = [n for n in s.get_neighbors(site, 2.7) if "O" in {sp.symbol for sp in n.species}]
        d = sorted(n.nn_distance for n in nn)
        cn = sum(1 for x in d if x < 2.5)
        key = (cn, round(float(np.mean(d[:cn])), 2), round(min(site.frac_coords[2], 1 - site.frac_coords[2]) if cn else 0, 2))
        # classify by CN and z of the site closest to layer positions
        seen.setdefault((cn, tuple(np.round(d[:cn], 3))), []).append(site.frac_coords)
    for (cn, ds), pos in seen.items():
        rows.append(dict(cif=tag, coordination_number=cn, M_O_distances_A=" ".join(f"{x:.3f}" for x in ds), mean_M_O_A=round(float(np.mean(ds)), 3), n_positions_in_cell=len(pos), example_frac_coords=" ".join(f"{x:.4f}" for x in pos[0])))
rows.sort(key=lambda r: (r["cif"], r["coordination_number"], r["mean_M_O_A"]))
with open(OUT + "/site_geometry.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(len(rows), "geometry rows")
