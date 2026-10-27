"""Compute coordination of metal sites for the W-type CIFs (pymatgen). Output printed as CSV rows."""
import warnings, sys, csv, os
warnings.filterwarnings("ignore")
import numpy as np
from pymatgen.io.cif import CifParser
LIT = r"D:/Work/Hexx/raw/theory_side/06_literature"
FILES = {
 "L01_Zn": LIT+"/L01_Morch2019/supporting/lt5018sup_cif/SrZn2Fe16O27/SrZn2Fe16O27_SINQ_HRPT.cif",
 "L02_Zn": LIT+"/L02_Morch2023/supporting/CIFs/SrZn2Fe16O27.cif",
 "L02_CoZn": LIT+"/L02_Morch2023/supporting/CIFs/SrCoZnFe16O27.cif",
}
def analyse(tag, path, block_hint):
    p = CifParser(path, occupancy_tolerance=1.0)
    out = []
    structs = p.parse_structures(primitive=False, on_error="ignore") if hasattr(p, "parse_structures") else p.get_structures(primitive=False)
    print(tag, "n structures:", len(structs), [s.composition.reduced_formula for s in structs], file=sys.stderr)
    return structs
if __name__ == "__main__":
    for tag, path in FILES.items():
        structs = analyse(tag, path, None)
        s = structs[0]
        print("==", tag, s.formula, s.lattice.abc)
        # metal sites
        sr = [x for x in s if "Sr" in [sp.symbol for sp in x.species]]
        for i, site in enumerate(s):
            els = [sp.symbol for sp in site.species]
            if "O" in els or "Sr" in els: continue
            nn = s.get_neighbors(site, 2.7)
            d = sorted(n.nn_distance for n in nn if "O" in [sp.symbol for sp in n.species])
            cn = sum(1 for x in d if x < 2.5)
            print(i, dict((sp.symbol, round(o,3)) for sp, o in site.species.items()), np.round(site.frac_coords,4), "CN(<2.5)=", cn, "d=", np.round(d[:7],3))
