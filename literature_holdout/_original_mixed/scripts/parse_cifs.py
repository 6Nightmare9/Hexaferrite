"""Parse W-type CIFs of L01/L02 with pymatgen; write structures/*.csv.
Run: D:/Work/Hexx/.venv/Scripts/python.exe -I literature/scripts/parse_cifs.py"""
import csv, re, glob, os, warnings, json
warnings.filterwarnings("ignore")
from pymatgen.io.cif import CifParser, CifBlock, CifFile
from pymatgen.core import Structure

ROOT = r"D:/Work/Hexx"
LIT = ROOT + "/raw/theory_side/06_literature"
OUT = ROOT + "/literature/structures"

def num(s):
    s = str(s).strip()
    m = re.match(r"^(-?[\d.]+)(?:\((\d+)\))?$", s)
    if not m: return None, None
    v = float(m.group(1)); u = None
    if m.group(2):
        dec = len(m.group(1).split(".")[1]) if "." in m.group(1) else 0
        u = int(m.group(2)) * 10 ** (-dec)
    return v, u

def blocks(path):
    cf = CifFile.from_file(path)
    return cf.data

def asym_rows(b):
    d = b.data
    lab = d.get("_atom_site_label")
    if not lab: return []
    rows = []
    for i, l in enumerate(lab):
        r = {"label": l}
        for k in ("type_symbol", "fract_x", "fract_y", "fract_z", "occupancy", "symmetry_multiplicity", "U_iso_or_equiv", "B_iso_or_equiv"):
            v = d.get("_atom_site_" + k)
            r[k] = v[i] if v else None
        rows.append(r)
    return rows

def cell(b):
    d = b.data
    return {k: d.get("_cell_length_" + k) for k in "abc"} | {"vol": d.get("_cell_volume")}

if __name__ == "__main__":
    out = []
    files = sorted(glob.glob(LIT + "/L01_Morch2019/supporting/lt5018sup_cif/*/*.cif")) + sorted(glob.glob(LIT + "/L02_Morch2023/supporting/CIFs/*.cif"))
    for f in files:
        base = os.path.basename(f)
        if "StructureFactors" in base or base.startswith("pdCIF") or "Magnetic" in base: continue
        try:
            for name, b in blocks(f).items():
                rows = asym_rows(b)
                if not rows: continue
                c = cell(b)
                sg = b.data.get("_symmetry_space_group_name_H-M")
                for r in rows:
                    out.append(dict(file=os.path.relpath(f, LIT).replace("\\", "/"), block=name, sg=sg, a=c["a"], b=c["b"], c=c["c"], **r))
        except Exception as e:
            print("ERR", f, e)
    keys = list(out[0].keys())
    with open(OUT + "/_raw_asym_units.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in out for k in r}))
        w.writeheader(); w.writerows(out)
    print(len(out), "rows")
