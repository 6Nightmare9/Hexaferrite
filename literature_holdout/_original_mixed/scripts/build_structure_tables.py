"""Build structures/*.csv from the raw asymmetric-unit rows extracted by parse_cifs.py.
Unified site numbering W1..W7 follows the label order of Morch 2019/2023 and Islam 2025 (see site_notation.md).
Run: D:/Work/Hexx/.venv/Scripts/python.exe -I literature/scripts/build_structure_tables.py"""
import csv, re, collections, os
OUT = r"D:/Work/Hexx/literature/structures"
rows = list(csv.DictReader(open(OUT + "/_raw_asym_units.csv", encoding="utf-8")))

def num(s):
    s = str(s).strip()
    m = re.match(r"^(-?[\d.]+)(?:\((\d+)\))?$", s)
    if not m: return None, None
    v = float(m.group(1)); u = None
    if m.group(2):
        dec = len(m.group(1).split(".")[1]) if "." in m.group(1) else 0
        u = int(m.group(2)) * 10 ** (-dec)
    return v, u

# hex-cell (Z=2) multiplicities of the 7 TM sites, atoms per formula unit
MULT_HEX = {1: 4, 2: 4, 3: 6, 4: 4, 5: 12, 6: 4, 7: 4}   # W7: 4 positions (4f), half occupied -> 2 atoms
SITE_NAME = {1: "4e tet(S)", 2: "4f_IV tet(S)", 3: "6g oct(S-S)", 4: "4f_VI oct(S)", 5: "12k oct(R-S)", 6: "4f_VI* oct(R)", 7: "2d bipyr(R) [4f half-occ.]"}
SPIN = {1: "down", 2: "down", 3: "up", 4: "up", 5: "up", 6: "down", 7: "up"}

def site_index(label):
    m = re.match(r"^(?:Fe|Co|Zn|Mg|Ni)(\d)", label)
    return int(m.group(1)) if m else None

SAMPLES = [  # (source, sample, file-suffix, block-name, hex?)
 ("L01", "SrMg2Fe16O27", "SrMg2Fe16O27_SINQ_HRPT.cif", "SrMg2Fe16O27"),
 ("L01", "SrCo2Fe16O27", "SrCo2Fe16O27_SINQ_HRPT.cif", "SrCo2Fe16O27"),
 ("L01", "SrNi2Fe16O27", "SrNi2Fe16O27_SINQ_HRPT.cif", "SrNi2Fe16O27"),
 ("L01", "SrZn2Fe16O27", "SrZn2Fe16O27_SINQ_HRPT.cif", "SrZn2Fe16O27"),
 ("L02", "SrCo2Fe16O27", "CIFs/SrCo2Fe16O27.cif", ""),
 ("L02", "SrCoZnFe16O27", "CIFs/SrCoZnFe16O27.cif", ""),
 ("L02", "SrZn2Fe16O27", "CIFs/SrZn2Fe16O27.cif", ""),
]
long_rows, lat_rows, comp_rows = [], [], []
for src, sample, suffix, block in SAMPLES:
    sel = [r for r in rows if r["file"].endswith(suffix) and r["block"] == block]
    assert sel, (src, sample)
    ortho = src == "L02" and sample != "SrZn2Fe16O27" or (src == "L02" and sel[0]["sg"] == "63")
    r0 = sel[0]
    lat_rows.append(dict(source=src, sample=sample, cif=r0["file"], space_group=r0["sg"], a_A=r0["a"], b_A=r0["b"], c_A=r0["c"]))
    # collect per site/species
    acc = collections.defaultdict(lambda: collections.defaultdict(float))     # site -> species -> weighted occ (atoms in cell)
    acc_u = collections.defaultdict(lambda: collections.defaultdict(float))
    wsum = collections.defaultdict(float)
    for r in sel:
        sp = re.match(r"([A-Z][a-z]?)", r["type_symbol"]).group(1)
        if sp in ("Sr", "O"): continue
        idx = site_index(r["label"])
        occ, u = num(r["occupancy"])
        if ortho:
            mult = {"8": 8, "16": 16, "4": 4}[str(r["symmetry_multiplicity"]).strip()] if r["symmetry_multiplicity"] else None
        else:
            mult = None
        if idx is None: continue
        # weights: hex -> MULT_HEX (position multiplicity); ortho -> given multiplicity (8f,8d,4b,16h)
        w = (mult if ortho else MULT_HEX[idx])
        acc[idx][sp] += w * occ
        acc_u[idx][sp] += (w * (u or 0)) ** 2
        if sp == "Fe": wsum[idx] += w * 1.0   # positions counted once per Fe row
    for idx in range(1, 8):
        # total positions in the cell for this site
        if ortho:
            tot = {1: 8, 2: 8, 3: 12, 4: 8, 5: 24, 6: 8, 7: 8}[idx]     # Cmcm cell has Z=4 -> 8/8/12/8/24/8/8 positions
            cell_f_units = 4
        else:
            tot = MULT_HEX[idx]; cell_f_units = 2
        for spc, a in acc[idx].items():
            frac = a / tot
            fu = a / cell_f_units
            unc = (acc_u[idx][spc] ** 0.5) / tot if acc_u[idx][spc] else None
            long_rows.append(dict(source=src, sample=sample, cif=r0["file"], site=f"W{idx}", site_name=SITE_NAME[idx], spin_gorter=SPIN[idx],
                                  species=spc, occupancy_fraction_of_site=round(frac, 4), occupancy_unc_approx=(round(unc, 4) if unc else ""), atoms_per_formula_unit=round(fu, 4)))
for fn, data in (("lattice_by_sample.csv", lat_rows), ("site_occupancy_by_sample.csv", long_rows)):
    with open(os.path.join(OUT, fn), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(data[0].keys())); w.writeheader(); w.writerows(data)

# composition & charge check (derived_by_us)
tab = collections.defaultdict(lambda: collections.defaultdict(float))
for r in long_rows:
    tab[(r["source"], r["sample"])][r["species"]] += r["atoms_per_formula_unit"]
out = []
for (src, sample), d in tab.items():
    m = sum(v for k, v in d.items() if k != "Fe")
    tot = sum(d.values())
    out.append(dict(source=src, sample=sample, **{f"n_{k}_per_fu": round(v, 3) for k, v in d.items()}, n_TM_total=round(tot, 3),
                    n_nonFe_total=round(m, 3), n_Fe2plus_from_charge_balance_2_minus_m=round(2 - m, 3)))
keys = sorted({k for o in out for k in o})
with open(os.path.join(OUT, "composition_charge_check.csv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(out)
for o in out: print(o)
