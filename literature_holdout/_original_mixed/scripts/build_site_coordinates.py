"""structures/site_coordinates.csv : every TM asymmetric-unit row of the W-phase block of the L01/L02 CIFs with the verbatim CIF line.
Run: D:/Work/Hexx/.venv/Scripts/python.exe -I literature/scripts/build_site_coordinates.py"""
import csv, re, os
LIT = r"D:/Work/Hexx/raw/theory_side/06_literature"
OUT = r"D:/Work/Hexx/literature/structures"
DOI = {"L01": "10.1107/S2052252519003130", "L02": "10.1107/S1600576723002133"}
SITE_NAME = {1: "4e tet(S)", 2: "4f_IV tet(S)", 3: "6g oct(S-S)", 4: "4f_VI oct(S)", 5: "12k oct(R-S)", 6: "4f_VI* oct(R)", 7: "2d bipyr(R) [4f half-occ.]"}
rows = list(csv.DictReader(open(OUT + "/_raw_asym_units.csv", encoding="utf-8")))
SEL = [("L01", "SrMg2Fe16O27", "SrMg2Fe16O27_SINQ_HRPT.cif", "SrMg2Fe16O27"), ("L01", "SrCo2Fe16O27", "SrCo2Fe16O27_SINQ_HRPT.cif", "SrCo2Fe16O27"),
       ("L01", "SrNi2Fe16O27", "SrNi2Fe16O27_SINQ_HRPT.cif", "SrNi2Fe16O27"), ("L01", "SrZn2Fe16O27", "SrZn2Fe16O27_SINQ_HRPT.cif", "SrZn2Fe16O27"),
       ("L02", "SrCo2Fe16O27", "CIFs/SrCo2Fe16O27.cif", ""), ("L02", "SrCoZnFe16O27", "CIFs/SrCoZnFe16O27.cif", ""), ("L02", "SrZn2Fe16O27", "CIFs/SrZn2Fe16O27.cif", "")]
out = []
for src, sample, suf, block in SEL:
    sel = [r for r in rows if r["file"].endswith(suf) and r["block"] == block]
    path = LIT + "/" + sel[0]["file"]
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    # restrict to the W-phase data block region: first occurrence only (first line matching label)
    for r in sel:
        sp = re.match(r"([A-Z][a-z]?)", r["type_symbol"]).group(1)
        m = re.match(r"^(?:Fe|Co|Zn|Mg|Ni)(\d)", r["label"])
        if sp in ("O",): continue
        idx = int(m.group(1)) if m else None
        if sp != "Sr" and idx is None: continue
        pat = re.compile(r"^\s*" + re.escape(r["label"]) + r"\s")
        hit = next((l.strip() for l in lines if pat.match(l)), "")
        out.append(dict(source=src, doi=DOI[src], cif_file=sel[0]["file"], sample_id_in_paper=sample, site=("Sr" if sp == "Sr" else f"W{idx}"),
                        site_name=("Sr (2b/4c)" if sp == "Sr" else SITE_NAME[idx]), cif_label=r["label"], species=sp, x=r["fract_x"], y=r["fract_y"], z=r["fract_z"],
                        occupancy=r["occupancy"], coordinate_unit="fractional (cell edges)", verbatim_cif_line=hit, space_group=r["sg"],
                        cell_a_A=r["a"], cell_b_A=r["b"], cell_c_A=r["c"], status="measured (Rietveld, X-ray+neutron combined)",
                        note="partially occupied CIF row; NOT a DFT configuration"))
with open(OUT + "/site_coordinates.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
print(len(out), "rows"); print(sum(1 for o in out if not o["verbatim_cif_line"]), "without verbatim line")
