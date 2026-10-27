import csv, collections, json
OUT = r"D:/Work/Hexx/literature/structures"
lat = list(csv.DictReader(open(OUT + "/lattice_by_sample.csv", encoding="utf-8")))
occ = list(csv.DictReader(open(OUT + "/site_occupancy_by_sample.csv", encoding="utf-8")))
chk = list(csv.DictReader(open(OUT + "/composition_charge_check.csv", encoding="utf-8")))
geo = list(csv.DictReader(open(OUT + "/site_geometry.csv", encoding="utf-8")))
recs = [json.loads(l) for l in open(r"D:/Work/Hexx/literature/records.jsonl", encoding="utf-8")]
L = []
L.append("# structures/ - сводка по CIF W-фазы (Mørch 2019, Mørch & Christensen 2023)\n")
L.append("Источники: `raw/theory_side/06_literature/L01_Morch2019/supporting/lt5018sup_cif/*` (DOI 10.1107/S2052252519003130) и `L02_Morch2023/supporting/CIFs/*` (DOI 10.1107/S1600576723002133). Разбор: pymatgen (`parse_cifs.py`, `build_site_coordinates.py`), сверка с текстом статей и SI (Table S4, S1-S3, Fig. 3).\n")
L.append("**Предупреждение (ТЗ §7).** Все CIF - частично заселённые усреднённые модели (Fe/Co/Zn/Mg/Ni на одной позиции с дробными заселённостями; W7 = половинная заселённость расщеплённой 4f). Это **не** конфигурации для DFT. Для DFT нужно перечислять симметрийно неэквивалентные упорядочения с целыми числами атомов и зарядовой нейтральностью (см. `composition_charge_check.csv`).\n")
L.append("## Файлы\n")
L.append("| файл | содержание |\n|---|---|\n| `site_coordinates.csv` | все строки _atom_site для TM и Sr W-фазы 7 образцов: координаты, заселённость, **дословная строка CIF**, ячейка, ПГ |\n| `lattice_by_sample.csv` | a, b, c (из CIF) по образцам |\n| `site_occupancy_by_sample.csv` | заселённость по единым позициям W1-W7 (веса по кратности; Cmcm-расщепление объединено) |\n| `composition_charge_check.csv` | состав на ф.е. из CIF и баланс заряда: n(Fe2+) = 2 - n(A) (derived_by_us) |\n| `site_geometry.csv` | КЧ и расстояния M-O по CIF L02 (pymatgen, derived_by_us) |\n| `_raw_asym_units.csv` | сырой разбор всех CIF (отладочный) |\n")
L.append("Не разбирались: `*_StructureFactors.cif` (L01), `pdCIF_*.cif` (L02, 3.4-4.2 MB, профили), `L01 *_Magnetic.cif` (P1-ячейка с повторами; моменты на тип позиции сверены: Zn 4e -3.541, 4f_IV -4.608, 6g 2.411, 4f_VI(S) 3.485, 12k 2.996, 4f_VI*(R) -4.591, 2d 3.530 muB на атом Fe), .mcif L02 (только сверка знаков/величин моментов).\n")
L.append("## Решётки (из CIF; значения статей - в records.yaml, Table S4 / Table 2)\n")
L.append("| источник | образец | ПГ в CIF | a, A | b, A | c, A |\n|---|---|---|---|---|---|")
for r in lat:
    L.append(f"| {r['source']} | {r['sample']} | {r['space_group']} | {r['a_A']} | {r['b_A']} | {r['c_A']} |")
L.append("\nL02: SrCo2W и SrCoZnW в ортогональной ячейке Cmcm (b = sqrt3 a_hex ~ 10.21 A, Z=4); в ячейке P63/mmc Z=2. Температура: L01 - 300 K (mCIF), L02 - комнатная.\n")
sites = ["W1", "W2", "W3", "W4", "W5", "W6", "W7"]
L.append("## Заселённость позиций не-Fe катионами (доля позиции; Fe = остаток)\n")
L.append("| источник | образец | " + " | ".join(sites) + " | всего A / ф.е. |\n|---|---|" + "---|" * 8)
by = collections.defaultdict(dict)
for r in occ:
    if r["species"] != "Fe":
        by[(r["source"], r["sample"])][r["site"]] = by[(r["source"], r["sample"])].get(r["site"], 0) + float(r["occupancy_fraction_of_site"])
tot = {(c["source"], c["sample"]): c["n_nonFe_total"] for c in chk}
for k in by:
    L.append(f"| {k[0]} | {k[1]} | " + " | ".join(f"{by[k].get(s, 0):.3f}" for s in sites) + f" | {tot[k]} |")
L.append("\nСтолбец W7 = 0 (по CIF Fe занимает бипирамиду; A не найден). В SrCoZnW значения - сумма Co+Zn. Сравнение L01 vs L02: см. `../CONFLICTS.md` (A).\n")
L.append("## Баланс заряда и Fe2+ (derived_by_us)\n")
L.append("| источник | образец | n(A) на ф.е. | n(Fe2+) = 2 - n(A) | замечание |\n|---|---|---|---|---|")
for c in chk:
    na = float(c["n_nonFe_total"]); fe2 = float(c["n_Fe2plus_from_charge_balance_2_minus_m"])
    rem = "нарушение: n(A) > 2 -> Fe2+ < 0" if fe2 < 0 else ""
    L.append(f"| {c['source']} | {c['sample']} | {na:.3f} | {fe2:.3f} | {rem} |")
L.append("\nМодель Sr2+, A2+, O27, только Fe2+/Fe3+ (условная). Состав из CIF = все позиции полностью заняты Fe+A (Fe + A = 18). Локализация Fe2+ по позициям в CIF **не задана** (везде Fe+3).\n")
L.append("## Координационные полиэдры (L02, 300 K; pymatgen)\n")
L.append("| CIF | КЧ | расстояния M-O, A | <M-O>, A | позиций в ячейке | пример дроб. коорд. |\n|---|---|---|---|---|---|")
for g in geo:
    L.append(f"| {g['cif']} | {g['coordination_number']} | {g['M_O_distances_A']} | {g['mean_M_O_A']} | {g['n_positions_in_cell']} | {g['example_frac_coords']} |")
L.append("\nОтнесение к W1-W7 см. `../site_notation.md` (по КЧ и z).\n")
# moments
L.append("## Моменты по позициям (muB)\n")
L.append("L01, Fig. 4 (усреднено по заселению; положительный = параллельно суммарному моменту; Co - планарный порядок):\n")
L.append("| позиция | Mg | Co | Ni | Zn |\n|---|---|---|---|---|")
d = collections.defaultdict(dict)
for r in recs:
    if r["lit_id"] == "L01" and r["name"].startswith("site_magnetic_moment_W"):
        d[r["name"][-2:]][r["ext"]["formula"][2:4]] = r["ext"]["value_original"]
for w in sites:
    if w in d:
        L.append(f"| {w} | " + " | ".join(d[w].get(m, "") for m in ("Mg", "Co", "Ni", "Zn")) + " |")
L.append("\nL02 (моменты Fe3+, Tables S1-S3; Co и Zn-строки в SI; компоненты My/Mz для планарных, Mz для SrZn2W) - см. `records.yaml` (`Fe_moment_*`). Определения и знаки различаются между L01 и L02 (CONFLICTS.md B.4-B.5).\n")
open(OUT + "/README.md", "w", encoding="utf-8").write("\n".join(L))
print("ok")
