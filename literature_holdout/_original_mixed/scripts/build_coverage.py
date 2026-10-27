"""coverage_map.csv: composition x property x conditions x source, with gaps and the action that closes them.
x = total A (Co/Zn/Mn) atoms per formula unit in SrFe(18-x)A(x)O27 (x = 2 <-> SrA2Fe16O27). NOTE: L02 uses x for the Zn fraction in SrCo(2-x)Zn(x)Fe16O27."""
import csv

OUT = r"D:/Work/Hexx/literature/coverage_map.csv"

PROPS = {  # property -> default gap action (RU, as requested)
 "lattice a,c": "нужен DFT (релаксация) ; для проверки позже - лабораторный РФА (недоступен как вход)",
 "cation distribution over 7 sites (occupancy)": "нужен DFT (энергии предпочтения позиций, Fe2+ локализация); нейтроны для Sr-Mn/Zn/Co x<2 - нужен эксперимент",
 "magnetic structure (axial/planar/cone), T_SR": "нужен DFT (MAE(T)/конфигурации) + нужен эксперимент (нейтроны при T)",
 "site moments": "нужен DFT",
 "Ms(T) (mass or volume)": "нужен DFT+модель; оцифровка M(T) если рисунок",
 "Curie temperature": "нужен DFT (J_ij -> MC/MFA)",
 "K1,K2 / HA (experimental)": "нужен эксперимент (монокристалл/торк) или оцифровка; для Sr-серии лит. данных нет",
 "MAE / K DFT": "нужен DFT (SOC, чувствительность к U)",
 "exchange J_ij (DFT)": "нужен DFT (своя конвенция H, знак, нормировка)",
 "g-factor, alpha, FMR linewidth": "нужен эксперимент (в корпусе нет ни одного значения для Sr W)",
 "natural FMR / mu'(f), mu''(f)": "нужна оцифровка (Ba-аналоги) / нужен эксперимент (Sr)",
 "eps'(f), eps''(f), sigma": "нужна оцифровка (L11 Fig. 6, Ba) / нужен эксперимент (Sr); в корпусе нет eps Sr W",
 "Mossbauer B_hf, IS, QS per site (literature)": "нужен эксперимент (литература по SrW закрыта: C027 1988; для Ba - C067,C074,C077 закрыты); корпус пуст",
 "synthesis window / phase purity / Fe2+ content": "нужен эксперимент (лит.) ; Mn - нет данных вообще",
}

COMPS = [  # (class, label, A content x)
 ("Sr end-member", "SrFe18O27 (x=0)", 0),
 ("Sr-Co", "SrFe17.5Co0.5O27 (x=0.5)", 0.5), ("Sr-Co", "SrFe17Co1O27 (x=1)", 1), ("Sr-Co", "SrFe16.5Co1.5O27 (x=1.5)", 1.5), ("Sr-Co", "SrCo2Fe16O27 (x=2)", 2),
 ("Sr-Zn", "SrFe17.5Zn0.5O27 (x=0.5)", 0.5), ("Sr-Zn", "SrFe17Zn1O27 (x=1)", 1), ("Sr-Zn", "SrFe16.5Zn1.5O27 (x=1.5)", 1.5), ("Sr-Zn", "SrZn2Fe16O27 (x=2)", 2),
 ("Sr-Mn", "SrFe17.5Mn0.5O27 (x=0.5)", 0.5), ("Sr-Mn", "SrFe17Mn1O27 (x=1)", 1), ("Sr-Mn", "SrFe16.5Mn1.5O27 (x=1.5)", 1.5), ("Sr-Mn", "SrMn2Fe16O27 (x=2)", 2),
 ("Sr mixed Co-Zn (x_total=2)", "SrCo1.5Zn0.5Fe16O27", 2), ("Sr mixed Co-Zn (x_total=2)", "SrCoZnFe16O27", 2), ("Sr mixed Co-Zn (x_total=2)", "SrCo0.5Zn1.5Fe16O27", 2),
 ("Sr mixed with Mn", "SrCo/Zn/Mn mixed (any)", None),
 ("Sr analogue (tier 2)", "SrMg2Fe16O27", 2), ("Sr analogue (tier 2)", "SrNi2Fe16O27", 2),
]

# covered cells: (composition label prefix or exact, property) -> (sources, conditions, tier, split, comment)
C = {}


def cov(comp, prop, src, cond, tier, split, comment="", status="covered"):
    C[(comp, prop)] = (status, src, cond, tier, split, comment)


E = "SrFe18O27 (x=0)"
cov(E, "lattice a,c", "L05 (DFT a=5.877, c=32.842 A, GGA relaxed); L05 cites exp. a 5.89-5.91 (secondary)", "DFT 0 K", 2, "literature_fit", "calculated only; experimental lattice from primary refs [50,51] not in corpus", "partial")
cov(E, "magnetic structure (axial/planar/cone), T_SR", "L05 (collinear Gorter + uniaxial MAE>0, DFT)", "DFT 0 K", 2, "literature_fit", "no experiment in corpus", "partial")
cov(E, "site moments", "L05 Table 2 (7 sites, spin+orbital)", "DFT GGA+U(3 eV)+SOC, 0 K", 2, "literature_fit", "")
cov(E, "Ms(T) (mass or volume)", "C004 Table 2 (81.0-81.1 A m2/kg, 3 syntheses 1448-1523 K); L05 mu0Ms=0.652 T (DFT)", "RT, 1/H^2 extrapolation; DFT 0 K", 2, "literature_fit", "Ms(T) curve absent; one RT point (+ ~2-3 wt % Fe3O4)")
cov(E, "MAE / K DFT", "L05 Table 3: Ku=0.148 MJ/m3, mu0Ha=0.569 T, kappa=0.661", "DFT 0 K, U=3 eV", 2, "literature_fit", "exp. Ku=0.344 MJ/m3 (300 K) only as secondary citation (Ueda) - primary needed")
cov(E, "synthesis window / phase purity / Fe2+ content", "C004 (sealed tube, T>=1448 K, 97-98 wt % W); L03 context", "1448-1523 K, 12 h, pO2~1e-8 atm (sealed)", 2, "literature_fit", "")
cov(E, "exchange J_ij (DFT)", "none (L05 defers J; L04 PRB 2024 not accessible)", "", 2, "", "L04 full text not obtained", "gap")
for x in ("SrCo2Fe16O27 (x=2)",):
    cov(x, "lattice a,c", "L01 Table S4 (5.8946(1), 32.7411(7)); L02 Table 2 (5.89575(2)/b 10.21042/32.76195); C015 composite (5.849/33.112, anomalous)", "RT, powder, 1300 C (L01) / 1200 C (L02)", 1, "L02 fit; L01 holdout; C015 holdout (low weight)", "L01 vs L02 differ in synthesis T and refinement software")
    cov(x, "cation distribution over 7 sites (occupancy)", "L01 CIF/Fig.3 (Co: 6g 33 %, 4e 14 %, 4f_IV 11 %, ...); L02 Table S1 (Co on 6g 47.7 % weighted, 4e 7 %)", "RT, NPD+XRD", 1, "L02 fit; L01 holdout", "ALTERNATIVES kept: see CONFLICTS (6g 33 % vs 48 %)")
    cov(x, "magnetic structure (axial/planar/cone), T_SR", "L01 (planar, 90 deg, R-factor scan); L02 (Cm'cm' planar + non-collinear canting); L02 TT=555 K (deriv.)/520 K (onset); C023 ~470 K; C015 cites 550 K (cone)", "RT diffraction; M(T) 100 Oe", 1, "L02 fit; C023 holdout", "transition nature (cone vs uniaxial) NOT resolved by diffraction")
    cov(x, "site moments", "L01 Fig.4 (avg per site); L02 Table S1 (Fe & Co moments, My/Mz)", "RT", 1, "L02 fit; L01 holdout", "Fe:Co moment ratio fixed 5:3 in both")
    cov(x, "Ms(T) (mass or volume)", "L01 Ms,VSM=78.5(1) A m2/kg (300 K, +-3 T); C044 cites 74.31 (secondary)", "300 K", 1, "L01 holdout", "no Ms(T); L02 M(T) only normalised (figure_only)", "partial")
    cov(x, "Curie temperature", "L01 TC=470(3) C = 743 K (TG); L02 790 K (M(T) 100 Oe, derivative) / abstract 780 K; C015 cites ~850 K", "powder; different methods", 1, "L02 fit; L01 holdout", "47 K spread between L01 and L02: see CONFLICTS")
    cov(x, "synthesis window / phase purity / Fe2+ content", "L03 (pO2 0.2-387 atm; n*(Co)=2.03 at 387 atm; flux crystal W-phase Co 1.79); L02 Table 1 (1473 K 2 h); L01 (1300 C 2 h); C023, C015", "1200-1300 C", 1, "L03 fit", "")
    cov(x, "K1,K2 / HA (experimental)", "none for Sr (Ba: C057 K1(Co0.62Zn1.38 Ba)=8.7e4 J/m3; L11 Ha(BaCo2W) formal)", "", 1, "", "Sr Co2W is easy-plane: K1<0; only DFT can fill", "gap")
    cov(x, "natural FMR / mu'(f), mu''(f)", "Ba analogue only: L12 (Co1.2Zn0.8 plane / Co0.6Zn1.4 axis, figures), C056/C036 (Ba, composite)", "RT", 3, "holdout", "Sr: none", "gap")
x = "SrZn2Fe16O27 (x=2)"
cov(x, "lattice a,c", "L01 Table S4 (5.90922(3), 32.8293(2)); L02 Table 2 (5.908047(10), 32.83467(6))", "RT", 1, "L02 fit; L01 holdout", "")
cov(x, "cation distribution over 7 sites (occupancy)", "L01 CIF (Zn: 4e 55 %, 4f_IV 24 %, 6g 0); L02 CIF (Zn: 4e 58 %, 4f_IV 14.5 %, 6g 9.7 %)", "RT, NPD+XRD", 1, "L02 fit; L01 holdout", "Zn on 4f_IV 24 vs 14.5 %; 6g 0 vs 9.7 % (same nominal composition, 1473 K 2 h both) - CONFLICTS")
cov(x, "magnetic structure (axial/planar/cone), T_SR", "L01, L02: uniaxial P63/mm'c', no transition up to TC", "RT; M(T) 300-900 K", 1, "L02 fit; L01 holdout", "")
cov(x, "site moments", "L01 Fig.4; L02 Table S3 (Fe moments Mz)", "RT", 1, "L02 fit; L01 holdout", "L02 Table S3 W4 sign vs mCIF")
cov(x, "Ms(T) (mass or volume)", "L01 Ms,VSM=76.2(1) raw, 80.4(3) impurity-corrected A m2/kg (300 K)", "300 K", 1, "L01 holdout", "", "partial")
cov(x, "Curie temperature", "L01 345(3) C = 618 K; L02 590 K", "TG vs VSM 50 Oe", 1, "L02 fit; L01 holdout", "28 K spread")
cov(x, "synthesis window / phase purity / Fe2+ content", "L03 (spinel impurity grows at 387 atm); L01 (ZnFe2O4 5.28 wt %); L02 (spinel 4.12, M 1.77, X 2.67 wt %)", "1200 C", 1, "L03 fit", "")
for lab in ("SrCo1.5Zn0.5Fe16O27", "SrCo0.5Zn1.5Fe16O27"):
    cov(lab, "Curie temperature", "L02 Table 3 (720 K / 650 K)", "M(T) 50 Oe heating 10 K/min", 1, "L02 fit", "no diffraction for these two compositions")
    cov(lab, "magnetic structure (axial/planar/cone), T_SR", "L02 Table 3 TT = 450 K / 265 K (derivative max)", "M(T)", 1, "L02 fit", "phase assignment not diffraction-verified", "partial")
cov("SrCoZnFe16O27", "lattice a,c", "L02 Table 2 (5.89874(3)/10.21562(5)/32.79282(9))", "RT", 1, "L02 fit", "")
cov("SrCoZnFe16O27", "cation distribution over 7 sites (occupancy)", "L02 Table S2 (all Co on 6g sub-sites: 8d 19 %, 4b 45 % -> 28 % of 6g; Zn on 4e 30 %, 4f_IV 10 %)", "RT", 1, "L02 fit", "")
cov("SrCoZnFe16O27", "magnetic structure (axial/planar/cone), T_SR", "L02 Cm'cm' planar; TT 405 K (deriv.)/360 K (onset)", "RT; M(T)", 1, "L02 fit", "")
cov("SrCoZnFe16O27", "site moments", "L02 Table S2", "RT", 1, "L02 fit", "")
cov("SrCoZnFe16O27", "Curie temperature", "L02 Table 3: 700 K (abstract 680 K)", "50 Oe", 1, "L02 fit", "")
for lab, comp in (("SrMg2Fe16O27", "Mg"), ("SrNi2Fe16O27", "Ni")):
    for p, s in (("lattice a,c", "L01 Table S4"), ("magnetic structure (axial/planar/cone), T_SR", "L01 (uniaxial)"), ("site moments", "L01 Fig. 4"), ("Curie temperature", "L01 Fig. 4")):
        cov(lab, p, s, "RT", 2, "L01 fit", "Ni: site occupancy not unambiguous (Fe/Ni contrast) per L01" if comp == "Ni" and p.startswith("site") else "")
    cov(lab, "cation distribution over 7 sites (occupancy)", "L01 CIF/Fig. 3", "RT", 2, "L01 fit", "Ni poorly determined" if comp == "Ni" else "")
    cov(lab, "Ms(T) (mass or volume)", "L01 Ms,VSM", "300 K", 2, "L01 fit", "")
# Mn explicit
for lab in [c[1] for c in COMPS if c[0] == "Sr-Mn"] + ["SrCo/Zn/Mn mixed (any)"]:
    for p in PROPS:
        pass
cov("SrMn2Fe16O27 (x=2)", "synthesis window / phase purity / Fe2+ content", "none for Sr; Ba: L11 (BaMn2W single-phase 1300 C 3 h, Ms 47.77 emu/g, Ha 11301 Oe); C056 (Ba, Mn0.4Zn0.4Co1.2); C053 (low-quality Ca/Sr-Mn-Zn)", "Ba analogue only", 3, "L11 holdout", "Mn2+ vs Mn3+ and Mn site preference unknown; Sr W with Mn not in corpus; C026 (Ce-Mn SrW) closed access", "gap")
cov("SrMn2Fe16O27 (x=2)", "Ms(T) (mass or volume)", "Ba analogue: L11 BaMn2W Ms=47.77 emu/g (300 K) - much lower than Co/Ni/Zn Ba analogues", "300 K", 3, "L11 holdout", "transfer Ba->Sr unjustified; Mn moment assumed 2.7 muB by authors (inconsistent with high-spin Mn2+)", "partial")
cov("SrMn2Fe16O27 (x=2)", "K1,K2 / HA (experimental)", "Ba analogue: L11 BaMn2W Ha=11301.48 Oe (approach-to-saturation, formal)", "300 K", 3, "L11 holdout", "", "partial")
cov("SrMn2Fe16O27 (x=2)", "natural FMR / mu'(f), mu''(f)", "Ba analogue (composite) C056: Mn0.4Zn0.4Co1.2 mu'' peaks 4.5 and 8 GHz", "RT, powder-paraffin 8:3", 3, "C056 holdout", "Mn only as 1/5 of A content; not a Sr", "partial")

# generic fills: Zn/Co intermediate x ; Ba analogue info
cov("SrFe17.5Co0.5O27 (x=0.5)", "synthesis window / phase purity / Fe2+ content", "L03 formalism: actual composition SrMe(2-d)Fe(16+d) (Fe2+); no samples with x<1.7 for Sr", "", 1, "", "x<2 is the user's main grid; only DFT/new experiment", "gap")
# Ba analogues of intermediate Co-Zn compositions (tier 3)
for lab in ("SrCo1.5Zn0.5Fe16O27", "SrCoZnFe16O27", "SrCo0.5Zn1.5Fe16O27"):
    cov(lab, "K1,K2 / HA (experimental)", "Ba analogue C057: BaCo0.62Zn1.38 K1=8.7e4 J/m3 (296 K) ; C060 SRT vs x ; L12 x=1.3 boundary", "", 3, "holdout", "Sr data absent", "partial")
    cov(lab, "natural FMR / mu'(f), mu''(f)", "Ba analogue L12 (figures)", "RT", 3, "L12 holdout", "", "partial")

rows = []
for cls, lab, x in COMPS:
    for p, act in PROPS.items():
        key = (lab, p)
        if key in C:
            st, src, cond, tier, split, cm = C[key]
            rows.append(dict(composition_class=cls, composition=lab, x_total_A=x, property=p, conditions=cond, status=st, sources=src, transferability_tier=tier, proposed_split=split, gap_action=("" if st == "covered" else act), comment=cm))
        else:
            comment = ""
            extra = ""
            if p.startswith("Mossbauer"):
                extra = "в корпусе нет литературных мёссбауэровских данных для W-фаз"
            if cls == "Sr-Mn" or cls == "Sr mixed with Mn":
                extra = (extra + "; " if extra else "") + "Sr-Mn-W: данных нет ни по структуре, ни по магнетизму"
            rows.append(dict(composition_class=cls, composition=lab, x_total_A=x, property=p, conditions="", status="gap", sources="", transferability_tier="", proposed_split="", gap_action=act, comment=extra))
# global gaps
for p, src, cm in (("g-factor, alpha, FMR linewidth", "none in corpus (Sr or Ba)", "L12 gives only model + figures; C056/C036 only mu'' peak positions of composites"),
                   ("Mossbauer B_hf, IS, QS per site (literature)", "none (L07 Rancourt/T01/T02 methodological only; primary Sr/Ba W Mossbauer papers closed-access)", "")):
    rows.append(dict(composition_class="ALL", composition="all W compositions", x_total_A="", property=p, conditions="", status="gap", sources=src, transferability_tier="", proposed_split="", gap_action=PROPS[p], comment=cm))
with open(OUT, "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
import collections
print(len(rows), collections.Counter(r["status"] for r in rows))
