"""Records from L01 (Morch 2019) and L02 (Morch & Christensen 2023)."""
import csv
from rc import add, vu

OUT = r"D:/Work/Hexx/literature/structures"
coords = list(csv.DictReader(open(OUT + "/site_coordinates.csv", encoding="utf-8")))

# ------------------------------------------------------------------ L01
SYN1 = ("sol-gel autocombustion: Sr(NO3)2, Me(NO3)2.6H2O, Fe(NO3)3.9H2O (Sigma-Aldrich technical, >98 %), [Sr]:[Me]:[Fe]=1:2:16, citric acid, "
        "NH4OH to neutral, dried 100 C, gel fired 350 C 30 min, crushed, then fired 1200 C (Mg, Zn) or 1300 C (Ni, Co), 2 h hold; "
        "furnace atmosphere and cooling rate not stated in text (heating scheme is a figure not reproduced in the XML)")
METH1 = "combined Rietveld refinement (FullProf) of 4 data sets: SLS-MS synchrotron XRD (0.778 A), Rigaku SmartLab Co-Ka XRD, SINQ-HRPT neutrons (1.89 and 2.45 A)"
CL1 = {"Mg": "MORCH_L01_MgNi", "Ni": "MORCH_L01_MgNi", "Co": "MORCH_L01_CoZn", "Zn": "MORCH_L01_CoZn"}
TIER1 = {"Mg": 2, "Ni": 2, "Co": 1, "Zn": 1}
SYN1T = {"Mg": 1200, "Zn": 1200, "Ni": 1300, "Co": 1300}
SAMP = {m: f"Sr{m}2Fe16O27" for m in ("Mg", "Co", "Ni", "Zn")}


def l01(name, me, **kw):
    kw.setdefault("sample", SAMP[me]); kw.setdefault("formula", SAMP[me]); kw.setdefault("cluster", CL1[me]); kw.setdefault("tier", TIER1[me])
    kw.setdefault("synth", SYN1 + f"; this sample fired at {SYN1T[me]} C"); kw.setdefault("method", METH1)
    kw.setdefault("group", f"L01_NPD_XRD_{SAMP[me]}")
    kw.setdefault("form", "powder"); kw.setdefault("texture", "none (powder)")
    return add("L01", name, **kw)


# lattice (Table S4)
A = {"Mg": "5.8955(1)", "Co": "5.8946(1)", "Ni": "5.8877(3)", "Zn": "5.90922(3)"}
C = {"Mg": "32.7308(6)", "Co": "32.7411(7)", "Ni": "32.702(2)", "Zn": "32.8293(2)"}
for me in ("Mg", "Co", "Ni", "Zn"):
    l01("lattice_a", me, vs=A[me], unit="angstrom", kind="length", status="measured", locator="Supporting information sup1.pdf, Table S4 (refined unit-cell parameters)",
        verbatim="a = b [A]   5.8955(1)   5.8946(1)   5.8877(3)   5.90922(3)  (columns SrMg2Fe16O27, SrCo2Fe16O27, SrNi2Fe16O27, SrZn2Fe16O27)", T={"K": 300, "note": "300 K from mCIF header _magn_diffrn_temperature; diffraction T not stated in main text"},
        definition="hexagonal a=b, space group P63/mmc", unc_note="(n) = standard uncertainty in last digit as printed (IUCr convention; not defined explicitly in the paper)")
    l01("lattice_c", me, vs=C[me], unit="angstrom", kind="length", status="measured", locator="Supporting information sup1.pdf, Table S4",
        verbatim="c [A]   32.7308(6)   32.7411(7)   32.702(2)   32.8293(2)", T={"K": 300, "note": "see lattice_a"}, definition="hexagonal c, space group P63/mmc",
        unc_note="(n) = standard uncertainty in last digit as printed")

# site occupancies from CIFs (non-Fe species), partially occupied CIF rows
for r in coords:
    if r["source"] != "L01" or r["species"] in ("Fe", "Sr"):
        continue
    me = r["species"]
    l01(f"site_occupancy_{me}_{r['site']}", me, vs=r["occupancy"], unit="1", kind="fraction", status="measured",
        locator=f"CIF {r['cif_file'].split('/')[-1]} (supplementary lt5018sup2.zip), _atom_site loop, row {r['cif_label']} ({r['site_name']})",
        verbatim=r["verbatim_cif_line"], definition=f"fraction of the {r['site_name']} crystallographic site occupied by {me} (all sites constrained fully occupied by Fe+{me}); partially occupied CIF, NOT a DFT configuration",
        T={"K": 300, "note": "RT; see lattice_a"}, unc_note="(n) as printed in CIF")

# Fig. 3 totals and the Zn cells that disagree with the CIF/Table S3
for me, tv in (("Mg", "1.69(1)"), ("Co", "2.2(1)"), ("Ni", "2.1(1)"), ("Zn", "1.9(1)")):
    l01(f"total_nonFe_content_per_fu_{me}", me, vs=tv, unit="1", kind="dimensionless", status="measured", locator="Fig. 3 (table embedded in figure), row 'Total Me' (SrMe2Fe16O27)",
        verbatim=f"Total Me  {tv}", definition="refined total Me per formula unit SrMe2Fe16O27 (nominal 2)", T={"K": 300, "note": "RT"})
l01("site_occupancy_Zn_W2_FIG3_alt", "Zn", vs="10.5(8)", unit="%", kind="fraction", status="measured", locator="Fig. 3 (table embedded in figure), row 4f_tet(S), column Zn (%)",
    verbatim="4f tet(S)   0   0   0   ... Zn (%) = 10.5(8)", definition="Zn on 4f_IV tet(S) as % of full site occupancy (Fig. 3)",
    note="CONFLICT inside L01: CIF row Zn2TS occupancy 0.236(9) and Table S3 0.236(9) give 23.6 %; Fig. 3 gives 10.5(8) %. Sum rule (total Zn 1.90) supports 23.6 %.",
    transfer="same sample; alternatives kept, not averaged")
l01("site_occupancy_Zn_W6_FIG3_alt", "Zn", vs="0", unit="%", kind="fraction", status="measured", locator="Fig. 3 (table embedded in figure), row 4f_oct(R), column Zn (%)",
    verbatim="4f oct(R)   ...  Zn (%) = 0", definition="Zn on 4f_VI* oct(R) as % of full site",
    note="CONFLICT inside L01: CIF Zn6OR occupancy 0.032(8) (3.2 %) vs Fig. 3 value 0.")

# magnetic moments: Fig. 4 table (average over occupants, incl. diamagnetic ions = 0)
FIG4 = {  # site -> Mg, Co, Ni, Zn
 "W7": ("3.6(4)", "2.4(1)", "3.7(5)", "3.5(6)"), "W4": ("3.6(3)", "3.2(1)", "3.8(4)", "3.4(5)"), "W3": ("2.0(6)", "2.9(1)", "3.1(5)", "2.4(5)"),
 "W1": ("-4.5(4)", "-4.2(1)", "-5.0(5)", "-1.6(1.3)"), "W2": ("-4.2(3)", "-3.3(1)", "-4.3(3)", "-3.5(6)"), "W6": ("-4.7(3)", "-3.8(1)", "-4.6(5)", "-4.4(4)"),
 "W5": ("3.5(2)", "3.5(1)", "3.5(3)", "2.9(3)")}
LAB = {"W7": "4f_1/2biV(R)", "W4": "4f_oct(S)", "W3": "6g_oct(S-S)", "W1": "4e_tet(S)", "W2": "4f_tet(S)", "W6": "4f_oct(R)", "W5": "12k_oct(R-S)"}
for w, vals in FIG4.items():
    for me, v in zip(("Mg", "Co", "Ni", "Zn"), vals):
        v2 = v.replace("(1.3)", "(13)") if "(1.3)" in v else v
        val, u = vu(v2) if "(13)" not in v2 else (float(v2.split("(")[0]), 1.3)
        l01(f"site_magnetic_moment_{w}", me, value=val, unc=u, unit="muB", kind="magnetic_moment", proposed_kind="magnetic_moment", status="measured",
            locator="Fig. 4 (table embedded in figure) row " + LAB[w],
            verbatim=f"{LAB[w]}  {v} (muB) [column {SAMP[me]}]", value_original=v,
            definition="average magnetic moment per site from joint X-ray/neutron refinement (NPD), averaged over site occupants (Fe and diamagnetic Zn/Mg counted); sign: positive = parallel to net moment, collinear model, sample-wide Gorter ferrimagnet (Co: planar, moments in ab plane)",
            T={"K": 300, "note": "RT"}, unc_note="(n) as printed; for Zn 4e '(1.3)' is 1.3 muB",
            note="moment constraint in refinement: mu(Fe3+):mu(Co2+)=5:3 (unquenched orbital moment ignored)" if me == "Co" else None, pending_unit="muB")
MS = {"Mg": ("41(3)", "60.3(1)", "63.7(2)", "435(2)"), "Co": ("59(3)", "78.5(1)", "78.4(2)", "470(3)"), "Ni": ("50(5)", "67.7(1)", "69.0(3)", "500(2)"), "Zn": ("57(5)", "76.2(1)", "80.4(3)", "345(3)")}
for me, (npd, vsm, corr, tc) in MS.items():
    l01("Ms_from_NPD_moments", me, vs=npd, unit="A*m^2/kg", kind="specific_magnetization", proposed_kind="specific_magnetization", status="calculated",
        locator="Fig. 4 (table) row Ms,NPD", verbatim=f"Ms,NPD (Am2/kg)  {npd}", definition="sum of refined site moments over the unit cell divided by cell mass (authors); neutron moments are lower than VSM for all samples",
        T={"K": 300, "note": "RT"}, pending_unit="A*m^2/kg", transformation="none (1 A m2/kg = 1 emu/g)")
    l01("Ms_VSM_raw", me, vs=vsm, unit="A*m^2/kg", kind="specific_magnetization", proposed_kind="specific_magnetization", status="measured",
        locator="Fig. 4 (table) row Ms,VSM; also Table S9/S10 (Supporting information)", verbatim=f"Ms,VSM (Am2/kg)  {vsm}",
        definition="Ms by approach to saturation (extrapolation), after graphical near-infinite-slope demagnetisation correction; includes impurity phases (spinel ~5-8 wt %)",
        T={"K": 300}, field="up to +-3 T, PPMS-VSM", form="powder 20-45 mg in cylindrical capsule (h 2.7-3.0 mm, d 3.5 mm)", texture="none", pending_unit="A*m^2/kg", transformation="none (1 A m2/kg = 1 emu/g)")
    l01("Ms_VSM_impurity_corrected", me, vs=corr, unit="A*m^2/kg", kind="specific_magnetization", proposed_kind="specific_magnetization", status="assumed_by_authors",
        locator="Fig. 4 (table) row 'Wt% Corrected Ms'; Supporting information S8, Table S10", verbatim=f"Wt% Corrected Ms (Am2/kg)  {corr}",
        definition="Ms,VSM corrected for impurity phases using refined wt% and LITERATURE Ms of the impurities (MgFe2O4 25, CoFe2O4 80, NiFe2O4 50, ZnFe2O4 1, X-type 63.29 A m2/kg)",
        T={"K": 300}, note="depends on assumed impurity Ms from other papers; conclusion section calls these 'Ms,VSM' (naming ambiguity)", pending_unit="A*m^2/kg")
    l01("Curie_temperature", me, value=float(vu(tc)[0]) + 273.15, unc=vu(tc)[1], unit="K", kind="temperature", status="measured", locator="Fig. 4 (table) row T_C (degC); Conclusions; Supporting S2 / Fig. S4",
        verbatim=f"T_C (C)  {tc}", value_original=f"{tc} degC", unit_original="degC", transformation="degC -> K: +273.15",
        method="thermogravimetry in external magnet (Netzsch STA 449 F3), pellets 60-120 mg, 10 C/min, Ar 50 ml/min; T_C = intercept of linear fits of the mass increase near T_C and the preceding mass signal",
        definition="Curie temperature, TG-magnet intercept method", T={"K": None}, field="external permanent magnet (strength not given)", form="pellet 60-120 mg", unc_note="(n) as printed, in degC")

IMP = {"Co": ("CoFe2O4", "5.27(5)"), "Mg": ("MgFe2O4", "8.24(5)"), "Ni": ("NiFe2O4", "7.08(6)"), "Zn": ("ZnFe2O4", "5.28(4)")}
for me, (ph, w) in IMP.items():
    l01(f"impurity_wt_{ph}", me, vs=w, unit="%", kind="fraction", status="measured", locator="Section 3.2 Purity and occupation, para 1",
        verbatim=f"small impurities of CoFe2O4 [5.27 (5) wt%], MgFe2O4 [8.24 (5) wt%], NiFe2O4 [7.08 (6) wt%] and ZnFe2O4 [5.28 (4) wt%] -> {ph} {w} wt%",
        definition="weight fraction of spinel impurity from SLS-MS Rietveld", transformation="wt % kept as %")
l01("impurity_wt_Sr2Mg2Fe28O46_Xtype", "Mg", vs="8.77(5)", unit="%", kind="fraction", status="measured", locator="Section 3.2 Purity and occupation, para 1",
    verbatim="The SrMg2Fe16O27 sample additionally contained a small impurity of X-type hexaferrite Sr2Mg2Fe28O46 [8.77 (5) wt%].", definition="X-type wt fraction")

# magnetic structure
for me in ("Mg", "Ni", "Zn"):
    l01("moment_angle_from_c_axis", me, value=0.0, unit="deg", kind="angle", status="measured", locator="Section 3.3 Magnetic structure; Fig. 1(b)",
        verbatim="For Mg, Ni and Zn there is no magnetic contribution to the (00l) reflections, indicating that the magnetic moment is aligned along the c axis.",
        definition="collinear ferrimagnet k=0, uniaxial easy axis c; magnetic space group P63/mm'c' (Supporting S5, OG 194.8.1501)", T={"K": 300},
        qualitative="uniaxial", note="no uncertainty given")
l01("moment_angle_from_c_axis", "Co", value=90.0, unit="deg", kind="angle", status="measured", locator="Section 3.3 Magnetic structure; Fig. 6 (R factors vs angle, 5 deg steps); Fig. 1(c)",
    verbatim="the decreases in R factors are small between conical and planar ordering, the lowest value is for 90 degrees and there is no evidence to suggest that the ordering should be conical rather than planar. SrCo2Fe16O27 must belong to a lower symmetry space group than P63/mmc",
    definition="angle between refined collinear moment direction and c; scan in 5 deg steps 0-90; planar ordering in ab plane at 300 K", T={"K": 300}, qualitative="planar",
    unc_note="no uncertainty given (scan step 5 deg; R-factor minimum shallow between conical and planar)",
    note="L01 says earlier BaCo2Fe16O27 work reported conical (69-71 deg) or planar order; Co content 1.87(7) vs 2.2(1) per f.u. suggested as reason")
l01("magnetic_structure_Co_Fig6_Rfactor_vs_angle", "Co", status="measured", figure_only=True, locator="Fig. 6",
    axes="x: angle of moment from c axis, 0-90 deg (5 deg steps); y: R factors (HRPT data)", verbatim="Figure 6: Resulting R factors from the refinements where the only difference is the constrained angle for the magnetic ordering with respect to the c axis.",
    definition="R-factor profile", note="needed only if conical-vs-planar question is re-opened; no digitisation done")
l01("hysteresis_M_H_300K", "Co", status="measured", figure_only=True, locator="Fig. 7 (all four samples)", axes="x: applied field H (+-3 T); y: mass magnetization (A m2/kg); VSM 300 K, powders",
    verbatim="Figure 7: Magnetic hysteresis for the four SrMe2Fe16O27 (Me = Mg, Co, Ni, Zn) samples. The samples were measured using an applied field of +-3 T with a PPMS-VSM.",
    definition="M(H) loops; Co sample approaches saturation differently (planar)", note="Hc between 5 and 35 kA/m and (BH)max 0.5-2.2 kJ/m3 for the series (text 3.4)")

# ------------------------------------------------------------------ L02
SYN2 = ("sol-gel autocombustion: Sr(NO3)2, Fe(NO3)3.9H2O, Co(NO3)2.6H2O, Zn(NO3)2.6H2O (Sigma-Aldrich ACS >=98 %), stoichiometric, citric acid (1:1 to nitrates), "
        "NH4OH neutralised, dried 373 K, gel fired 623 K 30 min, cooled in air, crushed, then calcined (Table 1)")
CALC = {"SrCo2Fe16O27": "1473 K, 2 h", "SrCo1.5Zn0.5Fe16O27": "1523 K, 4 h", "SrCoZnFe16O27": "1523 K, 2 h", "SrCo0.5Zn1.5Fe16O27": "1523 K, 1 h", "SrZn2Fe16O27": "1473 K, 2 h"}
METH2 = "combined Rietveld (TOPAS 6): SuperHRPD J-PARC neutron TOF (3 banks, 8 h) + Diamond I11 synchrotron XRD (0.826555 A, 1 h); magnetic symmetry-mode analysis (ISODISTORT)"


def l02(name, formula, **kw):
    kw.setdefault("sample", formula); kw.setdefault("formula", formula); kw.setdefault("cluster", "MORCH_L02"); kw.setdefault("tier", 1)
    kw.setdefault("synth", SYN2 + (f"; calcination {CALC[formula]} (atmosphere not stated)" if formula in CALC else "; calcination per Table 1 (series)"))
    kw.setdefault("method", METH2); kw.setdefault("group", f"L02_NPD_XRD_{formula}"); kw.setdefault("form", "powder"); kw.setdefault("texture", "none (powder)")
    return add("L02", name, **kw)


for f in CALC:
    h, m = CALC[f].split(", ")
    T = float(h.split()[0]); t = float(m.split()[0])
    l02("calcination_temperature", f, value=T, unit="K", kind="temperature", status="measured", locator="Table 1 (Calcination temperatures and times for the five compositions studied)",
        verbatim="Temperature (K) and time (h): 1473 K, 2 h | 1523 K, 4 h | 1523 K, 2 h | 1523 K, 1 h | 1473 K, 2 h (SrCo2, SrCo1.5Zn0.5, SrCoZn, SrCo0.5Zn1.5, SrZn2)",
        group=f"L02_synth_{f}", note="optimal sintering T 'found through multiple synthesis attempts'; diffraction only for x=0,1,2")
    l02("calcination_time", f, value=t, unit="h", kind="time", status="measured", locator="Table 1", verbatim=f"{CALC[f]}", group=f"L02_synth_{f}")

# lattice (Table 2): refined for the 3 diffraction samples
LAT2 = {"SrCo2Fe16O27": ("5.89575(2)", "10.21042(4)", "32.76195(7)", "Cm'cm' orthorhombic (Cmcm cell)"),
        "SrCoZnFe16O27": ("5.89874(3)", "10.21562(5)", "32.79282(9)", "Cm'cm' orthorhombic (Cmcm cell)"),
        "SrZn2Fe16O27": ("5.908047(10)", None, "32.83467(6)", "P63/mm'c' hexagonal")}
for f, (a, b, c, sg) in LAT2.items():
    for nm, v in (("lattice_a", a), ("lattice_b_orthorhombic", b), ("lattice_c", c)):
        if v is None:
            continue
        l02(nm, f, vs=v, unit="angstrom", kind="length", status="measured", locator="Table 2 (agreement factors and crystallographic parameters), rows a, b, c",
            verbatim=f"a (A) 5.89575 (2) | 5.89874 (3) | 5.908047 (10); b (A) 10.21042 (4) | 10.21562 (5) | 5.908047 (10)§; c (A) 32.76195 (7) | 32.79282 (9) | 32.83467 (6) -> {nm} = {v}",
            definition=f"{sg}; orthorhombic b = sqrt(3) a_hex (footnote: 5.908047 x 3^1/2 = 10.233038 A)", T={"K": 300, "note": "room temperature"},
            note="CIF cell values (e.g. SrCo2Fe16O27.cif a=5.89486(2)) differ slightly: SI S1 says lattice parameters were refined separately for neutron and X-ray data (temperature differences)")

COMP = {"SrCo2Fe16O27": [("Co", "1.71(3)")], "SrCoZnFe16O27": [("Co", "0.84(9)"), ("Zn", "0.81(13)")], "SrZn2Fe16O27": [("Zn", "1.799(11)")]}
for f, lst in COMP.items():
    for el, v in lst:
        l02(f"refined_content_{el}_per_fu", f, vs=v, unit="1", kind="dimensionless", status="measured", locator="Section 3.3 Occupancy and purity, last sentence",
            verbatim="The resulting compositions from the refinements are SrCo1.71 (3)Fe16.29 (3)O27, SrCo0.84 (9)Zn0.81 (13)Fe16.35 (9)O27 and SrZn1.799 (11)Fe16.201 (11)O27.",
            definition=f"refined {el} atoms per formula unit (all sites full; Fe+Co+Zn = 18)", note="Fe2+ not refined; charge balance in the 2-x model requires n(Fe2+) = 2 - n(non-Fe) (derived_by_us in structures/composition_charge_check.csv)")

# Table 3 / thermomagnetics
TC = {"SrCo2Fe16O27": ("790", "555"), "SrCo1.5Zn0.5Fe16O27": ("720", "450"), "SrCoZnFe16O27": ("700", "405"), "SrCo0.5Zn1.5Fe16O27": ("650", "265"), "SrZn2Fe16O27": ("590", None)}
THM = dict(method="isofield M(T), PPMS-VSM with oven, heating 300-900 K at 10 K/min; powder dispersed in zircar cement (QDS-4097-030) on Pt heater, Cu foil; pre-magnetised at 10 000 Oe",
           orient="powder, random", form="powder nugget in cement")
for f, (tc, tt) in TC.items():
    fld = "100 Oe (7.98 kA/m)" if f == "SrCo2Fe16O27" else "50 Oe (3.99 kA/m)"
    l02("Curie_temperature_Table3", f, value=float(tc), unit="K", kind="temperature", status="measured", locator="Table 3 (Curie temperatures and suggested transition temperatures, based on the maximum of the derivative)",
        verbatim="TC (K) 790 720 700 650 590 ; TT (K) 555 450 405 265 - (SrCo2, SrCo1.5Zn0.5, SrCoZn, SrCo0.5Zn1.5, SrZn2)", definition="T_C from maximum of |dM/dT| (derivative) on heating", field=fld, T={"K": None, "range": "300-900 K"},
        group=f"L02_MT_{f}", **THM, unc_note="no uncertainty given")
    if tt:
        l02("transition_temperature_Table3_derivative_max", f, value=float(tt), unit="K", kind="temperature", status="measured", locator="Table 3",
            verbatim="TT (K) 555 450 405 265 - ; 'suggested transition temperatures, based on the maximum of the derivative'", definition="planar(Cm'cm') -> higher-T uniaxial/other transition T from derivative maximum; nature of high-T phase (conical vs uniaxial) NOT determined by diffraction",
            field=fld, T={"K": None, "range": "300-900 K"}, group=f"L02_MT_{f}", **THM, unc_note="no uncertainty given")
l02("transition_temperature_onset_text", "SrCo2Fe16O27", value=520.0, unit="K", kind="temperature", status="measured", locator="Section 3.5 and Fig. 4 caption; Abstract",
    verbatim="This feature is observed around 520 K for SrCo2Fe16O27 and 360 K for SrCoZnFe16O27 but is not present in SrZn2Fe16O27.", definition="onset of magnetic transition (peak in derivative lower than TC) as read in the text",
    field="100 Oe", group="L02_MT_SrCo2Fe16O27", **THM, note="CONFLICT within L02: Table 3 derivative-max gives 555 K; text gives 520 K (different definition)")
l02("transition_temperature_onset_text", "SrCoZnFe16O27", value=360.0, unit="K", kind="temperature", status="measured", locator="Section 3.5 and Fig. 4 caption; Abstract",
    verbatim="This feature is observed around 520 K for SrCo2Fe16O27 and 360 K for SrCoZnFe16O27", definition="onset (text)", field="50 Oe", group="L02_MT_SrCoZnFe16O27", **THM,
    note="CONFLICT within L02: Table 3 gives 405 K")
l02("Curie_temperature_abstract", "SrCo2Fe16O27", value=780.0, unit="K", kind="temperature", status="measured", locator="Abstract",
    verbatim="magnetic transitions at 520 and 360 K for SrCo2Fe16O27 and SrCoZnFe16O27, and Curie temperatures of 780 and 680 K, respectively, while SrZn2Fe16O27 showed no transition but a Curie temperature at 590 K",
    definition="TC as in abstract", field="100 Oe", group="L02_MT_SrCo2Fe16O27", **THM, note="CONFLICT within L02: Table 3 gives 790 K")
l02("Curie_temperature_abstract", "SrCoZnFe16O27", value=680.0, unit="K", kind="temperature", status="measured", locator="Abstract",
    verbatim="Curie temperatures of 780 and 680 K, respectively", definition="TC as in abstract", field="50 Oe", group="L02_MT_SrCoZnFe16O27", **THM, note="CONFLICT within L02: Table 3 gives 700 K")
l02("thermomagnetic_M_T_300_900K", "SrCo2Fe16O27", status="measured", figure_only=True, locator="Fig. 4 (x = 0, 1, 2), Fig. S6 (x = 0-2 in steps of 0.5; a: 10-400 K, b: 325-875 K)",
    axes="x: temperature (325-875 K shown; 300-900 K measured; Fig. S6a 10-400 K); y: M/M(300 K) and its derivative, normalised", verbatim="Figure 4: Thermomagnetic measurements of SrCo2Fe16O27, SrCoZnFe16O27 and SrZn2Fe16O27. The applied fields were, respectively, 100, 50 and 50 Oe.",
    definition="isofield M(T)", group="L02_MT_all", **THM, note="digitisation not done; low-T (10-400 K) data in Fig. S6a could constrain low-T anisotropy changes")
l02("composition_of_planar_to_uniaxial_boundary_at_RT", "SrCo2-xZnxFe16O27", value=0.65, unit="1", kind="dimensionless", status="assumed_by_authors", locator="Section 3.5, last paragraph; Fig. 5 (tentative phase diagram)",
    verbatim="A sample with a composition around SrCo0.65Zn1.35Fe16O27 would have a magnetic phase transition between Cm'cm' and P63/mm'c' at room temperature",
    definition="Co content (per f.u.) at which authors extrapolate the Cm'cm'/P63/mm'c' boundary to 300 K (extrapolation from 3 neutron compositions + M(T) transitions; not measured by neutron)", group="L02_phase_diagram",
    note="the high-T phase assignment for x>0 is inferred; compare Ba: C060 x_Co=0.7-0.8 -> SRT 280-330 K and L12 x_Zn=1.3", T={"K": 300})
# impurities (Table 2)
IMP2 = {"SrCo2Fe16O27": [("spinel", "6.54", None), ("M-type", None, None), ("X-type", None, None)],
        "SrCoZnFe16O27": [("spinel", "7.54", "1.94 (I11)"), ("M-type", "0.19", "2.10 (I11)"), ("X-type", "0.52", "2.94 (I11)")],
        "SrZn2Fe16O27": [("spinel", "4.12", None), ("M-type", "1.77", None), ("X-type", "2.67", None)]}
for f, lst in IMP2.items():
    for ph, w, u in lst:
        if w is None:
            continue
        l02(f"impurity_wt_{ph}", f, value=float(w), unit="%", kind="fraction", status="measured", locator="Table 2, rows Wt% spinel / M-type / X-type hexaferrite",
            verbatim="Wt% spinel 6.54 | 7.54 (1.94) | 4.12 ; Wt% M-type hexaferrite - | 0.19 (2.10) | 1.77 ; Wt% X-type hexaferrite - | 0.52 (2.94) | 2.67 (parentheses = values from X-ray data (I11))",
            definition=f"{ph} weight fraction from combined refinement", unc_note=("parenthesis in Table 2 are I11-only phase fractions, not uncertainties" if u else None),
            note="for SrCoZnFe16O27 wt fractions from I11 and SHRPD were 'not reconcilable' (SI S1.1); spinel composition assumed from the W phase (Co0.5Zn0.5Fe2O4)" if f == "SrCoZnFe16O27" else None)

# magnetic structure descriptors
for f, msg, kind_ in (("SrCo2Fe16O27", "Cm'cm' (BNS 63.464)", "planar"), ("SrCoZnFe16O27", "Cm'cm'", "planar"), ("SrZn2Fe16O27", "P63/mm'c' (BNS 194.270, OG 194.8.1501)", "uniaxial")):
    l02("magnetic_structure_type", f, status="measured", qualitative=kind_, locator="Abstract; Section 3.2; Section 3.4",
        verbatim=("Using neutron powder diffraction on W-type hexaferrites, planar magnetic ordering described in Cm'cm' was found for SrCo2Fe16O27 and SrCoZnFe16O27, while SrZn2Fe16O27 had uniaxial ordering described in P63/mm'c'."),
        definition=f"magnetic space group {msg}; non-collinear symmetry-allowed canting modes active (mGM4+ for planar: canting toward c, AFM on 8f tet/8f oct/8f tbp; SrZn2: small AFM canting to ab on 6g and 12k)",
        T={"K": 300}, note="'truly collinear order not symmetry protected'; canting refined values in Tables S1-S3")

# per-site Fe moments (Tables S1-S3). row tuples: label, site, Mx, My, Mz
S1 = [("Fe1", "W1", None, "3.997(15)", "-1.90(4)"), ("Fe2", "W2", None, "3.901(14)", "-0.07(5)"), ("Fe3_1", "W3", "0.11(11)", "-4.11(7)", "-0.38(6)"), ("Fe3_2", "W3", None, "-3.96(13)", "-0.10(9)"),
      ("Fe4", "W4", None, "-3.791(13)", "-1.20(5)"), ("Fe5_1", "W5", None, "-3.04(4)", "-0.40(3)"), ("Fe5_2", "W5", None, "-3.56(7)", "-0.54(4)"), ("Fe6", "W6", None, "4.108(13)", "0.07(4)"), ("Fe7", "W7", None, "-3.236(17)", "-4.0(5)")]
S2 = [("Fe1", "W1", None, "4.32(5)", "-2.36(7)"), ("Fe2", "W2", None, "3.96(4)", "0.08(5)"), ("Fe3_1", "W3", "-0.8(2)", "-2.88(17)", "-0.19(5)"), ("Fe3_2", "W3", None, "-4.5(4)", "0.25(8)"),
      ("Fe4", "W4", None, "-3.547(16)", "-1.35(4)"), ("Fe5_1", "W5", "-0.02(4)", "-2.87(3)", "-0.41(3)"), ("Fe5_2", "W5", None, "-3.44(6)", "-0.32(5)"), ("Fe6", "W6", None, "4.081(16)", "0.13(5)"), ("Fe7", "W7", None, "-3.16(2)", "-4.0(4)")]
S3 = [("Fe1", "W1", None, None, "-3.99(4)"), ("Fe2", "W2", None, None, "-4.511(17)"), ("Fe3", "W3", "-0.159(9)", "-0.319(18)", "2.991(17)"), ("Fe4", "W4", None, None, "-3.554(13)"),
      ("Fe5", "W5", "0.162(7)", "0.324(13)", "3.091(5)"), ("Fe6", "W6", None, None, "-4.353(15)"), ("Fe7", "W7", None, None, "3.786(14)")]
for f, tab, tn in (("SrCo2Fe16O27", S1, "S1"), ("SrCoZnFe16O27", S2, "S2"), ("SrZn2Fe16O27", S3, "S3")):
    for lab, site, mx, my, mz in tab:
        for comp, v in (("Mx", mx), ("My", my), ("Mz", mz)):
            if v is None:
                continue
            note = None
            if f == "SrZn2Fe16O27" and site == "W4" and comp == "Mz":
                note = "SIGN CONFLICT: the same SI .mcif gives Fe4 the same sign as 6g/12k/2d (up-sublattice; magnitude 0.1089*32.829=3.575 muB) whereas Table S3 prints -3.554; sign here as printed"
            l02(f"Fe_moment_{comp}_{site}_{lab}", f, vs=v, unit="muB", kind="magnetic_moment", proposed_kind="magnetic_moment", status="measured",
                locator=f"Supporting information sup1.pdf, Table {tn} (Crystal and magnetic structure details), row {lab}",
                verbatim=f"Table {tn} row {lab} (Fe3+): {comp} = {v} muB  [row read from rendered page; text layer of SrCoZn/SrZn tables is scrambled]",
                definition=f"refined Fe3+ moment component {comp} along orthorhombic axes (a, b=sqrt3 a_hex, c) for Cmcm-based Cm'cm' (Co-containing) / hexagonal axes for P63/mm'c' (Zn); in the planar samples the in-plane component is My; Mz = canting toward c. Fe:Co moment ratio fixed 5:3",
                T={"K": 300, "note": "RT"}, note=note, pending_unit="muB")

# Derived weighted site occupancies for split sites (derived_by_us)
import collections
accum = {}
for r in coords:
    if r["source"] != "L02" or r["species"] in ("Fe", "Sr"):
        continue
    accum.setdefault((r["sample_id_in_paper"], r["site"], r["species"]), []).append(r)
for (f, site, sp), rows_ in accum.items():
    if site not in ("W3", "W5"):
        continue
    if len(rows_) < 2:
        continue
    mult = {"W3": (8, 4), "W5": (16, 8)}[site]
    tot = sum(mult)
    val = sum(m * float(vu(rw["occupancy"])[0]) for m, rw in zip(mult, rows_)) / tot
    l02(f"site_occupancy_{sp}_{site}_weighted_derived", f, value=round(val, 4), unit="1", kind="fraction", status="derived_by_us",
        locator="derived from CIF rows " + " + ".join(rw["cif_label"] for rw in rows_), verbatim=" | ".join(rw["verbatim_cif_line"] for rw in rows_),
        definition=f"{sp} fraction of unified site {site} = (8*occ_8d + 4*occ_4b)/12 for 6g; (16*occ_16h + 8*occ_8f)/24 for 12k (Cmcm splitting of P63/mmc sites)",
        transformation="multiplicity-weighted mean over the Cmcm sub-sites; uncertainty not propagated", note="no source uncertainty propagated")

# occupancy per CIF row (non-Fe)
for r in coords:
    if r["source"] != "L02" or r["species"] in ("Fe", "Sr"):
        continue
    f = r["sample_id_in_paper"]; sp = r["species"]
    l02(f"site_occupancy_{r['cif_label']}", f, vs=r["occupancy"], unit="1", kind="fraction", status="measured",
        locator=f"CIF {r['cif_file'].split('/')[-1]} (supplementary in5078sup2.zip) and SI Table {'S1' if f=='SrCo2Fe16O27' else 'S2' if f=='SrCoZnFe16O27' else 'S3'}, row {r['cif_label']} ({r['site_name']} = unified {r['site']})",
        verbatim=r["verbatim_cif_line"], definition=f"occupancy of {sp} at crystallographic (sub)site {r['cif_label']} (fraction of that sub-site); partially occupied CIF, NOT a DFT configuration; "
        "Co2+/Zn2+ refined on all sites with Fe3+ as remainder", T={"K": 300})
