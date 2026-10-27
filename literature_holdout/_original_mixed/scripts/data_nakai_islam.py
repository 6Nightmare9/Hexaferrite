"""Records from L03 (Nakai 2023), C004 (Nakai 2025), L05 (Islam 2025)."""
from rc import add

# ------------------------------------------------------------------ L03
SYN3 = ("solid-state: SrCO3 (3N, kept at 473 K before weighing), alpha-Fe2O3 (4N), Co3O4 (3N), ZnO (4N), NiO (4N), MgO (3N); stoichiometric SrMe2Fe16O27, agate mortar, pellets, pre-baked several times; "
        "pO2 = 0.2 atm (air) and 1 atm (O2 stream): 24 h, temperature NOT stated in text; 10 atm: pressure vessel furnace 1423-1473 K, 12 h; 387 atm: HIP O2-Ar (Ar:O2 = 4:1) ~1473 K, 2 h; cooling rate not stated")


def l03(name, me, **kw):
    kw.setdefault("sample", f"{me}2W"); kw.setdefault("formula", f"Sr{me}2Fe16O27"); kw.setdefault("cluster", "NAKAI"); kw.setdefault("tier", 1 if me in ("Co", "Zn") else 2)
    kw.setdefault("synth", SYN3); kw.setdefault("form", "sintered pellet (polycrystal), polished and C-coated for WDX"); kw.setdefault("texture", "not stated")
    kw.setdefault("T", {"K": 300, "note": "analysis at RT"})
    return add("L03", name, **kw)


for p, h, tt in ((0.2, 24, "air; T not stated"), (1.0, 24, "O2 stream; T not stated"), (10.0, 12, "1423-1473 K"), (387.0, 2, "~1473 K (HIP)")):
    l03("synthesis_oxygen_partial_pressure", "Me", value=p, unit="atm", kind="pressure", proposed_kind="pressure", pending_unit="atm", status="measured",
        locator="Section 2 Experimental Procedures, para 1; Abstract", verbatim="synthesized by solid state reactions under four different oxygen partial pressures, pO2 = 0.2, 1, 10, and 387 atm",
        definition=f"sintering oxygen partial pressure ({tt}; duration {h} h)", group="L03_synth", sample="Me2W (Me=Co,Zn,Ni,Mg)", formula="SrMe2Fe16O27",
        note="Mn not studied in this paper")
l03("synthesis_temperature_10atm_min", "Me", value=1423.0, unit="K", kind="temperature", status="measured", locator="Section 2, para 1",
    verbatim="The synthesis under pO2 = 10 atm was done in a high-pressure vessel furnace (Crystal Systems, HPF-1500-M) at 1423-1473 K for 12 hours.", sample="Me2W (Me=Co,Zn,Ni,Mg)", formula="SrMe2Fe16O27", group="L03_synth")
l03("synthesis_temperature_10atm_max", "Me", value=1473.0, unit="K", kind="temperature", status="measured", locator="Section 2, para 1",
    verbatim="at 1423-1473 K for 12 hours", sample="Me2W (Me=Co,Zn,Ni,Mg)", formula="SrMe2Fe16O27", group="L03_synth")
l03("synthesis_temperature_387atm", "Me", value=1473.0, unit="K", kind="temperature", status="measured", locator="Section 2, para 1",
    verbatim="The synthesis at pO2 = 387 atm were made in a hot isostatic pressing (HIP) equipment (Kobe Steel, O2-Dr.HIP) for 2 hours at around 1473 K under 2000 kgf/cm3 of oxygen-argon mixed gas (Ar:O2 = 4:1).",
    definition="'around 1473 K'; temperature gradient in HIP may be steep (authors)", sample="Me2W (Me=Co,Zn,Ni,Mg)", formula="SrMe2Fe16O27", group="L03_synth")
l03("Me_content_nstar_at_387atm", "Co", value=2.03, unit="1", kind="dimensionless", status="measured", locator="Section 3.2 SEM-WDX composition analysis (polycrystal), Fig. 3",
    verbatim="in Co2W, n*(Co) = 2.03 at pO2 = 387 atm, which is slightly higher than the ideal value", definition="n*(X) = 18 x n(X)/(n(Fe)+n(Me)): Co per formula unit (Sr:(Fe+Me)=1:18 always found); WDX, 3-4 spots averaged",
    method="SEM-WDX (Oxford INCA Wave 500), standards: Sr from SrF2, Fe/Co/Zn/Ni metals; SrFe12O19 crystal for intensities", sample="Co2W (polycrystal, pO2=387 atm)", formula="SrCo2Fe16O27", group="L03_WDX_Co2W_387atm",
    note="at lower pO2 the Fe content is higher and Co lower than 16:2 (Fig. 3: figure_only); delta in SrMe2-dFe16+dO27 attributed to Fe2+")
l03("Me_content_nstar_at_387atm", "Ni", value=1.93, unit="1", kind="dimensionless", status="measured", locator="Section 3.2, Fig. 3",
    verbatim="in Ni2W, n*(Ni) = 1.93 at pO2 = 387 atm, which is close to the ideal value", definition="n*(Ni) = 18 x n(Ni)/(n(Fe)+n(Ni))", method="SEM-WDX", sample="Ni2W (polycrystal, pO2=387 atm)", formula="SrNi2Fe16O27", group="L03_WDX_Ni2W_387atm")
l03("Me_and_Fe_content_vs_pO2", "Co", status="measured", figure_only=True, locator="Fig. 3",
    axes="x: pO2 (atm), log scale 0.1-1000; left y: n*(Fe) 15.6-16.4; right y: n*(Me) 1.6-2.4; four series Co, Ni, Zn, Mg (Fe in Me2W and Me in Me2W)",
    verbatim="Fig. 3 Fe and Me contents in Me2W versus partial oxygen pressure at the sintering process. Definition of n*(X) is given in the text.",
    definition="actual Fe and Me content per f.u. in the W phase vs sintering pO2; deviation from 16:2 order in air: Co < Ni < Zn < Mg; Zn2W and Mg2W do not reach ideal even at 387 atm", sample="Me2W", formula="SrMe2Fe16O27", group="L03_WDX_all",
    note="THE key evidence that nominal != actual composition and Fe2+ appears; needs digitisation (Stage 2 part 2) if used as a quantitative prior for Fe2+/vacancy content")
CRY = {"M-type (SrFe12O19)": (1.00, 11.93, 0.07), "X-type (Sr2MeFe28O46)": (1.96, 28.12, 1.88), "W-type (SrMeFe16O27)": (1.00, 16.21, 1.79), "spinel (MeFe2O4)": (0.00, 2.13, 0.87)}
for ph, (sr, fe, co) in CRY.items():
    for el, v in (("Sr", sr), ("Fe", fe), ("Co", co)):
        l03(f"Co2W_crystal_phase_content_{ph.split()[0]}_{el}", "Co", value=v, unit="1", kind="dimensionless", status="measured", locator="Table 1 (Sr, Fe, Co content of Co2W and impurities in the Co2W crystal)",
            verbatim=f"{ph}: Sr {sr:.2f}  Fe {fe:.2f}  Co {co:.2f}  (normalised so that Fe + Co = 12 (M), 30 (X), 18 (W), 3 (spinel))", definition=f"{el} content normalised to Fe+Co fixed by phase; WDX on flux-grown crystal",
            sample="Co2W flux-grown crystal", formula="SrCo2Fe16O27", form="single-crystal-like hexagonal prisms up to 5x5x5 mm, acid-treated; multiphase (not a perfect single crystal)",
            synth="flux growth: SrCO3, Fe2O3, Co3O4, B2O3 = Sr:Fe:Co:B = 3:14:2:2 in Pt crucible, 1613 K 20 h, cooled to 1223 K at -2.5 K/h, then -100 K/h; HCl flux removal", group="L03_Co2W_crystal",
            note="Co content of the W phase in the crystal (1.79) is lower than polycrystal (2.03 at 387 atm)")
l03("no_stacking_faults_Co2W_TEM", "Co", status="measured", qualitative="no stacking faults", locator="Section 3.4; Abstract",
    verbatim="No stacking faults were observed in a total of about 20 dark-field images taken from [1-100] or [11-20].", definition="STEM/TEM on ion-milled Co2W crystal", group="L03_Co2W_TEM", sample="Co2W crystal", formula="SrCo2Fe16O27")
l03("phase_purity_vs_pO2_qualitative", "Me", status="measured", qualitative="spinel impurity decreases with pO2 up to 10 atm; at 387 atm: Ni2W single phase, Co2W M-type+hematite (no spinel), Zn2W more spinel, Mg2W four phases", locator="Section 3.1, Fig. 2",
    verbatim="In Ni2W, the spinel ferrite disappears, and a single W phase is obtained. In Zn2W, the peaks of impurity spinel ferrite appear larger than those at pO2 = 10 atm.",
    definition="XRD (Cu Ka1, 2theta 20-50 deg); no quantitative phase fractions given", sample="Me2W", formula="SrMe2Fe16O27", group="L03_XRD_all", note="no wt% in the paper: a gap for phase-fraction priors (need digitisation of Fig. 2 or new data)")
l03("recommended_pO2_for_ideal_composition", "Me", value=100.0, unit="atm", kind="pressure", proposed_kind="pressure", pending_unit="atm", status="assumed_by_authors", locator="Conclusion",
    verbatim="To synthesize Me-stabilized W with the ideal composition a sufficiently high oxygen pressure (>=100 atm) is necessary to suppress the appearance of Fe2+.", definition="lower bound claimed by authors (inequality >= 100 atm); sample Me2W", sample="Me2W", formula="SrMe2Fe16O27", group="L03_conclusion")

# ------------------------------------------------------------------ C004
SYN4 = ("solid-state in sealed evacuated quartz tube (pO2 ~1e-8 atm at RT): SrFe12O19 + (8/3) alpha-Fe2O3 + (2/3) Fe (Fe reduced in H2) -> SrFe18O27, pellets, 1373-1523 K 12 h, cooled to ~500 K at 200 K/h then furnace cooled")


def c004(name, **kw):
    kw.setdefault("cluster", "NAKAI"); kw.setdefault("tier", 2); kw.setdefault("synth", SYN4); kw.setdefault("form", "polycrystalline pellet / powder (~100 um multidomain particles for VSM)")
    kw.setdefault("method", "XRD (X'Pert PRO, Cu Ka1; HighScore Plus phase ratios) + RT VSM (VSM-5SC-10HF, TOEI), high-field part up to 6.4 MA/m extrapolated to 1/H^2=0")
    kw.setdefault("T", {"K": 300, "note": "RT"}); kw.setdefault("texture", "unfixed particles freely rotating in field (random)")
    return add("C004", name, **kw)


for s, T, w, ms in (("S1448", 1448, "97.1", "81.1(1)"), ("S1473", 1473, "98.0", "81.1(1)"), ("S1523", 1523, "97.9", "81.0(1)")):
    c004("SrW_weight_fraction", vs=w, unit="%", kind="fraction", status="measured", locator="Table 2 (ratio of constituent phases and saturation magnetization at RT), p. 237",
         verbatim=f"{s}  SrW (wt %) {w}  SrM --  Fe3O4 (rest)  M (A m2/kg) {ms}", sample=s, formula="SrFe18O27", definition="SrW wt % by XRD intensity analysis (secondary phase: Fe3O4 2-3 wt %)", group=f"C004_{s}",
         note=f"synthesis T {T} K, 12 h", unc_note="no uncertainty given for wt %")
    c004("Ms_RT", vs=ms, unit="A*m^2/kg", kind="specific_magnetization", proposed_kind="specific_magnetization", pending_unit="A*m^2/kg", status="measured", locator="Table 2, p. 237",
         verbatim=f"{s} ... M (A m2/kg) {ms}", sample=s, formula="SrFe18O27", definition="saturation magnetization at RT from 1/H^2 extrapolation of high-field magnetisation (field roughly along easy axis for unfixed particles); includes 2-3 wt % Fe3O4",
         group=f"C004_{s}", transformation="none (1 A m2/kg = 1 emu/g)")
for s, Tm, ms in (("S1373", 1373, "78.9(1)"), ("S1423", 1423, "78.8(1)")):
    c004("Ms_RT_no_W_phase", vs=ms, unit="A*m^2/kg", kind="specific_magnetization", proposed_kind="specific_magnetization", pending_unit="A*m^2/kg", status="measured", locator="Table 2, p. 237",
         verbatim=f"{s}  SrW --  SrM 63-65 wt %  Fe3O4 33-35 wt %  alpha-Fe2O3 2 wt %  M {ms}", sample=s, formula="SrFe18O27 (nominal; actually SrM + Fe3O4 + a-Fe2O3)", definition="mixture below the W-phase stability boundary: not a W-phase value", group=f"C004_{s}")
c004("W_phase_lower_stability_limit_TL_min", value=1423.0, unit="K", kind="temperature", status="measured", locator="Section 3.2 and Fig. 5 caption; Abstract ('only above 1448 K')",
     verbatim="TL is estimated to be between 1423 and 1448 K", sample="SrFe18O27", formula="SrFe18O27", definition="redox-free decomposition SrW -> SrM + 2 Fe3O4 boundary (pO2-independent); SrW formed only above it", group="C004_phase_boundary", note="lower bound of TL")
c004("W_phase_lower_stability_limit_TL_max", value=1448.0, unit="K", kind="temperature", status="measured", locator="Section 3.2 and Fig. 5 caption",
     verbatim="TL is estimated to be between 1423 and 1448 K, and thus PL is estimated to be between PM(1423) and PM(1448), i.e. 2-5 x 10^-4 atm", sample="SrFe18O27", formula="SrFe18O27", definition="upper bound of TL", group="C004_phase_boundary")
c004("SrW_stable_pO2_window_at_1473K", status="assumed_by_authors", qualitative="SrW stable between pO2 = 10^-4.7 and 10^-2.7 atm at 1473 K (literature phase diagram, ref 25)", locator="Section 1, Fig. 2",
     verbatim="SrW is stable between PO2 = 10-4.7 and 10-2.7 atm at 1473 K", sample="SrFe18O27", formula="SrFe18O27", definition="from reported PO2-T diagram (secondary)", group="C004_intro", split="excluded_secondary_citation", cluster="SECONDARY")
c004("literature_Ms_RT_SrW", value=77.0, unc=None, unit="A*m^2/kg", kind="specific_magnetization", proposed_kind="specific_magnetization", pending_unit="A*m^2/kg", status="measured", locator="Section 1",
     verbatim="The value of magnetization of SrW at room temperature (RT) (77-81 A m2/kg 5-7)) is 10 % higher than that of M-type strontium ferrite", sample="SrFe18O27 (review of refs 5-7)", formula="SrFe18O27", definition="lower end of a quoted literature range 77-81 (primary sources not in corpus)",
     group="C004_intro", split="excluded_secondary_citation", cluster="SECONDARY", note="secondary citation: do not use as an independent datum")
c004("literature_magnetic_anisotropy_RT_SrW_low", value=1.30, unit="MA/m", kind="H_field", proposed_kind="anisotropy_field_or_unspecified", pending_unit="MA/m", status="measured", locator="Section 1",
     verbatim="The magnetic anisotropy of SrW at RT (1.30-1.55 MA/m 5-7)) is large and comparable to that of SrM (1.48 MA/m)", sample="SrFe18O27 (review of refs 5-7)", formula="SrFe18O27",
     definition="quantity not defined in the text (units MA/m => presumably anisotropy field H_A); lower end of range 1.30-1.55", group="C004_intro", split="excluded_secondary_citation", cluster="SECONDARY",
     note="definition ambiguous (HA vs mu0HA vs K); compare L05 DFT mu0Ha = 0.569 T and L11 Ba: 13.7 kOe = 1.09 MA/m")

# ------------------------------------------------------------------ L05 (DFT)
CMP = {"SrW": "SrFe18O27", "SrKW": "Sr0.5K0.5Fe18O27", "SrCaW": "Sr0.5Ca0.5Fe18O27", "SrBaW": "Sr0.5Ba0.5Fe18O27", "SrLaW": "Sr0.5La0.5Fe18O27"}
DFT = ("VASP, PAW, PBE GGA+U (Dudarev, U_eff = 3 eV on Fe 3d, from linear response), 92-atom hexagonal cell P63/mmc (SrW), gamma-centred 5x5x1 k-mesh (SCF), cutoff 600 eV, "
       "SCF 1e-6 eV, forces <0.01 eV/A; SOC noncollinear for MAE with 12x12x2 k-mesh and Bloechl tetrahedra; collinear Gorter ferrimagnet 4e,4f1,4f3 antiparallel to 6g,4f2,12k,2d; "
       "K/Ca/Ba/La replace half of the Sr atoms at 2b (symmetry P-6m2)")


def l05(name, c, **kw):
    kw.setdefault("sample", c); kw.setdefault("formula", CMP[c]); kw.setdefault("cluster", "ISLAM_L05"); kw.setdefault("tier", 2)
    kw.setdefault("synth", "n/a (DFT)"); kw.setdefault("method", DFT); kw.setdefault("status", "calculated"); kw.setdefault("T", {"K": 0, "note": "static DFT, no thermal effects"})
    kw.setdefault("form", "ideal crystal, 92-atom cell"); kw.setdefault("texture", "n/a"); kw.setdefault("density", "n/a"); kw.setdefault("orient", "magnetisation along c (easy axis) and a for MAE")
    kw.setdefault("group", f"L05_{c}_DFT")
    return add("L05", name, **kw)


l05("U_eff_Fe3d", "SrW", value=3.0, unit="eV", kind="energy", locator="Section 2 Computational details", verbatim="yielding a value of approximately U ~ 3 eV ... a Hubbard Ueff of 3 eV was applied specifically for Fe-3d orbitals",
    definition="U_eff = U - J (Dudarev), linear-response (Cococcioni) estimate, applied to W-type; other phases: 4.5 eV (M-type, ref 55), 3.8 eV (spinel, ref 56)", status="assumed_by_authors", note="'approximately 3 eV': no uncertainty; sensitivity to U not studied in this paper (needs our own check)")
l05("plane_wave_cutoff", "SrW", value=600.0, unit="eV", kind="energy", locator="Section 2", verbatim="a plane-wave basis set with a cutoff energy of 600 eV was employed", definition="PAW plane-wave cutoff", status="assumed_by_authors")
l05("SCF_energy_threshold", "SrW", value=1e-6, unit="eV", kind="energy", locator="Section 2", verbatim="total energy threshold of 10-6 eV ... A more stringent threshold of 10-8 eV was also tested", definition="SCF convergence; 1e-8 eV tested, energy differences ~1e-8 eV", status="assumed_by_authors")
l05("DFT_settings_summary", "SrW", qualitative="PBE+U PAW VASP; SOC noncollinear; k 5x5x1 (SCF), 12x12x2 (MAE); forces<0.01 eV/A; relaxed within GGA (no U) for lattice, GGA+U for energies", locator="Section 2",
    verbatim="the U parameter was estimated ... k-point mesh of 5 x 5 x 1 ... denser k-point mesh of 12 x 12 x 2 was utilized to ensure well-converged values for the MAE", definition="settings text", note="structures relaxed within GGA scheme (lattice parameters Table 1 are GGA); MAE with GGA+U+SOC")
T1 = {"SrW": ("5.877", "32.842", "982.59"), "SrKW": ("5.875", "32.875", "982.59"), "SrCaW": ("5.875", "32.762", "979.34"), "SrBaW": ("5.885", "32.918", "987.26"), "SrLaW": ("5.880", "32.753", "980.63")}
for c, (a, cc, v) in T1.items():
    l05("lattice_a_DFT_relaxed", c, vs=a, unit="angstrom", kind="length", locator="Table 1 (calculated optimised lattice parameters), p. 3",
        verbatim=f"{c}  a {a}  b {a}  c {cc}  V {v}", definition="a = b of relaxed hexagonal cell (GGA)", note="calculated a is slightly smaller than experiment (5.89-5.91 A); structure P-6m2 for substituted")
    l05("lattice_c_DFT_relaxed", c, vs=cc, unit="angstrom", kind="length", locator="Table 1, p. 3", verbatim=f"{c}  a {a}  b {a}  c {cc}  V {v}", definition="c of relaxed cell (GGA)")
    l05("cell_volume_DFT_relaxed", c, vs=v, unit="angstrom^3", kind="volume", proposed_kind="volume", pending_unit="angstrom^3", locator="Table 1, p. 3", verbatim=f"{c}  ...  V {v}", definition="volume per unit cell (A^3)")
l05("formation_enthalpy_GGA", "SrW", vs="-2.815", unit="eV/atom", kind="energy_per_atom", proposed_kind="energy_per_atom", pending_unit="eV/atom", locator="Table 1", verbatim="SrW ... Hf -2.815 (-2.208)",
    definition="formation enthalpy from elements, GGA (mixed GGA/GGA+U with O2 correction 1.36 eV: -2.208)", note="not benchmarked against experiment (none available)")
l05("decomposition_energy_vs_SrM_plus_2Fe3O4_GGAU", "SrW", vs="-8.511", unit="eV/f.u.", kind="energy_per_fu", proposed_kind="energy_per_formula_unit", pending_unit="eV/f.u.", locator="Table 1",
    verbatim="SrW ... Edec -0.035 (-8.511)  [GGA (GGA+U)]", definition="E(SrFe18O27) - E(SrFe12O19) - 2E(Fe3O4); GGA: -0.035 eV/f.u.; GGA+U (U=3/4.5/3.8 eV for W/M/spinel): -8.511 eV/f.u.",
    note="GGA and GGA+U values differ by 8.5 eV because different U per phase: cross-phase GGA+U energies NOT comparable; GGA value -0.035 (tiny) agrees with C004 'driving force small'")
l05("decomposition_energy_vs_SrM_plus_2Fe3O4_GGA", "SrW", vs="-0.035", unit="eV/f.u.", kind="energy_per_fu", proposed_kind="energy_per_formula_unit", pending_unit="eV/f.u.", locator="Table 1",
    verbatim="SrW ... Edec -0.035", definition="GGA decomposition energy")
# Table 2
SITES = ["4e", "4f1", "6g", "4f2", "12k", "4f3", "2d"]
WID = {"4e": "W1", "4f1": "W2", "6g": "W3", "4f2": "W4", "12k": "W5", "4f3": "W6", "2d": "W7"}
T2 = {"SrW": [("-3.950", "-0.017"), ("-3.952", "-0.017"), ("3.859", "0.025"), ("4.013", "0.025"), ("4.104", "0.016"), ("-3.984", "-0.014"), ("3.986", "0.016")],
      "SrKW": [("-3.948", "-0.016"), ("-3.946", "-0.017"), ("3.878", "0.026"), ("4.067", "0.021"), ("4.092", "0.016"), ("-3.968", "-0.014"), ("3.943", "0.016")],
      "SrCaW": [("-3.944", "-0.017"), ("-3.944", "-0.017"), ("3.811", "0.027"), ("4.068", "0.024"), ("4.103", "0.016"), ("-3.979", "-0.015"), ("3.986", "0.018")],
      "SrBaW": [("-3.947", "-0.017"), ("-3.946", "-0.017"), ("3.818", "0.028"), ("4.067", "0.022"), ("4.102", "0.016"), ("-3.989", "-0.015"), ("3.986", "0.019")],
      "SrLaW": [("-3.947", "-0.017"), ("-3.951", "-0.017"), ("3.793", "0.034"), ("4.041", "0.034"), ("4.089", "0.017"), ("-3.999", "-0.015"), ("3.982", "0.019")]}
for c, rows in T2.items():
    for s, (ms, ml) in zip(SITES, rows):
        l05(f"Fe_spin_moment_{WID[s]}_{s}", c, vs=ms, unit="muB", kind="magnetic_moment", proposed_kind="magnetic_moment", pending_unit="muB", locator="Table 2 (calculated average total spin ms and orbital ml moments), p. 5",
            verbatim=f"Fe({s})  ms = {ms}  ml = {ml}  [{c}]", definition="spin moment in PAW sphere (U=3 eV); sign: negative = antiparallel to net moment; orbital moment is along easy axis (c)",
            note="Table text layer misaligned; values read from rendered page")
        l05(f"Fe_orbital_moment_{WID[s]}_{s}", c, vs=ml, unit="muB", kind="magnetic_moment", proposed_kind="magnetic_moment", pending_unit="muB", locator="Table 2, p. 5",
            verbatim=f"Fe({s})  ms = {ms}  ml = {ml}  [{c}]", definition="orbital moment along the easy axis, GGA+U+SOC")
MT = {"SrW": ("54.991", "0.303", "0.652"), "SrKW": ("55.389", "0.307", "0.657"), "SrCaW": ("54.990", "0.324", "0.655"), "SrBaW": ("54.964", "0.322", "0.649"), "SrLaW": ("54.050", "0.410", "0.643")}
for c, (m, ml, ms) in MT.items():
    l05("total_spin_moment_per_cell", c, vs=m, unit="muB", kind="magnetic_moment", proposed_kind="magnetic_moment", pending_unit="muB", locator="Table 2, row m_tot/u.c.", verbatim=f"m_tot/u.c. {m} (ml {ml}); mu0Ms {ms} T",
         definition="total magnetic moment per unit cell (2 f.u.), 'ms' column (spin) ", note="does not equal the sum of the PAW-sphere site moments (48.9 muB from the table x multiplicities) - includes interstitial/O contributions")
    l05("mu0Ms_DFT", c, vs=ms, unit="T", kind="magnetic_polarization_mu0M", proposed_kind="magnetic_polarization_mu0M", pending_unit="T", locator="Table 2, row mu0Ms (T)", verbatim=f"mu0Ms (T) {ms}", definition="mu0*Ms (spontaneous magnetisation expressed in tesla, polarisation) from m_tot/V; NOT B inside the material",
         note="Ms = 0.652 T / mu0 = 5.19e5 A/m (derived_by_us); experimental mtot 56.60 muB/u.c. (ref 65, Ueda; secondary)")
T3 = {"SrW": ("0.148", "0.569", "0.661"), "SrKW": ("0.244", "0.932", "0.842"), "SrCaW": ("0.326", "1.254", "0.978"), "SrBaW": ("0.301", "1.167", "0.948"), "SrLaW": ("0.582", "2.276", "1.331")}
for c, (k, ha, ka) in T3.items():
    l05("Ku_uniaxial_MAE_constant", c, vs=k, unit="MJ/m^3", kind="energy_density", proposed_kind="energy_density", pending_unit="MJ/m^3", locator="Table 3 (calculated MAE constant Ku ...), p. 6; Fig. 4",
         verbatim=f"Ku (This work) {k} MJ/m3  [{c}]", definition="Ku = MAE/V, MAE = E(M||a) - E(M||c) from GGA+U+SOC total energies; Ku>0 uniaxial; equals effective first-order constant (K1+2K2 not separated); 3 methods in Fig. 4 agree in trend (DFT energy, SOC energy, orbital anisotropy/Bruno)",
         note="dipole-dipole (shape) contribution ~0.033 % of MAE and not included in Ku; at T = 0 K; experiment 0.344 MJ/m3 at 300 K quoted from Ueda (secondary): DFT underestimates by ~2x")
    l05("mu0Ha_anisotropy_field_DFT", c, vs=ha, unit="T", kind="mu0H", locator="Table 3, row mu0Ha (T)", verbatim=f"mu0Ha (T) {ha}  [{c}]",
         definition="mu0*Ha = 2Ku/Ms (uniaxial anisotropy field as mu0H; includes no demagnetisation)", status="calculated", transformation="authors' formula; our check: 2*0.148e6/5.19e5 = 0.570 T (matches 0.569)", note="Ha defined by an effective single K (K1) formula; not the same as HA from approach-to-saturation fits")
    l05("hardness_parameter_kappa", c, vs=ka, unit="1", kind="dimensionless", locator="Table 3, row kappa", verbatim=f"kappa {ka}  [{c}]", definition="magnetic hardness parameter kappa = sqrt(Ku/(mu0 Ms^2)) (formula text garbled in PDF; our check reproduces 0.661 for SrW)", status="calculated")
l05("MAE_site_resolved", "SrW", status="calculated", figure_only=True, locator="Fig. 5; Fig. 6 (d-orbital resolved); Fig. 7 (PDOS)", axes="Fig. 5: bars of sublattice MAE contribution (meV) for 4e,4f1,6g,4f2,12k,4f3,2d; error bars = std over atoms in sublattice",
    verbatim="in the pure SrW, the MAE contribution is primarily governed by the Fe(6g), Fe(4f2), Fe(12k), and Fe(2d) site with small contributions from Fe(4e), Fe(4f1), and Fe(4f3).",
    definition="site-resolved MAE (second-order perturbation, SOC-energy method); 2d (bipyramidal) gives the largest single-ion contribution per atom", group="L05_SrW_DFT")
l05("exchange_J_not_computed", "SrW", status="calculated", qualitative="no exchange constants J_ij in this paper", locator="Conclusion, last sentence",
    verbatim="Nonetheless, to fully understand the magnetic interactions in these systems, further calculations of the exchange coupling parameters are required", definition="J_ij explicitly deferred by authors; L04 (Islam 2024 PRB) full text not available -> J_ij for SrW not in corpus", note="GAP: exchange J for W-phase: need DFT (our own)")
# secondary exp values quoted in L05
l05("Ku_experimental_300K_cited", "SrW", value=0.344, unit="MJ/m^3", kind="energy_density", proposed_kind="energy_density", pending_unit="MJ/m^3", status="measured", locator="Section 3.4 and Table 3 footnote a (Ref. [65])",
    verbatim="which is roughly half of the experimental value of 0.344 MJ/m3 obtained at room temperature [65]", definition="experimental Ku of SrFe18O27 at 300 K quoted from Ueda et al. [65]; primary source not in corpus (anisotropy field decreases below ~100 K there)",
    split="excluded_secondary_citation", cluster="SECONDARY", note="secondary citation; locate primary (Ueda) before using as holdout")
l05("m_tot_experimental_cited", "SrW", value=56.60, unit="muB", kind="magnetic_moment", proposed_kind="magnetic_moment", pending_unit="muB", status="measured", locator="Table 2, last row (Expt.), footnote a (Ref. [65])",
    verbatim="mtot/u.c. (Expt.) 56.60 [Ref. 65]", definition="experimental total moment per unit cell from Ueda [65]; temperature not stated here", split="excluded_secondary_citation", cluster="SECONDARY", note="secondary citation")
