"""Records from tier-3 (Ba / other W) and tier-4 (methodology) sources: L11, L12, C023, C015, C036, C044, C052, C053, C056, C057, C060, C066, C112, C115."""
from rc import add

BA = "Ba -> Sr transfer NOT justified a priori: Ba2+ (1.35 A) vs Sr2+ (1.13-1.18 A) changes c/a, Fe2+ stability (Sr-W needs low pO2, ref C004) and site preference; use only as holdout/transfer test with inflated uncertainty"

# ------------------------------------------------------------------ L11 He 2024 (BaMe2Fe16O27, Me = Fe, Mn, Zn, Ni, Co)
SYN11 = ("solid-state: BaCO3 (99.0 %), Fe2O3 (99.0 %), ZnO (99.999 %), MnO (99.0 %), NiO (99.0 %), Co3O4 (99.0 %); planetary ball mill 200 rpm 2 h in deionised water; "
         "sintered 1300 C 3 h in muffle furnace (air assumed, not stated), except BaFe18O27 sintered in tube furnace under N2")
T11 = {"Fe": ("95.34", "375", "13739.57"), "Co": ("74.90", "71", "10750.90"), "Ni": ("69.49", "457", "10700.60"), "Zn": ("72.90", "439", "12206.23"), "Mn": ("47.77", "429", "11301.48")}


def l11(name, me, **kw):
    kw.setdefault("sample", f"BaMe2Fe16O27 Me={me}"); kw.setdefault("formula", f"BaFe16{me}2O27" if me != "Fe" else "BaFe18O27"); kw.setdefault("cluster", "HE_L11")
    kw.setdefault("tier", 3); kw.setdefault("synth", SYN11); kw.setdefault("form", "powder (hexagonal plates, mean grain size 5.7-12.8 um from SEM)"); kw.setdefault("texture", "BaFe18O27: (00l) intensity variation -> some texture; others random")
    kw.setdefault("T", {"K": 300}); kw.setdefault("transfer", BA); kw.setdefault("group", f"L11_{me}")
    return add("L11", name, **kw)


for me, (ms, hc, ha) in T11.items():
    l11("Ms_RT", me, vs=ms, unit="A*m^2/kg", kind="specific_magnetization", proposed_kind="specific_magnetization", pending_unit="A*m^2/kg", status="measured",
        locator="Table 2 (room temperature magnetic parameters), p. 066801-5; Fig. 5(b)-(f)", verbatim=f"{me}  Ms = {ms} emu/g  Hcj = {hc} Oe  Ha = {ha} Oe",
        definition="Ms from linear fit M vs 1/H^2 (law of approach to saturation, M = Ms(1 - B/H^2)), PPMS-VSM 300 K, +-30 kOe", field="+-30 kOe (fit at high field, 1/H^2 axis 2-7 x 1e-9 Oe^-2)",
        transformation="emu/g -> A m^2/kg: factor 1 (numerically identical)", method="VSM (Quantum Design DynaCool), law of approach to saturation")
    l11("Hcj_RT", me, vs=hc, unit="Oe", kind="H_field", status="measured", locator="Table 2, p. 066801-5", verbatim=f"{me}  Hcj = {hc} Oe", definition="intrinsic coercivity (H_cj) from M(H) loop; multidomain grains", method="VSM", field="+-30 kOe",
        transformation="none (CGS Oe kept; 1 Oe = 79.5775 A/m per paper)")
    l11("Ha_anisotropy_field_RT", me, vs=ha, unit="Oe", kind="H_field", status="measured", locator="Table 2, p. 066801-5; Eq. (3)", verbatim=f"{me}  Ha = {ha} Oe",
        definition="Ha obtained from the 1/H^2 coefficient B = Ha^2/15 = 4K1^2/(15 Ms^2) of the law of approach to saturation (uniaxial hexagonal, K2=0 assumed); NOT a direct HA measurement; Co2W is easy-plane so the formula (uniaxial) is not applicable there (authors note Co2W planar)",
        method="VSM, law of approach to saturation", transformation="none (Oe); 13739.57 Oe = 1.0934e6 A/m (derived_by_us)",
        note="Co sample: planar anisotropy -> value is formal; compare DFT mu0Ha(SrW)=0.569 T vs Ba Fe2W 1.374 T (13739.57 Oe as mu0Ha)" if me == "Co" else "compare DFT mu0Ha(SrW)=0.569 T (L05) vs 1.374 T here (13739.57 Oe as mu0H): a factor ~2.4 difference -> method/Ba-Sr/temperature differences, see CONFLICTS" if me == "Fe" else None)
l11("assumed_dopant_moments", "Mn", status="assumed_by_authors", qualitative="Ni2+ 3.0, Co2+ 4.0, Mn2+ 2.7 muB (vs Fe3+ 5 muB)", locator="Section 3.4",
    verbatim="the magnetic moments of the dopant ions Ni2+, Co2+, and Mn2+ (3.0 muB, 4.0 muB, and 2.7 muB, respectively) are smaller than Fe3+ ions (5 muB)", definition="textbook-type values used in the argument; Mn2+ 2.7 muB conflicts with high-spin Mn2+ (S=5/2, ~5 muB); not measured here",
    note="do NOT use as Mn parameter; flagged in CONFLICTS", group="L11_Mn")
l11("dielectric_permittivity_and_loss_vs_f", "Fe", status="measured", figure_only=True, locator="Fig. 6 (a) epsilon' and (b) dielectric loss",
    axes="x: frequency 100 Hz - 10 MHz (log); y(a): dielectric constant ~10^2-10^3 (log), y(b): dielectric loss 10^0-10^2 (log); five samples Fe, Mn, Ni, Zn, Co; Agilent 4980 LCR meter, RT",
    verbatim="the dielectric constant and dielectric loss decrease accompanying the frequency range from 100 Hz to 1 MHz. Above 1 MHz, the dielectric permittivity almost remains unchanged.",
    definition="eps'(f), tan-delta-like loss of pressed powder/pellet (geometry, density, electrodes not stated)", group="L11_dielectric", note="sample geometry/density not given -> intrinsic eps not extractable; Fe2W highest eps (mixed valence hopping)",
    form="pellet (not described)")

# ------------------------------------------------------------------ L12 Zhuravlev 2005
Z = dict(cluster="ZHURAVLEV", tier=3, transfer=BA, synth="polycrystalline BaCo2-xZnxFe16O27 (preparation not described in this paper)", form="polycrystal; toroids (<=1 GHz coax-waveguide), pins (cavity 1-27 GHz), thin plate (0.2-1.5 GHz microstrip resonator)",
         T={"K": None, "note": "temperature not stated (room temperature implied)"})
add("L12", "easy_plane_to_easy_axis_boundary_x_Zn", value=1.3, unit="1", kind="dimensionless", status="assumed_by_authors", locator="Section 3 Experimental, para 1",
    verbatim="At x < 1.3, an easy-plane state is realized; at x > 1.3, an easy-axis state is observed.", sample="BaCo2-xZnxW system", formula="BaCo2-xZnxFe16O27",
    definition="Zn content x of the easy-plane/easy-axis boundary (T not stated; anisotropy values taken from the authors' refs 7, 8)", group="L12_system", **Z,
    note="Ba system; compare L02 Sr: planar for x_Co=2, 1; boundary extrapolated at SrCo0.65Zn1.35 (x_Zn=1.35); C060: Ba x_Co=0.7-0.8 SRT near RT (x_Zn=1.2-1.3); consistent but Sr-Ba shift not quantified")
add("L12", "mu_spectrum_Co1.2Zn0.8W", status="measured", figure_only=True, locator="Fig. 1", axes="x: frequency 0-20 GHz; y: relative permeability mu' and mu'' 0-6 (one sample, solid lines = model sum DW + NFMR)",
    verbatim="Fig. 1. Spectra of the MP of easy-plane hexaferrites.", sample="Co1.2Zn0.8W (easy plane)", formula="BaCo1.2Zn0.8Fe16O27", definition="mu(f) of polycrystalline easy-plane hexaferrite; non-monotonic with additional maxima/steps",
    group="L12_Fig1", **Z, note="no numeric values in text; digitisation candidate (stage 2 part 2): provides mu'(f) at ~1-20 GHz for easy-plane W; fit parameters (damping, H_phi) not tabulated")
add("L12", "mu_spectrum_Co0.6Zn1.4W", status="measured", figure_only=True, locator="Fig. 2", axes="x: frequency 0-20 GHz; y: relative permeability mu', mu'' 0-5",
    verbatim="Fig. 2. Spectra of the MP of easy-axis hexaferrites.", sample="Co0.6Zn1.4W (easy axis)", formula="BaCo0.6Zn1.4Fe16O27", definition="mu(f) of polycrystalline easy-axis hexaferrite",
    group="L12_Fig2", **Z, note="authors: agreement of independent-grain model better for easy-plane than easy-axis (intergrain magnetic interaction)")
add("L12", "model_independent_grains_ellipsoid", status="assumed_by_authors", qualitative="mu_av = (1/3) int f(N_l)(mu_xx+mu_yy+mu_zz) dN_l, N_l from 1/3 to 1; resonance w_res = sqrt(w1 w2) with w1,w2 in Tables 1,2 (H_theta, H_phi for easy-plane; H_a = 2k1/M0 for easy-axis); M0 -> M0(1-p) porosity correction; fitted: damping constant and H_phi (easy-plane only)",
    locator="Section 2 Theory, Eqs. (1)-(3), Tables 1-2; Section 4", verbatim="The only parameters that varied when constructing theoretical curves were the damping constant and the value of the anisotropy field H_phi for easy-plane hexaferrites.",
    sample="polycrystalline hexaferrites", formula="generic", definition="domain-structure NFMR model; formulas for w1,w2 in Tables 1-2 are partly garbled in the PDF text (symbols lost): need image check before implementation",
    group="L12_model", **Z, note="gives definition of H_theta=2(k1+2k2+3k3-3k4)/M0 and H_phi=36 k4/M0; k_i anisotropy constants (units CGS implied); Kramers-Kronig smoothing of measured mu(f)")

# ------------------------------------------------------------------ C057 LoBue 2012 (Ba single crystal)
S57 = dict(cluster="SATIE", tier=3, transfer=BA, sample="BaCo0.62Zn1.38Fe16O27 single crystal", formula="BaCo0.62Zn1.38Fe16O27",
           synth="flux method: BaCO3, CoO, ZnO, Fe2O3, Na2CO3 in Pt crucible, 1350 C 24 h, cooled to 1000 C at 3 C/h then faster", form="single crystal, m = 31.81 mg, few mm; demagnetising factors N_par=0.4, N_perp=0.3 (estimated)",
           orient="field parallel and perpendicular to c", texture="single crystal", method="direct calorimetry (Peltier heat-flux sensors) of isothermal entropy change, field up to 1 T at 0.033 T/s, 77-300 K")
add("C057", "K1_first_anisotropy_constant_296K", value=8.7e4, unit="J/m^3", kind="energy_density", proposed_kind="energy_density", pending_unit="J/m^3", status="measured",
    locator="Abstract; Section III (p. 6-7)", verbatim="At 296 K one deduces a value of K1 = 8.7 10^4 J m-3 for the first anisotropy constant, in good agreement with literature.",
    definition="K1 in G = K1 sin^2(th) + K2 sin^4(th) + K3 sin^6(th) - mu0 Ms (H_par cos th + H_perp sin th); obtained from K_EP = K1+2K2+3K3 at 296 K assuming K2,K3 negligible at RT (from saturation field H* of entropy curves minus demagnetising field)",
    T={"K": 296}, field="up to 1 T", group="C057_entropy", **S57, note="uncertainty not given; Ba W with x_Co = 0.62 is uniaxial at RT (SRT 215 K)")
add("C057", "spin_reorientation_temperature_range_lo", value=200.0, unit="K", kind="temperature", status="measured", locator="Abstract; Introduction", verbatim="The spin reorientation transition temperature is estimated to take place between 200 and 220 K.",
    definition="intermediate T range where the SRT-related entropy change depends on T (easy-cone intermediate phase suggested: no latent heat in DSC); low bound", group="C057_entropy", T={"K": None}, **S57)
add("C057", "spin_reorientation_temperature_range_hi", value=220.0, unit="K", kind="temperature", status="measured", locator="Abstract; Introduction", verbatim="between 200 and 220 K ... the transition takes place at about 215 K",
    definition="upper bound; central 'Tsr ~ 215 K' from K(T) = 0.175 (T - 215) J/kg sign change", group="C057_entropy", T={"K": None}, **S57)
add("C057", "dMs_dT", value=802.5, unit="A/(m*K)", kind="dMs_dT", proposed_kind="magnetization_temperature_derivative", pending_unit="A/(m*K)", status="measured", locator="Section III, after Eq. (3)", verbatim="we obtain dMs/dT = 802.5 Am-1 (units as printed; per kelvin implied)",
    definition="fit of linear part of Delta s(Ha) at saturating field; unit as printed 'A m-1' (K^-1 missing in text)", group="C057_entropy", T={"K": None, "note": "temperature range of fit not stated"}, **S57, note="unit ambiguity: printed 'Am-1', must be A m-1 K-1")
add("C057", "max_entropy_change_at_SRT", value=0.18, unit="J/(kg*K)", kind="entropy_change", proposed_kind="specific_entropy", pending_unit="J/(kg*K)", status="measured", locator="Section II, last paragraph before Fig. 3", verbatim="The maximum entropy change at the reorientation is s = 0.18 JKg-1K-1.",
    definition="maximum isothermal entropy change associated with the field-driven SRT", group="C057_entropy", T={"K": None}, **S57)
add("C057", "Anisotropy_vs_T_Fig", status="measured", figure_only=True, locator="Fig. 1-3 (Delta s(Ha) curves; inset Fig. 3: K(T))", axes="Fig. 3: x temperature 80-300 K; y entropy change (J kg-1 K-1); inset: effective anisotropy K (J/kg) vs T from H* (parallel/perpendicular)",
    verbatim="K(T) = 0.175 (T - 215) J kg-1", definition="single effective anisotropy constant changing sign at 215 K reproduces calorimetry", group="C057_entropy", T={"K": None, "range": "77-300 K"}, **S57,
    note="K(T) curve (inset) is the only T-dependent anisotropy in this paper; digitisation candidate")

# ------------------------------------------------------------------ C060 Pasko 2011 (Ba, nano)
S60 = dict(cluster="SATIE", tier=3, transfer=BA, synth="sol-gel citrate precursor (Fe(NO3)3 -> Fe(OH)3 -> citrate; Ba(OH)2, CoO, ZnO); gel dried 150-170 C, 450 C 2 h, calcined up to 1300 C 2 h at 200 K/h (BaCo2W: W-ferrite major only at 1300 C)",
           form="loose powder (plate-like particles 50-2000 nm, laser-scattering mean ~700 nm for x=0.75)", texture="random loose powder; XRD at 1300 C shows plate texture in holder")
add("C060", "spin_reorientation_temperature_x0.8", value=55.0 + 273.15, unit="K", kind="temperature", status="measured", locator="Section 3.3, Fig. 5(a)", verbatim="both spin reorientation (at ~ 55 C) and ferrimagnetic-paramagnetic (at ~ 430 C) transitions are well resolved",
    sample="BaCoxZn2-xFe16O27 x=0.8", formula="BaCo0.8Zn1.2Fe16O27", definition="SRT from sharp drop of M(T) on first heating (TGA with permanent magnet); field 1-10 mT; hysteretic heating/cooling", value_original="~55 degC", unit_original="degC", transformation="degC -> K +273.15",
    method="thermogravimetric analyzer with magnet, loose powder, 1-10 mT", field="1-10 mT", T={"K": None, "range": "RT-600 C"}, group="C060_x0.8", **S60, unc_note="'~' approximate; not given")
add("C060", "Curie_temperature_x0.8", value=430.0 + 273.15, unit="K", kind="temperature", status="measured", locator="Section 3.3, Fig. 5(a)", verbatim="ferrimagnetic-paramagnetic (at ~ 430 C)",
    sample="BaCoxZn2-xFe16O27 x=0.8", formula="BaCo0.8Zn1.2Fe16O27", definition="TC from thermomagnetic curve", value_original="~430 degC", unit_original="degC", transformation="degC -> K +273.15", method="TGA with magnet", field="1-10 mT", T={"K": None}, group="C060_x0.8", **S60)
add("C060", "spin_reorientation_T_x0.7_chi_real_max", value=6.0 + 273.15, unit="K", kind="temperature", status="measured", locator="Section 3.3, Fig. 5(b)", verbatim="The real part of susceptibility has a maximum at ~ 6 C followed by a sharp drop attributed to SRT. Its imaginary part has a minimum at ~ 27 C which can also be taken as the spin reorientation temperature.",
    sample="BaCoxZn2-xFe16O27 x=0.7", formula="BaCo0.7Zn1.3Fe16O27", definition="AC susceptibility (PPMS): max of chi' at ~6 C; chi'' minimum at ~27 C (300 K) also taken as SRT (two definitions differ by 21 K)", value_original="~6 degC", unit_original="degC", transformation="degC -> K +273.15",
    method="AC susceptibility, PPMS", field="AC (amplitude/frequency not stated)", T={"K": None, "range": "about 23-300 K"}, group="C060_x0.7", **S60)
add("C060", "spin_reorientation_T_x0.7_chi_imag_min", value=27.0 + 273.15, unit="K", kind="temperature", status="measured", locator="Section 3.3, Fig. 5(b)", verbatim="Its imaginary part has a minimum at ~ 27 C",
    sample="BaCoxZn2-xFe16O27 x=0.7", formula="BaCo0.7Zn1.3Fe16O27", definition="alternative SRT definition", value_original="~27 degC", unit_original="degC", transformation="degC -> K +273.15", method="AC susceptibility, PPMS", T={"K": None}, group="C060_x0.7", **S60)
add("C060", "coercivity_BaCo2W_1300C", value=15.0, unit="mT", kind="mu0H", status="measured", locator="Section 3.2, after Fig. 3", verbatim="the coercivity decreases and reaches ~ 15 mT at 1300 C, a typical for this phase value",
    sample="BaCo2Fe16O27 powders A,B,C", formula="BaCo2Fe16O27", definition="coercivity (mu0Hc) of W-phase powder from VSM", method="VSM (LakeShore)", T={"K": 300}, group="C060_Hc", **S60, note="approximate '~'")

# ------------------------------------------------------------------ C056 Qin 2013 (Ba(MnZn)xCo2(1-x)W) -- Mn!
S56 = dict(cluster="QIN_C056", tier=3, transfer=BA + "; Mn-Zn co-substitution on Co2W", sample="Ba(MnZn)0.4Co1.2Fe16O27", formula="BaMn0.4Zn0.4Co1.2Fe16O27",
           synth="sol-gel: Fe3+, Co2+, Zn2+, Ba2+ nitrates + Mn acetate, citric acid 1:1.5, NH3 to pH 6-7, 80 C 4 h sol, dried 100 C, ashed 450 C, sintered 1235 C 4 h in electric furnace (air)",
           form="powder (hexagonal flakes, mean size 5 um); EM measured in toroid of powder:paraffin 8:3 by weight", texture="random (paraffin composite)")
add("C056", "Ms_RT_x0.4", value=67.320, unit="A*m^2/kg", kind="specific_magnetization", proposed_kind="specific_magnetization", pending_unit="A*m^2/kg", status="measured", locator="Section 3.1, Fig. 3", verbatim="The coercitity (Hc), magnetization (Ms) and retentivity (Mr) are 28.702 Oe, 67.320 emu/g, 3.7964 emu/g, respectively.",
    definition="'magnetization' from VSM hysteresis loop (saturation definition/field not stated)", method="VSM", T={"K": 300, "note": "RT assumed, not stated"}, group="C056_VSM", transformation="1 emu/g = 1 A m^2/kg", **S56, note="Mn-containing W: only Ba, Mn:Zn=1:1 co-doping")
add("C056", "Hc_x0.4", value=28.702, unit="Oe", kind="H_field", status="measured", locator="Section 3.1, Fig. 3", verbatim="Hc ... 28.702 Oe", definition="coercivity", method="VSM", T={"K": 300}, group="C056_VSM", **S56)
add("C056", "permittivity_real_range", value=5.5, unit="1", kind="dimensionless", status="measured", locator="Section 3.3, Fig. 6", verbatim="The imaginary part e'' of the permittivity is small, and the real part e' is in the range of 5.5 - 6.5 over the frequency range of 2 - 18 GHz.",
    definition="eps' lower bound of range 5.5-6.5 (2-18 GHz) of the powder:paraffin composite (NOT intrinsic ferrite eps)", method="AV3629 vector network analyser, transmission/reflection on annular sample", field="none", T={"K": 300}, group="C056_EM", **S56)
add("C056", "permittivity_real_range_hi", value=6.5, unit="1", kind="dimensionless", status="measured", locator="Section 3.3, Fig. 6", verbatim="e' is in the range of 5.5 - 6.5 over the frequency range of 2 - 18 GHz", definition="upper bound", method="VNA", T={"K": 300}, group="C056_EM", **S56)
add("C056", "mu_imag_peak_1", value=4.5, unit="GHz", kind="frequency", status="measured", locator="Section 3.3, Fig. 6", verbatim="Two peaks of the imaginary part m'' are located at 4.5 GHz and 8 GHz",
    definition="peaks of mu'' (domain-wall resonance ~4.5 GHz per authors, natural resonance ~8 GHz); composite, no demag/porosity correction", method="VNA 2-18 GHz", T={"K": 300}, group="C056_EM", **S56)
add("C056", "mu_imag_peak_2", value=8.0, unit="GHz", kind="frequency", status="measured", locator="Section 3.3, Fig. 6/7", verbatim="Two peaks of the imaginary part m'' are located at 4.5 GHz and 8 GHz",
    definition="second mu'' peak ('natural resonance' per authors)", method="VNA 2-18 GHz", T={"K": 300}, group="C056_EM", **S56)
add("C056", "reflection_loss_min", value=-40.7, unit="dB", kind="reflection_loss", proposed_kind="reflection_loss", pending_unit="dB", status="measured", locator="Abstract; Section 3.2 (d = 2.8 mm)", verbatim="the reflection loss peak was -40.7 dB at a frequency of 7.3 GHz, with a bandwidth of 6.6 GHz at a loss of less than -10 dB",
    definition="calculated from measured eps, mu (transmission-line formula) for 2.8 mm layer on metal; not an intrinsic material property", method="derived from VNA data", T={"K": 300}, group="C056_EM", **S56)
add("C056", "tan_delta_m_peak_values", value=0.82, unit="1", kind="dimensionless", status="measured", locator="Section 3.3, Fig. 7", verbatim="the first peak exhibiting a value of 0.78 at approximately 4.5 GHz and the second peak exhibiting a value of 0.82 at approximately 8 GHz",
    definition="magnetic loss tangent mu''/mu' peak at ~8 GHz (second peak; first peak 0.78 at 4.5 GHz)", method="VNA", T={"K": 300}, group="C056_EM", **S56)

# ------------------------------------------------------------------ C036 (Ba1-x(La0.5Sr0.5)xCo2W) only x=0 used
S36 = dict(cluster="LI_C036", tier=3, transfer=BA, sample="x = 0.0 (BaCo2Fe16O27)", formula="BaCo2Fe16O27",
           synth="solid-state: BaCO3, SrCO3, Fe2O3, Co3O4, La2O3; ball mill ZrO2 24 h; dried 100 C 5 h vacuum; calcined 1270 C 2 h; milled 16 h; PVA 8 wt %; pressed; sintered 1270 C 2 h (5 C/min)",
           form="sintered toroid (7 x 3.04 mm, 2 mm thick for coax measurement), VSM on powders", texture="random")
add("C036", "Ms_BaCo2W", value=73.1, unit="A*m^2/kg", kind="specific_magnetization", proposed_kind="specific_magnetization", pending_unit="A*m^2/kg", status="measured", locator="Section 3.3 text and Table 1", verbatim="X = 0.0   Ms (emu/g) 73.1   Hc (Oe) 130",
    definition="Ms of pure BaCo2W from VSM RT hysteresis", method="VSM LakeShore 7404", T={"K": 300, "note": "RT"}, group="C036_x0", transformation="1 emu/g = 1 A m^2/kg", **S36)
add("C036", "Hc_BaCo2W", value=130.0, unit="Oe", kind="H_field", status="measured", locator="Table 1", verbatim="X = 0.0  Ms 73.1  Hc 130", definition="coercivity", method="VSM", T={"K": 300}, group="C036_x0", **S36)
add("C036", "mu_imag_resonance_band_Co2W_to_x0.6", value=6.0, unit="GHz", kind="frequency", status="measured", locator="Section 3.4, Fig. 6(b)", verbatim="the resonance peaks of samples with x = 0.0-0.6 are located in the frequency range of 6-10 GHz, and the peak frequency shifts towards lower frequencies with increasing doping amount",
    definition="lower end of the range (6-10 GHz) of mu'' resonance peaks for x=0-0.6 (BaCo2W is the x=0 member; exact x=0 peak not given in text -> figure)", method="coaxial-line VNA 1-18 GHz (Agilent E5071C)", T={"K": 300}, group="C036_EM", **S36,
    note="range covers La-Sr co-doped members; value for x=0 alone needs Fig. 6 digitisation")
add("C036", "reflection_loss_min_Co2W", value=-29.01, unit="dB", kind="reflection_loss", proposed_kind="reflection_loss", pending_unit="dB", status="measured", locator="Section 3.4, Fig. 8", verbatim="undoped Co2W hexaferrite ... RLmin of -29.01 dB at a thickness of 1.1 mm and an effective absorption bandwidth reaching 6.28 GHz",
    definition="calculated from measured eps/mu with transmission-line theory", method="derived from VNA data", T={"K": 300}, group="C036_EM", **S36)

# ------------------------------------------------------------------ C044 Mahmood 2020 (BaMg2-xCoxW)
S44 = dict(cluster="MAHMOOD_C044", tier=3, transfer=BA, synth="high-energy ball milling (16 h, 250 rpm, acetone) of BaCO3, Fe2O3, CoO, MgO; pellets 1.2 cm x 2 mm, 5 t; sintered 1100/1200/1300 C 2 h in air (W phase single only at 1300 C)",
           form="sintered pellet; bulk density 4.72-4.82 g/cm3 (porosity 6.4-8.5 %)", texture="not stated", T={"K": 300})
for x, f, ms, hc, tc in ((0.0, "BaMg2Fe16O27", "59.44", "70", "452"), (1.0, "BaMgCoFe16O27", "65.22", "30", "473"), (2.0, "BaCo2Fe16O27", "72.56", "90", "483")):
    add("C044", f"Ms_1300C_x{x}", vs=ms, unit="A*m^2/kg", kind="specific_magnetization", proposed_kind="specific_magnetization", pending_unit="A*m^2/kg", status="measured", locator="Abstract; Section 3.3.1, Fig. 12",
        verbatim=f"x = {x}: saturation magnetization {ms} emu/g (59.44 -> 65.22 -> 72.56 emu/g for x = 0, 1, 2)", sample=f"BaMg2-xCoxFe16O27 x={x}", formula=f, definition="Ms from VSM hysteresis (up to 10 kOe) of 1300 C sample", method="VSM", group=f"C044_x{x}", **S44, transformation="1 emu/g = 1 A m^2/kg",
        note="x=1 almost saturated at 10 kOe -> low anisotropy")
    add("C044", f"Hc_1300C_x{x}", value=float(hc), unit="Oe", kind="H_field", status="measured", locator="Section 3.3.1 (text: 70 Oe for Mg2W, 30 Oe MgCo-W, 90 Oe Co2W)", verbatim=f"low coercivity values of 70 Oe for Mg2W, 30 Oe for MgCo-W, and 90 Oe for Co2W", sample=f"x={x}", formula=f,
        definition="coercivity (multidomain)", method="VSM", group=f"C044_x{x}", **S44)
    add("C044", f"Curie_temperature_x{x}", value=float(tc) + 273.15, unit="K", kind="temperature", status="measured", locator="Table 6; Section 3.3.2, Fig. 13", verbatim=f"Tc (C) 452 473 483 (x = 0.0, 1.0, 2.0)", sample=f"x={x}", formula=f,
        definition="dip of dM/dT in M(T) at 100 Oe", value_original=f"{tc} degC", unit_original="degC", transformation="degC -> K +273.15", method="thermomagnetic M(T), H = 100 Oe", field="100 Oe", group=f"C044_x{x}", **S44)
add("C044", "SRT_T_Co2W_cone_to_axis", value=304.0 + 273.15, unit="K", kind="temperature", status="assumed_by_authors", locator="Section 3.3.2", verbatim="The rise in the thermomagnetic curve of the sample with x = 2.0 into a peak at ~ 304 C is associated with the spin reorientation transition from conical to c-axis spin orientation.",
    sample="BaCo2Fe16O27", formula="BaCo2Fe16O27", definition="peak in M(T) at 100 Oe interpreted as cone->axis SRT (not verified by diffraction)", value_original="~304 degC", unit_original="degC", transformation="degC -> K +273.15", method="M(T), 100 Oe", field="100 Oe", group="C044_x2", **S44)
add("C044", "SRT_T_MgCoW", value=232.0 + 273.15, unit="K", kind="temperature", status="assumed_by_authors", locator="Section 3.3.2", verbatim="The weak peak at 232 C could be associated with spin reorientation transition, where the magnetic anisotropy changes from easy-cone to easy axis.",
    sample="BaMgCoFe16O27", formula="BaMgCoFe16O27", definition="weak derivative dip at 232 C ('could be associated')", value_original="232 degC", unit_original="degC", transformation="degC -> K +273.15", method="M(T), 100 Oe", field="100 Oe", group="C044_x1", **S44)
for x, f, a, c in ((0.0, "BaMg2Fe16O27", "5.9076", "32.9681"), (1.0, "BaMgCoFe16O27", "5.9061", "32.9565"), (2.0, "BaCo2Fe16O27", "5.9069", "32.9552")):
    add("C044", f"lattice_a_x{x}", vs=a, unit="angstrom", kind="length", status="measured", locator="Table 2 (Rietveld, FullProf)", verbatim=f"{f}  a = b {a}  c {c}", sample=f"x={x}", formula=f, definition="a=b, P63/mmc", method="XRD Rietveld (Shimadzu 7000)", group=f"C044_x{x}", **S44)
    add("C044", f"lattice_c_x{x}", vs=c, unit="angstrom", kind="length", status="measured", locator="Table 2", verbatim=f"{f}  a = b {a}  c {c}", sample=f"x={x}", formula=f, definition="c, P63/mmc", method="XRD Rietveld", group=f"C044_x{x}", **S44)
add("C044", "literature_Ms_SrCo2W_cited", value=74.31, unit="A*m^2/kg", kind="specific_magnetization", proposed_kind="specific_magnetization", pending_unit="A*m^2/kg", status="measured", locator="Section 3.3.1, last paragraph (their ref [55])", verbatim="in good agreement with the previously reported value of 71.18 emu/g for BaCo2Fe16O27 [41], and 74.31 emu/g for SrCo2Fe16O27 [55]",
    sample="SrCo2Fe16O27 (cited)", formula="SrCo2Fe16O27", definition="Ms of SrCo2W quoted from another paper (primary not in corpus)", cluster="SECONDARY", tier=1, split="excluded_secondary_citation", note="secondary citation; compare L01 78.5 and L02 (no Ms)", group="C044_cited")

# ------------------------------------------------------------------ C052 Koide 2015 (Fe2W endmember)
add("C052", "Curie_temperature_BaFe2W", value=503.0 + 273.15, unit="K", kind="temperature", status="measured", locator="Abstract", verbatim="the Curie temperatures of Fe2W and Fe2Y ferrites are 503 and 375 C, respectively",
    sample="Fe2W", formula="BaFe18O27", definition="TC of BaFe2+2Fe3+16O27", value_original="503 degC", unit_original="degC", transformation="degC -> K +273.15", cluster="KOIDE_C052", tier=3, transfer=BA,
    synth="BaFe2O4 + Fe3O4 + alpha-Fe2O3, ball mill 24 h, pellets sealed in evacuated silica tube (0.5 Pa), sintered 1300 C", form="toroid pellet", T={"K": None}, group="C052_Fe2W",
    note="authors state literature TC of Fe2W 'inconsistent' (refs 7-10); Ba-based Fe2+ endmember")

# ------------------------------------------------------------------ C023 Li 2014 (SrCo2W ceramic)
S23 = dict(cluster="LI_C023", tier=1, sample="WM (major phase SrCo2W, 70-90 wt %, CoFe2O4 impurity)", formula="SrCo2Fe16O27 (+ CoFe2O4)",
           synth="solid-state: SrCO3, Co2O3, Fe2O3 (99.5+ %) stoichiometric, calcined in air 1250 C 4 h, pellets 12x3 mm, sintered in OXYGEN 1250 C 8 h (WM)", form="ceramic disc/plate (10 x 1 mm polished), Ag electrodes", texture="not stated; random ceramic", density="not stated",
           method="VSM (field 100 Oe, T from top down)", field="100 Oe", T={"K": None, "range": "up to > 763 K"})
add("C023", "magnetocrystalline_anisotropy_transition_T", value=470.0, unit="K", kind="temperature", status="measured", locator="Abstract; Section III.B, Fig. 3(b); Conclusions", verbatim="For the WM specimen (SrCo2W), the M-T curve exhibits one drop at ~470 K (Fig. 3b).",
    definition="steep M decrease at 470-520 K assigned by authors to cone -> uniaxial (Cone1 <470 K, Cone2 470-520 K, uniaxial 520-763 K, PM >763 K), by analogy with BaCo2W (refs); NOT verified by diffraction; compare L02: planar(Cm'cm') with TT = 555 K (derivative) / 520 K (onset)", group="C023_MT", **S23,
    note="CONFLICT/consistency: C023 470 K vs L02 520-555 K for the nominally same SrCo2Fe16O27; sample 70-90 wt % W; different sintering atmosphere (O2 here)")
add("C023", "assumed_Curie_temperature_for_regions", value=763.0, unit="K", kind="temperature", status="assumed_by_authors", locator="Section III.B (text before Fig. 3)", verbatim="Below the Curie temperature (Tc,FM-PM = 763 K; 728 K), BaCo2Fe16O27 (BaCo2W) goes through the following major changes in MCA ... Region iv: PM phase (> 763 K)",
    definition="Ba TC (literature values 763 K; 728 K) used as boundary for the SrCo2W regions; NOT measured for the Sr specimen in this paper", group="C023_MT", **S23, note="not a Sr measurement; L02 gives TC(SrCo2W) = 790 K")
add("C023", "cone_angle_BaCo2W_300K_cited", value=70.0, unit="deg", kind="angle", status="measured", locator="Section III.B", verbatim="Cone 1 (~70 deg, 2-453 K; 68.5-70 deg, 300 K)", definition="cone half-angle of BaCo2W from refs 5,25,26 (secondary, Ba)", sample="BaCo2W (cited)", formula="BaCo2Fe16O27",
    cluster="SECONDARY", tier=3, split="excluded_secondary_citation", group="C023_cited", note="L01 found planar (90 deg) for SrCo2W; cone angle for Sr unknown")

# ------------------------------------------------------------------ C015 Sr Co2 W/Z composite
S15 = dict(cluster="JU_C015", tier=1, sample="W/Z composite (80 vol % W, 20 vol % Z), off-stoichiometric", formula="SrCo2Fe16O27 + Sr3Co2Fe24O41",
           synth="solid-state: SrCO3, Co3O4, Fe2O3 (>99.9 %), ball milled in ethanol, calcined in air 1100 C 16 h, pellets sintered 1200 C 16 h in air and quenched; stoichiometry varied by treating 1200-1400 C and under flowing O2", form="sintered pellet 5 mm x 0.5 mm, grain size 0.5 um (SEM)", texture="random",
           method="XRD (FullProf multiphase refinement), RT", T={"K": 300, "note": "RT"})
add("C015", "lattice_a_W_in_composite", vs="5.849", unit="angstrom", kind="length", status="measured", locator="Section III (first paragraph of results)", verbatim="lattice parameters a = b = 5.849, c = 33.112 A for W-type and a = b = 5.841, c = 51.949 A for Z-type systems",
    definition="W-type phase (P63/mmc) in 80:20 W:Z composite", group="C015_XRD", **S15, note="a is 0.045 A smaller and c 0.37 A larger than L01/L02 SrCo2W (5.895 / 32.74-32.76): anomalous; composite/off-stoichiometric (cation vacancies/oxygen excess claimed) -> low weight")
add("C015", "lattice_c_W_in_composite", vs="33.112", unit="angstrom", kind="length", status="measured", locator="Section III", verbatim="a = b = 5.849, c = 33.112 A for W-type", definition="W-type c", group="C015_XRD", **S15)
add("C015", "W_phase_volume_fraction", value=80.0, unit="%", kind="fraction", status="measured", locator="Section II and III", verbatim="The multiphase fitting yields the volume fraction of the respective phases to be 80% W-type and 20% Z-type.", definition="volume fraction from FullProf multiphase fit", group="C015_XRD", **S15)
add("C015", "cited_TC_and_cone_SrCo2W", status="assumed_by_authors", qualitative="TC ~ 850 K; conical magnetic structure below ~550 K for SrCo2Fe16O27 (their ref 14)", locator="Section I, para 7", verbatim="It exhibits ferromagnetic order below TC ~ 850 K and below T ~ 550 K, conical magnetic structure sets in.",
    sample="SrCo2Fe16O27 (cited)", formula="SrCo2Fe16O27", definition="quoted from ref [14], not measured here", cluster="SECONDARY", tier=1, split="excluded_secondary_citation", group="C015_cited", note="conflicts with L02 (TC 790 K; planar) - flagged in CONFLICTS")

# ------------------------------------------------------------------ C066 Suwa 1996 (qualitative only)
add("C066", "Sr_W_decomposition_vs_Ba", status="measured", qualitative="Fe2+ content of calcined compound increases with Sr; Sr[ZnLi0.5Fe0.5]W partially decomposes to M + spinel on annealing at 900-1100 C, Ba[ZnLi0.5Fe0.5]W stays single phase; K_A and H_A ~60 % of M-type; grain size 0.3-5.0 um; sinter 1175-1225 C",
    locator="Synopsis (English)", verbatim="W phase of Sr[ZnLi0.5Fe0.5]Fe16O27 partially decomposed into M phase and spinel phase. On the other hand, Ba[ZnLi0.5Fe0.5]Fe16O27 remained a single phase.",
    sample="Sr1-xBax[ZnLi0.5Fe0.5]Fe16O27", formula="Sr1-xBax(ZnLi0.5Fe0.5)Fe16O27", definition="all numeric dependences (as, Hca, Curie point vs x) are figures only in a scanned paper", cluster="SUWA_C066", tier=2, group="C066", figure_only=False,
    note="Li-Zn charge-compensated W (not Co/Zn/Mn); numeric content only in figures: Fig. 2-5; text layer partial")

# ------------------------------------------------------------------ C053 (low quality)
add("C053", "Ms_range_CaSr2_MnZn_W", status="measured", qualitative="Ms 13.97-40.45 emu/g across x = 0-1 (CaSr2(MnZn)x/2Fe16-xO27), sintered 1000 C 4 h", locator="Section on VSM / Table-4", verbatim="the saturation magnetization varies in the range of 13.97 emu/g to 40.45 emu/g",
    sample="CaSr2(MnZn)x/2Fe16-xO27 x=0-1", formula="CaSr2(MnZn)x/2Fe16-xO27", definition="formula as printed (Ca/Sr co-occupancy unclear; sintered at 1000 C - far below W-phase formation T of other sources); table text layer scrambled",
    cluster="KHOB_C053", tier=3, group="C053", note="EXCLUDED: internally inconsistent formula/synthesis, Ms far below any other W; listed only because it mentions Mn in a Sr-containing W")

# ------------------------------------------------------------------ methodology tier 4
M115 = dict(cluster="METHOD", tier=4, sample="Ba1-xLaxFe12O19 (M-type!) x = 0, 0.5, 1.0", formula="Ba1-xLaxFe12O19", synth="ceramic: BaCO3, La2O3, Fe2O3 + 3 wt % Bi2O3, calc 800 C 2 h, sintered 1050 C 2 h", T={"K": 0})
add("C115", "exchange_integrals_from_spin_configurations", status="calculated", qualitative="Heisenberg model; J_ij from total-energy differences between ferrimagnetic spin configurations (S_i = +-1 in the energy sum, n_i sublattice population, z_ij neighbours); neighbours beyond 4 A neglected; MFA gives upper and RPA lower bound of TC",
    locator="Section 'Calculation of exchange integrals' (Eqs. 5-9 in image; text only here)", verbatim="The exchange energy per unit cell in the complex system with N magnetic sublattices could be then written as ... where ni is the number of ith sublattice, zij is the number of neighbouring sites in jth sublattice to ith sublattice, Jij is the exchange integrals, Si and Sj represent the spins in the ith and jth sublattices, and are equal to +-1",
    definition="exchange-sign convention NOT recoverable from text (equation is an image) -> do not import numbers; methodology only", group="C115", **M115, note="M-type, Ba: numbers not transferable; sign/normalisation to be fixed by our own DFT")
add("C115", "Ueff_reproducing_TC_BaM", value=6.7, unit="eV", kind="energy", status="calculated", locator="Section 'Calculations of Curie temperature', Table 4, Conclusions", verbatim="the experimental Curie temperature is reproduced by calculations for the effective value Ueff ~ 6.7 eV, and the MFA and RPA estimations provide an upper and lower bound of the Curie temperature",
    definition="U_eff that matches exp TC with RPA for La-substituted BaM (GGA+U VASP PBE 500 eV, 11x11x1, forces<0.03 eV/A); methodological illustration that J(U) decays strongly with U", group="C115", **M115, note="tier 4: DO NOT use for W-phase; shows U-sensitivity")
add("C112", "exchange_J_from_22_spin_configurations", status="calculated", qualitative="SrM DFT+U: J_ij from up to 22 spin configurations, U varied 2-5 eV; J decays strongly with U; SOC changes largest J by ~0.2 %; expanded theoretical cell gives ~10 % smaller J than experimental cell; ES/TS/ED/TD static-vs-dynamic relaxation limits",
    locator="Section 'Exchange couplings' and Figs. 4-7", verbatim="The SFO exchange constants have been obtained from Eq. (2) employing up to 22 different spin configurations and systematically varying the U term in the 2-5 eV interval ... for the largest Js the difference was only ~0.2 %, indicating that SOC may be safely ignored",
    sample="SrFe12O19 (M-type)", formula="SrFe12O19", definition="method/sensitivity findings only (numbers not transferable to W)", cluster="METHOD", tier=4, T={"K": 0}, synth="n/a (DFT)", group="C112",
    note="supports our plan: J from spin-configuration fits, SOC-free; check cell (exp vs relaxed) and U sensitivity")
