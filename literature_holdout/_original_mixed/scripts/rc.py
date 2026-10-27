"""Record-building core for literature/records.yaml. Python source of truth; run build_records.py to emit YAML/JSONL.
Not part of project/; reads only literature inputs already extracted by hand/script."""
import re

REG = []          # all records (dicts)
_COUNT = {}

LIT = {  # lit_id -> bibliographic info
 "L01": dict(citation="Morch M.I., Ahlburg J.V., Saura-Muzquiz M., Eikeland A.Z., Christensen M. Structure and magnetic properties of W-type hexaferrites. IUCrJ 6 (2019) 492-499", doi="10.1107/S2052252519003130", year=2019, path="06_literature/L01_Morch2019", access="full_text"),
 "L02": dict(citation="Morch M.I., Christensen M. Controlling the magnetic structure in W-type hexaferrites. J. Appl. Cryst. 56 (2023)", doi="10.1107/S1600576723002133", year=2023, path="06_literature/L02_Morch2023", access="full_text"),
 "L03": dict(citation="Nakai S., Waki T., Tabata Y., Kato M., Ohta H., Nakamura H. Study on Factors Hindering the Single-Phase Formation of Divalent-Ion-Stabilized W-Type Ferrites. Mater. Trans. 64(5) (2023) 1072-1077", doi="10.2320/matertrans.MT-Y2022010", year=2023, path="06_literature/L03_Nakai2023", access="full_text"),
 "L05": dict(citation="Islam R., Madsen S.P., Christensen M. Enhancing magnetocrystalline anisotropy in W-type SrFe18O27 hexaferrite via tailored substitutions. Acta Mater. 288 (2025) 120828", doi="10.1016/j.actamat.2025.120828", year=2025, path="06_literature/L05_Islam2025", access="full_text", license="CC BY 4.0"),
 "L11": dict(citation="He S., Liu R., Liu X., Ye X., Wang L., Shen B. Influences of divalent ion substitution on the magnetic and dielectric properties of W-type barium ferrite. Chin. Phys. B 33 (2024) 066801", doi="10.1088/1674-1056/ad39d5", year=2024, path="06_literature/L11_He2024", access="full_text"),
 "L12": dict(citation="Zhuravlev V.A. Magnetic Permeability Spectra of Hexaferrites in the Microwave Range. Phys. Met. Metallogr. 100 Suppl. 1 (2005) S79-S81", doi=None, year=2005, path="06_literature/L12_Zhuravlev2005", access="full_text"),
 "C004": dict(citation="Nakai S., Waki T., Tabata Y., Nakamura H. Low oxygen potential synthesis of W-type strontium ferrite SrFe18O27 in a closed atmosphere. J. Ceram. Soc. Japan 133(6) (2025) 235-241", doi="10.2109/jcersj2.25033", year=2025, path=None, access="full_text", license="CC BY 4.0"),
 "C015": dict(citation="Electric-field-driven resistive transition in multiferroic SrCo2Fe16O27/Sr3Co2Fe24O41 composite (arXiv version of 2022 paper)", doi="10.1063/5.0087172", year=2022, path=None, access="full_text (arXiv)"),
 "C023": dict(citation="Li J., Zhang H.-F., Shao G.-Q. et al. Magnetic / dielectric anomalies and magnetodielectric / magnetoelectric effects in Z- and W-type hexaferrites (arXiv:1407.3605)", doi="10.48550/arxiv.1407.3605", year=2014, path=None, access="full_text (arXiv)"),
 "C036": dict(citation="Absorption Band Tunable La-Sr Co-Doped BaCo2-W Type Hexaferrites. Materials 16 (2023) 5897", doi="10.3390/ma16175897", year=2023, path=None, access="full_text (PMC XML)", license="CC BY 4.0"),
 "C044": dict(citation="Mahmood S.H., Al-Shiab Q., Bsoul I., Maswadeh Y., Awadallah A. Structural and magnetic properties of (Mg,Co)2W hexaferrites. Jordan J. Phys. 13(1) (2020)", doi="10.47011/13.1.1", year=2020, path=None, access="full_text (arXiv)"),
 "C052": dict(citation="Koide M., Kakizaki K., Kamishima K. Synthesis and magnetic properties of Fe2W and Fe2Y hexaferrites. J. Magn. Soc. Jpn. 39 (2015) 147-150", doi="10.3379/msjmag.1506r003", year=2015, path=None, access="full_text"),
 "C053": dict(citation="Khobaragade C.L. et al. Improvement of the structural, electrical and magnetic properties of calcium W-type hexaferrite nano sized powders using co-precipitation method. IJRBAT (2015)", doi="10.29369/ijrbat.2015.03.ii.0024", year=2015, path=None, access="full_text"),
 "C056": dict(citation="Qin X., Cheng Y., Zhou K., Huang S., Hui X. Microwave Absorbing Properties of W-Type Hexaferrite Ba(MnZn)xCo2(1-x)Fe16O27. J. Mater. Sci. Chem. Eng. 1 (2013) 8-13", doi="10.4236/msce.2013.14002", year=2013, path=None, access="full_text", license="CC BY"),
 "C057": dict(citation="LoBue M. et al. Direct calorimetric measurements of isothermal entropy change on single crystal W-type hexaferrites at the spin reorientation transition. J. Appl. Phys. 111 (2012) (arXiv:1210.0997)", doi="10.1063/1.3670062", year=2012, path=None, access="full_text (arXiv)"),
 "C060": dict(citation="Pasko A. et al. Magnetic and structural characterization of nanosized BaCoxZn2-xFe16O27 hexaferrite in the vicinity of spin reorientation transition. J. Phys.: Conf. Ser. 303 (2011) 012045", doi="10.1088/1742-6596/303/1/012045", year=2011, path=None, access="full_text (arXiv)"),
 "C066": dict(citation="Suwa K., Taguchi H., Takeishi T. Crystal Stability and Magnetic Properties of W-type Hexagonal Ferrite Sr1-xBax[ZnLi0.5Fe0.5]Fe16O27. J. Jpn. Soc. Powder Powder Metall. 43 (1996) 995", doi="10.2497/jjspm.43.995", year=1996, path=None, access="partial (scan, text partly garbled)"),
 "C112": dict(citation="An ab initio study of the magnetic properties of strontium hexaferrite. Sci. Rep. 11 (2021)", doi="10.1038/s41598-021-81028-7", year=2021, path=None, access="full_text (PMC XML)", license="CC BY"),
 "C115": dict(citation="Wu C. et al. Calculation of exchange integrals and Curie temperature for La-substituted barium hexaferrites. Sci. Rep. 6 (2016) 36200", doi="10.1038/srep36200", year=2016, path=None, access="full_text (PMC XML)", license="CC BY"),
}

# quantity kinds that exist in hexw.units (with the units we use); everything else goes to `pending`
SUPPORTED = {"length", "temperature", "fraction", "dimensionless", "mu0H", "H_field", "frequency", "time", "angle", "energy", "mass"}

SPLIT_ROLE = {"literature_fit": "literature_fit", "literature_holdout": "literature_holdout"}

CLUSTER_SPLIT = {   # cluster -> proposed split (see split_proposal.md)
 "MORCH_L02": "literature_fit",
 "MORCH_L01_MgNi": "literature_fit",
 "MORCH_L01_CoZn": "literature_holdout",
 "ISLAM_L05": "literature_fit",
 "NAKAI": "literature_fit",
 "LI_C023": "literature_holdout",
 "JU_C015": "literature_holdout",
 "HE_L11": "literature_holdout",
 "ZHURAVLEV": "literature_holdout",
 "SATIE": "literature_holdout",
 "MAHMOOD_C044": "literature_holdout",
 "QIN_C056": "literature_holdout",
 "LI_C036": "literature_holdout",
 "KOIDE_C052": "literature_holdout",
 "SUWA_C066": "literature_holdout",
 "KHOB_C053": "excluded_low_quality",
 "SECONDARY": "excluded_secondary_citation",
 "METHOD": "methodology_only",
}


def vu(s):
    """'5.8955(1)' -> (5.8955, 1e-4); '0.5' -> (0.5, None)."""
    s = str(s).strip().replace("−", "-").replace(",", ".")
    m = re.match(r"^(-?\d+\.?\d*)(?:\((\d+)\))?$", s)
    if not m:
        raise ValueError("bad value string " + s)
    v = float(m.group(1)); u = None
    if m.group(2):
        dec = len(m.group(1).split(".")[1]) if "." in m.group(1) else 0
        u = int(m.group(2)) * 10 ** (-dec)
    return v, u


def add(lit, name, *, vs=None, value=None, unc=None, unit=None, kind=None, status, locator, verbatim, sample, formula, cluster, tier,
        group=None, synth=None, T=None, field=None, orient=None, form=None, density=None, texture=None, method=None, definition=None,
        transformation=None, figure_only=False, axes=None, note=None, unit_original=None, value_original=None, proposed_kind=None,
        unc_note=None, qualitative=None, transfer=None, split=None, pending_unit=None):
    n = _COUNT.get(lit, 0) + 1
    _COUNT[lit] = n
    rid = f"{lit}-{n:03d}"
    if vs is not None:
        value, unc = vu(vs)
        value_original = value_original or str(vs)
    sp = split or CLUSTER_SPLIT[cluster]
    role = SPLIT_ROLE.get(sp)
    est = None
    pend = None
    numeric = value is not None and not figure_only
    if numeric and kind in SUPPORTED and role:
        est = dict(name=name, value=float(value), unit=unit, quantity_kind=kind, source_role=role,
                   provenance=dict(source_id=lit, locator=locator, note=(note or None),
                                   status="derived_by_us" if status == "derived_by_us" else "verified"),
                   uncertainty=(dict(kind="std", std=float(unc)) if unc is not None else dict(kind="unknown")), identifiable=True)
    elif numeric:
        pend = dict(name=name, value=float(value), unit=unit or pending_unit, proposed_quantity_kind=proposed_kind or kind,
                    reason=("schema Kind/unit not available" if role else f"split={sp}: kept as claim, not a parameter estimate"),
                    uncertainty_std=(float(unc) if unc is not None else None))
    r = dict(rec_id=rid, lit_id=lit, doi=LIT[lit]["doi"], name=name, estimate=est, pending=pend,
             ext=dict(locator=locator, sample_id_in_paper=sample, formula=formula, synthesis=synth,
                      measurement=dict(T=T, field=field, orientation=orient, form=form, density=density, texture=texture),
                      method=method, definition=definition, status=status, verbatim=verbatim, transformation=transformation,
                      measurement_group=group, figure_only=bool(figure_only), axes=axes, split=sp, split_cluster=cluster, tier=tier,
                      transfer_note=transfer, notes=note, value_original=value_original, unit_original=unit_original or unit,
                      uncertainty_note=unc_note, qualitative=qualitative))
    REG.append(r)
    return r
