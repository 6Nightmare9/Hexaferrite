"""Stage 3 round 1: leaks holdout -> fit (incident 2026-10-26), routing by source file, the final reference check,
derivatives of holdout sources, schema details. Synthetic data for the logic; the built literature/ for the content."""
from __future__ import annotations

import csv
import io
import json
from pathlib import Path

import pytest
import yaml

from hexw.ingest import literature_split as LS
from hexw.io import paths
from hexw.registry import schemas as S
from hexw.registry import versioning as V

REAL = Path(__file__).resolve().parents[2]
CIF = "L01_Morch2019/supporting/lt5018sup_cif/{s}/{s}_SINQ_HRPT.cif"
HDR = ["a", "block", "file", "label", "occupancy"]


def row(sample, block, label="Fe1T"):
    return ["8.37957(18)", block, CIF.format(s=sample), label, "0.69(5)"]


# ---------------------------------------------------------------- routing by source file
def test_rows_of_holdout_cifs_go_to_holdout_even_for_spinel_blocks():
    body = [row("SrCo2Fe16O27", "SrCo2Fe16O27"), row("SrCo2Fe16O27", "CoFe2O4"), row("SrZn2Fe16O27", "ZnFe2O4"),
            row("SrMg2Fe16O27", "SrMg2Fe16O27"), row("SrMg2Fe16O27", "MgFe2O4"), row("SrNi2Fe16O27", "NiFe2O4")]
    text = LS._csv_write([HDR] + body)
    fit, hold, n = LS.split_structure_csv(text)
    f, h = LS._csv_read(fit), LS._csv_read(hold)
    assert n == 3 and len(h) == 4 and len(f) == 4
    assert {r[1] for r in h[1:]} == {"SrCo2Fe16O27", "CoFe2O4", "ZnFe2O4"}
    assert {r[1] for r in f[1:]} == {"SrMg2Fe16O27", "MgFe2O4", "NiFe2O4"}
    assert not any("Co" in r[2] or "Zn" in r[2] for r in f[1:])


def test_l02_rows_of_the_same_formula_stay_in_fit():
    l02 = ["1", "SrCo2Fe16O27", "L02_Morch2023/supporting/CIFs/pdCIF_SrCo2Fe16O27.cif", "Fe1", "1"]
    fit, hold, n = LS.split_structure_csv(LS._csv_write([HDR, l02]))
    assert n == 0 and len(LS._csv_read(fit)) == 2


def test_other_column_layouts_route_by_any_cell():
    hdr = ["source", "sample", "cif", "site"]
    rows = [["L01", "SrZn2Fe16O27", CIF.format(s="SrZn2Fe16O27"), "4f"], ["L01", "SrMg2Fe16O27", "x", "4f"],
            ["L01", "", CIF.format(s="SrCo2Fe16O27"), "12k"]]
    fit, hold, n = LS.split_structure_csv(LS._csv_write([hdr] + rows))
    assert n == 2 and len(LS._csv_read(fit)) == 2


# ---------------------------------------------------------------- coverage map: fit only, no holdout mention
COV_HDR = ["composition_class", "composition", "x_total_A", "property", "conditions", "status", "sources",
           "transferability_tier", "proposed_split", "gap_action", "comment"]


def cov(composition, prop, sources, comment="", split="L03 fit"):
    return ["c", composition, "2", prop, "RT", "covered", sources, "1", split, "", comment]


def test_coverage_holdout_mentions_are_rewritten_without_ids_or_numbers():
    zn = cov("SrZn2Fe16O27 (x=2)", "synthesis window / phase purity / Fe2+ content",
             "L03 (spinel); L01 (ZnFe2O4 9.99 wt %); L02 (spinel 4.12 wt %)")
    co = cov("SrCo2Fe16O27 (x=2)", "synthesis window / phase purity / Fe2+ content", "L03 (x); L01 (1300 C 2 h); C023, C015")
    ba = cov("SrCo2Fe16O27 (x=2)", "K1,K2 / HA (experimental)", "none for Sr (Ba: C057 K1=8.7e4; L11 Ha formal)", split="")
    al = cov("all W compositions", "g-factor, alpha, FMR linewidth", "x", comment="L12 gives only model; C056/C036", split="")
    mg = cov("SrMg2Fe16O27", "lattice a,c", "L01 Table S4", split="L01 fit")
    text = LS._csv_write([COV_HDR, zn, co, ba, al, mg])
    out, kept, total = LS.filter_coverage(text)
    assert kept == total == 5
    rows = LS._csv_read(out)[1:]
    blob = " ".join(" ".join(r) for r in rows[:4])
    for bad in ("L01", "C057", "C023", "C015", "C056", "C036", "L12", "9.99", "8.7e4", "1300 C"):
        assert bad not in blob, bad
    assert "holdout" in blob and "L02 (spinel 4.12, M 1.77, X 2.67 wt %)" in blob
    assert rows[4][6] == "L01 Table S4"                        # Mg/Ni rows keep their L01 citation


def test_coverage_unknown_holdout_mention_aborts():
    bad = cov("SrZn2Fe16O27 (x=2)", "some new property", "L01 (something)")
    with pytest.raises(LS.LeakError):
        LS.filter_coverage(LS._csv_write([COV_HDR, bad]))
    with pytest.raises(LS.LeakError):
        LS.filter_coverage(LS._csv_write([COV_HDR, cov("SrNi2Fe16O27", "p", "C056 data")]))


# ---------------------------------------------------------------- L01 fit records: no hint at the withheld samples
def mk_claim(rec_id, sample, verbatim=None, synthesis=None, lit="L01", role="literature_fit"):
    src = S.LiteratureSource(locator=S.Locator(text="Table 1"), data_status="measured", tier=2,
                             sample_id_in_source=sample, formula_in_source=sample, verbatim_source=verbatim,
                             synthesis=synthesis)
    return S.LiteratureClaim(rec_id, lit, "x", S.SourceRole(role), src, qualitative="q")


def test_scrub_l01_fit_removes_other_samples_hints():
    syn = "gel fired; then fired 1200 C (Mg, Zn) or 1300 C (Ni, Co), 2 h hold; this sample fired at 1200 C"
    c = mk_claim("L01-088", "SrMg2Fe16O27", "For Mg, Ni and Zn there is no magnetic contribution", syn)
    assert LS.scrub_l01_fit(c)
    assert "Zn" not in c.source.verbatim_source and "withheld" not in c.source.verbatim_source
    assert "Zn" not in c.source.synthesis and "Co)" not in c.source.synthesis and "this sample fired at 1200 C" in c.source.synthesis


# ---------------------------------------------------------------- final check: fit refers to nothing of the holdout
@pytest.fixture
def sroot(tmp_path):
    for rel, txt in {"raw/theory_side/06_literature/L01_Morch2019/supporting/SrCo2Fe16O27/SrCo2Fe16O27_SINQ_HRPT.cif": "x",
                     "raw/theory_side/06_literature/L01_Morch2019/supporting/SrMg2Fe16O27/SrMg2Fe16O27_SINQ_HRPT.cif": "x",
                     "raw/theory_side/06_literature/L01_Morch2019/paper_fulltext.xml": "x",
                     "lit_search/downloads/C056/C056_qin_2013.pdf": "x",
                     "dft_feasibility/pilots/c_W92/Model_SrZn2.cif": "x"}.items():
        p = tmp_path / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(txt)
    srcs = [tmp_path / "raw/theory_side/06_literature/L01_Morch2019", tmp_path / "lit_search/downloads/C056",
            tmp_path / "dft_feasibility/pilots/c_W92", tmp_path / "dft_feasibility/pilots/make_bw.py"]
    return tmp_path, srcs


def find(sroot, texts, claims=()):
    root, srcs = sroot
    return LS.find_fit_references(root, srcs, list(claims), texts)


def test_reference_check_clean_fit_passes(sroot):
    texts = {"fit/a.csv": "L01_Morch2019/supporting/SrMg2Fe16O27/SrMg2Fe16O27_SINQ_HRPT.cif,SrMg2Fe16O27,5.89\n"
                          "L02,SrCo2Fe16O27,pdCIF_SrCo2Fe16O27.cif\n",
             "fit/README.md": "Source: raw/theory_side/06_literature/L01_Morch2019/supporting (Mg, Ni)\n"}
    assert find(sroot, texts, [mk_claim("L01-001", "SrMg2Fe16O27")]) == []


@pytest.mark.parametrize("text", [
    "x,L01_Morch2019/supporting/SrCo2Fe16O27/SrCo2Fe16O27_SINQ_HRPT.cif,y",         # path of a holdout-sample CIF
    "see SrZn2Fe16O27_SINQ_HRPT.cif",                                               # file name only (absent file: no term)
    "from lit_search\\downloads\\C056 folder",                                      # backslash spelling, whole holdout pub
    "see C056_qin_2013.pdf",                                                        # file name of a holdout source
    "built by dft_feasibility/pilots/make_bw.py",                                   # declared derivative (file)
    "| L01 | SrCo2Fe16O27 | 5.9 |",                                                  # publication next to a holdout sample
    "L01 SrZn2Fe16O27 spinel",
])
def test_reference_check_finds_leaks(sroot, text):
    root, srcs = sroot
    if text.startswith("see SrZn2"):
        # the CIF does not exist in the synthetic tree -> create it so that its name is a term
        p = root / "raw/theory_side/06_literature/L01_Morch2019/supporting/SrZn2Fe16O27/SrZn2Fe16O27_SINQ_HRPT.cif"
        p.parent.mkdir(parents=True)
        p.write_text("x")
    hits = find(sroot, {"fit/x.txt": text})
    assert hits, text


def test_reference_check_claims(sroot):
    assert find(sroot, {}, [mk_claim("L01-050", "SrCo2Fe16O27")])                       # fit claim of a holdout sample
    assert find(sroot, {}, [mk_claim("L01-051", "SrMg2Fe16O27", verbatim="the Zn sample SrZn2Fe16O27 differs")])
    assert find(sroot, {}, [mk_claim("L01-052", "x")]) == []
    # the same formula in a fit publication of another lit_id (L02 Co/Zn) is fine
    assert find(sroot, {}, [mk_claim("L02-001", "SrCo2Fe16O27", lit="L02")]) == []


def test_reference_check_messages_carry_no_values(sroot):
    hits = find(sroot, {"fit/x.csv": "L01,SrCo2Fe16O27,12.3456"})
    assert hits and not any("12.3456" in h for h in hits)


# ---------------------------------------------------------------- built literature/ (content)
def _need_real():
    if not (REAL / "literature" / "fit" / "records_fit.yaml").is_file():
        pytest.skip("literature split not built")


def test_built_fit_has_no_l01_co_zn_rows_and_no_holdout_ids_in_coverage():
    _need_real()
    fit = REAL / "literature" / "fit"
    raw = list(csv.reader(io.StringIO((fit / "structures" / "_raw_asym_units.csv").read_text(encoding="utf-8"), newline="")))
    h = raw[0]
    fcol = h.index("file")
    for r in raw[1:]:
        assert not ("L01" in r[fcol] and ("SrCo2" in r[fcol] or "SrZn2" in r[fcol])), r[fcol]
    for p in sorted((fit / "structures").glob("*.csv")):
        for r in csv.reader(io.StringIO(p.read_text(encoding="utf-8"), newline="")):
            line = ",".join(r)
            assert not ("L01" in line and ("SrCo2Fe16O27" in line or "SrZn2Fe16O27" in line)), p.name
    cov = (fit / "coverage_map.csv").read_text(encoding="utf-8")
    assert not LS._HOLDOUT_ID_RE.search(cov)
    for r in csv.DictReader(io.StringIO(cov, newline="")):
        if not r["composition"].startswith(("SrMg", "SrNi")):
            assert "L01" not in r["sources"] + r["comment"] + r["conditions"], r["composition"]


def test_built_fit_passes_the_final_reference_check():
    _need_real()
    cfg = paths.get_config()
    from hexw.params.literature import load_fit_records
    texts = {}
    for sub in ("fit", "excluded"):
        for q in sorted((cfg.literature / sub).rglob("*")):
            if q.is_file():
                texts[f"{sub}/{q.name}"] = q.read_text(encoding="utf-8", errors="ignore")
    assert LS.find_fit_references(cfg.root, cfg.holdout_sources, load_fit_records(), texts) == []


def test_l01_fit_claims_have_no_hint_at_holdout_samples():
    _need_real()
    doc = yaml.safe_load((REAL / "literature/fit/records_fit.yaml").read_text(encoding="utf-8"))
    for c in doc["claims"]:
        if c["lit_id"] != "L01":
            continue
        blob = json.dumps(c, ensure_ascii=False)
        assert "a further sample" not in blob and "SrCo2Fe16O27" not in blob and "SrZn2Fe16O27" not in blob, c["rec_id"]
        assert "(Mg, Zn)" not in blob and "(Ni, Co)" not in blob, c["rec_id"]
    by = {c["rec_id"]: c for c in doc["claims"]}
    for rid in ("L01-088", "L01-089"):
        assert "withheld" not in by[rid]["source"]["verbatim_source"]


def test_pilot_derivatives_of_the_l01_zn_cif_are_declared_holdout():
    _need_real()
    cfg = paths.get_config()
    pil = cfg.root / "dft_feasibility" / "pilots"
    if not pil.is_dir():
        pytest.skip("no dft_feasibility")
    for rel in ("c_W92/SrZn2Fe16O27_model.cif", "c_W92/W_SrZn2Fe16O27_k111.in", "c_W92/W_SrZn2Fe16O27_k221.in",
                "c_W92/W_k111_3it.in", "c_W92/W_k111_3it.out", "c_W92/SrFe18O27_model.cif", "make_bw.py"):
        if (pil / rel).exists():
            assert paths.classify(pil / rel) is paths.Zone.HOLDOUT, rel
    from _helpers import run_as
    s = run_as("hexw.dft._probe", "def go(p):\n    return open(p, 'rb').read(4)\n")
    for rel in ("c_W92/SrZn2Fe16O27_model.cif", "make_bw.py"):
        if (pil / rel).exists():
            with pytest.raises(paths.TheoryAccessError):
                s["go"](pil / rel)


def test_split_decision_has_incident_and_dft_rule():
    _need_real()
    t = (REAL / "literature/SPLIT_DECISION.md").read_text(encoding="utf-8")
    assert "Инцидент 2026-10-26" in t and "L02" in t and "pdCIF_SrZn2Fe16O27" in t
    assert "до этого момента никакой параметризации по fit не проводилось" in t


# ---------------------------------------------------------------- schema details
def pe(role):
    return S.ParameterEstimate("x", 1.0, "muB", "magnetic_moment", role, S.Provenance("L"))


def test_holdout_estimate_cannot_enter_a_model_version():
    h = pe("literature_holdout")                                # buildable: it is the data of a holdout claim
    assert h.source_role is S.SourceRole.LITERATURE_HOLDOUT
    with pytest.raises(V.LeakageError, match="holdout"):
        V.make_model_version([pe("literature_fit"), S.ParameterEstimate(
            "y", 2.0, "muB", "magnetic_moment", "literature_holdout", S.Provenance("L"))])
    V.make_model_version([pe("literature_fit")])                # fit is fine
    for role in ("excluded_secondary", "excluded_low_quality", "methodology_only"):
        with pytest.raises(ValueError):
            pe(role)
    pred = S.Prediction("p1", "mv-1", [S.DataRef("x", source_role="literature_holdout")])
    with pytest.raises(V.LeakageError):
        V.check_prediction_inputs_roles(pred)


def test_site_moment_eq_across_definitions_raises():
    fe = S.SiteMoment(3.0, "fe_ion")
    co = S.SiteMoment(3.0, "species_ion", species="Co")
    avg = S.SiteMoment(3.0, "site_average_incl_diamagnetic")
    for a, b in ((fe, co), (fe, avg), (co, avg), (co, S.SiteMoment(3.0, "species_ion", species="Zn"))):
        with pytest.raises(S.MomentDefinitionError):
            a == b
        with pytest.raises(S.MomentDefinitionError):
            a != b
    assert fe == S.SiteMoment(3.0, "fe_ion") and fe != S.SiteMoment(2.0, "fe_ion")
    assert fe == S.SiteMoment(3000.0, "fe_ion", unit="mmuB") if "mmuB" in S.UNITS[S.Kind.MAGNETIC_MOMENT] else True
    assert (fe == 3.0) is False


# ---------------------------------------------------------------- audit: future derivatives are caught
@pytest.fixture
def aroot(tmp_path):
    for rel in ("raw/theory_side/06_literature/L01_Morch2019", "raw/original_zip", "vault/lab_validation", "registry",
                "lit_search/downloads/C056", "literature/fit", "literature_holdout", "project/config", "pilots", "docs"):
        (tmp_path / rel).mkdir(parents=True)
    (tmp_path / "raw/theory_side/06_literature/L01_Morch2019/SrCo2Fe16O27_HRPT.cif").write_text("x")
    (tmp_path / "lit_search/downloads/C056/C056_qin.pdf").write_text("x")
    (tmp_path / "pilots/derived_declared.in").write_text("geometry from lit_search/downloads/C056/C056_qin.pdf")
    (tmp_path / "project/config/holdout_sources.json").write_text(json.dumps({"holdout_sources": [
        "raw/theory_side/06_literature/L01_Morch2019", "lit_search/downloads/C056", "pilots/derived_declared.in"]}))
    with paths.test_root(tmp_path):
        yield tmp_path


CFG_JSON = "project/config/holdout_sources.json"


def test_audit_finds_undeclared_derivatives_and_accepts_declared_or_allowlisted(aroot):
    from hexw import audit
    cfg = paths.get_config()
    (aroot / "pilots" / "new_pilot.in").write_text("cell from SrCo2Fe16O27_HRPT.cif of the holdout")
    (aroot / "pilots" / "other.py").write_text("p = r'raw" + chr(92) + "theory_side" + chr(92) + "06_literature" + chr(92) + "L01_Morch2019" + chr(92) + "x.cif'")
    (aroot / "pilots" / "clean.py").write_text("print('nothing to see: L02 pdCIF_SrZn2Fe16O27.cif')")
    (aroot / "docs" / "notes.md").write_text("see lit_search/downloads/C056/C056_qin.pdf")
    with paths.test_access():
        off, n, nt = audit.find_undeclared_holdout_references(cfg, allow=[CFG_JSON])
        assert sorted(off) == ["docs/notes.md", "pilots/new_pilot.in", "pilots/other.py"] and nt >= 3
        # the declared derivative is not an offender (classified HOLDOUT); allow-listed files are accepted
        off2, _n, _t = audit.find_undeclared_holdout_references(cfg, allow=[CFG_JSON, "docs/*", "pilots/new_pilot.in"])
        assert off2 == ["pilots/other.py"]
        assert paths.classify(aroot / "pilots" / "derived_declared.in") is paths.Zone.HOLDOUT


def test_real_project_has_no_undeclared_holdout_references():
    _need_real()
    from hexw import audit
    cfg = paths.get_config()
    with paths.test_access():
        off, n, nt = audit.find_undeclared_holdout_references(cfg)
    assert off == [] and n > 100 and nt >= 10


def test_audit_pinned_allowlist_accepts_only_the_pinned_bytes(aroot):
    """Stage 4 round 1: a script / document that names a holdout source is allowed only while its sha256 equals the pin."""
    import hashlib
    from hexw import audit
    cfg = paths.get_config()
    f = aroot / "docs" / "pipeline.py"
    f.write_text("INPUT = 'lit_search/downloads/C056/C056_qin.pdf'\n")
    pin = hashlib.sha256(f.read_bytes()).hexdigest()
    with paths.test_access():
        off, _n, _t = audit.find_undeclared_holdout_references(cfg, allow=[CFG_JSON], pinned={"docs/pipeline.py": pin})
        assert off == []
        f.write_text("INPUT = 'lit_search/downloads/C056/C056_qin.pdf'  # edited\n")
        off, _n, _t = audit.find_undeclared_holdout_references(cfg, allow=[CFG_JSON], pinned={"docs/pipeline.py": pin})
        assert off == ["docs/pipeline.py (sha256 differs from the pin)"]
        off, _n, _t = audit.find_undeclared_holdout_references(cfg, allow=[CFG_JSON], pinned={})
        assert off == ["docs/pipeline.py"]
