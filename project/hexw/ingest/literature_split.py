"""Reproducible fit / holdout split of the extracted literature base (decision: literature/SPLIT_DECISION.md).

Input : literature/{records.yaml, records.jsonl, coverage_map.csv, CONFLICTS.md, READING_LOG.md, split_proposal.md,
        site_notation.md, structures/, scripts/}  (first run)   or   literature_holdout/_original_mixed/  (re-run).
Output: literature/fit/            fit claims (LiteratureClaim), fit structure tables, notation, CONFLICTS (no holdout numbers)
        literature/excluded/       excluded citations / methodology of fit publications (never parametrization)
        literature_holdout/        holdout claims + structure tables + excluded records of holdout publications,
                                   SOURCE_FILES_SHA256.txt, SHA256SUMS, _original_mixed/ (originals, read-only)
        literature/SPLIT_DECISION.md, literature/split_manifest.json, literature/SHA256SUMS

Rules implemented (see DECISION below):
  * every parameter record becomes a validated `LiteratureClaim`; kinds that the schema lacked in stage 2 are filled
    now (pending -> ParameterEstimate; fraction -> fraction_site / fraction_mass / fraction_volume);
  * fit records of L01 (Mg, Ni) lose the Co/Zn columns of table rows from `verbatim_source`;
  * a leak scan compares numeric markers of the holdout side with the text of every fit output: any hit aborts the
    run before anything is written;
  * runs inside `paths.holdout_build()`: the only code allowed to write the HOLDOUT zone (hexw.ingest only);
  * refuses to overwrite finished outputs without --force.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import re
import shutil
import stat
import sys
from collections import Counter, OrderedDict
from pathlib import Path

import yaml

from hexw.io import paths
from hexw.registry import schemas as S

DATE = "2026-10-25"
AUTHOR = "руководитель разработки"
AUTHORED = Path(__file__).resolve().parent / "literature_authored"

# ---------------------------------------------------------------- the decision ("variant A with an amendment")
FIT_SOURCES = {"L02", "L05", "L03", "C004", "L11"}          # whole publications (secondary citations are excluded)
L01_FIT_SAMPLES = {"SrMg2Fe16O27", "SrNi2Fe16O27"}
L01_HOLDOUT_SAMPLES = {"SrCo2Fe16O27", "SrZn2Fe16O27"}
HOLDOUT_KIND = {
    "L01": "same_group_reproducibility",                    # samples Co and Zn only
    "C023": "independent_group", "C015": "independent_group",
    "L12": "transfer_Ba_to_Sr", "C057": "transfer_Ba_to_Sr", "C060": "transfer_Ba_to_Sr", "C044": "transfer_Ba_to_Sr",
    "C056": "transfer_Ba_to_Sr", "C036": "transfer_Ba_to_Sr", "C052": "transfer_Ba_to_Sr", "C066": "transfer_Ba_to_Sr",
}
HOLDOUT_PUBLICATIONS = set(HOLDOUT_KIND)                    # L01 is blocked as a whole publication
METHODOLOGY = {"C112", "C115"}
# DOIs that candidates.csv must give for the lit_search/downloads/C### folders that are declared as holdout sources
EXPECTED_DOI = {
    "C015": "10.1063/5.0087172", "C023": "10.48550/arxiv.1407.3605", "C036": "10.3390/ma16175897",
    "C044": "10.47011/13.1.1", "C052": "10.3379/msjmag.1506r003", "C056": "10.4236/msce.2013.14002",
    "C057": "10.1063/1.3670062", "C060": "10.1088/1742-6596/303/1/012045", "C066": "10.2497/jjspm.43.995",
}
FIT_TIER_OVERRIDE = {"L11": 3}
FIT_TRANSFER = {"L11": "Ba_to_Sr"}

INDEPENDENCE = (
    "L11 (He S., Liu R., Liu X., Ye X., Wang L., Shen B.; Chin. Phys. B 33 (2024) 066801; BaMe2W, solid-state synthesis) и "
    "C056 (Qin X., Cheng Y., Zhou K., Huang S., Hui X.; J. Mater. Sci. Chem. Eng. 1 (2013) 8; sol-gel composites "
    "Ba(MnZn)xCo2(1-x)W with paraffin, Central South University) - **разные группы авторов** (по спискам авторов и "
    "аффилиациям на первых страницах статей нет пересечений), разные годы, методы синтеза и серии образцов. "
    "Основание для объединения в одну выборку отсутствует, поэтому **L11 в fit, C056 в holdout** (как решено). "
    "Оговорка: C056 содержит Mn только в смешанной подрешётке A (Mn:Zn = 1:1 вместе с Co) и измеряет "
    "микроволновое поглощение композитов, поэтому проверяет Mn-ветку косвенно, не как чистый BaMn2W.")

SECONDARY_DOC = ("excluded_secondary_citation", "excluded_low_quality", "methodology_only")


class LiteratureSplitError(RuntimeError):
    pass


class OutputsExistError(LiteratureSplitError, PermissionError):
    pass


class LeakError(LiteratureSplitError):
    pass


# ---------------------------------------------------------------- helpers
def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha256_path(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _dump_yaml(obj) -> str:
    return yaml.safe_dump(obj, allow_unicode=True, sort_keys=False, width=120)


def _writable(p: Path):
    if p.exists():
        os.chmod(p, os.stat(p).st_mode | stat.S_IWRITE)


def _readonly(p: Path):
    os.chmod(p, os.stat(p).st_mode & ~(stat.S_IWRITE | stat.S_IWGRP | stat.S_IWOTH))


# ---------------------------------------------------------------- classification
def classify(rec: dict) -> dict:
    """-> {role, tier, transfer, holdout_kind}. Raises for anything the decision does not cover."""
    lit = rec["lit_id"]
    ext = rec["ext"]
    old = ext["split"]
    if old == "excluded_secondary_citation":
        return {"role": S.SourceRole.EXCLUDED_SECONDARY, "tier": ext["tier"], "transfer": None, "holdout_kind": None}
    if old == "excluded_low_quality":
        return {"role": S.SourceRole.EXCLUDED_LOW_QUALITY, "tier": ext["tier"], "transfer": None, "holdout_kind": None}
    if old == "methodology_only" or lit in METHODOLOGY:
        return {"role": S.SourceRole.METHODOLOGY_ONLY, "tier": ext["tier"], "transfer": None, "holdout_kind": None}
    if lit == "L01":
        sample = ext.get("sample_id_in_paper")
        if sample in L01_FIT_SAMPLES:
            return {"role": S.SourceRole.LITERATURE_FIT, "tier": ext["tier"], "transfer": None, "holdout_kind": None}
        if sample in L01_HOLDOUT_SAMPLES:
            return {"role": S.SourceRole.LITERATURE_HOLDOUT, "tier": ext["tier"], "transfer": None,
                    "holdout_kind": HOLDOUT_KIND["L01"]}
        raise LiteratureSplitError(f"{rec['rec_id']}: L01 record with unknown sample {sample!r}")
    if lit in FIT_SOURCES:
        return {"role": S.SourceRole.LITERATURE_FIT, "tier": FIT_TIER_OVERRIDE.get(lit, ext["tier"]),
                "transfer": FIT_TRANSFER.get(lit), "holdout_kind": None}
    if lit in HOLDOUT_KIND:
        return {"role": S.SourceRole.LITERATURE_HOLDOUT, "tier": ext["tier"], "transfer": None,
                "holdout_kind": HOLDOUT_KIND[lit]}
    raise LiteratureSplitError(f"{rec['rec_id']}: source {lit!r} is not covered by the split decision")


# ---------------------------------------------------------------- record -> claim
_PENDING_KINDS = {
    "magnetic_moment": "magnetic_moment", "specific_magnetization": "specific_magnetization",
    "energy_density": "energy_density", "pressure": "pressure", "volume": "volume",
    "magnetic_polarization_mu0M": "magnetic_polarization_mu0M", "energy_per_formula_unit": "energy_per_formula_unit",
    "energy_per_atom": "energy_per_atom", "reflection_loss": "reflection_loss",
}
_PAREN_RE = re.compile(r"^\s*[+-]?\d*\.?\d+\((\d+(?:\.\d+)?)\)")
_TABLE_RE = re.compile(r"\bTables?\s+(S?\d+[A-Za-z]?)", re.I)
_FIG_RE = re.compile(r"\bFig(?:ure)?s?\.?\s+(S?\d+[a-z]?)", re.I)
_PAGE_RE = re.compile(r"\bp(?:p|age)?\.?\s*(\d+)", re.I)


def _refine_kind(name: str, kind: str, unit: str) -> str:
    """Stage 3 refinement of the legacy `fraction` kind: explicit basis."""
    if kind != "fraction":
        return kind
    n = name.lower()
    if n.startswith("site_occupancy"):
        return "fraction_site"
    if "weight_fraction" in n or "impurity_wt" in n or "wt_" in n:
        return "fraction_mass"
    if "volume_fraction" in n:
        return "fraction_volume"
    return kind


def _convention(value_original, has_std: bool) -> str:
    if has_std and isinstance(value_original, str) and _PAREN_RE.match(value_original):
        return "parenthetical_se"
    return "none"


def _locator(text) -> dict:
    d = {"text": text, "page": None, "table": None, "figure": None, "supplementary": None}
    if not text:
        return d
    m = _TABLE_RE.search(text)
    if m:
        d["table"] = m.group(1)
    m = _FIG_RE.search(text)
    if m:
        d["figure"] = m.group(1)
    m = _PAGE_RE.search(text)
    if m:
        d["page"] = m.group(1)
    if re.search(r"supporting|supplement|\bSI\b|\bS\d+\b|sup1", text, re.I):
        d["supplementary"] = text
    return d


def to_claim(rec: dict, dec: dict) -> S.LiteratureClaim:
    ext = rec["ext"]
    role = dec["role"]
    est = rec.get("estimate")
    pend = rec.get("pending")
    migration = []
    estimate = None
    pending = None
    if role not in S.NON_PARAMETRIC_ROLES:
        if est:
            e = json.loads(json.dumps(est))
            k2 = _refine_kind(e["name"], e["quantity_kind"], e["unit"])
            if k2 != e["quantity_kind"]:
                migration.append(f"quantity_kind {e['quantity_kind']} -> {k2}")
                e["quantity_kind"] = k2
            e["source_role"] = role.value
            u = e.setdefault("uncertainty", {"kind": "unknown"})
            conv = _convention(ext.get("value_original"), u.get("kind") == "std")
            if conv != "none":
                u["convention"] = conv
            estimate = S.from_dict(S.ParameterEstimate, e)
        elif pend:
            pk = pend.get("proposed_quantity_kind")
            if pk in _PENDING_KINDS:
                std = pend.get("uncertainty_std")
                u = {"kind": "std", "std": std} if std is not None else {"kind": "unknown"}
                conv = _convention(ext.get("value_original"), std is not None)
                if conv != "none":
                    u["convention"] = conv
                e = {"name": pend["name"], "value": pend["value"], "unit": pend["unit"],
                     "quantity_kind": _PENDING_KINDS[pk], "source_role": role.value,
                     "provenance": {"source_id": rec["lit_id"], "locator": ext.get("locator"),
                                    "note": None, "status": "unverified"},
                     "uncertainty": u, "identifiable": True}
                e["quantity_kind"] = _refine_kind(e["name"], e["quantity_kind"], e["unit"])
                estimate = S.from_dict(S.ParameterEstimate, e)
                migration.append(f"pending ({pk}, {pend['unit']}) -> ParameterEstimate")
            else:
                pending = dict(pend)
    if role in S.NON_PARAMETRIC_ROLES and (est or pend):
        # excluded / methodology records carry no value that could be parametrized; keep it as `pending`
        # (documentation only) when a number was extracted
        if pend:
            pending = dict(pend)
        elif est:
            pending = {"name": est["name"], "value": est["value"], "unit": est["unit"],
                       "proposed_quantity_kind": est["quantity_kind"], "reason": "excluded record (no parametrization)",
                       "uncertainty_std": (est.get("uncertainty") or {}).get("std")}
    notes = " | ".join(x for x in (ext.get("notes"), ext.get("transfer_note")) if x) or None
    src = S.LiteratureSource(
        locator=S.Locator(**_locator(ext.get("locator") or f"{rec['lit_id']} (locator not recorded)")),
        data_status=ext["status"], tier=dec["tier"], doi=rec.get("doi") or None,
        sample_id_in_source=ext.get("sample_id_in_paper"), formula_in_source=ext.get("formula"),
        synthesis=ext.get("synthesis"), measurement_conditions=S.MeasurementConditions(**(ext.get("measurement") or {})),
        method=ext.get("method"), quantity_definition=ext.get("definition"), verbatim_source=ext.get("verbatim"),
        transformation=ext.get("transformation"), measurement_group=ext.get("measurement_group"),
        figure_only=bool(ext.get("figure_only")), axes=ext.get("axes"), split_cluster=ext.get("split_cluster"),
        transfer=dec["transfer"], holdout_kind=dec["holdout_kind"])
    if migration:
        notes = (notes + " | " if notes else "") + "schema migration: " + "; ".join(migration)
    return S.LiteratureClaim(
        rec_id=rec["rec_id"], lit_id=rec["lit_id"], name=rec["name"], role=role, source=src, estimate=estimate,
        pending=pending, qualitative=ext.get("qualitative"), value_original=ext.get("value_original"),
        unit_original=ext.get("unit_original"), uncertainty_note=ext.get("uncertainty_note"), notes=notes)


def claim_plain(c: S.LiteratureClaim) -> dict:
    """Plain dict for storage: the source trace is stored once (the loader re-attaches it to the estimate)."""
    d = S.to_plain(c)
    if d.get("estimate") and d["estimate"].get("provenance"):
        d["estimate"]["provenance"].pop("literature", None)
    return d


# ---------------------------------------------------------------- L01 fit: remove holdout columns
def scrub_l01_fit(c: S.LiteratureClaim) -> bool:
    """Remove from a fit record of L01 (Mg, Ni) everything that carries the Co/Zn samples' data."""
    s = c.source
    changed = False
    v = s.verbatim_source or ""
    own = c.value_original
    if re.search(r"\(columns [^)]*\)", v):
        label = re.split(r"\s{2,}", v.strip())[0]
        s.verbatim_source = f"{label}  {own}  [column {s.sample_id_in_source}; the columns of the other samples of this row are withheld: holdout]"
        changed = True
    elif re.match(r"^.*?\S\s{2,}(?:-?\d[\d.]*(?:\(\d+\))?\s{2,}){3}-?\d[\d.]*(?:\(\d+\))?\s*$", v):
        # a 4-column row in the sample order Mg, Co, Ni, Zn without a column legend
        label = re.split(r"\s{2,}", v.strip())[0]
        s.verbatim_source = (f"{label}  {own}  [column {s.sample_id_in_source}; the columns of the other samples of "
                             "this row are withheld: holdout]")
        changed = True
    elif "small impurities of" in v and "->" in v:
        s.verbatim_source = v.split("->", 1)[1].strip() + "  [impurity phases of the other samples are withheld: holdout]"
        changed = True
    elif re.search(r"\bMg, Ni and Zn\b", v):
        s.verbatim_source = v.replace("Mg, Ni and Zn", "Mg and Ni")          # no hint at the withheld sample
        changed = True
    if s.synthesis:
        # "fired 1200 C (Mg, Zn) or 1300 C (Ni, Co), 2 h hold ... this sample fired at T": keep this sample only
        syn = re.sub(r"fired 1200 C \(Mg, Zn\) or 1300 C \(Ni, Co\), 2 h hold", "fired, 2 h hold", s.synthesis)
        if syn != s.synthesis:
            s.synthesis, changed = syn, True
    if s.quantity_definition:
        d2 = re.sub(r"\s*\(Co: planar[^)]*\)", "", s.quantity_definition)
        if d2 != s.quantity_definition:
            s.quantity_definition, changed = d2, True
    if c.uncertainty_note:
        n2 = re.sub(r";\s*for Zn 4e '\(1\.3\)' is 1\.3 muB", "", c.uncertainty_note)
        if n2 != c.uncertainty_note:
            c.uncertainty_note, changed = n2, True
    if changed:
        c.notes = ((c.notes + " | ") if c.notes else "") + "holdout columns removed from the quoted source text (SPLIT_DECISION.md)"
    return changed


# ---------------------------------------------------------------- structure tables
def _csv_read(text: str):
    return list(csv.reader(io.StringIO(text, newline="")))


def _csv_write(rows) -> str:
    buf = io.StringIO(newline="")
    csv.writer(buf, lineterminator="\n").writerows(rows)
    return buf.getvalue()


def _is_l01_holdout_row(header, row) -> bool:
    """A structure row belongs to the holdout if it comes from a CIF of a holdout SAMPLE of L01. The routing is by the
    SOURCE FILE (CIF path) / sample, never by the phase block: the spinel blocks (CoFe2O4, ZnFe2O4) of a Co / Zn CIF
    are rows of that CIF (incident 2026-10-26). Any cell naming L01 and a holdout sample sends the whole row away."""
    if not any("L01" in c for c in row):
        return False
    return any(sm in c for c in row for sm in L01_HOLDOUT_SAMPLES)


def split_structure_csv(text: str):
    rows = _csv_read(text)
    header, body = rows[0], rows[1:]
    fit = [header] + [r for r in body if not _is_l01_holdout_row(header, r)]
    hold = [header] + [r for r in body if _is_l01_holdout_row(header, r)]
    return _csv_write(fit), _csv_write(hold), len(hold) - 1


def scrub_structure_readme(text: str) -> str:
    out, lines, i = [], text.replace("\r\n", "\n").split("\n"), 0
    hold_row = re.compile(r"^\|\s*L01\s*\|\s*(SrCo2Fe16O27|SrZn2Fe16O27)\s*\|")
    while i < len(lines):
        ln = lines[i]
        if hold_row.match(ln):
            i += 1
            continue
        if ln.startswith("| позиция |"):                     # moments table: keep the Mg and Ni columns only
            hdr = [c.strip() for c in ln.strip().strip("|").split("|")]
            keep = [j for j, c in enumerate(hdr) if c in ("позиция", "Mg", "Ni")]
            while i < len(lines) and lines[i].startswith("|"):
                cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                out.append("| " + " | ".join(cells[j] for j in keep) + " |")
                i += 1
            continue
        out.append(ln)
        i += 1
    t = "\n".join(out)
    t = t.replace("; Co - планарный порядок", "")
    t = re.sub(r"; моменты на тип позиции сверены:[^)]*\)", ")", t)
    t = t.replace(" Сравнение L01 vs L02: см. `../CONFLICTS.md` (A).", "")
    t = t.replace("L02 (моменты Fe3+, Tables S1-S3; Co и Zn-строки в SI;", "L02 (моменты Fe3+, Tables S1-S3;")
    head = ("> Редакция fit: образцы L01 Co и Zn (holdout, SPLIT_DECISION.md) удалены из таблиц; их CIF недоступны теории.\n\n")
    return head + t


_HOLDOUT_ID_RE = re.compile(r"(?<![A-Za-z0-9])(?:L12|C015|C023|C036|C044|C052|C056|C057|C060|C066)(?![A-Za-z0-9])")
_COVERAGE_TEXT_COLS = ("sources", "comment", "conditions", "property", "gap_action")
# The fit coverage map is built from fit records only: a cell that names a holdout publication (or L01 outside the
# Mg/Ni rows) is replaced by an authored neutral text keyed by (composition, property). L01 Co/Zn and the other
# holdout publications are mentioned only as "data exist in the holdout" (no numbers, no identifiers).
COVERAGE_REWRITE = {
    ("SrCo2Fe16O27 (x=2)", "synthesis window / phase purity / Fe2+ content", "sources"):
        "L03 (pO2 0.2-387 atm; n*(Co)=2.03 at 387 atm; flux crystal W-phase Co 1.79); L02 Table 1 (1473 K 2 h); "
        "further data exist in the holdout (not shown)",
    ("SrZn2Fe16O27 (x=2)", "synthesis window / phase purity / Fe2+ content", "sources"):
        "L03 (spinel impurity grows at 387 atm); L02 (spinel 4.12, M 1.77, X 2.67 wt %); "
        "further data exist in the holdout (not shown)",
    ("SrCo2Fe16O27 (x=2)", "K1,K2 / HA (experimental)", "sources"):
        "none for Sr (Ba: L11 Ha(BaCo2W) formal); further data exist in the holdout (not shown)",
    ("all W compositions", "*", "comment"): "further data exist in the holdout (not shown)",
}


def _coverage_cell_needs_rewrite(header, row, col) -> bool:
    d = dict(zip(header, row))
    t = d.get(col, "")
    if _HOLDOUT_ID_RE.search(t):
        return True
    mg_ni = d.get("composition", "").startswith(("SrMg", "SrNi"))
    return bool(re.search(r"(?<![A-Za-z0-9])L01(?![A-Za-z0-9])", t)) and not mg_ni


def scrub_coverage(header, body):
    """-> (rows, n_rewritten). Raises LeakError for a holdout mention without an authored rewrite."""
    out, n = [], 0
    for r in body:
        r = list(r)
        d = dict(zip(header, r))
        for col in _COVERAGE_TEXT_COLS:
            if col not in header or not _coverage_cell_needs_rewrite(header, r, col):
                continue
            new = (COVERAGE_REWRITE.get((d["composition"], d["property"], col))
                   or COVERAGE_REWRITE.get((d["composition"], "*", col)))
            if new is None:
                raise LeakError(f"coverage_map row {d.get('composition')!r}/{d.get('property')!r}: cell {col!r} names "
                                "a holdout publication and has no authored rewrite")
            r[header.index(col)] = new
            n += 1
        out.append(r)
    return out, n


def filter_coverage(text: str):
    rows = _csv_read(text)
    header = rows[0]
    k = header.index("proposed_split")
    fit_tokens = {"", "literature_fit", "L01 fit", "L02 fit", "L03 fit"}
    body, _n = scrub_coverage(header, [r for r in rows[1:] if r[k].strip() in fit_tokens])
    keep = [header] + body
    return _csv_write(keep), len(keep) - 1, len(rows) - 1


# ---------------------------------------------------------------- leak scan
_NUM_TOKEN = re.compile(r"\d+\.\d{3,}|\d+\.\d+\(\d+\)|\d{3,}\.\d+")


def holdout_tokens(holdout_claims, hold_struct_texts) -> set:
    toks = set()
    for c in holdout_claims:
        for t in (c.value_original, c.source.verbatim_source):
            if isinstance(t, str):
                toks |= set(_NUM_TOKEN.findall(t))
    for t in hold_struct_texts:
        toks |= set(_NUM_TOKEN.findall(t))
    return {t for t in toks if len(t) >= 5}


def scan_leaks(tokens: set, fit_texts: dict, allow=frozenset()) -> list:
    hits = []
    for tok in sorted(tokens - set(allow)):
        pat = re.compile(r"(?<![\d.])" + re.escape(tok) + r"(?!\d)")
        for name, text in fit_texts.items():
            if pat.search(text):
                hits.append((tok, name))
    return hits


# ---------------------------------------------------------------- final check: fit refers to nothing of the holdout
MIXED_SAMPLES = {"L01": L01_HOLDOUT_SAMPLES}      # publications that have fit AND holdout samples
_LIT_PREFIX = re.compile(r"^([LC]\d+)")
_MIN_PATH_TERM, _MIN_NAME_TERM = 6, 8


def _norm_ref(t: str) -> str:
    return t.replace("\\", "/").lower()


def holdout_reference_terms(root, sources, mixed=None, *, listing=None) -> dict:
    """{normalised path / file-name term: description} that no fit row or record may contain.

    * a holdout source that is a whole holdout publication (or a derivative of one): its relative path, its folder name
      and the relative path and name of every file inside;
    * a MIXED publication (L01: Mg/Ni fit, Co/Zn holdout): the folder itself may be cited (the fit rows quote the
      publication), but not the files of the holdout samples (the CIFs of the Co and Zn samples).

    `listing(path) -> files` enumerates a holdout source; by default the files are listed from this (hexw.ingest)
    frame, which the guard allows only inside `paths.holdout_build()`; hexw.audit passes its own lister."""
    mixed = MIXED_SAMPLES if mixed is None else mixed
    listing = _all_files if listing is None else listing
    root = Path(root)
    terms: dict = {}

    def add(t, why, n):
        t = _norm_ref(t)
        if len(t) >= n:
            terms.setdefault(t, why)

    for sp in sources:
        sp = Path(sp)
        try:
            rel = sp.relative_to(root).as_posix()
        except ValueError:
            rel = sp.as_posix()
        m = _LIT_PREFIX.match(sp.name)
        samples = mixed.get(m.group(1)) if m else None
        files = list(listing(sp)) if sp.exists() else []
        if samples is None:
            add(rel, "holdout source", _MIN_PATH_TERM)
            if not re.fullmatch(r"[LC]\d{2,3}", sp.name):
                add(sp.name, "holdout source name", _MIN_PATH_TERM)
            for q in files:
                add(q.relative_to(root).as_posix(), "file of a holdout source", _MIN_PATH_TERM)
                add(q.name, "file name of a holdout source", _MIN_NAME_TERM)
        else:
            for q in files:
                rq = q.relative_to(root).as_posix()
                if any(sm in rq for sm in samples):
                    add(rq, "file of a holdout sample", _MIN_PATH_TERM)
                    add(q.name, "file name of a holdout sample", _MIN_NAME_TERM)
    return terms


def find_fit_references(root, sources, fit_claims, fit_texts: dict, mixed=None, *, listing=None) -> list:
    """Final guard (incident 2026-10-26): -> list of findings (empty = clean). Nothing from the fit side may

    * contain the path or the file name of a holdout source file (see `holdout_reference_terms`);
    * (mixed publication, e.g. L01) put a holdout sample next to its publication in one line / row, or be a fit claim
      of a holdout sample (sample_id_in_source / formula_in_source), or mention the sample inside a fit claim of
      that publication.
    The findings contain locations and terms only, never values."""
    mixed = MIXED_SAMPLES if mixed is None else mixed
    terms = holdout_reference_terms(root, sources, mixed, listing=listing)
    hits = []
    for name, text in fit_texts.items():
        nt = _norm_ref(text)
        for t, why in terms.items():
            if t in nt:
                hits.append(f"{name}: {why} {t!r}")
    for lit, samples in mixed.items():
        lit_re = re.compile(r"(?<![A-Za-z0-9])" + re.escape(lit) + r"(?![0-9A-Za-z])")
        for name, text in fit_texts.items():
            for i, line in enumerate(text.splitlines(), 1):
                if lit_re.search(line) and any(sm in line for sm in samples):
                    hits.append(f"{name}:{i}: {lit} next to a holdout sample")
        for c in fit_claims:
            if c.lit_id != lit:
                continue
            src = c.source
            if src.sample_id_in_source in samples or src.formula_in_source in samples:
                hits.append(f"{c.rec_id}: fit claim of the holdout sample {src.sample_id_in_source!r}")
                continue
            blob = json.dumps(S.to_plain(c), ensure_ascii=False)
            if any(sm in blob for sm in samples):
                hits.append(f"{c.rec_id}: a holdout sample is named inside a fit claim of {lit}")
    return hits


# ---------------------------------------------------------------- main
def _lit_records(src: dict, lit_ids, role):
    out = []
    for lr in src.get("literature_records", []):
        if lr["lit_id"] in lit_ids:
            d = {k: lr.get(k) for k in ("lit_id", "citation", "doi", "year", "path", "access", "license")}
            d["role"] = role.value
            d["claims"] = []
            S.from_dict(S.LiteratureRecord, json.loads(json.dumps(d)))     # schema check
            out.append(d)
    return out


def _all_files(p: Path):
    if p.is_file():
        yield p
    elif p.is_dir():
        for q in sorted(p.rglob("*")):
            if q.is_file():
                yield q


def run(root=None, *, write: bool = True, force: bool = False) -> dict:
    cfg = paths.get_config() if root is None else paths.configure(root)
    with paths.holdout_build():
        return _run(cfg, write, force)


def _run(cfg: paths.Config, write: bool, force: bool) -> dict:
    lit, hold = cfg.literature, cfg.literature_holdout
    mixed = hold / "_original_mixed"
    first_run = (lit / "records.yaml").is_file()
    src_dir = lit if first_run else mixed
    if not (src_dir / "records.yaml").is_file():
        raise FileNotFoundError(f"no records.yaml in {lit} or {mixed}")
    done = (lit / "fit" / "records_fit.yaml").exists() or (hold / "records_holdout.yaml").exists()
    if write and done and not force:
        raise OutputsExistError("the split outputs already exist (literature/fit, literature_holdout); "
                                "re-run with --force to rebuild them from _original_mixed")

    raw = (src_dir / "records.yaml").read_bytes()
    doc = yaml.safe_load(raw)
    inputs = OrderedDict()
    for name in ("records.yaml", "records.jsonl", "coverage_map.csv", "CONFLICTS.md", "READING_LOG.md",
                 "split_proposal.md", "site_notation.md"):
        p = src_dir / name
        if p.is_file():
            inputs[name] = sha256_path(p)
    # (site_notation.md is common and lives in literature/fit afterwards; on a re-run it is read from fit/)
    notation_src = src_dir / "site_notation.md"
    if not notation_src.is_file() and (lit / "fit" / "site_notation.md").is_file():
        notation_src = lit / "fit" / "site_notation.md"
    struct_src = src_dir / "structures"
    if struct_src.is_dir():
        for q in sorted(struct_src.iterdir()):
            if q.is_file():
                inputs["structures/" + q.name] = sha256_path(q)

    # ---- claims
    claims: list[S.LiteratureClaim] = []
    for rec in doc["parameter_records"]:
        claims.append(to_claim(rec, classify(rec)))
    ids = [c.rec_id for c in claims]
    if len(set(ids)) != len(ids):
        raise LiteratureSplitError("duplicate rec_id in the input")
    for c in claims:
        if c.lit_id == "L01" and c.role is S.SourceRole.LITERATURE_FIT:
            scrub_l01_fit(c)
    fit = [c for c in claims if c.role is S.SourceRole.LITERATURE_FIT]
    hold_c = [c for c in claims if c.role is S.SourceRole.LITERATURE_HOLDOUT]
    excl = [c for c in claims if c.role in S.NON_PARAMETRIC_ROLES]
    excl_fit_side = [c for c in excl if c.lit_id not in HOLDOUT_PUBLICATIONS]
    excl_hold_side = [c for c in excl if c.lit_id in HOLDOUT_PUBLICATIONS]

    # ---- DOIs of the holdout downloads must match candidates.csv (C### <-> publication)
    cand_problems = _check_candidates(cfg)

    # ---- structure tables, notation, coverage
    files_fit: dict[str, str] = {}
    files_hold: dict[str, str] = {}
    n_hold_struct = 0
    if struct_src.is_dir():
        for q in sorted(struct_src.iterdir()):
            if not q.is_file():
                continue
            text = q.read_text(encoding="utf-8")
            if q.suffix == ".csv":
                f, h, n = split_structure_csv(text)
                files_fit["fit/structures/" + q.name] = f
                files_hold["structures_holdout/" + q.name] = h
                n_hold_struct += n
            elif q.name == "README.md":
                files_fit["fit/structures/README.md"] = scrub_structure_readme(text)
            else:
                raise LiteratureSplitError(f"structures/{q.name}: unknown file type, add a rule")
    if notation_src.is_file():
        files_fit["fit/site_notation.md"] = notation_src.read_text(encoding="utf-8")
    files_fit["fit/CONFLICTS.md"] = (AUTHORED / "CONFLICTS_fit.md").read_text(encoding="utf-8")
    cov_in = src_dir / "coverage_map.csv"
    cov_stats = None
    if cov_in.is_file():
        cov_text, kept, total = filter_coverage(cov_in.read_text(encoding="utf-8-sig"))
        files_fit["fit/coverage_map.csv"] = cov_text
        cov_stats = (kept, total)

    meta_common = {"decision": "literature/SPLIT_DECISION.md", "date": DATE, "author": AUTHOR,
                   "schema": "hexw.registry.schemas.LiteratureClaim (load with hexw.params.literature)"}
    fit_doc = {
        "meta": {**meta_common, "split": "fit", "n_claims": len(fit),
                 "n_estimates": sum(c.estimate is not None for c in fit),
                 "n_pending": sum(c.pending is not None for c in fit),
                 "n_qualitative": sum(c.estimate is None and c.pending is None for c in fit),
                 "note": "Only records with split == fit. L01: samples Mg and Ni only (Co/Zn columns removed from the "
                         "quoted source text). L11: Ba -> Sr transfer, tier 3. No holdout numbers."},
        "literature_records": _lit_records(doc, {c.lit_id for c in fit}, S.SourceRole.LITERATURE_FIT),
        "claims": [claim_plain(c) for c in fit],
    }
    hold_doc = {
        "meta": {**meta_common, "split": "holdout", "n_claims": len(hold_c),
                 "note": "HOLDOUT: not readable by theory code; read by hexw.validate after the freeze and hexw.audit."},
        "literature_records": _lit_records(doc, {c.lit_id for c in hold_c}, S.SourceRole.LITERATURE_HOLDOUT),
        "claims": [claim_plain(c) for c in hold_c],
    }
    excl_doc = {
        "meta": {**meta_common, "split": "excluded+methodology_only", "n_claims": len(excl_fit_side),
                 "note": "Never usable in a parametrization (SourceRole.NON_PARAMETRIC_ROLES)."},
        "claims": [claim_plain(c) for c in excl_fit_side],
    }
    excl_hold_doc = {
        "meta": {**meta_common, "split": "excluded (from holdout publications)", "n_claims": len(excl_hold_side)},
        "claims": [claim_plain(c) for c in excl_hold_side],
    }
    files_fit["fit/records_fit.yaml"] = _dump_yaml(fit_doc)
    files_fit["excluded/records_excluded.yaml"] = _dump_yaml(excl_doc)
    files_hold["records_holdout.yaml"] = _dump_yaml(hold_doc)
    files_hold["records_excluded_holdout_papers.yaml"] = _dump_yaml(excl_hold_doc)

    # ---- the fit side refers to no holdout source file and to no holdout sample (before anything is written)
    ref_hits = find_fit_references(cfg.root, cfg.holdout_sources, fit,
                                   {k: v for k, v in files_fit.items() if k.startswith(("fit/", "excluded/"))})
    if ref_hits:
        raise LeakError(f"fit outputs refer to holdout sources / samples (first 10): {ref_hits[:10]}")

    # ---- leak scan (before anything is written)
    toks = holdout_tokens(hold_c, [t for k, t in files_hold.items() if k.startswith("structures_holdout/")])
    # numbers that legitimately occur on the fit side too (coordinates, equal values of Mg/Ni rows, ...) are not
    # leaks: the structure CSVs are filtered row-by-row (deterministic), so they are the "legit" vocabulary, and the
    # scan looks for holdout-only numbers in the free text of fit outputs (README, quoted sources, notes).
    legit = set()
    for k, t in files_fit.items():
        if k.startswith("fit/structures/") and k.endswith(".csv"):
            legit |= set(_NUM_TOKEN.findall(t))
    for c in fit:
        legit |= set(_NUM_TOKEN.findall(c.value_original or ""))
        if c.estimate is not None:       # numbers of the fit estimates themselves (value, std, range)
            legit |= set(_NUM_TOKEN.findall(_dump_yaml(S.to_plain(c.estimate).get("uncertainty") or {})))
            legit |= set(_NUM_TOKEN.findall(str(c.estimate.value)))
    fit_texts = {k: v for k, v in files_fit.items()
                 if k.startswith(("fit/", "excluded/")) and not (k.startswith("fit/structures/") and k.endswith(".csv"))}
    hits = scan_leaks(toks, fit_texts, allow=legit)
    if hits:
        raise LeakError(f"holdout numbers found in fit outputs (first 10): {hits[:10]}")

    counts = _count_table(claims)
    summary = {
        "n_records": len(claims), "fit": len(fit), "holdout": len(hold_c), "excluded": len(excl),
        "excluded_fit_side": len(excl_fit_side), "excluded_holdout_side": len(excl_hold_side),
        "methodology_only": sum(c.role is S.SourceRole.METHODOLOGY_ONLY for c in claims),
        "by_source": counts["by_source"], "holdout_struct_rows": n_hold_struct, "coverage_rows_kept_total": cov_stats,
        "schema_estimates_fit": fit_doc["meta"]["n_estimates"], "leak_tokens_checked": len(toks),
        "candidate_problems": cand_problems, "first_run": first_run, "written": False,
    }
    if cand_problems:
        raise LiteratureSplitError(f"lit_search/candidates.csv does not match the declared holdout downloads: {cand_problems}")
    if not write:
        return summary

    # ---- write ----
    # 1. originals -> literature_holdout/_original_mixed (first run only)
    hold.mkdir(parents=True, exist_ok=True)
    if first_run:
        mixed.mkdir(parents=True, exist_ok=True)
        for name in ("records.yaml", "records.jsonl", "coverage_map.csv", "CONFLICTS.md", "READING_LOG.md",
                     "split_proposal.md", "structures", "scripts"):
            s_ = lit / name
            if s_.exists():
                shutil.move(str(s_), str(mixed / name))
        if (lit / "site_notation.md").is_file():          # common file: lives in fit/ now
            (lit / "site_notation.md").unlink()
        for q in _all_files(mixed):
            _readonly(q)
    # 2. fit / excluded / holdout files
    written = {}
    for rel, text in list(files_fit.items()) + [("holdout::" + k, v) for k, v in files_hold.items()]:
        base = hold if rel.startswith("holdout::") else lit
        rel2 = rel.replace("holdout::", "")
        p = base / rel2
        p.parent.mkdir(parents=True, exist_ok=True)
        _writable(p)
        p.write_bytes(text.replace("\r\n", "\n").encode("utf-8"))
        written[("holdout/" if base is hold else "") + rel2] = sha256_path(p)
    # 3. hashes of the holdout source files
    src_lines = []
    for sp in cfg.holdout_sources:
        rel = sp.relative_to(cfg.root).as_posix()
        if not sp.exists():
            src_lines.append(f"MISSING  {rel}")
            continue
        for q in _all_files(sp):
            src_lines.append(f"{sha256_path(q)}  {q.relative_to(cfg.root).as_posix()}")
    p = hold / "SOURCE_FILES_SHA256.txt"
    _writable(p)
    p.write_bytes(("\n".join(src_lines) + "\n").encode("utf-8"))
    # 4. SPLIT_DECISION.md + manifest + SHA256SUMS
    rec_sha = {
        "records_original_mixed_sha256": inputs["records.yaml"],
        "records_fit_sha256": written["fit/records_fit.yaml"],
        "records_holdout_sha256": written["holdout/records_holdout.yaml"],
        "records_excluded_sha256": written["excluded/records_excluded.yaml"],
    }
    decision = (AUTHORED / "SPLIT_DECISION.template.md").read_text(encoding="utf-8")
    decision = (decision.replace("{{COUNTS}}", counts["md"]).replace("{{INDEPENDENCE}}", INDEPENDENCE)
                .replace("{{HOLDOUT_PATHS}}", "\n".join(f"- `{sp.relative_to(cfg.root).as_posix()}`"
                                                        for sp in cfg.holdout_sources) or "- (none declared)")
                .replace("{{HASHES}}", "\n".join(f"{k}: {v}" for k, v in rec_sha.items())))
    (lit / "SPLIT_DECISION.md").write_bytes(decision.encode("utf-8"))
    manifest = {
        "date": DATE, "author": AUTHOR,
        "script": {"file": "hexw/ingest/literature_split.py", "sha256": sha256_path(Path(__file__))},
        "authored_inputs": {p.name: sha256_path(p) for p in sorted(AUTHORED.iterdir())},
        "inputs": dict(inputs), "hashes": rec_sha,
        "outputs": dict(sorted(written.items())),
        "holdout_sources": [sp.relative_to(cfg.root).as_posix() for sp in cfg.holdout_sources],
        "counts": {k: summary[k] for k in ("n_records", "fit", "holdout", "excluded", "methodology_only")},
        "by_source": summary["by_source"],
    }
    (lit / "split_manifest.json").write_bytes(json.dumps(manifest, ensure_ascii=False, indent=1).encode("utf-8"))
    _write_sums(lit, [q for q in sorted(lit.rglob("*")) if q.is_file() and ".git" not in q.relative_to(lit).parts
                      and q.name != "SHA256SUMS"])
    hold_files = [q for q in sorted(hold.rglob("*")) if q.is_file() and q.name != "SHA256SUMS"]
    _write_sums(hold, hold_files)
    for q in _all_files(hold):
        _readonly(q)
    summary["written"] = True
    summary["hashes"] = rec_sha
    return summary


def _write_sums(base: Path, files):
    lines = [f"{sha256_path(q)}  {q.relative_to(base).as_posix()}" for q in files]
    p = base / "SHA256SUMS"
    _writable(p)
    p.write_bytes(("\n".join(lines) + "\n").encode("utf-8"))


def _check_candidates(cfg: paths.Config) -> list:
    """C### <-> publication: DOIs in lit_search/results/candidates.csv must equal the declared ones."""
    f = cfg.root / "lit_search" / "results" / "candidates.csv"
    if not f.is_file():
        return []
    problems = []
    with open(f, encoding="utf-8-sig", newline="") as fh:
        rows = {r["id"]: r for r in csv.DictReader(fh)}
    for cid, doi in EXPECTED_DOI.items():
        r = rows.get(cid)
        if r is None:
            problems.append(f"{cid}: not in candidates.csv")
        elif (r.get("doi") or "").strip().lower() != doi.lower():
            problems.append(f"{cid}: candidates.csv has DOI {r.get('doi')!r}, expected {doi!r}")
    return problems


def _count_table(claims) -> dict:
    by = Counter()
    for c in claims:
        if c.lit_id == "L01":
            key = (c.lit_id + (" (Mg, Ni)" if c.role is S.SourceRole.LITERATURE_FIT else " (Co, Zn)"), c.role.value)
        else:
            key = (c.lit_id, c.role.value)
        by[key] += 1
    tot = Counter(c.role.value for c in claims)
    lines = ["Число записей по источникам:", "", "| Источник | роль | записей |", "|---|---|---|"]
    for (lid, role), n in sorted(by.items()):
        lines.append(f"| {lid} | {role} | {n} |")
    lines += ["", "Итого: " + ", ".join(f"{r} {n}" for r, n in sorted(tot.items())) + f"; всего {len(claims)}."]
    return {"md": "\n".join(lines), "by_source": {f"{a}|{b}": n for (a, b), n in sorted(by.items())}}


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="Split the literature base into fit / holdout (SPLIT_DECISION.md)")
    ap.add_argument("--force", action="store_true", help="rebuild outputs that already exist (from _original_mixed)")
    ap.add_argument("--dry-run", action="store_true", help="classify and scan for leaks, write nothing")
    ns = ap.parse_args(argv)
    try:
        res = run(write=not ns.dry_run, force=ns.force)
    except LiteratureSplitError as e:
        print(f"literature_split: {e}", file=sys.stderr)
        return 2
    print(json.dumps(res, ensure_ascii=False, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
