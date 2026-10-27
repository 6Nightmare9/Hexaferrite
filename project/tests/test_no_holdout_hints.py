"""Stage 4 round 1: the theory side carries no hint at the holdout content (the L01-vs-L02 Co/6g affair of review 4) and does not use the
full lit_search catalogue. Mentioning that a holdout EXISTS (isolation docs, the guard) is fine."""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
THEORY = [ROOT / "hexw" / p for p in ("params", "structures", "theory", "dft", "spin", "chem")]
DOCS = [ROOT / "docs" / n for n in ("STRUCTURE_MODEL.md", "ANISOTROPY_MODEL.md", "MAGNET_MODEL.md", "DFT.md")]

#: phrases that point at holdout samples / conflicts (case-insensitive)
FORBIDDEN = [
    r"L01[ -](?:Co|Zn)\b", r"L01\s*\(\s*Co", r"L01\s+Co\s*/\s*Zn", r"\bCo\s*/\s*Zn\s+(?:samples?\s+)?(?:of|in|from)\s+L01",
    r"L01\s*(?:vs\.?|versus|и|and|against)\s*L02", r"L01-vs-L02", r"L02\s*(?:vs\.?|versus)\s*L01",
    r"конфликт[^\n]{0,80}(?:holdout|L01)", r"(?:holdout|L01)[^\n]{0,80}конфликт",
    r"conflict[^\n]{0,80}(?:holdout|L01)", r"(?:holdout|L01)[^\n]{0,80}conflict",
    r"(?:проверяется|проверяют|checked)\s+(?:на|on)\s+holdout", r"holdout[- ]side", r"holdout\s+quantity",
    r"holdout-?(?:значени|образц)", r"L01\s*\(\s*Co,\s*Zn", r"Co,\s*Zn\s+(?:holdout|образц[^\n]{0,20}L01)",
]


def _files():
    for base in THEORY:
        if base.is_dir():
            yield from (q for q in base.rglob("*.py") if "__pycache__" not in q.parts)
    for q in DOCS:
        if q.is_file():
            yield q


def test_no_phrase_points_at_holdout_content():
    pats = [re.compile(p, re.I) for p in FORBIDDEN]
    bad = []
    files = list(_files())
    assert files
    for q in files:
        for i, line in enumerate(q.read_text(encoding="utf-8", errors="ignore").splitlines(), 1):
            for p in pats:
                if p.search(line):
                    bad.append(f"{q.relative_to(ROOT).as_posix()}:{i}: {line.strip()[:100]}")
    assert not bad, "\n".join(bad)


def test_the_forbidden_patterns_do_catch_the_old_wording():
    old = ["Conflict on Co at W3 (6g) between L01 and L02", "Co on W3: the L01-vs-L02 conflict is a HOLDOUT-side quantity",
           "Конфликт по Co на W3 (6g) между L01 и L02", "проверяется на holdout (L01 Co)", "Holdout (L01 Co/Zn и др.)"]
    pats = [re.compile(p, re.I) for p in FORBIDDEN]
    assert all(any(p.search(t) for p in pats) for t in old)
    ok = ["L02 Co/Zn, L01 Mg/Ni", "the holdout exists and is blocked by the guard"]
    assert not any(p.search(t) for p in pats for t in ok)


def test_the_neutral_wording_is_in_place():
    ens = (ROOT / "hexw" / "structures" / "ensemble.py").read_text(encoding="utf-8")
    occ = (ROOT / "hexw" / "params" / "occupancy.py").read_text(encoding="utf-8")
    doc = (ROOT / "docs" / "STRUCTURE_MODEL.md").read_text(encoding="utf-8")
    assert "the only source is L02; there is no cross-check in the fit" in ens
    assert "the only source is L02; there is no cross-check in the fit" in occ
    assert "единственный источник - L02; перекрёстной проверки в fit нет" in doc


def test_the_theory_side_does_not_use_the_full_search_catalogue():
    """candidates.csv (reason / topic_tags from abstracts) is holdout-derived; only candidates_meta.csv may be read, and only by
    ingest / audit / tests."""
    hits = []
    for base in THEORY + [ROOT / "hexw" / "util", ROOT / "hexw" / "validate", ROOT / "hexw" / "experiment"]:
        if not base.is_dir():
            continue
        for q in base.rglob("*.py"):
            t = q.read_text(encoding="utf-8", errors="ignore")
            if re.search(r"candidates\.csv|rows_final|rows_stage2|merged\.json|raw_hits", t):
                hits.append(q.relative_to(ROOT).as_posix())
    assert not hits, hits


def test_the_config_declares_the_search_results_as_holdout():
    cfg = json.loads((ROOT / "config" / "holdout_sources.json").read_text(encoding="utf-8"))
    src = set(cfg["holdout_sources"])
    for rel in ("rows_final.json", "rows_stage2.json", "merged.json", "raw_hits.json", "raw_hits2.json", "raw_hits3.json",
                "raw_hits4.json", "coverage_from_abstracts.csv", "coverage_matrix.csv", "review_list.txt", "merged_list.txt",
                "candidates.csv"):
        assert f"lit_search/results/{rel}" in src, rel
    assert "lit_search/results/candidates_meta.csv" not in src
    assert "lit_search/results/candidates.csv" not in cfg["reference_allowlist"]
    pins = cfg["reference_allowlist_sha256"]
    assert "lit_search/results/candidates_meta.csv" in pins
    assert all(re.fullmatch(r"[0-9a-f]{64}", v) for v in pins.values())
    for name in ("download", "download2", "epmc_xml", "final"):
        assert f"lit_search/scripts/{name}.py" in pins
        assert f"lit_search/scripts/{name}.py" not in cfg["reference_allowlist"]       # pinned, not free
