"""Stage 3: fit / holdout split of the literature base and protection of the HOLDOUT zone.

Synthetic root for the access rules; the real (already built) literature/ for the content checks (read only)."""
import hashlib
import json
import re
from pathlib import Path

import pytest
import yaml

from _helpers import run_as
from hexw.io import paths
from hexw.registry import schemas as S

THEORY_MODS = ("hexw.theory._probe", "hexw.params._probe", "hexw.dft._probe", "hexw.structures._probe",
               "hexw.spin._probe", "hexw.util._probe")
OPEN_SRC = "def go(p):\n    return open(p, 'rb').read()\n"


# ---------------------------------------------------------------- synthetic root
@pytest.fixture
def hroot(tmp_path):
    (tmp_path / "raw" / "theory_side" / "06_literature" / "L01").mkdir(parents=True)
    (tmp_path / "raw" / "theory_side" / "06_literature" / "L02").mkdir(parents=True)
    (tmp_path / "raw" / "original_zip").mkdir(parents=True)
    (tmp_path / "vault" / "lab_validation").mkdir(parents=True)
    (tmp_path / "registry").mkdir()
    (tmp_path / "lit_search" / "downloads" / "C056").mkdir(parents=True)
    (tmp_path / "literature" / "fit").mkdir(parents=True)
    (tmp_path / "literature_holdout").mkdir()
    (tmp_path / "project" / "config").mkdir(parents=True)
    (tmp_path / "project" / "config" / "holdout_sources.json").write_text(json.dumps(
        {"holdout_sources": ["raw/theory_side/06_literature/L01", "lit_search/downloads/C056"]}))
    files = {"raw/theory_side/06_literature/L01/a.cif": "cif", "raw/theory_side/06_literature/L01/paper.pdf": "pdf",
             "raw/theory_side/06_literature/L02/paper.pdf": "fit paper", "lit_search/downloads/C056/C056.pdf": "pdf",
             "literature_holdout/records_holdout.yaml": "claims: []", "literature/fit/records_fit.yaml": "claims: []"}
    for rel, txt in files.items():
        (tmp_path / rel).write_text(txt)
    with paths.test_root(tmp_path):
        yield tmp_path


@pytest.mark.parametrize("rel", ["literature_holdout/records_holdout.yaml",
                                 "raw/theory_side/06_literature/L01/paper.pdf",
                                 "raw/theory_side/06_literature/L01/a.cif",
                                 "lit_search/downloads/C056/C056.pdf"])
@pytest.mark.parametrize("mod", THEORY_MODS)
def test_theory_cannot_read_holdout(hroot, rel, mod):
    s = run_as(mod, OPEN_SRC)
    with pytest.raises(paths.TheoryAccessError):
        s["go"](hroot / rel)


def test_theory_cannot_list_or_use_open_data_on_holdout(hroot):
    s = run_as("hexw.theory._probe", "import os\nfrom hexw.io import paths\n"
               "def ls(p): return os.listdir(p)\n"
               "def od(p): return paths.open_data(p, 'rb', purpose='metadata')\n")
    with pytest.raises(paths.TheoryAccessError):
        s["ls"](hroot / "literature_holdout")
    with pytest.raises(paths.TheoryAccessError):
        s["od"](hroot / "literature_holdout" / "records_holdout.yaml")


def test_fit_files_and_fit_papers_stay_readable_for_theory(hroot):
    s = run_as("hexw.theory._probe", OPEN_SRC)
    assert s["go"](hroot / "literature" / "fit" / "records_fit.yaml") == b"claims: []"
    assert s["go"](hroot / "raw/theory_side/06_literature/L02/paper.pdf") == b"fit paper"


def test_zone_classification(hroot):
    assert paths.classify(hroot / "literature_holdout") is paths.Zone.HOLDOUT
    assert paths.classify(hroot / "raw/theory_side/06_literature/L01/a.cif") is paths.Zone.HOLDOUT
    assert paths.classify(hroot / "raw/theory_side/06_literature/L02/paper.pdf") is paths.Zone.RAW_THEORY
    assert paths.classify(hroot / "literature/fit/records_fit.yaml") is paths.Zone.OTHER


def test_holdout_is_read_only_and_needs_freeze_for_validate(hroot):
    with pytest.raises(paths.ReadOnlyZoneError):
        paths.open_data(hroot / "literature_holdout" / "x.yaml", "wb")
    s = run_as("hexw.validate._probe", "from hexw.io import paths\n"
               "def go(p): return paths.open_data(p, 'rb').read()\n")
    with pytest.raises(paths.FreezeRequiredError):
        s["go"](hroot / "literature_holdout" / "records_holdout.yaml")
    a = run_as("hexw.audit._probe", "from hexw.io import paths\n"
               "def go(p): return paths.open_data(p, 'rb', purpose='integrity').read()\n")
    assert a["go"](hroot / "literature_holdout" / "records_holdout.yaml") == b"claims: []"


def test_experiment_and_ingest_alone_cannot_read_holdout(hroot):
    for mod in ("hexw.experiment._probe", "hexw.ingest._probe", "hexw.io._probe"):
        s = run_as(mod, "from hexw.io import paths\ndef go(p): return paths.open_data(p, 'rb').read()\n")
        with pytest.raises(paths.AccessError):
            s["go"](hroot / "literature_holdout" / "records_holdout.yaml")


def test_holdout_build_only_for_ingest(hroot):
    th = run_as("hexw.theory._probe", "from hexw.io import paths\n"
                "def go():\n    with paths.holdout_build():\n        pass\n")
    with pytest.raises(paths.TheoryAccessError):
        th["go"]()
    with pytest.raises(paths.AccessError):                       # this test module is not hexw.ingest
        with paths.holdout_build():
            pass
    ing = run_as("hexw.ingest._probe", "from hexw.io import paths\n"
                 "def go(p):\n    with paths.holdout_build():\n        return open(p, 'rb').read()\n")
    assert ing["go"](hroot / "literature_holdout" / "records_holdout.yaml") == b"claims: []"


def test_undeclared_sources_with_holdout_dir_is_an_error(tmp_path):
    (tmp_path / "vault").mkdir()
    (tmp_path / "literature_holdout").mkdir()
    with pytest.raises(paths.RootNotFoundError):
        with paths.test_root(tmp_path):
            pass


# ---------------------------------------------------------------- real, already built literature/
REAL = Path(__file__).resolve().parents[2]


def _need_real():
    if not (REAL / "literature" / "fit" / "records_fit.yaml").is_file():
        pytest.skip("literature split not built")


def test_load_fit_records_from_theory_side():
    _need_real()
    s = run_as("hexw.theory._probe", "from hexw.params.literature import load_fit_records\n"
               "def go(): return load_fit_records()\n")
    recs = s["go"]()
    assert len(recs) > 300
    assert all(c.split == "fit" and c.role is S.SourceRole.LITERATURE_FIT for c in recs)
    assert {c.lit_id for c in recs} == {"L01", "L02", "L03", "L05", "L11", "C004"}


def test_real_holdout_zone_blocked_for_theory():
    _need_real()
    s = run_as("hexw.theory._probe", "def go(p):\n    return open(p, 'rb').read(4)\n")
    l01 = REAL / "raw/theory_side/06_literature/L01_Morch2019"
    with paths.test_access():                         # enumeration only (listing errors are swallowed by pathlib)
        found = ([q for q in l01.rglob("*") if q.is_file()][:3] + sorted((REAL / "dft_feasibility/pilots/cif").glob("*.cif"))
                 + sorted((REAL / "lit_search/downloads/C023").glob("*"))[:1])
    targets = [REAL / "literature_holdout" / "records_holdout.yaml",
               REAL / "literature_holdout" / "_original_mixed" / "records.yaml",
               *found]
    assert len(targets) >= 5
    for t in targets:
        with pytest.raises(paths.TheoryAccessError):
            s["go"](t)


def test_fit_content_rules():
    _need_real()
    text = (REAL / "literature/fit/records_fit.yaml").read_text(encoding="utf-8")
    doc = yaml.safe_load(text)
    assert all(c["role"] == "literature_fit" and c["split"] == "fit" for c in doc["claims"])
    l01 = [c for c in doc["claims"] if c["lit_id"] == "L01"]
    assert l01 and {c["source"]["sample_id_in_source"] for c in l01} <= {"SrMg2Fe16O27", "SrNi2Fe16O27"}
    l11 = [c for c in doc["claims"] if c["lit_id"] == "L11"]
    assert l11 and all(c["source"]["tier"] == 3 and c["source"]["transfer"] == "Ba_to_Sr" for c in l11)
    assert not re.search(r"32\.7411|32\.8293|5\.90922", text)       # Co/Zn lattice constants of L01 (holdout)


def test_holdout_assignment():
    _need_real()
    with paths.test_access():
        hd = yaml.safe_load((REAL / "literature_holdout/records_holdout.yaml").read_text(encoding="utf-8"))
    kinds = {}
    for c in hd["claims"]:
        kinds.setdefault(c["lit_id"], set()).add(c["source"]["holdout_kind"])
    assert kinds["L01"] == {"same_group_reproducibility"}
    assert kinds["C023"] == kinds["C015"] == {"independent_group"}
    for k in ("L12", "C057", "C060", "C044", "C056", "C036", "C052", "C066"):
        assert kinds[k] == {"transfer_Ba_to_Sr"}, k
    assert "L11" not in kinds


def test_excluded_cannot_become_parameters():
    _need_real()
    ex = yaml.safe_load((REAL / "literature/excluded/records_excluded.yaml").read_text(encoding="utf-8"))
    roles = {c["role"] for c in ex["claims"]}
    assert roles <= {"excluded_secondary", "excluded_low_quality", "methodology_only"}
    assert "methodology_only" in roles and "excluded_secondary" in roles
    assert all(not c.get("estimate") for c in ex["claims"])
    mo = {c["lit_id"] for c in ex["claims"] if c["role"] == "methodology_only"}
    assert {"C112", "C115"} <= mo


def test_split_manifest_hashes_match_files():
    _need_real()
    man = json.loads((REAL / "literature/split_manifest.json").read_text(encoding="utf-8"))

    def sha(p):
        return hashlib.sha256(p.read_bytes()).hexdigest()
    assert man["hashes"]["records_fit_sha256"] == sha(REAL / "literature/fit/records_fit.yaml")
    assert man["hashes"]["records_excluded_sha256"] == sha(REAL / "literature/excluded/records_excluded.yaml")
    dec = (REAL / "literature/SPLIT_DECISION.md").read_text(encoding="utf-8")
    assert man["hashes"]["records_fit_sha256"] in dec and "2026-10-25" in dec


def test_split_script_refuses_to_overwrite_without_force():
    _need_real()
    from hexw.ingest import literature_split as LS
    with pytest.raises(LS.OutputsExistError):
        LS.run(write=True, force=False)
