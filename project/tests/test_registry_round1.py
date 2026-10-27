import csv
import io
import json
import os
import stat
import subprocess

import pytest
import yaml

from hexw import audit
from hexw.io import paths
from hexw.registry import integrity as I
from hexw.registry import prior_use as PU
from hexw.registry.measurements import iter_measurements


def real_registry():
    reg = paths._default_root() / "registry"
    if not (reg / "observations.yaml").exists() or not (reg / ".git").exists():
        pytest.skip("real registry not populated")
    return reg


# ---- 13: iter_measurements
def test_iter_measurements_default_excludes_prior(tmp_path):
    obs = [{"observation_id": "OBS-XRD-001", "method": "XRD", "independent": True, "files": []},
           {"observation_id": "OBS-MOS-001", "method": "Mossbauer", "independent": True, "files": []},
           {"observation_id": "PRI-001", "method": "prior_analysis_artifact", "role": "prior_lab_analysis",
            "independent": False, "files": []},
           {"observation_id": "PRI-012", "method": "prior_analysis_artifact", "role": "prior_lab_analysis",
            "independent": False, "files": []}]
    (tmp_path / "observations.yaml").write_text(yaml.safe_dump({"observations": obs}), encoding="utf-8")
    assert [o["observation_id"] for o in iter_measurements(tmp_path)] == ["OBS-XRD-001", "OBS-MOS-001"]
    assert [o["observation_id"] for o in iter_measurements(tmp_path, methods=["XRD"])] == ["OBS-XRD-001"]
    assert len(list(iter_measurements(tmp_path, include_prior=True))) == 4


def test_real_registry_measurements_and_mtd_registration():
    reg = real_registry()
    ms = list(iter_measurements(reg))
    assert len(ms) == 13 and all(m["method"] in PU.LAB_METHODS and m.get("role") != "prior_lab_analysis" for m in ms)
    allobs = {o["observation_id"]: o for o in iter_measurements(reg, include_prior=True)}
    for pid in ("PRI-012", "PRI-013", "PRI-014"):
        o = allobs[pid]
        assert o["role"] == "prior_lab_analysis" and o["independent"] is False
        assert all(f["sha256"] for f in o["files"])
    assert not any(str(f.get("archive_path", "")).endswith(".mtd") for m in ms for f in m["files"])


def test_real_registry_user_reported_mossbauer_metadata():
    reg = real_registry()
    n = 0
    for o in iter_measurements(reg, methods=["Mossbauer"]):
        n += 1
        u = o["user_reported_metadata"]
        assert u["status"] == "user_reported" and u["reported_date"] == "2026-10-25"
        assert u["source_radionuclide"]["value"] == "Co-57" and u["source_radionuclide"]["matrix"] == "unknown"
        assert u["temperature"]["unit"] == "K"
        assert u["temperature"]["uncertainty"]["basis"].startswith("assumed_by_registry")
        assert u["absorber_area"]["status"] == "unknown"
        assert o["calibration"]["calibration_status"] == "pending"
    assert n == 5


# ---- 11: different_version is really checked
def test_different_version_statuses_in_real_manifest():
    reg = real_registry()
    rows = list(csv.DictReader(open(reg / "files_manifest_umnik.csv", encoding="utf-8", newline="")))
    dv = [r for r in rows if r["different_version_of"]]
    assert dv and all(r["status"] == "different_version" for r in dv)
    assert all(r["status"] != "different_version" for r in rows if not r["different_version_of"])


def _umnik_checks(mutator):
    """check_umnik on the real registry with an in-memory mutated files_manifest_umnik.csv."""
    cfg = paths.get_config()
    rows = list(csv.DictReader(open(cfg.registry / "files_manifest_umnik.csv", encoding="utf-8", newline="")))
    mutator(rows)
    orig = audit._rows

    def fake(path):
        return rows if os.path.basename(str(path)) == "files_manifest_umnik.csv" else orig(path)
    audit._rows = fake
    try:
        rep = audit.Report()
        audit.check_umnik(rep, cfg)
    finally:
        audit._rows = orig
    return {c["name"]: c["ok"] for c in rep.checks}


def test_different_version_check_fires():
    real_registry()
    key = ("umnik: status different_version <=> different_version_of, reference resolves to another file with "
           "different bytes")
    assert _umnik_checks(lambda rows: None)[key]

    def wrong_status(rows):
        next(r for r in rows if r["different_version_of"])["status"] = "new"

    def status_without_ref(rows):
        next(r for r in rows if not r["different_version_of"])["status"] = "different_version"

    def dangling(rows):
        next(r for r in rows if r["different_version_of"])["different_version_of"] = (
            "umnik:No/Such/file.png sha256=" + "0" * 64)
    for m in (wrong_status, status_without_ref, dangling):
        assert not _umnik_checks(m)[key], m.__name__


# ---- 14: audit writes nothing by default
def test_audit_default_writes_nothing_and_flag_writes_manifest():
    reg = real_registry()
    runs = reg / "runs"
    before = set(runs.iterdir()) if runs.exists() else set()
    buf = io.StringIO()
    assert audit.run_audit(None, out=buf) == 0, buf.getvalue()
    assert (set(runs.iterdir()) if runs.exists() else set()) == before
    buf = io.StringIO()
    assert audit.run_audit(None, write_manifest=True, out=buf) == 0, buf.getvalue()
    new = set(runs.iterdir()) - before
    assert len(new) == 1
    m = json.loads(next(iter(new)).read_text(encoding="utf-8"))
    assert m["git_head"]["project"] and m["git_head"]["registry"]
    assert m["checksums"] and "samples.yaml" in m["checksums"]
    assert I.git_dirty(reg) == []                       # runs/ is git-ignored: manifests do not dirty the tree
    for f in new:
        os.chmod(f, stat.S_IWRITE)
        f.unlink()


def test_cli_audit_flags():
    from hexw import cli
    assert cli.build_parser().parse_args(["audit"]).write_manifest is False
    assert cli.build_parser().parse_args(["audit", "--write-manifest"]).write_manifest is True


# ---- 15: registry in git
def test_sha256sums_detect_tampering_and_git_dirty(tmp_path):
    reg = tmp_path / "reg"
    reg.mkdir()
    (reg / "a.yaml").write_text("a: 1\n")
    (reg / "runs").mkdir()
    (reg / "runs" / "x.json").write_text("{}")
    I.write_sums(reg)
    assert I.verify_sums(reg) == [] and set(I.read_sums(reg)) == {"a.yaml"}
    (reg / "a.yaml").write_text("a: 2\n")
    assert I.verify_sums(reg) == ["hash mismatch a.yaml"]
    (reg / "b.yaml").write_text("b\n")
    assert "unlisted file b.yaml" in I.verify_sums(reg)
    subprocess.run(["git", "init", "-q", str(reg)], check=True)
    assert I.git_dirty(reg)                              # untracked files -> dirty
    assert I.git_head(reg) is None                       # no commit yet


def test_real_registry_is_git_clean_with_sums_and_readonly_keys():
    reg = real_registry()
    assert I.git_head(reg) and I.git_dirty(reg) == [] and I.verify_sums(reg) == []
    for n in I.KEY_READONLY:
        assert audit._is_readonly(reg / n)
    assert "runs/" in (reg / ".gitignore").read_text(encoding="utf-8")
    tracked = subprocess.run(["git", "-C", str(reg), "ls-files"], capture_output=True, text=True).stdout.split()
    assert not any(t.lower().endswith((".spc", ".mdi", ".mtd", ".jip", ".zip", ".png")) for t in tracked)
