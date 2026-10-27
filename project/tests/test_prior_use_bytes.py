"""prior_use.yaml is hand-edited: the generator must keep the existing bytes (line endings, BOM, missing final newline)
and only APPEND missing entries (stage 3 round 1; regression of the independent tester). Synthetic data only."""
from __future__ import annotations

import csv
import os
import stat
import subprocess

import pytest
import yaml

from test_umnik_prior_use import umnik_root, unlocked  # noqa: F401  (fixtures)
from hexw.ingest import registry_guard as G
from hexw.ingest import umnik

MANUAL = ("# header kept — ручной комментарий\n"
          "schema_version: 1\nindependent_lab_validation_available: false\n"
          "independent_lab_validation_note: set by hand\nentries:\n"
          "- observation_id: OBS-XRD-001\n  sample_id: SMP-001\n  method: XRD\n  prior_use_status: used_prior\n"
          "  validation_mode: retrospective\n"
          "  user_confirmation: {by: user, date: 2026-10-25}   # hand\n")
NEW = [{"observation_id": "OBS-XRD-002", "sample_id": "SMP-001", "method": "XRD",
        "prior_use_status": "unknown_prior_use", "validation_mode": "retrospective"}]
GEN = [{"observation_id": "OBS-XRD-001", "prior_use_status": "unknown_prior_use",
        "validation_mode": "retrospective"}] + NEW
VARIANTS = [("lf", b"\n", False, True), ("crlf", b"\r\n", False, True), ("lf_bom", b"\n", True, True),
            ("crlf_bom", b"\r\n", True, True), ("no_final_eol", b"\n", False, False),
            ("crlf_no_final_eol", b"\r\n", False, False)]


def make_raw(nl: bytes, bom: bool, trail: bool) -> bytes:
    raw = MANUAL.replace("\n", nl.decode()).encode("utf-8")
    if not trail:
        raw = raw.rstrip(b"\r\n")
    return (b"\xef\xbb\xbf" if bom else b"") + raw


def same_eol_style(tail: bytes, nl: bytes) -> bool:
    return tail.count(b"\r\n") == tail.count(b"\n") if nl == b"\r\n" else b"\r" not in tail


@pytest.mark.parametrize("name,nl,bom,trail", VARIANTS, ids=[v[0] for v in VARIANTS])
def test_merge_and_write_keep_existing_bytes(tmp_path, name, nl, bom, trail):
    raw = make_raw(nl, bom, trail)
    p = tmp_path / "prior_use.yaml"
    p.write_bytes(raw)
    text, rep = G.merge_prior_use_text(G.read_text_exact(p), GEN)
    assert rep["added"] == ["OBS-XRD-002"]
    G.write_text_keep_mode(p, text, newline="")
    new = p.read_bytes()
    assert new.startswith(raw) and len(new) > len(raw)                 # byte-for-byte prefix
    tail = new[len(raw):]
    assert same_eol_style(tail, nl) and b"OBS-XRD-002" in tail
    if not trail:
        assert tail.startswith(nl)                                      # only an EOL is added in front of the block
    doc = yaml.safe_load(new.decode("utf-8-sig"))
    assert [e["observation_id"] for e in doc["entries"]] == ["OBS-XRD-001", "OBS-XRD-002"]
    # second run: nothing to add, text identical, file left alone
    text2, rep2 = G.merge_prior_use_text(G.read_text_exact(p), GEN)
    assert rep2["added"] == [] and text2 == G.read_text_exact(p)


@pytest.mark.parametrize("nl", [b"\n", b"\r\n"])
def test_write_text_keep_mode_keeps_the_eol_style_of_the_replaced_file(tmp_path, nl):
    p = tmp_path / "x.yaml"
    p.write_bytes(b"a: 1" + nl + b"b: 2" + nl)
    G.write_text_keep_mode(p, "a: 3\nb: 4\n")
    assert p.read_bytes() == b"a: 3" + nl + b"b: 4" + nl
    q = tmp_path / "new.yaml"
    G.write_text_keep_mode(q, "a: 1\n")
    assert q.read_bytes() == b"a: 1\n"                                  # a new file is LF, never CRLF by platform


def test_read_only_flag_restored_and_bytes_exact(tmp_path):
    p = tmp_path / "ro.yaml"
    p.write_bytes(b"a: 1\r\n")
    os.chmod(p, stat.S_IREAD)
    G.write_text_keep_mode(p, "a: 1\r\nb: 2\r\n", newline="")
    assert p.read_bytes() == b"a: 1\r\nb: 2\r\n" and not os.access(p, os.W_OK)
    os.chmod(p, stat.S_IWRITE | stat.S_IREAD)


def git(reg, *a):
    return subprocess.run(["git", "-C", str(reg), "-c", "user.email=t@t", "-c", "user.name=t", *a], check=True,
                          capture_output=True)


@pytest.mark.parametrize("name,nl,bom,trail", VARIANTS, ids=[v[0] for v in VARIANTS])
def test_umnik_build_keeps_prior_use_prefix(umnik_root, name, nl, bom, trail):      # noqa: F811
    reg = umnik_root / "registry"
    raw = make_raw(nl, bom, trail)
    (reg / "prior_use.yaml").write_bytes(raw)
    obs = yaml.safe_load((reg / "observations.yaml").read_text(encoding="utf-8"))["observations"]
    obs.append({"observation_id": "OBS-XRD-002", "sample_id": "SMP-001", "method": "XRD", "independent": True,
                "files": [{"archive_path": "01_XRD/b.txt", "sha256": "00"}]})
    (reg / "observations.yaml").write_text(yaml.safe_dump({"observations": obs}), encoding="utf-8")
    (reg / "samples.yaml").write_text(yaml.safe_dump({"samples": [
        {"sample_id": "SMP-001", "nominal_formula": {"A": {"value": "Co"}, "x": {"value": 1.0}}}]}), encoding="utf-8")
    rows = list(csv.DictReader(open(reg / "files_manifest.csv", encoding="utf-8", newline="")))
    for r in rows:
        r["data_manifest_role"] = "raw"
    with open(reg / "files_manifest.csv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with unlocked(umnik_root):
        (umnik_root / "vault" / "lab_validation_umnik" / "Mn" / "a.jip").unlink()
    git(reg, "init", "-q")
    git(reg, "add", "-A")
    git(reg, "commit", "-q", "-m", "x")
    rep = umnik.build(force=True)
    assert "OBS-XRD-002" in rep["prior_use"]["added"]
    new = (reg / "prior_use.yaml").read_bytes()
    assert new.startswith(raw) and len(new) > len(raw), "existing prior_use bytes altered"
    assert same_eol_style(new[len(raw):], nl)
