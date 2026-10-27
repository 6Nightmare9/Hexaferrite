import contextlib
import csv
import hashlib
import io
import json

import pytest
import yaml

from hexw import audit
from hexw.io import paths
from hexw.registry import prior_use as PU
from hexw.registry import schemas as S
from hexw.registry.versioning import ValidationStore


# ---------------- ValidationResult carries the prior-use status ----------------
def _vr(**kw):
    return S.ValidationResult("r1", "mv", "p1", [S.DataRef("f", group="OBS-XRD-001")], S.Metric("e"), **kw)


def test_validation_result_defaults_to_retrospective():
    r = _vr()
    assert r.validation_mode is S.ValidationMode.RETROSPECTIVE and r.prior_use_status == "unknown_prior_use"


def test_validation_result_independent_needs_no_prior_use_found():
    with pytest.raises(ValueError):
        _vr(validation_mode="independent_candidate")
    with pytest.raises(ValueError):
        _vr(validation_mode="independent_candidate", prior_use_status="used_prior")
    with pytest.raises(ValueError):
        _vr(validation_mode="blind")
    with pytest.raises(ValueError):
        _vr(prior_use_status="whatever")
    ok = _vr(validation_mode="independent_candidate", prior_use_status="no_prior_use_found")
    assert ok.validation_mode is S.ValidationMode.INDEPENDENT_CANDIDATE


def test_validation_mode_survives_json_round_trip_and_store(tmp_path):
    r = _vr(validation_mode="retrospective", prior_use_status="used_prior")
    d = json.loads(S.dump_json(r))
    assert d["validation_mode"] == "retrospective" and d["prior_use_status"] == "used_prior"
    assert S.from_dict(S.ValidationResult, d) == r
    st = ValidationStore(tmp_path / "v.jsonl")
    st.append(r,S.Prediction("p1","mv",[S.DataRef("fit-input",source_role="literature_fit")]))
    rec = st.all()[0]["record"]
    assert rec["validation_mode"] == "retrospective" and rec["prior_use_status"] == "used_prior"


# ---------------- prior_use rules ----------------
E = {
    "OBS-A": {"prior_use_status": "used_prior", "validation_mode": "retrospective"},
    "OBS-M": {"prior_use_status": "unknown_prior_use", "validation_mode": "retrospective"},
    "OBS-N": {"prior_use_status": "no_prior_use_found", "validation_mode": "independent_candidate",
              "user_confirmation": {"by": "user"}},
}


def test_problems_detect_inconsistent_entries():
    assert PU.problems(E["OBS-A"]) == [] and PU.problems(E["OBS-N"]) == []
    assert PU.problems({"prior_use_status": "used_prior", "validation_mode": "independent_candidate"})
    nc = dict(E["OBS-N"], user_confirmation=None)
    assert any("user_confirmation" in p for p in PU.problems(nc))
    assert PU.problems({"prior_use_status": "unknown_prior_use", "validation_mode": "independent_candidate"})


def test_mode_is_independent_only_if_every_target_is():
    assert PU.mode_for_observations(E, ["OBS-N"]) == (S.ValidationMode.INDEPENDENT_CANDIDATE, "no_prior_use_found")
    assert PU.mode_for_observations(E, ["OBS-N", "OBS-A"]) == (S.ValidationMode.RETROSPECTIVE, "used_prior")
    assert PU.mode_for_observations(E, ["OBS-N", "OBS-M"]) == (S.ValidationMode.RETROSPECTIVE, "unknown_prior_use")
    # a missing record is never treated as independent
    assert PU.mode_for_observations(E, ["OBS-N", "OBS-?"])[0] is S.ValidationMode.RETROSPECTIVE
    assert PU.mode_for_observations(E, [])[0] is S.ValidationMode.RETROSPECTIVE


def test_apply_prior_use_sets_result_fields():
    r = S.ValidationResult("r", "mv", "p", [S.DataRef("x.txt", group="OBS-A")], S.Metric("e"))
    PU.apply_prior_use(r, E)
    assert (r.validation_mode, r.prior_use_status) == (S.ValidationMode.RETROSPECTIVE, "used_prior")
    r2 = S.ValidationResult("r2", "mv", "p", [S.DataRef("OBS-N")], S.Metric("e"))
    PU.apply_prior_use(r2, E)
    assert r2.validation_mode is S.ValidationMode.INDEPENDENT_CANDIDATE


def test_real_registry_has_prior_use_for_every_lab_observation():
    reg = paths._default_root() / "registry"
    if not (reg / "prior_use.yaml").exists():
        pytest.skip("real registry not populated")
    entries = PU.load_prior_use(reg)
    obs = yaml.safe_load((reg / "observations.yaml").read_text(encoding="utf-8"))["observations"]
    lab = {o["observation_id"]: o for o in obs if o["method"] in PU.LAB_METHODS}
    assert set(entries) == set(lab)
    for oid, e in entries.items():
        assert PU.problems(e) == []
        assert e["validation_mode"] == "retrospective"        # nothing is confirmed independent yet
        assert e["prior_use_status"] == "used_prior"
        if lab[oid]["method"] == "Mossbauer":                  # confirmed by the user on 2026-10-25
            assert e["user_confirmation"]["source"] == "user_confirmation"
            assert e["user_confirmation"]["date"] == "2026-10-25"
    assert PU.independent_validation_available(reg) is False   # no independent laboratory validation left


# ---------------- audit.check_umnik on a synthetic tree ----------------
def _sha(b):
    return hashlib.sha256(b).hexdigest()


@contextlib.contextmanager
def unlocked(root):
    """The audit hook forbids writes into vault/ of the configured root; tests build fixtures with it pointed away."""
    with paths.test_root(root / "_elsewhere", validate=False):
        yield


def _w(path, data: bytes):
    with unlocked(_ROOT[0]):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    return _sha(data)


_ROOT = [None]


def _csv(path, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


FM_COLS = ["zone", "archive_path", "sha256", "sample_ids", "observation_ids"]
UM_COLS = ["zone", "umnik_path", "bytes", "sha256", "role", "status", "matched_registered_paths",
           "name_variant_of_archive", "different_version_of", "sample_ids", "observation_ids",
           "byte_dup_group_umnik", "origin", "note"]


@pytest.fixture
def umnik_root(root):
    _ROOT[0] = root
    (root / "vault" / "lab_validation_umnik").mkdir()
    (root / "vault" / "legacy" / "prior_analyses_umnik").mkdir()
    a = _w(root / "vault" / "lab_validation" / "01_XRD" / "a.txt", b"measurement-a")
    jip = _w(root / "vault" / "lab_validation_umnik" / "Mn" / "a.jip", b"BMxx")
    _w(root / "vault" / "lab_validation_umnik" / "Mn" / "a.txt", b"measurement-a")      # copy of archive
    pr = _w(root / "vault" / "legacy" / "prior_analyses_umnik" / "Mn" / "XRD_fit" / "fit.csv", b"derived")
    reg = root / "registry"

    def um(zone, p, sha, role, status, obs="", matched="", dup=""):
        return {"zone": zone, "umnik_path": p, "bytes": 1, "sha256": sha, "role": role, "status": status,
                "matched_registered_paths": matched, "name_variant_of_archive": False,
                "different_version_of": "", "sample_ids": "SMP-001", "observation_ids": obs,
                "byte_dup_group_umnik": dup, "origin": "Mn", "note": ""}
    _csv(reg / "files_manifest.csv", [
        {"zone": "vault/lab_validation", "archive_path": "01_XRD/a.txt", "sha256": a,
         "sample_ids": "SMP-001", "observation_ids": "OBS-XRD-001"}])
    _csv(reg / "files_manifest_umnik.csv", [
        um("vault/lab_validation_umnik", "Mn/a.txt", a, "experimental_raw_or_export",
           "byte_identical_to_archive", "OBS-XRD-001", "01_XRD/a.txt"),
        um("vault/lab_validation_umnik", "Mn/a.jip", jip, "instrument_preview_image", "new", "OBS-XRD-001"),
        um("vault/legacy", "Mn/XRD_fit/fit.csv", pr, "prior_lab_analysis", "new", "PRI-001")])
    obs = [
        {"observation_id": "OBS-XRD-001", "sample_id": "SMP-001", "method": "XRD", "independent": True,
         "files": [{"archive_path": "01_XRD/a.txt", "sha256": a}]},
        {"observation_id": "PRI-001", "sample_id": None, "sample_ids": ["SMP-001"], "role": "prior_lab_analysis",
         "method": "prior_analysis_artifact", "independent": False,
         "files": [{"zone": "vault/legacy", "umnik_path": "x", "sha256": pr}]}]
    (reg / "observations.yaml").write_text(yaml.safe_dump({"observations": obs}), encoding="utf-8")
    (reg / "samples.yaml").write_text(yaml.safe_dump({"samples": [{"sample_id": "SMP-001"}]}), encoding="utf-8")
    pu = [{"observation_id": "OBS-XRD-001", "sample_id": "SMP-001", "method": "XRD",
           "prior_use_status": "used_prior", "validation_mode": "retrospective"}]
    (reg / "prior_use.yaml").write_text(
        yaml.safe_dump({"independent_lab_validation_available": False, "entries": pu}), encoding="utf-8")
    return root


def _run(root):
    rep = audit.Report()
    audit.check_umnik(rep, paths.get_config())
    return {c["name"]: c["ok"] for c in rep.checks}, rep


def test_check_umnik_passes_on_consistent_tree(umnik_root):
    res, rep = _run(umnik_root)
    assert rep.ok, [c for c in rep.checks if not c["ok"]]


def test_check_umnik_flags_cross_source_duplicate_registered_as_new_observation(umnik_root):
    reg = umnik_root / "registry"
    obs = yaml.safe_load((reg / "observations.yaml").read_text(encoding="utf-8"))["observations"]
    obs.insert(1, {"observation_id": "OBS-XRD-009", "sample_id": "SMP-001", "method": "XRD", "independent": True,
                   "files": [dict(obs[0]["files"][0])]})          # same bytes, second 'measurement'
    (reg / "observations.yaml").write_text(yaml.safe_dump({"observations": obs}), encoding="utf-8")
    res, rep = _run(umnik_root)
    assert not res["every lab file sha256 belongs to exactly one lab Observation (no duplicate observations)"]


def test_check_umnik_flags_dup_marked_new_or_linked_to_other_observation(umnik_root):
    reg = umnik_root / "registry"
    rows = list(csv.DictReader(open(reg / "files_manifest_umnik.csv", encoding="utf-8", newline="")))
    rows[0]["status"] = "new"
    _csv(reg / "files_manifest_umnik.csv", rows)
    res, _ = _run(umnik_root)
    assert not res["umnik: files byte-identical to registered archive files map to the same Observation (not new)"]
    rows[0]["status"], rows[0]["observation_ids"] = "byte_identical_to_archive", "OBS-XRD-777"
    _csv(reg / "files_manifest_umnik.csv", rows)
    res, _ = _run(umnik_root)
    assert not res["umnik: files byte-identical to registered archive files map to the same Observation (not new)"]


def test_check_umnik_flags_prior_analysis_outside_legacy(umnik_root):
    # a copy of a prior analysis in the laboratory zone must be reported
    with paths.test_access():
        data = (umnik_root / "vault" / "legacy" / "prior_analyses_umnik" / "Mn" / "XRD_fit" / "fit.csv").read_bytes()
    _w(umnik_root / "vault" / "lab_validation_umnik" / "Mn" / "fit.csv", data)
    reg = umnik_root / "registry"
    rows = list(csv.DictReader(open(reg / "files_manifest_umnik.csv", encoding="utf-8", newline="")))
    rows.append(dict(rows[-1], zone="vault/lab_validation_umnik", umnik_path="Mn/fit.csv",
                     byte_dup_group_umnik="UDUP-01"))
    rows[2]["byte_dup_group_umnik"] = "UDUP-01"
    _csv(reg / "files_manifest_umnik.csv", rows)
    res, _ = _run(umnik_root)
    assert not res["umnik: prior analyses (prior_lab_analysis) exist only in vault/legacy/prior_analyses_umnik (except the declared in-place .mtd projects)"]


def test_check_umnik_flags_missing_prior_use_and_unconfirmed_independence(umnik_root):
    reg = umnik_root / "registry"
    (reg / "prior_use.yaml").write_text(yaml.safe_dump({"entries": []}), encoding="utf-8")
    res, _ = _run(umnik_root)
    assert not res["prior_use: every lab Observation has a record (no extras)"]
    pu = [{"observation_id": "OBS-XRD-001", "sample_id": "SMP-001", "method": "XRD",
           "prior_use_status": "no_prior_use_found", "validation_mode": "independent_candidate"}]
    (reg / "prior_use.yaml").write_text(yaml.safe_dump({"entries": pu}), encoding="utf-8")
    res, _ = _run(umnik_root)
    assert not res["prior_use: statuses/modes consistent (independent_candidate needs no_prior_use_found + user confirmation)"]
    assert not res["prior_use: XRD with prior analyses -> used_prior/retrospective; Mossbauer -> unknown_prior_use until user confirms"]


def test_check_umnik_flags_hash_mismatch_and_extra_file(umnik_root):
    p = umnik_root / "vault" / "lab_validation_umnik" / "Mn" / "a.txt"
    with unlocked(umnik_root):
        p.write_bytes(b"tampered")
    _w(umnik_root / "vault" / "lab_validation_umnik" / "Mn" / "stray.bin", b"?")
    res, _ = _run(umnik_root)
    assert not res["umnik: all 3 files present, sha256 == files_manifest_umnik.csv"]
    assert not res["umnik: no unregistered files in vault/lab_validation_umnik or vault/legacy/prior_analyses_umnik"]


def test_jip_header_parser_reads_geometry_only():
    from hexw.ingest import umnik
    assert umnik.norm_name("SrW_Co3O4_x=1_5_T=1200.mdi") == umnik.norm_name("SrW_Co3O4_x=1.5_T=1200.mdi")
    assert umnik.norm_name("a#1##20260922_142035.txt") == umnik.norm_name("a#1##20260922-142035.txt")


def test_check_umnik_rejects_undeclared_prior_analysis_in_lab_zone(umnik_root):
    """Only the declared in-place .mtd (matched to 03_previous_XRD_projects) may be prior_lab_analysis in the lab zone;
    the positive case is covered by the audit of the real registry."""
    key = "umnik: prior analyses (prior_lab_analysis) exist only in vault/legacy/prior_analyses_umnik (except the declared in-place .mtd projects)"
    reg = umnik_root / "registry"
    rows = list(csv.DictReader(open(reg / "files_manifest_umnik.csv", encoding="utf-8", newline="")))
    rows[0]["role"] = "prior_lab_analysis"                          # a prior analysis in the lab zone, not declared
    _csv(reg / "files_manifest_umnik.csv", rows)
    res, _ = _run(umnik_root)
    assert not res[key]


def test_check_umnik_flags_inconsistent_independence_summary(umnik_root):
    key = "prior_use: summary independent_lab_validation_available is present and consistent with the entries"
    reg = umnik_root / "registry"
    pu = [{"observation_id": "OBS-XRD-001", "sample_id": "SMP-001", "method": "XRD",
           "prior_use_status": "used_prior", "validation_mode": "retrospective"}]
    for summary, ok in ((None, False), (True, False), (False, True)):
        doc = {"entries": pu}
        if summary is not None:
            doc["independent_lab_validation_available"] = summary
        (reg / "prior_use.yaml").write_text(yaml.safe_dump(doc), encoding="utf-8")
        res, _ = _run(umnik_root)
        assert res[key] is ok
