"""Freeze must bind exact forecasts, not just an opaque nonempty marker."""
import json
import pytest
from hexw.io import paths
from hexw.registry.freeze import freeze_run, verify_manifest, FreezeError
from hexw.workflow import atomic_json, core_fingerprint, file_hash


def completed(root):
    from hexw.registry.integrity import write_sums
    (root/"registry"/"metadata.json").write_text('{"synthetic":true}')
    write_sums(root/"registry")
    run = root / "runs" / "predictions" / "synthetic"
    run.mkdir(parents=True)
    for name in ("parameters.npz", "request.json", "summary.json"):
        (run / name).write_text("SYNTHETIC " + name)
    atomic_json(run / "checkpoint.json", {"status": "complete", "run_id": "synthetic",
        "source_fingerprint": core_fingerprint(root),
        "artifacts": {p.name: file_hash(p) for p in run.iterdir() if p.name != "checkpoint.json"}})
    return run


def test_freeze_binds_artifacts_registry_and_seal(root):
    run = completed(root)
    m = freeze_run(run)
    assert paths.is_frozen()
    assert verify_manifest(run / "freeze.json")["model_version"] == m["model_version"]
    assert freeze_run(run) == m
    (run / "summary.json").write_text("changed")
    assert not paths.is_frozen()
    with pytest.raises(FreezeError, match="artifact"):
        verify_manifest(m)


def test_freeze_rejects_new_code_or_fit_input(root):
    run = completed(root)
    freeze_run(run)
    p = root / "project" / "hexw" / "theory" / "new.py"
    p.parent.mkdir(parents=True)
    p.write_text("# changed physical source")
    assert not paths.is_frozen()
    with pytest.raises(FreezeError, match="changed"):
        freeze_run(run)


def test_freeze_refuses_incomplete_and_unsealed_manifest(root):
    run = completed(root)
    cp = json.loads((run / "checkpoint.json").read_text())
    cp["status"] = "running"
    atomic_json(run / "checkpoint.json", cp)
    with pytest.raises(FreezeError, match="completed"):
        freeze_run(run)
    with pytest.raises(FreezeError):
        verify_manifest({"schema": "hexw.freeze.v1", "seal_sha256": "wrong"})


def test_registry_changes_without_updating_checksum_file_are_detected(root):
    run=completed(root);freeze_run(run)
    (root/"registry"/"metadata.json").write_text('{"synthetic":false}')
    assert not paths.is_frozen()
    with pytest.raises(FreezeError,match="registry contents"):
        verify_manifest(run/"freeze.json")
