"""A frozen forecast is bound to code, fit corpus, numerical settings and exact artifacts."""
from __future__ import annotations

import json
from pathlib import Path

from hexw.io import paths
from hexw.workflow import atomic_json, core_fingerprint, digest, file_hash, utc

SCHEMA = "hexw.freeze.v1"


def _registry_valid(root):
    from hexw.registry.integrity import verify_sums
    errors=verify_sums(root/"registry")
    if errors:
        raise FreezeError("registry contents do not match SHA256SUMS: "+str(errors[:3]))


class FreezeError(ValueError):
    pass


def _inside(root, relative):
    p = (root / relative).resolve()
    if not p.is_relative_to(root.resolve()):
        raise FreezeError("manifest path escapes project")
    return p


def verify_manifest(manifest, root=None, *, check_sources=True):
    root = Path(root or paths.get_config().root)
    m = json.loads(Path(manifest).read_text(encoding="utf-8")) if not isinstance(manifest, dict) else manifest
    if m.get("schema") != SCHEMA:
        raise FreezeError("unsupported or legacy freeze marker")
    body = {k: v for k, v in m.items() if k != "seal_sha256"}
    if digest(body) != m.get("seal_sha256"):
        raise FreezeError("freeze manifest seal does not match")
    run = _inside(root, m["run_relative"])
    if not run.is_relative_to((root / "runs").resolve()):
        raise FreezeError("freeze does not reference a prediction run")
    if check_sources and core_fingerprint(root) != m["source_fingerprint"]:
        raise FreezeError("prediction code or fit corpus changed since freeze")
    if file_hash(root / "registry" / "SHA256SUMS") != m["registry_sums_sha256"]:
        raise FreezeError("sample registry changed since freeze")
    _registry_valid(root)
    for name, h in m["artifacts"].items():
        p = (run / name).resolve()
        if not p.is_relative_to(run) or file_hash(p) != h:
            raise FreezeError(f"frozen artifact changed: {name}")
    if file_hash(run / "parameters.npz") != m["parameters_sha256"]:
        raise FreezeError("frozen parameters changed")
    if file_hash(run / "request.json") != m["request_sha256"]:
        raise FreezeError("frozen inputs changed")
    return m


def freeze_run(run, *, activate=True):
    root = paths.get_config().root
    run = Path(run).resolve()
    if not run.is_relative_to((root / "runs").resolve()):
        raise FreezeError("only project prediction runs can be frozen")
    state = json.loads((run / "checkpoint.json").read_text(encoding="utf-8"))
    if state.get("status") != "complete":
        raise FreezeError("only a completed forecast can be frozen")
    if core_fingerprint(root) != state["source_fingerprint"]:
        raise FreezeError("code/corpus changed during or after this run; reproduce a new version before freezing")
    _registry_valid(root)
    for name, h in state["artifacts"].items():
        if file_hash(run / name) != h:
            raise FreezeError(f"run artifact changed: {name}")
    existing = run / "freeze.json"
    if existing.exists():
        m = verify_manifest(existing, root)
    else:
        m = {"schema": SCHEMA, "run_id": state["run_id"], "created_utc": utc(),
             "run_relative": run.relative_to(root).as_posix(), "source_fingerprint": state["source_fingerprint"],
             "artifacts": state["artifacts"], "parameters_sha256": file_hash(run / "parameters.npz"),
             "request_sha256": file_hash(run / "request.json"),
             "registry_sums_sha256": file_hash(root / "registry" / "SHA256SUMS"),
             "validation_status": "predictions fixed; no experimental accuracy claim"}
        m["model_version"] = "mv-" + digest({k: v for k, v in m.items() if k != "created_utc"})[:20]
        m["seal_sha256"] = digest(m)
        atomic_json(existing, m)
    if activate:
        atomic_json(root / "registry" / "FROZEN_MODEL.json", m)
    return m
