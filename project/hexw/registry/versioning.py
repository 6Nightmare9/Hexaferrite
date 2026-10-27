"""Model versioning by hash of the parameter set, overlap guards, append-only validation store."""
from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from hexw.registry.schemas import (
    LEAKY_ROLES, DataRef, ParameterEstimate, Prediction, SourceRole, ValidationResult,
    dump_json, to_plain,
)


class LeakageError(ValueError):
    """Validation targets leak into a model/prediction, or forbidden roles are used."""


class AppendOnlyError(RuntimeError):
    pass


def _canon(obj):
    """Canonical numbers for hashing: every int/float becomes a float (so 1 and 1.0 hash alike; json then
    prints repr(float)), -0.0 becomes 0.0. bool and str are left alone."""
    if isinstance(obj, bool):
        return obj
    if isinstance(obj, (int, float)):
        return float(obj) + 0.0
    if isinstance(obj, dict):
        return {k: _canon(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_canon(v) for v in obj]
    return obj


def parameter_set_hash(params: Iterable[ParameterEstimate]) -> str:
    items = sorted((_canon(to_plain(p)) for p in params), key=lambda d: d["name"])
    names = [d["name"] for d in items]
    if len(set(names)) != len(names):
        raise ValueError("duplicate parameter names in set")
    blob = json.dumps(items, sort_keys=True, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class ModelVersion:
    version_id: str
    param_hash: str
    n_params: int
    allows_calibration: bool = False


def make_model_version(params: Iterable[ParameterEstimate], *, allow_explicit_calibration: bool = False,
                       extra_components: dict | None = None) -> ModelVersion:
    """Version id = hash of parameters (+ optional code/corpus/DFT-settings digests).
    Laboratory and legacy roles are rejected; explicit_calibration only with an explicit flag."""
    params = list(params)
    for p in params:
        if p.source_role in LEAKY_ROLES:
            raise LeakageError(f"parameter {p.name!r} has role {p.source_role.value}: "
                               "laboratory / legacy / literature-holdout values cannot define a blind model version")
        if p.source_role == SourceRole.EXPLICIT_CALIBRATION and not allow_explicit_calibration:
            raise LeakageError(f"parameter {p.name!r} is an explicit calibration; "
                               "pass allow_explicit_calibration=True (new version, new validation)")
    h = parameter_set_hash(params)
    if extra_components:
        h = hashlib.sha256((h + json.dumps(extra_components, sort_keys=True)).encode()).hexdigest()
    return ModelVersion("mv-" + h[:16], h, len(params), allow_explicit_calibration)


def _keys(refs: Iterable[DataRef]):
    ids, shas, groups = set(), set(), set()
    for r in refs:
        ids.add(r.id)
        if r.sha256:
            shas.add(r.sha256)
        if r.group:
            groups.add(r.group)
    return ids, shas, groups


def find_overlap(inputs: Iterable[DataRef], targets: Iterable[DataRef]) -> dict:
    """Overlap by id, by file hash, or by semantic group (.mdi vs .txt of one observation)."""
    ia, sa, ga = _keys(inputs)
    ib, sb, gb = _keys(targets)
    # a ref's id may equal another ref's group
    return {
        "ids": sorted((ia & ib) | (ia & gb) | (ga & ib)),
        "sha256": sorted(sa & sb),
        "groups": sorted(ga & gb),
    }


def check_inputs_targets_disjoint(prediction: Prediction, result: ValidationResult) -> None:
    ov = find_overlap(prediction.inputs, result.targets)
    if any(ov.values()):
        raise LeakageError(f"Prediction.inputs and ValidationResult.targets overlap: {ov}")
    if prediction.model_version != result.model_version:
        raise LeakageError("prediction and validation refer to different model versions")


def check_prediction_inputs_roles(prediction: Prediction, allow_explicit_calibration: bool = False) -> None:
    for r in prediction.inputs:
        if r.source_role in LEAKY_ROLES:
            raise LeakageError(f"input {r.id!r} has role {r.source_role.value}")
        if r.source_role == SourceRole.EXPLICIT_CALIBRATION and not allow_explicit_calibration:
            raise LeakageError(f"input {r.id!r} is explicit calibration data")


class ValidationStore:
    """Append-only JSONL with a hash chain. No update/delete API; existing ids are never rewritten."""

    def __init__(self, path):
        self.path = Path(path)

    def _lines(self):
        if not self.path.exists():
            return []
        return [json.loads(l) for l in self.path.read_text(encoding="utf-8").splitlines() if l.strip()]

    def all(self) -> list[dict]:
        return self._lines()

    def append(self, result: ValidationResult, prediction: Prediction | None = None) -> str:
        if prediction is None:
            raise LeakageError("a validation record requires its frozen prediction")
        check_prediction_inputs_roles(prediction)
        check_inputs_targets_disjoint(prediction, result)
        if prediction.prediction_id != result.prediction_id:
            raise LeakageError("prediction id differs from the validation record")
        recs = self._lines()
        if any(r["record"]["result_id"] == result.result_id for r in recs):
            raise AppendOnlyError(f"result_id {result.result_id!r} already exists; results are never overwritten")
        prev = recs[-1]["chain"] if recs else "0" * 64
        body = json.loads(dump_json(result))
        chain = hashlib.sha256((prev + json.dumps(body, sort_keys=True, ensure_ascii=False)).encode()).hexdigest()
        line = json.dumps({"record": body, "prev": prev, "chain": chain}, ensure_ascii=False, sort_keys=True)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.path, "a", encoding="utf-8", newline="\n") as f:
            f.write(line + "\n")
        return chain

    def verify(self) -> bool:
        prev = "0" * 64
        for r in self._lines():
            chain = hashlib.sha256((prev + json.dumps(r["record"], sort_keys=True, ensure_ascii=False)).encode()).hexdigest()
            if r["prev"] != prev or r["chain"] != chain:
                return False
            prev = chain
        return True

    def for_model(self, version_id: str) -> list[dict]:
        return [r["record"] for r in self._lines() if r["record"]["model_version"] == version_id]
