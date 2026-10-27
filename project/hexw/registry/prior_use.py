"""Prior-use disclosure for laboratory observations (registry/prior_use.yaml).

A laboratory Observation that was already analysed (prior XRD fits, .mtd projects, the legacy model 04)
cannot serve as a blind test. Its status travels into ValidationResult.validation_mode.

Rules:
  * prior_use_status == used_prior        -> validation_mode retrospective
  * prior_use_status == unknown_prior_use -> validation_mode retrospective (conservative) until the user confirms
  * independent_candidate only with prior_use_status == no_prior_use_found AND a user_confirmation record.
"""
from __future__ import annotations

from pathlib import Path

import yaml

from hexw.registry.schemas import PRIOR_USE_STATUSES, ValidationMode, ValidationResult

ROLE_PRIOR_LAB_ANALYSIS = "prior_lab_analysis"
MODES = tuple(m.value for m in ValidationMode)
LAB_METHODS = ("XRD", "Mossbauer")


def load_prior_use(registry_dir) -> dict[str, dict]:
    doc = yaml.safe_load((Path(registry_dir) / "prior_use.yaml").read_text(encoding="utf-8"))
    return {e["observation_id"]: e for e in doc["entries"]}


def independent_validation_available(registry_dir) -> bool | None:
    """Summary field of prior_use.yaml (None if absent)."""
    doc = yaml.safe_load((Path(registry_dir) / "prior_use.yaml").read_text(encoding="utf-8"))
    return doc.get("independent_lab_validation_available")


def problems(entry: dict) -> list[str]:
    """Consistency problems of one prior_use entry (empty list = consistent)."""
    out = []
    st, mode = entry.get("prior_use_status"), entry.get("validation_mode")
    if st not in PRIOR_USE_STATUSES:
        out.append(f"bad prior_use_status {st!r}")
    if mode not in MODES:
        out.append(f"bad validation_mode {mode!r}")
    if mode == ValidationMode.INDEPENDENT_CANDIDATE.value:
        if st != "no_prior_use_found":
            out.append("independent_candidate requires no_prior_use_found")
        if not entry.get("user_confirmation"):
            out.append("independent_candidate requires user_confirmation")
    if st in ("used_prior", "unknown_prior_use") and mode != ValidationMode.RETROSPECTIVE.value:
        out.append(f"{st} must give validation_mode retrospective")
    return out


def mode_for_observations(entries: dict[str, dict], obs_ids) -> tuple[ValidationMode, str]:
    """Combined (validation_mode, prior_use_status) for a set of target observations.
    Independent only if EVERY target is; a missing record counts as unknown (never as independent)."""
    ids = list(obs_ids)
    if not ids:
        return ValidationMode.RETROSPECTIVE, "unknown_prior_use"
    sts = []
    for i in ids:
        e = entries.get(i)
        sts.append(e["prior_use_status"] if e and not problems(e) else "unknown_prior_use")
    if all(s == "no_prior_use_found" for s in sts):
        if all(entries[i]["validation_mode"] == ValidationMode.INDEPENDENT_CANDIDATE.value for i in ids):
            return ValidationMode.INDEPENDENT_CANDIDATE, "no_prior_use_found"
        return ValidationMode.RETROSPECTIVE, "no_prior_use_found"
    status = "used_prior" if "used_prior" in sts else "unknown_prior_use"
    return ValidationMode.RETROSPECTIVE, status


def apply_prior_use(result: ValidationResult, entries: dict[str, dict]) -> ValidationResult:
    """Set result.validation_mode / prior_use_status from the targets' Observation ids (DataRef.group or .id)."""
    obs_ids = [t.group or t.id for t in result.targets]
    result.validation_mode, result.prior_use_status = mode_for_observations(entries, obs_ids)
    return result
