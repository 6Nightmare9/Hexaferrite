"""Protection of the curated registry against being overwritten by generator scripts.

Rules (hexw.ingest.build_registry, hexw.ingest.umnik):
  * if any file the script would write already exists in registry/, the script refuses (RegistryOverwriteError);
  * overwriting is possible only with an explicit `--force` AND a clean `git status` of the registry repository
    (every manual edit is therefore already committed and recoverable); a registry that is not a git repository
    counts as not clean;
  * files that were made read-only (key registry files) are re-protected after the write.
Hand-made fields (user_confirmation, independent_lab_validation_available, ...) are never regenerated: see
`merge_prior_use_text` for prior_use.yaml (block-style and empty flow-style `entries: []` are appended to correctly; a non-empty flow list is refused; the merged text is verified before and after writing).
"""
from __future__ import annotations

import os
import stat
from pathlib import Path

from hexw.registry import integrity


class RegistryOverwriteError(PermissionError):
    """A generator tried to overwrite existing registry files without --force or with uncommitted changes."""


def existing_targets(reg, names) -> list[str]:
    reg = Path(reg)
    return [n for n in names if (reg / n).exists()]


def check_write_allowed(reg, names, *, force: bool = False) -> list[str]:
    """Raise unless writing `names` into registry `reg` is allowed. Returns the list of existing files that will
    be overwritten (empty: nothing existed, no force needed)."""
    reg = Path(reg)
    existing = existing_targets(reg, names)
    if not existing:
        return []
    if not force:
        raise RegistryOverwriteError(
            f"registry files already exist and would be overwritten: {existing}. Refusing. Re-run with --force "
            "(only with a clean `git status` of the registry) if you really want to regenerate them.")
    dirty = integrity.git_dirty(reg)
    if dirty is None:
        raise RegistryOverwriteError(
            f"--force refused: {reg} is not a git repository (or git is unavailable), so overwritten manual edits "
            "could not be recovered")
    if dirty:
        raise RegistryOverwriteError(
            f"--force refused: git status of the registry is not clean ({len(dirty)} changed/untracked paths, e.g. "
            f"{dirty[:3]}). Commit or stash manual edits first.")
    return existing


def read_text_exact(path, *, encoding: str = "utf-8") -> str:
    """Text of a file WITHOUT newline translation (CRLF stays CRLF); a UTF-8 BOM stays as U+FEFF in the text."""
    return Path(path).read_bytes().decode(encoding)


def detect_newline(text: str) -> str:
    """Dominant line ending of a text: CRLF only if CRLF is more frequent than bare LF; LF otherwise."""
    crlf = text.count("\r\n")
    return "\r\n" if crlf > text.count("\n") - crlf else "\n"


def write_text_keep_mode(path, text: str, *, encoding: str = "utf-8", newline=None) -> None:
    """Write a file that may be read-only; restore the read-only attribute afterwards. BYTE-exact: nothing is
    translated by the platform.

    newline=None: `text` is LF text; it is written with the line ending of the file that is being replaced (CRLF stays
    CRLF, LF stays LF) and with LF when the file is new. newline="": `text` is written exactly as given (used when the
    text already is the old bytes plus an appended block)."""
    path = Path(path)
    if newline is None:
        style = detect_newline(read_text_exact(path, encoding=encoding)) if path.exists() else "\n"
        text = text.replace("\r\n", "\n")
        if style != "\n":
            text = text.replace("\n", style)
    elif newline != "":
        text = text.replace("\r\n", "\n").replace("\n", newline)
    was_ro = path.exists() and not os.access(path, os.W_OK)
    if was_ro:
        os.chmod(path, os.stat(path).st_mode | stat.S_IWRITE)
    try:
        with open(path, "wb") as fh:
            fh.write(text.encode(encoding))
    finally:
        if was_ro:
            os.chmod(path, os.stat(path).st_mode & ~(stat.S_IWRITE | stat.S_IWGRP | stat.S_IWOTH))


def merge_prior_use_text(existing_text: str, generated_entries: list[dict]) -> tuple[str, dict]:
    """Merge generated prior_use entries into the existing prior_use.yaml text WITHOUT touching manual content.

    Existing entries (keyed by observation_id) are kept byte-for-byte: user_confirmation, hand-set statuses and the
    top-level fields (independent_lab_validation_available/_note, header comments) are never regenerated. Entries
    for observations that are not in the file yet are appended at the end. If `entries:` is not the last top-level
    key the text cannot be extended safely and the merge is refused.
    The text must come from `read_text_exact` (no newline translation); the result starts with `existing_text`
    EXACTLY (BOM, CRLF/LF, no trailing newline are all preserved: only bytes are appended, in the file's own line
    ending style) and must be written with `write_text_keep_mode(..., newline="")`.
    Returns (new_text, report); new_text == existing_text when nothing is added."""
    import yaml
    doc = yaml.safe_load(existing_text) or {}
    have = {e["observation_id"]: e for e in doc.get("entries", [])}
    new = [e for e in generated_entries if e["observation_id"] not in have]
    drift = sorted(e["observation_id"] for e in generated_entries
                   if e["observation_id"] in have and any(
                       have[e["observation_id"]].get(k) != e.get(k) for k in ("prior_use_status", "validation_mode")))
    report = {"kept": len(have), "added": [e["observation_id"] for e in new],
              "status_differs_from_generated_kept_manual": drift,
              "not_generated_but_present": sorted(set(have) - {e["observation_id"] for e in generated_entries})}
    if not new:
        return existing_text, report
    top_keys = list(doc.keys())
    if not top_keys or top_keys[-1] != "entries":
        raise RegistryOverwriteError("prior_use.yaml: `entries` is not the last top-level key; refusing to merge "
                                     "automatically (edit by hand)")
    nl = detect_newline(existing_text)
    style, line_end, head_end, tok_start, tok_end, item_indent = _entries_style(existing_text)
    if style == "flow_nonempty":
        raise RegistryOverwriteError("prior_use.yaml: `entries` is a non-empty flow-style list ([...] / {...}); a block "
                                     "cannot be appended to it safely - refusing (edit by hand or rewrite as a block list)")
    block = yaml.safe_dump(new, allow_unicode=True, sort_keys=False, width=120)
    if item_indent:
        block = "".join((" " * item_indent + ln) if ln.strip() else ln for ln in block.splitlines(True))
    block = block.replace("\n", nl)
    if style == "flow_empty":
        # `entries: []` -> `entries:` + block list. Only the `[]` token (and the blanks before it) change; a hand comment
        # on that line and everything after the line stay byte-for-byte.
        before = existing_text[:tok_start].rstrip(" \t")
        comment = existing_text[tok_end:line_end]
        new_text = before + (("  " + comment.strip()) if comment.strip() else "") + nl + block + existing_text[head_end:]
    else:
        tail = "" if existing_text.endswith("\n") else nl        # the old bytes stay a prefix; only an EOL is appended
        new_text = existing_text + tail + block
    _verify_prior_use_text(new_text, doc.get("entries") or [], [e["observation_id"] for e in new])
    return new_text, report


def _entries_style(text: str):
    """Style of the top-level `entries:` key. Returns (style, end of the line without EOL, end of the line with EOL,
    start and end of the flow token, indent of the existing list items); style is 'block' | 'flow_empty' |
    'flow_nonempty'."""
    import re
    m = re.search(r"(?m)^entries:[ \t]*(?P<rest>[^\r\n]*)", text)
    if m is None:
        raise RegistryOverwriteError("prior_use.yaml: no top-level `entries:` line found")
    line_end = m.end()
    eol = re.compile(r"\r?\n").match(text, line_end)
    head_end = eol.end() if eol else line_end
    rest = m.group("rest")
    body = rest.split("#", 1)[0].strip()
    if body == "":
        im = re.search(r"(?m)^([ \t]*)- ", text[head_end:])
        return "block", line_end, head_end, None, None, (len(im.group(1)) if im else 0)
    tok_start = m.start("rest") + (len(rest) - len(rest.lstrip()))
    if re.fullmatch(r"\[\s*\]", body):
        return "flow_empty", line_end, head_end, tok_start, tok_start + len(body), 0
    return "flow_nonempty", line_end, head_end, None, None, 0


def _verify_prior_use_text(new_text: str, old_entries: list, new_ids: list) -> None:
    """The merged text must parse, keep every existing entry unchanged and contain exactly the new ids appended."""
    import yaml
    try:
        doc = yaml.safe_load(new_text.lstrip("﻿"))
    except yaml.YAMLError as exc:
        raise RegistryOverwriteError(f"prior_use.yaml: the merged text is not valid YAML ({exc}); nothing written") from exc
    got = (doc or {}).get("entries")
    if not isinstance(got, list) or got[:len(old_entries)] != list(old_entries) or \
            [e.get("observation_id") for e in got[len(old_entries):]] != list(new_ids):
        raise RegistryOverwriteError("prior_use.yaml: the merged text does not keep the existing entries and append "
                                     "exactly the new ones; nothing written")


def write_prior_use_checked(path, new_text: str) -> None:
    """Write prior_use.yaml (old bytes + appended block, byte-exact), re-read and parse it; on any problem restore the
    previous bytes (remove a file that did not exist) and raise RegistryOverwriteError."""
    import yaml
    path = Path(path)
    old = path.read_bytes() if path.exists() else None
    write_text_keep_mode(path, new_text, newline="" if old is not None else None)
    try:
        yaml.safe_load(read_text_exact(path).lstrip("﻿"))
    except Exception as exc:                                   # noqa: BLE001 - any failure means: roll back
        if path.exists() and not os.access(path, os.W_OK):
            os.chmod(path, os.stat(path).st_mode | stat.S_IWRITE)
        if old is None:
            path.unlink()
        else:
            path.write_bytes(old)
        raise RegistryOverwriteError(f"{path.name}: unparsable after writing ({exc}); rolled back") from exc
