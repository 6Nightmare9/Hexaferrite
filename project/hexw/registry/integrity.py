"""Registry integrity: SHA256SUMS, git state (HEAD, clean tree). Reads only the registry, never the vault.

SHA256SUMS covers every tracked-style registry file except itself, `.git/`, `runs/` (audit manifests, git-ignored)
and `.gitignore`; format is the usual `<hex>  <relative/path>` with LF line ends.
"""
from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path

SUMS = "SHA256SUMS"
EXCLUDE_TOP = {".git", "runs", SUMS, ".gitignore", ".gitattributes", "FROZEN_MODEL.json"}
KEY_READONLY = ("samples.yaml", "observations.yaml", "prior_use.yaml")


def _sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def registry_files(reg) -> list[Path]:
    reg = Path(reg)
    return sorted((p for p in reg.rglob("*") if p.is_file()
                   and p.relative_to(reg).parts[0] not in EXCLUDE_TOP), key=lambda p: p.relative_to(reg).as_posix())


def compute_sums(reg) -> dict[str, str]:
    reg = Path(reg)
    return {p.relative_to(reg).as_posix(): _sha(p) for p in registry_files(reg)}


def write_sums(reg) -> dict[str, str]:
    reg = Path(reg)
    sums = compute_sums(reg)
    (reg / SUMS).write_text("".join(f"{h}  {n}\n" for n, h in sums.items()), encoding="utf-8", newline="\n")
    return sums


def read_sums(reg) -> dict[str, str]:
    out = {}
    for line in (Path(reg) / SUMS).read_text(encoding="utf-8").splitlines():
        if line.strip():
            h, n = line.split("  ", 1)
            out[n] = h
    return out


def verify_sums(reg) -> list[str]:
    """Problems (empty list = SHA256SUMS matches the files exactly)."""
    try:
        want = read_sums(reg)
    except OSError as e:
        return [f"{SUMS} unreadable: {e}"]
    have = compute_sums(reg)
    return ([f"missing file {n}" for n in sorted(set(want) - set(have))]
            + [f"unlisted file {n}" for n in sorted(set(have) - set(want))]
            + [f"hash mismatch {n}" for n in sorted(set(want) & set(have)) if want[n] != have[n]])


def _git(repo, *args) -> str | None:
    try:
        r = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    return r.stdout.strip() if r.returncode == 0 else None


def git_head(repo) -> str | None:
    return _git(repo, "rev-parse", "HEAD")


def git_dirty(repo) -> list[str] | None:
    """Lines of `git status --porcelain` (untracked included, ignored not); None if not a git repo."""
    out = _git(repo, "status", "--porcelain")
    return None if out is None else [l for l in out.splitlines() if l.strip()]
