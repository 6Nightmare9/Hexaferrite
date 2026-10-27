"""Copy the Umnik folder byte-for-byte into the vault (one-off, reproducible).

  python tools/copy_umnik.py <Umnik source dir> <vault dir>

Layout: raw measurements -> <vault>/lab_validation_umnik/<same relative path>;
prior analyses (Co/XRD_fit, Zn/XRD_fit, Co/Claude outputs) -> <vault>/legacy/prior_analyses_umnik/<same path>.
Verifies sha256(copy) == sha256(source) for every file, sets files read-only, never modifies the source
(source hashes are compared before/after) and refuses to overwrite existing files.
"""
import hashlib
import os
import shutil
import stat
import sys
from pathlib import Path


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def is_prior(rel: Path) -> bool:
    p = rel.parts
    return len(p) > 2 and ((p[0] in ("Co", "Zn") and p[1] == "XRD_fit") or (p[0] == "Co" and p[1] == "Claude outputs"))


def main(src, vault):
    src, vault = Path(src), Path(vault)
    lab, pri = vault / "lab_validation_umnik", vault / "legacy" / "prior_analyses_umnik"
    files = sorted(p for p in src.rglob("*") if p.is_file())
    before = {p: (sha(p), p.stat().st_size, p.stat().st_mtime_ns) for p in files}
    for p in files:
        rel = p.relative_to(src)
        dst = (pri if is_prior(rel) else lab) / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.exists():
            raise SystemExit(f"refusing to overwrite {dst}")
        shutil.copyfile(p, dst)
        shutil.copystat(p, dst)
        if sha(dst) != before[p][0]:
            raise SystemExit(f"hash mismatch after copy: {rel}")
        os.chmod(dst, stat.S_IREAD)
    after = {p: (sha(p), p.stat().st_size, p.stat().st_mtime_ns) for p in src.rglob("*") if p.is_file()}
    print(f"copied {len(files)} files; source unchanged: {before == after}")
    return 0 if before == after else 1


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:3]))
