"""Small orchestration primitives. No laboratory reads and no process launch on import."""
from __future__ import annotations

import dataclasses
import hashlib
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


def plain(obj):
    if dataclasses.is_dataclass(obj):
        return {f.name: plain(getattr(obj, f.name)) for f in dataclasses.fields(obj)}
    if isinstance(obj, dict):
        return {str(k): plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [plain(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return plain(obj.tolist())
    if isinstance(obj, (np.integer, np.floating, np.bool_)):
        return plain(obj.item())
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    if isinstance(obj, Path):
        return str(obj)
    if hasattr(obj, "value"):
        return plain(obj.value)
    return obj


def canonical(obj):
    return json.dumps(plain(obj), ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False)


def digest(obj):
    return hashlib.sha256(canonical(obj).encode("utf-8")).hexdigest()


def file_hash(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def atomic_json(path, obj):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(p.name + f".{os.getpid()}.tmp")
    with open(tmp, "w", encoding="utf-8", newline="\n") as f:
        json.dump(plain(obj), f, ensure_ascii=False, indent=2, allow_nan=False)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, p)


def utc():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def summary(values):
    a = np.asarray(values, float).ravel()
    valid = a[np.isfinite(a)]
    if not len(valid):
        return {"median": None, "lo95": None, "hi95": None, "n_valid": 0, "n_total": len(a)}
    q = np.quantile(valid, [0.025, 0.5, 0.975])
    return {"median": float(q[1]), "lo95": float(q[0]), "hi95": float(q[2]),
            "n_valid": len(valid), "n_total": len(a)}


def core_fingerprint(root):
    """Prediction code, numerical configuration and FIT data; never validation targets.

    Experimental analysis/UI may evolve without retrospectively changing frozen theoretical predictions.
    Their own source hashes are carried by analysis results instead.
    """
    root = Path(root)
    pkg = root / "project" / "hexw"
    files = [p for p in pkg.rglob("*.py") if p.relative_to(pkg).parts[0] not in {"experiment", "validate", "ui"}]
    files += list((root / "project" / "config" / "dft").glob("*.yaml"))
    files += [root/"project"/name for name in ("pyproject.toml","requirements.txt") if (root/"project"/name).is_file()]
    policy=root/"project/config/holdout_sources.json"
    if policy.is_file():
        files.append(policy)
    fit = root / "literature" / "fit"
    files += [p for p in fit.rglob("*") if p.is_file() and p.suffix != ".pyc" and ".git" not in p.parts]
    from hexw.structures.configs import CIF_REL
    template = root / CIF_REL
    if template.is_file():
        files.append(template)
    result={p.relative_to(root).as_posix(): file_hash(p) for p in sorted(files)}
    import sys
    import importlib.metadata
    result["runtime/python"]=sys.version
    for package in ("numpy","scipy","pymatgen","pymatgen-core","spglib","matplotlib","pyyaml"):
        try:
            result["runtime/"+package]=importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            result["runtime/"+package]="not-installed"
    return result


def friendly_resources():
    """Applied before numerical imports by entry points; never changes system-wide power settings."""
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
    os.environ.setdefault("MKL_NUM_THREADS", "2")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "2")
    if os.name == "nt":
        import ctypes
        kernel=ctypes.WinDLL("kernel32",use_last_error=True)
        kernel.GetCurrentProcess.restype=ctypes.c_void_p
        kernel.SetPriorityClass.argtypes=[ctypes.c_void_p,ctypes.c_uint]
        kernel.SetPriorityClass.restype=ctypes.c_int
        if not kernel.SetPriorityClass(kernel.GetCurrentProcess(),0x4000):
            raise OSError(ctypes.get_last_error(),"cannot set background calculation priority")
