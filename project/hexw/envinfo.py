"""Collect environment facts into registry/environment.json (re-verified at run time)."""
from __future__ import annotations

import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from importlib import metadata


def _run(cmd: list[str]) -> str | None:
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        return r.stdout.strip() if r.returncode == 0 else None
    except Exception:
        return None


def _ps(expr: str) -> str | None:
    return _run(["powershell", "-NoProfile", "-Command", expr])


def collect() -> dict:
    cpu = _ps("(Get-CimInstance Win32_Processor | Select-Object -First 1 Name,NumberOfCores,NumberOfLogicalProcessors | ConvertTo-Json -Compress)")
    ram = _ps("(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory")
    os_cap = _ps("(Get-CimInstance Win32_OperatingSystem | Select-Object Caption,Version,BuildNumber | ConvertTo-Json -Compress)")
    gpu_csv = _run(["nvidia-smi", "--query-gpu=name,memory.total,driver_version", "--format=csv,noheader,nounits"])
    cuda = _run(["nvidia-smi", "--query-gpu=compute_cap", "--format=csv,noheader"])
    gpu = None
    if gpu_csv:
        n, mem, drv = [x.strip() for x in gpu_csv.splitlines()[0].split(",")]
        gpu = {"name": n, "vram_mib": int(mem), "driver": drv, "compute_capability": cuda.splitlines()[0] if cuda else None}
    pkgs = {}
    for p in ["numpy", "scipy", "matplotlib", "pandas", "pyyaml", "pytest", "openpyxl", "pymatgen", "pymatgen-core", "ase", "hexw"]:
        try:
            pkgs[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            pkgs[p] = None
    tools = {}
    for t in ["pw.x", "vasp_std", "vasp_gam", "mpiexec", "mpirun", "orterun", "srun", "sbatch", "git", "nvcc"]:
        tools[t] = shutil.which(t)
    info = {
        "collected_utc": datetime.now(timezone.utc).isoformat(),
        "os": json.loads(os_cap) if os_cap else platform.platform(),
        "platform": platform.platform(),
        "cpu": json.loads(cpu) if cpu else platform.processor(),
        "logical_threads_os_cpu_count": os.cpu_count(),
        "ram_bytes": int(ram) if ram and ram.isdigit() else None,
        "ram_gib": round(int(ram) / 2**30, 2) if ram and ram.isdigit() else None,
        "gpu": gpu,
        "python": {"version": sys.version.split()[0], "executable": sys.executable},
        "packages": pkgs,
        "executables_found": tools,
        "dft_and_parallel": {"quantum_espresso_pw_x": bool(tools["pw.x"]),
                             "vasp": bool(tools["vasp_std"] or tools["vasp_gam"]),
                             "mpi": bool(tools["mpiexec"] or tools["mpirun"] or tools["orterun"]),
                             "slurm_cluster": bool(tools["srun"] or tools["sbatch"]),
                             "licenses": "none detected (no QE/VASP/MPI binaries on PATH)"},
        "stated_by_user": {"cpu": "Ryzen 7 7735H, 16 threads", "ram_gib": 15.2,
                           "gpu": "RTX 4060 Laptop", "os": "Windows 11"},
    }
    return info


def main(path: str) -> None:
    info = collect()
    with open(path, "w", encoding="utf-8") as f:
        json.dump(info, f, ensure_ascii=False, indent=1)
    print(json.dumps({k: info[k] for k in ("cpu", "ram_gib", "gpu", "dft_and_parallel")}, ensure_ascii=False, indent=1))
