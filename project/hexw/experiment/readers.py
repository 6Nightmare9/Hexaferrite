"""Conservative text readers. Instrument axes are never silently converted into physical units."""
from __future__ import annotations
import io
import re
import numpy as np
from hexw.io import paths


def text(path=None,*,raw=None):
    if raw is None:
        with paths.open_data(path, "rb", purpose="experiment") as f:
            raw = f.read()
    for enc in ("utf-8-sig", "utf-16", "cp1251"):
        try:
            s = raw.decode(enc)
            if "\x00" not in s:
                return s
        except UnicodeError:
            pass
    raise ValueError("unsupported binary format: export a text table from the instrument")


def xy(path=None,*,raw=None):
    """Two numeric columns, arbitrary text header; accepts whitespace/semicolon/comma delimiters."""
    rows = []
    for line in text(path,raw=raw).splitlines():
        line = line.strip()
        if not line or line.startswith(("#", ";")):
            continue
        if ";" in line:
            fields = [x.strip().replace(",", ".") for x in line.split(";")]
        else:
            fields = re.split(r"[\s,]+", line)
        try:
            v = [float(t) for t in fields]
        except ValueError:
            continue
        if len(v) == 2:
            rows.append(v)
    if len(rows) < 8:
        raise ValueError("expected >=8 rows with exactly two numeric columns; axis calibration must be supplied separately")
    a = np.array(rows)
    if not np.isfinite(a).all() or np.any(np.diff(a[:, 0]) <= 0):
        raise ValueError("axis must be finite and strictly increasing; duplicate/unfolded channels are not accepted")
    return a[:, 0], a[:, 1]


def xrd(path,*,raw=None):
    from pathlib import Path
    if Path(path).suffix.lower() != ".mdi":
        return (*xy(path,raw=raw), {})
    lines = text(path,raw=raw).splitlines()
    fields = lines[1].split()
    start, step, end, n = float(fields[0]), float(fields[1]), float(fields[-2]), int(fields[-1])
    vals = np.array([float(v) for line in lines[2:] for v in line.split()])
    # Recovered exports carry -1 followed by zero padding after the n intensities.
    # Only accept a bounded all-zero pad; never truncate arbitrary extra counts.
    trailer=[]
    tail=vals[n:]
    if 3<=len(tail)<=10 and tail[0]==-1 and np.all(tail[1:]==0):
        trailer=tail.tolist(); vals=vals[:n]
    if len(vals) != n or not np.isclose(start + step*(n-1), end, atol=step/2):
        raise ValueError("MDI point count/range mismatch")
    return start+step*np.arange(n), vals, {"wavelength_angstrom": float(fields[-3]), "header": lines[:2],"trailer":trailer}


def spc_channels(path,*,raw=None):
    """Return raw channel data; header Vstart/Vfinish is NOT a calibrated Doppler axis."""
    s = text(path,raw=raw)
    if "Data:" not in s:
        raise ValueError("SPC requires a textual Data: section")
    head, data = s.split("Data:", 1)
    meta = dict(re.findall(r"^\s*([^:\r\n]+):\s*([^\r\n]*)", head, re.M))
    rows = []
    for line in data.splitlines():
        try:
            rows.append([float(x) for x in re.split(r"[\s;]+", line.strip()) if x])
        except ValueError:
            continue
    if rows and all(len(r) == 2 for r in rows):
        a = np.array(rows)
        return a[:, 0], a[:, 1], meta
    values = np.array([v for row in rows for v in row])
    if len(values) < 8:
        raise ValueError("SPC has no usable counts")
    return np.arange(len(values)), values, meta
