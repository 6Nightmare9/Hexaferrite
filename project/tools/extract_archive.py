"""Stage-1 extractor (stdlib only). Never executes archive content.

Usage: python -I extract_archive.py <zip> <root>
Writes bytes unchanged into raw/theory_side, vault/lab_validation, vault/legacy.
Records raw stored name bytes into <root>/registry/zip_entries.json.
Nested Archive.zip is NOT extracted (its members are hashed in memory by hexw.ingest).
"""
import hashlib, io, json, os, shutil, sys, zipfile
from pathlib import PurePosixPath

ZONES = {
    "06_literature": "raw/theory_side", "07_software_documentation": "raw/theory_side",
    "metadata": "raw/theory_side",
    "01_XRD": "vault/lab_validation", "02_Mossbauer": "vault/lab_validation",
    "03_previous_XRD_projects": "vault/lab_validation",
    "05_original_duplicate_archive": "vault/lab_validation",
    "04_legacy_model": "vault/legacy",
}
ROOT_FILES_ZONE = "raw/theory_side"


def zone_of(name):
    parts = name.split("/")
    return ZONES.get(parts[0], ROOT_FILES_ZONE if len(parts) == 1 else None)


def raw_name_bytes(i):
    # zipfile decodes with utf-8 if flag 0x800 else cp437; invert to get stored bytes
    return i.orig_filename.encode("utf-8" if i.flag_bits & 0x800 else "cp437")


def main(zpath, root):
    zsha = hashlib.sha256(open(zpath, "rb").read()).hexdigest()
    os.makedirs(os.path.join(root, "raw", "original_zip"), exist_ok=True)
    dst = os.path.join(root, "raw", "original_zip", os.path.basename(zpath))
    if not os.path.exists(dst):
        shutil.copyfile(zpath, dst)
    assert hashlib.sha256(open(dst, "rb").read()).hexdigest() == zsha
    z = zipfile.ZipFile(zpath)
    bad = z.testzip()
    assert bad is None, bad
    entries = []
    for i in z.infolist():
        if i.is_dir():
            continue
        name = i.filename
        p = PurePosixPath(name)
        if p.is_absolute() or ".." in p.parts or ":" in name:
            raise SystemExit("unsafe path in zip: %r" % name)
        zone = zone_of(name)
        if zone is None:
            raise SystemExit("unmapped top-level: %r" % name)
        target = os.path.join(root, *zone.split("/"), *p.parts)
        os.makedirs(os.path.dirname(target), exist_ok=True)
        data = z.read(i)
        if os.path.exists(target):
            raise SystemExit("exists: %s" % target)
        with open(target, "wb") as f:
            f.write(data)
        raw = raw_name_bytes(i)
        entries.append({
            "archive_path": name, "zone": zone, "raw_name_hex": raw.hex(),
            "utf8_flag": bool(i.flag_bits & 0x800),
            "name_is_ascii": all(c < 128 for c in raw),
            "size": i.file_size, "crc32": "%08x" % i.CRC,
            "sha256": hashlib.sha256(data).hexdigest(),
            "date_time": list(i.date_time),
        })
    out = {"zip_sha256": zsha, "zip_bytes": os.path.getsize(zpath), "n_files": len(entries),
           "entries": entries}
    json.dump(out, open(os.path.join(root, "registry", "zip_entries.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("extracted", len(entries), "zip sha256", zsha)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
