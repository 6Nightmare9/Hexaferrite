"""Build registry/files_manifest.csv, samples.yaml, observations.yaml (preliminary, stage 1).

Uses ONLY: file names, header text lines (.spc up to 'Data:', .mdi first two lines), file sizes and
byte hashes. Never parses measurement channels / intensities / prediction tables.
"""
from __future__ import annotations

import csv
import io
import json
import re
import zipfile
from collections import defaultdict
from pathlib import Path, PurePosixPath

import yaml

from hexw.ingest import registry_guard as guard
from hexw.io import paths
from hexw.registry.schemas import PRELIMINARY, Observation, SampleRecord, to_plain

NESTED = "05_original_duplicate_archive/Archive.zip"
OUTPUTS = ("files_manifest.csv", "samples.yaml", "observations.yaml", "nested_archive_members.json")


def pf(value, unit=None, note=None):
    """Preliminary value derived from a file name."""
    d = {"value": value, "unit": unit, "status": PRELIMINARY}
    if note:
        d["note"] = note
    return d


def _num(s: str):
    return float(s.replace("_", ".").replace(",", "."))


# ---------------- header-only readers ----------------
def read_spc_header(path) -> dict:
    head = paths.read_prefix(path, 600, purpose="metadata")
    i = head.find(b"Data:")
    if i < 0:
        return {}
    text = head[:i].decode("latin-1")
    out = {}
    for line in text.replace("\r", "").split("\n"):
        if not line.strip():
            continue
        if line.startswith(("Date:", "Comment:")):
            k, _, v = line.partition(":")
            out[k] = v.strip()
            continue
        # 'Voltage: 608 Gain: 10 LBound: 23 HBound: 57' holds several fields
        for k, v in re.findall(r"([A-Za-z]+):\s*([^:]*?)(?=\s+[A-Za-z]+:|$)", line):
            out[k] = v.strip()
    return out


def read_mdi_header(path) -> dict:
    head = paths.read_prefix(path, 220, purpose="metadata")
    lines = head.decode("latin-1").replace("\r", "").split("\n")[:2]
    l1, l2 = lines[0], lines[1]
    m1 = re.match(r"\s*(\d\d)/(\d\d)/(\d\d)\s+(.*)$", l1)
    t = l2.split()
    out = {"line1_raw": l1.strip(), "line2_raw": l2.strip()}
    if m1:
        out["date_dd_mm_yy"] = f"{m1.group(1)}/{m1.group(2)}/{m1.group(3)}"
        rest = m1.group(4)
        out["setup_raw"] = rest.strip()
        mk = re.search(r"(\d+)KV/(\d+)mA", rest)
        if mk:
            out["tube_kV"], out["tube_mA"] = int(mk.group(1)), int(mk.group(2))
        ms = re.search(r"Slit:(\S+)", rest)
        if ms:
            out["slit"] = ms.group(1)
        out["monochromator"] = "ON" if "Monochromator: ON" in rest else None
    if len(t) >= 7:
        out.update({"two_theta_start_deg": float(t[0]), "step_deg": float(t[1]),
                    "anode": t[3].upper(), "wavelength_angstrom": float(t[4]),
                    "two_theta_end_deg": float(t[5]), "n_points": int(t[6])})
    return out


def count_lines(path) -> int:
    n = 0
    with paths.open_data(path, "rb", purpose="metadata") as f:
        while True:
            b = f.read(1 << 20)
            if not b:
                break
            n += b.count(b"\n")
    return n


# ---------------- name parsers ----------------
RE_SRMN = re.compile(r"^SrMn_(?P<prec>[A-Za-z0-9]+)_x=(?P<x>[\d.]+)_(?P<T>\d+)_(?P<d1>\d+h)_(?P<d2>\d+h)")
RE_SRW_XRD = re.compile(
    r"^SrW_(?P<prec>Co3O4|ZnO)_[xX]=(?P<x>[\d._]+)_T=(?P<T1>\d+)_(?P<t1>\d+)h_T=(?P<T2>\d+)_(?P<t2>\d+)h_(?P<date>\d{6})")
RE_MOS = re.compile(
    r"^(?P<mm>\d\d)\.(?P<yy>\d\d)_SrW_(?P<prec>Mn2O3|MnO)_x=(?P<x>[\d.,]+)_(?P<T>\d+)°C(?:_(?P<t>\d+)h)?_"
    r"(?P<dd>\d\d)\.(?P<m2>\d\d)\.(?P<y2>\d\d)_(?P<mass>\d+)mg\.spc$")


def xrd_group_key(name: str) -> str:
    stem = name.rsplit(".", 1)[0]
    stem = re.sub(r"#1##\d{8}[-_]\d{6}_100$", "", stem)
    stem = re.sub(r"(x=\d)_(\d)", r"\1.\2", stem)        # x=1_5 -> x=1.5 (display normalisation)
    stem = re.sub(r"_[xX]=", "_x=", stem)
    return stem


def sample_from_xrd(group_key: str) -> dict:
    m = RE_SRMN.match(group_key) or RE_SRW_XRD.match(group_key)
    if not m:
        raise ValueError(group_key)
    g = m.groupdict()
    x = _num(g["x"].replace("_", "."))
    out = {"x": x, "precursor": g["prec"], "series": group_key.split("_")[0]}
    if "T" in g and g.get("T"):                        # SrMn series
        out["synthesis"] = {
            "raw_tokens": f"{g['T']}_{g['d1']}_{g['d2']}",
            "stages": [{"T_C": pf(int(g["T"]), "C", "unit C not in name, assumed"),
                        "duration_raw": g["d1"] + "_" + g["d2"]}],
            "note": "duration tokens '0h_6h' are ambiguous; not interpreted"}
        out["flags"] = ["srmn_regime_tokens_ambiguous", "precursor_Mn_species_unspecified"]
    else:
        out["synthesis"] = {
            "raw_tokens": f"T={g['T1']}_{g['t1']}h_T={g['T2']}_{g['t2']}h",
            "stages": [{"T_C": pf(int(g["T1"]), "C"), "duration": pf(int(g["t1"]), "h")},
                       {"T_C": pf(int(g["T2"]), "C"), "duration": pf(int(g["t2"]), "h")}],
            "atmosphere": None, "cooling": None,
            "name_date_token_ddmmyy": pf(g["date"], None, "meaning of this date (synthesis/project) unknown")}
        out["flags"] = []
    return out


# ---------------- main builder ----------------
def _csv_rows(path):
    with open(path, encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def build(root=None, write=True, force=False) -> dict:
    """Preliminary registry from the zip entries. Refuses to overwrite existing registry files (see
    hexw.ingest.registry_guard): needs force=True AND a clean `git status` of the registry."""
    cfg = paths.get_config() if root is None else paths.configure(root)
    reg = cfg.registry
    if write:
        guard.check_write_allowed(reg, OUTPUTS, force=force)
    ze = json.loads((reg / "zip_entries.json").read_text(encoding="utf-8"))
    entries = ze["entries"]
    all_m = {r["path"]: r for r in _csv_rows(cfg.raw_theory / "ALL_FILES_MANIFEST.csv")}
    data_m = {r["archive_path"]: r for r in _csv_rows(cfg.raw_theory / "metadata" / "DATA_MANIFEST.csv")}

    def disk(e):
        return cfg.root / Path(*e["zone"].split("/")) / Path(*e["archive_path"].split("/"))

    rows, by_sha = [], defaultdict(list)
    for e in entries:
        p = disk(e)
        sha = paths.sha256_file(p)
        am, dm = all_m.get(e["archive_path"]), data_m.get(e["archive_path"])
        base = PurePosixPath(e["archive_path"]).name
        row = {
            "zone": e["zone"], "archive_path": e["archive_path"], "extracted": True,
            "raw_name_hex": e["raw_name_hex"], "utf8_flag": e["utf8_flag"],
            "bytes": p.stat().st_size, "sha256": sha,
            "sha256_matches_zip_member": sha == e["sha256"],
            "all_manifest_sha256": am["sha256"] if am else "",
            "all_manifest_match": (am["sha256"].lower() == sha) if am else "not_listed",
            "all_manifest_bytes_match": (int(am["bytes"]) == p.stat().st_size) if am else "not_listed",
            "data_manifest_sha256": dm["sha256"] if dm else "",
            "data_manifest_match": (dm["sha256"].lower() == sha) if dm else "not_listed",
            "data_manifest_original_name_match": (dm["original_name"] == base) if dm else "not_listed",
            "data_manifest_role": dm["role"] if dm else "",
            "byte_dup_group": "", "semantic_group": "", "sample_ids": "", "observation_ids": "", "note": "",
        }
        rows.append(row)
        by_sha[sha].append(row)

    # nested Archive.zip: members hashed in memory, NOT extracted
    nested_rows, nested_info = [], []
    with paths.open_data(cfg.vault_lab / NESTED, "rb", purpose="integrity") as f:
        zb = f.read()
    with zipfile.ZipFile(io.BytesIO(zb)) as nz:
        for i in nz.infolist():
            if i.is_dir():
                continue
            raw = i.orig_filename.encode("utf-8" if i.flag_bits & 0x800 else "cp437")
            sha = __import__("hashlib").sha256(nz.read(i)).hexdigest()
            nested_info.append({"member": i.filename, "utf8_flag": bool(i.flag_bits & 0x800),
                                "raw_hex": raw.hex(), "bytes": i.file_size, "sha256": sha})
            r = {"zone": "vault/lab_validation", "archive_path": f"{NESTED}::{i.filename}", "extracted": False,
                 "raw_name_hex": raw.hex(), "utf8_flag": bool(i.flag_bits & 0x800), "bytes": i.file_size,
                 "sha256": sha, "sha256_matches_zip_member": True, "all_manifest_sha256": "",
                 "all_manifest_match": "not_listed", "all_manifest_bytes_match": "not_listed",
                 "data_manifest_sha256": "", "data_manifest_match": "not_listed",
                 "data_manifest_original_name_match": "not_listed", "data_manifest_role": "nested_member",
                 "byte_dup_group": "", "semantic_group": "", "sample_ids": "", "observation_ids": "",
                 "note": "hashed in memory, not extracted"}
            nested_rows.append(r)
            by_sha[sha].append(r)
    rows += nested_rows

    # nested names vs DATA_MANIFEST duplicates_inside_Archive_zip
    dup_decl = {}
    for dm in data_m.values():
        if dm["duplicates_inside_Archive_zip"]:
            dup_decl[dm["archive_path"]] = dm["duplicates_inside_Archive_zip"]

    gid = 0
    for sha, rs in sorted(by_sha.items()):
        if len(rs) > 1:
            gid += 1
            for r in rs:
                r["byte_dup_group"] = f"DUP-{gid:02d}"
    row_by_path = {r["archive_path"]: r for r in rows}

    # ---- samples & observations ----
    samples, observations = [], []
    xrd_files = sorted((r for r in rows if r["archive_path"].startswith("01_XRD/")), key=lambda r: r["archive_path"])
    groups = defaultdict(list)
    for r in xrd_files:
        groups[xrd_group_key(PurePosixPath(r["archive_path"]).name)].append(r)
    n_smp = n_obs = 0
    mtd = [r for r in rows if r["archive_path"].startswith("03_previous_XRD_projects/")]
    for key in sorted(groups):
        rs = groups[key]
        info = sample_from_xrd(key)
        n_smp += 1
        n_obs += 1
        sid, bid, oid = f"SMP-{n_smp:03d}", f"BAT-{n_smp:03d}", f"OBS-XRD-{n_obs:03d}"
        flags = list(info["flags"])
        files, hdr = [], {}
        mdi = [r for r in rs if r["archive_path"].endswith(".mdi")]
        txt = [r for r in rs if r["archive_path"].endswith(".txt")]
        for r in rs:
            rep = "instrument_native_mdi" if r["archive_path"].endswith(".mdi") else "text_export_of_same_measurement"
            files.append({"archive_path": r["archive_path"], "zone": r["zone"], "sha256": r["sha256"],
                          "representation": rep})
            r["semantic_group"], r["sample_ids"], r["observation_ids"] = oid, sid, oid
        if mdi:
            hdr = read_mdi_header(cfg.root / "vault" / "lab_validation" / Path(*mdi[0]["archive_path"].split("/")))
        if txt:
            nl = count_lines(cfg.root / "vault" / "lab_validation" / Path(*txt[0]["archive_path"].split("/")))
            hdr["txt_export_line_count"] = nl
            if hdr.get("n_points") and nl not in (hdr["n_points"], hdr["n_points"] + 1):
                flags.append("txt_line_count_differs_from_mdi_n_points")
        if len(mdi) != 1 or len(txt) != 1:
            flags.append("incomplete_mdi_txt_pair")
        names = {PurePosixPath(r["archive_path"]).name.rsplit(".", 1)[0] for r in rs}
        pair_note = None
        if key.startswith("SrW_Co3O4_x=1.5"):
            flags.append("pair_names_differ_requires_confirmation")
            pair_note = ("mdi is named x=1.5 without timestamp, txt is named x=1_5 with #1## timestamp "
                         "(2026-09-22 14:20:35); same measurement is likely (header date 22/09/26, same stem) "
                         "but UNCONFIRMED; x=1.5 vs 1_5 may also be a typo for another sample")
        if "ZnO" in key and key.endswith("160926") and "x=1.5" in key:
            flags.append("previous_xrd_projects_exist_not_ground_truth")
        precursor = info["precursor"]
        A = {"Mn": "Mn", "Co3O4": "Co", "ZnO": "Zn"}[precursor]
        x = info["x"]
        fl = info["flags"]
        smp = SampleRecord(
            sample_id=sid, batch_id=bid, label=key,
            nominal_formula={"template": f"SrFe(18-x){A}(x)O27", "A": pf(A), "x": pf(x, "1"),
                             "note": "formula template assumed from series naming; x definition vs precursor "
                                     "stoichiometry unconfirmed"},
            precursors={"A_source": pf(precursor), "purity": None, "Sr_source": None, "Fe_source": None},
            synthesis=info["synthesis"], post_treatment=None, geometry=None, mass=None,
            metadata_status={"all_name_derived_fields": PRELIMINARY, "composition_measured": "unknown",
                             "atmosphere": "unknown", "cooling": "unknown", "density_porosity": "unknown"},
            files=files + [], observations=[oid],
            ambiguity_flags=flags + ["x_definition_vs_precursor_unconfirmed",
                                     "same_formula_does_not_imply_same_sample"],
            notes=pair_note)
        if "ZnO" in key and "x=1.5" in key:
            for r in mtd:
                smp.files.append({"archive_path": r["archive_path"], "zone": r["zone"], "sha256": r["sha256"],
                                  "representation": "previous_analysis_project_not_ground_truth",
                                  "link_status": PRELIMINARY})
                r["sample_ids"] = sid
                r["note"] = "previous analysis project; linked to sample by file name only"
        samples.append(smp)
        observations.append(Observation(
            observation_id=oid, sample_id=sid, method="XRD", files=files, semantic_group=oid,
            header_metadata={"source": "mdi first two lines + txt line count (metadata only)", **hdr},
            calibration=None, uncertainty=None, flags=flags + ["mdi_and_txt_are_one_observation"],
            independent=True))
        smp.metadata_status["observation_header"] = "from_instrument_file_header"

    # Moessbauer: one observation per .spc
    spc = sorted((r for r in rows if r["archive_path"].startswith("02_Mossbauer/") and r["extracted"]), key=lambda r: r["archive_path"])
    nested_by_sha = {r["sha256"]: r for r in nested_rows}
    for r in spc:
        name = PurePosixPath(r["archive_path"]).name
        m = RE_MOS.match(name)
        n_smp += 1
        sid, bid, oid = f"SMP-{n_smp:03d}", f"BAT-{n_smp:03d}", f"OBS-MOS-{n_smp - 8:03d}"
        g = m.groupdict()
        x = _num(g["x"])
        hdr = read_spc_header(cfg.vault_lab / Path(*r["archive_path"].split("/")))
        flags = ["x_definition_vs_precursor_unconfirmed", "same_formula_does_not_imply_same_sample",
                 "not_linked_to_xrd_samples_without_proof",
                 "A_site_valence_Mn3+_possible_for_Mn2O3_precursor_outside_Fe2+_Fe3+_model"]
        name_date = f"20{g['y2']}-{g['m2']}-{g['dd']}"
        hdate = hdr.get("Date", "")
        mh = re.match(r"(\d\d)\.(\d\d)\.(\d{4})", hdate)
        if mh and f"{mh.group(3)}-{mh.group(2)}-{mh.group(1)}" != name_date:
            flags.append("filename_date_differs_from_header_date")
        dup = nested_by_sha.get(r["sha256"])
        files = [{"archive_path": r["archive_path"], "zone": r["zone"], "sha256": r["sha256"],
                  "representation": "instrument_native_spc"}]
        if dup:
            files.append({"archive_path": dup["archive_path"], "zone": dup["zone"], "sha256": dup["sha256"],
                          "representation": "byte_duplicate_inside_Archive.zip_not_extracted",
                          "independent": False})
            dup["semantic_group"], dup["sample_ids"], dup["observation_ids"] = oid, sid, oid
        else:
            flags.append("no_byte_duplicate_found_in_Archive.zip")
        r["semantic_group"], r["sample_ids"], r["observation_ids"] = oid, sid, oid
        A = "Mn"
        smp = SampleRecord(
            sample_id=sid, batch_id=bid, label=name.rsplit(".", 1)[0],
            nominal_formula={"template": "SrFe(18-x)Mn(x)O27", "A": pf(A), "x": pf(x, "1"),
                             "note": "x from file name; Mn-series XRD (x=1.1-1.4) is a different set"},
            precursors={"A_source": pf(g["prec"]), "purity": None},
            synthesis={"raw_tokens": name,
                       "stages": [{"T_C": pf(int(g["T"]), "C"),
                                   "duration": pf(int(g["t"]), "h") if g["t"] else None}],
                       "name_series_token_mm_yy": pf(f"{g['mm']}.{g['yy']}"),
                       "name_date_ddmmyy": pf(name_date, None, "meaning of date (synthesis/mounting/start) unknown"),
                       "atmosphere": None, "cooling": None},
            mass={"value": int(g["mass"]), "unit": "mg", "status": PRELIMINARY,
                  "note": "role unconfirmed (sample vs absorber mass)"},
            metadata_status={"all_name_derived_fields": PRELIMINARY, "composition_measured": "unknown",
                             "observation_header": "from_instrument_file_header"},
            files=files, observations=[oid], ambiguity_flags=flags,
            notes=None)
        samples.append(smp)
        observations.append(Observation(
            observation_id=oid, sample_id=sid, method="Mossbauer", files=files, semantic_group=oid,
            header_metadata={"source": "spc header text up to 'Data:' (no channel counts read)",
                             "fields": hdr, "velocity_unit": "unspecified in file header (likely mm/s)"},
            calibration=None, uncertainty=None, flags=flags, independent=True))

    # semantic groups for dup rows / legacy / mtd already set; legacy dup png
    for r in rows:
        if r["byte_dup_group"] and not r["semantic_group"]:
            r["semantic_group"] = r["byte_dup_group"]

    # link unreferenced
    unlinked = [r["archive_path"] for r in rows
                if r["zone"].startswith("vault/lab") and not r["sample_ids"] and r["extracted"]]

    summary = {"n_files": len(entries), "n_rows": len(rows), "n_samples": len(samples),
               "n_observations": len(observations), "dup_groups": gid, "nested": nested_info,
               "unlinked_lab_files": unlinked}
    if write:
        reg.mkdir(parents=True, exist_ok=True)
        cols = list(rows[0].keys())
        buf = io.StringIO(newline="")
        w = csv.DictWriter(buf, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
        with open(reg / "files_manifest.csv", "w", encoding="utf-8", newline="") as f:
            f.write(buf.getvalue())
        hdr_txt = ("# PRELIMINARY sample registry (stage 1). Built ONLY from file names and header text.\n"
                   "# Every name-derived value has status preliminary_from_filename. null = unknown.\n"
                   "# Equal formula != same sample; distinct dates/regimes = distinct samples until proven.\n")
        guard.write_text_keep_mode(reg / "samples.yaml",
            hdr_txt + yaml.safe_dump({"samples": [to_plain(s) for s in samples]}, allow_unicode=True,
                                     sort_keys=False, width=120))
        guard.write_text_keep_mode(reg / "observations.yaml",
            "# Observations: one per measurement; .mdi+.txt of one diffractogram = ONE observation;\n"
            "# Archive.zip members are byte duplicates and not independent.\n"
            + yaml.safe_dump({"observations": [to_plain(o) for o in observations]}, allow_unicode=True,
                             sort_keys=False, width=120))
        guard.write_text_keep_mode(reg / "nested_archive_members.json",
            json.dumps(nested_info, ensure_ascii=False, indent=1))
    return summary


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="Build the PRELIMINARY registry (refuses to overwrite existing files)")
    ap.add_argument("--force", action="store_true",
                    help="overwrite existing registry files; allowed only if `git status` of registry/ is clean")
    ns = ap.parse_args(argv)
    try:
        print(json.dumps(build(force=ns.force), ensure_ascii=False, indent=1))
    except guard.RegistryOverwriteError as e:
        print(f"build_registry: {e}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
