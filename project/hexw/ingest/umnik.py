"""Register the second laboratory source (folder Umnik) in the isolated vault.

Input: vault/lab_validation_umnik/** (raw measurements) and vault/legacy/prior_analyses_umnik/** (prior
analyses), both already copied byte-for-byte by tools/copy_umnik.py.

Reads ONLY: file bytes for sha256, the 62-byte BMP header of .jip files, and (for prior-analysis CSV files)
the column that holds sample identifiers. Numerical content of measurements and prior analyses is never
read, printed or stored; .docx is not opened.

Writes (registry only): files_manifest_umnik.csv, prior_use.yaml; updates observations.yaml and samples.yaml
(existing SampleIDs / ObservationIDs are never changed; no new samples are created here).

Protection (hexw.ingest.registry_guard): the script refuses to run if any of these files exists, unless `--force`
is given AND `git status` of the registry is clean. prior_use.yaml is NEVER regenerated: hand-made content
(user_confirmation, independent_lab_validation_available, statuses set by the user) is kept byte-for-byte and only
entries of observations that are absent from the file are appended (guard.merge_prior_use_text).
"""
from __future__ import annotations

import csv
import io
import re
import struct
from collections import defaultdict
from pathlib import Path

import yaml

from hexw.ingest import registry_guard as guard
from hexw.io import paths
from hexw.registry.prior_use import LAB_METHODS, ROLE_PRIOR_LAB_ANALYSIS

LAB_ZONE = "vault/lab_validation_umnik"
PRIOR_ZONE = "vault/legacy"          # prior analyses live under vault/legacy/prior_analyses_umnik
PRIOR_SUBDIR = "prior_analyses_umnik"
MANIFEST_NAME = "files_manifest_umnik.csv"
# column that holds sample identifiers in prior-analysis CSV files (the ONLY column ever read)
ID_COLUMNS = {"SrW_XRD_lattice.csv": "файл", "w_check_summary.csv": "образец"}
MANIFEST_COLS = ["zone", "umnik_path", "bytes", "sha256", "role", "status", "matched_registered_paths",
                 "name_variant_of_archive", "different_version_of", "sample_ids", "observation_ids",
                 "byte_dup_group_umnik", "origin", "note"]

LEGACY04_NOTE_XRD = ("not_determined: legacy model 04 file names carry no sample references and its content "
                     "was not opened")
LEGACY04_NOTE_MOS = ("possible: legacy model 04 (FMR tables/methodology, Mn series mentioned in a title) may "
                     "have used Mossbauer data; not excluded")
USER_CONFIRMATION_NOTE = "pending: no user confirmation recorded"


def norm_name(name: str) -> str:
    """Name key tolerant to the archive-vs-desktop renames (x=1_5 vs x=1.5, '-' vs '_', degree sign)."""
    return re.sub(r"[^0-9A-Za-z]", "", name).lower()


def read_id_column(path, column: str) -> list[str]:
    """Return only the sample-identifier column of a ';' CSV (no other column is kept)."""
    with paths.open_data(path, "rb", purpose="metadata") as f:
        text = f.read().decode("utf-8-sig")
    return [r[column] for r in csv.DictReader(io.StringIO(text, newline=""), delimiter=";")]


def read_jip_header(path) -> dict:
    """.jip is a 1-bit Windows BMP preview (132x67); only its 54+8 byte header is parsed."""
    b = paths.read_prefix(path, 62, purpose="metadata")
    if b[:2] != b"BM":
        return {"container": "unknown", "contains_instrument_metadata": False}
    bf_size, = struct.unpack_from("<I", b, 2)
    off, = struct.unpack_from("<I", b, 10)
    w, h = struct.unpack_from("<ii", b, 18)
    bpp, = struct.unpack_from("<H", b, 28)
    return {"container": "Windows BMP (DIB) preview image", "width_px": w, "height_px": h,
            "bits_per_pixel": bpp, "declared_size_bytes": bf_size, "pixel_data_offset": off,
            "contains_instrument_metadata": False,
            "note": ("header holds image geometry only; instrument, radiation, range, step and date are "
                     "taken from the linked .mdi header")}


def _rows(path):
    with open(path, encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def scan(cfg: paths.Config) -> list[dict]:
    out = []
    for zone, base in ((LAB_ZONE, cfg.vault_lab_umnik), (PRIOR_ZONE, cfg.prior_analyses_umnik)):
        if not base.is_dir():
            continue
        for p in sorted(base.rglob("*")):
            if p.is_file():
                out.append({"zone": zone, "umnik_path": p.relative_to(base).as_posix(), "path": p,
                            "bytes": p.stat().st_size, "sha256": paths.sha256_file(p)})
    return out


def _is_prior(zone: str) -> bool:
    return zone == PRIOR_ZONE


OUTPUTS = (MANIFEST_NAME, "observations.yaml", "samples.yaml", "prior_use.yaml")


def build(root=None, write: bool = True, force: bool = False) -> dict:
    cfg = paths.get_config() if root is None else paths.configure(root)
    reg = cfg.registry
    if write:
        guard.check_write_allowed(reg, OUTPUTS, force=force)
    fm = [r for r in _rows(reg / "files_manifest.csv") if r["zone"].startswith("vault")]
    obs_doc = yaml.safe_load((reg / "observations.yaml").read_text(encoding="utf-8"))
    smp_doc = yaml.safe_load((reg / "samples.yaml").read_text(encoding="utf-8"))
    observations = [o for o in obs_doc["observations"] if o.get("role") != ROLE_PRIOR_LAB_ANALYSIS]
    samples = smp_doc["samples"]
    lab_obs = [o for o in observations if o["method"] in LAB_METHODS]

    files = scan(cfg)
    # ---- match with the registered archive/legacy files ----
    reg_by_sha = defaultdict(list)
    for r in fm:
        reg_by_sha[r["sha256"]].append(r)
    reg_by_norm = defaultdict(list)
    for r in fm:
        if "::" not in r["archive_path"]:
            reg_by_norm[norm_name(r["archive_path"].rsplit("/", 1)[-1])].append(r)
    obs_by_sha = {f["sha256"]: o["observation_id"] for o in lab_obs for f in o["files"]}
    sample_of_obs = {o["observation_id"]: o["sample_id"] for o in lab_obs}

    by_sha_umnik = defaultdict(list)
    for f in files:
        by_sha_umnik[f["sha256"]].append(f)
    dup_groups = {}
    for sha, v in sorted(by_sha_umnik.items(), key=lambda kv: kv[1][0]["umnik_path"]):
        if len(v) > 1:
            dup_groups[sha] = f"UDUP-{len(dup_groups) + 1:02d}"
    by_name_umnik = defaultdict(list)
    for f in files:
        by_name_umnik[(f["zone"], f["umnik_path"].rsplit("/", 1)[-1])].append(f)

    for f in files:
        name = f["umnik_path"].rsplit("/", 1)[-1]
        m = [r for r in reg_by_sha.get(f["sha256"], [])]
        f["matched"] = [r["archive_path"] for r in m]
        f["match_obs"] = sorted({r["observation_ids"] for r in m if r["observation_ids"]})
        f["match_samples"] = sorted({s for r in m for s in r["sample_ids"].split(";") if s})
        f["variant"] = bool(m) and all(r["archive_path"].rsplit("/", 1)[-1] != name
                                       for r in m if "::" not in r["archive_path"]) \
            and any("::" not in r["archive_path"] for r in m)
        diffs = [f"{r['archive_path']} sha256={r['sha256']}" for r in reg_by_norm.get(norm_name(name), [])
                 if r["sha256"] != f["sha256"]]
        # same-named prior-analysis files in other Umnik folders with other bytes
        for g in by_name_umnik[(f["zone"], name)]:
            if g["sha256"] != f["sha256"]:
                diffs.append(f"umnik:{g['umnik_path']} sha256={g['sha256']}")
        f["diffs"] = diffs
        f["dup_group"] = dup_groups.get(f["sha256"], "")
        ext = name.rsplit(".", 1)[-1].lower()
        f["ext"] = ext
        top = f["umnik_path"].split("/")[0]
        if _is_prior(f["zone"]):
            f["role"] = ROLE_PRIOR_LAB_ANALYSIS
            f["origin"] = "Claude outputs (prior AI session)" if "/Claude outputs/" in "/" + f["umnik_path"] \
                else "XRD_fit (prior analysis)"
        else:
            f["origin"] = top
            f["role"] = {"jip": "instrument_preview_image", "mtd": "previous_analysis_not_ground_truth"}.get(
                ext, "experimental_raw_or_export")

    # ---- .jip -> Observation (by file-name stem, same folder) ----
    lab_files = [f for f in files if not _is_prior(f["zone"])]
    stem_obs = defaultdict(set)
    for f in lab_files:
        if f["ext"] in ("mdi", "txt", "spc") and f["match_obs"]:
            d, _, n = f["umnik_path"].rpartition("/")
            stem_obs[(d, n.rsplit(".", 1)[0])].update(f["match_obs"])
    # fallback key: file name without the instrument time stamp suffix '#1##...' (mdi without txt-style name)
    def base_stem(n):
        return n.split("#1##")[0]
    base_obs = defaultdict(set)
    for (d, s), v in stem_obs.items():
        base_obs[(d, base_stem(s))].update(v)
    for f in lab_files:
        if f["ext"] != "jip":
            continue
        d, _, n = f["umnik_path"].rpartition("/")
        stem = n.rsplit(".", 1)[0]
        cands = stem_obs.get((d, stem)) or set()
        basis = "same_folder_same_file_name_stem"
        if len(cands) != 1:
            cands = base_obs.get((d, base_stem(stem)), set())
            basis = "same_folder_stem_without_timestamp"
        f["jip"] = read_jip_header(f["path"])
        if len(cands) == 1:
            f["link_obs"], f["link_basis"], f["link_flag"] = next(iter(cands)), basis, None
        else:
            f["link_obs"], f["link_basis"], f["link_flag"] = None, None, "ambiguous_or_missing_observation_link"
        f["match_obs"] = [f["link_obs"]] if f["link_obs"] else []
        f["match_samples"] = [sample_of_obs[f["link_obs"]]] if f["link_obs"] else []

    # .mtd -> sample via archive match (already in match_samples)

    # ---- prior analyses: sample linkage from names / identifier columns ----
    xrd_samples = {o["sample_id"] for o in lab_obs if o["method"] == "XRD"}
    smp_elem_x = {}
    for s in samples:
        nf = s["nominal_formula"]
        if s["sample_id"] in xrd_samples:
            smp_elem_x[(nf["A"]["value"], float(nf["x"]["value"]))] = s["sample_id"]
    smp_elem = defaultdict(set)
    for (el, _), sid in smp_elem_x.items():
        smp_elem[el].add(sid)
    txt_name_to_obs = {f["umnik_path"].rsplit("/", 1)[-1]: f["match_obs"] for f in lab_files if f["ext"] == "txt"}

    def ids_from_column(f):
        col = ID_COLUMNS[f["umnik_path"].rsplit("/", 1)[-1]]
        res = set()
        for v in read_id_column(f["path"], col):
            if v in txt_name_to_obs:
                for ob in txt_name_to_obs[v]:
                    res.add(sample_of_obs[ob])
            else:
                m = re.fullmatch(r"([A-Za-z]+)_([0-9.]+)", v.strip())
                if m and (m.group(1), float(m.group(2))) in smp_elem_x:
                    res.add(smp_elem_x[(m.group(1), float(m.group(2)))])
        return res

    prior = [f for f in files if _is_prior(f["zone"])]
    folder_ids = {}
    for f in prior:
        name = f["umnik_path"].rsplit("/", 1)[-1]
        if name in ID_COLUMNS:
            f["sample_set"] = ids_from_column(f)
            f["link_basis"] = "identifier_column_of_csv"
            folder_ids.setdefault(f["umnik_path"].rsplit("/", 1)[0], set()).update(f["sample_set"])
    for f in prior:
        if "sample_set" in f:
            continue
        d, _, name = f["umnik_path"].rpartition("/")
        sib = folder_ids.get(d, set())
        m = re.fullmatch(r"xrd_fits_([A-Za-z]+)\.png", name)
        if m and m.group(1) in smp_elem:
            f["sample_set"], f["link_basis"] = sib & smp_elem[m.group(1)], "filename_element_token"
        else:
            f["sample_set"], f["link_basis"] = set(sib), "sibling_csv_ids_in_same_folder_unverified"

    # ---- artifacts (one per distinct hash) ----
    art_by_sha, artifacts = {}, []
    for f in sorted(prior, key=lambda x: x["umnik_path"]):
        if f["sha256"] not in art_by_sha:
            art_by_sha[f["sha256"]] = {"observation_id": f"PRI-{len(artifacts) + 1:03d}", "copies": []}
            artifacts.append(art_by_sha[f["sha256"]])
        art_by_sha[f["sha256"]]["copies"].append(f)
    sample_ids_of_art = {}
    for a in artifacts:
        ss = sorted({s for c in a["copies"] for s in c["sample_set"]})
        sample_ids_of_art[a["observation_id"]] = ss
        for c in a["copies"]:
            c["art_id"] = a["observation_id"]
            c["match_samples"] = sorted(c["sample_set"])
    # ---- manifest rows ----
    rows = []
    for f in files:
        if f["matched"]:
            status = "byte_identical_to_archive"
        elif f["diffs"] and not _is_prior(f["zone"]):
            status = "different_version"
        else:
            status = "new"
        note = []
        if f["variant"]:
            note.append("same bytes as archive file under a different file name")
        if f.get("jip"):
            note.append("jip: " + f["jip"]["container"] + f" {f['jip']['width_px']}x{f['jip']['height_px']}"
                        f"x{f['jip']['bits_per_pixel']}bit; " + (f["link_basis"] or f["link_flag"]))
        if _is_prior(f["zone"]):
            note.append(f"{f['art_id']}; sample link: {f['link_basis']}; content not read")
        rows.append({
            "zone": f["zone"], "umnik_path": f["umnik_path"], "bytes": f["bytes"], "sha256": f["sha256"],
            "role": f["role"], "status": status, "matched_registered_paths": ";".join(f["matched"]),
            "name_variant_of_archive": f["variant"], "different_version_of": " | ".join(f["diffs"]),
            "sample_ids": ";".join(f["match_samples"]),
            "observation_ids": ";".join(f["match_obs"]) if not _is_prior(f["zone"]) else f["art_id"],
            "byte_dup_group_umnik": f["dup_group"], "origin": f["origin"], "note": "; ".join(note)})

    # ---- observations / samples ----
    for o in observations:
        o.pop("umnik_copies", None)
        o.pop("auxiliary_files", None)
        copies, aux = [], []
        for f in lab_files:
            if f["ext"] == "jip" and f["link_obs"] == o["observation_id"]:
                aux.append({"zone": f["zone"], "umnik_path": f["umnik_path"], "sha256": f["sha256"],
                            "representation": "instrument_preview_bitmap_jip", "link_basis": f["link_basis"],
                            "jip_header": f["jip"]})
            elif f["ext"] in ("mdi", "txt", "spc") and o["observation_id"] in f["match_obs"]:
                copies.append({"zone": f["zone"], "umnik_path": f["umnik_path"], "sha256": f["sha256"],
                               "same_observation_as_archive": True,
                               "renamed_vs_archive": f["variant"]})
        if copies:
            o["umnik_copies"] = copies
        if aux:
            o["auxiliary_files"] = aux
    pri_obs = []
    for a in artifacts:
        cs = a["copies"]
        pri_obs.append({
            "observation_id": a["observation_id"], "sample_id": None,
            "sample_ids": sample_ids_of_art[a["observation_id"]],
            "method": "prior_analysis_artifact", "role": ROLE_PRIOR_LAB_ANALYSIS,
            "files": [{"zone": c["zone"], "umnik_path": "prior_analyses_umnik/" + c["umnik_path"],
                       "sha256": c["sha256"], "origin": c["origin"]} for c in cs],
            "sample_link_basis": sorted({c["link_basis"] for c in cs}),
            "independent": False, "duplicate_of": None,
            "flags": ["derived_processing_result_not_reference", "not_an_experiment",
                      "numerical_content_not_used_by_registry"]
                     + (["prior_ai_session_output"] if any("Claude" in c["origin"] for c in cs) else [])
                     + (["same_name_other_bytes_exists"] if any(c["diffs"] for c in cs) else [])})
    # ---- samples ----
    for s in samples:
        s.pop("umnik_copies", None)
        s.pop("prior_analysis_artifacts", None)
        sid = s["sample_id"]
        cp = [{"zone": f["zone"], "umnik_path": f["umnik_path"], "sha256": f["sha256"], "role": f["role"]}
              for f in lab_files if sid in f["match_samples"]]
        if cp:
            s["umnik_copies"] = cp
        pa = [a["observation_id"] for a in artifacts if sid in sample_ids_of_art[a["observation_id"]]]
        if pa:
            s["prior_analysis_artifacts"] = pa

    # ---- prior_use ----
    mtd_by_sample = defaultdict(list)
    for r in fm:
        if r["data_manifest_role"] == "previous_analysis_not_ground_truth" and "::" not in r["archive_path"]:
            for sid in r["sample_ids"].split(";"):
                if sid:
                    mtd_by_sample[sid].append({"path": r["archive_path"], "sha256": r["sha256"]})
    entries = []
    for o in lab_obs:
        sid = o["sample_id"]
        arts = [a["observation_id"] for a in artifacts if sid in sample_ids_of_art[a["observation_id"]]]
        mtd = mtd_by_sample.get(sid, [])
        if o["method"] == "XRD" and (arts or mtd):
            status, basis = "used_prior", ("prior Umnik analyses ({}) name/identify this sample{}; default for XRD "
                                           "samples with prior analyses is retrospective").format(
                ", ".join(arts) or "none", "; .mtd projects exist" if mtd else "")
            legacy = LEGACY04_NOTE_XRD
        elif o["method"] == "XRD":
            status, basis, legacy = "unknown_prior_use", "no prior analysis found by names, not confirmed", \
                LEGACY04_NOTE_XRD
        else:
            status = "unknown_prior_use"
            basis = ("no prior analysis names this sample, but legacy model 04 may have used Mossbauer data; "
                     "stays unknown until the user confirms")
            legacy = LEGACY04_NOTE_MOS
        entries.append({
            "observation_id": o["observation_id"], "sample_id": sid, "method": o["method"],
            "prior_use_status": status, "validation_mode": "retrospective",
            "evidence": {"prior_analyses_umnik": arts,
                         "mtd_projects": mtd,
                         "legacy_model_04": legacy},
            "basis": basis,
            "user_confirmation": None})
    prior_header = (
        "# Prior-use disclosure for laboratory Observations (stage 1b).\n"
        "# prior_use_status: used_prior | unknown_prior_use | no_prior_use_found.\n"
        "# validation_mode: retrospective | independent_candidate. independent_candidate is allowed ONLY with\n"
        "# no_prior_use_found AND a user_confirmation record; everything else is retrospective.\n"
        "# This status is copied into every ValidationResult (hexw.registry.prior_use.apply_prior_use).\n"
        "# Evidence: file names and identifier columns. CSV files of prior analyses were read IN FULL into memory\n"
        "# when the identifier column was extracted; only that column of identifiers was output or used.\n")
    prior_path = reg / "prior_use.yaml"
    if prior_path.exists():          # never regenerate: keep manual fields, append only missing observations
        prior_doc, prior_report = guard.merge_prior_use_text(guard.read_text_exact(prior_path), entries)
    else:
        prior_doc = prior_header + yaml.safe_dump({"schema_version": 1, "entries": entries}, allow_unicode=True,
                                                  sort_keys=False, width=120)
        prior_report = {"kept": 0, "added": [e["observation_id"] for e in entries]}

    summary = {
        "n_files": len(files), "n_lab": len(lab_files), "n_prior": len(prior),
        "n_byte_identical_to_archive": sum(r["status"] == "byte_identical_to_archive" for r in rows),
        "new": [r["umnik_path"] for r in rows if r["status"] == "new" and r["zone"] == LAB_ZONE],
        "different_versions": [r["umnik_path"] for r in rows if r["status"] == "different_version"],
        "n_artifacts": len(artifacts), "prior_use": prior_report}
    if write:
        guard.check_write_allowed(reg, OUTPUTS, force=force)       # re-check right before the first write
        with open(reg / MANIFEST_NAME, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=MANIFEST_COLS)
            w.writeheader()
            w.writerows(rows)
        hdr = ("# Observations: one per measurement; .mdi+.txt of one diffractogram = ONE observation;\n"
               "# Archive.zip members are byte duplicates and not independent.\n"
               "# Umnik byte-identical copies are listed under umnik_copies (same Observation, not new\n"
               "# measurements); PRI-* entries are prior_lab_analysis artifacts (derived, not reference).\n")
        guard.write_text_keep_mode(reg / "observations.yaml",
            hdr + yaml.safe_dump({"observations": observations + pri_obs}, allow_unicode=True,
                                 sort_keys=False, width=120))
        old = (reg / "samples.yaml").read_text(encoding="utf-8")
        hdr_s = "".join(l + "\n" for l in old.splitlines() if l.startswith("#"))
        guard.write_text_keep_mode(reg / "samples.yaml",
            hdr_s + yaml.safe_dump({"samples": samples}, allow_unicode=True, sort_keys=False, width=120))
        if prior_path.exists() and prior_doc == guard.read_text_exact(prior_path):
            pass                                                   # merge added nothing: file untouched
        else:                                                      # old bytes + appended block (or a new file)
            guard.write_prior_use_checked(prior_path, prior_doc)      # parsed after writing, rolled back on failure
    return summary


def main(argv=None) -> int:
    import argparse
    import json
    ap = argparse.ArgumentParser(description="Register the Umnik lab source (refuses to overwrite the registry)")
    ap.add_argument("--force", action="store_true",
                    help="regenerate existing registry files; allowed only if `git status` of registry/ is clean. "
                         "prior_use.yaml is merged, never overwritten")
    ns = ap.parse_args(argv)
    try:
        print(json.dumps(build(force=ns.force), ensure_ascii=False, indent=1))
    except guard.RegistryOverwriteError as e:
        print(f"umnik: {e}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
