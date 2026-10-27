"""`hexw audit`: integrity (hashes vs manifests, duplicates, read-only) + leak-guard self-checks."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import os
import stat
import sys
import tempfile
import zipfile
from collections import defaultdict
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

import yaml

from hexw.io import paths
from hexw.registry.schemas import (
    DataRef, Metric, ParameterEstimate, Prediction, Provenance, RunManifest, RunState, SourceRole,
    Uncertainty, ValidationResult, dump_json,
)

NESTED_REL = "05_original_duplicate_archive/Archive.zip"
INPLACE_PRIOR_PREFIX = "03_previous_XRD_projects/"   # prior .mtd projects: role prior_lab_analysis, NOT moved
PKGS = ["numpy", "scipy", "matplotlib", "pandas", "pyyaml", "pytest", "openpyxl", "pymatgen", "ase"]


class Report:
    def __init__(self):
        self.checks: list[dict] = []

    def add(self, name, ok, detail=""):
        self.checks.append({"name": name, "ok": bool(ok), "detail": detail})

    @property
    def ok(self):
        return all(c["ok"] for c in self.checks)


def _rows(path):
    with open(path, encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def _is_readonly(p: Path) -> bool:
    st = p.stat()
    if hasattr(st, "st_file_attributes"):
        return bool(st.st_file_attributes & stat.FILE_ATTRIBUTE_READONLY)
    return not os.access(p, os.W_OK)


def check_integrity(rep: Report, cfg: paths.Config) -> dict:
    ze = json.loads((cfg.registry / "zip_entries.json").read_text(encoding="utf-8"))
    entries = ze["entries"]
    # original zip
    zips = list(cfg.original_zip_dir.glob("*.zip"))
    zsha = paths.sha256_file(zips[0]) if len(zips) == 1 else None
    rep.add("original zip copy present and sha256 == recorded", zsha == ze["zip_sha256"], f"sha256={zsha}")

    all_m = {r["path"]: r for r in _rows(cfg.raw_theory / "ALL_FILES_MANIFEST.csv")}
    data_m = {r["archive_path"]: r for r in _rows(cfg.raw_theory / "metadata" / "DATA_MANIFEST.csv")}
    missing, mism_all, mism_data, mism_zip, by_sha = [], [], [], [], defaultdict(list)
    zone_dirs = {"raw/theory_side": cfg.raw_theory, "vault/lab_validation": cfg.vault_lab,
                 "vault/legacy": cfg.vault_legacy}
    expected = set()
    for e in entries:
        p = zone_dirs[e["zone"]] / Path(*e["archive_path"].split("/"))
        expected.add(os.path.normcase(str(p)))
        if not p.is_file():
            missing.append(e["archive_path"])
            continue
        sha = paths.sha256_file(p)
        by_sha[sha].append(e["archive_path"])
        if sha != e["sha256"]:
            mism_zip.append(e["archive_path"])
        a = all_m.get(e["archive_path"])
        if a and a["sha256"].lower() != sha:
            mism_all.append(e["archive_path"])
        d = data_m.get(e["archive_path"])
        if d and d["sha256"].lower() != sha:
            mism_data.append(e["archive_path"])
    extra = []
    umnik_prior = os.path.normcase(str(cfg.prior_analyses_umnik)) + os.sep   # audited by check_umnik
    for d in zone_dirs.values():
        for q in d.rglob("*"):
            nq = os.path.normcase(str(q))
            if q.is_file() and nq not in expected and not nq.startswith(umnik_prior):
                extra.append(str(q))
    rep.add(f"all {len(entries)} files present on disk", not missing, f"missing={missing}")
    rep.add("no extra files in raw/theory_side, vault/*", not extra, f"extra={extra}")
    rep.add("sha256 == zip member content", not mism_zip, f"mismatches={mism_zip}")
    rep.add(f"sha256 == ALL_FILES_MANIFEST.csv ({len(all_m)} rows)", not mism_all, f"mismatches={mism_all}")
    rep.add(f"sha256 == DATA_MANIFEST.csv ({len(data_m)} rows)", not mism_data, f"mismatches={mism_data}")
    not_in_all = [e["archive_path"] for e in entries if e["archive_path"] not in all_m]
    rep.add("every file except ALL_FILES_MANIFEST.csv itself is listed in it",
            not_in_all == ["ALL_FILES_MANIFEST.csv"], f"unlisted={not_in_all}")

    # byte duplicates on disk
    dups = {s: v for s, v in by_sha.items() if len(v) > 1}
    png = sorted(v for v in dups.values() if all(x.endswith(".png") for x in v))
    rep.add("byte duplicate on disk: SrW_series_vs_x.png == SrW_series_vs_x_1.png",
            len(dups) == 1 and len(png) == 1, f"groups={list(dups.values())}")

    # nested archive: hashed in memory, not extracted
    with paths.open_data(cfg.vault_lab / NESTED_REL, "rb", purpose="integrity") as f:
        zb = f.read()
    nested = {}
    with zipfile.ZipFile(io.BytesIO(zb)) as nz:
        for i in nz.infolist():
            if not i.is_dir():
                nested[i.filename] = hashlib.sha256(nz.read(i)).hexdigest()
    spc = {p for p in by_sha.values() for p in p if p.startswith("02_Mossbauer/")}
    spc_sha = {s for s, v in by_sha.items() if v[0].startswith("02_Mossbauer/")}
    rep.add("Archive.zip not extracted (members only hashed in memory)",
            [x.name for x in (cfg.vault_lab / "05_original_duplicate_archive").iterdir()] == ["Archive.zip"])
    rep.add("Archive.zip: 5 members, each byte-identical to a distinct 02_Mossbauer/*.spc",
            len(nested) == 5 and set(nested.values()) == spc_sha and len(spc_sha) == 5,
            f"members={len(nested)}")
    return {"n_files": len(entries), "dups_on_disk": list(dups.values())}


def check_readonly(rep: Report, cfg: paths.Config):
    bad = []
    for base in (cfg.raw, cfg.vault):
        for q in base.rglob("*"):
            if q.is_file() and not _is_readonly(q):
                bad.append(str(q))
    rep.add("all files in raw/ and vault/ are read-only", not bad, f"writable={bad[:5]} (total {len(bad)})")


def check_registry(rep: Report, cfg: paths.Config):
    try:
        samples = yaml.safe_load((cfg.registry / "samples.yaml").read_text(encoding="utf-8"))["samples"]
        obs = yaml.safe_load((cfg.registry / "observations.yaml").read_text(encoding="utf-8"))["observations"]
        fm = _rows(cfg.registry / "files_manifest.csv")
    except Exception as e:  # noqa
        rep.add("registry files readable", False, repr(e))
        return
    rep.add("registry files readable", True, f"{len(samples)} samples, {len(obs)} observations")
    ids = [s["sample_id"] for s in samples]
    rep.add("SampleID unique, BatchID present", len(ids) == len(set(ids)) and all(s["batch_id"] for s in samples))
    rep.add("every sample has files and metadata_status",
            all(s["files"] and s["metadata_status"] for s in samples))
    sha_ok = {r["archive_path"]: r["sha256"] for r in fm}
    rep.add("sample file hashes agree with files_manifest.csv",
            all(sha_ok.get(f["archive_path"]) == f["sha256"] for s in samples for f in s["files"]))
    linked = {f["archive_path"] for s in samples for f in s["files"]}
    unlinked = [r["archive_path"] for r in fm if r["zone"] == "vault/lab_validation" and r["extracted"] == "True"
                and r["archive_path"] not in linked]
    rep.add("every lab file linked to a sample (except the Archive.zip container)",
            unlinked == [NESTED_REL], f"unlinked={unlinked}")
    rep.add("name-derived values carry status preliminary_from_filename",
            all(s["metadata_status"].get("all_name_derived_fields") == "preliminary_from_filename" for s in samples))
    xrd = [o for o in obs if o["method"] == "XRD"]
    rep.add(".mdi and .txt of one diffractogram form ONE observation",
            len(xrd) == 8 and all(len(o["files"]) == 2 for o in xrd), f"{len(xrd)} XRD observations")
    mos = [o for o in obs if o["method"] == "Mossbauer"]
    rep.add("Mossbauer: 5 observations, not linked to XRD samples",
            len(mos) == 5 and all("not_linked_to_xrd_samples_without_proof" in o["flags"] for o in mos))
    # distinct batches
    rep.add("equal formula does not merge samples (13 samples, 13 batches)",
            len({s["batch_id"] for s in samples}) == len(samples) == 13)


def check_umnik(rep: Report, cfg: paths.Config):
    """Second lab source (Umnik): hashes, cross-source duplicates, prior analyses only in vault, prior_use."""
    from hexw.registry import prior_use as pu_mod
    try:
        rows = _rows(cfg.registry / "files_manifest_umnik.csv")
        fm = _rows(cfg.registry / "files_manifest.csv")
        obs = yaml.safe_load((cfg.registry / "observations.yaml").read_text(encoding="utf-8"))["observations"]
        samples = yaml.safe_load((cfg.registry / "samples.yaml").read_text(encoding="utf-8"))["samples"]
        entries = pu_mod.load_prior_use(cfg.registry)
    except Exception as e:  # noqa
        rep.add("umnik registry files readable (files_manifest_umnik.csv, prior_use.yaml)", False, repr(e))
        return
    zone_dirs = {"vault/lab_validation_umnik": cfg.vault_lab_umnik, "vault/legacy": cfg.prior_analyses_umnik}
    missing, bad_sha, expected = [], [], set()
    for r in rows:
        base = zone_dirs.get(r["zone"])
        p = base / Path(*r["umnik_path"].split("/")) if base else None
        if p is None or not p.is_file():
            missing.append(r["umnik_path"])
            continue
        expected.add(os.path.normcase(str(p)))
        if paths.sha256_file(p) != r["sha256"]:
            bad_sha.append(r["umnik_path"])
    extra = [str(q) for b in zone_dirs.values() if b.is_dir() for q in b.rglob("*")
             if q.is_file() and os.path.normcase(str(q)) not in expected]
    rep.add(f"umnik: all {len(rows)} files present, sha256 == files_manifest_umnik.csv",
            bool(rows) and not missing and not bad_sha, f"missing={missing} mismatched={bad_sha}")
    rep.add("umnik: no unregistered files in vault/lab_validation_umnik or vault/legacy/prior_analyses_umnik",
            not extra, f"extra={extra}")

    # prior analyses live only under vault/legacy/prior_analyses_umnik
    prior_all = [r for r in rows if r["role"] == pu_mod.ROLE_PRIOR_LAB_ANALYSIS]
    prior_sha = {r["sha256"] for r in prior_all}

    def inplace(r):          # declared in-place prior artefact (.mtd of the archive), physically not moved
        mp = [x for x in r["matched_registered_paths"].split(";") if x]
        return r["zone"] == "vault/lab_validation_umnik" and bool(mp) and all(
            x.startswith(INPLACE_PRIOR_PREFIX) for x in mp)
    prior_rows = [r for r in prior_all if r["zone"] == "vault/legacy"]
    prior_inplace = [r for r in prior_all if inplace(r)]
    elsewhere = [r["umnik_path"] for r in rows if r["zone"] == "vault/lab_validation_umnik"
                 and (r["role"] == pu_mod.ROLE_PRIOR_LAB_ANALYSIS or r["sha256"] in prior_sha) and not inplace(r)]
    elsewhere += [r["archive_path"] for r in fm if r["zone"] != "vault/legacy" and r["sha256"] in prior_sha
                  and not r["archive_path"].startswith(INPLACE_PRIOR_PREFIX)]
    rep.add("umnik: prior analyses (prior_lab_analysis) exist only in vault/legacy/prior_analyses_umnik "
            "(except the declared in-place .mtd projects)",
            bool(prior_rows) and len(prior_rows) + len(prior_inplace) == len(prior_all) and not elsewhere,
            f"misplaced={elsewhere}")

    # cross-source duplicates are the SAME Observation, never a new one
    reg_by_sha = defaultdict(list)
    for r in fm:
        if r["zone"].startswith("vault"):
            reg_by_sha[r["sha256"]].append(r)
    lab_obs = [o for o in obs if o["method"] in pu_mod.LAB_METHODS]
    lab_ids = {o["observation_id"] for o in lab_obs}
    obs_of_sha = defaultdict(set)
    for o in lab_obs:
        for f in o["files"]:
            obs_of_sha[f["sha256"]].add(o["observation_id"])
    bad = []
    for r in rows:
        if r["zone"] != "vault/lab_validation_umnik" or r["role"] == "instrument_preview_image":
            continue
        m = reg_by_sha.get(r["sha256"], [])
        want = {x["observation_ids"] for x in m if x["observation_ids"]}
        got = set(filter(None, r["observation_ids"].split(";")))
        if m and (r["status"] != "byte_identical_to_archive" or (want and got != want)):
            bad.append(r["umnik_path"])
        if not m and r["status"] not in ("new", "different_version"):
            bad.append(r["umnik_path"])
    rep.add("umnik: files byte-identical to registered archive files map to the same Observation (not new)",
            not bad, f"inconsistent={bad}")
    rep.add("every lab file sha256 belongs to exactly one lab Observation (no duplicate observations)",
            all(len(v) == 1 for v in obs_of_sha.values()) and len(lab_ids) == len(lab_obs),
            f"shared={[s for s, v in obs_of_sha.items() if len(v) > 1]}")
    # different_version: status and different_version_of must agree, and the reference must resolve to a manifest row
    by_path = {r["umnik_path"]: r for r in rows}
    dv_bad = []
    for r in rows:
        ref, st = r["different_version_of"], r["status"]
        if (st == "different_version") != bool(ref):
            dv_bad.append(r["umnik_path"])
            continue
        if ref:
            for part in ref.split(" | "):
                m = re.fullmatch(r"umnik:(.+) sha256=([0-9a-f]{6,64})", part.strip())
                other = by_path.get(m.group(1)) if m else None
                if not other or not other["sha256"].startswith(m.group(2)) or other["sha256"] == r["sha256"]:
                    dv_bad.append(r["umnik_path"])
                    break
    rep.add("umnik: status different_version <=> different_version_of, reference resolves to another file with "
            "different bytes", not dv_bad, f"inconsistent={dv_bad}")
    umnik_sha = defaultdict(list)
    for r in rows:
        umnik_sha[r["sha256"]].append(r)
    undeclared = [[x["umnik_path"] for x in v] for v in umnik_sha.values()
                  if len(v) > 1 and not all(x["byte_dup_group_umnik"] for x in v)]
    rep.add("umnik: byte duplicates between its own files are declared (byte_dup_group_umnik)",
            not undeclared, f"undeclared={undeclared}")
    jips = [r for r in rows if r["umnik_path"].lower().endswith(".jip")]
    jip_sha = {r["sha256"] for r in jips}
    rep.add(".jip previews: each linked to exactly one Observation; not registered as measurement files",
            bool(jips) and all(r["observation_ids"] in lab_ids for r in jips)
            and not any(f["sha256"] in jip_sha for o in lab_obs for f in o["files"]), f"{len(jips)} jip")

    # prior-analysis artifacts are registered with their role
    pri = [o for o in obs if o.get("role") == pu_mod.ROLE_PRIOR_LAB_ANALYSIS]
    sid = {s["sample_id"] for s in samples}
    def pri_zone_ok(f):
        return f["zone"] == "vault/legacy" or (
            f["zone"] == "vault/lab_validation" and f.get("archive_path", "").startswith(INPLACE_PRIOR_PREFIX)) or (
            f["zone"] == "vault/lab_validation_umnik" and any(
                f["sha256"] == r["sha256"] for r in prior_inplace))
    rep.add("prior analyses (incl. in-place .mtd) registered as role prior_lab_analysis, not independent",
            bool(pri) and all(not o["independent"] and o["method"] not in pu_mod.LAB_METHODS
                              and all(pri_zone_ok(f) for f in o["files"])
                              and set(o["sample_ids"]) <= sid for o in pri)
            and {f["sha256"] for o in pri for f in o["files"]} == prior_sha, f"{len(pri)} artifacts")

    # prior_use
    probs = {k: pu_mod.problems(v) for k, v in entries.items() if pu_mod.problems(v)}
    rep.add("prior_use: every lab Observation has a record (no extras)", set(entries) == lab_ids,
            f"missing={sorted(lab_ids - set(entries))} extra={sorted(set(entries) - lab_ids)}")
    rep.add("prior_use: statuses/modes consistent (independent_candidate needs no_prior_use_found + user confirmation)",
            not probs, f"{probs}")
    with_prior = {s for o in pri for s in o["sample_ids"]}
    rep.add("prior_use: XRD with prior analyses -> used_prior/retrospective; Mossbauer -> unknown_prior_use "
            "until user confirms",
            all((e["prior_use_status"], e["validation_mode"]) == ("used_prior", "retrospective")
                for e in entries.values() if e["method"] == "XRD" and e["sample_id"] in with_prior)
            and all(e["prior_use_status"] == "unknown_prior_use" or e.get("user_confirmation")
                    for e in entries.values() if e["method"] == "Mossbauer"))
    summary = pu_mod.independent_validation_available(cfg.registry)
    any_indep = any(e.get("validation_mode") == "independent_candidate" for e in entries.values())
    rep.add("prior_use: summary independent_lab_validation_available is present and consistent with the entries",
            isinstance(summary, bool) and summary == any_indep, f"summary={summary!r} any_independent={any_indep}")


def check_registry_git(rep: Report, cfg: paths.Config) -> dict:
    """Registry under git: SHA256SUMS matches, tree clean, key files read-only; HEADs of both repositories."""
    from hexw.registry import integrity as I
    reg = cfg.registry
    head = I.git_head(reg)
    dirty = I.git_dirty(reg)
    rep.add("registry is a git repository with a commit", head is not None, f"HEAD={head}")
    rep.add("registry git tree is clean (nothing modified or untracked)", dirty == [], f"dirty={dirty}")
    probs = I.verify_sums(reg)
    rep.add("registry SHA256SUMS matches the files exactly", not probs, f"{probs}")
    writable = [n for n in I.KEY_READONLY if (reg / n).is_file() and not _is_readonly(reg / n)]
    rep.add("registry key files (samples, observations, prior_use) are read-only", not writable,
            f"writable={writable}")
    proj = Path(__file__).resolve().parents[1]
    proj_head = I.git_head(proj) if (proj / ".git").exists() else None
    proj_dirty = I.git_dirty(proj) if proj_head else None
    rep.add("project git repository has a commit (HEAD recorded)", proj_head is not None, f"HEAD={proj_head}")
    sums = I.read_sums(reg) if (reg / I.SUMS).is_file() else {}
    return {"git_head": {"project": proj_head, "registry": head},
            "project_tree_dirty": bool(proj_dirty), "registry_checksums": dict(sums)}


# ---------------- leak-guard self-checks ----------------
def _expect(exc, fn):
    try:
        fn()
    except exc:
        return True
    except Exception as e:  # noqa
        return False
    return False


def _probe(modname: str, src: str, **g):
    scope = {"__name__": modname, **g}
    exec(compile(src, f"<probe {modname}>", "exec"), scope)
    return scope


def _sha(p: Path, purpose=None) -> str:
    h = hashlib.sha256()
    f = paths.open_data(p, "rb", purpose=purpose) if purpose else open(p, "rb")
    with f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _check_sums(base: Path, purpose=None) -> list:
    """Verify base/SHA256SUMS; returns the list of problems."""
    sums = base / "SHA256SUMS"
    if not sums.is_file():
        return [f"{sums} missing"]
    bad = []
    listed = set()
    for ln in (paths.open_data(sums, "rb", purpose=purpose) if purpose else open(sums, "rb")).read().decode(
            "utf-8").splitlines():
        if not ln.strip():
            continue
        digest, rel = ln.split("  ", 1)
        listed.add(rel)
        q = base / rel
        if not q.is_file() or _sha(q, purpose) != digest:
            bad.append(rel)
    on_disk = {q.relative_to(base).as_posix() for q in base.rglob("*")
               if q.is_file() and q.name != "SHA256SUMS" and ".git" not in q.relative_to(base).parts}
    bad += [f"unlisted: {r}" for r in sorted(on_disk - listed)]
    return bad


def check_literature(rep: Report, cfg: paths.Config):
    """Stage 3: fit / holdout split of the literature base (literature/SPLIT_DECISION.md)."""
    import re as _re
    hold, lit = cfg.literature_holdout, cfg.literature
    # 1. HOLDOUT zone protected
    rep.add("HOLDOUT: literature_holdout/ exists and is classified as zone HOLDOUT",
            hold.is_dir() and paths.classify(hold / "records_holdout.yaml") is paths.Zone.HOLDOUT)
    declared = list(cfg.holdout_sources)
    rep.add("HOLDOUT: holdout sources are declared in project/config/holdout_sources.json",
            bool(declared) and cfg.holdout_config_file.is_file())
    missing = [str(q) for q in declared if not q.exists()]
    rep.add("HOLDOUT: every declared holdout source exists", not missing, f"missing: {missing}")
    wrong = [str(q) for q in declared if paths.classify(q) is not paths.Zone.HOLDOUT]
    rep.add("HOLDOUT: every declared holdout source is classified HOLDOUT", not wrong, str(wrong))
    first_src = next((f for q in declared if q.exists() for f in ([q] if q.is_file() else sorted(q.rglob("*")))
                      if f.is_file()), None)
    probe = _probe("hexw.theory._probe", "def go(p):\n    return open(p, 'rb')\n")
    targets = [hold / "records_holdout.yaml"] + ([first_src] if first_src else [])
    rep.add("HOLDOUT: hexw.theory open() on literature_holdout and on a holdout source is blocked",
            all(_expect(paths.TheoryAccessError, lambda t=t: probe["go"](t)) for t in targets))
    probe_p = _probe("hexw.params._probe", "from hexw.io import paths\ndef go(p):\n"
                     "    return paths.open_data(p, 'rb', purpose='metadata')\n")
    rep.add("HOLDOUT: hexw.params open_data on literature_holdout is blocked",
            _expect(paths.TheoryAccessError, lambda: probe_p["go"](hold / "records_holdout.yaml")))
    rep.add("HOLDOUT: write into the zone is refused",
            _expect(paths.ReadOnlyZoneError, lambda: paths.open_data(hold / "x.bin", "wb")))
    ro_bad = [q.name for q in hold.rglob("*") if q.is_file() and os.access(q, os.W_OK)]
    rep.add("HOLDOUT: files of literature_holdout are read-only", not ro_bad, str(ro_bad[:5]))
    rep.add("HOLDOUT: SHA256SUMS of literature_holdout verified (purpose integrity)",
            not (bad := _check_sums(hold, "integrity")), str(bad[:5]))

    # 2./3. fit part
    fit_file = lit / "fit" / "records_fit.yaml"
    fit_claims, schema_err = [], None
    try:
        with open(fit_file, encoding="utf-8") as f:
            fdoc = yaml.safe_load(f)
        from hexw.params.literature import load_fit_records
        fit_claims = load_fit_records()
        n_decl = fdoc["meta"]["n_claims"]
    except Exception as e:  # noqa
        schema_err, n_decl = repr(e), -1
    rep.add("literature/fit: every fit record is valid by schema (LiteratureClaim)",
            schema_err is None and len(fit_claims) == n_decl and n_decl > 0, schema_err or f"{len(fit_claims)} vs {n_decl}")
    rep.add("literature/fit: no record with split != fit",
            bool(fit_claims) and all(c.split == "fit" and c.role is SourceRole.LITERATURE_FIT for c in fit_claims))
    stray = [q.relative_to(lit).as_posix() for q in (lit / "fit").rglob("*") if q.is_file()
             and _re.search(r"holdout|_original_mixed", q.name)]
    rep.add("literature/fit: no holdout / original-mixed files inside", not stray, str(stray))
    ex = lit / "excluded" / "records_excluded.yaml"
    try:
        with open(ex, encoding="utf-8") as f:
            exdoc = yaml.safe_load(f)
        ok_ex = all(c["role"] in ("excluded_secondary", "excluded_low_quality", "methodology_only")
                    and not c.get("estimate") for c in exdoc["claims"])
    except Exception:  # noqa
        ok_ex = False
    rep.add("literature/excluded: only excluded_* / methodology_only roles, no parameter estimates", ok_ex)

    # 4. SPLIT_DECISION.md
    dec = lit / "SPLIT_DECISION.md"
    ok_dec, detail = dec.is_file(), "SPLIT_DECISION.md missing"
    if ok_dec:
        txt = dec.read_text(encoding="utf-8")
        m = _re.search(r"HASHES:BEGIN -->(.*?)<!-- HASHES:END", txt, _re.S)
        want = dict(ln.split(": ", 1) for ln in (m.group(1).strip().splitlines() if m else []) if ": " in ln)
        have = {"records_original_mixed_sha256": hold / "_original_mixed" / "records.yaml",
                "records_fit_sha256": fit_file, "records_excluded_sha256": ex,
                "records_holdout_sha256": hold / "records_holdout.yaml"}
        bad = [k for k, q in have.items() if want.get(k) != (_sha(q, "integrity") if hold in q.parents else _sha(q))]
        ok_dec, detail = (bool(want) and not bad), f"hash mismatch or missing: {bad or 'no HASHES block'}"
    rep.add("SPLIT_DECISION.md exists and the record hashes match the files", ok_dec, detail)
    rep.add("literature/ SHA256SUMS verified", not (bad := _check_sums(lit)), str(bad[:5]))

    # 5. agent-os is ignored (personal folder of the user)
    root_pkg = Path(__file__).resolve().parent
    hits = [q.name for q in root_pkg.rglob("*.py") if q.name != "audit.py"
            and _re.search(r"agent[-_]os", q.read_text(encoding="utf-8", errors="replace"), _re.I)]
    sums_text = "".join((b / "SHA256SUMS").read_text(encoding="utf-8") for b in (lit, hold) if (b / "SHA256SUMS").is_file())
    cfg_text = cfg.holdout_config_file.read_text(encoding="utf-8") if cfg.holdout_config_file.is_file() else ""
    rep.add("agent-os is ignored (no code, config, checksum list or zone refers to it)",
            not hits and "agent-os" not in sums_text and "agent-os" not in cfg_text
            and paths.classify(cfg.root / "agent-os" / "x") is paths.Zone.OTHER, str(hits))


_REF_SKIP_DIRS = frozenset({"vault", "literature_holdout", ".git", ".venv", "venv", "__pycache__", "node_modules",
                            "agent-os", "raw", ".pytest_cache", "hexw.egg-info"})
_REF_EXT = frozenset({".py", ".md", ".txt", ".csv", ".tsv", ".json", ".jsonl", ".yaml", ".yml", ".in", ".out", ".log",
                      ".cif", ".mcif", ".ini", ".cfg", ".toml", ".sh", ".ps1", ".bat", ".xml", ".html", ".tex", ".rst"})
_REF_MAX_BYTES = 8_000_000
# literature/fit and literature/excluded are policed by the stricter fit rules (find_fit_references), not by this scan
_REF_FIT_DIRS = (("literature", "fit"), ("literature", "excluded"))
# third-party installations / downloads (QE build, pseudopotentials, package closure: thousands of directories)
_REF_THIRD_PARTY = (("tools", "_dl"), ("tools", "pseudo"), ("tools", "qe_ucrt64"))


def _list_files(p: Path):
    """Files of a source (audit frame: the guard allows hexw.audit to enumerate the holdout zone)."""
    if p.is_file():
        return [p]
    return sorted(q for q in p.rglob("*") if q.is_file()) if p.is_dir() else []


def _pinned_allowlist(cfg: paths.Config) -> dict:
    """{relative path: sha256} of the files that may name a holdout source ONLY while their content is byte-identical to the pin
    (`reference_allowlist_sha256` of holdout_sources.json): the search-pipeline scripts and the documents that cite the catalogue."""
    try:
        d = json.loads(cfg.holdout_config_file.read_text(encoding="utf-8")).get("reference_allowlist_sha256", {})
    except (OSError, ValueError):
        return {}
    return {str(k): str(v).lower() for k, v in d.items()} if isinstance(d, dict) else {}


def _file_sha256(q: Path) -> str:
    h = hashlib.sha256()
    with open(q, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def find_undeclared_holdout_references(cfg: paths.Config, allow=None, pinned=None):
    """Files outside vault/ and literature_holdout/ (raw/ is immutable data: its manifests list paths by design, and
    it is guarded by the hook) whose TEXT contains the path or the file name of a holdout source. Such a file must be
    declared holdout (it is a derivative: classified HOLDOUT) or listed in `reference_allowlist` of
    holdout_sources.json. -> (offenders, n_files_scanned, n_terms). Catches derivatives made in the future, e.g. a DFT
    pilot built from a holdout CIF."""
    import fnmatch
    from hexw.ingest.literature_split import holdout_reference_terms
    if allow is None:
        try:
            allow = json.loads(cfg.holdout_config_file.read_text(encoding="utf-8")).get("reference_allowlist", [])
        except (OSError, ValueError):
            allow = []
    if pinned is None:
        pinned = _pinned_allowlist(cfg)
    # whole-publication view: the folder of a mixed publication (L01) counts too
    terms = holdout_reference_terms(cfg.root, cfg.holdout_sources, mixed={}, listing=_list_files)
    offenders, n = [], 0
    protected = []                       # identities (st_dev, st_ino) of the protected directories / declared sources
    for q in (cfg.literature_holdout, cfg.vault, cfg.original_zip_dir, *cfg.holdout_sources):
        try:
            protected.append(os.stat(q))
        except OSError:
            pass

    def is_protected(p: Path) -> bool:   # by identity: a junction / symlink to a protected directory is pruned too
        try:
            st = os.stat(p)
        except OSError:
            return False
        return any(os.path.samestat(st, b) for b in protected)

    for d, dirs, files in os.walk(cfg.root):
        dp = Path(d)
        rel_d = dp.relative_to(cfg.root).parts
        kept = []
        for x in dirs:
            if x in _REF_SKIP_DIRS or tuple(rel_d + (x,)) in _REF_FIT_DIRS + _REF_THIRD_PARTY:
                continue
            if is_protected(dp / x):                                                     # declared: not scanned
                continue
            kept.append(x)
        dirs[:] = kept
        for f in files:
            q = dp / f
            if q.suffix.lower() not in _REF_EXT:
                continue
            rel = q.relative_to(cfg.root).as_posix()
            try:
                if q.stat().st_size > _REF_MAX_BYTES:
                    continue
                text = q.read_text(encoding="utf-8", errors="ignore").replace("\\", "/").lower()
            except OSError:
                continue
            n += 1
            if not any(t in text for t in terms):
                continue
            if paths.classify(q) is paths.Zone.HOLDOUT or any(fnmatch.fnmatch(rel, pat) for pat in allow):
                continue
            if rel in pinned:                      # allowed only while the content equals the pin
                try:
                    if _file_sha256(q) == pinned[rel]:
                        continue
                except OSError:
                    pass
                offenders.append(rel + " (sha256 differs from the pin)")
                continue
            offenders.append(rel)
    return offenders, n, len(terms)


def check_holdout_references(rep: Report, cfg: paths.Config):
    """Stage 3 round 1: derivatives of holdout sources must be declared; the fit side must not refer to holdout."""
    from hexw.ingest.literature_split import find_fit_references
    from hexw.params.literature import load_fit_records
    off, n, nt = find_undeclared_holdout_references(cfg)
    rep.add(f"HOLDOUT: no file outside vault/ and literature_holdout/ names a holdout source (path or file name) unless "
            f"it is declared holdout or allow-listed ({n} text files, {nt} terms)", not off, f"offenders={off}")
    pins = _pinned_allowlist(cfg)
    bad_pins = [rel for rel, h in pins.items()
                if not (cfg.root / rel).is_file() or _file_sha256(cfg.root / rel) != h]
    rep.add(f"HOLDOUT: every pinned reference-allowlist file ({len(pins)}) exists and matches its sha256", not bad_pins,
            f"mismatch={bad_pins}")
    deriv = [cfg.root / "dft_feasibility" / "pilots" / "c_W92", cfg.root / "dft_feasibility" / "pilots" / "make_bw.py"]
    bad = [str(q) for q in deriv if q.exists() and paths.classify(q) is not paths.Zone.HOLDOUT]
    rep.add("HOLDOUT: DFT pilot derivatives of the L01 Zn HRPT CIF (pilots/c_W92, make_bw.py) are declared holdout",
            not bad, str(bad))
    fit_texts = {}
    for sub in ("fit", "excluded"):
        base = cfg.literature / sub
        for q in sorted(base.rglob("*")) if base.is_dir() else []:
            if q.is_file() and q.suffix.lower() in _REF_EXT:
                fit_texts[f"{sub}/" + q.relative_to(base).as_posix()] = q.read_text(encoding="utf-8", errors="ignore")
    try:
        hits = find_fit_references(cfg.root, cfg.holdout_sources, load_fit_records(), fit_texts, listing=_list_files)
    except Exception as e:  # noqa
        hits = [f"check crashed: {e!r}"]
    rep.add("literature/fit, excluded: no path / file name of a holdout source, no holdout sample next to its "
            f"publication ({len(fit_texts)} files)", not hits, f"{hits[:5]}")


def check_guards(rep: Report, cfg: paths.Config):
    from hexw.chem.validators import (ChargeError, ModelNotApplicableError, SiteOverflowError,
                                      check_charge_neutral, check_site_occupancy, fe2_fe3_model)
    from hexw.registry.schemas import SiteOccupancy, StructureCandidate
    from hexw.registry.versioning import (AppendOnlyError, LeakageError, ValidationStore,
                                          check_inputs_targets_disjoint, make_model_version)
    any_vault = next(p for p in sorted(cfg.vault_legacy.rglob("*")) if p.is_file())

    src = "from hexw.io import paths\ndef go(p):\n    return paths.open_data(p, 'rb', purpose='metadata')\n"
    for mod in ("hexw.theory._probe", "hexw.params._probe"):
        s = _probe(mod, src)
        rep.add(f"{mod}: open_data on vault raises TheoryAccessError",
                _expect(paths.TheoryAccessError, lambda: s["go"](any_vault)))
    s = _probe("hexw.theory._probe", "def go(p):\n    return open(p, 'rb')\n")
    rep.add("hexw.theory: builtin open() on vault blocked by audit hook",
            _expect(paths.TheoryAccessError, lambda: s["go"](any_vault)))
    s = _probe("hexw.theory._probe", "def go(p):\n    return open(p, 'rb')\n")
    zp = next(cfg.original_zip_dir.glob("*.zip"))
    rep.add("hexw.theory: original zip (contains lab bytes) is blocked too",
            _expect(paths.TheoryAccessError, lambda: s["go"](zp)))
    if not paths.is_frozen():
        s = _probe("hexw.experiment._probe",
                   "from hexw.io import paths\ndef go(p):\n    return paths.open_data(p, 'rb')\n")
        rep.add("hexw.experiment: vault access without freeze file raises FreezeRequiredError",
                _expect(paths.FreezeRequiredError, lambda: s["go"](any_vault)))
    else:
        rep.add("freeze file present (freeze-denial check skipped)", True)
    import concurrent.futures as cf
    s = _probe("hexw.dft._probe", "def go(p):\n    return open(p, 'rb')\n")
    rep.add("hexw.dft: builtin open() on vault blocked", _expect(paths.TheoryAccessError, lambda: s["go"](any_vault)))
    with cf.ThreadPoolExecutor(1) as ex:
        rep.add("thread-pool worker without an allowed frame is refused",
                _expect(paths.AccessError, lambda: ex.submit(open, any_vault, "rb").result()))
    rep.add("extended-prefix spelling of a vault path is RESTRICTED",
            paths.classify(chr(92) * 2 + "?" + chr(92) + str(any_vault)) is paths.Zone.RESTRICTED)
    rep.add("write into vault is refused",
            _expect(paths.ReadOnlyZoneError, lambda: paths.open_data(cfg.vault_lab / "x.bin", "wb")))

    # overlap
    pe = Prediction("p1", "mv-1", [DataRef("OBS-1", sha256="aa"), DataRef("passport:x=1")])
    ok_vr = ValidationResult("r1", "mv-1", "p1", [DataRef("OBS-2", sha256="bb")], Metric("rel_err"))
    bad_vr = ValidationResult("r2", "mv-1", "p1", [DataRef("f.txt", sha256="cc", group="OBS-1")], Metric("rel_err"))
    check_inputs_targets_disjoint(pe, ok_vr)
    rep.add("inputs/targets overlap raises LeakageError", _expect(LeakageError,
            lambda: check_inputs_targets_disjoint(pe, bad_vr)))

    # versioning
    def P(v):
        return ParameterEstimate("Ms", v, "A/m", "magnetization", SourceRole.DFT, Provenance("dft-1"))
    v1, v2 = make_model_version([P(1.0)]), make_model_version([P(1.0000001)])
    rep.add("changing a ParameterEstimate changes the model version", v1.version_id != v2.version_id)
    with tempfile.TemporaryDirectory() as td:
        st = ValidationStore(Path(td) / "v.jsonl")
        r = ValidationResult("r1", v1.version_id, "p1", [DataRef("OBS-2")], Metric("e"))
        st.append(r,Prediction(r.prediction_id,r.model_version,[DataRef("audit-fit-input")]))
        before = (Path(td) / "v.jsonl").read_text(encoding="utf-8")
        ok1 = _expect(AppendOnlyError, lambda: st.append(r,Prediction(r.prediction_id,r.model_version,[DataRef("audit-fit-input")])))
        st.append(ValidationResult("r2", v2.version_id, "p1", [DataRef("OBS-2")], Metric("e")),
                  Prediction("p1",v2.version_id,[DataRef("audit-fit-input")]))
        after = (Path(td) / "v.jsonl").read_text(encoding="utf-8")
        rep.add("validation store is append-only; old results untouched",
                ok1 and after.startswith(before) and st.verify() and len(st.for_model(v1.version_id)) == 1)

    # chemistry
    over = StructureCandidate("c", "x", [SiteOccupancy("12k", {"Fe": 0.8, "Co": 0.4})])
    rep.add("site overflow detected", _expect(SiteOverflowError, lambda: check_site_occupancy(over)))
    rep.add("charge violation detected", _expect(ChargeError, lambda: check_charge_neutral(
        {"Sr": 1, "Fe": 18, "O": 27}, {"Sr": 2, "Fe": 2.0, "O": -2})))
    m = fe2_fe3_model(1.0, "Co")
    rep.add("Fe2+/Fe3+ model neutral for x=1 (n_Fe2=1, n_Fe3=16)", m.n_fe2 == 1.0 and m.n_fe3 == 16.0)
    rep.add("x > 2 rejected as not applicable", _expect(ModelNotApplicableError, lambda: fe2_fe3_model(2.1, "Co")))


# ---------------------------------------------------------------------------------------------- subprocess policy
THEORY_PACKAGES = ("theory", "params", "structures", "dft", "spin")
_SUBPROCESS_CALLS = {"subprocess": {"run", "call", "check_call", "check_output", "Popen", "getoutput", "getstatusoutput"},
                     "os": {"system", "popen", "execl", "execle", "execlp", "execlpe", "execv", "execve", "execvp",
                            "execvpe", "spawnl", "spawnle", "spawnlp", "spawnlpe", "spawnv", "spawnve", "spawnvp",
                            "spawnvpe", "posix_spawn", "posix_spawnp", "startfile"},
                     "multiprocessing": {"Process", "Pool", "get_context", "set_start_method"},
                     "concurrent.futures": {"ProcessPoolExecutor"},
                     "asyncio": {"create_subprocess_exec", "create_subprocess_shell"},
                     "pty": {"spawn", "fork"}, "pexpect": {"spawn", "run"}}
_WHOLE_MODULES = ("sh", "plumbum")            # any call through these is a process launch
_POLICY_MESSAGE = ("start Python children only through hexw.util.pool.WorkerPool and other programs only through "
                   "hexw.util.external.run_external (white list)")


def find_subprocess_violations(src: str, filename: str = "<src>") -> list:
    """Static analysis of one source: [(line, text)] of calls that start a child process (subprocess.*, os.system /
    popen / exec* / spawn*, Popen, multiprocessing.Process / Pool, ProcessPoolExecutor, ...) in any import spelling
    (`import subprocess as sp`, `from subprocess import Popen as P`, `from os import system`), plus dynamic spellings
    (`__import__('subprocess')`, `importlib.import_module('subprocess')`, `getattr(os, 'system')`). A call that carries
    `sys.executable` or a 'python' string is reported as such. Calls of hexw.util.pool / hexw.util.external are the
    allowed way and are not reported. A syntax error is a violation too (the file cannot be checked)."""
    import ast
    try:
        tree = ast.parse(src, filename)
    except SyntaxError as e:
        return [(e.lineno or 0, f"cannot parse: {e.msg}")]
    mod_alias: dict = {}            # local name -> full module name (`import x.y` binds x; `import x.y as z` binds z)
    name_alias: dict = {}           # local name -> (module, attr) for `from x import y [as z]`
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            for a in n.names:
                if a.asname:
                    mod_alias[a.asname] = a.name
                else:
                    mod_alias[a.name.split(".")[0]] = a.name.split(".")[0]
                    mod_alias[a.name] = a.name
        elif isinstance(n, ast.ImportFrom) and n.module:
            for a in n.names:
                name_alias[a.asname or a.name] = (n.module, a.name)

    def dotted(node):
        parts = []
        while isinstance(node, ast.Attribute):
            parts.append(node.attr)
            node = node.value
        if isinstance(node, ast.Name):
            parts.append(node.id)
            return ".".join(reversed(parts))
        return None

    def flagged(mod, attr):
        return mod in _WHOLE_MODULES or attr in _SUBPROCESS_CALLS.get(mod, ())

    def resolve_module(head):
        """Full module name a dotted prefix refers to, or None."""
        first, _, rest = head.partition(".")
        if first in mod_alias:
            return mod_alias[first] + ("." + rest if rest else "")
        if first in name_alias:                       # `from concurrent import futures` -> futures.X
            m, a = name_alias[first]
            return f"{m}.{a}" + ("." + rest if rest else "")
        return head

    def mentions_python(call):
        for sub in ast.walk(call):
            if isinstance(sub, ast.Attribute) and dotted(sub) == "sys.executable":
                return True
            if isinstance(sub, ast.Constant) and isinstance(sub.value, str) and \
                    re.search(r"(^|[\\/\s])python[\d.]*(\.exe)?($|\s)", sub.value, re.I):
                return True
        return False

    out = []
    for n in ast.walk(tree):
        if not isinstance(n, ast.Call):
            continue
        f = n.func
        hit = None
        d = dotted(f)
        if d and "." in d:
            head, _, attr = d.rpartition(".")
            mod = resolve_module(head)
            if flagged(mod, attr):
                hit = f"{mod}.{attr}"
        elif isinstance(f, ast.Name) and f.id in name_alias:
            mod, attr = name_alias[f.id]
            if flagged(mod, attr):
                hit = f"{mod}.{attr}"
        if hit is None and isinstance(f, ast.Name) and f.id == "__import__" and n.args and \
                isinstance(n.args[0], ast.Constant) and str(n.args[0].value).split(".")[0] in ("subprocess", "pexpect"):
            hit = f"__import__({n.args[0].value!r})"
        if hit is None and d and d.endswith("import_module") and n.args and isinstance(n.args[0], ast.Constant) \
                and str(n.args[0].value).split(".")[0] in ("subprocess", "pexpect"):
            hit = f"import_module({n.args[0].value!r})"
        if hit is None and isinstance(f, ast.Name) and f.id == "getattr" and len(n.args) >= 2 and \
                isinstance(n.args[1], ast.Constant) and isinstance(n.args[1].value, str):
            tgt = dotted(n.args[0])
            tmod = resolve_module(tgt) if tgt else None
            if tmod and flagged(tmod, n.args[1].value):
                hit = f"getattr({tmod}, {n.args[1].value!r})"
        if hit:
            out.append((n.lineno, f"{hit}{' with sys.executable / python' if mentions_python(n) else ''}"))
    return sorted(set(out))


def scan_theory_sources(pkg_root: Path, packages=THEORY_PACKAGES) -> list:
    """[(relative path, line, text)] violations in the .py files of the theory-side packages."""
    bad = []
    for pkg in packages:
        base = Path(pkg_root) / pkg
        if not base.is_dir():
            continue
        for p in sorted(base.rglob("*.py")):
            for line, text in find_subprocess_violations(p.read_text(encoding="utf-8"), str(p)):
                bad.append((str(p.relative_to(pkg_root)), line, text))
    return bad


def check_subprocess_policy(rep: Report, cfg: paths.Config):
    import hexw
    from hexw.util import external
    pkg_root = Path(hexw.__file__).resolve().parent
    bad = scan_theory_sources(pkg_root)
    rep.add("subprocess policy: no subprocess / os.system / Popen / direct multiprocessing in "
            f"hexw.{'/'.join(THEORY_PACKAGES)} ({_POLICY_MESSAGE})", not bad, f"{bad[:5]}")
    probes = [(find_subprocess_violations("import subprocess, sys\nsubprocess.run([sys.executable, 'x.py'])\n"), True),
              (find_subprocess_violations("from os import system as s\ns('python x.py')\n"), True),
              (find_subprocess_violations("import subprocess as sp\nsp.Popen(['python', '-c', 'pass'])\n"), True),
              (find_subprocess_violations("from hexw.util.pool import WorkerPool\nWorkerPool(2)\n"), False)]
    rep.add("subprocess policy: the static check recognises aliased / from-imports and ignores WorkerPool",
            all(bool(r) == want for r, want in probes), f"{probes}")
    ok = True
    for bad_exe in ("python", sys.executable, "cmd", "powershell", "sh", "bash"):
        try:
            external.resolve_executable(bad_exe)
            ok = False
        except external.ExternalProgramError:
            pass
    rep.add("subprocess policy: run_external refuses python / sys.executable / shells (white list: "
            f"{sorted(external.ALLOWED_EXECUTABLES)})", ok)


def package_versions() -> dict:
    out = {"python": sys.version.split()[0]}
    for p in PKGS:
        try:
            out[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            out[p] = None
    return out


def run_audit(root=None, write_manifest=False, out=sys.stdout) -> int:
    """By default nothing is written. write_manifest=True stores a RunManifest in registry/runs (git-ignored)."""
    cfg = paths.get_config() if root is None else paths.configure(root)
    rep = Report()
    started = datetime.now(timezone.utc)
    info = {}
    for fn, args in ((check_integrity, ()), (check_readonly, ()), (check_registry, ()), (check_umnik, ()),
               (check_registry_git, ()), (check_guards, ()), (check_literature, ()),
               (check_holdout_references, ()), (check_subprocess_policy, ())):
        try:
            r = fn(rep, cfg, *args)
            if isinstance(r, dict):
                info.update(r)
        except Exception as e:  # noqa
            rep.add(f"{fn.__name__} crashed", False, repr(e))
    width = max(len(c["name"]) for c in rep.checks)
    for c in rep.checks:
        print(f"[{'OK' if c['ok'] else 'FAIL'}] {c['name']}" + ("" if c["ok"] else f"  -> {c['detail']}"), file=out)
    n_bad = sum(not c["ok"] for c in rep.checks)
    print(f"\naudit: {len(rep.checks) - n_bad}/{len(rep.checks)} checks passed; "
          f"{info.get('n_files', '?')} files; frozen_model={paths.is_frozen()}", file=out)
    if write_manifest:
        run = RunManifest(
            run_id="audit-" + started.strftime("%Y%m%dT%H%M%SZ"), command="audit",
            started_utc=started.isoformat(), state=RunState.COMPLETED if rep.ok else RunState.FAILED,
            settings={"root": str(cfg.root)}, package_versions=package_versions(),
            resources={"cpu_threads": os.cpu_count()},
            checksums=info.get("registry_checksums", {}), git_head=info.get("git_head", {}),
            finished_utc=datetime.now(timezone.utc).isoformat(),
            notes=json.dumps(rep.checks, ensure_ascii=False))
        d = cfg.registry / "runs"
        d.mkdir(parents=True, exist_ok=True)
        (d / (run.run_id + ".json")).write_text(dump_json(run, indent=1), encoding="utf-8")
    return 0 if rep.ok else 1
