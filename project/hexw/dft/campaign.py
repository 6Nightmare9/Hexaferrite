"""Campaign layer: plan -> generate -> (run) -> collect -> fit-j for one ordered structure.

Directory: <runs_root>/<campaign>/<composition>/{plan.json, design.json, manifest.json, journal.jsonl, collected.json,
fit_j.json, estimates.json, <config_id>/{pw.in, job.json, status.json, pw.out, result.json}}.

Nothing here invents numbers: energies come from parsed pw.x output, J from `hexw.dft.hamiltonian.fit_exchange`. 0 K energies are
used for exchange constants only; they are never converted to phase fractions (TZ section 8).
"""
from __future__ import annotations

import json
import math
import sys
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
from typing import Optional

import numpy as np

from hexw.dft import hamiltonian as H
from hexw.dft import magnet as M
from hexw.dft import parse as P
from hexw.dft import qe_input as Q
from hexw.dft import runner as R
from hexw.dft import settings as ST
from hexw.registry import schemas as S
from hexw.units import Kind


class CampaignError(RuntimeError):
    pass


# ---- names ------------------------------------------------------------------------------------------------------------
def composition_name(composition: dict, hypothesis=None, index: int = 0) -> str:
    x = sum(composition.values()) if composition else 0.0
    parts = "".join(f"{sp}{v:g}" for sp, v in sorted(composition.items()))
    name = f"Sr{parts}Fe{18 - x:g}O27"
    if hypothesis is not None and "Mn" in composition and getattr(hypothesis, "mn_valence", None):
        name += f"_Mn{hypothesis.mn_valence}"
    if hypothesis is not None and getattr(hypothesis, "fe2_scenario", None):
        name += f"_{hypothesis.fe2_scenario}"
    return f"{name}_cfg{index}"


def composition_dir(cfg: ST.DftConfig, name: str) -> Path:
    return cfg.runs_dir() / name


# ---- runtime estimate ---------------------------------------------------------------------------------------------------
def estimate_runtime(cfg: ST.DftConfig, kind: str, n_jobs: int, n_parallel: Optional[int] = None) -> dict:
    """Time / memory estimate from the FEASIBILITY pilots (config.systems[kind]); an estimate, not a measurement of this job."""
    s = cfg.system(kind)
    n_par = int(n_parallel or cfg.raw["engine"]["max_parallel"])
    cpu_h = float(s["seconds_per_scf_iteration"]) * float(s["scf_iterations"]) / 3600.0
    waves = math.ceil(n_jobs / n_par) if n_jobs else 0
    slow = float(cfg.raw["runtime"]["parallel_slowdown_factor"]) if n_par > 1 and n_jobs > 1 else 1.0
    ram = float(s["ram_GB"]) * min(n_par, max(n_jobs, 1))
    budget = float(cfg.raw["engine"].get("memory_budget_GB", 10.0))
    return {"system": kind, "n_jobs": n_jobs, "n_parallel": n_par, "cpu_hours_per_job": cpu_h,
            "cpu_hours_total": cpu_h * n_jobs, "waves": waves, "slowdown_factor": slow,
            "wall_hours_estimate": waves * cpu_h * slow, "ram_GB_peak_estimate": ram, "memory_budget_GB": budget,
            "memory_ok": ram <= budget,
            "basis": f"seconds/iteration {s['seconds_per_scf_iteration']} s x {s['scf_iterations']} iterations "
                     f"(config systems.{kind}; FEASIBILITY pilots); slowdown factor is an assumption",
            "timeout_hours": float(cfg.raw["runtime"]["timeout_hours"].get(kind, 12.0))}


# ---- plan -------------------------------------------------------------------------------------------------------------------
def make_design(cfg: ST.DftConfig, system: M.MagneticSystem, max_pairs: Optional[int] = None) -> M.Design:
    ex = cfg.exchange(system.kind)
    mp = max_pairs if max_pairs is not None else ex.get("max_pairs")
    model = M.build_model(system, float(ex["r_cut_A"]), ex["group_by"], mp, ex.get("far_r_max_A"))
    return M.select_design(model, int(ex["n_extra_fit"]), int(ex["n_holdout"]))


def plan_summary(cfg: ST.DftConfig, system: M.MagneticSystem, design: M.Design, n_parallel: Optional[int] = None) -> dict:
    m = design.model
    els = tuple(sorted({s.specie.symbol for s in system.structure}))
    prof = ST.load_profile(cfg, cfg.raw["pseudo_profile"], els)
    sssp = ST.load_profile(cfg, "sssp_efficiency", els)
    rt = estimate_runtime(cfg, system.kind, len(design.config_ids), n_parallel)
    return {
        "system": system.name, "kind": system.kind, "n_atoms": system.n_atoms,
        "formula": str(system.structure.composition.reduced_formula),
        "convention": H.CONVENTION_TEXT,
        "n_J": m.n_classes, "n_parameters_incl_E0": m.n_classes + 1,
        "n_fit": len(design.fit_rows), "n_holdout": len(design.holdout_rows), "n_configs": len(design.config_ids),
        "dof": len(design.fit_rows) - (m.n_classes + 1),
        "J_classes": [{k: v for k, v in c.items()} for c in m.class_info],
        "classes_not_in_the_fit": m.absorbed, "model_notes": m.notes + design.notes,
        "group_by": m.group_by, "r_cut_A": m.r_cut,
        "conditioning_fit_set": design.cond_fit, "conditioning_full_pool": design.pool_cond,
        "pool_size": int(len(m.pool_patterns)),
        "holdout_leverage": design.leverage_holdout,
        "nuisance_far_columns": design.model.far_info,
        "far_prior_sigma_meV": float(cfg.exchange(system.kind).get("far_prior_sigma_meV", 0.1)),
        "J_meaning": "effective J of the class = sum over ALL pairs of the two sublattices of J_ij S_i S_j / sum of S_i S_j of the "
                     "bonds shorter than r_cut (docs/DFT.md section 3)",
        "configs": [{"config_id": i, "role": r, "flipped_groups": f} for i, r, f in
                    zip(design.config_ids, design.roles, design.flips)],
        "runtime": rt,
        "settings": settings_record(cfg, prof),
        "settings_if_sssp_cutoffs": {"ecutwfc_Ry": sssp.ecutwfc, "ecutrho_Ry": sssp.ecutrho, "basis": sssp.cutoff_basis,
                                     "status": sssp.cutoff_status,
                                     "note": "cost grows roughly with the density grid (~(ecutrho ratio)^1.5, not measured): "
                                             "not feasible on the laptop"},
    }


def settings_record(cfg: ST.DftConfig, prof: ST.PseudoProfile) -> dict:
    d = cfg.raw["dft"]
    return {"config_sha256": cfg.sha256, "campaign": cfg.campaign, "functional": d["functional"],
            "hubbard": {"projector": d["hubbard"]["projector"], "U_eff_eV": d["hubbard"]["U_eff_eV"],
                        "status": d["hubbard"]["U_eff_status"]},
            "pseudo_profile": prof.name, "ecutwfc_Ry": prof.ecutwfc, "ecutrho_Ry": prof.ecutrho,
            "cutoff_status": prof.cutoff_status, "pseudo_files": prof.files,
            "occupations": d["occupations"], "scf": d["scf"], "starting_magnetization": d["starting_magnetization"],
            "fe2_localization": d.get("fe2_localization"),
            "literature": cfg.literature, "literature_deviations": cfg.raw["literature_basis"]["deviations_from_L05"]}


# ---- generate -----------------------------------------------------------------------------------------------------------------
def _atoms_record(system: M.MagneticSystem) -> list:
    out = []
    for i, s in enumerate(system.structure):
        try:
            ion = f"{s.specie.symbol}{int(round(s.specie.oxi_state))}+" if system.spin[i] > 0 else s.specie.symbol
        except Exception:                                            # noqa: BLE001 - no oxidation state
            ion = s.specie.symbol
        out.append({"i": i + 1, "species": s.specie.symbol, "ion": ion, "spin_S": float(system.spin[i]),
                    "group": int(system.group[i]), "base_sign": float(system.base_sign[i]),
                    "start_moment_muB": float(2 * system.spin[i])})
    return out


def _versions() -> dict:
    out = {"python": sys.version.split()[0]}
    for p in ("numpy", "scipy", "pymatgen", "PyYAML"):
        try:
            out[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            out[p] = None
    return out


def generate(cfg: ST.DftConfig, system: M.MagneticSystem, design: M.Design, comp: str, *, force: bool = False,
             n_parallel: Optional[int] = None) -> dict:
    """Write the input of every configuration (validated), design.json, plan.json and manifest.json. Existing inputs of
    finished / running jobs are never overwritten unless `force`."""
    d0 = composition_dir(cfg, comp)
    d0.mkdir(parents=True, exist_ok=True)
    elements = tuple(sorted({s.specie.symbol for s in system.structure}))
    prof = ST.load_profile(cfg, cfg.raw["pseudo_profile"], elements)
    journal = R.Journal(d0 / "journal.jsonl")
    plan = plan_summary(cfg, system, design, n_parallel)
    in_hash = {}
    written, kept = [], []
    for k, cid in enumerate(design.config_ids):
        dj = d0 / cid
        dj.mkdir(exist_ok=True)
        text = Q.build_input(system, design.patterns[k], cfg, prefix=f"p{k:02d}", config_id=cid, profile=prof)
        Q.assert_valid(text, cfg, system, design.patterns[k], prof)
        pin = dj / "pw.in"
        state = R.read_status(dj)["state"]
        if pin.is_file() and pin.read_text(encoding="utf-8") == text:
            kept.append(cid)
        elif pin.is_file() and state in ("running", "done") and not force:
            raise CampaignError(f"{cid}: input differs from the existing one and the job is {state}; use force")
        else:
            pin.write_text(text, encoding="utf-8")
            if state != "pending":
                (dj / "status.json").unlink(missing_ok=True)
            written.append(cid)
        in_hash[cid] = ST.sha256_file(pin)
        job = {"config_id": cid, "role": design.roles[k], "flipped_groups": design.flips[k],
               "pattern": [int(v) for v in design.patterns[k]], "system": system.name,
               "input_sha256": in_hash[cid], "config_sha256": cfg.sha256}
        (dj / "job.json").write_text(json.dumps(job, indent=1), encoding="utf-8")
        if not (dj / "status.json").is_file():
            R.write_status(dj, state="pending", config_sha256=cfg.sha256)
    dsg = {
        "system": system.name, "kind": system.kind, "meta": system.meta, "convention": H.CONVENTION_TEXT,
        "class_labels": design.model.class_labels, "class_info": design.model.class_info,
        "group_names": system.group_names, "config_ids": design.config_ids, "roles": design.roles,
        "flipped_groups": design.flips, "patterns": design.patterns.tolist(), "X": design.X.tolist(),
        "atoms": _atoms_record(system), "r_cut_A": design.model.r_cut, "group_by": design.model.group_by,
        "conditioning_fit_set": design.cond_fit, "config_sha256": cfg.sha256, "input_sha256": in_hash,
        "energy_noise_floor_meV": float(cfg.exchange(system.kind)["energy_noise_floor_meV"]),
        "max_condition_number": float(cfg.exchange(system.kind)["max_condition_number"]),
        "far_labels": design.model.far_labels, "far_info": design.model.far_info,
        "X_far": [] if design.X_far is None else design.X_far.tolist(),
        "far_prior_sigma_meV": float(cfg.exchange(system.kind).get("far_prior_sigma_meV", 0.1)),
        "holdout_rel_tol": float(cfg.exchange(system.kind).get("holdout_rel_tol", 0.02)),
        "max_promote_to_fit": int(cfg.exchange(system.kind).get("max_lost_fit_to_promote_holdout", 2)),
    }
    (d0 / "design.json").write_text(json.dumps(dsg, indent=1), encoding="utf-8")
    (d0 / "plan.json").write_text(json.dumps(plan, indent=1, default=float), encoding="utf-8")
    now = datetime.now(timezone.utc)
    man = S.RunManifest(
        run_id=f"dft-{cfg.campaign}-{comp}-{now.strftime('%Y%m%dT%H%M%SZ')}", command="dft generate",
        started_utc=now.isoformat(), state=S.RunState.COMPLETED,
        settings={"config_sha256": cfg.sha256, "campaign": cfg.campaign, "composition": comp, "kind": system.kind,
                  "pseudo_profile": prof.name, "ecutwfc_Ry": prof.ecutwfc, "ecutrho_Ry": prof.ecutrho,
                  "cutoff_status": prof.cutoff_status, "hubbard": cfg.raw["dft"]["hubbard"],
                  "input_sha256": in_hash, "qe_version_expected": cfg.raw["engine"]["version_expected"],
                  "convention": H.CONVENTION_TEXT},
        package_versions=_versions(), resources={"max_parallel": cfg.raw["engine"]["max_parallel"]},
        checksums={f: ST.sha256_file(prof.directory / f) for f in sorted(set(prof.files[e] for e in elements))},
        finished_utc=now.isoformat(),
        notes="checksums = sha256 of the pseudopotential files; stage: generate only, no pw.x was run")
    (d0 / "manifest.json").write_text(S.dump_json(man, indent=1), encoding="utf-8")
    journal.log("generate", composition=comp, written=written, kept=kept, config_sha256=cfg.sha256)
    return {"dir": str(d0), "written": written, "kept": kept, "plan": plan}


# ---- run helper -----------------------------------------------------------------------------------------------------------------
def job_dirs(comp_dir: Path) -> list:
    dsg = json.loads((Path(comp_dir) / "design.json").read_text(encoding="utf-8"))
    return [Path(comp_dir) / c for c in dsg["config_ids"]]


def status_table(comp_dir: Path) -> list:
    comp_dir = Path(comp_dir)
    dsg = json.loads((comp_dir / "design.json").read_text(encoding="utf-8"))
    rows = []
    for cid, role in zip(dsg["config_ids"], dsg["roles"]):
        st = R.read_status(comp_dir / cid)
        rows.append({"config_id": cid, "role": role, "state": st["state"], "attempts": st.get("attempts", 0),
                     "wall_s": st.get("wall_s"), "n_iterations": st.get("n_iterations"), "message": st.get("message", "")})
    return rows


# ---- collect ----------------------------------------------------------------------------------------------------------------------
def collect(cfg: ST.DftConfig, comp_dir: Path) -> dict:
    """Parse every finished job. Statuses: ok | config_not_preserved (a moment flipped / collapsed) | valence_changed (Fe2+ moved
    to other Fe atoms, total moment differs from sum(sign 2S) by more than 0.5 muB, Mn not 2+) | smearing_large (|TS| >= 1 meV) |
    not_converged / failed / timeout / pending (from the queue) | output_changed. Only `ok` goes into the fit."""
    comp_dir = Path(comp_dir)
    dsg = json.loads((comp_dir / "design.json").read_text(encoding="utf-8"))
    atoms = dsg["atoms"]
    spin = np.array([a["spin_S"] for a in atoms]); start = np.array([a["start_moment_muB"] for a in atoms])
    mc = cfg.raw["moments_check"]
    rows = []
    journal = R.Journal(comp_dir / "journal.jsonl")
    for k, (cid, role) in enumerate(zip(dsg["config_ids"], dsg["roles"])):
        dj = comp_dir / cid
        st = R.read_status(dj)
        row = {"config_id": cid, "role": role, "job_state": st["state"], "status": st["state"]}
        if st["state"] == "done":
            outp = dj / "pw.out"
            out_sha = ST.sha256_file(outp)
            if st.get("output_sha256") not in (None, out_sha):
                row.update(status="output_changed", note="pw.out differs from the recorded hash")
                rows.append(row)
                continue
            res = P.parse_pw_file(outp)
            pat = np.array(json.loads((dj / "job.json").read_text(encoding="utf-8"))["pattern"], float)
            sg = np.array([a["base_sign"] for a in atoms]) * np.array([pat[a["group"]] if a["group"] >= 0 else 0.0 for a in atoms])
            chk = P.check_configuration(res, sg, spin, min_fraction=float(mc["min_fraction_of_start"]), start_moment=start,
                                        require_sign=bool(mc["require_sign"]))
            val = P.check_valence(res, atoms, sg, cfg.raw.get("valence_check"))
            row.update(energy_Ry=res.energy_Ry, energy_meV=res.energy_meV, n_iterations=res.n_iterations,
                       hubbard_ns=res.hubbard_ns, hubbard_ns_start=res.hubbard_ns_start, valence_check=val,
                       total_magnetization=res.total_magnetization, absolute_magnetization=res.absolute_magnetization,
                       site_moments=res.site_moments, hubbard_d_moments=res.hubbard_d_moments,
                       smearing_contrib_Ry=res.smearing_contrib_Ry, wall_s=res.wall_s, cpu_s=res.cpu_s,
                       ram_estimate_MB=res.ram_estimate_MB, qe_version=res.qe_version, config_check=chk,
                       input_sha256=st.get("input_sha256"), output_sha256=out_sha, config_sha256=st.get("config_sha256"),
                       fermi_eV=res.fermi_eV, estimated_scf_accuracy_Ry=res.estimated_scf_accuracy_Ry)
            if not chk["ok"]:
                row["status"], row["note"] = "config_not_preserved", chk["reason"]
            elif not val["valence_ok"]:
                row["status"], row["note"] = "valence_changed", "; ".join(val["reasons"])
            elif not val["smearing_ok"]:
                row["status"], row["note"] = "smearing_large", val["smearing_reason"]
            else:
                row["status"] = "ok"
        (dj / "result.json").write_text(json.dumps(row, indent=1, default=float), encoding="utf-8")
        rows.append(row)
    out = {"composition_dir": str(comp_dir), "config_sha256": cfg.sha256, "rows": rows,
           "n_ok": sum(r["status"] == "ok" for r in rows),
           "n_config_not_preserved": sum(r["status"] == "config_not_preserved" for r in rows),
           "n_valence_changed": sum(r["status"] == "valence_changed" for r in rows),
           "n_smearing_large": sum(r["status"] == "smearing_large" for r in rows),
           "n_not_converged": sum(r["status"] == "not_converged" for r in rows)}
    (comp_dir / "collected.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    journal.log("collect", n_ok=out["n_ok"], n_config_not_preserved=out["n_config_not_preserved"],
                n_valence_changed=out["n_valence_changed"], n_smearing_large=out["n_smearing_large"])
    return out


# ---- fit ----------------------------------------------------------------------------------------------------------------------------
def fit_j(cfg: ST.DftConfig, comp_dir: Path, write: bool = True) -> dict:
    """OLS fit of the J classes on the configurations with status ok and role fit; hold-out configurations only predict."""
    comp_dir = Path(comp_dir)
    dsg = json.loads((comp_dir / "design.json").read_text(encoding="utf-8"))
    col = collect(cfg, comp_dir)
    rows = {r["config_id"]: r for r in col["rows"]}
    X_all = np.array(dsg["X"])
    names = dsg["class_labels"]
    p = len(names) + 1
    fit_ids = [c for c, r in zip(dsg["config_ids"], dsg["roles"]) if r == "fit" and rows[c]["status"] == "ok"]
    hold_ids = [c for c, r in zip(dsg["config_ids"], dsg["roles"]) if r == "holdout" and rows[c]["status"] == "ok"]
    excluded = {c: rows[c]["status"] for c in dsg["config_ids"] if rows[c]["status"] != "ok"}
    # If 1-2 planned fit configurations are lost, the hold-out configurations that finished are PROMOTED into the fit (with an
    # explicit mark) and the hold-out check is then unavailable: no hold-out configuration remains.
    n_fit_planned = list(dsg["roles"]).count("fit")
    lost = n_fit_planned - len(fit_ids)
    promoted, holdout_check = [], ("available" if hold_ids else "none_planned_or_finished")
    if 0 < lost <= int(dsg.get("max_promote_to_fit", 2)) and hold_ids:
        order_ix = {c: k for k, c in enumerate(dsg["config_ids"])}
        promoted = list(hold_ids)
        fit_ids = sorted(fit_ids + promoted, key=order_ix.get)
        hold_ids = []
        holdout_check = (f"unavailable: {lost} fit configuration(s) lost, hold-out configuration(s) {promoted} promoted into the "
                         "fit; the model was NOT tested on independent configurations")
    out = {"composition_dir": str(comp_dir), "convention": H.CONVENTION_TEXT, "class_labels": names,
           "excluded_configurations": excluded, "config_sha256": cfg.sha256, "holdout_promoted_to_fit": promoted,
           "holdout_check": holdout_check}
    if len(fit_ids) < p:
        out.update(status="incomplete", message=f"{len(fit_ids)} usable fit configurations for {p} parameters "
                                                f"(missing / not preserved: {excluded})")
        if write:
            (comp_dir / "fit_j.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
        return out
    idx = {c: k for k, c in enumerate(dsg["config_ids"])}
    E_fit = np.array([rows[c]["energy_meV"] for c in fit_ids])
    E_ref = float(E_fit.mean())
    Xf = X_all[[idx[c] for c in fit_ids]]
    Xh = X_all[[idx[c] for c in hold_ids]] if hold_ids else None
    Eh = np.array([rows[c]["energy_meV"] for c in hold_ids]) - E_ref if hold_ids else None
    Xfar = np.array(dsg.get("X_far") or []).reshape(len(X_all), -1)
    nuis = {}
    if Xfar.shape[1]:
        nuis = dict(X_nuis=Xfar[[idx[c] for c in fit_ids]], nuis_names=dsg["far_labels"],
                    nuis_prior_sigma_meV=float(dsg["far_prior_sigma_meV"]),
                    X_nuis_holdout=Xfar[[idx[c] for c in hold_ids]] if hold_ids else None)
    fr = H.fit_exchange(Xf, E_fit - E_ref, names, noise_floor_meV=float(dsg["energy_noise_floor_meV"]), X_holdout=Xh,
                        E_holdout_meV=Eh, max_cond=float(dsg["max_condition_number"]),
                        holdout_rel_tol=float(dsg.get("holdout_rel_tol", 0.02)), **nuis)
    inputs = {c: rows[c]["input_sha256"] for c in fit_ids + hold_ids}
    outputs = {c: rows[c]["output_sha256"] for c in fit_ids + hold_ids}
    prov_hash = ST.sha256_bytes(ST.canonical_json({"in": inputs, "out": outputs, "cfg": cfg.sha256}).encode())
    qe_versions = sorted({rows[c]["qe_version"] for c in fit_ids + hold_ids if rows[c].get("qe_version")})
    out.update(
        status="ok" if not fr.notes else "ok_with_notes", notes=fr.notes,
        E0_meV_relative_to_mean_fit_energy=fr.E0, E_ref_meV=E_ref, J_meV=fr.J.tolist(), sigma_J_stat_meV=fr.sigma_J.tolist(), sigma_J_meV=fr.sigma_J_total.tolist(),
        holdout_consistent=fr.holdout_consistent, sigma_model_meV=fr.sigma_model_meV,
        sigma_basis=fr.sigma_basis, sigma_used_meV=fr.sigma_used_meV, dof=fr.dof, rms_residual_meV=fr.rms_meV,
        s_resid_meV=None if math.isnan(fr.s_resid_meV) else fr.s_resid_meV,
        residuals_meV=dict(zip(fit_ids, fr.residuals.tolist())),
        loo_error_meV=dict(zip(fit_ids, [None if math.isnan(v) else v for v in fr.loo_residuals.tolist()])),
        holdout=None if fr.holdout is None else dict(fr.holdout, config_ids=hold_ids),
        max_condition_number=float(dsg["max_condition_number"]), conditioning=fr.cond, fit_configurations=fit_ids, holdout_configurations=hold_ids,
        nuisance=fr.nuisance, energy_range_meV=float(np.ptp(E_fit)), J_in_kelvin=[H.J_to_kelvin(v) for v in fr.J],
        class_info=dsg["class_info"], inputs_sha256=inputs, outputs_sha256=outputs, provenance_sha256=prov_hash,
        qe_versions=qe_versions)
    ests = parameter_estimates(cfg, comp_dir.name, out, fr)
    out["n_estimates"] = len(ests)
    if write:
        (comp_dir / "fit_j.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
        (comp_dir / "estimates.json").write_text("[" + ",\n".join(S.dump_json(e) for e in ests) + "]", encoding="utf-8")
    return out


def parameter_estimates(cfg: ST.DftConfig, comp: str, fit: dict, fr: H.FitResult) -> list:
    ok_cond = fr.cond["cond_scaled"] <= float(fit["max_condition_number"])
    ident = bool(ok_cond and fit["holdout_consistent"] is not False)          # a failed hold-out check: not identifiable
    ests = []
    for k, name in enumerate(fit["class_labels"]):
        info = fit["class_info"][k]
        note = {"role": "dft", "campaign": cfg.campaign, "composition": comp, "convention": H.CONVENTION_TEXT,
                "config_sha256": cfg.sha256, "inputs_outputs_sha256": fit["provenance_sha256"],
                "qe_versions": fit["qe_versions"], "sigma_basis": fit["sigma_basis"], "dof": fit["dof"], "holdout_consistent": fit["holdout_consistent"],
                "sigma_model_meV": fit["sigma_model_meV"], "holdout_check": fit["holdout_check"],
                "holdout_promoted_to_fit": fit["holdout_promoted_to_fit"],
                "nuisance_far_pairs": None if fit.get("nuisance") is None else
                {k: fit["nuisance"][k] for k in ("names", "prior_sigma_meV")},
                "J_meaning": "effective J of the class: sum over ALL pairs of the two sublattices of J_ij S_i S_j divided by the sum "
                             "of S_i S_j over the bonds shorter than r_cut (docs/DFT.md section 3)",
                "cutoff_status": cfg.raw["pseudo_profiles"][cfg.raw["pseudo_profile"]]["cutoff_status"],
                "U_status": cfg.raw["dft"]["hubbard"]["U_eff_status"],
                "temperature": "0 K, collinear, no relaxation", "combination": info.get("combination")}
        prov = S.Provenance(source_id=f"dft:{cfg.campaign}:{comp}", locator=str(Path(fit["composition_dir"]) / "fit_j.json"),
                            note=json.dumps(note, sort_keys=True), status="computed")
        unc = S.Uncertainty(kind="std", std=float(fr.sigma_J_total[k])) if ident else S.Uncertainty(kind="unknown")
        ests.append(S.ParameterEstimate(name=f"dft_J_{name}", value=float(fr.J[k]), unit="meV", quantity_kind=Kind.ENERGY,
                                        source_role=S.SourceRole.DFT, provenance=prov, uncertainty=unc,
                                        identifiable=ident))
    return ests
