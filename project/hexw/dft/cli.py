"""`hexw dft <plan|generate|run|status|collect|fit-j>` (also `python -m hexw.dft ...`).

plan      structure -> J classes, flip configurations, conditioning (SVD), time / memory estimate; writes nothing
generate  the same, and writes validated pw.x inputs + design/plan/manifest into <runs_root>/<campaign>/<composition>/
run       starts pw.x for the unfinished configurations; REFUSES without --i-accept-runtime when the estimate exceeds the
          threshold of the config; resumable (finished jobs are never touched)
status    states of all configurations
collect   parse outputs, check that the magnetic configuration was preserved
fit-j     least-squares J, hold-out check, ParameterEstimate (source_role=dft)
(generate --soc writes non-collinear SOC inputs with status not_feasible_on_laptop; they are never run)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional

EXIT_OK, EXIT_ERROR, EXIT_CONFIRM = 0, 1, 2


def _parse_composition(items) -> dict:
    out = {}
    for it in items or []:
        for part in str(it).split(","):
            if not part.strip():
                continue
            k, _, v = part.partition("=")
            out[k.strip()] = float(v)
    return out


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="hexw dft", description="Quantum ESPRESSO DFT+U adapter: exchange constants J from collinear configurations")
    ap.add_argument("--config", help="config yaml (default: config/dft/campaign_v1.yaml)")
    sub = ap.add_subparsers(dest="sub", required=True)

    def target(p, need_structure: bool):
        p.add_argument("--system", choices=["w", "hematite"], default="w", help="w: SrFe(18-x)A(x)O27 (92 atoms); hematite: smoke system")
        p.add_argument("--composition", action="append", metavar="A=x[,B=y]", help="substituents per formula unit, e.g. Mn=2 (empty = SrFe18O27)")
        p.add_argument("--mn-valence", type=int, choices=[2, 3], help="Mn hypothesis (default 2)")
        p.add_argument("--fe2-scenario", help="Fe2+ scenario (default uniform)")
        p.add_argument("--config-index", type=int, default=0, help="which structural configuration of generate_configurations (default 0 = most probable)")
        p.add_argument("--n-rep", type=int, default=2, help="number of structural configurations requested from generate_configurations")
        p.add_argument("--target", help="name of an existing composition directory (overrides the options above)")
        if need_structure:
            p.add_argument("--max-pairs", type=int, help="keep only the N leading (shortest distance) J classes")
        p.add_argument("--parallel", type=int, help="processes at a time (<= config engine.max_parallel)")
        p.add_argument("--kpoints", metavar="N1,N2,N3", help="override the k mesh of the system (k test); the config hash changes accordingly")
        p.add_argument("--fe2-start", choices=["seeded", "none"], help="override dft.fe2_localization.mode (second start of the k test: none)")
        p.add_argument("--tag", help="suffix of the composition directory name (variants of one structure, k test)")

    for name in ("plan", "generate"):
        p = sub.add_parser(name)
        target(p, True)
        p.add_argument("--json", action="store_true", help="print the plan as JSON")
        if name == "generate":
            p.add_argument("--force", action="store_true", help="overwrite inputs of finished/running jobs")
            p.add_argument("--soc", action="store_true", help="also write the SOC inputs (status not_feasible_on_laptop)")
    p = sub.add_parser("run")
    target(p, False)
    p.add_argument("--i-accept-runtime", dest="accept", action="store_true", help="confirm the estimated runtime")
    p.add_argument("--retry-failed", action="store_true")
    p.add_argument("--force", action="store_true")
    p.add_argument("--only", action="append", help="run only these config_ids")
    p.add_argument("--timeout-hours", type=float)
    for name in ("status", "collect", "fit-j"):
        target(sub.add_parser(name), False)
    return ap


def _load(ns):
    from hexw.dft import settings as ST
    return ST.load_config(ns.config)


def _apply_overrides(ns, cfg):
    """--kpoints / --fe2-start: a copy of the config with these values (new hash)."""
    import json as _json
    if not getattr(ns, "kpoints", None) and not getattr(ns, "fe2_start", None):
        return cfg
    raw = _json.loads(_json.dumps(cfg.raw))
    if ns.kpoints:
        k = [int(x) for x in ns.kpoints.split(",")]
        if len(k) != 3 or min(k) < 1:
            raise SystemExit("--kpoints needs three positive integers, e.g. 2,2,1")
        kind = "hematite" if getattr(ns, "system", "w") == "hematite" else "W92"
        raw["systems"][kind]["kpoints"] = k
        if k != [1, 1, 1]:
            raw["systems"][kind]["gamma_tricks"] = False
    if ns.fe2_start:
        raw["dft"].setdefault("fe2_localization", {})["mode"] = ns.fe2_start
    return cfg.with_overrides(**{kk: raw[kk] for kk in ("systems", "dft")})


def _comp_dir_and_system(ns, cfg, need_system: bool):
    """-> (composition name, comp_dir, builder -> (system) or None)"""
    from hexw.dft import campaign as C
    if getattr(ns, "target", None):
        return ns.target, C.composition_dir(cfg, ns.target), None
    if ns.system == "hematite":
        name = "hematite_alpha-Fe2O3"

        def build():
            from hexw.dft.hematite import hematite_system
            return hematite_system()
        return name, C.composition_dir(cfg, name), build
    comp = _parse_composition(ns.composition)
    from hexw.structures.ensemble import Hypothesis
    hyp = Hypothesis(mn_valence=(ns.mn_valence or 2) if "Mn" in comp else None, fe2_scenario=ns.fe2_scenario or "uniform")
    name = C.composition_name(comp, hyp, ns.config_index)

    def build():
        from hexw.dft.magnet import system_from_w_config
        from hexw.structures.configs import generate_configurations
        cs = generate_configurations(comp, hyp, n_rep=max(ns.n_rep, ns.config_index + 1))
        if ns.config_index >= len(cs.configs):
            raise SystemExit(f"only {len(cs.configs)} structural configurations exist for this composition")
        return system_from_w_config(cs, ns.config_index, name)
    return name, C.composition_dir(cfg, name), build


def _print_plan(plan: dict, out=None) -> None:
    p = (out or sys.stdout).write
    p(f"system {plan['system']}  ({plan['n_atoms']} atoms, {plan['formula']})\n")
    p(f"convention: {plan['convention']}\n")
    s = plan["settings"]
    p(f"settings: PBE+U ({s['hubbard']['projector']}), U_eff {s['hubbard']['U_eff_eV']} (status {s['hubbard']['status']}), "
      f"pseudo profile {s['pseudo_profile']} {s['ecutwfc_Ry']:g}/{s['ecutrho_Ry']:g} Ry [{s['cutoff_status']}], config sha256 {s['config_sha256'][:12]}\n")
    p(f"J classes: {plan['n_J']} ({plan['group_by']}, r_cut {plan['r_cut_A']} A); parameters incl. E0: {plan['n_parameters_incl_E0']}; "
      f"fit configs {plan['n_fit']} (dof {plan['dof']}), hold-out {plan['n_holdout']}, total {plan['n_configs']}\n")
    for c in plan["J_classes"]:
        extra = f"  [{c['combination']}]" if c.get("combination") else ""
        ang = f", Fe-O-Fe {c['mean_Fe_O_Fe_angle_corner_deg']:.0f} deg" if c.get("mean_Fe_O_Fe_angle_corner_deg") else ""
        p(f"  {c['label']:<22} bonds/cell {c['n_bonds']:>3}  d {c['d_min']:.2f}-{c['d_max']:.2f} A  paths {c['paths']}{ang}{extra}\n")
    if plan["classes_not_in_the_fit"]:
        p("not in the fit: " + "; ".join(f"{a['label']} ({a['reason']})" for a in plan["classes_not_in_the_fit"]) + "\n")
    for n in plan["model_notes"]:
        p(f"note: {n}\n")
    if plan.get("nuisance_far_columns"):
        p("nuisance columns (sublattice pairs without a bond < r_cut; prior sigma %g meV per far bond, ASSUMPTION): %s\n" % (
            plan["far_prior_sigma_meV"], ", ".join(f"{c['label']} ({c['n_bonds_far']} bonds {c['d_min']:.2f}-{c['d_max']:.2f} A)"
                                                    for c in plan["nuisance_far_columns"])))
    c = plan["conditioning_fit_set"]
    sv = c["singular_values_scaled"]
    p(f"conditioning of the fit set [1, X] (columns scaled): cond = {c['cond_scaled']:.3g}, rank {c['rank']}/{c['n_columns']}, "
      f"SVD min/max = {min(sv):.3g}/{max(sv):.3g}; full pool of {plan['pool_size']} patterns: cond = {plan['conditioning_full_pool']['cond_scaled']:.3g}\n")
    p("configurations: " + ", ".join(f"{c['config_id']}[{c['role'][0]}]" for c in plan["configs"]) + "\n")
    r = plan["runtime"]
    p(f"runtime estimate ({r['basis']}):\n  {r['cpu_hours_per_job']:.2f} CPU-h per job, {r['cpu_hours_total']:.1f} CPU-h total, "
      f"{r['n_parallel']} in parallel -> ~{r['wall_hours_estimate']:.1f} h wall; RAM ~{r['ram_GB_peak_estimate']:.1f} GB (budget {r['memory_budget_GB']:g} GB, "
      f"{'ok' if r['memory_ok'] else 'EXCEEDED'}); timeout {r['timeout_hours']:g} h per job\n")
    sp = plan["settings_if_sssp_cutoffs"]
    p(f"with the SSSP-recommended cutoffs {sp['ecutwfc_Ry']:g}/{sp['ecutrho_Ry']:g} Ry ({sp['basis']}): {sp['status']}\n")


def main(argv=None) -> int:
    ap = build_parser()
    ns = ap.parse_args(argv)
    from hexw.dft import campaign as C
    from hexw.dft import runner as R
    from hexw.dft import settings as ST
    try:
        cfg = _load(ns)
    except ST.DftConfigError as e:
        print(f"hexw dft: config error: {e}", file=sys.stderr)
        return EXIT_ERROR
    try:
        cfg = _apply_overrides(ns, cfg)
    except (SystemExit, ST.DftConfigError) as e:
        print(f"hexw dft: {e}", file=sys.stderr)
        return EXIT_ERROR
    name, comp_dir, build = _comp_dir_and_system(ns, cfg, ns.sub in ("plan", "generate"))
    if getattr(ns, "tag", None) and not getattr(ns, "target", None):
        name = f"{name}_{ns.tag}"
        comp_dir = C.composition_dir(cfg, name)
    try:
        if ns.sub in ("plan", "generate"):
            system = build()
            design = C.make_design(cfg, system, getattr(ns, "max_pairs", None))
            plan = C.plan_summary(cfg, system, design, ns.parallel)
            if ns.sub == "plan":
                if ns.json:
                    print(json.dumps(plan, indent=1, default=float))
                else:
                    _print_plan(plan)
                return EXIT_OK
            res = C.generate(cfg, system, design, name, force=ns.force, n_parallel=ns.parallel)
            _print_plan(plan)
            print(f"\ngenerated in {res['dir']}: written {len(res['written'])}, unchanged {len(res['kept'])}; inputs validated")
            if ns.soc:
                from hexw.dft import soc
                s = soc.generate_soc_inputs(system, cfg, comp_dir)
                print(f"SOC inputs written: status {s['state']} ({s['reference']})")
            return EXIT_OK
        if not (comp_dir / "design.json").is_file():
            print(f"hexw dft: {comp_dir} has no design.json: run `generate` first", file=sys.stderr)
            return EXIT_ERROR
        if ns.sub == "status":
            for r in C.status_table(comp_dir):
                print(f"{r['config_id']:<34} {r['role']:<8} {r['state']:<22} attempts {r['attempts']} wall {r['wall_s']} it {r['n_iterations']} {r['message']}")
            return EXIT_OK
        if ns.sub == "collect":
            col = C.collect(cfg, comp_dir)
            for r in col["rows"]:
                e = r.get("energy_meV")
                print(f"{r['config_id']:<34} {r['status']:<22} E = {e if e is None else format(e, '.4f')} meV  {r.get('note', '')}")
            print(f"ok {col['n_ok']}, config_not_preserved {col['n_config_not_preserved']}, valence_changed {col['n_valence_changed']}, "
                  f"smearing_large {col['n_smearing_large']}, not_converged {col['n_not_converged']}")
            return EXIT_OK
        if ns.sub == "fit-j":
            fit = C.fit_j(cfg, comp_dir)
            if fit["status"] == "incomplete":
                print("fit-j: " + fit["message"], file=sys.stderr)
                return EXIT_ERROR
            print(f"convention: {fit['convention']}")
            for lab, j, s, s0 in zip(fit["class_labels"], fit["J_meV"], fit["sigma_J_meV"], fit["sigma_J_stat_meV"]):
                print(f"  J[{lab:<24}] = {j:+9.3f} +- {s:.3f} meV  (fit only: +- {s0:.3f})")
            print(f"dof {fit['dof']}, rms residual {fit['rms_residual_meV']:.3f} meV, sigma used {fit['sigma_used_meV']:.3f} meV ({fit['sigma_basis']}), "
                  f"cond {fit['conditioning']['cond_scaled']:.3g}")
            if fit.get("holdout"):
                h = fit["holdout"]
                print("hold-out: " + ", ".join(f"{c}: {e:+.3f} meV (z {z:+.2f})" for c, e, z in zip(h["config_ids"], h["error_meV"], h["z"])))
            print(f"hold-out check: {fit['holdout_check']}")
            for n in fit["notes"]:
                print("note:", n)
            return EXIT_OK
        if ns.sub == "run":
            dirs = C.job_dirs(comp_dir)
            if ns.only:
                dirs = [d for d in dirs if d.name in ns.only]
            dsg = json.loads((comp_dir / "design.json").read_text(encoding="utf-8"))
            kind = dsg["kind"]
            todo = [d for d in dirs if R.classify(cfg, d, ns.retry_failed, ns.force) == "run"]
            est = C.estimate_runtime(cfg, kind, len(todo), ns.parallel)
            print(f"{len(todo)} of {len(dirs)} jobs to run; estimated {est['cpu_hours_total']:.2f} CPU-h, ~{est['wall_hours_estimate']:.2f} h wall "
                  f"on {est['n_parallel']} processes, RAM ~{est['ram_GB_peak_estimate']:.1f} GB ({est['basis']})")
            if not est["memory_ok"]:
                print("hexw dft: estimated memory exceeds the budget; lower --parallel", file=sys.stderr)
                return EXIT_ERROR
            if est["cpu_hours_total"] > float(cfg.raw["runtime"]["confirm_above_cpu_hours"]) and not ns.accept:
                print("hexw dft run: the estimated runtime is above the threshold; confirm with --i-accept-runtime. Nothing was started.",
                      file=sys.stderr)
                return EXIT_CONFIRM
            tmo = (ns.timeout_hours or est["timeout_hours"]) * 3600.0
            out = R.run_queue(cfg, dirs, max_parallel=ns.parallel, timeout_s=tmo, retry_failed=ns.retry_failed, force=ns.force,
                              journal=R.Journal(comp_dir / "journal.jsonl"), progress=lambda s: print(s, flush=True))
            bad = [k for k, v in out.items() if v in ("failed", "timeout", "not_converged")]
            return EXIT_ERROR if bad else EXIT_OK
    except (C.CampaignError, SystemExit, ValueError) as e:
        print(f"hexw dft: {e}", file=sys.stderr)
        return EXIT_ERROR
    return EXIT_OK


if __name__ == "__main__":
    raise SystemExit(main())
