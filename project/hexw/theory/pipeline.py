"""Composition -> joint material draws -> equilibria -> FMR/chi/mu, with resumable cases.

This is an ideal-crystal/equal-ground-domain MACROSPIN forecast. It does not silently
replace the missing atomistic spin dynamics, porosity or dielectric modules.
"""
from __future__ import annotations

import csv
import json
import os
import time
from collections import Counter
from pathlib import Path

import numpy as np

from hexw.io import paths
from hexw.params import anisotropy as A, moments as M
from hexw.params.literature import load_fit_records
from hexw.structures.ensemble import build_ensemble
from hexw.theory.request import ForecastRequest
from hexw.theory.magnet import (Gyro, MagnetModel, field_vector, find_equilibria, fwhm,
                               llg_frequency_Hz, susceptibility, passivity_margin,
                               sphere, thin_disk_perp_c, thin_disk_parallel_c, no_shape)
from hexw.workflow import atomic_json, canonical, core_fingerprint, digest, file_hash, plain, summary, utc

SHAPES = {"sphere": sphere, "disk_perp_c": thin_disk_perp_c,
          "disk_parallel_c": thin_disk_parallel_c, "no_shape": no_shape}


def _save_npz(path, **values):
    p = Path(path)
    tmp = p.with_name(p.name + ".tmp")
    with open(tmp, "wb") as f:
        np.savez_compressed(f, **values)
    os.replace(tmp, p)


def _csv(path, rows):
    if not rows:
        return
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _compute_material_draws(req: ForecastRequest):
    """Use the SAME joint draws in Ms, Tc and anisotropy; no median-parameter substitution."""
    d = A.joint_draws(req.n_draws, req.seed, req.n_prior, host=req.host,
                      co_form=req.co_form, l11_constraint=req.l11_constraint)
    an = A.anisotropy(req.composition, req.temperatures_K, draws=d)
    ms = M.ms_T(req.temperatures_K, req.composition, draws=d, seed=req.seed)
    ensemble = build_ensemble(req.composition, n_draws=req.n_draws, seed=req.seed)
    values = {"Ms_A_m": np.asarray(ms.A_per_m.samples), "Ms_A_m2_kg":np.asarray(ms.A_m2_per_kg.samples),
              "m_reduced": np.asarray(ms.m),
              "content_ratio": d.r, "Mn_valence": d.mn_valence,
              "Ms0_A_m": ms.ms0.A_per_m, "eta": ms.eta}
    for name, dist in zip(("K1", "K2", "K3", "K4"), an):
        if dist.samples is not None:
            values[name + "_J_m3"] = np.asarray(dist.samples)
    tc = np.asarray(A.C.tc_samples(req.composition, d.tc_params))
    values["Tc_K"] = tc
    records = load_fit_records()
    provenance = [{"record": c.rec_id, "name": c.name, "doi": c.source.doi,
                   "lit_id":c.lit_id,"measurement_group":c.source.measurement_group,
                   "locator": {k:getattr(c.source.locator,k) for k in ("page","table","figure")
                               if getattr(c.source.locator,k) is not None},
                   "full_locator_fit_record":c.rec_id,"locator_sha256":digest(plain(c.source.locator)),
                   "role": c.role.value,
                   "tier": c.source.tier, "status": c.source.data_status.value}
                  for c in records]
    meta = {"anisotropy_status": an.status, "anisotropy_notes": an.notes,
            "parameter_sources": provenance, "posterior_ess": d.posterior.ess,
            "moment_model": ms.A_per_m.notes, "structures": plain(ensemble.members),
            "structure_notes": ensemble.notes, "structure_ambiguity": ensemble.ambiguity,
            "scope": "ideal_crystal_macrospin_equal_ground_domain_average",
            "validity": "not independently validated; structural candidates use a literature template",
            "unimplemented": ["composition-specific atomistic exchange parameters/modes", "interacting/domain-wall dynamics",
                              "dielectric prediction from composition without independent dielectric parameters"],
            "optional_modules": ["coupled-sublattice LLG given explicit exchange matrix", "Debye/conductivity and spherical-pore mixing given explicit parameters"]}
    from hexw.params import sr_anisotropy as SR
    if req.anisotropy_model=="sr_single_crystal" or (req.anisotropy_model=="auto" and SR.supported(req.composition)):
        direct,details=SR.sample(req.composition,req.temperatures_K,d,values["Ms_A_m"],values["m_reduced"],seed=req.seed)
        values.update(direct)
        meta.update(anisotropy_status=details["status"],anisotropy_notes=details["note"],sr_single_crystal=details)
        meta["parameter_sources"].append(details["source"])
    return values, meta


def material_draws(req: ForecastRequest):
    """Content-addressed numerical cache; raw lab data and Python pickles are never used."""
    import inspect
    root = paths.get_config().root
    source = core_fingerprint(root)
    deps = {k: v for k, v in source.items() if k.startswith(("runtime/", "literature/fit/", "raw/")) or
            any(k.startswith("project/hexw/" + s) for s in ("params/", "structures/", "chem/", "units.py", "registry/schemas.py"))}
    settings = {k: getattr(req, k) for k in ("composition", "temperatures_K", "n_draws", "n_prior", "seed", "host", "co_form", "l11_constraint", "anisotropy_model")}
    key = digest({"settings": settings, "sources": deps, "implementation": inspect.getsource(_compute_material_draws)})
    cache = root / "runs" / "parameter_cache" / key
    seal = cache / "manifest.json"
    if seal.exists():
        saved = json.loads(seal.read_text(encoding="utf-8"))
        if any(file_hash(cache / k) != v for k, v in saved["files"].items()):
            raise ValueError("material parameter cache modified")
        with np.load(cache / "parameters.npz", allow_pickle=False) as a:
            values = {k: a[k] for k in a.files}
        return values, json.loads((cache / "material.json").read_text(encoding="utf-8"))
    values, meta = _compute_material_draws(req)
    cache.mkdir(parents=True, exist_ok=True)
    _save_npz(cache / "parameters.npz", **values)
    atomic_json(cache / "material.json", meta)
    atomic_json(seal, {"key": key, "settings": settings, "files": {
        k: file_hash(cache / k) for k in ("parameters.npz", "material.json")}})
    return values, meta


def compute_case(req, params, ti, hi):
    """All degenerate ground states are equally weighted (declared domain-population assumption)."""
    nf, nd = len(req.frequencies_GHz), req.n_draws
    frequencies = np.array(req.frequencies_GHz) * 1e9
    B = field_vector(req.fields_T[hi], np.deg2rad(req.theta_deg), np.deg2rad(req.phi_deg))
    gyro = Gyro.from_g(req.g)
    shape = SHAPES[req.shape]()
    chi_e = np.full((nd, nf, 3, 3), np.nan + 1j * np.nan)
    chi_i = chi_e.copy()
    lowest = np.full(nd, np.nan)
    modes, states = [], []
    failures = Counter()
    passivity = []
    for j in range(nd):
        Ms = float(params["Ms_A_m"][ti, j])
        if not np.isfinite(Ms) or Ms <= 1e-8 or params["m_reduced"][ti, j] <= 0:
            failures["outside_ferrimagnetic_phase"] += 1
            modes.append([])
            continue
        try:
            model = MagnetModel(Ms, *(float(params[f"K{k}_J_m3"][ti, j]) for k in range(1, 5)), N=shape)
            eqs = find_equilibria(model, B, n_grid=req.equilibrium_grid, verify=True)
            if eqs.ground is None or eqs.possibly_incomplete:
                failures["equilibrium_incomplete"] += 1
                modes.append([])
                continue
            ground = eqs.ground_degenerate or [eqs.ground]
            fs = [llg_frequency_Hz(model, B, e.m, gyro) for e in ground]
            if not np.isfinite(fs).all():
                raise ValueError("non-finite mode")
            ce = np.mean([susceptibility(model, B, e.m, frequencies, alpha=req.alpha,
                                          gyro=gyro, kind="external") for e in ground], axis=0)
            ci = np.mean([susceptibility(model, B, e.m, frequencies, alpha=req.alpha,
                                          gyro=gyro, kind="internal") for e in ground], axis=0)
            if not (np.isfinite(ce).all() and np.isfinite(ci).all()):
                raise ValueError("non-finite response")
            pm = passivity_margin(ce)
            if pm < -1e-7:
                raise ValueError("non-passive external response")
            chi_e[j], chi_i[j], lowest[j] = ce, ci, min(fs)
            modes.append(fs)
            states.extend(e.kind for e in ground)
            passivity.append(pm)
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as e:
            failures[type(e).__name__ + ":" + str(e)[:80]] += 1
            modes.append([])
    valid = np.isfinite(lowest)
    if not valid.any():
        return {"lowest_ground_mode_Hz": lowest, "chi_external": chi_e, "chi_internal": chi_i}, {
            "status": "not_predicted", "reason": "no numerically valid ferrimagnetic draws",
            "failures": dict(failures), "n_valid": 0, "n_total": nd}
    ce = chi_e[valid].mean(0)
    ci = chi_i[valid].mean(0)
    loss = ce[:, 0, 0].imag
    linewidth = fwhm(frequencies, loss, smooth=None, mode_x=float(np.median(lowest[valid])))
    peak_draws = np.full(nd, np.nan)
    width_draws = np.full(nd, np.nan)
    width_status = Counter()
    for j in np.flatnonzero(valid):
        curve = chi_e[j, :, 0, 0].imag
        # A maximum at the window boundary is a censored observation, not a resonance measurement.
        k = int(np.argmax(curve))
        if 0 < k < nf - 1:
            peak_draws[j] = frequencies[k]
        w = fwhm(frequencies, curve, smooth=None, mode_x=lowest[j])
        width_status[w.status] += 1
        if w.status == "valid":
            width_draws[j] = w.width
    result = {"status": "computed_unvalidated", "n_valid": int(valid.sum()), "n_total": nd,
              "valid_fraction": float(valid.mean()), "lowest_ground_mode_Hz": summary(lowest),
              "mode_definition": "lowest undamped eigenfrequency among degenerate ground states per draw",
              "all_ground_modes_Hz_by_draw": modes, "ground_state_counts": dict(Counter(states)),
              "failures": dict(failures), "external_passivity_min": min(passivity),
              "linewidth_of_mean_external_chi_xx": plain(linewidth),
              "sampled_absorption_peak_Hz": summary(peak_draws),
              "sampled_fwhm_Hz": summary(width_draws), "fwhm_status_counts": dict(width_status),
              "linewidth_note": "ensemble spectrum with specified alpha; not an inferred intrinsic damping or total sample linewidth",
              "statistics": "conditional on valid ferrimagnetic draws; invalid/above-Tc draws not filled with zero",
              "domain_population": "equal weights of all degenerate ground states; equilibrium, no magnetic-history fit"}
    values = {"lowest_ground_mode_Hz": lowest, "chi_external": chi_e, "chi_internal": chi_i,
              "mean_chi_external": ce, "mean_mu_internal": np.eye(3)[None] + ci,
              "peak_Hz": peak_draws, "fwhm_Hz": width_draws,
              "loss_xx_quantiles": np.quantile(chi_e[valid, :, 0, 0].imag, [.025, .5, .975], axis=0)}
    return values, result


def run_forecast(req: ForecastRequest, out=None, *, resume=True, progress=print):
    root = paths.get_config().root
    sources = core_fingerprint(root)
    ident = digest({"request": plain(req), "sources": sources})[:20]
    out = Path(out) if out is not None else root / "runs" / "predictions" / ident
    out = out.resolve()
    # Results cannot overwrite sources or the validation vault.
    if not out.is_relative_to((root / "runs").resolve()):
        raise ValueError("forecast outputs must stay below <project root>/runs")
    out.mkdir(parents=True, exist_ok=True)
    checkpoint = out / "checkpoint.json"
    if checkpoint.exists():
        state = json.loads(checkpoint.read_text(encoding="utf-8"))
        if state["run_id"] != ident:
            raise ValueError("output folder belongs to different inputs/code; choose a new run directory")
        if state["request"]!=plain(req) or json.loads((out/"request.json").read_text(encoding="utf-8"))!=plain(req):
            raise ValueError("input checkpoint modified")
        if not resume:
            raise FileExistsError("existing run is immutable; use resume or a different output directory")
        if state.get("status") == "complete":
            for name, h in state["artifacts"].items():
                if file_hash(out / name) != h:
                    raise ValueError(f"completed artifact changed: {name}")
            return out
    else:
        state = {"schema": "hexw.run.v1", "run_id": ident, "created_utc": utc(), "status": "running",
                 "source_fingerprint": sources, "request": plain(req), "completed_cases": {}, "artifacts": {}}
        atomic_json(out / "request.json", req)
        atomic_json(checkpoint, state)
    progress(f"run {ident}: parameter ensemble ({req.n_draws} draws)")
    if (out / "parameters.npz").is_file() and state.get("parameters_hash"):
        if file_hash(out / "parameters.npz") != state["parameters_hash"]:
            raise ValueError("parameter checkpoint modified")
        if file_hash(out/"material.json")!=state["material_hash"]:
            raise ValueError("material provenance checkpoint modified")
        with np.load(out / "parameters.npz", allow_pickle=False) as f:
            params = {k: f[k] for k in f.files}
        meta = json.loads((out / "material.json").read_text(encoding="utf-8"))
    else:
        params, meta = material_draws(req)
        _save_npz(out / "parameters.npz", **params)
        atomic_json(out / "material.json", meta)
        state["parameters_hash"] = file_hash(out / "parameters.npz")
        state["material_hash"] = file_hash(out/"material.json")
        atomic_json(checkpoint, state)
    if not all(f"K{k}_J_m3" in params for k in range(1, 5)):
        state["status"] = "unsupported_anisotropy"
        atomic_json(out / "summary.json", {"run_id": ident, "status": state["status"], "reason": meta["anisotropy_notes"],
                                           "request": plain(req), "available": {k: summary(v) for k, v in params.items()}})
        atomic_json(checkpoint, state)
        return out
    rows = []
    for ti, T in enumerate(req.temperatures_K):
        for hi, H in enumerate(req.fields_T):
            key = f"T{ti:03d}_H{hi:04d}"
            npz, js = out / f"{key}.npz", out / f"{key}.json"
            old = state["completed_cases"].get(key)
            if old:
                if file_hash(npz) != old["npz"] or file_hash(js) != old["json"]:
                    raise ValueError(f"case checkpoint modified: {key}")
                result = json.loads(js.read_text(encoding="utf-8"))
            else:
                progress(f"{key}: T={T:g} K, mu0 H={H:g} T")
                start = time.monotonic()
                data, result = compute_case(req, params, ti, hi)
                result.update(T_K=T, mu0H_T=H, wall_s=round(time.monotonic() - start, 3))
                _save_npz(npz, **data)
                atomic_json(js, result)
                state["completed_cases"][key] = {"npz": file_hash(npz), "json": file_hash(js)}
                state["updated_utc"] = utc()
                atomic_json(checkpoint, state)
            q = result.get("lowest_ground_mode_Hz", {})
            rows.append({"T_K": T, "mu0H_T": H, "status": result["status"], "n_valid": result["n_valid"],
                         "n_total": req.n_draws, "f_low_median_GHz": (q.get("median") / 1e9 if q.get("median") is not None else None),
                         "f_low_lo95_GHz": (q.get("lo95") / 1e9 if q.get("lo95") is not None else None),
                         "f_low_hi95_GHz": (q.get("hi95") / 1e9 if q.get("hi95") is not None else None)})
    _csv(out / "fmr_vs_field.csv", rows)
    parameters = [{"T_K": T, "Ms_A_m": summary(params["Ms_A_m"][i]),
                   "Ms_A_m2_kg":summary(params["Ms_A_m2_kg"][i]) if "Ms_A_m2_kg" in params else None,
                   **{f"K{k}_J_m3": summary(params[f"K{k}_J_m3"][i]) for k in range(1, 5)}}
                  for i, T in enumerate(req.temperatures_K)]
    summary_out = {"schema": "hexw.forecast.v1", "run_id": ident, "status": "complete_unvalidated",
                   "request": plain(req), "material_status": meta["anisotropy_status"],
                   "model_scope": meta["scope"], "Tc_K": summary(params["Tc_K"]),
                   "parameters": parameters, "field_sweep": rows,
                   "accuracy": {"target_relative_error": 0.1, "achieved": None, "independent_validation": "not_performed"},
                   "unimplemented": meta["unimplemented"], "created_utc": utc()}
    summary_out["electromagnetic_scope"]={"mu":"intrinsic internal-field tensor and external susceptibility computed separately",
        "epsilon":"user-specified passive dispersion" if req.dielectric_model else "not_predicted: dielectric parameters missing",
        "effective_medium":"isotropic trace approximation with spherical pores" if req.porosity_fraction is not None else "not_requested"}
    atomic_json(out / "summary.json", summary_out)
    export_structures(meta, out)
    export_spectra(req, out)
    make_plots(req, out, rows)
    state["status"] = "complete"
    state["finished_utc"] = utc()
    state["artifacts"] = {p.name: file_hash(p) for p in out.iterdir()
                          if p.is_file() and p.name not in {"checkpoint.json", "freeze.json"}}
    atomic_json(checkpoint, state)
    progress(f"complete: {out}")
    return out


def export_spectra(req, out):
    from hexw.theory.electromagnetic import dielectric, maxwell_garnett
    epsilon=None
    if req.dielectric_model:
        model={k:v for k,v in req.dielectric_model.items() if k!="provenance"}
        epsilon=dielectric(np.asarray(req.frequencies_GHz)*1e9,**model)
    for ti, T in enumerate(req.temperatures_K):
        for hi, H in enumerate(req.fields_T):
            key = f"T{ti:03d}_H{hi:04d}"
            with np.load(out / f"{key}.npz", allow_pickle=False) as a:
                if "mean_chi_external" not in a:
                    continue
                ce, mu = a["mean_chi_external"], a["mean_mu_internal"]
                rows = []
                for k, f in enumerate(req.frequencies_GHz):
                    row = {"f_GHz": f, "T_K": T, "mu0H_T": H}
                    if "loss_xx_quantiles" in a:
                        for q, val in zip(("lo95", "median", "hi95"), a["loss_xx_quantiles"][:, k]):
                            row[f"loss_xx_{q}"] = val
                    if epsilon is not None:
                        row.update(epsilon_real=epsilon[k].real,epsilon_imag=epsilon[k].imag)
                    if req.porosity_fraction is not None:
                        # Trace average is an EXPLICIT isotropic orientation approximation.
                        me=maxwell_garnett(np.trace(mu[k])/3,1,req.porosity_fraction)
                        row.update(mu_effective_isotropic_real=me.real,mu_effective_isotropic_imag=me.imag)
                        if epsilon is not None:
                            ee=maxwell_garnett(epsilon[k],1,req.porosity_fraction)
                            row.update(epsilon_effective_real=ee.real,epsilon_effective_imag=ee.imag)
                    for i in range(3):
                        for j in range(3):
                            for name, tensor in [("chi_external", ce), ("mu_internal", mu)]:
                                row[f"{name}_{i}{j}_real"] = tensor[k, i, j].real
                                row[f"{name}_{i}{j}_imag"] = tensor[k, i, j].imag
                    rows.append(row)
                _csv(out / f"{key}_spectrum.csv", rows)


def export_structures(meta, out):
    """Partial-occupancy CIF candidates on a cited template: explicitly NOT relaxed/DFT inputs."""
    if not meta.get("structures"):
        return
    from pymatgen.core import Species, Structure
    from pymatgen.io.cif import CifWriter
    from hexw.structures.configs import load_template
    t = load_template()
    for i, member in enumerate(meta["structures"]):
        c = member["candidate"]
        sites = {s.get("site_id") or s["site"]: s for s in c["sites"]}
        species, frac = [], []
        for si, xyz in zip(t.tm_site, t.tm_frac):
            s = sites[f"W{si+1}"]
            occ = {}
            for label, f in s["species"].items():
                if f <= 1e-12:
                    continue
                symbol = "Fe" if label.startswith("Fe") else label
                ox = c["oxidation_states"].get(label, 0)
                sp = Species(symbol, ox)
                occ[sp] = occ.get(sp, 0.) + f
            species.append(occ)
            frac.append(xyz)
        for sym, xyz in t.other_species:
            species.append({Species(sym, 2 if sym == "Sr" else -2): 1})
            frac.append(xyz)
        st = Structure(t.lattice, species, frac)
        text = "# CANDIDATE on L02 SrZn2W template; NOT relaxed; partial occupancy; NOT DFT input\n"
        text += str(CifWriter(st))
        (out / f"candidate_{i:02d}_template.cif").write_text(text, encoding="utf-8")
    atomic_json(out / "structural_hypotheses.json", {"template": t.source,
        "status": "unrelaxed_partial_occupancy_candidates", "hypotheses": meta["structures"],
        "notes": meta.get("structure_notes", [])})


def make_plots(req, out, rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for T in req.temperatures_K:
        rr = [r for r in rows if r["T_K"] == T and r["f_low_median_GHz"] is not None]
        if rr:
            x = [r["mu0H_T"] for r in rr]
            ax.plot(x, [r["f_low_median_GHz"] for r in rr], "o-", label=f"{T:g} K")
            ax.fill_between(x, [r["f_low_lo95_GHz"] for r in rr], [r["f_low_hi95_GHz"] for r in rr], alpha=.2)
    ax.set(xlabel="Applied field mu0 H (T)", ylabel="Lowest ground-state eigenfrequency (GHz)",
           title="Macrospin forecast; median and 95% scenario interval")
    if ax.lines:
        ax.legend()
    ax.grid(alpha=.25)
    fig.tight_layout()
    fig.savefig(out / "fmr_vs_field.png", dpi=160)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for hi, H in enumerate(req.fields_T):
        p = out / f"T000_H{hi:04d}_spectrum.csv"
        if p.exists():
            data = np.genfromtxt(p, delimiter=",", names=True, encoding="utf-8-sig")
            ax.plot(data["f_GHz"], data["chi_external_00_imag"], label=f"{H:g} T")
    ax.set(xlabel="Frequency (GHz)", ylabel="Im chi_external,xx (dimensionless)",
           title=f"Mean conditional spectrum, T={req.temperatures_K[0]:g} K; g={req.g:g}, alpha={req.alpha:g}")
    if ax.lines:
        ax.legend()
    ax.grid(alpha=.25)
    fig.tight_layout()
    fig.savefig(out / "fmr_spectra.png", dpi=160)
    plt.close(fig)
