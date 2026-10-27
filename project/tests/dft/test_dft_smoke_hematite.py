"""SMOKE TEST of the whole chain with the REAL pw.x on alpha-Fe2O3 (10 atoms, ~1 min per configuration, 5 configurations).

Not part of the ordinary pytest run: `HEXW_RUN_DFT=1 python -m pytest tests/dft/test_dft_smoke_hematite.py -m dft -s`. It writes to the real
campaign directory (dft_runs/campaign_v1/hematite_alpha-Fe2O3) and resumes finished configurations.
generate -> run (3 processes, run_external) -> collect -> fit-j -> hold-out check. The asserts are about the chain and about
physics that does not depend on the settings: finished, converged, magnetic configurations preserved, the antiferromagnetic
ground state of hematite is the lowest of the configurations, resume does not rerun finished jobs, the hold-out is predicted
with a finite error that is reported. No J value is asserted: only its sign structure and a loose order of magnitude.
"""
import json

import numpy as np
import pytest

from hexw.dft import campaign as C
from hexw.dft import runner as R
from hexw.dft.hematite import hematite_system


@pytest.mark.dft
def test_hematite_chain_with_real_pw(cfg, tmp_path):
    """Works in the real campaign directory <HEXW_ROOT>/dft_runs/campaign_v1/hematite_alpha-Fe2O3 (checkpoint/resume: finished
    configurations are reused, so a repeated run only re-collects; delete the directory for a fresh run)."""
    if not cfg.executable_path().is_file():
        pytest.skip("pw_big.exe not installed")
    sysm = hematite_system()
    design = C.make_design(cfg, sysm)
    comp = "hematite_alpha-Fe2O3"
    C.generate(cfg, sysm, design, comp)
    d0 = C.composition_dir(cfg, comp)
    out = R.run_queue(cfg, C.job_dirs(d0), timeout_s=900, journal=R.Journal(d0 / "journal.jsonl"))
    assert set(out.values()) <= {"done", "skip:done"}, out                  # skip:done = resumed from a previous run
    col = C.collect(cfg, d0)
    assert col["n_ok"] == 5 and col["n_config_not_preserved"] == 0
    rows = {r["config_id"]: r for r in col["rows"]}
    ground = rows["c00_ground"]
    assert all(ground["energy_meV"] < r["energy_meV"] for k, r in rows.items() if k != "c00_ground")      # AFM ground state
    assert abs(ground["total_magnetization"]) < 0.05                                                      # compensated
    assert all(3.0 < abs(m) < 4.5 for m in ground["site_moments"][:4])                                    # high-spin Fe3+
    fit = C.fit_j(cfg, d0)
    assert fit["status"] in ("ok", "ok_with_notes") and fit["dof"] == 0 and fit["holdout"] is not None
    J = np.array(fit["J_meV"])
    assert np.all(np.isfinite(J)) and np.all(np.abs(J) < 500) and np.any(J < 0)                           # antiferromagnetic coupling present
    # resume: nothing is re-run
    assert set(R.run_queue(cfg, C.job_dirs(d0), timeout_s=5).values()) == {"skip:done"}
    h = fit["holdout"]
    print("\nSMOKE hematite: J =", [f"{lab}: {j:+.2f}+-{s:.2f} meV" for lab, j, s in zip(fit["class_labels"], J, fit["sigma_J_meV"])])
    print("  hold-out", h["config_ids"], "error", h["error_meV"], "meV, sigma_pred", h["sigma_pred_meV"], "z", h["z"])
    print("  residuals", fit["residuals_meV"], "energy range", fit["energy_range_meV"], "meV")
    print("  wall per job", {k: r["wall_s"] for k, r in rows.items()}, "s")
    (tmp_path / "smoke_result.json").write_text(json.dumps(fit, default=float, indent=1))
