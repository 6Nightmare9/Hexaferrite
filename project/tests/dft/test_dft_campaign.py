"""generate -> run -> collect -> fit-j with a FAKE pw.x that returns energies of a known Heisenberg model (no QE needed)."""
import json
import subprocess

import numpy as np
import pytest

from hexw.dft import campaign as C
from hexw.dft import hamiltonian as H
from hexw.dft import runner as R
from hexw.dft.hematite import hematite_system
from hexw.registry import schemas as S
from _dftsynth import synthetic_pw_out

J_TRUE = np.array([-12.0, 3.0, -5.0])
E_BASE = -1000.0 * H.RY_TO_MEV / 1000.0          # arbitrary offset, meV


@pytest.fixture(scope="module")
def hem_design():
    from hexw.dft import settings as ST
    cfg = ST.load_config()
    sysm = hematite_system()
    return sysm, C.make_design(cfg, sysm)


def make_fake(comp_dir, flip_in=None, unconverged=(), total_mag_shift=None, ts_Ry=0.0, e_shift=None, ns_of=None, mom_of=None):
    """Fake pw.x: writes pw.out into stdout_path like the real run_external does. total_mag_shift / ts_Ry / ns_of / mom_of
    {config_id: ...} damage the output of single configurations."""
    dsg = json.loads((comp_dir / "design.json").read_text())
    X = np.array(dsg["X"])
    atoms = dsg["atoms"]

    def fake(exe, args=(), *, cwd=None, timeout=1, env=None, stdin_path=None, stdout_path=None, stderr_path=None, pid_callback=None):
        k = dsg["config_ids"].index(cwd.name)
        e_mev = E_BASE + 10.0 + X[k] @ J_TRUE[: X.shape[1]] + (e_shift or {}).get(cwd.name, 0.0)
        pat = np.array(json.loads((cwd / "job.json").read_text())["pattern"], float)
        sg = [a["base_sign"] * (pat[a["group"]] if a["group"] >= 0 else 0.0) for a in atoms]
        mom = [x * 3.8 for x in sg]
        if cwd.name in ((flip_in,) if isinstance(flip_in, str) else tuple(flip_in or ())):
            mom[0] = -mom[0]                        # one Fe turns over
        if mom_of and cwd.name in mom_of:
            mom = mom_of[cwd.name](mom)
        m_tot = sum(x * 2 * a["spin_S"] for x, a in zip(sg, atoms)) + (total_mag_shift or {}).get(cwd.name, 0.0)
        text = synthetic_pw_out(e_mev / H.RY_TO_MEV, mom, converged=cwd.name not in unconverged, total_mag=m_tot,
                                ts_Ry=(ts_Ry if isinstance(ts_Ry, float) else ts_Ry.get(cwd.name, 0.0)),
                                ns=(ns_of or {}).get(cwd.name))
        assert stdout_path is not None and pid_callback is not None
        pid_callback(4242)
        with open(stdout_path, "w", encoding="utf-8") as f:
            f.write(text)
        return subprocess.CompletedProcess([], 0, None, None)
    return fake


def test_full_chain_recovers_j_and_checks_holdout(tmp_cfg, hem_design, monkeypatch):
    cfg = tmp_cfg
    sysm, design = hem_design
    comp = "hematite_test"
    res = C.generate(cfg, sysm, design, comp)
    d0 = C.composition_dir(cfg, comp)
    assert len(res["written"]) == 5 and (d0 / "design.json").is_file() and (d0 / "manifest.json").is_file()
    man = json.loads((d0 / "manifest.json").read_text())
    assert man["settings"]["config_sha256"] == cfg.sha256 and man["command"] == "dft generate"       # config hash in the RunManifest
    assert man["settings"]["cutoff_status"] == "pilot" and len(man["checksums"]) == 2
    # a second generate is idempotent
    assert C.generate(cfg, sysm, design, comp)["written"] == []
    monkeypatch.setattr(R, "run_external", make_fake(d0))
    out = R.run_queue(cfg, C.job_dirs(d0), timeout_s=10)
    assert set(out.values()) == {"done"}
    col = C.collect(cfg, d0)
    assert col["n_ok"] == 5 and col["n_config_not_preserved"] == 0, [(r["config_id"], r["status"], r.get("note")) for r in col["rows"]]
    fit = C.fit_j(cfg, d0)
    assert fit["status"] == "ok_with_notes" and fit["dof"] == 0 and fit["sigma_basis"] == "noise_floor"
    assert fit["J_meV"] == pytest.approx(list(J_TRUE), abs=1e-3)
    assert max(abs(e) for e in fit["holdout"]["error_meV"]) < 1e-3 and len(fit["holdout"]["config_ids"]) == 1
    assert fit["sigma_basis"] == "noise_floor" and all(s > 0 for s in fit["sigma_J_meV"])
    est = json.loads((d0 / "estimates.json").read_text())
    assert len(est) == 3 and all(e["source_role"] == "dft" and e["unit"] == "meV" for e in est)
    e0 = S.from_dict(S.ParameterEstimate, est[0])
    prov = json.loads(e0.provenance.note)
    assert prov["config_sha256"] == cfg.sha256 and prov["qe_versions"] == ["7.6"] and len(prov["inputs_outputs_sha256"]) == 64
    assert prov["cutoff_status"] == "pilot" and prov["U_status"]["Co"] == "assumption" and "0 K" in prov["temperature"]
    # status table and resume
    assert {r["state"] for r in C.status_table(d0)} == {"done"}
    assert set(R.run_queue(cfg, C.job_dirs(d0), timeout_s=10).values()) == {"skip:done"}


def test_flipped_configuration_is_excluded_and_the_holdout_is_promoted(tmp_cfg, hem_design, monkeypatch):
    """One fit configuration lost (<= 2): the finished hold-out is promoted into the fit, explicitly marked; no hold-out check."""
    cfg = tmp_cfg
    sysm, design = hem_design
    comp = "hematite_flip"
    C.generate(cfg, sysm, design, comp)
    d0 = C.composition_dir(cfg, comp)
    victim = design.config_ids[2]
    monkeypatch.setattr(R, "run_external", make_fake(d0, flip_in=victim))
    R.run_queue(cfg, C.job_dirs(d0), timeout_s=10)
    col = C.collect(cfg, d0)
    row = next(r for r in col["rows"] if r["config_id"] == victim)
    assert row["status"] == "config_not_preserved" and row["config_check"]["flipped"] == [1]
    assert "config_not_preserved" in (d0 / victim / "result.json").read_text()
    fit = C.fit_j(cfg, d0)
    hold = design.config_ids[4]
    assert victim in fit["excluded_configurations"] and fit["holdout_promoted_to_fit"] == [hold]
    assert hold in fit["fit_configurations"] and fit["holdout"] is None and fit["holdout_consistent"] is None
    assert fit["holdout_check"].startswith("unavailable") and fit["status"] == "ok_with_notes"
    est = json.loads((d0 / "estimates.json").read_text())
    note = json.loads(est[0]["provenance"]["note"])
    assert note["holdout_promoted_to_fit"] == [hold] and note["holdout_check"].startswith("unavailable")


def test_more_than_two_lost_fit_configurations_is_incomplete(tmp_cfg, hem_design, monkeypatch):
    cfg = tmp_cfg
    sysm, design = hem_design
    comp = "hematite_flip3"
    C.generate(cfg, sysm, design, comp)
    d0 = C.composition_dir(cfg, comp)
    victims = design.config_ids[1:4]
    monkeypatch.setattr(R, "run_external", make_fake(d0, flip_in=victims))
    R.run_queue(cfg, C.job_dirs(d0), timeout_s=10)
    fit = C.fit_j(cfg, d0)
    assert fit["status"] == "incomplete" and set(victims) <= set(fit["excluded_configurations"])      # 3 lost: no promotion
    assert fit["holdout_promoted_to_fit"] == [] and "estimates.json" not in {p.name for p in d0.iterdir()}


def test_extra_fit_point_lets_the_fit_survive_one_bad_configuration(tmp_cfg, hem_design, monkeypatch):
    """With one spare configuration (n_extra_fit = 1) a flipped run is dropped; the hold-out is promoted, the fit keeps dof = 1."""
    from hexw.dft import magnet as M
    cfg = tmp_cfg
    sysm, _ = hem_design
    model = M.build_model(sysm, 3.5, "shell")
    design = M.select_design(model, 1, 1)                                                   # 2 J + E0, 4 fit rows
    comp = "hematite_spare"
    C.generate(cfg, sysm, design, comp)
    d0 = C.composition_dir(cfg, comp)
    victim = design.config_ids[2]
    monkeypatch.setattr(R, "run_external", make_fake(d0, flip_in=victim))
    R.run_queue(cfg, C.job_dirs(d0), timeout_s=10)
    fit = C.fit_j(cfg, d0)
    assert victim in fit["excluded_configurations"] and victim not in fit["fit_configurations"]
    assert fit["holdout_promoted_to_fit"] == [design.config_ids[-1]] and fit["dof"] == 1


def test_unconverged_run_never_enters_the_fit(tmp_cfg, hem_design, monkeypatch):
    cfg = tmp_cfg
    sysm, design = hem_design
    comp = "hematite_unconv"
    C.generate(cfg, sysm, design, comp)
    d0 = C.composition_dir(cfg, comp)
    monkeypatch.setattr(R, "run_external", make_fake(d0, unconverged=(design.config_ids[1],)))
    out = R.run_queue(cfg, C.job_dirs(d0), timeout_s=10)
    assert out[design.config_ids[1]] == "not_converged"
    col = C.collect(cfg, d0)
    assert col["n_not_converged"] == 1
    fit = C.fit_j(cfg, d0)
    assert fit["excluded_configurations"] == {design.config_ids[1]: "not_converged"}
    assert design.config_ids[1] not in fit["fit_configurations"]


@pytest.mark.parametrize("what,status", [
    ("mtot", "valence_changed"),            # |M_tot - sum(sign 2S)| = 1 muB > 0.5
    ("ts", "smearing_large"),               # |TS| = 2 meV >= 1 meV
])
def test_valence_and_smearing_checks_exclude_a_configuration(tmp_cfg, hem_design, monkeypatch, what, status):
    cfg = tmp_cfg
    sysm, design = hem_design
    comp = f"hematite_{what}"
    C.generate(cfg, sysm, design, comp)
    d0 = C.composition_dir(cfg, comp)
    victim = design.config_ids[1]
    kw = {"total_mag_shift": {victim: 1.0}} if what == "mtot" else {"ts_Ry": {victim: -2.0 / H.RY_TO_MEV}}
    monkeypatch.setattr(R, "run_external", make_fake(d0, **kw))
    R.run_queue(cfg, C.job_dirs(d0), timeout_s=10)
    col = C.collect(cfg, d0)
    row = next(r for r in col["rows"] if r["config_id"] == victim)
    assert row["status"] == status and row["note"]
    assert col["n_ok"] == 4 and victim in C.fit_j(cfg, d0)["excluded_configurations"]


def test_fe2_location_is_checked_against_tr_ns(tmp_cfg, monkeypatch):
    """Toy cell with one Fe2+ (atom 2): the Fe with the largest Tr[ns] must be atom 2, else valence_changed."""
    from pymatgen.core import Species
    cfg = tmp_cfg
    sysm = hematite_system()
    st = sysm.structure.copy()
    fe = [i for i, x in enumerate(st) if x.specie.symbol == "Fe"]
    st.replace(fe[1], Species("Fe", 2))                                                       # one Fe2+
    sysm.structure = st
    sysm.spin[fe[1]] = 2.0
    design = C.make_design(cfg, sysm)
    comp = "hematite_fe2"
    C.generate(cfg, sysm, design, comp)
    d0 = C.composition_dir(cfg, comp)
    assert "Fe4" in (d0 / design.config_ids[0] / "pw.in").read_text() or "Fe3" in (d0 / design.config_ids[0] / "pw.in").read_text()
    nd_ok = {a: (4.9, 0.4) for a in range(1, 5)}
    nd_ok[fe[1] + 1] = (4.9, 1.3)                                                             # the extra electron sits on atom 2
    nd_bad = dict(nd_ok)
    nd_bad[fe[1] + 1], nd_bad[fe[0] + 1] = (4.9, 0.4), (4.9, 1.3)                             # it moved to atom 1
    bad = design.config_ids[1]
    monkeypatch.setattr(R, "run_external", make_fake(d0, ns_of={c: (nd_bad if c == bad else nd_ok) for c in design.config_ids},
                                                     total_mag_shift={c: 0.0 for c in design.config_ids}))
    R.run_queue(cfg, C.job_dirs(d0), timeout_s=10)
    col = C.collect(cfg, d0)
    rows = {r["config_id"]: r for r in col["rows"]}
    assert rows[bad]["status"] == "valence_changed" and "Fe2+ moved" in rows[bad]["note"]
    assert rows[design.config_ids[0]]["status"] == "ok"
    assert rows[design.config_ids[0]]["valence_check"]["fe2_found"] == [fe[1] + 1]


def test_generate_refuses_to_overwrite_finished_input(tmp_cfg, hem_design, monkeypatch):
    cfg = tmp_cfg
    sysm, design = hem_design
    comp = "hematite_over"
    C.generate(cfg, sysm, design, comp)
    d0 = C.composition_dir(cfg, comp)
    monkeypatch.setattr(R, "run_external", make_fake(d0))
    R.run_queue(cfg, C.job_dirs(d0), timeout_s=10)
    cfg2 = cfg.with_overrides(dft={**cfg.raw["dft"], "scf": {**cfg.raw["dft"]["scf"], "mixing_beta": 0.2}})
    with pytest.raises(C.CampaignError, match="force"):
        C.generate(cfg2, sysm, design, comp)


def test_runtime_estimate_from_feasibility(cfg):
    e = C.estimate_runtime(cfg, "W92", 17)
    assert e["cpu_hours_per_job"] == pytest.approx(290 * 40 / 3600)                 # 3.2 h per SCF
    assert e["cpu_hours_total"] == pytest.approx(17 * 3.2222, rel=1e-3) and e["waves"] == 6
    assert e["ram_GB_peak_estimate"] == pytest.approx(3 * 3.3) and e["memory_ok"]
    assert e["wall_hours_estimate"] == pytest.approx(6 * 3.2222 * 1.2, rel=1e-3)
    assert not C.estimate_runtime(cfg, "W92", 17, n_parallel=4)["memory_ok"]        # 4 x 3.3 GB > 10 GB budget
    assert C.composition_name({"Mn": 2.0}, None, 0) == "SrMn2Fe16O27_cfg0"
