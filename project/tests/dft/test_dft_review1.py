"""Review round 1: keep-awake, orphan detection, Fe2+ species / starting_ns_eigenvalue, Tr[ns] parsing, valence check, relative
hold-out criterion, nuisance columns for far pairs."""
import os
import subprocess
import sys
import time
from datetime import datetime, timedelta, timezone

import numpy as np
import pytest

from hexw.dft import hamiltonian as H
from hexw.dft import parse as P
from hexw.dft import qe_input as Q
from hexw.dft import runner as R
from hexw.dft.hematite import hematite_system
from hexw.util import keepawake, proc
from _dftsynth import synthetic_pw_out


# ---- keep-awake -------------------------------------------------------------------------------------------------------
def test_keep_awake_sets_and_restores_the_flags(monkeypatch):
    calls = []
    monkeypatch.setattr(keepawake, "_set", lambda flags: calls.append(flags) or 1)
    monkeypatch.setattr(keepawake.sys, "platform", "win32")
    with keepawake.keep_awake() as st:
        assert st.active and calls == [keepawake.ES_CONTINUOUS | keepawake.ES_SYSTEM_REQUIRED]
        assert not (calls[0] & 0x2)                                   # ES_DISPLAY_REQUIRED is NOT requested
    assert calls[-1] == keepawake.ES_CONTINUOUS and not st.active     # restored: ES_CONTINUOUS alone clears the request


def test_keep_awake_failure_is_recorded_not_raised(monkeypatch):
    monkeypatch.setattr(keepawake, "_set", lambda flags: 0)
    monkeypatch.setattr(keepawake.sys, "platform", "win32")
    with keepawake.keep_awake() as st:
        assert not st.active and "returned 0" in st.error


def test_keep_awake_disabled_does_nothing(monkeypatch):
    monkeypatch.setattr(keepawake, "_set", lambda flags: pytest.fail("must not be called"))
    with keepawake.keep_awake(False) as st:
        assert not st.active


@pytest.mark.skipif(sys.platform != "win32", reason="Windows API")
def test_keep_awake_real_call_on_windows():
    with keepawake.keep_awake() as st:
        assert st.active and st.error is None


def test_queue_holds_keep_awake_for_its_duration(cfg, tmp_path, monkeypatch):
    events = []
    import contextlib

    @contextlib.contextmanager
    def fake_ka(enabled=True):
        events.append("on")
        yield keepawake.KeepAwakeState()
        events.append("off")
    monkeypatch.setattr(R, "keep_awake", fake_ka)
    d = tmp_path / "c00"
    d.mkdir()
    (d / "pw.in").write_text("&CONTROL\n/\n")

    def fake_run(exe, args=(), *, cwd=None, stdout_path=None, **kw):
        events.append("run")
        stdout_path.write_text(synthetic_pw_out(-5.0, [3.7, -3.7]))
        return subprocess.CompletedProcess([], 0, None, None)
    monkeypatch.setattr(R, "run_external", fake_run)
    j = R.Journal(tmp_path / "j.jsonl")
    R.run_queue(cfg, [d], timeout_s=5, journal=j)
    assert events == ["on", "run", "off"]
    assert '"event": "keep_awake"' in (tmp_path / "j.jsonl").read_text()


# ---- pid / orphan ---------------------------------------------------------------------------------------------------------
def test_pid_alive():
    assert proc.pid_alive(os.getpid())
    assert not proc.pid_alive(0) and not proc.pid_alive(-5) and not proc.pid_alive(None) and not proc.pid_alive("x")
    child = subprocess.Popen([sys.executable, "-c", "pass"])
    child.wait()
    deadline = time.time() + 5
    while proc.pid_alive(child.pid) and time.time() < deadline:
        time.sleep(0.05)
    assert not proc.pid_alive(child.pid)                              # finished and reaped


@pytest.fixture
def job(tmp_path):
    d = tmp_path / "c00"
    d.mkdir()
    (d / "pw.in").write_text("&CONTROL\n/ ! job\n")
    return d


def _stale(d, **kw):
    old = (datetime.now(timezone.utc) - timedelta(seconds=10_000)).isoformat()
    R.write_status(d, state="running", heartbeat_utc=old, **kw)


def test_stale_running_job_with_a_live_runner_is_orphan_and_not_restarted(cfg, job, monkeypatch):
    child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
    try:
        _stale(job, runner_pid=child.pid)
        calls = []
        monkeypatch.setattr(R, "run_external", lambda *a, **k: calls.append(1))
        assert R.classify(cfg, job) == "skip:orphan_running"
        st = R.read_status(job)
        assert st["state"] == "orphan_running" and f"runner pid {child.pid}" in st["message"]
        assert R.run_queue(cfg, [job], timeout_s=5)[job.name] == "skip:orphan_running" and not calls
        assert R.classify(cfg, job, retry_failed=True, force=True) == "skip:orphan_running"      # even --force
    finally:
        child.kill()
        child.wait()
    assert R.classify(cfg, job) == "run"                                # the process is gone: restart allowed


def test_live_pw_process_alone_also_blocks_a_restart(cfg, job):
    child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"])
    try:
        _stale(job, runner_pid=99999999, pw_pid=child.pid)              # runner dead, pw.x orphaned but alive
        assert R.classify(cfg, job) == "skip:orphan_running"
        assert f"pw.x pid {child.pid}" in R.read_status(job)["message"]
    finally:
        child.kill()
        child.wait()


def test_dead_runner_is_restarted_and_own_stale_record_does_not_block(cfg, job, monkeypatch):
    _stale(job, runner_pid=99999999)
    assert R.classify(cfg, job) == "run"
    _stale(job, runner_pid=os.getpid())                                 # this process, but it is not running that job
    assert R.classify(cfg, job) == "run"


def test_pid_is_recorded_while_running_and_cleared_at_the_end(cfg, job, monkeypatch):
    seen = {}

    def fake_run(exe, args=(), *, cwd=None, stdout_path=None, pid_callback=None, **kw):
        pid_callback(31337)
        seen.update(R.read_status(cwd))
        stdout_path.write_text(synthetic_pw_out(-5.0, [3.7, -3.7]))
        return subprocess.CompletedProcess([], 0, None, None)
    monkeypatch.setattr(R, "run_external", fake_run)
    R.run_queue(cfg, [job], timeout_s=5)
    assert seen["state"] == "running" and seen["pw_pid"] == 31337 and seen["runner_pid"] == os.getpid()
    assert R.read_status(job)["pw_pid"] is None


def test_heartbeat_cannot_overwrite_a_final_state(job):
    R.write_status(job, state="done")
    st = R.write_status(job, ("running", "orphan_running"), state="running", heartbeat_utc="x")
    assert st["state"] == "done" and R.read_status(job)["state"] == "done"


# ---- Fe2+ species and starting_ns_eigenvalue ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def hem_fe2():
    """Hematite cell with one spin-up Fe2+ (atom fe[0]) and one spin-down Fe2+ (atom fe[1]) in the ground state."""
    from pymatgen.core import Species
    sysm = hematite_system()
    st = sysm.structure.copy()
    fe = [i for i, x in enumerate(st) if x.specie.symbol == "Fe"]
    down = [i for i in fe if sysm.base_sign[i] < 0]
    up = [i for i in fe if sysm.base_sign[i] > 0]
    fe = [up[0], down[0], down[1], up[1]]
    st.replace(fe[0], Species("Fe", 2))
    st.replace(fe[1], Species("Fe", 2))
    sysm.structure = st
    sysm.spin[fe[0]] = sysm.spin[fe[1]] = 2.0
    return sysm, fe


def _nl(txt):
    return Q.parse_pw_input(txt)["namelists"]["system"]


def test_fe2_gets_its_own_species_magnetization_and_starting_ns(cfg, hem_fe2):
    sysm, fe = hem_fe2
    pat = np.ones(4)
    txt = Q.build_input(sysm, pat, cfg, prefix="t", config_id="x")
    assert Q.validate_input(txt, cfg, sysm, pat) == []
    P_ = Q.parse_pw_input(txt)
    labels = [s[0] for s in P_["species"]]
    assert labels == ["Fe1", "Fe2", "Fe3", "Fe4", "O"] and all(len(l) <= 3 for l in labels)   # Fe3 / Fe4: Fe2+ spin up / down
    pos = {i: lab for i, (lab, _) in enumerate(P_["positions"])}
    assert [pos[k] for k in fe] == ["Fe3", "Fe4", "Fe2", "Fe1"]
    nl = _nl(txt)
    assert float(nl["starting_magnetization(3)"]) == pytest.approx(0.32)                       # own value, sign of the target
    assert float(nl["starting_magnetization(4)"]) == pytest.approx(-0.32)
    assert float(nl["starting_magnetization(1)"]) == pytest.approx(0.4) and float(nl["starting_magnetization(2)"]) == pytest.approx(-0.4)
    # spin-up Fe2+ (species 3): the minority channel is spin 2; spin-down Fe2+ (species 4): spin 1; the 5th eigenvalue gets 1 electron
    assert float(nl["starting_ns_eigenvalue(5,2,3)"]) == 1.0 and float(nl["starting_ns_eigenvalue(5,1,4)"]) == 1.0
    assert sum(k.startswith("starting_ns_eigenvalue") for k in nl) == 2
    assert ("U", "Fe3-3d", 3.0) in P_["U"] and ("U", "Fe4-3d", 3.0) in P_["U"]
    # the same ions with the sublattices turned over: the species swap, the seeding follows the target spin
    pat2 = np.array([-1.0, -1, -1, -1])
    txt2 = Q.build_input(sysm, pat2, cfg, prefix="t")
    assert Q.validate_input(txt2, cfg, sysm, pat2) == []
    pos2 = {i: lab for i, (lab, _) in enumerate(Q.parse_pw_input(txt2)["positions"])}
    assert [pos2[k] for k in fe] == ["Fe4", "Fe3", "Fe1", "Fe2"]


def test_without_localization_fe2_is_an_ordinary_fe_species(cfg, hem_fe2):
    sysm, _ = hem_fe2
    cfg0 = cfg.with_overrides(dft={**cfg.raw["dft"], "fe2_localization": {"mode": "none"}})
    txt = Q.build_input(sysm, np.ones(4), cfg0, prefix="t")
    assert {s[0] for s in Q.parse_pw_input(txt)["species"]} == {"Fe1", "Fe2", "O"}
    assert "starting_ns_eigenvalue" not in txt and Q.validate_input(txt, cfg0, sysm, np.ones(4)) == []


def test_starting_ns_can_be_switched_off_and_is_validated(cfg, hem_fe2):
    sysm, _ = hem_fe2
    loc = {**cfg.raw["dft"]["fe2_localization"], "starting_ns_eigenvalue": {"use": False}}
    cfg1 = cfg.with_overrides(dft={**cfg.raw["dft"], "fe2_localization": loc})
    txt = Q.build_input(sysm, np.ones(4), cfg1, prefix="t")
    assert "starting_ns_eigenvalue" not in txt and "Fe4" in txt
    assert Q.validate_input(txt, cfg1, sysm, np.ones(4)) == []
    txt = Q.build_input(sysm, np.ones(4), cfg, prefix="t")
    old = "starting_ns_eigenvalue(5,1,4)=1"
    assert old in txt
    for new, frag in [("starting_ns_eigenvalue(6,1,4)=1", "bad indices"), ("starting_ns_eigenvalue(5,3,4)=1", "bad indices"),
                      ("starting_ns_eigenvalue(5,1,9)=1", "bad indices"), ("starting_ns_eigenvalue(5,1,4)=2", "bad indices"),
                      ("starting_ns_eigenvalue(5,1,5)=1", "no Hubbard U")]:                     # species 5 is O
        probs = Q.validate_input(txt.replace(old, new), cfg, sysm, np.ones(4))
        assert any(frag in p for p in probs), (new, probs)


def test_gamma_tricks_need_the_atomic_projector(cfg):
    sysm = hematite_system()
    k = {**cfg.raw["systems"]["hematite"], "kpoints": [1, 1, 1], "gamma_tricks": True}
    c_atomic = cfg.with_overrides(systems={**cfg.raw["systems"], "hematite": k})
    from hexw.dft import settings as ST
    txt = Q.build_input(sysm, np.ones(4), c_atomic, prefix="t")
    assert "K_POINTS gamma" in txt and "{atomic}" in txt and Q.validate_input(txt, c_atomic, sysm, np.ones(4)) == []
    hub = {**cfg.raw["dft"]["hubbard"], "projector": "ortho-atomic"}
    with pytest.raises(ST.DftConfigError, match="gamma_tricks"):
        c_atomic.with_overrides(dft={**cfg.raw["dft"], "hubbard": hub})                      # refused at config level
    raw = {**c_atomic.raw, "dft": {**cfg.raw["dft"], "hubbard": hub}, "systems": {**c_atomic.raw["systems"]}}
    raw["systems"]["W92"] = {**raw["systems"]["W92"], "gamma_tricks": False}
    c_bad = ST.DftConfig(raw, "x", None)                                                      # bypasses the config validation
    with pytest.raises(Q.InputError, match="Gamma-only"):
        Q.build_input(sysm, np.ones(4), c_bad, prefix="t")


# ---- Tr[ns] and the valence check ----------------------------------------------------------------------------------------
def test_full_tr_ns_is_parsed():
    txt = ("     Tr[ns(  1)] (up, down, total) =   5.00000  0.00000  5.00000\n"
           "     Tr[ns(  2)] (up, down, total) =   0.00000  5.00000  5.00000\n"
           "     Tr[ns(  1)] (up, down, total) =   4.96480  1.22242  6.18722\n"
           "     Tr[ns(  2)] (up, down, total) =   1.39240  4.91648  6.30888\n"
           "     total magnetization       =     1.00 Bohr mag/cell\n")
    r = P.parse_pw_output(txt)
    assert r.hubbard_ns[1] == {"up": 4.9648, "down": 1.22242, "total": 6.18722}            # the LAST block
    assert r.hubbard_ns_start[1] == {"up": 5.0, "down": 0.0, "total": 5.0} and r.total_magnetization == 1.0
    assert r.hubbard_d_moments[2] == pytest.approx(1.39240 - 4.91648)


def _atoms(ions, spins=None):
    out = []
    for i, ion in enumerate(ions):
        el = ion.rstrip("+0123456789")
        S = {"Fe2+": 2.0, "Fe3+": 2.5, "Mn2+": 2.5, "Mn3+": 2.0}.get(ion, 0.0) if spins is None else spins[i]
        out.append({"i": i + 1, "species": el, "ion": ion, "spin_S": S, "group": 0 if S else -1, "base_sign": 1.0})
    return out


def _res(moments, ns=None, M=None, ts=0.0):
    nat = len(moments)
    return P.parse_pw_output(synthetic_pw_out(-1.0, moments, total_mag=0.0 if M is None else M, ts_Ry=ts, ns=ns))


def test_valence_check_accepts_and_rejects():
    atoms = _atoms(["Fe3+", "Fe2+", "Fe3+", "O"])
    sg = np.array([1.0, 1.0, 1.0, 0.0])
    M = 5.0 + 4.0 + 5.0
    ok = P.check_valence(_res([3.8, 3.1, 3.8, 0.0], {1: (4.9, 0.4), 2: (4.9, 1.3), 3: (4.9, 0.4)}, M), atoms, sg)
    assert ok["ok"] and ok["fe2_found"] == [2] and ok["nd_gap"] == pytest.approx(0.9) and ok["M_expected_muB"] == 14.0
    moved = P.check_valence(_res([3.1, 3.8, 3.8, 0.0], {1: (4.9, 1.3), 2: (4.9, 0.4), 3: (4.9, 0.4)}, M), atoms, sg)
    assert not moved["valence_ok"] and "Fe2+ moved" in moved["reasons"][0] and moved["fe2_found"] == [1]
    m_off = P.check_valence(_res([3.8, 3.1, 3.8, 0.0], {1: (4.9, 0.4), 2: (4.9, 1.3), 3: (4.9, 0.4)}, M + 0.6), atoms, sg)
    assert not m_off["valence_ok"] and "total magnetization" in m_off["reasons"][0]
    edge = P.check_valence(_res([3.8, 3.1, 3.8, 0.0], {1: (4.9, 0.4), 2: (4.9, 1.3), 3: (4.9, 0.4)}, M + 0.5), atoms, sg)
    assert edge["valence_ok"]                                                              # exactly 0.5 is allowed (> 0.5 fails)
    ts = P.check_valence(_res([3.8, 3.1, 3.8, 0.0], {1: (4.9, 0.4), 2: (4.9, 1.3), 3: (4.9, 0.4)}, M, ts=-1.5 / H.RY_TO_MEV), atoms, sg)
    assert ts["valence_ok"] and not ts["smearing_ok"] and not ts["ok"] and ts["TS_meV"] == pytest.approx(1.5, rel=1e-4)


def test_valence_check_falls_back_to_the_smallest_moment_without_tr_ns():
    atoms = _atoms(["Fe3+", "Fe2+", "Fe3+", "O"])
    sg = np.array([1.0, 1.0, 1.0, 0.0])
    good = P.check_valence(_res([3.8, 3.1, 3.8, 0.0], None, 14.0), atoms, sg)
    assert good["ok"] and good["fe2_basis"] == "smallest |m|"
    bad = P.check_valence(_res([3.1, 3.8, 3.8, 0.0], None, 14.0), atoms, sg)
    assert not bad["valence_ok"]
    none = P.check_valence(P.parse_pw_output("!    total energy = -1.0 Ry\n total magnetization = 14 Bohr mag/cell\n"), atoms, sg)
    assert not none["valence_ok"] and any("cannot be checked" in r for r in none["reasons"])


def test_mn_must_stay_2plus():
    atoms = _atoms(["Fe3+", "Mn2+", "Fe3+", "O"])
    sg = np.array([1.0, -1.0, 1.0, 0.0])
    M = 5.0 - 5.0 + 5.0
    ok = P.check_valence(_res([3.8, -3.8, 3.8, 0.0], None, M), atoms, sg)
    assert ok["valence_ok"]
    mn3 = P.check_valence(_res([3.8, -3.0, 3.8, 0.0], None, M), atoms, sg)            # d4-like moment: 0.79 of the Fe3+ moment
    assert not mn3["valence_ok"] and "Mn2+ not preserved" in mn3["reasons"][0]
    m_off = P.check_valence(_res([3.8, -3.8, 3.8, 0.0], None, M + 1.0), atoms, sg)    # one electron moved: M_tot differs by 1
    assert not m_off["valence_ok"]


# ---- relative hold-out criterion --------------------------------------------------------------------------------------------
def _toy(n_hold_err):
    rng = np.random.default_rng(3)
    X = rng.normal(size=(8, 3)) * 10
    J = np.array([-10.0, 4.0, -2.0])
    E = 5.0 + X @ J
    Xh = rng.normal(size=(2, 3)) * 10
    Eh = 5.0 + Xh @ J + n_hold_err
    return X, E, Xh, Eh


def test_holdout_criterion_is_relative_to_the_energy_range():
    X, E, Xh, Eh = _toy(0.0)
    span = np.ptp(E)
    good = H.fit_exchange(X, E, list("abc"), noise_floor_meV=10.0, X_holdout=Xh, E_holdout_meV=Eh + 0.015 * span)
    assert good.holdout_consistent is True and good.holdout["rel_rms"] == pytest.approx(0.015, rel=1e-6)
    assert np.allclose(good.sigma_J_total, good.sigma_J)
    bad = H.fit_exchange(X, E, list("abc"), noise_floor_meV=10.0, X_holdout=Xh, E_holdout_meV=Eh + 0.03 * span)
    assert bad.holdout_consistent is False and bad.sigma_model_meV == pytest.approx(0.03 * span, rel=1e-6)
    assert np.all(bad.sigma_J_total > bad.sigma_J) and any("hold-out disagrees" in n for n in bad.notes)
    # a z-score against a small floor would have rejected the 1.5 % case; the relative criterion does not depend on the floor
    small = H.fit_exchange(X, E, list("abc"), noise_floor_meV=0.5, X_holdout=Xh, E_holdout_meV=Eh + 0.015 * span)
    assert small.holdout_consistent is True and max(abs(z) for z in small.holdout["z"]) > 3


def test_hematite_numbers_pass_the_relative_criterion():
    """7.3 meV on a range of 1430 meV (docs section 9) is 0.5 %: consistent; it was z = -13 against the old 0.5 meV floor."""
    assert 7.3 / 1430.0 < 0.02


# ---- nuisance columns -------------------------------------------------------------------------------------------------------------
def _alias_problem():
    rng = np.random.default_rng(11)
    n = 9
    X = rng.choice([-1.0, 1.0], size=(n, 3)) * 6
    N = 0.8 * X[:, :2] + rng.normal(size=(n, 2)) * 2.0                      # far-pair columns partly collinear with the classes
    return X, N


def test_nuisance_columns_inflate_sigma_and_remove_nothing_when_the_prior_is_tight():
    X, N = _alias_problem()
    E = 3.0 + X @ np.array([-5.0, 2.0, -1.0])
    base = H.fit_exchange(X, E, list("abc"), noise_floor_meV=10.0)
    wide = H.fit_exchange(X, E, list("abc"), noise_floor_meV=10.0, X_nuis=N, nuis_prior_sigma_meV=0.1, nuis_names=["f1", "f2"])
    assert np.all(wide.sigma_J > base.sigma_J)                              # aliasing is part of the uncertainty
    assert wide.nuisance["names"] == ["f1", "f2"] and np.all(np.array(wide.nuisance["sigma_inflation"]) > 1)
    assert np.allclose(wide.J, base.J, atol=1e-6)                           # data without far-pair energy: unbiased
    tight = H.fit_exchange(X, E, list("abc"), noise_floor_meV=10.0, X_nuis=N, nuis_prior_sigma_meV=1e-6)
    assert np.allclose(tight.sigma_J, base.sigma_J, rtol=1e-6)
    loose = H.fit_exchange(X, E, list("abc"), noise_floor_meV=10.0, X_nuis=N, nuis_prior_sigma_meV=1.0)
    assert np.all(loose.sigma_J >= wide.sigma_J)                            # monotone in the prior sigma
    assert wide.dof == base.dof                                             # nuisance parameters are not counted as fitted


def test_nuisance_matches_the_gaussian_marginal_formula():
    """sigma(J)^2 = sigma^2 [(A^T A)^-1 + (A^T A)^-1 A^T N (tau^-2 sigma^2 I + N^T P N)^-1 N^T A (A^T A)^-1]  (marginal of the joint
    posterior), checked against the direct inverse of the joint matrix."""
    X, N = _alias_problem()
    E = X @ np.array([-5.0, 2.0, -1.0])
    tau, floor = 0.1, 10.0
    fr = H.fit_exchange(X, E, list("abc"), noise_floor_meV=floor, X_nuis=N, nuis_prior_sigma_meV=tau)
    A = H.augmented(X)
    Bm = np.hstack([A, N])
    Mm = Bm.T @ Bm + np.diag([0] * 4 + [(floor / tau) ** 2] * 2)
    cov = floor ** 2 * np.linalg.inv(Mm)
    assert np.allclose(fr.sigma_J, np.sqrt(np.diag(cov))[1:4], rtol=1e-8)
    with pytest.raises(ValueError, match="prior sigma"):
        H.fit_exchange(X, E, list("abc"), noise_floor_meV=floor, X_nuis=N)


def test_nuisance_columns_enter_the_holdout_prediction():
    X, N = _alias_problem()
    g = np.array([0.0, 0.0])
    E = 3.0 + X @ np.array([-5.0, 2.0, -1.0]) + N @ g
    fr = H.fit_exchange(X[:7], E[:7], list("abc"), noise_floor_meV=10.0, X_nuis=N[:7], nuis_prior_sigma_meV=0.1,
                        X_holdout=X[7:], E_holdout_meV=E[7:], X_nuis_holdout=N[7:])
    assert fr.holdout_consistent and len(fr.holdout["error_meV"]) == 2
