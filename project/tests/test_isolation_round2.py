"""Isolation, round 2: copy2 / sqlite3 / child processes / multiprocessing / hook fast path.

Everything runs against a SYNTHETIC vault (pytest tmp_path via the `root` fixture). No test touches the real vault.
"""
from __future__ import annotations

import concurrent.futures as cf
import multiprocessing as mp
import os
import sqlite3
import subprocess
import sys
import textwrap
import time

import pytest

from hexw.io import paths
from test_isolation import NT, _mklink_junction, temp_module, vfile, vroot  # noqa: F401  (fixtures)

THEORY2 = textwrap.dedent('''
    import os, shutil, sqlite3, subprocess, sys

    def t_copy2(src, dst): return shutil.copy2(src, dst)
    def t_copy(src, dst): return shutil.copy(src, dst)
    def t_copyfile(src, dst): return shutil.copyfile(src, dst)
    def t_copytree(src, dst): return shutil.copytree(src, dst)
    def t_move(src, dst): return shutil.move(src, dst)
    def t_winapi_copyfile2(src, dst):
        import _winapi
        return _winapi.CopyFile2(src, dst, 0)
    def t_sqlite(db, uri=False):
        c = sqlite3.connect(db, uri=uri)
        try:
            return c.execute("select name from sqlite_master").fetchall()
        finally:
            c.close()
    def t_run(args, **kw): return subprocess.run(args, capture_output=True, text=True, **kw)
    def t_popen_env(env): return subprocess.run([sys.executable, "-c", "pass"], env=env, capture_output=True)
    def t_open(p, mode="rb"): return open(p, mode)
''')

PROBE_ENV = "import os; print(os.environ.get('HEXW_GUARD')); print(os.environ.get('HEXW_ROOT'))"


@pytest.fixture(params=["hexw.theory", "hexw.dft"])
def th(request, tmp_path):
    with temp_module(tmp_path, request.param, "probe_r2", THEORY2) as m:
        yield m


@pytest.fixture
def out(tmp_path):
    d = tmp_path / "out"
    d.mkdir()
    return d


# ---------------------------------------------------------------- 1. shutil.copy2 and friends
@pytest.mark.parametrize("fn", ["t_copy2", "t_copy", "t_copyfile", "t_move"])
def test_copy_family_from_theory_cannot_take_vault_file(vroot, th, out, fn):
    dst = out / "c.bin"
    with pytest.raises(paths.AccessError):
        getattr(th, fn)(vfile(vroot), str(dst))
    assert not dst.exists()
    assert os.path.exists(vfile(vroot))


@pytest.mark.skipif(not NT, reason="_winapi.CopyFile2 is Windows-only")
def test_native_copyfile2_event_is_guarded(vroot, th, out):
    dst = out / "c.bin"
    with pytest.raises(paths.AccessError):
        th.t_winapi_copyfile2(vfile(vroot), str(dst))
    assert not dst.exists()


def test_copytree_of_vault_dir_blocked(vroot, th, out):
    with pytest.raises(paths.AccessError):
        th.t_copytree(os.path.dirname(vfile(vroot)), str(out / "tree"))
    assert not (out / "tree").exists()


def test_copy_into_read_only_zones_blocked_for_everybody(vroot, th):
    src = vroot / "raw" / "theory_side" / "lit.txt"
    for dst in (vroot / "raw" / "theory_side" / "new.txt", vroot / "vault" / "new.txt"):
        with paths.test_access():
            for fn in (th.t_copy2, th.t_copyfile):
                with pytest.raises(paths.ReadOnlyZoneError):
                    fn(str(src), str(dst))
        assert not dst.exists()


def test_copy2_positive_controls(vroot, th, out):
    src = vroot / "raw" / "theory_side" / "lit.txt"
    th.t_copy2(str(src), str(out / "a.txt"))                   # theory-side file -> neutral place: fine
    assert (out / "a.txt").read_text() == "literature"
    with paths.test_access():                                   # test infrastructure may copy OUT of the vault
        import shutil
        shutil.copy2(vfile(vroot), out / "b.bin")
    assert (out / "b.bin").exists()


# ---------------------------------------------------------------- 2. sqlite3
def _db(path):
    c = sqlite3.connect(str(path))
    c.execute("create table t(x)")
    c.commit()
    c.close()


def _uri(p: str, quote: bool) -> str:
    from urllib.parse import quote as q
    p = p.replace("\\", "/")
    return "file:" + (q(p, safe="/:") if quote else p) + "?mode=ro"


def test_sqlite_plain_and_uri_forms_blocked_and_nothing_created(vroot, th):
    vdir = vroot / "vault" / "lab_validation" / "01_XRD"
    hashed = vdir / "S_x=1.1_6h#1##2026_100.txt"        # a '#' in the name, like the real files
    with paths.test_root(vroot / "_elsewhere", validate=False):       # build the synthetic file
        hashed.write_bytes(b"1 2")
    truncated = vdir / "S_x=1.1_6h"                     # what SQLite creates for an unescaped '#'
    def ls():
        with paths.test_access():
            return sorted(os.listdir(vdir))
    before = ls()
    p = str(hashed)
    cases = [
        (p, False),
        (_uri(p, quote=True), True),
        (_uri(p, quote=False), True),                      # unescaped '#': SQLite would open the truncated name
        ("file:///" + p.replace("\\", "/").replace("#", "%23"), True),
        ("file://localhost/" + p.replace("\\", "/").replace("#", "%23") + "?mode=rwc", True),
        ("file:" + p.replace("\\", "/").replace("#", "%23").replace("/", "%2F"), True),
        (str(truncated), False),
    ]
    for db, uri in cases:
        with pytest.raises(paths.AccessError):
            th.t_sqlite(db, uri)
    assert ls() == before and not truncated.exists()


def test_sqlite_original_zip_blocked(vroot, th):
    with pytest.raises(paths.AccessError):
        th.t_sqlite(str(vroot / "raw" / "original_zip" / "a.zip"))


def test_sqlite_write_into_theory_side_blocked_but_ro_uri_allowed(vroot, th, out):
    db = vroot / "raw" / "theory_side" / "ref.db"
    with paths.test_root(vroot / "_elsewhere", validate=False):
        _db(db)                                            # build a reference db while it is "elsewhere"
    with pytest.raises(paths.ReadOnlyZoneError):
        th.t_sqlite(str(db))                               # may create journal / modify: refused
    assert th.t_sqlite("file:" + str(db).replace("\\", "/") + "?mode=ro", True) == [("t",)]


def test_sqlite_positive_controls(vroot, th, out):
    _db(out / "n.db")
    assert th.t_sqlite(str(out / "n.db")) == [("t",)]
    assert th.t_sqlite(":memory:") == []
    assert th.t_sqlite("file:" + str(out / "n.db").replace("\\", "/") + "?mode=ro", True) == [("t",)]


def test_sqlite_whitelisted_caller_may_read_vault_db(vroot):
    db = vroot / "vault" / "lab_validation" / "x.db"
    with paths.test_root(vroot / "_elsewhere", validate=False):
        _db(db)
    with paths.test_access():
        c = sqlite3.connect("file:" + str(db).replace("\\", "/") + "?mode=ro", uri=True)
        assert c.execute("select name from sqlite_master").fetchall() == [("t",)]
        c.close()


def test_uri_decoding_helpers():
    t, ro = paths._sqlite_targets("file:///C:/a%20b/c%23d.db?mode=ro#frag")
    assert "C:/a b/c#d.db" in t and ro
    t, ro = paths._sqlite_targets("file:C:/v/x#1##y.txt")
    assert "C:/v/x" in t and "C:/v/x#1##y.txt" in t and not ro
    assert paths._sqlite_targets(":memory:") == ([], True)


# ---------------------------------------------------------------- 3. child processes
def test_theory_cannot_pass_other_hexw_root_in_env(vroot, th, tmp_path):
    decoy = tmp_path / "decoy"
    (decoy / "vault").mkdir(parents=True)
    for key in ("HEXW_ROOT", "hexw_root", "Hexw_Root"):
        env = dict(os.environ)
        env[key] = str(decoy)
        with pytest.raises(paths.TheoryAccessError):
            th.t_popen_env(env)


def test_theory_cannot_repoint_environ_root(vroot, th, tmp_path, monkeypatch):
    decoy = tmp_path / "decoy"
    (decoy / "vault").mkdir(parents=True)
    monkeypatch.setenv("HEXW_ROOT", str(decoy))
    with pytest.raises(paths.AccessError):
        th.t_run([sys.executable, "-c", "pass"])           # env=None: inherits the poisoned os.environ
    with pytest.raises(paths.AccessError):
        os.system("rem")                                    # external caller, no white-listed frame: refused too
    monkeypatch.undo()
    assert th.t_run([sys.executable, "-c", "pass"]).returncode == 0


def test_child_gets_guard_and_current_root(vroot, th):
    r = th.t_run([sys.executable, "-c", PROBE_ENV], env={"SYSTEMROOT": os.environ.get("SYSTEMROOT", "")})
    guard, root = r.stdout.split()
    assert guard == "1" and os.path.samefile(root, vroot)
    env = {"SYSTEMROOT": os.environ.get("SYSTEMROOT", ""), "hexw_guard": "0"}
    r = th.t_run([sys.executable, "-c", PROBE_ENV], env=env)
    assert r.stdout.split()[0] == "1" and env["HEXW_GUARD"] == "1"      # completed in place, forced to 1


def test_child_default_env_inherits_guard(vroot, th):
    r = th.t_run([sys.executable, "-c", PROBE_ENV])
    assert r.stdout.split()[0] == "1" and os.path.samefile(r.stdout.split()[1], os.environ["HEXW_ROOT"])


def test_child_importing_hexw_is_guarded_with_synthetic_root(vroot, th):
    code = ("import hexw, os, sys" + chr(10) + "try:" + chr(10) + "    open(sys.argv[1], 'rb')" + chr(10) +
            "except PermissionError as e:" + chr(10) + "    print('blocked', type(e).__name__)" + chr(10))
    r = th.t_run([sys.executable, "-c", code, str(vroot / "raw" / "theory_side" / "lit.txt")])
    assert r.stdout.strip() == "", r.stderr                  # theory_side stays readable in the child
    env = {k: v for k, v in os.environ.items() if k.upper() != "HEXW_ROOT"}      # completed with the CURRENT root
    r = th.t_run([sys.executable, "-c", code.replace("sys.argv[1]", "os.environ['T']")], env=dict(env, T=vfile(vroot)))
    assert r.stdout.strip().startswith("blocked"), (r.stdout, r.stderr)


def test_guarded_child_without_root_fails_loudly(monkeypatch):
    monkeypatch.delenv("HEXW_ROOT", raising=False)
    monkeypatch.setenv("HEXW_GUARD", "1")
    with pytest.raises(paths.RootNotFoundError, match="HEXW_GUARD"):
        paths._default_root()


# ---------------------------------------------------------------- 4. multiprocessing (spawn)
MP_WORKER = textwrap.dedent('''
    import sys
    def has_hexw(_=None): return int("hexw" in sys.modules)
    def read(p):
        with open(p, "rb") as f:
            return len(f.read())
    def read_name(p):
        try:
            return "LEAK %d" % read(p)
        except PermissionError as e:
            return "blocked " + type(e).__name__
''')


@pytest.fixture
def mp_worker(tmp_path, monkeypatch):
    d = tmp_path / "mpmods"
    d.mkdir()
    (d / "r2_mp_worker.py").write_text(MP_WORKER)
    monkeypatch.syspath_prepend(str(d))
    sys.modules.pop("r2_mp_worker", None)
    import r2_mp_worker
    yield r2_mp_worker
    sys.modules.pop("r2_mp_worker", None)


def test_multiprocessing_spawn_child_installs_hook_before_user_code(vroot, mp_worker):
    ctx = mp.get_context("spawn")
    with ctx.Pool(1) as pool:
        assert pool.apply(mp_worker.has_hexw) == 1                       # hexw imported in the child
        assert pool.apply(mp_worker.read_name, (vfile(vroot),)).startswith("blocked")
        assert pool.apply(mp_worker.read_name, (str(vroot / "raw" / "original_zip" / "a.zip"),)).startswith("blocked")
        assert pool.apply(mp_worker.read_name, (str(vroot / "raw" / "theory_side" / "lit.txt"),)).startswith("LEAK")


def test_process_pool_executor_spawn_child_is_guarded(vroot, mp_worker):
    with cf.ProcessPoolExecutor(1, mp_context=mp.get_context("spawn")) as ex:
        assert ex.submit(mp_worker.read_name, vfile(vroot)).result().startswith("blocked")


def test_mp_bootstrap_is_installed_once_and_documented():
    import multiprocessing.spawn as sp
    assert getattr(sp.get_command_line, "_hexw_wrapped", False)
    cmd = sp.get_command_line(parent_pid=1, pipe_handle=2)
    assert "import hexw" in cmd[cmd.index("-c") + 1] and "HEXW_GUARD" in cmd[cmd.index("-c") + 1]


# ---------------------------------------------------------------- 5. fast path of the hook
def test_fast_path_speed_outside_protected_zones(vroot):
    f = vroot / "raw" / "theory_side" / "lit.txt"
    p = str(f)
    for _ in range(3):
        open(p, "rb").close()
    best = 1e9
    n = 2000
    for _ in range(5):
        t = time.perf_counter()
        for _ in range(n):
            open(p, "rb").close()
        best = min(best, (time.perf_counter() - t) / n * 1e6)
    assert best < 400, f"{best:.0f} us/open (slow path is ~800 us, fast path ~50 us; loose bound for loaded CI)"


def test_fast_path_does_not_hide_write_into_read_only_zone(vroot):
    p = str(vroot / "raw" / "theory_side" / "lit.txt")
    for _ in range(3):
        open(p, "rb").close()
    assert paths._FAST.get(p) is not None
    for mode in ("wb", "ab", "r+b", "xb"):
        with pytest.raises(paths.ReadOnlyZoneError):
            open(p, mode)
    with pytest.raises(paths.ReadOnlyZoneError):
        os.open(p, os.O_WRONLY | getattr(os, "O_BINARY", 0))


@pytest.mark.skipif(not NT, reason="junctions are a Windows feature")
def test_fast_path_revalidates_when_a_directory_is_swapped_for_a_junction(vroot, th, tmp_path):
    d = tmp_path / "d"
    d.mkdir()
    (d / "s.txt").write_text("benign")
    p = str(d / "s.txt")
    for _ in range(3):
        th.t_open(p).close()                              # cached as OTHER
    assert paths._FAST.get(p) is not None
    d.rename(tmp_path / "d_old")
    _mklink_junction(d, vroot / "vault" / "lab_validation" / "01_XRD")
    assert os.path.exists(p)
    with pytest.raises(paths.AccessError):
        th.t_open(p)


def test_fast_path_revalidates_symlink_in_place_of_file(vroot, th, tmp_path):
    p = tmp_path / "f.txt"
    p.write_text("benign")
    for _ in range(3):
        th.t_open(str(p)).close()
    p.unlink()
    try:
        os.symlink(vfile(vroot), p)
    except OSError as e:
        pytest.skip(f"symlinks need privileges here: {e}")
    with pytest.raises(paths.AccessError):
        th.t_open(str(p))


def test_fast_path_cache_is_dropped_on_root_change(vroot, th):
    p = str(vroot / "raw" / "theory_side" / "lit.txt")
    for _ in range(3):
        open(p, "rb").close()
    assert p in paths._FAST
    with paths.test_root(vroot / "_elsewhere", validate=False):
        ent = paths._FAST.get(p)
        assert ent is not None and ent[1] is not paths._CFG        # not trusted under another config
