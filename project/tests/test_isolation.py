"""Data isolation against REAL bypasses (reviewer probes, stage 1 round 1).

Nothing here uses laboratory data: the "vault" is a synthetic tree in tmp_path. Modules of the theory side are
real modules imported from a temporary directory added to the `__path__` of the hexw subpackage, so the frames
are identified exactly as in production (name / __spec__ / __package__ / file)."""
from __future__ import annotations

import concurrent.futures as cf
import contextlib
import ctypes
import glob
import importlib
import os
import runpy
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from hexw.io import paths

NL = chr(10)
NT = os.name == "nt"
PKG_DIR = Path(paths.__file__).resolve().parents[1]

THEORY_SRC = textwrap.dedent('''
    import concurrent.futures as cf
    import glob, os, pathlib, subprocess, sys, threading
    import numpy as np
    import pandas as pd
    from hexw.io import paths

    def t_open(p): return open(p, "rb").read(3)
    def t_open_data(p): return paths.open_data(p, "rb", purpose="metadata").read(3)
    def t_sha(p): return paths.sha256_file(p)
    def t_os_open(p):
        fd = os.open(p, os.O_RDONLY | getattr(os, "O_BINARY", 0))
        try: return os.read(fd, 3)
        finally: os.close(fd)
    def t_os_open_write(p):
        return os.open(p, os.O_WRONLY | os.O_CREAT | getattr(os, "O_BINARY", 0))
    def t_pathlib(p): return pathlib.Path(p).read_bytes()
    def t_loadtxt(p): return np.loadtxt(p)
    def t_fromfile(p): return np.fromfile(p, dtype=np.uint8, count=3)
    def t_read_csv(p): return pd.read_csv(p, sep=" ", header=None)
    def t_listdir(p): return os.listdir(os.path.dirname(p))
    def t_scandir(p): return list(os.scandir(os.path.dirname(p)))
    def t_iterdir(p): return list(pathlib.Path(p).parent.iterdir())
    def t_glob(p): return glob.glob(os.path.dirname(p) + "/*")
    def t_popen(p): return subprocess.run([sys.executable, "-c", "pass", p], capture_output=True)
    def t_shell(p): return subprocess.run('type "' + p + '"', shell=True, capture_output=True)
    def t_pool_open(p):
        with cf.ThreadPoolExecutor(1) as ex:
            return ex.submit(open, p, "rb").result().read(3)
    def t_pool_fn(p):
        with cf.ThreadPoolExecutor(1) as ex:
            return ex.submit(t_open, p).result()
    def t_pool_open_data(p):
        with cf.ThreadPoolExecutor(1) as ex:
            return ex.submit(paths.open_data, p, "rb").result().read(3)
    def t_thread(p):
        box = []
        def w():
            try: box.append(open(p, "rb").read(3))
            except BaseException as e: box.append(e)
        t = threading.Thread(target=w); t.start(); t.join()
        if isinstance(box[0], BaseException): raise box[0]
        return box[0]
    def t_configure(p): return paths.configure(os.path.dirname(p))
''')

CHANNELS = ["t_open", "t_open_data", "t_sha", "t_os_open", "t_os_open_write", "t_pathlib", "t_loadtxt",
            "t_fromfile", "t_read_csv", "t_listdir", "t_scandir", "t_iterdir", "t_glob", "t_popen", "t_shell",
            "t_pool_open", "t_pool_fn", "t_pool_open_data", "t_thread"]
THEORY_PKGS = ["hexw.theory", "hexw.dft", "hexw.structures", "hexw.spin", "hexw.params"]


@pytest.fixture
def vroot(root):
    d = root / "vault" / "lab_validation" / "01_XRD"
    with paths.test_root(root / "_elsewhere", validate=False):      # fixtures may build the fake vault
        d.mkdir(parents=True)
        (d / "s.txt").write_text("1 2\n3 4\n")
    return root


def vfile(r) -> str:
    return str(r / "vault" / "lab_validation" / "01_XRD" / "s.txt")


@contextlib.contextmanager
def temp_module(tmp_path, pkg, name, source, *, do_import=True):
    parent = importlib.import_module(pkg)
    d = tmp_path / ("mods_" + pkg.replace(".", "_"))
    d.mkdir(exist_ok=True)
    (d / f"{name}.py").write_text(source, encoding="utf-8")
    parent.__path__.append(str(d))
    importlib.invalidate_caches()
    try:
        yield importlib.import_module(f"{pkg}.{name}") if do_import else None
    finally:
        parent.__path__.remove(str(d))
        sys.modules.pop(f"{pkg}.{name}", None)


@pytest.fixture(params=THEORY_PKGS)
def theory_mod(request, tmp_path):
    with temp_module(tmp_path, request.param, "probe_x", THEORY_SRC) as m:
        yield m


# ------------------------------------------------------------------ real modules, every channel
@pytest.mark.parametrize("ch", CHANNELS)
def test_real_theory_like_module_is_blocked_on_every_channel(vroot, theory_mod, ch):
    with pytest.raises(paths.AccessError):
        getattr(theory_mod, ch)(vfile(vroot))


def test_real_theory_module_blocked_even_inside_test_access(vroot, theory_mod):
    """The explicit test context does not whitewash a theory frame in the stack."""
    with paths.test_access():
        for ch in ("t_open", "t_open_data", "t_loadtxt", "t_read_csv", "t_pathlib"):
            with pytest.raises(paths.TheoryAccessError):
                getattr(theory_mod, ch)(vfile(vroot))


def test_theory_can_still_read_theory_side_and_registry(vroot, theory_mod):
    lit = vroot / "raw" / "theory_side" / "lit.txt"
    assert theory_mod.t_open(str(lit)) == b"lit"
    assert theory_mod.t_pathlib(str(lit)) == b"literature"
    assert "lit.txt" in theory_mod.t_listdir(str(lit))
    assert theory_mod.t_pool_open(str(lit)) == b"lit"


def test_theory_cannot_write_into_raw_even_theory_side(vroot, theory_mod):
    with pytest.raises(paths.ReadOnlyZoneError):
        theory_mod.t_os_open_write(str(vroot / "raw" / "theory_side" / "new.bin"))


def test_original_zip_is_restricted_for_theory(vroot, theory_mod):
    z = str(vroot / "raw" / "original_zip" / "a.zip")
    for ch in ("t_open", "t_pathlib", "t_listdir", "t_pool_open"):
        with pytest.raises(paths.AccessError):
            getattr(theory_mod, ch)(z)


# ------------------------------------------------------------------ white list: the test thread itself
def test_non_whitelisted_caller_is_denied_by_default(vroot):
    """Absence of a prohibition is not a permission: plain test code has no allowed frame."""
    with pytest.raises(paths.AccessError):
        open(vfile(vroot), "rb")
    with pytest.raises(paths.AccessError):
        os.listdir(os.path.dirname(vfile(vroot)))
    with pytest.raises(paths.AccessError):
        paths.open_data(vfile(vroot), "rb", purpose="metadata")


def test_test_access_context_is_explicit_and_per_thread(vroot):
    with paths.test_access():
        assert open(vfile(vroot), "rb").read(3) == b"1 2"
        assert os.listdir(os.path.dirname(vfile(vroot))) == ["s.txt"]
        with cf.ThreadPoolExecutor(1) as ex:                    # pool thread does NOT inherit it
            with pytest.raises(paths.AccessError):
                ex.submit(open, vfile(vroot), "rb").result()
    with pytest.raises(paths.AccessError):
        open(vfile(vroot), "rb")


def test_numpy_pandas_from_test_code_without_context_denied(vroot):
    import numpy as np
    import pandas as pd
    with pytest.raises(paths.AccessError):
        np.loadtxt(vfile(vroot))
    with pytest.raises(paths.AccessError):
        pd.read_csv(vfile(vroot), sep=" ", header=None)


def test_subprocess_with_vault_path_denied_for_non_whitelisted(vroot):
    with pytest.raises(paths.AccessError):
        subprocess.run([sys.executable, "-c", "pass", vfile(vroot)])
    with pytest.raises(paths.AccessError):
        subprocess.run("type " + vfile(vroot), shell=True)
    with pytest.raises(paths.AccessError):
        subprocess.run([sys.executable, "-c", "pass"], cwd=os.path.dirname(vfile(vroot)))
    with paths.test_access():
        subprocess.run([sys.executable, "-c", "pass", vfile(vroot)], check=True)


def test_listing_glob_denied_for_non_whitelisted(vroot):
    for fn in (lambda: glob.glob(os.path.dirname(vfile(vroot)) + "/*"),
               lambda: list(os.scandir(vroot / "vault")),
               lambda: glob.glob("*", root_dir=os.path.dirname(vfile(vroot)))):
        with pytest.raises(paths.AccessError):
            fn()
    assert list((vroot / "vault").rglob("*")) == []                  # pathlib swallows OSError: nothing listed
    assert glob.glob(str(vroot / "raw" / "theory_side") + "/*")        # theory_side stays open


# ------------------------------------------------------------------ frame identity: runpy / -m / file
def test_runpy_main_of_theory_module_is_denied(vroot, tmp_path):
    src = textwrap.dedent('''
        import os
        try:
            open(os.environ["PROBE_TARGET"], "rb").read(1)
            RESULT = "LEAK"
        except PermissionError as e:
            RESULT = "blocked:" + type(e).__name__
    ''')
    with temp_module(tmp_path, "hexw.theory", "probe_main", src, do_import=False):
        os.environ["PROBE_TARGET"] = vfile(vroot)
        try:
            g = runpy.run_module("hexw.theory.probe_main", run_name="__main__")
        finally:
            del os.environ["PROBE_TARGET"]
    assert g["RESULT"] == "blocked:TheoryAccessError"


def test_runpy_run_path_script_has_no_role(vroot, tmp_path):
    script = tmp_path / "loose_script.py"
    script.write_text(f"open({vfile(vroot)!r}, 'rb').read(1)\n")
    with pytest.raises(paths.AccessError):
        runpy.run_path(str(script), run_name="__main__")


def test_python_dash_m_subprocess(vroot, tmp_path):
    site = tmp_path / "site"
    site.mkdir()
    (site / "sitecustomize.py").write_text(textwrap.dedent('''
        import os
        d = os.environ.get("HEXW_TEST_THEORY_DIR")
        if d:
            import hexw.theory, hexw.dft
            hexw.theory.__path__.append(d); hexw.dft.__path__.append(d)
    '''))
    mods = tmp_path / "mods_cli"
    mods.mkdir()
    (mods / "probe_cli.py").write_text(textwrap.dedent('''
        import os, sys
        print("name=" + __name__)
        try:
            open(os.environ["PROBE_TARGET"], "rb").read(1)
            print("RESULT=LEAK")
        except PermissionError as e:
            print("RESULT=blocked:" + type(e).__name__)
    '''))
    env = dict(os.environ, PYTHONPATH=str(site), HEXW_ROOT=str(vroot), HEXW_TEST_THEORY_DIR=str(mods),
               PROBE_TARGET=vfile(vroot))
    for pkg in ("hexw.theory", "hexw.dft"):
        r = subprocess.run([sys.executable, "-m", f"{pkg}.probe_cli"], env=env, cwd=str(tmp_path),
                           capture_output=True, text=True)
        assert "name=__main__" in r.stdout, r.stderr
        assert "RESULT=blocked:TheoryAccessError" in r.stdout, (r.stdout, r.stderr)


def test_import_hexw_fails_loudly_without_root(tmp_path):
    empty = tmp_path / "emptyroot"
    empty.mkdir()
    with paths.test_access():       # passing a different HEXW_ROOT to a child is allowed only for test infrastructure
        r = subprocess.run([sys.executable, "-c", "import hexw"], env=dict(os.environ, HEXW_ROOT=str(empty)),
                           capture_output=True, text=True, cwd=str(tmp_path))
    assert r.returncode != 0 and "HEXW_ROOT" in r.stderr


def test_hook_is_active_after_plain_import_hexw(vroot, tmp_path):
    """`import hexw` alone (no explicit paths import) installs the hook in a fresh process."""
    code = ("import hexw, os"+NL+"try:"+NL+"    open(os.environ['PROBE_TARGET'], 'rb')"+NL+"except PermissionError:"+NL +
            "    print('blocked')"+NL)
    r = subprocess.run([sys.executable, "-c", code],
                       env=dict(os.environ, HEXW_ROOT=str(vroot), PROBE_TARGET=vfile(vroot)),
                       capture_output=True, text=True, cwd=str(tmp_path))
    assert r.stdout.strip() == "blocked", r.stderr


@pytest.mark.parametrize("sub", ["theory", "dft", "structures", "spin", "params"])
def test_frame_identified_by_file_name_inside_package(vroot, sub):
    """__name__ == '__main__' but the source file lives in hexw/<sub>/ -> identified by the file."""
    fn = str(PKG_DIR / sub / "fake_main.py")
    g = {"__name__": "__main__"}
    exec(compile("def go(p):\n    return open(p, 'rb').read(1)\n", fn, "exec"), g)
    with pytest.raises(paths.TheoryAccessError):
        g["go"](vfile(vroot))


def test_conflicting_name_and_file_is_denied(vroot):
    (vroot / "registry" / "FROZEN_MODEL.json").write_text("{}")
    fn = str(PKG_DIR / "theory" / "disguised.py")
    g = {"__name__": "hexw.experiment.disguised"}
    exec(compile("def go(p):\n    return open(p, 'rb').read(1)\n", fn, "exec"), g)
    with pytest.raises(paths.AccessError):
        g["go"](vfile(vroot))


def test_file_inside_whitelisted_package_counts_even_if_main(vroot):
    (vroot / "registry" / "FROZEN_MODEL.json").write_text("{}")
    fn = str(PKG_DIR / "experiment" / "some_tool.py")
    g = {"__name__": "__main__"}
    exec(compile("def go(p):\n    return open(p, 'rb').read(3)\n", fn, "exec"), g)
    assert g["go"](vfile(vroot)) == b"1 2"


def test_spec_name_alone_is_enough_to_deny(vroot):
    import types
    g = {"__name__": "__main__", "__spec__": types.SimpleNamespace(name="hexw.dft.runner")}
    exec(compile("def go(p):\n    return open(p, 'rb').read(1)\n", "<spec>", "exec"), g)
    with pytest.raises(paths.TheoryAccessError):
        g["go"](vfile(vroot))
    g = {"__name__": "__main__", "__package__": "hexw.spin"}
    exec(compile("def go(p):\n    return open(p, 'rb').read(1)\n", "<pkg>", "exec"), g)
    with pytest.raises(paths.TheoryAccessError):
        g["go"](vfile(vroot))


# ------------------------------------------------------------------ the allowed path
EXPERIMENT_SRC = textwrap.dedent('''
    import concurrent.futures as cf
    from hexw.io import paths
    def read_open_data(p): return paths.open_data(p, "rb").read(3)
    def read_builtin(p): return open(p, "rb").read(3)
    def read_in_pool(p):
        with cf.ThreadPoolExecutor(1) as ex:
            return ex.submit(read_builtin, p).result()
    def read_in_pool_bare(p):
        with cf.ThreadPoolExecutor(1) as ex:
            return ex.submit(open, p, "rb").result().read(3)
''')


@pytest.mark.parametrize("pkg", ["hexw.experiment", "hexw.validate"])
def test_experiment_access_requires_freeze_then_works(vroot, tmp_path, pkg):
    with temp_module(tmp_path, pkg, "reader", EXPERIMENT_SRC) as m:
        for fn in (m.read_open_data, m.read_builtin, m.read_in_pool):
            with pytest.raises(paths.FreezeRequiredError):
                fn(vfile(vroot))
        from test_freeze import completed
        from hexw.registry.freeze import freeze_run
        freeze_run(completed(vroot))  # production-valid proof also works in worker threads
        assert m.read_open_data(vfile(vroot)) == b"1 2"
        assert m.read_builtin(vfile(vroot)) == b"1 2"
        assert m.read_in_pool(vfile(vroot)) == b"1 2"
        with pytest.raises(paths.AccessError):                  # bare open() in a pool: no allowed frame
            m.read_in_pool_bare(vfile(vroot))


def test_experiment_helper_called_from_theory_is_denied(vroot, tmp_path):
    (vroot / "registry" / "FROZEN_MODEL.json").write_text("{}")
    with temp_module(tmp_path, "hexw.experiment", "reader", EXPERIMENT_SRC) as m, \
            temp_module(tmp_path, "hexw.theory", "caller", "def go(f, p):\n    return f(p)\n") as t:
        with pytest.raises(paths.TheoryAccessError):
            t.go(m.read_open_data, vfile(vroot))


def test_audit_and_ingest_keep_purpose_rule(vroot, tmp_path):
    src = "from hexw.io import paths\ndef go(p, purpose):\n    return paths.open_data(p, 'rb', purpose=purpose).read(3)\n"
    for pkg in ("hexw.ingest",):
        with temp_module(tmp_path, pkg, "r", src) as m:
            assert m.go(vfile(vroot), "integrity") == b"1 2"
            assert m.go(vfile(vroot), "metadata") == b"1 2"
            with pytest.raises(paths.AccessError):
                m.go(vfile(vroot), "analysis")
            with pytest.raises(paths.AccessError):
                m.go(vfile(vroot), None)


def test_hexw_io_paths_itself_does_not_grant_access_to_pool_threads(vroot):
    with cf.ThreadPoolExecutor(1) as ex:
        with pytest.raises(paths.AccessError):
            ex.submit(paths.open_data, vfile(vroot), "rb", purpose="metadata").result()
        with pytest.raises(paths.AccessError):
            ex.submit(paths.sha256_file, vfile(vroot)).result()


# ------------------------------------------------------------------ fail-closed hook
def test_hook_failure_is_fail_closed_only_for_restricted_paths(vroot, monkeypatch, tmp_path):
    other = tmp_path / "plain.txt"
    other.write_text("x")

    def boom(*a, **k):
        raise RuntimeError("simulated hook failure")
    with paths.test_access():
        monkeypatch.setattr(paths, "_scan_stack", boom)
        with pytest.raises(paths.AccessError, match="fail-closed"):
            open(vfile(vroot), "rb")
        assert open(other, "rb").read() == b"x"                # not a vault path: scan is never needed
        monkeypatch.undo()
        monkeypatch.setattr(paths, "_matched_bases", boom)      # classification itself fails
        with pytest.raises(paths.AccessError, match="fail-closed"):
            open(vfile(vroot), "rb")
        with pytest.raises(paths.AccessError, match="fail-closed"):
            os.listdir(os.path.dirname(vfile(vroot)))
        assert open(other, "rb").read() == b"x"                 # no vault marker in the arguments


def test_hook_failure_in_frame_identity_denies(vroot, monkeypatch):
    def boom(*a, **k):
        raise ValueError("x")
    with paths.test_access():
        monkeypatch.setattr(paths, "_frame_role", boom)
        with pytest.raises(paths.AccessError):
            open(vfile(vroot), "rb")


# ------------------------------------------------------------------ configure / root
def test_configure_other_root_forbidden_outside_test_context(vroot, tmp_path):
    other = tmp_path / "otherroot"
    (other / "vault").mkdir(parents=True)
    code = textwrap.dedent('''
        import os
        from hexw.io import paths
        paths.configure(os.environ["HEXW_ROOT"])             # same root: no-op
        try:
            paths.configure(os.environ["OTHER"])
            print("CHANGED")
        except paths.AccessError:
            print("refused")
    ''')
    r = subprocess.run([sys.executable, "-c", code], env=dict(os.environ, HEXW_ROOT=str(vroot), OTHER=str(other)),
                       capture_output=True, text=True, cwd=str(tmp_path))
    assert r.stdout.strip() == "refused", r.stderr


def test_configure_from_theory_denied(vroot, theory_mod, tmp_path):
    other = tmp_path / "o2"
    (other / "vault").mkdir(parents=True)
    with pytest.raises(paths.AccessError):
        theory_mod.t_configure(str(other / "vault" / "x"))
    src = "from hexw.io import paths\ndef go(r):\n    with paths.test_root(r):\n        pass\n"
    with temp_module(tmp_path, "hexw.theory", "sneaky", src) as m:
        with pytest.raises(paths.TheoryAccessError):
            m.go(other)
    assert paths.get_config().root == vroot


def test_root_without_vault_is_an_explicit_error(tmp_path):
    with pytest.raises(paths.RootNotFoundError, match="HEXW_ROOT"):
        with paths.test_root(tmp_path / "nowhere"):
            pass


# ------------------------------------------------------------------ path spellings (classification via identity)
def _variants(r) -> dict:
    p = vfile(r)
    v = {"plain": p, "forward": p.replace("\\", "/"), "upper": p.upper(), "lower": p.lower(),
         "double_sep": p.replace("\\", "\\\\"),
         "dotdot": str(r / "raw" / "theory_side" / ".." / ".." / "vault" / "lab_validation" / "01_XRD" / "s.txt"),
         "trailing_dot_dir": p.replace("vault", "vault.")}
    if NT:
        drive = os.path.splitdrive(p)[0]
        v["lowerdrive"] = p[0].lower() + p[1:]
        v["extended"] = "\\\\?\\" + p
        v["extended_upper"] = "\\\\?\\" + p.upper()
        v["unc_admin"] = "\\\\localhost\\" + drive[0] + "$" + p[2:]
        v["unc_ext"] = "\\\\?\\UNC\\localhost\\" + drive[0] + "$" + p[2:]
        v["ads"] = p + "::$DATA"
        buf = ctypes.create_unicode_buffer(1024)
        if ctypes.windll.kernel32.GetShortPathNameW(os.path.dirname(p), buf, 1024):
            v["short83_dir"] = os.path.join(buf.value, "s.txt")
    return v


def test_all_spellings_of_a_vault_file_are_restricted_and_denied(vroot, theory_mod):
    ran = 0
    for name, p in _variants(vroot).items():
        if not os.path.exists(p):
            continue                        # spelling not resolvable on this machine (e.g. admin share disabled)
        ran += 1
        assert paths.classify(p) is paths.Zone.RESTRICTED, name
        with pytest.raises(paths.AccessError):
            theory_mod.t_open(p)
        with pytest.raises(paths.AccessError):
            open(p, "rb")                    # non-whitelisted test code, builtin open
        with pytest.raises(paths.ReadOnlyZoneError):
            paths.open_data(p, "wb")
    assert ran >= 6


@pytest.mark.skipif(not NT, reason="Windows path syntax")
def test_extended_and_unc_prefixes_classified_without_existing_file(vroot):
    p = vfile(vroot)
    assert paths.classify("\\\\?\\" + p) is paths.Zone.RESTRICTED
    assert paths.classify("\\\\?\\" + p.upper()) is paths.Zone.RESTRICTED
    assert paths.classify("\\\\localhost\\" + p[0] + "$" + p[2:]) is paths.Zone.RESTRICTED
    assert paths.classify("\\\\?\\UNC\\localhost\\" + p[0] + "$" + p[2:]) is paths.Zone.RESTRICTED
    assert paths.classify(str(vroot / "raw" / "theory_side" / "lit.txt")) is paths.Zone.RAW_THEORY
    assert paths.classify(str(vroot / "registry" / "x")) is paths.Zone.REGISTRY


def test_relative_path_from_other_cwd(vroot, theory_mod, monkeypatch):
    monkeypatch.chdir(vroot / "raw" / "theory_side")
    rel = os.path.join("..", "..", "vault", "lab_validation", "01_XRD", "s.txt")
    with pytest.raises(paths.AccessError):
        theory_mod.t_open(rel)
    monkeypatch.chdir(vroot / "vault" / "lab_validation")
    with pytest.raises(paths.AccessError):
        theory_mod.t_open(os.path.join("01_XRD", "s.txt"))


def _mklink_junction(link: Path, target: Path):
    with paths.test_access():       # the helper process receives a vault path
        r = subprocess.run(["cmd", "/c", "mklink", "/J", str(link), str(target)], capture_output=True, text=True)
    if r.returncode != 0:
        pytest.skip(f"cannot create a junction here: {r.stdout.strip()} {r.stderr.strip()}")


@pytest.mark.skipif(not NT, reason="junctions are a Windows feature")
@pytest.mark.parametrize("target_rel", [("vault",), ("vault", "lab_validation"), ("vault", "lab_validation", "01_XRD")])
def test_junction_into_vault_is_restricted(vroot, theory_mod, tmp_path, target_rel):
    link = tmp_path / "jv"
    _mklink_junction(link, vroot.joinpath(*target_rel))
    inner = link.joinpath(*(("lab_validation", "01_XRD", "s.txt")[len(target_rel) - 1:]))
    assert os.path.exists(inner)
    assert paths.classify(inner) is paths.Zone.RESTRICTED
    with pytest.raises(paths.AccessError):
        theory_mod.t_open(str(inner))
    with pytest.raises(paths.AccessError):
        open(inner, "rb")
    with pytest.raises(paths.AccessError):
        os.listdir(inner.parent)
    with pytest.raises(paths.ReadOnlyZoneError):
        open(inner.parent / "new.txt", "wb")                   # write through the junction


def test_symlink_into_vault_is_restricted(vroot, theory_mod, tmp_path):
    link = tmp_path / "sl.txt"
    try:
        os.symlink(vfile(vroot), link)
    except OSError as e:
        pytest.skip(f"symlinks need privileges here: {e}")
    assert paths.classify(link) is paths.Zone.RESTRICTED
    with pytest.raises(paths.AccessError):
        theory_mod.t_open(str(link))


def test_write_modes_to_read_only_zones_denied_for_everybody(vroot):
    for p in (vfile(vroot), str(vroot / "raw" / "theory_side" / "lit.txt"), str(vroot / "vault" / "new.bin")):
        with paths.test_access():
            for mode in ("wb", "ab", "r+b", "xb"):
                with pytest.raises(paths.ReadOnlyZoneError):
                    open(p, mode)
            with pytest.raises(paths.ReadOnlyZoneError):
                os.open(p, os.O_WRONLY | getattr(os, "O_BINARY", 0))
