"""Stage 3 round 1: child-process arguments in EVERY spelling of a protected path (HOLDOUT and VAULT).

Regression of the independent tester (`cmd /c type P`, `python -c open(P)`, PowerShell Get-Content with
P = `...\\literature\\..\\lit_search\\...`, `\\\\localhost\\D$\\...`, `\\\\?\\UNC\\localhost\\D$\\...`, junction).
On Windows subprocess reports the ALREADY JOINED command line, so the hook scans the text for path fragments.
Everything is synthetic (tmp_path); the child is never started when the hook refuses."""
from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest

from _helpers import run_as
from hexw.io import paths

NT = os.name == "nt"
SRC = (
    "import os, subprocess, sys\n"
    "def run(a, **kw): return subprocess.run(a, capture_output=True, **kw)\n"
    "def chdir_run(d, a):\n"
    "    os.chdir(d)\n"
    "    return subprocess.run(a, capture_output=True)\n"
    "def system(c): return os.system(c)\n"
)


@pytest.fixture
def sroot(tmp_path):
    for rel in ("raw/theory_side/06_literature/L01", "raw/theory_side/06_literature/L02", "raw/original_zip",
                "vault/lab_validation/01_XRD", "registry", "lit_search/downloads/C056", "literature/fit",
                "literature_holdout", "project/config", "dft_feasibility/pilots/c_W92", "work"):
        (tmp_path / rel).mkdir(parents=True)
    (tmp_path / "project/config/holdout_sources.json").write_text(json.dumps({"holdout_sources": [
        "raw/theory_side/06_literature/L01", "lit_search/downloads/C056",
        "dft_feasibility/pilots/c_W92/SrZn2Fe16O27_model.cif"]}))
    for rel in ("raw/theory_side/06_literature/L01/a.cif", "raw/theory_side/06_literature/L02/p.pdf",
                "lit_search/downloads/C056/C056.pdf", "literature_holdout/records_holdout.yaml",
                "vault/lab_validation/01_XRD/s.txt", "dft_feasibility/pilots/c_W92/SrZn2Fe16O27_model.cif",
                "dft_feasibility/pilots/c_W92/ok.in", "literature/fit/records_fit.yaml"):
        (tmp_path / rel).write_text("x")
    with paths.test_root(tmp_path):
        yield tmp_path


TARGETS = {
    "holdout_src_file": "raw/theory_side/06_literature/L01/a.cif",
    "holdout_dl_file": "lit_search/downloads/C056/C056.pdf",
    "holdout_dir_file": "literature_holdout/records_holdout.yaml",
    "pilot_derivative": "dft_feasibility/pilots/c_W92/SrZn2Fe16O27_model.cif",
    "vault_file": "vault/lab_validation/01_XRD/s.txt",
    "vault_missing_file": "vault/lab_validation/01_XRD/not_there.txt",
}


def _unc_admin(p: str) -> str:
    return "\\\\localhost\\" + p[0] + "$" + p[2:]


def spellings(root, rel):
    p = str(root / rel).replace("/", "\\")
    d = p[0]
    out = {
        "plain": p,
        "forward": p.replace("\\", "/"),
        "lower": p.lower(),
        "dotdot": str(root / "literature" / ".." / rel).replace("/", "\\"),
        "dotdot_deep": str(root / "work" / ".." / "work" / ".." / rel).replace("/", "\\"),
        "extended": "\\\\?\\" + p,
        "unc_admin": _unc_admin(p),
        "unc_extended": "\\\\?\\UNC\\localhost\\" + d + "$" + p[2:],
        "rooted_no_drive": p[2:],
    }
    return out


def _skip_if_no_admin_share(root):
    probe = _unc_admin(str(root))
    try:
        os.stat(probe)
    except OSError:
        pytest.skip("administrative share of the drive is not reachable here")


def _cmd_forms(p):
    """(label, argv-or-string) for the typical readers."""
    return [
        ("cmd_type_list", ["cmd", "/c", "type", p]),
        ("cmd_type_str", 'cmd /c type "%s"' % p),
        ("python_c_literal", [sys.executable, "-c", "print(open(r'%s','rb').read(8))" % p]),
        ("python_c_escaped", [sys.executable, "-c", "print(open('%s','rb').read(8))" % p.replace("\\", "\\\\")]),
        ("powershell", ["powershell", "-NoProfile", "-Command", "Get-Content -LiteralPath '%s' -TotalCount 1" % p]),
        ("option_eq", [sys.executable, "-c", "pass", "--in=" + p]),
        ("option_colon_quoted", 'tool.exe /in:"%s" /x' % p),
        ("redirect", 'cmd /c "type %s > nul"' % p),
    ]


@pytest.mark.skipif(not NT, reason="Windows path spellings")
@pytest.mark.parametrize("target", sorted(TARGETS))
@pytest.mark.parametrize("spelling", ["plain", "forward", "lower", "dotdot", "dotdot_deep", "extended",
                                      "unc_admin", "unc_extended", "rooted_no_drive"])
def test_child_process_arguments_in_every_spelling_are_blocked(sroot, target, spelling):
    if spelling.startswith("unc"):
        _skip_if_no_admin_share(sroot)
    if spelling == "rooted_no_drive" and os.path.splitdrive(os.getcwd())[0].lower() != str(sroot)[:2].lower():
        pytest.skip("the drive-relative spelling needs the cwd on the same drive")
    p = spellings(sroot, TARGETS[target])[spelling]
    s = run_as("hexw.theory._probe", SRC)
    for label, cmd in _cmd_forms(p):
        with pytest.raises(paths.AccessError):
            s["run"](cmd, shell=isinstance(cmd, str))


@pytest.mark.skipif(not NT, reason="Windows path spellings")
@pytest.mark.parametrize("target", ["holdout_src_file", "vault_file", "pilot_derivative"])
def test_os_system_blocked(sroot, target):
    p = str(sroot / TARGETS[target])
    s = run_as("hexw.theory._probe", SRC)
    for q in (p, str(sroot / "literature" / ".." / TARGETS[target])):
        with pytest.raises(paths.AccessError):
            s["system"]('type "%s" > nul' % q)


@pytest.mark.parametrize("target", ["holdout_src_file", "holdout_dl_file", "holdout_dir_file", "vault_file",
                                    "vault_missing_file"])
def test_relative_paths_cwd_and_chdir_blocked(sroot, target):
    s = run_as("hexw.theory._probe", SRC)
    rel = TARGETS[target]
    here = sroot / "work"
    old = os.getcwd()
    try:
        os.chdir(here)
        up = os.path.join("..", *rel.split("/"))
        with pytest.raises(paths.AccessError):                               # relative path from the cwd
            s["run"](["cmd", "/c", "type", up])
        with pytest.raises(paths.AccessError):
            s["run"]([sys.executable, "-c", "print(open(r'%s','rb').read(8))" % up])
        with pytest.raises(paths.AccessError):                               # cwd keyword = protected directory
            s["run"]([sys.executable, "-c", "pass"], cwd=str((sroot / rel).parent))
        with pytest.raises(paths.AccessError):                               # relative cwd keyword
            s["run"]([sys.executable, "-c", "pass"], cwd=os.path.dirname(up))
        os.chdir(old)
        with pytest.raises(paths.AccessError):                               # the process cwd is protected itself
            s["chdir_run"](str((sroot / rel).parent), [sys.executable, "-c", "pass"])
    finally:
        os.chdir(old)


@pytest.mark.skipif(not NT, reason="mklink /J")
@pytest.mark.parametrize("target", ["holdout_src_file", "holdout_dl_file", "vault_file", "pilot_derivative"])
def test_junction_blocked(sroot, tmp_path_factory, target):
    link = tmp_path_factory.mktemp("jn") / "link"
    tdir = (sroot / TARGETS[target]).parent
    with paths.test_access():                      # building the fixture is test infrastructure
        r = subprocess.run(["cmd", "/c", "mklink", "/J", str(link), str(tdir)], capture_output=True, text=True)
    if r.returncode:
        pytest.skip("cannot create a junction: " + r.stderr)
    p = str(link / (sroot / TARGETS[target]).name)
    s = run_as("hexw.theory._probe", SRC)
    for label, cmd in _cmd_forms(p):
        with pytest.raises(paths.AccessError):
            s["run"](cmd, shell=isinstance(cmd, str))


def test_harmless_children_still_run(sroot):
    s = run_as("hexw.theory._probe", SRC)
    r = s["run"]([sys.executable, "-c", "print(1+1)"])
    assert r.returncode == 0 and r.stdout.strip() == b"2"
    ok = str(sroot / "dft_feasibility/pilots/c_W92/ok.in")
    assert s["run"]([sys.executable, "-c", "pass", ok, "--in=" + ok]).returncode == 0
    r = s["run"]([sys.executable, "-c", "print('a/b', 1/2)"], cwd=str(sroot / "work"))
    assert r.returncode == 0
    fit = str(sroot / "literature" / "fit" / "records_fit.yaml")
    assert s["run"]([sys.executable, "-c", "pass", fit]).returncode == 0
    assert s["run"]([sys.executable, "-c", "print('http://example.com/x', 'http://localhost/vault')"]).returncode == 0


def test_whitelisted_callers_and_test_access_still_pass(sroot):
    v = str(sroot / TARGETS["vault_file"])
    with paths.test_access():
        assert subprocess.run([sys.executable, "-c", "pass", v]).returncode == 0
    aud = run_as("hexw.audit._probe", SRC)
    assert aud["run"]([sys.executable, "-c", "pass", v]).returncode == 0           # audit: from_hook decision
    h = str(sroot / TARGETS["holdout_src_file"])
    assert aud["run"]([sys.executable, "-c", "pass", h]).returncode == 0
    th = run_as("hexw.theory._probe", SRC)
    with pytest.raises(paths.TheoryAccessError):
        th["run"]([sys.executable, "-c", "pass", h])
    with pytest.raises(paths.TheoryAccessError):
        th["run"]([sys.executable, "-c", "pass", v])


@pytest.mark.parametrize("pkg", ["theory", "params", "dft", "structures", "spin", "util"])
def test_all_theory_side_packages_blocked(sroot, pkg):
    s = run_as(f"hexw.{pkg}._probe", SRC)
    for t in ("holdout_src_file", "vault_file"):
        with pytest.raises(paths.TheoryAccessError):
            s["run"](["cmd", "/c", "type", str(sroot / "literature" / ".." / TARGETS[t])])


def test_too_many_fragments_fail_closed(sroot):
    s = run_as("hexw.theory._probe", SRC)
    junk = " ".join(f"a{i}/b{i}" for i in range(paths._PROC_MAX_CAND + 50))
    with pytest.raises(paths.AccessError):
        s["run"]([sys.executable, "-c", "pass", junk])
