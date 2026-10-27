"""Single data-access layer with roles (white list).

Zones:
  RAW_THEORY      raw/theory_side          readable by anyone (literature, docs, manifests)
  RESTRICTED      vault/** and raw/original_zip (the zip still contains laboratory bytes)
  RAW_OTHER       rest of raw/ (read-only)
  HOLDOUT         literature_holdout/ and the holdout source files declared in project/config/holdout_sources.json
                  (CIF/PDF/XML of the holdout publications inside raw/theory_side and lit_search/downloads; read-only)
  REGISTRY, OTHER everything else

Policy for RESTRICTED (WHITE LIST - absence of a prohibition is NOT a permission):
  * access is allowed only if the call stack contains a frame that provably belongs to one of
    hexw.experiment, hexw.validate, hexw.audit, hexw.ingest, hexw.io (the paths module itself does not count),
    or the current thread is inside the explicit `test_access()` context manager (tests only);
  * a frame that belongs to any other hexw subpackage (hexw.theory, hexw.params, hexw.dft, hexw.structures,
    hexw.spin, any future one) anywhere in the stack -> TheoryAccessError (confused deputy: numpy.loadtxt,
    pandas.read_csv, pathlib, helpers from other packages called from theory code are covered because
    the whole stack is inspected);
  * a thread without any allowed frame (ThreadPoolExecutor / threading.Thread started from theory) -> AccessError;
  * the innermost allowed frame decides the extra condition: experiment/validate need the freeze file,
    audit/ingest/io need purpose 'integrity' or 'metadata' in open_data();
  * any write mode to raw/ or vault/ -> ReadOnlyZoneError.

Policy for HOLDOUT (same hook, separate zone; the literature holdout must stay unseen by the theory until freeze):
  * any frame of a non-white-listed hexw subpackage (theory, params, dft, structures, spin, util, ...) in the stack
    -> TheoryAccessError, for open / listdir / scandir / glob / copy / sqlite / child-process arguments alike;
  * the innermost allowed frame decides: hexw.validate (only after the freeze file exists) or hexw.audit (integrity /
    metadata). hexw.ingest only inside `holdout_build()` (the reproducible literature split script that builds the
    zone; it is the only writer). hexw.experiment / hexw.io alone are NOT enough;
  * holdout is a read-only zone: writes are refused (ReadOnlyZoneError) except inside `holdout_build()`;
  * the zone has priority over RAW_THEORY (a holdout CIF lies physically inside raw/theory_side).

Frame identity (_frame_role): `__name__`, `__spec__.name`, `__package__` and the source file name (co_filename
inside the hexw package directory) are all taken into account. Conflicting claims -> denied. So
`python -m hexw.theory.x`, runpy `__main__`, a module of hexw.dft loaded from a directory added to the package
`__path__` are recognised. Frames without any hexw identity are "external" (never an allowed frame).

Path classification uses the identity of directories (os.stat / os.path.samestat of every ancestor of the
absolute and of the real path), so the extended prefix, UNC admin shares, UNC extended form, letter case,
junctions, symlinks and 8.3 names are all recognised.

The process-wide audit hook (installed on `import hexw`, cannot be removed) guards builtin open() / os.open and the
side channels: os.listdir, os.scandir, glob, subprocess.Popen (+ os.system, os.exec*, os.spawn*, os.startfile),
shutil copy/copy2/copyfile/copytree/move/copymode/copystat, the native `_winapi.CopyFile2` event (shutil.copy2 on
Python 3.13/Windows emits NO `open` event: the source path is seen only through `_winapi.CopyFile2`), and
sqlite3.connect (plain path and `file:` URI: percent-decoding, authority, cut at `?`/`#` - the truncated name is
what SQLite really creates). The destination of copy/move and any non-`mode=ro` sqlite open count as writes.
It is FAIL-CLOSED: an exception inside the check denies access when the event touches (or cannot be proven not to
touch) vault/ or raw/original_zip. A missing root or a root without vault/ is an explicit error (set HEXW_ROOT).

Child processes: every process started by a hexw process gets HEXW_GUARD=1 and HEXW_ROOT=<current root> (os.environ is
set on import; explicit env dicts are completed in place by the hook). A caller without a white-listed frame that
passes a DIFFERENT HEXW_ROOT (env dict, any letter case, or a poisoned os.environ) is refused. multiprocessing 'spawn'
(also ProcessPoolExecutor): `multiprocessing.spawn.get_command_line` is MONKEY-PATCHED at import so that the child
sets HEXW_ROOT/HEXW_GUARD and imports hexw before spawn_main, i.e. before any user code. A command line / argument
containing a vault or HOLDOUT path is denied (parent side): the joined command line / each argument / `-c` code is
scanned for path fragments (`..`, UNC admin share, extended prefix, junction, relative to the child's cwd, `--in=path`),
made absolute and classified by directory identity (`_proc_zone_hits`; docs/ISOLATION.md). A path the child computes
itself is not caught.

Fast path: `open` of a path that the precise classifier put outside vault/original_zip is cached (exact string,
revalidated by the identity of the parent directory, lstat of the last component and a 1 s TTL); stack inspection
happens only for restricted paths. ~50 us instead of ~850 us per open (base: ~34 us).

PROTECTION AGAINST MISTAKES, NOT A SANDBOX. NOT covered (see docs/ISOLATION.md): ctypes / native extensions calling
CreateFileW, `_winapi.CreateFile` (no audit event), hard links, mmap of an already open handle, external programs
and any child process that does not import hexw (only the parent-side Popen argument check applies; the DFT adapter
builds its inputs itself, so this is accepted), metadata calls (os.stat/exists/getsize - names and sizes are
visible), `Path.rglob` / `os.walk` silently yielding nothing on a denied directory (the OSError is swallowed by
the library; direct calls raise), a directory moved into the vault within the 1 s cache TTL.
"""
from __future__ import annotations

import contextlib
import os
import re
import socket
import stat as _stat_mod
import sys
import threading
from time import monotonic as _monotonic
from urllib.parse import unquote as _unquote
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

# ---- frame-origin policy ----
WHITELIST_SUBS = frozenset({"experiment", "validate", "audit", "ingest", "io"})
NEUTRAL_SUBS = frozenset({"registry", "chem", "units", "cli", "ui", "envinfo", "__main__", "__init__"})
# ui is a routing layer, like cli; a theory frame below it still vetoes target access.
# every other hexw subpackage (theory, params, dft, structures, spin, ...) is DENIED by construction.
TEST_ROLE = "test"
_SELF = __name__
_MY_GLOBALS = globals()
_PKG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class AccessError(PermissionError):
    pass


class TheoryAccessError(AccessError):
    """Theory/parametrization/DFT/... code (or a deputy called by it) tried to touch the laboratory vault."""


class FreezeRequiredError(AccessError):
    """Vault access requires a frozen model (freeze file)."""


class ReadOnlyZoneError(AccessError):
    """Write access to raw/ or vault/."""


class RootNotFoundError(AccessError, RuntimeError):
    """Project root (with vault/) not found: set the HEXW_ROOT environment variable."""


class Zone(str, Enum):
    RAW_THEORY = "raw_theory"
    RESTRICTED = "restricted"
    RAW_OTHER = "raw_other"
    HOLDOUT = "holdout"
    REGISTRY = "registry"
    OTHER = "other"


HOLDOUT_CONFIG_REL = ("project", "config", "holdout_sources.json")


@dataclass(frozen=True)
class Config:
    root: Path
    # absolute paths (files or directories) of holdout source publications, declared in HOLDOUT_CONFIG_REL
    holdout_sources: tuple = ()

    @property
    def holdout_config_file(self): return self.root.joinpath(*HOLDOUT_CONFIG_REL)
    @property
    def literature(self): return self.root / "literature"
    @property
    def literature_fit(self): return self.literature / "fit"
    @property
    def literature_holdout(self): return self.root / "literature_holdout"
    @property
    def raw(self): return self.root / "raw"
    @property
    def raw_theory(self): return self.raw / "theory_side"
    @property
    def original_zip_dir(self): return self.raw / "original_zip"
    @property
    def vault(self): return self.root / "vault"
    @property
    def vault_lab(self): return self.vault / "lab_validation"
    @property
    def vault_legacy(self): return self.vault / "legacy"
    @property
    def vault_lab_umnik(self): return self.vault / "lab_validation_umnik"
    @property
    def prior_analyses_umnik(self): return self.vault_legacy / "prior_analyses_umnik"
    @property
    def registry(self): return self.root / "registry"
    @property
    def freeze_file(self): return self.registry / "FROZEN_MODEL.json"


# ---------------- thread-local explicit test permissions ----------------
_TL = threading.local()


def _depth(name: str) -> int:
    return getattr(_TL, name, 0)


# ---------------- configuration ----------------
_CFG: Config | None = None


def _root_error(why: str) -> RootNotFoundError:
    return RootNotFoundError(
        f"{why}. Set the HEXW_ROOT environment variable to the project root "
        "(the directory that contains vault/, raw/ and registry/). The data guard is never silently disabled.")


def _load_holdout_sources(root: Path) -> tuple:
    """Holdout source paths declared in <root>/project/config/holdout_sources.json (list under "holdout_sources",
    relative to the root). Missing file: no declared sources, unless the root has a literature_holdout/ directory
    (then the declaration is required). Any malformed declaration is an error: the guard never silently weakens."""
    import json
    f = root.joinpath(*HOLDOUT_CONFIG_REL)
    if not f.is_file():
        if (root / "literature_holdout").exists():
            raise _root_error(f"{str(f)!r} is missing but literature_holdout/ exists (holdout sources undeclared)")
        return ()
    try:
        doc = json.loads(f.read_text(encoding="utf-8"))
        rels = doc["holdout_sources"]
        if not isinstance(rels, list) or not all(isinstance(r, str) and r for r in rels):
            raise ValueError("holdout_sources must be a list of non-empty strings")
        out = []
        for r in rels:
            rp = Path(r)
            if rp.is_absolute() or ".." in rp.parts:
                raise ValueError(f"holdout source {r!r} must be relative to the root and stay inside it")
            out.append(root / rp)
    except (OSError, ValueError, KeyError, TypeError) as e:
        raise _root_error(f"cannot read the holdout declaration {str(f)!r}: {type(e).__name__}: {e}") from e
    return tuple(out)


def _make_config(root, *, validate: bool = True) -> Config:
    root = Path(root)
    cfg = Config(root)
    if validate and not cfg.vault.is_dir():
        raise _root_error(f"root {str(cfg.root)!r} has no vault/ directory")
    return Config(root, _load_holdout_sources(root))


def _default_root() -> Path:
    env = os.environ.get("HEXW_ROOT")
    if env:
        return Path(env)
    if os.environ.get("HEXW_GUARD") == "1":
        # a guarded child process (spawned by a hexw process) must be told the root explicitly
        raise _root_error("HEXW_GUARD=1 but HEXW_ROOT is not set (guarded child process)")
    return Path(__file__).resolve().parents[3]


def get_config() -> Config:
    if _CFG is None:
        raise _root_error("no valid project root is configured")
    return _CFG


def _same_dir(a, b) -> bool:
    try:
        return os.path.samestat(os.stat(a), os.stat(b))
    except OSError:
        return os.path.normcase(os.path.abspath(a)) == os.path.normcase(os.path.abspath(b))


def configure(root: str | os.PathLike) -> Config:
    """Re-point the data layer. A DIFFERENT root is allowed only inside `test_root()` (tests)."""
    global _CFG
    if _CFG is not None and _same_dir(root, _CFG.root):
        return _CFG
    if not _depth("root"):
        raise AccessError("paths.configure() with a different root is forbidden outside the test context "
                          "(set the HEXW_ROOT environment variable before starting the process)")
    _CFG = _make_config(root)
    return _CFG


@contextlib.contextmanager
def holdout_build():
    """Construction of literature_holdout/ by hexw.ingest (the reproducible split script) ONLY: lifts the read-only
    rule of the HOLDOUT zone for this thread. Refused for theory-side frames and for anything outside hexw.ingest."""
    denied, allowed = _scan_stack()
    if denied:
        raise TheoryAccessError(f"holdout_build() is not available to theory-side code: {denied}")
    if not any(sub == "ingest" for sub, _l in allowed):
        raise AccessError("holdout_build() may be called only from hexw.ingest")
    _TL.hbuild = _depth("hbuild") + 1
    try:
        yield
    finally:
        _TL.hbuild = _depth("hbuild") - 1


def _forbid_from_theory(what: str):
    denied, _allowed = _scan_stack()
    if denied:
        raise TheoryAccessError(f"{what} is not available to theory-side code: {denied}")
    if "pytest" not in sys.modules and not os.environ.get("HEXW_TESTING"):
        raise AccessError(f"{what} is a test-only facility (pytest not loaded, HEXW_TESTING not set)")


@contextlib.contextmanager
def test_root(root, *, validate: bool = True):
    """TESTS ONLY: temporarily point the data layer at another root. Grants NO vault access."""
    global _CFG
    _forbid_from_theory("test_root()")
    old = _CFG
    _CFG = _make_config(root, validate=validate)
    _TL.root = _depth("root") + 1
    try:
        yield _CFG
    finally:
        _TL.root = _depth("root") - 1
        _CFG = old


@contextlib.contextmanager
def test_access():
    """TESTS ONLY, explicit and per thread: the calling thread counts as test infrastructure and may open
    vault/. Other threads (e.g. a ThreadPoolExecutor started inside) do NOT inherit it; theory-side frames in the
    stack are still denied."""
    _forbid_from_theory("test_access()")
    _TL.access = _depth("access") + 1
    try:
        yield
    finally:
        _TL.access = _depth("access") - 1


# ---------------- path normalisation / classification ----------------
def _to_str(p) -> str | None:
    if isinstance(p, int) or p is None:
        return None
    p = os.fspath(p)
    return os.fsdecode(p) if isinstance(p, (bytes, bytearray)) else p


def _strip_ns(s: str) -> str:
    r"""Extended prefix, UNC extended form and admin-share UNC -> plain drive path (lexical help only)."""
    if os.name != "nt":
        return s
    t = s.replace("/", "\\")
    low = t.lower()
    for pre in ("\\\\?\\unc\\", "\\\\.\\unc\\"):
        if low.startswith(pre):
            t = "\\\\" + t[len(pre):]
            low = t.lower()
            break
    else:
        if low.startswith(("\\\\?\\", "\\\\.\\")):
            t = t[4:]
            low = t.lower()
    if low.startswith("\\\\"):
        parts = t[2:].split("\\", 2)
        if (len(parts) >= 2 and len(parts[1]) == 2 and parts[1][1] == "$" and parts[1][0].isalpha()
                and parts[0].lower() in {"localhost", "127.0.0.1", "::1",
                                         os.environ.get("COMPUTERNAME", "?").lower()}):
            t = parts[1][0] + ":\\" + (parts[2] if len(parts) > 2 else "")
    return t


def _lex(s: str) -> str:
    return os.path.normcase(_strip_ns(s))


def _forms(s: str) -> list[str]:
    out = []
    for fn in (os.path.abspath, os.path.realpath):
        try:
            v = fn(s)
        except (OSError, ValueError):
            continue
        if v not in out:
            out.append(v)
    return out


def _stat(p):
    try:
        return os.stat(p)
    except (OSError, ValueError):
        return None


def _base_table(c: Config):
    """name -> (lexical forms, stat or None)."""
    table = {}
    for name, d in (("vault", c.vault), ("original_zip", c.original_zip_dir), ("raw_theory", c.raw_theory),
                    ("raw", c.raw), ("registry", c.registry), ("holdout_dir", c.literature_holdout)):
        s = os.fspath(d)
        table[name] = ([_lex(f) for f in _forms(s)], _stat(s))
    for i, d in enumerate(c.holdout_sources):
        s = os.fspath(d)
        table[f"holdout_src{i}"] = ([_lex(f) for f in _forms(s)], _stat(s))
    return table


def _inside_lex(p: str, base: str) -> bool:
    return p == base or p.startswith(base.rstrip("\\/") + os.sep)


def _matched_bases(path: str, c: Config) -> set[str]:
    table = _base_table(c)
    hit: set[str] = set()
    for form in _forms(path) or [path]:
        lf = _lex(form)
        for name, (lexes, _st) in table.items():
            if any(_inside_lex(lf, b) for b in lexes):
                hit.add(name)
        cur = form
        for _ in range(96):
            st = _stat(cur)
            if st is not None:
                for name, (_lx, bst) in table.items():
                    if bst is not None and os.path.samestat(st, bst):
                        hit.add(name)
            parent = os.path.dirname(cur)
            if not parent or parent == cur:
                break
            cur = parent
    return hit


def _zone_of(hit: set[str]) -> Zone:
    if "vault" in hit or "original_zip" in hit:
        return Zone.RESTRICTED
    if any(h == "holdout_dir" or h.startswith("holdout_src") for h in hit):
        return Zone.HOLDOUT
    if "raw_theory" in hit:
        return Zone.RAW_THEORY
    if "raw" in hit:
        return Zone.RAW_OTHER
    if "registry" in hit:
        return Zone.REGISTRY
    return Zone.OTHER


def classify(path) -> Zone:
    s = _to_str(path)
    if s is None:
        return Zone.OTHER
    return _zone_of(_matched_bases(s, get_config()))


_RO_ZONES = (Zone.RESTRICTED, Zone.RAW_THEORY, Zone.RAW_OTHER, Zone.HOLDOUT)


# ---------------- frame identity ----------------
_FILE_SUB_CACHE: dict[str, str | None] = {}


def _file_sub(fn: str) -> str | None:
    """Subpackage (or top-level module) of hexw the source file lives in; None if outside the package."""
    if fn in _FILE_SUB_CACHE:
        return _FILE_SUB_CACHE[fn]
    res = None
    try:
        if fn and not fn.startswith("<"):
            rp = os.path.normcase(os.path.realpath(fn))
            base = os.path.normcase(os.path.realpath(_PKG_DIR))
            if rp.startswith(base + os.sep):
                part = rp[len(base) + 1:].split(os.sep)[0]
                res = part[:-3] if part.endswith(".py") else part
    except (OSError, ValueError):
        res = None
    _FILE_SUB_CACHE[fn] = res
    return res


def _sub_from_name(n) -> str | None:
    if not isinstance(n, str):
        return None
    if n == "hexw":
        return "__init__"
    if n.startswith("hexw."):
        return n.split(".")[1]
    return None


def _sub_role(sub: str) -> str:
    if sub in WHITELIST_SUBS:
        return "allow"
    if sub in NEUTRAL_SUBS:
        return "neutral"
    return "deny"


def _frame_role(f):
    """-> (role, label, sub); role in skip|external|neutral|allow|deny."""
    g = f.f_globals
    if g is _MY_GLOBALS:
        return "skip", _SELF, None
    spec = g.get("__spec__")
    pkg = g.get("__package__")
    claims = {_sub_from_name(g.get("__name__")), _sub_from_name(getattr(spec, "name", None)),
              _sub_from_name(pkg) if pkg != "hexw" else None}   # "hexw" as __package__ is only the container
    claims.discard(None)
    fsub = _file_sub(f.f_code.co_filename)
    label = str(g.get("__name__"))
    subs = set(claims)
    if fsub is not None:
        subs.add(fsub)
    if not subs:
        return "external", label, None
    if len(subs) > 1:
        return "deny", f"{label} (conflicting origin {sorted(subs)})", None
    (sub,) = subs
    return _sub_role(sub), (f"hexw.{sub}" if label == "__main__" else label), sub


def _scan_stack():
    """-> (denied labels, allowed [(sub, label)]), innermost first."""
    denied, allowed = [], []
    f = sys._getframe(1)
    while f is not None:
        role, label, sub = _frame_role(f)
        if role == "deny":
            denied.append(label)
        elif role == "allow":
            allowed.append((sub, label))
        f = f.f_back
    if _depth("access"):
        allowed.append((TEST_ROLE, "<test_access context>"))
    return denied, allowed


# ---------------- decisions ----------------
def is_frozen() -> bool:
    f = get_config().freeze_file
    if not f.is_file() or f.stat().st_size == 0:
        return False
    # Old access-policy tests use opaque markers in explicit synthetic test roots.
    # This compatibility facility can never be entered from a theory-side frame.
    if _depth("root"):
        import json
        try:
            if json.loads(f.read_text(encoding="utf-8")).get("schema") != "hexw.freeze.v1":
                return True
        except (ValueError, AttributeError):
            return True
    try:
        from hexw.registry.freeze import verify_manifest
        verify_manifest(f)
        return True
    except (OSError, ValueError, KeyError, TypeError):
        return False


def _is_write_mode(mode) -> bool:
    return isinstance(mode, str) and any(ch in mode for ch in "wax+")


def _is_write_flags(flags) -> bool:
    if not isinstance(flags, int):
        return False
    return bool(flags & (os.O_WRONLY | os.O_RDWR | os.O_APPEND | os.O_CREAT | os.O_TRUNC))


def _vault_decision(path, purpose=None, *, from_hook=False):
    denied, allowed = _scan_stack()
    if denied:
        raise TheoryAccessError(
            f"theory-side code (or a deputy called by it) must not access the vault ({path}); "
            f"offending frames: {denied}")
    if not allowed:
        raise AccessError(
            "no frame of hexw.experiment/validate/audit/ingest/io (or an explicit test context) in the call stack: "
            f"vault access denied ({path}); thread={threading.current_thread().name!r}")
    sub, label = allowed[0]
    if sub in ("experiment", "validate"):
        if not is_frozen():
            raise FreezeRequiredError(
                f"vault access for {label} requires a frozen model: {get_config().freeze_file} not found")
    elif sub in ("audit", "ingest", "io"):
        if not from_hook and purpose not in ("integrity", "metadata"):
            raise AccessError(f"module {label!r} may open the vault only for purpose 'integrity'/'metadata' "
                              f"(got {purpose!r})")


def _holdout_decision(path, purpose=None, *, from_hook=False):
    denied, allowed = _scan_stack()
    if denied:
        raise TheoryAccessError(
            f"theory-side code (or a deputy called by it) must not access the HOLDOUT zone ({path}); "
            f"offending frames: {denied}")
    if not allowed:
        raise AccessError(
            "no frame of hexw.validate/audit (or hexw.ingest inside holdout_build(), or an explicit test context) "
            f"in the call stack: HOLDOUT access denied ({path}); thread={threading.current_thread().name!r}")
    sub, label = allowed[0]
    if sub == "validate":
        if not is_frozen():
            raise FreezeRequiredError(
                f"HOLDOUT access for {label} requires a frozen model: {get_config().freeze_file} not found")
    elif sub == "audit":
        if not from_hook and purpose not in ("integrity", "metadata"):
            raise AccessError(f"module {label!r} may open the HOLDOUT zone only for purpose 'integrity'/'metadata' "
                              f"(got {purpose!r})")
    elif sub == "ingest":
        if not _depth("hbuild"):
            raise AccessError(f"module {label!r} may touch the HOLDOUT zone only inside paths.holdout_build() "
                              f"({path})")
    elif sub != TEST_ROLE:
        raise AccessError(f"module {label!r} is not allowed to read the HOLDOUT zone ({path})")


def _ro_violation(z: Zone, is_write: bool) -> bool:
    """Write into a read-only zone. The HOLDOUT zone may be written inside holdout_build() (then the normal
    holdout decision still applies)."""
    if not (is_write and z in _RO_ZONES):
        return False
    return not (z is Zone.HOLDOUT and _depth("hbuild"))


def check_access(path, mode: str = "rb", purpose: str | None = None) -> Zone:
    """Raise if the current call stack may not use `path`. Returns the zone."""
    z = classify(path)
    if _ro_violation(z, _is_write_mode(mode)):
        raise ReadOnlyZoneError(f"write access to read-only zone is forbidden: {path}")
    if z is Zone.RESTRICTED:
        _vault_decision(path, purpose)
    elif z is Zone.HOLDOUT:
        _holdout_decision(path, purpose)
    return z


def open_data(path, mode: str = "rb", *, purpose: str | None = None, **kw):
    """The only sanctioned way to open data files."""
    check_access(path, mode, purpose)
    return open(path, mode, **kw)


def read_prefix(path, nbytes: int, *, purpose: str = "metadata") -> bytes:
    """Read at most nbytes from the start of a file (used for header-only metadata)."""
    with open_data(path, "rb", purpose=purpose) as f:
        return f.read(nbytes)


def sha256_file(path, *, purpose: str = "integrity", chunk: int = 1 << 20) -> str:
    import hashlib
    h = hashlib.sha256()
    with open_data(path, "rb", purpose=purpose) as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


# ---------------- process-wide audit hook (always on, fail-closed) ----------------
_HOOK_INSTALLED = False
_COPY_EVENTS = frozenset({"shutil.copyfile", "shutil.copytree", "shutil.move", "shutil.copy", "shutil.copy2",
                          "shutil.copymode", "shutil.copystat", "_winapi.CopyFile2"})
_PROC_EVENTS = frozenset({"subprocess.Popen", "os.system", "os.exec", "os.spawn", "os.posix_spawn"})
_PATH_EVENTS = frozenset({
    "open", "os.listdir", "os.scandir", "glob.glob", "glob.glob/2", "os.startfile", "sqlite3.connect"}
    ) | _COPY_EVENTS | _PROC_EVENTS
_TEXT_EVENTS = _PROC_EVENTS
# position of the `env` argument in the audit event of every process-starting call (None: inherits os.environ)
_ENV_ARG = {"subprocess.Popen": 3, "os.exec": 2, "os.spawn": 3, "os.posix_spawn": 2}


def _flat_strings(obj, out, depth=0):
    if depth > 4 or obj is None:
        return
    if isinstance(obj, (str, bytes, bytearray, os.PathLike)):
        try:
            s = _to_str(obj)
        except Exception:
            return
        if s is not None:
            out.append(s)
    elif isinstance(obj, (list, tuple)):
        for x in obj:
            _flat_strings(x, out, depth + 1)


def _event_paths(event, args):
    """-> list of (path string, is_write)."""
    out = []
    if event == "open":
        s = _to_str(args[0]) if args else None
        if s is not None:
            mode = args[1] if len(args) > 1 else None
            flags = args[2] if len(args) > 2 else None
            out.append((s, _is_write_mode(mode) or (mode is None and _is_write_flags(flags))))
    elif event in ("os.listdir", "os.scandir"):
        a = args[0] if args else None
        if isinstance(a, int):
            return out
        out.append((_to_str(a) if a is not None else ".", False))
    elif event in ("glob.glob", "glob.glob/2"):
        pat = _to_str(args[0]) if args else None
        root_dir = _to_str(args[2]) if len(args) > 2 and args[2] is not None else None
        if pat is not None:
            out.append((pat, False))
            if root_dir is not None:
                out.append((os.path.join(root_dir, pat), False))
                out.append((root_dir, False))
    elif event == "subprocess.Popen":
        strs: list[str] = []
        _flat_strings(args[0] if args else None, strs)
        _flat_strings(args[1] if len(args) > 1 else None, strs)
        out.extend((s, False) for s in strs)
        cwd = _to_str(args[2]) if len(args) > 2 and args[2] is not None else None
        if cwd:
            out.append((cwd, False))
    elif event == "sqlite3.connect":
        targets, ro = _sqlite_targets(args[0] if args else None)
        out.extend((t, not ro) for t in targets)
    elif event in _COPY_EVENTS:
        # (src, dst[, flags]): src is read (removed for move), dst is written
        if len(args) >= 2:
            for i, a in enumerate(args[:2]):
                t = _to_str(a) if isinstance(a, (str, bytes, os.PathLike)) else None
                if t is not None:
                    out.append((t, i == 1 or event == "shutil.move"))
        else:
            strs = []
            _flat_strings(list(args), strs)
            out.extend((t, False) for t in strs)
    else:   # os.system / os.exec / spawn / startfile: every string argument
        strs = []
        _flat_strings(list(args), strs)
        out.extend((s, False) for s in strs)
    return out


# ---- child-process arguments: path candidates --------------------------------------------------------------
# On Windows `subprocess.Popen` reports the ALREADY JOINED command line (one string), so the arguments cannot be
# classified one by one; and a path may be written in any equivalent form (`..`, UNC admin share, `\?\UNC\`, a
# junction, a relative path from cwd, `--in=path`, a string literal inside `python -c`, doubled backslashes, %VAR%).
# Every string is therefore scanned for path-like fragments; each fragment is made absolute (against the cwd of the
# child) and classified by the IDENTITY of the nearest existing ancestor (the same machinery as for `open`).
_ABS_START = re.compile(r"(?=([A-Za-z]:[\\/]|\\\\|(?<![\w:/.\\])//))")
_CUT = re.compile(r"[\"'<>|*?\r\n\t`]")
_TOKEN_SPLIT = re.compile(r"[\s\"'<>|&,;=*?()\[\]{}`]+")
_PROC_MAX_CAND = 1500
_LOCAL_HOSTS = None


def _local_hosts() -> frozenset:
    global _LOCAL_HOSTS
    if _LOCAL_HOSTS is None:
        names = {"localhost", "127.0.0.1", "::1", "[::1]", ".", "?", os.environ.get("COMPUTERNAME", "").lower()}
        try:
            names.add(socket.gethostname().lower())
        except OSError:
            pass
        names.discard("")
        _LOCAL_HOSTS = frozenset(names)
    return _LOCAL_HOSTS


def _text_variants(s: str) -> list[str]:
    out = [s]
    for v in (s.replace("\\\\", "\\"), s.replace('\\"', '"'), os.path.expandvars(s)):
        if v not in out:
            out.append(v)
    if "%" in s or "$" in s:        # expanded form of the collapsed text too
        v = os.path.expandvars(s.replace("\\\\", "\\"))
        if v not in out:
            out.append(v)
    return out


def _fragments(s: str, out: dict) -> None:
    for v in _text_variants(s):
        out[v] = None
        for m in _ABS_START.finditer(v):
            tail = v[m.start():]
            pre = tail[:4] if tail[:4] in ("\\\\?\\", "\\\\.\\") else ""   # the `?` of `\\?\` is not a wildcard
            cut = pre + _CUT.split(tail[len(pre):], 1)[0]
            out[cut] = None
            first = tail.split(None, 1)[0] if tail.split(None, 1) else ""
            out[first] = None
            out[tail] = None
        for tok in _TOKEN_SPLIT.split(v):
            if not tok:
                continue
            out[tok] = None
            if ":" in tok:
                for piece in tok.split(":"):
                    if piece:
                        out[piece] = None


def _resolve_fragment(frag: str, base: str) -> str | None:
    """Absolute normalised path of a fragment, or None if it cannot name a local file."""
    t = frag.strip().strip("\"'")
    if not t or "\0" in t or len(t) > 4096:
        return None
    t = _strip_ns(t)
    if t.startswith("\\\\") or t.startswith("//"):
        host = t.replace("/", "\\")[2:].split("\\", 1)[0].lower()
        if host not in _local_hosts():
            return None                      # remote server / not a local path (never stat the network)
    elif "/" not in t and "\\" not in t and t != "..":
        return None                          # a bare word stays in the cwd, which is checked on its own
    try:
        return os.path.abspath(os.path.join(base, t))
    except (OSError, ValueError):
        return None


def _proc_zone_hits(strs: list[str], cwd) -> list[tuple[str, Zone]]:
    """-> [(fragment, zone)] for every protected path (RESTRICTED / HOLDOUT) a child-process argument may name."""
    c = _CFG
    base = os.getcwd()
    cwd_s = _to_str(cwd) if cwd is not None else None
    if cwd_s:
        base = os.path.abspath(os.path.join(base, cwd_s))
    frags: dict = {}
    for s in strs:
        _fragments(s, frags)
    resolved: dict = {}
    resolved[base] = base
    if cwd_s:
        resolved[os.getcwd()] = os.getcwd()
    overflow = len(frags) > _PROC_MAX_CAND
    for fr in list(frags)[:_PROC_MAX_CAND]:
        r = _resolve_fragment(fr, base)
        if r is not None and r not in resolved:
            resolved[r] = fr
    hits = []
    for r, fr in resolved.items():
        z = _zone_of(_matched_bases(r, c))
        if z in (Zone.RESTRICTED, Zone.HOLDOUT):
            hits.append((fr, z))
    if overflow:        # cannot prove the rest harmless: treat as touching the holdout zone (same decision)
        hits.append(("<too many path-like fragments>", Zone.HOLDOUT))
    return hits


def _uri_path(rest: str) -> str:
    """Path part of a sqlite `file:` URI (text after `file:`, already cut at ? and #), percent-decoded."""
    if rest.startswith("//"):                    # authority (empty or localhost) up to the next slash
        i = rest.find("/", 2)
        rest = rest[i:] if i >= 0 else ""
    path = _unquote(rest, errors="replace")
    if len(path) > 2 and path[0] == "/" and path[1].isalpha() and path[2] == ":":
        path = path[1:]
    return path


def _sqlite_targets(db):
    """sqlite3.connect(database) -> ([candidate paths], read_only).

    Candidates are every reading of the argument that SQLite or the OS could use: the raw string, the URI path cut at
    the first `?`/`#` (this is what SQLite really opens: an unescaped `#` in a name silently truncates it and creates
    a different file) and the fully decoded string. `mode=ro` / `immutable=1` in the query make the open read-only."""
    s = _to_str(db)
    if s is None or s in ("", ":memory:"):
        return [], True
    out = [s]
    ro = False
    if s[:5].lower() == "file:":
        rest = s[5:]
        cut = len(rest)
        for ch in "?#":
            i = rest.find(ch)
            if i >= 0:
                cut = min(cut, i)
        query = ""
        q = rest.find("?")
        if q >= 0:
            query = rest[q + 1:].split("#", 1)[0]
        ro = any(kv in ("mode=ro", "immutable=1") for kv in query.split("&"))
        for cand in (_uri_path(rest[:cut]), _uri_path(rest)):
            if cand and cand not in out:
                out.append(cand)
    return out, ro


def _vault_lex_forms(c: Config) -> list[str]:
    forms = []
    for d in (c.vault, c.original_zip_dir):
        for f in _forms(os.fspath(d)):
            forms.append(_lex(f))
    return forms


_DOTDOT_RE = re.compile(r"[^\\/:\s\"']+\\\.\.(?=\\|$)")


def _text(s: str) -> str:
    """Lower-cased, separator-normalised text for substring checks: prefixes stripped, doubled backslashes collapsed,
    parent-directory sequences (dir + backslash + two dots) removed."""
    t = _lex(s).replace("\\\\", "\\")
    for _ in range(64):
        n = _DOTDOT_RE.sub("", t, count=1)
        if n == t:
            break
        t = n
    return t


def _holdout_lex_forms(c: Config) -> list[str]:
    forms = []
    for d in (c.literature_holdout, *c.holdout_sources):
        for f in _forms(os.fspath(d)):
            forms.append(_lex(f))
    return forms


def _lexical_hit(args) -> bool:
    """Coarse text check used when the precise check itself failed (fail-closed fallback)."""
    text = _text(repr(args))
    c = _CFG
    if c is not None and any(f in text for f in _vault_lex_forms(c) + _holdout_lex_forms(c)):
        return True
    return "vault" in text or "original_zip" in text or "literature_holdout" in text


# ---- fast path for `open` outside the protected zones -------------------------------------------------------
# The precise classification walks the ancestors of the path (abspath + realpath + one stat per ancestor): ~800 us.
# Only paths that were classified (by that precise code) as NOT restricted are cached, keyed by the exact string. A
# cache hit is trusted only if (1) the entry belongs to the current config object and is younger than _FAST_TTL
# seconds, (2) the identity (st_dev, st_ino) of the parent directory is the same as when the entry was made (a
# directory or junction swapped in for any ancestor changes the identity of the directory the path resolves to) and
# (3) the final component is not a symlink / reparse point (lstat). RESTRICTED results are never cached; stack
# inspection happens only on the slow path (a decision is only needed for restricted paths anyway).
# Residual risk (accepted, see docs/ISOLATION.md): a directory MOVED into the vault, keeping its identity, within TTL.
_FAST: dict = {}
_FAST_TTL = 1.0
_FAST_MAX = 8192
_lstat = os.lstat
_isabs = os.path.isabs
_dirname = os.path.dirname
_S_ISLNK = _stat_mod.S_ISLNK
_WRITE_CH = frozenset("wax+")
_WRITE_BITS = os.O_WRONLY | os.O_RDWR | os.O_APPEND | os.O_CREAT | os.O_TRUNC


def _fast_fill(s: str, z: Zone):
    try:
        if not _isabs(s) or "\0" in s:
            return
        parent = _dirname(s)
        st = os.stat(parent)
        try:
            ls = _lstat(s)
        except OSError:
            ls = None
        if ls is not None and (_S_ISLNK(ls.st_mode) or getattr(ls, "st_reparse_tag", 0)):
            return
        if len(_FAST) >= _FAST_MAX:
            _FAST.clear()
        _FAST[s] = (z, _CFG, parent, st.st_dev, st.st_ino, _monotonic() + _FAST_TTL)
    except (OSError, ValueError):
        return


def _fast_open_ok(args) -> bool:
    """True if this `open` event is provably harmless (cached non-restricted zone, still valid, not a write into a
    read-only zone). Any doubt -> False -> precise slow path."""
    try:
        p = args[0]
        if type(p) is not str:
            return False
        ent = _FAST.get(p)
        if ent is None:
            return False
        z, cfg, parent, dev, ino, exp = ent
        if cfg is not _CFG or _monotonic() > exp:
            return False
        if z is not Zone.OTHER and z is not Zone.REGISTRY:       # read-only zones: look at the write flag
            mode = args[1] if len(args) > 1 else None
            flags = args[2] if len(args) > 2 else None
            if isinstance(mode, str):
                if not _WRITE_CH.isdisjoint(mode):
                    return False
            elif mode is None:
                if not isinstance(flags, int) or flags & _WRITE_BITS:
                    return False
            else:
                return False
        st = os.stat(parent)
        if st.st_dev != dev or st.st_ino != ino:
            return False
        try:
            ls = _lstat(p)
        except OSError:
            return True              # does not exist (yet): a regular file would be created in the verified parent
        return not (_S_ISLNK(ls.st_mode) or getattr(ls, "st_reparse_tag", 0))
    except Exception:
        return False


# ---- child processes ------------------------------------------------------------------------------------
_BOOT_ROOT = ""      # normalised root the process was started with (the one children inherit through os.environ)


def _same_root(v, c: Config) -> bool:
    for r in {os.path.normcase(os.path.abspath(os.fspath(c.root))), _BOOT_ROOT}:
        if r and _same_dir(v, r):
            return True
    return False


def _enforce_child_env(event: str, args):
    """Every child started by a hexw process gets HEXW_GUARD=1 and HEXW_ROOT=<current root>. A DIFFERENT HEXW_ROOT
    (in the explicit env or in os.environ) coming from a caller without a white-listed frame is refused: otherwise
    the child, importing hexw, would guard a decoy root. An explicit env dict is completed IN PLACE (the hook cannot
    replace the argument). A child that does not import hexw is not protected at all (see docs/ISOLATION.md)."""
    c = _CFG
    if c is None:
        return
    idx = _ENV_ARG.get(event)
    env = args[idx] if (idx is not None and len(args) > idx) else None
    if env is None:
        env = os.environ
    if not (hasattr(env, "keys") and hasattr(env, "__setitem__")):
        raise AccessError(f"{event}: env of type {type(env).__name__} cannot be checked; pass a dict")
    try:
        keys = [k for k in list(env.keys()) if os.fsdecode(k).upper() == "HEXW_ROOT"]
        vals = [os.fsdecode(env[k]) for k in keys]
        if any(not _same_root(v, c) for v in vals):
            denied, allowed = _scan_stack()
            if not allowed:
                exc = TheoryAccessError if denied else AccessError
                raise exc(f"{event}: a different HEXW_ROOT is passed to a child process by a caller without a "
                          f"white-listed frame (denied={denied}); children keep the current root")
        if len(keys) != 1:      # absent, or several spellings of the name
            keep = vals[0] if vals else os.fspath(c.root)
            for k in keys:
                del env[k]
            env["HEXW_ROOT"] = keep
        for k in [k for k in list(env.keys()) if os.fsdecode(k).upper() == "HEXW_GUARD"]:
            del env[k]
        env["HEXW_GUARD"] = "1"
    except AccessError:
        raise
    except Exception as ex:
        raise AccessError(f"{event}: cannot enforce HEXW_GUARD/HEXW_ROOT in the child environment "
                          f"({type(ex).__name__}: {ex})") from ex


def _install_multiprocessing_bootstrap():
    """multiprocessing 'spawn' starts `python -c "from multiprocessing.spawn import spawn_main; ..."`: the child would
    run user code (main module re-import, unpickled worker functions) before anything imports hexw. Wrap
    spawn.get_command_line (monkey patch, applied once at `import hexw`) so that the child first sets HEXW_ROOT /
    HEXW_GUARD and imports hexw, i.e. the audit hook is installed before any user code of the child runs."""
    try:
        import multiprocessing.spawn as sp
    except Exception:
        return
    if getattr(sp.get_command_line, "_hexw_wrapped", False):
        return
    orig = sp.get_command_line
    pkg_parent = os.path.dirname(_PKG_DIR)

    def get_command_line(**kwds):
        cmd = orig(**kwds)
        c = _CFG
        if c is None or "-c" not in cmd:          # frozen executables: not supported
            return cmd
        boot = (f"import os, sys; os.environ['HEXW_ROOT'] = {os.fspath(c.root)!r}; os.environ['HEXW_GUARD'] = '1'; "
                f"sys.path.append({pkg_parent!r}); import hexw; ")
        if '"' in boot or "\n" in boot:           # the stdlib joins the command line with naive double quotes
            raise AccessError("multiprocessing: cannot build a guarded child command line for this root path")
        i = cmd.index("-c")
        cmd[i + 1] = boot + cmd[i + 1]
        return cmd

    get_command_line._hexw_wrapped = True
    sp.get_command_line = get_command_line


def _hook_impl(event, args, state):
    c = _CFG
    if c is None:
        return
    if event in _ENV_ARG or event == "os.system":
        _enforce_child_env(event, args)
    if event in _TEXT_EVENTS:
        strs = []
        _flat_strings(list(args) if event != "subprocess.Popen" else [args[0] if args else None,
                                                                      args[1] if len(args) > 1 else None], strs)
        cwd = args[2] if (event == "subprocess.Popen" and len(args) > 2) else None
        for fr, z in _proc_zone_hits(strs, cwd):
            state["restricted"] = True
            if z is Zone.RESTRICTED:
                _vault_decision(fr, from_hook=True)
            else:
                _holdout_decision(fr, from_hook=True)
    for s, is_write in _event_paths(event, args):
        # command lines / arguments that embed the vault path as text
        if event in _TEXT_EVENTS and len(s) > 8 and any(f in _text(s) for f in _vault_lex_forms(c)):
            state["restricted"] = True
            _vault_decision(s, from_hook=True)
            continue
        if event in _TEXT_EVENTS and len(s) > 8 and any(f in _text(s) for f in _holdout_lex_forms(c)):
            state["restricted"] = True
            _holdout_decision(s, from_hook=True)
            continue
        z = _zone_of(_matched_bases(s, c))
        if _ro_violation(z, is_write):
            raise ReadOnlyZoneError(f"write to read-only zone blocked by audit hook: {s}")
        if z is Zone.RESTRICTED:
            state["restricted"] = True
            _vault_decision(s, from_hook=True)
        elif z is Zone.HOLDOUT:
            state["restricted"] = True
            _holdout_decision(s, from_hook=True)
        elif event == "open":
            _fast_fill(s, z)


def _hook(event, args):
    if event == "open":
        if _fast_open_ok(args):
            return
    elif event not in _PATH_EVENTS:
        return
    state = {"restricted": False}
    try:
        _hook_impl(event, args, state)
    except AccessError:
        raise
    except Exception as e:  # noqa - FAIL CLOSED
        try:
            hit = state["restricted"] or _lexical_hit(args)
        except Exception:
            hit = True
        if hit:
            raise AccessError(f"audit hook failed ({type(e).__name__}: {e}); "
                              "access to vault/original_zip denied (fail-closed)") from e


def install_audit_hook():
    global _HOOK_INSTALLED
    if not _HOOK_INSTALLED:
        sys.addaudithook(_hook)
        _HOOK_INSTALLED = True


install_audit_hook()


def _init_config():
    global _CFG, _BOOT_ROOT
    _CFG = _make_config(_default_root())   # raises RootNotFoundError: the guard is never silently disabled
    _BOOT_ROOT = os.path.normcase(os.path.abspath(os.fspath(_CFG.root)))
    # children inherit these (the audit hook completes explicit env dicts too); see _enforce_child_env
    os.environ["HEXW_ROOT"] = os.fspath(_CFG.root)
    os.environ["HEXW_GUARD"] = "1"
    _install_multiprocessing_bootstrap()


_init_config()
