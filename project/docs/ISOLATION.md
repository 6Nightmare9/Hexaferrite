# Data isolation: what is protected and what is not

Goal: stop **accidental** leakage of laboratory data (`vault/`, `raw/original_zip/`) into theory-side code
(`hexw.theory`, `params`, `dft`, `structures`, `spin`, any new subpackage). This is protection against our own
mistakes, **not a sandbox against hostile code**. Implementation: `hexw/io/paths.py` (white list of frames +
process-wide audit hook installed on `import hexw`). Whitelisted: `hexw.experiment`, `validate`, `audit`, `ingest`, `io`
(plus the explicit test-only `test_access()` context).

## Protected (blocked, tested on a synthetic vault)

| Channel | How |
|---|---|
| `open`, `os.open`, `pathlib`, `numpy`/`pandas` readers, `zipfile`, `io.*`, threads / pools started from theory | `open` audit event + whole-stack inspection |
| `os.listdir`, `os.scandir`, `glob`, `Path.iterdir`, `shutil.copytree` source | audit events |
| `shutil.copy`, `copy2`, `copyfile`, `copytree`, `move`, `copymode`, `copystat` | `shutil.*` events. **`shutil.copy2` on Python 3.13/Windows goes through `_winapi.CopyFile2` and emits no `open` event**; the hook handles the `_winapi.CopyFile2` event (found with `sys.addaudithook`). No monkey patch is needed. Destination of copy/move counts as a write. |
| `sqlite3.connect` (plain path and `file:` URI) | `sqlite3.connect` event. URI is parsed: authority dropped, cut at the first `?`/`#` (SQLite really opens the truncated name and *creates* it), percent-decoded; the raw string, cut and uncut decoded forms are all checked. Connections that are not `mode=ro`/`immutable=1` count as writes (read-only zones refuse them). |
| junctions, symlinks, 8.3 names, UNC / extended prefixes, case | classification by directory identity (`os.stat` / `samestat` of every ancestor) |
| writes into `raw/` and `vault/` | `ReadOnlyZoneError` for everybody |
| subprocess / `os.system` / `os.exec*` / `os.spawn*` with a protected path (vault, original_zip, HOLDOUT) in the arguments, the command line, `-c` code or the cwd, in ANY spelling | **path fragments are extracted and classified by identity**, see "Child-process arguments" below |
| child process env | every child of a hexw process gets `HEXW_GUARD=1` and `HEXW_ROOT=<current root>` (os.environ set at import; explicit env dicts are completed in place, any letter case of the name). A caller without a whitelisted frame that passes a **different** `HEXW_ROOT` (or poisons `os.environ`) is refused. A guarded child without `HEXW_ROOT` fails loudly (no auto-discovery). |
| `multiprocessing` (spawn; `Pool`, `ProcessPoolExecutor`) | **Monkey patch**: `multiprocessing.spawn.get_command_line` is wrapped once at `import hexw`; the child command line becomes `import os,sys; os.environ[HEXW_ROOT/HEXW_GUARD]=...; sys.path.append(...); import hexw; <stock spawn_main>`. So the hook is installed in the child before the main module is re-imported or any pickled worker runs. If the child's main module itself touches the vault at import, the child dies (and `Pool` keeps respawning it: that is stock `Pool` behaviour). Not supported: frozen executables. |

## Child-process arguments (stage 3 round 1)

On Windows `subprocess.Popen` reports the already **joined command line as ONE string** (`cmd /c type P`), so the
arguments cannot be checked one by one, and a purely lexical comparison with the root path is bypassed by every
equivalent spelling of the same file. The hook therefore treats each string argument (and the joined line, `-c` code,
PowerShell text) as text that may contain paths:

1. fragments are extracted: every substring that starts like `X:\` or `X:/`, like `\\` (UNC, `\\?\`, `\\.\`) or `//`,
   cut at a quote / redirection / wildcard character and also taken as a whole up to the first blank; every token of the
   text split at blanks, quotes, `=`, `,`, `;`, `:`, brackets (so `--in=path`, `/in:"path"`, `open('path')`,
   `Get-Content -LiteralPath 'path'` work); the variants with doubled backslashes collapsed (paths written inside
   Python / JSON string literals), escaped quotes removed, and `%VAR%` / `$VAR` expanded;
2. a fragment that can name a file (contains a separator, or is `..`) is made absolute against the **cwd of the child**
   (`cwd=` or the current directory) with `os.path.abspath` (so `..`, `.`, forward slashes, letter case, rooted paths
   without a drive and relative paths are normalised); the extended prefix and the administrative-share UNC forms
   (`\\localhost\D$\..`, `\\?\UNC\localhost\D$\..`) are stripped first; UNC paths of a remote host are skipped (the
   network is never touched);
3. the absolute path is classified by the **identity of the nearest existing ancestor** (`os.stat` + `samestat`, the
   same machinery as for `open`), so a path that does not exist yet is classified by its closest existing parent and a
   junction / symlink / 8.3 name leads to the same directory identity;
4. the cwd of the child and the cwd of the process are classified too (a bare file name in a protected cwd);
5. any fragment in RESTRICTED (vault, original_zip) or HOLDOUT goes through the normal decision (theory-side frame in
   the stack -> `TheoryAccessError`; no allowed frame -> `AccessError`; `hexw.audit` / `validate` / test context pass as
   for `open`). **Fail-closed**: more than 1500 path-like fragments in one call cannot be checked and count as a hit.

The same code serves both zones; tests: `tests/test_isolation_subprocess.py` (synthetic root, HOLDOUT and VAULT, every
spelling above x `cmd /c type`, `python -c open(..)`, PowerShell `Get-Content`, `--in=`, `os.system`, cwd, chdir,
junction).

**Known limit (not caught)**: a path that the child process computes by itself - read from a config / data file,
assembled from pieces (`"D:" + "\\Work" + ...`), taken from an environment variable that is not expanded in the
command line, or produced by a program that expands wildcards itself - is not visible to the parent-side check. Only a
child that imports `hexw` is guarded by its own hook.

## Declared derivatives and references (stage 3 round 1)

`holdout_sources.json` also lists the derivatives of holdout sources (`dft_feasibility/pilots/c_W92`, `make_bw.py`:
built from the L01 Zn HRPT CIF) - the files stay where they are, but the theory side cannot read them. **Rule for the
DFT campaign: the geometry of the W phase is taken from the fit source (L02 CIFs), never from L01 Co/Zn.** `hexw audit`
scans the text files of the project (outside `vault/`, `literature_holdout/`, `raw/`, the declared holdout sources and the third-party trees `tools/_dl`, `tools/pseudo`, `tools/qe_ucrt64`): a file that names a holdout source
(path or file name) must be declared holdout or be in `reference_allowlist` of the same JSON (the declaration, the guard
and split code, tests, this document, name-only inventories). `literature/fit` and `literature/excluded` are checked by
the stricter rules of `hexw.ingest.literature_split.find_fit_references`.

### Search results and pinned references (stage 4 round 1)

The intermediate and final results of the literature search contain annotations, abstracts or reasons per paper, including
the holdout papers, so they are holdout-derived and are declared in `holdout_sources`: `lit_search/results/` -
`rows_final.json`, `rows_stage2.json`, `merged.json`, `raw_hits*.json`, `coverage_from_abstracts.csv`,
`coverage_matrix.csv`, `review_list.txt`, `merged_list.txt`, `build_out.txt` and the full catalogue `candidates.csv`
(`reason`, `topic_tags` come from abstracts). Allowed instead: `candidates_meta.csv` (id, doi, title, authors, year, venue,
relevance, oa_status, downloaded_file, sha256; no abstract-derived text; the file paths and hashes of the declared-holdout downloads are blanked). It still contains the ids of those papers (they are bibliography), so it cannot be made free of a holdout file name and is therefore pinned by sha256 instead of being allow-listed freely. Theory code does not use `candidates.csv`
(`tests/test_no_holdout_hints.py`).

Files that name a declared source but are not holdout themselves (the search scripts `lit_search/scripts/*.py`, the
documents `GUIDE/*.md` and `lit_search/SEARCH_LOG.md` that cite the catalogue, `candidates_meta.csv`) are in
`reference_allowlist_sha256`: allowed **only while their bytes equal the pinned sha256** (`hexw audit` checks every pin and
reports a changed file as an offender). A deliberate edit of such a file therefore needs a new pin in
`project/config/holdout_sources.json`.

## NOT protected (documented limits)

* `ctypes` / native extensions calling `CreateFileW`, and `_winapi.CreateFile` (no audit event at all).
* Hard links created before/outside the process; `mmap` of an already open handle.
* **External programs and any child that does not import hexw** (pw.x, `cmd`, a Python child run with a script that
  never imports hexw): only the parent-side check of the Popen arguments protects them (a protected path in the command
  line / cwd, in any spelling, is blocked; a path computed or read from a file by the child is not). Accepted for the DFT adapter: its inputs are
  built by our code from theory-side data.
* `hexw.util.keepawake` (`SetThreadExecutionState`) and `hexw.util.proc` (`OpenProcess`/`GetExitCodeProcess`, is a PID alive) call kernel32 through `ctypes`. `hexw audit` has no rule against `ctypes` (its static policy covers subprocess / os.system / Popen / multiprocessing only), so no exception was needed. Neither function reads a file or starts a process; they take an integer (flags / PID) and return an integer / bool. They live in `hexw.util`, outside the scanned theory packages, and are imported by `hexw.dft.runner`. `keepawake` changes no power setting (`powercfg` is not touched); closing the lid still suspends the machine.
* Metadata: `os.stat`, `os.path.exists/getsize`, `os.path.getmtime` work (names, sizes, times are visible).
* `Path.rglob` and `os.walk` (default `onerror=None`) **silently yield nothing** for a denied directory because the
  library swallows the `OSError`; direct `os.listdir`/`scandir`/`iterdir` raise.
* A directory *moved* into the vault keeping its identity within the 1 s fast-path cache TTL (see below).
* Code that replaces/unpatches the hook or edits `paths._CFG` directly (the hook cannot be removed through the
  audit API, but the process is the attacker's own).

## Performance (fast path)

`open` of a path outside the restricted zones is classified once by the precise (slow, ~800 us) code and cached by
exact string. A hit is trusted only when: same config object, entry younger than 1 s, identity `(st_dev, st_ino)` of the
parent directory unchanged (a junction/directory swapped in for any ancestor changes it), and the last component is not
a symlink / reparse point. Restricted paths are never cached, and stack inspection is done only for them. Writes into
read-only zones always take the slow path. Measured (10 000 open/close of a file in `raw/theory_side`, best of 5):
34 us without hook, 847 us before, ~50 us now.

## Zone HOLDOUT (stage 3: literature holdout)

Separate zone of the same hook (`hexw.io.paths.Zone.HOLDOUT`): `literature_holdout/` and the source files of the
holdout publications declared in `project/config/holdout_sources.json` (they stay where they are in `raw/theory_side`
and `lit_search/downloads`; the whole L01 publication, all its CIF/PDF/XML and the working CIF copies are included,
because the L01 fit data (Mg, Ni) reach the theory only through `literature/fit`). Decision and reasons:
`literature/SPLIT_DECISION.md`.

* any `hexw.theory/params/dft/structures/spin/util` frame in the stack -> `TheoryAccessError` (open, listdir, glob,
  copy, sqlite, child processes), same as for the vault;
* access only for `hexw.validate` after the freeze file exists and for `hexw.audit` (purpose integrity/metadata);
  `hexw.ingest` only inside `paths.holdout_build()` (the split script); `experiment`/`io` alone are refused;
* the zone is read-only (`ReadOnlyZoneError`); files are read-only + `SHA256SUMS`;
* a declared holdout directory without the declaration file is a configuration error (no silent weakening);
* known limit: `pathlib.rglob()/glob()` swallow `PermissionError`, so a listing from a plain (non-hexw) frame may look
  empty instead of raising; opening a file raises.
Theory loads the literature via `hexw.params.literature.load_fit_records()` (only `literature/fit/records_fit.yaml`).

## Subprocess policy (stage 4A)

Code of `hexw.theory`, `hexw.params`, `hexw.structures`, `hexw.dft`, `hexw.spin` starts child processes **only** through:

* **Python children: `hexw.util.pool.WorkerPool`** (spawn context; the data guard is installed in the child before any user
  code, start-up handshake, timeouts). No `subprocess`, `os.system`, `Popen`, `multiprocessing.Process/Pool`,
  `ProcessPoolExecutor` or `os.exec*/spawn*` in these packages, in particular nothing with `sys.executable` or a `"python"`
  string.
* **Other programs (the future DFT adapter): `hexw.util.external.run_external`**, which accepts only executables of a white
  list (`ALLOWED_EXECUTABLES`: Quantum ESPRESSO `pw.x`, `pw.exe`, `pw_big.exe` - the Windows build with a larger stack, needed even for a 10-atom cell), always without a shell, with a list of arguments (a command
  line as one string, `-c`, a missing timeout are refused). An interpreter, a shell (`python`, `cmd`, `powershell`, `sh`,
  `bash`, `sys.executable`) or any name outside the list raises `ExternalProgramError`. Adding an executable is a
  deliberate edit of the list (reviewed). Paths in the arguments/cwd are still checked by the audit hook (a theory-side frame
  is on the stack): the vault and HOLDOUT stay closed.

**Enforcement.** `hexw audit` runs a static analysis (`hexw.audit.find_subprocess_violations`, `scan_theory_sources`) over every
`.py` file of those five packages: `ast` walk, calls of `subprocess.*`, `os.system/popen/exec*/spawn*/startfile`,
`multiprocessing.Process/Pool/get_context`, `concurrent.futures.ProcessPoolExecutor`, `asyncio.create_subprocess_*`, `pty`,
`pexpect`, `sh`, `plumbum` in any import spelling (`import subprocess as sp`, `from subprocess import Popen as P`,
`from os import system`, dotted `concurrent.futures`), plus `__import__('subprocess')`, `importlib.import_module('subprocess')`
and `getattr(os, 'system')`. A call that carries `sys.executable` or a `python` string is reported as such; a file that cannot
be parsed is a violation. Calls of `WorkerPool`/`run_external` are not reported. `hexw.util` (the wrappers), `ingest`, `io`,
`audit` and `registry` are outside the scanned set. Checks: `tests/test_subprocess_policy.py`.

**Static limits**: a launch hidden behind `eval/exec`, a string assembled at run time, `ctypes`/`cffi` calls into the C
runtime, or a third-party library that spawns processes internally is not seen; the runtime hook (above) still blocks
protected paths in any child started from a theory-side frame.
