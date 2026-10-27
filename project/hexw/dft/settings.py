"""DFT campaign settings: loading, validation and hashing of config/dft/campaign_v2.yaml.

The settings are fixed from the fit literature (record L05, Islam 2025) and from physics BEFORE any comparison with
laboratory data. `load_config` checks the literature numbers against `hexw.params.literature` (fit records only), so a
config edited away from the literature value is refused. The config hash (sha256 of the canonical JSON of the parsed
document) goes into every manifest / status / provenance record.
"""
from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml

from hexw.io import paths

CONFIG_NAME = "campaign_v2.yaml"          # v1 is kept unchanged for the record (smoke run of the chain); a new hash is a new version
MAGNETIC_ELEMENTS = ("Fe", "Co", "Mn")


class DftConfigError(ValueError):
    pass


def default_config_path() -> Path:
    return Path(__file__).resolve().parents[2] / "config" / "dft" / CONFIG_NAME


def canonical_json(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


@dataclass
class DftConfig:
    raw: dict
    sha256: str                      # canonical JSON hash
    path: Optional[Path] = None
    literature: dict = field(default_factory=dict)      # checked fit-record values (rec_id, value, method)

    # ---- shortcuts ---------------------------------------------------------------------------------------------
    @property
    def campaign(self) -> str:
        return self.raw["campaign"]

    @property
    def root(self) -> Path:
        return paths.get_config().root

    def runs_dir(self) -> Path:
        return self.root / self.raw["runs_root"] / self.campaign

    def bin_dir(self) -> Path:
        return self.root / self.raw["engine"]["bin_dir"]

    def executable(self) -> str:
        return self.raw["engine"]["executable"]

    def executable_path(self) -> Path:
        return self.bin_dir() / self.executable()

    def system(self, name: str) -> dict:
        try:
            return self.raw["systems"][name]
        except KeyError:
            raise DftConfigError(f"unknown system {name!r}; known: {sorted(self.raw['systems'])}") from None

    def profile(self, name: Optional[str] = None) -> "PseudoProfile":
        return load_profile(self, name or self.raw["pseudo_profile"])

    def exchange(self, system: str) -> dict:
        """Global 'exchange' block with the overrides of the given system."""
        out = dict(self.raw["exchange"])
        out.update(self.system(system).get("exchange", {}))
        return out

    def hubbard_U(self, element: str) -> Optional[float]:
        u = self.raw["dft"]["hubbard"]["U_eff_eV"].get(element)
        return None if u is None else float(u)

    def with_overrides(self, **kw) -> "DftConfig":
        """A copy with top-level keys replaced (hash recomputed); used by tests."""
        raw = json.loads(json.dumps(self.raw))
        raw.update(kw)
        _validate(raw)
        return DftConfig(raw, sha256_bytes(canonical_json(raw).encode()), self.path)


@dataclass
class PseudoProfile:
    name: str
    directory: Path
    files: dict                      # element -> file name
    ecutwfc: float
    ecutrho: float
    cutoff_status: str
    note: str = ""
    cutoff_basis: str = ""


def _need(d: dict, *keys):
    cur = d
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            raise DftConfigError(f"config key {'.'.join(keys)} missing")
        cur = cur[k]
    return cur


def _validate(raw: dict) -> None:
    for k in ("schema", "campaign", "runs_root", "engine", "dft", "pseudo_profile", "pseudo_profiles", "systems",
              "exchange", "runtime", "moments_check", "soc", "literature_basis"):
        _need(raw, k)
    if raw["schema"] != 1:
        raise DftConfigError(f"unsupported schema {raw['schema']!r}")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", str(raw["campaign"])):
        raise DftConfigError("bad campaign name")
    eng = raw["engine"]
    if int(eng["omp_threads"]) != 1:
        raise DftConfigError("omp_threads must be 1 (OpenMP slows this QE build down, FEASIBILITY section 1)")
    if int(eng["max_parallel"]) < 1:
        raise DftConfigError("max_parallel must be >= 1")
    d = raw["dft"]
    if int(d["nspin"]) != 2:
        raise DftConfigError("only collinear nspin=2 is supported by the J campaign")
    hub = d["hubbard"]
    if hub["projector"] not in ("ortho-atomic", "atomic"):
        raise DftConfigError("hubbard.projector must be 'ortho-atomic' or 'atomic'")
    for el, u in hub["U_eff_eV"].items():
        if u is not None and not (isinstance(u, (int, float)) and 0.0 <= u <= 12.0):
            raise DftConfigError(f"U_eff for {el} out of range: {u!r}")
        if el not in hub["U_eff_status"]:
            raise DftConfigError(f"U_eff_status missing for {el}")
    for el in ("Co", "Mn"):
        if hub["U_eff_eV"].get(el) is not None and hub["U_eff_status"][el] != "assumption":
            raise DftConfigError(f"U on {el} is not in the fit corpus: its status must be 'assumption'")
    if raw["pseudo_profile"] not in raw["pseudo_profiles"]:
        raise DftConfigError("active pseudo_profile is not defined")
    for name, system in raw["systems"].items():
        for k in ("kpoints", "seconds_per_scf_iteration", "scf_iterations", "ram_GB"):
            _need(system, k)
        if len(system["kpoints"]) != 3 or any(int(n) < 1 for n in system["kpoints"]):
            raise DftConfigError(f"bad k mesh of system {name}")
        if system.get("gamma_tricks", False) and hub["projector"] != "atomic":
            raise DftConfigError(f"system {name}: gamma_tricks (K_POINTS gamma) works in QE 7.6 only with the 'atomic' projector")
    loc = d.get("fe2_localization")
    if loc is not None:
        if loc.get("mode", "none") not in ("none", "seeded"):
            raise DftConfigError("dft.fe2_localization.mode must be 'none' or 'seeded'")
        if loc.get("mode") == "seeded" and not (0.0 < float(loc.get("starting_magnetization_Fe2", 0)) <= 1.0):
            raise DftConfigError("dft.fe2_localization.starting_magnetization_Fe2 must be in (0, 1]")
        ns = loc.get("starting_ns_eigenvalue") or {}
        if ns.get("use") and not (1 <= int(ns.get("m_index", 0)) <= 5 and 0.0 <= float(ns.get("value", -1)) <= 1.0):
            raise DftConfigError("dft.fe2_localization.starting_ns_eigenvalue: m_index 1..5, value 0..1")
    vc = raw.get("valence_check")
    if vc is not None and any(float(vc[k]) <= 0 for k in vc):
        raise DftConfigError("valence_check thresholds must be positive")
    ex = raw["exchange"]
    if not (0 < float(ex["r_cut_A"]) < 8):
        raise DftConfigError("exchange.r_cut_A out of range")
    if ex["group_by"] not in ("sublattice_pair", "shell"):
        raise DftConfigError("exchange.group_by must be 'sublattice_pair' or 'shell'")
    if int(ex["n_holdout"]) < 1:
        raise DftConfigError("at least one hold-out configuration is mandatory (TZ 8.8)")
    if int(ex["n_extra_fit"]) < 0:
        raise DftConfigError("n_extra_fit must be >= 0")
    if float(ex.get("far_prior_sigma_meV", 0.1)) <= 0 or not (0.0 < float(ex.get("holdout_rel_tol", 0.02)) < 1.0):
        raise DftConfigError("exchange.far_prior_sigma_meV must be > 0 and holdout_rel_tol in (0, 1)")
    if raw["soc"].get("allow_run") not in (False, True):
        raise DftConfigError("soc.allow_run must be a boolean")


def _check_literature(raw: dict) -> dict:
    """Compare the literature numbers of the config with the fit records (hexw.params.literature)."""
    from hexw.params.literature import load_fit_records
    recs = {c.rec_id: c for c in load_fit_records()}
    out = {}
    for key, spec in raw["literature_basis"].items():
        if not isinstance(spec, dict) or "rec_id" not in spec:
            continue
        rid = spec["rec_id"]
        if rid not in recs:
            raise DftConfigError(f"literature_basis.{key}: record {rid} is not in the fit corpus")
        est = recs[rid].estimate
        if est is None or abs(float(est.value) - float(spec["value"])) > 1e-9 * max(1.0, abs(float(spec["value"]))):
            raise DftConfigError(f"literature_basis.{key}: config value {spec['value']} != record {rid} "
                                 f"value {None if est is None else est.value}")
        out[key] = {"rec_id": rid, "value": float(est.value), "unit": est.unit, "lit_id": recs[rid].lit_id,
                    "method": recs[rid].source.method}
    # the active U on Fe must be the L05 value
    ufe = raw["dft"]["hubbard"]["U_eff_eV"].get("Fe")
    if ufe is None or abs(float(ufe) - out["U_eff_Fe_eV"]["value"]) > 1e-9:
        raise DftConfigError(f"U_eff(Fe) = {ufe} differs from the literature value {out['U_eff_Fe_eV']['value']} (L05-001)")
    return out


def load_config(path: Optional[str | Path] = None, check_literature: bool = True) -> DftConfig:
    p = Path(path) if path else default_config_path()
    if not p.is_file():
        raise DftConfigError(f"config {p} not found")
    with open(p, encoding="utf-8") as f:
        raw = yaml.safe_load(f)
    if not isinstance(raw, dict):
        raise DftConfigError("config must be a mapping")
    _validate(raw)
    cfg = DftConfig(raw, sha256_bytes(canonical_json(raw).encode()), p)
    if check_literature:
        cfg.literature = _check_literature(raw)
    return cfg


# ---- pseudopotentials ----------------------------------------------------------------------------------------------
def upf_valence(p: Path) -> float:
    """Valence charge of a UPF file (v2 attribute z_valence or the v1 header line '... Z valence')."""
    head = p.read_text(errors="replace")[:20000]
    m = re.search(r'z_valence\s*=\s*"\s*([0-9.Ee+-]+)\s*"', head)
    if m:
        return float(m.group(1))
    m = re.search(r"^\s*([0-9.Ee+-]+)\s+Z valence", head, re.M)
    if m:
        return float(m.group(1))
    raise DftConfigError(f"cannot read the valence charge from {p.name}")


def sssp_recommendation(meta_path: Path, elements) -> tuple[float, float, dict]:
    with open(meta_path, encoding="utf-8") as f:
        meta = json.load(f)
    per = {}
    for el in elements:
        if el not in meta:
            raise DftConfigError(f"element {el} not in the SSSP metadata")
        per[el] = (float(meta[el]["cutoff_wfc"]), float(meta[el]["cutoff_rho"]), meta[el]["filename"])
    return max(v[0] for v in per.values()), max(v[1] for v in per.values()), per


def load_profile(cfg: DftConfig, name: str, elements=("Sr", "Fe", "O")) -> PseudoProfile:
    try:
        spec = cfg.raw["pseudo_profiles"][name]
    except KeyError:
        raise DftConfigError(f"unknown pseudopotential profile {name!r}") from None
    d = cfg.root / spec["dir"]
    basis = "config"
    if spec.get("cutoffs") == "from_sssp_metadata":
        ew, er, per = sssp_recommendation(cfg.root / spec["metadata"], elements)
        basis = "SSSP 1.3.0 PBE efficiency metadata: " + "; ".join(f"{k} {v[0]:g}/{v[1]:g}" for k, v in per.items())
    else:
        ew, er = float(spec["ecutwfc_Ry"]), float(spec["ecutrho_Ry"])
    return PseudoProfile(name, d, dict(spec["files"]), ew, er, spec["cutoff_status"], spec.get("note", ""), basis)
