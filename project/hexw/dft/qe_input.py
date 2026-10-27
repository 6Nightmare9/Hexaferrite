"""pw.x inputs for collinear DFT+U calculations of a magnetic configuration, and their validation.

Species of different spin orientation get different labels (Fe1 = up, Fe2 = down; QE labels are at most 3 characters):
the starting_magnetization of every species carries the sign of the target configuration. Fe2+ ions (oxidation state 2 in the
ordered structure; their positions come from `hexw.structures.configs.generate_configurations`, ONE placeholder realisation)
get species of their own, Fe3 = Fe2+ spin up, Fe4 = Fe2+ spin down, so that the extra electron can be localised on them:
(i) their own starting_magnetization (config dft.fe2_localization.starting_magnetization_Fe2) and (ii) optionally
`starting_ns_eigenvalue(m, ispin, ityp)` of the &SYSTEM namelist (QE 7.6 accepts it with the HUBBARD card and the `atomic`
projector, "Modifying starting ns matrices according to input values"; with the Fe2+ spin up the 5th (largest) eigenvalue of the
spin-down occupation matrix is set to 1, i.e. d6 = d5(up) + one minority electron). It only seeds the first SCF iteration:
the run is checked afterwards (`hexw.dft.campaign.collect`, status `valence_changed`). tot_magnetization is NOT fixed
(a fixed value would bias the configuration); instead the moments after the run are compared with the target
(`hexw.dft.parse.check_configuration`). Occupations: gaussian smearing with a small width (insulators, but a flip can
change the occupations). DFT+U: HUBBARD card (QE >= 7.1) with the projector of the config (default ortho-atomic).
"""
from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Optional

import numpy as np

from hexw.dft import settings as ST
from hexw.dft.magnet import MagneticSystem


class InputError(ValueError):
    pass


def species_label(element: str, sign: float, reduced: bool = False) -> str:
    """Fe1/Fe2 = spin up/down; reduced (Fe2+ when localisation is seeded): Fe3/Fe4. At most 3 characters (QE)."""
    if sign == 0:
        return element
    if reduced:
        return element + ("3" if sign > 0 else "4")
    return element + ("1" if sign > 0 else "2")


def is_reduced_site(site) -> bool:
    """Fe2+ (the only mixed-valence ion of the model)."""
    try:
        return site.specie.symbol == "Fe" and int(round(site.specie.oxi_state)) == 2
    except Exception:                                            # noqa: BLE001 - no oxidation state
        return False


def fe2_localization(cfg: ST.DftConfig) -> dict:
    """The dft.fe2_localization block with defaults (a config without it behaves as before: mode 'none')."""
    d = dict(cfg.raw["dft"].get("fe2_localization") or {})
    d.setdefault("mode", "none")
    if d["mode"] not in ("none", "seeded"):
        raise InputError(f"dft.fe2_localization.mode must be 'none' or 'seeded', not {d['mode']!r}")
    ns = dict(d.get("starting_ns_eigenvalue") or {})
    ns.setdefault("use", False)
    d["starting_ns_eigenvalue"] = ns
    return d


def build_input(system: MagneticSystem, pattern, cfg: ST.DftConfig, *, prefix: str, config_id: str = "",
                profile: Optional[ST.PseudoProfile] = None, outdir: str = "./tmp") -> str:
    """Text of pw.in for the flip pattern `pattern` (+-1 per flip group; +1 everywhere = the ground state)."""
    st = system.structure
    d = cfg.raw["dft"]
    sysspec = cfg.system(system.kind)
    elements = sorted({s.specie.symbol for s in st})
    prof = profile or ST.load_profile(cfg, cfg.raw["pseudo_profile"], tuple(elements))
    for el in elements:
        if el not in prof.files:
            raise InputError(f"no pseudopotential for {el} in profile {prof.name}")
    sigma = system.sigma(np.asarray(pattern, float))
    loc = fe2_localization(cfg)
    seeded = loc["mode"] == "seeded"
    labels = []
    for k, s in enumerate(st):
        el = s.specie.symbol
        sg = 0.0 if system.spin[k] <= 0 else float(np.sign(sigma[k]))
        red = seeded and sg != 0 and is_reduced_site(s)
        labels.append((el, species_label(el, sg, red), sg))
    species = []                                         # (label, element, sign) in a stable order
    for el in elements:
        for lab in sorted({l for e, l, _ in labels if e == el}):
            sg = next(s for e, l, s in labels if l == lab)
            species.append((lab, el, sg))
    reduced_labels = {species_label("Fe", +1, True), species_label("Fe", -1, True)}
    gamma_only = tuple(int(x) for x in sysspec["kpoints"]) == (1, 1, 1) and sysspec.get("gamma_tricks", False)
    if gamma_only and d["hubbard"]["projector"] != "atomic":
        raise InputError("K_POINTS gamma (gamma_tricks) is implemented in QE 7.6 only for the 'atomic' Hubbard projector "
                         "(ortho-atomic stops with 'Gamma-only calculation for this case not implemented')")
    nat = len(st)
    conv = float(d["scf"]["conv_thr_per_atom_Ry"]) * nat
    occ = d["occupations"]
    lines = [
        f"! hexw.dft  campaign={cfg.campaign}  config={config_id}  config_sha256={cfg.sha256}",
        f"! system={system.name} kind={system.kind} pseudo_profile={prof.name} cutoff_status={prof.cutoff_status}",
        "&CONTROL",
        f"  calculation='scf', prefix='{prefix}', outdir='{outdir}', pseudo_dir='{prof.directory.as_posix()}',",
        f"  verbosity='{d['verbosity']}', tprnfor=.false., tstress=.false., disk_io='{d['disk_io']}'",
        "/",
        "&SYSTEM",
        f"  ibrav=0, nat={nat}, ntyp={len(species)}, ecutwfc={prof.ecutwfc:g}, ecutrho={prof.ecutrho:g},",
        f"  nspin=2, occupations='smearing', smearing='{occ['smearing']}', degauss={float(occ['degauss_Ry']):g},",
    ]
    if occ.get("tot_magnetization") is not None:
        lines.append(f"  tot_magnetization={float(occ['tot_magnetization']):g},")
    nsv = loc["starting_ns_eigenvalue"]
    ns_lines = []
    for i, (lab, el, sg) in enumerate(species, start=1):
        if lab in reduced_labels:
            m0 = float(loc["starting_magnetization_Fe2"]) * sg
            if nsv["use"]:                                                # the minority-spin channel gets one more electron
                ispin = 2 if sg > 0 else 1
                ns_lines.append(f"  starting_ns_eigenvalue({int(nsv['m_index'])},{ispin},{i})={float(nsv['value']):g},")
        else:
            m0 = float(d["starting_magnetization"].get(el, 0.0)) * sg if sg != 0 else 0.0
        lines.append(f"  starting_magnetization({i})={m0:g},")
    lines += ns_lines
    lines[-1] = lines[-1].rstrip(",")
    lines += [
        "/",
        "&ELECTRONS",
        f"  electron_maxstep={int(d['scf']['electron_maxstep'])}, conv_thr={conv:.3e}, mixing_beta={float(d['scf']['mixing_beta']):g},",
        f"  mixing_mode='{d['scf']['mixing_mode']}', diagonalization='{d['scf']['diagonalization']}'",
        "/",
        "ATOMIC_SPECIES",
    ]
    for lab, el, _ in species:
        lines.append(f"{lab} {float(d["atomic_masses"][el]):g} {prof.files[el]}")
    lines.append("CELL_PARAMETERS angstrom")
    for v in st.lattice.matrix:
        lines.append("  %.10f %.10f %.10f" % tuple(v))
    lines.append("ATOMIC_POSITIONS crystal")
    for (el, lab, _), s in zip(labels, st):
        f = s.frac_coords
        lines.append("%s %.10f %.10f %.10f" % (lab, *f))
    k = sysspec["kpoints"]
    if tuple(int(x) for x in k) == (1, 1, 1) and sysspec.get("gamma_tricks", False):
        lines.append("K_POINTS gamma")
    else:
        lines += ["K_POINTS automatic", "%d %d %d 0 0 0" % tuple(int(x) for x in k)]
    hub = d["hubbard"]
    u_lines = []
    for lab, el, _ in species:
        u = cfg.hubbard_U(el)
        if u is not None and u > 0:
            u_lines.append(f"U {lab}-3d {u:g}")
    if u_lines:
        lines.append("HUBBARD {%s}" % hub["projector"])
        lines += u_lines
    return "\n".join(lines) + "\n"


# ---- parsing + validation -------------------------------------------------------------------------------------------
_NL = re.compile(r"&(\w+)(.*?)^\s*/", re.S | re.M | re.I)
_KV = re.compile(r"(\w+(?:\([\d,\s]+\))?)\s*=\s*('[^']*'|\"[^\"]*\"|[^,\s/]+)")


def parse_pw_input(text: str) -> dict:
    body = "\n".join(l for l in text.splitlines() if not l.strip().startswith("!"))
    out = {"namelists": {}, "species": [], "positions": [], "cell": [], "kpoints": None, "hubbard": None, "U": []}
    for m in _NL.finditer(body + "\n"):
        name = m.group(1).lower()
        kv = {}
        for k, v in _KV.findall(m.group(2)):
            kv[k.lower()] = v.strip("'\"")
        out["namelists"][name] = kv
    rest = _NL.sub("", body)
    mode = None
    for raw in rest.splitlines():
        line = raw.strip()
        if not line:
            continue
        up = line.split()[0].upper()
        if up in ("ATOMIC_SPECIES", "CELL_PARAMETERS", "ATOMIC_POSITIONS", "K_POINTS", "HUBBARD"):
            mode = up
            if up == "K_POINTS":
                out["kpoints"] = {"style": line.split()[1].lower() if len(line.split()) > 1 else "tpiba", "mesh": None}
            if up == "HUBBARD":
                mm = re.search(r"\{(.*?)\}", line)
                out["hubbard"] = mm.group(1) if mm else None
            continue
        t = line.split()
        if mode == "ATOMIC_SPECIES":
            out["species"].append((t[0], float(t[1]), t[2]))
        elif mode == "CELL_PARAMETERS":
            out["cell"].append([float(x) for x in t[:3]])
        elif mode == "ATOMIC_POSITIONS":
            out["positions"].append((t[0], [float(x) for x in t[1:4]]))
        elif mode == "K_POINTS":
            if out["kpoints"]["style"] == "automatic":
                out["kpoints"]["mesh"] = [int(x) for x in t[:6]]
        elif mode == "HUBBARD":
            out["U"].append((t[0], t[1], float(t[2])))
    return out


def validate_input(text: str, cfg: ST.DftConfig, system: Optional[MagneticSystem] = None, pattern=None,
                   profile: Optional[ST.PseudoProfile] = None) -> list[str]:
    """List of problems (empty = valid)."""
    bad: list[str] = []
    try:
        P = parse_pw_input(text)
    except Exception as e:                                       # noqa: BLE001
        return [f"cannot parse the input: {e!r}"]
    nl = P["namelists"]
    for need in ("control", "system", "electrons"):
        if need not in nl:
            bad.append(f"namelist &{need.upper()} missing")
    if bad:
        return bad
    c, s, e = nl["control"], nl["system"], nl["electrons"]
    if c.get("calculation") != "scf":
        bad.append("calculation must be 'scf' (no relaxation on the laptop)")
    if s.get("nspin") != "2":
        bad.append("nspin must be 2")
    if "noncolin" in s or "lspinorb" in s:
        bad.append("collinear input must not contain noncolin / lspinorb")
    if s.get("ibrav") != "0":
        bad.append("ibrav must be 0")
    try:
        nat, ntyp = int(s["nat"]), int(s["ntyp"])
    except (KeyError, ValueError):
        return bad + ["nat / ntyp missing"]
    if nat != len(P["positions"]):
        bad.append(f"nat={nat} but {len(P['positions'])} positions")
    if ntyp != len(P["species"]):
        bad.append(f"ntyp={ntyp} but {len(P['species'])} species")
    labels = [x[0] for x in P["species"]]
    if len(set(labels)) != len(labels):
        bad.append("duplicate species labels")
    if any(len(l) > 3 for l in labels):
        bad.append("species label longer than 3 characters")
    used = {p[0] for p in P["positions"]}
    if used - set(labels):
        bad.append(f"positions use undeclared species {sorted(used - set(labels))}")
    if set(labels) - used:
        bad.append(f"declared species without atoms {sorted(set(labels) - used)}")
    prof = profile or ST.load_profile(cfg, cfg.raw["pseudo_profile"], tuple(sorted({re.match(r"[A-Z][a-z]?", l).group(0) for l in labels})))
    pdir = Path(c.get("pseudo_dir", ""))
    zval = {}
    for lab, mass, fn in P["species"]:
        f = pdir / fn
        if not f.is_file():
            bad.append(f"pseudopotential {f} not found")
            continue
        el = re.match(r"[A-Z][a-z]?", lab).group(0)
        if not fn.lower().startswith(el.lower()):
            bad.append(f"pseudopotential {fn} does not belong to element {el} (species {lab})")
        if fn != prof.files.get(el):
            bad.append(f"pseudopotential {fn} for {el} differs from the profile ({prof.files.get(el)})")
        try:
            zval[lab] = ST.upf_valence(f)
        except ST.DftConfigError as ex:
            bad.append(str(ex))
    try:
        ew, er = float(s["ecutwfc"]), float(s["ecutrho"])
        if ew + 1e-9 < prof.ecutwfc or er + 1e-9 < prof.ecutrho:
            bad.append(f"cutoffs {ew:g}/{er:g} Ry below the profile {prof.ecutwfc:g}/{prof.ecutrho:g}")
        if er < 4 * ew - 1e-9:
            bad.append("ecutrho < 4 ecutwfc")
    except (KeyError, ValueError):
        bad.append("ecutwfc/ecutrho missing")
    # starting magnetizations
    sm = {}
    for i, lab in enumerate(labels, start=1):
        key = f"starting_magnetization({i})"
        if key not in s:
            bad.append(f"{key} missing")
            continue
        v = float(s[key])
        if not (-1.0 <= v <= 1.0) or not math.isfinite(v):
            bad.append(f"{key} = {v} outside [-1, 1]")
        sm[lab] = v
    for key, val in s.items():
        mm = re.fullmatch(r"starting_ns_eigenvalue\((\d+),(\d+),(\d+)\)", key.replace(" ", ""))
        if not mm:
            continue
        m_, isp, ityp = (int(x) for x in mm.groups())
        try:
            v_ok = 0.0 <= float(val) <= 1.0
        except ValueError:
            v_ok = False
        if not (1 <= m_ <= 5 and isp in (1, 2) and 1 <= ityp <= len(labels) and v_ok):
            bad.append(f"{key} = {val}: bad indices or value (m 1..5, spin 1..2, species 1..{len(labels)}, value 0..1)")
        elif not (cfg.hubbard_U(re.match(r"[A-Z][a-z]?", labels[ityp - 1]).group(0)) or 0) > 0:
            bad.append(f"{key}: species {labels[ityp - 1]} has no Hubbard U")
    if "tot_magnetization" in s and cfg.raw["dft"]["occupations"].get("tot_magnetization") is None:
        bad.append("tot_magnetization is fixed in the input but not in the config")
    try:
        if not float(e["conv_thr"]) > 0:
            bad.append("conv_thr must be positive")
        if int(e.get("electron_maxstep", 100)) < 20:
            bad.append("electron_maxstep < 20")
    except (KeyError, ValueError):
        bad.append("conv_thr missing")
    # geometry
    cell = np.array(P["cell"])
    if cell.shape != (3, 3) or abs(np.linalg.det(cell)) < 1.0:
        bad.append("bad cell")
    else:
        fr = np.array([p[1] for p in P["positions"]])
        d = fr[:, None, :] - fr[None, :, :]
        d -= np.round(d)
        dist = np.linalg.norm(d @ cell, axis=2) + np.eye(len(fr)) * 99
        if dist.min() < 1.0:
            bad.append(f"atoms closer than 1 A ({dist.min():.2f} A)")
    # electrons
    if zval and len(zval) == len(labels):
        nelec = sum(zval[p[0]] for p in P["positions"])
        if nelec <= 0 or abs(nelec - round(nelec)) > 1e-6:
            bad.append(f"non-integer or non-positive electron count {nelec}")
    # k points
    kp = P["kpoints"]
    if kp is None:
        bad.append("K_POINTS missing")
    elif kp["style"] == "automatic" and (not kp["mesh"] or len(kp["mesh"]) != 6 or min(kp["mesh"][:3]) < 1):
        bad.append("bad automatic k mesh")
    # Hubbard
    el_of = {lab: re.match(r"[A-Z][a-z]?", lab).group(0) for lab in labels}
    if P["hubbard"] is not None and P["hubbard"] != cfg.raw["dft"]["hubbard"]["projector"]:
        bad.append(f"Hubbard projector {P['hubbard']} differs from the config")
    have = {}
    for kind, lab, val in P["U"]:
        name = lab.rsplit("-", 1)[0]
        if kind != "U" or name not in labels:
            bad.append(f"Hubbard line for unknown species {lab}")
        have[name] = val
    for lab in labels:
        want = cfg.hubbard_U(el_of[lab])
        if (want or 0) > 0 and abs(have.get(lab, -1) - want) > 1e-9:
            bad.append(f"U for {lab} is {have.get(lab)} but the config says {want}")
        if not want and lab in have:
            bad.append(f"U given for {lab}, which has none in the config")
    # target pattern
    if system is not None and pattern is not None:
        sigma = system.sigma(np.asarray(pattern, float))
        if len(sigma) != len(P["positions"]):
            bad.append("atom number differs from the system")
        else:
            for k, (lab, _) in enumerate(P["positions"]):
                if system.spin[k] > 0:
                    if sm.get(lab, 0) * sigma[k] <= 0:
                        bad.append(f"atom {k + 1} ({lab}): starting moment sign disagrees with the target configuration")
                        break
                elif abs(sm.get(lab, 0.0)) > 0:
                    bad.append(f"non-magnetic atom {k + 1} ({lab}) has a starting moment")
                    break
    return bad


def assert_valid(text: str, cfg: ST.DftConfig, system=None, pattern=None, profile=None) -> None:
    bad = validate_input(text, cfg, system, pattern, profile)
    if bad:
        raise InputError("; ".join(bad))
