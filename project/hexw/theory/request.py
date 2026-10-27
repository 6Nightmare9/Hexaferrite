"""Explicit SI input contract for a material forecast; g and damping have no hidden material defaults."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import yaml

from hexw.params.curie import normalize_composition


def composition_from_formula(formula):
    """Accept the explicitly supported SrFe(18-x)A(x)O27 family only."""
    import re
    formula=formula.translate(str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")).replace(" ", "")
    tokens=list(re.finditer(r"([A-Z][a-z]?)(\d*(?:\.\d+)?)",formula))
    if "".join(t.group(0) for t in tokens)!=formula or not tokens:
        raise ValueError("use an explicit formula, for example SrCo1Zn1Fe16O27; variables/parentheses are not supported")
    counts={}
    for t in tokens:
        element=t.group(1)
        counts[element]=counts.get(element,0.)+float(t.group(2) or 1.)
    if set(counts)-{"Sr","Fe","O","Co","Zn","Mn","Mg","Ni"}:
        raise ValueError("supported elements: Sr,Fe,O and Co,Zn,Mn,Mg,Ni")
    comp={k:v for k,v in counts.items() if k in {"Co","Zn","Mn","Mg","Ni"}}
    if not np.isclose(counts.get("Sr",0),1) or not np.isclose(counts.get("O",0),27) or not np.isclose(counts.get("Fe",0),18-sum(comp.values())):
        raise ValueError("supported stoichiometry: SrFe(18-x)A(x)O27, total x<=2")
    return normalize_composition(comp)


@dataclass
class ForecastRequest:
    composition: dict
    temperatures_K: list
    fields_T: list
    frequencies_GHz: list
    g: float
    alpha: float
    theta_deg: float = 0.0
    phi_deg: float = 0.0
    shape: str = "sphere"
    n_draws: int = 128
    n_prior: int = 60000
    seed: int = 20261026
    equilibrium_grid: int = 80
    host: str = "dft"
    co_form: str = "ion_mf"
    l11_constraint: str = "none"
    label: str = ""
    assumptions: dict = field(default_factory=dict)
    dielectric_model: dict | None = None
    porosity_fraction: float | None = None
    anisotropy_model: str = "auto"

    def __post_init__(self):
        self.composition = normalize_composition(self.composition)
        for name, low in [("temperatures_K", 0), ("frequencies_GHz", 0)]:
            a = np.asarray(getattr(self, name), float)
            if a.ndim != 1 or not len(a) or not np.isfinite(a).all() or np.any(a < low):
                raise ValueError(f"{name}: finite nonnegative one-dimensional array required")
            if name == "frequencies_GHz" and (len(a) < 5 or np.any(np.diff(a) <= 0) or a[0] <= 0):
                raise ValueError("frequencies_GHz: strictly increasing positive grid with >=5 points")
            setattr(self, name, a.tolist())
        f = np.asarray(self.fields_T, float)
        if f.ndim != 1 or not len(f) or not np.isfinite(f).all() or (len(f) > 1 and np.any(np.diff(f) <= 0)):
            raise ValueError("fields_T: finite strictly increasing array required")
        self.fields_T = f.tolist()
        for key in ("g", "alpha"):
            x = float(getattr(self, key))
            if not np.isfinite(x) or x <= 0:
                raise ValueError(f"{key} must be finite and >0 for finite-width spectra")
            setattr(self, key, x)
        if not (2 <= self.n_draws <= 10000) or not (1000 <= self.n_prior <= 2000000):
            raise ValueError("n_draws in [2,10000], n_prior in [1000,2000000] required")
        if self.equilibrium_grid < 24:
            raise ValueError("equilibrium_grid must be >=24")
        if self.shape not in {"sphere", "disk_perp_c", "disk_parallel_c", "no_shape"}:
            raise ValueError("unsupported shape")
        if not np.isfinite([self.theta_deg, self.phi_deg]).all() or not 0 <= self.theta_deg <= 180:
            raise ValueError("finite angles, theta in [0,180] required")
        if self.host not in {"dft", "l11"} or self.co_form not in {"ion_mf", "cc_power"}:
            raise ValueError("invalid physical-model branch")
        if self.l11_constraint not in {"none", "co_zn_ratio"}:
            raise ValueError("invalid l11_constraint")
        if self.anisotropy_model not in {"auto","legacy","sr_single_crystal"}:
            raise ValueError("anisotropy_model: auto, legacy or sr_single_crystal required")
        if self.porosity_fraction is not None and not 0<=self.porosity_fraction<.5:
            raise ValueError("porosity_fraction must be in [0,0.5); dilute spherical-pore approximation")
        if self.dielectric_model is not None:
            from hexw.theory.electromagnetic import dielectric
            model=dict(self.dielectric_model)
            if not model.pop("provenance",None):
                raise ValueError("dielectric_model requires explicit parameter provenance")
            dielectric(np.asarray(self.frequencies_GHz)*1e9,**model)
        # These values are declared scenarios until a material-specific provenance is supplied.
        self.assumptions = {"g": "user-specified scenario, not inferred from laboratory FMR",
                            "alpha": "user-specified Gilbert damping; not a prediction of total linewidth",
                            "shape": "explicit ideal specimen geometry; not a measurement-cell model",
                            **self.assumptions}

    @classmethod
    def load(cls, path):
        with open(Path(path), encoding="utf-8") as f:
            data = yaml.safe_load(f)
        if "formula" in data:
            composition=composition_from_formula(data.pop("formula"))
            if "composition" in data and data["composition"]!=composition:
                raise ValueError("formula and composition disagree")
            data["composition"]=composition
        for key in ("fields_T", "frequencies_GHz", "temperatures_K"):
            if isinstance(data.get(key), dict):
                spec = data[key]
                data[key] = np.linspace(float(spec["start"]), float(spec["stop"]), int(spec["count"])).tolist()
        return cls(**data)
