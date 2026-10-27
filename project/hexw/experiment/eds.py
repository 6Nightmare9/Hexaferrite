"""Quantified EDS export -> atomic fractions; not raw-spectrum ZAF/standards correction."""
import numpy as np
from hexw.workflow import summary


def analyze(rows,*,normalize_element="Sr",seed=20261026):
    from pymatgen.core import Element
    symbols=[r["element"].strip() for r in rows]
    if not symbols or len(set(symbols))!=len(symbols):
        raise ValueError("unique chemical elements required in EDS table")
    mass=np.array([float(Element(s).atomic_mass) for s in symbols])
    weight=np.array([float(r["weight_percent"]) for r in rows])
    if not np.isfinite(weight).all() or np.any(weight<0) or weight.sum()<=0:
        raise ValueError("finite nonnegative quantified EDS weight percentages required")
    atoms=weight/mass; atomic=100*atoms/atoms.sum()
    if normalize_element not in symbols or atoms[symbols.index(normalize_element)]<=0:
        raise ValueError("normalization element missing/nonpositive; choose an observed element")
    ref=symbols.index(normalize_element)
    ratios=atoms/atoms[ref]
    intervals=None
    if all(r.get("std_weight_percent") not in (None,"") for r in rows):
        sd=np.array([float(r["std_weight_percent"]) for r in rows])
        if not np.isfinite(sd).all() or np.any(sd<0):
            raise ValueError("EDS standard uncertainties must be nonnegative")
        draws=np.maximum(np.random.default_rng(seed).normal(weight,sd,size=(10000,len(weight))),0)/mass
        good=draws[:,ref]>0
        draws=draws[good]/draws[good,ref,None]
        intervals={s:summary(draws[:,j]) for j,s in enumerate(symbols)}
    return {"status":"quantified_export_analysis","input_weight_sum_percent":float(weight.sum()),
        "atomic_percent":dict(zip(symbols,map(float,atomic))),"normalized_to":normalize_element+"=1",
        "atomic_ratios":dict(zip(symbols,map(float,ratios))),"ratio_scenario_intervals":intervals,
        "limits":["input must already contain instrument-quantified weight percentages; raw X-ray spectra not quantified here",
            "uncertainty propagation assumes independent reported Gaussian errors with nonnegative truncation",
            "oxygen and light-element quantification may be biased; oxygen stoichiometry is not imposed",
            "element ratios do not identify valence or crystallographic occupancy","no automatic update of frozen theory composition"]}
