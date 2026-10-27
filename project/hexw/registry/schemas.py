"""Registry entities (ТЗ section 5) as dataclasses with JSON/YAML (de)serialisation.
Units are mandatory wherever a number is stored."""
from __future__ import annotations

import json
import math
import re
from dataclasses import MISSING, dataclass, field, fields, is_dataclass
from enum import Enum
from typing import Any, Optional

import yaml

from hexw.units import UNITS, Kind, UnitError, check_unit

# Stage 3: fields added to existing dataclasses are omitted from the plain form (to_plain / dump_*) while they hold
# their default value, so the serialisation - and therefore parameter-set hashes / model versions - of objects built
# the old way does not change.
OMIT = {"omit_if_default": True}


class SourceRole(str, Enum):
    LITERATURE_FIT = "literature_fit"
    LITERATURE_HOLDOUT = "literature_holdout"
    DFT = "dft"
    LABORATORY_VALIDATION = "laboratory_validation"
    LEGACY_PREDICTION = "legacy_prediction"
    EXPLICIT_CALIBRATION = "explicit_calibration"
    # stage 3: records kept for traceability but NEVER usable as parameter values
    EXCLUDED_SECONDARY = "excluded_secondary"          # secondary citation of a value (primary source not in corpus)
    EXCLUDED_LOW_QUALITY = "excluded_low_quality"      # internally inconsistent / low-quality source
    METHODOLOGY_ONLY = "methodology_only"              # conventions and methods; numbers are not transferred
    # stage 4 round 1: a named model ASSUMPTION (not a literature or DFT number). Not leaky: it is the modeller's own
    # statement; every ParameterEstimate whose source_id names ASSUMPTION must carry this role (see ParameterEstimate).
    ASSUMPTION = "assumption"


# roles forbidden as inputs of a blind theoretical model version. literature_holdout: a ParameterEstimate with this
# role exists (it is the data of a holdout claim, `LiteratureClaim.estimate`), but it is a validation TARGET, never a
# parameter: make_model_version / check_prediction_inputs_roles refuse it (hexw.registry.versioning).
LEAKY_ROLES = frozenset({SourceRole.LABORATORY_VALIDATION, SourceRole.LEGACY_PREDICTION,
                         SourceRole.LITERATURE_HOLDOUT})
# roles that cannot be used in a parametrization at all (a ParameterEstimate cannot even be built with them)
NON_PARAMETRIC_ROLES = frozenset({SourceRole.EXCLUDED_SECONDARY, SourceRole.EXCLUDED_LOW_QUALITY,
                                  SourceRole.METHODOLOGY_ONLY})

PRELIMINARY = "preliminary_from_filename"


class Conclusion(str, Enum):
    PASS = "pass"
    FAIL = "fail"
    INCONCLUSIVE = "inconclusive"
    NOT_EVALUATED = "not_evaluated"


class ValidationMode(str, Enum):
    """How blind a comparison with laboratory data really is (see registry/prior_use.yaml)."""
    RETROSPECTIVE = "retrospective"                  # data were analysed/used before; not a blind test
    INDEPENDENT_CANDIDATE = "independent_candidate"  # no prior use found AND confirmed by the user


PRIOR_USE_STATUSES = ("used_prior", "unknown_prior_use", "no_prior_use_found")


class RunState(str, Enum):
    CREATED = "created"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    NOT_IMPLEMENTED = "not_implemented"


# A non-identifiable parameter (the data cannot determine it) must not carry a narrow uncertainty.
# "Narrow" = relative half-width < NONIDENTIFIABLE_MIN_REL_HALFWIDTH of |value|, where the half-width is
#   std            -> std (1 sigma)
#   interval       -> (upper - lower) / 2
#   distribution   -> sigma/std/scale (normal-like) or (high - low) / 2 / (upper - lower) / 2 (uniform-like)
# If the half-width cannot be expressed relative to the value (value == 0, or a distribution whose width is not
# recognisable) the uncertainty must be kind="unknown" or carry the explicit flag declared_wide=True. The flag
# never overrides a computable narrow width.
NONIDENTIFIABLE_MIN_REL_HALFWIDTH = 0.5


def _finite(x, what):
    if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x):
        raise ValueError(f"{what} must be a finite number, got {x!r}")


UNCERTAINTY_CONVENTIONS = ("parenthetical_se", "range", "sd", "none")
_INEQUALITY_ALIASES = {">=": "ge", "≥": "ge", "<=": "le", "≤": "le"}


@dataclass
class Uncertainty:
    kind: str = "unknown"          # std | interval | distribution | unknown
    std: Optional[float] = None
    lower: Optional[float] = None
    upper: Optional[float] = None
    distribution: Optional[str] = None   # e.g. normal, uniform
    params: dict = field(default_factory=dict)
    declared_wide: bool = False    # explicit declaration "this is a wide / undetermined interval"
    # stage 3 (omitted from the plain form while default):
    #   convention: how the SOURCE states the uncertainty. parenthetical_se = 5.8955(1) (standard uncertainty in the
    #   last digits, IUCr), sd = standard deviation, range = the value is a range lo-hi (kind interval, lower/upper
    #   are lo/hi, value = midpoint), none = not stated / not applicable.
    #   inequality: the reported value is only a bound ("ge": true value >= value, "le": true value <= value).
    convention: str = field(default="none", metadata=OMIT)
    inequality: Optional[str] = field(default=None, metadata=OMIT)

    def __post_init__(self):
        if self.kind not in ("std", "interval", "distribution", "unknown"):
            raise ValueError(f"bad uncertainty kind {self.kind}")
        if self.convention not in UNCERTAINTY_CONVENTIONS:
            raise ValueError(f"bad uncertainty convention {self.convention!r}; allowed {UNCERTAINTY_CONVENTIONS}")
        if self.inequality is not None:
            self.inequality = _INEQUALITY_ALIASES.get(self.inequality, self.inequality)
            if self.inequality not in ("ge", "le"):
                raise ValueError(f"bad inequality {self.inequality!r}; allowed 'ge' (>=) / 'le' (<=)")
        if self.convention in ("parenthetical_se", "sd") and self.kind != "std":
            raise ValueError(f"convention {self.convention!r} needs kind='std'")
        if self.convention == "range" and self.kind != "interval":
            raise ValueError("convention 'range' needs kind='interval' with lower=lo, upper=hi")
        for name in ("std", "lower", "upper"):
            v = getattr(self, name)
            if v is not None:
                _finite(v, f"uncertainty.{name}")
        for k, v in (self.params or {}).items():
            if isinstance(v, (int, float)) and not isinstance(v, bool) and not math.isfinite(v):
                raise ValueError(f"uncertainty.params[{k!r}] must be finite, got {v!r}")
        if self.kind == "std" and (self.std is None or not self.std >= 0):
            raise ValueError("std uncertainty needs std >= 0")
        if self.kind == "interval" and (self.lower is None or self.upper is None or not self.lower <= self.upper):
            raise ValueError("interval uncertainty needs lower <= upper")
        if self.kind == "distribution" and not self.distribution:
            raise ValueError("distribution uncertainty needs a distribution name")

    def half_width(self) -> Optional[float]:
        """Half-width of the uncertainty in the units of the value; None if it cannot be determined."""
        if self.kind == "std":
            return self.std
        if self.kind == "interval":
            return (self.upper - self.lower) / 2
        if self.kind == "distribution":
            p = self.params or {}
            for k in ("std", "sigma", "scale"):
                if isinstance(p.get(k), (int, float)) and not isinstance(p.get(k), bool):
                    return abs(p[k])
            for lo, hi in (("low", "high"), ("lower", "upper"), ("a", "b"), ("min", "max")):
                if lo in p and hi in p:
                    return (p[hi] - p[lo]) / 2
        return None


class ProvenanceStatus(str, Enum):
    """Verification status of the provenance of a value."""
    VERIFIED = "verified"                # checked against the primary source
    UNVERIFIED = "unverified"
    DERIVED_BY_US = "derived_by_us"      # computed by the project from sourced numbers
    ASSUMPTION = "assumption"            # no source: a modelling assumption

    @classmethod
    def _missing_(cls, value):
        # the stage-5 DFT adapter records its fit results as "computed": a value computed by the project
        if value == "computed":
            return cls.DERIVED_BY_US
        return None


@dataclass
class Provenance:
    source_id: str                      # LiteratureRecord id, DFT run id, file sha256, ...
    locator: Optional[str] = None       # page/table/section/file path
    note: Optional[str] = None
    status: ProvenanceStatus = ProvenanceStatus.UNVERIFIED     # verified | unverified | derived_by_us | assumption
    literature: Optional["LiteratureSource"] = field(default=None, metadata=OMIT)   # stage 3: full literature trace

    def __post_init__(self):
        self.status = ProvenanceStatus(self.status)


def _req_unit(obj):
    obj.quantity_kind = check_unit(obj.quantity_kind, obj.unit)


@dataclass
class ParameterEstimate:
    name: str
    value: float
    unit: str
    quantity_kind: Kind
    source_role: SourceRole
    provenance: Provenance
    uncertainty: Uncertainty = field(default_factory=Uncertainty)
    identifiable: bool = True

    def __post_init__(self):
        _req_unit(self)
        self.source_role = SourceRole(self.source_role)
        if self.source_role in NON_PARAMETRIC_ROLES:
            raise ValueError(f"parameter {self.name!r}: role {self.source_role.value} cannot be used in a "
                             "parametrization (excluded secondary citation / low quality / methodology only)")
        if not math.isfinite(self.value):
            raise ValueError("value must be finite")
        if not self.name:
            raise ValueError("name required")
        # stage 4 round 1: an assumption is declared as such, in both directions
        names_assumption = "assumption" in str(self.provenance.source_id).lower()
        if names_assumption and self.source_role is not SourceRole.ASSUMPTION:
            raise ValueError(f"parameter {self.name!r}: source_id {self.provenance.source_id!r} names an ASSUMPTION, so the "
                             f"role must be 'assumption' (got {self.source_role.value}); split mixed records into an "
                             "assumption estimate and a sourced estimate")
        if names_assumption:
            parts = [q.strip() for q in str(self.provenance.source_id).split(";") if q.strip()]
            if any("assumption" not in q.lower() for q in parts):
                raise ValueError(f"parameter {self.name!r}: mixed record {self.provenance.source_id!r} (ASSUMPTION together with "
                                 "sourced ids); make two estimates")
        if self.source_role is SourceRole.ASSUMPTION:
            if not names_assumption:
                raise ValueError(f"parameter {self.name!r}: role 'assumption' requires an 'ASSUMPTION' source_id")
            if self.provenance.status is not ProvenanceStatus.ASSUMPTION:
                raise ValueError(f"parameter {self.name!r}: role 'assumption' requires provenance.status 'assumption'")
        elif self.provenance.status is ProvenanceStatus.ASSUMPTION:
            raise ValueError(f"parameter {self.name!r}: provenance.status 'assumption' requires role 'assumption'")
        if not self.identifiable:
            self._check_nonidentifiable_uncertainty()

    def _check_nonidentifiable_uncertainty(self):
        u = self.uncertainty
        if u.kind == "unknown":
            return
        hw = u.half_width()
        if hw is None or self.value == 0:       # width not comparable with the value: needs the explicit flag
            if u.declared_wide:
                return
            raise ValueError("non-identifiable parameter: the uncertainty width cannot be compared with the value; "
                             "use kind='unknown' or declared_wide=True")
        rel = hw / abs(self.value)              # computable: the numbers decide, the flag cannot override them
        if not rel >= NONIDENTIFIABLE_MIN_REL_HALFWIDTH:
            raise ValueError(
                f"non-identifiable parameter {self.name!r} must not get a narrow uncertainty "
                f"(relative half-width {rel:.3g} < {NONIDENTIFIABLE_MIN_REL_HALFWIDTH})")


class DataStatus(str, Enum):
    """What kind of number a literature value is."""
    MEASURED = "measured"
    CALCULATED = "calculated"
    DIGITIZED = "digitized"
    ASSUMED_BY_AUTHORS = "assumed_by_authors"
    DERIVED_BY_US = "derived_by_us"


HOLDOUT_KINDS = ("same_group_reproducibility", "independent_group", "transfer_Ba_to_Sr")
TRANSFERS = ("Ba_to_Sr",)


@dataclass
class Locator:
    """Where in the source a value is: at least one of the fields must be given."""
    text: Optional[str] = None            # original free-text locator
    page: Optional[str] = None
    table: Optional[str] = None
    figure: Optional[str] = None
    supplementary: Optional[str] = None   # supplementary file/table/figure

    def __post_init__(self):
        if not any(getattr(self, f.name) for f in fields(self)):
            raise ValueError("locator needs at least one of text/page/table/figure/supplementary")


@dataclass
class MeasurementConditions:
    T: Optional[Any] = None               # temperature, e.g. {"K": 300, "note": "..."}
    field: Optional[Any] = None           # applied field
    orientation: Optional[Any] = None
    form: Optional[Any] = None            # powder / ceramic / single crystal / film / composite
    density: Optional[Any] = None
    texture: Optional[Any] = None


@dataclass
class LiteratureSource:
    """Full trace of a literature value (TZ sections 4 and 13)."""
    locator: Locator
    data_status: DataStatus
    tier: int                                         # 1 (direct, same system) .. 4 (methodology / indirect)
    doi: Optional[str] = None
    sample_id_in_source: Optional[str] = None
    formula_in_source: Optional[str] = None
    synthesis: Optional[str] = None
    measurement_conditions: MeasurementConditions = field(default_factory=MeasurementConditions)
    method: Optional[str] = None
    quantity_definition: Optional[str] = None
    verbatim_source: Optional[str] = None             # the exact string from the source
    transformation: Optional[str] = None              # what we did to get the stored number (null: nothing)
    measurement_group: Optional[str] = None           # table/figure/CIF of one measurement stay together
    figure_only: bool = False                         # value exists only in a figure
    axes: Optional[str] = None                        # axes description (mandatory if figure_only)
    split_cluster: Optional[str] = None               # unit of fit/holdout division (publication/series)
    transfer: Optional[str] = None                    # e.g. "Ba_to_Sr": value is from an analogue system
    holdout_kind: Optional[str] = None                # for holdout records only

    def __post_init__(self):
        if isinstance(self.locator, str):
            self.locator = Locator(text=self.locator)
        self.data_status = DataStatus(self.data_status)
        if isinstance(self.tier, bool) or not isinstance(self.tier, int) or not 1 <= self.tier <= 4:
            raise ValueError(f"tier must be an integer 1..4, got {self.tier!r}")
        if self.doi is not None and not re.match(r"^10\.\S+/\S+", self.doi):
            raise ValueError(f"bad DOI {self.doi!r}")
        if self.figure_only and not self.axes:
            raise ValueError("figure_only value needs `axes` (what the axes of the figure are)")
        if self.data_status is DataStatus.DIGITIZED and not self.figure_only:
            raise ValueError("a digitized value comes from a figure: set figure_only=True and axes")
        if self.transfer is not None:
            if self.transfer not in TRANSFERS:
                raise ValueError(f"bad transfer {self.transfer!r}; allowed {TRANSFERS}")
            if self.tier < 3:
                raise ValueError("transferred values (analogue system) must be tier 3 or 4")
        if self.holdout_kind is not None and self.holdout_kind not in HOLDOUT_KINDS:
            raise ValueError(f"bad holdout_kind {self.holdout_kind!r}; allowed {HOLDOUT_KINDS}")


_PENDING_KEYS = ("name", "value", "unit", "proposed_quantity_kind", "reason")
# role of a literature claim -> its split (literature/SPLIT_DECISION.md)
ROLE_TO_SPLIT = {
    SourceRole.LITERATURE_FIT: "fit",
    SourceRole.LITERATURE_HOLDOUT: "holdout",
    SourceRole.EXCLUDED_SECONDARY: "excluded",
    SourceRole.EXCLUDED_LOW_QUALITY: "excluded",
    SourceRole.METHODOLOGY_ONLY: "methodology_only",
}


@dataclass
class LiteratureClaim:
    """One extracted literature value: ParameterEstimate (when the schema has a Kind/unit for it), or a `pending`
    numeric value (no Kind yet), or a qualitative statement. Always with the full `source` trace."""
    rec_id: str
    lit_id: str
    name: str
    role: SourceRole
    source: LiteratureSource
    estimate: Optional[ParameterEstimate] = None
    pending: Optional[dict] = None
    qualitative: Optional[str] = None
    value_original: Optional[str] = None
    unit_original: Optional[str] = None
    uncertainty_note: Optional[str] = None
    notes: Optional[str] = None
    split: Optional[str] = None           # fit | holdout | excluded | methodology_only; must agree with `role`

    def __post_init__(self):
        self.role = SourceRole(self.role)
        if not self.rec_id or not self.lit_id or not self.name:
            raise ValueError("rec_id, lit_id and name are required")
        if self.role not in ROLE_TO_SPLIT:
            raise ValueError(f"{self.rec_id}: role {self.role.value} is not a literature role")
        want = ROLE_TO_SPLIT[self.role]
        if self.split is None:
            self.split = want
        elif self.split != want:
            raise ValueError(f"{self.rec_id}: split {self.split!r} contradicts role {self.role.value} ({want!r})")
        if self.role in NON_PARAMETRIC_ROLES and self.estimate is not None:
            raise ValueError(f"{self.rec_id}: role {self.role.value} cannot carry a ParameterEstimate")
        if self.estimate is not None:
            if self.estimate.source_role is not self.role:
                raise ValueError(f"{self.rec_id}: estimate role {self.estimate.source_role.value} != {self.role.value}")
            if self.pending is not None:
                raise ValueError(f"{self.rec_id}: estimate and pending are mutually exclusive")
            if self.estimate.provenance.literature is None:
                self.estimate.provenance.literature = self.source
        if self.pending is not None:
            miss = [k for k in _PENDING_KEYS if k not in self.pending]
            if miss:
                raise ValueError(f"{self.rec_id}: pending value lacks {miss}")
            v = self.pending["value"]
            if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
                raise ValueError(f"{self.rec_id}: pending value must be a finite number, got {v!r}")
        if self.estimate is None and self.pending is None and not (
                self.qualitative or self.notes or self.source.verbatim_source):
            raise ValueError(f"{self.rec_id}: a claim without a number needs a qualitative statement, notes or the "
                             "verbatim source text")


# ---- helpers for sources that print uncertainties in their own convention ----
_PAREN = re.compile(r"^\s*([+-]?\d*\.?\d+)\((\d+)\)\s*$")
_RANGE = re.compile(r"^\s*([+-]?\d*\.?\d+)\s*[-–—]\s*([+-]?\d*\.?\d+)\s*$")
_INEQ = re.compile(r"^\s*(>=|<=|≥|≤)\s*([+-]?\d*\.?\d+)\s*$")


def parse_reported(text: str) -> tuple[float, Uncertainty]:
    """'5.8955(1)' -> (5.8955, std 1e-4, parenthetical_se); '60-65' -> (62.5, interval [60, 65], range);
    '>= 5' / '<= 5' -> (5, unknown, inequality); plain number -> (x, unknown, none)."""
    m = _PAREN.match(text)
    if m:
        num, digits = m.group(1), m.group(2)
        dec = len(num.split(".")[1]) if "." in num else 0
        return float(num), Uncertainty(kind="std", std=int(digits) * 10.0 ** (-dec), convention="parenthetical_se")
    m = _INEQ.match(text)
    if m:
        return float(m.group(2)), Uncertainty(kind="unknown", inequality=m.group(1))
    m = _RANGE.match(text)
    if m:
        lo, hi = float(m.group(1)), float(m.group(2))
        if hi < lo:
            raise ValueError(f"range {text!r}: hi < lo")
        return (lo + hi) / 2, Uncertainty(kind="interval", lower=lo, upper=hi, convention="range")
    try:
        return float(text), Uncertainty()
    except ValueError:
        raise ValueError(f"cannot parse reported value {text!r}") from None


@dataclass
class DataRef:
    """Reference to an input or a target. Overlap is detected via id, sha256 or group."""
    id: str
    sha256: Optional[str] = None
    group: Optional[str] = None          # e.g. Observation id (so .mdi and .txt collide)
    source_role: Optional[SourceRole] = None

    def __post_init__(self):
        if self.source_role is not None:
            self.source_role = SourceRole(self.source_role)


@dataclass
class OutputValue:
    name: str
    value: float
    unit: str
    quantity_kind: Kind
    uncertainty: Uncertainty = field(default_factory=Uncertainty)

    def __post_init__(self):
        _req_unit(self)


@dataclass
class Prediction:
    prediction_id: str
    model_version: str
    inputs: list[DataRef]
    outputs: list[OutputValue] = field(default_factory=list)
    created_utc: Optional[str] = None
    notes: Optional[str] = None


@dataclass
class Metric:
    name: str
    value: Optional[float] = None
    unit: Optional[str] = None
    tolerance: Optional[float] = None
    tolerance_kind: Optional[str] = None   # relative | absolute

    def __post_init__(self):
        if self.value is not None and not self.unit:
            raise ValueError("metric value requires a unit (1 for dimensionless)")


@dataclass
class ValidationResult:
    result_id: str
    model_version: str
    prediction_id: str
    targets: list[DataRef]
    metric: Metric
    conclusion: Conclusion = Conclusion.NOT_EVALUATED
    created_utc: Optional[str] = None
    notes: Optional[str] = None
    # prior-use disclosure carried over from registry/prior_use.yaml; the conservative default is
    # "retrospective": a result is never presented as blind unless the registry says so.
    validation_mode: ValidationMode = ValidationMode.RETROSPECTIVE
    prior_use_status: str = "unknown_prior_use"

    def __post_init__(self):
        self.conclusion = Conclusion(self.conclusion)
        self.validation_mode = ValidationMode(self.validation_mode)
        if self.prior_use_status not in PRIOR_USE_STATUSES:
            raise ValueError(f"bad prior_use_status {self.prior_use_status!r}")
        if (self.validation_mode is ValidationMode.INDEPENDENT_CANDIDATE
                and self.prior_use_status != "no_prior_use_found"):
            raise ValueError("independent_candidate requires prior_use_status == no_prior_use_found")


@dataclass
class Claim:
    text: str
    value: Optional[float] = None
    unit: Optional[str] = None
    quantity_kind: Optional[Kind] = None
    locator: Optional[str] = None
    extraction_status: str = "unverified"

    def __post_init__(self):
        if self.value is not None:
            if not self.unit or self.quantity_kind is None:
                raise ValueError("numeric claim requires unit and quantity_kind")
            self.quantity_kind = check_unit(self.quantity_kind, self.unit)


@dataclass
class LiteratureRecord:
    lit_id: str
    citation: str
    doi: Optional[str] = None
    year: Optional[int] = None
    path: Optional[str] = None             # relative to raw/theory_side
    access: str = "unknown"                # full_text | partial | missing
    role: Optional[SourceRole] = None
    claims: list[Claim] = field(default_factory=list)
    license: Optional[str] = None


@dataclass
class SampleRecord:
    sample_id: str
    batch_id: str
    label: str
    nominal_formula: Optional[dict] = None
    measured_composition: Optional[dict] = None
    precursors: Optional[dict] = None
    synthesis: Optional[dict] = None
    post_treatment: Optional[dict] = None
    geometry: Optional[dict] = None
    mass: Optional[dict] = None
    metadata_status: dict = field(default_factory=dict)
    files: list[dict] = field(default_factory=list)
    observations: list[str] = field(default_factory=list)
    ambiguity_flags: list[str] = field(default_factory=list)
    notes: Optional[str] = None


@dataclass
class Observation:
    observation_id: str
    sample_id: Optional[str]
    method: str                          # XRD | Mossbauer | ...
    files: list[dict]                    # archive_path, zone, sha256, representation
    semantic_group: Optional[str] = None
    header_metadata: dict = field(default_factory=dict)
    calibration: Optional[dict] = None
    uncertainty: Optional[dict] = None
    flags: list[str] = field(default_factory=list)
    independent: bool = True
    duplicate_of: Optional[str] = None


# ---- unified transition-metal site IDs of the W phase (literature/site_notation.md; P63/mmc) ----
SITE_IDS = ("W1", "W2", "W3", "W4", "W5", "W6", "W7")
SITE_NOTATION: dict[str, dict] = {
    "W1": {"wyckoff": "4e", "block": "S", "env": "tetrahedron", "gorter_spin": "down", "atoms_per_fu": 2},
    "W2": {"wyckoff": "4f", "label": "4f_IV", "block": "S", "env": "tetrahedron", "gorter_spin": "down", "atoms_per_fu": 2},
    "W3": {"wyckoff": "6g", "block": "S-S", "env": "octahedron", "gorter_spin": "up", "atoms_per_fu": 3},
    "W4": {"wyckoff": "4f", "label": "4f_VI", "block": "S", "env": "octahedron", "gorter_spin": "up", "atoms_per_fu": 2},
    "W5": {"wyckoff": "12k", "block": "R-S", "env": "octahedron", "gorter_spin": "up", "atoms_per_fu": 6},
    "W6": {"wyckoff": "4f", "label": "4f_VI*", "block": "R", "env": "octahedron", "gorter_spin": "down", "atoms_per_fu": 2},
    "W7": {"wyckoff": "2d", "block": "R", "env": "trigonal bipyramid", "gorter_spin": "up", "atoms_per_fu": 1,
           "note": "half-occupied split 4f in the CIFs; one site of the pair is occupied"},
}
# author-specific labels -> unified ID. "4f" alone is ambiguous (W2, W4, W6) and is NOT mapped; Nakai's 4f2/4f3 are the
# reverse of Islam's (site_notation.md, warning 1): the scheme must be named.
SITE_ALIASES: dict[str, dict[str, str]] = {
    "unified": {w: w for w in SITE_IDS},
    "wyckoff_unique": {"4e": "W1", "6g": "W3", "12k": "W5", "2d": "W7"},
    "morch_cif": {"Fe1": "W1", "Fe2": "W2", "Fe3": "W3", "Fe4": "W4", "Fe5": "W5", "Fe6": "W6", "Fe7": "W7",
                  "1TS": "W1", "2TS": "W2", "3OSS": "W3", "4OS": "W4", "5ORS": "W5", "6OR": "W6", "7PR": "W7"},
    "islam2025": {"4e": "W1", "4f1": "W2", "6g": "W3", "4f2": "W4", "12k": "W5", "4f3": "W6", "2d": "W7"},
    "nakai": {"4e": "W1", "4f1": "W2", "6g": "W3", "4f3": "W4", "12k": "W5", "4f2": "W6", "2d": "W7"},
}


def site_id(label: str, scheme: str) -> str:
    """Unified W1..W7 ID of an author's site label. Raises KeyError for unknown/ambiguous labels."""
    try:
        return SITE_ALIASES[scheme][label]
    except KeyError:
        raise KeyError(f"label {label!r} is unknown or ambiguous in scheme {scheme!r} "
                       f"(schemes: {sorted(SITE_ALIASES)})") from None


class MomentDefinition(str, Enum):
    """What a quoted magnetic moment is the moment OF. Different definitions are never compared or mixed."""
    FE_ION = "fe_ion"                                     # moment of a Fe ion on that site
    SPECIES_ION = "species_ion"                           # moment of the named (non-Fe) ion, e.g. Co
    SITE_AVERAGE = "site_average_incl_diamagnetic"        # average over all occupants, Zn/Mg/vacancies count as 0


class MomentDefinitionError(ValueError):
    """Moments with different definitions were compared / subtracted."""


@dataclass
class SiteMoment:
    value: float
    definition: MomentDefinition
    unit: str = "muB"
    species: Optional[str] = None         # FE_ION -> "Fe"; SPECIES_ION -> the ion
    uncertainty: Optional[float] = None   # 1 sigma in `unit`

    def __post_init__(self):
        self.definition = MomentDefinition(self.definition)
        check_unit(Kind.MAGNETIC_MOMENT, self.unit)
        _finite(self.value, "moment value")
        if self.definition is MomentDefinition.FE_ION:
            if self.species not in (None, "Fe"):
                raise ValueError("fe_ion moment is the moment of Fe")
            self.species = "Fe"
        elif self.definition is MomentDefinition.SPECIES_ION:
            if not self.species or self.species == "Fe":
                raise ValueError("species_ion moment needs a non-Fe species (use fe_ion for Fe)")
        elif self.species is not None:
            raise ValueError("a site average has no single species")
        if self.uncertainty is not None:
            _finite(self.uncertainty, "moment uncertainty")

    def _muB(self) -> float:
        return self.value * UNITS[Kind.MAGNETIC_MOMENT][self.unit]

    def _check_comparable(self, other):
        if not isinstance(other, SiteMoment):
            raise MomentDefinitionError(f"cannot compare a site moment with {type(other).__name__}")
        if (self.definition, self.species) != (other.definition, other.species):
            raise MomentDefinitionError(
                f"moments with different definitions are not comparable: {self.definition.value}"
                f"({self.species}) vs {other.definition.value}({other.species}). Convert explicitly "
                "(to_site_average) and state the occupancy used.")

    def __eq__(self, other):
        """Equal only for the same definition (and species): comparing different definitions is an ERROR
        (MomentDefinitionError), not a quiet False. Other types: NotImplemented (-> False)."""
        if not isinstance(other, SiteMoment):
            return NotImplemented
        self._check_comparable(other)
        return (self._muB(), self._u_muB()) == (other._muB(), other._u_muB())

    def _u_muB(self):
        return None if self.uncertainty is None else self.uncertainty * UNITS[Kind.MAGNETIC_MOMENT][self.unit]

    def __sub__(self, other) -> float:
        """Difference in muB; refuses different definitions."""
        self._check_comparable(other)
        return self._muB() - other._muB()

    def __lt__(self, other): self._check_comparable(other); return self._muB() < other._muB()
    def __le__(self, other): self._check_comparable(other); return self._muB() <= other._muB()
    def __gt__(self, other): self._check_comparable(other); return self._muB() > other._muB()
    def __ge__(self, other): self._check_comparable(other); return self._muB() >= other._muB()

    def to_site_average(self, occupancy_of_species: float) -> "SiteMoment":
        """EXPLICIT conversion ion moment -> site average: moment * occupancy of that ion on the site."""
        if self.definition is MomentDefinition.SITE_AVERAGE:
            raise MomentDefinitionError("already a site average")
        if not 0.0 <= occupancy_of_species <= 1.0:
            raise ValueError("occupancy must be within [0, 1]")
        return SiteMoment(self.value * occupancy_of_species, MomentDefinition.SITE_AVERAGE, self.unit, None,
                          None if self.uncertainty is None else self.uncertainty * occupancy_of_species)


class PartialOccupancyError(ValueError):
    """A partially occupied (averaged) CIF structure was used where a DFT configuration is required."""


@dataclass
class SiteOccupancy:
    site: str
    species: dict[str, float]            # species label -> occupancy (fraction of the site)
    multiplicity: Optional[int] = None
    wyckoff: Optional[str] = None
    # stage 3 (omitted from the plain form while default)
    site_id: Optional[str] = field(default=None, metadata=OMIT)       # unified W1..W7
    space_group: Optional[str] = field(default=None, metadata=OMIT)   # group the `wyckoff` label refers to
    split_label: Optional[str] = field(default=None, metadata=OMIT)   # sub-site in a low-symmetry (magnetic) group,
    #                                                                   e.g. W3 -> "8d" / "4b" in Cmcm
    moments: list[SiteMoment] = field(default_factory=list, metadata=OMIT)

    def __post_init__(self):
        if self.site_id is not None and self.site_id not in SITE_IDS:
            raise ValueError(f"bad site_id {self.site_id!r}; allowed {SITE_IDS}")
        if self.split_label is not None and (self.site_id is None or not self.space_group):
            raise ValueError("a split sub-site needs site_id (parent W-site) and the space_group it belongs to")
        seen = set()
        for m in self.moments:
            key = (m.definition, m.species)
            if key in seen:
                raise ValueError(f"site {self.site}: duplicate moment of definition {key}")
            seen.add(key)

    def moment(self, definition: MomentDefinition, species: Optional[str] = None) -> Optional[SiteMoment]:
        for m in self.moments:
            if m.definition is MomentDefinition(definition) and (species is None or m.species == species):
                return m
        return None


@dataclass
class StructureCandidate:
    candidate_id: str
    formula: str
    sites: list[SiteOccupancy]
    oxidation_states: dict[str, float] = field(default_factory=dict)
    magnetic_hypothesis: Optional[str] = None
    source: Optional[str] = None
    status: str = "candidate"
    # True: this is a partially occupied (averaged) structure taken from a CIF - a description of the sample, NOT a
    # configuration that may be given to DFT (needs enumeration of integer-atom orderings, charge neutral).
    partial_occupancy_cif_not_dft_config: bool = field(default=False, metadata=OMIT)

    def require_dft_configuration(self) -> None:
        """Raise unless this candidate can be used as an ordered DFT configuration."""
        if self.partial_occupancy_cif_not_dft_config:
            raise PartialOccupancyError(f"{self.candidate_id}: flagged partial_occupancy_cif_not_dft_config")
        frac = [(s.site, sp, x) for s in self.sites for sp, x in s.species.items()
                if not (math.isclose(x, 0.0, abs_tol=1e-9) or math.isclose(x, 1.0, abs_tol=1e-9))]
        if frac:
            raise PartialOccupancyError(f"{self.candidate_id}: fractional occupancies {frac[:3]} - not an ordered "
                                        "DFT configuration")


@dataclass
class RunManifest:
    run_id: str
    command: str
    started_utc: str
    state: RunState = RunState.CREATED
    settings: dict = field(default_factory=dict)
    package_versions: dict = field(default_factory=dict)
    resources: dict = field(default_factory=dict)
    checksums: dict = field(default_factory=dict)     # registry file -> sha256 (from SHA256SUMS)
    git_head: dict = field(default_factory=dict)      # {"project": sha, "registry": sha}
    model_version: Optional[str] = None
    finished_utc: Optional[str] = None
    notes: Optional[str] = None

    def __post_init__(self):
        self.state = RunState(self.state)


# ---------------- serialisation ----------------
def _is_default(f, v) -> bool:
    if f.default is not MISSING:
        return v == f.default
    if f.default_factory is not MISSING:
        return v == f.default_factory()
    return False


def to_plain(obj: Any):
    if isinstance(obj, Enum):
        return obj.value
    if is_dataclass(obj) and not isinstance(obj, type):
        out = {}
        for f in fields(obj):
            v = getattr(obj, f.name)
            if f.metadata.get("omit_if_default") and _is_default(f, v):
                continue
            out[f.name] = to_plain(v)
        return out
    if isinstance(obj, dict):
        return {k: to_plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_plain(v) for v in obj]
    return obj


def dump_json(obj, **kw) -> str:
    return json.dumps(to_plain(obj), ensure_ascii=False, sort_keys=True, **kw)


def dump_yaml(obj) -> str:
    return yaml.safe_dump(to_plain(obj), allow_unicode=True, sort_keys=False)


_NESTED = {
    (ParameterEstimate, "provenance"): Provenance,
    (ParameterEstimate, "uncertainty"): Uncertainty,
    (OutputValue, "uncertainty"): Uncertainty,
    (ValidationResult, "metric"): Metric,
    (Provenance, "literature"): LiteratureSource,
    (LiteratureSource, "locator"): Locator,
    (LiteratureSource, "measurement_conditions"): MeasurementConditions,
    (LiteratureClaim, "source"): LiteratureSource,
    (LiteratureClaim, "estimate"): ParameterEstimate,
}
_NESTED_LIST = {
    (Prediction, "inputs"): DataRef, (Prediction, "outputs"): OutputValue,
    (ValidationResult, "targets"): DataRef, (StructureCandidate, "sites"): SiteOccupancy,
    (LiteratureRecord, "claims"): Claim, (SiteOccupancy, "moments"): SiteMoment,
}


def from_dict(cls, d: Optional[dict]):
    """Rebuild a dataclass from plain data (JSON/YAML round trip)."""
    if d is None:
        return None
    kwargs = {}
    for k, v in d.items():
        if (cls, k) in _NESTED and isinstance(v, dict):
            v = from_dict(_NESTED[(cls, k)], v)
        elif (cls, k) in _NESTED_LIST and isinstance(v, list):
            v = [from_dict(_NESTED_LIST[(cls, k)], x) if isinstance(x, dict) else x for x in v]
        kwargs[k] = v
    return cls(**kwargs)
