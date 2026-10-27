"""Typed physical quantities. Different field-like quantities are different kinds and never
convert implicitly. H [A/m], mu0*H [T], B (flux density) [T], B_hf [T], gamma [rad/(s*T)],
gamma/2pi [Hz/T] are separate kinds with explicit converters.

Conventions
  * CGS: 1 Oe <-> 1 G <-> 1e-4 T in vacuum (mu0*H and B). Hence 1 kOe <-> 0.1 T.
  * B_hf is traditionally quoted in tesla; a value in kOe is converted with 1 kOe = 0.1 T (vacuum-equivalent
    convention, no demagnetising/Lorentz corrections).
  * Magnetization M: SI A/m; CGS emu/cm3 with 1 emu/cm3 = 1000 A/m. In CGS the quantity that is compared with
    fields is 4*pi*M in GAUSS: 4*pi*M[G] = 4*pi*M[emu/cm3] (Kind.FOURPI_M, converters fourpiM_from_M / M_from_fourpiM).
    SI check: mu0*M[T] = 4*pi*M[G] * 1e-4.
  * Angular frequency omega [rad/s] = 2*pi*f [Hz]; gamma [rad/(s*T)]; gamma/2pi [Hz/T] (also MHz/Oe, GHz/kOe,
    both equal to 1e10 Hz/T per unit).
  * Constants come from scipy.constants (CODATA release recorded in CODATA_VERSION).

Stage 3 additions (only added, nothing renamed):
  * specific (mass) magnetization sigma: A*m^2/kg == emu/g (1 emu/g = 1 A*m^2/kg); M = sigma * density.
  * magnetic moment: muB (and J/T); energy density: J/m^3, MJ/m^3, erg/cm^3 (1 erg/cm^3 = 0.1 J/m^3);
    pressure: Pa, atm (101325 Pa); volume: angstrom^3; density: g/cm^3, kg/m^3.
  * energy per atom and energy per formula unit are DIFFERENT kinds (never mixed implicitly; explicit converter with
    the number of atoms per f.u.). Polarization mu0*M in tesla is its own kind (mu0M); reflection loss in dB.
  * fractions carry an explicit basis: mass / volume / mole / spectral_area / site (Kind.FRACTION stays for legacy data;
    kinds of different basis never add or compare).
  * kA/m and MA/m are valid for H_field and magnetization."""
from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from scipy import constants as _sc


def _codata_version() -> str:
    """CODATA release implemented by the installed scipy (private attribute, then the module docstring)."""
    try:
        from scipy.constants import _codata
        v = getattr(_codata, "_current_codata", None)
        if v:
            return str(v)
    except Exception:
        pass
    import re
    m = re.search(r"CODATA.{0,200}?(20[0-9][0-9])", _sc.__doc__ or "", re.S)
    return f"CODATA {m.group(1)}" if m else "unknown"


CODATA_VERSION = _codata_version()          # e.g. "CODATA 2022" for scipy >= 1.15
MU0 = _sc.mu_0                              # N/A^2 (NOT exactly 4*pi*1e-7 since the 2019 SI revision)
MU_B = _sc.physical_constants["Bohr magneton"][0]   # J/T
H_PLANCK = _sc.h                            # J s
HBAR = _sc.hbar                             # J s
E_CHARGE = _sc.e                            # C (exact)


class Kind(str, Enum):
    H_FIELD = "H_field"                    # magnetic field strength H
    MU0H = "mu0H"                          # mu0 * H
    B_FLUX = "B_flux"                      # magnetic flux density inside material
    B_HF = "B_hf"                          # hyperfine field
    GAMMA = "gamma"                        # gyromagnetic ratio gamma
    GAMMA_OVER_2PI = "gamma_over_2pi"      # gamma / 2pi
    FREQUENCY = "frequency"
    LENGTH = "length"
    ENERGY = "energy"
    TEMPERATURE = "temperature"
    MASS = "mass"
    VELOCITY = "velocity"                  # Moessbauer scale, isomer shift, etc.
    ANGLE = "angle"
    MAGNETIZATION = "magnetization"
    FRACTION = "fraction"
    TIME = "time"
    ANGULAR_FREQUENCY = "angular_frequency"  # omega, rad/s
    FOURPI_M = "fourpi_M"                  # CGS 4*pi*M, gauss
    DIMENSIONLESS = "dimensionless"
    # ---- stage 3 ----
    SPECIFIC_MAGNETIZATION = "specific_magnetization"   # sigma, A*m^2/kg == emu/g
    MAGNETIC_MOMENT = "magnetic_moment"                 # muB
    ENERGY_DENSITY = "energy_density"                   # J/m^3 (anisotropy constants K, ...)
    PRESSURE = "pressure"
    VOLUME = "volume"
    DENSITY = "density"                                 # mass density
    ENERGY_PER_ATOM = "energy_per_atom"
    ENERGY_PER_FORMULA_UNIT = "energy_per_formula_unit"
    MU0M = "magnetic_polarization_mu0M"                 # mu0*M in tesla
    REFLECTION_LOSS = "reflection_loss"                 # dB
    FRACTION_MASS = "fraction_mass"
    FRACTION_VOLUME = "fraction_volume"
    FRACTION_MOLE = "fraction_mole"
    FRACTION_SPECTRAL_AREA = "fraction_spectral_area"
    FRACTION_SITE = "fraction_site"                     # site occupancy


# kind -> {unit: factor to canonical}; first key is canonical
UNITS: dict[Kind, dict[str, float]] = {
    Kind.H_FIELD: {"A/m": 1.0, "Oe": 1000.0 / (4 * math.pi)},
    Kind.MU0H: {"T": 1.0, "mT": 1e-3, "G": 1e-4, "kG": 1e-1},
    Kind.B_FLUX: {"T": 1.0, "mT": 1e-3, "G": 1e-4},
    Kind.B_HF: {"T": 1.0, "mT": 1e-3, "kOe": 0.1},
    Kind.GAMMA: {"rad/(s*T)": 1.0, "rad/(s*Oe)": 1e4},
    Kind.GAMMA_OVER_2PI: {"Hz/T": 1.0, "GHz/T": 1e9, "MHz/T": 1e6, "MHz/Oe": 1e10, "GHz/kOe": 1e10},
    Kind.ANGULAR_FREQUENCY: {"rad/s": 1.0, "Mrad/s": 1e6, "Grad/s": 1e9},
    Kind.FOURPI_M: {"G": 1.0, "kG": 1e3},
    Kind.FREQUENCY: {"Hz": 1.0, "kHz": 1e3, "MHz": 1e6, "GHz": 1e9, "THz": 1e12},
    Kind.LENGTH: {"m": 1.0, "nm": 1e-9, "angstrom": 1e-10, "mm": 1e-3},
    Kind.ENERGY: {"J": 1.0, "eV": _sc.e, "meV": _sc.e * 1e-3},
    Kind.TEMPERATURE: {"K": 1.0},
    Kind.MASS: {"kg": 1.0, "g": 1e-3, "mg": 1e-6},
    Kind.VELOCITY: {"m/s": 1.0, "mm/s": 1e-3},
    Kind.ANGLE: {"rad": 1.0, "deg": math.pi / 180},
    Kind.MAGNETIZATION: {"A/m": 1.0, "kA/m": 1e3, "emu/cm3": 1e3},   # CGS: 1 emu/cm3 = 1000 A/m
    Kind.FRACTION: {"1": 1.0, "%": 0.01},
    Kind.TIME: {"s": 1.0, "h": 3600.0, "min": 60.0},
    Kind.DIMENSIONLESS: {"1": 1.0},
    # ---- stage 3 ----
    Kind.SPECIFIC_MAGNETIZATION: {"A*m^2/kg": 1.0, "emu/g": 1.0},      # 1 emu/g = 1e-3 A m^2 / 1e-3 kg
    Kind.MAGNETIC_MOMENT: {"muB": 1.0, "J/T": 1.0 / MU_B},
    Kind.ENERGY_DENSITY: {"J/m^3": 1.0, "kJ/m^3": 1e3, "MJ/m^3": 1e6, "erg/cm^3": 0.1},   # 1 erg/cm^3 = 0.1 J/m^3
    Kind.PRESSURE: {"Pa": 1.0, "kPa": 1e3, "MPa": 1e6, "GPa": 1e9, "bar": 1e5, "atm": 101325.0},
    Kind.VOLUME: {"angstrom^3": 1.0, "nm^3": 1e3},
    Kind.DENSITY: {"g/cm^3": 1.0, "kg/m^3": 1e-3},
    Kind.ENERGY_PER_ATOM: {"eV/atom": 1.0, "meV/atom": 1e-3},
    Kind.ENERGY_PER_FORMULA_UNIT: {"eV/f.u.": 1.0, "meV/f.u.": 1e-3},
    Kind.MU0M: {"T": 1.0, "mT": 1e-3},
    Kind.REFLECTION_LOSS: {"dB": 1.0},
    Kind.FRACTION_MASS: {"1": 1.0, "%": 0.01, "wt%": 0.01},
    Kind.FRACTION_VOLUME: {"1": 1.0, "%": 0.01, "vol%": 0.01},
    Kind.FRACTION_MOLE: {"1": 1.0, "%": 0.01, "mol%": 0.01},
    Kind.FRACTION_SPECTRAL_AREA: {"1": 1.0, "%": 0.01, "area%": 0.01},
    Kind.FRACTION_SITE: {"1": 1.0, "%": 0.01},
}
# kA/m and MA/m (added without touching the existing entries; the canonical unit stays the first key)
UNITS[Kind.H_FIELD].update({"kA/m": 1e3, "MA/m": 1e6, "kOe": 1e6 / (4 * math.pi)})
UNITS[Kind.MAGNETIZATION].update({"MA/m": 1e6})


class UnitError(ValueError):
    pass


def check_unit(kind: Kind | str, unit: str) -> Kind:
    k = Kind(kind)
    if not unit or not isinstance(unit, str):
        raise UnitError("unit is mandatory")
    if unit not in UNITS[k]:
        raise UnitError(f"unit {unit!r} is not valid for kind {k.value}; allowed: {list(UNITS[k])}")
    return k


@dataclass(frozen=True)
class Quantity:
    value: float
    unit: str
    kind: Kind

    def __post_init__(self):
        object.__setattr__(self, "kind", check_unit(self.kind, self.unit))
        if not math.isfinite(float(self.value)):
            raise UnitError("value must be finite")

    def canonical(self) -> "Quantity":
        canon = next(iter(UNITS[self.kind]))
        return Quantity(self.value * UNITS[self.kind][self.unit], canon, self.kind)

    def __add__(self, other):
        if not isinstance(other, Quantity) or other.kind != self.kind:
            raise UnitError(f"cannot add {getattr(other, 'kind', type(other))} to {self.kind}")
        o = other.canonical().value
        c = self.canonical()
        return Quantity(c.value + o, c.unit, self.kind)


def H(value, unit="A/m") -> Quantity: return Quantity(value, unit, Kind.H_FIELD)
def mu0H(value, unit="T") -> Quantity: return Quantity(value, unit, Kind.MU0H)
def B_flux(value, unit="T") -> Quantity: return Quantity(value, unit, Kind.B_FLUX)
def B_hf(value, unit="T") -> Quantity: return Quantity(value, unit, Kind.B_HF)
def gamma(value, unit="rad/(s*T)") -> Quantity: return Quantity(value, unit, Kind.GAMMA)
def gamma_over_2pi(value, unit="Hz/T") -> Quantity: return Quantity(value, unit, Kind.GAMMA_OVER_2PI)


def mu0H_from_H(h: Quantity) -> Quantity:
    """Explicit conversion mu0*H [T] from H [A/m]. NOT the flux density B inside matter."""
    if h.kind != Kind.H_FIELD:
        raise UnitError("expected H field")
    return Quantity(MU0 * h.canonical().value, "T", Kind.MU0H)


def H_from_mu0H(m: Quantity) -> Quantity:
    if m.kind != Kind.MU0H:
        raise UnitError("expected mu0H")
    return Quantity(m.canonical().value / MU0, "A/m", Kind.H_FIELD)


def gamma_over_2pi_from_gamma(g: Quantity) -> Quantity:
    if g.kind != Kind.GAMMA:
        raise UnitError("expected gamma")
    return Quantity(g.canonical().value / (2 * math.pi), "Hz/T", Kind.GAMMA_OVER_2PI)


def gamma_from_gamma_over_2pi(g: Quantity) -> Quantity:
    if g.kind != Kind.GAMMA_OVER_2PI:
        raise UnitError("expected gamma/2pi")
    return Quantity(g.canonical().value * 2 * math.pi, "rad/(s*T)", Kind.GAMMA)


def gamma_from_g(g: float) -> Quantity:
    """Gyromagnetic ratio gamma = g * mu_B / hbar [rad/(s*T)] (electron: pass |g|, e.g. 2.00231930436)."""
    return Quantity(g * MU_B / HBAR, "rad/(s*T)", Kind.GAMMA)


def angular_frequency_from_frequency(f: Quantity) -> Quantity:
    if f.kind != Kind.FREQUENCY:
        raise UnitError("expected frequency")
    return Quantity(2 * math.pi * f.canonical().value, "rad/s", Kind.ANGULAR_FREQUENCY)


def fourpiM_from_M(m: Quantity) -> Quantity:
    """CGS 4*pi*M in gauss from M (4*pi*M[G] = 4*pi*M[emu/cm3] = mu0*M[T]*1e4)."""
    if m.kind != Kind.MAGNETIZATION:
        raise UnitError("expected magnetization")
    return Quantity(4 * math.pi * m.canonical().value / 1e3, "G", Kind.FOURPI_M)


def M_from_fourpiM(q: Quantity) -> Quantity:
    if q.kind != Kind.FOURPI_M:
        raise UnitError("expected 4*pi*M")
    return Quantity(q.canonical().value * 1e3 / (4 * math.pi), "A/m", Kind.MAGNETIZATION)


# ---------------- stage 3 converters (all explicit; kinds never convert implicitly) ----------------
def specific_magnetization(value, unit="A*m^2/kg") -> Quantity: return Quantity(value, unit, Kind.SPECIFIC_MAGNETIZATION)
def magnetic_moment(value, unit="muB") -> Quantity: return Quantity(value, unit, Kind.MAGNETIC_MOMENT)
def energy_density(value, unit="J/m^3") -> Quantity: return Quantity(value, unit, Kind.ENERGY_DENSITY)
def mu0M(value, unit="T") -> Quantity: return Quantity(value, unit, Kind.MU0M)


def M_from_specific(sigma: Quantity, rho: Quantity) -> Quantity:
    """Volume magnetization M [A/m] = sigma [A*m^2/kg] * mass density [kg/m^3]."""
    if sigma.kind != Kind.SPECIFIC_MAGNETIZATION or rho.kind != Kind.DENSITY:
        raise UnitError("expected specific magnetization and density")
    return Quantity(sigma.canonical().value * rho.canonical().value * 1e3, "A/m", Kind.MAGNETIZATION)


def specific_from_M(m: Quantity, rho: Quantity) -> Quantity:
    if m.kind != Kind.MAGNETIZATION or rho.kind != Kind.DENSITY:
        raise UnitError("expected magnetization and density")
    return Quantity(m.canonical().value / (rho.canonical().value * 1e3), "A*m^2/kg", Kind.SPECIFIC_MAGNETIZATION)


def mu0M_from_M(m: Quantity) -> Quantity:
    """Polarization mu0*M [T] from M [A/m]."""
    if m.kind != Kind.MAGNETIZATION:
        raise UnitError("expected magnetization")
    return Quantity(MU0 * m.canonical().value, "T", Kind.MU0M)


def M_from_mu0M(q: Quantity) -> Quantity:
    if q.kind != Kind.MU0M:
        raise UnitError("expected mu0*M")
    return Quantity(q.canonical().value / MU0, "A/m", Kind.MAGNETIZATION)


def M_from_moment(moment_per_cell: Quantity, volume: Quantity) -> Quantity:
    """M [A/m] = (moment per cell in muB) * mu_B / (cell volume)."""
    if moment_per_cell.kind != Kind.MAGNETIC_MOMENT or volume.kind != Kind.VOLUME:
        raise UnitError("expected magnetic moment and volume")
    v_m3 = volume.canonical().value * 1e-30
    return Quantity(moment_per_cell.canonical().value * MU_B / v_m3, "A/m", Kind.MAGNETIZATION)


def energy_per_fu_from_per_atom(e: Quantity, atoms_per_fu: float) -> Quantity:
    """Explicit per-atom -> per-formula-unit conversion (SrFe18O27 has 46 atoms per f.u.)."""
    if e.kind != Kind.ENERGY_PER_ATOM:
        raise UnitError("expected energy per atom")
    if not atoms_per_fu > 0:
        raise UnitError("atoms per formula unit must be positive")
    return Quantity(e.canonical().value * atoms_per_fu, "eV/f.u.", Kind.ENERGY_PER_FORMULA_UNIT)


def energy_per_atom_from_per_fu(e: Quantity, atoms_per_fu: float) -> Quantity:
    if e.kind != Kind.ENERGY_PER_FORMULA_UNIT:
        raise UnitError("expected energy per formula unit")
    if not atoms_per_fu > 0:
        raise UnitError("atoms per formula unit must be positive")
    return Quantity(e.canonical().value / atoms_per_fu, "eV/atom", Kind.ENERGY_PER_ATOM)
