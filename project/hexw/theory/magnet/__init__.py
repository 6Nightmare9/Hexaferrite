"""hexw.theory.magnet: fast macrospin model (equilibria, FMR, susceptibility, line widths).
Pure computation: no file access, no process pools at import. See docs/MAGNET_MODEL.md for conventions."""
from .core import (MU0, Gyro, MagnetModel, angles, as_gyro, energy, field_vector, map_literature_constants,
                   resolve_gamma, sphere_derivs)
from .demag import ellipsoid, no_shape, sphere, thin_disk_parallel_c, thin_disk_perp_c
from .equilibrium import (Equilibrium, EquilibriumSet, cone_angle_analytic, find_equilibria, hysteresis_loop,
                          trace_branch)
from .fmr import (H_res, ModeSet, f_res, llg_frequency_Hz, llg_omega, natural_fmr, smit_beljers_omega,
                  sweep_branches)
from .linewidth import (derivative_peak_to_peak, find_spectrum_peaks, fwhm, lorentzian_fwhm_from_dpp,
                        width_field_to_freq, width_freq_to_field)
from .polycrystal import (LEVEL_ENSEMBLE, LEVEL_IDEAL, Texture, crystal_chi, effective_medium_mixing,
                          intergranular_interaction, polycrystal_chi, porosity_correction, rotation_from_c_axis,
                          sample_demag_correction, schloemann_porosity_linewidth)
from .susceptibility import (apparent_mu_ellipsoid, apparent_mu_shape, kramers_kronig_real, mu_r, passivity_margin,
                             susceptibility)
from .uncertainty import monte_carlo
