"""SYNTHETIC parameters for the magnet-model tests. They are NOT material data and carry no literature meaning."""
import numpy as np

from hexw.theory.magnet import Gyro

SYNTH_MS = 3.0e5            # A/m  (synthetic)
SYNTH_K1 = 1.0e5            # J/m^3 (synthetic)
SYNTH_K2 = 2.0e4            # J/m^3 (synthetic)
SYNTH_K3 = 5.0e3            # J/m^3 (synthetic)
SYNTH_K4 = -1.0e4           # J/m^3 (synthetic)
SYNTH_G = 2.0               # g used as an explicit test input only
GY = Gyro.from_g(SYNTH_G)
GAMMA_2PI = GY.gamma_over_2pi_Hz_per_T


def rand_rotation(rng):
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    a, b, c, d = q
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])
