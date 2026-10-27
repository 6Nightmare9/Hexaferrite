"""Explicit dielectric dispersion and isotropic spherical-pore mixing; e^(-i omega t)."""
import numpy as np

EPS0=8.8541878128e-12


def dielectric(frequencies_Hz,*,epsilon_inf,relaxations=(),sigma_dc_S_m=0.):
    f=np.asarray(frequencies_Hz,float)
    if f.ndim!=1 or np.any(f<=0) or not np.isfinite(f).all() or epsilon_inf<1 or sigma_dc_S_m<0:
        raise ValueError("finite positive f, epsilon_inf>=1 and sigma>=0 required")
    omega=2*np.pi*f
    eps=np.full(f.shape,complex(epsilon_inf))+1j*sigma_dc_S_m/(EPS0*omega)
    for term in relaxations:
        delta,tau=term["delta_epsilon"],term["tau_s"]
        if delta<0 or tau<=0 or not np.isfinite([delta,tau]).all():
            raise ValueError("passive Debye relaxation requires delta>=0,tau>0")
        eps+=delta/(1-1j*omega*tau)
    return eps


def maxwell_garnett(host,inclusion,fraction):
    """Scalar isotropic spherical inclusions; not componentwise tensor homogenization."""
    if not np.isfinite(fraction) or not 0<=fraction<1:
        raise ValueError("inclusion volume fraction must be in [0,1)")
    h,i=np.broadcast_arrays(np.asarray(host,complex),np.asarray(inclusion,complex))
    if not np.isfinite(h).all() or not np.isfinite(i).all() or np.any(h.imag<-1e-12) or np.any(i.imag<-1e-12):
        raise ValueError("finite passive scalar constitutive functions required")
    den=i+2*h-fraction*(i-h)
    if np.any(abs(den)<1e-20):
        raise ValueError("mixing pole in the requested range")
    return h*(i+2*h+2*fraction*(i-h))/den
