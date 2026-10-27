"""Peak metrology and a constrained structural-pattern screening fit (not full Rietveld)."""
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks, peak_widths
from hexw.workflow import summary


def analyze(two_theta, counts, *, wavelength_angstrom, instrument_fwhm_deg=None, prominence_fraction=.03):
    x, y = np.asarray(two_theta, float), np.asarray(counts, float)
    if len(x) != len(y) or len(x)<8 or np.any(np.diff(x)<=0) or not np.isfinite([x,y]).all():
        raise ValueError("finite monotonic two-theta/count arrays required")
    if wavelength_angstrom<=0 or np.any((x<=0)|(x>=180)) or np.any(y<0):
        raise ValueError("XRD requires positive wavelength, 0<2theta<180 deg, nonnegative counts")
    if instrument_fwhm_deg is not None and (not np.isfinite(instrument_fwhm_deg) or instrument_fwhm_deg<0):
        raise ValueError("instrumental FWHM must be finite and nonnegative")
    idx, props = find_peaks(y, prominence=prominence_fraction*np.ptp(y), distance=3)
    widths, _, left, right = peak_widths(y, idx, rel_height=.5)
    peaks=[]
    for k,p in enumerate(idx):
        theta=np.deg2rad(x[p]/2)
        width=np.interp(right[k],np.arange(len(x)),x)-np.interp(left[k],np.arange(len(x)),x)
        beta=None
        if instrument_fwhm_deg is not None and width>instrument_fwhm_deg:
            beta=np.deg2rad(np.sqrt(width**2-instrument_fwhm_deg**2))
        peaks.append({"two_theta_deg":float(x[p]),"d_angstrom":float(wavelength_angstrom/(2*np.sin(theta))),
            "intensity_counts":float(y[p]),"observed_fwhm_deg":float(width),
            "scherrer_D_nm":float(.9*wavelength_angstrom/(beta*np.cos(theta))/10) if beta else None})
    return {"method":"XRD_peak_screening","peaks":peaks,"n_peaks":len(peaks),
        "scherrer_note":"K=0.9 assumption; Gaussian instrumental-width subtraction; size/strain/overlap unresolved. Crystallites, not SEM grains.",
        "phase_fractions":None,"rietveld_status":"not_performed",
        "required_for_quantitative_refinement":["instrument profile/standard", "candidate structures of all phases", "background/texture model"]}


def screen_structure(x, y, structure, wavelength_angstrom):
    """Refine scale/background/zero/one peak width against a fixed crystal structure.

    Diffraction intensities are computed by pymatgen. Atomic coordinates/occupancy,
    phase mass fractions, strain and crystallite sizes are NOT refined here.
    """
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    pat=XRDCalculator(wavelength=wavelength_angstrom).get_pattern(structure, two_theta_range=(float(x[0])-1,float(x[-1])+1))
    if not len(pat.x):
        raise ValueError("structure has no reflections in range")
    x,y=np.asarray(x),np.asarray(y)
    z=(x-x.mean())/np.ptp(x)
    def curve(p):
        a,b,c,shift,w=p
        u=(x[:,None]-pat.x[None,:]-shift)/w
        profile=np.exp(-4*np.log(2)*u*u)@pat.y
        return a*profile+b+c*z
    sigma=np.sqrt(np.maximum(y,1))
    p0=[max(np.ptp(y)/100,.01),float(np.percentile(y,10)),0,0,.15]
    opt=least_squares(lambda p:(curve(p)-y)/sigma,p0,bounds=([0,0,-np.inf,-.5,.02],[np.inf,np.inf,np.inf,.5,2]),max_nfev=250)
    fit=curve(opt.x)
    rwp=float(np.sqrt(np.sum(((fit-y)/sigma)**2)/np.sum((y/sigma)**2)))
    return {"method":"fixed_structure_pattern_screening","success":bool(opt.success),"Rwp":rwp,
        "zero_shift_deg":float(opt.x[3]),"profile_fwhm_deg":float(opt.x[4]),
        "phase_purity_or_mass_fraction":None,"scope":"single fixed structure; fit quality is not an accuracy score or phase purity"}, fit
