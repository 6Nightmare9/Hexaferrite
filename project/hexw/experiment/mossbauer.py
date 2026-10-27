"""57Fe thin-absorber alternatives. No automatic identification of seven crystal sites."""
import numpy as np
from scipy.optimize import least_squares

# Nuclear Zeeman transition velocities: (-g_e*m_e+g_g*m_g)*mu_N*c/Egamma.
# Constants are an explicit approximate convention, not a velocity calibration.
GG, GE = .1806, -.1035
MUN=5.0507837393e-27
EGAMMA=14412.5*1.602176634e-19
TRANSITIONS=[(-1.5,-.5),(-.5,-.5),(.5,-.5),(-.5,.5),(.5,.5),(1.5,.5)]
COEFF=np.sort(np.array([(-GE*me+GG*mg)*MUN*299792458/EGAMMA*1000 for me,mg in TRANSITIONS]))


def sextet(v, Bhf_T, delta_mm_s, epsilon_mm_s, fwhm_mm_s, ratio24=2.):
    positions=Bhf_T*COEFF+delta_mm_s+epsilon_mm_s*np.array([1,-1,-1,-1,-1,1])
    weights=np.array([3,ratio24,1,1,ratio24,3.],float)
    weights/=weights.sum()
    half=fwhm_mm_s/2
    return np.sum(weights[None,:]*half/np.pi/((np.asarray(v)[:,None]-positions[None,:])**2+half**2),axis=1)


def fit_alternatives(v, counts, *, calibrated, n_components=(1,2,3), seed=20261026, starts=3):
    v,y=np.asarray(v,float),np.asarray(counts,float)
    if len(v)!=len(y) or np.any(np.diff(v)<=0) or not np.isfinite([v,y]).all() or np.any(y<=0):
        raise ValueError("finite monotonic velocity/counts with positive counts required")
    if not calibrated:
        return {"status":"calibration_required","physical_parameters":None,
                "reason":"alpha-Fe velocity calibration/zero and source convention are not confirmed"}, None
    rng=np.random.default_rng(seed)
    sigma=np.sqrt(y)
    alternatives=[]
    curves=[]
    for n in n_components:
        if not 1<=n<=7 or len(v)<=2+5*n:
            raise ValueError("1..7 components and enough observations required")
        def model(p):
            out=p[0]+p[1]*(v-v.mean())
            for a,b,d,e,w in p[2:].reshape(n,5):
                out=out-a*sextet(v,b,d,e,w)
            return out
        b0=float(np.percentile(y,95)); depth=max(float(np.ptp(y)),1.)
        lower=[0,-np.inf]+[0,5,-1.5,-1,.08]*n
        upper=[np.inf,np.inf]+[b0*10,65,2,1,3]*n
        best=None
        for _ in range(starts):
            p=[b0,0]
            for j in range(n):
                p.extend([depth/n,rng.uniform(35,53),rng.uniform(.15,.6),0,.4])
            fit=least_squares(lambda p:(model(p)-y)/sigma,p,bounds=(lower,upper),max_nfev=1200)
            if best is None or np.sum(fit.fun**2)<np.sum(best.fun**2):
                best=fit
        chi2=float(best.fun@best.fun); k=len(best.x); dof=len(v)-k
        sv=np.linalg.svd(best.jac,compute_uv=False)
        condition=float(sv[0]/max(sv[-1],1e-300))
        cov=np.linalg.pinv(best.jac.T@best.jac)*chi2/max(dof,1)
        pars=best.x[2:].reshape(n,5)
        areas=pars[:,0]/max(pars[:,0].sum(),1e-300)
        components=[dict(zip(("area_amplitude","Bhf_T","delta_mm_s","epsilon_mm_s","fwhm_mm_s"),map(float,p)),
                         spectral_area_fraction=float(a)) for p,a in zip(pars,areas)]
        alternatives.append({"n_sextets":n,"components":components,"chi2":chi2,"reduced_chi2":chi2/dof,
            "BIC":chi2+k*np.log(len(v)),"AIC":chi2+2*k,"success":bool(best.success),
            "jacobian_condition":condition,"local_parameter_std":np.sqrt(np.maximum(np.diag(cov),0)).tolist(),
            "identifiability":"poor" if condition>1e8 else "local_only"})
        curves.append(model(best.x))
    i=int(np.argmin([a["BIC"] for a in alternatives]))
    return {"status":"thin_absorber_fit","selected_by_BIC":i,"alternatives":alternatives,
        "conventions":{"g_ground":GG,"g_excited":GE,"gamma_energy_eV":14412.5,"intensity_ratio":"3:2:1:1:2:3 (random powder)"},
        "limits":["local covariance is not a global uncertainty interval", "thickness/relaxation/distributed hyperfine fields not fitted",
                  "spectral area is not atomic occupancy without recoil-free factors", "sextets are not automatically assigned to W1..W7",
                  "Bhf is a nuclear hyperfine field, not the FMR anisotropy field"]},curves[i]
