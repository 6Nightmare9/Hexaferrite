"""Uniform multi-sublattice LLG in SI, given EXPLICIT energy-density couplings.

E = -sum(i<j) J_ij m_i.m_j + sum K_i(1-m_iz^2) - B.sum M_i m_i
    + mu0/2 M_total.N.M_total.
J>0 is ferromagnetic; J<0 antiferromagnetic. M_i includes multiplicity.
This module does not invent exchange constants from a chemical formula.
"""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import block_diag
from hexw.theory.magnet import Gyro

MU0=4e-7*np.pi
J2=np.array([[0.,-1.],[1.,0.]])


@dataclass
class SpinSystem:
    M_A_m: np.ndarray
    J_J_m3: np.ndarray
    K_J_m3: np.ndarray
    g: np.ndarray
    alpha: np.ndarray
    N: np.ndarray
    provenance: dict

    def __post_init__(self):
        for name in ("M_A_m","J_J_m3","K_J_m3","g","alpha","N"):
            setattr(self,name,np.asarray(getattr(self,name),float))
        n=len(self.M_A_m)
        if not 1<=n<=36 or any(getattr(self,k).shape!=(n,) for k in ("M_A_m","K_J_m3","g","alpha")):
            raise ValueError("1..36 sublattices with matching M,K,g,alpha required")
        if self.J_J_m3.shape!=(n,n) or not np.allclose(self.J_J_m3,self.J_J_m3.T) or np.any(np.diag(self.J_J_m3)!=0):
            raise ValueError("symmetric zero-diagonal exchange matrix required; each bond counted once")
        if self.N.shape!=(3,3) or not np.allclose(self.N,self.N.T) or np.linalg.eigvalsh(self.N).min()<-1e-12:
            raise ValueError("symmetric positive-semidefinite demagnetizing tensor required")
        if any(not np.isfinite(getattr(self,k)).all() for k in ("M_A_m","J_J_m3","K_J_m3","g","alpha","N")):
            raise ValueError("all spin parameters must be finite")
        if np.any(self.M_A_m<=0) or np.any(self.g<=0) or np.any(self.alpha<0) or not self.provenance:
            raise ValueError("positive M,g, nonnegative alpha and parameter provenance required")

    def energy(self,m,B):
        total=np.sum(self.M_A_m[:,None]*m,axis=0)
        return (-.5*np.sum(self.J_J_m3*(m@m.T))+np.sum(self.K_J_m3*(1-m[:,2]**2))
                -np.dot(B,total)+.5*MU0*total@self.N@total)

    def gradient(self,m,B):
        total=np.sum(self.M_A_m[:,None]*m,axis=0)
        out=-self.J_J_m3@m-self.M_A_m[:,None]*np.asarray(B)[None,:]
        out[:,2]-=2*self.K_J_m3*m[:,2]
        return out+MU0*self.M_A_m[:,None]*(self.N@total)[None,:]


def tangent(m):
    frames=[]
    for z in m:
        axis=np.eye(3)[np.argmin(abs(z))]
        u=np.cross(axis,z); u/=np.linalg.norm(u)
        frames.append(np.column_stack([u,np.cross(z,u)]))
    return block_diag(*frames)


def equilibrium(system,B,*,starts=12,seed=20261026):
    n=len(system.M_A_m); rng=np.random.default_rng(seed)
    def unpack(a):
        a=a.reshape(n,2); th,ph=a[:,0],a[:,1]
        return np.column_stack([np.sin(th)*np.cos(ph),np.sin(th)*np.sin(ph),np.cos(th)])
    trials=[]
    scale=max(np.max(abs(system.J_J_m3)),np.max(abs(system.K_J_m3)),np.max(system.M_A_m)*np.linalg.norm(B),1.)
    for i in range(starts):
        if i<2:
            theta=np.zeros(n) if i==0 else np.where(np.arange(n)%2, np.pi,0)
        else:
            theta=np.arccos(rng.uniform(-1,1,n))
        initial=np.column_stack([theta,rng.uniform(0,2*np.pi,n)]).ravel()
        opt=minimize(lambda a:system.energy(unpack(a),B)/scale,initial,method="BFGS",options={"gtol":1e-8,"maxiter":1500})
        m=unpack(opt.x)
        residual=np.linalg.norm(np.cross(m,system.gradient(m,B)),axis=1).max()/scale
        if residual<2e-5:
            trials.append((system.energy(m,B),m,residual))
    if not trials:
        raise ValueError("no stationary spin state found")
    energy,m,residual=min(trials,key=lambda p:p[0])
    h=stiffness(system,m,B)
    if np.linalg.eigvalsh(h).min()<-scale*1e-6:
        raise ValueError("lowest found state is unstable; increase starts")
    return m,{"energy_J_m3":float(energy),"stationarity_relative":float(residual),
              "n_stationary_starts":len(trials),"scope":"lowest of multistart candidates; global minimum not guaranteed"}


def stiffness(system,m,B):
    """Exact constrained Cartesian Hessian projected onto local tangent planes."""
    n=len(m); t=tangent(m); moments=system.M_A_m
    h=np.kron(-system.J_J_m3,np.eye(3))
    for i in range(n):
        h[3*i+2,3*i+2]-=2*system.K_J_m3[i]
        for j in range(n):
            h[3*i:3*i+3,3*j:3*j+3]+=MU0*moments[i]*moments[j]*system.N
    lag=np.sum(m*system.gradient(m,B),axis=1)
    return t.T@(h-np.diag(np.repeat(lag,3)))@t


def response(system,m,B,frequencies_Hz):
    m=np.asarray(m,float); B=np.asarray(B,float)
    if m.shape!=(len(system.M_A_m),3) or not np.allclose(np.linalg.norm(m,axis=1),1,atol=1e-7):
        raise ValueError("one unit magnetization per sublattice required")
    if B.shape!=(3,) or not np.isfinite(B).all():
        raise ValueError("finite three-component applied B required")
    energy_scale=max(np.max(abs(system.J_J_m3)),np.max(abs(system.K_J_m3)),
                     np.max(system.M_A_m)*np.linalg.norm(B),MU0*system.M_A_m.sum()**2,1.)
    torque=np.linalg.norm(np.cross(m,system.gradient(m,B)),axis=1).max()
    if torque>2e-5*energy_scale:
        raise ValueError("magnetizations are not stationary; find equilibrium before computing a response")
    h=stiffness(system,m,B)
    scale=max(np.linalg.norm(h),1.)
    if np.linalg.eigvalsh(h).min()<-1e-7*scale:
        raise ValueError("response about an unstable state is not a passive prediction")
    t=tangent(m); n=len(m)
    blocks=[]; undamped=[]
    for M,g,a in zip(system.M_A_m,system.g,system.alpha):
        gamma=Gyro.from_g(float(g)).gamma_rad_per_s_T
        blocks.append(gamma/M/(1+a*a)*(J2-a*np.eye(2)))
        undamped.append(gamma/M*J2)
    L=block_diag(*blocks); A=L@h
    A0=block_diag(*undamped)@h
    eig=np.linalg.eigvals(A0)
    modes=np.sort(eig.imag[eig.imag>1e-3]/(2*np.pi))
    # Uniform applied h in A/m drives every sublattice; observe total M in A/m.
    C=np.hstack([M*t[3*i:3*i+3,2*i:2*i+2] for i,M in enumerate(system.M_A_m)])
    force=-MU0*L@C.T
    f=np.asarray(frequencies_Hz,float)
    if f.ndim!=1 or not np.isfinite(f).all() or np.any(f<=0):
        raise ValueError("positive finite frequency grid required")
    chi=np.array([C@np.linalg.solve(-2j*np.pi*freq*np.eye(2*n)-A,force) for freq in f])
    return {"undamped_modes_Hz":modes,"chi_external":chi,"dynamic_eigenvalues_s_inv":np.linalg.eigvals(A),
            "definition":"uniform coupled-sublattice LLG; exchange parameters and provenance supplied explicitly"}


def fit_exchange(orientations,energies_eV,pairs):
    """Map completed DFT energies to E0-sum J_ij m_i.m_j (J per listed bond, eV).

    Pairs may group equivalent bonds: each list contains (i,j) bonds whose sum is
    one column. Rank deficiency is refused; enough independent spin flips required.
    """
    m=np.asarray(orientations,float); energy=np.asarray(energies_eV,float)
    if m.ndim!=3 or m.shape[2]!=3 or len(m)!=len(energy) or not np.allclose(np.linalg.norm(m,axis=2),1):
        raise ValueError("unit spin configurations and matching DFT energies required")
    X=np.column_stack([np.ones(len(m))]+[-sum(np.sum(m[:,i]*m[:,j],axis=1) for i,j in group) for group in pairs])
    if np.linalg.matrix_rank(X)<X.shape[1] or len(energy)<=X.shape[1]:
        raise ValueError("exchange mapping is underdetermined or has no residual degrees of freedom")
    coef,_,_,_=np.linalg.lstsq(X,energy,rcond=None)
    residual=energy-X@coef
    covariance=np.linalg.pinv(X.T@X)*(residual@residual)/(len(energy)-len(coef))
    return {"E0_eV":float(coef[0]),"J_eV_per_grouped_bond":coef[1:].tolist(),
            "residual_rms_eV":float(np.sqrt(np.mean(residual**2))),"condition":float(np.linalg.cond(X)),
            "local_covariance":covariance.tolist(),"convention":"E0-sum J m_i.m_j; each listed bond counted once; classical unit spins"}
