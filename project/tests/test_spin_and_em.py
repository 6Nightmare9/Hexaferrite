import numpy as np
import pytest
from hexw.spin.multisublattice import SpinSystem,response,fit_exchange
from hexw.theory.magnet import Gyro
from hexw.theory.electromagnetic import dielectric,maxwell_garnett


def system(n,J,K,M):
    return SpinSystem(np.full(n,M),J,np.full(n,K),np.full(n,2.),np.full(n,.01),np.zeros((3,3)),{"source":"SYNTHETIC analytic test"})


def test_one_sublattice_reduces_to_kittel():
    s=system(1,np.zeros((1,1)),1e5,3e5)
    r=response(s,np.array([[0.,0.,1.]]),[0,0,.3],np.linspace(1e9,50e9,51))
    gamma=Gyro.from_g(2).gamma_over_2pi_Hz_per_T
    assert r["undamped_modes_Hz"]==pytest.approx([gamma*(.3+2e5/3e5)],rel=1e-12)
    assert np.min(r["chi_external"][:,0,0].imag)>=0


def test_antiferromagnet_optical_gap():
    J=-2e6; K=1e4; M=3e5
    s=system(2,np.array([[0.,J],[J,0.]]),K,M)
    r=response(s,np.array([[0.,0.,1.],[0.,0.,-1.]]),[0,0,0],[1e9,10e9])
    expected=Gyro.from_g(2).gamma_over_2pi_Hz_per_T*np.sqrt((2*K/M)*(2*K/M+2*abs(J)/M))
    assert r["undamped_modes_Hz"]==pytest.approx([expected,expected],rel=1e-12)


def test_response_refuses_a_nonstationary_state_even_with_positive_curvature():
    s=system(1,np.zeros((1,1)),1e5,3e5)
    with pytest.raises(ValueError,match="not stationary"):
        response(s,np.array([[.3,0,np.sqrt(.91)]]),[0,0,.3],[1e9])


def test_exchange_energy_mapping_and_rank_refusal():
    dot=np.array([1.,-1.,0.,.5])
    configs=np.array([[[0,0,1],[np.sqrt(1-z*z),0,z]] for z in dot])
    fit=fit_exchange(configs,2-.04*dot,[[(0,1)]])
    assert fit["J_eV_per_grouped_bond"]==pytest.approx([.04],abs=1e-14)
    with pytest.raises(ValueError,match="underdetermined"):
        fit_exchange(configs[:1],np.array([1.]),[[(0,1)]])


def test_debye_limits_passivity_and_porosity_endpoints():
    eps=dielectric([1e-6,1e18],epsilon_inf=3,relaxations=[{"delta_epsilon":7,"tau_s":1e-6}])
    assert eps.real==pytest.approx([10,3]) and np.all(eps.imag>=0)
    assert maxwell_garnett(eps,1,0)==pytest.approx(eps)
    assert maxwell_garnett(10,1,.5).real<10
