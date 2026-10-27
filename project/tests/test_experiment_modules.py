"""Synthetic standards with known geometry/spectral parameters; no laboratory labels."""
import numpy as np
import pytest
from hexw.experiment.xrd import analyze as xrd
from hexw.experiment.mossbauer import fit_alternatives, sextet, COEFF
from hexw.experiment.sem import analyze as sem
from hexw.experiment.eds import analyze as eds


def test_bragg_spacing_and_instrument_correction():
    x=np.linspace(20,40,2001); y=100+1000*np.exp(-4*np.log(2)*((x-30)/.2)**2)
    r=xrd(x,y,wavelength_angstrom=1.54,instrument_fwhm_deg=.1)
    assert len(r["peaks"])==1
    p=r["peaks"][0]
    assert p["d_angstrom"]==pytest.approx(1.54/(2*np.sin(np.pi/12)),rel=1e-5)
    assert p["observed_fwhm_deg"]==pytest.approx(.2,rel=.001)
    assert p["scherrer_D_nm"]>0 and r["phase_fractions"] is None


def test_missing_mossbauer_calibration_never_yields_Bhf():
    result,curve=fit_alternatives(np.linspace(-12,12,512),np.ones(512)*10000,calibrated=False)
    assert result["status"]=="calibration_required" and curve is None


def test_iron_sextet_recovers_synthetic_field_and_shift():
    v=np.linspace(-10,10,1001)
    # alpha-Fe ~33 T has outer lines near +-5.3 mm/s, an independent dimensional check.
    assert abs(COEFF[0]*33)==pytest.approx(5.3,rel=.02)
    y=100000-1000*sextet(v,45,.35,.03,.3)
    r,curve=fit_alternatives(v,y,calibrated=True,n_components=[1],starts=2)
    c=r["alternatives"][0]["components"][0]
    assert c["Bhf_T"]==pytest.approx(45,abs=.01)
    assert c["delta_mm_s"]==pytest.approx(.35,abs=.001)
    assert np.max(abs(curve-y))<.1


def test_calibrated_sem_disks_and_edge_censoring():
    yy,xx=np.mgrid[:160,:200]
    mask=(xx-50)**2+(yy-80)**2<=20**2
    mask|=(xx-140)**2+(yy-80)**2<=15**2
    mask|=(xx-1)**2+(yy-1)**2<=20**2
    r,out=sem(mask.astype(float),um_per_pixel=.1,mask=mask)
    assert r["n_grains"]==2
    assert sorted(g["equivalent_diameter_um"] for g in r["grains"])==pytest.approx([3,4],rel=.01)
    with pytest.raises(ValueError):
        sem(mask,um_per_pixel=0)


def test_quantified_eds_stoichiometry_preserves_oxygen_without_imposing_it():
    from pymatgen.core import Element
    stoich={"Sr":1,"Ni":2,"Fe":16,"O":26.5}
    total=sum(float(Element(s).atomic_mass)*n for s,n in stoich.items())
    rows=[{"element":s,"weight_percent":100*float(Element(s).atomic_mass)*n/total} for s,n in stoich.items()]
    r=eds(rows)
    assert r["atomic_ratios"]==pytest.approx(stoich)
    assert r["ratio_scenario_intervals"] is None
