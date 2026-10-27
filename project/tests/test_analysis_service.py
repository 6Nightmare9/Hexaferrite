import json
import numpy as np
import pytest
from hexw.experiment.service import analyze
from hexw.experiment.readers import xrd
from hexw.registry.freeze import freeze_run
from test_freeze import completed


def test_native_mdi_trailer_is_explicit_and_extra_counts_rejected(root):
    p=root/"synthetic.mdi"
    header="synthetic\n5 .02 1 cu 1.54184 5.18 10\n"
    p.write_text(header+"1 2 3 4 5 6 7 8 9 10\n-1 0 0\n")
    xx,yy,meta=xrd(p)
    assert len(xx)==len(yy)==10 and meta["trailer"]==[-1.,0.,0.]
    p.write_text(header+"1 2 3 4 5 6 7 8 9 10\n-1 0 0 0 0 0 0\n")
    assert len(xrd(p)[1])==10
    p.write_text(header+"1 2 3 4 5 6 7 8 9 10 11\n")
    with pytest.raises(ValueError,match="count"):
        xrd(p)


def test_analysis_resume_is_bound_to_code_data_and_outputs(root):
    freeze_run(completed(root))
    p=root/"synthetic_xrd.csv"
    xx=np.linspace(20,40,501);y=100+500*np.exp(-((xx-30)/.2)**2)
    np.savetxt(p,np.column_stack([xx,y]),delimiter=",")
    cfg={"method":"xrd","file":str(p),"options":{"wavelength_angstrom":1.54}}
    out=analyze(cfg)
    assert analyze(cfg)==out
    assert json.loads((out/"analysis.json").read_text())["theory_parameters_updated"] is False
    (out/"analysis.png").write_bytes(b"changed")
    with pytest.raises(ValueError,match="artifact changed"):
        analyze(cfg)


def test_measured_derivative_is_not_differentiated_again(root):
    freeze_run(completed(root))
    xx=np.linspace(-3,3,2401);gamma=.6
    yy=-2*xx/(gamma/2)**2/(1+(xx/(gamma/2))**2)**2
    p=root/"synthetic_fmr_derivative.csv";np.savetxt(p,np.column_stack([xx,yy]),delimiter=",")
    out=analyze({"method":"fmr","file":str(p),"options":{"axis_unit":"T","signal_kind":"derivative","baseline":0}})
    result=json.loads((out/"analysis.json").read_text())
    assert result["status"]=="valid" and result["delta_pp"]==pytest.approx(gamma/np.sqrt(3),rel=.01)
