import json
from types import SimpleNamespace
import numpy as np
import pytest
from hexw.params import sr_anisotropy as S


def test_direct_Sr_graph_branch_preserves_field_and_discloses_extrapolation(root,monkeypatch):
    p=root/"literature"/S.REL;p.parent.mkdir(parents=True)
    p.write_text(json.dumps({"host_Ha_T":2.,"series":{"Ni":[[1.,1.5]]},"digitization_std_T":.01,"source_url":"synthetic"}))
    n=20000
    d=SimpleNamespace(n=n,tc_params={},shape_s=np.ones(n)*.5,form="kuzmin")
    monkeypatch.setattr(S.C,"tc_samples",lambda *a:np.full(n,800.))
    monkeypatch.setattr(S.C,"reduced_magnetization",lambda *a:np.full(n,.8))
    ms=np.full((1,n),3e5); m=np.full((1,n),.8)
    values,info=S.sample({"Ni":1},[300],d,ms,m,seed=42)
    # K1/Ms cancels the uncertain moment scale exactly at the graph temperature.
    assert np.median(2*values["K1_J_m3"]/ms)==pytest.approx(1.5,rel=.01)
    assert info["status"]=="sr_graph_interpolation"
    assert info["source"]["status"]=="graph_digitized"
    _,extended=S.sample({"Ni":2},[300],d,ms,m,seed=42)
    assert extended["extrapolated_in_x"] and extended["status"]=="sr_graph_extrapolation"


def test_direct_branch_never_silently_mixes_metals():
    assert S.supported({"Mg":1}) and not S.supported({"Co":1,"Zn":1})
    assert not S.supported({"Mn":1,"Ni":1})
