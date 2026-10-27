"""End-to-end contracts with SYNTHETIC parameters and independent Kittel limits."""
import json
import numpy as np
import pytest

from hexw.io import paths
from hexw.theory.request import ForecastRequest
from hexw.theory import pipeline as P
from hexw.theory.magnet import Gyro


def request(**kwargs):
    return ForecastRequest(**({"composition": {"Co": 1, "Zn": 1}, "temperatures_K": [300],
                              "fields_T": [0, .3], "frequencies_GHz": np.linspace(.1, 40, 81).tolist(),
                              "g": 2, "alpha": .02, "n_draws": 4, "n_prior": 1000,
                              "equilibrium_grid": 40} | kwargs))


def synthetic(req):
    n = req.n_draws
    return {"Ms_A_m": np.full((1,n), 3e5), "m_reduced": np.ones((1,n)), "Tc_K": np.full(n, 700.),
            "K1_J_m3": np.full((1,n), 1e5), **{f"K{k}_J_m3": np.zeros((1,n)) for k in (2,3,4)}}, {
                "anisotropy_status": "SYNTHETIC", "scope": "SYNTHETIC", "unimplemented": []}


@pytest.mark.parametrize("kwargs", [{"alpha": 0}, {"g": float("nan")}, {"composition": {"Al": 1}},
                                   {"fields_T": [1,0]}, {"frequencies_GHz": [0,1,2]}, {"shape": "cell"}])
def test_bad_requests_rejected(kwargs):
    with pytest.raises(ValueError):
        request(**kwargs)


def test_synthetic_axial_case_matches_kittel_and_preserves_tensor_passivity():
    req = request()
    params, _ = synthetic(req)
    data, info = P.compute_case(req, params, 0, 1)
    expected = Gyro.from_g(2).gamma_over_2pi_Hz_per_T * (.3 + 2e5/3e5)
    assert np.allclose(data["lowest_ground_mode_Hz"], expected, rtol=1e-8)
    assert info["n_valid"] == 4 and info["external_passivity_min"] >= -1e-9
    assert np.all(data["mean_chi_external"][:,0,0].imag > 0)


def test_pipeline_resume_and_artifact_tampering(root, monkeypatch):
    monkeypatch.setattr(P, "material_draws", synthetic)
    out = P.run_forecast(request(), progress=lambda _: None)
    assert json.loads((out/"summary.json").read_text())["accuracy"]["achieved"] is None
    monkeypatch.setattr(P, "compute_case", lambda *a: pytest.fail("completed cases must be reused"))
    assert P.run_forecast(request(), out, progress=lambda _: None) == out
    (out/"fmr_vs_field.csv").write_text("tampered")
    with pytest.raises(ValueError, match="artifact changed"):
        P.run_forecast(request(), out)


def test_above_curie_is_not_silently_zero_fmr():
    req = request()
    params, _ = synthetic(req)
    params["m_reduced"][:] = 0
    data, info = P.compute_case(req, params, 0, 0)
    assert info["status"] == "not_predicted"
    assert np.isnan(data["lowest_ground_mode_Hz"]).all()


def test_resume_after_interruption_reuses_completed_case(root,monkeypatch):
    monkeypatch.setattr(P,"material_draws",synthetic)
    real=P.compute_case
    def interrupt(req,params,ti,hi):
        if hi==1:
            raise RuntimeError("synthetic power interruption")
        return real(req,params,ti,hi)
    monkeypatch.setattr(P,"compute_case",interrupt)
    with pytest.raises(RuntimeError,match="interruption"):
        P.run_forecast(request(),progress=lambda _:None)
    out=next((root/"runs"/"predictions").iterdir())
    calls=[]
    def resume(req,params,ti,hi):
        calls.append(hi)
        return real(req,params,ti,hi)
    monkeypatch.setattr(P,"compute_case",resume)
    P.run_forecast(request(),out,progress=lambda _:None)
    assert calls==[1]
    assert json.loads((out/"checkpoint.json").read_text())["status"]=="complete"


def test_material_cache_is_content_bound_and_tamper_checked(root,monkeypatch):
    calls=[]
    def compute(req):
        calls.append(1)
        return synthetic(req)
    monkeypatch.setattr(P,"_compute_material_draws",compute)
    a,_=P.material_draws(request())
    b,_=P.material_draws(request())
    assert len(calls)==1 and np.array_equal(a["Ms_A_m"],b["Ms_A_m"])
    cache=next((root/"runs"/"parameter_cache").iterdir())
    (cache/"material.json").write_text("modified")
    with pytest.raises(ValueError,match="cache modified"):
        P.material_draws(request())
