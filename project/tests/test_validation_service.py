import json
import numpy as np
import pytest
from hexw.validate.service import _match


def test_comparison_does_not_replace_peak_with_eigenfrequency(tmp_path):
    req={"composition":{"Co":1.,"Zn":1.},"temperatures_K":[300.],"fields_T":[0.],
         "theta_deg":0.,"phi_deg":0.,"shape":"sphere","g":2.,"alpha":.02}
    target={"composition":req["composition"],"T_K":300.,"mu0H_T":0.,"observable":"lowest_ground_mode_Hz",
            "unit":"Hz",**{k:req[k] for k in ("theta_deg","phi_deg","shape","g","alpha")},
            "quantity_definition":"external_chi_xx_absorption_peak"}
    with pytest.raises(ValueError,match="definition"):
        _match(target,req,{},tmp_path)
    target["quantity_definition"]="undamped_lowest_ground_eigenfrequency"
    np.savez(tmp_path/"T000_H0000.npz",lowest_ground_mode_Hz=[1e9,2e9,3e9])
    assert _match(target,req,{},tmp_path)["median"]==2e9
    target.pop("shape")
    with pytest.raises(ValueError,match="shape"):
        _match(target,req,{},tmp_path)


def test_signed_anisotropy_and_units_are_required(tmp_path):
    req={"composition":{"Ni":2.},"temperatures_K":[300.]}
    target={"composition":{"Ni":2.},"T_K":300.,"observable":"K1_J_m3","unit":"J/m^3"}
    p={"K1_J_m3":np.array([[1.,2.,3.]])}
    with pytest.raises(ValueError,match="signed"):
        _match(target,req,p,tmp_path)


def test_neutron_ordered_moment_and_impure_vsm_are_not_intrinsic_ms_targets():
    from hexw.validate.literature import admissible_magnetization
    assert admissible_magnetization("Ms_VSM_impurity_corrected","A*m^2/kg")
    assert not admissible_magnetization("Ms_NPD","A*m^2/kg")
    assert not admissible_magnetization("Ms_VSM","A*m^2/kg")


def test_legacy_predictions_cannot_be_relabelled_by_the_validation_entrypoint(tmp_path):
    from hexw.validate.service import check_target_origin
    doc={"source_role":"legacy_prediction","targets":[]}
    with pytest.raises(ValueError,match="source_role"):
        check_target_origin(doc,tmp_path/"targets.json",tmp_path)
    doc["source_role"]="laboratory_validation"
    with pytest.raises(ValueError,match="not laboratory"):
        check_target_origin(doc,tmp_path/"runs/predictions/a/targets.json",tmp_path)
