"""Sr-only single-crystal Ha(x,300 K) from a primary graph, with explicit transfer uncertainty.

The graph constrains a uniaxial anisotropy field, not damping, K2/K4 or atomic
exchange. K1=Ms*mu0Ha/2. Temperature extrapolation uses Ha~m^(p-1), p~N(3,.5).
"""
import json
from pathlib import Path
import numpy as np
from hexw.io import paths
from hexw.params import curie as C

REL=Path("fit/additions/Michioka_single_crystals_figure2.json")


def supported(composition):
    return len(composition)==1 and next(iter(composition)) in {"Mn","Mg","Ni"}


def sample(composition,temperatures_K,draws,Ms_A_m,m_reduced,*,seed):
    if not supported(composition):
        raise ValueError("Sr single-crystal branch supports one of Mn/Mg/Ni, total x<=2")
    doc=json.loads((paths.get_config().literature/REL).read_text(encoding="utf-8"))
    if doc.get("source_path"):
        from hexw.workflow import file_hash
        src=(paths.get_config().root/doc["source_path"]).resolve()
        if not src.is_relative_to(paths.get_config().root) or file_hash(src)!=doc["source_sha256"]:
            raise ValueError("Sr anisotropy primary source missing or checksum changed")
    sp,x=next(iter(composition.items()))
    points=sorted([[0,doc["host_Ha_T"]]]+doc["series"][sp])
    xp,ha=np.array(points).T
    central=float(np.interp(x,xp,ha))
    extrapolated=x<xp[0] or x>xp[-1]
    if x>xp[-1]:
        central=max(.05,ha[-1]+(x-xp[-1])*(ha[-1]-ha[-2])/(xp[-1]-xp[-2]))
    rng=np.random.default_rng([seed,913])
    # Sampling the drawing error does not account for a different specimen. Add
    # a separately declared material-transfer width; larger beyond the measured x.
    sigma_ln=.25 if not extrapolated else .5+.25*(x-xp[-1])
    H300=np.maximum(rng.normal(central,doc["digitization_std_T"],draws.n),.01)*np.exp(rng.normal(0,sigma_ln,draws.n))
    tc=np.maximum(C.tc_samples(composition,draws.tc_params),301.)
    m300=np.maximum(C.reduced_magnetization(300/tc,draws.shape_s,draws.form),1e-6)
    power=np.clip(rng.normal(3,.5,draws.n),1.1,5.)
    H=np.where(m_reduced>0,H300[None,:]*(np.maximum(m_reduced,0)/m300[None,:])**(power[None,:]-1),0.)
    K1=.5*Ms_A_m*H
    q2=np.clip(rng.normal(0,.1,draws.n),-.4,.4)
    q4=rng.normal(0,.01,draws.n)
    K2=K1*q2[None,:]*m_reduced**7
    K4=K1*q4[None,:]*m_reduced**18
    temps=np.asarray(temperatures_K)
    # Data are at room temperature only. Cover the unmeasured low-T high-order terms.
    extra=np.exp(rng.normal(0,.5,draws.n)[None,:]*np.minimum(abs(temps[:,None]-300)/300,2))
    values={"K1_J_m3":K1*extra,"K2_J_m3":K2*extra,"K3_J_m3":np.zeros_like(K1),"K4_J_m3":K4*extra,
            "mu0Ha_uniaxial_T":H*extra}
    status="sr_graph_extrapolation" if extrapolated or np.any(temps!=300) else "sr_graph_interpolation"
    source={"record":"SRSC-"+sp+"-Ha300","name":"mu0Ha_uniaxial","doi":None,
        "lit_id":"SRSC","source_sha256":doc.get("source_sha256"),"measurement_group":"SRSC_Ha300_"+sp,
        "url":doc["source_url"],"locator":{"page":"183 (PDF page 177)","figure":"2"},
        "role":"literature_fit","tier":2,"status":"graph_digitized",
        "basis":"nominal composition, Sr single crystals, room temperature; no laboratory targets used"}
    note=(f"Sr single-crystal primary graph; {sp} measured up to x={xp[-1]:g}. "
          f"Digitization std {doc['digitization_std_T']:g} T plus assumed specimen-transfer sigma_ln={sigma_ln:g}. "
          "K1=Ms*mu0Ha/2; K2/K4 and temperature law are priors, not measured. Not <=10% validated.")
    return values,{"status":status,"note":note,"source":source,"graph_range_x":[float(xp[0]),float(xp[-1])],
                   "temperature_data_K":300,"Ha300_central_T":central,"extrapolated_in_x":extrapolated}
