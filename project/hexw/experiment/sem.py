"""Calibrated 2-D segmentation with explicit ROI/mask and edge censoring."""
import numpy as np
from scipy import ndimage as ndi
from hexw.workflow import summary


def analyze(image, *, um_per_pixel, foreground="bright", mask=None, min_area_px=16):
    if um_per_pixel is None or not np.isfinite(um_per_pixel) or um_per_pixel<=0:
        raise ValueError("a confirmed positive um_per_pixel calibration is mandatory")
    a=np.asarray(image,float)
    if a.ndim==3:
        a=a[:,:,:3].mean(2)
    if a.ndim!=2 or not np.isfinite(a).all() or foreground not in {"bright","dark"}:
        raise ValueError("2-D image or RGB and explicit bright/dark foreground required")
    if mask is None:
        hist,edges=np.histogram(a,bins=256)
        p=hist/hist.sum(); centers=(edges[:-1]+edges[1:])/2
        w=np.cumsum(p); m=np.cumsum(p*centers)
        score=(m[-1]*w-m)**2/np.maximum(w*(1-w),1e-30)
        threshold=centers[int(np.argmax(score[:-1]))]
        mask=a>threshold if foreground=="bright" else a<threshold
        method="Otsu connected components"
    else:
        mask=np.asarray(mask,bool)
        if mask.shape!=a.shape:
            raise ValueError("mask and image dimensions differ")
        threshold=None; method="user supplied segmentation mask"
    lab,n=ndi.label(mask)
    edge=set(np.unique(np.concatenate([lab[0],lab[-1],lab[:,0],lab[:,-1]])))
    sizes=np.bincount(lab.ravel()); retained=[]; grains=[]
    for k in range(1,n+1):
        if k in edge or sizes[k]<min_area_px:
            continue
        pts=np.column_stack(np.where(lab==k))*um_per_pixel
        eig=np.linalg.eigvalsh(np.cov(pts.T))
        grains.append({"label":k,"area_um2":float(sizes[k]*um_per_pixel**2),
            "equivalent_diameter_um":float(2*np.sqrt(sizes[k]/np.pi)*um_per_pixel),
            "aspect_ratio_moment":float(np.sqrt(max(eig[-1],1e-30)/max(eig[0],1e-30)))})
        retained.append(k)
    output=np.isin(lab,retained)
    return {"status":"2D_segmentation_requires_visual_review","method":method,"threshold":threshold,
        "um_per_pixel":um_per_pixel,"n_grains":len(grains),"grains":grains,
        "equivalent_diameter_um":summary([g["equivalent_diameter_um"] for g in grains]),
        "area_fraction":float(output.mean()),
        "limits":["touching grains merge unless supplied mask separates them", "edge objects and objects below min_area_px excluded",
                  "2-D surface statistics are not a 3-D distribution or XRD coherent-domain size",
                  "scale bars/annotations must be cropped before analysis", "EDS composition cannot identify crystallographic occupancy"]},output
