import numpy as np
PP = {  # GBRV 1.x USPP (SSSP-efficiency 1.3.0 files for Sr, Co, Zn, Mn; GBRV site for Fe, O)
 'Fe': ('fe_pbe_v1.5.uspp.F.UPF', 55.847), 'O': ('o_pbe_v1.2.uspp.F.UPF', 15.999),
 'Sr': ('Sr_pbe_v1.uspp.F.UPF', 87.62), 'Zn': ('Zn_pbe_v1.uspp.F.UPF', 65.38),
 'Co': ('Co_pbe_v1.2.uspp.F.UPF', 58.933), 'Mn': ('mn_pbe_v1.5.uspp.F.UPF', 54.938)}
import os
PSEUDO_DIR = 'D:/Work/Hexx/tools/pseudo/use'
def write_pw(fname, struct, labels, spins, prefix, ecut, ecutrho, kgrid, U=4.0, maxstep=5, nbnd=None, extra_sys="", conv="1e-6", beta=0.3):
    """labels[i]: species label (e.g. Fe1), spins[label]: starting_magnetization (-1..1), elements from label[:2] or label[0]"""
    labs = []
    for l in labels:
        if l not in labs: labs.append(l)
    el = lambda l: ''.join(c for c in l if c.isalpha())
    nel = 0
    s = f"""&CONTROL
  calculation='scf', prefix='{prefix}', outdir='./tmp', pseudo_dir='{PSEUDO_DIR}',
  verbosity='low', tprnfor=.false., tstress=.false., disk_io='low'
/
&SYSTEM
  ibrav=0, nat={len(struct)}, ntyp={len(labs)}, ecutwfc={ecut}, ecutrho={ecutrho},
  nspin=2, occupations='smearing', smearing='gaussian', degauss=0.01,
  {extra_sys}
"""
    for i, l in enumerate(labs):
        s += f"  starting_magnetization({i+1})={spins.get(l,0.0)},\n"
    s += f"""/
&ELECTRONS
  electron_maxstep={maxstep}, conv_thr={conv}, mixing_beta={beta}, diagonalization='david'
/
ATOMIC_SPECIES
"""
    for l in labs:
        f, m = PP[el(l)]; s += f"{l} {m} {f}\n"
    s += "CELL_PARAMETERS angstrom\n"
    for v in struct.lattice.matrix: s += "  %.8f %.8f %.8f\n" % tuple(v)
    s += "ATOMIC_POSITIONS crystal\n"
    for l, site in zip(labels, struct): s += "%s %.8f %.8f %.8f\n" % ((l,) + tuple(site.frac_coords))
    s += f"K_POINTS automatic\n{kgrid[0]} {kgrid[1]} {kgrid[2]} 0 0 0\n"
    if U:
        s += "HUBBARD {atomic}\n"
        for l in labs:
            if el(l) == 'Fe': s += f"U {l}-3d {U}\n"
    open(fname, 'w').write(s)
