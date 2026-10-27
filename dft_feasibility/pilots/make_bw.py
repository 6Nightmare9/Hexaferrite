import sys,warnings,itertools; warnings.filterwarnings('ignore')
sys.path.insert(0, r'D:\Work\Hexx\dft_feasibility\pilots')
import numpy as np
from qe_common import write_pw
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

def mindist(s):
    d = s.distance_matrix + np.eye(len(s))*99
    return d.min()

# ---- b) M-type SrFe12O19, literature-type positions (resource pilot only)
best=None
for sr in ([2/3,1/3,.25],[1/3,2/3,.25]):
    s = Structure.from_spacegroup('P6_3/mmc', Lattice.hexagonal(5.8844, 23.0405),
        ['Sr','Fe','Fe','Fe','Fe','Fe','O','O','O','O','O'],
        [sr,[0,0,0],[0,0,.25],[1/3,2/3,.02715],[1/3,2/3,.19240],[.1688,.3376,-.10837],
         [0,0,.15117],[1/3,2/3,-.05477],[.1811,.3622,.25],[.1566,.3132,.0526],[.4975,.995,.1499]])
    print('M-type nat',len(s),s.composition, 'min dist',round(mindist(s),3))
    # Fe-O distances
    fe_o=[min(s.get_distance(i,j) for j,x in enumerate(s) if x.specie.symbol=='O') for i,x in enumerate(s) if x.specie.symbol=='Fe']
    print('  Fe-O nearest min/max',round(min(fe_o),2),round(max(fe_o),2))
    sc=mindist(s)+0
    if best is None or max(fe_o)<best[1]: best=(s,max(fe_o))
s=best[0]
# collinear Gorter: up: 2a,2b,12k ; down: 4f1,4f2
sga=SpacegroupAnalyzer(s,symprec=1e-3); ss=sga.get_symmetrized_structure()
labels=['']*len(s)
for eq,w in zip(ss.equivalent_indices,ss.wyckoff_symbols):
    for i in eq:
        el=s[i].specie.symbol
        labels[i] = el if el!='Fe' else ('Fe2' if w=='4f' else 'Fe1')
print('M-type wyckoff:',list(zip(ss.wyckoff_symbols,[len(e) for e in ss.equivalent_indices])))
s.to(filename='b_Mtype/SrFe12O19_model.cif')
for tag,ec,er,k in [('k111',40,240,(1,1,1)),('k221',40,240,(2,2,1))]:
    write_pw(f'b_Mtype/M_{tag}.in', s, labels, {'Fe1':0.4,'Fe2':-0.4}, 'M'+tag, ec, er, k, maxstep=6)
open('b_Mtype/nat.txt','w').write(str(len(s)))

# ---- c) W-type from Morch 2019 HRPT CIF positions (read-only source), ordered models
a,c=5.909219,32.82927
sites = [  # label, element, xyz, sign
 ('Sr','Sr',[0,0,.25],0), ('7PR','Fe',[1/3,2/3,.74525],+1), ('4OS','Fe',[1/3,2/3,-.07528],+1),
 ('3OSS','Fe',[.5,0,0],+1), ('1TS','Fe',[0,0,.056499],-1), ('2TS','Fe',[1/3,2/3,.094],-1),
 ('6OR','Fe',[1/3,2/3,.20856],-1), ('5ORS','Fe',[.83605,.67212,.151391],+1),
 ('O1','O',[.8223,.6447,.03517],0),('O2','O',[1/3,2/3,.03394],0),('O3','O',[.5123,.0244,.11054],0),
 ('O4','O',[0,0,.11473],0),('O5','O',[.16169,.32346,.18037],0),('O6','O',[1/3,2/3,-.18055],0),('O7','O',[.4841,.9683,.25],0)]
from pymatgen.symmetry.groups import SpaceGroup
sg=SpaceGroup('P6_3/mmc'); lat=Lattice.hexagonal(a,c)
pts=[]  # (label, element, sign, frac)
for lab,el,xyz,sg_ in sites:
    for f in sg.get_orbit(xyz, tol=1e-3):
        pts.append((lab,el,sg_,np.array(f)%1.0))
print('raw orbit total',len(pts))
# resolve split 7PR (half occupancy): keep one of each close pair
keep=[]
pr=[p for p in pts if p[0]=='7PR']
used=set()
for i,p in enumerate(pr):
    if i in used: continue
    for j in range(i+1,len(pr)):
        d=lat.get_distance_and_image(p[3],pr[j][3])[0]
        if j not in used and d<1.0: used.add(j); break
    keep.append(p)
pts=[p for p in pts if p[0]!='7PR']+keep
print('after 7PR resolve',len(pts))
def build(zn_on_1TS):
    el=[];fc=[];labels=[]
    for lab,e,sgn,f in pts:
        if e=='Fe' and zn_on_1TS and lab=='1TS': e2='Zn'; l='Zn'
        elif e=='Fe': e2='Fe'; l='Fe1' if sgn>0 else 'Fe2'
        else: e2=e; l=e
        el.append(e2); fc.append(f); labels.append(l)
    return Structure(lat,el,fc), labels
for name,zn in [('SrFe18O27',False),('SrZn2Fe16O27',True)]:
    s,labels=build(zn)
    d=s.distance_matrix+np.eye(len(s))*99
    print(name,len(s),s.composition,'min dist',round(d.min(),3),'nup',labels.count('Fe1'),'ndn',labels.count('Fe2'))
    s.to(filename=f'c_W92/{name}_model.cif')
    for tag,ec,er,k in [('k111',40,240,(1,1,1)),('k221',40,240,(2,2,1))]:
        write_pw(f'c_W92/W_{name}_{tag}.in', s, labels, {'Fe1':0.4,'Fe2':-0.4}, 'W'+name+tag, ec, er, k, maxstep=6)
