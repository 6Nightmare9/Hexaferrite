import sys,warnings; warnings.filterwarnings('ignore')
sys.path.insert(0, r'D:\Work\Hexx\dft_feasibility\pilots')
from qe_common import write_pw
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
# alpha-Fe2O3 R-3c (hex setting), literature-type parameters, only for resource pilot
s = Structure.from_spacegroup('R-3c', Lattice.hexagonal(5.038, 13.772), ['Fe','O'], [[0,0,0.3553],[0.3059,0,0.25]])
p = SpacegroupAnalyzer(s, symprec=1e-3).get_primitive_standard_structure()
print(len(p), p.composition)
# AFM: Fe along 3-fold axis ordered + - - +
fe = [i for i,x in enumerate(p) if x.specie.symbol=='Fe']
fe.sort(key=lambda i: p[i].coords @ (p.lattice.matrix.sum(0)/abs(p.lattice.matrix.sum(0)).sum()))  # along c(hex)
sp = {fe[0]:'Fe1',fe[1]:'Fe2',fe[2]:'Fe2',fe[3]:'Fe1'}
labels = [sp.get(i,'O') for i in range(len(p))]
p.to(filename='a_hematite/struct.cif')
for tag,ec,er,k in [('ec40_k2',40,200,(2,2,2)),('ec40_k3',40,200,(3,3,3))]:
    write_pw(f'a_hematite/hem_{tag}.in', p, labels, {'Fe1':0.4,'Fe2':-0.4}, 'hem'+tag, ec, er, k, maxstep=8)
