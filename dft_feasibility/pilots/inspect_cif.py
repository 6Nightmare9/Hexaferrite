import sys,warnings
warnings.filterwarnings('ignore')
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
p=sys.argv[1]
s=Structure.from_file(p)
print(p.split('/')[-1], s.composition.reduced_formula, 'natoms',len(s))
print(s.lattice.abc, s.lattice.angles, 'vol',round(s.volume,1))
sga=SpacegroupAnalyzer(s,symprec=1e-2)
print(sga.get_space_group_symbol(), 'prim natoms', len(sga.find_primitive()))
print('occupancy-disordered:', not s.is_ordered)
ss=sga.get_symmetrized_structure()
for eq,lab in zip(ss.equivalent_sites, ss.wyckoff_symbols):
    print(eq[0].species_string, lab, len(eq))
