import json, re
from collections import Counter

R = 'D:/Work/Hexx/lit_search/results/'
U = json.load(open(R + 'merged.json', encoding='utf-8'))
junk = re.compile(r'\.s\d{3}$|/review\d*$|/v\d+/|\.v\d+$|/decision|/response|\.t\d{3}$|\.g\d{3}$|\.f\d{3}$|10\.17760|10\.2139/ssrn')
noise = re.compile(r'BaV6O11|Laves|diffractometer|collimator|nanotube|sample environments|antipolar|exsolution|electric-dipole liquid|photocatal|desulfur|polyaniline|PANI|epoxy|acrylic|ferrofluid|humidity|spin liquid|Realising|Кристаллограф', re.I)
W = re.compile(r'w-type|w type|w-phase|fe18o27|fe16o27|\b(sr|ba|co|zn|mn|ni|mg|cu|fe)2-?w\b|\bw-(type )?ferrite|\bw hexa|\bw-hexa|[a-z]w-type|co2-w|fe2-w|fe2w', re.I)

manual4 = {
 '10.1103/physrevb.71.184433': 'exchange interactions in BaFe12O19 (M-type methodological analogue)',
 '10.1038/s41598-021-81028-7': 'ab initio magnetic properties of SrFe12O19 (M-type DFT analogue)',
 '10.1063/1.4936368': 'ab initio site preference in Ga/In-substituted SrM',
 '10.1063/1.4795821': 'ab initio magnetic anisotropy change in SrCoTi M-type',
 '10.1038/srep36200': 'exchange integrals and Tc calculation for La-substituted BaM',
 '2306.11952': 'first-principles Cr-substituted SrM',
 '10.3390/magnetism5020014': 'ab initio moments, MCA and Tc of hexaferrite',
 '10.1103/physrevmaterials.5.094415': 'anisotropy in La M-type via charge transfer (DFT)',
 '10.1103/physrevapplied.20.024016': 'Ce-substituted SrM magnetic anisotropy (DFT)',
 '10.1063/1.5084762': 'site preference Zn-Sn SrM (ab initio)',
 '1510.00513': 'electronic structure SrFe12O19 first-principles',
 '10.1016/j.jmmm.2013.08.006': 'theory of Zn-Sn pair substitution in SrM',
 '10.1063/1.362313': 'two-magnon FMR linewidth in uniaxial/planar hexagonal ferrite crystals',
 '10.1103/physrevb.79.094427': 'FMR linewidth mechanisms in polycrystalline ferrites',
 '10.1103/physrev.119.1523': 'anisotropic FMR linewidth in ferrites',
 '10.1063/1.365288': 'angle dependence of FMR linewidth, easy-axis/easy-plane hexagonal ferrites',
 '10.1023/a:1009444406700': 'FMR spectra in single-crystal hexagonal ferrites',
 '10.1016/j.jmmm.2006.10.1068': 'FMR equation of hexagonal ferrite in c-plane',
 '10.1063/1.3671793': 'FMR of hexagonal ferrite powders at mm waves',
 '10.1109/tmtt.1963.1125695': 'FMR linewidth, nonoriented polycrystalline hexagonal ferrites with large anisotropy',
 '10.1063/1.340778': 'FMR linewidth contribution in M-phase barium ferrite',
 '10.1063/1.358366': 'transverse susceptibility/anisotropy CoTi-doped BaM crystals',
 '10.1016/s0304-8853(01)00962-3': 'sublattice occupation by Mossbauer in Sr-La-Co M-type',
 '10.1016/j.ceramint.2017.05.264': 'Mossbauer of Zn-substituted SrCo Z-type',
 '10.1063/1.4922867': 'site occupancy in Al-substituted SrM',
 '10.1021/acs.chemmater.4c02205': 'site occupancy in Al-substituted M-type',
 '2312.12001': 'conical magnetic structure in BaM',
 '10.1007/s003390201556': 'neutron powder diffraction (Co,Ti)-BaM',
 '10.1016/0038-1098(79)90926-8': 'anisotropy of M, W, Y hexagonal ferrites with Fe2+ (title names W)',
 '2409.02625': 'FMR/magnetization near collinear-conical transition in hexaferrite',
 '10.1103/physrevb.110.014432': 'magnetic properties of diluted hexaferrites',
}
manual_rel2 = {'10.1016/0038-1098(79)90926-8'}


def tag_comp(t):
    c = []
    if re.search(r'\bSr|strontium', t): c.append('Sr')
    if re.search(r'\bBa(?![a-z])|barium', t): c.append('Ba')
    for el in ['Co', 'Zn', 'Mn', 'Ni', 'Mg', 'Cu']:
        if re.search(r'\b' + el + r'(?![a-z])', t): c.append(el)
    return c


def tag_prop(t):
    P = {'structure/occupancy': r'crystal structure|rietveld|site occup|cation distribution|refinement|lattice parameter|structural',
         'Ms': r'saturation magnetization|magnetization|magnetic properties',
         'anisotropy/HA': r'anisotrop|easy (axis|plane)|cone',
         'spin-reorientation': r'spin reorientation|spin-reorientation',
         'FMR/mu(f)': r'resonance|permeab|gyromagnetic|linewidth|microwave|absorb',
         'Mossbauer': r'm[öo]ssbauer',
         'dielectric': r'dielectric|permittivity|conductivity|hopping|impedance',
         'DFT': r'first.principles|\bdft\b|density functional|ab initio',
         'neutron': r'neutron'}
    return [k for k, v in P.items() if re.search(v, t, re.I)]


def tag_method(t):
    M = {'sol-gel/chemical': r'sol.gel|co-?precip|hydrothermal|auto.?combust', 'ceramic/solid-state': r'solid.state|ceramic|sinter|ball.mill',
         'XRD': r'xrd|diffraction', 'VSM': r'vsm|hysteresis', 'Mossbauer': r'm[öo]ssbauer', 'FMR': r'resonance',
         'DFT': r'first.principles|\bdft\b|density functional|ab initio', 'neutron': r'neutron', 'single-crystal': r'single.crystal'}
    return [k for k, v in M.items() if re.search(v, t, re.I)]


DROP = {'10.4028/www.scientific.net/jnanor.2.91', '10.1139/v11-013', '10.7566/jpsj.82.025003', '2209.03680', '10.1103/physrevx.4.011035'}
BOOST2 = {'10.1007/bf02061493', '10.1016/0304-8853(89)90098-x', '10.1109/20.42281', '10.1016/0304-8853(94)01406-x', '10.1038/170708a0', '10.1016/j.jmmm.2011.11.041'}
OVR3 = {'10.1007/10920503_386', '10.2497/jjspm.17.193', '10.13005/msri/110102', '10.3379/jmsjmag.22.s1_197', '10.21741/9781644902318-1'}
rows = []
screened = []
for u in U:
    d = u['doi'] or ''
    if junk.search(d): continue
    t = u['title'] + ' ' + u['abstract']
    key = d or u.get('arxiv', '')
    m4 = manual4.get(key)
    isW = bool(W.search(t))
    if re.search(r'Z-type|Y-type|U-type|X-type', u['title']) and not re.search(r'w-type', u['title'], re.I):
        isW = False
    if noise.search(u['title']): continue
    if not isW and not m4: continue
    if key in DROP: continue
    title = u['title']
    if isW:
        sr_t = bool(re.search(r'\bSr|strontium', title)) or (bool(re.search(r'strontium|\bSr[A-Z(0-9]', u['abstract'])) and not re.search(r'\bBa(Fe|Co|Zn)|barium', title, re.I))
        prop = tag_prop(t)
        role_sci = bool(re.search(r'anisotrop|resonance|m[öo]ssbauer|neutron|spin reorientation|cation|site|magnetic structure|crystal structure|magnetization', t, re.I))
        coznmn = bool(re.search(r'(Co|Zn|Mn)(?![a-z]{2})', t))
        if sr_t and not (re.search(r'Ba|barium', title) and not re.search(r'Sr|strontium', title)): rel = 1
        elif role_sci or coznmn: rel = 2
        else: rel = 3
        if not u['abstract'] and rel == 2 and not coznmn: rel = 3
        if re.search(r'first.principles|ab initio|dft', t, re.I) and rel == 2: rel = 3
        if re.search(r'dielectric|conductivity', title, re.I) and rel == 2: rel = 3
        reason = 'W-type, ' + ('Sr-based' if sr_t else 'non-Sr') + '; ' + ', '.join(prop or ['general']) + ' (по заголовку/реферату)'
    else:
        rel = 4; reason = 'методологический аналог: ' + m4
    if key in OVR3: rel = 3
    if key in BOOST2: rel = 2
    if key in manual_rel2:
        rel = 2; reason = 'W-структура названа в заголовке: ' + m4
    typ = u.get('type') or ''
    role = 'journal-article'
    if u.get('arxiv') and not d: role = 'preprint'
    if re.search(r'review|overview', title, re.I): role = 'review_navigation'
    if rel == 3 and not u['abstract'] and not re.search(r'dielectric|conductiv|first.principles|m[öo]ssbauer|NMR|resonance|gyromagnetic|permeab|anisotropy|Sr', title, re.I) and key not in BOOST2:
        screened.append(dict(key=key, title=title, year=u['year'], why='rel3, no abstract, synthesis/doping-only title'))
        continue
    rows.append(dict(u=u, rel=rel, reason=reason, role=role, comp=tag_comp(t), prop=tag_prop(t), meth=tag_method(t), key=key))

json.dump(screened, open(R + 'screened_out.json', 'w', encoding='utf-8'), ensure_ascii=False)
print(len(screened), 'screened;', len(rows), Counter(r['rel'] for r in rows))
json.dump(rows, open(R + 'rows_prelim.json', 'w', encoding='utf-8'), ensure_ascii=False)
for r in sorted(rows, key=lambda r: (r['rel'], r['u']['year'] or 0)):
    u = r['u']
    print(r['rel'], u['year'], u['doi'] or u.get('arxiv'), '|', u['title'][:85].replace('\n', ' '), '|', 'A' if u['abstract'] else '-', u.get('is_oa'))
