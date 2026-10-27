import json, re, csv, collections

R = 'D:/Work/Hexx/lit_search/results/'
rows = json.load(open(R + 'rows_final.json', encoding='utf-8'))
for r in rows:
    r['alt'] = ''
    u = r['u']
    if not u['doi']:
        r['alt'] = ('arxiv:' + u['arxiv']) if u.get('arxiv') else ('openalex:' + u.get('oa_id', '').split('/')[-1] if u.get('oa_id') else '')
clean = lambda s: re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', '', __import__('html').unescape(s or ''))).strip()

# ---- candidates.csv
cols = ['id', 'doi', 'title', 'authors', 'year', 'venue', 'type', 'relevance', 'topic_tags', 'oa_status', 'oa_url', 'downloaded_file', 'sha256',
        'retrieved_at', 'source_api', 'abstract_available', 'reason', 'alt_id', 'corpus_link', 'download_note']
with open(R + 'candidates.csv', 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(cols)
    for r in rows:
        u = r['u']
        tags = 'composition=%s; property=%s; method=%s' % (','.join(r['comp']) or '-', ','.join(r['prop']) or '-', ','.join(r['meth']) or '-')
        w.writerow([r['id'], u['doi'], clean(u['title']), clean(u['authors']), u['year'], clean(u['venue']), r['role'], r['rel'], tags, r['oa_status'],
                    r.get('oa_url_final', ''), r['downloaded_file'], r['sha256'], r['retrieved_at'], ';'.join(dict.fromkeys(u['apis'])),
                    'yes' if u['abstract'] else 'no', r['reason'], r['alt'], u['corpus_link'], r['dl_note']])

# ---- coverage
PROPS = [('структура/заселённость', r'crystal structure|rietveld|site occup|cation distribution|refinement|lattice parameter|structural|occupan|XRD|diffraction'),
         ('Ms', r'saturation magnetization|\bMs\b|magnetization'),
         ('анизотропия/HA', r'anisotrop|easy (axis|plane)|anisotropy field|\bHA\b'),
         ('Tспин-реориентации', r'spin.reorientation|reorientation|conical|cone'),
         ('ФМР/μ(f)', r'resonance|permeab|gyromagnetic|linewidth|\bFMR\b|absorption|reflection loss'),
         ('мёссбауэр', r'm[öo]ssbauer'),
         ('диэлектрика', r'dielectric|permittivity|conductivity|hopping|impedance'),
         ('DFT', r'first.principles|\bDFT\b|density functional|ab initio')]
ELS = ['Co', 'Zn', 'Mn', 'Ni', 'Mg', 'Cu']
cells = collections.defaultdict(list)
for r in rows:
    u = r['u']
    if r['rel'] == 4 and not re.search(r'W-type|W structure', u['title'] + u['abstract']): continue
    t = clean(u['title']) + ' ' + clean(u['abstract'])
    if not re.search(r'W-type|W type|Fe18O27|Fe16O27|[A-Za-z]2-?W\b|W structure|W-hexa|W hexa', t): continue
    tt = clean(u['title'])
    sr = bool(re.search(r'\bSr|strontium', tt)) or (not re.search(r'\bBa|barium', tt) and bool(re.search(r'strontium|\bSr[A-Z(0-9]', t)))
    ba = bool(re.search(r'\bBa(?![a-z])|barium', tt)) or (not sr and bool(re.search(r'\bBa(?![a-z])|barium', t)))
    systems = (['Sr'] if sr else []) + (['Ba'] if ba else [])
    if not systems: systems = ['(Ca/Pb/другое/не указано)']
    elems = [e for e in ELS if re.search(r'(?<![A-Za-z])' + e + r'(?![a-z])', tt)] or [e for e in ELS if re.search(r'(?<![A-Za-z])' + e + r'(?![a-z]{2})', t)]
    if re.search(r'Fe18O27|Fe2W|Fe2-W', t): elems.append('Fe18O27/Fe2W (без Me)')
    if not elems: elems = ['(Me не указан)']
    basis = 'реферат' if u['abstract'] else 'только заголовок'
    for s in systems:
        for e in elems:
            for pn, pr in PROPS:
                if re.search(pr, t, re.I):
                    cells[(s, e, pn)].append((r['id'], basis))
with open(R + 'coverage_from_abstracts.csv', 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['система', 'состав_Me', 'свойство', 'n_записей', 'id_записей (основание)', 'пометка'])
    for (s, e, p), v in sorted(cells.items()):
        w.writerow([s, e, p, len(v), '; '.join('%s(%s)' % x for x in v), 'по реферату/заголовку, не проверено по тексту'])
# matrix
keys = sorted({(s, e) for (s, e, p) in cells})
with open(R + 'coverage_matrix.csv', 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['система', 'состав_Me'] + [p for p, _ in PROPS] + ['пометка'])
    for s, e in keys:
        w.writerow([s, e] + [len(cells.get((s, e, p), [])) for p, _ in PROPS] + ['число записей; по реферату/заголовку, не проверено по тексту'])

# screened out
sc = json.load(open(R + 'screened_out.json', encoding='utf-8'))
with open(R + 'screened_out.csv', 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f); w.writerow(['key', 'title', 'year', 'why'])
    for s in sc: w.writerow([s['key'], clean(s['title']), s['year'], s['why']])
print(len(rows), collections.Counter(r['rel'] for r in rows))
print('downloaded', sum(1 for r in rows if r['downloaded_file']))
