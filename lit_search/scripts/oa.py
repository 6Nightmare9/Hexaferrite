import json, re, time, urllib.request, urllib.parse, hashlib, os, subprocess, datetime

R = 'D:/Work/Hexx/lit_search/results/'
DL = 'D:/Work/Hexx/lit_search/downloads/'
UA = {'User-Agent': 'hexw-litsearch/0.1'}
rows = json.load(open(R + 'rows_prelim.json', encoding='utf-8'))


def norm(t): return re.sub(r'[^a-z0-9]', '', re.sub(r'<[^>]+>', '', t.lower()))


# dedupe by DOI / title (prefer the record with DOI and abstract)
best = {}
for r in rows:
    u = r['u']
    k = norm(u['title'])
    if k in best:
        o = best[k]
        if (u['doi'] and not o['u']['doi']) or (len(u['abstract']) > len(o['u']['abstract']) and u['doi']):
            best[k] = r
    else:
        best[k] = r
rows = list(best.values())
# extra: arXiv version of same article as DOI 10.48550 -> keep one
seen = {}
out = []
for r in rows:
    u = r['u']
    t = norm(u['title'])[:60]
    if t in seen and (u['doi'].startswith('10.48550') or seen[t]['u']['doi'].startswith('10.48550')):
        if u['doi'].startswith('10.48550'): continue
        out.remove(seen[t])
    seen[t] = r
    out.append(r)
rows = sorted(out, key=lambda r: (r['rel'], -(r['u']['year'] or 0)))


def get(url, raw=False, tries=2):
    for i in range(tries):
        try:
            resp = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60)
            data = resp.read()
            return (data, resp.headers.get('Content-Type', ''), resp.geturl()) if raw else json.loads(data)
        except Exception as e:
            err = e; time.sleep(2)
    return None


# known-corpus links
CORP = {'10.1107/s2052252519003130': 'L01', '10.1107/s1600576723002133': 'L02', '10.2320/matertrans.mt-y2022010': 'L03',
        '10.1016/j.actamat.2025.120828': 'L05', '10.1103/physrevb.109.024414': 'L04', '10.1016/0304-8853(95)00406-8': 'L06',
        '10.1016/0168-583x(89)90428-x': 'L07', '10.1002/sca.21041': 'L08', '10.1088/1674-1056/ad39d5': 'L11',
        '10.1107/s1600576723011032': 'L10', '10.3390/magnetochemistry9010030': 'L13'}
CORP_FULL = {'L01', 'L02', 'L03', 'L05', 'L11', 'L10', 'L13'}
# add L07, L08 (not found by topic searches): metadata from Crossref
have = {r['u']['doi'] for r in rows}
for doi, rel, why in [('10.1016/0168-583x/89)90428-x', 4, ''), ]:
    pass
for doi, rel, why in [('10.1016/0168-583X(89)90428-X', 4, 'методология: точные заселённости позиций из мёссбауэровских спектров (L07 корпуса)'),
                      ('10.1002/sca.21041', 4, 'методология: количественность SEM/EDS (L08 корпуса)')]:
    if doi.lower() in have: continue
    m = get('https://api.crossref.org/works/' + urllib.parse.quote(doi))
    if not m: continue
    m = m['message']
    ti = re.sub(r'\s+', ' ', re.sub(r'<[^>]+>', '', (m.get('title') or [''])[0])).strip()
    rows.append(dict(u=dict(doi=doi.lower(), title=ti, authors='; '.join((a.get('given', '') + ' ' + a.get('family', '')).strip() for a in m.get('author', [])[:8]),
                            year=(m.get('issued', {}).get('date-parts') or [[None]])[0][0], venue=(m.get('container-title') or [''])[0], type='journal-article',
                            is_oa=None, oa_url='', pdf_url='', abstract='', apis=['crossref'], queries=['L07/L08 DOI lookup'], source_api='crossref'),
                     rel=rel, reason=why, role='journal-article', comp=[], prop=[], meth=['Mossbauer'] if 'L07' in why else [], key=doi.lower()))
    time.sleep(0.5)

# Crossref license info for all with DOI
now = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
for r in rows:
    u = r['u']
    u['corpus_link'] = CORP.get(u['doi'], '')
    u['cr_license'] = []
    u['cr_links'] = []
    if u['doi'] and not u['doi'].startswith('10.48550'):
        j = get('https://api.crossref.org/works/' + urllib.parse.quote(u['doi']))
        if j:
            m = j['message']
            u['cr_license'] = [l['URL'] for l in m.get('license', [])]
            u['cr_links'] = [l['URL'] for l in m.get('link', [])]
            if not u.get('abstract_clean') and m.get('abstract') and not u['abstract']:
                u['abstract'] = re.sub('<[^>]+>', '', m['abstract']).strip()
        time.sleep(0.3)

json.dump(rows, open(R + 'rows_stage2.json', 'w', encoding='utf-8'), ensure_ascii=False)
print(len(rows))
