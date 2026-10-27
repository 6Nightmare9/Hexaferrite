import json, re, time, urllib.request, urllib.parse, hashlib, os, subprocess, datetime, shutil, html

R = 'D:/Work/Hexx/lit_search/results/'
DL = 'D:/Work/Hexx/lit_search/downloads/'
UA = {'User-Agent': 'hexw-litsearch/0.1'}
rows = json.load(open(R + 'rows_stage2.json', encoding='utf-8'))

# drop duplicate arXiv listing without DOI
rows = [r for r in rows if not (not r['u']['doi'] and r['u'].get('oa_url', '').endswith('1407.3605'))]
CORP_FULL = {'L01', 'L02', 'L03', 'L05', 'L10', 'L11', 'L13'}
for i, r in enumerate(rows, 1):
    r['id'] = 'C%03d' % i


def fetch(url):
    req = urllib.request.Request(url, headers=UA)
    resp = urllib.request.urlopen(req, timeout=60)
    return resp.read(40_000_000), resp.headers.get('Content-Type', ''), resp.geturl()


def words(t):
    return set(w for w in re.findall(r'[a-z0-9]{4,}', re.sub(r'<[^>]+>', ' ', t.lower())))


def verify(path, title):
    try:
        txt = subprocess.run(['pdftotext', '-l', '2', path, '-'], capture_output=True, timeout=60).stdout.decode('utf-8', 'replace')
    except Exception as e:
        return False, 'pdftotext failed'
    tw = words(title) - {'type', 'with', 'from', 'their', 'study', 'properties', 'effect', 'the'}
    if not tw: return None, 'no title words'
    xw = words(txt)
    hit = len(tw & xw) / len(tw)
    return hit >= 0.6, 'title word overlap %.2f' % hit


now = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
for r in rows:
    u = r['u']
    r['retrieved_at'] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
    r['oa_status'] = 'closed_or_unknown'
    r['downloaded_file'] = ''; r['sha256'] = ''; r['dl_note'] = ''
    cc = [l for l in u['cr_license'] if 'creativecommons' in l]
    urls = []
    if u.get('arxiv'): urls.append('https://arxiv.org/pdf/' + re.sub('v\\d+$', '', u['arxiv']))
    elif u['doi'].startswith('10.48550/arxiv.'): urls.append('https://arxiv.org/pdf/' + u['doi'].split('arxiv.')[-1])
    for k in ('pdf_url', 'oa_url'):
        if u.get(k) and u[k] not in urls: urls.append(u[k])
    if u.get('arxiv') or u['doi'].startswith('10.48550'):
        r['oa_status'] = 'green_arxiv'
    elif cc:
        r['oa_status'] = 'gold_cc_license (' + cc[0].split('creativecommons.org/')[-1].strip('/') + ')'
    elif u.get('is_oa'):
        r['oa_status'] = 'oa_per_source_metadata (%s)' % u['apis'][0]
    if r['oa_status'] == 'closed_or_unknown': r['oa_url_final'] = ''; continue
    r['oa_url_final'] = (urls[0] if urls else (u['doi'] and 'https://doi.org/' + u['doi']))
    if r['rel'] == 3: r['dl_note'] = 'не качали (rel 3)'; continue
    if u['corpus_link'] in CORP_FULL: r['dl_note'] = 'уже в корпусе (%s), не качали' % u['corpus_link']; continue
    if not urls: r['dl_note'] = 'нет прямой ссылки'; continue
    # try download
    tried = []
    done = False
    for url in urls[:2]:
        cands = [url]
        for c in cands:
            try:
                data, ct, final = fetch(c)
                time.sleep(1.5)
                if data[:5] == b'%PDF-':
                    pdf = data
                else:
                    page = data[:2_000_000].decode('utf-8', 'replace')
                    m = re.search(r'<meta[^>]+name=["\']citation_pdf_url["\'][^>]+content=["\']([^"\']+)', page) or re.search(r'<meta[^>]+content=["\']([^"\']+)["\'][^>]+name=["\']citation_pdf_url', page)
                    if not m:
                        tried.append('html-stub:' + c[:60]); continue
                    pu = urllib.parse.urljoin(final, html.unescape(m.group(1)))
                    data2, ct2, final2 = fetch(pu)
                    time.sleep(1.5)
                    if data2[:5] != b'%PDF-':
                        tried.append('html-stub(pdf-meta not pdf):' + pu[:60]); continue
                    pdf = data2; final = final2
                d = DL + r['id'] + '/'
                os.makedirs(d, exist_ok=True)
                p = d + r['id'] + '.pdf'
                open(p, 'wb').write(pdf)
                ok, note = verify(p, u['title'])
                if ok is False:
                    shutil.rmtree(d)
                    tried.append('pdf-title-mismatch(%s):%s' % (note, final[:60])); continue
                r['downloaded_file'] = 'downloads/%s/%s.pdf' % (r['id'], r['id'])
                r['sha256'] = hashlib.sha256(pdf).hexdigest()
                r['oa_url_final'] = final
                r['dl_note'] = 'PDF проверен: ' + note
                done = True; break
            except Exception as e:
                tried.append('err:%s:%s' % (str(e)[:40], c[:50]))
        if done: break
    if not done: r['dl_note'] = 'не скачано: ' + '; '.join(tried)[:300]
    print(r['id'], u['doi'] or u.get('arxiv'), r['dl_note'][:100], flush=True)
json.dump(rows, open(R + 'rows_final.json', 'w', encoding='utf-8'), ensure_ascii=False)
json.dump(now, open(R + 'retrieved_at.json', 'w'))
print('done')
