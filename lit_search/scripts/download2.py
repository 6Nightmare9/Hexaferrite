import json, re, time, urllib.request, urllib.parse, hashlib, os, subprocess, shutil, html, http.client, datetime

R = 'D:/Work/Hexx/lit_search/results/'
DL = 'D:/Work/Hexx/lit_search/downloads/'
UA = {'User-Agent': 'hexw-litsearch/0.1'}
rows = json.load(open(R + 'rows_final.json', encoding='utf-8'))


def fetch(url, tries=3):
    err = None
    for i in range(tries):
        try:
            resp = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60)
            return resp.read(40_000_000), resp.headers.get('Content-Type', ''), resp.geturl()
        except http.client.IncompleteRead as e:
            err = e; time.sleep(3)
        except Exception as e:
            err = e; time.sleep(3)
            if 'HTTP Error 403' in str(e) or 'HTTP Error 404' in str(e): break
    raise err


def words(t): return set(re.findall(r'[a-z0-9]{4,}', re.sub(r'<[^>]+>', ' ', t.lower())))


def pdftext(p):
    return subprocess.run(['pdftotext', '-l', '2', p, '-'], capture_output=True, timeout=60).stdout.decode('utf-8', 'replace')


def verify(p, title, doi):
    txt = pdftext(p)
    tw = words(title) - {'type', 'with', 'from', 'their', 'study', 'properties', 'effect'}
    hit = len(tw & words(txt)) / max(1, len(tw))
    if hit >= 0.6: return True, 'title word overlap %.2f' % hit
    if doi and doi.lower() in txt.lower().replace(' ', ''): return True, 'DOI string found in first pages (title overlap %.2f)' % hit
    return False, 'title overlap %.2f' % hit


targets = [r for r in rows if r['rel'] != 3 and not r['downloaded_file'] and r['oa_status'] != 'closed_or_unknown' and 'уже в корпусе' not in r['dl_note']]
for r in targets:
    u = r['u']
    urls = []
    # Europe PMC route (PMC OA subset)
    if u['doi']:
        try:
            d, _, _ = fetch('https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:%22' + urllib.parse.quote(u['doi']) + '%22&format=json&resultType=lite')
            res = json.loads(d)['resultList']['result']
            for x in res:
                if x.get('pmcid') and x.get('isOpenAccess') == 'Y':
                    urls.append('https://europepmc.org/backend/ptpmcrender.fcgi?accid=%s&blobtype=pdf' % x['pmcid'])
        except Exception as e:
            pass
        time.sleep(0.5)
    for k in ('pdf_url', 'oa_url'):
        if u.get(k) and u[k] not in urls: urls.append(u[k])
    note = []
    for url in urls:
        try:
            data, ct, final = fetch(url)
            time.sleep(1.5)
            if data[:5] != b'%PDF-':
                page = data[:2_000_000].decode('utf-8', 'replace')
                m = re.search(r'<meta[^>]+name=["\']citation_pdf_url["\'][^>]+content=["\']([^"\']+)', page) or re.search(r'<meta[^>]+content=["\']([^"\']+)["\'][^>]+name=["\']citation_pdf_url', page)
                if not m: note.append('html-stub:' + url[:50]); continue
                data, ct, final = fetch(urllib.parse.urljoin(final, html.unescape(m.group(1))))
                time.sleep(1.5)
                if data[:5] != b'%PDF-': note.append('meta-not-pdf:' + url[:50]); continue
            d = DL + r['id'] + '/'
            os.makedirs(d, exist_ok=True)
            p = d + r['id'] + '.pdf'
            open(p, 'wb').write(data)
            ok, vn = verify(p, u['title'], u['doi'])
            if not ok:
                shutil.rmtree(d); note.append('mismatch(%s):%s' % (vn, url[:50])); continue
            r['downloaded_file'] = 'downloads/%s/%s.pdf' % (r['id'], r['id'])
            r['sha256'] = hashlib.sha256(data).hexdigest()
            r['oa_url_final'] = final
            r['dl_note'] = 'PDF проверен: ' + vn
            r['retrieved_at'] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
            break
        except Exception as e:
            note.append('err:%s:%s' % (str(e)[:40], url[:50]))
    else:
        pass
    if not r['downloaded_file']: r['dl_note'] = 'не скачано (повтор): ' + '; '.join(note)[:300]
    print(r['id'], u['doi'], r['dl_note'][:140], flush=True)
json.dump(rows, open(R + 'rows_final.json', 'w', encoding='utf-8'), ensure_ascii=False)
print('done')
