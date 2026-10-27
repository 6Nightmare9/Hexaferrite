import json, time, urllib.request, urllib.parse, sys, datetime
UA={'User-Agent':'hexw-litsearch/0.1'}
OUT='D:/Work/Hexx/lit_search/results/raw_hits.json'
Q={
1:["SrCo2W hexaferrite","SrZn2W hexaferrite","SrMn2W hexaferrite","SrFe18O27","Sr W-type hexaferrite","Sr(Co,Zn)2Fe16O27","SrMe2Fe16O27 W-type","strontium W-type ferrite cobalt zinc substitution"],
2:["W-type hexaferrite Ni2W","W-type hexaferrite Mg2W","W-type hexaferrite Cu2W","Fe2W hexaferrite","W-type ferrite BaFe2+2Fe16O27 substitution"],
3:["BaCo2W hexaferrite","BaZn2W hexaferrite","BaMn2W hexaferrite","BaFe18O27","Co-Zn W-type hexaferrite","spin reorientation W-type hexaferrite","cone anisotropy W-type hexaferrite","ferromagnetic resonance W-type hexaferrite","microwave permeability W-type hexaferrite natural resonance","Mossbauer W-type hexaferrite site occupancy","cation distribution W-type hexaferrite","neutron diffraction W-type hexaferrite"],
4:["first-principles SrFe18O27","first-principles BaFe18O27 W-type","exchange constants hexaferrite DFT+U","magnetocrystalline anisotropy hexaferrite spin-orbit first principles","SrFe12O19 density functional theory magnetic anisotropy","M-type hexaferrite exchange interactions first-principles"],
5:["W-type hexaferrite dielectric permittivity","W-type hexaferrite electrical conductivity Fe2+ Fe3+ hopping","W-type hexaferrite dielectric properties"],
6:["Mossbauer recoil-free fraction site populations ferrites","Smit-Beljers ferromagnetic resonance hexagonal anisotropy","polycrystalline ferromagnetic resonance linewidth Schlomann porosity","ferromagnetic resonance hexagonal ferrite single crystal anisotropy field"],
}
def get(url):
    for i in range(3):
        try:
            return json.load(urllib.request.urlopen(urllib.request.Request(url,headers=UA),timeout=40))
        except Exception as e:
            err=e; time.sleep(3)
    print('ERR',url[:100],err); return None
def abs_from_inv(d):
    if not d: return ''
    w=[]
    for k,v in d.items():
        for p in v: w.append((p,k))
    return ' '.join(k for p,k in sorted(w))
hits={}; log=[]
now=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
def add(rec,q,t,api):
    key=(rec['doi'] or rec['title'].lower())
    rec.update(topic_q=t,query=q,source_api=api)
    if key in hits:
        hits[key]['queries'].append(q); return
    rec['queries']=[q]; hits[key]=rec
for t,qs in Q.items():
    for q in qs:
        # OpenAlex
        u='https://api.openalex.org/works?'+urllib.parse.urlencode({'search':q,'per-page':25,'select':'id,doi,title,authorships,publication_year,primary_location,type,open_access,abstract_inverted_index,best_oa_location'})
        j=get(u); n=0
        if j:
            n=j['meta']['count']
            for w in j['results']:
                doi=(w.get('doi') or '').replace('https://doi.org/','').lower()
                loc=(w.get('primary_location') or {}).get('source') or {}
                oa=w.get('open_access') or {}
                bo=w.get('best_oa_location') or {}
                add(dict(doi=doi,title=w.get('title') or '',authors='; '.join(a['author']['display_name'] for a in w['authorships'][:8])+(' et al.' if len(w['authorships'])>8 else ''),year=w.get('publication_year'),venue=loc.get('display_name') or '',type=w.get('type'),is_oa=oa.get('is_oa'),oa_url=oa.get('oa_url') or '',pdf_url=bo.get('pdf_url') or '',abstract=abs_from_inv(w.get('abstract_inverted_index')),oa_id=w['id']),q,t,'openalex')
        log.append(('OpenAlex',q,n,now)); time.sleep(0.5)
        # Crossref
        u='https://api.crossref.org/works?'+urllib.parse.urlencode({'query.bibliographic':q,'rows':15,'select':'DOI,title,author,issued,container-title,type,abstract'})
        j=get(u); n=0
        if j:
            n=j['message'].get('total-results',0)
            for w in j['message']['items']:
                ti=(w.get('title') or [''])[0]
                add(dict(doi=w['DOI'].lower(),title=ti,authors='; '.join((a.get('given','')+' '+a.get('family','')).strip() for a in w.get('author',[])[:8]),year=(w.get('issued',{}).get('date-parts') or [[None]])[0][0],venue=(w.get('container-title') or [''])[0],type=w.get('type'),is_oa=None,oa_url='',pdf_url='',abstract=w.get('abstract','') or ''),q,t,'crossref')
        log.append(('Crossref',q,n,now)); time.sleep(0.5)
json.dump(dict(hits=list(hits.values()),log=log),open(OUT,'w',encoding='utf-8'),ensure_ascii=False)
print(len(hits))
