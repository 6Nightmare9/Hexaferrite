import json, time, urllib.request, urllib.parse, datetime, re
import xml.etree.ElementTree as ET
UA={'User-Agent':'hexw-litsearch/0.1'}
OUT='D:/Work/Hexx/lit_search/results/raw_hits3.json'
Q=json.load(open('D:/Work/Hexx/lit_search/scripts/queries3.json'))
def get(url,raw=False,tries=3,wait=4):
    for i in range(tries):
        try:
            r=urllib.request.urlopen(urllib.request.Request(url,headers=UA),timeout=40).read()
            return r if raw else json.loads(r)
        except Exception as e:
            err=e; time.sleep(wait*(i+1))
    print('ERR',url[:110],err,flush=True); return None
now=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
hits={}; log=[]
def add(rec,q,t,api):
    key=rec['doi'] or rec.get('arxiv') and 'arxiv:'+rec['arxiv'] or rec['title'].lower()
    rec.update(topic_q=t,source_api=api)
    if key in hits: hits[key].setdefault('queries',[]).append(q); return
    rec['queries']=[q]; hits[key]=rec
A='{http://www.w3.org/2005/Atom}'
for t,qs in Q.items():
    for q in qs:
        # Europe PMC
        j=get('https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode({'query':q,'format':'json','pageSize':25,'resultType':'core'}))
        n=0
        if j:
            n=j.get('hitCount',0)
            for w in j['resultList']['result']:
                add(dict(doi=(w.get('doi') or '').lower(),title=w.get('title',''),authors=w.get('authorString',''),year=int(w['pubYear']) if w.get('pubYear') else None,venue=w.get('journalTitle',''),type=w.get('pubType',''),is_oa=w.get('isOpenAccess')=='Y',oa_url='',pdf_url='',abstract=re.sub('<[^>]+>','',w.get('abstractText','') or ''),pmcid=w.get('pmcid','')),q,t,'europepmc')
        log.append(('EuropePMC',q,n,now)); time.sleep(0.5)
        # arXiv
        x=get('https://export.arxiv.org/api/query?'+urllib.parse.urlencode({'search_query':'all:'+q if ' ' not in q else 'all:"%s"'%q,'max_results':15}),raw=True)
        n=0
        if x:
            r=ET.fromstring(x)
            n=int(r.find('{http://a9.com/-/spec/opensearch/1.1/}totalResults').text)
            for e in r.findall(A+'entry'):
                aid=e.find(A+'id').text.split('/abs/')[-1]
                d=e.find('{http://arxiv.org/schemas/atom}doi')
                add(dict(doi=(d.text.lower() if d is not None else ''),arxiv=aid,title=' '.join(e.find(A+'title').text.split()),authors='; '.join(a.find(A+'name').text for a in e.findall(A+'author')[:8]),year=int(e.find(A+'published').text[:4]),venue='arXiv',type='preprint',is_oa=True,oa_url='https://arxiv.org/abs/'+aid,pdf_url='https://arxiv.org/pdf/'+aid,abstract=' '.join(e.find(A+'summary').text.split())),q,t,'arxiv')
        log.append(('arXiv',q,n,now)); time.sleep(3.5)
        # Crossref general query
        j=get('https://api.crossref.org/works?'+urllib.parse.urlencode({'query':q,'rows':15,'filter':'type:journal-article','select':'DOI,title,author,issued,container-title,type,abstract'}))
        n=0
        if j:
            n=j['message'].get('total-results',0)
            for w in j['message']['items']:
                add(dict(doi=w['DOI'].lower(),title=(w.get('title') or [''])[0],authors='; '.join((a.get('given','')+' '+a.get('family','')).strip() for a in w.get('author',[])[:8]),year=(w.get('issued',{}).get('date-parts') or [[None]])[0][0],venue=(w.get('container-title') or [''])[0],type=w.get('type'),is_oa=None,oa_url='',pdf_url='',abstract=re.sub('<[^>]+>','',w.get('abstract','') or '')),q,t,'crossref')
        log.append(('Crossref-query',q,n,now)); time.sleep(0.5)
        # DOAJ (W-type related only)
        if t in ('1','2','3','5'):
            j=get('https://doaj.org/api/search/articles/'+urllib.parse.quote('"%s"'%q if ' ' in q else q)+'?pageSize=15')
            n=0
            if j:
                n=j.get('total',0)
                for w in j['results']:
                    b=w['bibjson']; doi=''
                    for i in b.get('identifier',[]):
                        if i['type']=='doi': doi=i['id'].lower()
                    add(dict(doi=doi,title=b.get('title',''),authors='; '.join(a.get('name','') for a in b.get('author',[])[:8]),year=int(b['year']) if b.get('year') else None,venue=b.get('journal',{}).get('title',''),type='journal-article',is_oa=True,oa_url=(b.get('link') or [{}])[0].get('url',''),pdf_url='',abstract=b.get('abstract','') or ''),q,t,'doaj')
            log.append(('DOAJ',q,n,now)); time.sleep(0.5)
        j=get('https://api.semanticscholar.org/graph/v1/paper/search?'+urllib.parse.urlencode({'query':q,'limit':15,'fields':'title,year,externalIds,openAccessPdf,isOpenAccess,venue,authors,abstract'}),tries=2,wait=8)
        n=0
        if j:
            n=j.get('total',0)
            for w in j.get('data',[]):
                ex=w.get('externalIds') or {}
                oap=(w.get('openAccessPdf') or {}).get('url','')
                add(dict(doi=(ex.get('DOI') or '').lower(),arxiv=ex.get('ArXiv',''),title=w.get('title') or '',authors='; '.join(a['name'] for a in w.get('authors',[])[:8]),year=w.get('year'),venue=w.get('venue') or '',type='',is_oa=bool(w.get('isOpenAccess')),oa_url=oap,pdf_url=oap,abstract=w.get('abstract') or ''),q,t,'semanticscholar')
        log.append(('SemanticScholar',q,n,now)); time.sleep(6)
        print(t,q,flush=True)
json.dump(dict(hits=list(hits.values()),log=log),open(OUT,'w',encoding='utf-8'),ensure_ascii=False)
print(len(hits))
