import json,re
d1=json.load(open('D:/Work/Hexx/lit_search/results/raw_hits.json',encoding='utf-8'))
d2=json.load(open('D:/Work/Hexx/lit_search/results/raw_hits2.json',encoding='utf-8'));d3=json.load(open('D:/Work/Hexx/lit_search/results/raw_hits3.json',encoding='utf-8'));d4=json.load(open('D:/Work/Hexx/lit_search/results/raw_hits4.json',encoding='utf-8'))
M={}
def norm(t): return re.sub(r'[^a-z0-9]','',t.lower())
for r in d1['hits']+d2['hits']+d3['hits']+d4['hits']:
    r['abstract']=re.sub('<[^>]+>','',r.get('abstract') or '').strip()
    k=r['doi'] or ('arxiv:'+r['arxiv'] if r.get('arxiv') else norm(r['title']))
    k2=norm(r['title'])
    kk=M.get(('d',k)) or M.get(('t',k2))
    if kk is None:
        r['apis']=[r['source_api']]; r['queries']=list(r.get('queries',[])); kk=r
        M[('d',k)]=kk; M[('t',k2)]=kk
    else:
        for f in ('abstract','oa_url','pdf_url','venue','authors','pmcid','arxiv','oa_id'):
            if not kk.get(f) and r.get(f): kk[f]=r[f]
        if not kk['doi'] and r['doi']: kk['doi']=r['doi']
        if r.get('is_oa') and not kk.get('is_oa'): kk['is_oa']=True
        kk['apis'].append(r['source_api']); kk['queries']+=r.get('queries',[])
U={id(v):v for v in M.values()}.values()
U=list(U)
kw=re.compile(r'hexaferrit|hexagonal ferrit|w-type|w type|w-phase|magnetoplumbite|fe18o27|\bco2w\b|\bzn2w\b|\bmn2w\b|m-type|ferrite|ferrit',re.I)
U=[u for u in U if kw.search(u['title']+' '+u['abstract'])]
print(len(U))
json.dump(U,open('D:/Work/Hexx/lit_search/results/merged.json','w',encoding='utf-8'),ensure_ascii=False)
for i,u in enumerate([]):
    pass#print(i,u['year'],u['doi'] or u.get('arxiv'),'|',u['title'][:95],'|',len(u['abstract']))
