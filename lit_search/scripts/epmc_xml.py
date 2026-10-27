import json,urllib.request,urllib.parse,hashlib,os,re,time,datetime
R='D:/Work/Hexx/lit_search/results/';DL='D:/Work/Hexx/lit_search/downloads/'
UA={'User-Agent':'hexw-litsearch/0.1'}
rows=json.load(open(R+'rows_final.json',encoding='utf-8'))
def get(u):
    return urllib.request.urlopen(urllib.request.Request(u,headers=UA),timeout=60).read()
for r in rows:
    u=r['u']
    if r['downloaded_file'] or r['rel']==3 or r['oa_status']=='closed_or_unknown' or not u['doi'] or 'уже в корпусе' in r['dl_note']: continue
    try:
        res=json.loads(get('https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:%22'+urllib.parse.quote(u['doi'])+'%22&format=json&resultType=lite'))['resultList']['result']
    except Exception as e: print(r['id'],'search err',e); continue
    pm=[x['pmcid'] for x in res if x.get('pmcid') and x.get('isOpenAccess')=='Y']
    if not pm: print(r['id'],u['doi'],'no PMC OA'); continue
    try:
        x=get('https://www.ebi.ac.uk/europepmc/webservices/rest/%s/fullTextXML'%pm[0])
    except Exception as e: print(r['id'],pm[0],'xml err',e); continue
    txt=x[:200000].decode('utf-8','replace')
    t=re.sub(r'<[^>]+>',' ',txt.lower())
    tw=set(re.findall(r'[a-z0-9]{4,}',re.sub(r'<[^>]+>',' ',u['title'].lower())))-{'type','with','from','their','study','properties','effect'}
    hit=len(tw&set(re.findall(r'[a-z0-9]{4,}',t)))/max(1,len(tw))
    if not x.lstrip().startswith(b'<?xml') and b'<article' not in x[:2000] or hit<0.6: print(r['id'],'xml verify fail',hit); continue
    d=DL+r['id']+'/'; os.makedirs(d,exist_ok=True)
    open(d+r['id']+'.xml','wb').write(x)
    r['downloaded_file']='downloads/%s/%s.xml'%(r['id'],r['id']); r['sha256']=hashlib.sha256(x).hexdigest()
    r['oa_url_final']='https://europepmc.org/article/PMC/'+pm[0]; r['dl_note']='Europe PMC fullTextXML (%s) проверен: title overlap %.2f'%(pm[0],hit)
    r['retrieved_at']=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
    print(r['id'],u['doi'],'OK',pm[0]); time.sleep(1)
json.dump(rows,open(R+'rows_final.json','w',encoding='utf-8'),ensure_ascii=False)
