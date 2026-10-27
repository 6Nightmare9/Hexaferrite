import urllib.request,urllib.parse,re,json,time,datetime
import xml.etree.ElementTree as ET
A='{http://www.w3.org/2005/Atom}'
qs=['all:hexaferrite','all:"W-type" AND all:ferrite','all:"W-type" AND all:hexaferrite','all:hexagonal AND all:ferrite AND all:anisotropy','all:SrFe12O19','all:BaFe12O19','all:"M-type" AND all:hexaferrite AND all:DFT','all:Mossbauer AND all:hexaferrite','all:"ferromagnetic resonance" AND all:hexaferrite','all:"strontium hexaferrite"','all:"Fe18O27"','all:"first-principles" AND all:hexaferrite','all:"recoil-free" AND all:Mossbauer AND all:ferrite']
out={};log=[]
now=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
for q in qs:
    u='https://export.arxiv.org/api/query?'+urllib.parse.urlencode({'search_query':q,'max_results':120})
    try: x=urllib.request.urlopen(urllib.request.Request(u,headers={'User-Agent':'hexw-litsearch/0.1'}),timeout=60).read()
    except Exception as e: print('ERR',q,e); log.append(['arXiv',q,-1,now]); time.sleep(5); continue
    r=ET.fromstring(x); n=int(r.find('{http://a9.com/-/spec/opensearch/1.1/}totalResults').text)
    log.append(['arXiv',q,n,now])
    for e in r.findall(A+'entry'):
        aid=re.sub(r'v\d+$','',e.find(A+'id').text.split('/abs/')[-1])
        d=e.find('{http://arxiv.org/schemas/atom}doi')
        out[aid]=dict(doi=(d.text.lower() if d is not None else ''),arxiv=aid,title=' '.join(e.find(A+'title').text.split()),authors='; '.join(a.find(A+'name').text for a in e.findall(A+'author')[:8]),year=int(e.find(A+'published').text[:4]),venue='arXiv',type='preprint',is_oa=True,oa_url='https://arxiv.org/abs/'+aid,pdf_url='https://arxiv.org/pdf/'+aid,abstract=' '.join(e.find(A+'summary').text.split()),topic_q='arxiv_sweep',source_api='arxiv',queries=[q])
    time.sleep(4)
json.dump(dict(hits=list(out.values()),log=log),open('D:/Work/Hexx/lit_search/results/raw_hits4.json','w',encoding='utf-8'),ensure_ascii=False)
print(len(out)); 
for l in log: print(l)
