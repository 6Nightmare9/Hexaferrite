import urllib.request,urllib.parse,json,time,re
UA={'User-Agent':'hexw-litsearch/0.1'}
def g(u,raw=False):
    for i in range(4):
        try:
            d=urllib.request.urlopen(urllib.request.Request(u,headers=UA),timeout=40).read()
            return d.decode('utf-8','replace') if raw else json.loads(d)
        except Exception as e:
            err=str(e); time.sleep(10)
    return {'err':err}
for doi in ['10.1103/PhysRevB.109.024414','10.1016/0304-8853(95)00406-8','10.1016/0168-583X(89)90428-X','10.1002/sca.21041']:
    r=g('https://api.semanticscholar.org/graph/v1/paper/DOI:'+urllib.parse.quote(doi)+'?fields=title,isOpenAccess,openAccessPdf,externalIds')
    print(doi,'S2:',r.get('isOpenAccess'),r.get('openAccessPdf'),(r.get('externalIds') or {}).get('ArXiv'),r.get('err'))
    time.sleep(8)
for q in ['au:Islam_R AND au:Christensen_M AND all:hexaferrite','all:W-type AND all:SrFe18O27','au:Christensen_M AND all:hexaferrite','all:"magnetic properties of W-type"']:
    x=g('https://export.arxiv.org/api/query?search_query='+urllib.parse.quote(q,safe=':()')+'&max_results=5',raw=True)
    print(q,re.findall(r'<title>(.*?)</title>',x,re.S)[1:4] if isinstance(x,str) else x); time.sleep(4)
