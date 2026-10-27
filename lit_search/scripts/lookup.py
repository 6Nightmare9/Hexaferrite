import json,urllib.request,urllib.parse,sys,time,re
UA={'User-Agent':'hexw-litsearch/0.1'}
def g(u,raw=False):
    try:
        r=urllib.request.urlopen(urllib.request.Request(u,headers=UA),timeout=40).read()
        return r.decode('utf-8','replace') if raw else json.loads(r)
    except Exception as e: return {'err':str(e)} if not raw else 'ERR '+str(e)
for doi in ['10.1103/PhysRevB.109.024414','10.1016/0304-8853(95)00406-8','10.1016/0168-583X(89)90428-X','10.1002/sca.21041']:
    print('==',doi)
    c=g('https://api.crossref.org/works/'+urllib.parse.quote(doi))
    m=c.get('message',{})
    print('crossref title:',(m.get('title') or [''])[0],'| license:',[l['URL'] for l in m.get('license',[])][:3])
    e=g('https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:%22'+doi+'%22&format=json&resultType=core')
    for r in e.get('resultList',{}).get('result',[]): print('EPMC:',r.get('title'),r.get('isOpenAccess'),r.get('pmcid'),r.get('inEPMC'),r.get('inPMC'))
    time.sleep(1)
for q in ['all:SrFe18O27','ti:"W-type SrFe18O27"','all:(first AND principles AND SrFe18O27)','all:Islam AND all:Christensen AND all:hexaferrite','ti:"Zn2-xCoxW"','all:"Newbury" AND all:"Ritchie" AND all:"quantitative"']:
    x=g('https://export.arxiv.org/api/query?search_query='+urllib.parse.quote(q,safe=':()')+'&max_results=5',raw=True)
    print('arXiv',q,[ (t) for t in re.findall(r'<title>(.*?)</title>',x,re.S)[1:]][:5], re.findall(r'<id>http://arxiv.org/abs/(.*?)</id>',x)[:5]); time.sleep(3.5)
