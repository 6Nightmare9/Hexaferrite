import sys,re
src,dst=sys.argv[1],sys.argv[2]
t=open(src,encoding='utf8',errors='ignore').read()
i=t.index('data_SrZn2Fe16O27') if 'data_SrZn2Fe16O27' in t else t.index('data_')
j=t.find('\ndata_',i+5)
open(dst,'w').write(t[i:j] if j>0 else t[i:])
