import re,sys,json,os
for f in sys.argv[1:]:
    o=open(f,errors='ignore').read()
    g=lambda p:(re.findall(p,o) or [None])
    cpu=[float(x) for x in re.findall(r'total cpu time spent up to now is\s+([\d.]+)',o)]
    it=[b-a for a,b in zip(cpu[:-1],cpu[1:])]
    nel=g(r'number of electrons\s+=\s+([\d.]+)')[0]; nk=g(r'number of k points=\s+(\d+)')[0]
    nb=g(r'number of Kohn-Sham states=\s+(\d+)')[0]; ram=g(r'Estimated max dynamical RAM per process >\s+([\d.]+) MB')[0]
    fft=g(r'Dense  grid:\s+\d+ G-vectors\s+FFT dimensions: \(\s*(\d+),\s*(\d+),\s*(\d+)\)')[0]
    npw=g(r'Sum\s+\d+\s+\d+\s+\d+\s+(\d+)')[0]
    mag=g(r'total magnetization\s+=\s+([-\d.]+)')
    print(os.path.basename(f),'nel',nel,'nk',nk,'nbnd',nb,'FFT',fft,'ramEst_MB',ram,'setup_s',cpu[0] if cpu else None,
      'iter_s',[round(x,1) for x in it], 'mean_iter_s(2+)', round(sum(it[1:])/max(1,len(it)-1),1) if len(it)>1 else None, 'last_mag',mag[-1])
