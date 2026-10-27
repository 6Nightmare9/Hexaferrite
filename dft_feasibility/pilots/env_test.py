import subprocess,os,sys,time,re
Q=r"D:\Work\Hexx\tools\qe_ucrt64\ucrt64\bin"
inp=sys.argv[1]; to=float(sys.argv[2])
for cfg in sys.argv[3:]:
    env=dict(os.environ);env["PATH"]=Q+";"+env["PATH"];env["OMP_STACKSIZE"]="256M"
    for k in cfg.split(','):
        a,b=k.split('=');env[a]=b
    t=time.time()
    try:
        r=subprocess.run([Q+r"\pw_big.exe","-i",inp],env=env,capture_output=True,text=True,timeout=to,cwd=os.getcwd())
        o=r.stdout;st='ok'
    except subprocess.TimeoutExpired as e:
        o=(e.stdout or b'').decode(errors='ignore') if isinstance(e.stdout,bytes) else (e.stdout or '');st='TIMEOUT'
    cpu=re.findall(r'total cpu time spent up to now is\s+([\d.]+)',o)
    print(cfg,st,round(time.time()-t,1),'s; cpu-time stamps',cpu[:6],flush=True)
