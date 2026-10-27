import subprocess,os,sys
Q=r"D:\Work\Hexx\tools\qe_ucrt64\ucrt64\bin"
env=dict(os.environ);env["PATH"]=Q+";"+env["PATH"];env["OMP_NUM_THREADS"]=sys.argv[1]
for k in sys.argv[2:]:
    a,b=k.split('=');env[a]=b
r=subprocess.run([Q+r"\pw_big.exe","-i","hem_ec40_k2.in"],env=env,capture_output=True,text=True,timeout=300)
print(r.returncode);print(r.stdout[-1500:]);print(r.stderr[-800:])
