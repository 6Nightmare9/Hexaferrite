"""External test supervisor: no hexw import in the log-archiving parent."""
from pathlib import Path
import os,subprocess,sys,time,json,tempfile,shutil,stat
from datetime import datetime,timezone
import xml.etree.ElementTree as ET
if os.name=='nt':
    import ctypes
    k=ctypes.WinDLL('kernel32',use_last_error=True)
    k.GetCurrentProcess.restype=ctypes.c_void_p
    k.SetPriorityClass.argtypes=[ctypes.c_void_p,ctypes.c_uint]
    k.SetPriorityClass(k.GetCurrentProcess(),0x4000)
root=Path(__file__).resolve().parents[1];audit=root/'audit'
env=dict(os.environ,HEXW_ROOT=str(root),PYTHONDONTWRITEBYTECODE='1',PYTEST_DISABLE_PLUGIN_AUTOLOAD='1',HEXW_RUN_DFT='0',MPLCONFIGDIR=str(root/'tools/venv/mpl_config'),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2')
selection=sys.argv[1:] or ['tests'];start=time.monotonic()
workers=int(os.environ.get('HEXW_TEST_WORKERS','1'))
if workers not in (1,2):raise ValueError('shared-laptop test worker limit is 2')
if workers==2:
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
temp=Path(tempfile.mkdtemp(prefix='hexx-review-tests-'))
cmd=[sys.executable,'-B','-m','pytest',*selection,'-m','not slow and not dft','-q','-p','no:cacheprovider','--basetemp',str(temp/'pytest'),'--junitxml',str(temp/'results.xml'),'--durations=5']
if workers==2:cmd+=['-p','xdist.plugin','-n','2','--dist','loadscope']
running={'status':'running','pid':os.getpid(),'selection':selection,'log':str(temp/'pytest.log'),'started_utc':datetime.now(timezone.utc).isoformat()}
(audit/'current_test_run.json').write_text(json.dumps(running),encoding='utf-8');print(json.dumps(running),flush=True)
with (temp/'pytest.log').open('w',encoding='utf-8') as log:
    try:
        result=subprocess.run(cmd,cwd=root/'project',env=env,stdout=log,stderr=subprocess.STDOUT,timeout=2100)
        exitcode=result.returncode
    except subprocess.TimeoutExpired:
        exitcode=124
totals={}
if (temp/'results.xml').exists():
    for suite in ET.parse(temp/'results.xml').getroot().iter('testsuite'):
        for key in ('tests','failures','errors','skipped'):
            totals[key]=totals.get(key,0)+int(suite.attrib.get(key,0))
    totals['passed']=totals['tests']-totals['failures']-totals['errors']-totals['skipped']
record={'returncode':exitcode,'totals':totals,'workers':workers,'wall_s':round(time.monotonic()-start,2),'python':sys.executable,'selection':selection,'raw_log':str(temp/'pytest.log'),'raw_xml':str(temp/'results.xml')}
(audit/'last_test_run.json').write_text(json.dumps(record,indent=2),encoding='utf-8');print(json.dumps(record),flush=True)
try:
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    saved=root/'vault/recovery_audit/2026-10-26'/('tests-'+stamp);saved.mkdir(exist_ok=False)
    for name in ('pytest.log','results.xml'):
        source=temp/name
        if source.exists():
            target=saved/name;shutil.copyfile(source,target);target.chmod(stat.S_IREAD)
    record['archive']=str(saved)
except OSError as e:
    record['archive_error']=str(e)
(audit/'last_test_run.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
(audit/'current_test_run.json').write_text(json.dumps({'status':'complete',**record},indent=2),encoding='utf-8')
if exitcode:
    print((temp/'pytest.log').read_text(encoding='utf-8',errors='replace')[-20000:],flush=True)
print(json.dumps(record),flush=True)
sys.exit(exitcode)
