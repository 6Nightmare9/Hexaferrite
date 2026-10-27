"""Read-only check of a cloned HexW installation; never changes saved forecasts."""
from pathlib import Path
import os,sys,json
root=Path(__file__).resolve().parents[1]
os.environ['HEXW_ROOT']=str(root)
manifest=json.loads((root/'audit/GITHUB_EXPORT_MANIFEST.json').read_text(encoding='utf-8'))
missing=[];pointers=[]
for entry in manifest['source_files']:
 p=root/entry.get('export_path',entry['path'])
 if not p.is_file():missing.append(entry['path']);continue
 with p.open('rb') as f:
  if f.read(100).startswith(b'version https://git-lfs.github.com/spec/v1'):pointers.append(entry['path'])
if missing or pointers:
 print('ERROR: missing files:',len(missing),'unhydrated Git LFS pointers:',len(pointers))
 print('Run git lfs pull in the clone. First affected files:',(missing+pointers)[:5])
 raise SystemExit(1)
import hexw
if not Path(hexw.__file__).resolve().is_relative_to(root/'project'):
 print('ERROR: hexw imports from a different checkout:',hexw.__file__)
 raise SystemExit(1)
from hexw.registry.freeze import verify_manifest
try:
 m=verify_manifest(root/'registry/FROZEN_MODEL.json')
 print('Frozen reference verified:',m['model_version'])
except (ValueError,OSError) as e:
 print('NOTICE: the historical freeze is incompatible with this code/runtime:',e)
 print('Old results remain available. Run and freeze a new forecast before validation.')
from hexw.ui.server import Application
app=Application()
try:
 print('Saved forecasts:',len(app.runs()),'saved analyses:',len(app.analyses()))
finally:app.worker.shutdown()
print('Installation check complete. Start START_HEXW.cmd and open http://127.0.0.1:8765/')
