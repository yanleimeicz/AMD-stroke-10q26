import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import requests, sys, time, os
H={'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126.0 Safari/537.36','Accept':'text/html,application/xhtml+xml,application/json;q=0.9,*/*;q=0.8'}
out=f'{PROJ_ROOT}/derived/eye_eqtl/web'
os.makedirs(out,exist_ok=True)
for line in open(sys.argv[1],encoding='utf-8'):
    line=line.strip()
    if not line or line.startswith('#'): continue
    name,url=line.split('\t')
    try:
        s=requests.Session(); s.trust_env=False
        r=s.get(url,headers=H,timeout=120,allow_redirects=True)
        open(os.path.join(out,name),'wb').write(r.content)
        print(name,r.status_code,len(r.content),r.headers.get('content-type'),r.url)
    except Exception as e:
        print(name,'ERR',e)
    time.sleep(20)
