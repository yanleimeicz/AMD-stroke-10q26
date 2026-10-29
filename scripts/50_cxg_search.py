"""Search CELLxGENE Discover for human eye (retina/RPE/choroid) and brain-vascular datasets with H5AD size."""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import json, sys, requests
from pathlib import Path
w = lambda s: sys.stdout.buffer.write((s + "\n").encode("utf-8"))
cache = Path(f"{PROJ_ROOT}/derived/cxg_datasets.json")
if not cache.exists():
    r = requests.get("https://api.cellxgene.cziscience.com/curation/v1/datasets", timeout=600)
    cache.write_text(r.text, encoding="utf-8")
ds = json.loads(cache.read_text(encoding="utf-8"))
w(f"datasets: {len(ds)}")
KEY = ["retina", "pigment", "choroid", "macula", "eye", "fovea", "blood vessel", "vascul", "pericyte", "cerebral", "middle cerebral artery"]
for d in ds:
    org = [o["label"] for o in d.get("organism", [])]
    if "Homo sapiens" not in org: continue
    tis = [t["label"] for t in d.get("tissue", [])]
    title = d.get("title", "")
    txt = (" ".join(tis) + " " + title).lower()
    eye = any(k in txt for k in ["retina", "pigment", "choroid", "macula", "fovea", "eye"])
    vas = ("brain" in txt or "cerebr" in txt or "cortex" in txt) and any(k in txt for k in ["vascul", "vessel", "pericyte", "endothel"])
    if not (eye or vas): continue
    h5 = [a for a in d.get("assets", []) if a.get("filetype") == "H5AD"]
    size = h5[0].get("filesize", 0) / 1e6 if h5 else -1
    w(f"{'EYE' if eye else 'VAS'} | {size:8.0f} MB | cells={d.get('cell_count')} | {title[:90]} | tissues={';'.join(tis)[:80]} | {h5[0]['url'] if h5 else ''}")
