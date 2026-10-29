"""Decompress the eQTLGen SMR-formatted BESD trio (.besd/.epi/.esi .gz) in parallel."""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import gzip, shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
D = Path(f"{PROJ_ROOT}/data/raw/eqtl")
def gunzip(f):
    out = D / f.name[:-3].replace("cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense", "eqtlgen_besd")
    if out.exists(): return f"skip {out.name}"
    with gzip.open(f, "rb") as i, open(str(out) + ".tmp", "wb") as o: shutil.copyfileobj(i, o, 1 << 24)
    Path(str(out) + ".tmp").replace(out)
    return f"ok {out.name} {out.stat().st_size / 1e9:.2f} GB"
with ThreadPoolExecutor(3) as ex:
    for m in ex.map(gunzip, sorted(D.glob("*besd-dense.*.gz"))): print(m, flush=True)
