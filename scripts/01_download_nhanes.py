"""Download NHANES 2005-2006 (_D) and 2007-2008 (_E) XPT files needed for the PD x eye-disease analysis."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "raw" / "nhanes"
OUT.mkdir(parents=True, exist_ok=True)

CYCLES = {"D": "2005", "E": "2007"}
# Questionnaire / exam / lab components
COMPONENTS = [
    "DEMO", "RXQ_RX", "VIQ", "OPXRET", "VIX",      # demographics, Rx, vision Q, retinal imaging, vision exam
    "BMX", "BPX", "SMQ", "ALQ", "DIQ", "BPQ", "MCQ", "KIQ_U",
    "GHB", "GLU", "TCHOL", "HDL", "DPQ", "PAQ", "HIQ", "SLQ", "BIOPRO", "CRP",
]
URL = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{year}/DataFiles/{name}_{suf}.xpt"


def fetch(comp_suf):
    comp, suf = comp_suf
    dest = OUT / f"{comp}_{suf}.xpt"
    if dest.exists() and dest.stat().st_size > 0:
        return f"skip {dest.name}"
    r = requests.get(URL.format(year=CYCLES[suf], name=comp, suf=suf), timeout=300)
    if r.status_code != 200 or r.content[:6] != b"HEADER":
        return f"MISSING {comp}_{suf} (HTTP {r.status_code})"
    dest.write_bytes(r.content)
    return f"ok {dest.name} {len(r.content)/1e6:.1f} MB"


jobs = [(c, s) for s in CYCLES for c in COMPONENTS]
with ThreadPoolExecutor(16) as ex:
    for msg in ex.map(fetch, jobs):
        print(msg)
