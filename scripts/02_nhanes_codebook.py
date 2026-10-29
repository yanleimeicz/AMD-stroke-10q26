"""Download NHANES doc pages and parse variable name / SAS label / question text / value codes."""
from pathlib import Path
import re, requests
import pandas as pd
from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "data" / "raw" / "nhanes" / "docs"
DOC.mkdir(parents=True, exist_ok=True)
OUT = ROOT / "derived" / "nhanes_codebook.tsv"
CYCLES = {"D": "2005", "E": "2007"}
COMPS = ["VIQ", "VIX", "OPXRET", "MCQ", "DIQ", "BPQ", "SMQ", "ALQ", "PAQ", "DPQ", "KIQ_U", "HIQ", "DEMO", "SLQ"]

rows = []
for suf, year in CYCLES.items():
    for comp in COMPS:
        f = DOC / f"{comp}_{suf}.htm"
        if not f.exists():
            r = requests.get(f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{year}/DataFiles/{comp}_{suf}.htm", timeout=120)
            f.write_bytes(r.content)
        soup = BeautifulSoup(f.read_bytes(), "html.parser")
        for blk in soup.select("div.pagebreak"):
            info = {dt.get_text(strip=True).rstrip(":"): dd.get_text(" ", strip=True)
                    for dt, dd in zip(blk.select("dt"), blk.select("dd"))}
            tbl = blk.select_one("table")
            codes = ""
            if tbl:
                codes = " | ".join(
                    f"{tds[0].get_text(strip=True)}={tds[1].get_text(strip=True)} (n={tds[2].get_text(strip=True)})"
                    for tr in tbl.select("tbody tr") if len(tds := tr.select("td")) >= 3)
            rows.append(dict(cycle=suf, file=comp, var=info.get("Variable Name", ""),
                             label=info.get("SAS Label", ""), text=info.get("English Text", ""),
                             target=info.get("Target", ""), codes=codes))
pd.DataFrame(rows).to_csv(OUT, sep="\t", index=False)
print("rows:", len(rows))
