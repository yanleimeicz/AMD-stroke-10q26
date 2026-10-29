import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import openpyxl, csv
jobs = [(f"{PROJ_ROOT}/data/raw/pqtl/UKBPPP_Sun2023_SuppTables.xlsx", "ST9", 4, f"{PROJ_ROOT}/derived/pqtl/ukbppp_ST9_primary.tsv"),
        (f"{PROJ_ROOT}/data/raw/pqtl/UKBPPP_Sun2023_SuppTables.xlsx", "ST3", 2, f"{PROJ_ROOT}/derived/pqtl/ukbppp_ST3_proteins.tsv"),
        (f"{PROJ_ROOT}/data/raw/pqtl/deCODE_Ferkingstad2021_SuppTables.xlsx", "ST02", 2, f"{PROJ_ROOT}/derived/pqtl/decode_ST02_pqtl.tsv")]
import os; os.makedirs(f"{PROJ_ROOT}/derived/pqtl", exist_ok=True)
for f, sh, hdr, out in jobs:
    ws = openpyxl.load_workbook(f, read_only=True)[sh]
    with open(out, "w", newline="", encoding="utf-8") as fo:
        w = csv.writer(fo, delimiter="\t")
        for i, row in enumerate(ws.iter_rows(values_only=True)):
            if i < hdr: continue
            if i == hdr: row = [str(x).replace("\n", " ").strip() if x is not None else f"col{j}" for j, x in enumerate(row)]
            if all(x is None for x in row): continue
            w.writerow(["" if x is None else x for x in row])
    print("wrote", out, flush=True)
