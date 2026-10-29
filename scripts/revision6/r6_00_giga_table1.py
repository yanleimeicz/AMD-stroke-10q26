import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import pandas as pd
P=PROJ_ROOT
x=pd.read_excel(f"{P}/derived/revision6/giga_suppl_tables.xlsx", sheet_name="Suppl.Table1", header=None)
x.to_csv(f"{P}/derived/revision6/giga_suppl_table1.tsv", sep="\t", index=False, header=False)
t=pd.read_excel(f"{P}/derived/revision6/giga_suppl_tables.xlsx", sheet_name="Table.of.content", header=None)
t.to_csv(f"{P}/derived/revision6/giga_suppl_toc.tsv", sep="\t", index=False, header=False)
print(x.shape)
