"""Druggable-genome annotation for all genes tested in SMR (DGIdb GraphQL gene categories + drug interactions).
Finan et al. 2017 druggable genome = DGIdb category "DRUGGABLE GENOME" with source "HingoraniCasas".

Output -> results/smr_heidi/dgidb_gene_categories.tsv
"""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import requests, pandas as pd

OUT = Path(f"{PROJ_ROOT}/results/smr_heidi"); OUT.mkdir(parents=True, exist_ok=True)
genes = sorted(set(pd.read_csv(f"{PROJ_ROOT}/results/smr/smr_all.tsv", sep="\t", usecols=["GeneSymbol"]).GeneSymbol.dropna()))
extra = Path(f"{PROJ_ROOT}/results/smr_heidi/smr_genes_extra.txt")       # probes seen only in the BESD run (optional)
if extra.exists(): genes = sorted(set(genes) | {g.strip() for g in open(extra) if g.strip()})
Q = """query($names:[String!]) { genes(names:$names) { nodes { name
  geneCategoriesWithSources { name sourceNames } interactions { drug { name approved } } } } }"""


def batch(names):
    for attempt in range(6):
        try:
            r = requests.post("https://dgidb.org/api/graphql", json={"query": Q, "variables": {"names": names}}, timeout=180)
            rows = []
            for n in r.json()["data"]["genes"]["nodes"]:
                cats = {c["name"]: c["sourceNames"] for c in n["geneCategoriesWithSources"]}
                drugs = sorted({i["drug"]["name"] for i in n["interactions"] if i["drug"]["approved"]})
                rows.append(dict(GeneSymbol=n["name"],
                                 druggable_finan="HingoraniCasas" in cats.get("DRUGGABLE GENOME", []),
                                 druggable_any="DRUGGABLE GENOME" in cats,
                                 clinically_actionable="CLINICALLY ACTIONABLE" in cats,
                                 n_interactions=len(n["interactions"]), n_approved_drugs=len(drugs),
                                 approved_drugs=";".join(drugs[:15]), categories=";".join(sorted(cats))))
            return rows
        except Exception as e:
            print("retry", names[0], e, flush=True); time.sleep(10 * (attempt + 1))
    print("FAILED batch", names[0], flush=True); return []


chunks = [genes[i:i + 50] for i in range(0, len(genes), 50)]
with ThreadPoolExecutor(8) as ex:
    res = [r for rows in ex.map(batch, chunks) for r in rows]
d = pd.DataFrame(res).drop_duplicates("GeneSymbol")
d.to_csv(OUT / "dgidb_gene_categories.tsv", sep="\t", index=False)
print(f"queried {len(genes)} genes; found {len(d)}; druggable (Finan) {d.druggable_finan.sum()}, "
      f"druggable (any) {d.druggable_any.sum()}, with approved drugs {(d.n_approved_drugs > 0).sum()}", flush=True)
