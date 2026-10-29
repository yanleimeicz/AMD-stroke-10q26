"""Round 5: network fetches.
(1) GWAS Catalog associations for atrial fibrillation (Nielsen 2018, GCST006414; fallback Roselli 2018, GCST006061) -> AF instruments
(2) Ensembl protein-coding genes within +/-250 kb of the complement-annotated AMD instruments (gene-based pathway membership)
Outputs: results/revision5/af_gwascatalog_associations.tsv, results/revision5/complement_instrument_genes.tsv"""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import json, time, re
import requests, pandas as pd

O = f"{PROJ_ROOT}/results/revision5"
S = requests.Session(); S.trust_env = False

def get(url, params=None, tries=6):
    for i in range(tries):
        try:
            r = S.get(url, params=params, timeout=120, headers={"Accept": "application/json"})
            if r.status_code == 200: return r.json()
            print("status", r.status_code, url)
        except Exception as e:
            print("error", e)
        time.sleep(5 * (i + 1))
    raise RuntimeError(url)

# ---------------------------------------------------------------- (1) AF associations
rows = []
for acc in ("GCST006414", "GCST006061"):
    st = get(f"https://www.ebi.ac.uk/gwas/rest/api/studies/{acc}")
    print(acc, st.get("publicationInfo", {}).get("title"), st.get("initialSampleSize"))
    url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/{acc}/associations"
    js = get(url, params={"projection": "associationByStudy"})
    for a in js["_embedded"]["associations"]:
        for loc in a.get("loci", []):
            for ra in loc.get("strongestRiskAlleles", []):
                nm = ra.get("riskAlleleName", "")
                m = re.match(r"(rs\d+)-([ACGT?]+)", nm)
                if not m: continue
                rows.append(dict(study=acc, SNP=m.group(1), risk_allele=m.group(2), p=a.get("pvalue"),
                                 or_per_copy=a.get("orPerCopyNum"), beta=a.get("betaNum"), beta_dir=a.get("betaDirection"),
                                 ci=a.get("range"), raf=a.get("riskFrequency")))
    if rows:
        af = pd.DataFrame(rows); af.to_csv(f"{O}/af_gwascatalog_associations.tsv", sep="\t", index=False)
        print(acc, "associations:", len(af)); break

# ---------------------------------------------------------------- (2) genes near complement-annotated instruments (GRCh38)
inst = {"rs12354144": ("1", 195301905), "rs1329424": ("1", 196677046), "rs72744971": ("1", 197888054),
        "rs139779213": ("4", 109740713), "rs429608": ("6", 31962685), "rs2230199": ("19", 6718376),
        "rs11080055": ("17", 28322698), "rs6665824": ("1", 196851104), "rs114920504": ("1", 197941325),
        "rs11732475": ("4", 109732639), "rs635348": ("6", 31934512), "rs11569415": ("19", 6716268)}
COMP = {"CFH", "CFHR1", "CFHR2", "CFHR3", "CFHR4", "CFHR5", "CFI", "CFB", "C2", "C3", "C4A", "C4B", "C5", "C9", "VTN", "CD46", "CD55", "CR1", "CR2", "F13B"}
gr = []
for snp, (c, p) in inst.items():
    js = get(f"https://rest.ensembl.org/overlap/region/human/{c}:{max(1, p - 250000)}-{p + 250000}",
             params={"feature": "gene", "biotype": "protein_coding", "content-type": "application/json"})
    genes = sorted({g.get("external_name") for g in js if g.get("external_name")}, key=lambda x: x)
    near = None; dmin = 1e12
    for g in js:
        if not g.get("external_name"): continue
        d = 0 if g["start"] <= p <= g["end"] else min(abs(g["start"] - p), abs(g["end"] - p))
        if d < dmin: dmin, near = d, g["external_name"]
    comp = sorted(set(genes) & COMP)
    gr.append(dict(SNP=snp, chr=c, pos_grch38=p, nearest_protein_coding_gene=near, distance_bp=dmin,
                   complement_genes_within_250kb=";".join(comp), n_genes_within_250kb=len(genes), genes_within_250kb=";".join(genes),
                   complement_member_gene_based=bool(comp)))
    time.sleep(0.4)
g = pd.DataFrame(gr); g.to_csv(f"{O}/complement_instrument_genes.tsv", sep="\t", index=False)
print(g[["SNP", "nearest_protein_coding_gene", "distance_bp", "complement_genes_within_250kb", "complement_member_gene_based"]].to_string())
