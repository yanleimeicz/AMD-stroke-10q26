"""Round-3 follow-up: variants of the UKB IHD credible set (lead rs4237541) and the UKB AMD lead (rs3750846):
LD with the 10q26 key variants (UKB British LD) and Z across Pan-UKB traits. Output: panukb_extra_variants.tsv, panukb_ld_extra.tsv"""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import numpy as np, pandas as pd
from scipy import sparse
P = PROJ_ROOT; D = f"{P}/derived/panukb"; O = f"{P}/results/revision3/panukb"
ldv = pd.read_csv(f"{D}/ld/chr10_123000001_126000001.gz", sep=r"\s+")
M = sparse.load_npz(f"{D}/ld/chr10_123000001_126000001.npz").tocsr()
EXTRA = {"rs4237541": 124237947, "rs3750846": 124215565, "rs28595632": 124403545}
KEYP = {"rs72631113": 124213449, "rs10490924": 124214448, "rs3750848": 124215315, "rs11200638": 124220544, "rs2672592": 124230750}
allp = {**KEYP, **EXTRA}
idx = {}
for rs, pos in allp.items():
    x = ldv[(ldv.position == pos) & (ldv.rsid == rs)]
    if len(x): idx[rs] = int(x.index[0])
ii = list(idx.values()); A = M[ii][:, ii].toarray(); A = A + A.T; np.fill_diagonal(A, 1)
names = list(idx)
rows = []
for a in EXTRA:
    if a not in idx: continue
    for b in names:
        if b == a: continue
        rows.append(dict(snp_a=a, snp_b=b, r_UKB_allele2_coded=A[names.index(a), names.index(b)],
                         a_alleles=f"{ldv.allele1[idx[a]]}/{ldv.allele2[idx[a]]}", b_alleles=f"{ldv.allele1[idx[b]]}/{ldv.allele2[idx[b]]}"))
pd.DataFrame(rows).to_csv(f"{O}/panukb_ld_extra.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).to_string())
ph = pd.read_csv(f"{O}/panukb_phenotypes.tsv", sep="\t")
out = []
for t in ph.tag:
    d = pd.read_csv(f"{D}/region_{t}.tsv", sep="\t")
    for rs, pos in EXTRA.items():
        y = d[d.pos == pos]
        for _, r in y.iterrows():
            out.append(dict(trait=t, rsid=rs, pos=pos, ref=r.ref, alt=r.alt, beta_alt=r.beta_EUR, se=r.se_EUR, z_alt=r.beta_EUR / r.se_EUR,
                            neglog10p=r.neglog10_pval_EUR, af_alt=r.af_controls_EUR))
out = pd.DataFrame(out); out.to_csv(f"{O}/panukb_extra_variants.tsv", sep="\t", index=False)
print(out.to_string())
# IHD credible set members
z = pd.read_csv(f"{O}/susie_input/IHD_phe411_z.tsv", sep="\t")
print(z.reindex(z.z.abs().sort_values(ascending=False).index).head(12).to_string())
