"""Round-3 follow-up: a fourth 10q26 marker, rs4237541 (UKB CAD/IHD credible-set lead, ~7 kb downstream of rs2672592),
looked up in GIGASTROKE/MEGASTROKE/FinnGen and fitted jointly with rs3750848, rs2672592, rs72631113 using UKB British LD.
Output: results/revision3/panukb/fourvar_lookup.tsv, fourvar_joint_models.tsv"""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import gzip, numpy as np, pandas as pd
from scipy import sparse, stats
from joblib import Parallel, delayed
P = PROJ_ROOT; D = f"{P}/derived/panukb"; O = f"{P}/results/revision3/panukb"; RAW = f"{P}/data/raw/gwas"
V = ["rs3750848", "rs2672592", "rs72631113", "rs4237541"]
RISK = {"rs3750848": "G", "rs2672592": "T", "rs72631113": "T", "rs4237541": "C"}
P37 = {"rs3750848": 124215315, "rs2672592": 124230750, "rs72631113": 124213449, "rs4237541": 124237947}
P38 = {k: v - 1_759_516 for k, v in P37.items()}
p2 = lambda z: 2 * stats.norm.sf(abs(z))
# UKB LD, risk-coded
ldv = pd.read_csv(f"{D}/ld/chr10_123000001_126000001.gz", sep=r"\s+")
M = sparse.load_npz(f"{D}/ld/chr10_123000001_126000001.npz").tocsr()
ii, sg = [], []
for v in V:
    x = ldv[(ldv.position == P37[v]) & (ldv.rsid == v)].iloc[0]
    ii.append(int(x.name)); sg.append(1 if x.allele2 == RISK[v] else -1)
A = M[ii][:, ii].toarray(); A = A + A.T; np.fill_diagonal(A, 1); A = A * np.outer(sg, sg)
print(pd.DataFrame(A, index=V, columns=V).round(3))

def gwascat(path):
    out = {}
    with gzip.open(path, "rt") as f:
        h = f.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
        cr, cb, ce = ("hm_rsid", "hm_beta", "hm_effect_allele") if "hm_rsid" in ix else ("rsid", "beta", "effect_allele")
        for line in f:
            x = line.rstrip("\n").split("\t"); r = x[ix[cr]]
            if r in RISK:
                b = float(x[ix[cb]]); se = float(x[ix["standard_error"]])
                if x[ix[ce]].upper() != RISK[r]: b = -b
                out[r] = (b, se)
                if len(out) == len(V): break
    return out
def finngen(path):
    out = {}
    with gzip.open(path, "rt") as f:
        h = f.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
        for line in f:
            if not line.startswith("10\t"): continue
            x = line.rstrip("\n").split("\t"); pos = int(x[ix["pos"]])
            for r in V:
                if pos == P38[r] and r in x[ix["rsids"]].split(","):
                    b = float(x[ix["beta"]]); se = float(x[ix["sebeta"]])
                    if x[ix["alt"]] != RISK[r]: b = -b
                    out[r] = (b, se)
            if pos > 122_500_000: break
    return out
FILES = {"GIGA_EUR_AS": (f"{RAW}/GIGASTROKE_GCST90104539.h.tsv.gz", gwascat), "GIGA_EUR_IS": (f"{RAW}/GIGASTROKE_GCST90104540.h.tsv.gz", gwascat),
         "MEGA_EUR_AS": (f"{RAW}/MEGASTROKE_AS_EUR_GCST006906.h.tsv.gz", gwascat), "MEGA_EUR_IS": (f"{RAW}/MEGASTROKE_IS_EUR_GCST006908.h.tsv.gz", gwascat),
         "FinnGen_I9_STR": (f"{RAW}/finngen_R12_I9_STR.full.gz", finngen), "FinnGen_C_STROKE": (f"{RAW}/finngen_R12_C_STROKE.gz", finngen),
         "FinnGen_H7_AMD": (f"{RAW}/finngen_R12_H7_AMD.gz", finngen)}
res = Parallel(n_jobs=len(FILES))(delayed(fn)(p) for p, fn in FILES.values())
rows = []
for k, d in zip(FILES, res):
    for r, (b, se) in d.items(): rows.append(dict(dataset=k, rsid=r, beta_risk=b, se=se, z=b / se))
# Pan-UKB
kv = pd.read_csv(f"{O}/panukb_10q26_variants.tsv", sep="\t"); ex = pd.read_csv(f"{O}/panukb_extra_variants.tsv", sep="\t")
for t in kv.trait.unique():
    for r in V[:3]:
        y = kv[(kv.trait == t) & (kv.rsid == r)]
        if len(y): rows.append(dict(dataset="UKB_" + t, rsid=r, beta_risk=y.beta_risk.iloc[0], se=y.se.iloc[0], z=y.z_risk.iloc[0]))
    y = ex[(ex.trait == t) & (ex.rsid == "rs4237541")]
    if len(y):
        s = 1 if y.alt.iloc[0] == "C" else -1
        rows.append(dict(dataset="UKB_" + t, rsid="rs4237541", beta_risk=s * y.beta_alt.iloc[0], se=y.se.iloc[0], z=s * y.z_alt.iloc[0]))
lk = pd.DataFrame(rows); lk.to_csv(f"{O}/fourvar_lookup.tsv", sep="\t", index=False)

out = []
for ds in lk.dataset.unique():
    d = lk[lk.dataset == ds].set_index("rsid")
    if not all(v in d.index for v in V): continue
    z = d.loc[V, "z"].values.astype(float); Ri = np.linalg.inv(A)
    zj = (Ri @ z) / np.sqrt(np.diag(Ri)); full = z @ Ri @ z
    def cond(keep):   # chi2 of the variants not in `keep`, given `keep`
        k = [V.index(x) for x in keep]; return full - z[k] @ np.linalg.solve(A[np.ix_(k, k)], z[k])
    c_side = cond(["rs3750848", "rs4237541"]); c_lead = cond(["rs2672592", "rs72631113", "rs4237541"]); c_4 = cond(V[:3])
    c_side3 = cond(["rs3750848"])
    o = dict(dataset=ds, ld="UKB_British (risk-coded)", **{f"z_{v}": z[i] for i, v in enumerate(V)}, **{f"zjoint4_{v}": zj[i] for i, v in enumerate(V)},
             **{f"pjoint4_{v}": p2(zj[i]) for i, v in enumerate(V)},
             p_side_given_lead_and_rs4237541_2df=stats.chi2.sf(c_side, 2), p_lead_given_side_and_rs4237541_1df=stats.chi2.sf(c_lead, 1),
             p_rs4237541_given_three_1df=stats.chi2.sf(c_4, 1), p_side_and_rs4237541_given_lead_3df=stats.chi2.sf(c_side3, 3),
             note=("in-sample LD" if ds.startswith("UKB") else "UKB LD used as an approximation for this dataset"))
    out.append(o)
out = pd.DataFrame(out); out.to_csv(f"{O}/fourvar_joint_models.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(lk.pivot(index="dataset", columns="rsid", values="z").round(2).to_string())
print(out[["dataset", "pjoint4_rs3750848", "pjoint4_rs2672592", "pjoint4_rs72631113", "pjoint4_rs4237541", "p_side_given_lead_and_rs4237541_2df",
           "p_lead_given_side_and_rs4237541_1df", "p_rs4237541_given_three_1df"]].to_string())
