"""Round-3: 10q26 symmetric conditional analysis in Pan-UK Biobank (EUR) with UKB in-sample LD (UKB British, Weissbrod 2020).
Same models as scripts/revision2/rev2_10q26.py. Inputs: derived/panukb/region_*.tsv, derived/panukb/ld/*.npz/.gz.
Outputs (results/revision3/panukb): panukb_10q26_variants.tsv, panukb_ld_keypairs.tsv, panukb_symmetric_conditional.tsv,
panukb_r_sensitivity.tsv, panukb_proportionality.tsv, heterogeneity_with_other_datasets.tsv, susie_input/*.tsv(+ R matrix)"""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import gzip, json, os
import numpy as np, pandas as pd
from scipy import sparse, stats

P = PROJ_ROOT; D = f"{P}/derived/panukb"; O = f"{P}/results/revision3/panukb"
R10 = f"{P}/results/revision/10q26"; R2 = f"{P}/results/revision2"
os.makedirs(f"{O}/susie_input", exist_ok=True)
p2 = lambda z: 2 * stats.norm.sf(abs(z))
S = ["rs3750848", "rs2672592", "rs72631113"]
# GRCh37 positions (hg38 + 1,759,516 in this region) and AMD-risk alleles
KEY = {"rs61871745": (124210369, "A"), "rs72631113": (124213449, "T"), "rs10490924": (124214448, "T"), "rs36212733": (124215211, "C"),
       "rs3750848": (124215315, "G"), "rs11200638": (124220544, "A"), "rs2672592": (124230750, "T")}
ph = pd.read_csv(f"{O}/panukb_phenotypes.tsv", sep="\t")

def load(tag):
    d = pd.read_csv(f"{D}/region_{tag}.tsv", sep="\t", dtype={"chr": str})
    d = d[d.beta_EUR.notna() & d.se_EUR.notna()].copy()
    d["z"] = d.beta_EUR / d.se_EUR
    d["af_EUR"] = d.af_controls_EUR
    return d
REG = {t: load(t) for t in ph.tag}

# ---------------------------------------------------------------- LD (UKB British; PolyFun format)
ldv = pd.read_csv(f"{D}/ld/chr10_123000001_126000001.gz", sep=r"\s+")
ldv.columns = [c.lower() for c in ldv.columns]
M = sparse.load_npz(f"{D}/ld/chr10_123000001_126000001.npz").tocsr()
diag = M.diagonal()
def Rsub(idx):
    A = M[idx][:, idx].toarray(); A = A + A.T
    np.fill_diagonal(A, 1.0); return A
print("LD variants", len(ldv), "diag range", diag.min(), diag.max(), flush=True)
posmap = {}
for i, (pos, a1, a2) in enumerate(zip(ldv.position, ldv.allele1, ldv.allele2)):
    posmap.setdefault(int(pos), []).append((i, a1, a2))
def ld_index(pos, ref, alt):
    """returns (index, sign): sign=+1 if LD is for dosage of `alt` (allele2==alt), -1 if swapped"""
    for i, a1, a2 in posmap.get(int(pos), []):
        if a1 == ref and a2 == alt: return i, 1
        if a1 == alt and a2 == ref: return i, -1
    return None, 0

# key-variant table, coded to the AMD-risk allele
rows = []; kidx = {}
base = REG[ph.tag.iloc[0]]
for rs, (pos, risk) in KEY.items():
    x = base[base.pos == pos]
    if not len(x): print("missing", rs); continue
    ref, alt = x.ref.iloc[0], x.alt.iloc[0]
    sg = 1 if alt == risk else -1
    i, s_ld = ld_index(pos, ref, alt)
    kidx[rs] = (i, s_ld * sg)     # LD coded to risk allele
    for t, d in REG.items():
        y = d[(d.pos == pos) & (d.ref == ref) & (d.alt == alt)]
        if not len(y): continue
        y = y.iloc[0]
        rows.append(dict(rsid=rs, pos_grch37=pos, ref=ref, alt=alt, risk=risk, trait=t, beta_risk=sg * y.beta_EUR, se=y.se_EUR,
                         z_risk=sg * y.z, p=10 ** (-y.neglog10_pval_EUR), af_risk_controls=(y.af_EUR if sg == 1 else 1 - y.af_EUR),
                         low_confidence=y.low_confidence_EUR))
kv = pd.DataFrame(rows); kv.to_csv(f"{O}/panukb_10q26_variants.tsv", sep="\t", index=False)

# key-pair r: UKB vs 1000G EUR vs FinnGen SISu
kg = pd.read_csv(f"{R10}/ld_1000G_EUR_r_dprime.tsv", sep=r"\s+")
fg = json.load(open(f"{R10}/ld_finngen_sisu42_pairs.json"))
def kg_r(a, b):
    x = kg[((kg.SNP_A == a) & (kg.SNP_B == b)) | ((kg.SNP_A == b) & (kg.SNP_B == a))]
    return float(x.R.iloc[0]) if len(x) else np.nan
def fg_r(a, b):
    k = f"{a}|{b}" if f"{a}|{b}" in fg else f"{b}|{a}"
    return np.sqrt(fg[k]["r2"]) if k in fg else np.nan
def ukb_r(a, b):
    (i, sa), (j, sb) = kidx[a], kidx[b]
    if i is None or j is None: return np.nan
    A = Rsub([i, j]); return sa * sb * A[0, 1]
lp = []
names = [k for k in KEY if k in kidx]
for ii, a in enumerate(names):
    for b in names[ii + 1:]:
        lp.append(dict(snp_a=a, snp_b=b, r_UKB_British=ukb_r(a, b), r_1000G_EUR=kg_r(a, b), r_FinnGen_SISu=fg_r(a, b)))
lp = pd.DataFrame(lp); lp.to_csv(f"{O}/panukb_ld_keypairs.tsv", sep="\t", index=False)
print(lp.to_string(), flush=True)

pairs = [(S[0], S[1]), (S[0], S[2]), (S[1], S[2])]
LD = {"UKB_British": [ukb_r(*p) for p in pairs], "KG_EUR": [kg_r(*p) for p in pairs], "FinnGen_SISu": [fg_r(*p) for p in pairs]}
def Rmat(r): return np.array([[1, r[0], r[1]], [r[0], 1, r[2]], [r[1], r[2], 1]])

def models(z, R):
    out = {}
    for i, j in [(1, 0), (0, 1), (2, 0), (0, 2), (2, 1), (1, 2)]:
        r = R[i, j]; zc = (z[i] - r * z[j]) / np.sqrt(1 - r * r)
        out[f"z_{S[i]}|{S[j]}"] = zc; out[f"p_{S[i]}|{S[j]}"] = p2(zc)
    Ri = np.linalg.inv(R); zj = (Ri @ z) / np.sqrt(np.diag(Ri))
    for k in range(3): out[f"zjoint3_{S[k]}"] = zj[k]; out[f"pjoint3_{S[k]}"] = p2(zj[k])
    chi_full = z @ Ri @ z
    chi_side = chi_full - z[0] ** 2
    chi_lead = chi_full - z[1:] @ np.linalg.solve(R[1:, 1:], z[1:])
    out.update(chi2_side_given_lead_2df=chi_side, p_side_given_lead_2df=stats.chi2.sf(chi_side, 2),
               chi2_lead_given_side_1df=chi_lead, p_lead_given_side_1df=stats.chi2.sf(chi_lead, 1))
    for k in (1, 2):
        R2_ = R[np.ix_([0, k], [0, k])]; V = np.linalg.inv(R2_); b = V @ z[[0, k]]
        g = np.array([-b[1] / b[0] ** 2, 1 / b[0]])
        out[f"ratio_{S[k]}_to_lead"] = b[1] / b[0]; out[f"ratio_{S[k]}_to_lead_se"] = np.sqrt(g @ V @ g)
    return out
def zvec(t):
    d = kv[kv.trait == t].set_index("rsid")
    return np.array([d.loc[s, "z_risk"] for s in S], float) if all(s in d.index for s in S) else None

sym, sens = [], []
for t in ph.tag:
    z = zvec(t)
    if z is None: continue
    for ld in ("UKB_British", "KG_EUR", "FinnGen_SISu"):
        r = LD[ld]; o = dict(trait=t, ld_source=ld, ld_matched=(ld == "UKB_British"), r_lead_2672592=r[0], r_lead_72631113=r[1],
                             r_2672592_72631113=r[2], z_rs3750848=z[0], z_rs2672592=z[1], z_rs72631113=z[2])
        o.update(models(z, Rmat(r))); sym.append(o)
    r0 = np.array(LD["UKB_British"])
    for mode in ("both_lead_side", "lead_2672592_only"):
        for dl in np.round(np.arange(-0.10, 0.1001, 0.01), 3):
            r = r0.copy(); r[0] += dl
            if mode == "both_lead_side": r[1] += dl
            R = Rmat(r)
            if np.linalg.eigvalsh(R).min() < 1e-6 or r.max() >= 0.999: continue
            m = models(z, R)
            sens.append(dict(trait=t, ld_source="UKB_British", mode=mode, delta=dl, r_lead_2672592=r[0], r_lead_72631113=r[1],
                             z_rs2672592_given_lead=m["z_rs2672592|rs3750848"], z_lead_given_rs2672592=m["z_rs3750848|rs2672592"],
                             zjoint3_lead=m["zjoint3_rs3750848"], zjoint3_rs2672592=m["zjoint3_rs2672592"], zjoint3_rs72631113=m["zjoint3_rs72631113"],
                             p_side_given_lead_2df=m["p_side_given_lead_2df"], p_lead_given_side_1df=m["p_lead_given_side_1df"],
                             ratio_2672592_to_lead=m["ratio_rs2672592_to_lead"]))
sym = pd.DataFrame(sym); sym.to_csv(f"{O}/panukb_symmetric_conditional.tsv", sep="\t", index=False)
sens = pd.DataFrame(sens); sens.to_csv(f"{O}/panukb_r_sensitivity.tsv", sep="\t", index=False)

# proportionality (projected): UKB AMD joint effects (UKB LD) vs UKB stroke (UKB LD) -> fully in-sample; also FinnGen AMD (FinnGen LD) vs UKB stroke
def prop_test(zA, RA, zS, RS, ridge=0.0):
    RS = RS + ridge * np.eye(3); RA = RA + ridge * np.eye(3)
    m = RS @ np.linalg.solve(RA, zA); W = np.linalg.inv(RS)
    k = (m @ W @ zS) / (m @ W @ m); e = zS - k * m; chi = e @ W @ e
    return k, chi, stats.chi2.sf(chi, 2)
kvo = pd.read_csv(f"{R10}/key_variants_all_datasets.tsv", sep="\t")
def zother(t):
    d = kvo[kvo.trait == t].set_index("rsid"); return np.array([d.loc[s, "z_risk"] for s in S], float)
prop = []
amds = [("UKB " + t, zvec(t), "UKB_British") for t in ph.tag if t.startswith("AMD")] + [("FinnGen_H7_AMD", zother("FinnGen_H7_AMD"), "FinnGen_SISu")]
for an, zA, lda in amds:
    for t in ph.tag:
        if t.startswith("AMD"): continue
        zS = zvec(t)
        for ridge in (0, 0.02, 0.05):
            k, chi, p = prop_test(zA, Rmat(LD[lda]), zS, Rmat(LD["UKB_British"]), ridge)
            prop.append(dict(amd=an, amd_ld=lda, stroke="UKB " + t, stroke_ld="UKB_British", ridge=ridge, k=k, chi2=chi, p=p))
pd.DataFrame(prop).to_csv(f"{O}/panukb_proportionality.tsv", sep="\t", index=False)

# heterogeneity / power vs GIGASTROKE, MEGASTROKE, GIGA-minus-MEGA-FinnGen, FinnGen
h3 = pd.read_csv(f"{R2}/10q26_stroke_heterogeneity_3parts.tsv", sep="\t")
pw = lambda ez: stats.norm.cdf(abs(ez) - 1.96) + stats.norm.cdf(-abs(ez) - 1.96)
het = []
for ukt, sub in [("IS_I63", "IS"), ("IS_phe433.2", "IS"), ("CI_phe433.3", "IS"), ("AS_sr6150_3", "AS"), ("AS_sr20002_1081", "AS"), ("CVD_phe433", "AS")]:
    for s in S:
        u = kv[(kv.trait == ukt) & (kv.rsid == s)]
        if not len(u): continue
        bU, sU = float(u.beta_risk.iloc[0]), float(u.se.iloc[0])
        g = kvo[(kvo.trait == f"GIGA_EUR_{sub}") & (kvo.rsid == s)]
        bG, sG = float(g.beta_risk.iloc[0]), float(g.se.iloc[0])
        r3 = h3[(h3.subtype == sub) & (h3.rsid == s)].iloc[0]
        comp = {"GIGASTROKE_EUR": (bG, sG), "MEGASTROKE": (r3.beta_MEGA, r3.se_MEGA),
                "GIGA_minus_MEGA_minus_FinnGen": (r3.beta_GIGA_minus_MEGA_FinnGen, r3.se_GIGA_minus_MEGA_FinnGen),
                "FinnGen_I9_STR": (r3.beta_FinnGen_I9_STR, r3.se_FinnGen)}
        o = dict(ukb_trait=ukt, compared_subtype=sub, rsid=s, beta_ukb=bU, se_ukb=sU, z_ukb=bU / sU, p_ukb=p2(bU / sU))
        for k, (b, se) in comp.items():
            o[f"beta_{k}"] = b; o[f"se_{k}"] = se
            o[f"p_diff_ukb_vs_{k}"] = p2((bU - b) / np.hypot(sU, se))
            o[f"expected_z_ukb_if_{k}"] = b / sU; o[f"power_ukb_if_{k}"] = pw(b / sU)
        bs = np.array([comp["MEGASTROKE"][0], comp["GIGA_minus_MEGA_minus_FinnGen"][0], comp["FinnGen_I9_STR"][0], bU])
        ss = np.array([comp["MEGASTROKE"][1], comp["GIGA_minus_MEGA_minus_FinnGen"][1], comp["FinnGen_I9_STR"][1], sU])
        w = 1 / ss ** 2; bm = (w * bs).sum() / w.sum(); Q = (w * (bs - bm) ** 2).sum()
        o.update(Q_4parts=Q, p_Q_3df=stats.chi2.sf(Q, 3), I2=max(0, (Q - 3) / Q) if Q > 0 else 0, ratio_ukb_to_MEGA=bU / comp["MEGASTROKE"][0],
                 note="UKB contributed ~2,932 stroke cases to GIGASTROKE EUR (and to its remainder after subtracting MEGASTROKE); comparisons with GIGASTROKE are not fully independent")
        het.append(o)
het = pd.DataFrame(het); het.to_csv(f"{O}/heterogeneity_with_other_datasets.tsv", sep="\t", index=False)

# AMD effect comparison (UKB vs FinnGen / IAMDGC) for the side/lead ratio (two-variant model)
# SuSiE inputs: region 124.0-124.45 Mb, MAF >= 0.005, matched to LD by position+alleles
for t in ["AMD_phe362.29", "AMD_phe362.2", "AMD_sr6148_5", "IS_I63", "IS_phe433.2", "AS_sr6150_3", "CVD_phe433", "IHD_phe411"]:
    d = REG[t]; d = d[(d.pos >= 124_000_000) & (d.pos <= 124_450_000) & (d.af_EUR.between(0.005, 0.995))].copy()
    if "low_confidence_EUR" in d: d = d[d.low_confidence_EUR.astype(str) != "True"]
    ix = [ld_index(p, r, a) for p, r, a in zip(d.pos, d.ref, d.alt)]
    d["ld_i"] = [i for i, s in ix]; d["ld_s"] = [s for i, s in ix]
    d = d[d.ld_i.notna()].drop_duplicates("ld_i").sort_values("pos")
    idx = d.ld_i.astype(int).tolist(); sg = d.ld_s.values.astype(float)
    Rm = Rsub(idx) * np.outer(sg, sg)
    n = int(ph.loc[ph.tag == t, "n_cases_EUR"].iloc[0]) + int(ph.loc[ph.tag == t, "n_controls_EUR"].iloc[0])
    d[["pos", "ref", "alt", "af_EUR", "beta_EUR", "se_EUR", "z"]].assign(n=n).to_csv(f"{O}/susie_input/{t}_z.tsv", sep="\t", index=False)
    Rm.astype("<f8").tofile(f"{O}/susie_input/{t}_R.f64")
    print(t, "susie input", len(d), "n", n, flush=True)

print(sym[sym.ld_matched][["trait", "z_rs3750848", "z_rs2672592", "z_rs72631113", "p_rs2672592|rs3750848", "p_side_given_lead_2df",
                           "p_lead_given_side_1df", "pjoint3_rs3750848", "pjoint3_rs2672592", "pjoint3_rs72631113", "ratio_rs2672592_to_lead"]].to_string())
