"""Round 5, item 2: 10q26 stroke haplotype assignment.
Per dataset (two-variant joint model rs3750848 [lead, AMD haplotype] + rs2672592 [side]):
  contrast d = beta_side - beta_lead, var(d) = s_side^2 + s_lead^2 - 2 cov, cov = -r * s_lead * s_side
  (correlation of the two joint estimates is -r, r = LD between the markers in the dataset's matched panel).
Pooling: DerSimonian-Laird (z), REML (z and t_{k-1}), Hartung-Knapp-Sidik-Jonkman (HKSJ, t_{k-1}); prediction intervals.
Sets: main (MEGASTROKE, GIGASTROKE remainder, FinnGen I9_STR, UKB), no_remainder, full_GIGA (full GIGASTROKE EUR + FinnGen + UKB; overlapping).
Output: results/revision5/10q26_stroke_contrast_hksj.tsv"""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import json, numpy as np, pandas as pd
from scipy import stats, optimize

P = PROJ_ROOT; O = f"{P}/results/revision5"
R10 = f"{P}/results/revision/10q26"; R3 = f"{P}/results/revision3/panukb"; R4 = f"{P}/results/revision4"
S = ["rs3750848", "rs2672592", "rs72631113"]
p2 = lambda z: 2 * stats.norm.sf(abs(z))

fg = json.load(open(f"{R10}/ld_finngen_sisu42_pairs.json"))
def fg_r(a, b):
    k = f"{a}|{b}" if f"{a}|{b}" in fg else f"{b}|{a}"; return np.sqrt(fg[k]["r2"])
kg = pd.read_csv(f"{R10}/ld_1000G_EUR_r_dprime.tsv", sep=r"\s+")
def kg_r(a, b):
    x = kg[((kg.SNP_A == a) & (kg.SNP_B == b)) | ((kg.SNP_A == b) & (kg.SNP_B == a))]; return float(x.R.iloc[0])
ukl = pd.read_csv(f"{R3}/panukb_ld_keypairs.tsv", sep="\t")
def uk_r(a, b):
    x = ukl[((ukl.snp_a == a) & (ukl.snp_b == b)) | ((ukl.snp_a == b) & (ukl.snp_b == a))]; return float(x.r_UKB_British.iloc[0])
R = {"KG_EUR": kg_r(S[0], S[1]), "FinnGen_SISu": fg_r(S[0], S[1]), "UKB_British": uk_r(S[0], S[1])}
print("r(lead, rs2672592):", R)

ds = pd.read_csv(f"{R4}/10q26_stroke_datasets.tsv", sep="\t")
kv = pd.read_csv(f"{R10}/key_variants_all_datasets.tsv", sep="\t")

def joint2(b, se, r):
    z = b / se; Rm = np.array([[1, r], [r, 1]]); V = np.linalg.inv(Rm); bs = V @ z
    return bs * se, np.sqrt(np.diag(V)) * se

def per_dataset(sub):
    out = []
    for nm in ("MEGASTROKE", "GIGASTROKE_remainder", "FinnGen_I9_STR", "UKB"):
        d = ds[(ds.subtype == sub) & (ds.dataset == nm)]
        if d.empty: continue
        ld = d.ld.iloc[0]
        bl, sl = d[d.quantity == "joint_rs2672592_lead"][["beta", "se"]].values[0]
        bs, ss = d[d.quantity == "joint_rs2672592_side"][["beta", "se"]].values[0]
        out.append(dict(subtype=sub, dataset=nm, ld=ld, r=R[ld], b_lead=bl, s_lead=sl, b_side=bs, s_side=ss))
    # full GIGASTROKE EUR (overlaps MEGASTROKE, remainder, FinnGen share and UKB share)
    t = kv[kv.trait == f"GIGA_EUR_{sub}"].set_index("rsid")
    if len(t):
        b = t.loc[S[:2], "beta_risk"].values.astype(float); se = t.loc[S[:2], "se"].values.astype(float)
        bj, sj = joint2(b, se, R["KG_EUR"])
        out.append(dict(subtype=sub, dataset="GIGASTROKE_EUR_full", ld="KG_EUR", r=R["KG_EUR"], b_lead=bj[0], s_lead=sj[0], b_side=bj[1], s_side=sj[1]))
    x = pd.DataFrame(out)
    x["cov"] = -x.r * x.s_lead * x.s_side
    x["b_contrast"] = x.b_side - x.b_lead
    x["s_contrast"] = np.sqrt(x.s_side ** 2 + x.s_lead ** 2 - 2 * x["cov"])
    x["p_contrast"] = p2(x.b_contrast / x.s_contrast)
    x["p_lead"] = p2(x.b_lead / x.s_lead); x["p_side"] = p2(x.b_side / x.s_side)
    return x

def meta(y, s):
    y, s = np.asarray(y, float), np.asarray(s, float); k = len(y); w = 1 / s ** 2
    mu_f = (w * y).sum() / w.sum(); Q = (w * (y - mu_f) ** 2).sum(); C = w.sum() - (w ** 2).sum() / w.sum()
    t2_dl = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
    def nll(t2):
        v = s ** 2 + t2; wi = 1 / v; mu = (wi * y).sum() / wi.sum()
        return 0.5 * (np.log(v).sum() + np.log(wi.sum()) + (wi * (y - mu) ** 2).sum())
    t2_reml = optimize.minimize_scalar(nll, bounds=(0, max(1e-10, 10 * np.var(y) + 1e-6)), method="bounded").x
    res = {}
    tq = stats.t.ppf(0.975, k - 1)
    for lab, t2 in (("DL", t2_dl), ("REML", t2_reml)):
        wi = 1 / (s ** 2 + t2); mu = (wi * y).sum() / wi.sum(); se = 1 / np.sqrt(wi.sum())
        res[f"{lab}_est"] = mu; res[f"{lab}_se"] = se; res[f"{lab}_p_z"] = p2(mu / se)
        res[f"{lab}_lo"] = mu - 1.96 * se; res[f"{lab}_hi"] = mu + 1.96 * se; res[f"{lab}_tau2"] = t2
        if lab == "REML":
            res["REML_p_t"] = 2 * stats.t.sf(abs(mu / se), k - 1); res["REML_lo_t"] = mu - tq * se; res["REML_hi_t"] = mu + tq * se
        # HKSJ variance (based on this tau2)
        q = (wi * (y - mu) ** 2).sum() / (k - 1); se_h = np.sqrt(q / wi.sum())
        se_h_mod = max(se_h, se)          # Knapp-Hartung ad hoc modification (never narrower than the z-based SE)
        res[f"HKSJ_{lab}_se"] = se_h; res[f"HKSJ_{lab}_p"] = 2 * stats.t.sf(abs(mu / se_h), k - 1)
        res[f"HKSJ_{lab}_lo"] = mu - tq * se_h; res[f"HKSJ_{lab}_hi"] = mu + tq * se_h
        res[f"HKSJmod_{lab}_p"] = 2 * stats.t.sf(abs(mu / se_h_mod), k - 1)
        tp = stats.t.ppf(0.975, k - 2) if k > 2 else np.nan
        res[f"{lab}_pi_lo"] = mu - tp * np.sqrt(t2 + se ** 2) if k > 2 else np.nan
        res[f"{lab}_pi_hi"] = mu + tp * np.sqrt(t2 + se ** 2) if k > 2 else np.nan
    res["FE_est"] = mu_f; res["FE_p"] = p2(mu_f * np.sqrt(w.sum())); res["Q"] = Q
    res["Q_p"] = stats.chi2.sf(Q, k - 1) if k > 1 else np.nan; res["I2"] = max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0
    res["k"] = k; res["n_negative"] = int((y < 0).sum())
    return res

rows, dsrows = [], []
SETS = {"main: MEGASTROKE + GIGASTROKE remainder + FinnGen + UKB": ["MEGASTROKE", "GIGASTROKE_remainder", "FinnGen_I9_STR", "UKB"],
        "sensitivity: without GIGASTROKE remainder": ["MEGASTROKE", "FinnGen_I9_STR", "UKB"],
        "sensitivity: full GIGASTROKE EUR + FinnGen + UKB (overlapping)": ["GIGASTROKE_EUR_full", "FinnGen_I9_STR", "UKB"]}
for sub in ("IS", "AS"):
    x = per_dataset(sub); dsrows.append(x)
    for setname, members in SETS.items():
        xx = x[x.dataset.isin(members)]
        for q, (bcol, scol) in {"contrast_side_minus_lead": ("b_contrast", "s_contrast"), "side_rs2672592_given_lead": ("b_side", "s_side"),
                                "lead_rs3750848_given_side": ("b_lead", "s_lead")}.items():
            r = meta(xx[bcol], xx[scol]); r.update(subtype=sub, set=setname, quantity=q, datasets=";".join(xx.dataset)); rows.append(r)
res = pd.DataFrame(rows)
cols = ["subtype", "set", "quantity", "k", "datasets", "n_negative", "FE_est", "FE_p", "DL_est", "DL_se", "DL_lo", "DL_hi", "DL_p_z", "DL_tau2",
        "REML_est", "REML_se", "REML_p_z", "REML_p_t", "REML_lo_t", "REML_hi_t", "REML_tau2",
        "HKSJ_DL_p", "HKSJ_DL_lo", "HKSJ_DL_hi", "HKSJmod_DL_p", "HKSJ_REML_p", "HKSJ_REML_lo", "HKSJ_REML_hi", "HKSJmod_REML_p",
        "DL_pi_lo", "DL_pi_hi", "REML_pi_lo", "REML_pi_hi", "Q", "Q_p", "I2"]
res = res[cols]
dsx = pd.concat(dsrows, ignore_index=True)
with open(f"{O}/10q26_stroke_contrast_hksj.tsv", "w", encoding="utf-8") as f:
    res.to_csv(f, sep="\t", index=False)
dsx.to_csv(f"{O}/10q26_stroke_contrast_datasets.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(dsx.round(4).to_string())
print(res[["subtype", "set", "quantity", "k", "DL_est", "DL_p_z", "REML_p_t", "HKSJ_DL_p", "HKSJ_DL_lo", "HKSJ_DL_hi", "HKSJmod_DL_p", "HKSJ_REML_p", "DL_pi_lo", "DL_pi_hi", "I2", "Q_p"]].round(4).to_string())
