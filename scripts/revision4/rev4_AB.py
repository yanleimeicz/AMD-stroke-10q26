"""Round-4 reviewer response, items A and B.
A. 10q26: random-effects meta-analysis of two-variant joint effects (rs3750848 + rs2672592; rs3750848 + rs72631113) and of
   marginal effects across approximately independent European stroke datasets (MEGASTROKE; GIGASTROKE remainder after
   approximate inverse-variance subtraction of MEGASTROKE, FinnGen and UK Biobank; FinnGen I9_STR; UK Biobank), each with
   its matched LD; East Asian GIGASTROKE shown separately. DerSimonian-Laird and REML, prediction intervals, leave-one-out.
   Joint effects on the log-OR scale: b_std = R^-1 z (standardized joint effects, covariance R^-1), beta_joint_k = b_std_k * se_k,
   SE_joint_k = sqrt([R^-1]_kk) * se_k (summary-statistic approximation of COJO; se_k ~ 1/sqrt(2 p_k q_k N)).
B. Bound-free reporting of MR equivalence: smallest symmetric bound at which TOST P < 0.05; one-sided tests of the MR estimate
   against the value expected if the NHANES association were causal (Burgess & Labrecque: beta_exp = log(OR_obs) * p(1-p));
   TOST P as a function of the bound.
Outputs: results/revision4/10q26_stroke_random_effects.tsv, 10q26_stroke_loo.tsv, 10q26_stroke_datasets.tsv,
         equivalence_boundfree.tsv, causal_expectation_test.tsv, tost_curve.tsv"""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import json, os
import numpy as np, pandas as pd
from scipy import stats, optimize

P = PROJ_ROOT; O = f"{P}/results/revision4"; os.makedirs(O, exist_ok=True)
R10 = f"{P}/results/revision/10q26"; R2 = f"{P}/results/revision2"; R3 = f"{P}/results/revision3/panukb"
S = ["rs3750848", "rs2672592", "rs72631113"]
p2 = lambda z: 2 * stats.norm.sf(abs(z))

# ================================================================ A. data
kv = pd.concat([pd.read_csv(f"{R10}/key_variants_all_datasets.tsv", sep="\t"),
                pd.read_csv(f"{R2}/finngen_extra_endpoints_10q26.tsv", sep="\t")], ignore_index=True)
uk = pd.read_csv(f"{R3}/panukb_10q26_variants.tsv", sep="\t")
def get(t, src="kv"):
    d = (kv if src == "kv" else uk); d = d[d.trait == t].set_index("rsid")
    return d.loc[S, "beta_risk"].values.astype(float), d.loc[S, "se"].values.astype(float)

fg = json.load(open(f"{R10}/ld_finngen_sisu42_pairs.json"))
def fg_r(a, b):
    k = f"{a}|{b}" if f"{a}|{b}" in fg else f"{b}|{a}"; return np.sqrt(fg[k]["r2"])
kg = pd.read_csv(f"{R10}/ld_1000G_EUR_r_dprime.tsv", sep=r"\s+")
def kg_r(a, b):
    x = kg[((kg.SNP_A == a) & (kg.SNP_B == b)) | ((kg.SNP_A == b) & (kg.SNP_B == a))]; return float(x.R.iloc[0])
ukl = pd.read_csv(f"{R3}/panukb_ld_keypairs.tsv", sep="\t")
def uk_r(a, b):
    x = ukl[((ukl.snp_a == a) & (ukl.snp_b == b)) | ((ukl.snp_a == b) & (ukl.snp_b == a))]; return float(x.r_UKB_British.iloc[0])
pairs = [(S[0], S[1]), (S[0], S[2]), (S[1], S[2])]
LD = {"KG_EUR": [kg_r(*p) for p in pairs], "FinnGen_SISu": [fg_r(*p) for p in pairs], "UKB_British": [uk_r(*p) for p in pairs],
      "KG_EAS": [0.821883203381113, 0.678589714039345, 0.541761017423735]}
def Rmat(r): return np.array([[1, r[0], r[1]], [r[0], 1, r[2]], [r[1], r[2], 1]])

# UKB information share in GIGASTROKE EUR: 2,932 cases / 369,419 controls (Mishra 2022 Suppl. Table 1)
def neff(nc, nk): return 4 / (1 / nc + 1 / nk)
F_UKB = {"AS": neff(2932, 369419) / neff(73652, 1234808), "IS": neff(2932, 369419) / neff(62100, 1234808)}

def remainder(sub, fF=0.15, sub_ukb=True):
    bG, sG = get(f"GIGA_EUR_{sub}"); bM, sM = get(f"MEGA_EUR_{sub}"); bF, sF = get("FinnGen_I9_STR")
    IG, IM = 1 / sG ** 2, 1 / sM ** 2; IF = fF * IG
    num = IG * bG - IM * bM - IF * bF; IR = IG - IM - IF
    if sub_ukb and sub in F_UKB:
        bU, sU = get("IS_I63" if sub == "IS" else "AS_sr6150_3", "uk")
        IU = F_UKB[sub] * IG; num = num - IU * bU; IR = IR - IU
    return num / IR, 1 / np.sqrt(IR)

def datasets(sub, fF=0.15, sub_ukb=True):
    ds = {}
    ds["MEGASTROKE"] = (*get(f"MEGA_EUR_{sub}"), "KG_EUR")
    ds["GIGASTROKE_remainder"] = (*remainder(sub, fF, sub_ukb), "KG_EUR")
    if sub in ("IS", "AS"):
        ds["FinnGen_I9_STR"] = (*get("FinnGen_I9_STR"), "FinnGen_SISu")
        ds["UKB"] = (*get("IS_I63" if sub == "IS" else "AS_sr6150_3", "uk"), "UKB_British")
    return ds

def joint(b, se, r, k):
    """two-variant joint model of lead (index 0) and marker k (1 or 2)"""
    idx = [0, k]; z = b[idx] / se[idx]
    R2 = np.array([[1, r[k - 1]], [r[k - 1], 1]]); V = np.linalg.inv(R2); bs = V @ z
    bj = bs * se[idx]; sej = np.sqrt(np.diag(V)) * se[idx]
    return bj, sej

# ---------------------------------------------------------------- random-effects machinery
def re_meta(y, s):
    y, s = np.asarray(y, float), np.asarray(s, float); k = len(y); w = 1 / s ** 2
    mu_f = (w * y).sum() / w.sum(); Q = (w * (y - mu_f) ** 2).sum(); dfq = k - 1
    C = w.sum() - (w ** 2).sum() / w.sum()
    t2_dl = max(0.0, (Q - dfq) / C) if C > 0 else 0.0
    def nll(t2):                                      # REML negative log-likelihood
        v = s ** 2 + t2; wi = 1 / v; mu = (wi * y).sum() / wi.sum()
        return 0.5 * (np.log(v).sum() + np.log(wi.sum()) + (wi * (y - mu) ** 2).sum())
    t2_reml = optimize.minimize_scalar(nll, bounds=(0, max(1e-12, 10 * np.var(y) + 1e-6)), method="bounded").x if k > 1 else 0.0
    out = {}
    for lab, t2 in (("DL", t2_dl), ("REML", t2_reml)):
        wi = 1 / (s ** 2 + t2); mu = (wi * y).sum() / wi.sum(); se = 1 / np.sqrt(wi.sum())
        tq = stats.t.ppf(0.975, k - 2) if k > 2 else np.nan
        out[lab] = dict(est=mu, se=se, lo=mu - 1.96 * se, hi=mu + 1.96 * se, p=p2(mu / se), tau2=t2,
                        pi_lo=mu - tq * np.sqrt(t2 + se ** 2) if k > 2 else np.nan, pi_hi=mu + tq * np.sqrt(t2 + se ** 2) if k > 2 else np.nan)
    out["FE"] = dict(est=mu_f, se=1 / np.sqrt(w.sum()), p=p2(mu_f * np.sqrt(w.sum())))
    out["Q"] = Q; out["Q_p"] = stats.chi2.sf(Q, dfq) if dfq > 0 else np.nan; out["I2"] = max(0.0, (Q - dfq) / Q) if Q > 0 else 0.0
    return out

def flat(o, prefix=""):
    d = {f"{prefix}FE_est": o["FE"]["est"], f"{prefix}FE_p": o["FE"]["p"], f"{prefix}Q": o["Q"], f"{prefix}Q_p": o["Q_p"], f"{prefix}I2": o["I2"]}
    for lab in ("DL", "REML"):
        for k, v in o[lab].items(): d[f"{prefix}{lab}_{k}"] = v
    return d

rows, loo, dsrows = [], [], []
for sub in ("IS", "AS", "SVS"):
    for variant_label, fF, sub_ukb in (("main (remainder minus MEGASTROKE, FinnGen 0.15, UKB)", 0.15, True),
                                        ("sensitivity: FinnGen share 0.10", 0.10, True),
                                        ("sensitivity: FinnGen share 0.20", 0.20, True),
                                        ("sensitivity: UKB not subtracted", 0.15, False)):
        if sub == "SVS" and variant_label != "main (remainder minus MEGASTROKE, FinnGen 0.15, UKB)": continue
        ds = datasets(sub, fF, sub_ukb); names = list(ds)
        est = {}
        for nm, (b, se, ld) in ds.items():
            r = LD[ld]; e = {}
            for k, mk in ((1, "rs2672592"), (2, "rs72631113")):
                bj, sej = joint(b, se, r, k)
                e[f"joint_{mk}_lead"] = (bj[0], sej[0]); e[f"joint_{mk}_side"] = (bj[1], sej[1])
            for i, s_ in enumerate(S): e[f"marginal_{s_}"] = (b[i], se[i])
            est[nm] = e
            if variant_label.startswith("main"):
                for kk, (bb, ss) in e.items():
                    dsrows.append(dict(subtype=sub, dataset=nm, ld=ld, quantity=kk, beta=bb, se=ss, z=bb / ss, p=p2(bb / ss)))
        for q in est[names[0]]:
            y = [est[n][q][0] for n in names]; s_ = [est[n][q][1] for n in names]
            o = re_meta(y, s_)
            row = dict(subtype=sub, analysis=variant_label, quantity=q, k=len(names), datasets=";".join(names),
                       n_negative=int(sum(np.array(y) < 0)), n_positive=int(sum(np.array(y) > 0)),
                       n_nominal_negative=int(sum((np.array(y) < 0) & (p2(np.array(y) / np.array(s_)) < 0.05))))
            row.update(flat(o)); rows.append(row)
            if variant_label.startswith("main") and len(names) > 2:
                for drop in names:
                    keep = [n for n in names if n != drop]
                    o2 = re_meta([est[n][q][0] for n in keep], [est[n][q][1] for n in keep])
                    lr = dict(subtype=sub, quantity=q, dropped=drop, k=len(keep)); lr.update(flat(o2)); loo.append(lr)

# East Asian GIGASTROKE IS: shown separately (not pooled)
bE, sE = get("GIGA_EAS_IS")
for k, mk in ((1, "rs2672592"), (2, "rs72631113")):
    bj, sej = joint(bE, sE, LD["KG_EAS"], k)
    for lab, i in (("lead", 0), ("side", 1)):
        dsrows.append(dict(subtype="IS", dataset="GIGASTROKE_EAS (not pooled)", ld="KG_EAS", quantity=f"joint_{mk}_{lab}",
                           beta=bj[i], se=sej[i], z=bj[i] / sej[i], p=p2(bj[i] / sej[i])))
for i, s_ in enumerate(S):
    dsrows.append(dict(subtype="IS", dataset="GIGASTROKE_EAS (not pooled)", ld="KG_EAS", quantity=f"marginal_{s_}", beta=bE[i], se=sE[i],
                       z=bE[i] / sE[i], p=p2(bE[i] / sE[i])))

re_df = pd.DataFrame(rows); re_df.to_csv(f"{O}/10q26_stroke_random_effects.tsv", sep="\t", index=False)
pd.DataFrame(loo).to_csv(f"{O}/10q26_stroke_loo.tsv", sep="\t", index=False)
pd.DataFrame(dsrows).to_csv(f"{O}/10q26_stroke_datasets.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(pd.DataFrame(dsrows).round(4).to_string())
print(re_df[re_df.analysis.str.startswith("main")][["subtype", "quantity", "k", "n_negative", "FE_est", "FE_p", "DL_est", "DL_lo", "DL_hi", "DL_p",
      "DL_tau2", "REML_est", "REML_p", "REML_pi_lo", "REML_pi_hi", "Q_p", "I2"]].round(4).to_string())
print(pd.DataFrame(loo)[["subtype", "quantity", "dropped", "DL_est", "DL_p", "REML_est", "REML_p", "Q_p"]].round(4).to_string())
print("F_UKB", F_UKB, "LD", LD)

# ================================================================ B. equivalence without pre-set bounds
mr = pd.read_csv(f"{P}/results/revision/mr/mr_revised_all.tsv", sep="\t")
mr = mr[(mr.method == "IVW") & (mr.scheme == "A_dist") & mr.exclusion.isin(["main", "no_10q26", "no_major3", "no_all"])].copy()
mr["binary"] = mr.binary.astype(str).str.upper() == "TRUE"
z90 = stats.norm.isf(0.05)
def tost(b, se, d): return np.maximum(stats.norm.sf((b + d) / se), stats.norm.cdf((b - d) / se))
bf, ce, cv = [], [], []
PREV = [0.03, 0.05, 0.07, 0.10, 0.15, 0.20]
for _, r in mr.iterrows():
    b, se = r.b, r.se; dmin = abs(b) + z90 * se          # TOST P = 0.05 exactly at this symmetric bound
    o = dict(exposure=r.exposure, outcome=r.outcome, exclusion=r.exclusion, nsnp=r.nsnp, b=b, se=se, p=r.p, binary=r.binary,
             est=np.exp(b) if r.binary else b, lo95=np.exp(b - 1.96 * se) if r.binary else b - 1.96 * se,
             hi95=np.exp(b + 1.96 * se) if r.binary else b + 1.96 * se,
             lo90=np.exp(b - z90 * se) if r.binary else b - z90 * se, hi90=np.exp(b + z90 * se) if r.binary else b + z90 * se,
             delta_min_logscale=dmin, equivalence_bound_min=(np.exp(dmin) if r.binary else dmin),
             tost_p_at_5pct=tost(b, se, np.log(1.05) if r.binary else 0.05))
    # largest causal OR of manifest AMD excluded (upper 95% one-sided limit, translated)
    if r.binary:
        up = b + z90 * se
        for p in PREV:
            o[f"max_causal_OR_bin_excluded_above_p{p}"] = np.exp(up / (p * (1 - p)))
    bf.append(o)
    if r.binary:
        for orx, lab in ((2.08, "NHANES fully adjusted OR 2.08"), (1.20, "lower 95% limit of NHANES OR (1.20)")):
            for p in PREV:
                bexp = np.log(orx) * p * (1 - p); zz = (b - bexp) / se
                ce.append(dict(exposure=r.exposure, outcome=r.outcome, exclusion=r.exclusion, observational_OR=orx, label=lab, prevalence_p=p,
                               expected_MR_OR=np.exp(bexp), observed_MR_OR=np.exp(b), lo95=np.exp(b - 1.96 * se), hi95=np.exp(b + 1.96 * se),
                               z=zz, p_one_sided_less=stats.norm.cdf(zz), rejects_causal_expectation=stats.norm.cdf(zz) < 0.05))
    grid = np.round(np.arange(1.01, 1.1501, 0.005), 3) if r.binary else np.round(np.arange(0.01, 0.1501, 0.005), 3)
    for g in grid:
        d = np.log(g) if r.binary else g
        cv.append(dict(exposure=r.exposure, outcome=r.outcome, exclusion=r.exclusion, bound=g, bound_scale="OR per unit log-odds" if r.binary else "SD",
                       tost_p=tost(b, se, d)))
bf = pd.DataFrame(bf); bf.to_csv(f"{O}/equivalence_boundfree.tsv", sep="\t", index=False)
ce = pd.DataFrame(ce); ce.to_csv(f"{O}/causal_expectation_test.tsv", sep="\t", index=False)
pd.DataFrame(cv).to_csv(f"{O}/tost_curve.tsv", sep="\t", index=False)
k = bf[bf.outcome.isin(["AS", "IS", "WMH", "MEGA_AS", "MEGA_IS"])]
print(k[["exposure", "outcome", "exclusion", "est", "lo90", "hi90", "equivalence_bound_min", "max_causal_OR_bin_excluded_above_p0.05",
         "max_causal_OR_bin_excluded_above_p0.1"]].round(4).to_string())
c = ce[(ce.observational_OR == 2.08) & ce.outcome.isin(["AS", "IS", "MEGA_AS", "MEGA_IS"])]
print(c.pivot_table(index=["exposure", "outcome", "exclusion"], columns="prevalence_p", values="p_one_sided_less").round(4).to_string())
c = ce[(ce.observational_OR == 1.20) & ce.outcome.isin(["AS", "IS"])]
print(c.pivot_table(index=["exposure", "outcome", "exclusion"], columns="prevalence_p", values="p_one_sided_less").round(4).to_string())
