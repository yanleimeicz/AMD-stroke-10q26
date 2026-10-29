"""Round 6, items A(3), B, C, E, F.
A(3) 10q26 stroke: GLS pooling of two-variant joint effects (lead rs3750848 | side rs2672592) and their contrast across
     MEGASTROKE, GIGASTROKE remainder (GIGA - MEGA - FinnGen - UKB, proxy error added to the variance), FinnGen I9_STR and
     UK Biobank I63, with between-dataset sampling correlations from null-SNP z correlations (minus the polygenic baseline
     estimated from the non-overlapping FinnGen-UKB pair); fixed-effect GLS, random-effects GLS (REML tau2) and a
     Hartung-Knapp-type t(k-1) variance; strictly non-overlapping subsets.
B    in-sample-LD-only pool: FinnGen I9_STR (SISu LD) + UK Biobank I63 (UKB LD).
C    CAD: UK Biobank 411.4 (in-sample LD), FinnGen R12 I9_CORATHER / I9_CHD (SISu LD), CARDIoGRAMplusC4D 2015 (no UKB;
     1000G EUR LD); pooled non-overlapping; Aragam 2022 shown as consistency (includes UKB).
E    phenotype-matched MR table and benchmark differences (expected - observed) with 95% intervals by simulation.
F    heterogeneity main table."""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import json, numpy as np, pandas as pd
from scipy import stats, optimize
P = PROJ_ROOT; O = f"{P}/results/revision6"
R10 = f"{P}/results/revision/10q26"; R2 = f"{P}/results/revision2"; R3 = f"{P}/results/revision3/panukb"
S = ["rs3750848", "rs2672592", "rs72631113"]; RISK = {"rs3750848": "G", "rs2672592": "T", "rs72631113": "T"}
p2 = lambda z: 2 * stats.norm.sf(abs(z))
kv = pd.concat([pd.read_csv(f"{R10}/key_variants_all_datasets.tsv", sep="\t"), pd.read_csv(f"{R2}/finngen_extra_endpoints_10q26.tsv", sep="\t")], ignore_index=True)
uk = pd.read_csv(f"{R3}/panukb_10q26_variants.tsv", sep="\t")
def get(t, src="kv"):
    d = (kv if src == "kv" else uk); d = d[d.trait == t].drop_duplicates("rsid").set_index("rsid")
    return d.loc[S[:2], "beta_risk"].values.astype(float), d.loc[S[:2], "se"].values.astype(float)
fg = json.load(open(f"{R10}/ld_finngen_sisu42_pairs.json"))
kg = pd.read_csv(f"{R10}/ld_1000G_EUR_r_dprime.tsv", sep=r"\s+")
ukl = pd.read_csv(f"{R3}/panukb_ld_keypairs.tsv", sep="\t")
def fg_r(a, b): k = f"{a}|{b}" if f"{a}|{b}" in fg else f"{b}|{a}"; return np.sqrt(fg[k]["r2"])
def kg_r(a, b): x = kg[((kg.SNP_A == a) & (kg.SNP_B == b)) | ((kg.SNP_A == b) & (kg.SNP_B == a))]; return float(x.R.iloc[0])
def uk_r(a, b): x = ukl[((ukl.snp_a == a) & (ukl.snp_b == b)) | ((ukl.snp_a == b) & (ukl.snp_b == a))]; return float(x.r_UKB_British.iloc[0])
R = {"KG_EUR": kg_r(S[0], S[1]), "FinnGen_SISu": fg_r(S[0], S[1]), "UKB_British": uk_r(S[0], S[1])}
def neff(a, b): return 4 / (1 / a + 1 / b)
F_UKB = {"AS": neff(2932, 369419) / neff(73652, 1234808), "IS": neff(2211, 369419) / neff(62100, 1234808)}

def remainder(sub, fF=0.15):
    bG, sG = get(f"GIGA_EUR_{sub}"); bM, sM = get(f"MEGA_EUR_{sub}"); bF, sF = get("FinnGen_I9_STR")
    bU, sU = get("IS_I63" if sub == "IS" else "AS_sr6150_3", "uk")
    IG, IM = 1 / sG ** 2, 1 / sM ** 2; IF, IU = fF * IG, F_UKB[sub] * IG; IR = IG - IM - IF - IU
    b = (IG * bG - IM * bM - IF * bF - IU * bU) / IR
    v = 1 / IR + (IF / IR) ** 2 * np.maximum(0, 1 / IF - sF ** 2) + (IU / IR) ** 2 * np.maximum(0, 1 / IU - sU ** 2)
    return b, np.sqrt(v), 1 / np.sqrt(IR)

def joint(b, se, r):
    z = b / se; V = np.linalg.inv(np.array([[1, r], [r, 1]])); bs = V @ z
    bj = bs * se; sj = np.sqrt(np.diag(V)) * se; cov = V[0, 1] * se[0] * se[1]
    return dict(b_lead=bj[0], s_lead=sj[0], b_side=bj[1], s_side=sj[1], cov=cov, b_con=bj[1] - bj[0],
                s_con=np.sqrt(sj[0] ** 2 + sj[1] ** 2 - 2 * cov))

def stroke_datasets(sub):
    d = {}
    d["MEGASTROKE"] = (*get(f"MEGA_EUR_{sub}"), "KG_EUR", "MEGASTROKE EUR (Malik 2018; no UKB, no FinnGen)", {"IS": "34,217 / 406,111", "AS": "40,585 / 406,111"}[sub], "IS" if sub == "IS" else "any stroke", "clinic- and population-based; TOAST subtyping")
    bR, sR, sR0 = remainder(sub)
    d["GIGASTROKE_remainder"] = (bR, sR, "KG_EUR", "GIGASTROKE EUR minus MEGASTROKE, FinnGen and UKB (approximate)", {"IS": "~16,200 / ~265,000 (approx.: GIGASTROKE EUR minus MEGASTROKE, FinnGen/FINRISK and UKB rows)", "AS": "~19,100 / ~265,000 (approx.)"}[sub], "IS" if sub == "IS" else "any stroke", "biobank and population cohorts (HUNT, EstBB, Geisinger, CCHS, ...)")
    d["FinnGen_I9_STR"] = (*get("FinnGen_I9_STR"), "FinnGen_SISu", "FinnGen R12 I9_STR (stroke excl. SAH, incl. ICH)", "34,110 / 450,023", "stroke excl. SAH (includes ICH)", "hospital discharge and death registers (ICD)")
    d["UKB"] = (*get("IS_I63" if sub == "IS" else "AS_sr6150_3", "uk"), "UKB_British",
                "UK Biobank (Pan-UKB EUR)", "4,275 / 416,256" if sub == "IS" else "6,454 / 413,170", "ICD-10 I63" if sub == "IS" else "self-reported stroke", "hospital ICD-10 / self-report")
    out = {}
    for k, (b, se, ld, label, n, endpoint, ascert) in d.items():
        j = joint(b, se, R[ld]); j.update(ld=ld, r=R[ld], label=label, n_cases_controls=n, endpoint=endpoint, ascertainment=ascert,
                                       ancestry="Finnish" if k == "FinnGen_I9_STR" else "European"); out[k] = j
    return out

nz = pd.read_csv(f"{O}/null_z_correlations.tsv", sep="\t")
def rho_raw(a, b):
    x = nz[((nz.a == a) & (nz.b == b)) | ((nz.a == b) & (nz.b == a))]; return float(x.r_null_z.iloc[0]) if len(x) else 0.0
BASE = rho_raw("FinnGen_I9_STR", "UKB_I63")        # polygenic baseline: non-overlapping pair
NZ = {"IS": {"MEGASTROKE": "MEGA_IS", "GIGASTROKE_remainder": "GIGA_remainder2_IS", "FinnGen_I9_STR": "FinnGen_I9_STR", "UKB": "UKB_I63"},
      "AS": {"MEGASTROKE": "MEGA_AS", "GIGASTROKE_remainder": "GIGA_remainder2_AS", "FinnGen_I9_STR": "FinnGen_I9_STR", "UKB": "UKB_I63"}}
def rho(sub, a, b, mode):
    r = rho_raw(NZ[sub][a], NZ[sub][b])
    return r if mode == "raw" else max(0.0, r - BASE)

def gls(y, s, C):
    y, s = np.asarray(y, float), np.asarray(s, float); k = len(y); X = np.ones(k)
    Vs = np.outer(s, s) * C
    def fit(t2):
        W = np.linalg.inv(Vs + t2 * np.eye(k)); var = 1 / (X @ W @ X); mu = var * (X @ W @ y); res = y - mu
        return mu, var, W, res
    def nll(t2):
        V = Vs + t2 * np.eye(k); W = np.linalg.inv(V); var = 1 / (X @ W @ X); mu = var * (X @ W @ y); res = y - mu
        return 0.5 * (np.linalg.slogdet(V)[1] + np.log(X @ W @ X) + res @ W @ res)
    mu0, var0, W0, res0 = fit(0.0); Q = res0 @ W0 @ res0
    t2 = optimize.minimize_scalar(nll, bounds=(0, max(1e-10, 10 * np.var(y) + 1e-6)), method="bounded").x if k > 1 else 0.0
    mu, var, W, res = fit(t2); q = (res @ W @ res) / (k - 1) if k > 1 else np.nan
    se_hk = np.sqrt(q * var); tq = stats.t.ppf(0.975, k - 1) if k > 1 else np.nan
    pi = stats.t.ppf(0.975, k - 2) * np.sqrt(t2 + var) if k > 2 else np.nan
    return dict(k=k, FE_est=mu0, FE_se=np.sqrt(var0), FE_p=p2(mu0 / np.sqrt(var0)), RE_est=mu, RE_se=np.sqrt(var), RE_p=p2(mu / np.sqrt(var)),
                RE_lo=mu - 1.96 * np.sqrt(var), RE_hi=mu + 1.96 * np.sqrt(var), tau2=t2, Q=Q, Q_p=stats.chi2.sf(Q, k - 1) if k > 1 else np.nan,
                I2=max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0, HK_se=se_hk, HK_p=2 * stats.t.sf(abs(mu / se_hk), k - 1) if k > 1 else np.nan,
                HK_lo=mu - tq * se_hk, HK_hi=mu + tq * se_hk, HKmod_p=2 * stats.t.sf(abs(mu / max(se_hk, np.sqrt(var))), k - 1) if k > 1 else np.nan,
                PI_lo=mu - pi if k > 2 else np.nan, PI_hi=mu + pi if k > 2 else np.nan)

SETS = {"main: MEGASTROKE + GIGASTROKE remainder + FinnGen + UKB (overlap-corrected covariance)": ["MEGASTROKE", "GIGASTROKE_remainder", "FinnGen_I9_STR", "UKB"],
        "strict non-overlapping 1: MEGASTROKE + FinnGen + UKB": ["MEGASTROKE", "FinnGen_I9_STR", "UKB"],
        "strict non-overlapping 2: GIGASTROKE remainder + FinnGen + UKB": ["GIGASTROKE_remainder", "FinnGen_I9_STR", "UKB"],
        "strict non-overlapping 3: MEGASTROKE + FinnGen": ["MEGASTROKE", "FinnGen_I9_STR"],
        "in-sample LD only (item B): FinnGen + UKB": ["FinnGen_I9_STR", "UKB"]}
rows, het = [], []
for sub in ("IS", "AS"):
    ds = stroke_datasets(sub)
    for nm, j in ds.items():
        het.append(dict(subtype=sub, dataset=nm, **{k: j[k] for k in ("label", "endpoint", "ancestry", "n_cases_controls", "ascertainment", "ld", "r")},
                        b_lead=j["b_lead"], s_lead=j["s_lead"], p_lead=p2(j["b_lead"] / j["s_lead"]), b_side=j["b_side"], s_side=j["s_side"],
                        p_side=p2(j["b_side"] / j["s_side"]), b_contrast=j["b_con"], s_contrast=j["s_con"], p_contrast=p2(j["b_con"] / j["s_con"])))
    for mode in ("baseline-adjusted", "raw"):
        for setname, mem in SETS.items():
            if mode == "raw" and not setname.startswith("main"): continue
            C = np.array([[1.0 if a == b else rho(sub, a, b, mode) for b in mem] for a in mem])
            for q, (bc, sc) in {"side_rs2672592_given_lead": ("b_side", "s_side"), "lead_rs3750848_given_side": ("b_lead", "s_lead"),
                                "contrast_side_minus_lead": ("b_con", "s_con")}.items():
                g = gls([ds[m][bc] for m in mem], [ds[m][sc] for m in mem], C)
                g.update(subtype=sub, set=setname, covariance=mode, quantity=q, datasets=";".join(mem),
                         max_offdiag_rho=float(np.max(C - np.eye(len(mem)))) if len(mem) > 1 else 0.0); rows.append(g)
res = pd.DataFrame(rows)
cols = ["subtype", "set", "covariance", "quantity", "k", "datasets", "max_offdiag_rho", "FE_est", "FE_se", "FE_p", "RE_est", "RE_se", "RE_p", "RE_lo", "RE_hi",
        "tau2", "Q", "Q_p", "I2", "HK_se", "HK_p", "HK_lo", "HK_hi", "HKmod_p", "PI_lo", "PI_hi"]
res[cols].to_csv(f"{O}/10q26_stroke_gls_overlap.tsv", sep="\t", index=False)
res[res.set.str.startswith("in-sample")][cols].to_csv(f"{O}/insample_only_pool.tsv", sep="\t", index=False)
H = pd.DataFrame(het)
for sub in ("IS", "AS"):
    g = res[(res.subtype == sub) & res.set.str.startswith("main") & (res.covariance == "baseline-adjusted")]
    for q in g.quantity:
        x = g[g.quantity == q].iloc[0]
        H = pd.concat([H, pd.DataFrame([dict(subtype=sub, dataset=f"POOLED main ({q})", label="GLS random effects; HK t(3)",
                                            b_contrast=x.RE_est if q.startswith("contrast") else np.nan,
                                            b_side=x.RE_est if q.startswith("side") else np.nan, b_lead=x.RE_est if q.startswith("lead") else np.nan,
                                            Q=x.Q, Q_p=x.Q_p, I2=x.I2, tau2=x.tau2, HK_p=x.HK_p, HK_lo=x.HK_lo, HK_hi=x.HK_hi, PI_lo=x.PI_lo, PI_hi=x.PI_hi)])])
H.to_csv(f"{O}/heterogeneity_main_table.tsv", sep="\t", index=False)

# ---------------- C. CAD
cad = []
fgc = pd.read_csv(f"{P}/derived/revision6/finngen_cad_10q26_region.tsv", sep="\t")
for e in ("I9_CORATHER", "I9_CHD"):
    x = fgc[(fgc.endpoint == e) & fgc.rsids.isin(S[:2])].set_index("rsids")
    b = np.array([x.loc[s_, "beta"] * (1 if x.loc[s_, "alt"] == RISK[s_] else -1) for s_ in S[:2]]); se = x.loc[S[:2], "sebeta"].values
    j = joint(b, se, R["FinnGen_SISu"]); j.update(dataset=f"FinnGen R12 {e}", ld="FinnGen_SISu (in-sample)", independent_of_UKB=True); cad.append(j)
c4 = pd.read_csv(f"{P}/derived/revision6/c4d_10q26.tsv", sep="\t").set_index("markername")
# allele alignment by frequency: the C4D 2015 file labels rs2672592 effect allele "T" with frequency 0.63, whereas the T (risk)
# allele frequency is 0.36-0.39 in UKB/FinnGen and rs2672592-T is in r = 0.95 with rs72631113-T (C4D effect C, beta > 0);
# the reported beta therefore refers to the allele with frequency ~0.63 (G).  Align to the risk allele using frequency.
FR = {"rs3750848": 0.2144, "rs2672592": 0.3578}   # UKB EUR control frequency of the risk allele
b = np.array([c4.loc[s_, "beta"] * (1 if abs(c4.loc[s_, "effect_allele_freq"] - FR[s_]) < abs(c4.loc[s_, "effect_allele_freq"] - (1 - FR[s_])) else -1) for s_ in S[:2]])
se = c4.loc[S[:2], "se_dgc"].values
j = joint(b, se, R["KG_EUR"]); j.update(dataset="CARDIoGRAMplusC4D 2015 (Nikpay; 60,801 cases)", ld="1000G EUR", independent_of_UKB=True); cad.append(j)
bU, sU = get("CAD_phe411.4", "uk")
j = joint(bU, sU, R["UKB_British"]); j.update(dataset="UK Biobank coronary atherosclerosis (411.4)", ld="UKB (in-sample)", independent_of_UKB=False); cad.append(j)
cad = pd.DataFrame(cad)
for c in ("lead", "side", "con"):
    cad[f"p_{c}"] = p2(cad[f"b_{c}"] / cad[f"s_{c}"])
pool_rows = []
for nm, mem in {"independent replication: FinnGen I9_CORATHER + CARDIoGRAMplusC4D 2015": ["FinnGen R12 I9_CORATHER", "CARDIoGRAMplusC4D 2015 (Nikpay; 60,801 cases)"],
                "all non-overlapping: UKB + FinnGen I9_CORATHER + CARDIoGRAMplusC4D 2015": ["UK Biobank coronary atherosclerosis (411.4)", "FinnGen R12 I9_CORATHER", "CARDIoGRAMplusC4D 2015 (Nikpay; 60,801 cases)"]}.items():
    x = cad.set_index("dataset").loc[mem]
    for q, (bc, sc) in {"side": ("b_side", "s_side"), "lead": ("b_lead", "s_lead"), "contrast": ("b_con", "s_con")}.items():
        g = gls(x[bc].values, x[sc].values, np.eye(len(mem))); g.update(set=nm, quantity=q); pool_rows.append(g)
cadp = pd.DataFrame(pool_rows)
with open(f"{O}/cad_independent_validation.tsv", "w", encoding="utf-8") as f:
    cad.to_csv(f, sep="\t", index=False); f.write("\n"); cadp.to_csv(f, sep="\t", index=False)

# ---------------- E. phenotype-matched MR and benchmark differences
mr = pd.read_csv(f"{P}/results/revision/mr/mr_revised_all.tsv", sep="\t")
mr = mr[(((mr.scheme == "A_dist") & (mr.method == "IVW")) | ((mr.scheme == "B_52") & (mr.method == "IVW_correl"))) & mr.outcome.isin(["AS", "IS", "SVS"]) & mr.exclusion.isin(["main", "no_10q26", "no_all"]) & mr.scheme.isin(["A_dist", "B_52"])]
lab = {"AMD_FG": "registry AMD (FinnGen, any)", "DRYAMD_FG": "registry dry AMD (FinnGen)", "WETAMD_FG": "registry wet AMD (FinnGen)", "AMD_IAMDGC": "advanced AMD (IAMDGC)"}
E = mr[mr.exposure.isin(lab)][["exposure", "scheme", "method", "outcome", "exclusion", "nsnp", "b", "se", "OR", "lo", "hi", "p"]].copy()
E["phenotype"] = E.exposure.map(lab)
ea = pd.read_csv(f"{P}/results/revision5/earlyamd_exposure_mr.tsv", sep="\t")
ea = ea[ea.outcome.isin(["AS", "IS", "SVS"])].copy(); ea["scheme"] = "early-AMD instruments (P < 5e-8)"; ea["phenotype"] = "early AMD (Winkler 2020)"
E = pd.concat([E, ea[["exposure", "scheme", "outcome", "exclusion", "nsnp", "b", "se", "OR", "lo", "hi", "p", "phenotype"]]], ignore_index=True)
E.sort_values(["phenotype", "scheme", "outcome", "exclusion"]).to_csv(f"{O}/phenotype_matched_mr.tsv", sep="\t", index=False)
kap = pd.read_csv(f"{P}/results/revision5/kappa_earlyamd_scale.tsv", sep="\t")
rng = np.random.default_rng(7); ND = 200000
lo_obs, se_obs = np.log(2.08), (np.log(3.61) - np.log(1.20)) / (2 * 1.96)
brows = []
for _, x in E[(E.scheme.isin(["A_dist", "early-AMD instruments (P < 5e-8)"])) & E.outcome.isin(["AS", "IS"]) & E.exclusion.isin(["main", "no_10q26", "no_all"])].iterrows():
    if x.exposure == "EARLYAMD_W20": k_, ks = 1.0, 0.0
    else:
        kk = kap[(kap.exposure == x.exposure) & (kap.set == ("main" if x.exclusion == "main" else "no_10q26" if x.exclusion == "no_10q26" else "no_all"))]
        if kk.empty: kk = kap[(kap.exposure == x.exposure) & (kap.set == "main")]
        k_, ks = float(kk.kappa.iloc[0]), float(kk.se.iloc[0])
    for p in (0.03, 0.063, 0.10):
        lor = rng.normal(lo_obs, se_obs, ND); OR_ = np.exp(lor); kd = rng.normal(k_, ks, ND) if ks > 0 else np.full(ND, k_)
        g = p * (1 - p) * (OR_ - 1) / (1 + p * (OR_ - 1)); bexp = kd * g; bmr = rng.normal(x.b, x.se, ND); d = bexp - bmr
        brows.append(dict(phenotype=x.phenotype, exposure=x.exposure, outcome=x.outcome, exclusion=x.exclusion, prevalence=p, kappa=k_, kappa_se=ks,
                          expected_OR=np.exp(np.median(bexp)), observed_OR=x.OR, observed_lo=x.lo, observed_hi=x.hi,
                          diff_ratio_OR_median=np.exp(np.median(d)), diff_ratio_OR_lo=np.exp(np.quantile(d, 0.025)), diff_ratio_OR_hi=np.exp(np.quantile(d, 0.975)),
                          frac_draws_expected_gt_observed=(d > 0).mean()))
pd.DataFrame(brows).to_csv(f"{O}/benchmark_difference_ci.tsv", sep="\t", index=False)

pd.set_option("display.width", 260); pd.set_option("display.max_columns", 40)
print("R", R, "BASE", BASE)
print(H.round(4).to_string())
print(res[["subtype", "set", "covariance", "quantity", "k", "max_offdiag_rho", "FE_est", "FE_p", "RE_est", "RE_p", "HK_p", "HK_lo", "HK_hi", "Q_p", "I2"]].round(4).to_string())
print(cad.round(4).to_string()); print(cadp[["set", "quantity", "FE_est", "FE_p", "RE_p", "HK_p", "Q_p"]].round(4).to_string())
print(pd.DataFrame(brows)[lambda d: d.prevalence == 0.063].round(3).to_string())
