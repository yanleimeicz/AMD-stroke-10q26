"""Round-2 minor revision, items 1, 2 and 5 (10q26 symmetric conditional models; FinnGen stroke heterogeneity
and power; EPVS Z/FDR reconciliation).
Inputs: results/revision/10q26 (key-variant Z, LD from FinnGen SISu v4.2 / 1000G EUR / 1000G EAS / Ensembl),
raw MEGASTROKE, FinnGen R12 (I9_STR, C_STROKE, I9_TIA) and Duperron 2023 EPVS files.
Outputs: results/revision2/10q26_symmetric_conditional.tsv, 10q26_r_sensitivity.tsv, 10q26_proportionality_projected.tsv,
         finngen_heterogeneity_power.tsv, finngen_extra_endpoints_10q26.tsv, epvs_corrected.tsv
Usage: python rev2_10q26.py"""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import gzip, json, os
import numpy as np, pandas as pd
from scipy import stats, optimize
from joblib import Parallel, delayed

P = PROJ_ROOT
R10 = f"{P}/results/revision/10q26"; O = f"{P}/results/revision2"; RAW = f"{P}/data/raw/gwas"
os.makedirs(O, exist_ok=True)
S = ["rs3750848", "rs2672592", "rs72631113"]          # AMD lead (ARMS2 haplotype); stroke-side pair
RISK = {"rs3750848": "G", "rs2672592": "T", "rs72631113": "T", "rs10490924": "T", "rs11200638": "A"}
POS38 = {"rs72631113": 122453933, "rs10490924": 122454932, "rs3750848": 122455799, "rs11200638": 122461028, "rs2672592": 122471234}
p2 = lambda z: 2 * stats.norm.sf(abs(z))

# ---------------------------------------------------------------- LD sources (risk alleles all positively coupled, D' ~ 1)
fg = json.load(open(f"{R10}/ld_finngen_sisu42_pairs.json"))
def fg_r(a, b):
    k = f"{a}|{b}" if f"{a}|{b}" in fg else f"{b}|{a}"
    return np.sqrt(fg[k]["r2"])
kg = pd.read_csv(f"{R10}/ld_1000G_EUR_r_dprime.tsv", sep=r"\s+")
def kg_r(a, b):
    x = kg[((kg.SNP_A == a) & (kg.SNP_B == b)) | ((kg.SNP_A == b) & (kg.SNP_B == a))]
    return float(x.R.iloc[0])
ens = pd.read_csv(f"{R10}/ld_ensembl_EUR_FIN_EAS.tsv", sep="\t")
def ens_r(pop, a, b):
    x = ens[(ens["pop"] == pop) & (((ens.snp1 == a) & (ens.snp2 == b)) | ((ens.snp1 == b) & (ens.snp2 == a)))]
    return np.sqrt(float(x.r2.iloc[0]))
pairs = [(S[0], S[1]), (S[0], S[2]), (S[1], S[2])]
LD = {
    "FinnGen_SISu": [fg_r(*p) for p in pairs],
    "KG_EUR": [kg_r(*p) for p in pairs],
    "Ensembl_EUR": [ens_r("EUR", *p) for p in pairs],
    "Ensembl_FIN": [ens_r("FIN", *p) for p in pairs],
    "KG_EAS": [0.821883203381113, 0.678589714039345, 0.541761017423735],   # plink 1000G EAS (conditional_two_snp_matched_LD.tsv)
}
def Rmat(r):
    return np.array([[1, r[0], r[1]], [r[0], 1, r[2]], [r[1], r[2], 1]])

# ---------------------------------------------------------------- key-variant Z (risk-allele coded)
kv = pd.read_csv(f"{R10}/key_variants_all_datasets.tsv", sep="\t")

def read_finngen(path, rsids):
    out = {}
    with gzip.open(path, "rt") as f:
        h = f.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
        for line in f:
            if not line.startswith("10\t"): continue
            x = line.rstrip("\n").split("\t")
            pos = int(x[ix["pos"]])
            for r in rsids:
                if pos == POS38[r] and r in x[ix["rsids"]].split(","):
                    b = float(x[ix["beta"]]); se = float(x[ix["sebeta"]])
                    if x[ix["alt"]] != RISK[r]: b = -b
                    out[r] = dict(beta=b, se=se, p=float(x[ix["pval"]]), eaf=float(x[ix["af_alt"]]))
            if pos > 122500000: break
    return out

def read_gwascat_h(path, rsids):
    out = {}
    with gzip.open(path, "rt") as f:
        h = f.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
        for line in f:
            x = line.rstrip("\n").split("\t")
            r = x[ix["hm_rsid"]]
            if r in rsids:
                b = float(x[ix["hm_beta"]]); se = float(x[ix["standard_error"]])
                if x[ix["hm_effect_allele"]].upper() != RISK[r]: b = -b
                out[r] = dict(beta=b, se=se, p=float(x[ix["p_value"]]))
                if len(out) == len(rsids): break
    return out

extra_files = {
    "FinnGen_C_STROKE": (f"{RAW}/finngen_R12_C_STROKE.gz", read_finngen),
    "FinnGen_I9_TIA": (f"{RAW}/finngen_R12_I9_TIA.gz", read_finngen),
    "MEGA_EUR_AS": (f"{RAW}/MEGASTROKE_AS_EUR_GCST006906.h.tsv.gz", read_gwascat_h),
    "MEGA_EUR_IS": (f"{RAW}/MEGASTROKE_IS_EUR_GCST006908.h.tsv.gz", read_gwascat_h),
    "MEGA_EUR_SVS": (f"{RAW}/MEGASTROKE_SVS_EUR_GCST006909.h.tsv.gz", read_gwascat_h),
}
ext_rs = list(POS38)
avail = {k: v for k, v in extra_files.items() if os.path.exists(v[0])}
res = Parallel(n_jobs=len(avail))(delayed(fn)(path, ext_rs) for path, fn in avail.values())
ext_rows = []
for k, d in zip(avail, res):
    for r, v in d.items():
        ext_rows.append(dict(rsid=r, pos=POS38[r], risk=RISK[r], beta_risk=v["beta"], se=v["se"], z_risk=v["beta"] / v["se"], p=v["p"], trait=k))
ext = pd.DataFrame(ext_rows)
ext.to_csv(f"{O}/finngen_extra_endpoints_10q26.tsv", sep="\t", index=False)
allkv = pd.concat([kv[["rsid", "pos", "risk", "beta_risk", "se", "z_risk", "p", "trait"]], ext], ignore_index=True)

def zvec(trait):
    d = allkv[allkv.trait == trait].set_index("rsid")
    if not all(s in d.index for s in S): return None
    return np.array([d.loc[s, "z_risk"] for s in S], float)

# ---------------------------------------------------------------- item 1: symmetric conditional models
def models(z, R):
    out = {}
    for i, j in [(1, 0), (0, 1), (2, 0), (0, 2), (2, 1), (1, 2)]:
        r = R[i, j]; zc = (z[i] - r * z[j]) / np.sqrt(1 - r * r)
        out[f"z_{S[i]}|{S[j]}"] = zc; out[f"p_{S[i]}|{S[j]}"] = p2(zc)
    Ri = np.linalg.inv(R); zj = (Ri @ z) / np.sqrt(np.diag(Ri))
    for k in range(3): out[f"zjoint3_{S[k]}"] = zj[k]; out[f"pjoint3_{S[k]}"] = p2(zj[k])
    chi_full = z @ Ri @ z
    chi_side = chi_full - z[0] ** 2                                        # {rs2672592, rs72631113} | rs3750848, 2 df
    chi_lead = chi_full - z[1:] @ np.linalg.solve(R[1:, 1:], z[1:])        # rs3750848 | {rs2672592, rs72631113}, 1 df
    out.update(chi2_side_given_lead_2df=chi_side, p_side_given_lead_2df=stats.chi2.sf(chi_side, 2),
               chi2_lead_given_side_1df=chi_lead, p_lead_given_side_1df=stats.chi2.sf(chi_lead, 1))
    # ratio of joint standardized effects, side/lead, in two-variant models (delta-method SE)
    for k in (1, 2):
        R2 = R[np.ix_([0, k], [0, k])]; V = np.linalg.inv(R2); b = V @ z[[0, k]]
        ratio = b[1] / b[0]
        g = np.array([-b[1] / b[0] ** 2, 1 / b[0]]); se = np.sqrt(g @ V @ g)
        out[f"ratio_{S[k]}_to_lead"] = ratio; out[f"ratio_{S[k]}_to_lead_se"] = se
    return out

TRAITS = {  # trait -> list of LD sources (first = matched/primary)
    "FinnGen_H7_AMD": ["FinnGen_SISu", "Ensembl_FIN", "KG_EUR"], "FinnGen_DRY_AMD": ["FinnGen_SISu", "Ensembl_FIN", "KG_EUR"],
    "FinnGen_WET_AMD": ["FinnGen_SISu", "Ensembl_FIN", "KG_EUR"],
    "FinnGen_I9_STR": ["FinnGen_SISu", "Ensembl_FIN", "KG_EUR"], "FinnGen_C_STROKE": ["FinnGen_SISu", "Ensembl_FIN", "KG_EUR"],
    "FinnGen_I9_TIA": ["FinnGen_SISu", "Ensembl_FIN", "KG_EUR"],
    "GIGA_EUR_AS": ["KG_EUR", "Ensembl_EUR", "FinnGen_SISu"], "GIGA_EUR_IS": ["KG_EUR", "Ensembl_EUR", "FinnGen_SISu"],
    "GIGA_EUR_CES": ["KG_EUR", "Ensembl_EUR", "FinnGen_SISu"], "GIGA_EUR_LAS": ["KG_EUR", "Ensembl_EUR", "FinnGen_SISu"],
    "GIGA_EUR_SVS": ["KG_EUR", "Ensembl_EUR", "FinnGen_SISu"], "WMH": ["KG_EUR", "Ensembl_EUR", "FinnGen_SISu"],
    "MEGA_EUR_AS": ["KG_EUR", "Ensembl_EUR", "FinnGen_SISu"], "MEGA_EUR_IS": ["KG_EUR", "Ensembl_EUR", "FinnGen_SISu"],
    "MEGA_EUR_SVS": ["KG_EUR", "Ensembl_EUR", "FinnGen_SISu"],
    "GIGA_EAS_AS": ["KG_EAS"], "GIGA_EAS_IS": ["KG_EAS"], "GIGA_EAS_SVS": ["KG_EAS"],
}
rows = []
for t, lds in TRAITS.items():
    z = zvec(t)
    if z is None: continue
    for ld in lds:
        r = LD[ld]; o = dict(trait=t, ld_source=ld, ld_matched=(ld == lds[0]), r_lead_2672592=r[0], r_lead_72631113=r[1], r_2672592_72631113=r[2],
                            z_rs3750848=z[0], z_rs2672592=z[1], z_rs72631113=z[2])
        o.update(models(z, Rmat(r))); rows.append(o)
sym = pd.DataFrame(rows)
sym.to_csv(f"{O}/10q26_symmetric_conditional.tsv", sep="\t", index=False)

# r-sensitivity: shift both lead-side correlations jointly (and r(lead, rs2672592) alone) by delta
sens = []
for t, lds in TRAITS.items():
    z = zvec(t)
    if z is None: continue
    r0 = np.array(LD[lds[0]])
    for mode in ("both_lead_side", "lead_2672592_only"):
        for d in np.round(np.arange(-0.10, 0.1001, 0.01), 3):
            r = r0.copy(); r[0] += d
            if mode == "both_lead_side": r[1] += d
            R = Rmat(r)
            if np.linalg.eigvalsh(R).min() < 1e-6 or r.max() >= 0.999: continue
            m = models(z, R)
            sens.append(dict(trait=t, ld_source=lds[0], mode=mode, delta=d, r_lead_2672592=r[0], r_lead_72631113=r[1],
                             z_rs2672592_given_lead=m["z_rs2672592|rs3750848"], z_rs72631113_given_lead=m["z_rs72631113|rs3750848"],
                             z_lead_given_rs2672592=m["z_rs3750848|rs2672592"],
                             zjoint3_rs2672592=m["zjoint3_rs2672592"], zjoint3_rs72631113=m["zjoint3_rs72631113"], zjoint3_lead=m["zjoint3_rs3750848"],
                             p_side_given_lead_2df=m["p_side_given_lead_2df"], p_lead_given_side_1df=m["p_lead_given_side_1df"],
                             ratio_2672592_to_lead=m["ratio_rs2672592_to_lead"]))
sens = pd.DataFrame(sens)
sens.to_csv(f"{O}/10q26_r_sensitivity.tsv", sep="\t", index=False)

# proportionality test, projected: AMD joint effects (AMD LD) -> expected stroke Z under stroke LD; chi2 df = 2
def prop_test(zA, RA, zS, RS, ridge=0.0):
    RS = RS + ridge * np.eye(3); RA = RA + ridge * np.eye(3)
    m = RS @ np.linalg.solve(RA, zA); W = np.linalg.inv(RS)
    k = (m @ W @ zS) / (m @ W @ m); e = zS - k * m; chi = e @ W @ e
    return k, chi, stats.chi2.sf(chi, 2)
prop = []
for a in ["FinnGen_H7_AMD", "FinnGen_DRY_AMD", "FinnGen_WET_AMD"]:
    zA = zvec(a)
    for s, lds in TRAITS.items():
        if "AMD" in s: continue
        zS = zvec(s)
        if zS is None: continue
        for lds_s in lds[:1] + (["FinnGen_SISu"] if lds[0] != "FinnGen_SISu" and "EAS" not in s else []):
            for ridge in (0, 0.02, 0.05):
                k, chi, p = prop_test(zA, Rmat(LD["FinnGen_SISu"]), zS, Rmat(LD[lds_s]), ridge)
                W = np.linalg.inv(Rmat(LD[lds_s]) + ridge * np.eye(3))   # raw version (round-1 implementation): AMD Z used directly
                kr = (zA @ W @ zS) / (zA @ W @ zA); chir = (zS - kr * zA) @ W @ (zS - kr * zA)
                prop.append(dict(amd=a, amd_ld="FinnGen_SISu", stroke=s, stroke_ld=lds_s, ridge=ridge, k_projected=k, chi2_projected=chi,
                                 p_projected=p, k_raw=kr, chi2_raw=chir, p_raw=stats.chi2.sf(chir, 2),
                                 note=("EAS stroke vs Finnish AMD: different populations; invalid as a test of a shared variant" if "EAS" in s else "")))
pd.DataFrame(prop).to_csv(f"{O}/10q26_proportionality_projected.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- item 2: FinnGen heterogeneity, power, winner's curse
C_GWS = stats.norm.isf(2.5e-8)
def wc_mle(z, c=C_GWS):
    """Conditional MLE of the mean Z given selection |Z| > c (Ghosh, Zou & Wright 2008)."""
    nll = lambda mu: -(stats.norm.logpdf(z - mu) - np.log(stats.norm.cdf(-c - mu) + stats.norm.sf(c - mu)))
    lo, hi = (min(0, z) - 1, max(0, z) + 1)
    return optimize.minimize_scalar(nll, bounds=(lo, hi), method="bounded").x
def wc_mean(z, c=C_GWS):
    """Approximately unbiased average of naive and conditional MLE (Zhong & Prentice 2008, MSE-weighted proxy)."""
    return (z + wc_mle(z, c)) / 2

# GIGASTROKE EUR: 73,652 AS / 62,100 IS cases, 1,234,808 controls (Mishra 2022). FinnGen share of GIGASTROKE EUR
# (FINNGEN 9,803 cases + 174,426 controls; FINRISK 1,193 + 20,250; Mishra 2022 Suppl. Table 1)
neff = lambda ca, co: 4 / (1 / ca + 1 / co)
f_AS = (neff(9803, 174426) + neff(1193, 20250)) / neff(73652, 1234808)
f_IS = (neff(9803 * 62100 / 73652, 174426) + neff(1193 * 62100 / 73652, 20250)) / neff(62100, 1234808)   # IS cases pro rata
het = []
def g(trait, r):
    d = allkv[(allkv.trait == trait) & (allkv.rsid == r)]
    return (float(d.beta_risk.iloc[0]), float(d.se.iloc[0])) if len(d) else (np.nan, np.nan)
for fgt in ["FinnGen_I9_STR", "FinnGen_C_STROKE", "FinnGen_I9_TIA"]:
    for ref, f0 in [("GIGA_EUR_AS", f_AS), ("GIGA_EUR_IS", f_IS)]:
        for r in S:
            bF, sF = g(fgt, r); bG, sG = g(ref, r)
            if np.isnan(bF) or np.isnan(bG): continue
            bM, sM = g(ref.replace("GIGA", "MEGA"), r)
            for f in (0.10, f0, 0.20):
                bR = (bG - f * bF) / (1 - f); sR = sG / np.sqrt(1 - f)       # approximate GIGASTROKE-minus-FinnGen
                zd_naive = (bF - bG) / np.hypot(sF, sG); zd = (bF - bR) / np.hypot(sF, sR)
                zG = bG / sG; mle = wc_mle(zG) * sG; avg = wc_mean(zG) * sG
                ezR = bR / sF; ezM = bM / sF if not np.isnan(bM) else np.nan; ezW = mle / sF
                pw = lambda ez: stats.norm.cdf(abs(ez) - 1.96) + stats.norm.cdf(-abs(ez) - 1.96)
                het.append(dict(finngen_trait=fgt, reference=ref, rsid=r, f_finngen_share=round(f, 3), f_is_default=(f == f0),
                                beta_finngen=bF, se_finngen=sF, z_finngen=bF / sF, beta_ref=bG, se_ref=sG, z_ref=zG,
                                beta_ref_minus_finngen=bR, se_ref_minus_finngen=sR,
                                z_diff_naive=zd_naive, p_diff_naive=p2(zd_naive), z_diff_vs_ref_minus_finngen=zd, p_diff_vs_ref_minus_finngen=p2(zd),
                                beta_mega=bM, se_mega=sM, z_diff_vs_mega=(bF - bM) / np.hypot(sF, sM) if not np.isnan(bM) else np.nan,
                                p_diff_vs_mega=p2((bF - bM) / np.hypot(sF, sM)) if not np.isnan(bM) else np.nan,
                                beta_ref_wc_mle=mle, beta_ref_wc_avg=avg,
                                expected_z_finngen_refminusFG=ezR, power_refminusFG=pw(ezR), expected_z_finngen_mega=ezM,
                                power_mega=pw(ezM) if not np.isnan(ezM) else np.nan, expected_z_finngen_wc_mle=ezW, power_wc_mle=pw(ezW),
                                ratio_finngen_to_ref_minus_finngen=bF / bR))
het = pd.DataFrame(het)
het.to_csv(f"{O}/finngen_heterogeneity_power.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- item 5: EPVS (Duperron 2023): primary Z vs IVW beta/SE
fdr = pd.read_csv(f"{P}/results/revision/mr/fdr_cited_results.tsv", sep="\t")
fdr = fdr[fdr.item.str.contains("EPVS", na=False) & (fdr.layer == "Key-variant lookups")].copy()
fdr["rsid"] = fdr.item.str.split().str[0]; fdr["compartment"] = fdr.item.str.extract(r"(EPVS_\w+)$")[0]
KEEP = set(POS38) | {"rs10737680", "rs1269852", "rs3184504", "rs429358"}
def read_epvs(comp):
    fn = [f for f in os.listdir(RAW) if f.startswith(f"{comp}_Duperron2023")][0]
    out = []
    with gzip.open(f"{RAW}/{fn}", "rt") as f:
        h = f.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
        for line in f:
            x = line.rstrip("\n").split("\t")
            if x[ix["rsid"]] in KEEP:
                out.append(dict(compartment=comp, rsid=x[ix["rsid"]], effect_allele=x[ix["effect_allele"]], hm_code=int(x[ix["hm_code"]]),
                                beta_ivw=float(x[ix["beta"]]), se_ivw=float(x[ix["standard_error"]]), p_primary=float(x[ix["p_value"]]),
                                zscore_file=float(x[ix["Zscore"]]), n=int(x[ix["N"]]), n_event=int(x[ix["N_EVENT"]])))
    return out
ep = pd.DataFrame(sum(Parallel(n_jobs=3)(delayed(read_epvs)(c) for c in ["EPVS_WM", "EPVS_BG", "EPVS_HIP"]), []))
ep["z_primary_harmonised"] = np.where(ep.hm_code == 11, -ep.zscore_file, ep.zscore_file)   # hm_code 11: alleles flipped, Zscore left unflipped
ep["z_ivw_beta_over_se"] = ep.beta_ivw / ep.se_ivw
ep["z_from_p"] = np.sign(ep.z_primary_harmonised) * stats.norm.isf(ep.p_primary / 2)
ep["sign_check_beta_vs_z"] = np.sign(ep.beta_ivw) == np.sign(ep.z_primary_harmonised)
flip = ep.rsid.map(RISK).notna() & (ep.effect_allele != ep.rsid.map(RISK))
for c in ["beta_ivw", "z_primary_harmonised", "z_ivw_beta_over_se", "z_from_p"]:
    ep.loc[flip, c] = -ep.loc[flip, c]
ep["allele_reported"] = np.where(ep.rsid.map(RISK).notna(), ep.rsid.map(RISK) + " (AMD-risk)", ep.effect_allele)
ep = ep.merge(fdr[["rsid", "compartment", "p", "fdr", "n_tests", "status"]].rename(columns={"p": "p_in_fdr_file", "fdr": "fdr_76tests"}),
              on=["rsid", "compartment"], how="left")
ep.sort_values(["rsid", "compartment"]).to_csv(f"{O}/epvs_corrected.tsv", sep="\t", index=False)

print("== symmetric conditional (matched LD) ==")
print(sym[sym.ld_matched][["trait", "ld_source", "p_rs2672592|rs3750848", "p_rs72631113|rs3750848", "p_side_given_lead_2df",
                           "p_lead_given_side_1df", "pjoint3_rs3750848", "pjoint3_rs2672592", "pjoint3_rs72631113",
                           "ratio_rs2672592_to_lead", "ratio_rs2672592_to_lead_se"]].to_string())
print("== f shares", round(f_AS, 3), round(f_IS, 3))
print(het[het.f_is_default][["finngen_trait", "reference", "rsid", "z_finngen", "z_ref", "beta_ref_minus_finngen", "p_diff_naive",
                             "p_diff_vs_ref_minus_finngen", "p_diff_vs_mega", "expected_z_finngen_refminusFG", "power_refminusFG",
                             "expected_z_finngen_wc_mle", "power_wc_mle"]].to_string())
print(ep.to_string())
