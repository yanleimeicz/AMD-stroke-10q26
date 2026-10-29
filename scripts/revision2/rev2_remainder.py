"""Round-2 item 2 add-on: approximate inverse-variance subtraction of MEGASTROKE and/or FinnGen from GIGASTROKE EUR at the
three 10q26 markers, then the same symmetric conditional models (1000G EUR LD). MEGASTROKE information share is taken per
variant as SE_GIGA^2 / SE_MEGA^2; FinnGen share as 0.15 of GIGASTROKE information (case/control counts, Mishra 2022).
Output: results/revision2/10q26_giga_remainder_conditional.tsv"""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import numpy as np, pandas as pd
from scipy import stats
import importlib.util, sys
O = f"{PROJ_ROOT}/results/revision2"; R10 = f"{PROJ_ROOT}/results/revision/10q26"
S = ["rs3750848", "rs2672592", "rs72631113"]
kv = pd.concat([pd.read_csv(f"{R10}/key_variants_all_datasets.tsv", sep="\t"), pd.read_csv(f"{O}/finngen_extra_endpoints_10q26.tsv", sep="\t")])
kg = pd.read_csv(f"{R10}/ld_1000G_EUR_r_dprime.tsv", sep=r"\s+")
def kg_r(a, b):
    x = kg[((kg.SNP_A == a) & (kg.SNP_B == b)) | ((kg.SNP_A == b) & (kg.SNP_B == a))]; return float(x.R.iloc[0])
r = [kg_r(S[0], S[1]), kg_r(S[0], S[2]), kg_r(S[1], S[2])]
R = np.array([[1, r[0], r[1]], [r[0], 1, r[2]], [r[1], r[2], 1]])
p2 = lambda z: 2 * stats.norm.sf(abs(z))
def get(t):
    d = kv[kv.trait == t].set_index("rsid"); return d.loc[S, "beta_risk"].values.astype(float), d.loc[S, "se"].values.astype(float)
rows = []
for sub in ("AS", "IS", "SVS"):
    bG, sG = get(f"GIGA_EUR_{sub}"); bM, sM = get(f"MEGA_EUR_{sub}"); bF, sF = get("FinnGen_I9_STR")
    IG, IM = 1 / sG**2, 1 / sM**2
    for label, use_m, fF in (("GIGA_minus_MEGA", True, 0.0), ("GIGA_minus_MEGA_minus_FinnGen", True, 0.15), ("GIGA_minus_FinnGen", False, 0.15), ("MEGA", None, 0)):
        if use_m is None:
            b, se = bM, sM; info_share = IM / IG
        else:
            IF = fF * IG
            IR = IG - (IM if use_m else 0) - IF
            b = (IG * bG - (IM * bM if use_m else 0) - IF * bF) / IR; se = 1 / np.sqrt(IR); info_share = IR / IG
        z = b / se
        Ri = np.linalg.inv(R); full = z @ Ri @ z
        zc = lambda i, j: (z[i] - R[i, j] * z[j]) / np.sqrt(1 - R[i, j] ** 2)
        side = full - z[0] ** 2; lead = full - z[1:] @ np.linalg.solve(R[1:, 1:], z[1:])
        # proportionality with FinnGen AMD joint effects projected through 1000G EUR LD
        rows.append(dict(subtype=sub, dataset=label, info_share_of_GIGA=";".join(f"{x:.2f}" for x in info_share),
                         z_rs3750848=z[0], z_rs2672592=z[1], z_rs72631113=z[2],
                         p_rs2672592_given_lead=p2(zc(1, 0)), p_lead_given_rs2672592=p2(zc(0, 1)),
                         p_rs72631113_given_lead=p2(zc(2, 0)), p_lead_given_rs72631113=p2(zc(0, 2)),
                         p_side_given_lead_2df=stats.chi2.sf(side, 2), p_lead_given_side_1df=stats.chi2.sf(lead, 1)))
out = pd.DataFrame(rows)
out.to_csv(f"{O}/10q26_giga_remainder_conditional.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(out.round(5).to_string())

# per-variant heterogeneity: MEGASTROKE vs GIGASTROKE-minus-MEGASTROKE vs FinnGen I9_STR (approximately independent parts)
het = []
for sub in ("AS", "IS", "SVS"):
    bG, sG = get(f"GIGA_EUR_{sub}"); bM, sM = get(f"MEGA_EUR_{sub}"); bF, sF = get("FinnGen_I9_STR")
    IG, IM = 1 / sG**2, 1 / sM**2; IF = 0.15 * IG; IR = IG - IM - IF
    bR = (IG * bG - IM * bM - IF * bF) / IR; sR = 1 / np.sqrt(IR)
    for k, s in enumerate(S):
        b = np.array([bM[k], bR[k], bF[k]]); se = np.array([sM[k], sR[k], sF[k]]); w = 1 / se**2
        bm = (w * b).sum() / w.sum(); Q = (w * (b - bm) ** 2).sum()
        het.append(dict(subtype=sub, rsid=s, beta_MEGA=bM[k], se_MEGA=sM[k], beta_GIGA_minus_MEGA_FinnGen=bR[k], se_GIGA_minus_MEGA_FinnGen=sR[k],
                        beta_FinnGen_I9_STR=bF[k], se_FinnGen=sF[k],
                        p_diff_MEGA_vs_rest=p2((bM[k] - bR[k]) / np.hypot(sM[k], sR[k])),
                        Q_3parts=Q, p_Q_2df=stats.chi2.sf(Q, 2), I2=max(0, (Q - 2) / Q) if Q > 0 else 0))
het = pd.DataFrame(het); het.to_csv(f"{O}/10q26_stroke_heterogeneity_3parts.tsv", sep="\t", index=False)
print(het.round(5).to_string())
