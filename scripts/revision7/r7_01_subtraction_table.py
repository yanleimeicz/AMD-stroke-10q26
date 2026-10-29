"""Round 7: transparent table of the GIGASTROKE 'remainder' construction at 10q26 (per variant and joint model)."""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import importlib.util, numpy as np, pandas as pd, sys
P = PROJ_ROOT
spec = importlib.util.spec_from_file_location("pool", f"{P}/scripts/revision6/r6_06_pool.py")
src = open(f"{P}/scripts/revision6/r6_06_pool.py", encoding="utf-8").read()
cut = src.index("nz = pd.read_csv")          # only the definitions (data loading + remainder/joint), no pooling
ns = {}; exec(src[:cut], ns)
rows = []; jrows = []
from scipy.stats import norm
for sub in ("IS", "AS"):
    bG, sG = ns["get"](f"GIGA_EUR_{sub}"); bM, sM = ns["get"](f"MEGA_EUR_{sub}"); bF, sF = ns["get"]("FinnGen_I9_STR")
    bU, sU = ns["get"]("IS_I63" if sub == "IS" else "AS_sr6150_3", "uk")
    IG, IM = 1 / sG ** 2, 1 / sM ** 2; IF, IU = 0.15 * IG, ns["F_UKB"][sub] * IG; IR = IG - IM - IF - IU
    bR, sR, sR0 = ns["remainder"](sub)
    for i, snp in enumerate(ns["S"][:2]):
        rows.append(dict(outcome=sub, variant=snp, risk_allele=ns["RISK"][snp],
                         b_GIGA=bG[i], se_GIGA=sG[i], b_MEGA=bM[i], se_MEGA=sM[i], b_FinnGen=bF[i], se_FinnGen=sF[i],
                         b_UKB=bU[i], se_UKB=sU[i], info_share_MEGA=IM[i] / IG[i], info_share_FinnGen=0.15, info_share_UKB=ns["F_UKB"][sub],
                         info_share_remainder=IR[i] / IG[i], b_remainder=bR[i], se_remainder_nominal=sR0[i], se_remainder_with_proxy_error=sR[i]))
    for lab, b_, s_, ld in (("GIGASTROKE EUR (full, before subtraction)", bG, sG, "KG_EUR"), ("MEGASTROKE (subtracted)", bM, sM, "KG_EUR"),
                            ("FinnGen R12 (subtracted)", bF, sF, "FinnGen_SISu"), ("UK Biobank (subtracted)", bU, sU, "UKB_British"),
                            ("GIGASTROKE remainder (after subtraction, with proxy error)", bR, sR, "KG_EUR")):
        j = ns["joint"](np.asarray(b_[:2]), np.asarray(s_[:2]), ns["R"][ld])
        jrows.append(dict(outcome=sub, dataset=lab, LD_source=ld, r_LD=ns["R"][ld],
                          b_AMD_haplotype=j["b_lead"], se_AMD_haplotype=j["s_lead"], P_AMD_haplotype=2 * norm.sf(abs(j["b_lead"] / j["s_lead"])),
                          b_rs2672592_haplotype=j["b_side"], se_rs2672592_haplotype=j["s_side"], P_rs2672592_haplotype=2 * norm.sf(abs(j["b_side"] / j["s_side"])),
                          cov_joint=j["cov"], b_difference=j["b_con"], se_difference=j["s_con"], P_difference=2 * norm.sf(abs(j["b_con"] / j["s_con"]))))
pd.DataFrame(rows).to_csv(f"{P}/results/revision7/remainder_construction.tsv", sep="\t", index=False)
pd.DataFrame(jrows).to_csv(f"{P}/results/revision7/remainder_joint_before_after.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).to_string()); print(pd.DataFrame(jrows).to_string())
