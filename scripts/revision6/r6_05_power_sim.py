"""Round 6, item D: power of the 10q26 stroke haplotype tests.
(1) Minimum detectable contrast (80% power, two-sided alpha 0.05) from the observed SEs per dataset and pooled.
(2) Simulation of summary statistics for two linked markers (lead rs3750848, side rs2672592; r = 0.62 EUR) under
    H_shared  : the stroke-causal variant is the lead (AMD) haplotype (beta_L = b, beta_S = 0), and
    H_distinct: the stroke-causal variant is the side haplotype       (beta_L = 0, beta_S = b),
    with b = 0.031 log OR (pooled side-haplotype estimate), allele frequencies 0.23 (lead) and 0.43 (side), and effective
    sample sizes of each dataset. Tests: side|lead (rejects H_shared), lead|side (rejects H_distinct), contrast b_S - b_L,
    and the proportionality test against the AMD pattern (AMD causal = lead; equivalent to side|lead).
    Marginal z ~ MVN(sqrt(Neff/4) * R gamma, R), gamma_j = beta_j * sqrt(2 p_j q_j) (standardized effects), in-sample LD.
(3) Number of ischemic-stroke cases (1:10 case:control, in-sample LD) needed for 80% power to discriminate the two
    configurations (side|lead under H_distinct and lead|side under H_shared, each alpha 0.05)."""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import numpy as np, pandas as pd
from scipy import stats
P = PROJ_ROOT; O = f"{P}/results/revision6"
rng = np.random.default_rng(20261007)
r = 0.62; pL, pS = 0.23, 0.43; b = 0.031
R = np.array([[1, r], [r, 1]]); Ri = np.linalg.inv(R); Lc = np.linalg.cholesky(R)
sd = np.sqrt([2 * pL * (1 - pL), 2 * pS * (1 - pS)])
def neff(n1, n0): return 4 / (1 / n1 + 1 / n0)
DS = {"MEGASTROKE IS": neff(34217, 406111), "GIGASTROKE remainder IS (GIGA - MEGA - FinnGen)": neff(62100, 1234808) - neff(34217, 406111) - 0.15 * neff(62100, 1234808),
      "FinnGen I9_STR": neff(34110, 400000), "UKB I63": neff(4275, 416256)}
DS["pooled 4 datasets"] = sum(DS.values()); DS["GIGASTROKE EUR full IS"] = neff(62100, 1234808)

# (1) minimum detectable contrast from observed SEs (round 5 per-dataset estimates)
dd = pd.read_csv(f"{P}/results/revision5/10q26_stroke_contrast_datasets.tsv", sep="\t")
dd = dd[dd.subtype == "IS"]
rows = []
for _, x in dd.iterrows():
    rows.append(dict(dataset=x.dataset, se_contrast=x.s_contrast, mde80_contrast=2.8016 * x.s_contrast, observed_contrast=x.b_contrast,
                     se_side=x.s_side, mde80_side=2.8016 * x.s_side))
w = 1 / dd[dd.dataset != "GIGASTROKE_EUR_full"].s_contrast ** 2
rows.append(dict(dataset="pooled (4 datasets, fixed effect)", se_contrast=1 / np.sqrt(w.sum()), mde80_contrast=2.8016 / np.sqrt(w.sum())))
hk = pd.read_csv(f"{P}/results/revision5/10q26_stroke_contrast_hksj.tsv", sep="\t")
h = hk[(hk.subtype == "IS") & hk.set.str.startswith("main") & (hk.quantity == "contrast_side_minus_lead")].iloc[0]
se_hk = (h.HKSJ_DL_hi - h.HKSJ_DL_lo) / (2 * stats.t.ppf(0.975, 3))
rows.append(dict(dataset="pooled (4 datasets, Hartung-Knapp, t3)", se_contrast=se_hk,
                 mde80_contrast=se_hk * (stats.t.ppf(0.975, 3) + stats.t.ppf(0.8, 3)), observed_contrast=h.DL_est))
pd.DataFrame(rows).to_csv(f"{O}/contrast_power.tsv", sep="\t", index=False)

# (2) simulation
def sim(ne, truth, nsim=20000):
    beta = np.array([b, 0.0]) if truth == "shared (lead causal)" else np.array([0.0, b])
    mu = np.sqrt(ne / 4) * (R @ (beta * sd))
    z = mu + rng.standard_normal((nsim, 2)) @ Lc.T
    zj = z @ Ri.T; vj = np.diag(Ri)
    side_given_lead = 2 * stats.norm.sf(np.abs(zj[:, 1] / np.sqrt(vj[1])))
    lead_given_side = 2 * stats.norm.sf(np.abs(zj[:, 0] / np.sqrt(vj[0])))
    # contrast on log-OR scale: b_j = zj_j / (sqrt(ne/4) sd_j); var from Ri
    sc = 1 / (np.sqrt(ne / 4) * sd)
    bS, bL = zj[:, 1] * sc[1], zj[:, 0] * sc[0]
    vS, vL, cSL = Ri[1, 1] * sc[1] ** 2, Ri[0, 0] * sc[0] ** 2, Ri[0, 1] * sc[0] * sc[1]
    pc = 2 * stats.norm.sf(np.abs((bS - bL) / np.sqrt(vS + vL - 2 * cSL)))
    return dict(rej_side_given_lead=(side_given_lead < 0.05).mean(), rej_lead_given_side=(lead_given_side < 0.05).mean(),
                rej_contrast=(pc < 0.05).mean(), mean_z_side_given_lead=np.mean(zj[:, 1] / np.sqrt(vj[1])),
                mean_z_lead_given_side=np.mean(zj[:, 0] / np.sqrt(vj[0])))
srows = []
for nm, ne in DS.items():
    for truth in ("shared (lead causal)", "distinct (side causal)"):
        d = sim(ne, truth); d.update(dataset=nm, neff=ne, truth=truth); srows.append(d)
# (3) cases needed (1:10), in-sample LD
def power_pair(ncase):
    ne = neff(ncase, 10 * ncase)
    a = sim(ne, "distinct (side causal)", 8000)["rej_side_given_lead"]; c = sim(ne, "shared (lead causal)", 8000)["rej_lead_given_side"]
    return min(a, c), a, c
grid = [10000, 20000, 40000, 60000, 80000, 100000, 120000, 150000, 200000, 250000, 300000]
for n in grid:
    m, a, c = power_pair(n)
    srows.append(dict(dataset=f"hypothetical IS GWAS, {n} cases (1:10)", neff=neff(n, 10 * n), truth="both (min power)",
                      rej_side_given_lead=a, rej_lead_given_side=c, rej_contrast=np.nan))
pd.DataFrame(srows).to_csv(f"{O}/simulation_10q26.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 20)
print(pd.DataFrame(rows).round(4).to_string()); print(pd.DataFrame(srows).round(3).to_string())
