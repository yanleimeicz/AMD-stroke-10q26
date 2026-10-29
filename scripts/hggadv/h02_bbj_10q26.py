"""East Asian replication of the 10q26 haplotype partition in BioBank Japan (BBJ) summary statistics (pheweb.jp, GRCh37).
Outcomes: coronary artery disease (CAD; 32,512 cases), myocardial infarction (MI; subset of CAD) and ischemic stroke (IS; BBJ is
part of GIGASTROKE East Asian, so IS is a consistency check only). Models are those used for every other dataset:
two-variant joint model (rs3750848 = ARMS2 p.Ala69Ser haplotype, rs2672592 = linked haplotype) in both directions, three-variant
model with rs72631113 and the joint 2-df test of rs2672592 + rs72631113 conditional on rs3750848, with 1000 Genomes East Asian LD
(r +/- 0.03 sensitivity). The BBJ joint effects are compared with the European CAD datasets without UK Biobank participants."""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import json
import numpy as np, pandas as pd
from scipy import stats

P = PROJ_ROOT; D = f"{P}/derived/hggadv"; O = f"{P}/results/hggadv"
S = ["rs3750848", "rs2672592", "rs72631113"]; RISK = {"rs3750848": "G", "rs2672592": "T", "rs72631113": "T"}
# 1000 Genomes EAS correlations between risk alleles (plink, as used for GIGASTROKE East Asian):
# (rs3750848, rs2672592), (rs3750848, rs72631113), (rs2672592, rs72631113)
R_EAS = [0.821883203381113, 0.678589714039345, 0.541761017423735]
p2 = lambda z: 2 * stats.norm.sf(np.abs(z))


def load(ph):
    d = json.load(open(f"{D}/bbj_{ph}_b37.json"))["data"]
    t = pd.DataFrame({k: d[k] for k in ("rsid", "ref", "alt", "af", "beta", "sebeta", "pvalue", "position")})
    t = t[t.rsid.isin(S)].drop_duplicates("rsid").set_index("rsid").loc[S]
    sgn = np.where(t.alt == pd.Series(RISK).loc[S].values, 1, np.where(t.ref == pd.Series(RISK).loc[S].values, -1, np.nan))
    assert not np.isnan(sgn).any()
    t["beta_risk"] = t.beta * sgn; t["eaf_risk"] = np.where(sgn > 0, t.af, 1 - t.af)
    return t


def joint2(b, se, r):
    z = b / se; V = np.linalg.inv(np.array([[1, r], [r, 1]])); bs = V @ z
    bj = bs * se; sj = np.sqrt(np.diag(V)) * se; cov = V[0, 1] * se[0] * se[1]
    return dict(b_lead=bj[0], s_lead=sj[0], b_side=bj[1], s_side=sj[1], b_con=bj[1] - bj[0], s_con=np.sqrt(sj[0] ** 2 + sj[1] ** 2 - 2 * cov))


def test2df(z, r):
    """Joint test of variants 2,3 conditional on variant 1 (and 1-df test of variant 1 conditional on 2,3)."""
    R = np.array([[1, r[0], r[1]], [r[0], 1, r[2]], [r[1], r[2], 1]])
    zc = z[1:] - R[1:, :1] @ z[:1]; Rc = R[1:, 1:] - R[1:, :1] @ R[:1, 1:]
    chi = float(zc @ np.linalg.solve(Rc, zc)); p_side = stats.chi2.sf(chi, 2)
    z1 = z[0] - R[:1, 1:] @ np.linalg.solve(R[1:, 1:], z[1:]); v1 = 1 - R[:1, 1:] @ np.linalg.solve(R[1:, 1:], R[1:, :1])
    p_lead = p2(float(np.ravel(z1)[0] / np.sqrt(np.ravel(v1)[0])))
    bs = np.linalg.solve(R, z)
    return dict(p_side_2df=p_side, p_lead_given_both=float(p_lead), z_joint3_lead=bs[0] / np.sqrt(np.linalg.inv(R)[0, 0]),
                z_joint3_rs2672592=bs[1] / np.sqrt(np.linalg.inv(R)[1, 1]), z_joint3_rs72631113=bs[2] / np.sqrt(np.linalg.inv(R)[2, 2]))


rows, marg = [], []
N = {"CAD": "32,512 / 146,214", "MI": "14,992 / 146,214", "IS": "22,664 / 152,022"}
for ph in ("CAD", "MI", "IS"):
    t = load(ph)
    for s_ in S:
        x = t.loc[s_]; marg.append(dict(trait=f"BBJ_{ph}", rsid=s_, risk=RISK[s_], eaf_risk=x.eaf_risk, beta_risk=x.beta_risk, se=x.sebeta,
                                          z=x.beta_risk / x.sebeta, p=x.pvalue))
    b = t.beta_risk.values; se = t.sebeta.values; z = b / se
    for dr in (0.0, -0.03, 0.03):
        r = [min(0.999, v + dr) for v in R_EAS]
        j = joint2(b[:2], se[:2], r[0]); k = test2df(z, r)
        rows.append(dict(trait=f"BBJ_{ph}", cases_controls=N[ph], ld="1000G EAS", delta_r=dr, r_lead_side=r[0], **j,
                         p_lead=p2(j["b_lead"] / j["s_lead"]), p_side=p2(j["b_side"] / j["s_side"]), p_con=p2(j["b_con"] / j["s_con"]), **k))
res = pd.DataFrame(rows); mg = pd.DataFrame(marg)

# comparison with European CAD datasets without UK Biobank (joint effects from results/revision6)
eu = pd.read_csv(f"{P}/results/revision6/cad_independent_validation.tsv", sep="\t", nrows=4)
eu = eu[eu.dataset.isin(["FinnGen R12 I9_CORATHER", "CARDIoGRAMplusC4D 2015 (Nikpay; 60,801 cases)"])]
w = 1 / eu.s_side ** 2; b_eu = float((w * eu.b_side).sum() / w.sum()); s_eu = float(1 / np.sqrt(w.sum()))
wl = 1 / eu.s_lead ** 2; bl_eu = float((wl * eu.b_lead).sum() / wl.sum()); sl_eu = float(1 / np.sqrt(wl.sum()))
bb = res[(res.trait == "BBJ_CAD") & (res.delta_r == 0)].iloc[0]
cmp_ = dict(european_side_fixed=b_eu, european_side_se=s_eu, bbj_side=bb.b_side, bbj_side_se=bb.s_side,
            p_diff_side=float(p2((bb.b_side - b_eu) / np.hypot(bb.s_side, s_eu))),
            european_lead_fixed=bl_eu, european_lead_se=sl_eu, bbj_lead=bb.b_lead, bbj_lead_se=bb.s_lead,
            p_diff_lead=float(p2((bb.b_lead - bl_eu) / np.hypot(bb.s_lead, sl_eu))),
            power_side_at_european_effect=float(stats.norm.sf(1.96 - abs(b_eu) / bb.s_side) + stats.norm.cdf(-1.96 - abs(b_eu) / bb.s_side)))
# pooled independent CAD (FinnGen + C4D + BBJ), fixed effect, side and lead
allb = np.r_[eu.b_side.values, bb.b_side]; alls = np.r_[eu.s_side.values, bb.s_side]
wa = 1 / alls ** 2; bp = float((wa * allb).sum() / wa.sum()); sp = float(1 / np.sqrt(wa.sum()))
Q = float((wa * (allb - bp) ** 2).sum())
cmp_.update(pooled3_side=bp, pooled3_side_se=sp, pooled3_side_p=float(p2(bp / sp)), pooled3_Q_p=float(stats.chi2.sf(Q, 2)))
allbl = np.r_[eu.b_lead.values, bb.b_lead]; allsl = np.r_[eu.s_lead.values, bb.s_lead]; wl3 = 1 / allsl ** 2
bpl = float((wl3 * allbl).sum() / wl3.sum()); spl = float(1 / np.sqrt(wl3.sum()))
cmp_.update(pooled3_lead=bpl, pooled3_lead_se=spl, pooled3_lead_p=float(p2(bpl / spl)))

import os
os.makedirs(O, exist_ok=True)
mg.to_csv(f"{O}/bbj_10q26_marginal.tsv", sep="\t", index=False)
res.to_csv(f"{O}/bbj_10q26_conditional.tsv", sep="\t", index=False)
pd.Series(cmp_).to_csv(f"{O}/bbj_10q26_vs_european_cad.tsv", sep="\t", header=["value"])
pd.set_option("display.width", 250)
print(mg.round(4).to_string(index=False)); print()
print(res[["trait", "delta_r", "b_lead", "p_lead", "b_side", "p_side", "b_con", "p_con", "p_side_2df", "p_lead_given_both",
           "z_joint3_lead", "z_joint3_rs2672592", "z_joint3_rs72631113"]].round(4).to_string(index=False)); print()
print(pd.Series(cmp_).round(4).to_string())
