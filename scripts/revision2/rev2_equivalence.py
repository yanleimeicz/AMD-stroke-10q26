"""Round-2 minor revision, item 3: equivalence (TOST) anchored to the binary-exposure translation.
Burgess & Labrecque 2018: an MR estimate per unit log-odds of liability (beta_MR) relates to the causal log-OR of the
binary exposure (AMD itself) approximately as log(OR_bin) ~ beta_MR / [p(1-p)], p = AMD prevalence in the outcome population.
Outputs: results/revision2/equivalence_bounds.tsv (translation table) and equivalence_tost_by_exclusion.tsv
Usage: python rev2_equivalence.py"""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import numpy as np, pandas as pd
from scipy import stats

P = PROJ_ROOT; O = f"{P}/results/revision2"
PREV = [0.03, 0.05, 0.07, 0.10, 0.15, 0.20]
z90 = stats.norm.isf(0.05)

# ---------------------------------------------------------------- translation table
rows = []
for p in PREV:
    w = p * (1 - p)
    for bound in (1.025, 1.05, 1.10):
        rows.append(dict(kind="bound_to_causal_OR", prevalence_p=p, p1mp=w, mr_bound_upper=bound, mr_bound_lower=1 / bound if bound != 1.05 else 0.95,
                         causal_OR_upper=np.exp(np.log(bound) / w),
                         causal_OR_lower=np.exp(np.log(1 / bound if bound != 1.05 else 0.95) / w),
                         causal_OR_upper_linear=1 + (bound - 1) / w))                      # DA's linear form, MR ~ 1 + p(1-p)(OR-1)
    for orx, lab in ((2.08, "NHANES OR 2.08"), (1.97, "NHANES Model 3 OR 1.97")):
        rows.append(dict(kind="observational_OR_to_expected_MR", prevalence_p=p, p1mp=w, observational_OR=orx, label=lab,
                         expected_MR_OR_per_unit_logodds=np.exp(w * np.log(orx)),
                         expected_MR_OR_linear=1 + w * (orx - 1),
                         inside_pm5pct=np.exp(w * np.log(orx)) < 1.05))
tr = pd.DataFrame(rows)

# ---------------------------------------------------------------- TOST by exclusion set (IVW)
mr = pd.read_csv(f"{P}/results/revision/mr/mr_revised_all.tsv", sep="\t")
mr = mr[(mr.method == "IVW") & mr.exclusion.isin(["main", "no_10q26", "no_CFH", "no_major3", "no_all"])].copy()
mr["binary"] = mr.binary.astype(str).str.upper() == "TRUE"
def tost(b, se, bnd): return np.maximum(stats.norm.sf((b + bnd) / se), stats.norm.cdf((b - bnd) / se))
out = []
for _, r in mr.iterrows():
    b, se = r.b, r.se
    margin = abs(b) + z90 * se                     # smallest symmetric bound with TOST P < 0.05 (90% CI inside +/- margin)
    o = dict(exposure=r.exposure, scheme=r.scheme, outcome=r.outcome, exclusion=r.exclusion, nsnp=r.nsnp, b=b, se=se, p=r.p,
             binary=r.binary, OR=np.exp(b) if r.binary else np.nan,
             lo95=np.exp(b - 1.96 * se) if r.binary else b - 1.96 * se, hi95=np.exp(b + 1.96 * se) if r.binary else b + 1.96 * se,
             lo90=np.exp(b - z90 * se) if r.binary else b - z90 * se, hi90=np.exp(b + z90 * se) if r.binary else b + z90 * se)
    if r.binary:
        for bnd in (1.025, 1.05, 1.07, 1.10, 1.15):
            o[f"tost_p_{bnd}"] = tost(b, se, np.log(bnd))
        o["min_equiv_margin_OR"] = np.exp(margin)
        for p in PREV:   # largest causal OR of AMD itself that remains compatible (upper end of the 90% CI translated)
            o[f"causal_OR_excluded_above_p{p}"] = np.exp(margin / (p * (1 - p)))
    else:
        for bnd in (0.025, 0.05, 0.10):
            o[f"tost_p_{bnd}SD"] = tost(b, se, bnd)
        o["min_equiv_margin_SD"] = margin
    out.append(o)
eq = pd.DataFrame(out)
eq["exclusion"] = pd.Categorical(eq.exclusion, ["no_10q26", "main", "no_CFH", "no_major3", "no_all"])
eq["role"] = eq.exclusion.map({"no_10q26": "primary (10q26 removed)", "main": "all instruments (pleiotropy-biased)", "no_CFH": "sensitivity",
                               "no_major3": "sensitivity (10q26+CFH+C3 removed)", "no_all": "sensitivity (all major loci removed)"})
eq = eq.sort_values(["scheme", "exposure", "outcome", "exclusion"])
eq.to_csv(f"{O}/equivalence_tost_by_exclusion.tsv", sep="\t", index=False)

# condensed rows for the primary scheme added to the bounds file
key = eq[(eq.scheme == "A_dist") & eq.exposure.isin(["AMD_FG", "DRYAMD_FG", "WETAMD_FG"]) & eq.outcome.isin(["AS", "IS", "CES", "LAS", "SVS", "WMH"])]
pd.concat([tr, key.assign(kind="achieved_margin_A_dist")], ignore_index=True).to_csv(f"{O}/equivalence_bounds.tsv", sep="\t", index=False)

pd.set_option("display.width", 250)
print(tr.round(3).to_string())
cols = ["exposure", "outcome", "exclusion", "nsnp", "OR", "lo95", "hi95", "tost_p_1.05", "tost_p_0.05SD", "min_equiv_margin_OR", "min_equiv_margin_SD",
        "causal_OR_excluded_above_p0.05", "causal_OR_excluded_above_p0.1"]
print(key[cols].round(4).to_string())
