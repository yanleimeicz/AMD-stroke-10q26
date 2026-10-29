"""Round 5, item 5: two-variant conditional eQTL analyses at 10q26 in both directions.
AMD haplotype marker: rs3750848-G (EyeGEx: rs36212733-C, r = 1 with rs3750848 in 1000G EUR); vascular-side marker: rs2672592-T.
z_b|a = (z_b - r z_a) / sqrt(1 - r^2), r = 1000G EUR correlation (risk-coded), evaluated at r and r +/- 0.03.
Output: results/revision5/eqtl_conditional_both_directions.tsv"""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import numpy as np, pandas as pd
from scipy import stats

P = PROJ_ROOT
lk = pd.read_csv(f"{P}/results/tissue_eqtl/key_variant_eqtl_lookup.tsv", sep="\t")
kg = pd.read_csv(f"{P}/results/revision/10q26/ld_1000G_EUR_r_dprime.tsv", sep=r"\s+")
x = kg[((kg.SNP_A == "rs3750848") & (kg.SNP_B == "rs2672592")) | ((kg.SNP_A == "rs2672592") & (kg.SNP_B == "rs3750848"))]
r0 = float(x.R.iloc[0]); print("r(rs3750848, rs2672592) =", r0)
rows = []
for gene in ("HTRA1", "ARMS2", "PLEKHA1"):
    for dsn in sorted(lk.dataset.unique()):
        d = lk[(lk.gene == gene) & (lk.dataset == dsn)].set_index("SNP")
        lead = "rs3750848" if "rs3750848" in d.index else ("rs36212733" if "rs36212733" in d.index else None)
        if lead is None or "rs2672592" not in d.index: continue
        za = d.loc[lead, "b_allele"] / d.loc[lead, "se"]; zb = d.loc["rs2672592", "b_allele"] / d.loc["rs2672592", "se"]
        for dr in (-0.03, 0.0, 0.03):
            r = r0 + dr; den = np.sqrt(1 - r ** 2)
            zba = (zb - r * za) / den; zab = (za - r * zb) / den
            rows.append(dict(gene=gene, dataset=dsn, lead_marker=f"{lead}-{d.loc[lead, 'allele']}", side_marker="rs2672592-T",
                             n=d.loc[lead, "n"], r=r, r_shift=dr, z_lead_marginal=za, z_side_marginal=zb,
                             z_side_given_lead=zba, p_side_given_lead=2 * stats.norm.sf(abs(zba)),
                             z_lead_given_side=zab, p_lead_given_side=2 * stats.norm.sf(abs(zab))))
out = pd.DataFrame(rows); out.to_csv(f"{P}/results/revision5/eqtl_conditional_both_directions.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(out[out.r_shift == 0].round(4).to_string())
print(out[(out.gene == "HTRA1") & out.dataset.isin(["EyeGEx_Retina", "eQTLGen_Blood", "GTEx_Whole_Blood"])].round(4).to_string())
