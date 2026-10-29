"""Draft figure for the BBJ 10q26 replication: joint effects of the two haplotypes on CAD (independent datasets) and IS."""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import sys
sys.path.insert(0, f"{PROJ_ROOT}/scripts/hggadv/figures")
from hg_common import *
import pandas as pd
from pathlib import Path

O = Path(f"{PROJ_ROOT}/results/hggadv/BBJ_validation"); O.mkdir(parents=True, exist_ok=True)
eu = pd.read_csv(f"{PROJ_ROOT}/results/revision6/cad_independent_validation.tsv", sep="\t", nrows=4)
bbj = pd.read_csv(f"{PROJ_ROOT}/results/hggadv/bbj_10q26_conditional.tsv", sep="\t")
bbj = bbj[bbj.delta_r == 0].set_index("trait")
cmpv = pd.read_csv(f"{PROJ_ROOT}/results/hggadv/bbj_10q26_vs_european_cad.tsv", sep="\t", index_col=0)["value"]
rows = []
for lab, key in [("CAD · UK Biobank (in-sample LD; EUR)", "UK Biobank coronary atherosclerosis (411.4)"),
                 ("CAD · FinnGen (in-sample LD; FIN)", "FinnGen R12 I9_CORATHER"),
                 ("CAD · CARDIoGRAMplusC4D 2015 (EUR)", "CARDIoGRAMplusC4D 2015 (Nikpay; 60,801 cases)")]:
    x = eu[eu.dataset == key].iloc[0]; rows.append((lab, x.b_lead, x.s_lead, x.b_side, x.s_side))
x = bbj.loc["BBJ_CAD"]; rows.append(("CAD · BioBank Japan (EAS; new)", x.b_lead, x.s_lead, x.b_side, x.s_side))
rows.append(("CAD · pooled FinnGen + C4D + BBJ", np.nan, np.nan, cmpv["pooled3_side"], cmpv["pooled3_side_se"]))
rows[-1] = ("CAD · pooled FinnGen + C4D + BBJ", cmpv["pooled3_lead"], cmpv["pooled3_lead_se"], cmpv["pooled3_side"], cmpv["pooled3_side_se"])
x = bbj.loc["BBJ_MI"]; rows.append(("MI · BioBank Japan (EAS; subset of CAD)", x.b_lead, x.s_lead, x.b_side, x.s_side))
x = bbj.loc["BBJ_IS"]; rows.append(("IS · BioBank Japan (EAS; within GIGASTROKE EAS)", x.b_lead, x.s_lead, x.b_side, x.s_side))

fig = plt.figure(figsize=(W, 92 * MM)); ax = fig.add_axes([0.46, 0.14, 0.50, 0.74])
n = len(rows)
for i, (lab, bl, sl, bs, ss) in enumerate(rows):
    y = n - 1 - i
    for b, s, col, dy in ((bl, sl, VIOL, 0.15), (bs, ss, AQUA, -0.15)):
        if np.isnan(b): continue
        p = 2 * (1 - __import__("scipy").stats.norm.cdf(abs(b / s)))
        forest(ax, y + dy, b, b - 1.96 * s, b + 1.96 * s, color=col, filled=p < 0.05, log=False,
               marker="D" if "pooled" in lab else "o")
ax.axvline(0, color=GREY, ls="--", lw=0.6, zorder=0)
ax.set_yticks(range(n)); ax.set_yticklabels([r[0] for r in rows][::-1], fontsize=7)
ax.set_xlabel("Joint effect per allele, log OR (95% CI)")
ax.axhline(n - 5.5, color=LGREY, lw=0.6)
ax.legend(handles=[plt.Line2D([], [], marker="o", color=VIOL, ls="none", ms=3.4, label="AMD haplotype (rs3750848-G | rs2672592)"),
                   plt.Line2D([], [], marker="o", color=AQUA, ls="none", ms=3.4, label="rs2672592 haplotype (rs2672592-T | rs3750848)")],
          loc="lower center", bbox_to_anchor=(0.3, 1.0), ncol=1, fontsize=7)
save_supp_png(fig, "FigS14_BBJ_draft")
import shutil
shutil.copy2(f"{PROJ_ROOT}/scripts/hggadv/figures/out/supp/FigS14_BBJ_draft.png", O / "Figure_BBJ_10q26_draft.png")
