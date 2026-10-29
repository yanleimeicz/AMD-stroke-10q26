"""Figure 2 (HGG Advances): MR estimates per unit log-odds of AMD liability, one instrument per locus (IVW),
all loci / without 10q26 / without CFH / without 10q26, CFH and C3. The NHANES benchmark band of the
Human Genomics version is not shown (moved to the Supplemental Information)."""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
from hg_common import *

fig = plt.figure(figsize=(W, 112 * MM))
mr = tsv(RR / "mr/mr_revised_all.tsv")
mr = mr[(mr.scheme == "A_dist") & (mr.method == "IVW")]
fdr = tsv(RR / "mr/fdr_all_layers.tsv"); fdr = fdr[fdr.layer == "MR forward (distance-pruned instruments)"].set_index("item").fdr
EXP = [("AMD_FG", "AMD (FinnGen)"), ("DRYAMD_FG", "Dry AMD (FinnGen)"), ("WETAMD_FG", "Wet AMD (FinnGen)"), ("AMD_IAMDGC", "Advanced AMD (IAMDGC)")]
OUTC = [("AS", "Any stroke"), ("IS", "Ischemic stroke"), ("LAS", "Large-artery stroke"), ("CES", "Cardioembolic stroke"),
        ("SVS", "Small-vessel stroke"), ("WMH", "WMH volume (exp β)")]
SETS = [("main", INK, "o", 0.29), ("no_10q26", AQUA, "s", 0.10), ("no_CFH", VIOL, "^", -0.10), ("no_major3", GREY, "D", -0.29)]
x0, wb = 0.155, 0.198
for j, (e, el) in enumerate(EXP):
    ax = fig.add_axes([x0 + j * (wb + 0.014), 0.14, wb, 0.66])
    for i, (o, ol) in enumerate(OUTC):
        y = len(OUTC) - 1 - i
        for ex, col, mk, dy in SETS:
            r = mr[(mr.exposure == e) & (mr.outcome == o) & (mr.exclusion == ex)]
            if not len(r): continue
            r = r.iloc[0]
            est, lo, hi = np.exp(r.b), np.exp(r.b - 1.96 * r.se), np.exp(r.b + 1.96 * r.se)
            XL, XH = 0.80, 1.20
            forest(ax, y + dy, est, max(lo, XL), min(hi, XH), color=col, filled=r.p < 0.05, ms=2.8 if mk == "D" else 3.2, marker=mk, log=False)
            if hi > XH: ax.plot(XH, y + dy, ">", ms=3, color=col, clip_on=False, zorder=4)
            if lo < XL: ax.plot(XL, y + dy, "<", ms=3, color=col, clip_on=False, zorder=4)
            if ex == "main":
                q = fdr.get(f"{e} -> {o}", np.nan)
                if q < 0.05:
                    ax.text(hi + 0.006, y + dy, "FDR*", fontsize=7, va="center", color=INK, fontweight="bold")
    ax.axvline(1, color=GREY, lw=0.6, ls="--", zorder=0)
    ax.set_ylim(-0.6, len(OUTC) - 0.4); ax.set_xlim(0.80, 1.20)
    ax.set_xticks([0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15]); ax.set_xticklabels(["", "0.9", "", "1.0", "", "1.1", ""])
    ax.set_yticks(range(len(OUTC))); ax.set_yticklabels([o[1] for o in OUTC][::-1] if j == 0 else [])
    ax.set_title(el, fontsize=7, fontweight="bold", pad=3)
    if j == 0: ax.set_xlabel("OR per unit log-odds of AMD liability (95% CI)", loc="left")
hl = [plt.Line2D([], [], marker="o", ls="none", color=INK, ms=3.4, label="All loci"),
      plt.Line2D([], [], marker="s", ls="none", color=AQUA, ms=3.4, label="Without 10q26"),
      plt.Line2D([], [], marker="^", ls="none", color=VIOL, ms=3.6, label="Without CFH"),
      plt.Line2D([], [], marker="D", ls="none", color=GREY, ms=3.0, label="Without 10q26, CFH and C3"),
      plt.Line2D([], [], marker="o", ls="none", mfc="white", mec=INK, ms=3.4, label="Open: P ≥ 0.05")]
fig.legend(handles=hl, loc="lower center", ncol=5, bbox_to_anchor=(0.55, 0.885), handletextpad=0.3, columnspacing=0.9, fontsize=7)
fig.text(0.55, 0.875, "FDR*, FDR < 0.05 among the 24 primary forward tests (all loci)", ha="center", va="center", fontsize=7, color=INK)
save(fig, "Figure2")
