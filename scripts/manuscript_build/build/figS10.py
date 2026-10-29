"""Figure S10: LDSC genetic correlation (former Fig. 2c) and effect-direction concordance (former Fig. 3b)."""
from rev_common import *
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap

H = 145
fig = plt.figure(figsize=(W, H * MM))
F = lambda mm: mm / H

# ---------------- a: LDSC heatmap with ICH (former Fig. 2c)
l1 = tsv(R / "csvd/ldsc_rg.tsv"); l0 = tsv(R / "ldsc/ldsc_rg_all.tsv"); l2 = tsv(RR / "extra/ldsc_rg_ich_earlyamd.tsv")
ld = pd.concat([d[["trait1", "trait2", "rg", "se", "p"]] for d in (l2, l1, l0)]).drop_duplicates(["trait1", "trait2"])
lf = tsv(RR / "mr/fdr_all_layers.tsv"); lf = lf[lf.layer == "LDSC genetic correlation (FinnGen AMD)"].set_index("item").fdr
ROWS = [("AMD_FG", "AMD"), ("DRYAMD_FG", "Dry AMD"), ("WETAMD_FG", "Wet AMD")]
COLS = [("AS", "Any\nstroke"), ("IS", "Ischemic"), ("LAS", "Large-\nartery"), ("CES", "Cardio-\nembolic"), ("SVS", "Small-\nvessel"),
        ("LACUNAR", "Lacunar"), ("ICH_UKBWGS", "ICH"), ("WMH", "WMH"), ("WM_MD", "WM\nMD"), ("WM_FA", "WM\nFA"),
        ("EPVS_WM", "EPVS\nWM"), ("EPVS_BG", "EPVS\nBG"), ("EPVS_HIP", "EPVS\nHIP")]
M = np.full((len(ROWS), len(COLS)), np.nan); PM = M.copy(); FM = M.copy()
for i, (r_, _) in enumerate(ROWS):
    for j, (c_, _) in enumerate(COLS):
        d = ld[((ld.trait1 == r_) & (ld.trait2 == c_)) | ((ld.trait1 == c_) & (ld.trait2 == r_))]
        if len(d): M[i, j], PM[i, j] = d.rg.iloc[0], d.p.iloc[0]
        FM[i, j] = lf.get(f"{r_} ~ {c_}", np.nan)
axc = fig.add_axes([0.14, F(102), 0.725, F(22.6)])
cmap = LinearSegmentedColormap.from_list("div", [ANT, "#f0efec", CON])
im = axc.imshow(M, cmap=cmap, norm=TwoSlopeNorm(0, -0.4, 0.4), aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        if np.isnan(M[i, j]): continue
        mark = "*" if FM[i, j] < 0.05 else ("†" if PM[i, j] < 0.05 else "")
        axc.text(j, i, f"{M[i, j]:.2f}".replace("-", "−") + mark, ha="center", va="center", fontsize=7, color=INK)
axc.set_xticks(range(len(COLS))); axc.set_xticklabels([c[1] for c in COLS], fontsize=7)
axc.set_yticks(range(len(ROWS))); axc.set_yticklabels([r[1] for r in ROWS])
axc.tick_params(length=0); [s.set_visible(False) for s in axc.spines.values()]
axc.xaxis.tick_top()
cb = fig.add_axes([0.885, F(102), 0.012, F(22.6)]); cbar = fig.colorbar(im, cax=cb, orientation="vertical")
cbar.set_label("rg", fontsize=7); cbar.ax.tick_params(labelsize=7, length=2); cbar.outline.set_linewidth(0.4)
panel(fig, (0.005, F(H - 6)), "a")
fig.text(0.14, F(H - 7.5), "LDSC genetic correlation of FinnGen AMD with stroke and CSVD traits (†, P < 0.05 but FDR ≥ 0.05; no FDR-significant pair)",
         fontsize=7, color=INK2)

# ---------------- b: effect-direction concordance, one variant per locus (former Fig. 3b)
cc = tsv(RR / "mr/concordance_distance_pruned.tsv")
TR = [("IS", "Ischemic stroke"), ("AS", "Any stroke"), ("LAS", "Large-artery stroke"), ("SVS", "Small-vessel stroke"),
      ("LACUNAR", "Lacunar stroke"), ("WMH", "WMH volume"), ("WM_MD", "White-matter MD"), ("WM_FA", "White-matter FA"),
      ("EPVS_WM", "EPVS, white matter"), ("EPVS_BG", "EPVS, basal ganglia"), ("EPVS_HIP", "EPVS, hippocampus")]
axb = fig.add_axes([0.25, F(17), 0.42, F(62)])
for i, (t, tl) in enumerate(TR):
    y = len(TR) - 1 - i
    for ex, mk, fill, dy in [("AMD_FG", "o", True, 0.15), ("AMD_IAMDGC", "s", False, -0.15)]:
        r = cc[(cc.exposure == ex) & (cc.trait == t)].iloc[0]
        axb.plot(r.frac_same, y + dy, mk, ms=4, mfc=INK if fill else "white", mec=INK, mew=0.8)
    rf = cc[(cc.exposure == "AMD_FG") & (cc.trait == t)].iloc[0]; ri = cc[(cc.exposure == "AMD_IAMDGC") & (cc.trait == t)].iloc[0]
    axb.text(1.03, y, f"{int(rf.n_same)}/{int(rf.n_loci)}   {rf.fdr_perm_sumz:.2f}   |   {int(ri.n_same)}/{int(ri.n_loci)}   {ri.fdr_perm_sumz:.2f}",
             transform=axb.get_yaxis_transform(), va="center", fontsize=7)
axb.axvline(0.5, color=GREY, lw=0.6, ls="--", zorder=0)
axb.set_yticks(range(len(TR))); axb.set_yticklabels([t[1] for t in TR][::-1], fontsize=7)
axb.set_xlim(0.15, 0.85); axb.set_ylim(-0.7, len(TR) - 0.3)
axb.set_xlabel("Fraction of AMD loci with concordant effect\n(AMD risk allele raises the vascular trait)")
axb.text(1.03, len(TR) - 0.2, "FinnGen: n/N  FDR   |   IAMDGC: n/N  FDR", transform=axb.get_yaxis_transform(), fontsize=7, fontweight="bold")
axb.legend(handles=[plt.Line2D([], [], marker="o", ls="none", color=INK, ms=4, label="FinnGen AMD"),
                    plt.Line2D([], [], marker="s", ls="none", mfc="white", mec=INK, ms=4, label="IAMDGC advanced AMD")],
           loc="upper left", fontsize=7, bbox_to_anchor=(0.0, 1.10), ncol=2, handletextpad=0.3)
panel(fig, (0.005, F(90)), "b")
save_supp_png(fig, "FigS10")
