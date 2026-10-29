"""Figure S8 (former Fig. 5): the 3p14 (ADAMTS9) locus, panels a-c."""
from rev_common import *
from statistics import NormalDist
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap
from matplotlib.patches import Rectangle

H = 165                                   # mm; two columns: a (left), b over c (right)
fig = plt.figure(figsize=(W, H * MM))
F = lambda mm: mm / H
LEAD = "rs4611812"

# ---------------- a: regional association, 3p14 (AMD, T2D, aortic ADAMTS9 eQTL, IS)
g = tsv(FD0 / "region3p14_gwas.tsv"); t2 = tsv(FD / "region3p14_t2d.tsv"); eq = tsv(FD0 / "region3p14_aorta_eqtl_raw.tsv")
ld = tsv(FD0 / "region3p14_ld.tsv").set_index("SNP").r2
genes = tsv(P / "derived/adamts9/gene_symbols.tsv")
lead_pos = g[(g.SNP == LEAD) & (g.trait == "AMD_FG")].pos38.iloc[0]
X0, X1 = lead_pos - 300_000, lead_pos + 250_000
eq = eq[eq.gene_id == "ENSG00000163638"].rename(columns={"rsid": "SNP", "position": "pos38", "pvalue": "p"})
tracks = [("AMD (FinnGen)", g[g.trait == "AMD_FG"][["SNP", "pos38", "p"]]),
          ("Type 2 diabetes (Xue et al.)", t2[["SNP", "pos38", "p"]]),
          ("Aortic ADAMTS9 eQTL (GTEx)", eq[["SNP", "pos38", "p"]]),
          ("Ischemic stroke (GIGASTROKE)", g[g.trait == "IS"][["SNP", "pos38", "p"]])]
bins = [0, 0.2, 0.4, 0.6, 0.8, 1.01]; cols = ["#d9d8d4", "#b7d3f6", "#6da7ec", "#2a78d6", "#0d366b"]
cm = ListedColormap(cols); nm = BoundaryNorm(bins, cm.N)
top, hh, gap = F(160), F(29), F(2)
for i, (tl, d) in enumerate(tracks):
    ax = fig.add_axes([0.085, top - (i + 1) * hh - i * gap, 0.43, hh])
    d = d[d.pos38.between(X0, X1)].copy(); d["lp"] = -np.log10(d.p.clip(lower=1e-300)); d["r2"] = d.SNP.map(ld).fillna(0)
    d = d.sort_values("r2")
    ax.scatter(d.pos38 / 1e6, d.lp, s=6, c=d.r2, cmap=cm, norm=nm, edgecolors="none", zorder=3)
    s = d[d.SNP == LEAD]
    if len(s): ax.scatter(s.pos38 / 1e6, s.lp, s=34, marker="D", c=ANT, edgecolors=INK, linewidths=0.7, zorder=6)
    ax.set_xlim(X0 / 1e6, X1 / 1e6); ax.set_ylim(0, d.lp.max() * 1.2)
    ax.text(0.01, 0.95, tl, transform=ax.transAxes, va="top", fontsize=7, fontweight="bold")
    ax.set_ylabel("−log$_{10}$P", fontsize=7); ax.tick_params(labelbottom=False)
    if i == 0:
        cax = ax.inset_axes([0.74, 0.58, 0.23, 0.10])
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=nm, cmap=cm), cax=cax, orientation="horizontal", ticks=[0, 0.4, 0.8])
        cb.ax.tick_params(labelsize=7, length=1.5); cb.outline.set_linewidth(0.3)
        cax.set_title(f"r² with {LEAD}", fontsize=7, pad=2)
axg = fig.add_axes([0.085, top - 4 * hh - 3 * gap - F(10.75), 0.43, F(8.6)])
gs = genes[(genes.end38 >= X0) & (genes.start38 <= X1) & genes.gene.isin(["ADAMTS9", "ADAMTS9-AS2", "ADAMTS9-AS1", "PRICKLE2"])]
for _, r in gs.sort_values("start38").iterrows():
    yy = {"ADAMTS9": 0.6, "ADAMTS9-AS2": 0.6, "ADAMTS9-AS1": 0.12, "PRICKLE2": 0.12}[r.gene]
    s0, e0 = max(r.start38, X0) / 1e6, min(r.end38, X1) / 1e6
    axg.plot([s0, e0], [yy, yy], color=INK2, lw=2.2, solid_capstyle="butt")
    axg.text((s0 + e0) / 2, yy + 0.18, r.gene + (" →" if r.strand > 0 else " ←"), ha="center", fontsize=7, style="italic")
axg.set_xlim(X0 / 1e6, X1 / 1e6); axg.set_ylim(-0.1, 1.1); axg.set_yticks([]); axg.spines["left"].set_visible(False)
axg.set_xlabel("Chromosome 3 position (Mb, GRCh38)")
panel(fig, (0.005, 0.975), "a")

# ---------------- b: phenome-wide effects of rs4611812-C (z scores)
ke = tsv(RR / "adamts9/rs4611812C_effects_all_traits.tsv").set_index("trait")
bp = tsv(FD / "rs4611812C_bp_whr.tsv").set_index("trait")
zfromp = lambda p, sgn: sgn * NormalDist().inv_cdf(1 - p / 2)
ITEMS = [("AMD (FinnGen)", ke.loc["AMD_FG", "beta"] / ke.loc["AMD_FG", "se"], ke.loc["AMD_FG", "p"]),
         ("Type 2 diabetes", ke.loc["T2D", "beta"] / ke.loc["T2D", "se"], ke.loc["T2D", "p"]),
         ("HbA1c‡", zfromp(8e-13, 1), 8e-13),
         ("Diastolic BP", bp.loc["DBP", "beta_C"] / bp.loc["DBP", "se"], bp.loc["DBP", "p"]),
         ("Systolic BP", bp.loc["SBP", "beta_C"] / bp.loc["SBP", "se"], bp.loc["SBP", "p"]),
         ("WHRadjBMI", bp.loc["WHRadjBMI", "beta_C"] / bp.loc["WHRadjBMI", "se"], bp.loc["WHRadjBMI", "p"]),
         ("Carotid IMT‡", zfromp(4e-9, 1), 4e-9),
         ("Ischemic stroke", ke.loc["IS", "beta"] / ke.loc["IS", "se"], ke.loc["IS", "p"]),
         ("Large-artery stroke", ke.loc["LAS", "beta"] / ke.loc["LAS", "se"], ke.loc["LAS", "p"]),
         ("Coronary artery disease", ke.loc["CAD", "beta"] / ke.loc["CAD", "se"], ke.loc["CAD", "p"])]
axb = fig.add_axes([0.73, F(95), 0.18, F(65)])
for i, (lab, z, p) in enumerate(ITEMS):
    y = len(ITEMS) - 1 - i
    col = INK if i == 0 else (ANT if z > 0 else CON)
    axb.hlines(y, z - 1.96, z + 1.96, color=col, lw=0.8)
    axb.plot(z, y, "o", ms=3.8, color=col)
    axb.text(1.03, y, f"{p:.0e}".replace("e-0", "e−").replace("e-", "e−"), transform=axb.get_yaxis_transform(), va="center", fontsize=7)
axb.axvline(0, color=GREY, ls="--", lw=0.6, zorder=0)
axb.set_yticks(range(len(ITEMS))); axb.set_yticklabels([it[0] for it in ITEMS][::-1], fontsize=7)
axb.set_ylim(-0.6, len(ITEMS) - 0.3); axb.set_xlim(-10, 20)
axb.set_xlabel("z for rs4611812-C (± 1.96)")
axb.text(1.03, len(ITEMS) - 0.35, "P", transform=axb.get_yaxis_transform(), fontsize=7, fontweight="bold")
axb.text(0.5, -0.20, "Orange: opposite to AMD; blue: same as AMD", transform=axb.transAxes, ha="center", va="top", fontsize=7, color=INK2)
panel(fig, (0.555, 0.975), "b")

# ---------------- c: pairwise PP.H4 matrix + HyPrColoc main cluster
cp = tsv(RR / "adamts9/coloc_pairwise_3p14.tsv"); cp = cp[cp.p12 == 1e-5]
TT = [("AMD_FG", "AMD"), ("AORTA_ADAMTS9", "Aortic eQTL"), ("T2D", "T2D"), ("DBP", "DBP"),
      ("IS", "IS"), ("CIMT", "cIMT"), ("WHRadjBMI", "WHRadjBMI"), ("CAD", "CAD")]
n = len(TT); Mx = np.full((n, n), np.nan)
for i, (a, _) in enumerate(TT):
    for j, (b, _) in enumerate(TT):
        if i == j: continue
        q = cp[((cp.trait1 == a) & (cp.trait2 == b)) | ((cp.trait1 == b) & (cp.trait2 == a))]
        if len(q): Mx[i, j] = q["PP.H4"].iloc[0]
axc = fig.add_axes([0.64, F(18), 0.26, F(44)])
seq = LinearSegmentedColormap.from_list("seq", ["#f4f4f2", "#86b6ef", "#1c5cab", "#0d366b"])
axc.imshow(np.ma.masked_invalid(Mx), cmap=seq, vmin=0, vmax=1, aspect="equal")
for i in range(n):
    for j in range(n):
        if i == j:
            axc.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fc="white", ec="none"))
        elif not np.isnan(Mx[i, j]):
            v = Mx[i, j]; axc.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7, color="white" if v > 0.6 else INK)
axc.add_patch(Rectangle((-0.5, -0.5), 4, 4, fc="none", ec=ANT, lw=1.6))
axc.text(3.6, -0.75, "HyPrColoc cluster, posterior 0.95", ha="right", va="bottom", fontsize=7, color=ANT, fontweight="bold")
axc.set_xticks(range(n)); axc.set_xticklabels([t[1] for t in TT], rotation=45, ha="right", fontsize=7)
axc.set_yticks(range(n)); axc.set_yticklabels([t[1] for t in TT], fontsize=7)
axc.tick_params(length=0); [s.set_visible(False) for s in axc.spines.values()]
cbx = fig.add_axes([0.915, F(18), 0.010, F(44)])
cbb = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, 1), cmap=seq), cax=cbx)
cbb.set_label("Pairwise PP.H4 (coloc, p12 = $10^{-5}$)", fontsize=7); cbb.ax.tick_params(labelsize=7, length=1.5); cbb.outline.set_linewidth(0.3)
panel(fig, (0.555, F(74)), "c")

save_supp_png(fig, "FigS8")
