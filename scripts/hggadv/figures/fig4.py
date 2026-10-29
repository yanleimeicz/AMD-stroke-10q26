import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
from hg_common import *
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch

H = 222
fig = plt.figure(figsize=(W, H * MM))
F = lambda mm: mm / H
LA, LB, LC = "rs3750848", "rs2672592", "rs72631113"

# ---------------- a: regional association, AMD and IS
g = tsv(FD0 / "region10q26_gwas.tsv"); ld = tsv(FD0 / "region10q26_ld.tsv"); genes = tsv(FD0 / "genes_10q26_grch38.tsv")
ra = ld[ld.lead == LA].set_index("SNP").r2; rb = ld[ld.lead == LB].set_index("SNP").r2
X0, X1 = 122.30e6, 122.70e6
TR = [("AMD_FG", "AMD (FinnGen)"), ("IS", "Ischemic stroke (GIGASTROKE)")]
top, hh, gap = F(216), F(20), F(2.3)
axs = []
for i, (t, tl) in enumerate(TR):
    ax = fig.add_axes([0.085, top - (i + 1) * hh - i * gap, 0.455, hh]); axs.append(ax)
    d = g[(g.trait == t) & g.pos38.between(X0, X1)].copy()
    d["lp"] = -np.log10(d.p.clip(lower=1e-320))
    d["ra"] = d.SNP.map(ra).fillna(0); d["rb"] = d.SNP.map(rb).fillna(0)
    grp = np.where((d.ra >= 0.8) & (d.ra >= d.rb), "A", np.where(d.rb >= 0.8, "B", "N"))
    for k, c, z in [("N", LGREY, 2), ("A", VIOL, 3), ("B", AQUA, 4)]:
        s = d[grp == k]; ax.scatter(s.pos38 / 1e6, s.lp, s=6, c=c, edgecolors="none", zorder=z)
    for snp, mk, c in [(LA, "D", VIOL), (LB, "^", AQUA), (LC, "s", AQUA)]:
        s = d[d.SNP == snp]
        if len(s): ax.scatter(s.pos38 / 1e6, s.lp, s=34, marker=mk, c=c, edgecolors=INK, linewidths=0.7, zorder=6)
    ax.set_xlim(X0 / 1e6, X1 / 1e6); ax.set_ylim(0, d.lp.max() * (1.15 if i == 0 else 1.45))
    ax.text(0.01, 0.95, tl, transform=ax.transAxes, va="top", fontsize=7, fontweight="bold")
    ax.set_ylabel("−log$_{10}$P", fontsize=7); ax.tick_params(labelbottom=False)
axg = fig.add_axes([0.085, top - 2 * hh - 1 * gap - F(11), 0.455, F(8.6)])
for k, (_, r) in enumerate(genes[genes.end.between(X0, X1 + 2e5) | genes.start.between(X0, X1)].sort_values("start").iterrows()):
    yy = 0.6 if k % 2 == 0 else 0.15
    s0, e0 = max(r.start, X0) / 1e6, min(r.end, X1) / 1e6
    axg.plot([s0, e0], [yy, yy], color=INK2, lw=2.2, solid_capstyle="butt")
    axg.text((s0 + e0) / 2, yy + 0.18, r.gene + (" →" if r.strand > 0 else " ←"), ha="center", fontsize=7, style="italic")
axg.set_xlim(X0 / 1e6, X1 / 1e6); axg.set_ylim(-0.1, 1.1); axg.set_yticks([])
axg.spines["left"].set_visible(False); axg.set_xlabel("Chromosome 10 position (Mb, GRCh38)")
hl = [plt.Line2D([], [], marker="o", ls="none", color=VIOL, ms=4, label=f"r² ≥ 0.8 with {LA}"),
      plt.Line2D([], [], marker="o", ls="none", color=AQUA, ms=4, label=f"r² ≥ 0.8 with {LB}"),
      plt.Line2D([], [], marker="D", ls="none", mfc=VIOL, mec=INK, ms=4.5, label=f"{LA} (AMD haplotype)"),
      plt.Line2D([], [], marker="^", ls="none", mfc=AQUA, mec=INK, ms=5, label=f"{LB}"),
      plt.Line2D([], [], marker="s", ls="none", mfc=AQUA, mec=INK, ms=4.5, label=f"{LC} (proxy of {LB})")]
fig.legend(handles=hl, loc="center", bbox_to_anchor=(0.31, F(146)), ncol=2, fontsize=7, handletextpad=0.2, columnspacing=1.2, labelspacing=0.3)
panel(fig, (0.005, F(217)), "a")

# ---------------- b: in-sample SuSiE PIPs for AMD: FinnGen (any, dry, wet) and UK Biobank (phecode 362.29)
sk = tsv(RR / "10q26/finngen_susie_key_variants.tsv")
uk = tsv(P / "results/revision3/panukb/panukb_susie.tsv"); uk = uk[uk.trait == "AMD_phe362.29"].set_index("rsid")
VAR = [("rs10490924", "rs10490924 (p.Ala69Ser)"), ("rs3750848", "rs3750848"), ("rs36212733", "rs36212733"), ("rs11200638", "rs11200638"),
       ("rs374446560", "rs374446560"), ("rs541533723", "rs541533723"), ("rs2672592", "rs2672592"), ("rs72631113", "rs72631113")]
PMIN = 1.1e-4
axb = fig.add_axes([0.745, F(166), 0.235, F(40)])
for i, (v, vl) in enumerate(VAR):
    y = len(VAR) - 1 - i
    for t, mk, fill in [("H7_AMD", "o", True), ("DRY_AMD", "v", False), ("WET_AMD", "^", False)]:
        r = sk[(sk.trait == t) & (sk.rsid == v)].iloc[0]
        cs = int(r.cs); col = VIOL if cs == 1 else (GREY if cs == 2 else AQUA)
        if t == "H7_AMD":
            axb.hlines(y + 0.14, 1e-4, max(r.PIP, 1e-4), color=col, lw=2.4, alpha=0.85)
        axb.plot(max(r.PIP, PMIN), y + 0.14, mk, ms=3.6, mfc=col if fill else "white", mec=col, mew=0.8, zorder=4)
    if v in uk.index:
        u = uk.loc[v]; col = VIOL if str(u.cs) not in ("nan", "NA") else AQUA
        axb.plot(max(u.pip, PMIN), y - 0.22, "s", ms=3.4, mfc=col, mec=INK, mew=0.6, zorder=5, clip_on=False)
    else:
        axb.text(PMIN * 1.1, y - 0.22, "UKB: n/a", fontsize=7, va="center")
axb.set_xscale("log"); axb.set_xlim(1e-4, 1.2); axb.set_xticks([1e-4, 1e-3, 1e-2, 1e-1, 1])
axb.set_xticklabels(["≤$10^{-4}$", "$10^{-3}$", "0.01", "0.1", "1"]); axb.minorticks_off()
axb.set_yticks(range(len(VAR))); axb.set_yticklabels([v[1] for v in VAR][::-1], fontsize=7)
axb.set_ylim(-0.6, len(VAR) - 0.4); axb.set_xlabel("PIP (in-sample SuSiE)")
axb.legend(handles=[plt.Line2D([], [], color=VIOL, lw=2.4, label="CS 1"), plt.Line2D([], [], color=GREY, lw=2.4, label="CS 2"),
                    plt.Line2D([], [], color=AQUA, lw=2.4, label="No CS"),
                    plt.Line2D([], [], marker="o", ls="none", color=INK, ms=3.4, label="FinnGen AMD"),
                    plt.Line2D([], [], marker="v", ls="none", mfc="white", mec=INK, ms=3.4, label="Dry / wet"),
                    plt.Line2D([], [], marker="s", ls="none", mfc="white", mec=INK, ms=3.4, label="UK Biobank")],
           loc="lower center", bbox_to_anchor=(0.30, 1.0), fontsize=7, handlelength=1.0, borderaxespad=0.1, labelspacing=0.15, ncol=3, columnspacing=0.5, handletextpad=0.3)
panel(fig, (0.565, F(217)), "b")

# ---------------- c: joint effects of the two haplotypes (lead | side vs side | lead), log OR and 95% CI
hm = tsv(P / "results/story/haplotype_trait_matrix.tsv")
def hrow(trait, ds):
    q = hm[(hm.trait == trait) & (hm.dataset == ds)]
    assert len(q), (trait, ds)
    return q.iloc[0]
TOPR = [("AMD · UK Biobank (in-sample LD)", hrow("AMD", "UKB_AMD_phe362.29"), "o")]
R6 = P / "results/revision6"
ht = tsv(R6 / "heterogeneity_main_table.tsv"); ht = ht[ht.subtype == "IS"].set_index("dataset")
def ht_row(ds):   # per-dataset joint effects (lead, side) from the round-6 heterogeneity table
    r = ht.loc[ds]
    return [(r.b_lead, r.b_lead - 1.96 * r.s_lead, r.b_lead + 1.96 * r.s_lead, r.p_lead),
            (r.b_side, r.b_side - 1.96 * r.s_side, r.b_side + 1.96 * r.s_side, r.p_side)]
def pooled(q, col):   # GLS random effects with overlap-modelled covariance; Hartung-Knapp CI
    r = ht.loc[f"POOLED main ({q})"]; return r[col], r.HK_lo, r.HK_hi, r.HK_p
cv = tsv(R6 / "cad_independent_validation.tsv", nrows=4).set_index("dataset")   # first block: per-dataset joint models
def cad_row(ds):
    r = cv.loc[ds]
    return [(r.b_lead, r.b_lead - 1.96 * r.s_lead, r.b_lead + 1.96 * r.s_lead, r.p_lead),
            (r.b_side, r.b_side - 1.96 * r.s_side, r.b_side + 1.96 * r.s_side, r.p_side)]
def std(r):   # (lead, side) as (est, lo, hi, p) from a haplotype-matrix row
    return [(r.beta_lead_joint, r.beta_lead_joint - 1.96 * r.se_lead_joint, r.beta_lead_joint + 1.96 * r.se_lead_joint, r.p_lead_joint),
            (r.beta_side_joint, r.beta_side_joint - 1.96 * r.se_side_joint, r.beta_side_joint + 1.96 * r.se_side_joint, r.p_side_joint)]
ROWS = [("CAD · UK Biobank (in-sample LD)", std(hrow("CAD (coronary atherosclerosis)", "UKB_CAD_phe411.4")), "o"),
        ("CAD · FinnGen I9_CORATHER (in-sample LD)", cad_row("FinnGen R12 I9_CORATHER"), "o"),
        ("CAD · CARDIoGRAMplusC4D 2015", cad_row("CARDIoGRAMplusC4D 2015 (Nikpay; 60,801 cases)"), "o"),
        ("CAD · Aragam et al. (incl. UKB; consistency)", std(hrow("CAD", "Aragam2022")), "o"),
        ("WMH (SD)", std(hrow("WMH", "WMH")), "o"),
        ("IS · MEGASTROKE", ht_row("MEGASTROKE"), "o"),
        ("IS · GIGASTROKE remainder", ht_row("GIGASTROKE_remainder"), "o"),
        ("Stroke excl. SAH · FinnGen", ht_row("FinnGen_I9_STR"), "o"),
        ("IS · UK Biobank", ht_row("UKB"), "o"),
        ("IS · pooled (4 EUR, overlap-modeled; HK CI)", [pooled("lead_rs3750848_given_side", "b_lead"), pooled("side_rs2672592_given_lead", "b_side")], "D"),
        ("Difference, pooled (HK CI)", [None, pooled("contrast_side_minus_lead", "b_contrast")], "D"),
        ("IS · GIGASTROKE EUR (full; reference)", std(hrow("IS", "GIGA_EUR_IS")), "o"),
        ("IS · GIGASTROKE EAS", std(hrow("IS", "GIGA_EAS_IS")), "o")]
for lab, vals, _ in ROWS:
    if "pooled" in lab or "Difference" in lab: print(lab, [None if v is None else tuple(round(float(x), 3) for x in v[:3]) for v in vals])
def hplot(ax, rows):
    for k, (lab, vals, mk) in enumerate(rows):
        y = len(rows) - 1 - k
        if vals[0] is None:   # single pooled difference
            e, lo, hi, p = vals[1]
            forest(ax, y, e, lo, hi, color=INK, marker=mk, ms=4.2, filled=p < 0.05, log=False)
            continue
        for (e, lo, hi, p), col, dy in zip(vals, (VIOL, AQUA), (0.17, -0.17)):
            forest(ax, y + dy, e, lo, hi, color=col, marker=mk, ms=4.2 if mk == "D" else 3.4, filled=p < 0.05, log=False)
    ax.axvline(0, color=GREY, ls="--", lw=0.6, zorder=0)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows][::-1], fontsize=7, linespacing=1.0)
    ax.set_ylim(-0.6, len(rows) - 0.4); ax.tick_params(axis="y", length=0)
ax1 = fig.add_axes([0.335, F(121), 0.245, F(6)]); hplot(ax1, [(l, std(r), m) for l, r, m in TOPR])
ax1.set_xlim(-0.05, 0.50); ax1.set_xticks([0, 0.2, 0.4])
ax2 = fig.add_axes([0.335, F(70), 0.245, F(46)]); hplot(ax2, ROWS)
ax2.set_xlim(-0.24, 0.14); ax2.set_xticks([-0.2, -0.1, 0, 0.1])
ax2.set_xlabel("Joint effect per allele, log OR (95% CI)")
ax1.legend(handles=[plt.Line2D([], [], marker="o", color=VIOL, ms=3.4, lw=0.8, label=f"AMD haplotype ({LA}-G | {LB})"),
                    plt.Line2D([], [], marker="o", color=AQUA, ms=3.4, lw=0.8, label=f"rs2672592 haplotype ({LB}-T | {LA})")],
           loc="lower left", bbox_to_anchor=(-0.95, 1.0), fontsize=7, ncol=1, handletextpad=0.3, borderaxespad=0, labelspacing=0.2)
panel(fig, (0.005, F(133)), "c")

# ---------------- e: haplotype structure only (no tissue or expression arrows)
ax = fig.add_axes([0.615, F(66), 0.385, F(68)]); ax.set_xlim(0, 65); ax.set_ylim(0, 70); ax.axis("off")
ax.text(1, 66, "Haplotype structure", fontsize=7.5, fontweight="bold", va="center")
XA, XB, X0, X1 = 21, 32, 15, 37
for yy, a1, lab in [(48, "risk", "AMD risk\nhaplotype"), (31, "non-risk", "Other\nhaplotype")]:
    ax.plot([X0, X1], [yy, yy], color=INK2, lw=1.8, solid_capstyle="round", zorder=1)
    ax.scatter([XA], [yy], s=80, c=VIOL if a1 == "risk" else LGREY, edgecolors=INK, linewidths=0.6, zorder=3)
    ax.scatter([XB], [yy], s=80, c=AQUA, edgecolors=INK, linewidths=0.6, zorder=3)
    ax.text(X0 - 1.2, yy, lab, fontsize=7, fontweight="bold", va="center", ha="right", linespacing=1.0)
    ax.text(XA, yy - 4.5, "G (Ser)" if a1 == "risk" else "T (Ala)", ha="center", fontsize=7, va="center")
    ax.text(XB, yy - 4.5, "T", ha="center", fontsize=7, va="center")
ax.text(XA, 59.5, "rs3750848\n(p.Ala69Ser)", ha="center", fontsize=7, va="center", linespacing=1.0)
ax.text(XB + 1.5, 59.5, "rs2672592", ha="left", fontsize=7, va="center")
ax.add_patch(FancyBboxPatch((41, 45.5), 22, 10, boxstyle="round,pad=0,rounding_size=1.2", fc="#f1eff9", ec=VIOL, lw=0.8))
ax.text(52, 50.5, "AMD ↑\n(in-sample LD)", ha="center", va="center", fontsize=6.5, linespacing=1.1)
ax.add_patch(FancyBboxPatch((41, 31.5), 22, 10, boxstyle="round,pad=0,rounding_size=1.2", fc="#eaf7f2", ec=AQUA, lw=0.8))
ax.text(52, 36.5, "CAD ↓\n(conditional)", ha="center", va="center", fontsize=6.5, linespacing=1.1)
ax.annotate("", xy=(40.5, 51.5), xytext=(XA, 50.4), arrowprops=dict(arrowstyle="-", color=VIOL, lw=0.8, ls=(0, (2, 2)),
                                                                    connectionstyle="arc3,rad=-0.25"))
ax.plot([XB + 1.4, 40.5], [47.2, 37.5], color=AQUA, lw=0.8, ls=(0, (2, 2)))
ax.plot([XB + 1.4, 40.5], [31.6, 35.5], color=AQUA, lw=0.8, ls=(0, (2, 2)))
ax.add_patch(FancyBboxPatch((40, 18), 24.5, 10, boxstyle="round,pad=0,rounding_size=1.2", fc="white", ec=INK, lw=0.6, ls=(0, (2, 1.5))))
ax.text(52.25, 23, "Ischemic stroke ↓:\nhaplotype unresolved", ha="center", va="center", fontsize=6.2, linespacing=1.1)
ax.text(1, 8, "AMD risk alleles occur almost only on\nrs2672592-T chromosomes (D′ 0.93–0.97)", fontsize=7, va="center")
panel(fig, (0.60, F(133)), "e")

# ---------------- d: tissue eQTL z-scores (incl. RPE single-nucleus eQTL)
eq = tsv(R / "tissue_eqtl/key_variant_eqtl_lookup.tsv")
DS = [("Voigt_RPE_sn", "RPE (snRNA, Voigt)"), ("EyeGEx_Retina", "Retina (EyeGEx)"), ("ROSMAP_DLPFC", "Brain DLPFC (ROSMAP)"), ("BrainSeq_DLPFC", "Brain DLPFC (BrainSeq)"),
      ("GTEx_Brain_Cortex", "Brain cortex (GTEx)"), ("GTEx_Brain_Putamen", "Brain putamen (GTEx)"),
      ("GTEx_Artery_Aorta", "Aorta (GTEx)"), ("GTEx_Artery_Coronary", "Coronary artery (GTEx)"), ("GTEx_Artery_Tibial", "Tibial artery (GTEx)"),
      ("GTEx_Fibroblast", "Fibroblasts (GTEx)"), ("GTEx_Whole_Blood", "Whole blood (GTEx)"), ("eQTLGen_Blood", "Whole blood (eQTLGen)")]
COLS = [("HTRA1", LA), ("HTRA1", LB), ("ARMS2", LA), ("ARMS2", LB)]
Z = np.full((len(COLS), len(DS)), np.nan); PV = Z.copy(); NS = set()
vo = tsv(P / "results" / "revision3" / "eye_eqtl" / "sc_eqtl_10q26_3p14.tsv")
vo = vo[(vo.cell_type == "RPE") & (vo.gene == "HTRA1") & (vo.rsid == LA)].iloc[0]
for j, (gene, snp) in enumerate(COLS):
    for i, (ds, _) in enumerate(DS):
        if ds == "Voigt_RPE_sn":   # single-nucleus eQTL; only FDR < 0.05 results are public (t per AMD risk allele)
            if gene == "HTRA1" and snp == LA:
                Z[j, i] = vo.t * (1 if vo.effect_allele == vo.amd_risk_allele else -1); PV[j, i] = vo.p
            else:
                NS.add((j, i))
            continue
        q = eq[(eq.region == "10q26_HTRA1") & (eq.gene == gene) & (eq.dataset == ds) & (eq.SNP == snp)]
        if not len(q) and snp == LA:
            q = eq[(eq.region == "10q26_HTRA1") & (eq.gene == gene) & (eq.dataset == ds) & (eq.SNP == "rs36212733")]
        if len(q):
            Z[j, i] = q.b_allele.iloc[0] / q.se.iloc[0]; PV[j, i] = q.p.iloc[0]
axd = fig.add_axes([0.20, 7 / H, 0.62, 26 / H])
cmap = LinearSegmentedColormap.from_list("expr", ["#4a3aa7", "#f0efec", "#008300"])
im = axd.imshow(np.ma.masked_invalid(np.clip(Z, -8, 8)), cmap=cmap, norm=TwoSlopeNorm(0, -8, 8), aspect="auto")
for j in range(Z.shape[0]):
    for i in range(Z.shape[1]):
        if (j, i) in NS:
            axd.text(i, j, "n.s.", ha="center", va="center", fontsize=7, color=GREY)
        elif np.isnan(Z[j, i]):
            axd.text(i, j, "n/a", ha="center", va="center", fontsize=7, color=GREY)
        else:
            axd.text(i, j, (f"{Z[j, i]:.1f}" + ("*" if PV[j, i] < 0.05 else "")).replace("-", "−"), ha="center", va="center",
                     fontsize=7, color="white" if abs(Z[j, i]) > 6 else INK)
axd.set_xticks(range(len(DS))); axd.set_xticklabels([d[1] for d in DS], fontsize=7)
axd.set_yticks(range(len(COLS))); axd.set_yticklabels([f"{g_} · {s_}-{'G' if s_ == LA else 'T'}" for g_, s_ in COLS], fontsize=7)
axd.tick_params(length=0); [s.set_visible(False) for s in axd.spines.values()]
axd.xaxis.tick_top(); plt.setp(axd.get_xticklabels(), rotation=35, ha="left")
cb = fig.add_axes([0.835, 7 / H, 0.012, 26 / H]); cbar = fig.colorbar(im, cax=cb)
cbar.set_label("eQTL z or t (AMD risk allele)", fontsize=7); cbar.ax.tick_params(labelsize=7, length=2); cbar.outline.set_linewidth(0.4)
panel(fig, (0.005, F(62)), "d")
save(fig, "Figure4")
