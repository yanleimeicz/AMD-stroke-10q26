from rev_common import *
from matplotlib.patches import Patch

H = 162
fig = plt.figure(figsize=(W, H * MM))
F = lambda mm: mm / H

# ---------------- pathway-partitioned MR (FinnGen AMD filled; IAMDGC open)
pw = tsv(P / "results/story/pathway_mr.tsv")
pw = pw[pw.group != "complement"]
ch = tsv(P / "results/revision5/complement_heterogeneity.tsv")
ch = ch[(ch.definition == "genebased") & (ch.analysis == "all complement instruments")].assign(group="complement")
pw = pd.concat([pw, ch[["exposure", "outcome", "group", "nsnp", "b", "se", "p"]]], ignore_index=True)
ce = tsv(P / "results/revision5/causal_expectation_respecified.tsv")
ce = ce[(ce.exposure == "AMD_FG") & (ce.outcome == "IS") & (ce.exclusion == "no_10q26") & (ce.prevalence == 0.063) & (ce.g_form == "exact")
        & ce.kappa_source.str.startswith("kappa (early AMD)")].iloc[0]
EXP_PT, EXP_LO, EXP_HI = ce.expected_OR_point, ce.expected_OR_q025, ce.expected_OR_q975
GRP = [("all", "All loci", INK), ("no_10q26", "Without 10q26", INK), ("10q26", "10q26 (ARMS2–HTRA1)", ANT),
       ("complement", "Complement (gene-based)", VIOL), ("lipid", "Lipid", AQUA), ("other", "Other", GREY)]
OUTC = [("AS", "Any\nstroke", (0.80, 1.20), True, [0.9, 1.0, 1.1]), ("IS", "Ischemic\nstroke", (0.80, 1.20), True, [0.9, 1.0, 1.1]),
        ("SVS", "Small-vessel\nstroke", (0.70, 1.30), True, [0.8, 1.0, 1.2]), ("CAD", "Coronary\nartery disease", (0.55, 1.15), False, [0.6, 0.8, 1.0]),
        ("WMH", "WMH volume\n(exp β)", (0.90, 1.10), False, [0.95, 1.0, 1.05])]
x0, wd, gp = 0.19, 0.118, 0.012
for j, (o, ol, (xl, xh), band, xt) in enumerate(OUTC):
    ax = fig.add_axes([x0 + j * (wd + gp), F(76), wd, F(62)])
    if band:
        ax.axvspan(EXP_LO, EXP_HI, color=ANT, alpha=0.15, lw=0, zorder=0)
        ax.axvline(EXP_PT, color=ANT, lw=0.6, zorder=0)
    for k, (g, gl, col) in enumerate(GRP):
        y = len(GRP) - 1 - k
        for ex, dy, fill in [("AMD_FG", 0.16, True), ("AMD_IAMDGC", -0.16, False)]:
            q = pw[(pw.exposure == ex) & (pw.outcome == o) & (pw.group == g)]
            if not len(q): continue
            r = q.iloc[0]
            est, lo, hi = np.exp(r.b), np.exp(r.b - 1.96 * r.se), np.exp(r.b + 1.96 * r.se)
            forest(ax, y + dy, est, max(lo, xl), min(hi, xh), color=col, filled=fill, ms=3.4, log=False)
            if hi > xh: ax.plot(xh, y + dy, ">", ms=3, color=col, clip_on=False, zorder=4)
            if lo < xl: ax.plot(xl, y + dy, "<", ms=3, color=col, clip_on=False, zorder=4)
    ax.axvline(1, color=GREY, ls="--", lw=0.6, zorder=1)
    ax.set_xlim(xl, xh); ax.set_ylim(-0.6, len(GRP) - 0.4)
    ax.set_yticks(range(len(GRP))); ax.set_yticklabels([g[1] for g in GRP][::-1] if j == 0 else [])
    ax.set_title(ol, fontsize=7, fontweight="bold", pad=3, linespacing=1.1)
    ax.tick_params(axis="x", labelsize=7); ax.set_xticks(xt); ax.set_xticklabels([f"{v:g}" for v in xt])
fig.text(x0 + 2.5 * (wd + gp) - gp / 2, F(64), "OR per unit log-odds of AMD liability (95% CI)", ha="center", fontsize=8)
fig.legend(handles=[plt.Line2D([], [], marker="o", color=INK, ms=3.4, lw=0.8, label="FinnGen AMD"),
                    plt.Line2D([], [], marker="o", color=INK, mfc="white", ms=3.4, lw=0.8, label="IAMDGC advanced AMD"),
                    Patch(fc=ANT, alpha=0.15, ec="none", label="Expected if the NHANES association were causal (scaled to FinnGen AMD liability; 95% interval)")],
           loc="center", bbox_to_anchor=(0.55, F(153)), ncol=2, fontsize=7, handletextpad=0.4, columnspacing=1.0)
panel(fig, (0.005, F(156)), "a")

# ---------------- b: per-variant Wald ratios of complement variants on IS and CAD
cv = tsv(P / "results/revision5/complement_per_variant.tsv")
VARS = [("rs1329424", "CFH lead (rs1329424)"), ("rs1061170", "CFH Y402H (rs1061170)"), ("rs6677604", "CFHR1/3-deletion tag (rs6677604)"),
        ("rs2230199", "C3 p.Arg102Gly (rs2230199)"), ("rs641153", "CFB p.Arg32Gln (rs641153)"), ("rs429608", "C2/CFB lead (rs429608)")]
for j, (o, ol, xl, xh) in enumerate([("IS", "Ischemic stroke", 0.85, 1.25), ("CAD", "Coronary artery disease", 0.93, 1.12)]):
    ax = fig.add_axes([0.30 + j * 0.36, F(15), 0.25, F(36)])
    for k, (snp, lab) in enumerate(VARS):
        r = cv[(cv.SNP == snp) & (cv.outcome == o)].iloc[0]; y = len(VARS) - 1 - k
        est, lo, hi = np.exp(r.wald), np.exp(r.wald - 1.96 * r.se_wald), np.exp(r.wald + 1.96 * r.se_wald)
        forest(ax, y, est, max(lo, xl), min(hi, xh), color=VIOL, filled=r.p_out < 0.05, ms=3.6, log=False)
        ax.text(1.03, y, f"{est:.3f}" + ("*" if r.p_out < 0.05 else ""), transform=ax.get_yaxis_transform(), va="center", fontsize=7)
    ax.axvline(1, color=GREY, ls="--", lw=0.6, zorder=0)
    ax.set_yticks(range(len(VARS))); ax.set_yticklabels([v[1] for v in VARS][::-1] if j == 0 else [])
    ax.set_ylim(-0.6, len(VARS) - 0.4); ax.set_xlim(xl, xh)
    ax.set_title(ol, fontsize=7, fontweight="bold", pad=3)
    ax.set_xlabel("Wald ratio, OR per unit log-odds\nof AMD liability (95% CI)", fontsize=7)
fig.text(0.98, F(57), "Filled symbols and *, P < 0.05", ha="right", fontsize=7)
panel(fig, (0.005, F(57)), "b")
save(fig, "Fig5")
