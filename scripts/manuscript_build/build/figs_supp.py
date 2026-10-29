"""Revised supplementary figures S1-S6 -> submission/revision/figures/supp/FigS#.png and preview_supp_figures.pdf"""
import io, textwrap
from rev_common import *
from matplotlib.backends.backend_pdf import PdfPages
import fitz

SUP = OUT / "supp"; SUP.mkdir(parents=True, exist_ok=True)
A4 = (210 * MM, 297 * MM)


def load(path, page=0, zoom=4):
    path = Path(path)
    if path.suffix.lower() == ".pdf":
        doc = fitz.open(path); pix = doc[page].get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=False)
        return Image.open(io.BytesIO(pix.tobytes("png"))).convert("RGB")
    return Image.open(path).convert("RGB")


def to_img(f, dpi=300):
    grey_text_to_black(f)
    buf = io.BytesIO(); f.savefig(buf, format="png", dpi=dpi); plt.close(f); buf.seek(0)
    return Image.open(buf).convert("RGB")


def save_supp(name, fig):
    grey_text_to_black(fig)
    fig.savefig(SUP / f"{name}.png", dpi=300); plt.close(fig); print("saved", name, flush=True)


def compose(name, ims, ncol=1, letters=True, h=None):
    n = len(ims); nrow = int(np.ceil(n / ncol))
    fig = plt.figure(figsize=(W, (h or 90) * nrow * MM))
    for k, im in enumerate(ims):
        r_, c_ = divmod(k, ncol)
        ax = fig.add_axes([c_ / ncol + 0.02, 1 - (r_ + 1) / nrow + 0.01, 1 / ncol - 0.03, 1 / nrow - 0.03]); ax.imshow(im); ax.axis("off")
        if letters and n > 1:
            fig.text(c_ / ncol + 0.005, 1 - r_ / nrow - 0.025, "abcdefgh"[k], fontsize=10, fontweight="bold")
    save_supp(name, fig)


# ---------------- S1: NHANES flow + subgroups (unchanged content)
def flow_panel():
    fl = tsv(R / "nhanes/flowchart.tsv")
    f, ax = plt.subplots(figsize=(80 * MM, 95 * MM)); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, len(fl) * 2)
    for i, r in fl.iterrows():
        y = (len(fl) - 1 - i) * 2 + 0.3
        ax.add_patch(plt.Rectangle((0.3, y), 9.4, 1.4, fc="#eef3fb", ec=CON, lw=0.8))
        ax.text(5, y + 0.7, f"{r.step.replace('>=', '≥')}\nn = {int(r.n):,}", ha="center", va="center", fontsize=7)
        if i < len(fl) - 1:
            ax.annotate("", xy=(5, y - 0.55), xytext=(5, y), arrowprops=dict(arrowstyle="-|>", lw=0.8, color=INK2, mutation_scale=7))
    return to_img(f)
compose("FigS1", [flow_panel(), load(R / "nhanes/Fig_forest_subgroups.png")], ncol=1, h=95)

# ---------------- S2: MR with distance-defined instruments: all methods (a), leave-one-locus-out (b), locus contributions (c)
mr = tsv(RR / "mr/mr_revised_all.tsv"); mr = mr[(mr.scheme == "A_dist") & (mr.exclusion == "main")]
EXP = [("AMD_FG", "AMD (FinnGen)"), ("DRYAMD_FG", "Dry AMD"), ("WETAMD_FG", "Wet AMD"), ("AMD_IAMDGC", "AMD (IAMDGC)")]
OUTC = [("AS", "Any stroke"), ("IS", "Ischemic"), ("LAS", "Large-artery"), ("CES", "Cardioembolic"), ("SVS", "Small-vessel"), ("WMH", "WMH (exp β)")]
METH = [("IVW", "IVW"), ("Egger", "MR-Egger"), ("WeightedMedian", "Weighted median"), ("WeightedMode", "Weighted mode"), ("RAPS", "MR-RAPS")]
fig = plt.figure(figsize=(W, 225 * MM))
for i, (e, el) in enumerate(EXP):
    for j, (o, ol) in enumerate(OUTC):
        ax = fig.add_axes([0.14 + j * 0.142, 0.885 - i * 0.072, 0.125, 0.062])
        for k, (m, ml) in enumerate(METH):
            r = mr[(mr.exposure == e) & (mr.outcome == o) & (mr.method == m)]
            if not len(r): continue
            r = r.iloc[0]; y = len(METH) - 1 - k
            forest(ax, y, np.exp(r.b), np.exp(r.b - 1.96 * r.se), np.exp(r.b + 1.96 * r.se), color=INK if m == "IVW" else INK2,
                   filled=r.p < 0.05, ms=2.6, log=False, lw=0.6)
        ax.axvline(1, color=GREY, ls="--", lw=0.5, zorder=0)
        ax.set_xlim(0.75, 1.25); ax.set_xticks([0.8, 1.0, 1.2]); ax.set_ylim(-0.7, len(METH) - 0.3)
        ax.set_yticks(range(len(METH))); ax.set_yticklabels([m[1] for m in METH][::-1] if j == 0 else [], fontsize=7)
        ax.tick_params(labelbottom=(i == len(EXP) - 1))
        if i == 0: ax.set_title(ol, fontsize=7, fontweight="bold", pad=2)
        if j == len(OUTC) - 1: ax.text(1.04, 0.5, el, transform=ax.transAxes, rotation=270, va="center", fontsize=7, fontweight="bold")
fig.text(0.56, 0.625, "OR per unit log-odds of AMD liability (95% CI); filled symbols, P < 0.05", ha="center", fontsize=7)
fig.text(0.005, 0.965, "a", fontsize=10, fontweight="bold")
lo = tsv(RR / "mr/mr_leave_one_locus_out_wide.tsv"); lo = lo[lo.scheme == "A_dist"]
EXC = [("main", "All loci"), ("no_10q26", "−10q26"), ("no_ADAMTS9", "−ADAMTS9"), ("no_CFH", "−CFH"), ("no_C3", "−C3"),
       ("no_APOE", "−APOE"), ("no_all", "−all five")]
from statistics import NormalDist
def ci(b, p):
    z = NormalDist().inv_cdf(1 - max(p, 1e-300) / 2); se = abs(b) / z if z > 0 else np.nan
    return np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se)
for j, (o, ol) in enumerate([("AS", "Any stroke"), ("IS", "Ischemic stroke"), ("SVS", "Small-vessel stroke")]):
    ax = fig.add_axes([0.14 + j * 0.29, 0.355, 0.25, 0.19])
    for k, (e, el) in enumerate(EXP):
        r = lo[(lo.exposure == e) & (lo.outcome == o)].iloc[0]
        for x, (ex, xl) in enumerate(EXC):
            est, l_, h_ = ci(r[f"b_{ex}"], r[f"p_{ex}"])
            col = [INK, CON, AQUA, VIOL][k]
            xx = x + (k - 1.5) * 0.15
            ax.vlines(xx, l_, h_, color=col, lw=0.6)
            ax.plot(xx, est, "o", ms=2.6, mfc=col if r[f"p_{ex}"] < 0.05 else "white", mec=col, mew=0.6)
    ax.axhline(1, color=GREY, ls="--", lw=0.5, zorder=0)
    ax.set_xticks(range(len(EXC))); ax.set_xticklabels([x[1] for x in EXC], rotation=45, ha="right", fontsize=7)
    ax.set_title(ol, fontsize=7, fontweight="bold", pad=2); ax.set_ylim(0.82, 1.24)
    if j == 0: ax.set_ylabel("OR (95% CI)", fontsize=7)
fig.legend(handles=[plt.Line2D([], [], marker="o", ls="none", color=c, ms=3, label=l) for c, (_, l) in zip([INK, CON, AQUA, VIOL], EXP)],
           loc="center", bbox_to_anchor=(0.56, 0.585), ncol=4, fontsize=7, handletextpad=0.3)
fig.text(0.005, 0.59, "b", fontsize=10, fontweight="bold")
lc = tsv(RR / "mr/mr_locus_contributions.tsv"); lc = lc[(lc.exposure == "AMD_FG") & (lc.scheme == "A_dist")]
for j, (o, ol) in enumerate([("IS", "AMD (FinnGen) → ischemic stroke"), ("SVS", "AMD (FinnGen) → small-vessel stroke")]):
    ax = fig.add_axes([0.11 + j * 0.47, 0.05, 0.22, 0.19])
    d = lc[lc.outcome == o].sort_values("contrib")
    cols = [ANT if c < 0 else CON for c in d.contrib]
    ax.barh(range(len(d)), d.contrib, color=cols, height=0.65)
    for k, (_, r) in enumerate(d.iterrows()):
        ax.text(1.03, k, f"{int(r.n)} SNP{'s' if r.n > 1 else ''}, {r.weight_share * 100:.0f}%", transform=ax.get_yaxis_transform(),
                va="center", ha="left", fontsize=7, clip_on=False)
    ax.axvline(0, color=INK2, lw=0.6)
    ax.set_yticks(range(len(d))); ax.set_yticklabels(d.locus, fontsize=7)
    ax.set_xlabel("Contribution to IVW log-OR", fontsize=7); ax.set_title(ol, fontsize=7, fontweight="bold", pad=4, loc="left")
    xl_ = max(abs(d.contrib).max() * 1.15, 0.012); ax.set_xlim(-xl_, xl_)
fig.text(0.005, 0.275, "c", fontsize=10, fontweight="bold")
save_supp("FigS2", fig)

# ---------------- S3: GWAS-GWAS colocalization heatmap (a, redrawn) + prior sensitivity (b)
from matplotlib.colors import LinearSegmentedColormap as _LSC
cl = tsv(R / "coloc/coloc_gwas_gwas.tsv")
AMDS = ["AMD_FG", "DRYAMD_FG", "WETAMD_FG", "AMD_IAMDGC"]; STK = ["AS", "IS", "LAS", "CES", "SVS"]
cl["amd"] = np.where(cl.source.str.contains("AMD"), cl.source, cl.partner); cl["stk"] = np.where(cl.source.str.contains("AMD"), cl.partner, cl.source)
cl = cl[cl.amd.isin(AMDS) & cl.stk.isin(STK)]
cl["bin"] = (cl.pos38 // 2e6).astype(int)
best = cl.groupby(["chr", "bin", "amd", "stk"])["PP.H4"].max().reset_index()
keep = best[best["PP.H4"] > 0.5][["chr", "bin"]].drop_duplicates()
GENE = {(10, 61): "10q26 ARMS2–HTRA1", (3, 32): "3p14 ADAMTS9", (9, 66): "9q34 ABO", (7, 46): "7q21 CDK6", (1, 97): "1q31", (19, 5): "19p13 SMARCA4–LDLR"}
loci = [k for k in GENE if ((keep.chr == k[0]) & (keep.bin == k[1])).any()]
extra = [tuple(x) for x in keep.values.tolist() if tuple(x) not in GENE]
loci += extra
Mh = np.full((len(loci), len(AMDS) * len(STK)), np.nan)
for i, (c, bn) in enumerate(loci):
    for a_, am in enumerate(AMDS):
        for s_i, st in enumerate(STK):
            q = best[(best.chr == c) & (best.bin == bn) & (best.amd == am) & (best.stk == st)]
            if len(q): Mh[i, a_ * len(STK) + s_i] = q["PP.H4"].iloc[0]
f = plt.figure(figsize=(W, 70 * MM))
ax = f.add_axes([0.20, 0.12, 0.72, 0.62])
cmh = _LSC.from_list("h4", ["#f4f4f2", "#86b6ef", "#1c5cab", "#0d366b"])
ax.imshow(np.ma.masked_invalid(Mh), cmap=cmh, vmin=0, vmax=1, aspect="auto")
for i in range(Mh.shape[0]):
    for j in range(Mh.shape[1]):
        if not np.isnan(Mh[i, j]):
            ax.text(j, i, f"{Mh[i, j]:.2f}"[1:] if Mh[i, j] < 1 else "1.0", ha="center", va="center", fontsize=7, color="white" if Mh[i, j] > 0.6 else INK)
ax.set_yticks(range(len(loci))); ax.set_yticklabels([GENE.get(l, f"chr{l[0]}") for l in loci], fontsize=7)
ax.set_xticks(range(Mh.shape[1])); ax.set_xticklabels(STK * len(AMDS), fontsize=7)
ax.tick_params(length=0); [sp.set_visible(False) for sp in ax.spines.values()]
for a_, am in enumerate(["AMD (FinnGen)", "Dry AMD (FinnGen)", "Wet AMD (FinnGen)", "AMD (IAMDGC)"]):
    ax.text(a_ * 5 + 2, -0.95, am, ha="center", va="bottom", fontsize=7, fontweight="bold")
    if a_: ax.axvline(a_ * 5 - 0.5, color="white", lw=2.5)
cbx = f.add_axes([0.935, 0.12, 0.012, 0.62]); cb = f.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, 1), cmap=cmh), cax=cbx)
cb.set_label("PP.H4", fontsize=7); cb.ax.tick_params(labelsize=7, length=1.5); cb.outline.set_linewidth(0.3)
s3a = to_img(f)
p12 = tsv(RR / "10q26/coloc_abf_p12_sensitivity.tsv"); p12 = p12[p12.amd == "H7_AMD"]
cp = tsv(RR / "adamts9/coloc_pairwise_3p14.tsv")
f = plt.figure(figsize=(W, 62 * MM))
ax0 = f.add_axes([0.08, 0.2, 0.25, 0.68]); ax1 = f.add_axes([0.55, 0.2, 0.25, 0.68])
TS = [("GIGA_EUR_SVS", "SVS", CON), ("GIGA_EUR_AS", "Any stroke", INK), ("LACUNAR", "Lacunar", VIOL), ("GIGA_EUR_IS", "IS", AQUA),
      ("GIGA_EUR_LAS", "LAS", GREY), ("GIGA_EAS_IS", "IS, East Asian", ANT)]
for t, tl, c in TS:
    d = p12[p12.stroke == t].sort_values("p12")
    ax0.plot(d.p12, d["PP.H4.abf"], "-o", color=c, ms=3, lw=0.9, label=tl)
for (t1, t2, tl, c) in [("AMD_FG", "T2D", "AMD–T2D", ANT), ("AMD_FG", "AORTA_ADAMTS9", "AMD–aortic eQTL", CON), ("AMD_FG", "CIMT", "AMD–cIMT", VIOL),
                        ("AMD_FG", "IS", "AMD–IS", AQUA), ("AMD_FG", "LAS", "AMD–LAS", GREY)]:
    d = cp[((cp.trait1 == t1) & (cp.trait2 == t2)) | ((cp.trait1 == t2) & (cp.trait2 == t1))].sort_values("p12")
    ax1.plot(d.p12, d["PP.H4"], "-o", color=c, ms=3, lw=0.9, label=tl)
for ax_, t in [(ax0, "10q26: AMD (FinnGen) with stroke traits"), (ax1, "3p14: AMD (FinnGen) with other traits")]:
    ax_.set_xscale("log"); ax_.set_xlabel("Prior probability of a shared variant (p12)"); ax_.set_ylabel("PP.H4")
    ax_.set_ylim(0, 1.02); ax_.set_title(t, fontsize=7, fontweight="bold")
    ax_.legend(fontsize=7, loc="center left", bbox_to_anchor=(1.02, 0.5), handlelength=1.4)
s3b = to_img(f)
compose("FigS3", [s3a, s3b], ncol=1, h=70)

# ---------------- S4: 12q24 (unchanged)
# ---------------- S4: 12q24 (SH2B3) regional plot, redrawn from the data of scripts/49_cojo_susie.R (rev_s4_export.R)
s4 = tsv(FD / "s4_region12q24.tsv")
S4P = [("DRYAMD_FG", "Dry AMD (FinnGen)"), ("IS", "Ischemic stroke (GIGASTROKE)"), ("AS", "Any stroke (GIGASTROKE)"), ("AMD_IAMDGC", "Advanced AMD (IAMDGC)")]
S4P += [(p_, "SH2B3 eQTL, " + p_.replace("SH2B3 ", "").replace("GTEx_Artery_Coronary", "coronary artery (GTEx)")
         .replace("GTEx_Brain_Hippocampus", "hippocampus (GTEx)").replace("eQTLGen_Blood", "whole blood (eQTLGen)")
         .replace("EyeGEx_Retina", "retina (EyeGEx)")) for p_ in s4.panel.unique() if p_.startswith("SH2B3")]
R2C = [(-1, 0.2, "#a6a6a6", "< 0.2"), (0.2, 0.4, "#86b3d1", "0.2–0.4"), (0.4, 0.6, "#2c7fb8", "0.4–0.6"), (0.6, 0.8, "#f39c12", "0.6–0.8"), (0.8, 1.01, "#d7191c", "≥ 0.8")]
f = plt.figure(figsize=(W, 205 * MM)); n4 = len(S4P); hh4 = 0.86 / n4
for i, (pk, pl) in enumerate(S4P):
    ax = f.add_axes([0.09, 0.955 - (i + 1) * hh4 + 0.008, 0.70, hh4 - 0.016])
    d = s4[s4.panel == pk].copy(); d["lp"] = -np.log10(d.p.clip(lower=1e-300)); d["r2"] = d.r2.fillna(-0.5)
    d = d.sort_values("r2")
    for lo_, hi_, c_, _ in R2C:
        q = d[(d.r2 >= lo_) & (d.r2 < hi_)] if lo_ > -1 else d[d.r2 < hi_]
        ax.scatter(q.pos38 / 1e6, q.lp, s=4, c=c_, edgecolors="none", zorder=2)
    for snp in ["rs3184504", "rs587914"]:
        q = d[d.SNP == snp]
        if len(q):
            ax.scatter(q.pos38 / 1e6, q.lp, s=28, marker="D", c=VIOL, edgecolors=INK, linewidths=0.6, zorder=4)
            ax.annotate(snp, (q.pos38.iloc[0] / 1e6, q.lp.iloc[0]), xytext=(-4, 4) if snp == "rs3184504" else (4, 4),
                        textcoords="offset points", ha="right" if snp == "rs3184504" else "left", fontsize=7)
    ax.set_xlim(111.1, 111.8); ax.set_ylim(0, d.lp.max() * 1.25)
    ax.text(0.005, 0.97, pl, transform=ax.transAxes, va="top", fontsize=7, fontweight="bold")
    ax.set_ylabel("−log$_{10}$P", fontsize=7); ax.tick_params(labelbottom=(i == n4 - 1))
ax.set_xlabel("Chromosome 12 position (Mb, GRCh38)")
f.legend(handles=[plt.Line2D([], [], marker="o", ls="none", color=c_, ms=4, label=l_) for _, _, c_, l_ in R2C] +
         [plt.Line2D([], [], marker="D", ls="none", mfc=VIOL, mec=INK, ms=4.5, label="rs3184504 / rs587914")],
         title="r² with rs3184504", loc="upper left", bbox_to_anchor=(0.80, 0.95), fontsize=7, title_fontsize=7)
save_supp("FigS4", f)

# ---------------- S5: SMR-HEIDI genes shared under the relaxed criterion (no drug annotation)
sm = tsv(R / "smr_heidi/smr_heidi_all.tsv")
AMDT = ["AMD_FG", "AMD_IAMDGC", "DRYAMD_FG", "WETAMD_FG"]; STKT = ["AS", "IS", "LAS", "CES", "SVS", "WMH"]
a = sm[sm.outcome.isin(AMDT) & (sm.p_fdr < 0.05) & sm.heidi_pass]
s_ = sm[sm.outcome.isin(STKT) & (sm.p_SMR < 0.05) & sm.heidi_pass]
genes = sorted(set(a.GeneSymbol) & set(s_.GeneSymbol))
rows = []
for gname in genes:
    ra = a[a.GeneSymbol == gname].sort_values("p_SMR").iloc[0]; rs_ = s_[s_.GeneSymbol == gname].sort_values("p_SMR").iloc[0]
    rows.append((gname, ra.outcome, ra.b_SMR, ra.se_SMR, rs_.outcome, rs_.b_SMR, rs_.se_SMR, rs_.p_fdr))
dd = pd.DataFrame(rows, columns=["g", "oa", "ba", "sa", "os", "bs", "ss", "fdr_s"])
dd["opp"] = np.sign(dd.ba) != np.sign(dd.bs); dd = dd.sort_values(["opp", "g"], ascending=[False, True]).reset_index(drop=True)
LAB = {"AMD_FG": "AMD", "AMD_IAMDGC": "AMD (IAMDGC)", "DRYAMD_FG": "dry AMD", "WETAMD_FG": "wet AMD", "AS": "AS", "IS": "IS",
       "LAS": "LAS", "CES": "CES", "SVS": "SVS", "WMH": "WMH"}
f = plt.figure(figsize=(W, (40 + 5.2 * len(dd)) * MM))
axL = f.add_axes([0.17, 0.08, 0.25, 0.86]); axR = f.add_axes([0.56, 0.08, 0.25, 0.86])
for k, r in dd.iterrows():
    y = len(dd) - 1 - k; col = ANT if r.opp else CON
    for ax_, b, s, o in [(axL, r.ba, r.sa, r.oa), (axR, r.bs, r.ss, r.os)]:
        ax_.hlines(y, np.exp(b - 1.96 * s), np.exp(b + 1.96 * s), color=col, lw=0.8)
        ax_.plot(np.exp(b), y, "o", ms=3.4, color=col)
        ax_.text(1.02, y, LAB[o] + (f"  FDR {r.fdr_s:.2f}" if ax_ is axR else ""), transform=ax_.get_yaxis_transform(), va="center", fontsize=7)
for ax_, t in [(axL, "AMD (FDR < 0.05, HEIDI passed)"), (axR, "Stroke or WMH (P < 0.05, HEIDI passed)")]:
    ax_.axvline(1, color=GREY, ls="--", lw=0.6, zorder=0); ax_.set_xscale("log"); ax_.set_ylim(-0.7, len(dd) - 0.3)
    ax_.set_title(t, fontsize=7, fontweight="bold"); ax_.set_xlabel("SMR OR per SD expression (95% CI)")
    ax_.minorticks_off(); ax_.set_xticks([0.25, 0.5, 1, 2]); ax_.set_xticklabels(["0.25", "0.5", "1", "2"])
axL.set_yticks(range(len(dd))); axL.set_yticklabels(list(dd.g)[::-1], fontsize=7, style="italic"); axR.set_yticks([])
f.legend(handles=[plt.Line2D([], [], marker="o", color=ANT, ms=3.4, label="Opposite directions"),
                  plt.Line2D([], [], marker="o", color=CON, ms=3.4, label="Same direction")], loc="lower center", ncol=2, fontsize=7,
         bbox_to_anchor=(0.5, 0.0))
save_supp("FigS5", f)
print("S5 genes:", len(dd), "opposite:", int(dd.opp.sum()))

# ---------------- S6: plasma cis-pQTL MR, proteins shared by AMD and stroke or WMH (redrawn in the style of S5)
pq = tsv(R / "pqtl_mr/pqtl_mr_shared.tsv"); pq = pq[pq.source == "UKB-PPP"]
rows = []
for gname, sub in pq.groupby("gene"):
    ra = sub[sub.group == "AMD"].sort_values("pval").iloc[0]; rs_ = sub[sub.group == "STROKE_WMH"].sort_values("pval").iloc[0]
    rows.append((gname, ra.outcome, ra.b, ra.se, ra.fdr, ra.b_dec, ra.se_dec, rs_.outcome, rs_.b, rs_.se, rs_.fdr, rs_.b_dec, rs_.se_dec, ra.tier))
pp = pd.DataFrame(rows, columns=["g", "oa", "ba", "sa", "fa", "bda", "sda", "os", "bs", "ss", "fs", "bds", "sds", "tier"])
pp["opp"] = np.sign(pp.ba) != np.sign(pp.bs)
pp["both"] = pp.tier.astype(str).str.contains("both")
pp = pp.sort_values(["both", "opp", "g"], ascending=[False, False, True]).reset_index(drop=True)
f = plt.figure(figsize=(W, (45 + 5.2 * len(pp)) * MM))
axL = f.add_axes([0.15, 0.11, 0.24, 0.83]); axR = f.add_axes([0.58, 0.11, 0.24, 0.83])
XL, XH = 0.25, 4.0
clip = lambda v: min(max(v, XL * 1.04), XH / 1.04)
LAB6 = dict(LAB, AMD_IAMDGC="IAMDGC")
for k, r in pp.iterrows():
    y = len(pp) - 1 - k; col = ANT if r.opp else CON
    for ax_, bb, ss_, fd, bd, sd, o in [(axL, r.ba, r.sa, r.fa, r.bda, r.sda, r.oa), (axR, r.bs, r.ss, r.fs, r.bds, r.sds, r.os)]:
        ax_.hlines(y + 0.12, clip(np.exp(bb - 1.96 * ss_)), clip(np.exp(bb + 1.96 * ss_)), color=col, lw=0.8)
        ax_.plot(clip(np.exp(bb)), y + 0.12, "o", ms=3.4, color=col)
        if pd.notna(bd) and pd.notna(sd):
            e_ = np.exp(bd); inside = XL < e_ < XH
            ax_.hlines(y - 0.22, clip(np.exp(bd - 1.96 * sd)), clip(np.exp(bd + 1.96 * sd)), color=GREY, lw=0.6)
            ax_.plot(clip(e_), y - 0.22, "D" if inside else ("<" if e_ <= XL else ">"), ms=2.8 if inside else 3.6,
                     mfc="white" if inside else INK2, mec=INK2, mew=0.6)
        ax_.text(1.02, y, f"{LAB6.get(o, o)}  FDR {fd:.1g}".replace("e-0", "e−").replace("e-", "e−"), transform=ax_.get_yaxis_transform(), va="center", fontsize=7)
for ax_, t in [(axL, "AMD (most significant outcome)"), (axR, "Stroke or WMH (most significant outcome)")]:
    ax_.axvline(1, color=GREY, ls="--", lw=0.6, zorder=0); ax_.set_xscale("log"); ax_.set_ylim(-0.7, len(pp) - 0.3)
    ax_.set_title(t, fontsize=7, fontweight="bold"); ax_.set_xlabel("OR per SD of plasma protein (95% CI)"); ax_.minorticks_off()
    ax_.set_xlim(XL, XH); ax_.set_xticks([0.25, 0.5, 1, 2, 4]); ax_.set_xticklabels(["0.25", "0.5", "1", "2", "4"])
axL.set_yticks(range(len(pp))); axL.set_yticklabels([g_ + (" †" if b_ else "") for g_, b_ in zip(pp.g, pp.both)][::-1], fontsize=7, style="italic")
axR.set_yticks([])
f.legend(handles=[plt.Line2D([], [], marker="o", color=ANT, ms=3.4, label="Opposite directions (UKB-PPP)"),
                  plt.Line2D([], [], marker="o", color=CON, ms=3.4, label="Same direction (UKB-PPP)"),
                  plt.Line2D([], [], marker="D", color=GREY, mfc="white", mec=INK2, ms=3, label="deCODE (triangle: beyond axis range)")],
         loc="lower center", ncol=3, fontsize=7, bbox_to_anchor=(0.5, 0.0))
save_supp("FigS6", f)
print("S6 proteins:", len(pp), "both:", int(pp.both.sum()))

# ---------------- S7: single-cell expression of ADAMTS9 / ADAMTS9-AS2 (former Fig. 5d; unchanged content)
from matplotlib.colors import LinearSegmentedColormap
DS7 = [("rpe_choroid", "RPE–choroid (n cells)"), ("brain_vascular_multiomics", "Brain vasculature (n cells)")]
rows7 = []
for ds, dl in DS7:
    t = tsv(R / f"adamts9/scrna/{ds}_celltype_expr.tsv")
    t = t[t.n_cells >= 200].sort_values("ADAMTS9_frac", ascending=False).head(9 if ds == "rpe_choroid" else 8)
    for _, r in t.iterrows():
        rows7.append((dl, f"{r.cell_type} ({int(r.n_cells):,})", r.ADAMTS9_mean, r.ADAMTS9_frac, r["ADAMTS9-AS2_mean"], r["ADAMTS9-AS2_frac"]))
cdf = pd.DataFrame(rows7, columns=["ds", "ct", "m1", "f1", "m2", "f2"])
f = plt.figure(figsize=(120 * MM, 112 * MM))
axd = f.add_axes([0.53, 0.13, 0.17, 0.83])
seq2 = LinearSegmentedColormap.from_list("seq2", ["#cde2fb", "#3987e5", "#0d366b"])
ys = []; y = 0; labels = []; seps = []
for ds_, sub in list(cdf.groupby("ds", sort=False))[::-1]:
    for _, r in sub.iloc[::-1].iterrows():
        for x, m, fr in [(0, r.m1, r.f1), (1, r.m2, r.f2)]:
            axd.scatter(x, y, s=4 + 70 * fr, c=[m], cmap=seq2, vmin=0, vmax=2.3, edgecolors=INK2, linewidths=0.3, zorder=3)
        labels.append(r.ct); ys.append(y); y += 1
    seps.append((y - 0.5, ds_)); y += 1.7
axd.set_yticks(ys); axd.set_yticklabels(labels, fontsize=7)
for yb, dl in seps:
    axd.text(-0.72, yb + 0.25, dl, fontsize=7, fontweight="bold", ha="right", va="bottom")
axd.set_xticks([0, 1]); axd.set_xticklabels(["ADAMTS9", "ADAMTS9-AS2"], style="italic", rotation=30, ha="right"); axd.set_xlim(-0.6, 1.6)
axd.set_ylim(-0.8, y - 0.3); axd.tick_params(axis="y", length=0)
cax = f.add_axes([0.80, 0.62, 0.025, 0.22])
cb = f.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, 2.3), cmap=seq2), cax=cax)
cb.set_label("Mean expression", fontsize=7); cb.ax.tick_params(labelsize=7, length=1.5); cb.outline.set_linewidth(0.3)
f.text(0.78, 0.53, "Expressing cells", fontsize=7)
for fr, yv in [(0.25, 0.48), (0.5, 0.42), (0.75, 0.36)]:
    f.add_artist(plt.Line2D([0.80], [yv], marker="o", ms=np.sqrt(4 + 70 * fr), mfc=LGREY, mec=INK2, mew=0.3, transform=f.transFigure))
    f.text(0.83, yv, f"{int(fr * 100)}%", va="center", fontsize=7)
save_supp("FigS7", f)

# ---------------- preview PDF with legends
LEG = {}
txt = open(P / "submission/revision/figure_legends_final.md", encoding="utf-8").read()
for block in txt.split("\n### ")[1:]:
    head, *body = block.split("\n")
    if head.startswith("Figure S"):
        key = head.split()[1].rstrip(".")
        LEG[key] = (head.split(". ", 1)[1] if ". " in head else head, " ".join(l.strip() for l in body if l.strip()).replace("**", ""))
with PdfPages(P / "submission/revision/preview_supp_figures.pdf") as pdf:
    for k in ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11"]:
        im = Image.open(SUP / f"Fig{k}.png")
        title, legend = LEG.get(k, (k, ""))
        fig = plt.figure(figsize=A4)
        fig.text(0.06, 0.965, f"Figure {k}", fontsize=10, fontweight="bold")
        fig.text(0.06, 0.945, title, fontsize=8.5, fontweight="bold", va="top", family=["Arial", "DejaVu Sans"])
        w_, h_ = im.size; avail_h = 0.62; ar = (h_ / w_) * (0.88 * A4[0]) / A4[1]
        hh = min(avail_h, ar); ww = 0.88 * hh / ar if ar > 0 else 0.88
        ax = fig.add_axes([0.06 + (0.88 - ww) / 2, 0.93 - hh, ww, hh]); ax.imshow(im); ax.axis("off")
        fig.text(0.06, 0.93 - hh - 0.015, "\n".join(textwrap.wrap(legend, 125)), fontsize=7.5, va="top", linespacing=1.35, family=["Arial", "DejaVu Sans"])
        grey_text_to_black(fig); pdf.savefig(fig); plt.close(fig)
print("preview done")
