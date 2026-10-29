from rev_common import *
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap

fig = plt.figure(figsize=(W, 215 * MM))

# ---------------- a: NHANES (unchanged)
mm_ = tsv(R / "nhanes/main_models.tsv"); se_ = tsv(R / "nhanes/sensitivity.tsv")
rows = []
lab = {"Model 1": "M1: unadjusted", "Model 2": "M2: + age, sex, race/ethnicity",
       "Model 3": "M3: + education, PIR, smoking, alcohol, BMI", "Model 4": "M4: + hypertension, diabetes, lipids, CVD, CKD"}
for m in ["Model 1", "Model 2", "Model 3", "Model 4"]:
    r = mm_[(mm_.model == m) & (mm_.term == "amd_any")].iloc[0]
    rows.append((lab[m], r.OR, r.lo, r.hi, r.p, "main"))
sens_lab = {"S3 Excluding retinopathy": "Excluding retinopathy (M3)", "S4 Age >=50 years": "Age ≥50 years (M3)",
            "S5 Full model + log CRP": "M4 + log CRP", "S7 Multiple imputation (m=20), Model 4": "Multiple imputation (M4)",
            "S2 Excluding diabetes": "Excluding diabetes (M3)", "S1 Self-reported AMD (questionnaire)": "Self-reported AMD (M3)"}
for k, v in sens_lab.items():
    r = se_[se_.model == k].iloc[0]
    rows.append((v, r.OR, r.lo, r.hi, r.p, "sens"))
axa = fig.add_axes([0.36, 0.735, 0.40, 0.235])
n = len(rows); ys = np.arange(n)[::-1].astype(float); ys[4:] -= 1.0
for (l, o, lo, hi, p, kind), y in zip(rows, ys):
    forest(axa, y, o, lo, hi, color=INK if kind == "main" else INK2, filled=(p < 0.05))
    axa.text(1.03, y, f"{o:.2f} ({lo:.2f}–{hi:.2f})", transform=axa.get_yaxis_transform(), va="center", fontsize=7)
axa.axvline(1, color=GREY, lw=0.6, ls="--", zorder=0)
axa.set_yticks(ys); axa.set_yticklabels([r[0] for r in rows], fontsize=7); axa.set_ylim(-1.3, n - 0.2)
axa.set_xlim(0.5, 9); axa.set_xticks([0.5, 1, 2, 4, 8]); axa.set_xticklabels(["0.5", "1", "2", "4", "8"])
axa.set_xlabel("OR for stroke, any AMD vs none (95% CI)")
axa.text(1.03, ys[0] + 0.85, "OR (95% CI)", transform=axa.get_yaxis_transform(), fontsize=7, fontweight="bold")
axa.text(-0.86, ys[4] + 0.95, "Sensitivity analyses", transform=axa.get_yaxis_transform(), fontsize=7, color=INK2, fontweight="bold")
axa.text(-0.86, ys[0] + 0.85, "Adjustment models (any AMD)", transform=axa.get_yaxis_transform(), fontsize=7, color=INK2, fontweight="bold")
panel(fig, (0.005, 0.978), "a")

# ---------------- b: distance-defined instruments, IVW; all / without 10q26 / without CFH
mr = tsv(RR / "mr/mr_revised_all.tsv")
ce = tsv(P / "results/revision5/causal_expectation_respecified.tsv")
ce = ce[(ce.exposure == "AMD_FG") & (ce.outcome == "IS") & (ce.exclusion == "no_10q26") & (ce.prevalence == 0.063) & (ce.g_form == "exact")
        & ce.kappa_source.str.startswith("kappa (early AMD)")].iloc[0]
EXP_PT, EXP_LO, EXP_HI = ce.expected_OR_point, ce.expected_OR_q025, ce.expected_OR_q975
print(f"expected OR {EXP_PT:.3f} ({EXP_LO:.3f}-{EXP_HI:.3f}), kappa {ce.kappa:.2f}")
mr = mr[(mr.scheme == "A_dist") & (mr.method == "IVW")]
fdr = tsv(RR / "mr/fdr_all_layers.tsv"); fdr = fdr[fdr.layer == "MR forward (distance-pruned instruments)"].set_index("item").fdr
EXP = [("AMD_FG", "AMD (FinnGen)"), ("DRYAMD_FG", "Dry AMD (FinnGen)"), ("WETAMD_FG", "Wet AMD (FinnGen)"), ("AMD_IAMDGC", "Advanced AMD (IAMDGC)")]
OUTC = [("AS", "Any stroke"), ("IS", "Ischemic stroke"), ("LAS", "Large-artery stroke"), ("CES", "Cardioembolic stroke"),
        ("SVS", "Small-vessel stroke"), ("WMH", "WMH volume (exp β)")]
SETS = [("main", INK, "o", 0.29), ("no_10q26", AQUA, "s", 0.10), ("no_CFH", VIOL, "^", -0.10), ("no_major3", GREY, "D", -0.29)]
x0, wb = 0.155, 0.198
for j, (e, el) in enumerate(EXP):
    ax = fig.add_axes([x0 + j * (wb + 0.014), 0.30, wb, 0.325])
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
    # expected MR effect if the NHANES association were causal, kappa-scaled to FinnGen AMD liability (p = 0.063; 95% interval)
    ax.fill_betweenx([0.5, len(OUTC) - 0.4], EXP_LO, EXP_HI, color=ANT, alpha=0.15, lw=0, zorder=0)
    ax.vlines(EXP_PT, 0.5, len(OUTC) - 0.4, color=ANT, lw=0.6, zorder=0)
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
fig.legend(handles=hl, loc="lower center", ncol=5, bbox_to_anchor=(0.55, 0.667), handletextpad=0.3, columnspacing=0.9, fontsize=7)
from matplotlib.patches import Patch
fig.legend(handles=[Patch(fc=ANT, alpha=0.15, ec="none", label="Expected if the NHANES association were causal (scaled to FinnGen AMD liability; 95% interval)")],
           loc="center", bbox_to_anchor=(0.55, 0.658), fontsize=7, handlelength=1.6, handletextpad=0.4, borderaxespad=0, borderpad=0)
fig.text(0.55, 0.646, "FDR*, FDR < 0.05 (all loci)",
         ha="center", va="center", fontsize=7, color=INK)
panel(fig, (0.005, 0.655), "b")

# ---------------- c: positive controls (same MR pipeline)
pc = tsv(P / "results" / "story" / "positive_controls.tsv")
af = tsv(P / "results/revision5/positive_control_af.tsv"); af = af[(af.set == "all") & af.outcome.isin(["IS", "CES"])].copy()
af["exposure"] = "AF"
pc = pd.concat([pc, af[["exposure", "outcome", "b", "se", "p", "OR", "lo", "hi"]]], ignore_index=True)
SUBS = [("AMD instruments and atrial fibrillation liability", "OR per unit log-odds of liability",
         [(("AMD_FG", "AMD_IAMDGC"), "FinnGen AMD instruments → IAMDGC AMD"), (("AMD_IAMDGC", "AMD_FG"), "IAMDGC instruments → FinnGen AMD"),
          (("AF", "IS"), "Atrial fibrillation → ischemic stroke"), (("AF", "CES"), "Atrial fibrillation → cardioembolic stroke")],
         [0.31, 0.05, 0.18, 0.12], (0.8, 6.5), True),
        ("Systolic blood pressure (436 instruments)", "OR per mmHg",
         [(("SBP", "IS"), "SBP → ischemic stroke"), (("SBP", "SVS"), "SBP → small-vessel stroke"), (("SBP", "CAD"), "SBP → coronary artery disease")],
         [0.76, 0.05, 0.18, 0.12], (0.995, 1.06), False)]
for title, xl, items, pos, xlim, logx in SUBS:
    ax = fig.add_axes(pos)
    for k, ((e, o), lab) in enumerate(items):
        r = pc[(pc.exposure == e) & (pc.outcome == o)].iloc[0]; y = len(items) - 1 - k
        forest(ax, y, r.OR, r.lo, r.hi, color=INK, filled=r.p < 0.05, ms=3.6, log=logx)
        ax.text(1.03, y, f"{r.OR:.2f}" if r.OR > 1.2 else f"{r.OR:.3f}", transform=ax.get_yaxis_transform(), va="center", fontsize=7)
    ax.axvline(1, color=GREY, ls="--", lw=0.6, zorder=0)
    ax.set_yticks(range(len(items))); ax.set_yticklabels([it[1] for it in items][::-1], fontsize=7)
    ax.set_ylim(-0.6, len(items) - 0.4); ax.set_xlim(*xlim)
    if logx:
        ax.set_xticks([1, 2, 4]); ax.set_xticklabels(["1", "2", "4"]); ax.minorticks_off()
    else:
        ax.set_xticks([1.0, 1.02, 1.04, 1.06]); ax.set_xticklabels(["1.00", "1.02", "1.04", "1.06"])
    ax.set_xlabel(xl, fontsize=7); ax.set_title(title, fontsize=7, fontweight="bold", pad=3)
panel(fig, (0.005, 0.205), "c")
save(fig, "Fig2")
