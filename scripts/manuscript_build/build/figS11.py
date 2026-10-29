"""Figure S11 (former Fig. 4c): symmetric conditional analyses at 10q26 in each AMD, stroke and coronary dataset."""
from rev_common import *

fig = plt.figure(figsize=(120 * MM, 130 * MM))
# ---------------- c: symmetric conditional tests (same models for AMD and stroke; matched LD; bars Δr ± 0.03)
R2 = P / "results" / "revision2"
R3 = P / "results" / "revision3" / "panukb"
sc = pd.concat([tsv(R2 / "10q26_symmetric_conditional.tsv"), tsv(R3 / "panukb_symmetric_conditional.tsv")], ignore_index=True)
sc = sc[sc.ld_matched.astype(str) == "True"].set_index("trait")
rs = pd.concat([tsv(R2 / "10q26_r_sensitivity.tsv"), tsv(R3 / "panukb_r_sensitivity.tsv")], ignore_index=True)
rs = rs[(rs["mode"] == "both_lead_side") & (rs.delta.abs() <= 0.0301)]
GRP = [("FinnGen AMD (FinnGen LD)", [("FinnGen_H7_AMD", "AMD"), ("FinnGen_DRY_AMD", "Dry AMD"), ("FinnGen_WET_AMD", "Wet AMD")]),
       ("UK Biobank AMD (in-sample LD)", [("AMD_phe362.29", "AMD")]),
       ("GIGASTROKE EUR and WMH", [("GIGA_EUR_IS", "IS"), ("GIGA_EUR_AS", "AS"), ("GIGA_EUR_SVS", "SVS"), ("WMH", "WMH")]),
       ("MEGASTROKE EUR", [("MEGA_EUR_AS", "AS"), ("MEGA_EUR_IS", "IS"), ("MEGA_EUR_SVS", "SVS")]),
       ("FinnGen stroke", [("FinnGen_I9_STR", "I9_STR"), ("FinnGen_C_STROKE", "C_STROKE"), ("FinnGen_I9_TIA", "I9_TIA")]),
       ("UK Biobank (in-sample LD)", [("IS_I63", "I63"), ("CAD_phe411.4", "CAD")]),
       ("GIGASTROKE EAS", [("GIGA_EAS_IS", "IS")])]
axc = fig.add_axes([0.36, 0.085, 0.58, 0.885])
lp = lambda p: -np.log10(max(float(p), 1e-320))
XMAX = 12
yt, yl = [], []
y = 0.0
rows_c = []
for gname, items in GRP:
    y -= 1.0; rows_c.append(("H", gname, y))
    note = {"FinnGen AMD (FinnGen LD)": "LD-sensitive (Z ≈ 47)", "UK Biobank AMD (in-sample LD)": "stable (Z ≈ 12)"}.get(gname)
    if note:
        y -= 0.85; rows_c.append(("N", note, y))
    for t, tl in items:
        y -= 1.0; rows_c.append(("T", (t, tl), y))
for kind, obj, yy in rows_c:
    if kind == "H":
        ht = axc.text(-0.45, yy, obj, transform=axc.get_yaxis_transform(), fontsize=7, fontweight="bold", va="center", bbox=dict(fc="white", ec="none", pad=0.4), zorder=5)
        continue
    if kind == "N":
        axc.text(-0.45, yy, obj, transform=axc.get_yaxis_transform(), fontsize=7, color=INK2, va="center", bbox=dict(fc="white", ec="none", pad=0.4), zorder=5)
        continue
    t, tl = obj; m = sc.loc[t]; sub = rs[rs.trait == t]
    for col_p, colr, dy in [("p_side_given_lead_2df", AQUA, 0.17), ("p_lead_given_side_1df", VIOL, -0.17)]:
        v = lp(m[col_p]); vals = [lp(p) for p in sub[col_p]] if len(sub) else [v]
        lo_, hi_ = min(vals), max(vals)
        if v > XMAX:
            axc.plot(XMAX, yy + dy, ">", ms=4, color=colr, clip_on=False, zorder=4)
            pv = float(m[col_p]); e = int(np.floor(np.log10(pv))); mant = pv / 10 ** e
            lab = f"10$^{{{e}}}$" if mant < 1.5 else f"{mant:.0f} × 10$^{{{e}}}$"
            axc.text(XMAX - 0.35, yy + dy, lab, fontsize=7, ha="right", va="center", color=colr)
            continue
        axc.hlines(yy + dy, lo_, min(hi_, XMAX), color=colr, lw=0.9, zorder=2)
        axc.plot(v, yy + dy, "o", ms=3.4, color=colr, zorder=3)
    yt.append(yy); yl.append(tl)
axc.axvline(-np.log10(0.05), color=GREY, ls="--", lw=0.6, zorder=0)
axc.set_xlim(0, XMAX); axc.set_xticks([0, 2, 4, 6, 8, 10, 12])
axc.set_yticks(yt); axc.set_yticklabels(yl, fontsize=7); axc.set_ylim(y - 0.7, 2.1)
axc.tick_params(axis="y", length=0)
axc.set_xlabel("−log$_{10}$P (matched LD)")
axc.legend(handles=[plt.Line2D([], [], marker="o", color=AQUA, ms=3.4, lw=0.9, label="side | lead (2 df)"),
                    plt.Line2D([], [], marker="o", color=VIOL, ms=3.4, lw=0.9, label="lead | side (1 df)")],
           loc="upper right", fontsize=7, borderaxespad=0.1, handletextpad=0.3, handlelength=1.4, ncol=1, labelspacing=0.2, frameon=False)

save_supp_png(fig, "FigS11")
