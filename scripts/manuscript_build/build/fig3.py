from rev_common import *

H = 200
fig = plt.figure(figsize=(W, H * MM))
F = lambda mm: mm / H
ST = P / "results" / "story"

# ---------------- a: per-locus effects on AMD (x) and IS (y), FinnGen AMD instruments; IVW slopes
pl = tsv(ST / "per_locus_effects.tsv"); pl = pl[(pl.exposure == "AMD_FG") & (pl.outcome == "IS")]
pw = tsv(ST / "pathway_mr.tsv"); pw = pw[(pw.exposure == "AMD_FG") & (pw.outcome == "IS")].set_index("group")
PCOL = {"complement": VIOL, "10q26": ANT, "lipid": AQUA, "other": GREY}
PLAB = {"complement": "Complement", "10q26": "10q26", "lipid": "Lipid", "other": "Other"}
axa = fig.add_axes([0.10, F(122), 0.40, F(68)])
for pwy in ["other", "lipid", "complement", "10q26"]:
    d = pl[pl.pathway == pwy]
    axa.errorbar(d.beta_amd, d.beta_out, xerr=1.96 * d.se_amd, yerr=1.96 * d.se_out, fmt="o", ms=3.6, color=PCOL[pwy],
                 ecolor=PCOL[pwy], elinewidth=0.6, capsize=0, zorder=3 if pwy != "other" else 2)
xx = np.array([0, 0.75])
axa.plot(xx, pw.loc["all", "b"] * xx, color=INK, lw=0.9, zorder=1, label=f"IVW slope, all loci (OR {pw.loc['all', 'OR']:.3f})")[0]
axa.plot(xx, pw.loc["no_10q26", "b"] * xx, color=INK, lw=0.9, ls=(0, (4, 2)), zorder=1,
         label=f"IVW slope, without 10q26 (OR {pw.loc['no_10q26', 'OR']:.3f})")
axa.axhline(0, color=GREY, lw=0.5, ls=":", zorder=0)
LBL = {"rs61871744": ("10q26", (6, -10)), "rs1329424": ("CFH", (6, 4)), "rs2230199": ("C3", (5, 6)), "rs429358": ("APOE", (3, -14))}
for snp, (lab, off) in LBL.items():
    r = pl[pl.SNP == snp].iloc[0]
    axa.annotate(lab, (r.beta_amd, r.beta_out), xytext=off, textcoords="offset points", fontsize=7, fontweight="bold")
axa.set_xlim(0, 0.75); axa.set_ylim(-0.10, 0.075)
axa.set_xlabel("Effect on AMD, log OR (FinnGen)"); axa.set_ylabel("Effect on IS, log OR (GIGASTROKE)")
lg = axa.legend(handles=[plt.Line2D([], [], marker="o", ls="none", color=PCOL[k], ms=3.6, label=PLAB[k]) for k in ["complement", "10q26", "lipid", "other"]],
                loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=4, fontsize=7, handletextpad=0.2, columnspacing=1.0, borderaxespad=0.1)
axa.add_artist(lg)
axa.legend(handles=[l for l in axa.get_lines() if l.get_label().startswith("IVW")], loc="lower left", fontsize=7, handlelength=2.2,
           handletextpad=0.4, borderaxespad=0.3, labelspacing=0.3)
panel(fig, (0.005, F(195)), "a")

# ---------------- b: leave-one-locus-out MR for SVS (FinnGen AMD)
mr = tsv(RR / "mr/mr_revised_all.tsv")
mr = mr[(mr.scheme == "A_dist") & (mr.method == "IVW") & (mr.exposure == "AMD_FG") & (mr.outcome == "SVS")].set_index("exclusion")
EXC = [("main", "All loci"), ("no_10q26", "Without 10q26"), ("no_ADAMTS9", "Without ADAMTS9"), ("no_CFH", "Without CFH"),
       ("no_C3", "Without C3"), ("no_APOE", "Without APOE"), ("no_major3", "Without 10q26, CFH and C3"), ("no_all", "Without all major loci")]
axb = fig.add_axes([0.76, F(128), 0.14, F(56)])
for k, (ex, lab) in enumerate(EXC):
    r = mr.loc[ex]; y = len(EXC) - 1 - k
    est, lo, hi = np.exp(r.b), np.exp(r.b - 1.96 * r.se), np.exp(r.b + 1.96 * r.se)
    col = ANT if ex == "main" else (AQUA if ex == "no_10q26" else INK)
    forest(axb, y, est, lo, hi, color=col, filled=r.p < 0.05, ms=3.6, log=False)
    axb.text(1.04, y, f"{est:.2f}", transform=axb.get_yaxis_transform(), va="center", fontsize=7)
axb.axvline(1, color=GREY, ls="--", lw=0.6, zorder=0)
axb.set_yticks(range(len(EXC))); axb.set_yticklabels([e[1] for e in EXC][::-1], fontsize=7)
axb.set_ylim(-0.6, len(EXC) - 0.4); axb.set_xlim(0.78, 1.22); axb.set_xticks([0.8, 1.0, 1.2])
axb.set_xlabel("OR for SVS per unit\nlog-odds of AMD liability")
axb.set_title("Small-vessel stroke (FinnGen AMD)", fontsize=7, fontweight="bold", pad=3)
panel(fig, (0.53, F(195)), "b")

# ---------------- c: 11 unique LAVA blocks (best pair per block) + MEGASTROKE replication at 10q26
ub = tsv(RR / "mr/lava_unique_blocks.tsv")
bv = tsv(R / "csvd/lava_bivariate_fdr05.tsv")
sens = tsv(RR / "mr/lava_10q26_sensitivity.tsv")
LOC = {1597: "10q26 ARMS2–HTRA1", 1596: "10q26 FGFR2–TACC2 (adjacent)", 2453: "21q22 near ERG", 2302: "18q23 near MBP–GALR1",
       1613: "11p15 near KCNQ1–CDKN1C", 2159: "16q23 near MAF–WWOX", 2163: "16q23 near CDH13", 1845: "12q24 near MED13L",
       1852: "12q24 near SCARB1–NCOR2", 699: "4q25 near PITX2", 1873: "13q12 near FLT1"}
TN = {"IS": "IS", "SVS": "SVS", "LACUNAR": "lacunar", "WMH": "WMH", "WM_FA": "WM FA†", "AMD_FG": "AMD",
      "DRYAMD_FG": "dry AMD", "WETAMD_FG": "wet AMD", "AMD_IAMDGC": "AMD (IAMDGC)"}
rows = []
for _, u in ub.iterrows():
    out, amd = u.best_pair.split()
    b = bv[(bv.locus == u.locus) & (bv.phen1 == out) & (bv.phen2 == amd)].iloc[0]
    s = -1 if out == "WM_FA" else 1
    lo, hi = (b["rho.lower"], b["rho.upper"]) if s > 0 else (-b["rho.upper"], -b["rho.lower"])
    rows.append(dict(locus=u.locus, lab=f"{LOC[u.locus]}  ·  {TN[amd]} – {TN[out]}", rho=s * b.rho, lo=lo, hi=hi, p=b.p, n=int(u.n_fdr_within_trait)))
d = pd.DataFrame(rows).sort_values("rho").reset_index(drop=True)
axa = fig.add_axes([0.44, F(22), 0.37, F(77)])
for i, r in d.iterrows():
    y = len(d) - 1 - i
    forest(axa, y, r.rho, r.lo, r.hi, color=ANT if r.rho < 0 else CON, log=False, ms=3.8)
    axa.text(1.04, y, f"{r.p:.1e}".replace("e-0", "e−").replace("e-", "e−") + f"  ({r.n})", transform=axa.get_yaxis_transform(), va="center", fontsize=7)
    if r.locus in (1597, 1596):
        m = sens[(sens.phen1 == "AMD_FG") & (sens.phen2 == "MEGA_SVS" if r.locus == 1597 else sens.phen2 == "MEGA_IS")]
        m = m[(m.locus.astype(str) == str(r.locus)) & (m.scenario == "ldsc_intercept")].iloc[0]
        yy = y - 0.36
        axa.hlines(yy, m["rho.lower"], m["rho.upper"], color=GREY, lw=0.8)
        axa.plot(m.rho, yy, "D", ms=3.4, mfc="white", mec=INK, mew=0.8)
        axa.text(1.04, yy, f"{m.p:.1e}".replace("e-0", "e−").replace("e-", "e−") + "  rep.", transform=axa.get_yaxis_transform(),
                 va="center", fontsize=7, color=INK2)
axa.axvline(0, color=GREY, lw=0.6, ls="--", zorder=0)
axa.set_yticks(range(len(d))); axa.set_yticklabels(list(d.lab)[::-1], fontsize=7)
axa.set_xlim(-1.05, 1.05); axa.set_ylim(-0.8, len(d) - 0.3)
axa.set_xlabel("Local genetic correlation, ρ (95% CI)\n(oriented so that ρ > 0 = same direction of disease risk)")
axa.text(1.04, len(d) - 0.25, "P  (n pairs)", transform=axa.get_yaxis_transform(), fontsize=7, fontweight="bold")
axa.text(-0.98, len(d) + 0.15, "Opposing (4 blocks)", color=ANT, fontsize=7, fontweight="bold")
axa.text(0.98, len(d) + 0.15, "Concordant (7 blocks)", color=CON, fontsize=7, fontweight="bold", ha="right")
axa.legend(handles=[plt.Line2D([], [], marker="D", ls="none", mfc="white", mec=INK, ms=3.6,
                               label="FinnGen AMD with MEGASTROKE (no FinnGen)")],
           loc="lower left", fontsize=7, borderaxespad=0.2, handletextpad=0.3)
panel(fig, (0.005, F(107)), "c")

save(fig, "Fig3")
