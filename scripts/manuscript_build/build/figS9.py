"""Figure S9 (former Fig. 6a): key shared variants x stroke subtypes and CSVD markers."""
from rev_common import *
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch

fig = plt.figure(figsize=(W, 128 * MM))

# ---------------- a: key variants x traits (z for the AMD risk allele), evidence tier from Table 2
kv = tsv(FD0 / "fig6_key_variants_raw.tsv")
kv = pd.concat([kv[kv.SNP != "rs920065566"], tsv(FD0 / "fig6_abo_rs947073006.tsv")], ignore_index=True)
ex = tsv(RR / "extra/key_variants_ich_earlyamd_eas.tsv")
ep = tsv(P / "results" / "revision2" / "epvs_corrected.tsv")  # primary (sample-size) meta-analysis Z
from statistics import NormalDist
amd = kv[kv.trait == "AMD_FG"].set_index("SNP")
ROWS = [("rs3750848", "10q26 · rs3750848 (AMD haplotype)", "Primary"),
        ("rs2672592", "10q26 · rs2672592 (vascular signal)", "Primary (heterogeneous for stroke)"),
        ("rs4611812", "3p14 · rs4611812 (ADAMTS9)", "Primary (AMD–T2D); exploratory (stroke)"),
        ("rs429358", "19q13 · rs429358 (APOE)", "Primary (positive control)"),
        ("rs1269852", "6p21 · rs1269852 (C4A/C4B)", "Exploratory"),
        ("rs947073006", "9q34 · rs947073006 (ABO)", "Exploratory"),
        ("rs3184504", "12q24 · rs3184504 (SH2B3)", "Exploratory")]
COLS = [("IS", "IS"), ("EAS_IS", "IS\nEAS"), ("LAS", "LAS"), ("CES", "CES"), ("SVS", "SVS"), ("LACUNAR", "Lacu-\nnar"), ("ICH_UKBWGS", "ICH"),
        ("WMH", "WMH"), ("WM_MD", "WM\nMD"), ("WM_FA", "WM\nFA†"), ("EPVS_WM", "EPVS\nWM"), ("EPVS_BG", "EPVS\nBG"), ("EPVS_HIP", "EPVS\nHIP")]
Z = np.full((len(ROWS), len(COLS)), np.nan); PV = Z.copy(); AZ = []
for i, (s, _, _) in enumerate(ROWS):
    a = amd.loc[s]; risk = a.ea if a.beta > 0 else a.oa; AZ.append(abs(a.beta / a.se))
    for j, (t, _) in enumerate(COLS):
        if t in ("EAS_IS", "ICH_UKBWGS"):
            q = ex[(ex.SNP == s) & (ex.trait == t)]
            if len(q): Z[i, j], PV[i, j] = q.z_risk.iloc[0], q.p.iloc[0]
            continue
        if t.startswith("EPVS"):
            q = ep[(ep.rsid == s) & (ep.compartment == t)]
            if len(q):
                q = q.iloc[0]; sign = 1 if q.effect_allele == risk else -1
                Z[i, j], PV[i, j] = sign * q.z_primary_harmonised, q.p_primary
                continue
            q = kv[(kv.SNP == s) & (kv.trait == t)]
            if not len(q): continue
            q = q.iloc[0]; sign = 1 if q.ea == risk else (-1 if q.oa == risk else np.nan)
            Z[i, j] = sign * np.sign(q.beta) * NormalDist().inv_cdf(1 - q.p / 2); PV[i, j] = q.p
            continue
        q = kv[(kv.SNP == s) & (kv.trait == t)]
        if not len(q): continue
        q = q.iloc[0]
        sign = 1 if q.ea == risk else (-1 if q.oa == risk else np.nan)
        Z[i, j] = sign * q.beta / q.se * (-1 if t == "WM_FA" else 1); PV[i, j] = q.p
axa = fig.add_axes([0.33, 0.21, 0.62, 0.598])
cmap = LinearSegmentedColormap.from_list("div", [ANT, "#f0efec", CON])
im = axa.imshow(np.ma.masked_invalid(np.clip(Z, -6, 6)), cmap=cmap, norm=TwoSlopeNorm(0, -6, 6), aspect="auto")
for i in range(Z.shape[0]):
    for j in range(Z.shape[1]):
        if np.isnan(Z[i, j]):
            axa.text(j, i, "n/a", ha="center", va="center", fontsize=7, color=GREY)
        else:
            st = "*" if PV[i, j] < 0.05 else ""
            txt = f"{Z[i, j]:.1f}" if round(Z[i, j], 1) != 0 else "0.0"
            axa.text(j, i, (txt + st).replace("-", "−"), ha="center", va="center", fontsize=7,
                     fontweight="bold" if PV[i, j] < 1e-4 else "normal", color="white" if abs(Z[i, j]) > 4.5 else INK)
axa.set_xticks(range(len(COLS))); axa.set_xticklabels([c[1] for c in COLS], fontsize=7); axa.xaxis.tick_top()
axa.set_yticks(range(len(ROWS))); axa.set_yticklabels([f"{lab}\nAMD z = {z:.0f}\n{tier}" for (_, lab, tier), z in zip(ROWS, AZ)], fontsize=7, linespacing=1.1)
axa.tick_params(length=0); [sp.set_visible(False) for sp in axa.spines.values()]
for xb in [6.5]:
    axa.axvline(xb, color="white", lw=2.5)
axa.text(3.0, -1.6, "Stroke", ha="center", fontsize=7, fontweight="bold", color=INK2)
axa.text(9.5, -1.6, "CSVD markers", ha="center", fontsize=7, fontweight="bold", color=INK2)
cb = fig.add_axes([0.33, 0.135, 0.30, 0.016]); cbar = fig.colorbar(im, cax=cb, orientation="horizontal")
cbar.set_label("z for the AMD risk allele (> 0 concordant, < 0 opposing); * P < 0.05, bold P < 10$^{-4}$", fontsize=7)
cbar.ax.tick_params(labelsize=7, length=2); cbar.outline.set_linewidth(0.4)

save_supp_png(fig, "FigS9")
