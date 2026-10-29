from rev_common import *
from matplotlib.patches import FancyBboxPatch

H = 128
fig = plt.figure(figsize=(W, H * MM))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 170); ax.set_ylim(0, H); ax.axis("off")

def box(x, y, w, h, fc, ec, lw=0.8):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.5", fc=fc, ec=ec, lw=lw))

def arrow(x, y0, y1):
    ax.annotate("", xy=(x, y1), xytext=(x, y0), arrowprops=dict(arrowstyle="-|>", lw=0.8, color=INK2, mutation_scale=7))

# clinical observation
box(1.5, 112, 167, 12, "#f3f2ef", INK, lw=1.0)
ax.text(85, 118, "Clinical observation: AMD patients have more strokes (NHANES 2005–2008, OR 2.08; 95% CI 1.20–3.61)",
        ha="center", va="center", fontsize=8, fontweight="bold")
arrow(85, 111.6, 107.6)

ax.text(21.5, 104, "Question", ha="center", fontsize=7.5, fontweight="bold")
ax.text(87, 104, "Main data and analyses", ha="center", fontsize=7.5, fontweight="bold")
ax.text(150.5, 104, "Answer", ha="center", fontsize=7.5, fontweight="bold")

rows = [
    ("Question 1", "Are MR estimates compatible\nwith the effect implied by\nthe NHANES association?",
     "MR, one instrument per AMD locus: FinnGen and IAMDGC AMD →\nGIGASTROKE, MEGASTROKE, CAD (Aragam et al.), WMH\n"
     "NHANES benchmark scaled to early-AMD liability\n"
     "Positive controls: independent AMD GWAS; systolic blood\npressure; atrial fibrillation",
     "MR does not explain\nthe observational\nassociation; a meaningful\ncausal effect is\nnot excluded"),
    ("Question 2", "Which loci produce the\nconflicting MR estimates?",
     "Per-locus Wald ratios and leave-one-locus-out MR\nLAVA local genetic correlation (448 bivariate tests)",
     "10q26 alone produces\nthe apparent protection"),
    ("Question 3", "Do AMD and vascular\ndisease map to the\nsame haplotype?",
     "In-sample fine-mapping and conditional analysis: FinnGen, UK Biobank\n"
     "Stroke: GIGASTROKE EUR/EAS, MEGASTROKE, FinnGen, UK Biobank;\nCAD: UK Biobank, FinnGen, CARDIoGRAMplusC4D 2015, Aragam et al.\n"
     "Tissue eQTL incl. RPE single-nucleus eQTL; MR partitioned\nby locus group (complement, 10q26, lipid, other)",
     "AMD: ARMS2 p.Ala69Ser\nhaplotype; CAD: preferentially\nrs2672592 haplotype;\nstroke: unresolved;\ncomplement: no robust\nassociation (limited power)"),
]
HS = [32, 22, 34]
y = 101
for i, ((q, qt, data, ans), h) in enumerate(zip(rows, HS)):
    box(1.5, y - h, 40, h, "#eef3fb", CON)
    ax.text(21.5, y - 3.2, q, ha="center", va="top", fontsize=7.5, fontweight="bold")
    ax.text(21.5, y - h / 2 - 2, qt, ha="center", va="center", fontsize=7, linespacing=1.2)
    box(44, y - h, 86, h, "white", GREY)
    ax.text(46, y - h / 2, data, va="center", fontsize=7, linespacing=1.3)
    box(132.5, y - h, 36, h, "#fdf0ea", ANT)
    ax.text(150.5, y - h / 2, ans, ha="center", va="center", fontsize=7, linespacing=1.25)
    if i < len(rows) - 1:
        arrow(21.5, y - h - 0.3, y - h - 3.6)
    y -= h + 4
save(fig, "Fig1")
