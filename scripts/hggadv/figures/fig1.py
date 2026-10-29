"""Figure 1 (HGG Advances): study design, two questions."""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
from hg_common import *
from matplotlib.patches import FancyBboxPatch

H = 112
fig = plt.figure(figsize=(W, H * MM))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 170); ax.set_ylim(0, H); ax.axis("off")

def box(x, y, w, h, fc, ec, lw=0.8):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.5", fc=fc, ec=ec, lw=lw))

def arrow(x, y0, y1):
    ax.annotate("", xy=(x, y1), xytext=(x, y0), arrowprops=dict(arrowstyle="-|>", lw=0.8, color=INK2, mutation_scale=7))

box(1.5, 96, 167, 14, "#f3f2ef", INK, lw=1.0)
ax.text(85, 103, "AMD and stroke co-occur in observational studies, yet published MR suggests that\n"
        "genetic liability to AMD protects against small-vessel stroke",
        ha="center", va="center", fontsize=8, fontweight="bold", linespacing=1.3)
arrow(85, 95.6, 91.6)
ax.text(21.5, 88, "Question", ha="center", fontsize=7.5, fontweight="bold")
ax.text(87, 88, "Data and analyses", ha="center", fontsize=7.5, fontweight="bold")
ax.text(150.5, 88, "Answer", ha="center", fontsize=7.5, fontweight="bold")
rows = [
    ("Question 1", "Does the protective\nMR estimate reflect\nAMD liability or\na single locus?",
     "MR with one instrument per AMD locus: FinnGen and IAMDGC AMD →\nGIGASTROKE, MEGASTROKE (no FinnGen), WMH volume\n"
     "IVW with MR-Egger, weighted median, MR-RAPS, MR-PRESSO, CAUSE\n"
     "Per-locus Wald ratios, leave-one-locus-out MR, MR-Clust\nGenome-wide (LDSC) and local (LAVA) genetic correlation",
     "The protective estimate\narises from 10q26;\nwithout 10q26, estimates\nare close to 1"),
    ("Question 2", "At 10q26, do AMD\nand vascular disease\nmap to the same\nhaplotype?",
     "In-sample LD: FinnGen fine-mapping; UK Biobank conditional\nanalysis and fine-mapping\n"
     "Coronary disease: UK Biobank, FinnGen, CARDIoGRAMplusC4D 2015,\nAragam et al., BioBank Japan; ischemic stroke: MEGASTROKE,\n"
     "GIGASTROKE remainder, FinnGen, UK Biobank, GIGASTROKE East Asian\nExpression: RPE single-nucleus, retina, blood and artery eQTLs",
     "AMD: ARMS2 p.Ala69Ser\nhaplotype\nCoronary disease: rs2672592\nhaplotype (conditional;\nnot in BioBank Japan)\nIschemic stroke:\nunresolved"),
]
HS = [36, 40]
y = 85
for i, ((q, qt, data, ans), h) in enumerate(zip(rows, HS)):
    box(1.5, y - h, 40, h, "#eef3fb", CON)
    ax.text(21.5, y - 3.2, q, ha="center", va="top", fontsize=7.5, fontweight="bold")
    ax.text(21.5, y - h / 2 - 2, qt, ha="center", va="center", fontsize=7, linespacing=1.2)
    box(44, y - h, 86, h, "white", GREY)
    ax.text(46, y - h / 2, data, va="center", fontsize=7, linespacing=1.3)
    box(132.5, y - h, 36, h, "#fdf0ea", ANT)
    t = ax.text(150.5, y - h / 2, ans, ha="center", va="center", fontsize=7, linespacing=1.25)
    if i < len(rows) - 1:
        arrow(21.5, y - h - 0.3, y - h - 3.6)
    y -= h + 4
save(fig, "Figure1")
