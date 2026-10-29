"""Shared style for Human Genomics (BMC) figures: 170 mm width, Arial, >=7 pt text, bold lowercase panel labels."""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

P = Path(PROJ_ROOT)
R = P / "results"
FD0 = P / "submission" / "build" / "figdata"          # original export (regional GWAS, LD, genes)
FD = P / "submission" / "revision" / "build" / "figdata"
RR = P / "results" / "revision"
OUT = P / "submission" / "revision" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
MM = 1 / 25.4
W = 170 * MM

# color roles (dataviz reference palette; first three slots validated all-pairs, light mode)
CON = "#2a78d6"      # concordant (same direction for AMD and stroke/CSVD)
ANT = "#eb6834"      # opposing direction
AQUA = "#1baf7a"     # rs2672592-side LD group
VIOL = "#4a3aa7"     # rs3750848 / AMD-haplotype LD group
GREY = "#8a8985"
LGREY = "#c9c8c4"
INK = "#0b0b0b"
INK2 = "#52514e"

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "legend.title_fontsize": 7,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "lines.linewidth": 0.8,
    "axes.edgecolor": INK2, "axes.labelcolor": "#000000", "xtick.color": INK2, "ytick.color": INK2,
    "xtick.labelcolor": "#000000", "ytick.labelcolor": "#000000",   # tick labels black (tick marks stay grey)
    "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "savefig.facecolor": "white", "figure.facecolor": "white", "legend.frameon": False,
})


def tsv(path, **kw):
    return pd.read_csv(path, sep="\t", **kw)


def panel(fig, ax_or_xy, letter, dx=-0.0, dy=0.0):
    """Bold lowercase panel label at the top-left of an axes (or at figure coords)."""
    if isinstance(ax_or_xy, tuple):
        x, y = ax_or_xy
    else:
        bb = ax_or_xy.get_position()
        x, y = bb.x0 + dx, bb.y1 + dy
    fig.text(x, y, letter, fontsize=10, fontweight="bold", ha="left", va="bottom")


def grey_text_to_black(fig):
    """Recolour grey/achromatic text (annotations, labels, legends, 'n.s.', 'n/a', ...) to pure black (#000000).
    Coloured text (blue/orange/green/purple emphasis) and white text inside heatmap cells are left unchanged."""
    from matplotlib.text import Text
    from matplotlib.colors import to_rgba
    fig.canvas.draw()                                   # make sure tick labels exist
    for t in fig.findobj(Text):
        r, g, b, _ = to_rgba(t.get_color())
        if max(r, g, b) - min(r, g, b) < 0.04 and max(r, g, b) < 0.95:
            t.set_color("#000000")


def save(fig, name):
    grey_text_to_black(fig)
    for ext in ("pdf",):
        fig.savefig(OUT / f"{name}.{ext}")
    fig.savefig(OUT / f"{name}.png", dpi=200)
    fig.savefig(OUT / f"{name}_tmp.png", dpi=600)
    im = Image.open(OUT / f"{name}_tmp.png").convert("RGB")
    im.save(OUT / f"{name}.tif", compression="tiff_lzw", dpi=(600, 600))
    (OUT / f"{name}_tmp.png").unlink()
    plt.close(fig)
    print("saved", name, f"{fig.get_size_inches()[0] / MM:.0f} x {fig.get_size_inches()[1] / MM:.0f} mm", flush=True)


def forest(ax, y, est, lo, hi, color=INK, marker="o", ms=3.5, log=True, filled=True, lw=0.8, z=3):
    ax.hlines(y, lo, hi, color=color, lw=lw, zorder=z)
    ax.plot(est, y, marker, ms=ms, mfc=color if filled else "white", mec=color, mew=0.8, ls="none", zorder=z + 1)
    if log:
        ax.set_xscale("log")


def save_supp_png(fig, name, dpi=300):
    """Supplementary figure: PNG (300 dpi) in figures/supp, grey text recoloured to black."""
    grey_text_to_black(fig)
    (OUT / "supp").mkdir(exist_ok=True)
    fig.savefig(OUT / "supp" / f"{name}.png", dpi=dpi)
    plt.close(fig)
    print("saved", name, f"{fig.get_size_inches()[0] / MM:.0f} x {fig.get_size_inches()[1] / MM:.0f} mm", flush=True)
