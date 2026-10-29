"""Figure S12 (HGG Advances supplement): (A) NHANES association of photograph-graded AMD with prevalent stroke;
(B) instrument relevance and positive controls. Panels a and c of Fig. 2 in the Human Genomics version."""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
from hg_common import *
src = open(f"{PROJ_ROOT}/scripts/manuscript_build/build/fig2.py", encoding="utf-8").read()
a_code = src.split("# ---------------- a: NHANES (unchanged)")[1].split("# ---------------- b:")[0]
c_code = src.split("# ---------------- c: positive controls (same MR pipeline)")[1].split('save(fig, "Fig2")')[0]
fig = plt.figure(figsize=(W, 150 * MM))
# shift panel a (designed at y 0.735-0.97 of a 215-mm figure) and panel c (y 0.05-0.17) into a 150-mm figure
a_code = a_code.replace("[0.36, 0.735, 0.40, 0.235]", "[0.36, 0.585, 0.40, 0.34]").replace('panel(fig, (0.005, 0.978), "a")', 'panel(fig, (0.005, 0.965), "A")')
c_code = (c_code.replace("[0.31, 0.05, 0.18, 0.12]", "[0.31, 0.08, 0.18, 0.30]").replace("[0.76, 0.05, 0.18, 0.12]", "[0.76, 0.08, 0.18, 0.30]")
          .replace('panel(fig, (0.005, 0.205), "c")', 'panel(fig, (0.005, 0.44), "B")'))
exec(a_code); exec(c_code)
save(fig, "FigureS12")
