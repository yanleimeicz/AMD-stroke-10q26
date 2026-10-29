"""HGG Advances figure style: reuses the Human Genomics revision style (rev_common) with uppercase panel labels
and a separate output folder."""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import sys
sys.path.insert(0, f"{PROJ_ROOT}/scripts/manuscript_build/build")
import rev_common as _rc
from rev_common import *          # noqa: F401,F403
OUT = Path(f"{PROJ_ROOT}/scripts/hggadv/figures/out")
OUT.mkdir(parents=True, exist_ok=True)
_rc.OUT = OUT


def panel(fig, ax_or_xy, letter, dx=-0.0, dy=0.0):
    """Bold uppercase panel label (HGG Advances: uppercase letters, upper-left corner)."""
    return _rc.panel(fig, ax_or_xy, letter.upper(), dx, dy)
