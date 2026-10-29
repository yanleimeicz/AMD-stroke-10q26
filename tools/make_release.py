"""Copy the manuscript-relevant scripts into code_release/scripts/ and replace hard-coded project paths.

Run from anywhere:  python code_release/tools/make_release.py
Only files under code_release/ are written; the original project folders are read-only inputs.
Every replacement is logged to code_release/PATH_CHANGES.md (and tools/path_changes.tsv).
"""
import re, shutil, sys
from pathlib import Path

PROJ = Path(__file__).resolve().parents[2]          # original project (G:\Comorbidity)
REL = PROJ / "code_release"
DST = REL / "scripts"

TOP = """01_download_nhanes.py 02_nhanes_codebook.py 10_nhanes_stroke_amd.R 11_nhanes_forest_plots.R
20_gwascatalog_lookup.py 21_check_sources.py 30_prepare_gwas.R 30b_check_iamdgc_sign.R 31_mr_bidirectional.R
32_coloc_gwas.R 33_eqtl_smr.R 34_ldsc_rg.R 35_coloc_post.R 36_mr_steiger.R 37_mr_coloc_figures.R
38_gunzip_besd.py 38a_make_ma.R 38b_ld_subset.R 38_smr_heidi_druggable.R 38c_dgidb_categories.py 39_pqtl_mr.R
45_csvd_gwas_prep.R 46_csvd_mr_ldsc.R 47_lava_local_rg.R 48a_fetch_tissue_eqtl.R 48b_gwas_regions.R 48_tissue_eqtl.R
49_cojo_susie.R 50_cxg_search.py 50_scrna.py 52_scrna_adamts9.py 74a_adamts9_data.R 74b_adamts9_cimt.R 74_adamts9.R
75_c4_conditional.R 75a_rev_adamts9_fetch.R 75b_rev_adamts9_finngen_susie.R 75c_rev_adamts9_mvmr.R
75d_rev_adamts9_locus.R 76_rev_extra.R rev_mr_00_prep.R rev_mr_01_mr.R rev_mr_01b_raps.R
rev_mr_02_fdr_concordance.R rev_mr_03_robust.R""".split()
SUB = """revision/R1_10q26_insample.R
revision2/rev2_10q26.py revision2/rev2_concordant_blocks.R revision2/rev2_equivalence.py revision2/rev2_remainder.py
revision3/panukb_01_fetch.py revision3/panukb_02_analysis.py revision3/panukb_03_susie.R revision3/panukb_04_followup.py
revision3/panukb_05_fourvar.py revision3/r3_fetch.py revision3/r3_pdftext.py revision3/r3_eye_build.py
revision3/urls1.tsv revision3/urls2.tsv revision3/urls3.tsv
revision4/rev4_AB.py revision4/rev4_C_coloc.R revision4/rev4_C_extra.py revision4/rev4_C_fetch.py
revision5/r5_00_fetch.py revision5/r5_01_kappa_complement_af.R revision5/r5_02_contrast.py revision5/r5_04_corivw.R
revision5/r5_05_eqtl.py story/story_01_mr.R story/story_02_tables.R""".split()
if (PROJ / "results" / "revision6" / "SUMMARY.md").exists():     # round-6 scripts only once round 6 is finalised (SUMMARY.md written)
    SUB += [f"revision6/{f.name}" for f in sorted((PROJ / "scripts" / "revision6").iterdir()) if f.suffix in (".py", ".R")]
SUB += [f"revision7/{f.name}" for f in sorted((PROJ / "scripts" / "revision7").iterdir()) if f.suffix in (".py", ".R")]

COPIES = [(PROJ / "scripts" / f, DST / f) for f in TOP + SUB]
COPIES += [
    (PROJ / "tmp/revmr/cause_fix.R", DST / "helpers/cause_fix.R"),
    (PROJ / "tmp/rev_adamts9/phewas_gc.py", DST / "helpers/phewas_gc.py"),
    (PROJ / "tmp/pqtl/xlsx2tsv.py", DST / "helpers/pqtl_xlsx2tsv.py"),
    (PROJ / "submission/revision/build_docx.py", DST / "manuscript_build/build_docx.py"),
    (PROJ / "submission/revision/build_af1.py", DST / "manuscript_build/build_af1.py"),
    (PROJ / "submission/build/fig_export.R", DST / "manuscript_build/build_v1/fig_export.R"),
    (PROJ / "submission/build/fig_export_abo.R", DST / "manuscript_build/build_v1/fig_export_abo.R"),
    (PROJ / "results/integration/scrna_genes.txt", REL / "config/scrna_genes.txt"),
    (PROJ / "tmp/rev_adamts9/finngen_counts.tsv", REL / "config/finngen_counts.tsv"),
    (PROJ / "derived/eye_eqtl/ld_snps.txt", REL / "config/eye_eqtl/ld_snps.txt"),
    (PROJ / "derived/eye_eqtl/ld2_snps.txt", REL / "config/eye_eqtl/ld2_snps.txt"),
]
for f in """fig1.py fig2.py fig3.py fig4.py fig5.py figS8.py figS9.py figS10.py figS11.py figs_supp.py rev_common.py
build_supp_tables_rev.py build_strobe_mr_rev.py check_text_colors.py rev_fig_export.R rev_s4_export.R""".split():
    COPIES.append((PROJ / "submission/revision/build" / f, DST / "manuscript_build/build" / f))

# special relocations (file moved into the release) -> replacement target
RELOC = {"tmp/revmr/cause_fix.R": "scripts/helpers/cause_fix.R"}

R_HDR = """# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
"""
PY_HDR = """import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\\\", "/")
"""

LOG = []
ROOT_RX = r"(?:G:[\\/]+Comorbidity|/workspace/Comorbidity)"


def r_sub(line, fname, ln):
    def rep(m):
        rest = m.group(2).replace("\\", "/").lstrip("/")
        worker = "clusterEvalQ" in line
        base = 'Sys.getenv("COMORBIDITY_ROOT")' if worker else "PROJ_ROOT"
        if rest in RELOC: rest = RELOC[rest]
        new = base if not rest else f'file.path({base}, "{rest}")'
        LOG.append((fname, ln, m.group(0), new))
        return new
    return re.sub(r'(["\'])' + ROOT_RX + r'([^"\']*)\1', rep, line)


def py_sub(line, fname, ln):
    def rep(m):
        pre, q, rest = m.group(1) or "", m.group(2), m.group(3)
        rest = rest.replace("\\\\", "/").replace("\\", "/").lstrip("/")
        if not rest:
            new = "PROJ_ROOT"
        else:
            if "{" in rest and "f" not in pre.lower():
                rest = rest.replace("{", "{{").replace("}", "}}")
            new = f'f{q}{{PROJ_ROOT}}/{rest}{q}'
        LOG.append((fname, ln, m.group(0), new))
        return new
    return re.sub(r'(?<![A-Za-z0-9_])([rRfF]{0,2})(["\'])' + ROOT_RX + r'([^"\']*)\2', rep, line)


def insert_py_header(src):
    lines = src.splitlines(keepends=True)
    i = 0
    # skip shebang/encoding comments and module docstring and __future__ imports
    while i < len(lines) and lines[i].startswith("#"): i += 1
    if i < len(lines) and lines[i].lstrip().startswith(('"""', "'''", 'r"""')):
        q = '"""' if '"""' in lines[i] else "'''"
        if lines[i].count(q) >= 2: i += 1
        else:
            i += 1
            while i < len(lines) and q not in lines[i]: i += 1
            i += 1
    while i < len(lines) and lines[i].startswith("from __future__"): i += 1
    return "".join(lines[:i]) + PY_HDR + "".join(lines[i:])


def insert_r_header(src):
    lines = src.splitlines(keepends=True)
    i = 0
    while i < len(lines) and lines[i].lstrip().startswith("#"): i += 1
    return "".join(lines[:i]) + R_HDR + "".join(lines[i:])


def main():
    if DST.exists(): shutil.rmtree(DST)
    for s, d in COPIES:
        if not s.exists(): sys.exit(f"missing {s}")
        d.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(s, d)
        rel = d.relative_to(REL).as_posix()
        if d.suffix not in (".R", ".py"): continue
        src = d.read_text(encoding="utf-8")
        out, changed = [], False
        for ln, line in enumerate(src.splitlines(keepends=True), 1):
            new = (r_sub if d.suffix == ".R" else py_sub)(line, rel, ln)
            changed |= new != line; out.append(new)
        if changed:
            txt = "".join(out)
            txt = insert_r_header(txt) if d.suffix == ".R" else insert_py_header(txt)
            d.write_text(txt, encoding="utf-8", newline="")
    # leftovers check
    left = [(p, i, l.strip()) for p in DST.rglob("*") if p.suffix in (".R", ".py")
            for i, l in enumerate(p.read_text(encoding="utf-8").splitlines(), 1) if re.search(ROOT_RX, l)]
    with open(REL / "tools/path_changes.tsv", "w", encoding="utf-8") as f:
        f.write("file\tline_in_original\told\tnew\n")
        for r in LOG: f.write("\t".join(map(str, r)) + "\n")
    print(f"copied {len(COPIES)} files; {len(LOG)} replacements in {len({r[0] for r in LOG})} files; leftovers: {len(left)}")
    for x in left: print("LEFT", x)


if __name__ == "__main__":
    main()
