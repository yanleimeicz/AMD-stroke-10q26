# PATH_CHANGES — hard-coded paths replaced in the copied scripts

Only the copies in `code_release/scripts/` were edited; the originals in the project folder are unchanged. The edits are mechanical (string-literal substitution plus a short header); no analysis logic was changed. They were made by `tools/make_release.py`; the line-by-line log (file, line number in the original, old literal, new expression) is `tools/path_changes.tsv`.

**Total: 320 replacements in 72 files.**

## 1. Header added to every edited file

R (inserted after the leading comment block):

```r
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
```

Python (inserted after the module docstring):

```python
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
```

`PROJ_ROOT` uses forward slashes on all platforms. With `Rscript scripts/xx.R` or `python scripts/xx.py` the default resolves to the repository root (the parent of `scripts/`), wherever the command is started from. Both headers were tested (R 4.5.2 on Windows; Python 3.12 on Linux), with and without `COMORBIDITY_ROOT` set.

## 2. Substitution rules

| Language | Original literal | Replacement |
|---|---|---|
| R | `"G:/Comorbidity"` | `PROJ_ROOT` |
| R | `"G:/Comorbidity/<sub/path>"` | `file.path(PROJ_ROOT, "<sub/path>")` |
| R, inside `clusterEvalQ(...)` (evaluated on parallel workers, where `PROJ_ROOT` is not defined) | `"G:/Comorbidity/Rlib"` | `file.path(Sys.getenv("COMORBIDITY_ROOT"), "Rlib")` |
| R (helper moved into the release) | `"G:/Comorbidity/tmp/revmr/cause_fix.R"` | `file.path(PROJ_ROOT, "scripts/helpers/cause_fix.R")` |
| Python | `r"G:\Comorbidity"`, `"G:/Comorbidity"` | `PROJ_ROOT` (e.g. `Path(r"G:\Comorbidity")` → `Path(PROJ_ROOT)`) |
| Python | `r"G:\Comorbidity\<sub\path>"`, `"G:/Comorbidity/<sub/path>"`, `'G:/Comorbidity/<sub/path>'` | `f"{PROJ_ROOT}/<sub/path>"` (backslashes converted to `/`) |
| Python | `'/workspace/Comorbidity'` (Linux sandbox mount; `revision3/r3_eye_build.py`) | `PROJ_ROOT` |

## 3. Pattern counts

| File type | Original pattern (root + first sub-folder) | Context | n |
|---|---|---|---|
| .R | `"G:/Comorbidity/derived…"` | main process | 80 |
| .R | `"G:/Comorbidity/data…"` | main process | 65 |
| .R | `"G:/Comorbidity/results…"` | main process | 50 |
| .R | `"G:/Comorbidity/Rlib…"` | main process | 47 |
| .R | `"G:/Comorbidity/tmp…"` | main process | 18 |
| .R | `"G:/Comorbidity/Rlib…"` | parallel worker (clusterEvalQ) | 14 |
| .R | `"G:/Comorbidity"` | main process | 4 |
| .R | `"G:/Comorbidity/scripts…"` | main process | 2 |
| .R | `"G:/Comorbidity/Rlib_hypr…"` | main process | 1 |
| .R | `"G:/Comorbidity/tmp…"` | relocated helper | 1 |
| .R | `"G:/Comorbidity/submission…"` | main process | 1 |
| .py | `"G:/Comorbidity"` | main process | 11 |
| .py | `r"G:\Comorbidity\data…"` | main process | 6 |
| .py | `r"G:\Comorbidity\results…"` | main process | 6 |
| .py | `r"G:\Comorbidity\derived…"` | main process | 5 |
| .py | `"G:/Comorbidity/results…"` | main process | 3 |
| .py | `r"G:\Comorbidity\submission…"` | main process | 3 |
| .py | `r'G:\Comorbidity\derived…'` | main process | 1 |
| .py | `'/workspace/Comorbidity'` | main process | 1 |
| .py | `r"G:\Comorbidity"` | main process | 1 |

## 4. Replacements per file

| File (under code_release/) | n | Lines in original |
|---|---|---|
| `scripts/10_nhanes_stroke_amd.R` | 2 | 4, 8 |
| `scripts/11_nhanes_forest_plots.R` | 2 | 2, 4 |
| `scripts/30_prepare_gwas.R` | 4 | 4, 7, 56 |
| `scripts/30b_check_iamdgc_sign.R` | 3 | 2, 3 |
| `scripts/31_mr_bidirectional.R` | 5 | 5, 8, 9 |
| `scripts/32_coloc_gwas.R` | 6 | 4, 7, 8, 67 |
| `scripts/33_eqtl_smr.R` | 7 | 5, 8, 9, 10, 13 |
| `scripts/34_ldsc_rg.R` | 6 | 2, 4, 5, 16 |
| `scripts/35_coloc_post.R` | 7 | 3, 6, 7, 15, 49 |
| `scripts/36_mr_steiger.R` | 2 | 2, 4 |
| `scripts/37_mr_coloc_figures.R` | 3 | 2, 4 |
| `scripts/38_gunzip_besd.py` | 1 | 5 |
| `scripts/38_smr_heidi_druggable.R` | 12 | 7, 10, 11, 12, 13, 20, 24, 25, 59, 128 |
| `scripts/38a_make_ma.R` | 4 | 3, 6, 7 |
| `scripts/38b_ld_subset.R` | 5 | 2, 4, 5, 6 |
| `scripts/38c_dgidb_categories.py` | 3 | 11, 12, 13 |
| `scripts/39_pqtl_mr.R` | 10 | 8, 11, 13, 14, 165, 206, 218 |
| `scripts/45_csvd_gwas_prep.R` | 4 | 8, 11, 44 |
| `scripts/46_csvd_mr_ldsc.R` | 13 | 4, 7, 8, 100, 144, 145, 175, 189, 193 |
| `scripts/47_lava_local_rg.R` | 7 | 6, 9, 12, 13 |
| `scripts/48_tissue_eqtl.R` | 6 | 6, 8, 14, 31, 78 |
| `scripts/48a_fetch_tissue_eqtl.R` | 5 | 5, 7, 32, 61, 62 |
| `scripts/48b_gwas_regions.R` | 5 | 3, 5, 13, 18 |
| `scripts/49_cojo_susie.R` | 9 | 6, 8, 10, 11, 12, 17, 30, 155 |
| `scripts/50_cxg_search.py` | 1 | 5 |
| `scripts/50_scrna.py` | 2 | 17 |
| `scripts/52_scrna_adamts9.py` | 2 | 15 |
| `scripts/74_adamts9.R` | 8 | 9, 11, 13, 14, 29, 97, 132 |
| `scripts/74a_adamts9_data.R` | 9 | 6, 8, 10, 11, 16, 43, 48, 64 |
| `scripts/74b_adamts9_cimt.R` | 2 | 5, 7 |
| `scripts/75_c4_conditional.R` | 6 | 4, 6, 7, 8 |
| `scripts/75a_rev_adamts9_fetch.R` | 4 | 6, 8, 12, 57 |
| `scripts/75b_rev_adamts9_finngen_susie.R` | 3 | 3, 5, 30 |
| `scripts/75c_rev_adamts9_mvmr.R` | 6 | 7, 10, 11, 12 |
| `scripts/75d_rev_adamts9_locus.R` | 9 | 6, 9, 10, 11, 163, 231 |
| `scripts/76_rev_extra.R` | 10 | 7, 10, 61, 95, 125, 126 |
| `scripts/helpers/pqtl_xlsx2tsv.py` | 7 | 2, 3, 4, 5 |
| `scripts/manuscript_build/build/build_strobe_mr_rev.py` | 1 | 15 |
| `scripts/manuscript_build/build/build_supp_tables_rev.py` | 3 | 14, 15, 525 |
| `scripts/manuscript_build/build/rev_common.py` | 1 | 9 |
| `scripts/manuscript_build/build/rev_fig_export.R` | 2 | 2, 4 |
| `scripts/manuscript_build/build/rev_s4_export.R` | 3 | 2, 4 |
| `scripts/manuscript_build/build_v1/fig_export.R` | 3 | 2, 4, 49 |
| `scripts/manuscript_build/build_v1/fig_export_abo.R` | 3 | 2, 4, 6 |
| `scripts/rev_mr_00_prep.R` | 4 | 4, 7, 22 |
| `scripts/rev_mr_01_mr.R` | 8 | 7, 10, 11, 12, 13, 56 |
| `scripts/rev_mr_01b_raps.R` | 2 | 3, 28 |
| `scripts/rev_mr_02_fdr_concordance.R` | 4 | 4, 7 |
| `scripts/rev_mr_03_robust.R` | 13 | 4, 7, 8, 11, 54, 62, 75, 92, 103 |
| `scripts/revision/R1_10q26_insample.R` | 7 | 9, 11, 55, 62, 63, 187 |
| `scripts/revision2/rev2_10q26.py` | 1 | 13 |
| `scripts/revision2/rev2_concordant_blocks.R` | 6 | 8, 11, 12, 26 |
| `scripts/revision2/rev2_equivalence.py` | 1 | 9 |
| `scripts/revision2/rev2_remainder.py` | 2 | 8 |
| `scripts/revision3/panukb_01_fetch.py` | 1 | 8 |
| `scripts/revision3/panukb_02_analysis.py` | 1 | 9 |
| `scripts/revision3/panukb_03_susie.R` | 2 | 2, 4 |
| `scripts/revision3/panukb_04_followup.py` | 1 | 5 |
| `scripts/revision3/panukb_05_fourvar.py` | 1 | 7 |
| `scripts/revision3/r3_eye_build.py` | 1 | 5 |
| `scripts/revision3/r3_fetch.py` | 1 | 3 |
| `scripts/revision4/rev4_AB.py` | 1 | 17 |
| `scripts/revision4/rev4_C_coloc.R` | 3 | 3, 5, 9 |
| `scripts/revision4/rev4_C_extra.py` | 1 | 10 |
| `scripts/revision4/rev4_C_fetch.py` | 1 | 10 |
| `scripts/revision5/r5_00_fetch.py` | 1 | 8 |
| `scripts/revision5/r5_01_kappa_complement_af.R` | 13 | 8, 11, 12, 13, 14, 33, 38, 45, 50, 120 |
| `scripts/revision5/r5_02_contrast.py` | 1 | 11 |
| `scripts/revision5/r5_04_corivw.R` | 8 | 4, 6, 7, 8 |
| `scripts/revision5/r5_05_eqtl.py` | 1 | 8 |
| `scripts/story/story_01_mr.R` | 13 | 3, 6, 7, 8, 9, 13, 95, 175, 179, 181 |
| `scripts/story/story_02_tables.R` | 5 | 2, 4, 33, 46, 75 |

## 5. Copied scripts without hard-coded paths (not edited)

These already use relative paths (`Path(__file__).resolve().parents[1]` in `01_download_nhanes.py` and `02_nhanes_codebook.py`; `Path(__file__).parent` in `build_docx.py` and `build_af1.py`, which expect the manuscript markdown parts next to them), take paths as command-line arguments, or import their paths from `rev_common.py`:

- `scripts/01_download_nhanes.py`
- `scripts/02_nhanes_codebook.py`
- `scripts/20_gwascatalog_lookup.py`
- `scripts/21_check_sources.py`
- `scripts/helpers/cause_fix.R`
- `scripts/helpers/phewas_gc.py`
- `scripts/manuscript_build/build/check_text_colors.py`
- `scripts/manuscript_build/build/fig1.py`
- `scripts/manuscript_build/build/fig2.py`
- `scripts/manuscript_build/build/fig3.py`
- `scripts/manuscript_build/build/fig4.py`
- `scripts/manuscript_build/build/fig5.py`
- `scripts/manuscript_build/build/figS10.py`
- `scripts/manuscript_build/build/figS11.py`
- `scripts/manuscript_build/build/figS8.py`
- `scripts/manuscript_build/build/figS9.py`
- `scripts/manuscript_build/build/figs_supp.py`
- `scripts/manuscript_build/build_af1.py`
- `scripts/manuscript_build/build_docx.py`
- `scripts/revision3/r3_pdftext.py`

## 6. Deliberately not changed

- Executable names. PLINK, SMR and GCTA are called as Windows binaries under `<root>/data/raw/tools/` (`plink/plink.exe`, `smr-1.3.1-win-x86_64/smr-1.3.1-win.exe`, `gcta*.exe`). On Linux or macOS, put the corresponding binaries there under the same names or edit the `PLINK`/`SMR`/`GCTA` variable definitions (see `environment/external_tools.md`).
- Relative sub-paths below the root (`data/raw/...`, `derived/...`, `results/...`, `tmp/...`, `submission/...`). Outputs are written to the same relative locations as in the original project.
- `.libPaths(c(file.path(PROJ_ROOT, "Rlib"), ...))`. The project-local R library is still prepended; if `<root>/Rlib` does not exist R ignores it and uses the default libraries.
