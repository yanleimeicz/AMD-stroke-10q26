# Code for "Fine-mapping of the ARMS2–HTRA1 locus explains the apparent protective effect of genetic liability to age-related macular degeneration on stroke"

This repository contains the analysis code for the study above (submitted to *HGG Advances*). It reproduces the
population analysis in NHANES, the Mendelian randomization (MR), genetic-correlation, fine-mapping, conditional,
colocalization, QTL and single-cell analyses, and the figures and supplementary tables of the manuscript.
No data are included; all inputs are public and listed in [DATA_SOURCES.md](DATA_SOURCES.md).

- **Citation:** Zhang Y, Zhao M, Song W, Yan L. Fine-mapping of the ARMS2–HTRA1 locus explains the apparent protective effect of genetic liability to age-related macular degeneration on stroke. *HGG Advances* (submitted). Please cite the published article once available.
- **License:** MIT (see [LICENSE](LICENSE)). The license covers the code only; the downloaded data remain under the terms of their providers (FinnGen, UK Biobank/Pan-UKB, GWAS Catalog, NHANES, CELLxGENE, etc.).
- **Contact:** Lei Yan, Department of Neurology, The Affiliated Brain Hospital of Nanjing Medical University (leiyan@njmu.edu.cn).
- **Use of AI tools:** as stated in the manuscript (Acknowledgments), an artificial intelligence assistant (Claude Code, Anthropic) was used to help write the analysis code. The authors reviewed the code and are responsible for all analytical decisions (data sets, instruments, LD reference panels, colocalization priors, sensitivity analyses) and for the results.

## Contents

| Path | Contents |
|---|---|
| `scripts/` | Analysis scripts (R and Python), numbered roughly in run order; `revision*/` and `story/` hold analyses added during peer review |
| `scripts/helpers/` | Small helpers: patched CAUSE function (`cause_fix.R`), GWAS Catalog PheWAS look-up (`phewas_gc.py`), pQTL supplementary-table converter (`pqtl_xlsx2tsv.py`) |
| `scripts/manuscript_build/` | Figure, supplementary-table, STROBE-MR and DOCX builders (`build/`, `build_v1/`, `build_docx.py`, `build_af1.py`) |
| `config/` | Small input lists: genes for the single-cell look-up (`scrna_genes.txt`), FinnGen R12 endpoint case/control counts (`finngen_counts.tsv`), rsID lists for the ocular-eQTL LD (`eye_eqtl/`) |
| `environment/` | `python_requirements.txt`, `R_packages.tsv`, `external_tools.md` |
| `DATA_SOURCES.md` | Every public dataset with accession, download URL and expected local path |
| `PATH_CHANGES.md` | How hard-coded paths of the original scripts were replaced by `COMORBIDITY_ROOT` |
| `tools/` | `make_release.py` (builds this folder from the project), `path_changes.tsv` (line-level log) |

## Setup

1. **Software.** R 4.5.2 with the packages in `environment/R_packages.tsv`; Python 3.14 with
   `environment/python_requirements.txt`; a separate Python 3.11 environment with Scanpy 1.11.5 for the two
   single-cell scripts; PLINK 1.9, SMR 1.3.1 and GCTA 1.94.1 (`environment/external_tools.md`). The analyses were run
   on Windows 11; the scripts call the Windows executables `plink.exe`, `smr-1.3.1-win.exe` and `gcta*.exe` under
   `data/raw/tools/` (see `environment/external_tools.md` for other systems).
2. **Project root.** All scripts resolve paths from one root directory, read from the environment variable
   `COMORBIDITY_ROOT`. If it is not set, the root defaults to the parent of the `scripts/` folder (this repository).
   The root must contain (or link to) the data folders:

   ```
   <root>/
     scripts/            # this repository's scripts (the default root is the repository itself)
     data/raw/{nhanes,gwas,gwas/cardiometab,eqtl,pqtl,sc,ref,tools}/   # downloads, see DATA_SOURCES.md
     derived/            # intermediate files (created by the scripts; ~56 GB in the original run)
     results/            # result tables and figures (created by the scripts)
     tmp/                # scratch files (created by the scripts)
     submission/         # figure/table outputs of scripts/manuscript_build (created)
   ```

   ```bash
   # Linux / macOS
   export COMORBIDITY_ROOT=/path/to/project
   # Windows (cmd)
   set COMORBIDITY_ROOT=D:\path\to\project
   ```
   Run R scripts with `Rscript scripts/<name>.R` and Python scripts with `python scripts/<name>.py`. Optional: a
   project-local R library in `<root>/Rlib` is put first on `.libPaths()` if it exists.
3. **Small inputs that the scripts expect under `tmp/` or `results/`.** Copy `config/finngen_counts.tsv` to
   `<root>/tmp/rev_adamts9/finngen_counts.tsv`; save the eQTL Catalogue `tabix_ftp_paths.tsv` as
   `<root>/tmp/tissue/tabix_paths.tsv`; see "Known gaps" below for the small 10q26 LD/fine-mapping input tables.
4. Keep temporary files on the data disk: point `TMPDIR`/`TEMP`/`TMP`, `MPLCONFIGDIR` and `JOBLIB_TEMP_FOLDER`
   to `<root>/tmp`.

## Run order and correspondence to the manuscript

Steps are listed in dependency order. Within a step, run the scripts in the order given. "S" numbers refer to
Additional file 1 (Figures S1–S11, Supplementary Methods/Notes) and Additional file 2 (Tables S1–S18).

| Step | Scripts | Main outputs | Manuscript |
|---|---|---|---|
| **0. Data download** | `01_download_nhanes.py`, `02_nhanes_codebook.py`; GWAS/QTL/single-cell files per `DATA_SOURCES.md` (helpers: `20_gwascatalog_lookup.py`, `21_check_sources.py`); `38_gunzip_besd.py` (eQTLGen BESD); `helpers/pqtl_xlsx2tsv.py` (UKB-PPP, deCODE); `50_cxg_search.py` (CELLxGENE metadata) | `data/raw/`, `derived/pqtl/` | Table 1; Table S1 |
| **1. NHANES** | `10_nhanes_stroke_amd.R` → `11_nhanes_forest_plots.R` | `results/nhanes/` | Results §1; Fig. 2a; Fig. S1; Table S2 |
| **2. GWAS harmonization** | `30_prepare_gwas.R`, `30b_check_iamdgc_sign.R` (check), `45_csvd_gwas_prep.R`, `rev_mr_00_prep.R` (IAMDGC P < 1e-300 recovery, MEGASTROKE) | `derived/gwas/`, `derived/revision/` | Methods "GWAS summary statistics" |
| **3. First-round MR, colocalization, LDSC, LAVA** (inputs to later steps) | `31_mr_bidirectional.R`, `36_mr_steiger.R`, `37_mr_coloc_figures.R`; `32_coloc_gwas.R`, `35_coloc_post.R`; `34_ldsc_rg.R`; `46_csvd_mr_ldsc.R`; `47_lava_local_rg.R` | `results/mr/`, `results/coloc/`, `results/ldsc/`, `results/csvd/` | Tables S6–S8; Fig. 3c; Fig. S3; Fig. S10 |
| **4. Revised MR (primary)** | `rev_mr_01_mr.R` → `rev_mr_01b_raps.R` → `rev_mr_02_fdr_concordance.R` → `rev_mr_03_robust.R` (MR-Clust, MRlap, CAUSE with `helpers/cause_fix.R`, LAVA sensitivity at 10q26) | `results/revision/mr/` | Results §1–2; Fig. 2b; Fig. 3a–c; Fig. S2, S10; Tables S3–S5, S8, S9, S16 |
| **5. Expected effect under causality, equivalence, positive controls, locus groups** | `revision2/rev2_equivalence.py`; `story/story_01_mr.R`, `story/story_02_tables.R`; `revision5/r5_00_fetch.py` → `revision5/r5_01_kappa_complement_af.R` → `revision5/r5_04_corivw.R`; `revision4/rev4_AB.py` (part B) | `results/revision2/`, `results/story/`, `results/revision5/`, `results/revision4/` | Results §1 and §5; Fig. 2b–c; Fig. 3a; Fig. 5; Table S5, S17 (complement pQTL); Supplementary Note 1 |
| **6. Extra outcomes and loci** | `76_rev_extra.R` (ICH, early AMD, East Asian key variants; needs `46_csvd_mr_ldsc.R`), `75_c4_conditional.R` (C4) | `results/revision/extra/` | Table 2; Tables S6, S10, S11; Fig. S9, S10 |
| **7. 10q26 locus** | `revision/R1_10q26_insample.R` (FinnGen in-sample SuSiE, coloc.bf_bf, conditional/proportionality, prior sensitivity) → `revision2/rev2_10q26.py` → `revision2/rev2_remainder.py` → `revision3/panukb_01_fetch.py` → `panukb_02_analysis.py` → `panukb_03_susie.R` → `panukb_04_followup.py` → `panukb_05_fourvar.py` → `revision4/rev4_AB.py` (part A, random-effects meta-analysis) → `revision5/r5_02_contrast.py` (Hartung–Knapp, haplotype contrast) | `results/revision/10q26/`, `results/revision2/`, `results/revision3/panukb/`, `results/revision4/`, `results/revision5/` | Results §3; Table 2; Fig. 4a–c, 4e; Fig. S3, S11; Tables S11, S12; Supplementary Note 2 |
| **8. Expression and protein QTL** | Tissue eQTL: `48a_fetch_tissue_eqtl.R` → `48b_gwas_regions.R` → `48_tissue_eqtl.R` → `49_cojo_susie.R`. Ocular eQTL: `revision3/r3_fetch.py` with `revision3/urls1.tsv`–`urls3.tsv`, `revision3/r3_pdftext.py`, 1000G LD with PLINK (below), `revision3/r3_eye_build.py`; `revision5/r5_05_eqtl.py` (conditional eQTL in both directions). SMR–HEIDI: `33_eqtl_smr.R` → `38a_make_ma.R` → `38b_ld_subset.R` → `38_smr_heidi_druggable.R` → `38c_dgidb_categories.py`. pQTL MR: `39_pqtl_mr.R` | `results/tissue_eqtl/`, `results/revision3/eye_eqtl/`, `results/smr/`, `results/smr_heidi/`, `results/pqtl_mr/` | Results §4–5; Fig. 4d; Fig. S4 (12q24), S5, S6; Tables S11 (12q24), S13, S17; Supplementary Note 3 |
| **9. 3p14 (ADAMTS9) locus** | `74a_adamts9_data.R` → `74b_adamts9_cimt.R` → `74_adamts9.R`; `75a_rev_adamts9_fetch.R` → `75b_rev_adamts9_finngen_susie.R` → `75c_rev_adamts9_mvmr.R` → `75d_rev_adamts9_locus.R` (HyPrColoc; needs `Rlib_hypr` or hyprcoloc installed); `helpers/phewas_gc.py`; `revision4/rev4_C_fetch.py` → `rev4_C_extra.py` → `rev4_C_coloc.R` (Pan-UKB replication) | `results/adamts9/`, `results/revision/adamts9/`, `results/revision4/` | Results §6; Table 2; Fig. S8; Tables S14, S15; Supplementary Note 4 |
| **10. LAVA concordant blocks follow-up** | `revision2/rev2_concordant_blocks.R` | `results/revision2/concordant_blocks_*.tsv` | Results §2; Table S8 |
| **11. Single cell** (Python 3.11 / Scanpy environment) | `50_scrna.py config/scrna_genes.txt`; `52_scrna_adamts9.py` | `results/scrna/`, `results/adamts9/scrna/` | Methods "Single-cell localization"; Fig. S7; Table S18; Supplementary Methods S6 |
| **12. Figures, tables, documents** | Plot extracts: `manuscript_build/build_v1/fig_export.R`, `fig_export_abo.R`, `manuscript_build/build/rev_fig_export.R`, `rev_s4_export.R`. Figures: `build/fig1.py` (schematic), `fig2.py`, `fig3.py`, `fig4.py`, `fig5.py`, `figs_supp.py` (Figs. S1–S7), `figS8.py`–`figS11.py` (they import `rev_common.py` from their own folder). Tables: `build/build_supp_tables_rev.py` (Additional file 2), `build/build_strobe_mr_rev.py` (Additional file 3). QA: `build/check_text_colors.py`. Documents: `build_af1.py`, `build_docx.py` (need the manuscript markdown parts next to them) | `submission/revision/figures/`, `submission/revision/additional_files/` | Figs. 1–5; Figs. S1–S11; Additional files 1–3 |

PLINK commands for the ocular eQTL LD (step 8), run in `<root>/derived/eye_eqtl` after copying the SNP lists
`config/eye_eqtl/ld_snps.txt` and `ld2_snps.txt` there, with `REF=<root>/data/raw/ref/g1000_eur/g1000_eur`:

```
plink --bfile $REF --extract ld_snps.txt  --r --matrix --out ld_eur_r --allow-no-sex
plink --bfile $REF --extract ld_snps.txt  --r2 dprime-signed --ld-window-r2 0 --ld-window 99999 --ld-window-kb 99999 --out ld_eur_pairs --allow-no-sex
plink --bfile $REF --extract ld_snps.txt  --freq --out ld_eur_frq --allow-no-sex
plink --bfile $REF --extract ld2_snps.txt --r2 dprime-signed --ld-window-r2 0 --ld-window 99999 --ld-window-kb 99999 --out ld2_eur_pairs --allow-no-sex
plink --bfile $REF --extract ld2_snps.txt --r --ld-window-r2 0 --ld-window 99999 --ld-window-kb 99999 --out ld2_eur_r --allow-no-sex
plink --bfile $REF --extract ld_snps.txt  --r --ld-window-r2 0 --ld-window 99999 --ld-window-kb 99999 --out ld_eur_rsigned --allow-no-sex
```

## Hardware and runtime

The analyses were run on a workstation with a 16-core/32-thread CPU (AMD Ryzen 9 9950X3D), 128 GB RAM and an NVIDIA
RTX 5090 GPU (the GPU is not needed for the reported analyses). Raw downloads take about 32 GB (GWAS 18 GB, eQTL
6.6 GB, reference panels 3.9 GB, single-cell 3.5 GB) and intermediate files about 56 GB. Several R scripts start
PSOCK clusters with up to 30 workers (`makeCluster(...)`); reduce these numbers on smaller machines. Most scripts
finish in minutes. The longest steps are the genome-wide preparation and LD clumping (step 2–3), transcriptome-wide
SMR–HEIDI and eQTL colocalization (step 8), LAVA across 2,495 blocks (step 3–4), CAUSE (step 4) and the remote
tabix/range downloads (steps 7–9), each of which can take from tens of minutes to several hours depending on CPU
and network speed. A full re-run takes roughly one to two days on comparable hardware.

## Known gaps

The following small input tables used by the 10q26 and 3p14 analyses were produced by interactive commands during
revision, and no script for them is included. Their provenance is described in `DATA_SOURCES.md`, and the values are
reported in Additional file 2:

- `results/revision/10q26/ld_finngen_sisu42_pairs.json` (FinnGen LD server, SISu v4.2), `ld_ensembl_EUR_FIN_EAS.tsv`
  (Ensembl REST LD), `ld_1000G_EUR_r_dprime.tsv` (PLINK `--r dprime` on the 1000G EUR panel);
- `results/revision/10q26/finngen_susie_key_variants.tsv` and `finngen_susie_credsets_10q26.tsv` (extracts of the
  FinnGen R12 SuSiE files for the AMD endpoints at 10q26), `htra1_pqtl_mr_decode.tsv` (Wald ratio of the deCODE HTRA1
  cis-pQTL rs61871680);
- `results/revision/mr/input/iamdgc_52_independent_gwascatalog.tsv` (GWAS Catalog associations of GCST003219);
- `results/revision4/3p14_amd_replication.tsv` (comparison of the Pan-UKB, FinnGen and IAMDGC rs4611812 estimates);
- `derived/eye_eqtl/web/eyeeqtl/*.html` (eye-eqtl.com query pages).

## Analyses not reported in the manuscript

An earlier version of this work included exploratory analyses that were removed from the manuscript because they
did not bear on its conclusions (manuscript Table 3, footnote). Their scripts are **not** part of this release:

| Analysis | Scripts (original project `scripts/`) |
|---|---|
| Transcriptome (GEO) differential expression | `40_geo_parse.R`, `40_geo_probe.py`, `41_geo_deg.R` |
| Machine-learning classification | `60_ml.py` |
| Gene integration and protein–protein interaction hubs | `70_integration.py`, `71_ppi_hub.py`, `73_disease_gene_intersection.py` |
| Receiver operating characteristic (single-gene ROC) | `72_roc_singlegene.py` |
| Molecular docking | `80_drug_targets.py`, `81_docking.py`, `82_docking_curated.py`, `83_redock_rmsd.py`, `84_fix_cdk6_loops.py` |
| Molecular dynamics and MM/GBSA | `90_md.py`, `91_md_analysis.py`, `92_mmgbsa.py` |
| Single-cell high- versus low-expression differential expression | `51_scrna_highlow.py` |

Also not included: project-scoping and literature-screening scripts (`00_extract_pdf_text.py`,
`03_screen_nhanes_neuro_eye.py`, `04_pubmed_saturation.py`, `05_finngen_endpoints.py`, `06_pubmed_titles.py`),
installation and diagnostic utilities (`install_r_pkgs*.R`, `install_r_tissue.R`, `check_r_env.R`, `debug_sc*.py`,
`dl.py`, `dl_seg.py`, `gpu_test.py`, `inspect_*.py`, `sizes.py`, `speedtest.py`, `revision2/rev2_peek*.py`,
`revision3/r3_check.py`) and superseded figure scripts of the first submission.

## Rebuilding this folder

`python tools/make_release.py` (run inside the original project) recopies the scripts and reapplies the path
replacements; `python tools/make_path_changes_md.py` (run from `code_release/`) regenerates `PATH_CHANGES.md`.

## Round-6 additions (sample-overlap audit, independent coronary data, power)

| Step | Scripts (`scripts/revision6/`) | Output (`results/revision6/`) | Manuscript |
|---|---|---|---|
| Overlap audit of GIGASTROKE, MEGASTROKE, FinnGen and UK Biobank | `r6_00_giga_table1.py`, `r6_01_pdownload.py`, `r6_02_peek.R` | `overlap_matrix.tsv` | Methods; Table S12 |
| Null-variant Z-statistic correlations between datasets | `r6_03_overlap_corr.R` | `null_z_correlations.tsv` | Methods; Table 4; Table S12 |
| FinnGen coronary endpoints at 10q26 (in-sample LD) | `r6_04_fg_cad_region.R` | `cad_independent_validation.tsv` | Results (10q26); Fig. 4c |
| Power and simulation of the two-haplotype model | `r6_05_power_sim.py` | `contrast_power.tsv`, `simulation_10q26.tsv` | Results; Supplementary Note 5 |
| GLS meta-analysis modeling overlap; heterogeneity table; phenotype-matched MR; benchmark differences; sign-concordance inference | `r6_06_pool.py` | `10q26_stroke_gls_overlap.tsv`, `heterogeneity_main_table.tsv`, `insample_only_pool.tsv`, `phenotype_matched_mr.tsv`, `benchmark_difference_ci.tsv`, `sign_concordance_inference.tsv` | Table 4; Table S5; Table S12 |

## Round-7 additions (GIGASTROKE remainder construction)

| Step | Scripts (`scripts/revision7/`) | Output (`results/revision7/`) | Manuscript |
|---|---|---|---|
| Per-variant inverse-variance subtraction of MEGASTROKE, FinnGen and UK Biobank from GIGASTROKE European (information shares, estimates before and after, proxy error) and the two-haplotype joint model of each component | `r7_01_subtraction_table.py` (reuses the definitions of `revision6/r6_06_pool.py`) | `remainder_construction.tsv`, `remainder_joint_before_after.tsv` | Supplementary Methods S4; Supplementary Note 2; Table S12 |
| Loci shared by AMD and stroke/CSVD outside 10q26 (moved from main Table 3) | — (transcribed from the former Table 3) | `other_shared_loci_summary.tsv` | Table S10 |


## HGG Advances version (October 2026)

`scripts/hggadv/` contains the code added for the HGG Advances version: `h01_mr_diagnostics.R` (instrument diagnostics for Table 2 and Table S19: F statistics, Cochran's Q, I², MR-Egger intercept, I²_GX and MR-PRESSO, run in parallel), `figures/` (Figures 1–4 and Figures S12–S13 with uppercase panel labels; they reuse `scripts/manuscript_build/build/rev_common.py`) and `h02_bbj_10q26.py` and `h03_bbj_figure.py` (East Asian replication in BioBank Japan, Table S20 and Figure S14). The scripts that assemble the manuscript files are not included. All other analyses are unchanged.
