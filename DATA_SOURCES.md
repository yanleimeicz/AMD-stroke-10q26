# Data sources

All data used in the study are public summary-level or de-identified data. **No raw or derived data are included in
this repository.** Download the files listed below into the paths given (relative to the project root,
`COMORBIDITY_ROOT`); several scripts download their inputs themselves (column "Fetched by").
Sample sizes, phenotype definitions, genome builds and PMIDs for every dataset are given in Additional file 2:
Table S1 of the manuscript; Table 1 of the manuscript summarises the main sources.

GWAS Catalog harmonised files follow the pattern
`http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/<GCST range>/<GCST>/harmonised/<file>.h.tsv.gz`
(the range folder is `GCSTxxxxx001-GCSTxxxx000`; `scripts/20_gwascatalog_lookup.py` lists the files of an accession).

## 1. Observational data

| Dataset | Content used | URL | Local path | Fetched by |
|---|---|---|---|---|
| NHANES 2005–2006 (cycle D) and 2007–2008 (cycle E) | DEMO, OPXRET (retinal imaging, AMD grading), VIQ, VIX, MCQ (stroke), BMX, BPX, BPQ, SMQ, ALQ, DIQ, KIQ_U, GHB, GLU, TCHOL, HDL, CRP, BIOPRO, RXQ_RX, DPQ, PAQ, HIQ, SLQ (XPT files) and their documentation pages | https://wwwn.cdc.gov/nchs/nhanes/ (files: `https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{2005,2007}/DataFiles/<COMP>_{D,E}.xpt`) | `data/raw/nhanes/` | `01_download_nhanes.py`, `02_nhanes_codebook.py` |

## 2. GWAS summary statistics

### AMD

| Trait | Study | Accession / file | URL | Local path |
|---|---|---|---|---|
| Any AMD; dry AMD; wet AMD | FinnGen release 12 | endpoints `H7_AMD`, `DRY_AMD`, `WET_AMD` | https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_<ENDPOINT>.gz ; access page https://www.finngen.fi/en/access_results | `data/raw/gwas/finngen_R12_<ENDPOINT>.gz` |
| Advanced AMD | IAMDGC (Fritsche et al. 2016) | GCST003219 | …/GCST003001-GCST004000/GCST003219/harmonised/26691988-GCST003219-EFO_0001365.h.tsv.gz | `data/raw/gwas/IAMDGC_GCST003219.h.tsv.gz` |
| 52 conditionally independent AMD variants | IAMDGC | GWAS Catalog associations of GCST003219 | GWAS Catalog REST API | `results/revision/mr/input/iamdgc_52_independent_gwascatalog.tsv` (small input table) |
| Early AMD | Winkler et al. 2020 | GCST010723 | …/GCST010001-GCST011000/GCST010723/harmonised/32843070-GCST010723-EFO_0001365.h.tsv.gz | `data/raw/gwas/EARLYAMD_Winkler2020_GCST010723.h.tsv.gz` |

### Stroke

| Trait | Study | Accession | Local path |
|---|---|---|---|
| Any, ischemic, cardioembolic, large-artery, small-vessel stroke (European) | GIGASTROKE (Mishra et al. 2022) | GCST90104539 (AS), GCST90104540 (IS), GCST90104541 (CES), GCST90104542 (LAS), GCST90104543 (SVS) | `data/raw/gwas/GIGASTROKE_<GCST>.h.tsv.gz` |
| Any, ischemic, large-artery, small-vessel stroke (East Asian) | GIGASTROKE | GCST90104544 (AS), GCST90104545 (IS), GCST90104547 (LAS), GCST90104548 (SVS) | `data/raw/gwas/GIGASTROKE_EAS_<trait>_<GCST>.h.tsv.gz` |
| Any, large-artery, ischemic, small-vessel, cardioembolic stroke (European; FinnGen-free) | MEGASTROKE (Malik et al. 2018) | GCST006906 (AS), GCST006907 (LAS), GCST006908 (IS), GCST006909 (SVS), GCST006910 (CES) | `data/raw/gwas/MEGASTROKE_<trait>_EUR_<GCST>.h.tsv.gz` (also copied to `data/raw/gwas/cardiometab/MEGASTROKE_<trait>.h.tsv.gz` for the 3p14 scripts) |
| Stroke endpoints | FinnGen R12 | `I9_STR` (stroke excl. SAH), `C_STROKE`, `I9_TIA`, `I9_ICH`, `I9_SAH`, `I9_STR_EMBOLIC` (fine-mapping only) | `data/raw/gwas/finngen_R12_<ENDPOINT>.gz` |

### Cerebral small vessel disease and related traits

| Trait | Study | Accession / file | Local path |
|---|---|---|---|
| White matter hyperintensity volume | UK Biobank BIG40 (Smith et al. 2021), IDP 1437 | https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/1437.txt.gz | `data/raw/gwas/BIG40_1437.txt.gz` |
| Lacunar stroke | Traylor et al. 2021 | GCST90014122 | `data/raw/gwas/LACUNAR_Traylor2021_GCST90014122.h.tsv.gz` |
| Enlarged perivascular spaces (WM, BG, HIP) | Duperron et al. 2023 | GCST90244151, GCST90244153, GCST90244155 | `data/raw/gwas/EPVS_<region>_Duperron2023_<GCST>.h.tsv.gz` |
| White-matter mean diffusivity, fractional anisotropy | Persyn et al. 2020 | GCST010103 (`GWAS_UKB_PC1_MD.zip`), GCST010102 (`GWAS_UKB_PC1_FA.zip`) | `data/raw/gwas/WM_{MD,FA}_Persyn2020_<GCST>.zip` |
| Intracerebral hemorrhage | UK Biobank whole-genome sequencing | GCST90473584 | `data/raw/gwas/ICH_UKBWGS_GCST90473584.h.tsv.gz` |

### Cardiometabolic and vascular traits (locus-level analyses, positive controls)

| Trait | Study | Accession / file | Local path / fetched by |
|---|---|---|---|
| Coronary artery disease | Aragam et al. 2022 | GCST90132314 (harmonised; remote tabix and full download) | `derived/story/GCST90132314.h.tsv.gz` (+ `.tbi` in `tmp/rev_adamts9/`); `75a_rev_adamts9_fetch.R`, `story/story_01_mr.R` |
| Coronary artery disease (earlier) | van der Harst & Verweij 2018 | GCST005195 | `data/raw/gwas/cardiometab/CAD_vdHarst2018.h.tsv.gz` |
| Type 2 diabetes | Xue et al. 2018 | GCST006867 | `data/raw/gwas/cardiometab/T2D_Xue2018.h.tsv.gz` |
| WHR adjusted for BMI | Pulit et al. 2019 | GCST008999 (`fat-distn.giant.ukbb.meta-analysis.whradjbmi.combined.tbl.gz`) | `data/raw/gwas/cardiometab/WHRadjBMI_Pulit2019.tbl.gz` |
| Systolic / diastolic blood pressure | Evangelou et al. 2018 | GCST006624 (`Evangelou_30224653_SBP.txt.gz`), GCST006630 (`Evangelou_30224653_DBP.txt.gz`) | `data/raw/gwas/cardiometab/{SBP,DBP}_Evangelou2018.txt.gz` |
| HbA1c | Chen et al. 2021 (MAGIC) | GCST90002244 | `data/raw/gwas/cardiometab/HbA1c_Chen2021.h.tsv.gz` |
| LDL / HDL cholesterol | GLGC (Graham et al. 2021), European | GCST90239658 (LDL), GCST90239652 (HDL) | remote tabix in `75a_rev_adamts9_fetch.R` |
| cIMT; carotid plaque | Franceschini et al. 2018 | GCST007436, GCST007435 | `74b_adamts9_cimt.R` (remote) |
| Atrial fibrillation (positive control) | Nielsen et al. 2018 | GCST006414 (GWAS Catalog associations; fallback Roselli et al. 2018, GCST006061) | `revision5/r5_00_fetch.py` |
| FinnGen R12 endpoints for 3p14 and 10q26 look-ups | FinnGen R12 | `T2D`, `I9_STR`, `C_STROKE`, `I9_ICH`, `I9_CHD`, `I9_CORATHER`, `I9_HYPTENS`, `I9_MI_STRICT`, `I9_TIA`, `I9_SAH` (case/control counts in `config/finngen_counts.tsv`; copy it to `tmp/rev_adamts9/finngen_counts.tsv`) | `75a_rev_adamts9_fetch.R` (remote tabix) |

### FinnGen fine-mapping and LD

| Resource | URL | Fetched by |
|---|---|---|
| FinnGen R12 SuSiE fine-mapping (`.SUSIE.snp.bgz`, `.SUSIE.cred.bgz`) for `H7_AMD`, `DRY_AMD`, `WET_AMD`, `T2D`, `E4_DM2`, `T2D_WIDE`, `KELA_DIAB_INSUL_EXMORE`, `I9_STR`, `I9_STR_EMBOLIC`, `C_STROKE`, `I9_CORATHER`, `I9_CHD`, `I9_HYPTENS`, `I9_HYPTENSESS`, `I9_MI_STRICT`, `I9_REVASC` | https://storage.googleapis.com/finngen-public-data-r12/finemap/full/susie/finngen_R12_<ENDPOINT>.SUSIE.{snp,cred}.bgz | `75b_rev_adamts9_finngen_susie.R`, `revision/R1_10q26_insample.R` |
| FinnGen in-sample LD (SISu v4.2 reference) for key 10q26 variant pairs | FinnGen LD server (https://ld.finngen.fi) | queried interactively; stored as `results/revision/10q26/ld_finngen_sisu42_pairs.json` (input to `revision2/rev2_10q26.py`, `revision3/panukb_02_analysis.py`, `revision4/rev4_AB.py`, `revision5/r5_02_contrast.py`) |

### Pan-UK Biobank (in-sample LD analyses at 10q26 and 3p14)

| Resource | Details | URL | Fetched by |
|---|---|---|---|
| Pan-UKB per-phenotype summary statistics (European ancestry, GRCh37) | phecode 362.29 (AMD, primary), 362.2, self-reported 6148-5; ICD-10 I63 (ischemic stroke, primary), phecodes 433.2, 433.3, 433, 430.2, ICD-10 G45, self-reported 6150-3 and 20002-1081; phecode 411.4 (coronary atherosclerosis), 411; phecode 250.2 and ICD-10 E11 (type 2 diabetes) | https://pan-ukb-us-east-1.s3.amazonaws.com (`sumstats_flat_files/…`, `phenotype_manifest.tsv.bgz`); https://pan.ukbb.broadinstitute.org | `revision3/panukb_01_fetch.py` (10q26), `revision4/rev4_C_fetch.py` (3p14); region-only range requests |
| UK Biobank in-sample LD matrices (Weissbrod et al. 2020; British ancestry, n ≈ 337,000) | blocks `chr10_123000001_126000001`, `chr3_64000001_67000001` (`.npz` + `.gz`) | https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/ | same scripts |

## 3. QTL data

| Resource | Content used | URL | Local path / fetched by |
|---|---|---|---|
| eQTLGen phase I | whole-blood cis-eQTL (significant, FDR < 0.05), allele frequencies, SMR-formatted BESD (`cis-eQTL-SMR_20191212.tar.gz`) | https://www.eqtlgen.org/ ; https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/ (`2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz`, `2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz`, `SMR_formatted/cis-eQTL-SMR_20191212.tar.gz`) | `data/raw/eqtl/` (`eqtlgen_cis_sig.txt.gz`, `eqtlgen_besd.*`; `38_gunzip_besd.py`) |
| GTEx v8 (via eQTL Catalogue) | gene-expression QTL summary statistics, datasets QTD000131 (aorta), QTD000136 (coronary artery), QTD000141 (tibial artery), QTD000171/176/156/196/181/166 (brain), QTD000356 (whole blood), QTD000216 (fibroblasts); aortic SuSiE credible sets (QTD000131) | https://www.ebi.ac.uk/eqtl/ ; dataset paths from `tabix_ftp_paths.tsv` in https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources (save as `tmp/tissue/tabix_paths.tsv`) | remote tabix in `48a_fetch_tissue_eqtl.R`, `74a_adamts9_data.R`; aortic credible sets in `data/raw/eqtl/susie_aorta/` |
| ROSMAP and BrainSeq DLPFC (via eQTL Catalogue) | QTD000434 (ROSMAP), QTD000051 (BrainSeq) | as above | as above |
| EyeGEx neural retina (Ratnapriya et al. 2019) | nominal cis-eQTL by chromosome and permutation results | Zenodo record 18331341, e.g. https://zenodo.org/api/records/18331341/files/new_nominal.eQTL.chr10.with_thresholds.txt.gz/content | `data/raw/eqtl/eyegex/retina_nominal_chr<N>.txt.gz`, `retina_perms_full.txt` |
| Retina / RPE / choroid single-nucleus eQTL (Voigt et al. 2026) | Supplementary Table 3 (all eSNPs with FDR < 0.05, 30 cell types) | bioRxiv doi:10.64898/2026.03.30.714946 (v2), supplementary file `media-1.csv` | `revision3/r3_fetch.py` with `revision3/urls*.tsv` → `derived/eye_eqtl/web/` |
| Retina and RPE–choroid bulk eQTL (Orozco et al. 2020) | significant eQTL look-ups at 10q26 | https://eye-eqtl.com (gene and variant query pages) | query pages saved as `derived/eye_eqtl/web/eyeeqtl/<query>.html` (saved outside the released scripts; parsed by `revision3/r3_eye_build.py`) |
| UKB-PPP plasma cis-pQTL (Sun et al. 2023) | Supplementary Tables 3 and 9 | https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06592-6/MediaObjects/41586_2023_6592_MOESM3_ESM.xlsx | `data/raw/pqtl/UKBPPP_Sun2023_SuppTables.xlsx` → `helpers/pqtl_xlsx2tsv.py` |
| deCODE plasma pQTL (Ferkingstad et al. 2021) | Supplementary Table 2 | https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-021-00978-w/MediaObjects/41588_2021_978_MOESM4_ESM.xlsx | `data/raw/pqtl/deCODE_Ferkingstad2021_SuppTables.xlsx` → `helpers/pqtl_xlsx2tsv.py` |
| DGIdb gene categories | druggable-genome annotation of SMR genes | https://dgidb.org/api/graphql | `38c_dgidb_categories.py` |

## 4. Single-cell data (CZ CELLxGENE Discover)

| Dataset | CELLxGENE dataset ID (Table S1) | h5ad download used | Local path |
|---|---|---|---|
| RPE–choroid atlas | 5ed90bc0-12f3-4354-b314-0ce0dbb96458 | https://datasets.cellxgene.cziscience.com/68493594-1bde-4b8b-8916-473a10caf035.h5ad | `data/raw/sc/rpe_choroid_scrna_all.h5ad` |
| Human retina, periphery (Cowan et al. 2020) | ccfdcad0-7104-46b9-addf-fd66a2a15907 | https://datasets.cellxgene.cziscience.com/052ea638-dead-4657-a01b-3c0278bea59f.h5ad | `data/raw/sc/retina_rpe_periphery_cowan.h5ad` |
| Human Brain Cell Atlas v1.0, vascular supercluster (Siletti et al. 2023) | 0d7f4c06-a6bd-47d2-a42b-3a7196704f77 | https://datasets.cellxgene.cziscience.com/3df231b4-4c8e-4d7b-b854-9c27223d0738.h5ad | `data/raw/sc/brain_vascular_siletti.h5ad` |
| Brain vascular multi-omics (Reid et al. 2025) | 203025fe-fa99-4d57-81da-458ed8f0c334 | https://datasets.cellxgene.cziscience.com/fd50a30c-3d7e-47c2-9689-fab95ea6f20a.h5ad | `data/raw/sc/brain_vascular_multiomics_snrna.h5ad` |

Dataset metadata were searched with the CELLxGENE curation API (https://api.cellxgene.cziscience.com/curation/v1/datasets; `50_cxg_search.py`).
The h5ad asset URLs above are version-specific; if they have changed, download the same dataset IDs from https://cellxgene.cziscience.com/.

## 5. Reference panels

| Resource | URL | Local path |
|---|---|---|
| 1000 Genomes phase 3 European reference (n = 503; PLINK format) | https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe/download (MAGMA `g1000_eur.zip`); project page https://www.internationalgenome.org/ | `data/raw/ref/g1000_eur/` |
| 1000 Genomes EUR/FIN/EAS pairwise LD for the key 10q26 variants | Ensembl REST API (https://rest.ensembl.org) | `results/revision/10q26/ld_ensembl_EUR_FIN_EAS.tsv`, `ld_1000G_EUR_r_dprime.tsv` (small input tables) |
| LD scores (European, HapMap3 SNP list) | https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/eur_w_ld_chr.tar.bz2, `w_hm3.snplist.bz2` | `data/raw/ref/eur_w_ld_chr/` |
| LAVA LD blocks (2,495 blocks, GRCh37) | https://github.com/josefin-werme/LAVA (`support_data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile`) | `data/raw/ref/lava/` |
| UK Biobank LD matrices (Weissbrod et al. 2020) | https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/ | fetched by the Pan-UKB scripts |

## 6. Small result tables used directly by the figure scripts

The figure and table builders (`scripts/manuscript_build/`) read the result TSVs written by the analysis scripts
under `results/` (for example `results/revision/mr/mr_revised_all.tsv`, `results/revision5/causal_expectation_respecified.tsv`,
`results/story/pathway_mr.tsv`, `results/revision5/10q26_stroke_contrast_hksj.tsv`) and plotting extracts in
`submission/build/figdata/` and `submission/revision/build/figdata/` (written by the `fig_export*.R`,
`rev_fig_export.R` and `rev_s4_export.R` scripts). These files are not included here; they are reproduced by
running the pipeline, and every number shown in the figures is also in Additional file 2.
