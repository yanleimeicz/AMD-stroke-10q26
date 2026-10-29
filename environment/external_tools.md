# External command-line tools and reference files

All tools were run on Windows 11 (x86_64). The scripts expect them at the paths below, relative to the project
root (`COMORBIDITY_ROOT`). On Linux or macOS, install the native builds and either place them under the same
names or edit the `PLINK`, `SMR` and `GCTA` variable definitions at the top of the scripts listed.

| Tool | Version used | Expected location (relative to root) | Download | Used by |
|---|---|---|---|---|
| PLINK 1.9 | v1.9.0-b.7.11 64-bit (19 Aug 2025) | `data/raw/tools/plink/plink.exe` | https://www.cog-genomics.org/plink/1.9/ (Windows build `plink_win64_20250819.zip` from https://s3.amazonaws.com/plink1-assets/) | clumping, LD matrices, allele frequencies, proxies: 31, 32, 33, 35, 38a, 38b, 38_smr_heidi, 39, 46, 47, 49, 74, 74a, 75_c4, 75c, rev_mr_01, rev_mr_03, r5_01, r5_04, story_01, fig_export.R; eye-eQTL LD (commands in README step 8) |
| SMR | 1.3.1 (Windows x86_64) | `data/raw/tools/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe` | https://yanglab.westlake.edu.cn/software/smr/ | 38_smr_heidi_druggable.R, 48a_fetch_tissue_eqtl.R, 74_adamts9.R |
| GCTA (COJO) | 1.94.1 (Windows x86_64) | `data/raw/tools/gcta-1.94.1-win-x86_64/…/gcta*.exe` (found by pattern) | https://github.com/jianyangqt/gcta/releases/tag/v1.94.1 | 49_cojo_susie.R |
| micromamba | latest (2026) | `data/raw/tools/micromamba.exe` | https://github.com/mamba-org/micromamba-releases | creates the single-cell Python environment (see python_requirements.txt) |
| R | 4.5.2 | on PATH or `RSCRIPT` | https://cran.r-project.org | all `.R` scripts |
| Python | 3.14.5 (general) and 3.11.16 (single-cell env) | on PATH | https://www.python.org ; conda-forge | all `.py` scripts |

Remote tabix queries (eQTL Catalogue, FinnGen, Pan-UKB, GWAS Catalog harmonised files) are done with
`Rsamtools` (R) or HTTP range requests (Python); no separate `tabix`/`htslib` binary is required.

## Reference data (not tools, but needed by the tools)

| Resource | Expected location | Source |
|---|---|---|
| 1000 Genomes phase 3 European reference (PLINK binary, n = 503) | `data/raw/ref/g1000_eur/g1000_eur.{bed,bim,fam}` | MAGMA auxiliary files, https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe/download (`g1000_eur.zip`) |
| 1000 Genomes EUR/FIN/EAS pairwise LD for the key 10q26 variants | `results/revision/10q26/ld_ensembl_EUR_FIN_EAS.tsv` (small input table) | Ensembl REST API LD endpoint (https://rest.ensembl.org/ld/human/pairwise/…; populations 1000GENOMES:phase_3:EUR/FIN/EAS) |
| LD scores (European, HapMap3) | `data/raw/ref/eur_w_ld_chr/`, `w_hm3.snplist` | https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/eur_w_ld_chr.tar.bz2 (mirror used: https://ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr/) |
| LAVA LD-block locus file | `data/raw/ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile` | https://github.com/josefin-werme/LAVA (support_data) |
| UK Biobank in-sample LD matrices (Weissbrod et al. 2020, British ancestry, n ≈ 337k) | downloaded on the fly by `revision3/panukb_01_fetch.py`, `revision4/rev4_C_fetch.py` | https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/ (`chr10_123000001_126000001.npz/.gz`, chr3 block) |

## Notes for non-Windows systems

- `38b_ld_subset.R` calls `plink.exe` with an explicit path; the other scripts use the `PLINK` variable.
- `R` parallelism uses PSOCK clusters (`parallel::makeCluster`), which works on all platforms. The scripts were run with
  up to 30 workers on a 16-core/32-thread CPU with 128 GB RAM; reduce the worker counts in `makeCluster(...)`
  calls on smaller machines.
