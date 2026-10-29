# Standardise additional cerebral small vessel disease (CSVD) GWAS to the format of 30_prepare_gwas.R and save as RDS.
# Columns: SNP chr pos ea oa eaf beta se p n [ncase ncontrol]  (pos = GRCh38 for GWAS-Catalog harmonised files)
#   LACUNAR  Traylor 2021 Lancet Neurol, European lacunar stroke (GCST90014122; 6,030 cases / 219,389 controls per README)
#   EPVS_WM / EPVS_BG / EPVS_HIP  Duperron 2023 Nat Med, extensive perivascular space burden, European (GCST90244151/53/55;
#            beta/se from the secondary IVW meta-analysis, P from the primary sample-size meta-analysis)
#   WM_MD / WM_FA  Persyn 2020 Nat Commun, UK Biobank white matter microstructure (first PC of mean diffusivity / FA)
# Usage: Rscript 45_csvd_gwas_prep.R [trait ...]
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table) })
setDTthreads(30)
RAW <- file.path(PROJ_ROOT, "data/raw/gwas"); OUT <- file.path(PROJ_ROOT, "derived/gwas")

TRAITS <- list(
  LACUNAR  = list(file = "LACUNAR_Traylor2021_GCST90014122.h.tsv.gz", fmt = "gwascat", ncase = 6030, ncontrol = 219389, label = "Lacunar stroke (Traylor 2021)", gcst = "GCST90014122"),
  EPVS_WM  = list(file = "EPVS_WM_Duperron2023_GCST90244151.h.tsv.gz", fmt = "gwascat", ncase = 9317, ncontrol = 29281, label = "Extensive WM-PVS (Duperron 2023)", gcst = "GCST90244151"),
  EPVS_BG  = list(file = "EPVS_BG_Duperron2023_GCST90244153.h.tsv.gz", fmt = "gwascat", ncase = 8950, ncontrol = 29953, label = "Extensive BG-PVS (Duperron 2023)", gcst = "GCST90244153"),
  EPVS_HIP = list(file = "EPVS_HIP_Duperron2023_GCST90244155.h.tsv.gz", fmt = "gwascat", ncase = 9163, ncontrol = 29708, label = "Extensive HIP-PVS (Duperron 2023)", gcst = "GCST90244155"),
  WM_MD    = list(file = "WM_MD_Persyn2020_GCST010103.zip", fmt = "persyn", n = 17467, label = "White matter MD, PC1 (Persyn 2020)", gcst = "GCST010103"),
  WM_FA    = list(file = "WM_FA_Persyn2020_GCST010102.zip", fmt = "persyn", n = 17663, label = "White matter FA, PC1 (Persyn 2020)", gcst = "GCST010102")
)
saveRDS(TRAITS, file.path(OUT, "csvd_trait_meta.rds"))

pickcol <- function(d, cands) { x <- intersect(cands, names(d)); if (length(x)) x[1] else NA_character_ }

std_gwascat <- function(f) {   # same logic as 30_prepare_gwas.R (harmonised columns preferred)
  d <- fread(f)
  m <- list(SNP = pickcol(d, c("hm_rsid", "rsid", "variant_id")),
            chr = pickcol(d, c("hm_chrom", "chromosome")), pos = pickcol(d, c("hm_pos", "base_pair_location")),
            ea = pickcol(d, c("hm_effect_allele", "effect_allele")), oa = pickcol(d, c("hm_other_allele", "other_allele")),
            beta = pickcol(d, c("hm_beta", "beta")), or = pickcol(d, c("hm_odds_ratio", "odds_ratio")),
            se = pickcol(d, c("standard_error", "se")), p = pickcol(d, c("p_value", "pvalue")),
            eaf = pickcol(d, c("hm_effect_allele_frequency", "effect_allele_frequency")))
  message("  column map: ", paste(names(m), unlist(m), sep = "=", collapse = " "))
  o <- data.table(SNP = d[[m$SNP]], chr = as.character(d[[m$chr]]), pos = as.numeric(d[[m$pos]]),
                  ea = toupper(d[[m$ea]]), oa = toupper(d[[m$oa]]),
                  beta = if (!is.na(m$beta)) as.numeric(d[[m$beta]]) else log(as.numeric(d[[m$or]])),
                  se = as.numeric(d[[m$se]]), p = as.numeric(d[[m$p]]),
                  eaf = if (!is.na(m$eaf)) as.numeric(d[[m$eaf]]) else NA_real_)
  if (!is.na(m$or)) o[is.na(beta), beta := log(as.numeric(d[[m$or]]))[is.na(o$beta)]]
  if ("N" %in% names(d)) o[, n_snp := as.numeric(d$N)]
  o
}
std_persyn <- function(f) {   # zip with one whitespace-delimited table; column names checked at run time
  tmp <- file.path(file.path(PROJ_ROOT, "tmp/csvd"), sub("\\.zip$", "", basename(f))); dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  unzip(f, exdir = tmp); ff <- list.files(tmp, full.names = TRUE, recursive = TRUE)
  ff <- ff[file.info(ff)$size == max(file.info(ff)$size)][1]
  d <- fread(ff); message("  columns: ", paste(names(d), collapse = " "))
  m <- list(SNP = pickcol(d, c("SNP", "rsid", "RSID", "ID", "MarkerName")), chr = pickcol(d, c("CHR", "chr", "Chr", "CHROM")),
            pos = pickcol(d, c("BP", "POS", "pos", "bp", "Position")),
            ea = pickcol(d, c("A1", "ALLELE1", "Allele1", "EA", "effect_allele")), oa = pickcol(d, c("AX", "A2", "ALLELE0", "Allele2", "OA", "other_allele", "ALLELE2")),
            beta = pickcol(d, c("BETA", "beta", "Effect", "b")), se = pickcol(d, c("SE", "se", "StdErr")),
            p = pickcol(d, c("P", "p", "P_BOLT_LMM_INF", "P_BOLT_LMM", "pval", "P-value")),
            eaf = pickcol(d, c("A1_FREQ", "A1FREQ", "EAF", "Freq1", "MAF", "freq")))
  message("  column map: ", paste(names(m), unlist(m), sep = "=", collapse = " "))
  if ("INFO" %in% names(d)) d <- d[INFO >= 0.8]                          # imputation quality (positions GRCh37)
  o <- data.table(SNP = d[[m$SNP]], chr = as.character(d[[m$chr]]), pos = as.numeric(d[[m$pos]]), ea = toupper(d[[m$ea]]), oa = toupper(d[[m$oa]]),
                  beta = as.numeric(d[[m$beta]]), se = as.numeric(d[[m$se]]), p = as.numeric(d[[m$p]]),
                  eaf = if (!is.na(m$eaf)) as.numeric(d[[m$eaf]]) else NA_real_)
  o[is.na(eaf) | (eaf >= 0.01 & eaf <= 0.99)]
}

args <- commandArgs(TRUE); todo <- if (length(args)) args else names(TRAITS)
for (k in todo) {
  tr <- TRAITS[[k]]; f <- file.path(RAW, tr$file); o <- file.path(OUT, paste0(k, ".rds"))
  if (!file.exists(f)) { message("missing ", k); next }
  if (file.exists(o) && file.mtime(o) > file.mtime(f)) { message("up-to-date ", k); next }
  t0 <- Sys.time(); message("== ", k)
  d <- switch(tr$fmt, gwascat = std_gwascat(f), persyn = std_persyn(f))
  d <- d[!is.na(beta) & !is.na(se) & se > 0 & !is.na(p) & grepl("^rs[0-9]+$", SNP) & nchar(ea) == 1 & nchar(oa) == 1]
  if ("n_snp" %in% names(d)) {   # drop variants with low per-SNP sample size (< 50% of the maximum)
    d <- d[is.na(n_snp) | n_snp >= 0.5 * max(n_snp, na.rm = TRUE)]; d[, n_snp := NULL]
  }
  d <- d[order(p)][!duplicated(SNP)]
  if (!is.null(tr$n)) d[, n := tr$n] else { d[, ncase := tr$ncase]; d[, ncontrol := tr$ncontrol]; d[, n := tr$ncase + tr$ncontrol] }
  saveRDS(d, o, compress = FALSE)
  message(sprintf("   %s: %d SNPs, %d with p<5e-8, %.0fs", k, nrow(d), sum(d$p < 5e-8), as.numeric(Sys.time() - t0, units = "secs")))
  print(head(d[order(p)], 3))
}
message("CSVD prep done")
