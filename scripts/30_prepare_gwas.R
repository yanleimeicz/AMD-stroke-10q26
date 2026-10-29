# Standardise GWAS summary statistics to a common format and save as RDS.
# Columns: SNP chr pos ea oa eaf beta se p n ncase ncontrol  (pos = GRCh38 for FinnGen / GWAS-Catalog harmonised; GRCh37 for BIG40)
# Usage: Rscript 30_prepare_gwas.R [trait ...]   (default: all available)
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
setDTthreads(16)
RAW <- file.path(PROJ_ROOT, "data/raw/gwas"); OUT <- file.path(PROJ_ROOT, "derived/gwas"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TRAITS <- list(
  AMD_IAMDGC = list(file = "IAMDGC_GCST003219.h.tsv.gz", fmt = "iamdgc", ncase = 16144, ncontrol = 17832, label = "Advanced AMD (IAMDGC)"),
  AMD_FG     = list(file = "finngen_R12_H7_AMD.gz",  fmt = "finngen", ncase = 12495, ncontrol = 461686, label = "AMD (FinnGen R12)"),
  DRYAMD_FG  = list(file = "finngen_R12_DRY_AMD.gz", fmt = "finngen", ncase = 8570,  ncontrol = 329258, label = "Dry AMD (FinnGen R12)"),
  WETAMD_FG  = list(file = "finngen_R12_WET_AMD.gz", fmt = "finngen", ncase = 6699,  ncontrol = 331070, label = "Wet AMD (FinnGen R12)"),
  AS  = list(file = "GIGASTROKE_GCST90104539.h.tsv.gz", fmt = "gwascat", ncase = 73652, ncontrol = 1234808, label = "Any stroke (GIGASTROKE)"),
  IS  = list(file = "GIGASTROKE_GCST90104540.h.tsv.gz", fmt = "gwascat", ncase = 62100, ncontrol = 1234808, label = "Ischemic stroke (GIGASTROKE)"),
  CES = list(file = "GIGASTROKE_GCST90104541.h.tsv.gz", fmt = "gwascat", ncase = 10804, ncontrol = 1234808, label = "Cardioembolic stroke (GIGASTROKE)"),
  LAS = list(file = "GIGASTROKE_GCST90104542.h.tsv.gz", fmt = "gwascat", ncase = 6399,  ncontrol = 1234808, label = "Large artery stroke (GIGASTROKE)"),
  SVS = list(file = "GIGASTROKE_GCST90104543.h.tsv.gz", fmt = "gwascat", ncase = 6811,  ncontrol = 1234808, label = "Small vessel stroke (GIGASTROKE)"),
  WMH = list(file = "BIG40_1437.txt.gz", fmt = "big40", n = 33224, label = "WMH volume (UKB BIG40)")
)
saveRDS(TRAITS, file.path(OUT, "trait_meta.rds"))

pickcol <- function(d, cands) { x <- intersect(cands, names(d)); if (length(x)) x[1] else NA_character_ }

std_finngen <- function(f) {
  d <- fread(f, select = c("#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"))
  setnames(d, c("chr", "pos", "oa", "ea", "SNP", "p", "beta", "se", "eaf"))
  d[, SNP := sub(",.*", "", SNP)]
  d[, chr := as.character(chr)][chr == "23", chr := "X"]
  d
}
std_gwascat <- function(f) {
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
  # if beta missing but OR present (some harmonised rows)
  if (!is.na(m$or)) o[is.na(beta), beta := log(as.numeric(d[[m$or]]))[is.na(o$beta)]]
  o
}
std_iamdgc <- function(f) {
  # Fritsche 2016 provides only P and overall direction ("+"/"-") relative to effect_allele.
  # Recover z, then approximate log-OR: se = 1 / sqrt(2 f (1 - f) N s (1 - s)), beta = z * se (f from 1000G EUR).
  d <- fread(f, select = c("variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele", "p_value", "overall", "hm_chrom", "hm_pos"))
  d[, p_value := suppressWarnings(as.numeric(p_value))]
  d <- d[overall %in% c("+", "-") & !is.na(p_value) & p_value > 0]
  frq <- fread(file.path(PROJ_ROOT, "derived/g1000_eur.frq"))[, .(SNP, A1, A2, MAF)]
  d <- merge(d, frq, by.x = "variant_id", by.y = "SNP")
  comp <- c(A = "T", T = "A", C = "G", G = "C")
  ea <- toupper(d$effect_allele)
  d[, f := fifelse(ea == A1, MAF, fifelse(ea == A2, 1 - MAF, fifelse(comp[ea] == A1, MAF, fifelse(comp[ea] == A2, 1 - MAF, NA_real_))))]
  d <- d[!is.na(f) & f > 0.005 & f < 0.995]
  N <- 16144 + 17832; s <- 16144 / N
  z <- ifelse(d$overall == "+", 1, -1) * qnorm(pmax(d$p_value, 1e-300) / 2, lower.tail = FALSE)
  se <- 1 / sqrt(2 * d$f * (1 - d$f) * N * s * (1 - s))
  data.table(SNP = d$variant_id, chr = as.character(fifelse(is.na(d$hm_chrom), as.character(d$chromosome), as.character(d$hm_chrom))),
             pos = as.numeric(fifelse(is.na(d$hm_pos), NA_real_, as.numeric(d$hm_pos))),
             ea = toupper(d$effect_allele), oa = toupper(d$other_allele), beta = z * se, se = se, p = d$p_value, eaf = d$f)
}
std_big40 <- function(f) {
  d <- fread(f)
  pcol <- grep("pval", names(d), value = TRUE)[1]
  o <- data.table(SNP = d$rsid, chr = as.character(d$chr), pos = d$pos, ea = d$a2, oa = d$a1,   # BIG40: effect allele = a2
                  beta = d$beta, se = d$se,
                  p = if (grepl("log", pcol)) 10^(-d[[pcol]]) else d[[pcol]], eaf = NA_real_)
  o[chr == "0X" | chr == "23", chr := "X"][, chr := sub("^0", "", chr)]
  o
}

args <- commandArgs(TRUE); todo <- if (length(args)) args else names(TRAITS)
for (k in todo) {
  tr <- TRAITS[[k]]; f <- file.path(RAW, tr$file); o <- file.path(OUT, paste0(k, ".rds"))
  if (!file.exists(f)) { message("missing ", k); next }
  if (file.exists(o) && file.mtime(o) > file.mtime(f)) { message("up-to-date ", k); next }
  t0 <- Sys.time(); message("== ", k)
  d <- switch(tr$fmt, finngen = std_finngen(f), gwascat = std_gwascat(f), iamdgc = std_iamdgc(f), big40 = std_big40(f))
  d <- d[!is.na(beta) & !is.na(se) & se > 0 & !is.na(p) & grepl("^rs[0-9]+$", SNP) & nchar(ea) == 1 & nchar(oa) == 1]
  d <- d[order(p)][!duplicated(SNP)]
  if (!is.null(tr$n)) d[, n := tr$n] else { d[, ncase := tr$ncase]; d[, ncontrol := tr$ncontrol]; d[, n := tr$ncase + tr$ncontrol] }
  saveRDS(d, o, compress = FALSE)
  message(sprintf("   %s: %d SNPs, %d with p<5e-8, %.0fs", k, nrow(d), sum(d$p < 5e-8), as.numeric(Sys.time() - t0, units = "secs")))
  print(head(d[order(p)], 3))
}
