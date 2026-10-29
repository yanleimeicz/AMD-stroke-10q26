# Revision C, step 0: (i) re-derive IAMDGC with P values below 1e-300 kept (log-scale z), (ii) MEGASTROKE EUR outcomes.
# Original 30_prepare_gwas.R parsed P with as.numeric(), which underflows to 0 for P < ~1e-308 and dropped those rows
# (202 variants, including the 10q26 and CFH leads). Here z is recovered from the mantissa/exponent string on the log scale.
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages(library(data.table))
setDTthreads(30)
RAW <- file.path(PROJ_ROOT, "data/raw/gwas"); OUT <- file.path(PROJ_ROOT, "derived/revision"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

log10p_from_string <- function(s) {          # "1.2E-734" -> -733.92 ; "0.003" -> log10(0.003)
  s <- toupper(trimws(s)); m <- regmatches(s, regexec("^([0-9.]+)E([-+]?[0-9]+)$", s))
  out <- suppressWarnings(log10(as.numeric(s)))
  hit <- lengths(m) == 3
  out[hit] <- log10(as.numeric(vapply(m[hit], `[`, "", 2))) + as.numeric(vapply(m[hit], `[`, "", 3))
  out
}

# ---- IAMDGC (fixed)
d <- fread(file.path(RAW, "IAMDGC_GCST003219.h.tsv.gz"), colClasses = list(character = "p_value"),
           select = c("variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele", "p_value", "overall", "hm_chrom", "hm_pos"))
d[, l10p := log10p_from_string(p_value)]
d <- d[overall %in% c("+", "-") & is.finite(l10p)]
frq <- fread(file.path(PROJ_ROOT, "derived/g1000_eur.frq"))[, .(SNP, A1, A2, MAF)]
d <- merge(d, frq, by.x = "variant_id", by.y = "SNP")
comp <- c(A = "T", T = "A", C = "G", G = "C"); ea <- toupper(d$effect_allele)
d[, f := fifelse(ea == A1, MAF, fifelse(ea == A2, 1 - MAF, fifelse(comp[ea] == A1, MAF, fifelse(comp[ea] == A2, 1 - MAF, NA_real_))))]
d <- d[!is.na(f) & f > 0.005 & f < 0.995]
N <- 16144 + 17832; s <- 16144 / N
lp <- d$l10p * log(10) - log(2)                                    # log of P/2
z <- ifelse(d$overall == "+", 1, -1) * qnorm(lp, lower.tail = FALSE, log.p = TRUE)
se <- 1 / sqrt(2 * d$f * (1 - d$f) * N * s * (1 - s))
o <- data.table(SNP = d$variant_id, chr = as.character(fifelse(is.na(d$hm_chrom), as.character(d$chromosome), as.character(d$hm_chrom))),
                pos = as.numeric(d$hm_pos), ea = toupper(d$effect_allele), oa = toupper(d$other_allele),
                beta = z * se, se = se, p = 10^pmax(d$l10p, -300), log10p = d$l10p, eaf = d$f)
o <- o[grepl("^rs[0-9]+$", SNP) & nchar(ea) == 1 & nchar(oa) == 1][order(log10p)][!duplicated(SNP)]
o[, `:=`(ncase = 16144, ncontrol = 17832, n = N)]
saveRDS(o, file.path(OUT, "AMD_IAMDGC_fix.rds"), compress = FALSE)
message(sprintf("IAMDGC_fix: %d SNPs; %d with P < 1e-300 (previously dropped)", nrow(o), sum(o$log10p < -300)))
print(o[SNP %in% c("rs3750848", "rs2672592", "rs10490924", "rs3750846", "rs570618", "rs3184504", "rs4611812")])

# ---- MEGASTROKE EUR (Malik et al. 2018; GWAS Catalog harmonised files; no FinnGen / FINRISK)
std_gwascat <- function(f) {
  x <- fread(f)
  pick <- function(c) { y <- intersect(c, names(x)); if (length(y)) y[1] else NA_character_ }
  m <- list(SNP = pick(c("hm_rsid", "rsid", "variant_id")), chr = pick(c("hm_chrom", "chromosome")), pos = pick(c("hm_pos", "base_pair_location")),
            ea = pick(c("hm_effect_allele", "effect_allele")), oa = pick(c("hm_other_allele", "other_allele")),
            beta = pick(c("hm_beta", "beta")), or = pick(c("hm_odds_ratio", "odds_ratio")), se = pick(c("standard_error", "se")),
            p = pick(c("p_value", "pvalue")), eaf = pick(c("hm_effect_allele_frequency", "effect_allele_frequency")))
  o <- data.table(SNP = x[[m$SNP]], chr = as.character(x[[m$chr]]), pos = as.numeric(x[[m$pos]]), ea = toupper(x[[m$ea]]), oa = toupper(x[[m$oa]]),
                  beta = suppressWarnings(as.numeric(x[[m$beta]])), se = suppressWarnings(as.numeric(x[[m$se]])), p = suppressWarnings(as.numeric(x[[m$p]])),
                  eaf = if (!is.na(m$eaf)) suppressWarnings(as.numeric(x[[m$eaf]])) else NA_real_)
  if (!is.na(m$or)) { orv <- suppressWarnings(as.numeric(x[[m$or]])); o[is.na(beta) & !is.na(orv), beta := log(orv[is.na(o$beta) & !is.na(orv)])] }
  o[!is.na(beta) & !is.na(se) & se > 0 & !is.na(p) & grepl("^rs[0-9]+$", SNP) & nchar(ea) == 1 & nchar(oa) == 1][order(p)][!duplicated(SNP)]
}
MEGA <- list(MEGA_AS = list("MEGASTROKE_AS_EUR_GCST006906.h.tsv.gz", 40585, 406111), MEGA_IS = list("MEGASTROKE_IS_EUR_GCST006908.h.tsv.gz", 34217, 406111),
             MEGA_LAS = list("MEGASTROKE_LAS_EUR_GCST006907.h.tsv.gz", 4373, 297290), MEGA_SVS = list("MEGASTROKE_SVS_EUR_GCST006909.h.tsv.gz", 5386, 343560),
             MEGA_CES = list("MEGASTROKE_CES_EUR_GCST006910.h.tsv.gz", 7193, 355468))
for (k in names(MEGA)) {
  f <- file.path(RAW, MEGA[[k]][[1]]); if (!file.exists(f)) { message("missing ", f); next }
  o <- std_gwascat(f); o[, `:=`(ncase = MEGA[[k]][[2]], ncontrol = MEGA[[k]][[3]], n = MEGA[[k]][[2]] + MEGA[[k]][[3]])]
  saveRDS(o, file.path(OUT, paste0(k, ".rds")), compress = FALSE)
  message(sprintf("%s: %d SNPs, %d at P < 5e-8", k, nrow(o), sum(o$p < 5e-8)))
}
message("prep done")
