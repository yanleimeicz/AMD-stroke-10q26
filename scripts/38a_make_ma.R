# GWAS summary statistics -> SMR .ma format (SNP A1 A2 freq b se p n; A1 = effect allele).
# freq = GWAS effect-allele frequency where available, otherwise 1000G EUR frequency of the effect allele.
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
DER <- file.path(PROJ_ROOT, "derived/gwas"); OUT <- file.path(PROJ_ROOT, "derived/smr_ma"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
frq <- fread(file.path(PROJ_ROOT, "derived/g1000_eur.frq"))[, .(SNP, fA1 = A1, fA2 = A2, MAF)]
OUTCOMES <- c("AMD_FG", "AMD_IAMDGC", "DRYAMD_FG", "WETAMD_FG", "AS", "IS", "LAS", "CES", "SVS", "WMH")
args <- commandArgs(TRUE); if (length(args)) OUTCOMES <- args

for (k in OUTCOMES) {
  o <- file.path(OUT, paste0(k, ".ma"))
  if (file.exists(o)) { message("skip ", k); next }
  d <- readRDS(file.path(DER, paste0(k, ".rds")))[, .(SNP, A1 = ea, A2 = oa, freq = eaf, b = beta, se, p, n)]
  d <- merge(d, frq, by = "SNP", all.x = TRUE)
  d[, f1kg := fifelse(A1 == fA1, MAF, fifelse(A1 == fA2, 1 - MAF, NA_real_))]
  n_gwas <- d[!is.na(freq), .N]
  d[is.na(freq) | freq <= 0 | freq >= 1, freq := f1kg]
  d <- d[!is.na(freq) & freq > 0 & freq < 1]
  d <- unique(d, by = "SNP")
  fwrite(d[, .(SNP, A1, A2, freq, b, se, p, n)], o, sep = " ")
  message(sprintf("%s: %d SNPs written (%d with GWAS freq)", k, nrow(d), n_gwas))
}
message("ma done")
