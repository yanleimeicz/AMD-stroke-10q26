# Export extra plotting data for the revised figures -> submission/revision/build/figdata/
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressMessages(library(data.table)); setDTthreads(16)
P <- PROJ_ROOT; OUT <- file.path(P, "submission/revision/build/figdata"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CM <- file.path(P, "data/raw/gwas/cardiometab")
rs <- fread(file.path(P, "derived/adamts9/region_snps.tsv"))          # SNP pos38 pos37 A1 A2

# T2D (Xue 2018) regional association at 3p14
x <- fread(file.path(CM, "T2D_Xue2018.h.tsv.gz"), select = c("hm_rsid", "hm_effect_allele", "hm_other_allele", "hm_beta", "standard_error", "p_value"))
setnames(x, c("SNP", "ea", "oa", "beta", "se", "p"))
x <- merge(x[SNP %in% rs$SNP], rs[, .(SNP, pos38, pos37)], by = "SNP")
x[, trait := "T2D"]
fwrite(x[, .(SNP, trait, pos37, pos38, ea, oa, beta, se, p)], file.path(OUT, "region3p14_t2d.tsv"), sep = "\t")
cat("T2D region rows:", nrow(x), "\n")

# rs4611812 effects (C allele) for DBP, SBP (Evangelou 2018, GRCh37 position 3:64699445) and WHRadjBMI (Pulit 2019)
alignC <- function(d, trait) {
  d <- d[1]; s <- if (toupper(d$ea) == "C") 1 else if (toupper(d$oa) == "C") -1 else NA_real_
  data.table(trait = trait, ea = d$ea, oa = d$oa, beta_C = s * d$beta, se = d$se, p = d$p)
}
res <- list()
for (bp in c("SBP", "DBP")) {
  b <- fread(file.path(CM, paste0(bp, "_Evangelou2018.txt.gz")), select = c("MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P"))
  setnames(b, c("marker", "ea", "oa", "beta", "se", "p"))
  b <- b[startsWith(marker, "3:64699445")]
  print(b); if (nrow(b)) res[[bp]] <- alignC(b, bp)
}
w <- fread(file.path(CM, "WHRadjBMI_Pulit2019.tbl.gz"), select = c("SNP", "Allele1", "Allele2", "Effect", "StdErr", "P-value"))
setnames(w, c("SNP", "ea", "oa", "beta", "se", "p")); w <- w[sub(":.*", "", SNP) == "rs4611812"]
print(w); if (nrow(w)) res[["WHRadjBMI"]] <- alignC(w, "WHRadjBMI")
res <- rbindlist(res); print(res)
fwrite(res, file.path(OUT, "rs4611812C_bp_whr.tsv"), sep = "\t")
cat("EXPORT DONE\n")
