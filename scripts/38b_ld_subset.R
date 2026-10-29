# 1000G EUR LD reference restricted to eQTLGen SNPs (SMR reads .bim slowly on Windows)
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
LD <- file.path(PROJ_ROOT, "derived/ld_eqtlgen/g1000_eur_eqtlgen")
writeLines(fread(file.path(PROJ_ROOT, "data/raw/eqtl/eqtlgen_besd.esi"), header = FALSE, select = 2)[[1]], paste0(LD, "_snps.txt"))
system2(file.path(PROJ_ROOT, "data/raw/tools/plink/plink.exe"), c("--bfile", file.path(PROJ_ROOT, "data/raw/ref/g1000_eur/g1000_eur"),
        "--extract", paste0(LD, "_snps.txt"), "--make-bed", "--threads", "30", "--out", LD))
