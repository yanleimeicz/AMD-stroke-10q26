# Lookup of the ABO GWAS-GWAS colocalization SNP rs947073006 across traits (for Fig. 6a)
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressMessages({ library(data.table); library(parallel) })
P <- PROJ_ROOT; OUT <- file.path(P, "submission/build/figdata")
TR <- c("AMD_FG", "DRYAMD_FG", "WETAMD_FG", "AS", "IS", "LAS", "CES", "SVS", "LACUNAR", "WMH", "WM_MD", "WM_FA", "EPVS_WM", "EPVS_BG", "EPVS_HIP")
cl <- makeCluster(8); clusterExport(cl, "P"); clusterEvalQ(cl, { .libPaths(c(file.path(Sys.getenv("COMORBIDITY_ROOT"), "Rlib"), .libPaths())); library(data.table) })
kv <- rbindlist(parLapplyLB(cl, TR, function(k) {
  d <- as.data.table(readRDS(file.path(P, "derived/gwas", paste0(k, ".rds"))))[SNP == "rs947073006"]
  d[, trait := k]; d[, .(SNP, trait, ea, oa, beta, se, p)]
}), fill = TRUE)
stopCluster(cl)
kv[, label := "9q34 rs947073006 (ABO)"]
print(kv)
fwrite(kv, file.path(OUT, "fig6_abo_rs947073006.tsv"), sep = "\t")
cat("ABO DONE\n")
