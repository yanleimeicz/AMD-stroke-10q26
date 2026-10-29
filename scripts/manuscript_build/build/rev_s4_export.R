# Export the 12q24 (SH2B3) regional data used for Figure S4 (same selection as scripts/49_cojo_susie.R plot_region)
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressMessages(library(data.table))
D <- file.path(PROJ_ROOT, "derived/tissue_eqtl"); FD <- file.path(PROJ_ROOT, "submission/revision/build/figdata")
rsnps <- fread(file.path(D, "region_snps.tsv"))[!is.na(pos37)]
gw <- readRDS(file.path(D, "gwas_regions.rds")); eq <- readRDS(file.path(D, "eqtl_harmonised.rds"))
pos38 <- rsnps[, .(SNP, pos38)]
ld <- fread(file.path(D, "cojo", "ld_rs3184504.ld"))[, .(SNP = SNP_B, r2 = R2)]
panel <- function(d, lab) { d <- merge(d[, .(SNP, p)], pos38, by = "SNP"); d[, panel := lab]; d }
reg <- "12q24_SH2B3"
P <- rbind(rbindlist(lapply(c("AMD_IAMDGC", "DRYAMD_FG", "IS", "AS"), function(k) panel(gw[trait == k & region == reg], k))),
           rbindlist(lapply(unique(eq[region == reg & gene == "SH2B3"]$dataset), function(ds) {
             e <- eq[region == reg & gene == "SH2B3" & dataset == ds]; if (min(e$p) > 1e-4) return(NULL); panel(e, paste("SH2B3", ds)) })))
P <- P[pos38 >= 111.1e6 & pos38 <= 111.8e6]
P <- merge(P, ld, by = "SNP", all.x = TRUE); P[SNP == "rs3184504", r2 := 1]
fwrite(P, file.path(FD, "s4_region12q24.tsv"), sep = "\t")
print(P[, .N, by = panel])
