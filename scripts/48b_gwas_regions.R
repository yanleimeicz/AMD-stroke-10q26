# Extract GWAS summary statistics for the 7 shared-locus regions (SNP set defined by rsID from FinnGen GRCh38 windows,
# so GRCh37 traits are matched by rsID). Output -> derived/tissue_eqtl/gwas_regions.rds, region_snps_bim.tsv
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
OUT <- file.path(PROJ_ROOT, "derived/tissue_eqtl"); G <- file.path(PROJ_ROOT, "derived/gwas")
REG <- fread(file.path(OUT, "regions_hg38.tsv"), colClasses = c(chr = "character"))
TRAITS <- c("AMD_FG", "AMD_IAMDGC", "DRYAMD_FG", "WETAMD_FG", "AS", "IS", "LAS", "SVS", "WMH")

fg <- as.data.table(readRDS(file.path(G, "AMD_FG.rds")))
snps <- rbindlist(lapply(seq_len(nrow(REG)), function(j)
  fg[chr == REG$chr[j] & pos >= REG$start[j] & pos <= REG$end[j] & grepl("^rs", SNP), .(SNP, region = REG$region[j], pos38 = pos)]))
snps <- unique(snps, by = "SNP"); rm(fg); gc()
bim <- fread(file.path(PROJ_ROOT, "data/raw/ref/g1000_eur/g1000_eur.bim"), header = FALSE, select = c(1, 2, 4, 5, 6),
             col.names = c("chr", "SNP", "pos37", "A1", "A2"))[SNP %in% snps$SNP]
snps <- merge(snps, bim[, .(SNP, pos37, A1, A2)], by = "SNP", all.x = TRUE)
fwrite(snps, file.path(OUT, "region_snps.tsv"), sep = "\t")

cl <- makeCluster(5); clusterExport(cl, c("G", "snps")); clusterEvalQ(cl, { .libPaths(c(file.path(Sys.getenv("COMORBIDITY_ROOT"), "Rlib"), .libPaths())); library(data.table) })
res <- parLapply(cl, TRAITS, function(k) {
  d <- as.data.table(readRDS(file.path(G, paste0(k, ".rds"))))[SNP %in% snps$SNP]
  d[, trait := k]; d
})
stopCluster(cl)
gw <- rbindlist(res, fill = TRUE)
gw <- merge(gw, snps[, .(SNP, region, pos38, pos37)], by = "SNP")
saveRDS(gw, file.path(OUT, "gwas_regions.rds"))
print(gw[, .N, by = .(trait, region)][order(region, trait)])
