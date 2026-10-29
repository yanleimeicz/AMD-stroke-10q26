# Export plotting data for the Human Genomics figures -> submission/build/figdata/*.tsv
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
P <- PROJ_ROOT; OUT <- file.path(P, "submission/build/figdata"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
PLINK <- file.path(P, "data/raw/tools/plink/plink.exe"); BFILE <- file.path(P, "data/raw/ref/g1000_eur/g1000_eur")
TMP <- file.path(P, "tmp/figs"); dir.create(TMP, showWarnings = FALSE, recursive = TRUE)

# ---- 10q26 regional GWAS (AMD_FG, IS, SVS from gwas_regions.rds; LACUNAR from full file)
gw <- as.data.table(readRDS(file.path(P, "derived/tissue_eqtl/gwas_regions.rds")))
print(names(gw)); print(gw[, .N, by = .(trait, region)][region == "10q26_HTRA1"])
r10 <- gw[region == "10q26_HTRA1" & trait %in% c("AMD_FG", "IS", "SVS")]
lac <- as.data.table(readRDS(file.path(P, "derived/gwas/LACUNAR.rds")))[SNP %in% r10$SNP]
lac <- merge(lac, unique(r10[, .(SNP, region, pos38, pos37)]), by = "SNP"); lac[, trait := "LACUNAR"]
r10 <- rbindlist(list(r10, lac), fill = TRUE)
fwrite(r10[, .(SNP, trait, pos37, pos38, ea, oa, beta, se, p)], file.path(OUT, "region10q26_gwas.tsv"), sep = "\t")

# LD (r2) with rs3750848 and rs2672592 in the 10q26 window (GRCh37 1000G EUR)
ld_with <- function(lead, chr, from, to, tag) {
  o <- file.path(TMP, tag)
  system2(PLINK, c("--bfile", BFILE, "--chr", chr, "--from-bp", from, "--to-bp", to, "--r2", "--ld-snp", lead,
                   "--ld-window-kb", 2000, "--ld-window", 99999, "--ld-window-r2", 0, "--out", o), stdout = FALSE)
  d <- fread(paste0(o, ".ld")); d[, .(SNP = SNP_B, r2 = R2, lead = lead)]
}
p37 <- range(r10$pos37, na.rm = TRUE)
ld10 <- rbind(ld_with("rs3750848", 10, p37[1], p37[2], "ld10_a"), ld_with("rs2672592", 10, p37[1], p37[2], "ld10_b"))
fwrite(ld10, file.path(OUT, "region10q26_ld.tsv"), sep = "\t")

# ---- 3p14 regional (AMD_FG, IS, LAS, cIMT, aorta ADAMTS9 eQTL) + LD with rs4611812
g3 <- as.data.table(readRDS(file.path(P, "derived/adamts9/gwas_region.rds")))
print(names(g3)); print(unique(g3$trait))
g3 <- g3[trait %in% c("AMD_FG", "IS", "LAS")]
ci <- as.data.table(readRDS(file.path(P, "derived/adamts9/gwas_cimt.rds"))); print(names(ci)); print(head(ci, 2))
fwrite(g3[, .(SNP, trait, pos37, pos38, ea, oa, beta, se, p)], file.path(OUT, "region3p14_gwas.tsv"), sep = "\t")
fwrite(ci, file.path(OUT, "region3p14_cimt_raw.tsv"), sep = "\t")
eq <- fread(file.path(P, "derived/adamts9/eqtl/GTEx_artery_aorta.tsv.gz")); print(names(eq)); print(head(eq, 2))
fwrite(eq, file.path(OUT, "region3p14_aorta_eqtl_raw.tsv"), sep = "\t")
p37b <- range(g3$pos37, na.rm = TRUE)
fwrite(ld_with("rs4611812", 3, p37b[1], p37b[2], "ld3"), file.path(OUT, "region3p14_ld.tsv"), sep = "\t")

# ---- Fig 6 key-variant lookup across traits (full GWAS files)
abo <- fread(file.path(P, "results/tissue_eqtl/coloc_tissue_eqtl_gwas.tsv"))[region == "9q34_ABO" & dataset == "GTEx_Artery_Tibial" & trait == "LAS"]
print(abo); abo_snp <- abo$coloc_snp[1]
KEY <- data.table(SNP = c("rs3750848", "rs2672592", "rs4611812", "rs429358", "rs1269852", abo_snp, "rs3184504", "rs10737680"),
                  label = c("10q26 rs3750848 (AMD haplotype)", "10q26 rs2672592 (HTRA1, Malik 2024)", "3p14 rs4611812 (ADAMTS9)",
                            "19q13 rs429358 (APOE e4)", "6p21 rs1269852 (C4A eQTL)", paste0("9q34 ", abo_snp, " (ABO)"),
                            "12q24 rs3184504 (SH2B3 R262W)", "1q31 rs10737680 (CFH)"))
TR <- c("AMD_FG", "DRYAMD_FG", "WETAMD_FG", "AS", "IS", "LAS", "CES", "SVS", "LACUNAR", "WMH", "WM_MD", "WM_FA", "EPVS_WM", "EPVS_BG", "EPVS_HIP")
cl <- makeCluster(8); clusterExport(cl, c("P", "KEY"))
clusterEvalQ(cl, { .libPaths(c(file.path(Sys.getenv("COMORBIDITY_ROOT"), "Rlib"), .libPaths())); library(data.table) })
kv <- rbindlist(parLapplyLB(cl, TR, function(k) {
  d <- as.data.table(readRDS(file.path(P, "derived/gwas", paste0(k, ".rds"))))[SNP %in% KEY$SNP]
  d[, trait := k]; d[, .(SNP, trait, ea, oa, beta, se, p)]
}), fill = TRUE)
stopCluster(cl)
kv <- merge(kv, KEY, by = "SNP")
fwrite(kv, file.path(OUT, "fig6_key_variants_raw.tsv"), sep = "\t")
print(kv[, .N, by = trait])
cat("EXPORT DONE\n")
