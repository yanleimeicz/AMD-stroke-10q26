# Post-process GWAS-GWAS colocalisation: independent shared loci, gene annotation of the top shared SNP,
# allele-aligned effects of the shared SNP across all traits, and regional association plots (LD from 1000G EUR).
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
setDTthreads(16)
DER <- file.path(PROJ_ROOT, "derived/gwas"); OUT <- file.path(PROJ_ROOT, "results/coloc"); dir.create(file.path(OUT, "regional"), showWarnings = FALSE)
PLINK <- file.path(PROJ_ROOT, "data/raw/tools/plink/plink.exe"); BFILE <- file.path(PROJ_ROOT, "data/raw/ref/g1000_eur/g1000_eur")
cl <- fread(file.path(OUT, "coloc_gwas_gwas.tsv"))[PP.H4 > 0.5]
cl[, region := paste0(chr, ":", floor(pos38 / 2e6))]
# one row per region x (AMD trait, stroke trait) pair with max PP.H4
cl[, amd := fifelse(source %like% "AMD", source, partner)][, stk := fifelse(source %like% "AMD", partner, source)]
loci <- cl[order(-PP.H4)][!duplicated(paste(region, amd, stk))]

# nearest genes of the top shared SNP (FinnGen annotation)
ng <- fread(file.path(PROJ_ROOT, "data/raw/gwas/finngen_R12_H7_AMD.gz"), select = c("rsids", "#chrom", "pos", "nearest_genes"))
setnames(ng, c("rsids", "chr", "pos38", "nearest_genes")); ng[, SNP := sub(",.*", "", rsids)]
loci <- merge(loci[, !c("nearest_genes")], ng[, .(top_snp_H4 = SNP, top_pos38 = pos38, top_nearest = nearest_genes)], by = "top_snp_H4", all.x = TRUE)

# allele-aligned effects of each top SNP in every trait (effect allele = AMD_FG effect allele)
traits <- c("AMD_FG", "AMD_IAMDGC", "DRYAMD_FG", "WETAMD_FG", "AS", "IS", "LAS", "CES", "SVS", "WMH")
G <- lapply(setNames(traits, traits), function(k) readRDS(file.path(DER, paste0(k, ".rds")))[SNP %in% unique(loci$top_snp_H4)])
ref <- G$AMD_FG[, .(SNP, ref_ea = ea, ref_oa = oa)]
eff <- rbindlist(lapply(traits, function(k) {
  m <- merge(ref, G[[k]][, .(SNP, ea, oa, beta, se, p)], by = "SNP")
  m[, b := fifelse(ea == ref_ea, beta, fifelse(ea == ref_oa, -beta, NA_real_))]
  m[, .(SNP, trait = k, ea = ref_ea, b, se, p)]
}))
effw <- dcast(eff, SNP + ea ~ trait, value.var = c("b", "p"))
fwrite(eff, file.path(OUT, "shared_snp_effects_long.tsv"), sep = "\t")
loci <- merge(loci, effw, by.x = "top_snp_H4", by.y = "SNP", all.x = TRUE)
loci[, direction := {
  ba <- mapply(function(s, a) eff[SNP == s & trait == a, b][1], top_snp_H4, amd)
  bs <- mapply(function(s, t) eff[SNP == s & trait == t, b][1], top_snp_H4, stk)
  fifelse(sign(ba) == sign(bs), "concordant (same risk allele)", "discordant (opposite)")
}]
setorder(loci, -PP.H4)
fwrite(loci, file.path(OUT, "coloc_shared_loci_annotated.tsv"), sep = "\t")
print(loci[, .(region, amd, stk, PP.H4 = round(PP.H4, 3), top_snp_H4, top_nearest, direction)])

# regional plots for each independent region (best pair)
reg_best <- loci[!duplicated(region)]
for (i in seq_len(nrow(reg_best))) {
  r <- reg_best[i]; snp <- r$top_snp_H4; ch <- as.character(r$chr); c0 <- r$top_pos38
  win <- G$AMD_FG[0]
  d1 <- readRDS(file.path(DER, paste0(r$amd, ".rds"))); d2 <- readRDS(file.path(DER, paste0(r$stk, ".rds")))
  pos <- readRDS(file.path(DER, "AMD_FG.rds"))[chr == ch & pos >= c0 - 5e5 & pos <= c0 + 5e5, .(SNP, pos)]
  m <- merge(pos, d1[, .(SNP, p1 = p)], by = "SNP"); m <- merge(m, d2[, .(SNP, p2 = p)], by = "SNP")
  rm(d1, d2); gc()
  tmp <- tempfile(tmpdir = file.path(PROJ_ROOT, "tmp")); writeLines(m$SNP, paste0(tmp, ".snps"))
  system2(PLINK, c("--bfile", BFILE, "--extract", paste0(tmp, ".snps"), "--r2", "--ld-snp", snp, "--ld-window-kb", "1000",
                   "--ld-window", "99999", "--ld-window-r2", "0", "--out", tmp), stdout = FALSE, stderr = FALSE)
  ld <- if (file.exists(paste0(tmp, ".ld"))) fread(paste0(tmp, ".ld"))[, .(SNP = SNP_B, r2 = R2)] else data.table(SNP = character(), r2 = numeric())
  m <- merge(m, ld, by = "SNP", all.x = TRUE); m[is.na(r2), r2 := 0]
  long <- rbind(m[, .(SNP, pos, logp = -log10(p1), r2, trait = r$amd)], m[, .(SNP, pos, logp = -log10(p2), r2, trait = r$stk)])
  long[, trait := factor(trait, c(r$amd, r$stk))]
  p <- ggplot(long, aes(pos / 1e6, logp, colour = r2)) + geom_point(size = .8) +
    geom_point(data = long[SNP == snp], shape = 18, size = 3.5, colour = "purple") +
    scale_colour_gradientn(colours = c("navy", "lightblue", "green", "orange", "red"), limits = c(0, 1), name = expression(r^2)) +
    facet_wrap(~trait, ncol = 1, scales = "free_y") + theme_bw(base_size = 9) +
    labs(x = paste0("chr", ch, " position (Mb, GRCh38)"), y = expression(-log[10](P)),
         title = sprintf("%s  %s ~ %s  PP.H4 = %.2f  (%s)", snp, r$amd, r$stk, r$PP.H4, r$top_nearest))
  ggsave(file.path(OUT, "regional", sprintf("regional_chr%s_%s_%s_%s.pdf", ch, snp, r$amd, r$stk)), p, width = 6.5, height = 5)
}
message("coloc post done")
