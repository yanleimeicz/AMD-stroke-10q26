# MHC class III: is the C4-region AMD / WMH association independent of the established CFB/C2 AMD signal?
# Two-SNP conditional models (z_cond(b|a) = (z_b - r z_a)/sqrt(1-r^2)), r from 1000G EUR (signed, A1-coded), r +/- 0.1 sensitivity.
# Inputs: derived/tissue_eqtl/gwas_regions.rds (region 6p21_C4). Output: results/revision/extra/c4_conditional.tsv, c4_leads.tsv
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
D <- file.path(PROJ_ROOT, "derived/tissue_eqtl"); OUT <- file.path(PROJ_ROOT, "results/revision/extra")
W <- file.path(PROJ_ROOT, "tmp/rev_c4"); dir.create(W, recursive = TRUE, showWarnings = FALSE)
PLINK <- file.path(PROJ_ROOT, "data/raw/tools/plink/plink.exe"); BFILE <- file.path(PROJ_ROOT, "data/raw/ref/g1000_eur/g1000_eur")
gw <- as.data.table(readRDS(file.path(D, "gwas_regions.rds")))
cat("columns:", paste(names(gw), collapse = ","), "\n"); print(gw[, .N, by = .(region, trait)][grepl("C4", region)])
g <- gw[grepl("C4", region)]
pcol <- if ("p" %in% names(g)) "p" else grep("^p", names(g), value = TRUE)[1]
g[, z := beta / se]
TR <- intersect(c("AMD_FG", "DRYAMD_FG", "WETAMD_FG", "AMD_IAMDGC", "WMH"), unique(g$trait))

# CFB/C2 (and adjacent NELFE/SKIV2L) window, GRCh38 chr6:31.89-31.97 Mb; candidate known AMD variants
cand <- c("rs641153", "rs9332739", "rs429608", "rs4151667", "rs116503776", "rs12614")
win <- g[pos >= 31890000 & pos <= 31970000]
leads <- rbindlist(lapply(c("AMD_FG", "AMD_IAMDGC"), function(k) win[trait == k][order(-abs(z))][1:3][, rank := 1:3]))
known <- g[SNP %in% cand & trait %in% c("AMD_FG", "AMD_IAMDGC")]
fwrite(rbind(leads[, .(trait, SNP, pos, ea, beta, se, z, set = paste0("window_top", rank))],
             known[, .(trait, SNP, pos, ea, beta, se, z, set = "known_candidate")]), file.path(OUT, "c4_leads.tsv"), sep = "\t")
print(leads[, .(trait, SNP, pos, z)]); print(known[, .(trait, SNP, z)])

bim <- fread(paste0(BFILE, ".bim"), header = FALSE, col.names = c("chr", "SNP", "cm", "bp", "A1", "A2"))[chr == 6]
pair_r <- function(a, b) {
  o <- file.path(W, sprintf("r_%s_%s", a, b))
  fwrite(data.table(c(a, b)), paste0(o, ".snps"), col.names = FALSE)
  system2(PLINK, c("--bfile", BFILE, "--chr", "6", "--extract", paste0(o, ".snps"), "--r", "--keep-allele-order",
                   "--ld-window-r2", "0", "--out", o), stdout = FALSE, stderr = FALSE)
  f <- paste0(o, ".ld"); if (!file.exists(f)) return(NA_real_); x <- fread(f); if (!nrow(x)) NA_real_ else x$R[1]
}
zA1 <- function(k, s) { x <- merge(bim[SNP == s, .(SNP, A1, A2)], g[trait == k & SNP == s], by = "SNP")
  if (!nrow(x)) return(NA_real_); with(x[1], fifelse(toupper(ea) == A1, z, fifelse(toupper(ea) == A2, -z, NA_real_))) }

c4snps <- c("rs1269852", "rs1270942")
conditioners <- unique(c(leads$SNP, known$SNP))
res <- rbindlist(lapply(c4snps, function(b) rbindlist(lapply(conditioners, function(a) {
  r <- pair_r(a, b); if (is.na(r)) return(NULL)
  rbindlist(lapply(TR, function(k) {
    za <- zA1(k, a); zb <- zA1(k, b); if (is.na(za) || is.na(zb)) return(NULL)
    rbindlist(lapply(c(r - 0.1, r, r + 0.1), function(rr) { rr <- max(min(rr, 0.99), -0.99)
      data.table(trait = k, c4_snp = b, conditioned_on = a, r_1000G = round(r, 3), r_used = round(rr, 3),
                 z_c4 = zb, z_cond_snp = za, p_c4 = 2 * pnorm(-abs(zb)),
                 pC_c4_given_cond = 2 * pnorm(-abs((zb - rr * za) / sqrt(1 - rr^2))),
                 pC_cond_given_c4 = 2 * pnorm(-abs((za - rr * zb) / sqrt(1 - rr^2)))) }))
  }))
}))))
fwrite(res, file.path(OUT, "c4_conditional.tsv"), sep = "\t")
print(res[abs(r_used - r_1000G) < 1e-9])
cat("C4 conditional done\n")
