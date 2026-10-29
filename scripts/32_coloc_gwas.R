# GWAS-GWAS colocalisation between AMD and stroke subtypes / WMH at genome-wide significant loci of either trait.
# Locus = lead SNP (clumped, p < 5e-8) +/- 500 kb (GRCh38 positions of the AMD dataset; matched by rsID).
# coloc.abf with default priors (p1 = p2 = 1e-4, p12 = 1e-5); PP.H4 > 0.8 = strong, 0.5-0.8 = suggestive.
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(coloc); library(ieugwasr) })
setDTthreads(16)
DER <- file.path(PROJ_ROOT, "derived/gwas"); OUT <- file.path(PROJ_ROOT, "results/coloc"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PLINK <- file.path(PROJ_ROOT, "data/raw/tools/plink/plink.exe"); BFILE <- file.path(PROJ_ROOT, "data/raw/ref/g1000_eur/g1000_eur")
avail <- setdiff(sub("\\.rds$", "", list.files(DER, "\\.rds$")), "trait_meta")
A_TRAITS <- intersect(c("AMD_IAMDGC", "AMD_FG", "DRYAMD_FG", "WETAMD_FG"), avail)
S_TRAITS <- intersect(c("AS", "IS", "LAS", "CES", "SVS", "WMH"), avail)
L <- lapply(setNames(c(A_TRAITS, S_TRAITS), c(A_TRAITS, S_TRAITS)), function(k) readRDS(file.path(DER, paste0(k, ".rds"))))
CC <- c(A_TRAITS, setdiff(S_TRAITS, "WMH"))
POS_REF <- if ("AMD_FG" %in% avail) L$AMD_FG[, .(SNP, chr38 = chr, pos38 = pos)] else L[[A_TRAITS[1]]][, .(SNP, chr38 = chr, pos38 = pos)]

lead_loci <- function(k) {
  d <- L[[k]][p < 5e-8]
  if (!nrow(d)) return(NULL)
  cl <- ld_clump(dplyr::tibble(rsid = d$SNP, pval = d$p), clump_kb = 1000, clump_r2 = 0.01, clump_p = 1, plink_bin = PLINK, bfile = BFILE)
  x <- merge(d[SNP %in% cl$rsid, .(SNP, p)], POS_REF, by = "SNP")
  x[, trait := k]
  x
}

ds <- function(k, d) {
  eaf <- d$eaf; eaf[is.na(eaf)] <- 0.3
  maf <- pmin(eaf, 1 - eaf); maf[maf <= 0] <- 1e-4
  x <- list(beta = d$beta, varbeta = d$se^2, snp = d$SNP, MAF = maf, N = d$n[1])
  if (k %in% CC) { x$type <- "cc"; x$s <- d$ncase[1] / d$n[1] } else { x$type <- "quant"; x$sdY <- 1 }
  x
}

loci <- rbindlist(lapply(c(A_TRAITS, S_TRAITS), lead_loci), fill = TRUE)
fwrite(loci, file.path(OUT, "lead_loci.tsv"), sep = "\t")
message("lead loci: ", paste(names(table(loci$trait)), table(loci$trait), collapse = ", "))

res <- list()
for (i in seq_len(nrow(loci))) {
  lc <- loci[i]; src_is_amd <- lc$trait %in% A_TRAITS
  partners <- if (src_is_amd) S_TRAITS else A_TRAITS
  win <- POS_REF[chr38 == lc$chr38 & pos38 >= lc$pos38 - 5e5 & pos38 <= lc$pos38 + 5e5, SNP]
  d1 <- L[[lc$trait]][SNP %in% win]
  for (k2 in partners) {
    d2 <- L[[k2]][SNP %in% win]
    m <- merge(d1, d2, by = "SNP", suffixes = c(".1", ".2"))
    # align alleles
    flip <- m$ea.1 == m$oa.2 & m$oa.1 == m$ea.2
    keep <- (m$ea.1 == m$ea.2 & m$oa.1 == m$oa.2) | flip
    m[flip, beta.2 := -beta.2]; m <- m[keep]
    if (nrow(m) < 50) next
    a <- m[, .(SNP, beta = beta.1, se = se.1, eaf = eaf.1, n = n.1, ncase = if ("ncase.1" %in% names(m)) ncase.1 else NA)]
    b <- m[, .(SNP, beta = beta.2, se = se.2, eaf = ifelse(is.na(eaf.2), eaf.1, eaf.2), n = n.2, ncase = if ("ncase.2" %in% names(m)) ncase.2 else NA)]
    cr <- tryCatch(coloc.abf(ds(lc$trait, a), ds(k2, b)), error = function(e) NULL)
    if (is.null(cr)) next
    top <- cr$results[order(-cr$results$SNP.PP.H4), ][1, ]
    res[[length(res) + 1]] <- data.table(source = lc$trait, lead_snp = lc$SNP, chr = lc$chr38, pos38 = lc$pos38, partner = k2,
                                         nsnps = cr$summary["nsnps"], PP.H0 = cr$summary["PP.H0.abf"], PP.H1 = cr$summary["PP.H1.abf"],
                                         PP.H2 = cr$summary["PP.H2.abf"], PP.H3 = cr$summary["PP.H3.abf"], PP.H4 = cr$summary["PP.H4.abf"],
                                         top_snp_H4 = top$snp, top_snp_PP = top$SNP.PP.H4,
                                         min_p_source = min(m$p.1), min_p_partner = min(m$p.2))
  }
  if (i %% 20 == 0) message(sprintf("  %d/%d loci", i, nrow(loci)))
}
res <- rbindlist(res)
# annotate nearest gene from FinnGen (if available)
if ("AMD_FG" %in% avail) {
  fgf <- file.path(PROJ_ROOT, "data/raw/gwas/finngen_R12_H7_AMD.gz")
  ng <- fread(fgf, select = c("rsids", "nearest_genes"))[, .(SNP = sub(",.*", "", rsids), nearest_genes)]
  res <- merge(res, ng[!duplicated(SNP)], by.x = "lead_snp", by.y = "SNP", all.x = TRUE)
}
setorder(res, -PP.H4)
fwrite(res, file.path(OUT, "coloc_gwas_gwas.tsv"), sep = "\t")
print(res[PP.H4 > 0.5, .(source, lead_snp, chr, pos38, partner, nsnps, PP.H3 = round(PP.H3, 3), PP.H4 = round(PP.H4, 3), top_snp_H4, nearest_genes)])
message("coloc done")
