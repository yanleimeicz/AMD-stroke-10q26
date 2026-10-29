# Transcriptome-wide SMR/MR using eQTLGen whole-blood cis-eQTLs (FDR < 0.05 file) for AMD and stroke outcomes.
# Primary instrument per gene = top cis-eQTL (p < 5e-8) present in the outcome GWAS (SMR approach, Zhu et al. 2016).
# beta_eQTL from z: b = z / sqrt(2p(1-p)(N + z^2)), p = 1000G EUR frequency of the assessed allele.
# For genes significant in both an AMD and a stroke/WMH outcome: multi-SNP IVW (clumped r2 < 0.1) as sensitivity.
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(ieugwasr) })
setDTthreads(30)
DER <- file.path(PROJ_ROOT, "derived/gwas"); OUT <- file.path(PROJ_ROOT, "results/smr"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PLINK <- file.path(PROJ_ROOT, "data/raw/tools/plink/plink.exe"); BFILE <- file.path(PROJ_ROOT, "data/raw/ref/g1000_eur/g1000_eur")
EQTL <- file.path(PROJ_ROOT, "data/raw/eqtl/eqtlgen_cis_sig.txt.gz")

# 1000G EUR allele frequencies
frq_file <- file.path(PROJ_ROOT, "derived/g1000_eur.frq")
if (!file.exists(frq_file)) system2(PLINK, c("--bfile", BFILE, "--freq", "--out", sub("\\.frq$", "", frq_file)), stdout = FALSE)
frq <- fread(frq_file)[, .(SNP, A1, A2, MAF)]

eq <- fread(EQTL, select = c("Pvalue", "SNP", "SNPChr", "SNPPos", "AssessedAllele", "OtherAllele", "Zscore", "Gene", "GeneSymbol", "NrSamples"))
eq <- eq[Pvalue < 5e-8]
eq <- merge(eq, frq, by = "SNP")
eq[, p_ass := fifelse(AssessedAllele == A1, MAF, fifelse(AssessedAllele == A2, 1 - MAF, NA_real_))]
eq <- eq[!is.na(p_ass) & p_ass > 0.01 & p_ass < 0.99]
eq[, `:=`(bx = Zscore / sqrt(2 * p_ass * (1 - p_ass) * (NrSamples + Zscore^2)),
          sx = 1 / sqrt(2 * p_ass * (1 - p_ass) * (NrSamples + Zscore^2)))]
message(sprintf("eQTLGen: %d SNP-gene pairs (p<5e-8), %d genes", nrow(eq), uniqueN(eq$Gene)))

avail <- setdiff(sub("\\.rds$", "", list.files(DER, "\\.rds$")), "trait_meta")
OUTCOMES <- intersect(c("AMD_IAMDGC", "AMD_FG", "DRYAMD_FG", "WETAMD_FG", "AS", "IS", "LAS", "CES", "SVS", "WMH"), avail)

smr_one <- function(k) {
  g <- readRDS(file.path(DER, paste0(k, ".rds")))[SNP %in% eq$SNP, .(SNP, ea, oa, by = beta, sy = se, py = p)]
  m <- merge(eq, g, by = "SNP")
  m[, by_al := fifelse(ea == AssessedAllele & oa == OtherAllele, by, fifelse(ea == OtherAllele & oa == AssessedAllele, -by, NA_real_))]
  m <- m[!is.na(by_al)]
  setorder(m, Gene, Pvalue)
  top <- m[, .SD[1], by = Gene]
  top[, `:=`(zx = bx / sx, zy = by_al / sy)]
  top[, b_smr := by_al / bx]
  top[, T_smr := zx^2 * zy^2 / (zx^2 + zy^2)]
  top[, p_smr := pchisq(T_smr, 1, lower.tail = FALSE)]
  top[, se_smr := abs(b_smr) / sqrt(T_smr)]
  top[, `:=`(outcome = k, p_fdr = p.adjust(p_smr, "BH"))]
  top[, .(outcome, Gene, GeneSymbol, topSNP = SNP, chr = SNPChr, pos37 = SNPPos, p_eqtl = Pvalue, bx, sx, by = by_al, sy, py,
          b_smr, se_smr, OR = exp(b_smr), p_smr, p_fdr)]
}
smr <- rbindlist(lapply(OUTCOMES, function(k) { message("SMR ", k); smr_one(k) }))
fwrite(smr, file.path(OUT, "smr_all.tsv"), sep = "\t")

AMDS <- intersect(c("AMD_IAMDGC", "AMD_FG", "DRYAMD_FG", "WETAMD_FG"), OUTCOMES)
STRK <- intersect(c("AS", "IS", "LAS", "CES", "SVS", "WMH"), OUTCOMES)
sig_a <- smr[outcome %in% AMDS & p_fdr < 0.05, unique(Gene)]
sig_s <- smr[outcome %in% STRK & p_fdr < 0.05, unique(Gene)]
nom_s <- smr[outcome %in% STRK & p_smr < 0.05, unique(Gene)]
shared_strict <- intersect(sig_a, sig_s)
shared_loose <- intersect(sig_a, nom_s)
message(sprintf("AMD FDR genes: %d; stroke/WMH FDR genes: %d; shared (FDR both): %d; shared (AMD FDR & stroke p<0.05): %d",
                length(sig_a), length(sig_s), length(shared_strict), length(shared_loose)))
wide <- dcast(smr[Gene %in% shared_loose], Gene + GeneSymbol ~ outcome, value.var = c("OR", "p_smr", "p_fdr"))
wide[, tier := fifelse(Gene %in% shared_strict, "FDR in both", "AMD FDR + stroke nominal")]
fwrite(wide, file.path(OUT, "smr_shared_genes.tsv"), sep = "\t")

# multi-SNP IVW sensitivity for shared genes
ivw_gene <- function(gene, k) {
  e <- eq[Gene == gene]
  g <- readRDS(file.path(DER, paste0(k, ".rds")))[SNP %in% e$SNP, .(SNP, ea, oa, by = beta, sy = se)]
  m <- merge(e, g, by = "SNP")
  m[, by_al := fifelse(ea == AssessedAllele & oa == OtherAllele, by, fifelse(ea == OtherAllele & oa == AssessedAllele, -by, NA_real_))]
  m <- m[!is.na(by_al)]
  if (nrow(m) > 1) {
    cl <- tryCatch(ld_clump(dplyr::tibble(rsid = m$SNP, pval = m$Pvalue), clump_kb = 1000, clump_r2 = 0.1, clump_p = 1,
                            plink_bin = PLINK, bfile = BFILE), error = function(e) NULL)
    if (!is.null(cl)) m <- m[SNP %in% cl$rsid]
  }
  if (!nrow(m)) return(NULL)
  w <- m$bx^2 / m$sy^2; b <- sum(w * m$by_al / m$bx) / sum(w); se <- 1 / sqrt(sum(w))
  if (nrow(m) >= 2) {   # multiplicative random effects
    q <- sum(w * (m$by_al / m$bx - b)^2); phi <- max(1, q / (nrow(m) - 1)); se <- se * sqrt(phi)
  }
  data.table(Gene = gene, outcome = k, nsnp = nrow(m), b_ivw = b, se_ivw = se, OR_ivw = exp(b), p_ivw = 2 * pnorm(-abs(b / se)))
}
if (length(shared_loose)) {
  iv <- rbindlist(lapply(shared_loose, function(gn) rbindlist(lapply(OUTCOMES, function(k) ivw_gene(gn, k)))))
  iv <- merge(iv, unique(eq[, .(Gene, GeneSymbol)]), by = "Gene")
  fwrite(iv, file.path(OUT, "smr_shared_genes_multiSNP_ivw.tsv"), sep = "\t")
}
print(wide[order(tier)][, 1:min(12, ncol(wide))])
message("SMR done")
