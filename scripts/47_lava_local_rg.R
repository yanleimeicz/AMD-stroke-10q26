# Local genetic correlation (LAVA; Werme et al. 2022) between AMD and cerebral small vessel disease / stroke traits.
# 2,495 semi-independent LD blocks (LAVA blocks_s2500_m25_f1_w200, GRCh37); LD reference 1000G EUR (MAF >= 0.5%);
# sample overlap from the LDSC cross-trait intercept matrix (46_csvd_mr_ldsc.R). Univariate filter P < 1e-4 for both traits,
# bivariate local rg with FDR (BH) per AMD target. Run in 22 parallel workers (one per chromosome).
# Usage: Rscript 47_lava_local_rg.R [prep|run|summ]...
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(parallel); library(ggplot2) })
setDTthreads(30)
DER <- file.path(PROJ_ROOT, "derived/gwas"); W <- file.path(PROJ_ROOT, "derived/lava"); OUT <- file.path(PROJ_ROOT, "results/csvd")
dir.create(file.path(W, "ss"), recursive = TRUE, showWarnings = FALSE); dir.create(file.path(W, "ref"), showWarnings = FALSE)
dir.create(file.path(W, "res"), showWarnings = FALSE); dir.create(OUT, showWarnings = FALSE)
PLINK <- file.path(PROJ_ROOT, "data/raw/tools/plink/plink.exe"); BFILE <- file.path(PROJ_ROOT, "data/raw/ref/g1000_eur/g1000_eur")
LOCF <- file.path(PROJ_ROOT, "data/raw/ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
STEP <- commandArgs(TRUE); if (!length(STEP)) STEP <- c("prep", "run", "summ")

AMD <- c("AMD_FG", "DRYAMD_FG", "WETAMD_FG", "AMD_IAMDGC")
OTHER <- c("IS", "SVS", "LACUNAR", "WMH", "EPVS_WM", "EPVS_BG", "EPVS_HIP", "WM_MD", "WM_FA")
PH <- c(AMD, OTHER); PH <- PH[file.exists(file.path(DER, paste0(PH, ".rds")))]
ov <- as.matrix(read.table(file.path(W, "sample_overlap.txt"), check.names = FALSE)); PH <- intersect(PH, rownames(ov))
message("phenotypes: ", paste(PH, collapse = ", "))

if ("prep" %in% STEP) {
  cl <- makeCluster(22); clusterExport(cl, c("PLINK", "BFILE", "W"))
  parLapply(cl, 1:22, function(c) {
    o <- file.path(W, "ref", paste0("g1000_eur_chr", c))
    if (!file.exists(paste0(o, ".bim"))) system2(PLINK, c("--bfile", BFILE, "--chr", c, "--maf", "0.005", "--snps-only", "just-acgt", "--memory", "3000", "--threads", "1", "--make-bed", "--out", o), stdout = FALSE)
  }); stopCluster(cl)
  bim <- rbindlist(lapply(1:22, function(c) fread(file.path(W, "ref", paste0("g1000_eur_chr", c, ".bim")), header = FALSE, select = c(1, 2), col.names = c("CHR", "SNP"))))
  info <- list()
  for (k in PH) {
    d <- readRDS(file.path(DER, paste0(k, ".rds")))
    if (!"ncase" %in% names(d)) d[, `:=`(ncase = NA_real_, ncontrol = NA_real_)]
    d <- d[, .(SNP, A1 = ea, A2 = oa, Z = beta / se, N = n, ncase, ncontrol)]
    d <- merge(d, bim, by = "SNP")
    for (c in 1:22) fwrite(d[CHR == c, .(SNP, A1, A2, Z, N)], file.path(W, "ss", paste0(k, "_chr", c, ".txt")), sep = "\t")
    info[[k]] <- data.table(phenotype = k, cases = d$ncase[1], controls = d$ncontrol[1])
    message(sprintf("  %s: %d SNPs in reference", k, nrow(d))); rm(d); gc()
  }
  fwrite(rbindlist(info), file.path(W, "pheno_info.tsv"), sep = "\t")
}

if ("run" %in% STEP) {
  info <- fread(file.path(W, "pheno_info.tsv"))[phenotype %in% PH]
  loci_all <- fread(LOCF)
  worker <- function(c) {
    suppressPackageStartupMessages({ library(LAVA); library(data.table) })
    f <- file.path(W, "res", paste0("chr", c, ".rds")); if (file.exists(f)) return(f)
    inf <- copy(info)[, filename := file.path(W, "ss", paste0(phenotype, "_chr", c, ".txt"))]
    inf[is.na(cases), `:=`(cases = NA, controls = NA)]
    fi <- file.path(W, "res", paste0("info_chr", c, ".txt")); fwrite(inf, fi, sep = "\t", na = "NA")
    input <- process.input(input.info.file = fi, sample.overlap.file = file.path(W, "sample_overlap.txt"),
                           ref.prefix = file.path(W, "ref", paste0("g1000_eur_chr", c)), phenos = inf$phenotype)
    loci <- loci_all[CHR == c]
    U <- list(); B <- list()
    for (i in seq_len(nrow(loci))) {
      locus <- tryCatch(process.locus(loci[i], input), error = function(e) NULL)
      if (is.null(locus)) next
      u <- tryCatch(run.univ(locus), error = function(e) NULL); if (is.null(u)) next
      u <- as.data.table(u); u[, `:=`(locus = loci$LOC[i], chr = c, start = loci$START[i], stop = loci$STOP[i], n_snps = locus$n.snps)]
      U[[length(U) + 1]] <- u
      pass <- u[p < 1e-4, phen]
      for (t in intersect(AMD, pass)) {
        oth <- setdiff(intersect(OTHER, pass), t); if (!length(oth)) next
        b <- tryCatch(run.bivar(locus, phenos = c(t, oth), target = t), error = function(e) NULL)
        if (!is.null(b)) { b <- as.data.table(b); b[, `:=`(locus = loci$LOC[i], chr = c, start = loci$START[i], stop = loci$STOP[i])]; B[[length(B) + 1]] <- b }
      }
    }
    saveRDS(list(univ = rbindlist(U, fill = TRUE), bivar = rbindlist(B, fill = TRUE)), f); f
  }
  cl <- makeCluster(22); clusterExport(cl, c("W", "info", "loci_all", "AMD", "OTHER", "PH"))
  ord <- order(-table(loci_all$CHR)[as.character(1:22)])
  fs <- parLapplyLB(cl, (1:22)[ord], worker); stopCluster(cl)
  message("LAVA run done: ", length(fs), " chromosomes")
}

if ("summ" %in% STEP) {
  rr <- lapply(1:22, function(c) readRDS(file.path(W, "res", paste0("chr", c, ".rds"))))
  univ <- rbindlist(lapply(rr, `[[`, "univ"), fill = TRUE); biv <- rbindlist(lapply(rr, `[[`, "bivar"), fill = TRUE)
  biv[, p_fdr := p.adjust(p, "BH"), by = phen1]
  # annotate key regions (GRCh37)
  REG <- data.table(region = c("1q31 CFH", "10q26 ARMS2/HTRA1", "19p13 C3", "19q13 APOE", "6p21 MHC/C4", "12q24 SH2B3", "3p14 ADAMTS9", "9q34 ABO", "7q21 CDK6"),
                    chr = c(1, 10, 19, 19, 6, 12, 3, 9, 7), s = c(196621008, 124134094, 6677846, 45409011, 28477797, 111843752, 64501333, 136125788, 92234235),
                    e = c(196716634, 124274424, 6730573, 45412650, 33448354, 111889427, 64673679, 136150630, 92465908))
  ann <- function(d) { d[, region := NA_character_]
    for (i in seq_len(nrow(REG))) d[chr == REG$chr[i] & start <= REG$e[i] & stop >= REG$s[i], region := REG$region[i]]; d }
  univ <- ann(univ); biv <- ann(biv)
  fwrite(univ, file.path(OUT, "lava_univariate.tsv"), sep = "\t"); fwrite(biv, file.path(OUT, "lava_bivariate.tsv"), sep = "\t")
  sig <- biv[p_fdr < 0.05][order(p)]
  fwrite(sig, file.path(OUT, "lava_bivariate_fdr05.tsv"), sep = "\t")
  cnt <- biv[, .(n_tested = .N, n_nominal = sum(p < 0.05), n_fdr = sum(p_fdr < 0.05), n_fdr_neg = sum(p_fdr < 0.05 & rho < 0), n_fdr_pos = sum(p_fdr < 0.05 & rho > 0)), by = .(phen1, phen2)]
  fwrite(cnt, file.path(OUT, "lava_counts.tsv"), sep = "\t")
  print(cnt); print(sig[, .(phen1, phen2, locus, chr, start, stop, region, rho = round(rho, 3), rho.lower, rho.upper, p = signif(p, 3), p_fdr = signif(p_fdr, 3))])
  # genome plot of bivariate local rg
  cs <- cumsum(c(0, tapply(biv$stop, biv$chr, max)[as.character(1:21)])); names(cs) <- 1:22
  biv[, gpos := (start + stop) / 2 + cs[as.character(chr)]]
  biv[, lab := ifelse(p_fdr < 0.05, ifelse(is.na(region), paste0("chr", chr, ":", round(start / 1e6, 1), "M"), region), NA)]
  g <- ggplot(biv, aes(gpos, -log10(p) * sign(rho))) + geom_hline(yintercept = 0) +
    geom_point(aes(colour = phen2, size = p_fdr < 0.05), alpha = 0.8) + scale_size_manual(values = c(`FALSE` = 1.2, `TRUE` = 2.6)) +
    ggrepel::geom_text_repel(aes(label = lab), size = 2.4, max.overlaps = 40, na.rm = TRUE) +
    facet_wrap(~phen1, ncol = 1) + labs(x = "genomic position (chr 1-22)", y = "sign(local rg) x -log10 P", colour = NULL, size = "FDR<0.05") +
    theme_bw(9) + theme(axis.text.x = element_blank())
  ggsave(file.path(OUT, "Fig_LAVA_local_rg.pdf"), g, width = 12, height = 10); ggsave(file.path(OUT, "Fig_LAVA_local_rg.png"), g, width = 12, height = 10, dpi = 200)
  if (nrow(sig)) {
    sig[, key := paste0(ifelse(is.na(region), paste0("chr", chr, ":", round(start / 1e6, 2), "M"), region), " | ", phen1)]
    g <- ggplot(sig, aes(rho, reorder(key, rho), colour = phen2)) + geom_vline(xintercept = 0, lty = 2) +
      geom_pointrange(aes(xmin = pmax(-1, rho.lower), xmax = pmin(1, rho.upper)), position = position_dodge(width = 0.6), size = 0.3) +
      labs(x = "local genetic correlation (95% CI)", y = NULL, colour = NULL) + theme_bw(9)
    ggsave(file.path(OUT, "Fig_LAVA_sig_loci.pdf"), g, width = 8, height = max(3, 0.25 * nrow(sig) + 1.5))
    ggsave(file.path(OUT, "Fig_LAVA_sig_loci.png"), g, width = 8, height = max(3, 0.25 * nrow(sig) + 1.5), dpi = 200)
  }
  message("LAVA summary done")
}
