# AMD <-> cerebral small vessel disease (CSVD): bidirectional MR, LDSC genetic correlation (incl. the cross-trait intercept
# matrix used by LAVA for sample overlap), and direction concordance of genome-wide significant AMD loci in CSVD traits.
# Instruments / methods as in 31_mr_bidirectional.R; sensitivity analyses drop major AMD loci (CFH, ARMS2/HTRA1, C3) and APOE.
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(TwoSampleMR); library(ieugwasr); library(MRPRESSO); library(ggplot2); library(GenomicSEM) })
setDTthreads(30); set.seed(20261006)
DER <- file.path(PROJ_ROOT, "derived/gwas"); OUT <- file.path(PROJ_ROOT, "results/csvd"); dir.create(file.path(OUT, "plots"), recursive = TRUE, showWarnings = FALSE)
PLINK <- file.path(PROJ_ROOT, "data/raw/tools/plink/plink.exe"); BFILE <- file.path(PROJ_ROOT, "data/raw/ref/g1000_eur/g1000_eur")
STEP <- commandArgs(TRUE); if (!length(STEP)) STEP <- c("mr", "ldsc", "concord")
cache <- new.env()
load_trait <- function(k) { if (is.null(cache[[k]])) cache[[k]] <- readRDS(file.path(DER, paste0(k, ".rds"))); cache[[k]] }

AMD <- c("AMD_FG", "DRYAMD_FG", "WETAMD_FG", "AMD_IAMDGC")
CSVD <- c("LACUNAR", "SVS", "WMH", "EPVS_WM", "EPVS_BG", "EPVS_HIP", "WM_MD", "WM_FA")
CSVD <- CSVD[file.exists(file.path(DER, paste0(CSVD, ".rds")))]
message("CSVD traits available: ", paste(CSVD, collapse = ", "))

# GRCh38, +/-500 kb
MAJOR <- data.table(gene = c("CFH", "ARMS2/HTRA1", "C3", "APOE"), chr = c("1", "10", "19", "19"),
                    start = c(196651878, 122454932, 6677704, 44905796), end = c(196747504, 122514908, 6720662, 44909393))
in_regions <- function(d, reg, flank = 5e5) {
  hit <- rep(FALSE, nrow(d))
  for (i in seq_len(nrow(reg))) hit <- hit | (d$chr == reg$chr[i] & d$pos >= reg$start[i] - flank & d$pos <= reg$end[i] + flank)
  hit
}
clump_local <- function(d, r2 = 0.001, kb = 10000) {
  x <- tryCatch(ld_clump(dplyr::tibble(rsid = d$SNP, pval = d$p), clump_kb = kb, clump_r2 = r2, clump_p = 1,
                         plink_bin = PLINK, bfile = BFILE), error = function(e) { message("clump error: ", conditionMessage(e)); NULL })
  if (is.null(x)) d[0] else d[SNP %in% x$rsid]
}
POSMAP <- NULL
get_iv <- function(k, pthr = 5e-8) {
  d <- copy(load_trait(k)[p < pthr])
  if (!nrow(d)) return(d)
  if (anyNA(d$pos)) {        # fill missing GRCh38 positions from FinnGen by rsID
    if (is.null(POSMAP)) POSMAP <<- load_trait("AMD_FG")[, .(SNP, chr38 = chr, pos38 = pos)]
    d <- merge(d, POSMAP, by = "SNP", all.x = TRUE)
    d[is.na(pos) & !is.na(pos38), `:=`(chr = chr38, pos = pos38)][, c("chr38", "pos38") := NULL]
  }
  d <- clump_local(d)
  d[, F := (beta / se)^2]
  d[F > 10]
}
fmt <- function(d, k, type) {
  x <- as.data.frame(d); x$Phenotype <- k
  format_data(x, type = type, snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "eaf",
              effect_allele_col = "ea", other_allele_col = "oa", pval_col = "p", samplesize_col = "n",
              ncase_col = if ("ncase" %in% names(x)) "ncase" else "ncase_", ncontrol_col = if ("ncontrol" %in% names(x)) "ncontrol" else "ncontrol_",
              chr_col = "chr", pos_col = "pos", phenotype_col = "Phenotype")
}
PREV <- c(AMD_IAMDGC = 0.02, AMD_FG = 0.03, DRYAMD_FG = 0.02, WETAMD_FG = 0.015, SVS = 0.005, LACUNAR = 0.005,
          EPVS_WM = 0.24, EPVS_BG = 0.23, EPVS_HIP = 0.24)
add_r2 <- function(dat) {
  for (side in c("exposure", "outcome")) {
    k <- dat[[side]][1]; b <- dat[[paste0("beta.", side)]]; eaf <- dat[[paste0("eaf.", side)]]
    if (all(is.na(eaf))) eaf <- dat[[paste0("eaf.", ifelse(side == "exposure", "outcome", "exposure"))]]
    eaf[is.na(eaf)] <- 0.3
    r <- if (k %in% names(PREV)) get_r_from_lor(b, eaf, dat[[paste0("ncase.", side)]][1], dat[[paste0("ncontrol.", side)]][1], PREV[[k]])
         else get_r_from_bsen(b, dat[[paste0("se.", side)]], dat[[paste0("samplesize.", side)]])
    dat[[paste0("rsq.", side)]] <- r^2; dat[[paste0("units.", side)]] <- ifelse(k %in% names(PREV), "log odds", "SD")
  }
  dat
}
run_pair <- function(ek, ok, iv, tag = "main", presso = TRUE) {
  if (nrow(iv) < 1) return(NULL)
  out <- load_trait(ok)[SNP %in% iv$SNP]; if (!nrow(out)) return(NULL)
  dat <- harmonise_data(fmt(iv, ek, "exposure"), fmt(out, ok, "outcome"), action = 2); dat <- dat[dat$mr_keep, ]
  if (nrow(dat) < 1) return(NULL)
  dat <- add_r2(dat)
  st <- tryCatch(directionality_test(dat), error = function(e) NULL)
  meth <- if (nrow(dat) >= 3) c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode")
          else if (nrow(dat) == 2) "mr_ivw" else "mr_wald_ratio"
  res <- mr(dat, method_list = meth)
  het <- if (nrow(dat) >= 3) mr_heterogeneity(dat, method_list = "mr_ivw") else NULL
  ple <- if (nrow(dat) >= 3) mr_pleiotropy_test(dat) else NULL
  pr <- NULL
  if (presso && nrow(dat) >= 4) pr <- tryCatch({
    p <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure",
                   OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = max(1000, 20 * nrow(dat)), SignifThreshold = 0.05)
    m <- p$`Main MR results`
    data.frame(presso_global_p = p$`MR-PRESSO results`$`Global Test`$Pvalue,
               presso_n_outlier = length(p$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`),
               presso_corr_b = m$`Causal Estimate`[2], presso_corr_se = m$Sd[2], presso_corr_p = m$`P-value`[2])
  }, error = function(e) NULL)
  res$exposure <- ek; res$outcome <- ok; res$analysis <- tag
  res$OR <- exp(res$b); res$lo <- exp(res$b - 1.96 * res$se); res$hi <- exp(res$b + 1.96 * res$se)
  summ <- data.frame(exposure = ek, outcome = ok, analysis = tag, nsnp = nrow(dat),
                     F_mean = mean((dat$beta.exposure / dat$se.exposure)^2), F_min = min((dat$beta.exposure / dat$se.exposure)^2),
                     Q = if (!is.null(het)) het$Q else NA, Q_p = if (!is.null(het)) het$Q_pval else NA,
                     egger_int = if (!is.null(ple)) ple$egger_intercept else NA, egger_int_p = if (!is.null(ple)) ple$pval else NA,
                     steiger_correct = if (!is.null(st)) st$correct_causal_direction else NA, steiger_p = if (!is.null(st)) st$steiger_pval else NA)
  if (!is.null(pr)) summ <- cbind(summ, pr)
  list(res = res, summ = summ, dat = dat)
}

# ------------------------------------------------------------------ MR
if ("mr" %in% STEP) {
  ivs <- list(); ivlog <- list()
  for (k in c(AMD, CSVD)) {
    f <- file.path(file.path(PROJ_ROOT, "results/mr"), paste0("instruments_", k, ".tsv"))
    if (file.exists(f) && k %in% c(AMD, "SVS", "WMH")) { iv <- fread(f); thr <- if (k == "SVS") 5e-6 else 5e-8 }   # reuse 31_ instruments
    else { iv <- get_iv(k, 5e-8); thr <- 5e-8; if (nrow(iv) < 3) { iv <- get_iv(k, 5e-6); thr <- 5e-6 } }
    iv[, chr := as.character(chr)]
    ivs[[k]] <- iv; ivlog[[k]] <- data.frame(trait = k, p_threshold = thr, n_iv = nrow(iv), F_min = if (nrow(iv)) min(iv$F) else NA)
    fwrite(iv, file.path(OUT, paste0("instruments_", k, ".tsv")), sep = "\t"); message(sprintf("IV %s: %d (p<%g)", k, nrow(iv), thr))
  }
  fwrite(rbindlist(ivlog), file.path(OUT, "mr_instrument_summary.tsv"), sep = "\t")
  pairs <- rbind(expand.grid(e = AMD, o = CSVD, dir = "AMD->CSVD", stringsAsFactors = FALSE),
                 expand.grid(e = CSVD, o = AMD, dir = "CSVD->AMD", stringsAsFactors = FALSE))
  R <- list(); S <- list(); D <- list()
  for (i in seq_len(nrow(pairs))) {
    ek <- pairs$e[i]; ok <- pairs$o[i]; key <- paste(ek, ok, sep = "__")
    message(sprintf("[%d/%d] %s -> %s", i, nrow(pairs), ek, ok))
    r <- run_pair(ek, ok, ivs[[ek]]); if (is.null(r)) next
    R[[key]] <- r$res; S[[key]] <- r$summ; D[[key]] <- r$dat
    if (ek %in% AMD) for (sv in list(list("excl_CFH_HTRA1_C3", MAJOR[gene != "APOE"]), list("excl_APOE", MAJOR[gene == "APOE"]),
                                     list("excl_CFH_HTRA1_C3_APOE", MAJOR))) {
      r2 <- run_pair(ek, ok, ivs[[ek]][!in_regions(ivs[[ek]], sv[[2]])], sv[[1]], presso = FALSE)
      if (!is.null(r2)) { R[[paste(key, sv[[1]])]] <- r2$res; S[[paste(key, sv[[1]])]] <- r2$summ }
    }
  }
  res <- rbindlist(R, fill = TRUE); summ <- rbindlist(S, fill = TRUE)
  ivw <- res[method %in% c("Inverse variance weighted", "Wald ratio")]
  ivw[analysis == "main", p_fdr := p.adjust(pval, "BH")]
  fwrite(res, file.path(OUT, "mr_all_methods.tsv"), sep = "\t"); fwrite(summ, file.path(OUT, "mr_diagnostics.tsv"), sep = "\t")
  fwrite(ivw, file.path(OUT, "mr_ivw.tsv"), sep = "\t"); saveRDS(D, file.path(OUT, "mr_harmonised.rds"))
  # forest plot: AMD -> CSVD, main vs excluding major loci
  fp <- ivw[exposure %in% AMD & analysis %in% c("main", "excl_CFH_HTRA1_C3_APOE")]
  fp[, outcome := factor(outcome, levels = rev(CSVD))]
  g <- ggplot(fp, aes(OR, outcome, colour = analysis)) + geom_vline(xintercept = 1, lty = 2) +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), size = 0.3) +
    facet_wrap(~exposure, nrow = 1, scales = "free_x") + scale_colour_manual(values = c(main = "#2b5c8a", excl_CFH_HTRA1_C3_APOE = "#c0392b")) +
    labs(x = "OR / beta-scale effect per log-odds AMD liability (IVW)", y = NULL, colour = NULL) + theme_bw(10) + theme(legend.position = "bottom")
  ggsave(file.path(OUT, "Fig_MR_AMD_to_CSVD_forest.pdf"), g, width = 12, height = 4.5); ggsave(file.path(OUT, "Fig_MR_AMD_to_CSVD_forest.png"), g, width = 12, height = 4.5, dpi = 200)
  fr <- ivw[outcome %in% AMD & analysis == "main"]; fr[, exposure := factor(exposure, levels = rev(CSVD))]
  g <- ggplot(fr, aes(OR, exposure)) + geom_vline(xintercept = 1, lty = 2) + geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.3, colour = "#2b5c8a") +
    facet_wrap(~outcome, nrow = 1, scales = "free_x") + labs(x = "OR of AMD per unit CSVD liability/trait (IVW)", y = NULL) + theme_bw(10)
  ggsave(file.path(OUT, "Fig_MR_CSVD_to_AMD_forest.pdf"), g, width = 12, height = 4); ggsave(file.path(OUT, "Fig_MR_CSVD_to_AMD_forest.png"), g, width = 12, height = 4, dpi = 200)
  print(ivw[, .(exposure, outcome, analysis, nsnp, OR = round(OR, 3), pval = signif(pval, 3), p_fdr = signif(p_fdr, 3))])
}

# ------------------------------------------------------------------ LDSC (also provides the intercept matrix for LAVA)
if ("ldsc" %in% STEP) {
  REF <- file.path(PROJ_ROOT, "data/raw/ref"); LD <- file.path(REF, "eur_w_ld_chr")
  OLD <- file.path(PROJ_ROOT, "derived/ldsc"); W <- file.path(PROJ_ROOT, "derived/ldsc_csvd"); dir.create(W, showWarnings = FALSE); setwd(W)
  HM3 <- file.path(OLD, "w_hm3.snplist"); if (!file.exists(HM3)) HM3 <- file.path(LD, "w_hm3.snplist")
  lt <- c("AMD_FG", "DRYAMD_FG", "WETAMD_FG", "AMD_IAMDGC", "IS", CSVD)
  neff <- function(d) if ("ncase" %in% names(d)) 4 / (1 / d$ncase[1] + 1 / d$ncontrol[1]) else d$n[1]
  todo <- lt[!file.exists(file.path(OLD, paste0(lt, ".sumstats.gz"))) & !file.exists(file.path(W, paste0(lt, ".sumstats.gz")))]
  if (length(todo)) {
    fs <- c(); Ns <- c()
    for (k in todo) { d <- load_trait(k); f <- file.path(W, paste0(k, ".txt")); fwrite(d[, .(SNP, A1 = ea, A2 = oa, effect = beta, P = p)], f, sep = "\t")
                      fs <- c(fs, f); Ns <- c(Ns, neff(d)); rm(d); rm(list = k, envir = cache); gc() }
    munge(files = fs, hm3 = HM3, trait.names = todo, N = Ns, parallel = TRUE, cores = min(length(todo), 10))
  }
  ss <- ifelse(file.exists(file.path(OLD, paste0(lt, ".sumstats.gz"))), file.path(OLD, paste0(lt, ".sumstats.gz")), file.path(W, paste0(lt, ".sumstats.gz")))
  binary <- c(AMD_FG = .03, DRYAMD_FG = .02, WETAMD_FG = .015, AMD_IAMDGC = .02, IS = .025, SVS = .005, LACUNAR = .005, EPVS_WM = .24, EPVS_BG = .23, EPVS_HIP = .24)
  sp <- ifelse(lt %in% c("EPVS_WM", "EPVS_BG", "EPVS_HIP"), 0.5, 0.5); pp <- binary[lt]   # munged with Neff -> sample prevalence 0.5
  sp[!lt %in% names(binary)] <- NA; pp[!lt %in% names(binary)] <- NA
  res <- ldsc(traits = ss, sample.prev = sp, population.prev = as.numeric(pp), ld = LD, wld = LD, trait.names = lt)
  S <- res$S; V <- res$V; I <- res$I; k <- nrow(S); dimnames(I) <- list(lt, lt)
  sdg <- sqrt(diag(S)); Rg <- S / outer(sdg, sdg)
  idx <- which(lower.tri(S, diag = TRUE), arr.ind = TRUE); vidx <- function(i, j) which((idx[, 1] == max(i, j)) & (idx[, 2] == min(i, j)))
  out <- list()
  for (i in 1:k) for (j in 1:k) if (i < j) {
    a <- vidx(i, j); b <- vidx(i, i); c <- vidx(j, j)
    gg <- c(1 / sqrt(S[i, i] * S[j, j]), -0.5 * S[i, j] / (S[i, i]^1.5 * sqrt(S[j, j])), -0.5 * S[i, j] / (S[j, j]^1.5 * sqrt(S[i, i])))
    se <- sqrt(as.numeric(t(gg) %*% V[c(a, b, c), c(a, b, c)] %*% gg))
    out[[length(out) + 1]] <- data.table(trait1 = lt[i], trait2 = lt[j], rg = Rg[i, j], se = se, p = 2 * pnorm(-abs(Rg[i, j] / se)), gcov_int = I[i, j])
  }
  out <- rbindlist(out)
  h2 <- data.table(trait = lt, h2 = diag(S), se = sqrt(diag(V)[sapply(1:k, function(i) vidx(i, i))]), intercept = diag(I))
  fwrite(out, file.path(OUT, "ldsc_rg.tsv"), sep = "\t"); fwrite(h2, file.path(OUT, "ldsc_h2.tsv"), sep = "\t")
  ov <- round(cov2cor(I), 5); fwrite(as.data.table(ov, keep.rownames = "trait"), file.path(OUT, "ldsc_intercept_matrix.tsv"), sep = "\t")
  dir.create(file.path(PROJ_ROOT, "derived/lava"), showWarnings = FALSE); write.table(ov, file.path(file.path(PROJ_ROOT, "derived/lava"), "sample_overlap.txt"), quote = FALSE)
  am <- out[(trait1 %in% AMD) != (trait2 %in% AMD)]; am[, amd := ifelse(trait1 %in% AMD, trait1, trait2)]; am[, other := ifelse(trait1 %in% AMD, trait2, trait1)]
  am[, other := factor(other, levels = c("IS", CSVD))]; am[, lab := sprintf("%.2f%s", rg, ifelse(p < 0.05, "*", ""))]
  g <- ggplot(am, aes(other, amd, fill = rg)) + geom_tile(colour = "white") + geom_text(aes(label = lab), size = 3) +
    scale_fill_gradient2(low = "#2166ac", high = "#b2182b", limits = c(-1, 1)) + labs(x = NULL, y = NULL, fill = "rg") + theme_minimal(10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(OUT, "Fig_LDSC_rg_heatmap.pdf"), g, width = 8, height = 3.2); ggsave(file.path(OUT, "Fig_LDSC_rg_heatmap.png"), g, width = 8, height = 3.2, dpi = 200)
  print(h2); print(am[order(p)][, .(amd, other, rg = round(rg, 3), se = round(se, 3), p = signif(p, 3))])
}

# ------------------------------------------------------------------ direction concordance of AMD loci + key variants
if ("concord" %in% STEP) {
  key <- data.table(SNP = c("rs3750848", "rs61871745", "rs2672592", "rs429358", "rs7412", "rs3184504", "rs10737680", "rs2230199"),
                    locus = c("10q26 ARMS2 (AMD lead)", "10q26 ARMS2", "10q26 HTRA1 (Malik 2024)", "APOE e4", "APOE e2", "12q24 SH2B3 R262W", "1q31 CFH", "19p13 C3 R102G"))
  hz <- fread(file.path(PROJ_ROOT, "results/smr_heidi/smr_heidi_all.tsv"))
  c4 <- unique(hz[GeneSymbol %in% c("C4A", "C4B") & outcome == "WMH", .(SNP = topSNP, locus = paste0("6p21 MHC (", GeneSymbol, " eQTL top SNP)"))])
  key <- unique(rbind(key, c4), by = "SNP")
  for (amd in c("AMD_FG", "AMD_IAMDGC")) {
    iv <- fread(file.path(file.path(PROJ_ROOT, "results/mr"), paste0("instruments_", amd, ".tsv")))   # clumped GW-significant AMD loci
    iv <- iv[, .(SNP, chr = as.character(chr), pos, ea, oa, b_amd = beta, se_amd = se, p_amd = p)]
    iv[b_amd < 0, `:=`(ea = oa, oa = ea, b_amd = -b_amd)]                                 # align to AMD risk allele
    iv[, major := in_regions(iv, MAJOR)]
    rows <- list(); cons <- list()
    for (k in CSVD) {
      d <- load_trait(k)[SNP %in% c(iv$SNP, key$SNP), .(SNP, ea2 = ea, oa2 = oa, b = beta, se, p)]
      m <- merge(iv, d, by = "SNP")
      m[, b_al := fifelse(ea2 == ea & oa2 == oa, b, fifelse(ea2 == oa & oa2 == ea, -b, NA_real_))]
      m <- m[!is.na(b_al)]; m[, trait := k]; rows[[k]] <- m
      for (sub in c("all", "excl_major")) {
        mm <- if (sub == "all") m else m[major == FALSE]
        for (thr in c(1, 0.05)) {
          x <- mm[p < thr]; nconc <- sum(x$b_al > 0); n <- nrow(x)
          cons[[length(cons) + 1]] <- data.table(amd = amd, trait = k, loci = sub, p_threshold = thr, n = n, n_same_dir = nconc,
                                                 frac_same = if (n) nconc / n else NA, binom_p = if (n) binom.test(nconc, n)$p.value else NA)
        }
      }
      # key variants, aligned to the AMD risk allele as given in the AMD GWAS
      a <- load_trait(amd)[SNP %in% key$SNP, .(SNP, ea, oa, b_amd = beta, se_amd = se, p_amd = p)]
      a[b_amd < 0, `:=`(ea = oa, oa = ea, b_amd = -b_amd)]
      kk <- merge(merge(a, d, by = "SNP"), key, by = "SNP")
      kk[, b_al := fifelse(ea2 == ea & oa2 == oa, b, fifelse(ea2 == oa & oa2 == ea, -b, NA_real_))]
      rows[[paste0(k, "_key")]] <- kk[, .(SNP, locus, amd, trait = k, risk_allele = ea, b_amd, p_amd, b_csvd = b_al, se_csvd = se, p_csvd = p)]
    }
    cons <- rbindlist(cons); fwrite(cons, file.path(OUT, paste0("concordance_", amd, ".tsv")), sep = "\t")
    allm <- rbindlist(rows[CSVD], fill = TRUE); fwrite(allm, file.path(OUT, paste0("concordance_loci_", amd, ".tsv")), sep = "\t")
    keyt <- rbindlist(rows[paste0(CSVD, "_key")], fill = TRUE); fwrite(keyt, file.path(OUT, paste0("key_variants_", amd, ".tsv")), sep = "\t")
    allm[, sig := p < 0.05]
    g <- ggplot(allm, aes(b_amd, b_al / se)) + geom_hline(yintercept = 0, lty = 2) + geom_hline(yintercept = c(-1.96, 1.96), lty = 3, colour = "grey60") +
      geom_point(aes(colour = major, shape = sig), size = 1.8, alpha = 0.8) + facet_wrap(~trait, scales = "free_y", nrow = 2) +
      scale_colour_manual(values = c(`FALSE` = "#2b5c8a", `TRUE` = "#c0392b"), labels = c("other loci", "CFH/HTRA1/C3/APOE")) +
      labs(x = paste0("AMD log-OR (risk allele, ", amd, ")"), y = "CSVD z-score for the AMD risk allele", colour = NULL, shape = "CSVD P<0.05") + theme_bw(9)
    ggsave(file.path(OUT, paste0("Fig_concordance_", amd, ".pdf")), g, width = 11, height = 5.5); ggsave(file.path(OUT, paste0("Fig_concordance_", amd, ".png")), g, width = 11, height = 5.5, dpi = 200)
    keyt[, `:=`(lo = b_csvd - 1.96 * se_csvd, hi = b_csvd + 1.96 * se_csvd)]
    g <- ggplot(keyt, aes(b_csvd, trait)) + geom_vline(xintercept = 0, lty = 2) + geom_pointrange(aes(xmin = lo, xmax = hi, colour = p_csvd < 0.05), size = 0.25) +
      facet_wrap(~locus, scales = "free_x", nrow = 2) + scale_colour_manual(values = c(`FALSE` = "grey50", `TRUE` = "#c0392b")) +
      labs(x = "CSVD effect of the AMD risk allele (beta / log-OR)", y = NULL, colour = "P<0.05") + theme_bw(9)
    ggsave(file.path(OUT, paste0("Fig_key_variants_", amd, ".pdf")), g, width = 13, height = 6); ggsave(file.path(OUT, paste0("Fig_key_variants_", amd, ".png")), g, width = 13, height = 6, dpi = 200)
    print(cons[loci == "all"]); print(keyt[, .(SNP, locus, trait, b_csvd = round(b_csvd, 4), p_csvd = signif(p_csvd, 3))])
  }
}
message("CSVD MR/LDSC done")
