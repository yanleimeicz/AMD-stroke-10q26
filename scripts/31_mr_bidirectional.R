# Bidirectional two-sample MR: AMD <-> stroke (GIGASTROKE subtypes) and AMD <-> WMH.
# Instruments: p < 5e-8 (fallback 5e-6 when < 3 SNPs), LD clump r2 < 0.001 / 10,000 kb (1000G EUR, local plink), F > 10.
# Methods: IVW (MRE), MR-Egger, weighted median, simple mode, weighted mode; Cochran Q, Egger intercept,
# MR-PRESSO, Steiger, leave-one-out; sensitivity excluding major AMD (CFH, ARMS2/HTRA1, C3) and lipid loci.
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(TwoSampleMR); library(ieugwasr); library(MRPRESSO); library(ggplot2) })
setDTthreads(16); set.seed(20261005)
DER <- file.path(PROJ_ROOT, "derived/gwas"); OUT <- file.path(PROJ_ROOT, "results/mr"); dir.create(file.path(OUT, "plots"), recursive = TRUE, showWarnings = FALSE)
PLINK <- file.path(PROJ_ROOT, "data/raw/tools/plink/plink.exe"); BFILE <- file.path(PROJ_ROOT, "data/raw/ref/g1000_eur/g1000_eur")
META <- readRDS(file.path(DER, "trait_meta.rds"))
cache <- new.env()
load_trait <- function(k) { if (is.null(cache[[k]])) cache[[k]] <- readRDS(file.path(DER, paste0(k, ".rds"))); cache[[k]] }

# Regions (GRCh38, +/-500 kb) for sensitivity analyses when AMD is the exposure
AMD_MAJOR <- data.table(gene = c("CFH", "ARMS2/HTRA1", "C3"), chr = c("1", "10", "19"), start = c(196651878, 122454932, 6677704), end = c(196747504, 122514908, 6720662))
LIPID <- data.table(gene = c("APOE", "CETP", "LIPC", "ABCA1", "LPL"), chr = c("19", "16", "15", "9", "8"),
                    start = c(44905796, 56961923, 58410569, 104781006, 19901717), end = c(44909393, 56983845, 58569844, 104928155, 19967259))
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
  if (anyNA(d$pos) && "AMD_FG" %in% names(META)) {        # fill missing GRCh38 positions from FinnGen by rsID
    if (is.null(POSMAP)) POSMAP <<- load_trait("AMD_FG")[, .(SNP, chr38 = chr, pos38 = pos)]
    d <- merge(d, POSMAP, by = "SNP", all.x = TRUE)
    d[is.na(pos) & !is.na(pos38), `:=`(chr = chr38, pos = pos38)][, c("chr38", "pos38") := NULL]
  }
  d <- clump_local(d)
  d[, F := (beta / se)^2]
  d[F > 10]
}

to_exp <- function(d, k) {
  x <- as.data.frame(d); x$Phenotype <- k
  format_data(x, type = "exposure", snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "eaf",
              effect_allele_col = "ea", other_allele_col = "oa", pval_col = "p", samplesize_col = "n",
              ncase_col = if ("ncase" %in% names(x)) "ncase" else "ncase_", ncontrol_col = if ("ncontrol" %in% names(x)) "ncontrol" else "ncontrol_",
              chr_col = "chr", pos_col = "pos", phenotype_col = "Phenotype")
}
to_out <- function(k, snps) {
  d <- load_trait(k)[SNP %in% snps]; x <- as.data.frame(d); x$Phenotype <- k
  format_data(x, type = "outcome", snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "eaf",
              effect_allele_col = "ea", other_allele_col = "oa", pval_col = "p", samplesize_col = "n",
              ncase_col = if ("ncase" %in% names(x)) "ncase" else "ncase_", ncontrol_col = if ("ncontrol" %in% names(x)) "ncontrol" else "ncontrol_",
              chr_col = "chr", pos_col = "pos", phenotype_col = "Phenotype")
}

PREV <- c(AMD_IAMDGC = 0.02, AMD_FG = 0.03, DRYAMD_FG = 0.02, WETAMD_FG = 0.015, AS = 0.03, IS = 0.025, CES = 0.005, LAS = 0.005, SVS = 0.005)
add_r2 <- function(dat) {   # variance explained for Steiger (binary: liability scale via get_r_from_lor)
  for (side in c("exposure", "outcome")) {
    k <- dat[[side]][1]; b <- dat[[paste0("beta.", side)]]; eaf <- dat[[paste0("eaf.", side)]]
    if (all(is.na(eaf))) eaf <- dat[[paste0("eaf.", ifelse(side == "exposure", "outcome", "exposure"))]]
    eaf[is.na(eaf)] <- 0.3
    if (k %in% names(PREV)) {
      r <- get_r_from_lor(b, eaf, dat[[paste0("ncase.", side)]][1], dat[[paste0("ncontrol.", side)]][1], PREV[[k]])
    } else {
      r <- get_r_from_bsen(b, dat[[paste0("se.", side)]], dat[[paste0("samplesize.", side)]])
    }
    dat[[paste0("rsq.", side)]] <- r^2; dat[[paste0("units.", side)]] <- ifelse(k %in% names(PREV), "log odds", "SD")
  }
  dat
}

run_pair <- function(ek, ok, iv, tag = "main", presso = TRUE) {
  if (nrow(iv) < 1) return(NULL)
  exp <- to_exp(iv, ek); out <- to_out(ok, iv$SNP)
  if (!nrow(out)) return(NULL)
  dat <- harmonise_data(exp, out, action = 2); dat <- dat[dat$mr_keep, ]
  if (nrow(dat) < 1) return(NULL)
  dat <- add_r2(dat)
  st <- tryCatch(directionality_test(dat), error = function(e) NULL)
  meth <- if (nrow(dat) >= 3) c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode")
          else if (nrow(dat) == 2) c("mr_ivw") else c("mr_wald_ratio")
  res <- mr(dat, method_list = meth)
  het <- if (nrow(dat) >= 3) mr_heterogeneity(dat, method_list = "mr_ivw") else NULL
  ple <- if (nrow(dat) >= 3) mr_pleiotropy_test(dat) else NULL
  pr <- NULL
  if (presso && nrow(dat) >= 4) {
    pr <- tryCatch({
      p <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure",
                     OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = max(1000, 20 * nrow(dat)), SignifThreshold = 0.05)
      mr <- p$`Main MR results`
      data.frame(presso_global_p = p$`MR-PRESSO results`$`Global Test`$Pvalue,
                 presso_n_outlier = length(p$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`),
                 presso_corr_b = mr$`Causal Estimate`[2], presso_corr_se = mr$Sd[2], presso_corr_p = mr$`P-value`[2])
    }, error = function(e) NULL)
  }
  loo <- if (nrow(dat) >= 3) mr_leaveoneout(dat) else NULL
  res$exposure <- ek; res$outcome <- ok; res$analysis <- tag
  res$OR <- exp(res$b); res$lo <- exp(res$b - 1.96 * res$se); res$hi <- exp(res$b + 1.96 * res$se)
  summ <- data.frame(exposure = ek, outcome = ok, analysis = tag, nsnp = nrow(dat),
                     F_mean = mean((dat$beta.exposure / dat$se.exposure)^2), F_min = min((dat$beta.exposure / dat$se.exposure)^2),
                     Q = if (!is.null(het)) het$Q else NA, Q_p = if (!is.null(het)) het$Q_pval else NA,
                     egger_int = if (!is.null(ple)) ple$egger_intercept else NA, egger_int_p = if (!is.null(ple)) ple$pval else NA,
                     steiger_correct = if (!is.null(st)) st$correct_causal_direction else NA,
                     steiger_p = if (!is.null(st)) st$steiger_pval else NA)
  if (!is.null(pr)) summ <- cbind(summ, pr)
  list(res = res, summ = summ, dat = dat, loo = loo)
}

AMD <- c("AMD_IAMDGC", "AMD_FG", "DRYAMD_FG", "WETAMD_FG")
STROKE <- c("AS", "IS", "LAS", "CES", "SVS")
avail <- sub("\\.rds$", "", list.files(DER, "\\.rds$")); avail <- setdiff(avail, "trait_meta")
AMD <- intersect(AMD, avail); STROKE <- intersect(STROKE, avail); has_wmh <- "WMH" %in% avail
message("available: ", paste(avail, collapse = ", "))

ivs <- list(); iv_log <- list()
for (k in c(AMD, STROKE, if (has_wmh) "WMH")) {
  iv <- get_iv(k, 5e-8); thr <- 5e-8
  if (nrow(iv) < 3) { iv <- get_iv(k, 5e-6); thr <- 5e-6 }
  ivs[[k]] <- iv; iv_log[[k]] <- data.frame(trait = k, p_threshold = thr, n_iv = nrow(iv), F_mean = mean(iv$F), F_min = min(iv$F))
  fwrite(iv, file.path(OUT, paste0("instruments_", k, ".tsv")), sep = "\t")
  message(sprintf("IV %s: %d SNPs (p<%g)", k, nrow(iv), thr))
}
fwrite(rbindlist(iv_log), file.path(OUT, "instrument_summary.tsv"), sep = "\t")

pairs <- rbind(expand.grid(e = AMD, o = c(STROKE, if (has_wmh) "WMH"), dir = "AMD->stroke", stringsAsFactors = FALSE),
               expand.grid(e = c(STROKE, if (has_wmh) "WMH"), o = AMD, dir = "stroke->AMD", stringsAsFactors = FALSE))
ALL_RES <- list(); ALL_SUM <- list(); ALL_LOO <- list(); ALL_DAT <- list()
for (i in seq_len(nrow(pairs))) {
  ek <- pairs$e[i]; ok <- pairs$o[i]
  message(sprintf("[%d/%d] %s -> %s", i, nrow(pairs), ek, ok))
  r <- run_pair(ek, ok, ivs[[ek]], "main")
  if (is.null(r)) next
  key <- paste(ek, ok, sep = "__")
  ALL_RES[[key]] <- r$res; ALL_SUM[[key]] <- r$summ; ALL_DAT[[key]] <- r$dat
  if (!is.null(r$loo)) { r$loo$exposure <- ek; r$loo$outcome <- ok; ALL_LOO[[key]] <- r$loo }
  if (ek %in% AMD) {   # sensitivity: drop major AMD loci; drop lipid loci
    for (sv in list(list("excl_major_AMD_loci", AMD_MAJOR), list("excl_lipid_loci", LIPID), list("excl_major_and_lipid", rbind(AMD_MAJOR, LIPID)))) {
      iv2 <- ivs[[ek]][!in_regions(ivs[[ek]], sv[[2]])]
      r2 <- run_pair(ek, ok, iv2, sv[[1]], presso = FALSE)
      if (!is.null(r2)) { ALL_RES[[paste(key, sv[[1]])]] <- r2$res; ALL_SUM[[paste(key, sv[[1]])]] <- r2$summ }
    }
  }
}
res <- rbindlist(ALL_RES, fill = TRUE); summ <- rbindlist(ALL_SUM, fill = TRUE)
ivw <- res[method %in% c("Inverse variance weighted", "Wald ratio") & analysis == "main"]
ivw[, p_fdr := p.adjust(pval, "BH")]; ivw[, p_bonf := pmin(1, pval * .N)]
fwrite(res, file.path(OUT, "mr_all_methods.tsv"), sep = "\t")
fwrite(summ, file.path(OUT, "mr_diagnostics.tsv"), sep = "\t")
fwrite(ivw, file.path(OUT, "mr_ivw_main.tsv"), sep = "\t")
if (length(ALL_LOO)) fwrite(rbindlist(ALL_LOO, fill = TRUE), file.path(OUT, "mr_leaveoneout.tsv"), sep = "\t")
saveRDS(ALL_DAT, file.path(OUT, "harmonised_data.rds"))

# Plots for pairs with nominal IVW significance
sig <- ivw[pval < 0.05]
for (j in seq_len(nrow(sig))) {
  key <- paste(sig$exposure[j], sig$outcome[j], sep = "__"); dat <- ALL_DAT[[key]]; rr <- ALL_RES[[key]]
  if (is.null(dat) || nrow(dat) < 3) next
  ggsave(file.path(OUT, "plots", paste0("scatter_", key, ".pdf")), mr_scatter_plot(rr, dat)[[1]], width = 6, height = 5)
  ggsave(file.path(OUT, "plots", paste0("funnel_", key, ".pdf")), mr_funnel_plot(mr_singlesnp(dat))[[1]], width = 6, height = 5)
  ggsave(file.path(OUT, "plots", paste0("loo_", key, ".pdf")), mr_leaveoneout_plot(mr_leaveoneout(dat))[[1]], width = 6, height = max(4, nrow(dat) * 0.12))
}
print(ivw[, .(exposure, outcome, nsnp, OR = round(OR, 3), lo = round(lo, 3), hi = round(hi, 3), pval = signif(pval, 3), p_fdr = signif(p_fdr, 3))])
message("MR done")
