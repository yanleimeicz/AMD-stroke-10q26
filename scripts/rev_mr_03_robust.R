# Revision C, step 3: MR-Clust, MRlap (sample-overlap correction), CAUSE (if installed), and LAVA sensitivity at 10q26
# (FinnGen-free MEGASTROKE outcomes, corrected IAMDGC, IAMDGC overlap scenarios, merged/extended blocks).
# Usage: Rscript rev_mr_03_robust.R [clust] [mrlap] [cause] [lava]
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
setDTthreads(30); set.seed(20261007)
DER <- file.path(PROJ_ROOT, "derived/gwas"); DRV <- file.path(PROJ_ROOT, "derived/revision"); R <- file.path(PROJ_ROOT, "results")
OUT <- file.path(R, "revision", "mr"); TMP <- file.path(PROJ_ROOT, "tmp/revmr")
path_of <- function(k) { f <- file.path(DRV, paste0(k, ".rds")); if (file.exists(f)) f else file.path(DER, paste0(k, ".rds")) }
STEP <- commandArgs(TRUE); if (!length(STEP)) STEP <- c("clust", "mrlap", "cause", "lava")
LD <- file.path(PROJ_ROOT, "data/raw/ref/eur_w_ld_chr"); HM3 <- file.path(PROJ_ROOT, "derived/ldsc/w_hm3.snplist")

# ------------------------------------------------------------------ MR-Clust on harmonised instruments
if ("clust" %in% STEP && requireNamespace("mrclust", quietly = TRUE)) {
  inst <- fread(file.path(OUT, "instruments_all_schemes.tsv"))
  CL <- list()
  for (e in c("AMD_FG", "WETAMD_FG", "AMD_IAMDGC")) for (s in c("A_dist", "C_orig")) {
    ex <- inst[exposure == e & scheme == s]
    for (o in c("IS", "SVS", "AS", "LAS")) {
      od <- readRDS(path_of(o))[SNP %in% ex$SNP, .(SNP, eo = ea, oo = oa, by = beta, sy = se)]
      m <- merge(ex, od, by = "SNP"); m[, fl := fifelse(eo == ea & oo == oa, 1, fifelse(eo == oa & oo == ea, -1, NA_real_))]
      m <- m[!is.na(fl)][, by := by * fl]; if (nrow(m) < 6) next
      th <- m$by / m$beta; thse <- abs(m$sy / m$beta)
      r <- tryCatch(mrclust::mr_clust_em(theta = th, theta_se = thse, bx = m$beta, by = m$by, bxse = m$se, byse = m$sy, obs_names = m$SNP),
                    error = function(e) NULL)
      if (is.null(r)) next
      a <- as.data.table(r$results$best)
      a <- merge(a, m[, .(observation = SNP, locus)], by = "observation", all.x = TRUE)
      CL[[paste(e, s, o)]] <- a[, `:=`(exposure = e, scheme = s, outcome = o)]
    }
  }
  cl <- rbindlist(CL, fill = TRUE); fwrite(cl, file.path(OUT, "mrclust_assignments.tsv"), sep = "\t")
  sm <- cl[probability >= 0.8, .(n = .N, loci = paste(sort(unique(locus)), collapse = ","), mean_theta = mean(theta)), by = .(exposure, scheme, outcome, cluster_class)]
  fwrite(sm, file.path(OUT, "mrclust_summary.tsv"), sep = "\t"); print(sm)
}

# ------------------------------------------------------------------ MRlap (LDSC-based correction for sample overlap, weak instruments, winner's curse)
if ("mrlap" %in% STEP && requireNamespace("MRlap", quietly = TRUE)) {
  prep <- function(k) { d <- readRDS(path_of(k)); data.frame(rsid = d$SNP, chr = d$chr, pos = d$pos, alt = d$ea, ref = d$oa, beta = d$beta, se = d$se, N = d$n) }
  pairs <- data.table(e = c("AMD_FG", "AMD_FG", "AMD_FG", "DRYAMD_FG", "DRYAMD_FG", "WETAMD_FG"), o = c("IS", "SVS", "AS", "AS", "IS", "IS"))
  run1 <- function(i) {
    suppressPackageStartupMessages(library(MRlap))
    e <- pairs$e[i]; o <- pairs$o[i]; wd <- file.path(TMP, paste0("mrlap_", e, "_", o)); dir.create(wd, showWarnings = FALSE); setwd(wd)
    r <- tryCatch(MRlap(exposure = prep(e), exposure_name = e, outcome = prep(o), outcome_name = o, ld = LD, hm3 = HM3,
                        MR_threshold = 5e-8, MR_pruning_dist = 1000, MR_pruning_LD = 0, MR_reverse = 1e-3, save_logfiles = FALSE, verbose = FALSE),
                  error = function(err) list(error = conditionMessage(err)))
    if (!is.null(r$error)) return(data.table(exposure = e, outcome = o, error = r$error))
    data.table(exposure = e, outcome = o, observed_b = r$MRcorrection$observed_effect, observed_se = r$MRcorrection$observed_effect_se,
               observed_p = r$MRcorrection$observed_effect_p, corrected_b = r$MRcorrection$corrected_effect, corrected_se = r$MRcorrection$corrected_effect_se,
               corrected_p = r$MRcorrection$corrected_effect_p, p_difference = r$MRcorrection$p_difference, m_IVs = r$MRcorrection$m_IVs,
               int_crosstrait = r$LDSC$int_crosstrait, int_crosstrait_se = r$LDSC$int_crosstrait_se, h2_exp = r$LDSC$h2_exp, h2_out = r$LDSC$h2_out, rg = r$LDSC$rg)
  }
  cl <- makeCluster(nrow(pairs)); clusterExport(cl, c("pairs", "prep", "path_of", "DRV", "DER", "TMP", "LD", "HM3"))
  clusterEvalQ(cl, { .libPaths(c(file.path(Sys.getenv("COMORBIDITY_ROOT"), "Rlib"), .libPaths())); library(data.table) })
  res <- rbindlist(parLapply(cl, seq_len(nrow(pairs)), run1), fill = TRUE); stopCluster(cl)
  fwrite(res, file.path(OUT, "mrlap_overlap_corrected.tsv"), sep = "\t"); print(res)
}

# ------------------------------------------------------------------ CAUSE (correlated pleiotropy), AMD_FG -> IS / SVS
if ("cause" %in% STEP && requireNamespace("cause", quietly = TRUE)) {
  suppressPackageStartupMessages({ library(cause); library(loo); library(dplyr) })
  fx <- new.env(parent = asNamespace("cause")); sys.source(file.path(PROJ_ROOT, "scripts/helpers/cause_fix.R"), envir = fx)
  assignInNamespace("in_sample_elpd_loo", fx$in_sample_elpd_loo, ns = "cause")
  CA <- list()
  for (o in c("IS", "SVS")) {
    x1 <- readRDS(path_of("AMD_FG")); x2 <- readRDS(path_of(o))
    X <- gwas_merge(as.data.frame(x1[, .(snp = SNP, beta_hat = beta, se = se, A1 = ea, A2 = oa)]), as.data.frame(x2[, .(snp = SNP, beta_hat = beta, se = se, A1 = ea, A2 = oa)]),
                    snp_name_cols = c("snp", "snp"), beta_hat_cols = c("beta_hat", "beta_hat"), se_cols = c("se", "se"), A1_cols = c("A1", "A1"), A2_cols = c("A2", "A2"))
    rm(x1, x2); gc()
    vs <- sample(X$snp, min(1e6, nrow(X)))
    params <- est_cause_params(X, vs)
    # LD pruning of candidate variants (P < 1e-3 in AMD) with local PLINK, r2 < 0.01 / 1 Mb, prioritising AMD P values
    X$p1 <- 2 * pnorm(-abs(X$beta_hat_1 / X$seb1)); cand <- X[X$p1 < 1e-3, ]
    fwrite(data.table(SNP = cand$snp, P = pmax(cand$p1, 1e-300)), file.path(TMP, paste0("cause_cand_", o, ".txt")), sep = "\t")
    system2(file.path(PROJ_ROOT, "data/raw/tools/plink/plink.exe"), c("--bfile", file.path(PROJ_ROOT, "data/raw/ref/g1000_eur/g1000_eur"), "--clump", file.path(TMP, paste0("cause_cand_", o, ".txt")),
            "--clump-p1", "1e-3", "--clump-p2", "1", "--clump-r2", "0.01", "--clump-kb", "1000", "--threads", "30", "--out", file.path(TMP, paste0("cause_clump_", o))), stdout = FALSE)
    top <- fread(file.path(TMP, paste0("cause_clump_", o, ".clumped")))$SNP
    res <- cause(X = X, variants = top, param_ests = params)
    s <- summary(res)
    CA[[o]] <- data.table(exposure = "AMD_FG", outcome = o, n_variants = length(top), p_causal_vs_sharing = s$p,
                          q_causal = paste(signif(s$quants[[2]][, "q"], 3), collapse = ";"), gamma_causal = paste(signif(s$quants[[2]][, "gamma"], 3), collapse = ";"),
                          eta_sharing = paste(signif(s$quants[[1]][, "eta"], 3), collapse = ";"))
    saveRDS(res, file.path(TMP, paste0("cause_", o, ".rds")))
  }
  ca <- rbindlist(CA); fwrite(ca, file.path(OUT, "cause_results.tsv"), sep = "\t"); print(ca)
}

# ------------------------------------------------------------------ LAVA sensitivity at 10q26 (chr10)
if ("lava" %in% STEP) {
  suppressPackageStartupMessages(library(LAVA))
  W <- file.path(DRV, "lava10"); dir.create(W, showWarnings = FALSE, recursive = TRUE)
  REF <- file.path(PROJ_ROOT, "derived/lava/ref/g1000_eur_chr10")
  bim <- fread(paste0(REF, ".bim"), header = FALSE, select = c(1, 2, 4), col.names = c("CHR", "SNP", "BP"))[BP > 120e6 & BP < 127e6]
  PH <- c("AMD_FG", "DRYAMD_FG", "WETAMD_FG", "AMD_IAMDGC", "AMD_IAMDGC_fix", "IS", "SVS", "AS", "MEGA_IS", "MEGA_SVS", "MEGA_AS")
  PH <- PH[file.exists(vapply(PH, path_of, ""))]
  info <- rbindlist(lapply(PH, function(k) {
    d <- readRDS(path_of(k))[SNP %in% bim$SNP, .(SNP, A1 = ea, A2 = oa, Z = beta / se, N = n)]
    m <- readRDS(path_of(k))[1]
    f <- file.path(W, paste0(k, "_chr10.txt")); fwrite(d, f, sep = "\t")
    data.table(phenotype = k, cases = if ("ncase" %in% names(m)) m$ncase else NA, controls = if ("ncontrol" %in% names(m)) m$ncontrol else NA, filename = f)
  }))
  fwrite(info, file.path(W, "info.txt"), sep = "\t", na = "NA")
  base <- as.matrix(read.table(file.path(PROJ_ROOT, "derived/lava/sample_overlap.txt"), check.names = FALSE))
  ov_ldsc <- fread(file.path(R, "csvd", "ldsc_intercept_matrix.tsv"))
  ldsc_val <- function(a, b) { if (!a %in% ov_ldsc$trait || !b %in% names(ov_ldsc)) return(NA_real_); ov_ldsc[trait == a][[b]] }
  mk_ov <- function(iamdgc_val) {     # pairwise use only; AMD vs MEGASTROKE = 0 (no FinnGen in MEGASTROKE; IAMDGC has no UKB)
    M <- diag(length(PH)); dimnames(M) <- list(PH, PH)
    for (a in PH) for (b in PH) if (a != b) {
      aa <- sub("_fix$", "", a); bb <- sub("_fix$", "", b); bb2 <- if (bb == "AS") "IS" else bb; aa2 <- if (aa == "AS") "IS" else aa
      v <- if (aa2 %in% rownames(base) && bb2 %in% colnames(base)) base[aa2, bb2] else 0
      if (grepl("^MEGA_", a) || grepl("^MEGA_", b)) v <- 0
      if (grepl("AMD_IAMDGC", a) || grepl("AMD_IAMDGC", b)) {
        other <- if (grepl("AMD_IAMDGC", a)) b else a
        if (other %in% c("IS", "SVS", "AS")) v <- if (is.character(iamdgc_val)) { x <- ldsc_val("AMD_IAMDGC", if (other == "AS") "IS" else other); ifelse(is.na(x), 0, x) } else iamdgc_val
        if (grepl("^AMD|^DRY|^WET", other)) v <- 0
      }
      M[a, b] <- v
    }
    M
  }
  LOC <- data.table(LOC = c("1596", "1597", "1596_1597", "1595_1598"), CHR = 10,
                    START = c(122803507, 123856185, 122803507, 121874446), STOP = c(123856184, 124894142, 124894142, 125798319))
  OUTS <- intersect(c("IS", "SVS", "AS", "MEGA_IS", "MEGA_SVS", "MEGA_AS"), PH)
  AMDS <- intersect(c("AMD_FG", "DRYAMD_FG", "WETAMD_FG", "AMD_IAMDGC", "AMD_IAMDGC_fix"), PH)
  RES <- list()
  for (sc in list(list("ldsc_intercept", "ldsc"), list("iamdgc_0", 0), list("iamdgc_0.1", 0.1), list("iamdgc_0.2", 0.2))) {
    M <- mk_ov(sc[[2]]); fo <- file.path(W, paste0("overlap_", sc[[1]], ".txt")); write.table(round(M, 5), fo, quote = FALSE)
    for (a in AMDS) for (o in OUTS) {
      if (sc[[1]] != "ldsc_intercept" && !grepl("AMD_IAMDGC", a)) next          # overlap scenarios only matter for IAMDGC pairs
      input <- tryCatch(process.input(input.info.file = file.path(W, "info.txt"), sample.overlap.file = fo, ref.prefix = REF, phenos = c(a, o)),
                        error = function(e) NULL); if (is.null(input)) next          # pairwise input: SNPs shared by the two traits only
      for (i in seq_len(nrow(LOC))) {
        locus <- tryCatch(process.locus(LOC[i], input), error = function(e) NULL); if (is.null(locus)) next
        u <- as.data.table(run.univ(locus))
        b <- tryCatch(as.data.table(run.bivar(locus, phenos = c(a, o))), error = function(e) NULL); if (is.null(b) || !nrow(b)) next
        RES[[length(RES) + 1]] <- b[, `:=`(scenario = sc[[1]], locus = LOC$LOC[i], n_snps = locus$n.snps,
                                           h2_p_amd = u[phen == a, p], h2_p_out = u[phen == o, p])]
      }
    }
  }
  lv <- rbindlist(RES, fill = TRUE); fwrite(lv, file.path(OUT, "lava_10q26_sensitivity.tsv"), sep = "\t")
  print(lv[, .(scenario, locus, phen1, phen2, rho = round(rho, 3), lo = round(rho.lower, 2), hi = round(rho.upper, 2), p = signif(p, 3), h2_p_out = signif(h2_p_out, 2))])
}
message("rev_mr_03 done: ", paste(STEP, collapse = ","))
