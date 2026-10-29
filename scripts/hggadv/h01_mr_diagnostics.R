# Revision C, step 1: rebuilt instruments, leave-one-locus-out MR, locus contributions, power and equivalence.
# Schemes: A_dist  = one lead per +/-1 Mb (P < 5e-8) in each exposure;
#          B_52    = the 52 conditionally independent IAMDGC variants (Fritsche 2016; GWAS Catalog GCST003219), correlated IVW;
#          C_orig  = original 1000G r2 < 0.001 / 10 Mb clumped instruments (results/mr/instruments_*.tsv), for comparison.
# Outcome SNPs missing are replaced by the best 1000G EUR proxy (r2 > 0.8, +/-500 kb; allele mapping from phased LD).
# Locus exclusion windows are +/-2 Mb (GRCh38). Outcomes: GIGASTROKE EUR (contains FinnGen), BIG40 WMH, MEGASTROKE EUR (no FinnGen).
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(TwoSampleMR); library(MendelianRandomization); library(ggplot2) })
setDTthreads(30); set.seed(20261007)
DER <- file.path(PROJ_ROOT, "derived/gwas"); DRV <- file.path(PROJ_ROOT, "derived/revision")
OUT <- file.path(PROJ_ROOT, "results/revision/mr"); OUT2 <- file.path(PROJ_ROOT, "results/hggadv"); dir.create(OUT2, showWarnings = FALSE)
PLINK <- file.path(PROJ_ROOT, "data/raw/tools/plink/plink.exe"); BFILE <- file.path(PROJ_ROOT, "data/raw/ref/g1000_eur/g1000_eur")
TMP <- file.path(PROJ_ROOT, "tmp/revmr"); dir.create(TMP, showWarnings = FALSE, recursive = TRUE)
path_of <- function(k) { f <- file.path(DRV, paste0(k, ".rds")); if (file.exists(f)) f else file.path(DER, paste0(k, ".rds")) }

EXPO <- c(AMD_FG = "AMD_FG", DRYAMD_FG = "DRYAMD_FG", WETAMD_FG = "WETAMD_FG", AMD_IAMDGC = "AMD_IAMDGC_fix")
OUTC <- c("AS", "IS", "LAS", "CES", "SVS", "WMH", "MEGA_AS", "MEGA_IS", "MEGA_LAS", "MEGA_SVS", "MEGA_CES")
OUTC <- OUTC[file.exists(vapply(OUTC, path_of, ""))]
BINARY <- setdiff(OUTC, "WMH")

LOCI <- data.table(locus = c("CFH", "10q26", "C3", "APOE", "ADAMTS9", "CETP", "LIPC", "ABCA1", "LPL"),
                   chr = c("1", "10", "19", "19", "3", "16", "15", "9", "8"),
                   start = c(196651878, 122454932, 6677704, 44905796, 64515654, 56961923, 58410569, 104781006, 19901717),
                   end = c(196747504, 122514908, 6720662, 44909393, 64687889, 56983845, 58569844, 104928155, 19967259))
FLANK <- 2e6
locus_of <- function(chr, pos) {
  lab <- rep("other", length(pos))
  for (i in seq_len(nrow(LOCI))) lab[!is.na(pos) & chr == LOCI$chr[i] & pos >= LOCI$start[i] - FLANK & pos <= LOCI$end[i] + FLANK] <- LOCI$locus[i]
  lab
}
LIPID <- c("CETP", "LIPC", "ABCA1", "LPL"); MAJOR3 <- c("CFH", "10q26", "C3")
EXCL <- list(main = character(0), no_10q26 = "10q26", no_ADAMTS9 = "ADAMTS9", no_CFH = "CFH", no_C3 = "C3", no_APOE = "APOE",
             no_major3 = MAJOR3, no_major3_APOE = c(MAJOR3, "APOE"), no_major3_lipid = c(MAJOR3, "APOE", LIPID),
             no_all = c(MAJOR3, "APOE", "ADAMTS9", LIPID))

# ------------------------------------------------------------------ 1. exposures and instruments
b52 <- fread(file.path(OUT, "input", "iamdgc_52_independent_gwascatalog.tsv"))
dist_prune <- function(d, win = 1e6) {
  d <- d[order(if ("log10p" %in% names(d)) log10p else log10(p))]; keep <- logical(nrow(d))
  taken <- data.table(chr = character(0), pos = numeric(0))
  for (i in seq_len(nrow(d))) {
    if (is.na(d$pos[i])) next
    if (!nrow(taken[chr == d$chr[i] & abs(pos - d$pos[i]) < win])) { keep[i] <- TRUE; taken <- rbind(taken, data.table(chr = d$chr[i], pos = d$pos[i])) }
  }
  d[keep]
}
EX <- list(); INST <- list()
posmap <- NULL
for (e in names(EXPO)) {
  d <- readRDS(path_of(EXPO[[e]]))
  if (is.null(posmap) && e == "AMD_FG") posmap <- d[, .(SNP, chr38 = chr, pos38 = pos)]
  if (anyNA(d$pos)) { d <- merge(d, posmap, by = "SNP", all.x = TRUE); d[is.na(pos) & !is.na(pos38), `:=`(chr = chr38, pos = pos38)]; d[, c("chr38", "pos38") := NULL] }
  sig <- if ("log10p" %in% names(d)) d[log10p < log10(5e-8)] else d[p < 5e-8]
  A <- dist_prune(sig)[, F := (beta / se)^2][F > 10]
  B <- d[SNP %in% b52$SNP][, F := (beta / se)^2]
  orig_f <- file.path(file.path(PROJ_ROOT, "results/mr"), paste0("instruments_", e, ".tsv"))
  C <- if (file.exists(orig_f)) d[SNP %in% fread(orig_f)$SNP][, F := (beta / se)^2] else d[0]
  INST[[e]] <- list(A_dist = A, B_52 = B[F > 10], C_orig = C)
  EX[[e]] <- d[SNP %in% unique(c(A$SNP, B$SNP, C$SNP))]
  message(sprintf("%s: A_dist %d, B_52 %d (F>10: %d), C_orig %d", e, nrow(A), nrow(B), sum(B$F > 10), nrow(C)))
  rm(d); gc()
}
inst_tab <- rbindlist(lapply(names(INST), function(e) rbindlist(lapply(names(INST[[e]]), function(s) {
  x <- INST[[e]][[s]]; if (!nrow(x)) return(NULL)
  x[, .(exposure = e, scheme = s, SNP, chr, pos, ea, oa, beta, se, p, F, locus = locus_of(chr, pos))] }))))

# ------------------------------------------------------------------ 2. proxies (1000G EUR, r2 > 0.8, +/-500 kb, phased)
U <- unique(inst_tab$SNP)
writeLines(U, file.path(TMP, "inst_snps.txt"))
system2(PLINK, c("--bfile", BFILE, "--r2", "in-phase", "with-freqs", "--ld-snp-list", file.path(TMP, "inst_snps.txt"),
                 "--ld-window-kb", "500", "--ld-window", "99999", "--ld-window-r2", "0.8", "--threads", "30", "--out", file.path(TMP, "proxies")), stdout = FALSE, stderr = FALSE)
PX <- fread(file.path(TMP, "proxies.ld"))[SNP_A != SNP_B]
PX[, c("hapA", "hapB") := tstrsplit(PHASE, "/", fixed = TRUE)]   # e.g. "AG/CT": A(snpA)-G(snpB) and C(snpA)-T(snpB)
PX[, `:=`(a1 = substr(hapA, 1, 1), b1 = substr(hapA, 2, 2), a2 = substr(hapB, 1, 1), b2 = substr(hapB, 2, 2))]
PX <- PX[, .(target = SNP_A, proxy = SNP_B, r2 = R2, a1, b1, a2, b2)][order(target, -r2)]

# ------------------------------------------------------------------ 3. outcomes (instrument SNPs + proxies), with proxy substitution
OD <- list(); attr_log <- list()
for (o in OUTC) {
  d <- readRDS(path_of(o))[SNP %in% c(U, PX$proxy)]
  have <- d[SNP %in% U]; miss <- setdiff(U, have$SNP); sub <- list()
  for (t in miss) {
    cand <- PX[target == t & proxy %in% d$SNP]; if (!nrow(cand)) next
    pr <- cand[1]; r <- d[SNP == pr$proxy][1]
    # map outcome effect allele (proxy allele) to the coupled target allele
    if (r$ea == pr$b1 && r$oa == pr$b2) { tea <- pr$a1; toa <- pr$a2 } else if (r$ea == pr$b2 && r$oa == pr$b1) { tea <- pr$a2; toa <- pr$a1 } else next
    sub[[t]] <- data.table(SNP = t, chr = r$chr, pos = NA_real_, ea = tea, oa = toa, beta = r$beta, se = r$se, p = r$p, eaf = NA_real_,
                           proxy = pr$proxy, proxy_r2 = pr$r2)
  }
  have[, `:=`(proxy = NA_character_, proxy_r2 = NA_real_)]
  keep <- c("SNP", "chr", "pos", "ea", "oa", "beta", "se", "p", "eaf", "proxy", "proxy_r2")
  if (!"eaf" %in% names(have)) have[, eaf := NA_real_]
  dd <- rbind(have[, ..keep], rbindlist(sub), fill = TRUE)
  m <- readRDS(path_of(o))[1]; dd[, n := m$n]; if ("ncase" %in% names(m)) dd[, `:=`(ncase = m$ncase, ncontrol = m$ncontrol)]
  OD[[o]] <- dd
  attr_log[[o]] <- data.table(outcome = o, n_instr_union = length(U), direct = nrow(have), proxied = length(sub), missing = length(miss) - length(sub))
  message(sprintf("outcome %s: %d direct, %d proxied, %d missing", o, nrow(have), length(sub), length(miss) - length(sub)))
}
attr_tab <- rbindlist(attr_log)
# per-exposure/scheme attrition
attr2 <- rbindlist(lapply(names(INST), function(e) rbindlist(lapply(names(INST[[e]]), function(s) rbindlist(lapply(OUTC, function(o) {
  snps <- INST[[e]][[s]]$SNP; x <- OD[[o]][SNP %in% snps]
  data.table(exposure = e, scheme = s, outcome = o, n_instr = length(snps), n_direct = sum(is.na(x$proxy)), n_proxy = sum(!is.na(x$proxy)), n_lost = length(snps) - nrow(x)) }))))))
harm <- function(ex, oc) {
  m <- merge(ex[, .(SNP, chr, pos, ea, oa, bx = beta, sx = se, eafx = if ("eaf" %in% names(ex)) eaf else NA_real_)],
             oc[, .(SNP, eo = ea, oo = oa, by = beta, sy = se, eafo = eaf, proxy)], by = "SNP")
  comp <- c(A = "T", T = "A", C = "G", G = "C")
  m[, flip := fifelse(eo == ea & oo == oa, 1, fifelse(eo == oa & oo == ea, -1,
              fifelse(comp[eo] == ea & comp[oo] == oa, 1, fifelse(comp[eo] == oa & comp[oo] == ea, -1, NA_real_))))]
  m[, pal := paste0(ea, oa) %in% c("AT", "TA", "CG", "GC")]
  m <- m[!is.na(flip) & !(pal & (is.na(eafx) | (eafx > 0.42 & eafx < 0.58)))]
  m[, by := by * flip]; m[]
}

# ------------------------------------------------------------------ diagnostics for the main (A_dist) MR: F, Cochran's Q, I2, Egger intercept, I2_GX, MR-PRESSO
# combinations are run in parallel (30 workers, PSOCK cluster on Windows)
library(parallel)
JOBS <- list()
for (e in names(INST)) {
  ex <- INST[[e]][["A_dist"]]; ex[, locus := locus_of(chr, pos)]
  for (o in c("AS", "IS", "LAS", "CES", "SVS", "WMH", "MEGA_IS", "MEGA_SVS")) for (xn in c("main", "no_10q26")) {
    m <- merge(harm(ex, OD[[o]]), ex[, .(SNP, locus, F)], by = "SNP")[!locus %in% EXCL[[xn]]]
    JOBS[[length(JOBS) + 1]] <- list(e = e, o = o, xn = xn, m = as.data.frame(m))
  }
}
one <- function(J) {
  suppressPackageStartupMessages({ library(data.table); library(TwoSampleMR); library(MRPRESSO) })
  m <- J$m
  w <- m$bx^2 / m$sy^2; r <- m$by / m$bx; b <- sum(w * r) / sum(w)
  Q <- sum(w * (r - b)^2); df <- nrow(m) - 1
  eg <- mr_egger_regression(m$bx, m$by, m$sx, m$sy, default_parameters())
  bxa <- abs(m$bx); wq <- 1 / m$sx^2; mu <- sum(wq * bxa) / sum(wq); Qx <- sum(wq * (bxa - mu)^2); I2GX <- max(0, (Qx - df) / Qx)
  pr <- tryCatch(mr_presso(BetaOutcome = "by", BetaExposure = "bx", SdOutcome = "sy", SdExposure = "sx", OUTLIERtest = TRUE,
                           DISTORTIONtest = FALSE, data = m, NbDistribution = 5000, SignifThreshold = 0.05, seed = 1),
                 error = function(err) NULL)
  gp <- if (!is.null(pr)) pr$`MR-PRESSO results`$`Global Test`$Pvalue else NA
  ot <- if (!is.null(pr)) pr$`MR-PRESSO results`$`Outlier Test` else NULL
  nout <- if (!is.null(ot) && is.data.frame(ot)) sum(suppressWarnings(as.numeric(sub("<", "", as.character(ot$Pvalue)))) < 0.05, na.rm = TRUE) else NA
  data.table(exposure = J$e, outcome = J$o, exclusion = J$xn, nsnp = nrow(m), F_mean = mean(m$F), F_min = min(m$F),
             Q = Q, Q_df = df, Q_p = pchisq(Q, df, lower.tail = FALSE), I2 = max(0, (Q - df) / Q),
             egger_intercept = eg$b_i, egger_intercept_se = eg$se_i, egger_intercept_p = eg$pval_i, I2GX = I2GX,
             presso_global_p = suppressWarnings(as.numeric(sub("<", "", as.character(gp)))), presso_n_outliers = nout)
}
cl <- makeCluster(30)
clusterEvalQ(cl, .libPaths(c(file.path(Sys.getenv("COMORBIDITY_ROOT"), "Rlib"), .libPaths())))
res <- parLapply(cl, JOBS, one)
stopCluster(cl)
dx <- rbindlist(res)
fwrite(dx, file.path(OUT2, "mr_diagnostics_A_dist.tsv"), sep = "	")
print(dx[outcome %in% c("IS", "SVS")], digits = 3)
