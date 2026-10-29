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
OUT <- file.path(PROJ_ROOT, "results/revision/mr"); dir.create(file.path(OUT, "plots"), recursive = TRUE, showWarnings = FALSE)
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
fwrite(inst_tab, file.path(OUT, "instruments_all_schemes.tsv"), sep = "\t")
fwrite(dcast(inst_tab[, .N, by = .(exposure, scheme, locus)], exposure + scheme ~ locus, value.var = "N", fill = 0),
       file.path(OUT, "instruments_per_locus.tsv"), sep = "\t")

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
fwrite(attr2, file.path(OUT, "instrument_attrition.tsv"), sep = "\t")

# ------------------------------------------------------------------ 4. LD matrix for scheme B (signed r relative to bim A1)
writeLines(b52$SNP, file.path(TMP, "b52.txt"))
system2(PLINK, c("--bfile", BFILE, "--extract", file.path(TMP, "b52.txt"), "--r", "square", "spaces", "--write-snplist", "--out", file.path(TMP, "b52ld")), stdout = FALSE, stderr = FALSE)
b52_snps <- readLines(file.path(TMP, "b52ld.snplist")); R52 <- as.matrix(read.table(file.path(TMP, "b52ld.ld"))); dimnames(R52) <- list(b52_snps, b52_snps)
bim52 <- fread(paste0(BFILE, ".bim"), header = FALSE)[V2 %in% b52_snps, .(SNP = V2, A1 = V5, A2 = V6)]

# ------------------------------------------------------------------ 5. MR engine
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
ivw_mre <- function(m) {
  w <- m$bx^2 / m$sy^2; b <- sum(w * m$by / m$bx) / sum(w); se <- 1 / sqrt(sum(w))
  if (nrow(m) >= 3) { phi <- sum(w * (m$by / m$bx - b)^2) / (nrow(m) - 1); se <- se * sqrt(max(1, phi)) }
  c(b = b, se = se, p = 2 * pnorm(-abs(b / se)), n = nrow(m))
}
run_methods <- function(m) {
  out <- list(IVW = ivw_mre(m))
  if (nrow(m) >= 3) {
    x <- data.frame(beta.exposure = m$bx, se.exposure = m$sx, beta.outcome = m$by, se.outcome = m$sy)
    eg <- tryCatch(mr_egger_regression(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome, default_parameters()), error = function(e) NULL)
    wm <- tryCatch(mr_weighted_median(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome, list(nboot = 1000)), error = function(e) NULL)
    wmo <- tryCatch(mr_weighted_mode(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome, default_parameters()), error = function(e) NULL)
    ra <- tryCatch({ r <- mr.raps::mr.raps.overdispersed.robust(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome, loss.function = "huber", suppress.warning = TRUE)
                     list(b = r$beta.hat, se = r$beta.se, pval = 2 * pnorm(-abs(r$beta.hat / r$beta.se))) }, error = function(e) NULL)
    for (nm in c("Egger", "WeightedMedian", "WeightedMode", "RAPS")) {
      r <- switch(nm, Egger = eg, WeightedMedian = wm, WeightedMode = wmo, RAPS = ra)
      if (!is.null(r)) out[[nm]] <- c(b = r$b, se = r$se, p = r$pval, n = nrow(m))
    }
    if (!is.null(eg)) out$Egger_intercept <- c(b = eg$b_i, se = eg$se_i, p = eg$pval_i, n = nrow(m))
  }
  out
}
corr_ivw <- function(m) {     # correlated-instrument IVW (scheme B) with signed LD aligned to exposure effect alleles
  s <- intersect(m$SNP, rownames(R52)); m <- m[match(s, SNP)]
  if (nrow(m) < 3) return(NULL)
  a1 <- bim52[match(m$SNP, SNP)]; sgn <- ifelse(m$ea == a1$A1, 1, ifelse(m$ea == a1$A2, -1, NA)); sgn[is.na(sgn)] <- 1
  rho <- R52[s, s] * outer(sgn, sgn)
  inp <- mr_input(bx = m$bx, bxse = m$sx, by = m$by, byse = m$sy, correlation = rho, snps = m$SNP)
  r <- tryCatch(mr_ivw(inp, correl = TRUE, model = "random"), error = function(e) NULL); if (is.null(r)) return(NULL)
  eg <- tryCatch(mr_egger(inp, correl = TRUE), error = function(e) NULL)
  list(IVW_correl = c(b = r@Estimate, se = r@StdError, p = r@Pvalue, n = nrow(m)),
       Egger_correl = if (!is.null(eg)) c(b = eg@Estimate, se = eg@StdError.Est, p = eg@Causal.pval, n = nrow(m)) else NULL)
}

RES <- list(); CONTR <- list(); k <- 0
for (e in names(INST)) for (s in names(INST[[e]])) {
  ex <- INST[[e]][[s]]; if (!nrow(ex)) next
  ex[, locus := locus_of(chr, pos)]
  for (o in OUTC) {
    m0 <- harm(ex, OD[[o]]); m0 <- merge(m0, ex[, .(SNP, locus)], by = "SNP")
    for (xn in names(EXCL)) {
      m <- m0[!locus %in% EXCL[[xn]]]; if (nrow(m) < 2) next
      if (s == "B_52") { r <- corr_ivw(m); if (is.null(r)) next; rr <- r } else rr <- run_methods(m)
      k <- k + 1
      RES[[k]] <- rbindlist(lapply(names(rr), function(nm) if (!is.null(rr[[nm]])) data.table(exposure = e, scheme = s, outcome = o, exclusion = xn,
                     method = nm, b = rr[[nm]][["b"]], se = rr[[nm]][["se"]], p = rr[[nm]][["p"]], nsnp = rr[[nm]][["n"]],
                     n_removed = nrow(m0) - nrow(m))))
      if (xn == "main" && s != "B_52") {   # locus contributions to the IVW estimate: sum_j w_j r_j / sum w
        w <- m$bx^2 / m$sy^2; m[, contrib := w * (by / bx) / sum(w)]
        CONTR[[length(CONTR) + 1]] <- m[, .(contrib = sum(contrib), n = .N, weight_share = sum(bx^2 / sy^2) / sum(w)), by = locus][, `:=`(exposure = e, scheme = s, outcome = o)]
      }
    }
  }
}
res <- rbindlist(RES)
res[, binary := outcome %in% BINARY]
res[, `:=`(OR = ifelse(binary, exp(b), NA), lo = ifelse(binary, exp(b - 1.96 * se), NA), hi = ifelse(binary, exp(b + 1.96 * se), NA))]
# power (minimum detectable effect at 80% power, alpha 0.05, from the IVW SE) and TOST equivalence
EQ_BIN <- log(1.05); EQ_Q <- 0.05
res[method %in% c("IVW", "IVW_correl"), `:=`(
  mde_logOR = (qnorm(0.975) + qnorm(0.8)) * se,
  tost_p = { bnd <- ifelse(binary, EQ_BIN, EQ_Q); pmax(1 - pnorm((b + bnd) / se), pnorm((b - bnd) / se)) })]
res[method %in% c("IVW", "IVW_correl"), `:=`(mde_OR_up = exp(mde_logOR), mde_OR_down = exp(-mde_logOR),
                                             equivalent = tost_p < 0.05)]
fwrite(res, file.path(OUT, "mr_revised_all.tsv"), sep = "\t")
contr <- rbindlist(CONTR); fwrite(contr, file.path(OUT, "mr_locus_contributions.tsv"), sep = "\t")
fwrite(attr_tab, file.path(OUT, "outcome_proxy_summary.tsv"), sep = "\t")

# ------------------------------------------------------------------ 6. key tables and figures
key <- res[method %in% c("IVW", "IVW_correl") & outcome %in% c("AS", "IS", "LAS", "SVS", "CES", "WMH", "MEGA_AS", "MEGA_IS", "MEGA_SVS", "MEGA_LAS")]
wide <- dcast(key[exclusion %in% c("main", "no_10q26", "no_ADAMTS9", "no_CFH", "no_C3", "no_APOE", "no_major3", "no_all")],
              exposure + scheme + outcome ~ exclusion, value.var = c("b", "p"))
fwrite(wide, file.path(OUT, "mr_leave_one_locus_out_wide.tsv"), sep = "\t")
pd <- key[scheme == "A_dist" & outcome %in% c("AS", "IS", "LAS", "SVS", "MEGA_IS", "MEGA_SVS") & binary]
pd[, exclusion := factor(exclusion, levels = names(EXCL))]
g <- ggplot(pd, aes(exp(b), exclusion, colour = p < 0.05)) + geom_vline(xintercept = 1, lty = 2) +
  geom_pointrange(aes(xmin = exp(b - 1.96 * se), xmax = exp(b + 1.96 * se)), size = 0.25) +
  facet_grid(exposure ~ outcome) + scale_colour_manual(values = c(`FALSE` = "grey40", `TRUE` = "#d4502a")) +
  labs(x = "OR per unit log-odds of AMD liability (distance-pruned instruments)", y = NULL, colour = "P < 0.05") + theme_bw(8)
ggsave(file.path(OUT, "plots", "Fig_LOLO_forest.pdf"), g, width = 12, height = 7); ggsave(file.path(OUT, "plots", "Fig_LOLO_forest.png"), g, width = 12, height = 7, dpi = 200)
cp <- contr[scheme == "A_dist" & outcome %in% c("AS", "IS", "LAS", "SVS")]
g <- ggplot(cp, aes(outcome, contrib, fill = locus)) + geom_col() + geom_hline(yintercept = 0) + facet_wrap(~exposure, nrow = 1) +
  labs(y = "contribution to IVW log-OR (sum = IVW estimate)", x = NULL) + theme_bw(9)
ggsave(file.path(OUT, "plots", "Fig_locus_contributions.pdf"), g, width = 11, height = 4); ggsave(file.path(OUT, "plots", "Fig_locus_contributions.png"), g, width = 11, height = 4, dpi = 200)
print(key[scheme == "A_dist" & exclusion %in% c("main", "no_10q26", "no_ADAMTS9", "no_CFH", "no_major3", "no_all") & outcome %in% c("IS", "SVS", "AS"),
          .(exposure, outcome, exclusion, nsnp, OR = round(exp(b), 3), p = signif(p, 2), mde_OR_down = round(mde_OR_down, 3), equivalent)])
message("rev_mr_01 done")
