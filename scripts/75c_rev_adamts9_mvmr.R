# Revision R4 (ADAMTS9 cardiometabolic alternative), genome-wide MR step:
#   causal effects of T2D liability (Xue 2018), WHRadjBMI (Pulit 2019) and SBP (Evangelou 2018) on stroke (GIGASTROKE EUR)
#   by univariable IVW and multivariable MR (MVMR-IVW), instruments excluding the 3p14 ADAMTS9 region (+/-1 Mb),
#   then the effect of rs4611812-C on stroke that would be expected if it acted only through these risk factors.
# Output -> results/revision/adamts9/mr_cardiometabolic_to_stroke.tsv, mvmr_cardiometabolic_to_stroke.tsv,
#           expected_vs_observed_rs4611812.tsv
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressMessages({ library(data.table) })
setDTthreads(30)
RAW <- file.path(PROJ_ROOT, "data/raw/gwas"); CM <- file.path(RAW, "cardiometab")
OUT <- file.path(PROJ_ROOT, "results/revision/adamts9"); TMP <- file.path(PROJ_ROOT, "tmp/rev_adamts9")
PLINK <- file.path(PROJ_ROOT, "data/raw/tools/plink/plink.exe"); BFILE <- file.path(PROJ_ROOT, "data/raw/ref/g1000_eur/g1000_eur")
KEY <- "rs4611812"; KEY_C <- "C"
BP <- Sys.getenv("BP", "SBP"); SUF <- if (BP == "SBP") "" else paste0("_", BP)   # blood-pressure exposure: SBP (default) or DBP
bim <- fread(paste0(BFILE, ".bim"), header = FALSE, select = c(1, 2, 4, 5, 6), col.names = c("chr", "SNP", "pos37", "a1", "a2"))
bim[, chr := as.character(chr)]

# ---------------- exposures (columns: SNP ea oa beta se p eaf n)
t2d <- fread(file.path(CM, "T2D_Xue2018.h.tsv.gz"), select = c("hm_rsid", "hm_effect_allele", "hm_other_allele", "hm_beta",
                                                              "standard_error", "p_value", "hm_effect_allele_frequency", "n"))
setnames(t2d, c("SNP", "ea", "oa", "beta", "se", "p", "eaf", "n")); t2d <- t2d[grepl("^rs", SNP) & !is.na(beta)]
whr <- fread(file.path(CM, "WHRadjBMI_Pulit2019.tbl.gz"), select = c("SNP", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value"))
setnames(whr, c("SNP", "ea", "oa", "eaf", "beta", "se", "p")); whr[, SNP := sub(":.*", "", SNP)]
whr[, `:=`(ea = toupper(ea), oa = toupper(oa), n = 694649)]; whr <- whr[grepl("^rs", SNP)]
sbp <- fread(file.path(CM, paste0(BP, "_Evangelou2018.txt.gz")), select = c("MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P", "TotalSampleSize"))
setnames(sbp, c("marker", "ea", "oa", "eaf", "beta", "se", "p", "n"))
sbp[, c("chr", "pos37") := tstrsplit(marker, ":", keep = 1:2)][, pos37 := as.integer(pos37)]
sbp[, `:=`(ea = toupper(ea), oa = toupper(oa))]
sbp <- merge(sbp, bim[, .(chr, pos37, SNP, a1, a2)], by = c("chr", "pos37"))
sbp <- sbp[(ea == a1 & oa == a2) | (ea == a2 & oa == a1), .(SNP, ea, oa, beta, se, p, eaf, n)]
EXP <- list(T2D = t2d, WHRadjBMI = whr, BP = sbp); names(EXP)[3] <- BP
for (k in names(EXP)) { EXP[[k]] <- unique(EXP[[k]], by = "SNP"); message(k, ": ", nrow(EXP[[k]]), " variants") }

# ---------------- outcomes (GIGASTROKE EUR, harmonised)
OUTC <- c(AS = "GIGASTROKE_GCST90104539.h.tsv.gz", IS = "GIGASTROKE_GCST90104540.h.tsv.gz", CES = "GIGASTROKE_GCST90104541.h.tsv.gz",
          LAS = "GIGASTROKE_GCST90104542.h.tsv.gz", SVS = "GIGASTROKE_GCST90104543.h.tsv.gz")
YO <- lapply(OUTC, function(f) {
  d <- fread(file.path(RAW, f), select = c("rsid", "effect_allele", "other_allele", "beta", "standard_error", "p_value"))
  setnames(d, c("SNP", "ea", "oa", "beta", "se", "p")); d[, SNP := trimws(SNP)]; unique(d[grepl("^rs", SNP)], by = "SNP")
})

# ---------------- instruments: P < 5e-8, outside 3p14 +/-1 Mb, clumped r2 < 0.001 / 10 Mb (1000G EUR)
reg3 <- bim[chr == "3" & pos37 >= 63699445 & pos37 <= 65699445, SNP]
clump <- function(d, tag) {
  f <- file.path(TMP, paste0("clump_", tag, SUF)); fwrite(d[, .(SNP, P = p)], paste0(f, ".txt"), sep = "\t")
  system2(PLINK, c("--bfile", BFILE, "--clump", paste0(f, ".txt"), "--clump-p1", "1", "--clump-p2", "1", "--clump-r2", "0.001",
                   "--clump-kb", "10000", "--out", f), stdout = FALSE, stderr = FALSE)
  cl <- fread(paste0(f, ".clumped")); cl$SNP
}
IV <- lapply(names(EXP), function(k) {
  d <- EXP[[k]][p < 5e-8 & !SNP %in% reg3 & SNP %in% bim$SNP]
  s <- clump(d, k); message(k, " instruments: ", length(s)); s
}); names(IV) <- names(EXP)

pal <- function(a, b) paste0(a, b) %in% c("AT", "TA", "CG", "GC")
harm <- function(x, y) {        # align outcome y to exposure x; returns merged with by aligned
  m <- merge(x[, .(SNP, ea, oa, bx = beta, sx = se, eafx = eaf)], y[, .(SNP, ea2 = ea, oa2 = oa, by = beta, sy = se)], by = "SNP")
  m[, by := fifelse(ea == ea2 & oa == oa2, by, fifelse(ea == oa2 & oa == ea2, -by, NA_real_))]
  m[!is.na(by) & !(pal(ea, oa) & !is.na(eafx) & eafx > 0.42 & eafx < 0.58)]
}
ivw <- function(b, s, bx) {
  w <- 1 / s^2; fit <- lm(b ~ 0 + bx, weights = w)
  est <- unname(coef(fit)[1]); se <- unname(summary(fit)$coef[1, 2]) / min(1, summary(fit)$sigma); p <- 2 * pnorm(-abs(est / se))
  c(est = unname(est), se = se, p = p, nsnp = length(b))
}
uni <- rbindlist(lapply(names(EXP), function(k) rbindlist(lapply(names(YO), function(o) {
  m <- harm(EXP[[k]][SNP %in% IV[[k]]], YO[[o]]); r <- ivw(m$by, m$sy, m$bx)
  data.table(exposure = k, outcome = o, nsnp = r["nsnp"], theta = r["est"], se = r["se"], p = r["p"])
}))))
fwrite(uni, file.path(OUT, paste0("mr_cardiometabolic_to_stroke", SUF, ".tsv")), sep = "\t"); print(uni)

# ---------------- MVMR: union of instruments, re-clumped on the minimum P across exposures
U <- unique(unlist(IV))
base <- EXP$T2D[SNP %in% U, .(SNP, ea, oa, eaf)]
for (k in names(EXP)) {
  m <- merge(base[, .(SNP, ea0 = ea, oa0 = oa)], EXP[[k]][, .(SNP, ea, oa, beta, se, p)], by = "SNP")
  m[, beta := fifelse(ea == ea0 & oa == oa0, beta, fifelse(ea == oa0 & oa == ea0, -beta, NA_real_))]
  base <- merge(base, m[!is.na(beta), .(SNP, b = beta, s = se, pp = p)], by = "SNP")
  setnames(base, c("b", "s", "pp"), paste0(c("b_", "s_", "p_"), k))
}
base[, pmin := pmin(p_T2D, p_WHRadjBMI, get(paste0("p_", BP)))]
keep <- clump(base[, .(SNP, p = pmin)], "union")
base <- base[SNP %in% keep & !(pal(ea, oa) & eaf > 0.42 & eaf < 0.58)]
message("MVMR instruments: ", nrow(base))
mv <- rbindlist(lapply(names(YO), function(o) {
  y <- YO[[o]]; m <- merge(base, y[, .(SNP, ea2 = ea, oa2 = oa, by = beta, sy = se)], by = "SNP")
  m[, by := fifelse(ea == ea2 & oa == oa2, by, fifelse(ea == oa2 & oa == ea2, -by, NA_real_))]; m <- m[!is.na(by)]
  fit <- lm(as.formula(paste0("by ~ 0 + b_T2D + b_WHRadjBMI + b_", BP)), data = m, weights = 1 / sy^2)
  sm <- summary(fit)$coef; sc <- 1 / min(1, summary(fit)$sigma)
  V <- vcov(fit) * sc^2
  data.table(outcome = o, exposure = rownames(sm), theta = sm[, 1], se = sm[, 2] * sc, p = 2 * pnorm(-abs(sm[, 1] / (sm[, 2] * sc))),
             nsnp = nrow(m), vcov = list(V))
}))
fwrite(mv[, !"vcov"], file.path(OUT, paste0("mvmr_cardiometabolic_to_stroke", SUF, ".tsv")), sep = "\t"); print(mv[, !"vcov"])

# ---------------- expected effect of rs4611812-C on stroke through T2D, WHRadjBMI and SBP
al <- function(d) {     # effect of the C allele
  r <- d[SNP == KEY][1]; if (is.na(r$SNP)) return(c(NA, NA))
  s <- if (r$ea == KEY_C) 1 else if (r$oa == KEY_C) -1 else NA; c(s * r$beta, r$se)
}
bk <- sapply(EXP, al)                       # 2 x 3: beta, se of rs4611812-C on each risk factor
obs <- sapply(YO, al)
res <- rbindlist(lapply(names(YO), function(o) {
  th <- mv[outcome == o]; theta <- setNames(th$theta, sub("^b_", "", th$exposure)); V <- th$vcov[[1]]
  b <- bk[1, names(theta)]; sb <- bk[2, names(theta)]
  e_mv <- sum(b * theta); v_mv <- as.numeric(t(b) %*% V %*% b) + sum(theta^2 * sb^2)
  u <- uni[outcome == o]; e_t2d <- bk[1, "T2D"] * u[exposure == "T2D", theta]
  v_t2d <- bk[1, "T2D"]^2 * u[exposure == "T2D", se]^2 + u[exposure == "T2D", theta]^2 * bk[2, "T2D"]^2
  e_sum <- sum(sapply(names(EXP), function(k) bk[1, k] * u[exposure == k, theta]))
  data.table(outcome = o, obs_beta = obs[1, o], obs_se = obs[2, o], obs_p = 2 * pnorm(-abs(obs[1, o] / obs[2, o])),
             exp_mvmr = e_mv, exp_mvmr_se = sqrt(v_mv), exp_T2D_only = e_t2d, exp_T2D_only_se = sqrt(v_t2d),
             exp_sum_univariable = e_sum, prop_explained_mvmr = e_mv / obs[1, o],
             residual_beta = obs[1, o] - e_mv, residual_se = sqrt(obs[2, o]^2 + v_mv),
             residual_p = 2 * pnorm(-abs((obs[1, o] - e_mv) / sqrt(obs[2, o]^2 + v_mv))),
             b_T2D = bk[1, "T2D"], b_WHRadjBMI = bk[1, "WHRadjBMI"], b_BP = bk[1, BP], bp_trait = BP)
}))
fwrite(res, file.path(OUT, paste0("expected_vs_observed_rs4611812", SUF, ".tsv")), sep = "\t"); print(res)
message("mvmr step done")
