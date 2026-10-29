# Round 5 (review round 5, plan A): items 1, 3, 4
# 1. Re-specified causal-expectation test: kappa (early-AMD liability per unit exposure liability), exact derivative,
#    propagated uncertainty of log OR(NHANES) and kappa; early-AMD liability as exposure (matches the NHANES phenotype).
# 3. Complement per-variant Wald ratios (functional variants and leads), within-group heterogeneity and leave-one-out,
#    IAMDGC conditionally independent complement variants with correlated IVW, power at the re-specified benchmark.
# 4. Architecture-matched positive control: atrial fibrillation liability (GWAS Catalog, Nielsen 2018) -> cardioembolic and
#    ischemic stroke, with the same causal-expectation framework.
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(MendelianRandomization) })
setDTthreads(30); set.seed(20261007)
DER <- file.path(PROJ_ROOT, "derived/gwas"); DRV <- file.path(PROJ_ROOT, "derived/revision"); SD <- file.path(PROJ_ROOT, "derived/story")
OUT <- file.path(PROJ_ROOT, "results/revision5"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
PLINK <- file.path(PROJ_ROOT, "data/raw/tools/plink/plink.exe"); BFILE <- file.path(PROJ_ROOT, "data/raw/ref/g1000_eur/g1000_eur")
TMP <- file.path(PROJ_ROOT, "tmp/r5"); dir.create(TMP, showWarnings = FALSE, recursive = TRUE)
path_of <- function(k) { f <- file.path(DRV, paste0(k, ".rds")); if (file.exists(f)) f else file.path(DER, paste0(k, ".rds")) }
NSIM <- 200000
LOR <- log(2.08); SE_LOR <- (log(3.61) - log(1.20)) / (2 * qnorm(0.975)); PREV <- c(0.03, 0.063, 0.10)
g_exact <- function(OR, p) p * (1 - p) * (OR - 1) / (1 + p * (OR - 1))     # d log-risk / d logit(prevalence)
g_bl <- function(OR, p) log(OR) * p * (1 - p)                               # Burgess-Labrecque approximation

LOCI <- data.table(locus = c("CFH", "10q26", "C3", "APOE", "ADAMTS9", "CETP", "LIPC", "ABCA1", "LPL"),
                   chr = c("1", "10", "19", "19", "3", "16", "15", "9", "8"),
                   start = c(196651878, 122454932, 6677704, 44905796, 64515654, 56961923, 58410569, 104781006, 19901717),
                   end = c(196747504, 122514908, 6720662, 44909393, 64687889, 56983845, 58569844, 104928155, 19967259))
locus_of <- function(chr, pos) { lab <- rep("other", length(pos))
  for (i in seq_len(nrow(LOCI))) lab[!is.na(pos) & as.character(chr) == LOCI$chr[i] & pos >= LOCI$start[i] - 2e6 & pos <= LOCI$end[i] + 2e6] <- LOCI$locus[i]
  lab }
LIPID <- c("CETP", "LIPC", "ABCA1", "LPL"); MAJOR3 <- c("CFH", "10q26", "C3")
EXCL <- list(main = character(0), no_10q26 = "10q26", no_major3 = MAJOR3, no_all = c(MAJOR3, "APOE", "ADAMTS9", LIPID))

# ------------------------------------------------------------------ exposures, instruments, proxies
EXPO <- c(AMD_FG = "AMD_FG", DRYAMD_FG = "DRYAMD_FG", WETAMD_FG = "WETAMD_FG", AMD_IAMDGC = "AMD_IAMDGC_fix")
inst <- fread(file.path(PROJ_ROOT, "results/revision/mr/instruments_all_schemes.tsv"))[scheme == "A_dist" & exposure %in% names(EXPO)]
EXD <- list()
for (e in names(EXPO)) { d <- readRDS(path_of(EXPO[[e]])); if (!"eaf" %in% names(d)) d[, eaf := NA_real_]; EXD[[e]] <- d
  inst[exposure == e, eafx := d$eaf[match(SNP, d$SNP)]] }
inst[, chr := as.character(chr)]
pw_map <- unique(fread(file.path(PROJ_ROOT, "results/story/per_locus_effects.tsv"))[, .(SNP, pathway, gene)], by = "SNP")
inst <- merge(inst, pw_map, by = "SNP", all.x = TRUE); inst[is.na(pathway), pathway := "other"]
# gene-based complement membership (protein-coding complement gene within +/-250 kb; Ensembl; results/revision5/complement_instrument_genes.tsv)
GENEB <- c("rs1329424", "rs139779213", "rs429608", "rs2230199", "rs6665824", "rs11732475", "rs635348", "rs11569415", "rs11080055")
inst[, pathway_gene := fifelse(SNP %in% GENEB, "complement", fifelse(pathway == "complement", "other", pathway))]
early <- readRDS(path_of("EARLYAMD_W20")); if (!"eaf" %in% names(early)) early[, eaf := NA_real_]

einst <- fread(file.path(PROJ_ROOT, "results/revision/extra/instruments_EARLYAMD_W20.tsv"))
einst[, `:=`(chr = as.character(chr), eafx = eaf, locus = locus_of(chr, pos), exposure = "EARLYAMD_W20")]

COMPV <- c("rs1329424", "rs1061170", "rs6677604", "rs2230199", "rs641153", "rs429608", "rs10033900", "rs139779213", "rs11080055",
           "rs62358361", "rs12354144", "rs72744971")
b52 <- fread(file.path(PROJ_ROOT, "results/revision/mr/input/iamdgc_52_independent_gwascatalog.tsv"))
CGRP <- list(CFH_region = c("rs570618", "rs10922109", "rs61818925", "rs35292876", "rs148553336", "rs121913059", "rs187328863", "rs191281603"),
             C2_CFB = c("rs116503776", "rs181705462", "rs144629244"), C3 = c("rs2230199", "rs147859257"),
             CFI = c("rs10033900", "rs141853578"), C9 = "rs62358361")
U <- unique(c(inst$SNP, einst$SNP, COMPV, unlist(CGRP))); writeLines(U, file.path(TMP, "snps.txt"))
system2(PLINK, c("--bfile", BFILE, "--r2", "in-phase", "with-freqs", "--ld-snp-list", file.path(TMP, "snps.txt"),
                 "--ld-window-kb", "500", "--ld-window", "99999", "--ld-window-r2", "0.8", "--threads", "30", "--out", file.path(TMP, "proxies")), stdout = FALSE, stderr = FALSE)
PX <- fread(file.path(TMP, "proxies.ld"))[SNP_A != SNP_B]
PX[, c("hapA", "hapB") := tstrsplit(PHASE, "/", fixed = TRUE)]
PX[, `:=`(a1 = substr(hapA, 1, 1), b1 = substr(hapA, 2, 2), a2 = substr(hapB, 1, 1), b2 = substr(hapB, 2, 2))]
PX <- PX[, .(target = SNP_A, proxy = SNP_B, r2 = R2, a1, b1, a2, b2)][order(target, -r2)]

get_outcome <- function(d, snps = U, use_proxy = TRUE) {
  d <- as.data.table(d)[SNP %in% c(snps, PX$proxy)]
  if (!"eaf" %in% names(d)) d[, eaf := NA_real_]
  have <- d[SNP %in% snps]; sub <- list()
  if (use_proxy) for (t in setdiff(snps, have$SNP)) {
    cand <- PX[target == t & proxy %in% d$SNP]; if (!nrow(cand)) next
    pr <- cand[1]; r <- d[SNP == pr$proxy][1]
    if (r$ea == pr$b1 && r$oa == pr$b2) { tea <- pr$a1; toa <- pr$a2 } else if (r$ea == pr$b2 && r$oa == pr$b1) { tea <- pr$a2; toa <- pr$a1 } else next
    sub[[t]] <- data.table(SNP = t, ea = tea, oa = toa, beta = r$beta, se = r$se, p = r$p, eaf = NA_real_, proxy = pr$proxy)
  }
  have[, proxy := NA_character_]
  rbind(have[, .(SNP, ea = toupper(ea), oa = toupper(oa), beta, se, p, eaf, proxy)], rbindlist(sub), fill = TRUE)
}
# CAD (Aragam 2022) from the local full harmonised file
cadx <- fread(file.path(SD, "GCST90132314.h.tsv.gz"), select = c("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "p_value", "rsid"))
setnames(cadx, c("chr", "pos", "ea", "oa", "beta", "se", "eaf", "p", "SNP")); cadx[, SNP := sub(",.*", "", as.character(SNP))]
cadx <- unique(cadx[grepl("^rs", SNP) & !is.na(beta) & !is.na(se) & se > 0], by = "SNP")
OUTC <- c("AS", "IS", "SVS", "CES", "WMH", "MEGA_AS", "MEGA_IS", "MEGA_CES")
OD <- lapply(setNames(OUTC, OUTC), function(o) get_outcome(readRDS(path_of(o))))
OD$CAD <- get_outcome(cadx[SNP %in% c(U, PX$proxy)])
OD$EARLY <- get_outcome(early)
OD$AMD_FG <- get_outcome(EXD$AMD_FG); OD$AMD_IAMDGC <- get_outcome(EXD$AMD_IAMDGC)

comp <- c(A = "T", T = "A", C = "G", G = "C")
harm <- function(ex, oc) {
  m <- merge(ex[, .(SNP, chr, pos, ea = toupper(ea), oa = toupper(oa), bx = beta, sx = se, eafx, locus,
                    pathway = if ("pathway" %in% names(ex)) pathway else NA_character_)],
             oc[, .(SNP, eo = ea, oo = oa, by = beta, sy = se, proxy)], by = "SNP")
  m[, flip := fifelse(eo == ea & oo == oa, 1, fifelse(eo == oa & oo == ea, -1,
              fifelse(comp[eo] == ea & comp[oo] == oa, 1, fifelse(comp[eo] == oa & comp[oo] == ea, -1, NA_real_))))]
  m[, pal := paste0(ea, oa) %in% c("AT", "TA", "CG", "GC")]
  m <- m[!is.na(flip) & !(pal & (is.na(eafx) | (eafx > 0.42 & eafx < 0.58)))]
  m[, by := by * flip]
  m[bx < 0, `:=`(bx = -bx, by = -by)]            # align to the exposure-increasing (AMD risk) allele
  m[]
}
ivw_mre <- function(m) {
  if (nrow(m) == 1) return(c(b = m$by / m$bx, se = m$sy / abs(m$bx), se_fe = m$sy / abs(m$bx), p = 2 * pnorm(-abs(m$by / m$sy)), n = 1, Q = 0, Qp = NA))
  w <- m$bx^2 / m$sy^2; b <- sum(w * m$by / m$bx) / sum(w); se0 <- 1 / sqrt(sum(w)); Q <- sum(w * (m$by / m$bx - b)^2)
  se <- if (nrow(m) >= 3) se0 * sqrt(max(1, Q / (nrow(m) - 1))) else se0
  c(b = b, se = se, se_fe = se0, p = 2 * pnorm(-abs(b / se)), n = nrow(m), Q = Q, Qp = pchisq(Q, nrow(m) - 1, lower.tail = FALSE))
}

# ================================================================== 1a. kappa
kap <- list()
for (e in names(EXPO)) for (x in names(EXCL)) {
  m <- harm(inst[exposure == e & !locus %in% EXCL[[x]]], OD$EARLY); r <- ivw_mre(m)
  kap[[length(kap) + 1]] <- data.table(exposure = e, set = x, nsnp = r[["n"]], kappa = r[["b"]], se = r[["se"]], se_fixed = r[["se_fe"]], p = r[["p"]], Q_p = r[["Qp"]],
                                        n_proxy = sum(!is.na(m$proxy)))
}
for (e in c("AMD_FG", "AMD_IAMDGC")) for (g in c("complement", "complement_genebased", "10q26", "lipid", "other")) {
  m <- harm(if (g == "complement_genebased") inst[exposure == e & pathway_gene == "complement"] else inst[exposure == e & pathway == g], OD$EARLY); if (!nrow(m)) next; r <- ivw_mre(m)
  kap[[length(kap) + 1]] <- data.table(exposure = e, set = paste0("pathway_", g), nsnp = r[["n"]], kappa = r[["b"]], se = r[["se"]], se_fixed = r[["se_fe"]], p = r[["p"]], Q_p = r[["Qp"]],
                                        n_proxy = sum(!is.na(m$proxy)))
}
kap <- rbindlist(kap); fwrite(kap, file.path(OUT, "kappa_earlyamd_scale.tsv"), sep = "\t"); print(kap)

# ================================================================== 1b-c. expectation test
mr <- fread(file.path(PROJ_ROOT, "results/revision/mr/mr_revised_all.tsv"))[scheme == "A_dist" & method == "IVW" & exposure %in% names(EXPO) &
        outcome %in% c("AS", "IS", "SVS", "MEGA_AS", "MEGA_IS") & exclusion %in% names(EXCL)]
lor_s <- rnorm(NSIM, LOR, SE_LOR); z0 <- rnorm(NSIM)
expect_rows <- function(e, o, x, b, se, k, sk, label) {
  ks <- k + sk * z0[sample.int(NSIM)]
  rbindlist(lapply(PREV, function(pp) rbindlist(lapply(c("exact", "BL"), function(gf) {
    gfun <- if (gf == "exact") g_exact else g_bl
    bexp_s <- ks * gfun(exp(lor_s), pp); bexp0 <- k * gfun(2.08, pp)
    p_sim <- mean(pnorm((b - bexp_s) / se))
    bup <- b + qnorm(0.95) * se; cc <- bup / k
    ormax <- if (cc <= 0) 1 else if (gf == "exact") { if (cc < 1 - pp) 1 + cc / (pp * (1 - pp - cc)) else Inf } else exp(cc / (pp * (1 - pp)))
    data.table(exposure = e, outcome = o, exclusion = x, kappa_source = label, kappa = k, kappa_se = sk, prevalence = pp, g_form = gf,
               b = b, se = se, MR_OR = exp(b), MR_lo = exp(b - 1.96 * se), MR_hi = exp(b + 1.96 * se),
               expected_OR_point = exp(bexp0), expected_OR_q025 = exp(quantile(bexp_s, 0.025)), expected_OR_q975 = exp(quantile(bexp_s, 0.975)),
               p_one_sided_point = pnorm((b - bexp0) / se), p_one_sided_propagated = p_sim,
               max_causal_OR_NHANES_type = ormax, fraction_of_observed_excess = ifelse(is.finite(ormax), log(ormax) / LOR, NA_real_))
  }))))
}
ex1 <- list()
for (i in seq_len(nrow(mr))) {
  r <- mr[i]; kk <- kap[exposure == r$exposure & set == r$exclusion]
  ex1[[length(ex1) + 1]] <- expect_rows(r$exposure, r$outcome, r$exclusion, r$b, r$se, kk$kappa, kk$se, "kappa (early AMD), same exclusion set")
  ex1[[length(ex1) + 1]] <- expect_rows(r$exposure, r$outcome, r$exclusion, r$b, r$se, 1, 0, "kappa = 1 (as in round 4)")
}
ex1 <- rbindlist(ex1); fwrite(ex1, file.path(OUT, "causal_expectation_respecified.tsv"), sep = "\t")
print(ex1[kappa_source != "kappa = 1 (as in round 4)" & prevalence == 0.063 & g_form == "exact" & outcome %in% c("AS", "IS"),
          .(exposure, outcome, exclusion, kappa = round(kappa, 3), MR_OR = round(MR_OR, 3), expOR = round(expected_OR_point, 3), p_point = signif(p_one_sided_point, 2),
            p_prop = signif(p_one_sided_propagated, 2), ORmax = round(max_causal_OR_NHANES_type, 2), frac = round(fraction_of_observed_excess, 2))])

# ================================================================== 1d. early-AMD liability as exposure
emr <- list()
for (o in c("AS", "IS", "SVS", "MEGA_AS", "MEGA_IS")) {
  m0 <- harm(einst, OD[[o]])
  for (x in c("main", "no_10q26", "no_major3")) {
    m <- m0[!locus %in% EXCL[[x]]]; if (!nrow(m)) next; r <- ivw_mre(m)
    row <- data.table(exposure = "EARLYAMD_W20", outcome = o, exclusion = x, nsnp = r[["n"]], b = r[["b"]], se = r[["se"]], p = r[["p"]],
                      OR = exp(r[["b"]]), lo = exp(r[["b"]] - 1.96 * r[["se"]]), hi = exp(r[["b"]] + 1.96 * r[["se"]]))
    for (pp in PREV) { bexp_s <- g_exact(exp(lor_s), pp); row[, paste0("p_below_expected_exact_p", pp) := mean(pnorm((r[["b"]] - bexp_s) / r[["se"]]))]
                       row[, paste0("expected_OR_exact_p", pp) := exp(g_exact(2.08, pp))] }
    emr[[length(emr) + 1]] <- row
  }
}
emr <- rbindlist(emr, fill = TRUE); fwrite(emr, file.path(OUT, "earlyamd_exposure_mr.tsv"), sep = "\t"); print(emr)

# ================================================================== 3. complement
# 3a. per-variant Wald ratios (FinnGen AMD exposure; aligned to AMD-risk allele)
fg <- EXD$AMD_FG[SNP %in% COMPV]
fgx <- fg[, .(SNP, chr = as.character(chr), pos, ea, oa, beta, se, eafx = eaf, locus = locus_of(chr, pos))]
lab <- c(rs1329424 = "CFH lead (A_dist)", rs1061170 = "CFH Y402H", rs6677604 = "CFHR1/3 deletion tag", rs2230199 = "C3 R102G",
         rs641153 = "CFB R32Q", rs429608 = "C2/CFB lead (A_dist)", rs10033900 = "CFI (IAMDGC)", rs139779213 = "CFI lead (A_dist, rare)",
         rs11080055 = "TMEM97/VTN", rs62358361 = "C9 (IAMDGC, rare)", rs12354144 = "1q31 (KCNT2), A_dist 'CFH region'",
         rs72744971 = "1q32 (CFHR/F13B region), A_dist 'CFH region'")
pv <- rbindlist(lapply(c("AS", "IS", "SVS", "CES", "WMH", "CAD"), function(o) {
  m <- harm(fgx, OD[[o]]); if (!nrow(m)) return(NULL)
  m[, .(SNP, label = lab[SNP], outcome = o, beta_amd = bx, se_amd = sx, p_amd = 2 * pnorm(-abs(bx / sx)), beta_out = by, se_out = sy,
        p_out = 2 * pnorm(-abs(by / sy)), wald = by / bx, se_wald = sy / bx, OR_per_logodds = ifelse(o == "WMH", NA, exp(by / bx)),
        lo = ifelse(o == "WMH", NA, exp(by / bx - 1.96 * sy / bx)), hi = ifelse(o == "WMH", NA, exp(by / bx + 1.96 * sy / bx)), proxy)]
}))
fwrite(pv, file.path(OUT, "complement_per_variant.tsv"), sep = "\t")
print(pv[outcome %in% c("IS", "CAD"), .(label, outcome, bx = round(beta_amd, 3), wald = round(wald, 4), se = round(se_wald, 4), p = signif(p_out, 2))])

# 3b. heterogeneity / leave-one-out within the A_dist complement group (FinnGen and IAMDGC)
het <- list()
for (def in c("annotation", "genebased")) for (e in c("AMD_FG", "AMD_IAMDGC")) for (o in c("AS", "IS", "SVS", "CAD", "WMH")) {
  m <- harm(if (def == "genebased") inst[exposure == e & pathway_gene == "complement"] else inst[exposure == e & pathway == "complement"], OD[[o]]); if (nrow(m) < 2) next
  r <- ivw_mre(m); w <- m$bx^2 / m$sy^2
  kk <- kap[exposure == e & set == ifelse(def == "genebased", "pathway_complement_genebased", "pathway_complement")]
  bexp <- if (nrow(kk)) kk$kappa * g_exact(2.08, 0.063) else NA_real_
  het[[length(het) + 1]] <- data.table(definition = def, exposure = e, outcome = o, analysis = "all complement instruments", dropped = NA_character_, nsnp = r[["n"]],
    b = r[["b"]], se = r[["se"]], p = r[["p"]], Q = r[["Q"]], Q_p = r[["Qp"]], I2 = max(0, (r[["Q"]] - (r[["n"]] - 1)) / r[["Q"]]),
    snps = paste(m$SNP, collapse = ";"), weight_share = paste(sprintf("%s:%.2f", m$SNP, w / sum(w)), collapse = ";"),
    expected_b_kappa = bexp, power_two_sided_at_expected = ifelse(o == "WMH", NA, pnorm(abs(bexp) / r[["se"]] - qnorm(0.975))),
    power_one_sided_reject_expected_if_null = ifelse(o == "WMH", NA, pnorm(bexp / r[["se"]] - qnorm(0.95))),
    p_below_expected_point = ifelse(o %in% c("AS", "IS", "SVS"), pnorm((r[["b"]] - bexp) / r[["se"]]), NA_real_),
    p_below_expected_propagated = ifelse(o %in% c("AS", "IS", "SVS") && nrow(kk), mean(pnorm((r[["b"]] - (kk$kappa + kk$se * z0) * g_exact(exp(lor_s), 0.063)) / r[["se"]])), NA_real_))
  for (s in m$SNP) { rr <- ivw_mre(m[SNP != s])
    het[[length(het) + 1]] <- data.table(definition = def, exposure = e, outcome = o, analysis = "leave-one-out", dropped = s, nsnp = rr[["n"]], b = rr[["b"]], se = rr[["se"]], p = rr[["p"]],
                                         Q = rr[["Q"]], Q_p = rr[["Qp"]]) }
  # complement without C2/CFB (MHC class III)
  mm <- m[!SNP %in% c("rs429608", "rs635348")]; if (nrow(mm) >= 1 && nrow(mm) < nrow(m)) { rr <- ivw_mre(mm)
    het[[length(het) + 1]] <- data.table(definition = def, exposure = e, outcome = o, analysis = "without C2/CFB (MHC III)", dropped = "C2/CFB", nsnp = rr[["n"]], b = rr[["b"]], se = rr[["se"]], p = rr[["p"]], Q = rr[["Q"]], Q_p = rr[["Qp"]]) }
}
het <- rbindlist(het, fill = TRUE); het[!outcome %in% "WMH", `:=`(OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se))]
fwrite(het, file.path(OUT, "complement_heterogeneity.tsv"), sep = "\t"); print(het[analysis != "leave-one-out"])

# 3c. IAMDGC conditionally independent complement variants, correlated IVW
csn <- unlist(CGRP); writeLines(csn, file.path(TMP, "comp52.txt"))
system2(PLINK, c("--bfile", BFILE, "--extract", file.path(TMP, "comp52.txt"), "--r", "square", "spaces", "--write-snplist", "--out", file.path(TMP, "comp52ld")), stdout = FALSE, stderr = FALSE)
ls_ <- readLines(file.path(TMP, "comp52ld.snplist")); RL <- as.matrix(read.table(file.path(TMP, "comp52ld.ld"))); dimnames(RL) <- list(ls_, ls_)
bimc <- fread(paste0(BFILE, ".bim"), header = FALSE)[V2 %in% ls_, .(SNP = V2, A1 = V5, A2 = V6)]
iam <- EXD$AMD_IAMDGC[SNP %in% csn, .(SNP, chr = as.character(chr), pos, ea, oa, beta, se, eafx = eaf, locus = locus_of(chr, pos))]
civ <- list()
for (o in c("IS", "AS", "SVS", "CAD")) {
  m <- harm(iam, get_outcome(if (o == "CAD") cadx[SNP %in% csn] else readRDS(path_of(o)), snps = csn, use_proxy = FALSE))
  grp_of <- setNames(rep(names(CGRP), lengths(CGRP)), unlist(CGRP)); m[, group := grp_of[SNP]]
  sets <- list(all_complement = m$SNP, without_C2_CFB = m[group != "C2_CFB", SNP], CFH_region = m[group == "CFH_region", SNP])
  for (sn in names(sets)) {
    mm <- m[SNP %in% intersect(sets[[sn]], rownames(RL))]; if (nrow(mm) < 2) next
    a1 <- bimc[match(mm$SNP, SNP)]; sgn <- ifelse(mm$ea == a1$A1, 1, ifelse(mm$ea == a1$A2, -1, NA)); sgn[is.na(sgn)] <- 1
    rho <- RL[mm$SNP, mm$SNP, drop = FALSE] * outer(sgn, sgn)
    inp <- mr_input(bx = mm$bx, bxse = mm$sx, by = mm$by, byse = mm$sy, correlation = rho, snps = mm$SNP)
    rr <- tryCatch(mr_ivw(inp, correl = TRUE, model = if (nrow(mm) >= 3) "random" else "fixed"), error = function(e) NULL); if (is.null(rr)) next
    civ[[length(civ) + 1]] <- data.table(exposure = "AMD_IAMDGC (conditionally independent complement variants)", outcome = o, set = sn, nsnp = nrow(mm),
      snps = paste(mm$SNP, collapse = ";"), b = rr@Estimate, se = rr@StdError, p = rr@Pvalue,
      OR = exp(rr@Estimate), lo = exp(rr@Estimate - 1.96 * rr@StdError), hi = exp(rr@Estimate + 1.96 * rr@StdError),
      heterogeneity_p = tryCatch(rr@Heter.Stat[2], error = function(e) NA_real_))
  }
}
civ <- rbindlist(civ); fwrite(civ, file.path(OUT, "complement_iamdgc_cor_ivw.tsv"), sep = "\t"); print(civ)

# ================================================================== 4. AF positive control
aff <- file.path(OUT, "af_gwascatalog_associations.tsv")
if (file.exists(aff)) {
  af <- fread(aff)[!is.na(p) & p < 5e-8 & risk_allele %in% c("A", "C", "G", "T")]
  af <- af[!is.na(or_per_copy) & or_per_copy > 0]
  af[, b := log(or_per_copy)]
  af[, se := abs(b) / qnorm(p / 2, lower.tail = FALSE)]          # SE from the reported P value (catalog CIs are rounded)
  af[, se_ci := (log(ci_upper) - log(ci_lower)) / (2 * 1.96)]
  af[b < 0, `:=`(b = -b)]   # risk allele reported with OR < 1 should not occur; guard
  pos <- unique(readRDS(path_of("IS"))[, .(SNP, chr = as.character(chr), pos, ea_o = ea, oa_o = oa)], by = "SNP")
  af <- merge(af, pos, by = "SNP")
  af[, oa := fifelse(risk_allele == ea_o, oa_o, fifelse(risk_allele == oa_o, ea_o, NA_character_))]
  af <- af[!is.na(oa) & !(paste0(risk_allele, oa) %in% c("AT", "TA", "CG", "GC"))]
  af <- af[order(p)]; keep <- logical(nrow(af)); tk <- data.table(chr = character(0), pos = numeric(0))
  for (i in seq_len(nrow(af))) if (!nrow(tk[chr == af$chr[i] & abs(pos - af$pos[i]) < 1e6])) { keep[i] <- TRUE; tk <- rbind(tk, data.table(chr = af$chr[i], pos = af$pos[i])) }
  af <- af[keep]
  afi <- af[, .(SNP, chr, pos, ea = risk_allele, oa, beta = b, se, eafx = NA_real_,
                locus = fifelse(chr == "4" & pos > 109.6e6 & pos < 111.7e6, "PITX2", fifelse(chr == "16" & pos > 71.8e6 & pos < 74.1e6, "ZFHX3", "other")))]
  fwrite(afi, file.path(OUT, "af_instruments.tsv"), sep = "\t")
  ODAF <- lapply(setNames(c("CES", "IS", "AS", "MEGA_CES", "MEGA_IS"), c("CES", "IS", "AS", "MEGA_CES", "MEGA_IS")),
                 function(o) get_outcome(readRDS(path_of(o)), snps = afi$SNP, use_proxy = FALSE))
  pcr <- list()
  for (o in names(ODAF)) {
    m0 <- harm(afi, ODAF[[o]])
    for (x in c("all", "without_PITX2", "without_PITX2_ZFHX3", "PITX2_only")) {
      m <- switch(x, all = m0, without_PITX2 = m0[locus != "PITX2"], without_PITX2_ZFHX3 = m0[!locus %in% c("PITX2", "ZFHX3")], PITX2_only = m0[locus == "PITX2"])
      if (!nrow(m)) next; r <- ivw_mre(m)
      w <- m$bx^2 / m$sy^2
      row <- data.table(exposure = "Atrial fibrillation liability (GWAS Catalog associations)", outcome = o, set = x, nsnp = r[["n"]], b = r[["b"]], se = r[["se"]], p = r[["p"]],
                        OR = exp(r[["b"]]), lo = exp(r[["b"]] - 1.96 * r[["se"]]), hi = exp(r[["b"]] + 1.96 * r[["se"]]),
                        pitx2_weight_share = sum(w[m$locus == "PITX2"]) / sum(w))
      for (orx in c(5, 2.5)) for (pp in c(0.02, 0.03, 0.05)) {
        bexp_e <- g_exact(orx, pp); bexp_b <- g_bl(orx, pp)
        row[, paste0("expOR_exact_OR", orx, "_p", pp) := exp(bexp_e)]
        row[, paste0("p_below_expected_exact_OR", orx, "_p", pp) := pnorm((r[["b"]] - bexp_e) / r[["se"]])]
        row[, paste0("expOR_BL_OR", orx, "_p", pp) := exp(bexp_b)]
      }
      pcr[[length(pcr) + 1]] <- row
    }
  }
  pcr <- rbindlist(pcr, fill = TRUE); fwrite(pcr, file.path(OUT, "positive_control_af.tsv"), sep = "\t")
  print(pcr[, .(outcome, set, nsnp, OR = round(OR, 3), lo = round(lo, 3), hi = round(hi, 3), p = signif(p, 2), pitx2_w = round(pitx2_weight_share, 2),
                expOR = round(expOR_exact_OR5_p0.03, 3), p_below = signif(p_below_expected_exact_OR5_p0.03, 2))])
} else message("AF associations file missing")
message("r5_01 done")
