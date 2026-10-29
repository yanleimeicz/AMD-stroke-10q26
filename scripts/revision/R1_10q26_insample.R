# Revision R1 (10q26): is the AMD/stroke antagonism at ARMS2-HTRA1 one signal or two?
# Uses (i) FinnGen R12 in-sample SuSiE fine-mapping of AMD (log Bayes factors per single effect),
#      (ii) FinnGen SISu4.2 in-sample LD (api.finngen.fi) and 1000G EUR / EAS LD,
#      (iii) GIGASTROKE EUR and EAS, FinnGen I9_STR (same population as FinnGen AMD), lacunar stroke and WMH.
# Analyses: key-variant table; two-SNP conditional in both directions with matched LD (+/-0.1 sensitivity);
# three-SNP joint model; proportional-colocalisation test (single shared causal configuration);
# coloc.bf_bf of AMD SuSiE single effects vs stroke; coloc.abf FinnGen AMD vs FinnGen stroke with p12 sensitivity;
# GIGASTROKE per-variant effective N; prior dependence of the earlier PP.H4.
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(coloc) })
T <- file.path(PROJ_ROOT, "tmp/rev10q26"); OUT <- file.path(PROJ_ROOT, "results/revision/10q26")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
lo <- 121954932; hi <- 122954932                     # GRCh38, +/-500 kb around ARMS2-HTRA1

KEY <- data.table(rsid = c("rs10490924", "rs3750848", "rs36212733", "rs11200638", "rs2672592", "rs72631113"),
                  pos = c(122454932, 122455799, 122455695, 122461028, 122471234, 122453933),
                  risk = c("T", "G", "C", "A", "T", "T"))       # AMD-risk alleles (FinnGen H7_AMD positive beta)

# ---------------------------------------------------------------- loaders -> data.table(pos, ea, oa, beta, se, eaf, p, rsid)
std <- function(d, label, n = NA_real_, ncase = NA_real_) {
  d <- d[pos >= lo & pos <= hi & is.finite(beta) & is.finite(se) & se > 0]
  d[, `:=`(z = beta / se, trait = label, N = n, s = if (is.na(ncase)) NA_real_ else ncase / n)]
  d[, id := paste(pos, pmin(ea, oa), pmax(ea, oa), sep = ":")]
  unique(d, by = "id")
}
fin_amd <- function(e) {
  x <- fread(file.path(T, paste0(e, ".10q26.tsv")))
  d <- x[, .(pos = position, oa = allele1, ea = allele2, beta, se, eaf = maf, p, rsid = NA_character_)]
  list(ss = d, susie = x)
}
gw <- function(f) {
  x <- fread(file.path(T, f))
  x[, .(pos = base_pair_location, ea = effect_allele, oa = other_allele, beta, se = standard_error,
        eaf = effect_allele_frequency, p = p_value, rsid)]
}
NC <- list(H7_AMD = c(12495, 461686), DRY_AMD = c(8570, 329258), WET_AMD = c(6699, 331070))
amd <- lapply(names(NC), function(e) { a <- fin_amd(e); a$ss <- std(a$ss, paste0("FinnGen_", e), sum(NC[[e]]), NC[[e]][1]); a })
names(amd) <- names(NC)

stroke <- list(
  GIGA_EUR_AS  = std(gw("eur_539.10q26.tsv"), "GIGA_EUR_AS",  73652 + 1234808, 73652),
  GIGA_EUR_IS  = std(gw("eur_540.10q26.tsv"), "GIGA_EUR_IS",  62100 + 1234808, 62100),
  GIGA_EUR_CES = std(gw("eur_541.10q26.tsv"), "GIGA_EUR_CES", 10804 + 1234808, 10804),
  GIGA_EUR_LAS = std(gw("eur_542.10q26.tsv"), "GIGA_EUR_LAS", 6399 + 1234808, 6399),
  GIGA_EUR_SVS = std(gw("eur_543.10q26.tsv"), "GIGA_EUR_SVS", 6811 + 1234808, 6811),
  GIGA_EAS_AS  = std(gw("eas_AS_GCST90104544.10q26.tsv"),  "GIGA_EAS_AS",  27413 + 237242, 27413),
  GIGA_EAS_IS  = std(gw("eas_IS_GCST90104545.10q26.tsv"),  "GIGA_EAS_IS",  19032 + 237242, 19032),
  GIGA_EAS_SVS = std(gw("eas_SVS_GCST90104548.10q26.tsv"), "GIGA_EAS_SVS", 5532 + 237242, 5532))
fgs <- file.path(T, "fg_I9_STR_10q26.tsv")
if (file.exists(fgs)) {
  x <- fread(fgs); setnames(x, 1, "chrom")
  stroke$FinnGen_I9_STR <- std(x[, .(pos, oa = ref, ea = alt, beta, se = sebeta, eaf = af_alt, p = pval, rsid = rsids)],
                               "FinnGen_I9_STR", 34110 + 450023, 34110)
}
gr <- readRDS(file.path(PROJ_ROOT, "derived/tissue_eqtl/gwas_regions.rds"))
for (tr in intersect(c("LACUNAR", "WMH"), unique(gr$trait))) {
  y <- gr[trait == tr & region == "10q26_HTRA1", .(pos = pos38, ea, oa, beta, se, eaf, p, rsid = SNP)]
  nn <- suppressWarnings(gr[trait == tr, max(n, na.rm = TRUE)]); nc <- suppressWarnings(gr[trait == tr, max(ncase, na.rm = TRUE)])
  if (tr == "LACUNAR") { nn <- 6030 + 219389; nc <- 6030 }
  stroke[[tr]] <- std(y, tr, nn, if (is.finite(nc)) nc else NA_real_)
}
if (!"LACUNAR" %in% names(stroke) && file.exists(file.path(PROJ_ROOT, "derived/gwas/LACUNAR.rds"))) {
  y <- as.data.table(readRDS(file.path(PROJ_ROOT, "derived/gwas/LACUNAR.rds")))
  if ("pos38" %in% names(y)) y[, pos := pos38]
  y <- y[chr %in% c("10", 10) & pos >= lo & pos <= hi, .(pos, ea, oa, beta, se, eaf, p, rsid = SNP)]
  stroke$LACUNAR <- std(y, "LACUNAR", 6030 + 219389, 6030)
}
cat("datasets:", names(stroke), "\n")

# ---------------------------------------------------------------- key-variant table aligned to the AMD-risk allele
align_key <- function(d) {
  m <- merge(KEY, d, by = "pos")
  m[, sgn := fifelse(ea == risk, 1, fifelse(oa == risk, -1, NA_real_))]
  m <- m[!is.na(sgn)]
  m[, .(rsid = rsid.x, pos, risk, beta_risk = beta * sgn, se, z_risk = z * sgn, p, eaf_risk = fifelse(sgn == 1, eaf, 1 - eaf),
        Neff_rel = 1 / (se^2 * 2 * eaf * (1 - eaf)))]
}
keytab <- rbindlist(c(lapply(amd, function(a) align_key(a$ss)[, trait := a$ss$trait[1]]),
                      lapply(names(stroke), function(k) align_key(stroke[[k]])[, trait := k])), fill = TRUE)
keytab[, Neff_ratio_to_max := Neff_rel / max(Neff_rel), by = trait]
fwrite(keytab, file.path(OUT, "key_variants_all_datasets.tsv"), sep = "\t")

# ---------------------------------------------------------------- LD sources (signed r for AMD-risk alleles; all positive)
LDsrc <- list(
  FinnGen_SISu = c(A_B = sqrt(0.440), A_C = sqrt(0.456), B_C = sqrt(0.955)),
  KG_EUR       = c(A_B = 0.618701, A_C = 0.644931, B_C = 0.950453),
  KG_EAS       = c(A_B = sqrt(0.675492), A_C = sqrt(0.460484), B_C = sqrt(0.293505)))
ld_for <- function(k) if (grepl("^FinnGen", k)) "FinnGen_SISu" else if (grepl("EAS", k)) "KG_EAS" else "KG_EUR"

zget <- function(k, d) { t <- align_key(d); setNames(t$z_risk, t$rsid) }
cond2 <- function(za, zb, r) c(zC_b_given_a = (zb - r * za) / sqrt(1 - r^2), zC_a_given_b = (za - r * zb) / sqrt(1 - r^2))
pairs <- list(A_B = c("rs3750848", "rs2672592"), A_C = c("rs3750848", "rs72631113"), B_C = c("rs2672592", "rs72631113"))
allsets <- c(setNames(lapply(amd, `[[`, "ss"), sapply(amd, function(a) a$ss$trait[1])), stroke)
cres <- rbindlist(lapply(names(allsets), function(k) {
  z <- zget(k, allsets[[k]]); src <- ld_for(k)
  rbindlist(lapply(names(pairs), function(pn) {
    pr <- pairs[[pn]]; if (!all(pr %in% names(z))) return(NULL)
    r0 <- LDsrc[[src]][[pn]]
    rbindlist(lapply(c(-0.1, 0, 0.1), function(dr) {
      r <- min(0.99, r0 + dr); cz <- cond2(z[pr[1]], z[pr[2]], r)
      data.table(trait = k, ld_source = src, snp_a = pr[1], snp_b = pr[2], r = r, r_shift = dr,
                 z_a = z[pr[1]], z_b = z[pr[2]], p_a = 2 * pnorm(-abs(z[pr[1]])), p_b = 2 * pnorm(-abs(z[pr[2]])),
                 zC_b_given_a = cz[1], pC_b_given_a = 2 * pnorm(-abs(cz[1])),
                 zC_a_given_b = cz[2], pC_a_given_b = 2 * pnorm(-abs(cz[2])))
    }))
  }))
}))
fwrite(cres, file.path(OUT, "conditional_two_snp_matched_LD.tsv"), sep = "\t")

# three-SNP joint model (rs3750848, rs2672592, rs72631113)
jres <- rbindlist(lapply(names(allsets), function(k) {
  z <- zget(k, allsets[[k]]); s3 <- c("rs3750848", "rs2672592", "rs72631113"); if (!all(s3 %in% names(z))) return(NULL)
  L <- LDsrc[[ld_for(k)]]; R <- matrix(c(1, L["A_B"], L["A_C"], L["A_B"], 1, L["B_C"], L["A_C"], L["B_C"], 1), 3)
  Ri <- solve(R); bj <- Ri %*% z[s3]; zj <- as.vector(bj / sqrt(diag(Ri)))
  data.table(trait = k, ld_source = ld_for(k), snp = s3, z_marginal = z[s3], z_joint = zj, p_joint = 2 * pnorm(-abs(zj)))
}))
fwrite(jres, file.path(OUT, "joint_three_snp_model.tsv"), sep = "\t")

# ---------------------------------------------------------------- proportional colocalisation test
# H0: one shared causal configuration -> E[z_stroke] = k * E[z_AMD] across markers. chi2 = min_k e' (R_s + k^2 R_a)^-1 e
prop_test <- function(zA, zS, RA, RS, ridge = 0) {
  RA <- RA + diag(ridge, nrow(RA)); RS <- RS + diag(ridge, nrow(RS))
  f <- function(k) { e <- zS - k * zA; as.numeric(t(e) %*% solve(RS + k^2 * RA) %*% e) }
  o <- optimize(f, c(-2, 2)); list(k = o$minimum, chi2 = o$objective, df = length(zA) - 1,
                                   p = pchisq(o$objective, length(zA) - 1, lower.tail = FALSE))
}
mk <- function(src, s) { L <- LDsrc[[src]]; m <- c(rs3750848 = "A", rs2672592 = "B", rs72631113 = "C")
  R <- diag(length(s)); for (i in seq_along(s)) for (j in seq_along(s)) if (i != j) R[i, j] <- L[[paste(sort(m[c(s[i], s[j])]), collapse = "_")]]; R }
s3 <- c("rs3750848", "rs2672592", "rs72631113")
pres <- rbindlist(lapply(names(amd), function(e) {
  zA <- zget(e, amd[[e]]$ss)[s3]
  rbindlist(lapply(names(stroke), function(k) {
    zS <- zget(k, stroke[[k]]); if (!all(s3 %in% names(zS)) || any(is.na(zA))) return(NULL); zS <- zS[s3]
    rbindlist(lapply(c(0, 0.02, 0.05), function(rg) {
      o <- prop_test(zA, zS, mk("FinnGen_SISu", s3), mk(ld_for(k), s3), rg)
      data.table(amd = e, stroke = k, stroke_ld = ld_for(k), ridge = rg, k = o$k, chi2 = o$chi2, df = o$df, p = o$p,
                 zAMD = paste(round(zA, 2), collapse = "/"), zStroke = paste(round(zS, 2), collapse = "/"))
    }))
  }))
}))
fwrite(pres, file.path(OUT, "proportionality_test.tsv"), sep = "\t")

# ---------------------------------------------------------------- coloc.bf_bf: AMD SuSiE single effects vs stroke (ABF)
abf_vec <- function(d) {
  if (all(is.na(d$s))) return(NULL)
  lbf <- coloc:::approx.bf.estimates(z = d$z, V = d$se^2, type = "cc")   # sd.prior 0.2 for case-control
  setNames(lbf$lABF, d$id)
}
bres <- rbindlist(lapply(names(amd), function(e) {
  sx <- amd[[e]]$susie[position >= lo & position <= hi]
  sx[, id := paste(position, pmin(allele1, allele2), pmax(allele1, allele2), sep = ":")]
  sx <- unique(sx, by = "id")
  css <- sort(unique(sx[cs > 0]$cs)); if (!length(css)) return(NULL)
  bf1 <- t(as.matrix(sx[, paste0("lbf_variable", css), with = FALSE])); colnames(bf1) <- sx$id; rownames(bf1) <- paste0("AMD_CS", css)
  rbindlist(lapply(names(stroke), function(k) {
    b2 <- abf_vec(stroke[[k]]); if (is.null(b2)) return(NULL)
    common <- intersect(colnames(bf1), names(b2)); if (length(common) < 50) return(NULL)
    bf2 <- matrix(b2[common], nrow = 1, dimnames = list(k, common))
    rbindlist(lapply(c(1e-5, 5e-6, 1e-6), function(p12) {
      r <- coloc.bf_bf(bf1[, common, drop = FALSE], bf2, p1 = 1e-4, p2 = 1e-4, p12 = p12)$summary
      data.table(amd = e, stroke = k, p12 = p12, n_snps = length(common), r)
    }))
  }))
}), fill = TRUE)
fwrite(bres, file.path(OUT, "coloc_bf_bf_AMDsusie_vs_stroke.tsv"), sep = "\t")

# ---------------------------------------------------------------- coloc.abf with p12 sensitivity (single-variant)
cabf <- rbindlist(lapply(names(amd), function(e) {
  a <- amd[[e]]$ss
  rbindlist(lapply(names(stroke), function(k) {
    b <- stroke[[k]]; m <- merge(a, b, by = "id", suffixes = c(".a", ".b"))
    m <- m[is.finite(eaf.a) & eaf.a > 0 & eaf.a < 1 & (!is.na(s.b[1]) | (is.finite(eaf.b) & eaf.b > 0 & eaf.b < 1))]
    if (nrow(m) < 50) return(NULL)
    d1 <- list(beta = m$beta.a, varbeta = m$se.a^2, type = "cc", s = m$s.a[1], N = m$N.a[1], snp = m$id, MAF = pmin(m$eaf.a, 1 - m$eaf.a))
    d2 <- if (is.na(m$s.b[1])) list(beta = m$beta.b, varbeta = m$se.b^2, type = "quant", N = m$N.b[1], snp = m$id, MAF = pmin(m$eaf.b, 1 - m$eaf.b)) else
      list(beta = m$beta.b, varbeta = m$se.b^2, type = "cc", s = m$s.b[1], N = m$N.b[1], snp = m$id)
    rbindlist(lapply(c(1e-5, 5e-6, 2e-6, 1e-6), function(p12) {
      invisible(capture.output(r <- suppressMessages(coloc.abf(d1, d2, p12 = p12))))
      top <- r$results[which.max(r$results$SNP.PP.H4), ]
      data.table(amd = e, stroke = k, p12 = p12, n_snps = nrow(m), t(r$summary[2:6]), top_snp_H4 = top$snp, top_SNP.PP.H4 = top$SNP.PP.H4)
    }))
  }))
}))
fwrite(cabf, file.path(OUT, "coloc_abf_p12_sensitivity.tsv"), sep = "\t")

# ---------------------------------------------------------------- prior dependence of the reported PP.H4 = 0.95 (AMD-SVS)
cg <- fread(file.path(PROJ_ROOT, "results/coloc/coloc_gwas_gwas.tsv"))
if (all(c("PP.H3", "PP.H4") %in% names(cg))) {
  cg[, ratio := PP.H4 / PP.H3][, f_shared := (ratio / 1000) / (1 + ratio / 1000)]
  fwrite(cg[PP.H4 > 0.5][order(-PP.H4)], file.path(OUT, "earlier_coloc_prior_dependence.tsv"), sep = "\t")
}

# ---------------------------------------------------------------- effective N across the region (GIGASTROKE)
neff <- rbindlist(lapply(names(stroke), function(k) {
  d <- stroke[[k]][eaf > 0.05 & eaf < 0.95]; d[, Neff_rel := 1 / (se^2 * 2 * eaf * (1 - eaf))]
  data.table(trait = k, n_var = nrow(d), Neff_median = median(d$Neff_rel), Neff_p05 = quantile(d$Neff_rel, 0.05),
             Neff_p95 = quantile(d$Neff_rel, 0.95))
}))
fwrite(neff, file.path(OUT, "effective_N_region.tsv"), sep = "\t")
cat("done\n")
