# Revision (editorial items S1, S5, S6): intracerebral hemorrhage (ICH), early AMD and East Asian stroke look-ups.
#   prep    : standardise new GWAS (ICH UKB-WGS GCST90473584; early AMD Winkler 2020 GCST010723; GIGASTROKE EAS AS/IS/LAS/SVS)
#   lookup  : key variants aligned to the FinnGen AMD risk allele
#   mr      : AMD (FinnGen any/dry/wet, IAMDGC) -> ICH; early AMD -> AS/IS/LAS/SVS/WMH/ICH (main + major-locus exclusions)
#   ldsc    : rg of AMD / early AMD with ICH, IS, SVS
# MR helpers are taken verbatim from 46_csvd_mr_ldsc.R (lines up to the "MR" step).
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table) })
setDTthreads(30)
RAW <- file.path(PROJ_ROOT, "data/raw/gwas"); DER <- file.path(PROJ_ROOT, "derived/gwas"); OUTX <- file.path(PROJ_ROOT, "results/revision/extra")
dir.create(OUTX, recursive = TRUE, showWarnings = FALSE)
STEPX <- commandArgs(TRUE); if (!length(STEPX)) STEPX <- c("prep", "lookup", "mr", "ldsc")

NEW <- list(
  ICH_UKBWGS   = list(file = "ICH_UKBWGS_GCST90473584.h.tsv.gz", ncase = 2407, ncontrol = 456033, label = "ICH (UK Biobank WGS, non-Finnish European; GCST90473584)"),
  EARLYAMD_W20 = list(file = "EARLYAMD_Winkler2020_GCST010723.h.tsv.gz", ncase = 14034, ncontrol = 91214, label = "Early AMD (Winkler 2020; GCST010723)"),
  EAS_AS  = list(file = "GIGASTROKE_EAS_AS_GCST90104544.h.tsv.gz",  ncase = 27413, ncontrol = 237242, label = "Any stroke, East Asian (GIGASTROKE)"),
  EAS_IS  = list(file = "GIGASTROKE_EAS_IS_GCST90104545.h.tsv.gz",  ncase = 19032, ncontrol = 237242, label = "Ischemic stroke, East Asian (GIGASTROKE)"),
  EAS_LAS = list(file = "GIGASTROKE_EAS_LAS_GCST90104547.h.tsv.gz", ncase = 1735,  ncontrol = 237242, label = "Large-artery stroke, East Asian (GIGASTROKE)"),
  EAS_SVS = list(file = "GIGASTROKE_EAS_SVS_GCST90104548.h.tsv.gz", ncase = 5532,  ncontrol = 237242, label = "Small-vessel stroke, East Asian (GIGASTROKE)")
)

pickcol <- function(d, cands) { x <- intersect(cands, names(d)); if (length(x)) x[1] else NA_character_ }
std_gwascat <- function(f) {
  hdr <- names(fread(f, nrows = 0))
  want <- intersect(c("hm_rsid", "rsid", "variant_id", "hm_chrom", "chromosome", "hm_pos", "base_pair_location", "hm_effect_allele",
                      "effect_allele", "hm_other_allele", "other_allele", "hm_beta", "beta", "hm_odds_ratio", "odds_ratio",
                      "standard_error", "se", "p_value", "pvalue", "hm_effect_allele_frequency", "effect_allele_frequency", "n"), hdr)
  d <- fread(f, select = want)
  m <- list(SNP = pickcol(d, c("hm_rsid", "rsid", "variant_id")), chr = pickcol(d, c("hm_chrom", "chromosome")),
            pos = pickcol(d, c("hm_pos", "base_pair_location")), ea = pickcol(d, c("hm_effect_allele", "effect_allele")),
            oa = pickcol(d, c("hm_other_allele", "other_allele")), beta = pickcol(d, c("hm_beta", "beta")),
            or = pickcol(d, c("hm_odds_ratio", "odds_ratio")), se = pickcol(d, c("standard_error", "se")),
            p = pickcol(d, c("p_value", "pvalue")), eaf = pickcol(d, c("hm_effect_allele_frequency", "effect_allele_frequency")))
  message("  column map: ", paste(names(m), unlist(m), sep = "=", collapse = " "))
  o <- data.table(SNP = d[[m$SNP]], chr = as.character(d[[m$chr]]), pos = as.numeric(d[[m$pos]]),
                  ea = toupper(d[[m$ea]]), oa = toupper(d[[m$oa]]),
                  beta = if (!is.na(m$beta)) as.numeric(d[[m$beta]]) else log(as.numeric(d[[m$or]])),
                  se = as.numeric(d[[m$se]]), p = as.numeric(d[[m$p]]),
                  eaf = if (!is.na(m$eaf)) as.numeric(d[[m$eaf]]) else NA_real_)
  if (!is.na(m$or)) o[is.na(beta), beta := log(as.numeric(d[[m$or]]))[is.na(o$beta)]]
  o
}

if ("prep" %in% STEPX) for (k in names(NEW)) {
  tr <- NEW[[k]]; f <- file.path(RAW, tr$file); o <- file.path(DER, paste0(k, ".rds"))
  if (!file.exists(f)) { message("missing ", k); next }
  if (file.exists(o) && file.mtime(o) > file.mtime(f)) { message("up-to-date ", k); next }
  t0 <- Sys.time(); message("== ", k)
  d <- std_gwascat(f)
  d <- d[!is.na(beta) & !is.na(se) & se > 0 & !is.na(p) & grepl("^rs[0-9]+$", SNP) & nchar(ea) == 1 & nchar(oa) == 1]
  if (k == "ICH_UKBWGS") d <- d[is.na(eaf) | (eaf > 0.005 & eaf < 0.995)]
  d <- d[order(p)][!duplicated(SNP)]
  d[, ncase := tr$ncase]; d[, ncontrol := tr$ncontrol]; d[, n := tr$ncase + tr$ncontrol]
  saveRDS(d, o, compress = FALSE)
  message(sprintf("   %s: %d SNPs, %d with p<5e-8, %.0fs", k, nrow(d), sum(d$p < 5e-8), as.numeric(Sys.time() - t0, units = "secs")))
  print(head(d[order(p)], 3)); rm(d); gc()
}

# ---- helpers from 46_csvd_mr_ldsc.R (definitions only)
src <- readLines(file.path(PROJ_ROOT, "scripts/46_csvd_mr_ldsc.R"))
stop_at <- grep('^if \\("mr" %in% STEP\\)', src)[1] - 1
src <- src[seq_len(stop_at)]; src <- src[!grepl("^STEP <- commandArgs", src)]
eval(parse(text = src))
OUT <- OUTX
PREV <- c(PREV, ICH_UKBWGS = 0.005, EARLYAMD_W20 = 0.10, AS = 0.03, IS = 0.025, LAS = 0.004, CES = 0.006)

# ---- look-ups
KEY <- data.table(SNP = c("rs3750848", "rs2672592", "rs4611812", "rs6795735", "rs429358", "rs7412", "rs1269852", "rs1270942",
                          "rs947073006", "rs3184504", "rs10737680"),
                  locus = c("10q26 ARMS2-HTRA1 (AMD haplotype)", "10q26 HTRA1 (stroke-side)", "3p14 ADAMTS9 (aortic eQTL lead)",
                            "3p14 ADAMTS9 (AMD lead)", "19q13 APOE e4 (rs429358)", "19q13 APOE e2 (rs7412)", "6p21 C4 (rs1269852)",
                            "6p21 C4 (rs1270942)", "9q34 ABO", "12q24 SH2B3 R262W", "1q31 CFH"))
if ("lookup" %in% STEPX) {
  amd <- load_trait("AMD_FG")[SNP %in% KEY$SNP, .(SNP, ea_amd = ea, oa_amd = oa, b_amd = beta, p_amd = p)]
  amd[, risk := fifelse(b_amd > 0, ea_amd, oa_amd)]
  tr <- intersect(c("ICH_UKBWGS", "EARLYAMD_W20", "EAS_AS", "EAS_IS", "EAS_LAS", "EAS_SVS", "AS", "IS", "LAS", "SVS"),
                  sub("\\.rds$", "", list.files(DER, "\\.rds$")))
  L <- rbindlist(lapply(tr, function(k) {
    x <- merge(load_trait(k)[SNP %in% KEY$SNP, .(SNP, ea, oa, beta, se, p, eaf)], amd[, .(SNP, risk)], by = "SNP")
    if (!nrow(x)) return(NULL)
    x[, b_risk := fifelse(ea == risk, beta, fifelse(oa == risk, -beta, NA_real_))]
    x[, .(trait = k, SNP, amd_risk_allele = risk, b_risk, se, z_risk = b_risk / se, OR_risk = exp(b_risk), p, eaf)]
  }))
  L <- merge(KEY, L, by = "SNP", all.x = TRUE)
  L <- merge(L, amd[, .(SNP, b_amd_FG = abs(b_amd), p_amd_FG = p_amd)], by = "SNP", all.x = TRUE)
  fwrite(L[order(locus, trait)], file.path(OUTX, "key_variants_ich_earlyamd_eas.tsv"), sep = "\t")
  print(L[order(locus, trait), .(locus, trait, SNP, OR_risk = round(OR_risk, 3), p = signif(p, 3))])
}

# ---- MR
if ("mr" %in% STEPX) {
  AMDX <- c("AMD_FG", "DRYAMD_FG", "WETAMD_FG", "AMD_IAMDGC")
  ivs <- list()
  for (k in AMDX) { iv <- fread(file.path(file.path(PROJ_ROOT, "results/mr"), paste0("instruments_", k, ".tsv"))); iv[, chr := as.character(chr)]; ivs[[k]] <- iv }
  if (file.exists(file.path(DER, "EARLYAMD_W20.rds"))) {
    iv <- get_iv("EARLYAMD_W20", 5e-8); iv[, chr := as.character(chr)]; ivs[["EARLYAMD_W20"]] <- iv
    fwrite(iv, file.path(OUTX, "instruments_EARLYAMD_W20.tsv"), sep = "\t"); message("early AMD IVs: ", nrow(iv))
  }
  outs_amd <- intersect("ICH_UKBWGS", sub("\\.rds$", "", list.files(DER, "\\.rds$")))
  outs_early <- intersect(c("AS", "IS", "LAS", "SVS", "WMH", "LACUNAR", "ICH_UKBWGS"), sub("\\.rds$", "", list.files(DER, "\\.rds$")))
  pairs <- rbind(if (length(outs_amd)) expand.grid(e = AMDX, o = outs_amd, stringsAsFactors = FALSE),
                 if (!is.null(ivs[["EARLYAMD_W20"]])) expand.grid(e = "EARLYAMD_W20", o = outs_early, stringsAsFactors = FALSE))
  R <- list(); S <- list()
  for (i in seq_len(nrow(pairs))) {
    ek <- pairs$e[i]; ok <- pairs$o[i]; key <- paste(ek, ok, sep = "__"); message(sprintf("[%d/%d] %s -> %s", i, nrow(pairs), ek, ok))
    r <- run_pair(ek, ok, ivs[[ek]]); if (is.null(r)) next
    R[[key]] <- r$res; S[[key]] <- r$summ
    for (sv in list(list("excl_CFH_HTRA1_C3", MAJOR[gene != "APOE"]), list("excl_CFH_HTRA1_C3_APOE", MAJOR))) {
      r2 <- run_pair(ek, ok, ivs[[ek]][!in_regions(ivs[[ek]], sv[[2]])], sv[[1]], presso = FALSE)
      if (!is.null(r2)) { R[[paste(key, sv[[1]])]] <- r2$res; S[[paste(key, sv[[1]])]] <- r2$summ }
    }
  }
  res <- rbindlist(R, fill = TRUE); summ <- rbindlist(S, fill = TRUE)
  ivw <- res[method %in% c("Inverse variance weighted", "Wald ratio")]
  ivw[analysis == "main", p_fdr := p.adjust(pval, "BH")]
  fwrite(res, file.path(OUTX, "mr_ich_earlyamd_all_methods.tsv"), sep = "\t")
  fwrite(ivw, file.path(OUTX, "mr_ich_earlyamd_ivw.tsv"), sep = "\t"); fwrite(summ, file.path(OUTX, "mr_ich_earlyamd_diagnostics.tsv"), sep = "\t")
  print(ivw[, .(exposure, outcome, analysis, nsnp, OR = round(OR, 3), lo = round(lo, 3), hi = round(hi, 3), p = signif(pval, 3), p_fdr = signif(p_fdr, 3))])
}

# ---- LDSC
if ("ldsc" %in% STEPX) {
  suppressPackageStartupMessages(library(GenomicSEM))
  REF <- file.path(PROJ_ROOT, "data/raw/ref"); LD <- file.path(REF, "eur_w_ld_chr")
  OLD <- c(file.path(PROJ_ROOT, "derived/ldsc"), file.path(PROJ_ROOT, "derived/ldsc_csvd")); W <- file.path(PROJ_ROOT, "derived/ldsc_rev")
  dir.create(W, showWarnings = FALSE); setwd(W)
  HM3 <- file.path(OLD[1], "w_hm3.snplist"); if (!file.exists(HM3)) HM3 <- file.path(LD, "w_hm3.snplist")
  lt <- intersect(c("AMD_FG", "DRYAMD_FG", "WETAMD_FG", "EARLYAMD_W20", "ICH_UKBWGS", "IS", "SVS", "WMH"), sub("\\.rds$", "", list.files(DER, "\\.rds$")))
  findss <- function(k) { f <- c(file.path(OLD, paste0(k, ".sumstats.gz")), file.path(W, paste0(k, ".sumstats.gz"))); f[file.exists(f)][1] }
  neff <- function(d) if ("ncase" %in% names(d)) 4 / (1 / d$ncase[1] + 1 / d$ncontrol[1]) else d$n[1]
  todo <- lt[is.na(sapply(lt, findss))]
  if (length(todo)) {
    fs <- c(); Ns <- c()
    for (k in todo) { d <- load_trait(k); f <- file.path(W, paste0(k, ".txt")); fwrite(d[, .(SNP, A1 = ea, A2 = oa, effect = beta, P = p)], f, sep = "\t")
                      fs <- c(fs, f); Ns <- c(Ns, neff(d)); rm(d); rm(list = k, envir = cache); gc() }
    munge(files = fs, hm3 = HM3, trait.names = todo, N = Ns, parallel = TRUE, cores = min(length(todo), 10))
  }
  ss <- sapply(lt, findss)
  binary <- c(AMD_FG = .03, DRYAMD_FG = .02, WETAMD_FG = .015, EARLYAMD_W20 = .10, ICH_UKBWGS = .005, IS = .025, SVS = .005)
  sp <- ifelse(lt %in% names(binary), 0.5, NA); pp <- ifelse(lt %in% names(binary), binary[lt], NA)
  res <- ldsc(traits = ss, sample.prev = sp, population.prev = as.numeric(pp), ld = LD, wld = LD, trait.names = lt)
  S <- res$S; V <- res$V; I <- res$I; k <- nrow(S)
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
  fwrite(out, file.path(OUTX, "ldsc_rg_ich_earlyamd.tsv"), sep = "\t"); fwrite(h2, file.path(OUTX, "ldsc_h2_ich_earlyamd.tsv"), sep = "\t")
  print(h2); print(out[order(p)][, .(trait1, trait2, rg = round(rg, 3), se = round(se, 3), p = signif(p, 3))])
}
message("76_rev_extra done: ", paste(STEPX, collapse = ","))
