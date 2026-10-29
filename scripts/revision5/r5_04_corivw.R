# Round 5, item 3c (revised): conditionally independent complement signals (IAMDGC 52-variant set; Fritsche 2016) with
# FinnGen AMD effect sizes (the corrected IAMDGC file lacks most of these variants) and correlated-instrument IVW (1000G EUR LD).
# Outcomes: GIGASTROKE IS, AS, SVS; CAD (Aragam 2022). Missing outcome variants replaced by 1000G proxies (r2 > 0.8).
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
DER <- file.path(PROJ_ROOT, "derived/gwas"); DRV <- file.path(PROJ_ROOT, "derived/revision"); SD <- file.path(PROJ_ROOT, "derived/story")
OUT <- file.path(PROJ_ROOT, "results/revision5"); TMP <- file.path(PROJ_ROOT, "tmp/r5")
PLINK <- file.path(PROJ_ROOT, "data/raw/tools/plink/plink.exe"); BFILE <- file.path(PROJ_ROOT, "data/raw/ref/g1000_eur/g1000_eur")
path_of <- function(k) { f <- file.path(DRV, paste0(k, ".rds")); if (file.exists(f)) f else file.path(DER, paste0(k, ".rds")) }
CGRP <- list(CFH_region = c("rs570618", "rs10922109", "rs61818925", "rs35292876", "rs148553336", "rs121913059", "rs187328863", "rs191281603"),
             C2_CFB = c("rs116503776", "rs181705462", "rs144629244"), C3 = c("rs2230199", "rs147859257"),
             CFI = c("rs10033900", "rs141853578"), C9 = "rs62358361")
grp_of <- setNames(rep(names(CGRP), lengths(CGRP)), unlist(CGRP)); csn <- unlist(CGRP)
PX <- fread(file.path(TMP, "proxies.ld"))[SNP_A != SNP_B & SNP_A %in% csn]
PX[, c("hapA", "hapB") := tstrsplit(PHASE, "/", fixed = TRUE)]
PX[, `:=`(a1 = substr(hapA, 1, 1), b1 = substr(hapA, 2, 2), a2 = substr(hapB, 1, 1), b2 = substr(hapB, 2, 2))]
PX <- PX[, .(target = SNP_A, proxy = SNP_B, r2 = R2, a1, b1, a2, b2)][order(target, -r2)]
getx <- function(d, snps) {
  d <- as.data.table(d)[SNP %in% c(snps, PX$proxy)]; have <- d[SNP %in% snps]; sub <- list()
  for (t in setdiff(snps, have$SNP)) { cand <- PX[target == t & proxy %in% d$SNP]; if (!nrow(cand)) next
    pr <- cand[1]; r <- d[SNP == pr$proxy][1]
    if (r$ea == pr$b1 && r$oa == pr$b2) { tea <- pr$a1; toa <- pr$a2 } else if (r$ea == pr$b2 && r$oa == pr$b1) { tea <- pr$a2; toa <- pr$a1 } else next
    sub[[t]] <- data.table(SNP = t, ea = tea, oa = toa, beta = r$beta, se = r$se, proxy = pr$proxy) }
  rbind(have[, .(SNP, ea = toupper(ea), oa = toupper(oa), beta, se, proxy = NA_character_)], rbindlist(sub), fill = TRUE)
}
fg <- readRDS(path_of("AMD_FG")); EXP <- getx(fg, csn); EXP <- merge(EXP, fg[, .(SNP, eaf)], by = "SNP", all.x = TRUE)
message("FinnGen exposure rows: ", nrow(EXP), " (proxied ", sum(!is.na(EXP$proxy)), ")")
cadx <- fread(file.path(SD, "GCST90132314.h.tsv.gz"), select = c("effect_allele", "other_allele", "beta", "standard_error", "rsid"))
setnames(cadx, c("ea", "oa", "beta", "se", "SNP")); cadx[, SNP := sub(",.*", "", as.character(SNP))]
cadx <- unique(cadx[SNP %in% c(csn, PX$proxy) & !is.na(beta)], by = "SNP")
writeLines(csn, file.path(TMP, "comp52.txt"))
system2(PLINK, c("--bfile", BFILE, "--extract", file.path(TMP, "comp52.txt"), "--r", "square", "spaces", "--write-snplist", "--out", file.path(TMP, "comp52ld")), stdout = FALSE, stderr = FALSE)
ls_ <- readLines(file.path(TMP, "comp52ld.snplist")); RL <- as.matrix(read.table(file.path(TMP, "comp52ld.ld"))); dimnames(RL) <- list(ls_, ls_)
bimc <- fread(paste0(BFILE, ".bim"), header = FALSE)[V2 %in% ls_, .(SNP = V2, A1 = V5, A2 = V6)]
comp <- c(A = "T", T = "A", C = "G", G = "C")
res <- list(); per <- list()
for (o in c("IS", "AS", "SVS", "CAD")) {
  oc <- getx(if (o == "CAD") cadx else readRDS(path_of(o)), csn)
  m <- merge(EXP[, .(SNP, ea, oa, bx = beta, sx = se, eafx = eaf, px = proxy)], oc[, .(SNP, eo = ea, oo = oa, by = beta, sy = se, py = proxy)], by = "SNP")
  m[, flip := fifelse(eo == ea & oo == oa, 1, fifelse(eo == oa & oo == ea, -1, fifelse(comp[eo] == ea & comp[oo] == oa, 1, fifelse(comp[eo] == oa & comp[oo] == ea, -1, NA_real_))))]
  m <- m[!is.na(flip) & !(paste0(ea, oa) %in% c("AT", "TA", "CG", "GC") & (is.na(eafx) | (eafx > 0.42 & eafx < 0.58)))]
  m[, by := by * flip]; m[, group := grp_of[SNP]]
  per[[o]] <- m[, .(outcome = o, SNP, group, ea, bx, sx, by, sy, wald = by / bx, se_wald = sy / abs(bx), px, py)]
  sets <- list(all_complement = m$SNP, without_C2_CFB = m[group != "C2_CFB", SNP], CFH_region = m[group == "CFH_region", SNP])
  for (sn in names(sets)) {
    mm <- m[SNP %in% intersect(sets[[sn]], rownames(RL))]; if (nrow(mm) < 2) next
    a1 <- bimc[match(mm$SNP, SNP)]; sgn <- ifelse(mm$ea == a1$A1, 1, ifelse(mm$ea == a1$A2, -1, NA)); sgn[is.na(sgn)] <- 1
    rho <- RL[mm$SNP, mm$SNP, drop = FALSE] * outer(sgn, sgn)
    inp <- mr_input(bx = mm$bx, bxse = mm$sx, by = mm$by, byse = mm$sy, correlation = rho, snps = mm$SNP)
    rr <- tryCatch(mr_ivw(inp, correl = TRUE, model = if (nrow(mm) >= 3) "random" else "fixed"), error = function(e) { message(e$message); NULL }); if (is.null(rr)) next
    res[[length(res) + 1]] <- data.table(exposure = "FinnGen AMD at IAMDGC conditionally independent complement variants", outcome = o, set = sn, nsnp = nrow(mm),
      snps = paste(mm$SNP, collapse = ";"), groups = paste(mm$group, collapse = ";"), b = rr@Estimate, se = rr@StdError, p = rr@Pvalue,
      OR = exp(rr@Estimate), lo = exp(rr@Estimate - 1.96 * rr@StdError), hi = exp(rr@Estimate + 1.96 * rr@StdError),
      heterogeneity_p = tryCatch(rr@Heter.Stat[2], error = function(e) NA_real_))
  }
}
res <- rbindlist(res); fwrite(res, file.path(OUT, "complement_iamdgc_cor_ivw.tsv"), sep = "\t"); print(res)
fwrite(rbindlist(per), file.path(OUT, "complement_iamdgc_variants_harmonised.tsv"), sep = "\t")
message("r5_04 done")
