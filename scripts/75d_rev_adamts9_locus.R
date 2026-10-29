# Revision R4/R10/S6 (3p14 ADAMTS9 locus): cardiometabolic alternative, multi-trait colocalization, subtype heterogeneity,
# coronary disease, dry vs wet AMD, replication, PheWAS and tissue-scan disclosure.
# Inputs: derived/adamts9 (region GWAS, aorta eQTL), derived/revision_adamts9 (tabix extracts, FinnGen SuSiE),
#         data/raw/gwas/cardiometab (T2D, WHRadjBMI, SBP, DBP, MEGASTROKE), tmp/rev_adamts9 (PheWAS, GIGASTROKE key SNPs)
# Output -> results/revision/adamts9/
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib_hypr"), file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressMessages({ library(data.table); library(coloc); library(hyprcoloc) })
setDTthreads(16)
D0 <- file.path(PROJ_ROOT, "derived/adamts9"); D1 <- file.path(PROJ_ROOT, "derived/revision_adamts9")
CM <- file.path(PROJ_ROOT, "data/raw/gwas/cardiometab"); TMP <- file.path(PROJ_ROOT, "tmp/rev_adamts9")
OUT <- file.path(PROJ_ROOT, "results/revision/adamts9"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
KEY <- "rs4611812"; KEY_C <- "C"; POS38 <- 64713769; W <- 5e5
rs <- fread(file.path(D0, "region_snps.tsv"))                     # SNP pos38 pos37 A1 A2
inwin <- rs[pos38 >= POS38 - W & pos38 <= POS38 + W, SNP]

# ---------------- assemble region data: list of data.tables (SNP ea oa beta se p eaf), with type / N / s
TR <- list(); META <- list()
add <- function(name, d, type, N, s = NA_real_, label = name) {
  d[, SNP := as.character(SNP)]
  d <- unique(d[grepl("^rs", SNP) & !is.na(beta) & !is.na(se) & se > 0], by = "SNP")
  TR[[name]] <<- d; META[[name]] <<- data.table(trait = name, label = label, type = type, N = N, s = s, n_snps = nrow(d))
}
gw <- as.data.table(readRDS(file.path(D0, "gwas_region.rds")))
for (k in c("AMD_FG", "DRYAMD_FG", "WETAMD_FG", "AMD_IAMDGC", "AS", "IS", "LAS", "CES", "SVS", "LACUNAR", "WMH")) {
  d <- gw[trait == k]; nc <- d$ncase[1]; nt <- d$ncontrol[1]; nn <- if (!is.na(nc)) nc + nt else d$n[1]
  add(k, d[, .(SNP, ea, oa, beta, se, p, eaf)], if (is.na(nc)) "quant" else "cc", nn, if (is.na(nc)) NA else nc / (nc + nt))
}
ci <- as.data.table(readRDS(file.path(D0, "gwas_cimt.rds")))
add("CIMT", ci[trait == "CIMT", .(SNP, ea, oa, beta, se, p, eaf)], "quant", 71128)
add("PLAQUE", ci[trait == "PLAQUE", .(SNP, ea, oa, beta, se, p, eaf)], "cc", 48434, 21540 / 48434)
aq <- fread(file.path(D0, "eqtl", "GTEx_artery_aorta.tsv.gz"))[gene_id == "ENSG00000163638"]
add("AORTA_ADAMTS9", aq[, .(SNP = rsid, ea = alt, oa = ref, beta, se, p = pvalue, eaf = ifelse(is.na(an), maf, ac / an))], "quant", 387,
    label = "Aorta ADAMTS9 eQTL (GTEx)")
fgload <- function(f) { d <- fread(file.path(D1, f)); d[, .(SNP, ea, oa, beta, se, p, eaf)] }
fgN <- fread(file.path(TMP, "finngen_counts.tsv"))
for (e in c("T2D", "I9_STR", "C_STROKE", "I9_ICH", "I9_CHD", "I9_CORATHER", "I9_HYPTENS")) {
  r <- fgN[phenocode == e]; add(paste0("FG_", e), fgload(paste0("FG_", e, ".tsv.gz")), "cc", r$n_case + r$n_control, r$n_case / (r$n_case + r$n_control))
}
add("CAD", fgload("CAD_Aragam2022.tsv.gz"), "cc", 181522 + 984168, 181522 / (181522 + 984168), "CAD (Aragam 2022)")
add("LDL", fgload("LDL_GLGC2021.tsv.gz"), "quant", 1320016); add("HDL", fgload("HDL_GLGC2021.tsv.gz"), "quant", 1320016)
for (e in c("AS", "IS", "LAS", "SVS")) {
  cc <- c(AS = 27413, IS = 19032, LAS = 1735, SVS = 5532)[e]
  add(paste0("EAS_", e), fgload(paste0("EAS_", e, ".tsv.gz")), "cc", cc + 237242, cc / (cc + 237242))
}
# downloaded full files: keep region SNPs only
rsel <- function(d) d[SNP %in% rs$SNP]
x <- fread(file.path(CM, "T2D_Xue2018.h.tsv.gz"), select = c("hm_rsid", "hm_effect_allele", "hm_other_allele", "hm_beta", "standard_error", "p_value", "hm_effect_allele_frequency"))
setnames(x, c("SNP", "ea", "oa", "beta", "se", "p", "eaf")); add("T2D", rsel(x), "cc", 62892 + 596424, 62892 / (62892 + 596424), "T2D (Xue 2018)")
x <- fread(file.path(CM, "WHRadjBMI_Pulit2019.tbl.gz"), select = c("SNP", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value"))
setnames(x, c("SNP", "ea", "oa", "eaf", "beta", "se", "p")); x[, SNP := sub(":.*", "", SNP)]; x[, `:=`(ea = toupper(ea), oa = toupper(oa))]
add("WHRadjBMI", rsel(x), "quant", 694649)
for (bp in c("SBP", "DBP")) {
  x <- fread(file.path(CM, paste0(bp, "_Evangelou2018.txt.gz")), select = c("MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P"))
  setnames(x, c("marker", "ea", "oa", "eaf", "beta", "se", "p"))
  x <- x[startsWith(marker, "3:")][, pos37 := as.integer(tstrsplit(marker, ":")[[2]])][, `:=`(ea = toupper(ea), oa = toupper(oa))]
  x <- merge(x, rs[!is.na(pos37), .(SNP, pos37, A1, A2)], by = "pos37")[(ea == A1 & oa == A2) | (ea == A2 & oa == A1)]
  add(bp, x[, .(SNP, ea, oa, beta, se, p, eaf)], "quant", 757601)
}
MEGA <- list()
for (e in c("AS", "IS", "LAS", "SVS", "CES")) {
  f <- file.path(CM, paste0("MEGASTROKE_", e, ".h.tsv.gz")); if (!file.exists(f)) next
  x <- fread(f, select = c("hm_rsid", "hm_effect_allele", "hm_other_allele", "hm_beta", "standard_error", "p_value", "hm_effect_allele_frequency"))
  setnames(x, c("SNP", "ea", "oa", "beta", "se", "p", "eaf")); MEGA[[e]] <- unique(x[SNP %in% c(KEY, "rs6795735", "rs17676309")], by = "SNP")
}
meta <- rbindlist(META); fwrite(meta, file.path(OUT, "traits_used.tsv"), sep = "\t"); print(meta)

# ---------------- 1. key-variant effects aligned to the rs4611812-C allele (raises aortic ADAMTS9)
alignC <- function(d, snp = KEY, al = KEY_C) {
  r <- d[SNP == snp][1]; if (is.na(r$SNP)) return(data.table(beta = NA_real_, se = NA_real_, p = NA_real_))
  s <- if (r$ea == al) 1 else if (r$oa == al) -1 else NA_real_
  data.table(beta = s * r$beta, se = r$se, p = r$p)
}
keyeff <- rbindlist(lapply(names(TR), function(k) cbind(trait = k, label = META[[k]]$label, type = META[[k]]$type, alignC(TR[[k]]))))
keyeff[, OR := fifelse(type == "cc", exp(beta), NA_real_)]
fwrite(keyeff, file.path(OUT, "rs4611812C_effects_all_traits.tsv"), sep = "\t")

# ---------------- 2. pairwise coloc.abf (+/-500 kb) with prior sensitivity
mkds <- function(k) {
  d <- TR[[k]][SNP %in% inwin]; m <- META[[k]]
  maf <- pmin(d$eaf, 1 - d$eaf)
  l <- list(beta = d$beta, varbeta = d$se^2, snp = d$SNP, type = m$type, N = m$N)
  if (m$type == "cc") l$s <- m$s
  if (m$type == "quant") { ok <- !is.na(maf) & maf > 0; l$MAF <- ifelse(is.na(maf) | maf <= 0, 0.01, maf) }
  l
}
align_pair <- function(a, b) {   # restrict to shared SNPs with matching alleles (coloc.abf uses Z, sign irrelevant, alleles only to drop mismatches)
  m <- merge(TR[[a]][SNP %in% inwin, .(SNP, ea, oa)], TR[[b]][SNP %in% inwin, .(SNP, ea2 = ea, oa2 = oa)], by = "SNP")
  m[(ea == ea2 & oa == oa2) | (ea == oa2 & oa == ea2), SNP]
}
SET <- c("AMD_FG", "WETAMD_FG", "DRYAMD_FG", "AMD_IAMDGC", "AORTA_ADAMTS9", "IS", "LAS", "AS", "CIMT", "PLAQUE",
         "T2D", "FG_T2D", "WHRadjBMI", "SBP", "DBP", "LDL", "HDL", "CAD", "FG_I9_CORATHER", "FG_I9_STR")
SET <- intersect(SET, names(TR))
pairs <- combn(SET, 2, simplify = FALSE)
cres <- rbindlist(lapply(pairs, function(pr) {
  sn <- align_pair(pr[1], pr[2]); if (length(sn) < 50) return(NULL)
  d1 <- mkds(pr[1]); d2 <- mkds(pr[2]); f <- function(l) { i <- match(sn, l$snp); for (z in c("beta", "varbeta", "snp", "MAF")) if (!is.null(l[[z]])) l[[z]] <- l[[z]][i]; l }
  d1 <- f(d1); d2 <- f(d2)
  rbindlist(lapply(c(1e-5, 5e-6, 1e-6), function(p12) {
    r <- suppressMessages(coloc.abf(d1, d2, p12 = p12))
    top <- r$results[which.max(r$results$SNP.PP.H4), ]
    data.table(trait1 = pr[1], trait2 = pr[2], p12 = p12, nsnps = r$summary["nsnps"], PP.H0 = r$summary["PP.H0.abf"],
               PP.H1 = r$summary["PP.H1.abf"], PP.H2 = r$summary["PP.H2.abf"], PP.H3 = r$summary["PP.H3.abf"],
               PP.H4 = r$summary["PP.H4.abf"], top_snp_H4 = top$snp, top_snp_PP = top$SNP.PP.H4)
  }))
}))
fwrite(cres, file.path(OUT, "coloc_pairwise_3p14.tsv"), sep = "\t")
message("pairwise coloc: ", nrow(cres) / 3, " pairs")

# ---------------- 3. HyPrColoc (multi-trait), alleles aligned to AMD_FG
hyp_run <- function(traits, tag, ...) {
  ref <- TR[["AMD_FG"]][SNP %in% inwin, .(SNP, ea0 = ea, oa0 = oa)]
  B <- ref[, .(SNP)]; S <- ref[, .(SNP)]
  for (k in traits) {
    m <- merge(ref, TR[[k]][, .(SNP, ea, oa, beta, se)], by = "SNP")
    m[, beta := fifelse(ea == ea0 & oa == oa0, beta, fifelse(ea == oa0 & oa == ea0, -beta, NA_real_))]
    m <- m[!is.na(beta)]
    B <- merge(B, m[, .(SNP, b = beta)], by = "SNP"); S <- merge(S, m[, .(SNP, s = se)], by = "SNP")
    setnames(B, "b", k); setnames(S, "s", k)
  }
  bm <- as.matrix(B[, -1]); sm <- as.matrix(S[, -1]); rownames(bm) <- rownames(sm) <- B$SNP
  bin <- as.integer(sapply(traits, function(k) META[[k]]$type == "cc"))
  r <- hyprcoloc(bm, sm, trait.names = traits, snp.id = B$SNP, binary.outcomes = bin, ...)
  out <- as.data.table(r$results); out[, `:=`(analysis = tag, nsnps = nrow(bm))]; out
}
HT <- list(
  core = c("AMD_FG", "AORTA_ADAMTS9", "IS", "CIMT", "T2D", "WHRadjBMI"),
  extended = c("AMD_FG", "AORTA_ADAMTS9", "IS", "LAS", "CIMT", "T2D", "WHRadjBMI", "SBP", "DBP", "CAD"),
  amd_cardiometabolic = c("AMD_FG", "WETAMD_FG", "T2D", "FG_T2D", "WHRadjBMI", "SBP"),
  stroke_cardiometabolic = c("IS", "LAS", "CIMT", "T2D", "WHRadjBMI", "SBP", "DBP"))
hres <- rbindlist(lapply(names(HT), function(t) {
  tr <- intersect(HT[[t]], names(TR))
  rbindlist(list(
    tryCatch(hyp_run(tr, paste0(t, "_default")), error = function(e) data.table(analysis = t, error = conditionMessage(e))),
    tryCatch(hyp_run(tr, paste0(t, "_prior.c0.05"), prior.c = 0.05), error = function(e) NULL),
    tryCatch(hyp_run(tr, paste0(t, "_reg0.7_align0.7"), reg.thresh = 0.7, align.thresh = 0.7), error = function(e) NULL)), fill = TRUE)
}), fill = TRUE)
fwrite(hres, file.path(OUT, "hyprcoloc_3p14.tsv"), sep = "\t"); print(hres)

# ---------------- 4. FinnGen in-sample SuSiE: coloc.bf_bf and credible-set overlap; aorta eQTL credible set
sus <- lapply(c(AMD = "H7_AMD", WETAMD = "WET_AMD", T2D = "T2D", INSULIN = "KELA_DIAB_INSUL_EXMORE"), function(t) {
  f <- file.path(D1, paste0("susie_fg_", t, ".tsv.gz")); if (file.exists(f)) fread(f) else NULL })
sus <- sus[!sapply(sus, is.null)]
lbfm <- function(d) {      # rows = effects with a credible set, cols = variants (id chr:pos:ref:alt)
  L <- sort(unique(d$cs[d$cs > 0])); if (!length(L)) return(NULL)
  m <- t(as.matrix(d[, paste0("lbf_variable", L), with = FALSE])); colnames(m) <- d$v; rownames(m) <- paste0("L", L); m
}
aq2 <- aq[, .(v = paste(chromosome, position, ref, alt, sep = ":"), z = beta / se, se, maf)]
aq2[, lbf := { V <- se^2; W <- 0.15^2; r <- W / (V + W); 0.5 * (log(1 - r) + r * z^2) }]   # Wakefield ABF, prior sd 0.15
bfb <- list()
for (a in names(sus)) {
  ma <- lbfm(sus[[a]]); if (is.null(ma)) next
  for (b in setdiff(names(sus), a)) {
    if (a > b) next
    mb <- lbfm(sus[[b]]); if (is.null(mb)) next
    r <- tryCatch(coloc.bf_bf(ma, mb), error = function(e) NULL); if (is.null(r)) next
    bfb[[paste(a, b)]] <- data.table(trait1 = paste0("FG_", a), trait2 = paste0("FG_", b), as.data.table(r$summary))
  }
  ae <- matrix(aq2$lbf, nrow = 1, dimnames = list("L1", aq2$v))
  r <- tryCatch(coloc.bf_bf(ma, ae), error = function(e) NULL)
  if (!is.null(r)) bfb[[paste(a, "aorta")]] <- data.table(trait1 = paste0("FG_", a), trait2 = "AORTA_ADAMTS9 (single-effect ABF)", as.data.table(r$summary))
}
bfb <- rbindlist(bfb, fill = TRUE); fwrite(bfb, file.path(OUT, "coloc_bf_bf_finngen_insample.tsv"), sep = "\t"); print(bfb)
cs_aorta <- fread(file.path(PROJ_ROOT, "data/raw/eqtl/susie_aorta/aorta_ADAMTS9_cs.tsv"))[gene_id == "ENSG00000163638"]
cs_aorta[, v := sub("^chr", "", gsub("_", ":", variant))]
csl <- c(list(AORTA_ADAMTS9 = cs_aorta$v), do.call(c, lapply(names(sus), function(a) {
  d <- sus[[a]][cs > 0]; split(d$v, paste0("FG_", a, "_CS", d$cs)) })))
ov <- rbindlist(lapply(names(csl), function(x) rbindlist(lapply(names(csl), function(y)
  data.table(set1 = x, set2 = y, n1 = length(csl[[x]]), n2 = length(csl[[y]]), n_shared = length(intersect(csl[[x]], csl[[y]])),
             shared = paste(intersect(csl[[x]], csl[[y]]), collapse = ";"))))))
fwrite(ov[set1 < set2], file.path(OUT, "credible_set_overlap.tsv"), sep = "\t")
cs_members <- rbindlist(lapply(names(csl), function(x) data.table(set = x, variant = csl[[x]])))
fwrite(cs_members, file.path(OUT, "credible_set_members.tsv"), sep = "\t")

# ---------------- 5. subtype heterogeneity (GIGASTROKE EUR; shared controls) and stroke vs CAD
gk <- fread(file.path(TMP, "giga_keysnps.tsv"))
ncase <- c(AS = 73652, IS = 62100, CES = 10804, LAS = 6399, SVS = 6811); N0 <- 1234808
rho <- function(a, b) sqrt(ncase[a] * ncase[b] / ((ncase[a] + N0) * (ncase[b] + N0)))   # distinct cases, fully shared controls
het <- rbindlist(lapply(unique(gk$SNP), function(s) {
  d <- gk[SNP == s]; setkey(d, trait); sub <- c("LAS", "SVS", "CES")
  pw <- rbindlist(lapply(combn(sub, 2, simplify = FALSE), function(pr) {
    b1 <- d[pr[1], beta]; b2 <- d[pr[2], beta]; s1 <- d[pr[1], se]; s2 <- d[pr[2], se]; r <- rho(pr[1], pr[2])
    se_d <- sqrt(s1^2 + s2^2 - 2 * r * s1 * s2); z <- (b1 - b2) / se_d
    data.table(SNP = s, effect_allele = d$ea[1], comparison = paste(pr, collapse = " vs "), beta1 = b1, beta2 = b2, diff = b1 - b2,
               se_diff = se_d, z = z, p = 2 * pnorm(-abs(z)))
  }))
  b <- d[sub, beta]; se <- d[sub, se]
  Sg <- diag(se^2); for (i in 1:3) for (j in 1:3) if (i != j) Sg[i, j] <- rho(sub[i], sub[j]) * se[i] * se[j]
  Si <- solve(Sg); one <- rep(1, 3); Q <- as.numeric(t(b) %*% Si %*% b - (t(one) %*% Si %*% b)^2 / (t(one) %*% Si %*% one))
  rbind(pw, data.table(SNP = s, effect_allele = d$ea[1], comparison = "Q (LAS, SVS, CES)", diff = Q, p = pchisq(Q, 2, lower.tail = FALSE)), fill = TRUE)
}), fill = TRUE)
fwrite(het, file.path(OUT, "subtype_heterogeneity.tsv"), sep = "\t"); print(het)
# stroke vs coronary disease at the key variants (same effect allele as GIGASTROKE rows)
cadcmp <- rbindlist(lapply(unique(gk$SNP), function(s) {
  d <- gk[SNP == s & trait == "IS"]; out <- list()
  for (k in c("CAD", "FG_I9_CORATHER", "FG_I9_CHD")) {
    r <- TR[[k]][SNP == s][1]; if (is.na(r$SNP)) next
    sg <- if (r$ea == d$ea) 1 else if (r$oa == d$ea) -1 else NA
    bc <- sg * r$beta; z <- (d$beta - bc) / sqrt(d$se^2 + r$se^2)
    out[[k]] <- data.table(SNP = s, effect_allele = d$ea, IS_beta = d$beta, IS_se = d$se, IS_p = d$p, cad_trait = k, CAD_beta = bc,
                           CAD_se = r$se, CAD_p = r$p, diff_z = z, diff_p = 2 * pnorm(-abs(z)), same_direction = sign(d$beta) == sign(bc))
  }
  rbindlist(out)
}))
fwrite(cadcmp, file.path(OUT, "stroke_vs_cad_key_variants.tsv"), sep = "\t"); print(cadcmp)
# dry vs wet AMD (FinnGen; overlapping cases and shared controls -> LDSC intercept 0.361 as correlation)
dw <- rbindlist(lapply(c(KEY, "rs6795735", "rs3750848", "rs2672592"), function(s) {
  a <- TR$DRYAMD_FG[SNP == s][1]; b <- TR$WETAMD_FG[SNP == s][1]; if (is.na(a$SNP) || is.na(b$SNP)) return(NULL)
  sg <- if (b$ea == a$ea) 1 else -1; bb <- sg * b$beta; r <- 0.36089
  se_d <- sqrt(a$se^2 + b$se^2 - 2 * r * a$se * b$se); z <- (a$beta - bb) / se_d
  data.table(SNP = s, effect_allele = a$ea, dry_beta = a$beta, dry_se = a$se, wet_beta = bb, wet_se = b$se, rho = r, diff_z = z, p = 2 * pnorm(-abs(z)))
}))
fwrite(dw, file.path(OUT, "dry_vs_wet_amd.tsv"), sep = "\t"); print(dw)

# ---------------- 6. replication of the stroke association at rs4611812 (C allele)
rep <- list()
gC <- function(tr) { d <- gk[SNP == KEY & trait == tr]; s <- if (d$ea == KEY_C) 1 else -1; c(s * d$beta, d$se) }
for (e in names(MEGA)) {
  m <- alignC(MEGA[[e]]); g <- gC(e)
  wG <- 1 / g[2]^2; wM <- 1 / m$se^2
  if (wG > wM) { bN <- (g[1] * wG - m$beta * wM) / (wG - wM); sN <- sqrt(1 / (wG - wM)) } else { bN <- NA; sN <- NA }
  rep[[e]] <- data.table(dataset = c(paste0("GIGASTROKE EUR ", e), paste0("MEGASTROKE EUR ", e), paste0("GIGASTROKE minus MEGASTROKE ", e, " (approx.)")),
                         beta = c(g[1], m$beta, bN), se = c(g[2], m$se, sN))
}
for (k in c("FG_I9_STR", "FG_C_STROKE", "FG_I9_ICH", "EAS_AS", "EAS_IS", "EAS_LAS", "EAS_SVS", "FG_I9_CORATHER", "FG_I9_CHD", "CAD")) {
  if (is.null(TR[[k]])) next; m <- alignC(TR[[k]]); rep[[k]] <- data.table(dataset = k, beta = m$beta, se = m$se)
}
rep <- rbindlist(rep); rep[, `:=`(OR = exp(beta), p = 2 * pnorm(-abs(beta / se)))]
fwrite(rep, file.path(OUT, "replication_rs4611812C.tsv"), sep = "\t"); print(rep)

# ---------------- 7. tissue-scan disclosure
cov <- fread(file.path(PROJ_ROOT, "results/adamts9/eqtl_coverage_3p14.tsv"))
ca <- cov[gene == "ADAMTS9"][order(min_p)]; ca[, rank := .I]
ca[, p_bonf_datasets := pmin(1, min_p * .N)]
fwrite(ca, file.path(OUT, "tissue_scan_ADAMTS9.tsv"), sep = "\t")
message("datasets scanned: ", uniqueN(cov$dataset), "; ADAMTS9 tested in ", nrow(ca), "; aorta rank ", ca[dataset == "GTEx_artery_aorta", rank])
message("locus step done")
