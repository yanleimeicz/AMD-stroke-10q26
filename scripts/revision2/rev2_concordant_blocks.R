# Round-2 minor revision, item 4 (DA N5): follow-up of the concordant LAVA blocks.
# For every concordant block/pair with within-trait FDR < 0.05 (lava_bivariate_with_global_fdr.tsv, rho > 0):
#   lead variants of each trait in the block and their association with the other trait; genome-wide significant hits;
#   coloc.abf (default priors p12 = 1e-5, and 5e-6); direction concordance of variants with P < 1e-4 / 1e-3;
#   presence of AMD MR instruments (A_dist) in the block. Blocks are GRCh37 (1000G EUR bim); traits are matched by rsID.
# Output: results/revision2/concordant_blocks_followup.tsv, concordant_blocks_leads.tsv
# Usage: Rscript rev2_concordant_blocks.R
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(coloc); library(parallel) })
setDTthreads(8)
DER <- file.path(PROJ_ROOT, "derived/gwas"); REF <- file.path(PROJ_ROOT, "derived/lava/ref"); O <- file.path(PROJ_ROOT, "results/revision2")
lava <- fread(file.path(PROJ_ROOT, "results/revision/mr/lava_bivariate_with_global_fdr.tsv"))
LABEL <- c(`699` = "4q25 (downstream of PITX2)", `1596` = "10q26.13 FGFR2-TACC2 (adjacent to ARMS2/HTRA1 block 1597)",
           `1845` = "12q24.21 MED13L", `1873` = "13q12.3 FLT1/POMP", `2159` = "16q23.2 WWOX-MAF",
           `2163` = "16q23.3 CDH13", `2453` = "21q22.2 ERG")
pairs <- lava[locus %in% as.integer(names(LABEL)) & rho > 0 & p_fdr < 0.05, .(locus, chr, start, stop, partner = phen1, amd = phen2, rho, rho.lower, rho.upper, p, p_fdr, fdr_global)]
print(pairs)
traits <- unique(c(pairs$partner, pairs$amd))
L <- lapply(setNames(traits, traits), function(k) {
  d <- readRDS(file.path(DER, paste0(k, ".rds")))
  if (!"ncase" %in% names(d)) d[, ncase := NA_real_]
  if (!"eaf" %in% names(d)) d[, eaf := NA_real_]
  d[, .(SNP, chr, pos, ea, oa, beta, se, eaf, n, ncase, p)]
})
CC <- c("AMD_FG", "DRYAMD_FG", "WETAMD_FG", "AMD_IAMDGC", "AS", "IS", "SVS", "LAS", "CES", "LACUNAR")
inst <- fread(file.path(PROJ_ROOT, "results/revision/mr/instruments_all_schemes.tsv"))[scheme == "A_dist"]

ds <- function(k, d) {
  eaf <- d$eaf; eaf[is.na(eaf)] <- 0.3
  maf <- pmin(eaf, 1 - eaf); maf[maf <= 0] <- 1e-4
  x <- list(beta = d$beta, varbeta = d$se^2, snp = d$SNP, MAF = maf, N = d$n[1])
  if (k %in% CC) { x$type <- "cc"; x$s <- d$ncase[1] / d$n[1] } else { x$type <- "quant"; x$sdY <- 1 }
  x
}

run_pair <- function(i) {
  pr <- pairs[i]
  bim <- fread(file.path(REF, paste0("g1000_eur_chr", pr$chr, ".bim")), header = FALSE, select = c(2, 4), col.names = c("SNP", "pos19"))
  snps <- bim[pos19 >= pr$start & pos19 <= pr$stop]
  a <- L[[pr$amd]][SNP %in% snps$SNP]; b <- L[[pr$partner]][SNP %in% snps$SNP]
  m <- merge(a, b, by = "SNP", suffixes = c(".a", ".b"))
  flip <- m$ea.a == m$oa.b & m$oa.a == m$ea.b
  keep <- (m$ea.a == m$ea.b & m$oa.a == m$oa.b) | flip
  m[flip, beta.b := -beta.b]; m <- m[keep]
  m <- merge(m, snps, by = "SNP")
  m[, `:=`(z.a = beta.a / se.a, z.b = beta.b / se.b)]
  A <- m[, .(SNP, beta = beta.a, se = se.a, eaf = eaf.a, n = n.a, ncase = ncase.a)]
  B <- m[, .(SNP, beta = beta.b, se = se.b, eaf = ifelse(is.na(eaf.b), eaf.a, eaf.b), n = n.b, ncase = ncase.b)]
  c1 <- suppressMessages(coloc.abf(ds(pr$amd, A), ds(pr$partner, B)))
  c2 <- suppressMessages(coloc.abf(ds(pr$amd, A), ds(pr$partner, B), p12 = 5e-6))
  top <- as.data.table(c1$results)[order(-SNP.PP.H4)][1]
  la <- m[which.min(p.a)]; lb <- m[which.min(p.b)]
  conc <- function(th) { x <- m[p.a < th]; if (!nrow(x)) return(c(NA, 0)); c(mean(sign(x$z.a) == sign(x$z.b)), nrow(x)) }
  c4 <- conc(1e-4); c3 <- conc(1e-3)
  zz <- m[p.a < 1e-3 | p.b < 1e-3]
  ins <- inst[exposure == pr$amd & SNP %in% m$SNP]
  data.table(locus = pr$locus, region = LABEL[as.character(pr$locus)], chr = pr$chr, start_hg19 = pr$start, stop_hg19 = pr$stop,
             amd = pr$amd, partner = pr$partner, lava_rho = pr$rho, lava_rho_lo = pr$rho.lower, lava_rho_hi = pr$rho.upper,
             lava_p = pr$p, lava_fdr_within = pr$p_fdr, lava_fdr_global = pr$fdr_global, nsnps = nrow(m),
             amd_lead = la$SNP, amd_lead_pos19 = la$pos19, amd_lead_p = la$p.a, amd_lead_z = la$z.a, amd_lead_p_in_partner = la$p.b, amd_lead_z_in_partner = la$z.b,
             partner_lead = lb$SNP, partner_lead_pos19 = lb$pos19, partner_lead_p = lb$p.b, partner_lead_z = lb$z.b, partner_lead_p_in_amd = lb$p.a, partner_lead_z_in_amd = lb$z.a,
             n_gws_amd = sum(m$p.a < 5e-8), n_gws_partner = sum(m$p.b < 5e-8),
             n_amd_instruments_A_dist = nrow(ins), amd_instruments = paste(ins$SNP, collapse = ";"),
             concord_frac_pA_1e4 = c4[1], n_pA_1e4 = c4[2], concord_frac_pA_1e3 = c3[1], n_pA_1e3 = c3[2],
             zcor_p1e3_either = if (nrow(zz) > 2) cor(zz$z.a, zz$z.b) else NA_real_, n_p1e3_either = nrow(zz),
             PP.H0 = c1$summary[["PP.H0.abf"]], PP.H1 = c1$summary[["PP.H1.abf"]], PP.H2 = c1$summary[["PP.H2.abf"]],
             PP.H3 = c1$summary[["PP.H3.abf"]], PP.H4 = c1$summary[["PP.H4.abf"]], PP.H4_p12_5e6 = c2$summary[["PP.H4.abf"]],
             PP.H3_p12_5e6 = c2$summary[["PP.H3.abf"]], top_snp_H4 = top$snp, top_snp_PP = top$SNP.PP.H4)
}
res <- rbindlist(lapply(seq_len(nrow(pairs)), run_pair), fill = TRUE)
fwrite(res, file.path(O, "concordant_blocks_followup.tsv"), sep = "\t")

# top 5 variants per trait per block (for gene annotation)
leads <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  pr <- pairs[i]
  bim <- fread(file.path(REF, paste0("g1000_eur_chr", pr$chr, ".bim")), header = FALSE, select = c(2, 4), col.names = c("SNP", "pos19"))
  snps <- bim[pos19 >= pr$start & pos19 <= pr$stop]
  rbindlist(lapply(c(pr$amd, pr$partner), function(k) {
    d <- merge(L[[k]][SNP %in% snps$SNP, .(SNP, chr, pos, ea, oa, beta, se, p)], snps, by = "SNP")[order(p)][1:5]
    d[, `:=`(locus = pr$locus, pair = paste(pr$amd, pr$partner), trait = k)]
  }))
}))
fwrite(unique(leads), file.path(O, "concordant_blocks_leads.tsv"), sep = "\t")
print(res[, .(locus, region, amd, partner, lava_rho = round(lava_rho, 2), nsnps, amd_lead_p = signif(amd_lead_p, 2), amd_lead_p_in_partner = signif(amd_lead_p_in_partner, 2),
              partner_lead_p = signif(partner_lead_p, 2), partner_lead_p_in_amd = signif(partner_lead_p_in_amd, 2), n_gws_amd, n_gws_partner,
              n_amd_instruments_A_dist, concord_frac_pA_1e4, n_pA_1e4, H3 = round(PP.H3, 3), H4 = round(PP.H4, 3))])
