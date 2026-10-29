# Revision C, step 2: (i) effect-direction concordance with one variant per locus (distance-pruned AMD leads) and an
# LD-preserving block sign-flip null; (ii) LAVA results by unique LD block; (iii) one FDR family per analytical layer
# and a table of every result quoted in the manuscript with its family-wise FDR.
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
setDTthreads(30); set.seed(20261007)
DER <- file.path(PROJ_ROOT, "derived/gwas"); DRV <- file.path(PROJ_ROOT, "derived/revision"); R <- file.path(PROJ_ROOT, "results")
OUT <- file.path(R, "revision", "mr"); dir.create(file.path(OUT, "plots"), recursive = TRUE, showWarnings = FALSE)
path_of <- function(k) { f <- file.path(DRV, paste0(k, ".rds")); if (file.exists(f)) f else file.path(DER, paste0(k, ".rds")) }

# ------------------------------------------------------------------ (i) concordance
inst <- fread(file.path(OUT, "instruments_all_schemes.tsv"))[scheme == "A_dist" & exposure %in% c("AMD_FG", "AMD_IAMDGC")]
inst[, `:=`(risk = ifelse(beta > 0, ea, oa), nonrisk = ifelse(beta > 0, oa, ea), z_amd = abs(beta / se))]
TR <- c("IS", "SVS", "LAS", "AS", "LACUNAR", "WMH", "WM_MD", "WM_FA", "EPVS_WM", "EPVS_BG", "EPVS_HIP")
NPERM <- 20000
conc <- list(); vlev <- list()
for (t in TR) {
  d <- readRDS(path_of(t))[SNP %in% inst$SNP, .(SNP, ea_o = ea, oa_o = oa, b_o = beta, se_o = se, p_o = p)]
  for (e in unique(inst$exposure)) {
    x <- merge(inst[exposure == e], d, by = "SNP")
    x[, al := fifelse(ea_o == risk & oa_o == nonrisk, 1, fifelse(ea_o == nonrisk & oa_o == risk, -1, NA_real_))]
    x <- x[!is.na(al)]; x[, z_o := al * b_o / se_o]
    if (t == "WM_FA") x[, z_o := -z_o]                         # higher FA = healthier; reverse so that > 0 = worse vascular trait
    x <- x[order(chr, pos)]
    # blocks: leads closer than 5 Mb on the same chromosome share a sign flip (residual long-range LD, e.g. CFH, MHC)
    x[, block := cumsum(c(1, as.integer(diff(pos) > 5e6 | chr[-1] != chr[-.N])))]
    obs_frac <- mean(x$z_o > 0); obs_S <- sum(x$z_o)
    nb <- max(x$block); fl <- matrix(sample(c(-1, 1), nb * NPERM, TRUE), nb, NPERM)
    S_null <- colSums(x$z_o * fl[x$block, , drop = FALSE])
    F_null <- colMeans((x$z_o * fl[x$block, , drop = FALSE]) > 0)
    conc[[paste(e, t)]] <- data.table(exposure = e, trait = t, n_loci = nrow(x), n_blocks = nb, n_same = sum(x$z_o > 0), frac_same = obs_frac,
      binom_p = binom.test(sum(x$z_o > 0), nrow(x))$p.value,
      perm_p_frac = (sum(abs(F_null - 0.5) >= abs(obs_frac - 0.5)) + 1) / (NPERM + 1),
      sum_z = obs_S, perm_p_sumz = (sum(abs(S_null) >= abs(obs_S)) + 1) / (NPERM + 1))
    vlev[[paste(e, t)]] <- x[, .(exposure = e, trait = t, SNP, chr, pos, z_amd, z_o, block)]
  }
  rm(d); gc()
}
conc <- rbindlist(conc); conc[, fdr_perm_sumz := p.adjust(perm_p_sumz, "BH")]
fwrite(conc, file.path(OUT, "concordance_distance_pruned.tsv"), sep = "\t"); fwrite(rbindlist(vlev), file.path(OUT, "concordance_variant_level.tsv"), sep = "\t")
print(conc)

# ------------------------------------------------------------------ (ii) LAVA by unique block
lv <- fread(file.path(R, "csvd", "lava_bivariate.tsv"))
lv[, fdr_global := p.adjust(p, "BH")]                       # one family over all bivariate tests
blk <- lv[, .(n_tests = .N, n_fdr_within_trait = sum(p_fdr < 0.05), n_fdr_global = sum(fdr_global < 0.05),
              n_neg_sig = sum(p_fdr < 0.05 & rho < 0), n_pos_sig = sum(p_fdr < 0.05 & rho > 0),
              best_pair = paste(phen1, phen2)[which.min(p)], best_rho = rho[which.min(p)], best_p = min(p),
              sign_consistent_sig = if (any(p_fdr < 0.05)) length(unique(sign(rho[p_fdr < 0.05]))) == 1 else NA),
          by = .(locus, chr, start, stop, region)][n_fdr_within_trait > 0 | n_fdr_global > 0][order(best_p)]
fwrite(lv, file.path(OUT, "lava_bivariate_with_global_fdr.tsv"), sep = "\t"); fwrite(blk, file.path(OUT, "lava_unique_blocks.tsv"), sep = "\t")
print(blk)

# ------------------------------------------------------------------ (iii) FDR layers
lay <- list()
mr0 <- fread(file.path(R, "mr", "mr_ivw_main.tsv"))
fw <- mr0[exposure %in% c("AMD_FG", "DRYAMD_FG", "WETAMD_FG", "AMD_IAMDGC")]; rv <- mr0[!exposure %in% c("AMD_FG", "DRYAMD_FG", "WETAMD_FG", "AMD_IAMDGC")]
lay$mr_forward_orig <- fw[, .(layer = "MR forward (original instruments)", tier = "primary", item = paste(exposure, "->", outcome), estimate = OR, p = pval)]
lay$mr_reverse_orig <- rv[, .(layer = "MR reverse (original instruments)", tier = "primary", item = paste(exposure, "->", outcome), estimate = OR, p = pval)]
rr <- fread(file.path(OUT, "mr_revised_all.tsv"))
a <- rr[scheme == "A_dist" & method == "IVW" & outcome %in% c("AS", "IS", "LAS", "CES", "SVS", "WMH")]
lay$mr_forward_rev <- a[exclusion == "main", .(layer = "MR forward (distance-pruned instruments)", tier = "primary", item = paste(exposure, "->", outcome), estimate = exp(b), p)]
lay$mr_sens <- a[exclusion != "main", .(layer = "MR locus-exclusion sensitivity", tier = "exploratory", item = paste(exposure, "->", outcome, "|", exclusion), estimate = exp(b), p)]
lay$mr_mega <- rr[scheme == "A_dist" & method == "IVW" & grepl("^MEGA_", outcome) & exclusion %in% c("main", "no_major3", "no_major3_lipid"),
                  .(layer = "MR replication, MEGASTROKE (no FinnGen)", tier = "replication", item = paste(exposure, "->", outcome, "|", exclusion), estimate = exp(b), p)]
cs <- fread(file.path(R, "csvd", "mr_ivw.tsv"))[analysis == "main" & !outcome %in% c("SVS", "IS", "WMH")]
lay$mr_csvd <- cs[, .(layer = "MR to CSVD extension traits", tier = "exploratory", item = paste(exposure, "->", outcome), estimate = OR, p = pval)]
l1 <- fread(file.path(R, "ldsc", "ldsc_rg_all.tsv")); l2 <- fread(file.path(R, "csvd", "ldsc_rg.tsv"))
ld <- unique(rbind(l2[, .(trait1, trait2, rg, p)], l1[, .(trait1, trait2, rg, p)]), by = c("trait1", "trait2"))
ld <- ld[trait1 %in% c("AMD_FG", "DRYAMD_FG", "WETAMD_FG") & !trait2 %in% c("AMD_FG", "DRYAMD_FG", "WETAMD_FG", "AMD_IAMDGC")]
lay$ldsc <- ld[, .(layer = "LDSC genetic correlation (FinnGen AMD)", tier = "primary", item = paste(trait1, "~", trait2), estimate = rg, p)]
lay$lava <- lv[, .(layer = "LAVA bivariate local rg", tier = "primary", item = paste0(phen2, " ~ ", phen1, " | block ", locus, ifelse(is.na(region), "", paste0(" ", region))), estimate = rho, p)]
kv <- fread(file.path(R, "csvd", "key_variants_AMD_FG.tsv"))
lay$keyvar <- kv[, .(layer = "Key-variant lookups", tier = "exploratory", item = paste(SNP, locus, trait), estimate = b_csvd, p = p_csvd)]
sm <- fread(file.path(R, "smr_heidi", "smr_heidi_all.tsv"))
lay$smr <- sm[, .(layer = paste0("SMR eQTLGen (", outcome, ")"), tier = "exploratory", item = paste(GeneSymbol, "->", outcome), estimate = OR_SMR, p = p_SMR)]
pq <- fread(file.path(R, "pqtl_mr", "pqtl_mr_all.tsv"))[source == "UKB-PPP"]
lay$pqtl <- pq[, .(layer = paste0("pQTL MR UKB-PPP (", outcome, ")"), tier = "exploratory", item = paste(gene, "->", outcome), estimate = OR, p = pval)]
L <- rbindlist(lay, fill = TRUE)
L[, `:=`(n_tests = .N, fdr = p.adjust(p, "BH")), by = layer]
fwrite(L, file.path(OUT, "fdr_all_layers.tsv"), sep = "\t")

# results quoted in the manuscript (Results, Discussion, Table 2, Abstract)
pick <- function(lyr, pat) L[grepl(lyr, layer, fixed = TRUE) & grepl(pat, item)]
cited <- rbindlist(list(
  pick("MR forward (original instruments)", "^(AMD_FG|WETAMD_FG|DRYAMD_FG|AMD_IAMDGC) -> (IS|SVS|LAS|AS)$"),
  pick("MR reverse (original instruments)", "^CES -> "),
  pick("MR forward (distance-pruned instruments)", "."),
  pick("MR locus-exclusion sensitivity", "^DRYAMD_FG -> (AS|IS|SVS) \\| (no_major3|no_major3_lipid)$"),
  pick("MR locus-exclusion sensitivity", "^(AMD_FG|WETAMD_FG) -> (IS|SVS) \\| no_major3$"),
  pick("MR locus-exclusion sensitivity", "-> WMH \\| no_major3$"),
  pick("LDSC genetic correlation", "WMH|WM_MD"),
  L[layer == "LAVA bivariate local rg" & p < 0.05 & fdr < 0.2],
  pick("Key-variant lookups", "^(rs3184504|rs2672592|rs3750848|rs429358|rs1269852|rs10737680) "),
  L[grepl("^SMR eQTLGen", layer) & grepl("^(SLC12A3|HSPH1|CSK|EDEM2|PLTP|SH2B3|C4A|C4B) ->", item)],
  L[grepl("^pQTL MR", layer) & grepl("^(APOE|MMP12|WARS1|AGER) ->", item)]), fill = TRUE)
cited[, status := fifelse(fdr < 0.05, "FDR < 0.05", fifelse(p < 0.05, "nominal only (exploratory)", "not significant"))]
fwrite(cited, file.path(OUT, "fdr_cited_results.tsv"), sep = "\t")
print(cited[status != "FDR < 0.05" & p < 0.05, .(layer, item, estimate = signif(estimate, 3), p = signif(p, 3), n_tests, fdr = signif(fdr, 3))])
message("rev_mr_02 done")
