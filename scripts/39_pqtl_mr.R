# Plasma proteome cis-pQTL MR on AMD, stroke subtypes and WMH (Funct Integr Genomics 2026 xQTL design).
# Exposures: UKB-PPP (Sun 2023 Nature, ST9 discovery primary pQTLs, p < 1.7e-11, Olink, EUR n = 34,557) = discovery;
#            deCODE (Ferkingstad 2021 Nat Genet, ST2 primary cis signals, SomaScan, n = 35,559) = replication.
# Instruments: cis only (variant within 1 Mb of the encoding gene, as flagged by the source tables), F > 10;
#   proteins with > 1 cis lead are LD-clumped (r2 < 0.001, 1000G EUR). Missing outcome SNPs -> 1000G EUR proxy (r2 > 0.8, in-phase alleles).
# MR: Wald ratio (1 IV) / IVW multiplicative random effects (>= 2 IVs); BH-FDR within source x outcome; Steiger directionality.
# Outputs -> results/pqtl_mr/
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(TwoSampleMR); library(ggplot2); library(parallel) })
setDTthreads(16); set.seed(20261006)
DER <- file.path(PROJ_ROOT, "derived/gwas"); PQ <- file.path(PROJ_ROOT, "derived/pqtl"); OUT <- file.path(PROJ_ROOT, "results/pqtl_mr")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PLINK <- file.path(PROJ_ROOT, "data/raw/tools/plink/plink.exe"); BFILE <- file.path(PROJ_ROOT, "data/raw/ref/g1000_eur/g1000_eur")
TMP <- file.path(PROJ_ROOT, "tmp/pqtl"); dir.create(TMP, recursive = TRUE, showWarnings = FALSE)
AMD <- c("AMD_FG", "AMD_IAMDGC", "DRYAMD_FG", "WETAMD_FG"); STK <- c("AS", "IS", "LAS", "CES", "SVS", "WMH")
TRAITS <- c(AMD, STK)
PREV <- c(AMD_IAMDGC = 0.02, AMD_FG = 0.03, DRYAMD_FG = 0.02, WETAMD_FG = 0.015, AS = 0.03, IS = 0.025, CES = 0.005, LAS = 0.005, SVS = 0.005)
N_UKB <- 34557; N_DEC <- 35559

# Regions of interest (GRCh38) for annotation
REG <- data.table(region = c("xMHC", "10q26 ARMS2/HTRA1", "12q24 SH2B3", "1q31 CFH", "19p13 C3", "9q34 ABO", "19q13 APOE"),
                  chr = c("6", "10", "12", "1", "19", "9", "19"),
                  start = c(25500000, 121454932, 110884000, 195651878, 5677704, 132136000, 43905796),
                  end = c(33500000, 123514908, 112946000, 197747504, 7720662, 134160000, 45909393))
flag_region <- function(vchr, vpos) {
  vapply(seq_along(vchr), function(i) {
    h <- REG$region[REG$chr == vchr[i] & vpos[i] >= REG$start & vpos[i] <= REG$end]; if (length(h)) h[1] else ""
  }, character(1))
}

# ---------------- 1. instruments
u <- fread(file.path(PQ, "ukbppp_ST9_primary.tsv"), sep = "\t", quote = "")
setnames(u, c("vid", "chr", "pos", "region_id", "rstart", "rend", "mhc", "protein_id", "target", "uniprot", "SNP", "eaf",
              "beta", "se", "logp", paste0("rep_", 1:4), "cistrans", "cis_gene", "annot_gene", names(u)[23:ncol(u)]))
u <- u[cistrans == "cis" & grepl("^rs[0-9]+$", SNP)]
al <- tstrsplit(u$vid, ":", fixed = TRUE)
u[, `:=`(oa = al[[3]], ea = al[[4]], chr = as.character(chr))]
u[, gene := sub(":.*", "", protein_id)]
u[, `:=`(n = N_UKB, F = (beta / se)^2, p = 10^(-logp))]
u <- u[F > 10 & logp >= -log10(5e-8) & ea %in% c("A", "C", "G", "T") & oa %in% c("A", "C", "G", "T")]

dc <- fread(file.path(PQ, "decode_ST02_pqtl.tsv"), sep = "\t", quote = "")
cn <- names(dc)
pick <- function(pat) cn[grepl(pat, cn, fixed = TRUE)][1]
dc <- dc[, .(gene = get(pick("gene")), seqid = SeqId, short = get(pick("shortname")), SNP = variant,
             chr = sub("chr", "", get(pick("chr (var.)"))), pos = get(pick("pos (var.)")),
             ea = Amin, oa = Amaj, maf = as.numeric(get(pick("MAF"))) / 100, cistrans = get(pick("cis/")),
             beta = as.numeric(get(pick("beta (adj.)"))), logp = as.numeric(get(pick("-Log10(P) (adj.)"))),
             rank = as.integer(get(pick("Rank (cond.")))) ]
dc <- dc[cistrans == "cis" & rank == 1 & grepl("^rs[0-9]+$", SNP) & ea %in% c("A", "C", "G", "T") & oa %in% c("A", "C", "G", "T") &
            is.finite(beta) & beta != 0]   # deCODE codes deletions as "!"
dc[, z := qnorm(-logp * log(10) - log(2), lower.tail = FALSE, log.p = TRUE)]    # two-sided P -> |z| without underflow
dc[, `:=`(se = abs(beta) / z, eaf = maf, n = N_DEC, p = 10^(-logp), protein_id = paste0(short, "_", seqid))]
dc[, F := (beta / se)^2]
dc <- dc[F > 10 & logp >= -log10(5e-8)]

clump_multi <- function(d) {     # proteins with > 1 cis lead: keep LD-independent leads (r2 < 0.001)
  multi <- d[, .N, by = protein_id][N > 1, protein_id]
  if (!length(multi)) return(d)
  keep <- rbindlist(lapply(multi, function(pid) {
    x <- d[protein_id == pid]
    cl <- tryCatch(ieugwasr::ld_clump(dplyr::tibble(rsid = x$SNP, pval = x$p), clump_kb = 10000, clump_r2 = 0.001, clump_p = 1,
                                      plink_bin = PLINK, bfile = BFILE), error = function(e) NULL)
    if (is.null(cl)) x[which.min(p)] else x[SNP %in% cl$rsid]
  }))
  rbind(d[!protein_id %in% multi], keep)
}
u <- clump_multi(u); dc <- clump_multi(dc)
u[, region := flag_region(chr, pos)]; dc[, region := flag_region(chr, pos)]
message(sprintf("UKB-PPP: %d proteins, %d cis SNPs; deCODE: %d proteins, %d cis SNPs",
                uniqueN(u$protein_id), nrow(u), uniqueN(dc$protein_id), nrow(dc)))

# allele coding check between the two platforms (same gene, same lead SNP)
cmp <- merge(u[, .(gene, SNP, ea_u = ea, b_u = beta)], dc[, .(gene, SNP, ea_d = ea, oa_d = oa, b_d = beta)], by = c("gene", "SNP"))
cmp[, b_d_al := fifelse(ea_d == ea_u, b_d, fifelse(oa_d == ea_u, -b_d, NA_real_))]
message(sprintf("platform check: %d shared gene-SNP leads, sign concordance %.1f%%", nrow(cmp), 100 * mean(sign(cmp$b_u) == sign(cmp$b_d_al), na.rm = TRUE)))

# ---------------- 2. proxies (1000G EUR, r2 > 0.8, phased alleles) for every instrument SNP
iv_snps <- unique(c(u$SNP, dc$SNP))
writeLines(iv_snps, file.path(TMP, "iv_snps.txt"))
system2(PLINK, c("--bfile", BFILE, "--r2", "in-phase", "--ld-snp-list", file.path(TMP, "iv_snps.txt"), "--ld-window-kb", "500",
                 "--ld-window", "99999", "--ld-window-r2", "0.8", "--threads", "30", "--out", file.path(TMP, "proxies")), stdout = FALSE, stderr = FALSE)
px <- fread(file.path(TMP, "proxies.ld"))[SNP_A != SNP_B]
ph <- tstrsplit(px$PHASE, "/", fixed = TRUE)          # "A1B1/A2B2": target allele + proxy allele in phase
px[, `:=`(tA1 = substr(ph[[1]], 1, 1), pB1 = substr(ph[[1]], 2, 2), tA2 = substr(ph[[2]], 1, 1), pB2 = substr(ph[[2]], 2, 2))]
px <- px[nchar(ph[[1]]) == 2 & nchar(ph[[2]]) == 2][order(SNP_A, -R2)]
message(sprintf("proxies: %d target SNPs with >= 1 proxy", uniqueN(px$SNP_A)))

# ---------------- 3. outcome lookups (parallel, one worker per GWAS)
need <- unique(c(iv_snps, px$SNP_B))
cache_f <- file.path(TMP, "outcome_lookup_cache.rds")
if (file.exists(cache_f) && all(need %in% readRDS(cache_f)$need)) {   # cached lookup covers all SNPs needed
  OUTC <- readRDS(cache_f)$OUTC
} else {
  cl <- makeCluster(length(TRAITS)); clusterExport(cl, c("DER", "need"))
  OUTC <- parLapply(cl, TRAITS, function(k) {
    suppressPackageStartupMessages(library(data.table))
    x <- readRDS(file.path(DER, paste0(k, ".rds")))[SNP %in% need]; x[, trait := k]; x
  })
  stopCluster(cl); names(OUTC) <- TRAITS
  saveRDS(list(need = need, OUTC = OUTC), cache_f)
}

outcome_for <- function(k, snps) {           # direct hits + proxies re-coded to the target SNP alleles
  o <- OUTC[[k]]; hit <- o[SNP %in% snps]; hit[, `:=`(proxy = "", proxy_r2 = 1)]
  miss <- setdiff(snps, hit$SNP)
  if (length(miss)) {
    cand <- px[SNP_A %in% miss & SNP_B %in% o$SNP][, .SD[1], by = SNP_A]
    if (nrow(cand)) {
      po <- merge(cand, o, by.x = "SNP_B", by.y = "SNP")
      po[, ok := (ea == pB1 & oa == pB2) | (ea == pB2 & oa == pB1)]
      po <- po[ok == TRUE]
      po[, `:=`(ea_t = fifelse(ea == pB1, tA1, tA2), oa_t = fifelse(ea == pB1, tA2, tA1))]
      hit <- rbind(hit, po[, .(SNP = SNP_A, chr, pos, ea = ea_t, oa = oa_t, eaf, beta, se, p, n,
                              ncase = if ("ncase" %in% names(po)) ncase else NA_real_,
                              ncontrol = if ("ncontrol" %in% names(po)) ncontrol else NA_real_,
                              trait, proxy = SNP_B, proxy_r2 = R2)], fill = TRUE)
    }
  }
  hit
}

# ---------------- 4. MR per source x outcome (parallel over outcomes)
add_r2 <- function(dat, k) {   # variance explained (binary outcomes on the liability scale) for Steiger
  eaf <- dat$eaf.exposure; eaf[is.na(eaf)] <- 0.3
  re <- get_r_from_bsen(dat$beta.exposure, dat$se.exposure, dat$samplesize.exposure)
  if (k %in% names(PREV)) {
    eo <- dat$eaf.outcome; eo[is.na(eo)] <- eaf[is.na(eo)]
    ro <- get_r_from_lor(dat$beta.outcome, eo, dat$ncase.outcome[1], dat$ncontrol.outcome[1], PREV[[k]])
  } else ro <- get_r_from_bsen(dat$beta.outcome, dat$se.outcome, dat$samplesize.outcome)
  dat$r.exposure <- re; dat$r.outcome <- ro; dat$rsq.exposure <- re^2; dat$rsq.outcome <- ro^2
  dat
}

mr_one <- function(k, exp, src) {      # one outcome: harmonise all proteins, Wald ratio / IVW-MRE, Steiger
  o <- as.data.frame(outcome_for(k, unique(exp$SNP))); o$Phenotype <- rep(k, nrow(o))
  if (!"ncase" %in% names(o)) { o$ncase <- NA_real_; o$ncontrol <- NA_real_ }
  out <- format_data(o, type = "outcome", snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "eaf",
                     effect_allele_col = "ea", other_allele_col = "oa", pval_col = "p", samplesize_col = "n",
                     ncase_col = "ncase", ncontrol_col = "ncontrol", phenotype_col = "Phenotype")
  out$id.outcome <- k
  prox <- data.table(SNP = o$SNP, proxy = o$proxy)
  dat <- suppressMessages(harmonise_data(exp, out, action = 2)); dat <- dat[dat$mr_keep, ]
  if (!nrow(dat)) return(NULL)
  r <- suppressMessages(mr(dat, method_list = c("mr_wald_ratio", "mr_ivw_mre")))
  r <- as.data.table(r)[(nsnp == 1 & method == "Wald ratio") | (nsnp >= 2 & grepl("Inverse variance", method))]
  st <- tryCatch(as.data.table(suppressMessages(directionality_test(add_r2(dat, k))))[
                   , .(id.exposure, steiger_dir = correct_causal_direction, steiger_p = steiger_pval)],
                 error = function(e) data.table(id.exposure = character(), steiger_dir = logical(), steiger_p = numeric()))
  sn <- merge(as.data.table(dat)[, .(id.exposure, SNP)], prox, by = "SNP", all.x = TRUE)
  sn <- sn[, .(snps = paste(SNP, collapse = ";"), proxies = paste(unique(proxy[!is.na(proxy) & proxy != ""]), collapse = ";")), by = id.exposure]
  r <- merge(merge(r, st, by = "id.exposure", all.x = TRUE), sn, by = "id.exposure", all.x = TRUE)
  r[, `:=`(source = src, outcome = k)]
  r
}

run_source <- function(src, iv) {
  ex <- as.data.frame(iv); ex$Phenotype <- ex$protein_id
  exp <- format_data(ex, type = "exposure", snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "eaf",
                     effect_allele_col = "ea", other_allele_col = "oa", pval_col = "p", samplesize_col = "n",
                     phenotype_col = "Phenotype", chr_col = "chr", pos_col = "pos")
  exp$id.exposure <- exp$exposure
  if (nzchar(Sys.getenv("PQTL_SEQ"))) return(rbindlist(lapply(TRAITS, mr_one, exp = exp, src = src), fill = TRUE))
  cl <- makeCluster(length(TRAITS))
  clusterEvalQ(cl, { .libPaths(c(file.path(Sys.getenv("COMORBIDITY_ROOT"), "Rlib"), .libPaths())); suppressPackageStartupMessages({ library(data.table); library(TwoSampleMR) }) })
  clusterExport(cl, c("outcome_for", "OUTC", "px", "add_r2", "PREV", "mr_one"))
  res <- parLapply(cl, TRAITS, mr_one, exp = exp, src = src)
  stopCluster(cl)
  rbindlist(res, fill = TRUE)
}

ann_u <- u[, .(gene = gene[1], region = paste(unique(region[region != ""]), collapse = ";"), F_min = min(F)), by = protein_id]
ann_d <- dc[, .(gene = gene[1], region = paste(unique(region[region != ""]), collapse = ";"), F_min = min(F)), by = protein_id]
res <- rbind(merge(run_source("UKB-PPP", u), ann_u, by.x = "id.exposure", by.y = "protein_id"),
             merge(run_source("deCODE", dc), ann_d, by.x = "id.exposure", by.y = "protein_id"))
setnames(res, "id.exposure", "protein_id")
res[, `:=`(OR = exp(b), lci = exp(b - 1.96 * se), uci = exp(b + 1.96 * se))]
res[, fdr := p.adjust(pval, "BH"), by = .(source, outcome)]
res[, group := fifelse(outcome %in% AMD, "AMD", "STROKE_WMH")]
res <- res[, .(source, protein_id, gene, outcome, group, method, nsnp, b, se, pval, fdr, OR, lci, uci,
               steiger_dir, steiger_p, region, F_min, snps, proxies)][order(source, pval)]
fwrite(res, file.path(OUT, "pqtl_mr_all.tsv"), sep = "\t")
message("MR rows: ", nrow(res))
print(res[fdr < 0.05, .N, by = .(source, outcome)][order(source, outcome)])

# ---------------- 5. shared proteins (UKB-PPP discovery) + deCODE replication
ukb <- res[source == "UKB-PPP"]
sig_g <- function(d, grp, thr_fdr = TRUE) unique(d[group == grp & (if (thr_fdr) fdr < 0.05 else pval < 0.05), gene])
a_f <- sig_g(ukb, "AMD"); s_f <- sig_g(ukb, "STROKE_WMH"); a_p <- sig_g(ukb, "AMD", FALSE); s_p <- sig_g(ukb, "STROKE_WMH", FALSE)
strict <- intersect(a_f, s_f)
loose <- setdiff(union(intersect(a_f, s_p), intersect(s_f, a_p)), strict)
message(sprintf("UKB-PPP FDR<0.05: AMD %d proteins, stroke/WMH %d proteins; shared strict %d, loose %d",
                length(a_f), length(s_f), length(strict), length(loose)))

dec <- res[source == "deCODE", .(gene, outcome, b_dec = b, se_dec = se, p_dec = pval, nsnp_dec = nsnp)]
dec <- dec[order(p_dec)][!duplicated(paste(gene, outcome))]
rep <- merge(ukb[fdr < 0.05, .(gene, protein_id, outcome, group, b, se, pval, fdr, OR, steiger_dir, region)], dec,
             by = c("gene", "outcome"), all.x = TRUE)
rep[, replicated := !is.na(p_dec) & sign(b_dec) == sign(b) & p_dec < 0.05]
rep[, OR_dec := exp(b_dec)]
fwrite(rep[order(pval)], file.path(OUT, "pqtl_replication_decode.tsv"), sep = "\t")
message(sprintf("replication: %d / %d UKB-PPP FDR-significant pairs tested in deCODE, %d replicated",
                sum(!is.na(rep$p_dec)), nrow(rep), sum(rep$replicated)))

# eQTL (eQTLGen single-SNP SMR) direction for the same gene x outcome
smr <- fread(file.path(PROJ_ROOT, "results/smr/smr_all.tsv"))[, .(gene = GeneSymbol, outcome, b_smr, p_smr, fdr_smr = p_fdr)]
smr <- smr[order(p_smr)][!duplicated(paste(gene, outcome))]

sh_genes <- c(strict, loose)
shared <- merge(ukb[gene %in% sh_genes], dec, by = c("gene", "outcome"), all.x = TRUE)
shared <- merge(shared, smr, by = c("gene", "outcome"), all.x = TRUE)
shared[, `:=`(tier = fifelse(gene %in% strict, "FDR in both", "FDR in one + P<0.05 in other"),
              replicated_decode = !is.na(p_dec) & sign(b_dec) == sign(b) & p_dec < 0.05,
              eqtl_same_dir = fifelse(is.na(b_smr), NA, sign(b_smr) == sign(b)), OR_dec = exp(b_dec))]
fwrite(shared[order(tier, gene, outcome)], file.path(OUT, "pqtl_mr_shared.tsv"), sep = "\t")

# ---------------- 6. Venn: AMD eQTL genes (SMR FDR<0.05) x stroke/WMH pQTL proteins (UKB-PPP FDR<0.05), and the reverse
smr_all <- fread(file.path(PROJ_ROOT, "results/smr/smr_all.tsv"))
amd_e <- unique(smr_all[outcome %in% AMD & p_fdr < 0.05, GeneSymbol]); stk_e <- unique(smr_all[outcome %in% STK & p_fdr < 0.05, GeneSymbol])
x1 <- intersect(amd_e, s_f); x2 <- intersect(stk_e, a_f)
writeLines(c(paste0("# AMD eQTL genes (n=", length(amd_e), ") x stroke/WMH pQTL proteins (n=", length(s_f), "): ", length(x1)), x1, "",
             paste0("# stroke/WMH eQTL genes (n=", length(stk_e), ") x AMD pQTL proteins (n=", length(a_f), "): ", length(x2)), x2, "",
             paste0("# AMD pQTL x stroke/WMH pQTL (strict shared): ", length(strict)), strict),
           file.path(OUT, "xqtl_intersections.txt"))
suppressPackageStartupMessages(library(VennDiagram)); futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
pdf(file.path(OUT, "Fig_venn_xQTL.pdf"), width = 10, height = 4)
grid::grid.newpage()
vp <- function(x) grid::viewport(x = x, width = 0.33, height = 0.9)
for (i in 1:3) {
  s <- list(list(amd_e, s_f, c("AMD eQTL genes", "Stroke/WMH pQTL proteins")),
            list(stk_e, a_f, c("Stroke/WMH eQTL genes", "AMD pQTL proteins")),
            list(a_f, s_f, c("AMD pQTL proteins", "Stroke/WMH pQTL proteins")))[[i]]
  grid::pushViewport(vp(c(1, 3, 5)[i] / 6))
  tryCatch(grid::grid.draw(draw.pairwise.venn(length(s[[1]]), length(s[[2]]), length(intersect(s[[1]], s[[2]])), category = s[[3]],
                                     fill = c("#4C78A8", "#E45756"), alpha = 0.5, cex = 1.2, cat.cex = 0.8, cat.pos = c(-20, 20), ind = FALSE)),
           error = function(e) grid::grid.text(paste(s[[3]], collapse = " x ")))
  grid::popViewport()
}
dev.off()

# ---------------- 7. forest plot of shared proteins across all outcomes (both platforms)
if (length(sh_genes)) {
  fp <- res[gene %in% sh_genes]
  fp[, outcome := factor(outcome, levels = rev(TRAITS))]
  fp[, sig := fifelse(fdr < 0.05, "FDR < 0.05", fifelse(pval < 0.05, "P < 0.05", "NS"))]
  g <- ggplot(fp, aes(x = OR, y = outcome, colour = source, shape = sig)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_pointrange(aes(xmin = lci, xmax = uci), position = position_dodge(width = 0.6), size = 0.3) +
    scale_x_log10() + scale_shape_manual(values = c("FDR < 0.05" = 16, "P < 0.05" = 17, "NS" = 1)) +
    scale_colour_manual(values = c("UKB-PPP" = "#2b5c8a", "deCODE" = "#c0392b")) +
    facet_wrap(~ gene, scales = "free_x") + labs(x = "OR per 1-SD higher genetically predicted protein (95% CI)", y = NULL) +
    theme_bw(base_size = 9)
  nf <- length(unique(fp$gene)); nc <- min(4, nf)
  ggsave(file.path(OUT, "Fig_forest_shared_proteins.pdf"), g, width = 2.8 * nc, height = 2.6 * ceiling(nf / nc) + 0.6, limitsize = FALSE)
  ggsave(file.path(OUT, "Fig_forest_shared_proteins.png"), g, width = 2.8 * nc, height = 2.6 * ceiling(nf / nc) + 0.6, dpi = 200, limitsize = FALSE)
}
message("pQTL MR done")
