# Transcriptome-wide SMR + HEIDI (SMR 1.3.1, eQTLGen whole-blood cis-eQTL BESD) for AMD and stroke/WMH outcomes,
# restricted to the druggable genome (Finan 2017 via DGIdb), followed by eQTL-GWAS colocalisation (coloc.abf).
# Following the "druggable-gene MR + colocalisation" design (Adv Sci 2026, 2.4).
# SMR defaults: top cis-eQTL p < 5e-8, HEIDI on up to 20 SNPs with eQTL p < 1.57e-3 (>= 3 SNPs required), LD from 1000G EUR.
# Significance: SMR FDR < 0.05 within outcome; pleiotropy/linkage filter P_HEIDI >= 0.05.
# coloc: eQTL (quant, N = 31,684) vs GWAS (cc with case fraction s, or quant for WMH), default priors; PP.H4 > 0.8 strong, 0.5-0.8 moderate.
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(coloc); library(ggplot2); library(parallel) })
setDTthreads(30)
SMR <- file.path(PROJ_ROOT, "data/raw/tools/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe")
PLINK <- file.path(PROJ_ROOT, "data/raw/tools/plink/plink.exe")
BFILE <- file.path(PROJ_ROOT, "data/raw/ref/g1000_eur/g1000_eur"); BESD <- file.path(PROJ_ROOT, "data/raw/eqtl/eqtlgen_besd")
MA <- file.path(PROJ_ROOT, "derived/smr_ma"); OUT <- file.path(PROJ_ROOT, "results/smr_heidi"); RAW <- file.path(OUT, "raw")
dir.create(RAW, recursive = TRUE, showWarnings = FALSE)
AMDS <- c("AMD_FG", "AMD_IAMDGC", "DRYAMD_FG", "WETAMD_FG"); STRK <- c("AS", "IS", "LAS", "CES", "SVS", "WMH")
OUTCOMES <- c(AMDS, STRK)
SHARED13 <- c("SH2B3", "ACAD10", "FAM109A", "MAPKAPK5-AS1", "C4A", "C4B", "CYP21A2", "CYP21A1P", "TNXA", "AGER", "PBX2",
              "BTN3A2", "RP3-473L9.4")
N_EQTL <- 31684
meta <- readRDS(file.path(PROJ_ROOT, "derived/gwas/trait_meta.rds"))

# ---- 1. SMR + HEIDI per outcome (SMR's text parser is slow on Windows: keep only SNPs present in both the eQTL
#         and LD-reference files, then run the 10 outcomes as concurrent processes with 3 threads each)
MAF_DIR <- file.path(PROJ_ROOT, "derived/smr_ma_filtered"); dir.create(MAF_DIR, showWarnings = FALSE)
LD <- file.path(PROJ_ROOT, "derived/ld_eqtlgen/g1000_eur_eqtlgen")            # 1000G EUR restricted to eQTLGen SNPs
if (!file.exists(paste0(LD, ".bim"))) {
  writeLines(fread(paste0(BESD, ".esi"), header = FALSE, select = 2)[[1]], paste0(LD, "_snps.txt"))
  system2(PLINK, c("--bfile", BFILE, "--extract", paste0(LD, "_snps.txt"), "--make-bed", "--threads", "30", "--out", LD), stdout = FALSE)
}
todo <- OUTCOMES[!file.exists(file.path(RAW, paste0(OUTCOMES, ".smr")))]
if (length(todo)) {
  keep <- fread(paste0(LD, ".bim"), header = FALSE, select = 2)[[1]]
  message("SNPs in eQTL & LD reference: ", length(keep))
  for (k in todo) {
    f <- file.path(MAF_DIR, paste0(k, ".ma"))
    if (!file.exists(f)) fwrite(fread(file.path(MA, paste0(k, ".ma")))[SNP %chin% keep], f, sep = " ")
  }
  ps <- lapply(todo, function(k) {
    o <- file.path(RAW, k)
    system2(SMR, c("--bfile", LD, "--gwas-summary", file.path(MAF_DIR, paste0(k, ".ma")), "--beqtl-summary", BESD,
                   "--diff-freq", "0.2", "--diff-freq-prop", "0.2", "--thread-num", "3", "--out", o),
            stdout = paste0(o, ".stdout.log"), stderr = paste0(o, ".stdout.log"), wait = FALSE)
  })
  repeat {                                   # wait until every SMR process has exited (.smr is written incrementally)
    Sys.sleep(60)
    alive <- length(grep("smr-1.3.1", system("tasklist /FI \"IMAGENAME eq smr-1.3.1-win.exe\" /NH", intern = TRUE)))
    if (alive == 0) break
  }
  miss <- todo[!file.exists(file.path(RAW, paste0(todo, ".smr")))]
  if (length(miss)) stop("SMR failed for ", paste(miss, collapse = ", "), " (see raw/*.stdout.log)")
}
smr <- rbindlist(lapply(OUTCOMES, function(k) fread(file.path(RAW, paste0(k, ".smr")))[, outcome := k]), fill = TRUE)
smr[, p_fdr := p.adjust(p_SMR, "BH"), by = outcome]
smr[, p_bonf := pmin(1, p_SMR * .N), by = outcome]
smr[, heidi_pass := !is.na(p_HEIDI) & p_HEIDI >= 0.05]
smr[, sig := p_fdr < 0.05 & heidi_pass]
smr[, OR_SMR := exp(b_SMR)]
# eQTLGen BESD stores Ensembl IDs in the gene column: map to symbols from the eQTLGen cis-eQTL table
gmap <- unique(fread(file.path(PROJ_ROOT, "data/raw/eqtl/eqtlgen_cis_sig.txt.gz"), select = c("Gene", "GeneSymbol")), by = "Gene")
smr[, Gene := NULL]
smr <- merge(smr, gmap, by.x = "probeID", by.y = "Gene", all.x = TRUE)
smr[is.na(GeneSymbol), GeneSymbol := probeID]
fwrite(smr, file.path(OUT, "smr_heidi_all.tsv"), sep = "\t")

# ---- 2. druggable genome annotation
cat_f <- file.path(OUT, "dgidb_gene_categories.tsv")
dg <- if (file.exists(cat_f)) fread(cat_f) else data.table(GeneSymbol = character(), druggable_finan = logical())
smr <- merge(smr, dg[, .(GeneSymbol, druggable_finan, druggable_any, clinically_actionable, n_approved_drugs, approved_drugs)],
             by = "GeneSymbol", all.x = TRUE)
for (v in c("druggable_finan", "druggable_any", "clinically_actionable")) smr[is.na(get(v)), (v) := FALSE]
smr[is.na(n_approved_drugs), n_approved_drugs := 0L]
smr[, druggable := druggable_finan | n_approved_drugs > 0]

per_outcome <- smr[, .(n_tested = .N, n_fdr = sum(p_fdr < 0.05), n_bonf = sum(p_bonf < 0.05),
                       n_heidi_tested = sum(p_fdr < 0.05 & !is.na(p_HEIDI)), n_fdr_heidi_pass = sum(sig),
                       n_druggable_tested = sum(druggable), n_druggable_sig = sum(sig & druggable)), by = outcome]
fwrite(per_outcome, file.path(OUT, "smr_heidi_counts.tsv"), sep = "\t"); print(per_outcome)

# shared genes: FDR + HEIDI in an AMD outcome and FDR + HEIDI in a stroke/WMH outcome (strict),
# or FDR + HEIDI in AMD and nominal SMR p < 0.05 + HEIDI in stroke/WMH (loose)
g_a <- smr[outcome %in% AMDS & sig, unique(probeID)]
g_s <- smr[outcome %in% STRK & sig, unique(probeID)]
g_s_nom <- smr[outcome %in% STRK & p_SMR < 0.05 & heidi_pass, unique(probeID)]
shared_strict <- intersect(g_a, g_s); shared_loose <- intersect(g_a, g_s_nom)
sym <- unique(smr[, .(probeID, GeneSymbol, druggable, druggable_finan, n_approved_drugs, approved_drugs)])
message(sprintf("AMD sig (FDR+HEIDI) %d; stroke/WMH sig %d; shared strict %d; shared loose %d",
                length(g_a), length(g_s), length(shared_strict), length(shared_loose)))

# ---- 3. eQTL-GWAS colocalisation for candidate probes x all outcomes
cand <- unique(c(shared_loose, smr[GeneSymbol %in% SHARED13, unique(probeID)],
                 smr[sig & druggable, unique(probeID)]))
cand <- cand[cand %in% smr[sig | GeneSymbol %in% SHARED13 | probeID %in% shared_loose, probeID]]
message("coloc candidate probes: ", length(cand))
if (!length(cand)) stop("no candidate probes for colocalisation")
qf <- file.path(RAW, "coloc_probes.txt"); writeLines(cand, qf)
qo <- file.path(RAW, "coloc_query")
if (!file.exists(paste0(qo, ".txt")))
  system2(SMR, c("--beqtl-summary", BESD, "--extract-probe", qf, "--query", "1", "--out", qo),
          stdout = paste0(qo, ".log"), stderr = paste0(qo, ".log"))
eq <- fread(paste0(qo, ".txt"))   # SNP Chr BP A1 A2 Freq Probe Probe_Chr Probe_bp Gene Orientation b SE p
eq <- eq[is.finite(b) & is.finite(SE) & SE > 0 & Freq > 0 & Freq < 1]
message("eQTL query rows: ", nrow(eq))
gw_snps <- unique(eq$SNP)
gw <- rbindlist(lapply(OUTCOMES, function(k)
  fread(file.path(MA, paste0(k, ".ma")))[SNP %chin% gw_snps][, outcome := k]))

coloc_one <- function(job) {
  e <- eq[Probe == job$probe]; g <- gw[outcome == job$outcome]
  m <- merge(e, g, by = "SNP", suffixes = c(".e", ".g"))
  m[, by := fifelse(A1.g == A1.e & A2.g == A2.e, b.g, fifelse(A1.g == A2.e & A2.g == A1.e, -b.g, NA_real_))]
  m <- m[!is.na(by) & is.finite(se) & se > 0]
  m <- unique(m, by = "SNP")
  if (nrow(m) < 50) return(NULL)
  maf <- pmin(m$Freq, 1 - m$Freq)
  d1 <- list(beta = m$b.e, varbeta = m$SE^2, snp = m$SNP, type = "quant", N = N_EQTL, MAF = maf)
  tm <- meta[[job$outcome]]
  d2 <- if (!is.null(tm$ncase)) list(beta = m$by, varbeta = m$se^2, snp = m$SNP, type = "cc", s = tm$ncase / (tm$ncase + tm$ncontrol))
        else list(beta = m$by, varbeta = m$se^2, snp = m$SNP, type = "quant", N = median(m$n), MAF = maf)
  cr <- tryCatch(suppressMessages(coloc.abf(d1, d2)), error = function(e) NULL)
  if (is.null(cr)) return(NULL)
  top <- cr$results[order(-cr$results$SNP.PP.H4), ][1, ]
  data.table(probeID = job$probe, outcome = job$outcome, nsnps = nrow(m),
             PP.H0 = cr$summary[["PP.H0.abf"]], PP.H1 = cr$summary[["PP.H1.abf"]], PP.H2 = cr$summary[["PP.H2.abf"]],
             PP.H3 = cr$summary[["PP.H3.abf"]], PP.H4 = cr$summary[["PP.H4.abf"]], coloc_topSNP = top$snp, SNP.PP.H4 = top$SNP.PP.H4)
}
jobs <- CJ(probe = unique(eq$Probe), outcome = OUTCOMES)
cl <- makeCluster(30)
clusterEvalQ(cl, { .libPaths(c(file.path(Sys.getenv("COMORBIDITY_ROOT"), "Rlib"), .libPaths())); suppressPackageStartupMessages({ library(data.table); library(coloc) }); setDTthreads(1) })
clusterExport(cl, c("eq", "gw", "meta", "N_EQTL"))
cres <- rbindlist(parLapply(cl, split(jobs, seq_len(nrow(jobs))), coloc_one))
stopCluster(cl)
fwrite(cres, file.path(OUT, "coloc_eqtl_gwas.tsv"), sep = "\t")
message("coloc done: ", nrow(cres), " gene-outcome pairs; PP.H4 > 0.8: ", sum(cres$PP.H4 > 0.8))

# ---- 4. summary tables and figure
res <- merge(smr, cres, by = c("probeID", "outcome"), all.x = TRUE)
res[, coloc_class := fcase(is.na(PP.H4), NA_character_, PP.H4 > 0.8, "strong", PP.H4 > 0.5, "moderate", default = "no")]
fwrite(res[order(p_SMR)][sig | GeneSymbol %in% SHARED13 | !is.na(PP.H4)], file.path(OUT, "smr_heidi_coloc_candidates.tsv"), sep = "\t")

long <- res[probeID %in% unique(c(shared_loose, smr[GeneSymbol %in% SHARED13, probeID])),
            .(GeneSymbol, probeID, outcome, b_SMR, se_SMR, OR_SMR, p_SMR, p_fdr, p_HEIDI, nsnp_HEIDI, PP.H4, coloc_class,
              druggable, druggable_finan, n_approved_drugs, approved_drugs)]
long[, `:=`(group = fifelse(outcome %in% AMDS, "AMD", "Stroke/WMH"),
            tier = fifelse(probeID %in% shared_strict, "FDR+HEIDI in both", fifelse(probeID %in% shared_loose, "AMD FDR+HEIDI, stroke nominal", "13 shared genes (previous SMR)")))]
fwrite(long[order(tier, GeneSymbol, outcome)], file.path(OUT, "druggable_shared_summary.tsv"), sep = "\t")

wide <- dcast(long, GeneSymbol + probeID + tier + druggable + druggable_finan + n_approved_drugs + approved_drugs ~ outcome,
              value.var = c("OR_SMR", "p_SMR", "p_HEIDI", "PP.H4"))
fwrite(wide, file.path(OUT, "druggable_shared_wide.tsv"), sep = "\t")

fp <- long[druggable == TRUE & tier != "13 shared genes (previous SMR)"]
if (!nrow(fp)) fp <- long[tier != "13 shared genes (previous SMR)"]
if (nrow(fp)) {
  fp[, `:=`(lo = exp(b_SMR - 1.96 * se_SMR), hi = exp(b_SMR + 1.96 * se_SMR),
            mark = fcase(p_fdr < 0.05 & !is.na(p_HEIDI) & p_HEIDI >= 0.05, "FDR<0.05 & HEIDI pass", p_SMR < 0.05, "P<0.05", default = "n.s."))]
  fp[, outcome := factor(outcome, levels = rev(OUTCOMES))]
  p <- ggplot(fp, aes(OR_SMR, outcome, colour = group, shape = mark)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) + geom_point(size = 2) +
    scale_x_log10() + scale_shape_manual(values = c("FDR<0.05 & HEIDI pass" = 16, "P<0.05" = 1, "n.s." = 4)) +
    facet_wrap(~ GeneSymbol, scales = "free_x") +
    labs(x = "OR per 1-SD increase in genetically predicted expression (SMR)", y = NULL, colour = NULL, shape = NULL) +
    theme_bw(base_size = 9) + theme(legend.position = "bottom")
  ggsave(file.path(OUT, "Fig_SMR_druggable_shared_forest.pdf"), p, width = 10, height = 2.2 + 1.8 * ceiling(uniqueN(fp$GeneSymbol) / 3))
  ggsave(file.path(OUT, "Fig_SMR_druggable_shared_forest.png"), p, width = 10, height = 2.2 + 1.8 * ceiling(uniqueN(fp$GeneSymbol) / 3), dpi = 200)
}
message("SMR-HEIDI-druggable done")
