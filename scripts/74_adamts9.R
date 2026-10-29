# 3p14 ADAMTS9 deep-dive (report section 4.2 follow-up): is the opposite-direction AMD / ischaemic-stroke signal mediated
# by the same arterial ADAMTS9 cis-eQTL?
#   1) multi-tissue SMR + HEIDI with the SMR 1.3.1 binary (region BESD built from eQTL Catalogue / EyeGEx summary data,
#      LD = 1000G EUR); top cis-eQTL P < 5e-8, and P < 5e-6 as a relaxed second pass
#   2) coloc.abf for every tissue x gene x trait (default priors; p12 = 5e-6 sensitivity), GWAS-GWAS coloc in the window
#   3) three-way consistency AMD ~ IS ~ aortic ADAMTS9 eQTL (credible sets, joint three-trait ABF), LD among the leads
#   4) key-SNP effects across tissues / traits, forest plot and LocusZoom-style region plot
# Inputs: derived/adamts9/ (74a_adamts9_data.R). Output -> results/adamts9/
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressMessages({ library(data.table); library(parallel); library(coloc) })
D <- file.path(PROJ_ROOT, "derived/adamts9"); OUT <- file.path(PROJ_ROOT, "results/adamts9"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
TMP <- file.path(D, "smr"); dir.create(TMP, showWarnings = FALSE)
SMR <- file.path(PROJ_ROOT, "data/raw/tools/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe")
PLINK <- file.path(PROJ_ROOT, "data/raw/tools/plink/plink.exe"); BREG <- file.path(D, "g1000_eur_3p14")
OFF <- 14324
TRAITS <- c("AMD_FG", "DRYAMD_FG", "WETAMD_FG", "AMD_IAMDGC", "AS", "IS", "LAS", "CES", "SVS", "WMH",
            "LACUNAR", "EPVS_WM", "EPVS_BG", "EPVS_HIP", "WM_MD", "WM_FA", "CIMT", "PLAQUE")
QUANT <- c("WMH", "WM_MD", "WM_FA", "CIMT")
KEYSNP <- c("rs4611812", "rs6795735", "rs17676309", "rs6445425")   # GWAS-GWAS coloc lead, AMD lead, Franceschini 2018 cIMT/eQTL SNP
AMD <- TRAITS[1:4]
FOCUS <- c(ADAMTS9 = "ENSG00000163638", "ADAMTS9-AS2" = "ENSG00000241684", MAGI1 = "ENSG00000151276")
setDTthreads(8)

# ---------------- eQTL table
eq <- rbindlist(lapply(list.files(file.path(D, "eqtl"), pattern = "^(GTEx|ROSMAP|BrainSeq).*\\.tsv\\.gz$", full.names = TRUE), function(f) {
  d <- fread(f); if (!nrow(d)) return(NULL)
  d[, .(dataset, gene_id, SNP = rsid, ea = alt, oa = ref, beta, se, p = pvalue, maf, af_ea = ac / an, n, pos38 = position)]
}))
frq <- fread(file.path(PROJ_ROOT, "derived/g1000_eur.frq"))[, .(SNP, fA1 = A1, fA2 = A2, MAF)]
frq <- frq[SNP %in% fread(file.path(D, "region_snps.tsv"))$SNP]      # region only
CACHE <- Sys.getenv("ADAMTS9_RERUN") != "1"                               # reuse SMR / coloc tables if already written
ey <- fread(file.path(D, "eqtl", "EyeGEx_Retina.tsv.gz"))
ey <- merge(ey, frq, by = "SNP", all.x = TRUE)
ey[, af_ea := fifelse(toupper(ea) == fA1, MAF, fifelse(toupper(ea) == fA2, 1 - MAF, NA_real_))]
eq <- rbind(eq, ey[is.finite(se) & se > 0, .(dataset = "EyeGEx_Retina", gene_id, SNP, ea, oa, beta, se, p, maf = pmin(af_ea, 1 - af_ea),
                                               af_ea, n = 406, pos38 = vpos)])
eq <- unique(eq[!is.na(SNP) & grepl("^rs", SNP) & is.finite(se) & se > 0], by = c("dataset", "gene_id", "SNP"))

# gene symbols / TSS (Ensembl REST, cached)
gsf <- file.path(D, "gene_symbols.tsv")
if (!file.exists(gsf)) {
  ids <- unique(eq$gene_id); res <- list()
  for (ch in split(ids, ceiling(seq_along(ids) / 200))) {
    r <- httr::POST("https://rest.ensembl.org/lookup/id", httr::content_type_json(), httr::accept_json(),
                    body = jsonlite::toJSON(list(ids = ch)))
    j <- jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"), simplifyVector = FALSE)
    res[[length(res) + 1]] <- rbindlist(lapply(j[!vapply(j, is.null, TRUE)], function(x)
      data.table(gene_id = x$id, gene = if (is.null(x$display_name)) x$id else x$display_name, biotype = x$biotype,
                 start38 = x$start, end38 = x$end, strand = x$strand)))
    Sys.sleep(1)
  }
  fwrite(rbindlist(res), gsf, sep = "\t")
}
gs <- fread(gsf); gs[, tss38 := fifelse(strand == 1, start38, end38)]
eq <- merge(eq, gs[, .(gene_id, gene)], by = "gene_id", all.x = TRUE); eq[is.na(gene), gene := gene_id]
cov <- eq[, .(n_snps = .N, min_p = min(p), top_snp = SNP[which.min(p)], n = n[1]), by = .(gene, gene_id, dataset)]
fwrite(cov[order(gene, min_p)], file.path(OUT, "eqtl_coverage_3p14.tsv"), sep = "\t")
sel_genes <- unique(c(FOCUS, cov[min_p < 1e-5, gene_id]))
eqs <- eq[gene_id %in% sel_genes]
message("eQTL rows ", nrow(eq), "; datasets ", uniqueN(eq$dataset), "; genes analysed: ",
        paste(unique(gs[gene_id %in% sel_genes, gene]), collapse = ", "))

# ---------------- GWAS region + .ma files (region only)
gw <- readRDS(file.path(D, "gwas_region.rds"))
if (file.exists(file.path(D, "gwas_cimt.rds"))) gw <- rbind(gw, readRDS(file.path(D, "gwas_cimt.rds")), fill = TRUE)   # 74b
gw[, `:=`(ea = toupper(ea), oa = toupper(oa))]
for (k in TRAITS) {
  d <- merge(gw[trait == k, .(SNP, A1 = ea, A2 = oa, freq = eaf, b = beta, se, p, n)], frq, by = "SNP", all.x = TRUE)
  d[, f1kg := fifelse(A1 == fA1, MAF, fifelse(A1 == fA2, 1 - MAF, NA_real_))]
  d[is.na(freq) | freq <= 0 | freq >= 1, freq := f1kg]
  d <- unique(d[!is.na(freq) & freq > 0 & freq < 1 & is.finite(se) & se > 0], by = "SNP")
  if (all(is.na(d$n))) d[, n := 10000]
  fwrite(d[, .(SNP, A1, A2, freq, b, se, p, n)], file.path(TMP, paste0(k, ".ma")), sep = " ")
}

# ---------------- region BESD per eQTL dataset (query format -> --make-besd); positions GRCh37
rs <- fread(file.path(D, "region_snps.tsv"))
eqs <- merge(eqs, rs[, .(SNP, pos37)], by = "SNP", all.x = TRUE)
eqs[is.na(pos37), pos37 := pos38 - OFF]
eqs <- merge(eqs, gs[, .(gene_id, tss38, strand)], by = "gene_id", all.x = TRUE)
for (ds in unique(eqs$dataset)) {
  q <- eqs[dataset == ds & !is.na(af_ea) & af_ea > 0 & af_ea < 1 & !is.na(tss38),
           .(SNP, Chr = 3, BP = pos37, A1 = toupper(ea), A2 = toupper(oa), Freq = af_ea, Probe = gene_id, Probe_Chr = 3,
             Probe_bp = tss38 - OFF, Gene = gene, Orientation = fifelse(strand == 1, "+", "-"), b = beta, se, p)]
  if (!nrow(q)) next
  qf <- file.path(TMP, paste0(ds, ".query")); fwrite(q, qf, sep = "\t")
  system2(SMR, c("--qfile", qf, "--make-besd", "--out", file.path(TMP, ds)), stdout = FALSE, stderr = FALSE)
}
besds <- sub("\\.besd$", "", list.files(TMP, pattern = "\\.besd$"))

# ---------------- SMR + HEIDI (two eQTL thresholds)
smr_f <- file.path(OUT, "smr_heidi_multitissue.tsv")
if (CACHE && file.exists(smr_f)) { smr_res <- fread(smr_f); message("SMR: cached") } else {
jobs <- CJ(ds = besds, trait = TRAITS, thr = c("5e-8", "5e-6"))
cl <- makeCluster(30)
clusterExport(cl, c("jobs", "TMP", "SMR", "BREG"))
clusterEvalQ(cl, { .libPaths(c(file.path(Sys.getenv("COMORBIDITY_ROOT"), "Rlib"), .libPaths())); suppressMessages(library(data.table)) })
smr_res <- rbindlist(parLapplyLB(cl, seq_len(nrow(jobs)), function(i) {
  j <- jobs[i]; o <- file.path(TMP, sprintf("res_%s_%s_%s", j$ds, j$trait, j$thr))
  system2(SMR, c("--bfile", BREG, "--gwas-summary", file.path(TMP, paste0(j$trait, ".ma")), "--beqtl-summary", file.path(TMP, j$ds),
                 "--peqtl-smr", j$thr, "--diff-freq", "0.2", "--diff-freq-prop", "0.5", "--thread-num", "1", "--out", o),
          stdout = FALSE, stderr = FALSE)
  f <- paste0(o, ".smr"); if (!file.exists(f)) return(NULL)
  d <- data.table::fread(f); if (!nrow(d)) return(NULL)
  d[, `:=`(dataset = j$ds, trait = j$trait, peqtl_threshold = j$thr)]; d
}), fill = TRUE)
stopCluster(cl)
if (nrow(smr_res)) {
  smr_res <- smr_res[order(peqtl_threshold == "5e-6")]                 # keep the strict-threshold row when both exist
  smr_res <- unique(smr_res, by = c("dataset", "trait", "probeID"))
  smr_res[, OR_SMR := exp(b_SMR)]
  smr_res[, p_fdr_trait := p.adjust(p_SMR, "BH"), by = trait]
  smr_res[, heidi_pass := !is.na(p_HEIDI) & p_HEIDI >= 0.05]
}
fwrite(smr_res, smr_f, sep = "\t")
message("SMR rows: ", nrow(smr_res))
}

# ---------------- coloc.abf: tissue eQTL x trait; GWAS x GWAS
gwl <- split(gw, by = "trait")
meta <- gw[, .(ncase = ncase[1], ncontrol = ncontrol[1], nmax = max(n, na.rm = TRUE)), by = trait]
mk_gwas <- function(m, tr, b, se, maf) {
  mt <- meta[trait == tr]
  if (tr %in% QUANT) list(beta = b, varbeta = se^2, snp = m$SNP, type = "quant", N = mt$nmax, MAF = maf)
  else list(beta = b, varbeta = se^2, snp = m$SNP, type = "cc", s = mt$ncase / (mt$ncase + mt$ncontrol), MAF = maf)
}
co_f <- file.path(OUT, "coloc_multitissue.tsv")
if (CACHE && file.exists(co_f)) { co <- fread(co_f); message("coloc: cached") } else {
ejobs <- CJ(key = unique(eqs[, paste(gene, dataset, sep = "|")]), trait = TRAITS)
cl <- makeCluster(30)
clusterExport(cl, c("eqs", "gwl", "ejobs", "meta", "QUANT", "mk_gwas"))
clusterEvalQ(cl, { .libPaths(c(file.path(Sys.getenv("COMORBIDITY_ROOT"), "Rlib"), .libPaths())); suppressMessages({ library(data.table); library(coloc) }) })
co <- rbindlist(parLapplyLB(cl, seq_len(nrow(ejobs)), function(i) {
  k <- strsplit(ejobs$key[i], "\\|")[[1]]; tr <- ejobs$trait[i]
  e <- eqs[gene == k[1] & dataset == k[2]]; g <- gwl[[tr]]
  m <- merge(e, g[, .(SNP, gea = ea, goa = oa, gb = beta, gse = se, gp = p)], by = "SNP")
  m <- m[(toupper(ea) == gea & toupper(oa) == goa) | (toupper(ea) == goa & toupper(oa) == gea)]
  m[, gb_al := fifelse(toupper(ea) == gea, gb, -gb)]
  m <- m[is.finite(gb) & is.finite(gse) & gse > 0 & !is.na(maf) & maf > 0.01 & maf < 0.99]
  if (nrow(m) < 50) return(NULL)
  d1 <- list(beta = m$beta, varbeta = m$se^2, snp = m$SNP, type = "quant", N = m$n[1], MAF = m$maf)
  d2 <- mk_gwas(m, tr, m$gb_al, m$gse, m$maf)
  r <- tryCatch(suppressWarnings(suppressMessages(coloc.abf(d1, d2))), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  r2 <- tryCatch(suppressWarnings(suppressMessages(coloc.abf(d1, d2, p12 = 5e-6))), error = function(e) NULL)
  pp <- r$results[order(-r$results$SNP.PP.H4), c("snp", "SNP.PP.H4")]
  cs95 <- pp$snp[seq_len(which(cumsum(pp$SNP.PP.H4) >= 0.95)[1])]
  top <- pp$snp[1]; mt <- m[SNP == top]
  data.table(gene = k[1], dataset = k[2], trait = tr, nsnps = nrow(m), min_p_eqtl = min(m$p), min_p_gwas = min(m$gp),
             top_snp_eqtl = m$SNP[which.min(m$p)], top_snp_gwas = m$SNP[which.min(m$gp)],
             PP.H3 = r$summary[["PP.H3.abf"]], PP.H4 = r$summary[["PP.H4.abf"]],
             PP.H4_p12_5e6 = if (is.null(r2)) NA_real_ else r2$summary[["PP.H4.abf"]],
             coloc_snp = top, coloc_snp_PP = pp$SNP.PP.H4[1], cs95_size = length(cs95), cs95 = paste(cs95, collapse = ","),
             dir_expr_to_risk = sign(mt$beta * mt$gb_al))
}))
gg <- CJ(t1 = AMD, t2 = setdiff(TRAITS, AMD))
cog <- rbindlist(lapply(seq_len(nrow(gg)), function(i) {
  a <- gwl[[gg$t1[i]]]; b <- gwl[[gg$t2[i]]]
  m <- merge(a[, .(SNP, ea, oa, b1 = beta, s1 = se, p1 = p)], b[, .(SNP, ea2 = ea, oa2 = oa, b2 = beta, s2 = se, p2 = p)], by = "SNP")
  m <- m[(ea == ea2 & oa == oa2) | (ea == oa2 & oa == ea2)]; m[, b2 := fifelse(ea == ea2, b2, -b2)]
  m <- merge(m, frq[, .(SNP, MAF)], by = "SNP")[MAF > 0.01 & is.finite(s1) & is.finite(s2) & s1 > 0 & s2 > 0]
  if (nrow(m) < 50) return(NULL)
  r <- suppressWarnings(suppressMessages(coloc.abf(mk_gwas(m, gg$t1[i], m$b1, m$s1, m$MAF), mk_gwas(m, gg$t2[i], m$b2, m$s2, m$MAF))))
  pp <- r$results[order(-r$results$SNP.PP.H4), c("snp", "SNP.PP.H4")]
  cs95 <- pp$snp[seq_len(which(cumsum(pp$SNP.PP.H4) >= 0.95)[1])]
  top <- pp$snp[1]; mt <- m[SNP == top]
  data.table(gene = "GWAS", dataset = gg$t1[i], trait = gg$t2[i], nsnps = nrow(m), min_p_eqtl = min(m$p1), min_p_gwas = min(m$p2),
             top_snp_eqtl = m$SNP[which.min(m$p1)], top_snp_gwas = m$SNP[which.min(m$p2)],
             PP.H3 = r$summary[["PP.H3.abf"]], PP.H4 = r$summary[["PP.H4.abf"]], PP.H4_p12_5e6 = NA_real_,
             coloc_snp = top, coloc_snp_PP = pp$SNP.PP.H4[1], cs95_size = length(cs95), cs95 = paste(cs95, collapse = ","),
             dir_expr_to_risk = sign(mt$b1 * mt$b2))
}))
stopCluster(cl)
co <- rbind(co, cog)
co[, coloc_class := fifelse(PP.H4 > 0.8, "strong", fifelse(PP.H4 > 0.5, "moderate", fifelse(PP.H3 > 0.8, "distinct variants", "no/uncertain")))]
fwrite(co[order(gene, trait, -PP.H4)], co_f, sep = "\t")
message("coloc rows: ", nrow(co), "; strong: ", sum(co$PP.H4 > 0.8))
}

# ---------------- single-trait fine-mapping (95% credible sets) and three-way consistency
cs_single <- function(d) {   # 95% credible set conditional on one causal variant (null posterior reported separately)
  f <- finemap.abf(d); pnull <- f$SNP.PP[f$snp == "null"]; f <- f[f$snp != "null", ]
  f$pp <- f$SNP.PP / sum(f$SNP.PP); f <- f[order(-f$pp), ]
  structure(f$snp[seq_len(which(cumsum(f$pp) >= 0.95)[1])], pp_null = pnull)
}
aorta <- eqs[dataset == "GTEx_artery_aorta" & gene == "ADAMTS9"]
three <- list()
if (nrow(aorta)) {
  m <- Reduce(function(x, y) merge(x, y, by = "SNP"), list(
    aorta[, .(SNP, ea, oa, be = beta, se_e = se, maf, ne = n)],
    gwl[["AMD_FG"]][, .(SNP, a_ea = ea, a_oa = oa, ba = beta, sa = se)],
    gwl[["IS"]][, .(SNP, i_ea = ea, i_oa = oa, bi = beta, si = se)]))
  m <- m[(toupper(ea) == a_ea | toupper(ea) == a_oa) & (toupper(ea) == i_ea | toupper(ea) == i_oa) & maf > 0.01 & maf < 0.99]
  m[, ba := fifelse(toupper(ea) == a_ea, ba, -ba)]; m[, bi := fifelse(toupper(ea) == i_ea, bi, -bi)]
  d_e <- list(beta = m$be, varbeta = m$se_e^2, snp = m$SNP, type = "quant", N = m$ne[1], MAF = m$maf)
  d_a <- mk_gwas(m, "AMD_FG", m$ba, m$sa, m$maf); d_i <- mk_gwas(m, "IS", m$bi, m$si, m$maf)
  cs <- list(aorta_ADAMTS9 = cs_single(d_e), AMD_FG = cs_single(d_a), IS = cs_single(d_i))
  three$credible_sets <- rbindlist(lapply(names(cs), function(k) data.table(set = k, size = length(cs[[k]]), pp_null = attr(cs[[k]], "pp_null"),
                                                                           snps = paste(cs[[k]], collapse = ","))))
  three$overlap <- data.table(pair = c("aorta&AMD", "aorta&IS", "AMD&IS", "all3"),
                              n_shared = c(length(intersect(cs[[1]], cs[[2]])), length(intersect(cs[[1]], cs[[3]])),
                                           length(intersect(cs[[2]], cs[[3]])), length(Reduce(intersect, cs))),
                              shared = c(paste(intersect(cs[[1]], cs[[2]]), collapse = ","), paste(intersect(cs[[1]], cs[[3]]), collapse = ","),
                                         paste(intersect(cs[[2]], cs[[3]]), collapse = ","), paste(Reduce(intersect, cs), collapse = ",")))
  # joint posterior of each SNP being the single causal variant shared by all three traits (normalised product of the
  # per-trait Wakefield ABFs; HyPrColoc could not be compiled against the installed RcppEigen)
  lab <- function(d) { f <- finemap.abf(d); f <- f[f$snp != "null", ]; setNames(f$lABF., f$snp) }
  L <- lab(d_e) + lab(d_a)[names(lab(d_e))] + lab(d_i)[names(lab(d_e))]
  pj <- exp(L - max(L)); pj <- sort(pj / sum(pj), decreasing = TRUE)
  csj <- names(pj)[seq_len(which(cumsum(pj) >= 0.95)[1])]
  three$joint_abf <- data.table(SNP = names(pj)[1:min(20, length(pj))], joint_PP = pj[1:min(20, length(pj))])
  three$credible_sets <- rbind(three$credible_sets, data.table(set = "joint_3trait", size = length(csj), pp_null = NA_real_, snps = paste(csj, collapse = ",")))
  for (k in names(three)) fwrite(three[[k]], file.path(OUT, sprintf("three_way_%s.tsv", k)), sep = "\t")
}

# ---------------- LD among lead SNPs (GWAS leads per trait, eQTL leads per dataset for ADAMTS9 / ADAMTS9-AS2)
leads <- unique(c(KEYSNP, "rs7428936", "rs62247658", "rs13096863",
                  gw[, .SD[which.min(p)], by = trait]$SNP,
                  eqs[gene %in% c("ADAMTS9", "ADAMTS9-AS2"), .SD[which.min(p)], by = .(gene, dataset)]$SNP))
lf <- file.path(TMP, "leads.txt"); writeLines(leads, lf)
system2(PLINK, c("--bfile", BREG, "--extract", lf, "--r2", "square", "--write-snplist", "--out", file.path(TMP, "leads_ld")), stdout = FALSE)
if (file.exists(file.path(TMP, "leads_ld.ld"))) {
  ld <- as.matrix(fread(file.path(TMP, "leads_ld.ld"), header = FALSE)); sn <- readLines(file.path(TMP, "leads_ld.snplist"))
  dimnames(ld) <- list(sn, sn)
  ldt <- data.table(SNP = sn, ld)
  fwrite(ldt, file.path(OUT, "ld_r2_leads.tsv"), sep = "\t")
}
lead_tab <- rbind(gw[, .SD[which.min(p)], by = trait][, .(source = trait, type = "GWAS", lead = SNP, p)],
                  eqs[gene %in% c("ADAMTS9", "ADAMTS9-AS2"), .SD[which.min(p)], by = .(gene, dataset)][, .(source = paste(gene, dataset), type = "eQTL", lead = SNP, p)])
if (exists("ld")) lead_tab[, r2_rs4611812 := ifelse(lead %in% rownames(ld) & "rs4611812" %in% rownames(ld), ld[cbind(match(lead, rownames(ld)), match("rs4611812", rownames(ld)))], NA_real_)]
fwrite(lead_tab[order(type, p)], file.path(OUT, "lead_snps_ld.tsv"), sep = "\t")

# ---------------- key-SNP effects, aligned to the allele that raises aortic ADAMTS9 expression at rs4611812
key <- unique(c(KEYSNP, co[gene == "ADAMTS9" & dataset == "GTEx_artery_aorta" & trait %in% c("AMD_FG", "IS"), coloc_snp]))
ref <- eqs[dataset == "GTEx_artery_aorta" & gene == "ADAMTS9" & SNP %in% key]
ref[, up_allele := fifelse(beta > 0, toupper(ea), toupper(oa))]
upa <- setNames(ref$up_allele, ref$SNP)
ke <- eqs[SNP %in% names(upa) & gene %in% c("ADAMTS9", "ADAMTS9-AS2", "MAGI1")]
ke[, b_up := fifelse(toupper(ea) == upa[SNP], beta, -beta)]
kg <- gw[SNP %in% names(upa)]; kg[, b_up := fifelse(ea == upa[SNP], beta, -beta)]
kse <- rbind(ke[, .(SNP, up_allele = upa[SNP], type = "eQTL", source = paste(gene, dataset), b_up, se, p)],
             kg[, .(SNP, up_allele = upa[SNP], type = "GWAS", source = trait, b_up, se, p)])
fwrite(kse[order(SNP, type, p)], file.path(OUT, "key_snp_effects.tsv"), sep = "\t")

# ---------------- forest plot (rs4611812, per aortic-ADAMTS9-raising allele) + tissue eQTL panel
fp <- kse[SNP == "rs4611812"]
if (nrow(fp)) {
  g <- fp[type == "GWAS"][match(TRAITS, source)][!is.na(source)]
  e <- fp[type == "eQTL"][order(-abs(b_up / se))][1:min(.N, 20)]
  png(file.path(OUT, "Fig_forest_ADAMTS9.png"), width = 11, height = 6.5, units = "in", res = 200); pdf_also <- TRUE
  draw <- function() {
    par(mfrow = c(1, 2), mar = c(4, 11, 3, 1))
    q <- qnorm(0.975); y <- rev(seq_len(nrow(g)))
    plot(exp(g$b_up), y, xlim = range(exp(c(g$b_up - q * g$se, g$b_up + q * g$se)), na.rm = TRUE), log = "x", yaxt = "n",
         pch = ifelse(g$p < 0.05, 19, 1), col = ifelse(g$source %in% AMD, "#c0392b", "#2b5c8a"),
         xlab = "OR (or exp(beta)) per aortic ADAMTS9-raising allele", ylab = "", main = sprintf("rs4611812 (%s)", upa["rs4611812"]))
    segments(exp(g$b_up - q * g$se), y, exp(g$b_up + q * g$se), y, col = ifelse(g$source %in% AMD, "#c0392b", "#2b5c8a"))
    abline(v = 1, lty = 2); axis(2, y, sprintf("%s  P=%.1e", g$source, g$p), las = 1, cex.axis = 0.75)
    y2 <- rev(seq_len(nrow(e)))
    plot(e$b_up, y2, xlim = range(c(e$b_up - q * e$se, e$b_up + q * e$se)), yaxt = "n", pch = ifelse(e$p < 0.05, 19, 1),
         xlab = "eQTL beta (same allele)", ylab = "", main = "Tissue eQTL")
    segments(e$b_up - q * e$se, y2, e$b_up + q * e$se, y2); abline(v = 0, lty = 2)
    axis(2, y2, sprintf("%s  P=%.1e", sub("GTEx_", "", e$source), e$p), las = 1, cex.axis = 0.65)
  }
  draw(); dev.off(); pdf(file.path(OUT, "Fig_forest_ADAMTS9.pdf"), width = 11, height = 6.5); draw(); dev.off()
}

# ---------------- LocusZoom-style region plot (GRCh37), colour = r2 with rs4611812
system2(PLINK, c("--bfile", BREG, "--r2", "--ld-snp", "rs4611812", "--ld-window", "99999", "--ld-window-kb", "2000",
                 "--ld-window-r2", "0", "--out", file.path(TMP, "ld_rs4611812")), stdout = FALSE)
ldr <- fread(file.path(TMP, "ld_rs4611812.ld"))[, .(SNP = SNP_B, r2 = R2)]
tracks <- list(AMD_FG = gwl[["AMD_FG"]][, .(SNP, pos = pos37, p)], IS = gwl[["IS"]][, .(SNP, pos = pos37, p)],
               LACUNAR = gwl[["LACUNAR"]][, .(SNP, pos = pos37, p)],
               "ADAMTS9 eQTL, aorta" = eqs[dataset == "GTEx_artery_aorta" & gene == "ADAMTS9", .(SNP, pos = pos37, p)],
               "ADAMTS9 eQTL, fibroblast" = eqs[dataset == "GTEx_fibroblast" & gene == "ADAMTS9", .(SNP, pos = pos37, p)],
               "ADAMTS9-AS2 eQTL, fibroblast" = eqs[dataset == "GTEx_fibroblast" & gene == "ADAMTS9-AS2", .(SNP, pos = pos37, p)])
xl <- c(64.40e6, 64.90e6)
genes_draw <- gs[gene %in% c("ADAMTS9", "ADAMTS9-AS2", "MAGI1", "PRICKLE2")][, .(gene, s = start38 - OFF, e = end38 - OFF, strand)]
draw_region <- function() {
  par(mfrow = c(length(tracks) + 1, 1), mar = c(1, 5, 1.6, 1), oma = c(3, 0, 0, 0))
  pal <- colorRampPalette(c("#2b3a8a", "#36a2c4", "#2ca25f", "#f39c12", "#c0392b"))(5)
  for (k in names(tracks)) {
    d <- merge(tracks[[k]], ldr, by = "SNP", all.x = TRUE)[!is.na(pos) & pos >= xl[1] & pos <= xl[2]]
    d[is.na(r2), r2 := 0]; d <- d[order(r2)]
    plot(d$pos / 1e6, -log10(d$p), pch = 19, cex = 0.5, col = pal[pmin(5, floor(d$r2 * 5) + 1)], xlim = xl / 1e6,
         ylab = "-log10 P", xlab = "", main = k, cex.main = 0.9, xaxt = "n")
    lk <- d[SNP == "rs4611812"]; if (nrow(lk)) points(lk$pos / 1e6, -log10(lk$p), pch = 23, bg = "purple", cex = 1.3)
  }
  plot(NA, xlim = xl / 1e6, ylim = c(0, nrow(genes_draw) + 1), yaxt = "n", ylab = "", xlab = "chr3 (Mb, GRCh37)")
  for (i in seq_len(nrow(genes_draw))) {
    segments(genes_draw$s[i] / 1e6, i, genes_draw$e[i] / 1e6, i, lwd = 4, col = "grey40")
    text(max(xl[1], min(genes_draw$s[i], genes_draw$e[i])) / 1e6, i + 0.35, paste0(genes_draw$gene[i], ifelse(genes_draw$strand[i] == 1, " ->", " <-")),
         adj = 0, cex = 0.8)
  }
  mtext("chr3 (Mb, GRCh37); colour = r2 with rs4611812 (purple diamond)", side = 1, outer = TRUE, line = 1.5, cex = 0.8)
}
png(file.path(OUT, "Fig_region_3p14.png"), width = 8, height = 11, units = "in", res = 200); draw_region(); dev.off()
pdf(file.path(OUT, "Fig_region_3p14.pdf"), width = 8, height = 11); draw_region(); dev.off()
message("adamts9 done")
