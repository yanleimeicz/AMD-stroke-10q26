# Multi-tissue eQTL evidence for the AMD-stroke shared loci (report section 4 / 5.1 follow-up):
#   key-variant lookup across tissues, coloc.abf (eQTL vs AMD / stroke / WMH GWAS), SMR (top cis-eQTL Wald ratio).
# Sources: eQTL Catalogue (GTEx v8 arteries, brain, blood, fibroblast; ROSMAP / BrainSeq DLPFC), eQTLGen blood (BESD),
#          EyeGEx retina (Ratnapriya 2019, Zenodo 18331341; nominal QTLtools output, SE back-calculated from P).
# Inputs from 48a_fetch_tissue_eqtl.R and 48b_gwas_regions.R. Output -> results/tissue_eqtl/
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressMessages({ library(data.table); library(parallel) })
D <- file.path(PROJ_ROOT, "derived/tissue_eqtl"); OUT <- file.path(PROJ_ROOT, "results/tissue_eqtl")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
REG <- fread(file.path(D, "regions_hg38.tsv"), colClasses = c(chr = "character"))
GENES <- fread(file.path(D, "target_genes.tsv"))
rsnps <- fread(file.path(D, "region_snps.tsv"))
gw <- readRDS(file.path(D, "gwas_regions.rds"))
frq <- fread(file.path(PROJ_ROOT, "derived/g1000_eur.frq"))[, .(SNP, fA1 = A1, fA2 = A2, MAF)]

# ---------------- harmonised eQTL table: dataset, gene, region, SNP, ea, oa, beta, se, p, maf, n
eq <- list()
ds <- fread(file.path(D, "datasets.tsv"))
for (f in list.files(D, pattern = "__.*\\.tsv\\.gz$", full.names = TRUE)) {
  d <- fread(f); if (!nrow(d)) next
  reg <- sub(".*__(.*)\\.tsv\\.gz$", "\\1", basename(f))
  eq[[f]] <- d[, .(dataset, gene_id, region = reg, SNP = rsid, ea = alt, oa = ref, beta, se, p = pvalue, maf, n, pos38 = position)]
}
for (f in list.files(D, pattern = "^eQTLGen_Blood__.*\\.txt$", full.names = TRUE)) {
  d <- fread(f); if (!nrow(d)) next
  parts <- strsplit(sub("\\.txt$", "", basename(f)), "__")[[1]]
  eq[[f]] <- d[, .(dataset = "eQTLGen_Blood", gene_id = Probe, region = parts[2], SNP, ea = A1, oa = A2, beta = b, se = SE, p,
                   maf = pmin(Freq, 1 - Freq), n = 31684, pos38 = NA_real_)]
}
# EyeGEx retina: QTLtools nominal columns (README); variant IDs = rsID[:pos37:ref:alt]; positions GRCh38
ey_dir <- file.path(PROJ_ROOT, "data/raw/eqtl/eyegex")
for (j in seq_len(nrow(REG))) {
  f <- file.path(ey_dir, sprintf("retina_nominal_chr%s.txt.gz", REG$chr[j]))
  if (!file.exists(f)) { message("EyeGEx missing chr", REG$chr[j]); next }
  d <- fread(f, header = TRUE, select = c(1, 8, 9, 10, 12, 13, 14, 15),
             col.names = c("gene_id", "vid", "vchr", "vpos", "ea", "oa", "p", "beta"))
  d[, gene_id := sub("\\..*", "", gene_id)]
  d <- d[gene_id %in% GENES[region == REG$region[j]]$gene_id & vpos >= REG$start[j] & vpos <= REG$end[j]]
  d[, SNP := sub(":.*", "", vid)]; d <- d[grepl("^rs", SNP)]
  d[, se := abs(beta) / qnorm(p / 2, lower.tail = FALSE)]
  d <- merge(d, frq[, .(SNP, MAF)], by = "SNP", all.x = TRUE)
  eq[[f]] <- d[is.finite(se) & se > 0, .(dataset = "EyeGEx_Retina", gene_id, region = REG$region[j], SNP, ea, oa, beta, se, p,
                                          maf = MAF, n = 406, pos38 = vpos)]
}
eq <- rbindlist(eq, fill = TRUE)
eq <- merge(eq, GENES[, .(gene_id, gene)], by = "gene_id")
eq <- eq[!is.na(se) & se > 0 & !is.na(SNP) & SNP != ""]
eq <- unique(eq, by = c("dataset", "gene_id", "SNP"))
saveRDS(eq, file.path(D, "eqtl_harmonised.rds"))
cov <- eq[, .(n_snps = .N, min_p = min(p), top_snp = SNP[which.min(p)], n = n[1]), by = .(region, gene, dataset)]
fwrite(cov[order(region, gene, min_p)], file.path(OUT, "eqtl_dataset_gene_summary.tsv"), sep = "\t")
message("eQTL rows: ", nrow(eq), "; datasets: ", paste(unique(eq$dataset), collapse = ", "))

# ---------------- key-variant lookup (betas aligned to a stated allele)
fg <- gw[trait == "AMD_FG"]
risk <- function(s) { r <- fg[SNP == s]; if (!nrow(r)) NA_character_ else if (r$beta > 0) r$ea else r$oa }
KEY <- data.table(SNP = c("rs3750848", "rs36212733", "rs10490924", "rs11200638", "rs61871745", "rs2672592", "rs3184504",
                          "rs429358", "rs7412", "rs4611812", "rs42036"))
KEY[, allele := sapply(SNP, risk)]
KEY[SNP == "rs3184504", allele := "T"]; KEY[SNP == "rs429358", allele := "C"]; KEY[SNP == "rs7412", allele := "T"]
KEY[, note := fifelse(SNP %in% c("rs3184504"), "SH2B3 p.R262W: T = W262",
              fifelse(SNP == "rs429358", "C = APOE e4", fifelse(SNP == "rs7412", "T = APOE e2", "allele = AMD risk allele (FinnGen)")))]
align <- function(d) {
  d <- merge(d, KEY[, .(SNP, allele)], by = "SNP")
  d[, b_allele := fifelse(toupper(ea) == allele, beta, fifelse(toupper(oa) == allele, -beta, NA_real_))]; d
}
lk_eq <- align(eq[SNP %in% KEY$SNP])[, .(region, gene, dataset, SNP, allele, b_allele, se, p, n)]
lk_gw <- align(gw[SNP %in% KEY$SNP, .(SNP, trait, ea, oa, beta, se, p)])[, .(SNP, allele, trait, b_allele, se, p)]
fwrite(KEY, file.path(OUT, "key_variants.tsv"), sep = "\t")
fwrite(lk_eq[order(SNP, gene, p)], file.path(OUT, "key_variant_eqtl_lookup.tsv"), sep = "\t")
fwrite(lk_gw[order(SNP, trait)], file.path(OUT, "key_variant_gwas_lookup.tsv"), sep = "\t")

# ---------------- coloc.abf: every (eQTL dataset x target gene) vs every trait in the same region
TRAITS <- c("AMD_FG", "AMD_IAMDGC", "DRYAMD_FG", "WETAMD_FG", "AS", "IS", "LAS", "SVS", "WMH")
jobs <- CJ(key = unique(eq[, paste(region, gene, dataset, sep = "|")]), trait = TRAITS)
cl <- makeCluster(30)
clusterExport(cl, c("eq", "gw", "jobs"))
clusterEvalQ(cl, { .libPaths(c(file.path(Sys.getenv("COMORBIDITY_ROOT"), "Rlib"), .libPaths())); suppressMessages({ library(data.table); library(coloc) }) })
res <- parLapply(cl, seq_len(nrow(jobs)), function(i) {
  k <- strsplit(jobs$key[i], "\\|")[[1]]; tr <- jobs$trait[i]
  e <- eq[region == k[1] & gene == k[2] & dataset == k[3]]
  g <- gw[trait == tr & region == k[1]]
  m <- merge(e, g[, .(SNP, gea = toupper(ea), goa = toupper(oa), gb = beta, gse = se, gp = p, geaf = eaf, gn = n, ncase, ncontrol)], by = "SNP")
  m <- m[(toupper(ea) == gea & toupper(oa) == goa) | (toupper(ea) == goa & toupper(oa) == gea)]
  m[, gb_al := fifelse(toupper(ea) == gea, gb, -gb)]
  m <- m[is.finite(gb) & is.finite(gse) & gse > 0 & !is.na(maf) & maf > 0.01 & maf < 0.99]
  if (nrow(m) < 50) return(NULL)
  d1 <- list(beta = m$beta, varbeta = m$se^2, snp = m$SNP, type = "quant", N = m$n[1], MAF = m$maf)
  if (tr == "WMH") d2 <- list(beta = m$gb_al, varbeta = m$gse^2, snp = m$SNP, type = "quant", N = max(m$gn, na.rm = TRUE), MAF = m$maf)
  else { s <- m$ncase[1] / (m$ncase[1] + m$ncontrol[1]); d2 <- list(beta = m$gb_al, varbeta = m$gse^2, snp = m$SNP, type = "cc", s = s, MAF = m$maf) }
  r <- tryCatch(suppressWarnings(suppressMessages(coloc.abf(d1, d2))), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  top <- r$results[which.max(r$results$SNP.PP.H4), "snp"]; mt <- m[SNP == top]
  data.table(region = k[1], gene = k[2], dataset = k[3], trait = tr, nsnps = nrow(m),
             min_p_eqtl = min(m$p), min_p_gwas = min(m$gp), top_snp_eqtl = m$SNP[which.min(m$p)], top_snp_gwas = m$SNP[which.min(m$gp)],
             PP.H3 = r$summary[["PP.H3.abf"]], PP.H4 = r$summary[["PP.H4.abf"]], coloc_snp = top,
             coloc_snp_PP = max(r$results$SNP.PP.H4),
             dir_expr_to_risk = sign(mt$beta * mt$gb_al))   # +1: allele raising expression raises risk
})
stopCluster(cl)
co <- rbindlist(res)
co[, coloc_class := fifelse(PP.H4 > 0.8, "strong", fifelse(PP.H4 > 0.5, "moderate", fifelse(PP.H3 > 0.8, "distinct variants", "no/uncertain")))]
fwrite(co[order(region, gene, trait, -PP.H4)], file.path(OUT, "coloc_tissue_eqtl_gwas.tsv"), sep = "\t")

# ---------------- SMR (Wald ratio at the top cis-eQTL, P_eQTL < 1e-5; chi2_SMR = zG^2 zE^2 / (zG^2 + zE^2))
smr <- rbindlist(lapply(split(eq, by = c("region", "gene", "dataset")), function(e) {
  t <- e[which.min(p)]; if (t$p > 1e-5) return(NULL)
  g <- gw[SNP == t$SNP & region == t$region]
  if (!nrow(g)) return(NULL)
  g[, gb_al := fifelse(toupper(ea) == toupper(t$ea), beta, fifelse(toupper(oa) == toupper(t$ea), -beta, NA_real_))]
  g <- g[!is.na(gb_al)]
  zE <- t$beta / t$se; zG <- g$gb_al / g$se
  chi <- zG^2 * zE^2 / (zG^2 + zE^2); b <- g$gb_al / t$beta
  data.table(region = t$region, gene = t$gene, dataset = t$dataset, top_eqtl = t$SNP, p_eqtl = t$p, trait = g$trait,
             b_SMR = b, se_SMR = abs(b) / sqrt(chi), OR_SMR = exp(b), p_SMR = pchisq(chi, 1, lower.tail = FALSE))
}))
smr[, p_fdr := p.adjust(p_SMR, "BH"), by = trait]
fwrite(smr[order(region, gene, trait, p_SMR)], file.path(OUT, "smr_tissue_topeqtl.tsv"), sep = "\t")
message("coloc pairs ", nrow(co), "; strong ", sum(co$PP.H4 > 0.8), "; SMR rows ", nrow(smr))
message("tissue eQTL done")
