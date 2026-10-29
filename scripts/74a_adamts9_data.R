# 3p14 ADAMTS9 deep-dive, data step (report section 4.2 follow-up):
#   region cis-eQTL summary statistics (all genes) from all 49 GTEx v8 tissues + ROSMAP / BrainSeq DLPFC (eQTL Catalogue,
#   remote tabix, GRCh38) and EyeGEx retina; GWAS summary statistics of 16 AMD / stroke / CSVD traits in the region;
#   region-restricted 1000G EUR plink files for LD / SMR.
# Output -> derived/adamts9/
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressMessages({ library(data.table); library(Rsamtools); library(parallel) })
OUT <- file.path(PROJ_ROOT, "derived/adamts9"); dir.create(file.path(OUT, "eqtl"), recursive = TRUE, showWarnings = FALSE)
CHR <- "3"; S38 <- 63900000; E38 <- 65400000; OFF <- 14324          # GRCh38 window; pos38 = pos37 + 14324 in this region
PLINK <- file.path(PROJ_ROOT, "data/raw/tools/plink/plink.exe"); BFILE <- file.path(PROJ_ROOT, "data/raw/ref/g1000_eur/g1000_eur")
G <- file.path(PROJ_ROOT, "derived/gwas")
TRAITS <- c("AMD_FG", "DRYAMD_FG", "WETAMD_FG", "AMD_IAMDGC", "AS", "IS", "LAS", "CES", "SVS", "WMH",
            "LACUNAR", "EPVS_WM", "EPVS_BG", "EPVS_HIP", "WM_MD", "WM_FA")

# ---------------- eQTL Catalogue: all GTEx tissues + DLPFC cohorts, all genes in the window
paths <- fread(file.path(PROJ_ROOT, "tmp/tissue/tabix_paths.tsv"))[quant_method == "ge"]
DS <- paths[study_label == "GTEx" | dataset_id %in% c("QTD000434", "QTD000051"),
            .(dataset = dataset_id, label = fifelse(study_label == "GTEx", paste0("GTEx_", sample_group), paste0(study_label, "_DLPFC")),
              n = sample_size, url = sub("^ftp://", "http://", ftp_path))]
DS[label == "ROSMAP_DLPFC" | grepl("ROSMAP", label), label := "ROSMAP_DLPFC"]
DS[grepl("BrainSeq", label), label := "BrainSeq_DLPFC"]
fwrite(DS, file.path(OUT, "eqtl_datasets.tsv"), sep = "\t")
COLS <- c("molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "ma_samples", "maf", "pvalue", "beta", "se",
          "type", "ac", "an", "r2", "molecular_trait_object_id", "gene_id", "median_tpm", "rsid")
fetch <- function(i) {
  f <- file.path(OUT, "eqtl", paste0(DS$label[i], ".tsv.gz"))
  if (file.exists(f)) return(paste("skip", DS$label[i]))
  for (attempt in 1:6) {
    res <- tryCatch({
      x <- scanTabix(TabixFile(DS$url[i]), param = GRanges(CHR, IRanges(S38, E38)))[[1]]
      if (length(x)) fread(text = paste(x, collapse = "\n"), header = FALSE, col.names = COLS, sep = "\t")[type == "SNP"] else data.table()
    }, error = function(e) NULL)
    if (!is.null(res)) break
    Sys.sleep(20 * attempt)
  }
  if (is.null(res)) return(paste("FAILED", DS$label[i]))
  res[, `:=`(dataset = DS$label[i], n = DS$n[i])]
  fwrite(res, f, sep = "\t")
  sprintf("%s %s: %d rows, %d genes", format(Sys.time(), "%H:%M:%S"), DS$label[i], nrow(res), uniqueN(res$gene_id))
}
cl <- makeCluster(8)
clusterExport(cl, c("DS", "OUT", "COLS", "CHR", "S38", "E38"))
clusterEvalQ(cl, { .libPaths(c(file.path(Sys.getenv("COMORBIDITY_ROOT"), "Rlib"), .libPaths())); suppressMessages({ library(data.table); library(Rsamtools) }) })
for (m in parLapplyLB(cl, seq_len(nrow(DS)), fetch)) message(m)
stopCluster(cl)

# ---------------- EyeGEx retina (QTLtools nominal; SE back-calculated from P), all genes in the window
ey <- fread(file.path(PROJ_ROOT, "data/raw/eqtl/eyegex/retina_nominal_chr3.txt.gz"), header = TRUE, select = c(1, 8, 9, 10, 12, 13, 14, 15),
            col.names = c("gene_id", "vid", "vchr", "vpos", "ea", "oa", "p", "beta"))
ey <- ey[vpos >= S38 & vpos <= E38]
ey[, gene_id := sub("\\..*", "", gene_id)]; ey[, SNP := sub(":.*", "", vid)]; ey <- ey[grepl("^rs", SNP)]
ey[, se := abs(beta) / qnorm(p / 2, lower.tail = FALSE)]
fwrite(ey, file.path(OUT, "eqtl", "EyeGEx_Retina.tsv.gz"), sep = "\t")
message("EyeGEx rows ", nrow(ey), ", genes ", uniqueN(ey$gene_id))

# ---------------- GWAS region (SNP set: FinnGen rsIDs in the GRCh38 window; GRCh37 traits matched by rsID)
fg <- as.data.table(readRDS(file.path(G, "AMD_FG.rds")))
snps <- unique(fg[chr == CHR & pos >= S38 & pos <= E38 & grepl("^rs", SNP), .(SNP, pos38 = pos)], by = "SNP"); rm(fg); gc()
bim <- fread(paste0(BFILE, ".bim"), header = FALSE, select = c(1, 2, 4, 5, 6), col.names = c("chr", "SNP", "pos37", "A1", "A2"))
bim <- bim[chr == 3 & pos37 >= S38 - OFF - 50000 & pos37 <= E38 - OFF + 50000]
snps <- merge(snps, bim[, .(SNP, pos37, A1, A2)], by = "SNP", all.x = TRUE)
fwrite(snps, file.path(OUT, "region_snps.tsv"), sep = "\t")
cl <- makeCluster(8); clusterExport(cl, c("G", "snps"))
clusterEvalQ(cl, { .libPaths(c(file.path(Sys.getenv("COMORBIDITY_ROOT"), "Rlib"), .libPaths())); library(data.table) })
gw <- rbindlist(parLapplyLB(cl, TRAITS, function(k) {
  d <- as.data.table(readRDS(file.path(G, paste0(k, ".rds"))))[SNP %in% snps$SNP]
  d[, trait := k]; d
}), fill = TRUE)
stopCluster(cl)
gw <- merge(gw, snps[, .(SNP, pos38, pos37)], by = "SNP")
saveRDS(gw, file.path(OUT, "gwas_region.rds"))
print(gw[, .(n_snps = .N, min_p = min(p)), by = trait])

# ---------------- region plink files (GRCh37)
system2(PLINK, c("--bfile", BFILE, "--chr", "3", "--from-bp", S38 - OFF - 500000, "--to-bp", E38 - OFF + 500000,
                 "--make-bed", "--out", file.path(OUT, "g1000_eur_3p14")), stdout = FALSE)
message("data step done")
