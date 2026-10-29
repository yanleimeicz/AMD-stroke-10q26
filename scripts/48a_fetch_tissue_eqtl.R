# Fetch region-level cis-eQTL summary statistics for the shared loci from:
#   eQTL Catalogue (GTEx v8 arteries / brain / blood / fibroblast, ROSMAP & BrainSeq DLPFC) via remote tabix (GRCh38),
#   eQTLGen whole blood via the local BESD (smr --query, GRCh37 positions, rsIDs).
# EyeGEx retina (Zenodo) is read in 48_tissue_eqtl.R. Output -> derived/tissue_eqtl/
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressMessages({ library(data.table); library(Rsamtools) })
OUT <- file.path(PROJ_ROOT, "derived/tissue_eqtl"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

REGIONS <- data.table(   # GRCh38 windows
  region = c("10q26_HTRA1", "12q24_SH2B3", "19q13_APOE", "6p21_C4", "3p14_ADAMTS9", "9q34_ABO", "7q21_CDK6"),
  chr    = c("10", "12", "19", "6", "3", "9", "7"),
  start  = c(121955799, 110946804, 44658684, 31732539, 64350000, 133000000, 92450000),
  end    = c(122971234, 112005200, 45158684, 32232539, 64850000, 133520000, 92980000))
TARGETS <- list(
  "10q26_HTRA1"  = c(HTRA1 = "ENSG00000166033", ARMS2 = "ENSG00000254636", PLEKHA1 = "ENSG00000107679"),
  "12q24_SH2B3"  = c(SH2B3 = "ENSG00000111252", ATXN2 = "ENSG00000204842", ALDH2 = "ENSG00000111275",
                     BRAP = "ENSG00000089234", ACAD10 = "ENSG00000111271"),
  "19q13_APOE"   = c(APOE = "ENSG00000130203", APOC1 = "ENSG00000130208", TOMM40 = "ENSG00000130204", NECTIN2 = "ENSG00000130202"),
  "6p21_C4"      = c(C4A = "ENSG00000244731", C4B = "ENSG00000224389", CYP21A2 = "ENSG00000231852", TNXB = "ENSG00000168477"),
  "3p14_ADAMTS9" = c(ADAMTS9 = "ENSG00000163638", "ADAMTS9-AS2" = "ENSG00000241684"),
  "9q34_ABO"     = c(ABO = "ENSG00000175164"),
  "7q21_CDK6"    = c(CDK6 = "ENSG00000105810"))
gene_map <- rbindlist(lapply(names(TARGETS), function(r) data.table(region = r, gene = names(TARGETS[[r]]), gene_id = TARGETS[[r]])))
fwrite(REGIONS, file.path(OUT, "regions_hg38.tsv"), sep = "\t"); fwrite(gene_map, file.path(OUT, "target_genes.tsv"), sep = "\t")

DS <- data.table(
  dataset = c("QTD000131", "QTD000136", "QTD000141", "QTD000171", "QTD000176", "QTD000156", "QTD000196", "QTD000181",
              "QTD000166", "QTD000356", "QTD000216", "QTD000434", "QTD000051"),
  label = c("GTEx_Artery_Aorta", "GTEx_Artery_Coronary", "GTEx_Artery_Tibial", "GTEx_Brain_Cortex", "GTEx_Brain_FrontalCortex",
            "GTEx_Brain_Caudate", "GTEx_Brain_Putamen", "GTEx_Brain_Hippocampus", "GTEx_Brain_Cerebellum", "GTEx_Whole_Blood",
            "GTEx_Fibroblast", "ROSMAP_DLPFC", "BrainSeq_DLPFC"))
paths <- fread(file.path(PROJ_ROOT, "tmp/tissue/tabix_paths.tsv"))[quant_method == "ge"]
DS <- merge(DS, paths[, .(dataset = dataset_id, n = sample_size, ftp_path)], by = "dataset", sort = FALSE)
DS[, url := sub("^ftp://", "http://", ftp_path)]
fwrite(DS[, .(dataset, label, n, url)], file.path(OUT, "datasets.tsv"), sep = "\t")
COLS <- c("molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "ma_samples", "maf", "pvalue", "beta", "se",
          "type", "ac", "an", "r2", "molecular_trait_object_id", "gene_id", "median_tpm", "rsid")

for (i in seq_len(nrow(DS))) for (j in seq_len(nrow(REGIONS))) {
  f <- file.path(OUT, sprintf("%s__%s.tsv.gz", DS$label[i], REGIONS$region[j]))
  if (file.exists(f)) next
  ids <- TARGETS[[REGIONS$region[j]]]
  for (attempt in 1:5) {
    res <- tryCatch({
      x <- scanTabix(TabixFile(DS$url[i]), param = GRanges(REGIONS$chr[j], IRanges(REGIONS$start[j], REGIONS$end[j])))[[1]]
      d <- if (length(x)) fread(text = paste(x, collapse = "\n"), header = FALSE, col.names = COLS, sep = "\t") else data.table()
      if (nrow(d)) d <- d[gene_id %in% ids & type == "SNP"]
      d
    }, error = function(e) { message("retry ", attempt, " ", DS$label[i], " ", REGIONS$region[j], ": ", conditionMessage(e)); NULL })
    if (!is.null(res)) break
    Sys.sleep(30 * attempt)
  }
  if (is.null(res)) { message("FAILED ", f); next }
  res[, `:=`(dataset = DS$label[i], n = DS$n[i])]
  fwrite(res, f, sep = "\t")
  message(format(Sys.time(), "%H:%M:%S"), " ", DS$label[i], " ", REGIONS$region[j], ": ", nrow(res), " rows, genes ", uniqueN(res$gene_id))
  Sys.sleep(2)
}

# eQTLGen whole blood: full cis summary statistics of each target probe from the BESD
SMR <- file.path(PROJ_ROOT, "data/raw/tools/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe")
BESD <- file.path(PROJ_ROOT, "data/raw/eqtl/eqtlgen_besd")
epi <- fread(paste0(BESD, ".epi"), header = FALSE)
gm <- gene_map[gene_id %in% epi$V2]
todo <- gm[!file.exists(file.path(OUT, sprintf("eQTLGen_Blood__%s__%s.txt", region, gene)))]
if (nrow(todo)) {    # one pass over the BESD for all target probes
  pf <- file.path(OUT, "eqtlgen_probes.txt"); fwrite(data.table(todo$gene_id), pf, col.names = FALSE)
  system2(SMR, c("--beqtl-summary", BESD, "--query", "1", "--extract-probe", pf, "--out", file.path(OUT, "eqtlgen_targets")),
          stdout = FALSE, stderr = FALSE)
  q <- fread(file.path(OUT, "eqtlgen_targets.txt"))
  for (k in seq_len(nrow(todo))) {
    fwrite(q[Probe == todo$gene_id[k]], file.path(OUT, sprintf("eQTLGen_Blood__%s__%s.txt", todo$region[k], todo$gene[k])), sep = "\t")
    message("eQTLGen ", todo$gene[k], ": ", nrow(q[Probe == todo$gene_id[k]]))
  }
}
message("fetch done")
