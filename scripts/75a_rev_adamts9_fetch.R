# Revision (editorial items R4/R10/S6), data step: region summary statistics by remote tabix (GRCh38)
#   3p14 (ADAMTS9) chr3:63.9-65.4 Mb and 10q26 chr10:121.9-123.0 Mb
#   FinnGen R12 endpoints (T2D, stroke, ICH, CHD, coronary atherosclerosis, hypertension) and
#   GWAS Catalog harmonised files with tabix index (CAD Aragam 2022, GLGC 2021 LDL/HDL, GIGASTROKE East Asian).
# Output -> derived/revision_adamts9/<source>.tsv.gz (SNP chr pos ea oa beta se p eaf n ncase)
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
OUT <- file.path(PROJ_ROOT, "derived/revision_adamts9"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
REG <- list(c("3", 63900000, 65400000), c("10", 121900000, 123000000))
FG <- "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_"
GC <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
fgn <- fread(file.path(PROJ_ROOT, "tmp/rev_adamts9/finngen_counts.tsv"))       # phenocode n_case n_control (from PheWeb API)
SRC <- rbind(
  data.table(name = paste0("FG_", fgn$phenocode), url = paste0(FG, fgn$phenocode, ".gz"), fmt = "finngen",
             ncase = fgn$n_case, nctrl = fgn$n_control),
  data.table(name = c("CAD_Aragam2022", "LDL_GLGC2021", "HDL_GLGC2021", "EAS_AS", "EAS_IS", "EAS_LAS", "EAS_SVS"),
             url = paste0(GC, c("GCST90132001-GCST90133000/GCST90132314/harmonised/GCST90132314.h.tsv.gz",
                                "GCST90239001-GCST90240000/GCST90239658/harmonised/GCST90239658.h.tsv.gz",
                                "GCST90239001-GCST90240000/GCST90239652/harmonised/GCST90239652.h.tsv.gz",
                                "GCST90104001-GCST90105000/GCST90104544/harmonised/GCST90104544.h.tsv.gz",
                                "GCST90104001-GCST90105000/GCST90104545/harmonised/GCST90104545.h.tsv.gz",
                                "GCST90104001-GCST90105000/GCST90104547/harmonised/GCST90104547.h.tsv.gz",
                                "GCST90104001-GCST90105000/GCST90104548/harmonised/GCST90104548.h.tsv.gz")),
             fmt = c("cad", "glgc", "glgc_hdl", "gc", "gc", "gc", "gc"),
             ncase = c(181522, NA, NA, 27413, 19032, 1735, 5532), nctrl = c(984168, NA, NA, 237242, 237242, 237242, 237242)))
HDR <- list(finngen = c("chr", "pos", "oa", "ea", "SNP", "genes", "p", "mlogp", "beta", "se", "eaf", "eaf_case", "eaf_ctrl"),
            cad = c("chr", "pos", "ea", "oa", "beta", "se", "eaf", "p", "variant_id", "markername", "freqse", "minfreq", "maxfreq",
                    "direction", "hetisq", "hetchisq", "hetdf", "hetpval", "cases", "effective_cases", "n", "meta", "conv", "or", "hm_code", "SNP"),
            glgc = c("chr", "pos", "ea", "oa", "beta", "se", "eaf", "p", "variant_id", "n", "N_studies", "conv", "hm_code", "SNP"),
            glgc_hdl = c("chr", "pos", "ea", "oa", "beta", "se", "eaf", "p", "SNP", "n", "nstudy", "mlog10p", "mlog10p_gc", "p_gc", "conv", "hm_code", "variant_id"),
            gc = c("chr", "pos", "ea", "oa", "beta", "se", "eaf", "p", "variant_id", "ci_upper", "ci_lower", "conv", "or", "hm_code", "SNP"))
fetch <- function(i) {
  s <- SRC[i]; f <- file.path(OUT, paste0(s$name, ".tsv.gz"))
  if (file.exists(f)) return(paste("skip", s$name))
  res <- NULL
  for (attempt in 1:6) {
    res <- tryCatch({
      tf <- TabixFile(s$url)
      x <- unlist(lapply(REG, function(r) scanTabix(tf, param = GRanges(r[1], IRanges(as.numeric(r[2]), as.numeric(r[3]))))[[1]]))
      d <- fread(text = paste(x, collapse = "\n"), header = FALSE, sep = "\t", na.strings = c("NA", ""))
      setnames(d, HDR[[s$fmt]][seq_len(ncol(d))]); d
    }, error = function(e) NULL)
    if (!is.null(res)) break
    Sys.sleep(15 * attempt)
  }
  if (is.null(res)) return(paste("FAILED", s$name))
  res[, SNP := sub(",.*", "", as.character(SNP))]
  if (!"n" %in% names(res)) res[, n := s$ncase + s$nctrl]
  res[, `:=`(chr = as.character(chr), ncase = s$ncase, source = s$name)]
  keep <- intersect(c("SNP", "chr", "pos", "ea", "oa", "beta", "se", "p", "eaf", "n", "ncase", "source"), names(res))
  res <- res[!is.na(beta) & !is.na(se) & se > 0, ..keep]
  fwrite(res, f, sep = "\t")
  sprintf("%s %s: %d rows", format(Sys.time(), "%H:%M:%S"), s$name, nrow(res))
}
cl <- makeCluster(8)
clusterExport(cl, c("SRC", "OUT", "REG", "HDR"))
clusterEvalQ(cl, { .libPaths(c(file.path(Sys.getenv("COMORBIDITY_ROOT"), "Rlib"), .libPaths())); suppressMessages({ library(data.table); library(Rsamtools) }) })
for (m in parLapplyLB(cl, seq_len(nrow(SRC)), fetch)) message(m)
stopCluster(cl)
message("fetch done")
