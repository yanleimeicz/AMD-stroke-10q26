# Round 6, item C: FinnGen R12 CAD endpoints (I9_CHD, I9_CORATHER) at 10q26 by remote tabix (GRCh38)
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
suppressMessages({ library(data.table); library(Rsamtools) })
P <- PROJ_ROOT; O <- file.path(P, "derived/revision6")
FG <- "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_"
out <- rbindlist(lapply(c("I9_CHD", "I9_CORATHER"), function(e) {
  tf <- TabixFile(paste0(FG, e, ".gz"), index = paste0(FG, e, ".gz.tbi"))
  x <- scanTabix(tf, param = GRanges("10", IRanges(122440000, 122480000)))[[1]]
  d <- fread(text = paste(x, collapse = "\n"), header = FALSE)
  setnames(d, c("chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "mlogp", "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls"))
  d[, endpoint := e]; d }))
fwrite(out, file.path(O, "finngen_cad_10q26_region.tsv"), sep = "\t")
print(out[rsids %in% c("rs3750848", "rs2672592", "rs72631113")])
