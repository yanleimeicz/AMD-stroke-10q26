# Revision R4: FinnGen R12 in-sample SuSiE fine-mapping (lbf_variable1-10) at 3p14 for AMD and cardiometabolic / vascular
# endpoints, by remote tabix of the public SUSIE.snp.bgz files. Output -> derived/revision_adamts9/susie_fg_<trait>.tsv.gz
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
OUT <- file.path(PROJ_ROOT, "derived/revision_adamts9")
U <- "https://storage.googleapis.com/finngen-public-data-r12/finemap/full/susie/finngen_R12_"
TR <- c("H7_AMD", "WET_AMD", "DRY_AMD", "T2D", "I9_STR", "C_STROKE", "I9_CORATHER", "I9_CHD", "I9_HYPTENS", "I9_MI_STRICT",
        "E4_DM2", "T2D_WIDE", "KELA_DIAB_INSUL_EXMORE", "I9_REVASC", "I9_HYPTENSESS")
HDR <- c("trait", "region", "v", "rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se", "p", "mean", "sd",
         "prob", "cs", "cs_specific_prob", "low_purity", "lead_r2", "mean_99", "sd_99", "prob_99", "cs_99", "cs_specific_prob_99",
         "low_purity_99", "lead_r2_99", paste0("alpha", 1:10), paste0("mean", 1:10), paste0("sd", 1:10), paste0("lbf_variable", 1:10))
f1 <- function(t) {
  f <- file.path(OUT, paste0("susie_fg_", t, ".tsv.gz"))
  if (file.exists(f)) return(paste("skip", t))
  for (seqn in c("chr3", "3")) {
    x <- tryCatch(scanTabix(TabixFile(paste0(U, t, ".SUSIE.snp.bgz")), param = GRanges(seqn, IRanges(64213769, 65213769)))[[1]],
                  error = function(e) NULL)
    if (!is.null(x)) break
  }
  if (is.null(x)) return(paste("no file / error", t))
  if (!length(x)) return(paste("no fine-mapped region", t))
  d <- fread(text = paste(x, collapse = "\n"), header = FALSE, sep = "\t", na.strings = c("NA", ""))
  setnames(d, HDR[seq_len(ncol(d))])
  fwrite(d, f, sep = "\t")
  sprintf("%s: %d rows, regions %s, CS %s", t, nrow(d), paste(unique(d$region), collapse = ";"),
          paste(sort(unique(d$cs[d$cs > 0])), collapse = ","))
}
cl <- makeCluster(8)
clusterExport(cl, c("OUT", "U", "HDR"))
clusterEvalQ(cl, { .libPaths(c(file.path(Sys.getenv("COMORBIDITY_ROOT"), "Rlib"), .libPaths())); suppressMessages({ library(data.table); library(Rsamtools) }) })
for (m in parLapplyLB(cl, TR, f1)) message(m)
stopCluster(cl)
message("susie fetch done")
