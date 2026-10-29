# Post-hoc Steiger directionality for every MR pair using liability-scale r for binary traits.
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(TwoSampleMR) })
OUT <- file.path(PROJ_ROOT, "results/mr")
DAT <- readRDS(file.path(OUT, "harmonised_data.rds"))
PREV <- c(AMD_IAMDGC = 0.02, AMD_FG = 0.03, DRYAMD_FG = 0.02, WETAMD_FG = 0.015, AS = 0.03, IS = 0.025, CES = 0.005, LAS = 0.005, SVS = 0.005)
r_of <- function(dat, side) {
  k <- dat[[side]][1]; b <- dat[[paste0("beta.", side)]]; se <- dat[[paste0("se.", side)]]
  eaf <- dat[[paste0("eaf.", side)]]; other <- dat[[paste0("eaf.", ifelse(side == "exposure", "outcome", "exposure"))]]
  eaf[is.na(eaf)] <- other[is.na(eaf)]; eaf[is.na(eaf)] <- 0.3
  if (k %in% names(PREV)) get_r_from_lor(b, eaf, dat[[paste0("ncase.", side)]][1], dat[[paste0("ncontrol.", side)]][1], PREV[[k]])
  else get_r_from_bsen(b, se, dat[[paste0("samplesize.", side)]])
}
res <- rbindlist(lapply(names(DAT), function(key) {
  d <- DAT[[key]]; if (is.null(d) || !nrow(d)) return(NULL)
  rx <- r_of(d, "exposure"); ry <- r_of(d, "outcome")
  st <- mr_steiger(p_exp = d$pval.exposure, p_out = d$pval.outcome, n_exp = d$samplesize.exposure, n_out = d$samplesize.outcome,
                   r_exp = rx, r_out = ry)
  data.table(exposure = d$exposure[1], outcome = d$outcome[1], nsnp = nrow(d), r2_exposure = st$r2_exp, r2_outcome = st$r2_out,
             correct_direction = st$correct_causal_direction, steiger_p = st$steiger_test)
}))
fwrite(res, file.path(OUT, "mr_steiger_liability.tsv"), sep = "\t"); print(res)
