# Revision C, step 1b: MR-RAPS (over-dispersed, Huber loss) for schemes A_dist and C_orig, appended to mr_revised_all.tsv.
# Re-uses the data preparation of rev_mr_01_mr.R (sections 1-5) so that instruments, proxies and harmonisation are identical.
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
src <- readLines(file.path(PROJ_ROOT, "scripts/rev_mr_01_mr.R"))
eval(parse(text = src[seq_len(grep("^RES <- list\\(\\)", src)[1] - 1)]))
suppressPackageStartupMessages(library(mr.raps))
RA <- list()
for (e in names(INST)) for (s in c("A_dist", "C_orig")) {
  ex <- INST[[e]][[s]]; if (!nrow(ex)) next
  ex[, locus := locus_of(chr, pos)]
  for (o in OUTC) {
    m0 <- harm(ex, OD[[o]]); m0 <- merge(m0, ex[, .(SNP, locus)], by = "SNP")
    for (xn in names(EXCL)) {
      m <- m0[!locus %in% EXCL[[xn]]]; if (nrow(m) < 3) next
      r <- tryCatch(mr.raps.overdispersed.robust(m$bx, m$by, m$sx, m$sy, loss.function = "huber", suppress.warning = TRUE), error = function(e) NULL)
      if (is.null(r)) { message("RAPS failed: ", e, " ", s, " ", o, " ", xn); next }
      RA[[length(RA) + 1]] <- data.table(exposure = e, scheme = s, outcome = o, exclusion = xn, method = "RAPS", b = r$beta.hat, se = r$beta.se,
                                         p = 2 * pnorm(-abs(r$beta.hat / r$beta.se)), nsnp = nrow(m), n_removed = nrow(m0) - nrow(m))
    }
  }
}
ra <- rbindlist(RA)
ra[, binary := outcome %in% BINARY]
ra[, `:=`(OR = ifelse(binary, exp(b), NA), lo = ifelse(binary, exp(b - 1.96 * se), NA), hi = ifelse(binary, exp(b + 1.96 * se), NA))]
res <- fread(file.path(OUT, "mr_revised_all.tsv"))[method != "RAPS"]
res <- rbind(res, ra, fill = TRUE)
fwrite(res, file.path(OUT, "mr_revised_all.tsv"), sep = "\t")
# consistency check against the backup made before the rerun of rev_mr_01
bk <- file.path(PROJ_ROOT, "tmp/revmr/bak/mr_revised_all.tsv")
if (file.exists(bk)) { a <- fread(bk); k <- c("exposure", "scheme", "outcome", "exclusion", "method"); mm <- merge(a, res, by = k)
  message(sprintf("rows matched with backup: %d; max |delta b| = %.3g", nrow(mm), max(abs(mm$b.x - mm$b.y)))) }
iv <- res[method == "IVW", .(exposure, scheme, outcome, exclusion, OR_ivw = round(exp(b), 3), p_ivw = signif(p, 2))]
print(merge(ra[scheme == "A_dist" & exclusion %in% c("main", "no_10q26", "no_CFH", "no_major3") & outcome %in% c("AS", "IS", "SVS", "MEGA_IS", "MEGA_SVS"),
               .(exposure, scheme, outcome, exclusion, nsnp, OR_raps = round(exp(b), 3), lo = round(exp(b - 1.96 * se), 3), hi = round(exp(b + 1.96 * se), 3), p_raps = signif(p, 2))],
            iv, by = c("exposure", "scheme", "outcome", "exclusion")), nrows = 500)
message("rev_mr_01b done")
