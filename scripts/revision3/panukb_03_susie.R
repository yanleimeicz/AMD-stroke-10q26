# Round-3: SuSiE-RSS on Pan-UKB EUR 10q26 with UKB in-sample (British) LD; kriging_rss LD-mismatch diagnostics.
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressMessages(library(susieR))
O <- file.path(PROJ_ROOT, "results/revision3/panukb"); I <- file.path(O, "susie_input")
POS <- c(rs61871745 = 124210369, rs72631113 = 124213449, rs10490924 = 124214448, rs36212733 = 124215211,
         rs3750848 = 124215315, rs11200638 = 124220544, rs2672592 = 124230750)
out <- list(); diag <- list(); csm <- list()
for (f in list.files(I, pattern = "_z.tsv$")) {
  t <- sub("_z.tsv$", "", f)
  d <- read.delim(file.path(I, f)); p <- nrow(d)
  R <- matrix(readBin(file.path(I, paste0(t, "_R.f64")), "double", p * p, endian = "little"), p, p)
  n <- d$n[1]
  kr <- tryCatch(kriging_rss(d$z, R, n = n), error = function(e) NULL)
  nout <- if (is.null(kr)) NA else sum(kr$conditional_dist$logLR > 2 & abs(kr$conditional_dist$z) > 2)
  fit <- tryCatch(susie_rss(z = d$z, R = R, n = n, L = 10, max_iter = 500, estimate_residual_variance = FALSE),
                  error = function(e) { message(t, ": ", conditionMessage(e)); NULL })
  if (is.null(fit)) next
  cs <- fit$sets$cs; pip <- fit$pip
  csid <- rep(NA_integer_, p)
  if (length(cs)) for (k in seq_along(cs)) csid[cs[[k]]] <- as.integer(sub("L", "", names(cs)[k]))
  maxz <- which.max(abs(d$z))
  if (length(cs)) for (k in seq_along(cs)) csm[[length(csm) + 1]] <- data.frame(trait = t, cs = names(cs)[k], pos = d$pos[cs[[k]]], ref = d$ref[cs[[k]]],
      alt = d$alt[cs[[k]]], z = d$z[cs[[k]]], pip = pip[cs[[k]]])
  for (rs in names(POS)) {
    i <- which(d$pos == POS[rs])
    out[[length(out) + 1]] <- data.frame(trait = t, rsid = rs, pos_grch37 = POS[rs], in_input = length(i) > 0,
      z = if (length(i)) d$z[i[1]] else NA, pip = if (length(i)) pip[i[1]] else NA, cs = if (length(i)) csid[i[1]] else NA)
  }
  csum <- if (length(cs)) paste(sapply(seq_along(cs), function(k) {
      m <- cs[[k]]; top <- m[which.max(pip[m])]
      sprintf("%s:size=%d,lead_pos=%d,lead_pip=%.3f,purity_min_r=%.2f", names(cs)[k], length(m), d$pos[top], pip[top], fit$sets$purity[k, "min.abs.corr"])
    }), collapse = "; ") else "none"
  diag[[length(diag) + 1]] <- data.frame(trait = t, n_variants = p, n = n, max_abs_z = max(abs(d$z)), pos_max_z = d$pos[maxz],
      n_credible_sets = length(cs), credible_sets = csum, kriging_n_outliers_logLR2_absz2 = nout)
  message(t, " CS=", length(cs), " ", csum)
}
write.table(do.call(rbind, out), file.path(O, "panukb_susie.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, csm), file.path(O, "panukb_susie_cs_members.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, diag), file.path(O, "panukb_susie_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
