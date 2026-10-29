# Round-4 item C: 3p14 in Pan-UK Biobank (EUR) with UKB in-sample LD. SuSiE-RSS (L = 10) per trait, coloc.susie and coloc.abf
# for AMD-T2D, AMD-IS, T2D-IS, AMD-CAD, T2D-CAD; overlap of UKB credible sets with the FinnGen AMD/T2D 10-variant credible set.
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressMessages({library(susieR); library(coloc)})
O <- file.path(PROJ_ROOT, "results/revision4"); I <- file.path(O, "c_input")
OFF <- 14324L                                                   # GRCh38 - GRCh37 at 3p14
KEY <- c(rs4611812 = 64699445L, rs6795735 = 64705365L, rs6445425 = 64719091L, rs7428936 = 64725174L - OFF, rs62247658 = 64729479L - OFF,
         rs17676309 = NA_integer_)
fgcs <- read.delim(file.path(PROJ_ROOT, "results/revision/adamts9/credible_set_members.tsv"))
fg_amd <- fgcs[fgcs$set == "FG_AMD_CS1", ]
fg_pos37 <- as.integer(sapply(strsplit(fg_amd$variant, ":"), `[`, 2)) - OFF
tags <- sub("_z.tsv$", "", list.files(I, pattern = "_z.tsv$"))
d0 <- read.delim(file.path(I, paste0(tags[1], "_z.tsv"))); p <- nrow(d0)
R <- matrix(readBin(file.path(I, "R.f64"), "double", p * p, endian = "little"), p, p)
fits <- list(); dat <- list(); summ <- list(); vtab <- list()
for (t in tags) {
  d <- read.delim(file.path(I, paste0(t, "_z.tsv"))); stopifnot(nrow(d) == p, all(d$key == d0$key))
  n <- d$n[1]
  kr <- tryCatch(kriging_rss(d$z, R, n = n), error = function(e) NULL)
  nout <- if (is.null(kr)) NA else sum(kr$conditional_dist$logLR > 2 & abs(kr$conditional_dist$z) > 2)
  f <- susie_rss(z = d$z, R = R, n = n, L = 10, max_iter = 500, estimate_residual_variance = FALSE, coverage = 0.95)
  names(f$pip) <- d$key
  fits[[t]] <- f
  dat[[t]] <- list(beta = setNames(d$beta, d$key), varbeta = setNames(d$se^2, d$key), type = "cc", s = d$s[1], N = n,
                   MAF = setNames(pmin(d$af, 1 - d$af), d$key), snp = d$key, LD = `dimnames<-`(R, list(d$key, d$key)))
  cs <- f$sets$cs
  csum <- if (length(cs)) paste(sapply(seq_along(cs), function(k) { m <- cs[[k]]; top <- m[which.max(f$pip[m])]
      sprintf("%s:size=%d,lead=%s,lead_pip=%.3f,purity=%.2f,n_in_FinnGen_AMD_CS=%d", names(cs)[k], length(m), d$key[top], f$pip[top],
              f$sets$purity[k, "min.abs.corr"], sum(d$pos[m] %in% fg_pos37)) }), collapse = "; ") else "none"
  summ[[t]] <- data.frame(trait = t, n = n, n_cases = round(d$s[1] * n), n_variants = p, max_abs_z = max(abs(d$z)), key_max_z = d$key[which.max(abs(d$z))],
                          n_cs = length(cs), credible_sets = csum, kriging_outliers = nout)
  for (rs in names(KEY)) { if (is.na(KEY[rs])) next
    i <- which(d$pos == KEY[rs]); if (!length(i)) next
    csid <- NA; if (length(cs)) for (k in seq_along(cs)) if (i[1] %in% cs[[k]]) csid <- names(cs)[k]
    vtab[[length(vtab) + 1]] <- data.frame(trait = t, rsid = rs, key = d$key[i[1]], beta = d$beta[i[1]], se = d$se[i[1]], z = d$z[i[1]],
                                           pip = f$pip[i[1]], cs = csid) }
  # all FinnGen AMD credible-set variants
  for (q in which(d$pos %in% fg_pos37)) vtab[[length(vtab) + 1]] <- data.frame(trait = t, rsid = paste0("FinnGen_AMD_CS_member@", d$pos[q]),
      key = d$key[q], beta = d$beta[q], se = d$se[q], z = d$z[q], pip = f$pip[q], cs = NA)
  message(t, " CS=", length(cs), " ", csum)
}
pairs <- list(c("AMD_phe362.29", "T2D_phe250.2"), c("AMD_phe362.29", "T2D_E11"), c("AMD_phe362.29", "IS_I63"), c("T2D_phe250.2", "IS_I63"),
              c("AMD_phe362.29", "CAD_phe411.4"), c("T2D_phe250.2", "CAD_phe411.4"), c("IS_I63", "CAD_phe411.4"))
co <- list()
for (pr in pairs) {
  a <- pr[1]; b <- pr[2]
  for (p12 in c(1e-5, 5e-6, 1e-6)) {
    ab <- coloc.abf(dat[[a]][c("beta", "varbeta", "type", "s", "N", "MAF", "snp")], dat[[b]][c("beta", "varbeta", "type", "s", "N", "MAF", "snp")], p12 = p12)
    co[[length(co) + 1]] <- data.frame(trait1 = a, trait2 = b, method = "coloc.abf", p12 = p12, hit1 = NA, hit2 = NA,
      PP.H0 = ab$summary["PP.H0.abf"], PP.H1 = ab$summary["PP.H1.abf"], PP.H2 = ab$summary["PP.H2.abf"], PP.H3 = ab$summary["PP.H3.abf"], PP.H4 = ab$summary["PP.H4.abf"],
      lead_snp_H4 = ab$results$snp[which.max(ab$results$SNP.PP.H4)])
    if (length(fits[[a]]$sets$cs) && length(fits[[b]]$sets$cs)) {
      cs <- tryCatch(coloc.susie(fits[[a]], fits[[b]], p12 = p12), error = function(e) NULL)
      if (!is.null(cs) && !is.null(cs$summary)) for (j in seq_len(nrow(cs$summary))) { s <- cs$summary[j, ]
        co[[length(co) + 1]] <- data.frame(trait1 = a, trait2 = b, method = "coloc.susie", p12 = p12, hit1 = s$hit1, hit2 = s$hit2,
          PP.H0 = s$PP.H0.abf, PP.H1 = s$PP.H1.abf, PP.H2 = s$PP.H2.abf, PP.H3 = s$PP.H3.abf, PP.H4 = s$PP.H4.abf, lead_snp_H4 = NA) }
    }
  }
}
co <- do.call(rbind, co); rownames(co) <- NULL
write.table(co, file.path(O, "3p14_panukb_coloc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, vtab), file.path(O, "3p14_panukb_susie.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, summ), file.path(O, "3p14_panukb_susie_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(do.call(rbind, summ)); print(co[co$p12 == 1e-5, ]); print(do.call(rbind, vtab)[1:30, ])
