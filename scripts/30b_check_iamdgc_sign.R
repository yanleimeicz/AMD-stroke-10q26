# Validate IAMDGC z->beta conversion: sign concordance with FinnGen AMD at genome-wide significant SNPs.
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths())); suppressPackageStartupMessages(library(data.table))
a <- readRDS(file.path(PROJ_ROOT, "derived/gwas/AMD_IAMDGC.rds")); f <- readRDS(file.path(PROJ_ROOT, "derived/gwas/AMD_FG.rds"))
m <- merge(a[p < 5e-8], f[, .(SNP, ea2 = ea, oa2 = oa, b2 = beta, p2 = p)], by = "SNP")
comp <- c(A = "T", T = "A", C = "G", G = "C")
pal <- (m$ea == comp[m$oa])
same <- m$ea == m$ea2 | comp[m$ea] == m$ea2
m[, b2a := fifelse(same, b2, -b2)]
m <- m[!pal]
cat(sprintf("IAMDGC p<5e-8 SNPs in FinnGen: %d (non-palindromic); sign concordance: %.1f%%; among FinnGen p<5e-8: %.1f%%; cor(beta) = %.2f\n",
            nrow(m), 100 * mean(sign(m$beta) == sign(m$b2a)), 100 * m[p2 < 5e-8, mean(sign(beta) == sign(b2a))], cor(m$beta, m$b2a)))
print(m[order(p)][1:8, .(SNP, chr, ea, oa, beta = round(beta, 3), p, b_FG_aligned = round(b2a, 3), p2)])
