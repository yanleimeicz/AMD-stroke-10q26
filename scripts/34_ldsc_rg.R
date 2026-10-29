# Genome-wide genetic correlation (LDSC via GenomicSEM) between AMD and stroke subtypes / WMH.
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(GenomicSEM) })
DER <- file.path(PROJ_ROOT, "derived/gwas"); REF <- file.path(PROJ_ROOT, "data/raw/ref")
W <- file.path(PROJ_ROOT, "derived/ldsc"); OUT <- file.path(PROJ_ROOT, "results/ldsc")
dir.create(W, recursive = TRUE, showWarnings = FALSE); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(W)

# LD scores: original eur_w_ld_chr (GenomicSEM layout <ld>/<chr>.l2.ldscore.gz, <chr>.l2.M_5_50)
LD <- file.path(REF, "eur_w_ld_chr")
HM3 <- file.path(LD, "w_hm3.snplist")
if (!file.exists(HM3)) {   # rebuild HapMap3 SNP list: SNPs of the LD-score files + alleles from 1000G EUR .bim
  HM3 <- file.path(W, "w_hm3.snplist")
  if (!file.exists(HM3)) {
    snps <- rbindlist(lapply(1:22, function(c) fread(file.path(LD, paste0(c, ".l2.ldscore.gz")), select = "SNP")))
    bim <- fread(file.path(PROJ_ROOT, "data/raw/ref/g1000_eur/g1000_eur.bim"), header = FALSE, select = c(2, 5, 6), col.names = c("SNP", "A1", "A2"))
    hm <- merge(snps, bim[!duplicated(SNP)], by = "SNP")
    fwrite(hm, HM3, sep = "\t"); message("built HM3 list: ", nrow(hm), " SNPs")
  }
}

traits <- c("AMD_FG", "AMD_IAMDGC", "DRYAMD_FG", "WETAMD_FG", "AS", "IS", "LAS", "CES", "SVS", "WMH")
files <- c(); Ns <- c()
for (k in traits) {
  f <- file.path(W, paste0(k, ".txt"))
  d <- readRDS(file.path(DER, paste0(k, ".rds")))
  neff <- if ("ncase" %in% names(d)) 4 / (1 / d$ncase[1] + 1 / d$ncontrol[1]) else d$n[1]
  if (!file.exists(f)) fwrite(d[, .(SNP, A1 = ea, A2 = oa, effect = beta, P = p)], f, sep = "\t")
  files <- c(files, f); Ns <- c(Ns, neff); rm(d); gc()
}
if (!all(file.exists(paste0(traits, ".sumstats.gz"))))
  munge(files = files, hm3 = HM3, trait.names = traits, N = Ns, parallel = TRUE, cores = 10)
res <- ldsc(traits = paste0(traits, ".sumstats.gz"), sample.prev = rep(0.5, length(traits)) * c(rep(1, 9), NA),
            population.prev = c(0.03, 0.02, 0.02, 0.015, 0.03, 0.025, 0.005, 0.005, 0.005, NA), ld = LD, wld = LD, trait.names = traits)
S <- res$S; V <- res$V; k <- nrow(S)
sdg <- sqrt(diag(S)); R <- S / outer(sdg, sdg)
# delta-method SE for rg from the sampling covariance of S
idx <- which(lower.tri(S, diag = TRUE), arr.ind = TRUE); vidx <- function(i, j) which((idx[, 1] == max(i, j)) & (idx[, 2] == min(i, j)))
out <- list()
for (i in 1:k) for (j in 1:k) if (i < j) {
  a <- vidx(i, j); b <- vidx(i, i); c <- vidx(j, j)
  g <- c(1 / sqrt(S[i, i] * S[j, j]), -0.5 * S[i, j] / (S[i, i]^1.5 * sqrt(S[j, j])), -0.5 * S[i, j] / (S[j, j]^1.5 * sqrt(S[i, i])))
  vv <- V[c(a, b, c), c(a, b, c)]; se <- sqrt(as.numeric(t(g) %*% vv %*% g))
  out[[length(out) + 1]] <- data.table(trait1 = traits[i], trait2 = traits[j], rg = R[i, j], se = se, p = 2 * pnorm(-abs(R[i, j] / se)))
}
out <- rbindlist(out)
h2 <- data.table(trait = traits, h2_obs_or_liab = diag(S), se = sqrt(diag(V)[sapply(1:k, function(i) vidx(i, i))]), intercept = diag(res$I))
fwrite(out, file.path(OUT, "ldsc_rg_all.tsv"), sep = "\t"); fwrite(h2, file.path(OUT, "ldsc_h2.tsv"), sep = "\t")
print(h2); print(out[(trait1 %like% "AMD") != (trait2 %like% "AMD")][order(p)])
message("LDSC done")
