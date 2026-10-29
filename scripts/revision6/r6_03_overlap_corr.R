# Round 6, item A(2): genome-wide null-SNP z-statistic correlations between stroke GWAS datasets (sample-overlap check)
# HapMap3 SNPs (LDSC w_hm3), aligned to the w_hm3 A1 allele; SNPs with |z| < 1.96 in both datasets; thinned to 1 SNP / 250 kb
# to reduce LD dependence; jackknife SE over 22 chromosomes. Also constructs the per-SNP GIGASTROKE remainder
# (GIGA - MEGA - f * FinnGen) and checks its correlation with MEGASTROKE, FinnGen and UK Biobank.
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
suppressMessages({ library(data.table) })
P <- PROJ_ROOT; O <- file.path(P, "results/revision6"); dir.create(O, showWarnings = FALSE)
hm3 <- fread(file.path(P, "derived/ldsc/w_hm3.snplist"))
bim <- rbindlist(lapply(1:22, function(c) fread(file.path(P, sprintf("derived/lava/ref/g1000_eur_chr%d.bim", c)),
                         col.names = c("chr", "SNP", "cm", "pos37", "b1", "b2"))[SNP %in% hm3$SNP]))
bim <- bim[!duplicated(SNP)]
hm <- merge(hm3, bim[, .(SNP, chr = as.integer(chr), pos37)], by = "SNP")
pal <- function(a, b) (a == "A" & b == "T") | (a == "T" & b == "A") | (a == "C" & b == "G") | (a == "G" & b == "C")
hm <- hm[!pal(A1, A2)]
cat("HapMap3 SNPs with hg19 positions:", nrow(hm), "\n")
align <- function(d) {               # d: SNP, ea, oa, beta, se
  d <- merge(hm[, .(SNP, A1, A2, chr, pos37)], d, by = "SNP")
  s <- ifelse(d$ea == d$A1 & d$oa == d$A2, 1, ifelse(d$ea == d$A2 & d$oa == d$A1, -1, NA))
  d[, `:=`(beta = s * beta, z = s * beta / se)][!is.na(s) & is.finite(z), .(SNP, chr, pos37, beta, se, z)]
}
rd <- function(f) { x <- readRDS(f); setDT(x); x[, .(SNP, ea = toupper(ea), oa = toupper(oa), beta, se)] }
D <- list()
D$GIGA_IS <- align(rd(file.path(P, "derived/gwas/IS.rds")))
D$GIGA_AS <- align(rd(file.path(P, "derived/gwas/AS.rds")))
D$MEGA_IS <- align(rd(file.path(P, "derived/revision/MEGA_IS.rds")))
D$MEGA_AS <- align(rd(file.path(P, "derived/revision/MEGA_AS.rds")))
gz <- function(f) {                  # decompress (gz or bgz) to project tmp, then fread
  out <- file.path(P, "tmp", paste0(basename(f), ".txt"))
  if (!file.exists(out)) { con <- gzfile(f, "rb"); oc <- file(out, "wb")
    repeat { b <- readBin(con, "raw", 1e8); if (!length(b)) break; writeBin(b, oc) }; close(con); close(oc) }
  out }
fg <- fread(gz(file.path(P, "derived/revision6/finngen_R12_I9_STR.gz")), select = c("rsids", "ref", "alt", "beta", "sebeta"), showProgress = FALSE)
fg <- fg[!grepl(",", rsids)]
D$FinnGen_I9_STR <- align(fg[, .(SNP = rsids, ea = alt, oa = ref, beta, se = sebeta)])
uk <- fread(gz(file.path(P, "derived/revision6/panukb_icd10_I63.tsv.bgz")), select = c("chr", "pos", "ref", "alt", "beta_EUR", "se_EUR", "low_confidence_EUR"),
            showProgress = FALSE)
uk <- uk[low_confidence_EUR %in% c(FALSE, "false") & is.finite(beta_EUR)]
uk[, chr := suppressWarnings(as.integer(chr))]
uk <- merge(uk, hm[, .(SNP, chr, pos37)], by.x = c("chr", "pos"), by.y = c("chr", "pos37"))
D$UKB_I63 <- align(uk[, .(SNP, ea = alt, oa = ref, beta = beta_EUR, se = se_EUR)])
cat(sapply(D, nrow), "\n")

# per-SNP GIGASTROKE remainder (IS, AS): GIGA - MEGA - f * FinnGen (UKB inside MEGASTROKE)
mk_rem <- function(g, m, fgd, f = 0.15) {
  x <- merge(merge(g[, .(SNP, chr, pos37, bG = beta, sG = se)], m[, .(SNP, bM = beta, sM = se)], by = "SNP"), fgd[, .(SNP, bF = beta, sF = se)], by = "SNP")
  x[, `:=`(IG = 1 / sG^2, IM = 1 / sM^2)]; x[, IF := f * IG]; x[, IR := IG - IM - IF]
  x <- x[IR > 0.02 * IG]
  x[, bR := (IG * bG - IM * bM - IF * bF) / IR]
  x[, sR := sqrt(1 / IR + (IF / IR)^2 * pmax(0, 1 / IF - sF^2))]   # proxy error: Var(b_Fcomp - b_F) = 1/IF - 1/I_F(R12)
  x[, .(SNP, chr, pos37, beta = bR, se = sR, z = bR / sR)]
}
D$GIGA_remainder_IS <- mk_rem(D$GIGA_IS, D$MEGA_IS, D$FinnGen_I9_STR)
D$GIGA_remainder_AS <- mk_rem(D$GIGA_AS, D$MEGA_AS, D$FinnGen_I9_STR)
# remainder v2: additionally subtract the UK Biobank component (2,932 AS / 2,211 IS incident cases, 369,419 controls;
# not part of the public MEGASTROKE EUR summary statistics), using Pan-UKB I63 as proxy (round-4 specification)
ne <- function(a, b) 4 / (1 / a + 1 / b)
fU <- c(IS = ne(2211, 369419) / ne(62100, 1234808), AS = ne(2932, 369419) / ne(73652, 1234808))
mk_rem2 <- function(g, m, fgd, ukd, fu, f = 0.15) {
  x <- merge(merge(merge(g[, .(SNP, chr, pos37, bG = beta, sG = se)], m[, .(SNP, bM = beta, sM = se)], by = "SNP"),
                   fgd[, .(SNP, bF = beta, sF = se)], by = "SNP"), ukd[, .(SNP, bU = beta, sU = se)], by = "SNP")
  x[, `:=`(IG = 1 / sG^2, IM = 1 / sM^2)]; x[, `:=`(IF = f * IG, IU = fu * IG)]; x[, IR := IG - IM - IF - IU]
  x <- x[IR > 0.02 * IG]
  x[, bR := (IG * bG - IM * bM - IF * bF - IU * bU) / IR]
  x[, sR := sqrt(1 / IR + (IF / IR)^2 * pmax(0, 1 / IF - sF^2) + (IU / IR)^2 * pmax(0, 1 / IU - sU^2))]
  x[, .(SNP, chr, pos37, beta = bR, se = sR, z = bR / sR)]
}
D$GIGA_remainder2_IS <- mk_rem2(D$GIGA_IS, D$MEGA_IS, D$FinnGen_I9_STR, D$UKB_I63, fU["IS"])
D$GIGA_remainder2_AS <- mk_rem2(D$GIGA_AS, D$MEGA_AS, D$FinnGen_I9_STR, D$UKB_I63, fU["AS"])

thin <- function(x) { x <- x[order(chr, pos37)]; x[, bin := floor(pos37 / 250000)]; x[!duplicated(paste(chr, bin))] }
pair <- function(a, b) {
  x <- merge(D[[a]][, .(SNP, chr, pos37, za = z)], D[[b]][, .(SNP, zb = z)], by = "SNP")
  x <- x[abs(za) < 1.96 & abs(zb) < 1.96]
  set.seed(1); x <- x[sample(.N)]; x <- thin(x)
  r <- cor(x$za, x$zb)
  jk <- sapply(1:22, function(c) { y <- x[chr != c]; cor(y$za, y$zb) })
  se <- sqrt((22 - 1) / 22 * sum((jk - mean(jk))^2))
  data.table(a = a, b = b, n_snp = nrow(x), r_null_z = r, se_jackknife = se)
}
prs <- list(c("GIGA_IS", "MEGA_IS"), c("GIGA_IS", "FinnGen_I9_STR"), c("GIGA_IS", "UKB_I63"), c("MEGA_IS", "FinnGen_I9_STR"),
            c("MEGA_IS", "UKB_I63"), c("FinnGen_I9_STR", "UKB_I63"), c("GIGA_remainder_IS", "MEGA_IS"), c("GIGA_remainder_IS", "FinnGen_I9_STR"),
            c("GIGA_remainder_IS", "UKB_I63"), c("GIGA_AS", "MEGA_AS"), c("GIGA_AS", "FinnGen_I9_STR"), c("MEGA_AS", "FinnGen_I9_STR"),
            c("MEGA_AS", "UKB_I63"), c("GIGA_remainder_AS", "MEGA_AS"), c("GIGA_remainder_AS", "FinnGen_I9_STR"), c("GIGA_remainder_AS", "UKB_I63"), c("GIGA_remainder2_IS", "MEGA_IS"), c("GIGA_remainder2_IS", "FinnGen_I9_STR"), c("GIGA_remainder2_IS", "UKB_I63"), c("GIGA_remainder2_AS", "MEGA_AS"), c("GIGA_remainder2_AS", "FinnGen_I9_STR"), c("GIGA_remainder2_AS", "UKB_I63"))
res <- rbindlist(lapply(prs, function(p) pair(p[1], p[2])))
# analytic expectations under nesting: corr(GIGA, MEGA) = sqrt(Neff_M / Neff_G); corr(GIGA, FinnGen comp) ~ sqrt(share)
ne <- function(a, b) 4 / (1 / a + 1 / b)
exp_tab <- data.table(a = c("GIGA_IS", "GIGA_AS", "GIGA_IS"), b = c("MEGA_IS", "MEGA_AS", "FinnGen_I9_STR"),
                      expected = c(sqrt(ne(34217, 406111) / ne(62100, 1234808)), sqrt(ne(40585, 406111) / ne(73652, 1234808)), NA))
res <- merge(res, exp_tab, by = c("a", "b"), all.x = TRUE)
fwrite(res, file.path(O, "null_z_correlations.tsv"), sep = "\t")
print(res)
