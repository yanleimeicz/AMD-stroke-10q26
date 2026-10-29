# 3p14 ADAMTS9 deep-dive: carotid intima-media thickness (cIMT, GCST007436, up to 71,128) and carotid plaque
# (GCST007435, 21,540 cases / 26,894 controls) from Franceschini et al. 2018 Nat Commun (PMID 30510157), European ancestry.
# GWAS Catalog harmonised files (GRCh38) are read for the 3p14 window only via remote tabix.
# Output -> derived/adamts9/gwas_cimt.rds (same columns as gwas_region.rds)
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressMessages({ library(data.table); library(Rsamtools) })
D <- file.path(PROJ_ROOT, "derived/adamts9")
S38 <- 63900000; E38 <- 65400000
BASE <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000"
STUDIES <- data.table(trait = c("CIMT", "PLAQUE"), gcst = c("GCST007436", "GCST007435"),
                      ncase = c(NA, 21540), ncontrol = c(NA, 26894), n = c(71128, 21540 + 26894))
rs <- fread(file.path(D, "region_snps.tsv"))
out <- rbindlist(lapply(seq_len(nrow(STUDIES)), function(i) {
  url <- sprintf("%s/%s/harmonised/%s.h.tsv.gz", BASE, STUDIES$gcst[i], STUDIES$gcst[i])
  hdr <- strsplit(readLines(gzcon(url(url)), n = 1), "\t")[[1]]
  for (attempt in 1:5) {
    x <- tryCatch(scanTabix(TabixFile(url), param = GRanges("3", IRanges(S38, E38)))[[1]], error = function(e) NULL)
    if (!is.null(x)) break
    Sys.sleep(20 * attempt)
  }
  d <- fread(text = paste(x, collapse = "\n"), header = FALSE, col.names = hdr, sep = "\t")
  message(STUDIES$trait[i], ": ", nrow(d), " rows; columns: ", paste(hdr, collapse = ","))
  # GWAS-SSF harmonised format (chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error,
  # effect_allele_frequency, p_value, rsid, n) or the legacy hm_* format
  pick <- function(...) { for (k in c(...)) if (k %in% hdr) return(d[[k]]); rep(NA, nrow(d)) }
  b <- as.numeric(pick("hm_beta", "beta"))
  orr <- as.numeric(pick("hm_odds_ratio", "odds_ratio")); b[is.na(b) & is.finite(orr)] <- log(orr[is.na(b) & is.finite(orr)])
  snp <- as.character(pick("hm_rsid", "rsid", "rs_id")); snp2 <- as.character(pick("rs_id"))
  snp[!grepl("^rs", snp) & grepl("^rs", snp2)] <- snp2[!grepl("^rs", snp) & grepl("^rs", snp2)]
  nn <- as.numeric(pick("n")); nn[is.na(nn)] <- STUDIES$n[i]
  r <- data.table(SNP = snp, chr = "3", pos = as.numeric(pick("hm_pos", "base_pair_location")),
                  ea = toupper(pick("hm_effect_allele", "effect_allele")), oa = toupper(pick("hm_other_allele", "other_allele")),
                  beta = b, se = as.numeric(pick("standard_error")), p = as.numeric(pick("p_value")),
                  eaf = as.numeric(pick("hm_effect_allele_frequency", "effect_allele_frequency")),
                  ncase = STUDIES$ncase[i], ncontrol = STUDIES$ncontrol[i], n = nn, trait = STUDIES$trait[i])
  r[grepl("^rs", SNP) & is.finite(beta) & is.finite(se) & se > 0]
}))
out <- merge(out, rs[, .(SNP, pos38, pos37)], by = "SNP", all.x = TRUE)
out[is.na(pos38), pos38 := pos]; out[is.na(pos37), pos37 := pos - 14324]
saveRDS(out, file.path(D, "gwas_cimt.rds"))
print(out[, .(n_snps = .N, min_p = min(p), top = SNP[which.min(p)]), by = trait])
print(out[SNP %in% c("rs4611812", "rs17676309", "rs6795735", "rs6445425")])
message("cimt done")
