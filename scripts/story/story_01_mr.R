# Story revision: positive controls, per-locus effects, pathway-partitioned MR, AMD -> CAD MR.
# Outputs -> results/story/: positive_controls.tsv, per_locus_effects.tsv, pathway_mr.tsv, amd_cad_mr.tsv
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(Rsamtools) })
setDTthreads(30); set.seed(20261007)
DER <- file.path(PROJ_ROOT, "derived/gwas"); DRV <- file.path(PROJ_ROOT, "derived/revision"); SD <- file.path(PROJ_ROOT, "derived/story")
OUT <- file.path(PROJ_ROOT, "results/story"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE); dir.create(SD, showWarnings = FALSE, recursive = TRUE)
PLINK <- file.path(PROJ_ROOT, "data/raw/tools/plink/plink.exe"); BFILE <- file.path(PROJ_ROOT, "data/raw/ref/g1000_eur/g1000_eur")
TMP <- file.path(PROJ_ROOT, "tmp/story"); dir.create(TMP, showWarnings = FALSE, recursive = TRUE)
path_of <- function(k) { f <- file.path(DRV, paste0(k, ".rds")); if (file.exists(f)) f else file.path(DER, paste0(k, ".rds")) }

# ------------------------------------------------------------------ instruments (A_dist) with eaf from exposure files
inst <- fread(file.path(PROJ_ROOT, "results/revision/mr/instruments_all_schemes.tsv"))[scheme == "A_dist" & exposure %in% c("AMD_FG", "AMD_IAMDGC")]
EXF <- c(AMD_FG = "AMD_FG", AMD_IAMDGC = "AMD_IAMDGC_fix")
for (e in names(EXF)) {
  d <- readRDS(path_of(EXF[[e]]))
  if (!"eaf" %in% names(d)) d[, eaf := NA_real_]
  inst[exposure == e, eafx := d$eaf[match(SNP, d$SNP)]]
}
inst[, chr := as.character(chr)]

# gene / pathway annotation (GRCh38 positions; nearest established AMD locus gene, Fritsche 2016)
ann <- fread(text = "SNP\tgene\tpathway
rs10801195\t1q31.2 (193.4 Mb)\tother
rs12354144\tCFH region (KCNT2)\tcomplement
rs1329424\tCFH\tcomplement
rs72744971\tCFH/CFHR region\tcomplement
rs182465192\t1q32.1 (199.0 Mb)\tother
rs12566266\t1q32.1 (200.0 Mb)\tother
rs13411964\t2q31.2 (179.9 Mb)\tother
rs6445425\tADAMTS9-AS2\tother
rs116675086\tCOL8A1\tother
rs3775221\t4q22.1 (87.1 Mb)\tother
rs542286634\t4q24 (103.6 Mb)\tother
rs1436556086\t4q24 (107.2 Mb)\tother
rs11097986\t4q25 (108.4 Mb)\tother
rs139779213\tCFI\tcomplement
rs189631141\t5p13.2 (35.6 Mb)\tother
rs429608\tC2/CFB\tcomplement
rs7758685\tVEGFA\tother
rs3812112\tCOL10A1\tother
rs10808141\tKMT2E/SRPK2\tother
rs12413042\t10q26 (121.4 Mb)\t10q26
rs61871744\tARMS2-HTRA1\t10q26
rs567542414\t10q26 (123.8 Mb)\t10q26
rs3138142\tRDH5/CD63\tother
rs8000179\tB3GALTL\tother
rs1956526\tRAD51B\tother
rs759283\t14q24.2 (72.9 Mb)\tother
rs2043085\tLIPC\tlipid
rs74026308\t15q24.1 (73.9 Mb)\tother
rs3825991\tABHD2\tother
rs1800775\tCETP\tlipid
rs2230199\tC3\tcomplement
rs429358\tAPOE\tlipid
rs5749492\tTIMP3/SYN3\tother
rs6665824\tCFH\tcomplement
rs114920504\tCFH/CFHR region\tcomplement
rs11884770\tCOL4A3\tother
rs7428936\tADAMTS9-AS2\tother
rs140647181\tCOL8A1\tother
rs11732475\tCFI\tcomplement
rs3131857\tMHC class I (29.6 Mb)\tother
rs3131064\tMHC class I (30.8 Mb)\tother
rs635348\tC2/CFB\tcomplement
rs943080\tVEGFA\tother
rs7803454\tPILRB/PILRA\tother
rs1142\tKMT2E/SRPK2\tother
rs13278062\tTNFRSF10A\tother
rs12357257\tARHGAP21\tother
rs3750846\tARMS2-HTRA1\t10q26
rs3138141\tRDH5/CD63\tother
rs61941274\tACAD10\tother
rs61985136\tRAD51B\tother
rs11080055\tTMEM97/VTN\tcomplement
rs6565597\tNPLOC4/TSPAN10\tother
rs10422209\tCNN2\tother
rs11569415\tC3\tcomplement
rs8135665\tSLC16A8\tother")
inst <- merge(inst, ann, by = "SNP", all.x = TRUE)
if (anyNA(inst$pathway)) { message("unannotated: ", paste(inst[is.na(pathway), SNP], collapse = ",")); inst[is.na(pathway), `:=`(gene = "unannotated", pathway = "other")] }

# ------------------------------------------------------------------ proxies (1000G EUR, r2 > 0.8, +/-500 kb, phased)
U <- unique(inst$SNP); writeLines(U, file.path(TMP, "inst_snps.txt"))
system2(PLINK, c("--bfile", BFILE, "--r2", "in-phase", "with-freqs", "--ld-snp-list", file.path(TMP, "inst_snps.txt"),
                 "--ld-window-kb", "500", "--ld-window", "99999", "--ld-window-r2", "0.8", "--threads", "30", "--out", file.path(TMP, "proxies")), stdout = FALSE, stderr = FALSE)
PX <- fread(file.path(TMP, "proxies.ld"))[SNP_A != SNP_B]
PX[, c("hapA", "hapB") := tstrsplit(PHASE, "/", fixed = TRUE)]
PX[, `:=`(a1 = substr(hapA, 1, 1), b1 = substr(hapA, 2, 2), a2 = substr(hapB, 1, 1), b2 = substr(hapB, 2, 2))]
PX <- PX[, .(target = SNP_A, proxy = SNP_B, r2 = R2, a1, b1, a2, b2)][order(target, -r2)]

# ------------------------------------------------------------------ CAD (Aragam 2022): remote tabix point queries (instruments + proxies) with the local index
CADF <- file.path(SD, "CAD_Aragam2022_points.tsv.gz")
CADURL <- "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132314/harmonised/GCST90132314.h.tsv.gz"
CADIDX <- file.path(PROJ_ROOT, "tmp/rev_adamts9/GCST90132314.h.tsv.gz.tbi")
CADHDR <- c("chr", "pos", "ea", "oa", "beta", "se", "eaf", "p", "variant_id", "markername", "freqse", "minfreq", "maxfreq",
            "direction", "hetisq", "hetchisq", "hetdf", "hetpval", "cases", "effective_cases", "n", "meta", "conv", "or", "hm_code", "SNP")
tabix_fetch <- function(url, gr, hdr) {
  for (attempt in 1:6) {
    res <- tryCatch({
      tf <- TabixFile(url, index = CADIDX); x <- unlist(scanTabix(tf, param = gr))
      if (!length(x)) return(data.table())
      d <- fread(text = paste(x, collapse = "\n"), header = FALSE, sep = "\t", na.strings = c("NA", ""))
      setnames(d, hdr[seq_len(ncol(d))]); d }, error = function(e) { message("tabix error: ", conditionMessage(e)); NULL })
    if (!is.null(res)) return(res)
    Sys.sleep(15 * attempt)
  }
  stop("tabix failed: ", url)
}
pos38map <- unique(rbind(readRDS(path_of("AS"))[, .(SNP, chr = as.character(chr), pos)], readRDS(path_of("AMD_FG"))[, .(SNP, chr = as.character(chr), pos)]), by = "SNP")
CADLOCAL <- file.path(SD, "GCST90132314.h.tsv.gz")    # full file (parallel ranged download)
cad_all <- NULL
fetch_points <- function(snps, outf) {
  if (file.exists(outf)) return(fread(outf))
  if (is.null(cad_all)) {
    x <- fread(CADLOCAL, select = c("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "p_value", "rsid"))
    setnames(x, c("chr", "pos", "ea", "oa", "beta", "se", "eaf", "p", "SNP"))
    x[, SNP := sub(",.*", "", as.character(SNP))]
    cad_all <<- unique(x[grepl("^rs", SNP) & !is.na(beta) & !is.na(se) & se > 0], by = "SNP")
    message("CAD full file rows: ", nrow(cad_all))
  }
  x <- cad_all[SNP %in% snps, .(SNP, chr = as.character(chr), pos, ea = toupper(ea), oa = toupper(oa), beta, se, p, eaf)]
  fwrite(x, outf, sep = "\t"); x
}
cad <- fetch_points(unique(c(U, PX$proxy)), CADF); cad[, chr := as.character(chr)]
message("CAD rows: ", nrow(cad))

get_outcome <- function(o) {
  d <- if (o == "CAD") cad else readRDS(path_of(o))
  d <- as.data.table(d)[SNP %in% c(U, PX$proxy)]
  if (!"eaf" %in% names(d)) d[, eaf := NA_real_]
  have <- d[SNP %in% U]; miss <- setdiff(U, have$SNP); sub <- list()
  for (t in miss) {
    cand <- PX[target == t & proxy %in% d$SNP]; if (!nrow(cand)) next
    pr <- cand[1]; r <- d[SNP == pr$proxy][1]
    if (r$ea == pr$b1 && r$oa == pr$b2) { tea <- pr$a1; toa <- pr$a2 } else if (r$ea == pr$b2 && r$oa == pr$b1) { tea <- pr$a2; toa <- pr$a1 } else next
    sub[[t]] <- data.table(SNP = t, ea = tea, oa = toa, beta = r$beta, se = r$se, p = r$p, eaf = NA_real_, proxy = pr$proxy, proxy_r2 = pr$r2)
  }
  have[, `:=`(proxy = NA_character_, proxy_r2 = NA_real_)]
  rbind(have[, .(SNP, ea, oa, beta, se, p, eaf, proxy, proxy_r2)], rbindlist(sub), fill = TRUE)
}
OUTC <- c("AS", "IS", "SVS", "LAS", "WMH", "CAD", "AMD_FG", "AMD_IAMDGC_fix")
OD <- setNames(lapply(OUTC, get_outcome), OUTC)
for (o in OUTC) message(sprintf("%s: %d instr rows (%d proxied)", o, nrow(OD[[o]]), sum(!is.na(OD[[o]]$proxy))))

comp <- c(A = "T", T = "A", C = "G", G = "C")
harm <- function(ex, oc) {
  m <- merge(ex[, .(SNP, chr, pos, ea, oa, bx = beta, sx = se, eafx, gene, pathway, locus)],
             oc[, .(SNP, eo = ea, oo = oa, by = beta, sy = se, proxy)], by = "SNP")
  m[, flip := fifelse(eo == ea & oo == oa, 1, fifelse(eo == oa & oo == ea, -1,
              fifelse(comp[eo] == ea & comp[oo] == oa, 1, fifelse(comp[eo] == oa & comp[oo] == ea, -1, NA_real_))))]
  m[, pal := paste0(ea, oa) %in% c("AT", "TA", "CG", "GC")]
  m <- m[!is.na(flip) & !(pal & (is.na(eafx) | (eafx > 0.42 & eafx < 0.58)))]
  m[, by := by * flip]
  m[bx < 0, `:=`(bx = -bx, by = -by)]           # align to the AMD risk allele
  m[]
}
ivw_mre <- function(m) {
  if (nrow(m) == 1) return(c(b = m$by / m$bx, se = m$sy / abs(m$bx), p = 2 * pnorm(-abs(m$by / m$sy)), n = 1))
  w <- m$bx^2 / m$sy^2; b <- sum(w * m$by / m$bx) / sum(w); se <- 1 / sqrt(sum(w))
  if (nrow(m) >= 3) { phi <- sum(w * (m$by / m$bx - b)^2) / (nrow(m) - 1); se <- se * sqrt(max(1, phi)) }
  c(b = b, se = se, p = 2 * pnorm(-abs(b / se)), n = nrow(m))
}

# ------------------------------------------------------------------ 1a. positive control: AMD instruments in the other AMD dataset
pc <- list()
for (pr in list(c("AMD_FG", "AMD_IAMDGC_fix"), c("AMD_IAMDGC", "AMD_FG"))) {
  m <- harm(inst[exposure == pr[1]], OD[[pr[2]]]); r <- ivw_mre(m)
  pc[[length(pc) + 1]] <- data.table(control = "AMD instruments in an independent AMD GWAS", exposure = pr[1], outcome = sub("_fix", "", pr[2]),
                                      unit = "per unit log-odds AMD liability", nsnp = r[["n"]], b = r[["b"]], se = r[["se"]], p = r[["p"]],
                                      OR = exp(r[["b"]]), lo = exp(r[["b"]] - 1.96 * r[["se"]]), hi = exp(r[["b"]] + 1.96 * r[["se"]]),
                                      concordant_sign = mean(sign(m$by) == 1))
}
# ------------------------------------------------------------------ 1b. positive control: SBP -> stroke (existing) and SBP -> CAD (new)
sbp_st <- fread(file.path(PROJ_ROOT, "results/revision/adamts9/mr_cardiometabolic_to_stroke.tsv"))[exposure == "SBP"]
for (i in seq_len(nrow(sbp_st))) pc[[length(pc) + 1]] <- data.table(control = "Systolic blood pressure (Evangelou 2018) -> stroke (from revision/adamts9)",
  exposure = "SBP", outcome = sbp_st$outcome[i], unit = "per 1 mmHg", nsnp = sbp_st$nsnp[i], b = sbp_st$theta[i], se = sbp_st$se[i], p = sbp_st$p[i],
  OR = exp(sbp_st$theta[i]), lo = exp(sbp_st$theta[i] - 1.96 * sbp_st$se[i]), hi = exp(sbp_st$theta[i] + 1.96 * sbp_st$se[i]), concordant_sign = NA_real_)
sbp_iv <- fread(file.path(PROJ_ROOT, "tmp/rev_adamts9/clump_SBP.clumped"))$SNP
bim <- fread(paste0(BFILE, ".bim"), header = FALSE, select = c(1, 2, 4, 5, 6), col.names = c("chr", "SNP", "pos37", "a1", "a2")); bim[, chr := as.character(chr)]
sbp <- fread(file.path(PROJ_ROOT, "data/raw/gwas/cardiometab/SBP_Evangelou2018.txt.gz"), select = c("MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P"))
setnames(sbp, c("marker", "ea", "oa", "eaf", "beta", "se", "p")); sbp[, c("chr", "pos37") := tstrsplit(marker, ":", keep = 1:2)][, pos37 := as.integer(pos37)]
sbp[, `:=`(ea = toupper(ea), oa = toupper(oa))]
sbp <- merge(sbp, bim[SNP %in% sbp_iv, .(chr, pos37, SNP, a1, a2)], by = c("chr", "pos37"))
sbp <- unique(sbp[(ea == a1 & oa == a2) | (ea == a2 & oa == a1), .(SNP, ea, oa, beta, se, p, eaf)], by = "SNP")
SBPCAD <- file.path(SD, "CAD_Aragam2022_SBPinstr.tsv.gz")
invisible(fetch_points(sbp$SNP, SBPCAD))
cs <- fread(SBPCAD)
m <- merge(sbp[, .(SNP, ea, oa, bx = beta, sx = se, eafx = eaf)], cs[, .(SNP, eo = ea, oo = oa, by = beta, sy = se)], by = "SNP")
m[, by := fifelse(eo == ea & oo == oa, by, fifelse(eo == oa & oo == ea, -by, NA_real_))]
m <- m[!is.na(by) & !(paste0(ea, oa) %in% c("AT", "TA", "CG", "GC") & eafx > 0.42 & eafx < 0.58)]
r <- ivw_mre(m)
pc[[length(pc) + 1]] <- data.table(control = "Systolic blood pressure (Evangelou 2018) -> CAD (Aragam 2022)", exposure = "SBP", outcome = "CAD",
  unit = "per 1 mmHg", nsnp = r[["n"]], b = r[["b"]], se = r[["se"]], p = r[["p"]], OR = exp(r[["b"]]), lo = exp(r[["b"]] - 1.96 * r[["se"]]),
  hi = exp(r[["b"]] + 1.96 * r[["se"]]), concordant_sign = NA_real_)
pcd <- rbindlist(pc, fill = TRUE); fwrite(pcd, file.path(OUT, "positive_controls.tsv"), sep = "\t"); print(pcd)

# ------------------------------------------------------------------ 2. per-locus Wald ratios (aligned to AMD risk allele)
PO <- c("AS", "IS", "SVS", "WMH", "CAD")
pl <- rbindlist(lapply(c("AMD_FG", "AMD_IAMDGC"), function(e) rbindlist(lapply(PO, function(o) {
  m <- harm(inst[exposure == e], OD[[o]])
  m[, .(exposure = e, locus = fifelse(locus == "other", gene, locus), gene, pathway, SNP, chr, pos, beta_amd = bx, se_amd = sx,
        outcome = o, beta_out = by, se_out = sy, p_out = 2 * pnorm(-abs(by / sy)), wald = by / bx, se_wald = sy / bx,
        proxy = proxy)]
}))))
fwrite(pl, file.path(OUT, "per_locus_effects.tsv"), sep = "\t")

# ------------------------------------------------------------------ 3. pathway-partitioned MR + test against causal expectation
GROUPS <- list(all = c("complement", "10q26", "lipid", "other"), no_10q26 = c("complement", "lipid", "other"),
               complement = "complement", `10q26` = "10q26", lipid = "lipid", other = "other")
ORobs <- 2.08
pw <- rbindlist(lapply(c("AMD_FG", "AMD_IAMDGC"), function(e) rbindlist(lapply(c("AS", "IS", "SVS", "LAS", "WMH", "CAD"), function(o) {
  m0 <- harm(inst[exposure == e], OD[[o]])
  rbindlist(lapply(names(GROUPS), function(g) {
    m <- m0[pathway %in% GROUPS[[g]]]; if (!nrow(m)) return(NULL)
    r <- ivw_mre(m); bin <- o != "WMH"
    out <- data.table(exposure = e, outcome = o, group = g, nsnp = r[["n"]], method = ifelse(r[["n"]] == 1, "Wald ratio", "IVW (MRE)"),
                      b = r[["b"]], se = r[["se"]], p = r[["p"]],
                      OR = ifelse(bin, exp(r[["b"]]), NA), lo = ifelse(bin, exp(r[["b"]] - 1.96 * r[["se"]]), NA), hi = ifelse(bin, exp(r[["b"]] + 1.96 * r[["se"]]), NA),
                      mde80_logOR = (qnorm(0.975) + qnorm(0.8)) * r[["se"]])
    out[, mde80_OR := ifelse(bin, exp(mde80_logOR), NA)]
    if (o %in% c("AS", "IS", "SVS", "LAS")) for (pp in c(0.05, 0.063)) {
      bexp <- log(ORobs) * pp * (1 - pp); z <- (r[["b"]] - bexp) / r[["se"]]
      out[, paste0(c("expected_OR_p", "p_below_expected_p"), pp) := list(exp(bexp), pnorm(z))]
    }
    out
  }), fill = TRUE)
}), fill = TRUE)), fill = TRUE)
fwrite(pw, file.path(OUT, "pathway_mr.tsv"), sep = "\t")
print(pw[exposure == "AMD_FG", .(outcome, group, nsnp, OR = round(OR, 3), lo = round(lo, 3), hi = round(hi, 3), b = round(b, 4), p = signif(p, 2))])

# ------------------------------------------------------------------ 6. AMD -> CAD (A_dist IVW): all, no_10q26, no_major3, no_all
LIPID <- c("CETP", "LIPC", "ABCA1", "LPL"); MAJOR3 <- c("CFH", "10q26", "C3")
EXCL <- list(main = character(0), no_10q26 = "10q26", no_major3 = MAJOR3, no_all = c(MAJOR3, "APOE", "ADAMTS9", LIPID))
cadmr <- rbindlist(lapply(c("AMD_FG", "AMD_IAMDGC"), function(e) {
  m0 <- harm(inst[exposure == e], OD[["CAD"]])
  rbindlist(lapply(names(EXCL), function(x) {
    m <- m0[!locus %in% EXCL[[x]]]; r <- ivw_mre(m)
    # weighted median (simple implementation, bootstrap SE)
    wm <- function(bx, by, sy) { b <- by / bx; w <- (bx / sy)^2; o <- order(b); b <- b[o]; w <- w[o] / sum(w); cw <- cumsum(w) - 0.5 * w
                                  i <- max(which(cw < 0.5)); b[i] + (b[i + 1] - b[i]) * (0.5 - cw[i]) / (cw[i + 1] - cw[i]) }
    wmb <- if (nrow(m) >= 3) wm(m$bx, m$by, m$sy) else NA_real_
    wmse <- if (nrow(m) >= 3) sd(replicate(1000, { bxs <- rnorm(nrow(m), m$bx, m$sx); bys <- rnorm(nrow(m), m$by, m$sy); wm(bxs, bys, m$sy) })) else NA_real_
    data.table(exposure = e, outcome = "CAD (Aragam 2022)", exclusion = x, nsnp = r[["n"]], b = r[["b"]], se = r[["se"]], p = r[["p"]],
               OR = exp(r[["b"]]), lo = exp(r[["b"]] - 1.96 * r[["se"]]), hi = exp(r[["b"]] + 1.96 * r[["se"]]),
               wm_OR = exp(wmb), wm_p = 2 * pnorm(-abs(wmb / wmse)),
               n_proxy = sum(!is.na(m$proxy)))
  }))
}))
fwrite(cadmr, file.path(OUT, "amd_cad_mr.tsv"), sep = "\t"); print(cadmr)
message("story_01 done")
