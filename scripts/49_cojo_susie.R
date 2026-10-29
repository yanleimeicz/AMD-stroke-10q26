# Conditional analysis (GCTA-COJO) and multi-signal fine-mapping colocalisation (SuSiE-RSS + coloc.susie) at 10q26 / 12q24.
# Question (report section 4): is the 10q26 AMD signal (rs3750848 / rs36212733 haplotype) the same variant as the stroke signal
# (rs2672592, r2 ~ 0.36)? Which variant drives HTRA1 expression in retina / artery / brain / blood? 12q24: rs3184504 vs SH2B3 eQTL.
# LD reference: 1000G EUR (n = 503) -> LD mismatch with FinnGen / GIGASTROKE is expected; kriging_rss diagnostics + outlier-removed refit.
# Inputs: derived/smr_ma/*.ma, derived/tissue_eqtl/{region_snps.tsv, gwas_regions.rds, eqtl_harmonised.rds}. Output -> results/tissue_eqtl/
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressMessages({ library(data.table); library(parallel); library(susieR); library(coloc); library(ggplot2) })
D <- file.path(PROJ_ROOT, "derived/tissue_eqtl"); OUT <- file.path(PROJ_ROOT, "results/tissue_eqtl"); W <- file.path(D, "cojo")
dir.create(W, recursive = TRUE, showWarnings = FALSE)
BFILE <- file.path(PROJ_ROOT, "data/raw/ref/g1000_eur/g1000_eur")
PLINK <- file.path(PROJ_ROOT, "data/raw/tools/plink/plink.exe")
GCTA <- list.files(file.path(PROJ_ROOT, "data/raw/tools"), pattern = "^gcta.*\\.exe$", recursive = TRUE, full.names = TRUE)[1]
stopifnot(!is.na(GCTA))
rsnps <- fread(file.path(D, "region_snps.tsv"))[!is.na(pos37)]
gw <- readRDS(file.path(D, "gwas_regions.rds"))
eq <- readRDS(file.path(D, "eqtl_harmonised.rds"))
frq <- fread(file.path(PROJ_ROOT, "derived/g1000_eur.frq"))[, .(SNP, fA1 = A1, MAF)]
bim <- fread(paste0(BFILE, ".bim"), header = FALSE, select = c(1, 2, 4, 5, 6), col.names = c("chr", "SNP", "pos37", "A1", "A2"))[SNP %in% rsnps$SNP]

# ---------------- COJO inputs: GWAS (.ma from 38a_make_ma.R) and tissue eQTL (HTRA1 / ARMS2 / SH2B3) per region
CJ_REG <- c("10q26_HTRA1" = "10", "12q24_SH2B3" = "12")
TRAITS <- c("AMD_FG", "AMD_IAMDGC", "DRYAMD_FG", "WETAMD_FG", "AS", "IS", "LAS", "SVS")
for (r in names(CJ_REG)) {   # small region-level reference panels (GCTA otherwise re-reads the 22M-variant bim per job)
  fwrite(data.table(rsnps[region == r]$SNP), file.path(W, paste0(r, ".snps")), col.names = FALSE)
  system2(PLINK, c("--bfile", BFILE, "--extract", file.path(W, paste0(r, ".snps")), "--make-bed", "--out", file.path(W, paste0("ref_", r))),
          stdout = FALSE, stderr = FALSE)
}
inputs <- list()
for (k in TRAITS) {
  ma <- fread(sprintf(file.path(PROJ_ROOT, "derived/smr_ma/%s.ma"), k))
  for (r in names(CJ_REG)) {
    f <- file.path(W, sprintf("%s__%s.ma", k, r)); fwrite(ma[SNP %in% rsnps[region == r]$SNP], f, sep = " ")
    inputs[[f]] <- data.table(source = "GWAS", name = k, region = r, ma = f)
  }
}
EQ_GENES <- c("HTRA1", "ARMS2", "SH2B3")
for (key in unique(eq[gene %in% EQ_GENES, paste(region, gene, dataset, sep = "|")])) {
  k <- strsplit(key, "\\|")[[1]]; if (!k[1] %in% names(CJ_REG)) next
  e <- eq[region == k[1] & gene == k[2] & dataset == k[3]]
  if (min(e$p) > 1e-6) next
  e <- merge(e, frq, by = "SNP"); e[, freq := fifelse(toupper(ea) == fA1, MAF, 1 - MAF)]
  f <- file.path(W, sprintf("eQTL_%s_%s__%s.ma", k[2], k[3], k[1]))
  fwrite(e[, .(SNP, A1 = toupper(ea), A2 = toupper(oa), freq, b = beta, se, p, n)], f, sep = " ")
  inputs[[f]] <- data.table(source = "eQTL", name = paste(k[2], k[3], sep = "_"), region = k[1], ma = f)
}
inputs <- rbindlist(inputs)
COND <- list("10q26_HTRA1" = list(rs3750848 = "rs3750848", rs2672592 = "rs2672592", both = c("rs3750848", "rs2672592")),
             "12q24_SH2B3" = list(rs3184504 = "rs3184504", rs587914 = "rs587914", both = c("rs3184504", "rs587914")))
for (r in names(COND)) for (cn in names(COND[[r]])) fwrite(data.table(COND[[r]][[cn]]), file.path(W, sprintf("cond_%s_%s.txt", r, cn)), col.names = FALSE)
jobs <- rbind(inputs[, .(source, name, region, ma, mode = "slct", cond = NA_character_)],
              inputs[, CJ(cond = names(COND[[region[1]]])), by = .(source, name, region, ma)][, mode := "cond"], fill = TRUE)
jobs[, out := file.path(W, sprintf("%s__%s__%s%s", make.names(name), region, mode, fifelse(is.na(cond), "", paste0("_", cond))))]
run_one <- function(i) {
  j <- jobs[i]
  args <- c("--bfile", file.path(W, paste0("ref_", j$region)),
            "--cojo-file", j$ma, "--cojo-collinear", "0.9", "--cojo-wind", "2000", "--thread-num", "1", "--out", j$out)
  args <- if (j$mode == "slct") c(args, "--cojo-slct", "--cojo-p", "5e-8") else
    c(args, "--cojo-cond", file.path(W, sprintf("cond_%s_%s.txt", j$region, j$cond)))
  system2(GCTA, args, stdout = paste0(j$out, ".stdout"), stderr = paste0(j$out, ".stdout"))
  j$out
}
cl <- makeCluster(30); clusterExport(cl, c("jobs", "GCTA", "W", "CJ_REG")); clusterEvalQ(cl, library(data.table))
invisible(parLapply(cl, seq_len(nrow(jobs)), run_one)); stopCluster(cl)

KEYS <- c("rs3750848", "rs36212733", "rs2672592", "rs3184504", "rs587914")
slct <- rbindlist(lapply(which(jobs$mode == "slct"), function(i) {
  f <- paste0(jobs$out[i], ".jma.cojo"); if (!file.exists(f)) return(data.table(source = jobs$source[i], name = jobs$name[i], region = jobs$region[i], n_signals = 0L))
  d <- fread(f); data.table(source = jobs$source[i], name = jobs$name[i], region = jobs$region[i], n_signals = nrow(d),
                            signals = paste(sprintf("%s(p=%.1e,pJ=%.1e)", d$SNP, d$p, d$pJ), collapse = "; "))
}), fill = TRUE)
cond <- rbindlist(lapply(which(jobs$mode == "cond"), function(i) {
  f <- paste0(jobs$out[i], ".cma.cojo"); if (!file.exists(f)) return(NULL)
  d <- fread(f); top <- d[which.min(pC)]
  key <- d[SNP %in% KEYS, .(SNP, p, pC)]
  data.table(source = jobs$source[i], name = jobs$name[i], region = jobs$region[i], conditioned_on = jobs$cond[i],
             min_p_before = min(d$p), top_after = top$SNP, top_pC = top$pC, n_pC_5e8 = sum(d$pC < 5e-8, na.rm = TRUE),
             key_snps = paste(sprintf("%s:p=%.1e->pC=%.1e", key$SNP, key$p, key$pC), collapse = "; "))
}))
fwrite(slct[order(region, source, name)], file.path(OUT, "cojo_slct.tsv"), sep = "\t")
fwrite(cond[order(region, source, name, conditioned_on)], file.path(OUT, "cojo_cond.tsv"), sep = "\t")
message("COJO done: ", nrow(slct), " slct, ", nrow(cond), " cond")

# ---------------- two-SNP conditional model (robust to regional LD mismatch: depends only on the pairwise r)
# z_cond(b | a) = (z_b - r z_a) / sqrt(1 - r^2), z aligned to the 1000G A1 of each SNP, r = signed allelic correlation (A1-coded).
pair_r <- function(a, b, chr) {
  o <- file.path(W, sprintf("r_%s_%s", a, b))
  fwrite(data.table(c(a, b)), paste0(o, ".snps"), col.names = FALSE)
  system2(PLINK, c("--bfile", BFILE, "--chr", chr, "--extract", paste0(o, ".snps"), "--r", "--keep-allele-order",
                   "--ld-window-r2", "0", "--out", o), stdout = FALSE, stderr = FALSE)
  fread(paste0(o, ".ld"))$R[1]
}
PAIRS2 <- data.table(region = c("10q26_HTRA1", "12q24_SH2B3"), a = c("rs3750848", "rs3184504"), b = c("rs2672592", "rs587914"), chr = c("10", "12"))
PAIRS2[, r := mapply(pair_r, a, b, chr)]
zA1 <- function(d, s) { x <- merge(bim[SNP == s, .(SNP, A1, A2)], d[SNP == s], by = "SNP")
  if (!nrow(x)) return(NA_real_); with(x[1], fifelse(toupper(ea) == A1, beta / se, fifelse(toupper(ea) == A2, -beta / se, NA_real_))) }
two <- rbindlist(lapply(seq_len(nrow(PAIRS2)), function(i) {
  pr <- PAIRS2[i]
  srcs <- c(lapply(setNames(c(TRAITS, "WMH"), c(TRAITS, "WMH")), function(k) gw[trait == k, .(SNP, ea, beta, se)]),
            lapply(split(eq[region == pr$region & gene %in% EQ_GENES], by = c("gene", "dataset")), function(e) e[, .(SNP, ea, beta, se)]))
  rbindlist(lapply(names(srcs), function(nm) {
    za <- zA1(srcs[[nm]], pr$a); zb <- zA1(srcs[[nm]], pr$b)
    if (is.na(za) || is.na(zb)) return(NULL)
    rbindlist(lapply(c(pr$r - 0.1, pr$r, pr$r + 0.1), function(r) data.table(
      region = pr$region, source = nm, snp_a = pr$a, snp_b = pr$b, r_1000G = pr$r, r_used = round(r, 3), z_a = za, z_b = zb,
      p_a = 2 * pnorm(-abs(za)), p_b = 2 * pnorm(-abs(zb)),
      pC_b_given_a = 2 * pnorm(-abs((zb - r * za) / sqrt(1 - r^2))), pC_a_given_b = 2 * pnorm(-abs((za - r * zb) / sqrt(1 - r^2))),
      zC_b_given_a = (zb - r * za) / sqrt(1 - r^2), zC_a_given_b = (za - r * zb) / sqrt(1 - r^2))))
  }))
}))
fwrite(two, file.path(OUT, "conditional_two_snp.tsv"), sep = "\t")

# ---------------- LD with key SNPs (for plots / credible-set annotation)
ldr2 <- function(snp, chr) {
  o <- file.path(W, paste0("ld_", snp))
  system2(PLINK, c("--bfile", BFILE, "--chr", chr, "--r2", "--ld-snp", snp, "--ld-window-kb", "1500", "--ld-window", "999999",
                   "--ld-window-r2", "0", "--out", o), stdout = FALSE, stderr = FALSE)
  d <- fread(paste0(o, ".ld")); d[, .(SNP = SNP_B, r2 = R2)]
}
LD <- list(rs3750848 = ldr2("rs3750848", "10"), rs2672592 = ldr2("rs2672592", "10"), rs3184504 = ldr2("rs3184504", "12"))

# ---------------- SuSiE-RSS + coloc.susie at 10q26 (HTRA1 +/- 250 kb, GRCh38 122.21-122.72 Mb)
win <- rsnps[region == "10q26_HTRA1" & pos38 >= 122210000 & pos38 <= 122720000 & SNP %in% bim$SNP]
b10 <- bim[SNP %in% win$SNP][order(pos37)]
fwrite(b10[, .(SNP)], file.path(W, "susie10.snps"), col.names = FALSE)
system2(PLINK, c("--bfile", BFILE, "--extract", file.path(W, "susie10.snps"), "--r", "square", "--keep-allele-order",
                 "--out", file.path(W, "susie10")), stdout = FALSE, stderr = FALSE)
R <- as.matrix(fread(file.path(W, "susie10.ld"), header = FALSE)); dimnames(R) <- list(b10$SNP, b10$SNP)
ok <- which(colSums(is.na(R)) == 0); R <- R[ok, ok]; b10 <- b10[ok]
zalign <- function(d) {      # z aligned to the 1000G A1 (plink --keep-allele-order -> A1 is the counted allele)
  d <- merge(b10[, .(SNP, A1, A2)], d, by = "SNP")
  d[, z := fifelse(toupper(ea) == A1, beta / se, fifelse(toupper(ea) == A2, -beta / se, NA_real_))]
  d <- d[is.finite(z)]; d[match(b10$SNP, SNP)][!is.na(SNP)]
}
SETS <- list()
for (k in c("AMD_FG", "AMD_IAMDGC", "DRYAMD_FG", "WETAMD_FG", "IS", "SVS", "AS")) {
  d <- zalign(gw[trait == k, .(SNP, ea, oa, beta, se, n = fifelse(is.na(ncase), n, ncase + ncontrol))]); SETS[[k]] <- d
}
for (key in unique(eq[region == "10q26_HTRA1" & gene %in% c("HTRA1", "ARMS2", "PLEKHA1"), paste(gene, dataset, sep = "|")])) {
  k <- strsplit(key, "\\|")[[1]]; e <- eq[region == "10q26_HTRA1" & gene == k[1] & dataset == k[2]]
  if (min(e$p) > 1e-6) next
  SETS[[paste0("eQTL_", k[1], "_", k[2])]] <- zalign(e[, .(SNP, ea, oa, beta, se, n)])
}
fit_one <- function(nm) {
  d <- SETS[[nm]]; Rs <- R[d$SNP, d$SNP]; n <- max(d$n, na.rm = TRUE)
  kr <- tryCatch(kriging_rss(d$z, Rs, n = n)$conditional_dist, error = function(e) NULL)
  outl <- if (!is.null(kr)) d$SNP[kr$logLR > 2 & abs(kr$z) > 2] else character()
  fit <- function(dd) {
    s <- tryCatch(susie_rss(z = dd$z, R = R[dd$SNP, dd$SNP], n = n, L = 10, estimate_residual_variance = FALSE,
                            max_iter = 1000, coverage = 0.95, min_abs_corr = 0.5), error = function(e) NULL)
    if (!is.null(s)) { names(s$pip) <- dd$SNP; colnames(s$alpha) <- dd$SNP }
    s
  }
  list(name = nm, n = n, nsnp = nrow(d), outliers = outl, full = fit(d), clean = if (length(outl)) fit(d[!SNP %in% outl]) else NULL)
}
cl <- makeCluster(min(30, length(SETS))); clusterExport(cl, c("SETS", "R")); clusterEvalQ(cl, { .libPaths(c(file.path(Sys.getenv("COMORBIDITY_ROOT"), "Rlib"), .libPaths())); library(susieR); library(data.table) })
FITS <- parLapply(cl, names(SETS), fit_one); names(FITS) <- names(SETS); stopCluster(cl)
r2of <- function(s, ref) { v <- LD[[ref]][SNP == s]$r2; if (length(v)) v[1] else if (s == ref) 1 else NA_real_ }
cs_tab <- rbindlist(lapply(FITS, function(f) rbindlist(lapply(c("full", "clean"), function(v) {
  s <- f[[v]]; if (is.null(s) || is.null(s$sets$cs)) return(data.table(dataset = f$name, version = v, n_cs = 0L, n_outliers = length(f$outliers)))
  rbindlist(lapply(names(s$sets$cs), function(cn) {
    idx <- s$sets$cs[[cn]]; snp <- names(s$pip)[idx]
    lead <- snp[which.max(s$pip[idx])]
    data.table(dataset = f$name, version = v, n_outliers = length(f$outliers), cs = cn, cs_size = length(idx), lead = lead,
               lead_pip = max(s$pip[idx]), r2_lead_rs3750848 = r2of(lead, "rs3750848"), r2_lead_rs2672592 = r2of(lead, "rs2672592"),
               has_rs3750848 = "rs3750848" %in% snp, has_rs36212733 = "rs36212733" %in% snp, has_rs2672592 = "rs2672592" %in% snp,
               snps = paste(head(snp[order(-s$pip[idx])], 15), collapse = ","))
  }))
}), fill = TRUE)), fill = TRUE)
fwrite(cs_tab, file.path(OUT, "susie_10q26_credible_sets.tsv"), sep = "\t")

PAIRS <- rbind(CJ(a = c("AMD_FG", "AMD_IAMDGC", "DRYAMD_FG", "WETAMD_FG"), b = c("IS", "SVS", "AS")),
               CJ(a = c("AMD_FG", "AMD_IAMDGC", "IS", "SVS", "AS"), b = grep("^eQTL_", names(SETS), value = TRUE)))
cs_res <- rbindlist(lapply(seq_len(nrow(PAIRS)), function(i) rbindlist(lapply(c("full", "clean"), function(v) {
  s1 <- FITS[[PAIRS$a[i]]][[v]] %||% FITS[[PAIRS$a[i]]]$full; s2 <- FITS[[PAIRS$b[i]]][[v]] %||% FITS[[PAIRS$b[i]]]$full
  if (is.null(s1) || is.null(s2) || is.null(s1$sets$cs) || is.null(s2$sets$cs)) return(data.table(a = PAIRS$a[i], b = PAIRS$b[i], version = v, note = "no credible set in one trait"))
  r <- tryCatch(coloc.susie(s1, s2), error = function(e) NULL)
  if (is.null(r) || is.null(r$summary)) return(data.table(a = PAIRS$a[i], b = PAIRS$b[i], version = v, note = "coloc.susie failed / no overlap"))
  sm <- as.data.table(r$summary)
  sm[, .(a = PAIRS$a[i], b = PAIRS$b[i], version = v, hit1, hit2, nsnps, PP.H3 = PP.H3.abf, PP.H4 = PP.H4.abf,
         hit1_r2_rs3750848 = sapply(hit1, r2of, "rs3750848"), hit2_r2_rs3750848 = sapply(hit2, r2of, "rs3750848"),
         hit1_r2_rs2672592 = sapply(hit1, r2of, "rs2672592"), hit2_r2_rs2672592 = sapply(hit2, r2of, "rs2672592"))]
}), fill = TRUE)), fill = TRUE)
fwrite(cs_res, file.path(OUT, "coloc_susie_10q26.tsv"), sep = "\t")
diag <- rbindlist(lapply(FITS, function(f) data.table(dataset = f$name, n = f$n, nsnp = f$nsnp, n_kriging_outliers = length(f$outliers),
                                                      outliers = paste(head(f$outliers, 20), collapse = ","))))
fwrite(diag, file.path(OUT, "susie_10q26_ld_mismatch_diagnostics.tsv"), sep = "\t")
message("SuSiE done")

# ---------------- LocusZoom-style comparison plots
pos38 <- rsnps[, .(SNP, pos38)]
panel <- function(d, lab) { d <- merge(d[, .(SNP, p)], pos38, by = "SNP"); d[, panel := lab]; d }
plot_region <- function(reg, lo, hi, gw_tr, eq_gene, ref, marks, file) {
  P <- rbind(rbindlist(lapply(gw_tr, function(k) panel(gw[trait == k & region == reg], k))),
             rbindlist(lapply(unique(eq[region == reg & gene == eq_gene]$dataset), function(ds) {
               e <- eq[region == reg & gene == eq_gene & dataset == ds]; if (min(e$p) > 1e-4) return(NULL); panel(e, paste(eq_gene, ds)) })))
  P <- P[pos38 >= lo & pos38 <= hi]
  P <- merge(P, LD[[ref]], by = "SNP", all.x = TRUE); P[SNP == ref, r2 := 1]
  P[, lp := pmin(-log10(p), 300)]; P[, panel := factor(panel, levels = unique(panel))]
  P[, r2bin := cut(r2, c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf), labels = c("<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", ">0.8"))]
  M <- P[SNP %in% marks]
  g <- ggplot(P[order(r2, na.last = FALSE)], aes(pos38 / 1e6, lp, colour = r2bin)) + geom_point(size = 0.8) +
    geom_point(data = M, aes(pos38 / 1e6, lp), inherit.aes = FALSE, shape = 23, size = 2.4, fill = "purple", colour = "black") +
    geom_text(data = M, aes(pos38 / 1e6, lp, label = SNP), inherit.aes = FALSE, size = 2.3, vjust = -0.9) +
    facet_wrap(~panel, ncol = 1, scales = "free_y", strip.position = "right") +
    scale_colour_manual(values = c("<0.2" = "grey65", "0.2-0.4" = "#86b3d1", "0.4-0.6" = "#2c7fb8", "0.6-0.8" = "#f39c12", ">0.8" = "#d7191c"),
                        na.value = "grey85", name = paste0("r2 with ", ref)) +
    labs(x = "chr position (GRCh38, Mb)", y = "-log10 P") + theme_bw(base_size = 8) + theme(strip.text.y = element_text(angle = 0))
  ggsave(file, g, width = 7.5, height = 1.15 * nlevels(P$panel) + 0.8, limitsize = FALSE)
  ggsave(sub("\\.pdf$", ".png", file), g, width = 7.5, height = 1.15 * nlevels(P$panel) + 0.8, dpi = 200, limitsize = FALSE)
}
plot_region("10q26_HTRA1", 122.25e6, 122.65e6, c("AMD_FG", "IS", "SVS"), "HTRA1", "rs3750848",
            c("rs3750848", "rs2672592"), file.path(OUT, "Fig_region_10q26_HTRA1.pdf"))
plot_region("12q24_SH2B3", 111.1e6, 111.8e6, c("AMD_IAMDGC", "DRYAMD_FG", "IS", "AS"), "SH2B3", "rs3184504",
            c("rs3184504", "rs587914"), file.path(OUT, "Fig_region_12q24_SH2B3.pdf"))
message("49 done")
