# Story revision: 10q26 haplotype x trait matrix (joint two-variant effects) and complement pQTL summary.
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages(library(data.table))
R <- file.path(PROJ_ROOT, "results"); OUT <- file.path(R, "story"); dir.create(OUT, showWarnings = FALSE)
LEAD <- "rs3750848"; SIDE <- "rs2672592"; RISK <- c(rs3750848 = "G", rs2672592 = "T")
RV <- c(KG_EUR = 0.618701, FinnGen_SISu = 0.6636264009214823, UKB_British = 0.659085214138031, KG_EAS = 0.821883203381113)

# joint two-variant model from marginal (beta, se) on the risk-allele scale (COJO approximation):
# z_joint = R^-1 z ; beta_joint_i = z_joint_i * se_i ; se_joint_i = se_i * sqrt((R^-1)_ii)
joint2 <- function(b1, s1, b2, s2, r) {
  z <- c(b1 / s1, b2 / s2); Ri <- solve(matrix(c(1, r, r, 1), 2)); zj <- as.vector(Ri %*% z)
  bj <- zj * c(s1, s2); sj <- c(s1, s2) * sqrt(diag(Ri)); data.table(bl = bj[1], sl = sj[1], bs = bj[2], ss = sj[2])
}
rows <- list()
add <- function(trait, dataset, ld, b1, s1, b2, s2, scale = "log-OR") {
  j <- joint2(b1, s1, b2, s2, RV[[ld]])
  rows[[length(rows) + 1]] <<- data.table(trait, dataset, LD = ld, scale,
    beta_lead_marginal = b1, se_lead_marginal = s1, beta_side_marginal = b2, se_side_marginal = s2,
    beta_lead_joint = j$bl, se_lead_joint = j$sl, p_lead_joint = 2 * pnorm(-abs(j$bl / j$sl)),
    beta_side_joint = j$bs, se_side_joint = j$ss, p_side_joint = 2 * pnorm(-abs(j$bs / j$ss)))
}
kv <- fread(file.path(R, "revision/10q26/key_variants_all_datasets.tsv"))
getkv <- function(tr) { x <- kv[trait == tr & rsid %in% c(LEAD, SIDE)]; x <- x[match(c(LEAD, SIDE), rsid)]; if (anyNA(x$beta_risk)) NULL else x }
for (tr in list(c("AMD", "FinnGen_H7_AMD", "FinnGen_SISu"), c("Dry AMD", "FinnGen_DRY_AMD", "FinnGen_SISu"), c("Wet AMD", "FinnGen_WET_AMD", "FinnGen_SISu"),
                c("SVS", "GIGA_EUR_SVS", "KG_EUR"), c("AS", "GIGA_EUR_AS", "KG_EUR"), c("IS", "GIGA_EUR_IS", "KG_EUR"), c("LAS", "GIGA_EUR_LAS", "KG_EUR"),
                c("CES", "GIGA_EUR_CES", "KG_EUR"), c("Lacunar stroke", "LACUNAR", "KG_EUR"), c("WMH", "WMH", "KG_EUR"),
                c("IS", "GIGA_EAS_IS", "KG_EAS"), c("SVS", "GIGA_EAS_SVS", "KG_EAS"), c("AS", "GIGA_EAS_AS", "KG_EAS"))) {
  x <- getkv(tr[2]); if (is.null(x)) { message("missing ", tr[2]); next }
  add(tr[1], tr[2], tr[3], x$beta_risk[1], x$se[1], x$beta_risk[2], x$se[2], ifelse(tr[1] == "WMH", "SD", "log-OR"))
}
# MEGASTROKE
for (o in c("AS", "IS", "SVS", "LAS", "CES")) {
  f <- file.path(file.path(PROJ_ROOT, "derived/revision"), paste0("MEGA_", o, ".rds")); if (!file.exists(f)) next
  d <- readRDS(f)[SNP %in% c(LEAD, SIDE)]; if (nrow(d) < 2) next
  d[, b := ifelse(ea == RISK[SNP], beta, -beta)]; d <- d[match(c(LEAD, SIDE), SNP)]
  add(o, paste0("MEGASTROKE_", o), "KG_EUR", d$b[1], d$se[1], d$b[2], d$se[2])
}
# UK Biobank (Pan-UKB, in-sample LD)
pu <- fread(file.path(R, "revision3/panukb/panukb_10q26_variants.tsv"))
for (tr in list(c("AMD", "AMD_phe362.29"), c("CAD (coronary atherosclerosis)", "CAD_phe411.4"), c("IHD", "IHD_phe411"), c("IS", "IS_I63"),
                c("AS (self-report)", "AS_sr6150_3"), c("TIA", "TIA_G45"), c("ICH", "ICH_phe430.2"))) {
  x <- pu[trait == tr[2] & rsid %in% c(LEAD, SIDE)]; x <- x[match(c(LEAD, SIDE), rsid)]; if (anyNA(x$beta_risk)) next
  add(tr[1], paste0("UKB_", tr[2]), "UKB_British", x$beta_risk[1], x$se[1], x$beta_risk[2], x$se[2])
}
# CAD (Aragam 2022; region file)
cd <- fread(file.path(PROJ_ROOT, "derived/revision_adamts9/CAD_Aragam2022.tsv.gz"))[SNP %in% c(LEAD, SIDE)]
if (nrow(cd) == 2) { cd[, b := ifelse(ea == RISK[SNP], beta, -beta)]; cd <- cd[match(c(LEAD, SIDE), SNP)]
  add("CAD", "Aragam2022", "KG_EUR", cd$b[1], cd$se[1], cd$b[2], cd$se[2]) } else message("CAD Aragam: key variants missing")
mat <- rbindlist(rows)
# FinnGen I9_STR, GIGASTROKE remainder and random-effects pooled IS (from revision4)
ds <- fread(file.path(R, "revision4/10q26_stroke_datasets.tsv"))
for (s in c("IS", "AS")) for (dd in unique(ds[subtype == s]$dataset)) {
  a <- ds[subtype == s & dataset == dd & quantity == "joint_rs2672592_lead"]; b <- ds[subtype == s & dataset == dd & quantity == "joint_rs2672592_side"]
  if (!nrow(a) || !nrow(b)) next
  mat <- rbind(mat, data.table(trait = s, dataset = paste0(dd, " (revision4)"), LD = a$ld[1], scale = "log-OR",
    beta_lead_joint = a$beta[1], se_lead_joint = a$se[1], p_lead_joint = a$p[1], beta_side_joint = b$beta[1], se_side_joint = b$se[1], p_side_joint = b$p[1]), fill = TRUE)
}
re <- fread(file.path(R, "revision4/10q26_stroke_random_effects.tsv"))
for (s in c("IS", "AS")) {
  x <- re[subtype == s & grepl("^main", analysis)]
  a <- x[quantity == "joint_rs2672592_lead"]; b <- x[quantity == "joint_rs2672592_side"]
  if (nrow(a) && nrow(b)) mat <- rbind(mat, data.table(trait = s, dataset = "Random-effects meta-analysis (4 European datasets; DL)", LD = "matched", scale = "log-OR",
    beta_lead_joint = a$DL_est[1], se_lead_joint = a$DL_se[1], p_lead_joint = a$DL_p[1], beta_side_joint = b$DL_est[1], se_side_joint = b$DL_se[1], p_side_joint = b$DL_p[1]), fill = TRUE)
}
fwrite(mat, file.path(OUT, "haplotype_trait_matrix.tsv"), sep = "\t")
print(mat[, .(trait, dataset, LD, bl = round(beta_lead_joint, 4), pl = signif(p_lead_joint, 2), bs = round(beta_side_joint, 4), ps = signif(p_side_joint, 2))])

# ------------------------------------------------------------------ 4. complement pQTL MR summary (no new computation)
pq <- fread(file.path(R, "pqtl_mr/pqtl_mr_all.tsv"))
G <- c("CFH", "CFB", "CFI", "C3", "C5"); O <- c("AMD_FG", "DRYAMD_FG", "WETAMD_FG", "AMD_IAMDGC", "AS", "IS", "SVS", "LAS", "WMH")
cp <- pq[gene %in% G & outcome %in% O, .(source, gene, protein_id, outcome, group, method, instrument = snps, nsnp, b, se, pval, fdr, OR, lci, uci, region)]
missense <- c(rs2230199 = "C3 p.Arg102Gly (missense)", rs641153 = "CFB p.Arg32Gln (missense)", rs1061170 = "CFH p.Tyr402His (missense)",
              rs800292 = "CFH p.Val62Ile (missense)", rs4151667 = "CFB p.Leu9His (missense)", rs147859257 = "C3 p.Lys155Gln (missense)")
cp[, coding := fifelse(instrument %in% names(missense), missense[instrument], "not coding / not determined")]
amd <- readRDS(file.path(PROJ_ROOT, "derived/gwas/AMD_FG.rds"))[SNP %in% unique(cp$instrument), .(instrument = SNP, amd_gwas_p = p)]
cp <- merge(cp, amd, by = "instrument", all.x = TRUE)
cp[, in_major_amd_locus := region != "" & !is.na(region) | (!is.na(amd_gwas_p) & amd_gwas_p < 5e-8)]
dir <- dcast(cp, gene + outcome ~ source, value.var = "b", fun.aggregate = function(x) x[1])
if (all(c("UKB-PPP", "deCODE") %in% names(dir))) { dir[, platform_direction_agree := sign(`UKB-PPP`) == sign(deCODE)]
  cp <- merge(cp, dir[, .(gene, outcome, platform_direction_agree)], by = c("gene", "outcome"), all.x = TRUE) }
setorder(cp, gene, outcome, source)
fwrite(cp, file.path(OUT, "complement_pqtl_summary.tsv"), sep = "\t")
print(cp[, .(gene, source, outcome, instrument, OR = round(OR, 3), p = signif(pval, 2), fdr = signif(fdr, 2), coding, amd_gwas_p = signif(amd_gwas_p, 2), platform_direction_agree)])
message("story_02 done")
