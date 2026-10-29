# Summary figures: MR forest plots (main vs excluding major AMD loci) and colocalisation heatmap.
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
MR <- file.path(PROJ_ROOT, "results/mr"); CO <- file.path(PROJ_ROOT, "results/coloc")
lab <- c(AMD_IAMDGC = "AMD (IAMDGC)", AMD_FG = "AMD (FinnGen)", DRYAMD_FG = "Dry AMD", WETAMD_FG = "Wet AMD", AS = "Any stroke",
         IS = "Ischemic stroke", LAS = "Large artery", CES = "Cardioembolic", SVS = "Small vessel", WMH = "WMH volume")
r <- fread(file.path(MR, "mr_all_methods.tsv"))[method %in% c("Inverse variance weighted", "Wald ratio")]
r[, `:=`(Exposure = lab[exposure], Outcome = factor(lab[outcome], rev(lab)))]
r[, Analysis := factor(fifelse(analysis == "main", "All instruments", fifelse(analysis == "excl_major_AMD_loci", "Excl. CFH/ARMS2-HTRA1/C3",
                       fifelse(analysis == "excl_lipid_loci", "Excl. lipid loci", "Excl. major + lipid"))),
                       c("All instruments", "Excl. CFH/ARMS2-HTRA1/C3", "Excl. lipid loci", "Excl. major + lipid"))]
fwd <- r[exposure %like% "AMD" & analysis %in% c("main", "excl_major_AMD_loci")]
p1 <- ggplot(fwd, aes(OR, Outcome, colour = Analysis)) + geom_vline(xintercept = 1, lty = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = .6), size = .25) +
  facet_wrap(~Exposure, nrow = 1) + scale_x_log10() + theme_bw(base_size = 9) + scale_colour_manual(values = c("#2b5c8a", "#c0392b")) +
  labs(x = "OR per unit log-odds of genetic liability (IVW)", y = NULL, title = "Forward MR: genetic liability to AMD -> stroke subtypes / WMH") +
  theme(legend.position = "bottom")
ggsave(file.path(MR, "Fig_MR_forward.pdf"), p1, width = 11, height = 4); ggsave(file.path(MR, "Fig_MR_forward.png"), p1, width = 11, height = 4, dpi = 300)
rev_ <- r[!exposure %like% "AMD" & analysis == "main"]
rev_[, Exposure := factor(Exposure, lab)][, Outcome := lab[outcome]]
p2 <- ggplot(rev_, aes(OR, Exposure)) + geom_vline(xintercept = 1, lty = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = .25, colour = "#2b5c8a") + facet_wrap(~Outcome, nrow = 1) + scale_x_log10() +
  theme_bw(base_size = 9) + labs(x = "OR (IVW)", y = NULL, title = "Reverse MR: stroke subtypes / WMH -> AMD")
ggsave(file.path(MR, "Fig_MR_reverse.pdf"), p2, width = 11, height = 3.4); ggsave(file.path(MR, "Fig_MR_reverse.png"), p2, width = 11, height = 3.4, dpi = 300)

# all five methods for main analyses (supplementary forest)
m5 <- fread(file.path(MR, "mr_all_methods.tsv"))[analysis == "main"]
m5[, pair := paste(lab[exposure], "->", lab[outcome])]
p3 <- ggplot(m5, aes(OR, method)) + geom_vline(xintercept = 1, lty = 2, colour = "grey50") + geom_pointrange(aes(xmin = lo, xmax = hi), size = .15) +
  facet_wrap(~pair, ncol = 6, scales = "free_x") + theme_bw(base_size = 6) + labs(x = "OR", y = NULL)
ggsave(file.path(MR, "FigS_MR_all_methods.pdf"), p3, width = 14, height = 16)

# colocalisation heatmap (max PP.H4 per region x AMD trait x stroke trait)
cl <- fread(file.path(CO, "coloc_gwas_gwas.tsv"))
cl[, amd := fifelse(source %like% "AMD", source, partner)][, stk := fifelse(source %like% "AMD", partner, source)]
cl[, region := paste0("chr", chr, ":", round(pos38 / 1e6, 1), "Mb")]
ann <- fread(file.path(CO, "coloc_shared_loci_annotated.tsv"))[, .(region2 = region, top_nearest)]
best <- cl[, .(PP.H4 = max(PP.H4)), by = .(chr, bin = floor(pos38 / 2e6), amd, stk)]
keep <- best[PP.H4 > 0.5, unique(paste(chr, bin))]
best <- best[paste(chr, bin) %in% keep]
gene <- c("10 61" = "ARMS2/HTRA1 (10q26)", "3 32" = "ADAMTS9 (3p14)", "9 66" = "ABO (9q34)", "7 46" = "CDK6 (7q21)", "1 97" = "1q31", "19 5" = "SMARCA4/LDLR (19p13)")
best[, locus := gene[paste(chr, bin)]][is.na(locus), locus := paste0("chr", chr, " bin", bin)]
best[, pair := paste(lab[amd], "~", lab[stk])]
p4 <- ggplot(best, aes(pair, locus, fill = PP.H4)) + geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.2f", PP.H4)), size = 2.2) +
  scale_fill_gradient(low = "white", high = "#c0392b", limits = c(0, 1)) + theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + labs(x = NULL, y = NULL, title = "GWAS-GWAS colocalisation (posterior probability of a shared causal variant)")
ggsave(file.path(CO, "Fig_coloc_heatmap.pdf"), p4, width = 11, height = 3.8); ggsave(file.path(CO, "Fig_coloc_heatmap.png"), p4, width = 11, height = 3.8, dpi = 300)
message("figures done")
