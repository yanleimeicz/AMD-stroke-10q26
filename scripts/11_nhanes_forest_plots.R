# Forest plots for NHANES stroke x AMD results (main/sensitivity and subgroups)
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(forestploter); library(grid) })
OUT <- file.path(PROJ_ROOT, "results/nhanes")
fmt_ci <- function(or, lo, hi) ifelse(is.na(or), "", sprintf("%.2f (%.2f-%.2f)", or, lo, hi))
fmt_p  <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
th <- forest_theme(base_size = 9, ci_pch = 15, ci_col = "#2b5c8a", refline_gp = gpar(lty = "dashed", col = "grey40"),
                   core = list(bg_params = list(fill = c("white"))))

# ---- Figure A: sequential models + sensitivity
m <- read.delim(file.path(OUT, "main_models.tsv")); s <- read.delim(file.path(OUT, "sensitivity.tsv"))
lab_term <- c(amd_any = "Any AMD", "amd_catEarly AMD" = "  Early AMD", "amd_catLate AMD" = "  Late AMD",
              amd_self = "Self-reported AMD", stroke = "Stroke (reverse)")
m$label <- paste0(m$model, ": ", lab_term[m$term])
s$label <- s$model
a <- rbind(
  data.frame(label = "Sequential models (outcome: stroke)", n = NA, OR = NA, lo = NA, hi = NA, p = NA),
  m[, c("label", "n", "OR", "lo", "hi", "p")],
  data.frame(label = "Sensitivity analyses", n = NA, OR = NA, lo = NA, hi = NA, p = NA),
  s[, c("label", "n", "OR", "lo", "hi", "p")])
a$` ` <- paste(rep(" ", 22), collapse = " ")
a$`OR (95% CI)` <- fmt_ci(a$OR, a$lo, a$hi); a$P <- fmt_p(a$p); a$N <- ifelse(is.na(a$n), "", a$n)
names(a)[names(a) == "label"] <- "Model / analysis"
p1 <- forest(a[, c("Model / analysis", "N", " ", "OR (95% CI)", "P")], est = a$OR, lower = a$lo, upper = a$hi,
             ci_column = 3, ref_line = 1, x_trans = "log", xlim = c(0.5, 30), ticks_at = c(0.5, 1, 2, 5, 10, 25), theme = th)
p1 <- edit_plot(p1, row = which(is.na(a$OR)), gp = gpar(fontface = "bold"))
ggplot2::ggsave(file.path(OUT, "Fig_forest_models_sensitivity.pdf"), p1, width = 9.5, height = 6.2)
png(file.path(OUT, "Fig_forest_models_sensitivity.png"), width = 9.5, height = 6.2, units = "in", res = 300); plot(p1); dev.off()

# ---- Figure B: subgroups
g <- read.delim(file.path(OUT, "subgroups.tsv"))
glab <- c(age_grp = "Age, years", sex = "Sex", race2 = "Race/ethnicity", smoke2 = "Smoking", obese = "BMI, kg/m2",
          htn_f = "Hypertension", dm_f = "Diabetes")
rows <- list()
for (k in unique(g$model)) {
  gi <- g[g$model == k, ]
  rows[[k]] <- rbind(data.frame(label = glab[k], ev = "", OR = NA, lo = NA, hi = NA, p = NA, pint = gi$p_interaction[1]),
                     data.frame(label = paste0("   ", gi$level), ev = paste0(gi$events_amd, "/", gi$n_amd),
                                OR = gi$OR, lo = gi$lo, hi = gi$hi, p = gi$p, pint = NA))
}
b <- do.call(rbind, rows)
b$` ` <- paste(rep(" ", 22), collapse = " ")
b$`OR (95% CI)` <- fmt_ci(b$OR, b$lo, b$hi); b$P <- fmt_p(b$p); b$`P interaction` <- fmt_p(b$pint)
names(b)[names(b) == "ev"] <- "Stroke/AMD n"
names(b)[names(b) == "label"] <- "Subgroup"
p2 <- forest(b[, c("Subgroup", "Stroke/AMD n", " ", "OR (95% CI)", "P", "P interaction")], est = b$OR, lower = b$lo, upper = b$hi,
             ci_column = 3, ref_line = 1, x_trans = "log", xlim = c(0.25, 12), ticks_at = c(0.25, 0.5, 1, 2, 4, 8), theme = th)
p2 <- edit_plot(p2, row = which(is.na(b$OR)), gp = gpar(fontface = "bold"))
ggplot2::ggsave(file.path(OUT, "Fig_forest_subgroups.pdf"), p2, width = 9, height = 6)
png(file.path(OUT, "Fig_forest_subgroups.png"), width = 9, height = 6, units = "in", res = 300); plot(p2); dev.off()
cat("figures written\n")
