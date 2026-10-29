# NHANES 2005-2008: age-related macular degeneration (retinal photograph grading) and stroke
# Survey-weighted logistic regression (Taylor linearization), sequential models, sensitivity,
# subgroup and interaction analyses. Outputs -> results/nhanes/
# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
.libPaths(c(file.path(PROJ_ROOT, "Rlib"), .libPaths()))
suppressPackageStartupMessages({ library(haven); library(dplyr); library(survey); library(mice); library(mitools) })
options(survey.lonely.psu = "adjust", width = 200)

ROOT <- PROJ_ROOT
RAW  <- file.path(ROOT, "data/raw/nhanes")
OUT  <- file.path(ROOT, "results/nhanes"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(20261005)

rd <- function(f) { p <- file.path(RAW, paste0(f, ".xpt")); if (file.exists(p)) read_xpt(p) else NULL }
pick <- function(df, cols) df[, intersect(c("SEQN", cols), names(df)), drop = FALSE]
yn <- function(x) ifelse(x == 1, 1, ifelse(x == 2, 0, NA))

# ---------------------------------------------------------------- 1. Build data
build_cycle <- function(s) {
  d <- pick(rd(paste0("DEMO_", s)), c("RIDAGEYR", "RIAGENDR", "RIDRETH1", "DMDEDUC2", "DMDMARTL", "INDFMPIR",
                                       "RIDEXPRG", "SDMVSTRA", "SDMVPSU", "WTMEC2YR"))
  add <- list(
    OPXRET = c("OPDUARM", "OPDUGA", "OPDUEXU", "OPDUARMA", "OPDURL4"),
    VIQ    = c("VIQ310", "VIQ090", "VIQ071"),
    MCQ    = c("MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F", "MCQ180F", "MCQ220"),
    SMQ    = c("SMQ020", "SMQ040"),
    ALQ    = c("ALQ101", "ALQ110", "ALQ120Q"),
    BMX    = c("BMXBMI"),
    BPQ    = c("BPQ020", "BPQ050A", "BPQ080", "BPQ090D"),
    BPX    = c(paste0("BPXSY", 1:4), paste0("BPXDI", 1:4)),
    DIQ    = c("DIQ010", "DIQ050", "DID070"),
    GHB    = c("LBXGH"), GLU = c("LBXGLU"), TCHOL = c("LBXTC"),
    BIOPRO = c("LBXSCR"), CRP = c("LBXCRP"))
  for (f in names(add)) d <- left_join(d, pick(rd(paste0(f, "_", s)), add[[f]]), by = "SEQN")
  d$cycle <- s
  d
}
raw <- bind_rows(build_cycle("D"), build_cycle("E"))

ckd_epi21 <- function(scr, age, female) {   # CKD-EPI 2021 (race-free)
  k <- ifelse(female, 0.7, 0.9); a <- ifelse(female, -0.241, -0.302)
  142 * pmin(scr / k, 1)^a * pmax(scr / k, 1)^(-1.200) * 0.9938^age * ifelse(female, 1.012, 1)
}

d <- raw %>% mutate(
  WT = WTMEC2YR / 2,
  age = RIDAGEYR,
  female = as.integer(RIAGENDR == 2),
  sex = factor(ifelse(female == 1, "Female", "Male"), levels = c("Male", "Female")),
  race = factor(RIDRETH1, 1:5, c("Mexican American", "Other Hispanic", "Non-Hispanic White",
                                 "Non-Hispanic Black", "Other/Multiracial")),
  race = relevel(race, "Non-Hispanic White"),
  edu = factor(case_when(DMDEDUC2 %in% 1:2 ~ "<High school", DMDEDUC2 == 3 ~ "High school",
                         DMDEDUC2 %in% 4:5 ~ ">High school"), levels = c(">High school", "High school", "<High school")),
  pir_cat = factor(case_when(INDFMPIR < 1.3 ~ "<1.3", INDFMPIR < 3.5 ~ "1.3-3.5", INDFMPIR >= 3.5 ~ ">=3.5"),
                   levels = c(">=3.5", "1.3-3.5", "<1.3")),
  married = case_when(DMDMARTL %in% c(1, 6) ~ 1, DMDMARTL %in% 2:5 ~ 0),
  smoke = factor(case_when(SMQ020 == 2 ~ "Never", SMQ020 == 1 & SMQ040 == 3 ~ "Former",
                           SMQ020 == 1 & SMQ040 %in% 1:2 ~ "Current"), levels = c("Never", "Former", "Current")),
  ever_drink = case_when(ALQ101 == 1 ~ 1, ALQ101 == 2 & ALQ110 == 1 ~ 1, ALQ101 == 2 & ALQ110 == 2 ~ 0),
  alcohol = factor(case_when(ever_drink == 0 ~ "Never",
                             ever_drink == 1 & ALQ120Q == 0 ~ "Former",
                             ever_drink == 1 & ALQ120Q >= 1 & ALQ120Q <= 365 ~ "Current"),
                   levels = c("Never", "Former", "Current")),
  bmi = BMXBMI,
  bmi_cat = factor(case_when(bmi < 25 ~ "<25", bmi < 30 ~ "25-30", bmi >= 30 ~ ">=30"), levels = c("<25", "25-30", ">=30")),
  sbp = rowMeans(across(num_range("BPXSY", 1:4)), na.rm = TRUE),
  dbp = rowMeans(across(num_range("BPXDI", 1:4), ~ ifelse(.x == 0, NA, .x)), na.rm = TRUE),
  htn = case_when(BPQ020 == 1 | BPQ050A == 1 | sbp >= 140 | dbp >= 90 ~ 1,
                  BPQ020 == 2 | (!is.na(sbp) & sbp < 140 & !is.na(dbp) & dbp < 90) ~ 0),
  dm = case_when(DIQ010 == 1 | DIQ050 == 1 | DID070 %in% 1 | LBXGH >= 6.5 | LBXGLU >= 126 ~ 1,
                 DIQ010 %in% 2:3 | (!is.na(LBXGH) & LBXGH < 6.5) ~ 0),
  hchol = case_when(BPQ080 == 1 | BPQ090D == 1 | LBXTC >= 240 ~ 1,
                    BPQ080 == 2 | (!is.na(LBXTC) & LBXTC < 240) ~ 0),
  cvd = case_when(MCQ160B == 1 | MCQ160C == 1 | MCQ160D == 1 | MCQ160E == 1 ~ 1,
                  MCQ160B == 2 & MCQ160C == 2 & MCQ160D == 2 & MCQ160E == 2 ~ 0),
  egfr = ckd_epi21(LBXSCR, age, female == 1),
  ckd = as.integer(egfr < 60),
  log_crp = log(LBXCRP + 0.01),
  # outcome / exposure
  stroke = yn(MCQ160F),
  amd_sev = OPDUARM,                                    # 0 none, 1 early, 2 late
  amd_any = as.integer(OPDUARM >= 1),
  amd_cat = factor(OPDUARM, 0:2, c("No AMD", "Early AMD", "Late AMD")),
  amd_self = yn(VIQ310),
  retinopathy = OPDUARMA,
  # subgroup variables
  age_grp = cut(age, c(39, 59, 69, Inf), c("40-59", "60-69", ">=70")),
  race2 = factor(ifelse(race == "Non-Hispanic White", "Non-Hispanic White", "Other"), c("Non-Hispanic White", "Other")),
  smoke2 = factor(ifelse(smoke == "Never", "Never", "Ever"), c("Never", "Ever")),
  obese = factor(ifelse(bmi >= 30, ">=30", "<30"), c("<30", ">=30")),
  htn_f = factor(htn, 0:1, c("No", "Yes")), dm_f = factor(dm, 0:1, c("No", "Yes"))
)

# ---------------------------------------------------------------- 2. Analytic sample + flowchart
core_cov <- c("age", "sex", "race", "edu", "pir_cat", "smoke", "alcohol", "bmi")
full_cov <- c(core_cov, "htn", "dm", "hchol", "cvd", "ckd")
flow <- list()
s0 <- rep(TRUE, nrow(d));                         flow[["NHANES 2005-2008 participants"]] <- sum(s0)
s1 <- s0 & d$age >= 40;                           flow[["Aged >=40 years"]] <- sum(s1)
s2 <- s1 & !is.na(d$amd_sev);                     flow[["With gradable retinal photographs (AMD graded)"]] <- sum(s2)
s3 <- s2 & !is.na(d$stroke);                      flow[["With valid stroke history"]] <- sum(s3)
s4 <- s3 & complete.cases(d[, core_cov]);         flow[["Complete core covariates (analytic sample)"]] <- sum(s4)
s5 <- s4 & complete.cases(d[, full_cov]);         flow[["Complete full-model covariates"]] <- sum(s5)
d$in_core <- s4; d$in_full <- s5; d$in_mi <- s3
flow_df <- data.frame(step = names(flow), n = unlist(flow))
write.table(flow_df, file.path(OUT, "flowchart.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(flow_df)

# questionnaire-based AMD sample for sensitivity S1 (no retinal-image requirement)
d$in_s1 <- d$age >= 40 & !is.na(d$amd_self) & !is.na(d$stroke) & complete.cases(d[, core_cov])

des_all <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WT, nest = TRUE, data = d)
des  <- subset(des_all, in_core)
desF <- subset(des_all, in_full)

# ---------------------------------------------------------------- 3. Table 1 (weighted, by stroke status)
tab1_vars_cont <- c("age", "bmi", "sbp", "egfr", "log_crp")
tab1_vars_cat  <- c("sex", "race", "edu", "pir_cat", "smoke", "alcohol", "htn_f", "dm_f", "hchol", "cvd", "ckd", "amd_cat")
t1 <- list()
for (v in tab1_vars_cont) {
  f <- as.formula(paste0("~", v))
  m <- svyby(f, ~stroke, des, svymean, na.rm = TRUE); mt <- svymean(f, des, na.rm = TRUE)
  p <- tryCatch(svyttest(as.formula(paste0(v, "~stroke")), des)$p.value, error = function(e) NA)
  t1[[v]] <- data.frame(variable = v, level = "mean (SE)",
                        total = sprintf("%.1f (%.2f)", coef(mt), SE(mt)),
                        no_stroke = sprintf("%.1f (%.2f)", m[m$stroke == 0, v], m[m$stroke == 0, "se"]),
                        stroke = sprintf("%.1f (%.2f)", m[m$stroke == 1, v], m[m$stroke == 1, "se"]), p = p)
}
for (v in tab1_vars_cat) {
  f <- as.formula(paste0("~factor(", v, ")"))
  tt <- svyby(f, ~stroke, des, svymean, na.rm = TRUE); tot <- svymean(f, des, na.rm = TRUE)
  p <- tryCatch(svychisq(as.formula(paste0("~", v, "+stroke")), des)$p.value, error = function(e) NA)
  lv <- sub(paste0("factor\\(", v, "\\)"), "", names(coef(tot)))
  un <- table(des$variables[[v]], des$variables$stroke)
  for (i in seq_along(lv)) {
    cn <- names(coef(tot))[i]
    t1[[paste(v, lv[i])]] <- data.frame(variable = v, level = lv[i],
      total = sprintf("%.1f%%", 100 * coef(tot)[i]),
      no_stroke = sprintf("%d (%.1f%%)", un[lv[i], "0"], 100 * tt[tt$stroke == 0, cn]),
      stroke = sprintf("%d (%.1f%%)", un[lv[i], "1"], 100 * tt[tt$stroke == 1, cn]), p = if (i == 1) p else NA)
  }
}
t1 <- bind_rows(t1)
nn <- table(des$variables$stroke)
t1 <- rbind(data.frame(variable = "N (unweighted)", level = "", total = sum(nn), no_stroke = nn["0"], stroke = nn["1"], p = NA), t1)
write.table(t1, file.path(OUT, "table1_by_stroke.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------- 4. Model helpers
fit_or <- function(design, outcome, expo, covs, label = "") {
  f <- as.formula(paste(outcome, "~", paste(c(expo, covs), collapse = " + ")))
  m <- svyglm(f, design = design, family = quasibinomial())
  cf <- summary(m)$coefficients
  rows <- grep(paste0("^", expo), rownames(cf))
  tq <- qt(0.975, m$df.residual)
  data.frame(model = label, term = rownames(cf)[rows], n = nobs(m),
             OR = exp(cf[rows, 1]), lo = exp(cf[rows, 1] - tq * cf[rows, 2]), hi = exp(cf[rows, 1] + tq * cf[rows, 2]),
             p = cf[rows, 4], row.names = NULL)
}
p_trend <- function(design, outcome, covs) {
  f <- as.formula(paste(outcome, "~", paste(c("amd_sev", covs), collapse = " + ")))
  summary(svyglm(f, design = design, family = quasibinomial()))$coefficients["amd_sev", 4]
}
M <- list(
  "Model 1" = character(0),
  "Model 2" = c("age", "sex", "race"),
  "Model 3" = core_cov,
  "Model 4" = full_cov)

# ---------------------------------------------------------------- 5. Main sequential models (AMD -> stroke)
main <- list()
for (k in names(M)) {
  dz <- if (k == "Model 4") desF else des
  main[[paste(k, "any")]] <- fit_or(dz, "stroke", "amd_any", M[[k]], k)
  sev <- fit_or(dz, "stroke", "amd_cat", M[[k]], k)
  sev$p_trend <- p_trend(dz, "stroke", M[[k]])
  main[[paste(k, "sev")]] <- sev
}
main <- bind_rows(main)
# weighted prevalence of stroke by AMD status
prev <- svyby(~stroke, ~amd_cat, des, svymean, na.rm = TRUE)
write.table(prev, file.path(OUT, "stroke_prevalence_by_amd.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(main, file.path(OUT, "main_models.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n== Weighted stroke prevalence by AMD\n"); print(prev)
cat("\n== Main models (exposure AMD, outcome stroke)\n"); print(main, digits = 3)

# ---------------------------------------------------------------- 6. Sensitivity analyses (core Model 3 unless noted)
sens <- list()
# S1 self-reported AMD as exposure (questionnaire sample, no retinal-image requirement)
sens$S1 <- fit_or(subset(des_all, in_s1), "stroke", "amd_self", core_cov, "S1 Self-reported AMD (questionnaire)")
# S2 exclude diabetes
sens$S2 <- fit_or(subset(des, dm %in% 0), "stroke", "amd_any", setdiff(core_cov, character(0)), "S2 Excluding diabetes")
# S3 exclude any retinopathy
sens$S3 <- fit_or(subset(des, retinopathy %in% 0), "stroke", "amd_any", core_cov, "S3 Excluding retinopathy")
# S4 age >= 50
sens$S4 <- fit_or(subset(des, age >= 50), "stroke", "amd_any", core_cov, "S4 Age >=50 years")
# S5 full model + CRP
sens$S5 <- fit_or(subset(desF, !is.na(log_crp)), "stroke", "amd_any", c(full_cov, "log_crp"), "S5 Full model + log CRP")
# S6 reverse direction: stroke -> AMD
sens$S6 <- fit_or(des, "amd_any", "stroke", core_cov, "S6 Reverse: stroke -> any AMD (Model 3)")
sens$S6b <- fit_or(desF, "amd_any", "stroke", full_cov, "S6 Reverse: stroke -> any AMD (Model 4)")
# S7 multiple imputation of covariates (m = 20), full model
mi_vars <- c("SEQN", "SDMVPSU", "SDMVSTRA", "WT", "stroke", "amd_any", full_cov)
dm_in <- d[d$in_mi, mi_vars]
dm_in$htn <- factor(dm_in$htn); dm_in$dm <- factor(dm_in$dm); dm_in$hchol <- factor(dm_in$hchol)
dm_in$cvd <- factor(dm_in$cvd); dm_in$ckd <- factor(dm_in$ckd)
pm <- make.predictorMatrix(dm_in); pm[, c("SEQN", "SDMVPSU", "SDMVSTRA")] <- 0
imp <- mice(dm_in, m = 20, predictorMatrix = pm, printFlag = FALSE, seed = 20261005)
imps <- imputationList(lapply(1:20, function(i) complete(imp, i)))
des_mi <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WT, nest = TRUE, data = imps)
fmi <- with(des_mi, svyglm(as.formula(paste("stroke ~ amd_any +", paste(full_cov, collapse = "+"))), family = quasibinomial()))
mi <- summary(MIcombine(fmi))["amd_any", ]
sens$S7 <- data.frame(model = "S7 Multiple imputation (m=20), Model 4", term = "amd_any", n = nrow(dm_in),
                      OR = exp(mi$results), lo = exp(mi$`(lower`), hi = exp(mi$`upper)`),
                      p = 2 * pnorm(-abs(mi$results / mi$se)))
sens <- bind_rows(sens)
write.table(sens, file.path(OUT, "sensitivity.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n== Sensitivity\n"); print(sens, digits = 3)

# ---------------------------------------------------------------- 7. Subgroups + interaction (core Model 3)
subs <- list(age_grp = "age", sex = "sex", race2 = "race", smoke2 = "smoke", obese = "bmi", htn_f = NULL, dm_f = NULL)
sg <- list()
for (g in names(subs)) {
  covs <- setdiff(core_cov, subs[[g]])
  dz <- subset(des, !is.na(des$variables[[g]]))
  lv <- levels(droplevels(dz$variables[[g]]))
  for (l in lv) {
    r <- tryCatch(fit_or(subset(dz, dz$variables[[g]] == l), "stroke", "amd_any", covs, g),
                  error = function(e) data.frame(model = g, term = "amd_any", n = NA, OR = NA, lo = NA, hi = NA, p = NA))
    r$level <- l
    tb <- table(subset(dz, dz$variables[[g]] == l)$variables[, c("amd_any", "stroke")])
    r$events_amd <- tb["1", "1"]; r$n_amd <- sum(tb["1", ])
    sg[[paste(g, l)]] <- r
  }
  fi <- svyglm(as.formula(paste("stroke ~ amd_any *", g, "+", paste(covs, collapse = "+"))), design = dz, family = quasibinomial())
  sg[[paste(g, l)]]$p_interaction <- NA
  sg[[paste(g, lv[1])]]$p_interaction <- regTermTest(fi, as.formula(paste0("~amd_any:", g)))$p
}
sg <- bind_rows(sg)
write.table(sg, file.path(OUT, "subgroups.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n== Subgroups\n"); print(sg, digits = 3)
saveRDS(d, file.path(ROOT, "derived/nhanes_stroke_amd_analytic.rds"))
cat("\nDone.\n")
