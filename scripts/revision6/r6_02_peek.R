# ---- project root: set the environment variable COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT <- local({
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(a)) dirname(normalizePath(a[1])) else getwd()
  while (basename(d) != "scripts" && dirname(d) != d) d <- dirname(d)
  gsub("\\\\", "/", Sys.getenv("COMORBIDITY_ROOT", unset = normalizePath(if (basename(d) == "scripts") dirname(d) else "..", winslash = "/", mustWork = FALSE)))
})
Sys.setenv(COMORBIDITY_ROOT = PROJ_ROOT)   # inherited by parallel (PSOCK) workers
suppressMessages(library(data.table))
for (f in c(file.path(PROJ_ROOT, "derived/gwas/IS.rds"),file.path(PROJ_ROOT, "derived/revision/MEGA_IS.rds"),file.path(PROJ_ROOT, "derived/gwas/AS.rds"))) {
  x <- readRDS(f); cat(f, nrow(x), "\n"); print(head(x, 3)) }
