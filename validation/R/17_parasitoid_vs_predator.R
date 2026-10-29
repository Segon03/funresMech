# ============================================================================
# 17_parasitoid_vs_predator.R
#
# Comparison of the parasitoid process (attacked hosts remain available and
# can be found again) and the predator process (attacked hosts are removed;
# for a parasitoid, perfect and cost-free discrimination of parasitized hosts)
# in D2, D2 without the atypical trial, and D3 (S1 Text, section 11; S2 Fig).
# Exploratory preset, two seeds per combination (12 fits): the direction of
# the AIC difference is robust, its magnitude is approximate.
#
# Usage, from the root of the package repository with funresMech installed:
#   QUICK <- TRUE; source("validation/R/17_parasitoid_vs_predator.R")   # ~5 min
#   rm(list = ls()); source("validation/R/17_parasitoid_vs_predator.R")  # ~2 h
# Results: validation/tables/process_comparison_fits.csv,
#          validation/tables/process_comparison_profiles.csv
# ============================================================================

library(funresMech)
stopifnot(packageVersion("funresMech") >= "1.1.0.9000")
source(file.path("validation", "R", "helpers.R"))
if (!exists("QUICK")) QUICK <- FALSE
DIR_RES <- file.path("validation", "results", if (QUICK) "process_comparison_quick" else "process_comparison")
DIR_TAB <- file.path("validation", "tables")
dir.create(DIR_RES, recursive = TRUE, showWarnings = FALSE); dir.create(DIR_TAB, showWarnings = FALSE)
SUFFIX <- if (QUICK) "_quick" else ""
GRID <- if (QUICK) c(1, 2, 3) else seq(0.5, 5, by = 0.5)

jobs <- expand.grid(data = c("D2", "D2_without_atypical", "D3"), process = c("parasitoid", "predator"),
                    seed_run = 1:2, stringsAsFactors = FALSE)
if (QUICK) jobs <- jobs[jobs$data == "D3" & jobs$seed_run == 1, ]
jobs$file <- file.path(DIR_RES, sprintf("%s__%s__s%d.rds", jobs$data, jobs$process, jobs$seed_run))
pending <- jobs[!file.exists(jobs$file), ]
cat(sprintf("Fits: %d in total, %d pending.\n", nrow(jobs), nrow(pending)))

run_one <- function(job) {
  library(funresMech)
  future::plan(future::sequential)
  d <- read_published(job$data)
  ctrl <- if (QUICK) funres_control(n_sim = 300, itermax = 8, NP = 16, n_reeval = 2, bb = FALSE)
          else funres_control(n_sim = 2000, itermax = 50, n_reeval = 3, bb = FALSE)
  set.seed(fit_seed("17", "process_comparison", job$data, job$seed_run, job$process))
  t0 <- proc.time()[["elapsed"]]
  f <- suppressWarnings(fit_funres(d, T = T_EXP, model = funres_model(process = job$process),
                                   z_grid = GRID, control = ctrl, extend_z = FALSE, time_unit = "h"))
  saveRDS(list(job = job, fit = f, seconds = proc.time()[["elapsed"]] - t0), job$file)
  TRUE
}
run_jobs(pending, run_one,
         globals = mget(c("run_one", "pending", "read_published", "fit_seed", "SEEDS_FIT", "T_EXP", "GRID", "QUICK")),
         workers = min(max(1L, nrow(pending)), max(1L, parallel::detectCores() - 1L)),
         log_file = file.path(DIR_RES, "errors.log"))

# ---- tables ----------------------------------------------------------------------
O <- lapply(jobs$file[file.exists(jobs$file)], readRDS)
R <- do.call(rbind, lapply(O, function(o) {
  f <- o$fit; p <- f$par
  data.frame(data = o$job$data, process = o$job$process, seed_run = o$job$seed_run,
             trials = nrow(f$data), z_hat = p[["z"]], z_low = f$ci$z_low, z_high = f$ci$z_high,
             a = p[["a"]], h = p[["h"]], k = p[["k"]], s = p[["s"]], nll = f$nll, nll_se = f$nll_se,
             AIC = f$aic$AIC_full, reliability = attr(reliability(f), "overall"), hours = o$seconds / 3600)
}))
write.csv(R, file.path(DIR_TAB, paste0("process_comparison_fits", SUFFIX, ".csv")), row.names = FALSE)
P <- do.call(rbind, lapply(O, function(o)
  data.frame(data = o$job$data, process = o$job$process, seed_run = o$job$seed_run,
             o$fit$profile[, c("z", "nll", "nll_se", "a", "h", "k", "s", "n_floor")])))
write.csv(P, file.path(DIR_TAB, paste0("process_comparison_profiles", SUFFIX, ".csv")), row.names = FALSE)

# DeltaAIC = AIC(predator) - AIC(parasitoid): positive values favor the parasitoid process
w <- merge(R[R$process == "parasitoid", c("data", "seed_run", "AIC", "z_hat", "h", "k")],
           R[R$process == "predator", c("data", "seed_run", "AIC", "z_hat", "h", "k")],
           by = c("data", "seed_run"), suffixes = c("_parasitoid", "_predator"))
w$dAIC_predator_minus_parasitoid <- w$AIC_predator - w$AIC_parasitoid
print(w, digits = 3, row.names = FALSE)
