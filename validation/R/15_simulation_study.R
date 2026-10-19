# ============================================================================
# 15_simulation_study.R
#
# Simulation study of funresMech 1.2.0 with known parameters (main text,
# Table 2 and Fig 2; S1 Text, section 7).
#
#   core_exploratory  z {1, 1.75, 2.75} x k {0.2, 1, 5} x design {large: 8
#                     densities, 64 trials; small: 5 densities, 30 trials};
#                     8 datasets per scenario; exploratory preset (144 fits).
#   core_standard     The same scenarios, datasets 1 to 4, standard preset
#                     (72 fits on the same data).
#   dense_design      z x k with 13 densities and 5 trials each (65 trials);
#                     8 datasets per scenario; exploratory preset (72 fits).
#
# Every fit uses a z grid from 0.5 to 5 in steps of 0.5 (with refinement and
# extension) and is followed by verify_fit(). The NLL at the estimate and at
# the true parameters is evaluated precisely (5 x 20,000 simulated trials per
# density); a fit with a precise NLL more than one unit above that of the
# true parameters is an optimization failure.
#
# Usage, from the root of the package repository with funresMech installed:
#   QUICK <- TRUE;  source("validation/R/15_simulation_study.R")  # check, ~15 min
#   rm(list = ls()); source("validation/R/15_simulation_study.R")   # complete run
# The complete run took about one day on 23 cores. It can be interrupted and
# resumed. Results: validation/tables/simulation_study_fits.csv and
# validation/tables/simulation_study_checks.csv.
# ============================================================================

library(funresMech)
stopifnot(packageVersion("funresMech") >= "1.2.0")
source(file.path("validation", "R", "helpers.R"))

if (!exists("QUICK")) QUICK <- FALSE
WORKERS <- if (QUICK) 2L else max(1L, parallel::detectCores() - 1L)
DIR_RES <- file.path("validation", "results", if (QUICK) "simulation_study_quick" else "simulation_study")
DIR_TAB <- file.path("validation", "tables")
dir.create(DIR_RES, recursive = TRUE, showWarnings = FALSE); dir.create(DIR_TAB, showWarnings = FALSE)
SUFFIX <- if (QUICK) "_quick" else ""
GRID <- if (QUICK) c(0.5, 1, 2, 3) else seq(0.5, 5, by = 0.5)
N_PRECISE <- if (QUICK) 2000L else 20000L

control_for <- function(preset, ...) {
  if (QUICK) {
    p <- switch(preset, exploratory = list(n_sim = 200, itermax = 5), standard = list(n_sim = 300, itermax = 8))
    funres_control(preset, n_sim = p$n_sim, itermax = p$itermax, NP = 12, n_reeval = 2, ...)
  } else funres_control(preset, ...)
}

Z <- c(1, 1.75, 2.75); K <- c(0.2, 1, 5)
SCEN_CORE <- do.call(rbind, apply(expand.grid(z = Z, k = K, design = c("large", "small"),
                                              stringsAsFactors = FALSE), 1, function(r)
  scenario(z = as.numeric(r[["z"]]), k = as.numeric(r[["k"]]), design = r[["design"]], block = "core")))
SCEN_DENSE <- do.call(rbind, apply(expand.grid(z = Z, k = K), 1, function(r)
  scenario(z = r[["z"]], k = r[["k"]], design = "dense", block = "dense")))
if (QUICK) { SCEN_CORE <- SCEN_CORE[c(2, 14), ]; SCEN_DENSE <- SCEN_DENSE[5, ] }
SCEN_CORE$id <- scenario_id(SCEN_CORE); SCEN_DENSE$id <- scenario_id(SCEN_DENSE)
SCEN <- rbind(SCEN_CORE, SCEN_DENSE)

jobs <- rbind(
  do.call(rbind, lapply(SCEN_CORE$id, function(id) rbind(
    data.frame(arm = "core_exploratory", scenario = id, dataset = seq_len(if (QUICK) 1 else 8),
               variant = "exploratory"),
    data.frame(arm = "core_standard", scenario = id, dataset = seq_len(if (QUICK) 1 else 4),
               variant = "standard")))),
  do.call(rbind, lapply(SCEN_DENSE$id, function(id)
    data.frame(arm = "dense_design", scenario = id, dataset = seq_len(if (QUICK) 1 else 8),
               variant = "exploratory"))))
jobs$file <- file.path(DIR_RES, sprintf("%s__%s__d%d.rds", jobs$arm, jobs$scenario, jobs$dataset))
pending <- jobs[!file.exists(jobs$file), ]
cat(sprintf("Fits: %d in total, %d pending; %d workers.\n", nrow(jobs), nrow(pending), WORKERS))

run_one <- function(job) {
  library(funresMech)
  future::plan(future::sequential)
  t0 <- proc.time()[["elapsed"]]
  e <- SCEN[SCEN$id == job$scenario, ]
  d <- simulate_dataset(e, job$dataset); p0 <- true_parameters(e); m <- funres_model()
  set.seed(fit_seed("15", job$arm, job$scenario, job$dataset, job$variant))
  f <- suppressWarnings(fit_funres(d, T_EXP, m, z_grid = GRID, control = control_for(job$variant, bb = TRUE)))
  f <- verify_fit(f)
  res <- fit_row(f, e$z)
  res$nll_precise_hat <- precise_nll(f$par, m, d)
  res$nll_precise_true <- precise_nll(p0, m, d)
  out <- list(job = job[, c("arm", "scenario", "dataset", "variant")], res = res, checks = checks_long(f),
              fit = lighten(f), seconds = proc.time()[["elapsed"]] - t0)
  saveRDS(out, job$file)
  TRUE
}

run_jobs(pending, run_one,
         globals = mget(c("run_one", "pending", "SCEN", "DESIGNS", "T_EXP", "GRID", "N_PRECISE", "QUICK",
                          "control_for", "true_parameters", "simulate_dataset", "generating_model",
                          "data_seed", "fit_seed", "SEEDS_DATA", "SEEDS_FIT", "precise_nll", "lighten",
                          "fit_row", "checks_long")),
         workers = WORKERS, log_file = file.path(DIR_RES, "errors.log"))

# ---- tables ------------------------------------------------------------------
O <- lapply(jobs$file[file.exists(jobs$file)], readRDS)
if (length(O)) {
  G <- do.call(rbind, lapply(O, function(o)
    cbind(o$job, o$res, seconds = o$seconds,
          SCEN[SCEN$id == o$job$scenario, c("z", "k", "s_over_h", "h_over_T", "E", "design", "block")])))
  G$deficit <- G$nll_precise_hat - G$nll_precise_true
  G$optimization_failure <- G$deficit > 1
  G$open <- G$low_open | G$high_open
  G$z_error <- G$z_hat - G$z
  write.csv(G, file.path(DIR_TAB, paste0("simulation_study_fits", SUFFIX, ".csv")), row.names = FALSE)
  C <- do.call(rbind, lapply(O, function(o) {
    k <- merge(o$job, o$checks)
    k$covers <- o$res$covers
    k$optimization_failure <- (o$res$nll_precise_hat - o$res$nll_precise_true) > 1
    k$z_error <- o$res$z_hat - SCEN$z[SCEN$id == o$job$scenario]
    k }))
  write.csv(C, file.path(DIR_TAB, paste0("simulation_study_checks", SUFFIX, ".csv")), row.names = FALSE)
  for (a in unique(G$arm)) {
    w <- wilson(G$covers[G$arm == a])
    cat(sprintf("Coverage, %s: %.1f%% [%.1f, %.1f] (n = %d)\n", a, 100 * w[["p"]], 100 * w[["lower"]],
                100 * w[["upper"]], w[["n"]]))
  }
  cat("Summaries: source(\"validation/R/16_simulation_study_analysis.R\")\n")
}
