# ============================================================================
# 18_extended_range.R
#
# Simulation study of funresMech 1.2.0 outside the base scenario (S1 Text,
# section 7). Each block changes one factor of the base scenario (h/T = 0.025,
# s/h = 0.5, E = 20 expected encounters per trial at the geometric mean
# density) and keeps the large design (8 densities, 64 trials):
#
#   s/h = 1.5   highly variable handling times
#   h/T = 0.1   long handling, close to saturation at the highest densities
#   E = 5       few encounters per trial
#
# Each block crosses z {1, 1.75, 2.75} with k {0.2, 1, 5}; 8 datasets per
# scenario; exploratory preset followed by verify_fit() (216 fits). The
# datasets are the same as in the earlier run of these scenarios with the
# pre-release code, so the two runs can be compared pair by pair.
#
# Usage, from the root of the package repository with funresMech installed:
#   QUICK <- TRUE;  source("validation/R/18_extended_range.R")  # check, ~15 min
#   rm(list = ls()); source("validation/R/18_extended_range.R")   # complete run
# The complete run needs about 200 core-hours (about 9 h on 23 cores). It can
# be interrupted and resumed. Results: validation/tables/extended_range_fits.csv,
# extended_range_checks.csv and extended_range_summary.csv.
# ============================================================================

library(funresMech)
stopifnot(packageVersion("funresMech") >= "1.1.0.9000")
source(file.path("validation", "R", "helpers.R"))

if (!exists("QUICK")) QUICK <- FALSE
WORKERS <- if (QUICK) 2L else max(1L, parallel::detectCores() - 1L)
DIR_RES <- file.path("validation", "results", if (QUICK) "extended_range_quick" else "extended_range")
DIR_TAB <- file.path("validation", "tables")
dir.create(DIR_RES, recursive = TRUE, showWarnings = FALSE); dir.create(DIR_TAB, showWarnings = FALSE)
SUFFIX <- if (QUICK) "_quick" else ""
GRID <- if (QUICK) c(0.5, 1, 2, 3) else seq(0.5, 5, by = 0.5)
N_PRECISE <- if (QUICK) 2000L else 20000L

control_for <- function(preset, ...) {
  if (QUICK) effort_control(preset, n_sim = 200, itermax = 5, NP = 12, n_reeval = 2, ...)
  else effort_control(preset, ...)
}

Z <- c(1, 1.75, 2.75); K <- c(0.2, 1, 5)
BLOCKS <- list("s/h = 1.5" = list(s_over_h = 1.5), "h/T = 0.1" = list(h_over_T = 0.1),
               "E = 5" = list(E = 5))
SCEN <- do.call(rbind, lapply(names(BLOCKS), function(b)
  do.call(rbind, apply(expand.grid(z = Z, k = K), 1, function(r)
    do.call(scenario, c(list(z = r[["z"]], k = r[["k"]], block = b), BLOCKS[[b]]))))))
if (QUICK) SCEN <- SCEN[c(5, 14, 23), ]
SCEN$id <- scenario_id(SCEN)

jobs <- do.call(rbind, lapply(SCEN$id, function(id)
  data.frame(arm = "extended_range", scenario = id, dataset = seq_len(if (QUICK) 1 else 8),
             variant = "exploratory")))
jobs$file <- file.path(DIR_RES, sprintf("%s__%s__d%d.rds", jobs$arm, jobs$scenario, jobs$dataset))
pending <- jobs[!file.exists(jobs$file), ]
cat(sprintf("Fits: %d in total, %d pending; %d workers.\n", nrow(jobs), nrow(pending), WORKERS))

run_one <- function(job) {
  library(funresMech)
  future::plan(future::sequential)
  t0 <- proc.time()[["elapsed"]]
  e <- SCEN[SCEN$id == job$scenario, ]
  d <- simulate_dataset(e, job$dataset); p0 <- true_parameters(e); m <- funres_model()
  set.seed(fit_seed("18", job$arm, job$scenario, job$dataset, job$variant))
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
                          "control_for", "EFFORT", "effort_control", "true_parameters", "simulate_dataset", "generating_model",
                          "data_seed", "fit_seed", "SEEDS_DATA", "SEEDS_FIT", "precise_nll", "lighten",
                          "fit_row", "checks_long")),
         workers = WORKERS, log_file = file.path(DIR_RES, "errors.log"))

# ---- tables ------------------------------------------------------------------
O <- lapply(jobs$file[file.exists(jobs$file)], readRDS)
if (length(O)) {
  G <- do.call(rbind, lapply(O, function(o)
    cbind(o$job, o$res, seconds = o$seconds,
          SCEN[SCEN$id == o$job$scenario, c("z", "k", "s_over_h", "h_over_T", "E", "design", "block")])))
  G$excludes_1 <- excludes_one(G$z_low, G$z_high)
  G$deficit <- G$nll_precise_hat - G$nll_precise_true
  G$optimization_failure <- G$deficit > 1
  G$open <- G$low_open | G$high_open
  G$z_error <- G$z_hat - G$z
  G$h_ratio <- G$h_hat / (G$h_over_T * T_EXP)
  G$k_ratio <- G$k_hat / G$k
  write.csv(G, file.path(DIR_TAB, paste0("extended_range_fits", SUFFIX, ".csv")), row.names = FALSE)

  C <- do.call(rbind, lapply(O, function(o) {
    k <- merge(o$job, o$checks)
    k$covers <- o$res$covers
    k$optimization_failure <- (o$res$nll_precise_hat - o$res$nll_precise_true) > 1
    k }))
  write.csv(C, file.path(DIR_TAB, paste0("extended_range_checks", SUFFIX, ".csv")), row.names = FALSE)

  summarise <- function(g, label) {
    w <- wilson(g$covers)
    data.frame(group = label, n = nrow(g), coverage = w[["p"]], coverage_low = w[["lower"]],
               coverage_high = w[["upper"]],
               type_I_error = if (any(g$z == 1)) mean(g$excludes_1[g$z == 1]) else NA,
               power = if (any(g$z > 1)) mean(g$excludes_1[g$z > 1]) else NA,
               open_interval = mean(g$open), optimization_failure = mean(g$optimization_failure),
               median_z_error = stats::median(g$z_error), median_ci_width = stats::median(g$z_high - g$z_low,
                                                                                         na.rm = TRUE),
               median_h_ratio = stats::median(g$h_ratio), median_k_ratio = stats::median(g$k_ratio),
               median_hours = stats::median(g$seconds) / 3600)
  }
  S <- rbind(summarise(G, "all"),
             do.call(rbind, lapply(split(G, G$block), function(g) summarise(g, g$block[1]))),
             do.call(rbind, lapply(split(G, list(G$block, G$z), drop = TRUE), function(g)
               summarise(g, sprintf("%s, z = %g", g$block[1], g$z[1])))),
             do.call(rbind, lapply(split(G, list(G$block, G$k), drop = TRUE), function(g)
               summarise(g, sprintf("%s, k = %g", g$block[1], g$k[1])))))
  write.csv(S, file.path(DIR_TAB, paste0("extended_range_summary", SUFFIX, ".csv")), row.names = FALSE)
  for (i in seq_len(1 + length(unique(G$block))))
    cat(sprintf("Coverage, %s: %.1f%% [%.1f, %.1f] (n = %d)\n", S$group[i], 100 * S$coverage[i],
                100 * S$coverage_low[i], 100 * S$coverage_high[i], S$n[i]))
}
