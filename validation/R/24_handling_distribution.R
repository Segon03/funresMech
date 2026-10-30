# ============================================================================
# 24_handling_distribution.R
#
# Which family should be the default for handling times? Datasets are
# simulated with lognormal, Gamma or normal (truncated at zero) handling
# times and each one is fitted twice, with lognormal and with Gamma handling
# (the two built-in families that are candidates for the default). For each
# fit: recovery of z, h and the SD of handling times (s), coverage of the
# interval of z, and the likelihood, so that the two families can be compared
# by AIC on the same data (both have five free parameters).
#
# Handling times matter most when they take a large share of the trial, so
# all scenarios use h/T = 0.1 (the long-handling block of 18_extended_range.R).
#   process     parasitoid, predator
#   generator   lognormal, Gamma, truncated normal
#   s/h         0.25 (nearly symmetric for the three families), 0.6 (skewed)
#   z           1, 1.75
#   k = 1, E = 20, large design (8 densities, 64 trials); 6 datasets each.
# 24 scenarios, 144 datasets, 288 fits. The truncated normal cannot have a
# CV above 0.76, which is why the larger s/h is 0.6 and not 1.
#
# Usage, from the root of the package repository with funresMech installed:
#   QUICK <- TRUE; source("validation/R/24_handling_distribution.R")  # check, a few minutes
#   rm(list = ls()); source("validation/R/24_handling_distribution.R")
# Approximate cost: 100 core-hours (about 4.5 h on 23 cores). The run can be
# interrupted and resumed. Results:
#   validation/tables/handling_fits.csv       one row per fit
#   validation/tables/handling_summary.csv    recovery and coverage by group
#   validation/tables/handling_selection.csv  Gamma vs lognormal by AIC
# Optional: FIT_TNORM <- TRUE also fits the truncated normal (an R sampler,
# several times slower).
# ============================================================================

library(funresMech)
stopifnot(packageVersion("funresMech") >= "1.1.0.9000")
source(file.path("validation", "R", "helpers.R"))

if (!exists("QUICK")) QUICK <- FALSE
if (!exists("FIT_TNORM")) FIT_TNORM <- FALSE
WORKERS <- if (QUICK) 2L else max(1L, parallel::detectCores() - 1L)
DIR_RES <- file.path("validation", "results", if (QUICK) "handling_quick" else "handling")
DIR_TAB <- file.path("validation", "tables")
dir.create(DIR_RES, recursive = TRUE, showWarnings = FALSE); dir.create(DIR_TAB, showWarnings = FALSE)
SUFFIX <- if (QUICK) "_quick" else ""
GRID <- if (QUICK) c(0.5, 1, 2, 3) else seq(0.5, 5, by = 0.5)
N_PRECISE <- if (QUICK) 2000L else 20000L
N_DATA <- if (QUICK) 1L else 6L
FAMILIES <- c("lognormal", "gamma", if (FIT_TNORM) "tnorm")

control_for <- function(level, ...) {
  if (QUICK) effort_control(level, n_sim = 200, itermax = 5, NP = 12, n_reeval = 2, ...)
  else effort_control(level, ...)
}

# ---- scenarios -------------------------------------------------------------------
GEN <- data.frame(process = rep(c("parasitoid", "predator"), each = 3),
                  generator = rep(c("lognormal", "gamma", "tnorm"), 2),
                  generating_model = c("correct", "handling_gamma", "handling_tnorm",
                                       "predator", "predator_handling_gamma", "predator_handling_tnorm"),
                  stringsAsFactors = FALSE)
SCEN <- do.call(rbind, lapply(seq_len(nrow(GEN)), function(i)
  do.call(rbind, apply(expand.grid(z = c(1, 1.75), s_over_h = c(0.25, 0.6)), 1, function(r)
    cbind(scenario(z = r[["z"]], k = 1, s_over_h = r[["s_over_h"]], h_over_T = 0.1,
                   generating_model = GEN$generating_model[i], block = "handling"),
          process = GEN$process[i], generator = GEN$generator[i], stringsAsFactors = FALSE)))))
SCEN$id <- scenario_id(SCEN)
stopifnot(!anyDuplicated(SCEN$id))

# Seeds: added once to the seed tables, then always read from them. Scenarios
# already in the table (lognormal generators) keep their dataset seeds.
add_seeds <- function() {
  sd <- utils::read.csv(file.path("validation", "data", "seeds_datasets.csv"), stringsAsFactors = FALSE)
  sf <- utils::read.csv(file.path("validation", "data", "seeds_fits.csv"), stringsAsFactors = FALSE)
  new_d <- setdiff(SCEN$id, sd$scenario)
  if (length(new_d)) {
    base <- vapply(new_d, function(i) sum(utf8ToInt(i)) + 30000, 0)
    while (anyDuplicated(base) || any(base %in% sd$seed_base)) {
      j <- which(duplicated(base) | base %in% sd$seed_base)[1]; base[j] <- base[j] + 1
    }
    sd <- rbind(sd, data.frame(scenario = new_d, seed_base = unname(base)))
  }
  want <- expand.grid(scenario = SCEN$id, dataset = 1:6, variant = paste0("fit_", c("lognormal", "gamma", "tnorm")),
                      stringsAsFactors = FALSE)
  have <- paste(sf$scenario, sf$dataset, sf$variant)[sf$script == "24"]
  want <- want[!paste(want$scenario, want$dataset, want$variant) %in% have, ]
  if (nrow(want))
    sf <- rbind(sf, data.frame(script = "24", arm = "handling", scenario = want$scenario, dataset = want$dataset,
                               variant = want$variant,
                               seed = vapply(seq_len(nrow(want)), function(i)
                                 sum(utf8ToInt(paste("24", want$scenario[i], want$variant[i]))) * 7 +
                                   want$dataset[i] * 101, 0)))
  if (!length(new_d) && !nrow(want)) return(invisible(FALSE))
  utils::write.csv(sd, file.path("validation", "data", "seeds_datasets.csv"), row.names = FALSE)
  utils::write.csv(sf, file.path("validation", "data", "seeds_fits.csv"), row.names = FALSE)
  invisible(TRUE)
}
if (add_seeds()) {
  SEEDS_DATA <- utils::read.csv(file.path("validation", "data", "seeds_datasets.csv"), stringsAsFactors = FALSE)
  SEEDS_FIT  <- utils::read.csv(file.path("validation", "data", "seeds_fits.csv"), stringsAsFactors = FALSE)
}
if (QUICK) SCEN <- SCEN[c(1, 6, 11, 16, 21), ]

jobs <- expand.grid(scenario = SCEN$id, dataset = seq_len(N_DATA), family = FAMILIES, stringsAsFactors = FALSE)
jobs$arm <- "handling"; jobs$variant <- paste0("fit_", jobs$family)
jobs$file <- file.path(DIR_RES, sprintf("%s_d%d_%s.rds", jobs$scenario, jobs$dataset, jobs$family))
pending <- jobs[!file.exists(jobs$file), ]
cat(sprintf("Handling study: %d fits in total, %d pending; %d workers.\n", nrow(jobs), nrow(pending), WORKERS))

fit_model <- function(family, process)
  funres_model(handling = if (family == "tnorm") tnorm_sampler() else family, process = process)

run_one <- function(job) {
  library(funresMech)
  future::plan(future::sequential)
  t0 <- proc.time()[["elapsed"]]
  e <- SCEN[SCEN$id == job$scenario, ]
  d <- simulate_dataset(e, job$dataset); p0 <- true_parameters(e)
  m <- fit_model(job$family, e$process)
  set.seed(fit_seed("24", job$arm, job$scenario, job$dataset, job$variant))
  f <- suppressWarnings(fit_funres(d, T_EXP, m, z_grid = GRID, control = control_for("exploratory", bb = TRUE)))
  f <- verify_fit(f)
  res <- fit_row(f, e$z)
  res$nll_precise_hat <- precise_nll(f$par, m, d)
  # likelihood of the true parameters under the generating model (same for both fits of a dataset)
  res$nll_precise_true <- precise_nll(p0, generating_model(e), d)
  out <- list(job = job[, c("arm", "scenario", "dataset", "family", "variant")], res = res,
              checks = checks_long(f), fit = lighten(f), seconds = proc.time()[["elapsed"]] - t0)
  saveRDS(out, job$file)
  TRUE
}

run_jobs(pending, run_one,
         globals = mget(c("run_one", "pending", "SCEN", "DESIGNS", "T_EXP", "GRID", "N_PRECISE", "QUICK",
                          "control_for", "EFFORT", "effort_control", "true_parameters", "simulate_dataset",
                          "generating_model", "fit_model", "tnorm_sampler", "tnorm_mean", "tnorm_cv", "tnorm_draw",
                          "data_seed", "fit_seed", "SEEDS_DATA", "SEEDS_FIT", "precise_nll", "lighten",
                          "fit_row", "excludes_one", "checks_long")),
         workers = WORKERS, log_file = file.path(DIR_RES, "errors.log"))

# ---- tables ------------------------------------------------------------------
O <- lapply(jobs$file[file.exists(jobs$file)], readRDS)
if (length(O)) {
  info <- c("z", "k", "s_over_h", "h_over_T", "E", "design", "process", "generator", "generating_model")
  G <- do.call(rbind, lapply(O, function(o)
    data.frame(o$job, o$res, seconds = o$seconds, SCEN[SCEN$id == o$job$scenario, info], row.names = NULL)))
  G$excludes_1 <- excludes_one(G$z_low, G$z_high)
  G$deficit <- G$nll_precise_hat - G$nll_precise_true
  G$optimization_failure <- G$deficit > 1
  G$open <- G$low_open | G$high_open
  G$z_error <- G$z_hat - G$z
  G$h_ratio <- G$h_hat / (G$h_over_T * T_EXP)
  # s is the SD of handling times for the built-in families (a shape for the truncated normal)
  G$s_ratio <- ifelse(G$family == "tnorm", NA_real_, G$s_hat / (G$s_over_h * G$h_over_T * T_EXP))
  G$AIC <- 2 * G$nll_precise_hat + 2 * 5
  utils::write.csv(G, file.path(DIR_TAB, paste0("handling_fits", SUFFIX, ".csv")), row.names = FALSE)

  summarise <- function(g) {
    w <- wilson(g$covers)
    data.frame(n = nrow(g), coverage = w[["p"]], coverage_low = w[["lower"]], coverage_high = w[["upper"]],
               type_I_error = if (any(g$z == 1)) mean(g$excludes_1[g$z == 1]) else NA,
               power = if (any(g$z > 1)) mean(g$excludes_1[g$z > 1]) else NA,
               median_z_error = stats::median(g$z_error), median_h_ratio = stats::median(g$h_ratio),
               median_s_ratio = stats::median(g$s_ratio), open_interval = mean(g$open),
               optimization_failure = mean(g$optimization_failure))
  }
  S <- do.call(rbind, lapply(split(G, list(G$process, G$generator, G$s_over_h, G$family), drop = TRUE), function(g)
    cbind(g[1, c("process", "generator", "s_over_h", "family")], summarise(g), row.names = NULL)))
  S <- rbind(S, do.call(rbind, lapply(split(G, list(G$generator, G$family), drop = TRUE), function(g)
    cbind(process = "both", g[1, "generator", drop = FALSE], s_over_h = NA, g[1, "family", drop = FALSE],
          summarise(g), row.names = NULL))))
  utils::write.csv(S, file.path(DIR_TAB, paste0("handling_summary", SUFFIX, ".csv")), row.names = FALSE)

  # Gamma vs lognormal on the same dataset: dAIC = AIC(lognormal) - AIC(Gamma) (> 0 favors Gamma)
  W <- reshape(G[G$family %in% c("lognormal", "gamma"), c("scenario", "dataset", "process", "generator", "s_over_h",
                                                           "z", "family", "AIC")],
               idvar = c("scenario", "dataset", "process", "generator", "s_over_h", "z"), timevar = "family",
               direction = "wide")
  if (all(c("AIC.lognormal", "AIC.gamma") %in% names(W))) {
    W$dAIC <- W$AIC.lognormal - W$AIC.gamma
    sel <- function(g) data.frame(n = nrow(g), median_dAIC = stats::median(g$dAIC, na.rm = TRUE),
                                  gamma_better = mean(g$dAIC > 0, na.rm = TRUE),
                                  gamma_better_by_2 = mean(g$dAIC > 2, na.rm = TRUE),
                                  lognormal_better_by_2 = mean(g$dAIC < -2, na.rm = TRUE))
    P <- rbind(
      do.call(rbind, lapply(split(W, list(W$process, W$generator, W$s_over_h), drop = TRUE), function(g)
        cbind(g[1, c("process", "generator", "s_over_h")], sel(g), row.names = NULL))),
      do.call(rbind, lapply(split(W, W$generator), function(g)
        cbind(process = "both", generator = g$generator[1], s_over_h = NA, sel(g), row.names = NULL))))
    utils::write.csv(P, file.path(DIR_TAB, paste0("handling_selection", SUFFIX, ".csv")), row.names = FALSE)
    print(P[P$process == "both", ], digits = 3)
  }
  print(S[S$process == "both", c("generator", "family", "n", "coverage", "median_z_error", "median_h_ratio",
                                  "median_s_ratio")], digits = 3)
}
