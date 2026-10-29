# ============================================================================
# 12_extended_validation.R
#
# Extended validation of funresMech 1.2.0 (S1 Text, sections 9 to 11).
#
#   misspecification  Data generated with lognormal or Weibull search times,
#                     Gamma handling times (s/h = 1), the predator process, or
#                     individual variation in a (lognormal, sigma = 0.5), and
#                     with the default model; k = 0.3, z = 1 or 2.25, 8
#                     densities with 64 trials; 8 datasets per combination.
#                     Each dataset is fitted with the default model and, when
#                     the package can represent it, with the generating model
#                     (exploratory preset; 160 fits).
#   budget            Datasets 1 to 4 of the nine core scenarios with 8
#                     densities of 15_simulation_study.R (36 datasets), fitted
#                     with three independent exploratory runs (the first one
#                     followed by verify_fit()), with the best of the three
#                     runs at each profile point, and with the standard preset.
#   D2_D3_final       D2, D2 without the atypical trial (2 of 10 hosts at
#                     density 10), and D3 with the final preset (grid from 0.5
#                     to 5 in steps of 0.25, with extension and refinement),
#                     followed by verify_fit().
#
# The misspecification and budget fits use a z grid from 0.5 to 5 in steps of
# 0.5. In the runs reported in S1 Text, they were made before the local
# refinement of the profile and the shape-preserving interpolation of the
# interval limits were added to the package; the results of a new run can
# therefore differ slightly in the interval metrics.
#
# Usage, from the root of the package repository with funresMech installed:
#   QUICK <- TRUE;  source("validation/R/12_extended_validation.R")   # check, ~30 min
#   rm(list = ls()); source("validation/R/12_extended_validation.R")   # complete run
# The complete run takes about two days on 23 cores; it can be interrupted
# and resumed. Results: validation/tables/misspecification_fits.csv,
# budget_fits.csv, D2_D3_final_preset_fits.csv, D2_D3_final_preset_profiles.csv.
# ============================================================================

library(funresMech)
stopifnot(packageVersion("funresMech") >= "1.1.0.9000")
source(file.path("validation", "R", "helpers.R"))

if (!exists("QUICK")) QUICK <- FALSE
WORKERS <- if (QUICK) 2L else max(1L, parallel::detectCores() - 1L)
DIR_RES <- file.path("validation", "results", if (QUICK) "extended_quick" else "extended")
DIR_TAB <- file.path("validation", "tables")
dir.create(DIR_RES, recursive = TRUE, showWarnings = FALSE); dir.create(DIR_TAB, showWarnings = FALSE)
SUFFIX <- if (QUICK) "_quick" else ""
GRID <- if (QUICK) c(0.5, 1, 2, 3) else seq(0.5, 5, by = 0.5)
GRID_FINAL <- if (QUICK) c(0.5, 1, 2, 3) else seq(0.5, 5, by = 0.25)
N_PRECISE <- if (QUICK) 2000L else 20000L
WORKERS_FINAL <- 2L      # each D2/D3 fit evaluates its profile points with 2 processes

control_for <- function(preset, ...) {
  if (QUICK) {
    p <- switch(preset, exploratory = list(n_sim = 200, itermax = 5),
                standard = list(n_sim = 300, itermax = 8), final = list(n_sim = 400, itermax = 10))
    effort_control(preset, n_sim = p$n_sim, itermax = p$itermax, NP = 12, n_reeval = 2, ...)
  } else effort_control(preset, ...)
}

# ---- scenarios -----------------------------------------------------------------
GEN <- c("correct", "search_weibull", "search_lognormal", "handling_gamma", "predator", "heterogeneity_a")
SCEN_MIS <- do.call(rbind, lapply(GEN, function(g) do.call(rbind, lapply(c(1, 2.25), function(z)
  scenario(z = z, k = 0.3, s_over_h = if (g == "handling_gamma") 1 else 0.5, generating_model = g,
           sigma_a = if (g == "heterogeneity_a") 0.5 else 0, block = "misspecification")))))
SCEN_BUD <- do.call(rbind, apply(expand.grid(z = c(1, 1.75, 2.75), k = c(0.2, 1, 5)), 1, function(r)
  scenario(z = r[["z"]], k = r[["k"]], design = "large", block = "core")))
if (QUICK) {
  SCEN_MIS <- SCEN_MIS[SCEN_MIS$generating_model %in% c("search_weibull", "predator") & SCEN_MIS$z == 2.25, ]
  SCEN_BUD <- SCEN_BUD[2, ]
}
SCEN_MIS$id <- scenario_id(SCEN_MIS); SCEN_BUD$id <- scenario_id(SCEN_BUD)
SCEN <- rbind(SCEN_MIS, SCEN_BUD)
N_DATA <- if (QUICK) 1L else 8L
N_BUDGET <- if (QUICK) 1L else 4L

jobs <- list()
add <- function(arm, scen, dataset, variant)
  jobs[[length(jobs) + 1]] <<- data.frame(arm = arm, scenario = scen, dataset = dataset, variant = variant,
                                          stringsAsFactors = FALSE)
for (i in seq_len(nrow(SCEN_MIS))) for (j in seq_len(N_DATA)) {
  add("misspecification", SCEN_MIS$id[i], j, "default_fit")
  if (!SCEN_MIS$generating_model[i] %in% c("correct", "heterogeneity_a"))
    add("misspecification", SCEN_MIS$id[i], j, "correct_fit")
}
for (i in seq_len(nrow(SCEN_BUD))) for (j in seq_len(N_BUDGET))
  for (v in c("exploratory", "exploratory_run2", "exploratory_run3", "standard")) add("budget", SCEN_BUD$id[i], j, v)
for (nm in c("D2", "D2_without_atypical", "D3")) add("D2_D3_final", nm, NA, "final")
jobs <- do.call(rbind, jobs)
jobs$file <- file.path(DIR_RES, sprintf("%s__%s__d%s__%s.rds", jobs$arm, jobs$scenario, jobs$dataset, jobs$variant))
jobs <- jobs[order(jobs$arm != "D2_D3_final", jobs$variant != "standard"), ]     # long jobs first
pending <- jobs[!file.exists(jobs$file), ]
cat(sprintf("Fits: %d in total, %d pending; %d workers.\n", nrow(jobs), nrow(pending), WORKERS))

run_one <- function(job) {
  library(funresMech)
  future::plan(future::sequential)
  t0 <- proc.time()[["elapsed"]]
  seed <- fit_seed("12", job$arm, job$scenario, job$dataset, job$variant)
  if (job$arm == "D2_D3_final") {
    options(parallelly.maxWorkers.localhost = Inf)
    future::plan(future::multisession, workers = WORKERS_FINAL)
    on.exit(future::plan(future::sequential), add = TRUE)
    d <- read_published(sub("_without_atypical", "", job$scenario))
    scr <- screen_outliers(d, T_EXP)
    removed <- integer(0)
    if (grepl("_without_atypical", job$scenario)) {
      removed <- scr$row[scr$flagged]
      stopifnot(identical(removed, 13L), d$dens[13] == 10, d$par[13] == 2)
      d <- d[-removed, ]
    }
    checkpoint <- sub("\\.rds$", "_fit.rds", job$file)
    if (file.exists(checkpoint)) f <- readRDS(checkpoint) else {
      set.seed(seed)
      f <- fit_funres(d, T_EXP, z_grid = GRID_FINAL, control = control_for("final"))
      saveRDS(f, checkpoint)                   # checkpoint before verify_fit()
    }
    f <- verify_fit(f)
    res <- cbind(data.frame(data = job$scenario, trials = nrow(d), rows_removed = paste(removed, collapse = ",")),
                 fit_row(f), f$aic[, c("AIC_full", "AIC_restricted", "AIC_BB")])
    out <- list(res = res, checks = checks_long(f), fit = lighten(f), screening = scr)
  } else {
    e <- SCEN[SCEN$id == job$scenario, ]
    d <- simulate_dataset(e, job$dataset); p0 <- true_parameters(e)
    m <- if (job$variant == "correct_fit") generating_model(e) else funres_model()
    preset <- if (job$variant == "standard") "standard" else "exploratory"
    set.seed(seed)
    f <- suppressWarnings(fit_funres(d, T_EXP, m, z_grid = GRID,
                                     control = control_for(preset, bb = job$arm == "misspecification" ||
                                                             job$variant == "exploratory")))
    if (job$arm == "budget" && job$variant == "exploratory") f <- verify_fit(f)
    res <- fit_row(f, e$z)
    if (job$arm == "budget") {
      res$nll_precise_hat <- precise_nll(f$par, m, d)
      res$nll_precise_true <- precise_nll(p0, m, d)
    }
    out <- list(res = res, checks = checks_long(f), fit = lighten(f))
  }
  out$seconds <- proc.time()[["elapsed"]] - t0
  out$job <- job[, c("arm", "scenario", "dataset", "variant")]
  saveRDS(out, job$file)
  TRUE
}

n_final <- sum(pending$arm == "D2_D3_final")
run_jobs(pending, run_one,
         globals = mget(c("run_one", "pending", "SCEN", "DESIGNS", "T_EXP", "GRID", "GRID_FINAL", "N_PRECISE",
                          "WORKERS_FINAL", "QUICK", "control_for", "EFFORT", "effort_control", "true_parameters", "simulate_dataset",
                          "generating_model", "data_seed", "fit_seed", "SEEDS_DATA", "SEEDS_FIT", "read_published",
                          "precise_nll", "lighten", "fit_row", "checks_long")),
         workers = max(1L, WORKERS - n_final * (WORKERS_FINAL - 1L)),
         log_file = file.path(DIR_RES, "errors.log"))

# ---- tables ------------------------------------------------------------------------
load_arm <- function(arm) lapply(jobs$file[jobs$arm == arm & file.exists(jobs$file)], readRDS)
row_of <- function(o) cbind(o$job, o$res, seconds = o$seconds,
                            SCEN[SCEN$id == o$job$scenario, c("z", "k", "s_over_h", "h_over_T", "E", "design",
                                                               "generating_model", "block")])
write_table <- function(x, name) if (!is.null(x) && nrow(x))
  write.csv(x, file.path(DIR_TAB, paste0(name, SUFFIX, ".csv")), row.names = FALSE)

OM <- load_arm("misspecification")
if (length(OM)) write_table(do.call(rbind, lapply(OM, row_of)), "misspecification_fits")

OB <- load_arm("budget")
if (length(OB)) {
  B <- do.call(rbind, lapply(OB, row_of))
  # Best of three: at each profile point, the best of the three exploratory runs,
  # re-evaluated with new simulations (no selection bias).
  b3 <- list()
  for (key in unique(paste(B$scenario, B$dataset))) {
    oo <- OB[vapply(OB, function(o) paste(o$job$scenario, o$job$dataset) == key &&
                      grepl("exploratory", o$job$variant), logical(1))]
    if (length(oo) < 3) next
    e <- SCEN[SCEN$id == oo[[1]]$job$scenario, ]; d <- oo[[1]]$fit$data
    P <- do.call(rbind, lapply(oo, function(o) o$fit$profile))
    P <- do.call(rbind, lapply(split(P, round(P$z, 8)), function(q) q[which.min(q$nll), ]))
    ctl <- oo[[1]]$fit$control
    set.seed(fit_seed("12", "budget", e$id, oo[[1]]$job$dataset, "best_of_3_reevaluation"))
    P$nll <- vapply(seq_len(nrow(P)), function(i) mean(vapply(seq_len(ctl$n_reeval), function(r)
      funresMech:::.nll_sim(c(a = P$a[i], h = P$h[i], z = P$z[i], k = P$k[i], s = P$s[i]),
                            funres_model(), d, T_EXP, ctl$n_sim), 0)), 0)
    P <- P[order(P$z), ]
    ci <- funresMech:::.profile_ci(P$z, P$nll)
    ib <- which.min(P$nll)
    par <- c(a = P$a[ib], h = P$h[ib], z = P$z[ib], k = P$k[ib], s = P$s[ib])
    b3[[key]] <- data.frame(arm = "budget", scenario = e$id, dataset = oo[[1]]$job$dataset, variant = "best_of_3",
                            z_hat = par[["z"]], z_low = ci$z_low, z_high = ci$z_high,
                            low_open = ci$low_censored, high_open = ci$high_censored,
                            covers = (is.na(ci$z_low) || ci$z_low <= e$z) && (is.na(ci$z_high) || e$z <= ci$z_high),
                            excludes_1 = !is.na(ci$z_low) && ci$z_low > 1,
                            nll_precise_hat = precise_nll(par, funres_model(), d),
                            nll_precise_true = oo[[1]]$res$nll_precise_true,
                            seconds = sum(vapply(oo, function(o) o$seconds, 0)),
                            e[, c("z", "k", "s_over_h", "h_over_T", "E", "design", "generating_model", "block")])
  }
  common <- c("arm", "scenario", "dataset", "variant", "z_hat", "z_low", "z_high", "low_open", "high_open", "covers",
              "excludes_1", "nll_precise_hat", "nll_precise_true", "seconds", "z", "k", "s_over_h", "h_over_T",
              "E", "design", "generating_model", "block")
  B <- B[, common]
  if (length(b3)) B <- rbind(B, do.call(rbind, b3)[, common])
  B$deficit <- B$nll_precise_hat - B$nll_precise_true
  best <- tapply(B$nll_precise_hat, paste(B$scenario, B$dataset), min)
  B$gap <- B$nll_precise_hat - best[paste(B$scenario, B$dataset)]
  B$success_0.5 <- B$gap < 0.5; B$success_1 <- B$gap < 1
  write_table(B, "budget_fits")
}

OF <- load_arm("D2_D3_final")
if (length(OF)) {
  write_table(do.call(rbind, lapply(OF, function(o) cbind(o$res, seconds = o$seconds))), "D2_D3_final_preset_fits")
  write_table(do.call(rbind, lapply(OF, function(o) data.frame(data = o$job$scenario, o$checks))),
              "D2_D3_final_preset_checks")
  write_table(do.call(rbind, lapply(OF, function(o) data.frame(data = o$job$scenario, o$fit$profile))),
              "D2_D3_final_preset_profiles")
}
cat("Summaries: source(\"validation/R/20_figures_and_tables.R\")\n")
