# ============================================================================
# 11_optimizer_comparison.R
#
# Comparison of optimization strategies and of two ways of estimating z
# (main text, Results; S1 Fig; S1 Text, section 8).
#
#   strategies  Four strategies with the same budget of likelihood
#               evaluations (exploratory preset: 2,000 simulated trials per
#               density, 50 generations of 40 members), at the true z:
#               DEoptim | CMA-ES | CRS2 | DEoptim with common random numbers
#               and a Nelder-Mead polish. Five datasets per scenario, three
#               independent runs per strategy and dataset (135 fits per
#               strategy). The optimum of each run is re-evaluated precisely
#               (5 x 20,000 simulated trials per density); the gap of a run
#               is its precise NLL minus the best one found for the dataset.
#   z_method    Likelihood profile on a grid (0.5 to 5, step 0.5) against
#               joint estimation of z followed by an adaptive profile for the
#               interval (three datasets per scenario, 27 fits each).
#
# Scenarios: z {1, 2, 3} x k {0.1, 1, 5}; h = 0.6 h, s = 0.5 h, T = 24 h, and a
# such that 20 encounters per trial are expected at the geometric mean
# density. Design: densities 2, 4, 8, 12, 16, 24, 32, 48 with 8 trials each.
#
# Usage, from the root of the package repository with funresMech installed
# (packages cmaes and nloptr are also needed):
#   QUICK <- TRUE; source("validation/R/11_optimizer_comparison.R")    # check
#   rm(list = ls()); source("validation/R/11_optimizer_comparison.R")   # ~3 h on 22 cores
# Results: validation/tables/optimizer_comparison_fits.csv,
#          validation/tables/joint_vs_profile_fits.csv
# ============================================================================

library(funresMech)
library(future)
stopifnot(packageVersion("funresMech") >= "1.1.0.9000")
for (p in c("cmaes", "nloptr")) if (!requireNamespace(p, quietly = TRUE)) stop("Install package ", p)

if (!exists("QUICK")) QUICK <- FALSE
WORKERS <- if (QUICK) 2L else max(1L, parallel::detectCores() - 1L)
DIR_RES <- file.path("validation", "results", if (QUICK) "optimizers_quick" else "optimizers")
DIR_TAB <- file.path("validation", "tables")
dir.create(DIR_RES, recursive = TRUE, showWarnings = FALSE); dir.create(DIR_TAB, showWarnings = FALSE)
SUFFIX <- if (QUICK) "_quick" else ""

T_EXP <- 24
DENS <- rep(c(2, 4, 8, 12, 16, 24, 32, 48), each = 8)
X_REF <- exp(mean(log(unique(DENS))))
SCEN <- expand.grid(z = c(1, 2, 3), k = c(0.1, 1, 5))
if (QUICK) SCEN <- SCEN[c(2, 5), ]
SCEN$scenario <- sprintf("z%g_k%g", SCEN$z, SCEN$k)
SCEN$a <- 20 / (T_EXP * X_REF^SCEN$z)
SCEN$h <- 0.6; SCEN$s <- 0.5
N_DATA_A <- if (QUICK) 1 else 5
N_DATA_B <- if (QUICK) 1 else 3
N_RUNS <- if (QUICK) 1 else 3

model <- funres_model()
ctrl_base <- if (QUICK) {
  funres_control(n_sim = 300, itermax = 6, NP = 12, n_reeval = 2, bb = FALSE)
} else funres_control(n_sim = 2000, itermax = 50, n_reeval = 3, bb = FALSE)
STRATEGIES <- list(
  DEoptim            = list(method = "DEoptim"),
  CMAES              = list(method = "CMAES"),
  CRS2               = list(method = "CRS2"),
  DEoptim_CRN_polish = list(method = "DEoptim", crn = TRUE, polish = 200))

# ---- simulated datasets (fixed seeds) ------------------------------------------
dataset_of <- function(i_scen, j) {
  e <- SCEN[i_scen, ]
  set.seed(1000 * i_scen + j)
  s <- simulate_response(model, dens = unique(DENS), T = T_EXP,
                         par = c(a = e$a, h = e$h, z = e$z, k = e$k, s = e$s), n_sim = 8)
  data.frame(dens = s$dens, par = s$response)
}

# ---- jobs ------------------------------------------------------------------------
jobs_A <- do.call(rbind, lapply(seq_len(nrow(SCEN)), function(i)
  expand.grid(part = "strategies", i_scen = i, dataset = seq_len(N_DATA_A), strategy = names(STRATEGIES),
              run = seq_len(N_RUNS), stringsAsFactors = FALSE)))
jobs_B <- do.call(rbind, lapply(seq_len(nrow(SCEN)), function(i)
  expand.grid(part = "z_method", i_scen = i, dataset = seq_len(N_DATA_B), strategy = c("profile", "joint"),
              run = 1, stringsAsFactors = FALSE)))
jobs <- rbind(jobs_A, jobs_B)
jobs$file <- file.path(DIR_RES, sprintf("%s_%s_d%d_%s_r%d.rds", jobs$part, SCEN$scenario[jobs$i_scen],
                                        jobs$dataset, jobs$strategy, jobs$run))
pending <- jobs[!file.exists(jobs$file), ]
cat(sprintf("Jobs: %d in total, %d pending, %d workers.\n", nrow(jobs), nrow(pending), WORKERS))

run_one <- function(job) {
  library(funresMech)
  future::plan(future::sequential)
  e <- SCEN[job$i_scen, ]
  d <- dataset_of(job$i_scen, job$dataset)
  t0 <- proc.time()[["elapsed"]]
  if (job$part == "strategies") {
    st <- STRATEGIES[[job$strategy]]
    ctrl <- do.call(funres_control, c(list("exploratory", n_sim = ctrl_base$n_sim, itermax = ctrl_base$itermax,
                                           NP = ctrl_base$NP, n_reeval = ctrl_base$n_reeval, bb = FALSE), st))
    set.seed(job$run * 7919 + job$i_scen * 101 + job$dataset)
    pt <- funresMech:::.fit_point(e$z, d, T_EXP, model, ctrl, funresMech:::.make_bounds(list(), T_EXP),
                                  funresMech:::.h_ref(d))
    sec <- proc.time()[["elapsed"]] - t0
    nll_precise <- mean(vapply(1:5, function(i)
      funresMech:::.nll_sim(pt$par, model, d, T_EXP, if (QUICK) 2000 else 20000), 0))
    out <- data.frame(scenario = e$scenario, z_true = e$z, k_true = e$k, dataset = job$dataset,
                      strategy = job$strategy, run = job$run, nll_opt = pt$nll_opt, nll_reeval = pt$nll_reeval,
                      nll_precise = nll_precise, a_ratio = pt$par[["a"]] / e$a, h_ratio = pt$par[["h"]] / e$h,
                      k_ratio = pt$par[["k"]] / e$k, s_ratio = pt$par[["s"]] / e$s,
                      on_bound = any(pt$at_bound), n_evaluations = pt$nfeval, seconds = sec)
  } else {
    set.seed(job$i_scen * 31 + job$dataset)
    f <- suppressWarnings(fit_funres(
      d, T_EXP, model, z_grid = if (QUICK) c(1, 2, 3) else seq(0.5, 5, by = 0.5),
      control = ctrl_base, extend_z = TRUE, refine = FALSE, z_method = job$strategy, z_step = 0.25))
    sec <- proc.time()[["elapsed"]] - t0
    ci <- f$ci
    out <- data.frame(scenario = e$scenario, z_true = e$z, k_true = e$k, dataset = job$dataset,
                      strategy = job$strategy, run = 1, z_hat = f$par[["z"]], z_low = ci$z_low, z_high = ci$z_high,
                      covers = (is.na(ci$z_low) || ci$z_low <= e$z) && (is.na(ci$z_high) || e$z <= ci$z_high),
                      open = ci$low_censored || ci$high_censored,
                      excludes_1 = (!is.na(ci$z_low) && ci$z_low > 1) || (!is.na(ci$z_high) && ci$z_high < 1),
                      n_points = nrow(f$profile), reliability = attr(reliability(f), "overall"), seconds = sec)
  }
  saveRDS(out, job$file)
  TRUE
}

if (nrow(pending)) {
  plan(multisession, workers = WORKERS)
  t_start <- Sys.time()
  fs <- lapply(seq_len(nrow(pending)), function(i)
    future(run_one(pending[i, ]), seed = TRUE,
           globals = list(run_one = run_one, pending = pending, i = i, SCEN = SCEN, dataset_of = dataset_of,
                          STRATEGIES = STRATEGIES, model = model, ctrl_base = ctrl_base, T_EXP = T_EXP,
                          DENS = DENS, QUICK = QUICK)))
  for (i in seq_along(fs)) {
    tryCatch(value(fs[[i]]), error = function(e) message("Job ", i, " failed: ", conditionMessage(e)))
    if (i %% 10 == 0 || i == length(fs)) cat(sprintf("[%s] %d/%d jobs done\n", format(Sys.time(), "%H:%M"), i, length(fs)))
  }
  plan(sequential)
  cat("Total time:", format(Sys.time() - t_start), "\n")
}

# ---- tables --------------------------------------------------------------------------
read_part <- function(part) {
  f <- jobs$file[jobs$part == part & file.exists(jobs$file)]
  do.call(rbind, lapply(f, readRDS))
}
A <- read_part("strategies"); B <- read_part("z_method")
if (!is.null(A)) {
  best <- tapply(A$nll_precise, paste(A$scenario, A$dataset), min)
  A$gap <- A$nll_precise - best[paste(A$scenario, A$dataset)]
  A$success <- A$gap < 0.5
  write.csv(A, file.path(DIR_TAB, paste0("optimizer_comparison_fits", SUFFIX, ".csv")), row.names = FALSE)
  print(aggregate(cbind(success, gap, on_bound, seconds) ~ strategy, A, mean), digits = 3)
}
if (!is.null(B)) {
  write.csv(B, file.path(DIR_TAB, paste0("joint_vs_profile_fits", SUFFIX, ".csv")), row.names = FALSE)
  print(aggregate(cbind(covers, open, n_points, seconds) ~ strategy, B, mean), digits = 3)
}
