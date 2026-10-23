# ============================================================================
# helpers.R
#
# Shared definitions of the validation scripts of funresMech 1.2.0: designs,
# scenarios, true parameters, simulation of datasets, seeds, and summaries of
# fits. Sourced by 11, 12, 15, 17 and 18. Run all scripts from the root of the
# package repository.
#
# Seeds. Every simulated dataset and every fit uses a fixed seed, read from
# validation/data/seeds_datasets.csv and validation/data/seeds_fits.csv.
# These tables reproduce exactly the seeds of the runs reported in the article
# and in S1 Text.
# ============================================================================

T_EXP <- 24                                   # trial duration (h)

DESIGNS <- list(
  large = rep(c(2, 4, 8, 12, 16, 24, 32, 48), each = 8),                 # 8 densities, 64 trials
  small = rep(c(2, 6, 12, 24, 48), each = 6),                            # 5 densities, 30 trials
  dense = rep(c(2, 3, 4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48), each = 5)) # 13 densities, 65 trials

# Base scenario: h/T = 0.025 (h = 0.6 h), s/h = 0.5, E = 20 expected encounters
# per trial at the geometric mean density, default model.
BASE_SCENARIO <- list(z = 1, k = 1, s_over_h = 0.5, h_over_T = 0.025, E = 20, design = "large",
                      generating_model = "correct", sigma_a = 0)

scenario <- function(...) {
  e <- utils::modifyList(BASE_SCENARIO, list(...))
  as.data.frame(e, stringsAsFactors = FALSE)
}
scenario_id <- function(e) sprintf("%s_z%g_k%g_sh%g_hT%g_E%g_%s", e$generating_model, e$z, e$k,
                                   e$s_over_h, e$h_over_T, e$E, e$design)

true_parameters <- function(e) {
  dens <- DESIGNS[[e$design]]
  x_ref <- exp(mean(log(unique(dens))))
  c(a = e$E / (T_EXP * x_ref^e$z), h = e$h_over_T * T_EXP, z = e$z, k = e$k,
    s = e$s_over_h * e$h_over_T * T_EXP)
}

generating_model <- function(e) switch(e$generating_model,
  search_weibull   = funres_model(search = "weibull"),
  search_lognormal = funres_model(search = "lognormal"),
  handling_gamma   = funres_model(handling = "gamma"),
  predator         = funres_model(process = "predator"),
  funres_model())

SEEDS_DATA <- utils::read.csv(file.path("validation", "data", "seeds_datasets.csv"), stringsAsFactors = FALSE)
SEEDS_FIT  <- utils::read.csv(file.path("validation", "data", "seeds_fits.csv"), stringsAsFactors = FALSE)

data_seed <- function(id, j) {
  b <- SEEDS_DATA$seed_base[SEEDS_DATA$scenario == id]
  if (length(b) != 1) stop("No seed for scenario ", id)
  b * 1000 + j
}
fit_seed <- function(script, arm, scen, dataset, variant) {
  i <- which(SEEDS_FIT$script == script & SEEDS_FIT$arm == arm & SEEDS_FIT$scenario == scen &
             SEEDS_FIT$variant == variant &
             (if (is.na(dataset)) is.na(SEEDS_FIT$dataset) else SEEDS_FIT$dataset %in% dataset))
  if (length(i) != 1) stop("No seed for ", paste(script, arm, scen, dataset, variant))
  SEEDS_FIT$seed[i]
}

# One simulated dataset (one row per trial). With individual variation in a,
# each trial draws its own a from a lognormal distribution with the same mean.
simulate_dataset <- function(e, j) {
  dens <- DESIGNS[[e$design]]
  p <- true_parameters(e)
  set.seed(data_seed(e$id, j))
  if (e$sigma_a > 0) {
    response <- vapply(dens, function(x) {
      pi <- p; pi[["a"]] <- p[["a"]] * exp(stats::rnorm(1, -e$sigma_a^2 / 2, e$sigma_a))
      simulate_response(funres_model(), x, T_EXP, pi, n_sim = 1)$response
    }, numeric(1))
    return(data.frame(dens = dens, par = response))
  }
  out <- do.call(rbind, lapply(unique(dens), function(x)
    simulate_response(generating_model(e), x, T_EXP, p, n_sim = sum(dens == x))))
  data.frame(dens = out$dens, par = out$response)
}

read_published <- function(name) {
  d <- utils::read.csv(file.path("validation", "data", paste0(sub("_without_atypical", "", name), ".csv")))
  if (grepl("_without_atypical", name)) {
    stopifnot(d$dens[13] == 10, d$par[13] == 2)     # the atypical trial (2 of 10 hosts)
    d <- d[-13, ]
  }
  d
}

# Precise NLL: mean of five evaluations with N_PRECISE simulated trials per density.
precise_nll <- function(par, model, d, n = N_PRECISE) {
  mean(vapply(1:5, function(i) funresMech:::.nll_sim(par, model, d, T_EXP, n), 0))
}

# Fitted objects without the verification runs, to save disk space.
lighten <- function(f) { f$verification_runs <- NULL; f$joint <- NULL; f }

fit_row <- function(f, z_true = NA_real_) {
  ci <- f$ci; p <- f$par; pr <- f$profile
  data.frame(z_hat = p[["z"]], a_hat = p[["a"]], h_hat = p[["h"]], k_hat = p[["k"]],
             s_hat = p[["s"]], z_low = ci$z_low, z_high = ci$z_high,
             low_open = ci$low_censored, high_open = ci$high_censored,
             covers = (is.na(ci$z_low) || ci$z_low <= z_true) && (is.na(ci$z_high) || z_true <= ci$z_high),
             excludes_1 = !is.na(ci$z_low) && ci$z_low > 1,
             nll = f$nll, nll_se = f$nll_se, dAIC_z1 = f$aic$delta_AIC,
             dAIC_BB = f$aic$delta_AIC_BB, n_points = nrow(pr),
             floor_points = sum(pr$n_floor > 0, na.rm = TRUE),
             floor_retries = if (is.null(pr$n_retry)) NA_integer_ else sum(pr$n_retry, na.rm = TRUE),
             extended = f$extended, refined = isTRUE(f$refined),
             widened_bounds = paste(f$widened, collapse = ","),
             reliability = attr(reliability(f), "overall"))
}

checks_long <- function(f) {
  r <- as.data.frame(reliability(f))
  data.frame(check = r$check, status = r$status, value = suppressWarnings(as.numeric(r$value)))
}

wilson <- function(x) {
  n <- length(x); p <- mean(x); z <- 1.96
  c(p = p, lower = (p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n),
    upper = (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / (1 + z^2 / n), n = n)
}

# Run a list of jobs in parallel. Each job saves its own file, so an
# interrupted run resumes where it stopped.
run_jobs <- function(pending, fun, globals, workers, log_file) {
  if (!nrow(pending)) return(invisible(0))
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  t0 <- Sys.time()
  fs <- vector("list", nrow(pending))
  for (i in seq_len(nrow(pending))) {
    g <- globals; g$fun <- fun; g$pending <- pending; g$i <- i
    fs[[i]] <- future::future(fun(pending[i, ]), seed = TRUE, globals = g)
    cat(sprintf("[%s] started %d/%d\n", format(Sys.time(), "%d %H:%M"), i, nrow(pending)))
  }
  errors <- 0
  for (i in seq_along(fs)) tryCatch(future::value(fs[[i]]), error = function(e) {
    errors <<- errors + 1
    msg <- sprintf("%s | job %d | %s", format(Sys.time()), i, conditionMessage(e))
    message("Failed: ", msg)
    cat(msg, "\n", file = log_file, append = TRUE)
  })
  cat("Total time:", format(Sys.time() - t0), "; errors:", errors, "\n")
  invisible(errors)
}
