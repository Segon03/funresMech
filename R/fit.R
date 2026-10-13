############################################################
# Fitting the mechanistic model by simulated maximum likelihood
############################################################

#' @importFrom DEoptim DEoptim DEoptim.control
#' @importFrom future future value
#' @importFrom stats sd median
#' @noRd
NULL

#' Settings of the simulated-likelihood fit
#'
#' Collects the settings of the optimiser and of the Monte Carlo
#' approximation of the likelihood. Three presets are offered. **The
#' exploratory preset is not meant for results to be reported**: it is a
#' quick first look. Results must be checked with the reliability
#' diagnostics returned by [fit_funres()] whatever the preset.
#'
#' | preset | n_sim | itermax | NP | n_reeval |
#' |---|---|---|---|---|
#' | exploratory | 2000 | 50 | 40 | 3 |
#' | standard | 5000 | 100 | 40 | 5 |
#' | final | 10000 | 200 | 40 | 5 |
#'
#' The final preset is the configuration of Okuyama (2026).
#'
#' @param preset `"standard"` (default), `"exploratory"` or `"final"`.
#' @param n_sim Simulated trials per host density for each evaluation of the
#'   likelihood. Larger values reduce Monte Carlo noise.
#' @param itermax Generations of differential evolution (DEoptim).
#' @param NP Population size of differential evolution.
#' @param reltol,steptol Early stopping of DEoptim: the search stops when the
#'   best value has not improved by a relative amount `reltol` in `steptol`
#'   generations. By default `steptol = itermax`, i.e. no early stopping,
#'   because a noisy objective can stall temporarily.
#' @param n_reeval Number of fresh evaluations of the likelihood at each
#'   optimum. Their mean (an unbiased estimate, unlike the minimum found by
#'   the optimiser) is used for the profile, the estimate of z and its
#'   interval; their standard error measures Monte Carlo noise.
#' @param method Optimisation method. Currently `"DEoptim"`.
#' @param bounds Named list of search bounds `c(lower, upper)` replacing the
#'   defaults: `lambda` (encounter rate at the reference density, per unit
#'   time), `h` (mean handling time, time units), `k`, `cv_h` (s / h) and `z`.
#'   See Details of [fit_funres()].
#' @param cut Stop each simulated trial once all hosts have been attacked
#'   (exact; only saves time).
#' @return An object of class `funres_control`.
#' @examples
#' funres_control()
#' funres_control("final")
#' funres_control("exploratory", n_sim = 3000)
#' @export
funres_control <- function(preset = c("standard", "exploratory", "final"),
                           n_sim, itermax, NP, reltol, steptol, n_reeval,
                           method = "DEoptim", bounds = list(), cut = TRUE) {
  preset <- match.arg(preset)
  p <- switch(preset,
              exploratory = list(n_sim = 2000L, itermax = 50L, NP = 40L, n_reeval = 3L),
              standard    = list(n_sim = 5000L, itermax = 100L, NP = 40L, n_reeval = 5L),
              final       = list(n_sim = 10000L, itermax = 200L, NP = 40L, n_reeval = 5L))
  if (!missing(n_sim))    p$n_sim <- as.integer(n_sim)
  if (!missing(itermax))  p$itermax <- as.integer(itermax)
  if (!missing(NP))       p$NP <- as.integer(NP)
  if (!missing(n_reeval)) p$n_reeval <- as.integer(n_reeval)
  p$reltol  <- if (missing(reltol)) 1e-8 else reltol
  p$steptol <- if (missing(steptol)) p$itermax else as.integer(steptol)
  method <- match.arg(method, "DEoptim")
  if (p$n_sim < 100) stop("n_sim must be >= 100.", call. = FALSE)
  if (p$NP < 4) stop("NP must be >= 4.", call. = FALSE)
  if (p$n_reeval < 1) stop("n_reeval must be >= 1.", call. = FALSE)
  structure(c(list(preset = preset), p,
              list(method = method, bounds = bounds, cut = isTRUE(cut))),
            class = "funres_control")
}

#' @export
print.funres_control <- function(x, ...) {
  cat("funres_control: preset '", x$preset, "'", if (x$preset == "exploratory")
    " (exploratory: not for reporting results)", "\n", sep = "")
  cat("  n_sim = ", x$n_sim, ", itermax = ", x$itermax, ", NP = ", x$NP,
      ", n_reeval = ", x$n_reeval, ", method = ", x$method, "\n", sep = "")
  if (length(x$bounds))
    cat("  user bounds: ", paste(names(x$bounds), collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' Check the data of one species
#' @noRd
.check_data_spp <- function(data_spp) {
  if (!all(c("dens", "par") %in% names(data_spp)))
    stop("data_spp needs columns 'dens' and 'par'.", call. = FALSE)
  if (anyNA(data_spp$dens) || anyNA(data_spp$par))
    stop("Missing values in host density or parasitism.", call. = FALSE)
  if (any(data_spp$dens < 1))
    stop("Host densities must be >= 1.", call. = FALSE)
  if (any(data_spp$dens != round(data_spp$dens)) || any(data_spp$par != round(data_spp$par)))
    stop("Host densities and responses must be whole numbers.", call. = FALSE)
  if (any(data_spp$par < 0 | data_spp$par > data_spp$dens))
    stop("Parasitism must be between 0 and the host density.", call. = FALSE)
  invisible(TRUE)
}

# .fit_point
#
# One optimisation: z fixed (profile) or z in the parameter vector (joint,
# z_free = TRUE). Differential evolution on the dimensionless theta scale
# (R/reparam.R), then n_reeval fresh evaluations of the NLL at the optimum.
# The optimiser's best value (nll_opt) is the minimum of many noisy
# evaluations and is biased downwards; the mean of the re-evaluations
# (nll_reeval) is not, and is the value used downstream.
.fit_point <- function(z, data_spp, T_exp, model, ctrl, bounds, H_ref,
                       z_free = FALSE, initialpop = NULL) {
  ptab <- .par_table(model, bounds, z_free)
  fn <- function(theta) {
    .nll_sim(.theta_to_nat(theta, ptab, model, z, H_ref, T_exp), model,
             data_spp, T_exp, ctrl$n_sim, ctrl$cut)
  }
  res <- DEoptim(
    fn = fn, lower = ptab$lower, upper = ptab$upper,
    control = DEoptim.control(itermax = ctrl$itermax, NP = ctrl$NP,
                              reltol = ctrl$reltol, steptol = ctrl$steptol,
                              trace = FALSE, initialpop = initialpop))

  theta <- res$optim$bestmem
  names(theta) <- ptab$name
  par <- .theta_to_nat(theta, ptab, model, z, H_ref, T_exp)
  reevals <- vapply(seq_len(ctrl$n_reeval), function(i) fn(theta), numeric(1))

  rng <- ptab$upper - ptab$lower
  at_lower <- (theta <= ptab$lower + 1e-3 * rng) & !ptab$natural_lower
  at_upper <- theta >= ptab$upper - 1e-3 * rng
  at_bound <- setNames(at_lower | at_upper, ptab$name)

  tr <- res$member$bestvalit
  it <- res$optim$iter
  late_gain <- if (it >= 4) tr[ceiling(0.75 * it)] - tr[it] else NA_real_
  pop <- res$member$pop
  spread <- if (!is.null(pop)) apply(pop, 2, stats::sd) / rng else rep(NA_real_, length(rng))

  list(z = par[["z"]], par = par, theta = theta, nll_opt = res$optim$bestval,
       nll_reeval = mean(reevals),
       nll_reeval_se = if (length(reevals) > 1) stats::sd(reevals) / sqrt(length(reevals))
                       else NA_real_,
       at_bound = at_bound, iter = it, nfeval = res$optim$nfeval,
       late_gain = late_gain, spread = setNames(spread, ptab$name), pop = pop)
}

# TRUE when the current future plan can run this package's code: either the
# plan is sequential (same R session) or the workers have the INSTALLED
# package with the same code for the point fit.
.workers_ok <- function() {
  if (inherits(future::plan(), "sequential")) return(TRUE)
  mine <- paste(deparse(body(.fit_point)), collapse = "\n")
  theirs <- tryCatch(
    value(future({
      f <- tryCatch(utils::getFromNamespace(".fit_point", "funresMech"),
                    error = function(e) NULL)
      if (is.null(f)) NA_character_ else paste(deparse(body(f)), collapse = "\n")
    })),
    error = function(e) NA_character_)
  isTRUE(identical(mine, theirs))
}

.uses_cpp_sampler <- function(model) {
  model$search$code == .user_cpp_code || model$handling$code == .user_cpp_code
}

# Evaluates several values of z, in parallel when the user's future plan
# allows it. Random streams are derived from the caller's RNG state
# (future's seed = TRUE), so set.seed() makes the whole profile
# reproducible. Falls back to a sequential plan when the workers cannot run
# this version of the package or when a C++ sampler is used (external
# pointers cannot be sent to other R processes).
.eval_profile_points <- function(zs, model, ...) {
  seq_reason <- NULL
  if (!inherits(future::plan(), "sequential")) {
    if (.uses_cpp_sampler(model)) {
      seq_reason <- "a C++ sampler cannot be sent to parallel workers"
    } else if (!.workers_ok()) {
      seq_reason <- paste("the parallel workers do not have this version of",
                          "funresMech installed (install the package to use them)")
    }
  }
  if (!is.null(seq_reason)) {
    message("Evaluating the profile sequentially: ", seq_reason, ".")
    old_plan <- future::plan(future::sequential)
    on.exit(future::plan(old_plan), add = TRUE)
  }
  fs <- lapply(zs, function(z)
    future(utils::getFromNamespace(".fit_point", "funresMech")(z, model = model, ...),
           seed = TRUE))
  lapply(fs, value)
}

.points_to_df <- function(pts) {
  do.call(rbind, lapply(pts, function(p) {
    data.frame(z = p$z, nll = p$nll_reeval, nll_se = p$nll_reeval_se,
               nll_opt = p$nll_opt,
               a = p$par[["a"]], h = p$par[["h"]], k = p$par[["k"]], s = p$par[["s"]],
               at_bound = any(p$at_bound),
               bound_names = paste(names(p$at_bound)[p$at_bound], collapse = ","),
               iter = p$iter, late_gain = p$late_gain,
               max_spread = suppressWarnings(max(p$spread, na.rm = TRUE)),
               stringsAsFactors = FALSE)
  }))
}

# Text for an upper limit of z that the profile does not cross.
.high_open_text <- function(z_max, z_cap = 20) {
  if (z_max >= z_cap) {
    paste0("upper limit of z not identified (profile flat up to z = ", z_max, ")")
  } else {
    paste0("upper limit of z not reached within the grid (z >= ", z_max,
           "); extend the z grid")
  }
}

# Extra values of z inside a bracket where the profile crosses the 95%
# threshold, when that bracket is wider than the base grid step.
.refine_points <- function(z, z_low, z_high, step0, max_new = 12) {
  bracket_points <- function(zc) {
    if (is.na(zc) || !any(z < zc) || !any(z >= zc)) return(numeric(0))
    za <- max(z[z < zc]); zb <- min(z[z >= zc])
    if (zb - za <= 1.5 * step0) return(numeric(0))
    cand <- seq(za + step0, zb - 0.5 * step0, by = step0)
    cand[cand > za + 0.5 * step0]
  }
  new <- sort(unique(round(c(bracket_points(z_low), bracket_points(z_high)), 10)))
  new <- new[!new %in% round(z, 10)]
  utils::head(new, max_new)
}

#' Fit the mechanistic functional response model
#'
#' Estimates the parameters of a [funres_model()] by maximum likelihood, with
#' the likelihood of the observed counts approximated by simulation of the
#' foraging process (Okuyama 2012, 2026). The density-scaling exponent `z` is
#' estimated by likelihood profiling: it is fixed at each value of a grid,
#' the other parameters are optimised, and the estimate is the minimum of the
#' profile. The 95% interval contains the values of `z` whose negative
#' log-likelihood (NLL) is no more than `qchisq(0.95, 1) / 2 = 1.92` above
#' the minimum.
#'
#' **Time units.** `T`, `h` and `s` share the time unit of the experiment,
#' and `a` is per unit time. Internally time is measured in units of `T`, so
#' the default search bounds do not depend on the unit: the expected number of
#' encounters per trial at the reference density (geometric mean of the
#' densities) is searched in \eqn{[10^{-3}, 10^4]}, \eqn{h/T} in
#' \eqn{[10^{-6}, 1]}, \eqn{k} in \eqn{[10^{-5}, 100]} and \eqn{s/h} in
#' \eqn{[0, 10]}. They can be changed with `funres_control(bounds = ...)`.
#'
#' **Reliability.** A simulated likelihood is noisy, and an optimiser can stop
#' at a suboptimal point; a returned estimate is not by itself evidence of a
#' reliable fit. Check the `notes` of the result, which flag parameters on a
#' search bound, large Monte Carlo noise, an optimiser still improving at the
#' end of the run and open intervals, and act on them.
#'
#' @param data Data frame with one row per trial.
#' @param T Trial duration (time units of the experiment).
#' @param model A [funres_model()].
#' @param z_grid Values of z for the profile (z = 1 is added when the grid
#'   spans it). Ignored when `z` is fixed in `model`.
#' @param control A [funres_control()].
#' @param extend_z If the profile does not cross the 95% threshold above the
#'   estimate, evaluate the values in `z_extend` as well.
#' @param z_extend Extra values of z for `extend_z`.
#' @param refine Add grid points inside brackets wider than the grid step
#'   where the profile crosses the threshold.
#' @param dens,response Names of the columns with host density and number of
#'   hosts parasitized (prey consumed).
#' @param time_unit Label of the time unit, used in printed output.
#' @return An object of class `funres_fit`.
#' @seealso [funres_model()], [funres_control()], [simulate_response()]
#' @examples
#' \donttest{
#' d <- data.frame(dens = rep(c(2, 5, 10, 20), each = 4),
#'                 par  = c(0, 1, 1, 2, 2, 3, 3, 4, 5, 6, 6, 7, 8, 9, 10, 12))
#' set.seed(1)
#' fit <- fit_funres(d, T = 24, z_grid = c(0.5, 1, 2),
#'                   control = funres_control("exploratory", n_sim = 300,
#'                                            itermax = 10, NP = 20),
#'                   extend_z = FALSE, time_unit = "h")
#' fit
#' }
#' @export
fit_funres <- function(data, T, model = funres_model(),
                       z_grid = seq(0.5, 5, by = 0.25), control = funres_control(),
                       extend_z = TRUE, z_extend = c(8, 10, 15, 20), refine = TRUE,
                       dens = "dens", response = "par", time_unit = "time units") {
  cl <- match.call()
  if (!inherits(model, "funres_model")) stop("'model' must be a funres_model().", call. = FALSE)
  if (!inherits(control, "funres_control")) stop("'control' must be a funres_control().", call. = FALSE)
  if (!(is.numeric(T) && length(T) == 1 && is.finite(T) && T > 0))
    stop("T must be one positive number.", call. = FALSE)
  if (!all(c(dens, response) %in% names(data)))
    stop("Columns '", dens, "' and '", response, "' not found in 'data'.", call. = FALSE)
  data_spp <- data.frame(dens = data[[dens]], par = data[[response]])
  .check_data_spp(data_spp)
  seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL

  H_ref  <- .h_ref(data_spp)
  bounds <- .make_bounds(control$bounds, T)
  args <- list(data_spp = data_spp, T_exp = T, ctrl = control, bounds = bounds,
               H_ref = H_ref)

  z_fixed <- !is.null(model$fix$z)
  if (z_fixed) {
    z_grid <- model$fix$z
    extend_z <- FALSE; refine <- FALSE
  } else {
    z_grid <- sort(unique(round(z_grid, 10)))
    if (min(z_grid) <= 1 && max(z_grid) >= 1) z_grid <- sort(unique(c(z_grid, 1)))
  }
  step0 <- if (length(z_grid) > 1) stats::median(diff(z_grid)) else 1

  run <- function(zs) .points_to_df(do.call(.eval_profile_points,
                                            c(list(zs, model = model), args)))
  prof <- run(z_grid)
  ci <- if (!z_fixed && nrow(prof) > 1) .profile_ci(prof$z, prof$nll) else NULL

  extended <- FALSE
  if (isTRUE(extend_z) && !is.null(ci) && ci$high_censored) {
    z_new <- z_extend[z_extend > max(prof$z)]
    if (length(z_new) > 0) {
      prof <- rbind(prof, run(z_new))
      prof <- prof[order(prof$z), ]
      ci <- .profile_ci(prof$z, prof$nll)
      extended <- TRUE
    }
  }
  refined <- FALSE
  if (isTRUE(refine) && !is.null(ci)) {
    z_new <- .refine_points(prof$z, ci$z_low, ci$z_high, step0)
    if (length(z_new) > 0) {
      prof <- rbind(prof, run(z_new))
      prof <- prof[order(prof$z), ]
      ci <- .profile_ci(prof$z, prof$nll)
      refined <- TRUE
    }
  }
  rownames(prof) <- NULL

  i_best <- which.min(prof$nll)
  par <- c(a = prof$a[i_best], h = prof$h[i_best], z = prof$z[i_best],
           k = prof$k[i_best], s = prof$s[i_best])
  n_free <- .n_free(model)
  aic_full <- 2 * prof$nll[i_best] + 2 * n_free
  nll_1 <- if (!z_fixed) .nll_at_z(prof$z, prof$nll, 1) else NA_real_
  aic_restricted <- if (is.na(nll_1)) NA_real_ else 2 * nll_1 + 2 * (n_free - 1)

  fit <- structure(list(
    call = cl, model = model, control = control, data = data_spp, T = T,
    time_unit = time_unit, H_ref = H_ref, bounds = bounds, profile = prof,
    par = par, nll = prof$nll[i_best], nll_se = prof$nll_se[i_best],
    n_par = n_free, ci = ci,
    aic = data.frame(AIC_full = aic_full, AIC_restricted = aic_restricted,
                     delta_AIC = aic_restricted - aic_full),
    extended = extended, refined = refined, seed = seed,
    version = as.character(utils::packageVersion("funresMech"))),
    class = "funres_fit")
  fit$notes <- .fit_notes(fit, i_best, z_extend)
  fit
}

# Diagnostic notes of a fit (to be extended by the reliability checks)
.fit_notes <- function(fit, i_best, z_extend) {
  prof <- fit$profile; ci <- fit$ci
  notes <- character(0)
  if (!fit$model$search$user && is.null(fit$model$fix$k) && prof$k[i_best] < 0.1)
    notes <- c(notes, paste0(
      "k is very small (", signif(prof$k[i_best], 2), "): search times are ",
      "extremely variable and z is weakly identified from above (ridge ",
      "between z and k)."))
  if (!is.null(ci) && ci$high_censored)
    notes <- c(notes, paste0(
      "The profile does not cross the 95% threshold above z_hat: ",
      .high_open_text(max(prof$z), max(z_extend)), "."))
  if (!is.null(ci) && ci$low_censored)
    notes <- c(notes, paste0(
      "The profile does not cross the 95% threshold below z_hat: the lower ",
      "limit of z is <= ", min(prof$z), "."))
  if (prof$at_bound[i_best])
    notes <- c(notes, paste0(
      "At z_hat the optimiser ended on a search bound (",
      prof$bound_names[i_best], "): widen that bound and refit; the estimate ",
      "of that parameter is only a limit."))
  gap <- prof$nll[i_best] - prof$nll_opt[i_best]
  if (is.finite(gap) && gap > 2)
    notes <- c(notes, paste0(
      "The re-evaluated NLL at the optimum exceeds the optimiser's by ",
      round(gap, 2), " units: Monte Carlo noise is large; increase n_sim."))
  if (is.finite(prof$nll_se[i_best]) && prof$nll_se[i_best] > 0.25)
    notes <- c(notes, paste0(
      "Monte Carlo standard error of the NLL at the optimum is ",
      signif(prof$nll_se[i_best], 2), ", large relative to the 1.92 threshold ",
      "of the interval; increase n_sim."))
  if (is.finite(prof$late_gain[i_best]) && prof$late_gain[i_best] > 1)
    notes <- c(notes, paste0(
      "The optimiser was still improving in the last quarter of the run (by ",
      round(prof$late_gain[i_best], 2), " NLL units): increase itermax."))
  if (fit$control$preset == "exploratory")
    notes <- c(notes, "Exploratory settings: do not report these results without a fit with the standard or final settings.")
  notes
}

#' @export
print.funres_fit <- function(x, digits = 4, ...) {
  print(x$model)
  u <- x$time_unit
  cat("\nEstimates (T = ", x$T, " ", u, "; ", nrow(x$data), " trials):\n", sep = "")
  p <- x$par
  lab <- c(a = paste0("encounter-rate coefficient (1/", u, ")"),
           h = paste0("mean handling time (", u, ")"),
           z = "density-scaling exponent",
           k = if (x$model$search$user) "shape of search sampler" else "search dispersion (CV^2 = 1/k)",
           s = if (x$model$handling$user) "shape of handling sampler" else paste0("SD of handling time (", u, ")"))
  fixed <- names(x$model$fix)
  for (nm in names(p))
    cat(sprintf("  %-2s = %-12s %s%s\n", nm, format(signif(p[[nm]], digits)),
                lab[[nm]], if (nm %in% fixed) " [fixed]" else ""))
  if (!is.null(x$ci)) {
    f <- function(v) if (is.na(v)) "open" else format(round(v, 2))
    cat("  95% interval of z: [", f(x$ci$z_low), ", ", f(x$ci$z_high), "]\n", sep = "")
  }
  cat(sprintf("\nNLL = %.3f (Monte Carlo SE %.3f); %d free parameters; AIC = %.2f\n",
              x$nll, x$nll_se, x$n_par, x$aic$AIC_full))
  if (!is.na(x$aic$delta_AIC))
    cat(sprintf("AIC(z = 1) - AIC(z free) = %.2f\n", x$aic$delta_AIC))
  if (length(x$notes)) {
    cat("\nNotes:\n")
    for (n in x$notes) cat(" * ", strwrap(n, exdent = 3, width = 76), sep = "\n")
  }
  invisible(x)
}

#' @export
coef.funres_fit <- function(object, ...) object$par

#' @export
logLik.funres_fit <- function(object, ...) {
  structure(-object$nll, df = object$n_par, nobs = nrow(object$data), class = "logLik")
}

# Wrapper kept for the shiny application (version 1.1.0 interface)
fit_profile <- function(data_spp, T_exp, z_grid, n_sim = 3000, itermax = 50,
                        NP = 40, reltol = 1e-2, n_reeval = 3L,
                        extend_z = TRUE, z_extend = c(8, 10, 15, 20),
                        refine = TRUE, model = funres_model(), bounds = list()) {
  ctrl <- funres_control("exploratory", n_sim = n_sim, itermax = itermax, NP = NP,
                         n_reeval = n_reeval, bounds = bounds)
  f <- fit_funres(data_spp, T_exp, model = model, z_grid = z_grid, control = ctrl,
                  extend_z = extend_z, z_extend = z_extend, refine = refine)
  list(profile = f$profile, par = f$par, nll = f$nll, aic = f$aic, ci = f$ci,
       extended = f$extended, refined = f$refined, notes = f$notes)
}
