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
#' Collects the settings of the optimizer and of the Monte Carlo
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
#' The final preset reproduces the settings of the published analyses of this method (Okuyama 2026).
#'
#' **Default optimizer.** Differential evolution with common random numbers
#' and a Nelder-Mead polish was chosen in a comparison of strategies with
#' equal budgets on simulated data with known parameters (validation/,
#' phase 4): it reached the best optimum found within 0.5 NLL units in 80% of
#' the runs and never missed it by more than 2.2 units, whereas plain
#' differential evolution (the method of version 1.1.0) reached it in 43% of
#' the runs, CMA-ES in 78% (with occasional failures of up to 68 units) and
#' CRS2 in 67%. Even the best strategy missed the best optimum in one run out
#' of five with the exploratory budget: check every fit with [reliability()]
#' and [verify_fit()].
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
#'   the optimizer) is used for the profile, the estimate of z and its
#'   interval; their standard error measures Monte Carlo noise.
#' @param method Optimization method: `"DEoptim"` (differential evolution,
#'   default), `"CMAES"` (covariance matrix adaptation evolution strategy,
#'   package 'cmaes') or `"CRS2"` (controlled random search with local
#'   mutation, package 'nloptr'). All use a budget of about `itermax * NP`
#'   evaluations of the likelihood.
#' @param crn Use common random numbers during the search (default): every
#'   evaluation uses the same random seed, so the objective becomes
#'   deterministic in the parameters (but piecewise constant). Re-evaluations
#'   at the optimum always use fresh random numbers.
#' @param polish Number of Nelder-Mead evaluations used to refine the best
#'   point found by the global method (default 200; 0 = none).
#' @param bounds Named list of search bounds `c(lower, upper)` replacing the
#'   defaults: `lambda` (encounter rate at the reference density, per unit
#'   time), `h` (mean handling time, time units), `k`, `cv_h` (s / h) and `z`.
#'   See Details of [fit_funres()].
#' @param widen Maximum number of automatic widenings of a search bound on
#'   which the optimum ends (each widening multiplies the bound by 10 and
#'   refits the affected points of the profile). 0 disables it.
#' @param bb Fit the reference Beta-Binomial model with the same expected
#'   value (Okuyama 2026) and report its AIC.
#' @param cut Stop each simulated trial once all hosts have been attacked
#'   (exact; only saves time).
#' @param floor_retry Maximum number of times the re-evaluations at an
#'   optimum are repeated with four times more simulated trials when an
#'   observed outcome was never produced by the simulation (probability
#'   floor, which adds 708 NLL units). Such a zero is usually Monte Carlo
#'   noise in the tail of a distribution; more simulations resolve it. A point
#'   where it persists is flagged and left out of the interval of z.
#' @return An object of class `funres_control`.
#' @examples
#' funres_control()
#' funres_control("final")
#' funres_control("exploratory", n_sim = 3000)
#' @export
funres_control <- function(preset = c("standard", "exploratory", "final"),
                           n_sim, itermax, NP, reltol, steptol, n_reeval,
                           method = "DEoptim", bounds = list(), widen = 2L,
                           bb = TRUE, crn = TRUE, polish = 200L, cut = TRUE,
                           floor_retry = 2L) {
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
  method <- match.arg(method, c("DEoptim", "CMAES", "CRS2"))
  if (p$n_sim < 100) stop("n_sim must be >= 100.", call. = FALSE)
  if (p$NP < 4) stop("NP must be >= 4.", call. = FALSE)
  if (p$n_reeval < 1) stop("n_reeval must be >= 1.", call. = FALSE)
  structure(c(list(preset = preset), p,
              list(method = method, bounds = bounds, widen = as.integer(widen),
                   bb = isTRUE(bb), crn = isTRUE(crn), polish = as.integer(polish),
                   cut = isTRUE(cut), floor_retry = max(0L, as.integer(floor_retry)))),
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

# .optimize
#
# Minimizes fn over the box [ptab$lower, ptab$upper] with the method of the
# control object and the same budget of evaluations (about itermax * NP).
# Returns theta, value (best value found, biased low for a noisy objective),
# trace (best value per generation/iteration), pop (final population, when
# the method has one) and nfeval. An optional Nelder-Mead polish (ctrl$polish
# evaluations) refines the best point.
.optimize <- function(fn, ptab, ctrl, initialpop = NULL) {
  lo <- ptab$lower; up <- ptab$upper; d <- length(lo)
  rng <- up - lo
  to_theta <- function(u) lo + u * rng                 # unit cube -> theta
  fn_u <- function(u) if (any(u < 0 | u > 1)) 1e10 else fn(to_theta(u))
  out <- switch(ctrl$method,
    DEoptim = {
      # DEoptim advises NP >= 10 x number of parameters; the population is a
      # setting of funres_control() (same budget for every method), so that
      # advisory warning is muffled. Other warnings pass through.
      r <- withCallingHandlers(
        DEoptim(fn = fn, lower = lo, upper = up,
                control = DEoptim.control(itermax = ctrl$itermax, NP = ctrl$NP,
                                          reltol = ctrl$reltol, steptol = ctrl$steptol,
                                          trace = FALSE, initialpop = initialpop)),
        warning = function(w) {
          if (grepl("ten times the length of the parameter vector", conditionMessage(w)))
            invokeRestart("muffleWarning")
        })
      list(theta = unname(r$optim$bestmem), value = r$optim$bestval,
           trace = r$member$bestvalit, pop = r$member$pop, nfeval = r$optim$nfeval)
    },
    CMAES = {
      if (!requireNamespace("cmaes", quietly = TRUE))
        stop("method 'CMAES' needs the package 'cmaes'.", call. = FALSE)
      r <- cmaes::cma_es(stats::runif(d), fn_u, lower = rep(0, d), upper = rep(1, d),
                         control = list(maxit = ctrl$itermax, lambda = ctrl$NP,
                                        sigma = 0.3, diag.value = TRUE,
                                        stopfitness = -Inf))
      tr <- if (!is.null(r$diagnostic$value)) cummin(apply(r$diagnostic$value, 1, min)) else numeric(0)
      list(theta = to_theta(r$par), value = r$value, trace = tr, pop = NULL,
           nfeval = r$counts[[1]])
    },
    CRS2 = {
      if (!requireNamespace("nloptr", quietly = TRUE))
        stop("method 'CRS2' needs the package 'nloptr'.", call. = FALSE)
      r <- nloptr::nloptr(stats::runif(d), fn_u, lb = rep(0, d), ub = rep(1, d),
                          opts = list(algorithm = "NLOPT_GN_CRS2_LM",
                                      maxeval = ctrl$itermax * ctrl$NP,
                                      population = ctrl$NP, xtol_rel = 0,
                                      ftol_rel = 0, ftol_abs = 0))
      list(theta = to_theta(r$solution), value = r$objective, trace = numeric(0),
           pop = NULL, nfeval = r$iterations)
    })
  if (isTRUE(ctrl$polish > 0)) {
    u0 <- (out$theta - lo) / rng
    pl <- stats::optim(u0, fn_u, method = "Nelder-Mead",
                       control = list(maxit = ctrl$polish, reltol = 0))
    if (pl$value < out$value) {
      out$theta <- to_theta(pl$par); out$value <- pl$value
    }
    out$nfeval <- out$nfeval + pl$counts[[1]]
  }
  out
}

# .fit_point
#
# One optimization: z fixed (profile) or z in the parameter vector (joint,
# z_free = TRUE). Differential evolution on the dimensionless theta scale
# (R/reparam.R), then n_reeval fresh evaluations of the NLL at the optimum.
# The optimizer's best value (nll_opt) is the minimum of many noisy
# evaluations and is biased downwards; the mean of the re-evaluations
# (nll_reeval) is not, and is the value used downstream.
.fit_point <- function(z, data_spp, T_exp, model, ctrl, bounds, H_ref,
                       z_free = FALSE, initialpop = NULL) {
  ptab <- .par_table(model, bounds, z_free)
  # Common random numbers (optional): every evaluation during the search uses
  # the same seed, which makes the objective deterministic in theta; the
  # caller's random stream is restored after each evaluation.
  crn_seed <- if (isTRUE(ctrl$crn)) sample.int(.Machine$integer.max, 1L) else NULL
  fn <- function(theta) {
    if (!is.null(crn_seed)) {
      old <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old, envir = globalenv()))
      set.seed(crn_seed)
    }
    .nll_sim(.theta_to_nat(theta, ptab, model, z, H_ref, T_exp), model,
             data_spp, T_exp, ctrl$n_sim, ctrl$cut)
  }
  res <- .optimize(fn, ptab, ctrl, initialpop)

  theta <- res$theta
  names(theta) <- ptab$name
  par <- .theta_to_nat(theta, ptab, model, z, H_ref, T_exp)
  re <- lapply(seq_len(ctrl$n_reeval), function(i)
    .nll_sim(par, model, data_spp, T_exp, ctrl$n_sim, ctrl$cut, detail = TRUE))
  reevals <- vapply(re, function(r) if (is.list(r)) r$nll else r, numeric(1))
  n_floor <- max(vapply(re, function(r) if (is.list(r)) r$n_floor else NA_integer_, 0L))
  # An observed outcome that the simulation never produced adds 708 NLL units.
  # Usually it is Monte Carlo noise in a tail: repeat the re-evaluations with
  # four times more trials (at most floor_retry times).
  n_retry <- 0L
  n_sim_re <- ctrl$n_sim
  while (is.finite(n_floor) && n_floor > 0 && n_retry < (if (is.null(ctrl$floor_retry)) 0L else ctrl$floor_retry)) {
    n_retry <- n_retry + 1L
    n_sim_re <- n_sim_re * 4L
    re <- lapply(seq_len(ctrl$n_reeval), function(i)
      .nll_sim(par, model, data_spp, T_exp, n_sim_re, ctrl$cut, detail = TRUE))
    reevals <- vapply(re, function(r) if (is.list(r)) r$nll else r, numeric(1))
    n_floor <- max(vapply(re, function(r) if (is.list(r)) r$n_floor else NA_integer_, 0L))
  }

  rng <- ptab$upper - ptab$lower
  at_lower <- (theta <= ptab$lower + 0.01 * rng) & !ptab$natural_lower
  at_upper <- theta >= ptab$upper - 0.01 * rng
  at_bound <- stats::setNames(at_lower | at_upper, ptab$name)
  bound_side <- stats::setNames(ifelse(at_lower, "lower", ifelse(at_upper, "upper", "")),
                         ptab$name)

  tr <- res$trace
  it <- length(tr)
  late_gain <- if (it >= 4) tr[ceiling(0.75 * it)] - tr[it] else NA_real_
  pop <- res$pop
  spread <- if (!is.null(pop)) apply(pop, 2, stats::sd) / rng else rep(NA_real_, length(rng))

  list(z = par[["z"]], par = par, theta = theta, nll_opt = res$value,
       nll_reeval = mean(reevals),
       nll_reeval_se = if (length(reevals) > 1) stats::sd(reevals) / sqrt(length(reevals))
                       else NA_real_,
       n_floor = n_floor, n_retry = n_retry,
       at_bound = at_bound, bound_side = bound_side, iter = it, nfeval = res$nfeval,
       late_gain = late_gain, spread = stats::setNames(spread, ptab$name), pop = pop)
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
    }, seed = NULL)),   # seed = NULL: loading packages on the worker may draw random numbers
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
.eval_profile_points <- function(zs, model, ..., progress = NULL) {
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
  out <- vector("list", length(fs))
  for (j in seq_along(fs)) {
    out[[j]] <- value(fs[[j]])
    if (is.function(progress)) progress(zs[j])
  }
  out
}

.points_to_df <- function(pts) {
  do.call(rbind, lapply(pts, function(p) {
    data.frame(z = p$z, nll = p$nll_reeval, nll_se = p$nll_reeval_se,
               nll_opt = p$nll_opt,
               a = p$par[["a"]], h = p$par[["h"]], k = p$par[["k"]], s = p$par[["s"]],
               at_bound = any(p$at_bound),
               bound_names = if (any(p$at_bound))
                 paste0(names(p$at_bound)[p$at_bound], "(", p$bound_side[p$at_bound],
                        ")", collapse = ",") else "",
               n_floor = p$n_floor, n_retry = if (is.null(p$n_retry)) 0L else p$n_retry, iter = p$iter, late_gain = p$late_gain,
               max_spread = suppressWarnings(max(p$spread, na.rm = TRUE)),
               stringsAsFactors = FALSE)
  }))
}

# Adaptive profile around a joint estimate: points at z_hat -/+ step, 2 step,
# ... (the step doubles after every two points on a side) until the NLL
# exceeds the minimum by more than the 95% threshold plus a margin, or the
# bounds of z are reached. Both sides are evaluated together (parallel).
.adaptive_profile <- function(joint, run, step, zb, include_one = TRUE,
                              margin = 0.5, max_points = 40) {
  prof <- joint
  z0 <- joint$z
  thr <- .ci95_threshold() + margin
  sides <- list(low = list(open = TRUE, last = z0, step = step, n = 0),
                high = list(open = TRUE, last = z0, step = step, n = 0))
  first <- if (include_one && abs(z0 - 1) > 1e-8 && 1 > zb[1] && 1 < zb[2]) 1 else numeric(0)
  while ((sides$low$open || sides$high$open) && nrow(prof) < max_points) {
    zs <- first; first <- numeric(0)
    for (nm in names(sides)) {
      sd_ <- sides[[nm]]
      if (!sd_$open) next
      znew <- if (nm == "low") sd_$last - sd_$step else sd_$last + sd_$step
      znew <- round(znew, 6)
      if (nm == "low" && znew <= zb[1]) { znew <- zb[1]; sides[[nm]]$open <- FALSE }
      if (nm == "high" && znew >= zb[2]) { znew <- zb[2]; sides[[nm]]$open <- FALSE }
      sides[[nm]]$last <- znew
      sides[[nm]]$n <- sd_$n + 1
      if (sides[[nm]]$n %% 2 == 0) sides[[nm]]$step <- sd_$step * 2
      zs <- c(zs, znew)
    }
    zs <- setdiff(zs, round(prof$z, 6))
    if (length(zs)) prof <- rbind(prof, run(zs))
    m <- min(prof$nll)
    for (nm in names(sides)) {
      zl <- sides[[nm]]$last
      v <- prof$nll[which.min(abs(prof$z - zl))]
      if (v - m > thr) sides[[nm]]$open <- FALSE
    }
  }
  prof <- prof[order(prof$z), ]
  rownames(prof) <- NULL
  prof
}

# Widen the dimensionless bounds that were hit, by a factor of 10 on the
# natural scale (for cv_h, the upper bound is multiplied by 10; its lower
# bound 0 is natural). Bounds of user samplers are not changed here.
.widen_bounds <- function(bounds, hits) {
  key <- c(lambda = "lambda", h = "h", k = "k", s = "cv_h", z = "z")
  for (hh in hits) {
    nm <- key[[hh$name]]
    if (hh$side == "lower") bounds[[nm]][1] <- bounds[[nm]][1] / 10
    if (hh$side == "upper") bounds[[nm]][2] <- bounds[[nm]][2] * 10
  }
  bounds
}

# Bound hits of a profile row ("h(upper),k(lower)" -> list of name/side)
.parse_hits <- function(txt) {
  if (is.na(txt) || !nzchar(txt)) return(list())
  lapply(strsplit(txt, ",")[[1]], function(x) {
    list(name = sub("\\(.*", "", x), side = sub(".*\\((.*)\\)", "\\1", x))
  })
}

# Text for a lower limit of z that the profile does not cross.
.low_open_text <- function(z_min, z_floor = 0.1) {
  if (z_min <= z_floor) {
    paste0("lower limit of z not identified (profile flat down to z = ", z_min, ")")
  } else {
    paste0("lower limit of z not reached within the grid (z <= ", z_min,
           "); extend the z grid downwards")
  }
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

# Extra values of z to resolve the interval. Inside each bracket where the
# profile crosses the 95% threshold, and in the two grid cells next to the
# minimum, points are added so that the spacing is at most `step`, which is
# the base grid step or, for a narrow interval, a quarter of its width (never
# below 0.05). The linear or spline interpolation between distant points is
# the main source of error of the limits.
.refine_points <- function(z, z_low, z_high, step0, max_new = 12, z_hat = NA_real_) {
  width <- if (is.finite(z_low) && is.finite(z_high)) z_high - z_low else NA_real_
  step <- if (is.finite(width)) max(0.05, min(step0, round(width / 4 / 0.05) * 0.05)) else step0
  fill <- function(za, zb) {
    if (zb - za <= 1.5 * step) return(numeric(0))
    cand <- seq(za + step, zb - 0.5 * step, by = step)
    cand[cand > za + 0.5 * step]
  }
  bracket_points <- function(zc) {
    if (is.na(zc) || !any(z < zc) || !any(z >= zc)) return(numeric(0))
    fill(max(z[z < zc]), min(z[z >= zc]))
  }
  around_hat <- function() {
    if (!is.finite(z_hat)) return(numeric(0))
    c(if (any(z < z_hat)) fill(max(z[z < z_hat]), z_hat),
      if (any(z > z_hat)) fill(z_hat, min(z[z > z_hat])))
  }
  new <- sort(unique(round(c(bracket_points(z_low), bracket_points(z_high), around_hat()), 10)))
  new <- new[!new %in% round(z, 10)]
  utils::head(new, max_new)
}

# Interval of z from the profile points that are usable: points where an
# observed outcome was never simulated (probability floor) even after the
# retries are left out, because their NLL is inflated by 708 units. The
# limits are found on a monotone spline through the points.
.ci_profile <- function(prof) {
  ok <- is.na(prof$n_floor) | prof$n_floor == 0
  if (sum(ok) < 2) ok <- rep(TRUE, nrow(prof))
  .profile_ci(prof$z[ok], prof$nll[ok], interpolate = "spline")
}

#' Fit the mechanistic functional response model
#'
#' Estimates the parameters of a [funres_model()] by maximum likelihood, with
#' the likelihood of the observed counts approximated by simulation of the
#' foraging process (Okuyama 2012, 2026). The density-scaling exponent `z` is
#' estimated by likelihood profiling: it is fixed at each value of a grid,
#' the other parameters are optimized, and the estimate is the minimum of the
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
#' **Reliability.** A simulated likelihood is noisy, and an optimizer can stop
#' at a suboptimal point; a returned estimate is not by itself evidence of a
#' reliable fit. Check the `notes` of the result, which flag parameters on a
#' search bound, large Monte Carlo noise, an optimizer still improving at the
#' end of the run and open intervals, and act on them.
#'
#' @param data Data frame with one row per trial.
#' @param T Trial duration (time units of the experiment).
#' @param model A [funres_model()].
#' @param z_grid Values of z for the profile (z = 1 is added when the grid
#'   spans it). Ignored when `z` is fixed in `model`.
#' @param control A [funres_control()].
#' @param extend_z If the profile does not cross the 95% threshold above the
#'   estimate, evaluate the values of `z_extend` above the grid; if it does
#'   not cross it below the estimate, evaluate the values of `z_extend` below
#'   the grid (within the search bounds of z, by default \eqn{[0.1, 20]}).
#' @param z_extend Extra values of z for `extend_z`. Values above the grid are
#'   used for an open upper limit and values below it for an open lower limit.
#' @param refine Add grid points inside brackets wider than the grid step
#'   where the profile crosses the threshold.
#' @param z_method `"profile"` (default): z is fixed at each value of
#'   `z_grid` and the minimum of the profile is the estimate. `"joint"`: z is
#'   estimated together with the other parameters, and the profile is then
#'   computed only around the estimate, in steps of `z_step` (growing when the
#'   profile is flat), until the 95% threshold is crossed on each side or the
#'   search bounds of z are reached. The joint estimate is kept only if no
#'   profile point has a lower NLL. **Experimental:** in simulations with
#'   known parameters (validation/, phase 4) the joint method gave intervals
#'   that covered the true z in 67% of the datasets, against 96% for the
#'   profile, because with five free parameters the neighboring profile
#'   points are less well optimized. Confirm its results with the profile.
#' @param z_step Initial step of the adaptive profile of the joint method.
#' @param dens,response Names of the columns with host density and number of
#'   hosts parasitized (prey consumed).
#' @param time_unit Label of the time unit, used in printed output.
#' @param progress Optional function called with the value of z after each
#'   point of the profile is fitted (used by the shiny application).
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
                       extend_z = TRUE, z_extend = c(0.1, 0.25, 8, 10, 15, 20),
                       refine = TRUE, z_method = c("profile", "joint"), z_step = 0.25,
                       dens = "dens", response = "par", time_unit = "time units",
                       progress = NULL) {
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
                                            c(list(zs, model = model, progress = progress), args)))
  z_method <- match.arg(z_method)
  joint <- NULL
  if (z_method == "joint" && !z_fixed) {
    # z estimated with the other parameters; the profile is then computed
    # only around the estimate, as far as needed to find the 95% limits.
    joint <- .points_to_df(do.call(.eval_profile_points,
                                   c(list(NA_real_, model = model, z_free = TRUE,
                                          progress = progress), args)))
    prof <- .adaptive_profile(joint, run, z_step, bounds$z, include_one = TRUE)
    extend_z <- FALSE
  } else {
    prof <- run(z_grid)
  }
  ci <- if (!z_fixed && nrow(prof) > 1) .ci_profile(prof) else NULL

  extended <- FALSE
  z_extend <- sort(unique(z_extend[z_extend >= bounds$z[1] & z_extend <= bounds$z[2]]))
  if (isTRUE(extend_z) && !is.null(ci) && (ci$high_censored || ci$low_censored)) {
    z_new <- c(if (ci$low_censored) z_extend[z_extend < min(prof$z)],
               if (ci$high_censored) z_extend[z_extend > max(prof$z)])
    if (length(z_new) > 0) {
      prof <- rbind(prof, run(z_new))
      prof <- prof[order(prof$z), ]
      ci <- .ci_profile(prof)
      extended <- TRUE
    }
  }
  refined <- FALSE
  for (round in seq_len(if (isTRUE(refine)) 2L else 0L)) {
    if (is.null(ci)) break
    z_new <- .refine_points(prof$z, ci$z_low, ci$z_high, step0, z_hat = ci$z_hat)
    if (length(z_new) == 0) break
    prof <- rbind(prof, run(z_new))
    prof <- prof[order(prof$z), ]
    ci <- .ci_profile(prof)
    refined <- TRUE
  }
  rownames(prof) <- NULL

  # Automatic widening of search bounds hit at the optimum (at most twice).
  # Only the points whose optimum ended on a bound are refitted.
  widened <- character(0)
  user_par <- c(if (model$search$user) "k", if (model$handling$user) "s")
  for (round in seq_len(max(0L, as.integer(control$widen)))) {
    i_best <- which.min(prof$nll)
    hits <- Filter(function(hh) !(hh$name %in% user_par),
                   .parse_hits(prof$bound_names[i_best]))
    if (!length(hits)) break
    bounds <- .widen_bounds(bounds, hits)
    args$bounds <- bounds
    widened <- c(widened, vapply(hits, function(hh) paste0(hh$name, "(", hh$side, ")"), ""))
    redo <- prof$at_bound
    prof <- rbind(prof[!redo, , drop = FALSE], run(prof$z[redo]))
    prof <- prof[order(prof$z), ]
    rownames(prof) <- NULL
    if (!is.null(ci)) ci <- .ci_profile(prof)
  }

  i_best <- which.min(prof$nll)
  par <- c(a = prof$a[i_best], h = prof$h[i_best], z = prof$z[i_best],
           k = prof$k[i_best], s = prof$s[i_best])
  n_free <- .n_free(model)
  aic_full <- 2 * prof$nll[i_best] + 2 * n_free
  nll_1 <- if (!z_fixed) .nll_at_z(prof$z, prof$nll, 1) else NA_real_
  aic_restricted <- if (is.na(nll_1)) NA_real_ else 2 * nll_1 + 2 * (n_free - 1)

  bb <- if (isTRUE(control$bb)) tryCatch(.fit_bb_model(data_spp, model),
                                          error = function(e) NULL) else NULL
  aic_bb <- if (is.null(bb)) NA_real_ else 2 * bb$nll + 2 * bb$n_par

  fit <- structure(list(
    call = cl, model = model, control = control, data = data_spp, T = T,
    time_unit = time_unit, H_ref = H_ref, bounds = bounds, profile = prof,
    par = par, nll = prof$nll[i_best], nll_se = prof$nll_se[i_best],
    n_par = n_free, ci = ci,
    aic = data.frame(AIC_full = aic_full, AIC_restricted = aic_restricted,
                     delta_AIC = aic_restricted - aic_full,
                     AIC_BB = aic_bb, delta_AIC_BB = aic_bb - aic_full),
    bb = bb,
    extended = extended, refined = refined, widened = widened, z_method = z_method,
    joint = joint, seed = seed,
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
  n_drop <- sum(prof$n_floor > 0, na.rm = TRUE)
  if (n_drop > 0 && !(is.finite(prof$n_floor[i_best]) && prof$n_floor[i_best] > 0))
    notes <- c(notes, paste0(
      n_drop, " profile point(s) had observed outcomes never produced by the ",
      "simulation, even after repeating the re-evaluation with more trials; ",
      "they are left out of the interval of z."))
  if (!is.null(ci) && ci$high_censored)
    notes <- c(notes, paste0(
      "The profile does not cross the 95% threshold above z_hat: ",
      .high_open_text(max(prof$z), max(c(z_extend, fit$bounds$z[2]))), "."))
  if (!is.null(ci) && ci$low_censored)
    notes <- c(notes, paste0(
      "The profile does not cross the 95% threshold below z_hat: ",
      .low_open_text(min(prof$z), fit$bounds$z[1]), "."))
  if (is.finite(prof$n_floor[i_best]) && prof$n_floor[i_best] > 0)
    notes <- c(notes, paste0(
      prof$n_floor[i_best], " observed outcome(s) were never produced by the ",
      "simulation at the optimum (their probability was set to the numerical ",
      "floor, adding 708 NLL units each). The fit is not reliable: increase ",
      "n_sim, check the atypical-trial screening, and check that the model ",
      "can produce these outcomes."))
  if (length(fit$widened))
    notes <- c(notes, paste0(
      "Search bounds were widened automatically after the optimum ended on ",
      "them: ", paste(fit$widened, collapse = ", "), "."))
  if (prof$at_bound[i_best])
    notes <- c(notes, paste0(
      "At z_hat the optimizer ended on a search bound (",
      prof$bound_names[i_best], "): the estimate of that parameter is only a ",
      "limit. Widen the bound (funres_control(bounds = ...)) and refit, or ",
      "check whether the data can identify it."))
  gap <- prof$nll[i_best] - prof$nll_opt[i_best]
  if (is.finite(gap) && gap > 2)
    notes <- c(notes, paste0(
      "The re-evaluated NLL at the optimum exceeds the optimizer's by ",
      round(gap, 2), " units: Monte Carlo noise is large; increase n_sim."))
  if (is.finite(prof$nll_se[i_best]) && prof$nll_se[i_best] > 0.25)
    notes <- c(notes, paste0(
      "Monte Carlo standard error of the NLL at the optimum is ",
      signif(prof$nll_se[i_best], 2), ", large relative to the 1.92 threshold ",
      "of the interval; increase n_sim."))
  if (is.finite(prof$late_gain[i_best]) && prof$late_gain[i_best] > 1)
    notes <- c(notes, paste0(
      "The optimizer was still improving in the last quarter of the run (by ",
      round(prof$late_gain[i_best], 2), " NLL units): increase itermax."))
  if (fit$control$preset == "exploratory")
    notes <- c(notes, "Exploratory settings: do not report these results without a fit with the standard or final settings.")
  notes
}

# Density unit of a: the encounter rate is a x^z, so a is per unit of time
# per (host or prey)^z, and its units depend on the value of z
.a_density_unit <- function(fit) {
  who <- if (identical(fit$model$process, "predator")) "prey" else "host"
  z <- fit$par[["z"]]
  if (isTRUE(all.equal(z, 1))) who else paste0(who, "^", format(signif(z, 3)))
}

#' @export
print.funres_fit <- function(x, digits = 4, ...) {
  print(x$model)
  u <- x$time_unit
  cat("\nEstimates (T = ", x$T, " ", u, "; ", nrow(x$data), " trials):\n", sep = "")
  p <- x$par
  lab <- c(a = paste0("encounter-rate coefficient (1/", u, " per ", .a_density_unit(x), ")"),
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
  if (!is.na(x$aic$delta_AIC_BB))
    cat(sprintf("AIC(Beta-Binomial, same mean) - AIC(mechanistic) = %.2f (positive favors the mechanistic model)\n",
                x$aic$delta_AIC_BB))
  rel <- reliability(x)
  cat("\nReliability: ", toupper(attr(rel, "overall")),
      if (is.null(x[["verification"]])) " (cheap checks only; run verify_fit() before reporting)",
      ". See reliability().\n", sep = "")
  bad <- rel[rel$status != "ok", , drop = FALSE]
  for (j in seq_len(nrow(bad)))
    cat(" * [", bad$status[j], "] ", bad$check[j], ": ", bad$detail[j], "\n", sep = "")
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
                        extend_z = TRUE, z_extend = c(0.1, 0.25, 8, 10, 15, 20),
                        refine = TRUE, model = funres_model(), bounds = list()) {
  ctrl <- funres_control("exploratory", n_sim = n_sim, itermax = itermax, NP = NP,
                         n_reeval = n_reeval, bounds = bounds)
  f <- fit_funres(data_spp, T_exp, model = model, z_grid = z_grid, control = ctrl,
                  extend_z = extend_z, z_extend = z_extend, refine = refine)
  list(profile = f$profile, par = f$par, nll = f$nll, aic = f$aic, ci = f$ci,
       extended = f$extended, refined = f$refined, notes = f$notes)
}
