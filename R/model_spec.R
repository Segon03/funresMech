############################################################
# Model specification: distributions of search and handling times,
# foraging process and fixed parameters
############################################################

# Family codes shared with src/engine.cpp
.search_codes   <- c(gamma = 0L, lognormal = 1L, weibull = 2L)
.handling_codes <- c(lognormal = 0L, gamma = 1L, weibull = 2L, constant = 3L)
.user_r_code    <- 9L
.user_cpp_code  <- 10L

#' User-defined distribution of search or handling times
#'
#' Wraps a sampler written by the user so that it can be used as the
#' distribution of search times or handling times in [funres_model()].
#'
#' The sampler must return non-negative random variates with **mean 1**; the
#' simulation engine multiplies them by the required mean (the expected search
#' time \eqn{1/\lambda}, with \eqn{\lambda = a x^z}, or the mean handling time
#' \eqn{h}). Its single shape parameter is estimated together with the other
#' parameters, within `lower` and `upper`, unless `fixed` is given. In fitted
#' models the shape parameter is reported as `k` (search) or `s` (handling).
#'
#' Two kinds of sampler are accepted:
#' * an R function `function(n, shape)` returning `n` values. It is called in
#'   blocks of 256 values and is therefore slower than the built-in families;
#'   use [test_sampler()] to check it and measure its cost;
#' * an external pointer to a compiled C++ function `double f(double shape)`
#'   returning one value (see the example below). It is as fast as the
#'   built-in families. It must draw random numbers with R's generator (for
#'   example `R::rgamma()`), so that results are reproducible with
#'   [set.seed()].
#'
#' @param fun R function `function(n, shape)` or external pointer to a C++
#'   sampler (see Details).
#' @param lower,upper Bounds of the shape parameter (positive). The search is
#'   carried out on the log scale.
#' @param fixed Optional fixed value of the shape parameter; if given, the
#'   parameter is not estimated.
#' @param name Label used in printed output.
#' @return An object of class `funres_sampler`.
#' @seealso [funres_model()], [test_sampler()]
#' @examples
#' # Inverse Gaussian search times with mean 1 and shape parameter `shape`
#' # (written in R; requires no extra package)
#' r_invgauss1 <- function(n, shape) {
#'   mu <- 1
#'   y <- rnorm(n)^2
#'   x <- mu + mu^2 * y / (2 * shape) -
#'        mu / (2 * shape) * sqrt(4 * mu * shape * y + mu^2 * y^2)
#'   ifelse(runif(n) <= mu / (mu + x), x, mu^2 / x)
#' }
#' ig <- funres_sampler(r_invgauss1, lower = 0.01, upper = 100,
#'                      name = "inverse Gaussian")
#' m <- funres_model(search = ig)
#' m
#'
#' \dontrun{
#' # The same idea in C++ (compiled with Rcpp); the function must return a
#' # mean-1 variate and use R's random number generator.
#' Rcpp::sourceCpp(code = '
#'   #include <Rcpp.h>
#'   typedef double (*funres_rng_fun)(double);
#'   double gamma_mean1(double shape) { return R::rgamma(shape, 1.0 / shape); }
#'   // [[Rcpp::export]]
#'   SEXP gamma_mean1_ptr() {
#'     return Rcpp::XPtr<funres_rng_fun>(new funres_rng_fun(&gamma_mean1), true);
#'   }')
#' m2 <- funres_model(search = funres_sampler(gamma_mean1_ptr(), 1e-3, 100,
#'                                            name = "Gamma (C++)"))
#' }
#' @export
funres_sampler <- function(fun, lower = 1e-3, upper = 1e3, fixed = NULL,
                           name = "user-defined") {
  is_r   <- is.function(fun)
  is_cpp <- typeof(fun) == "externalptr"
  if (!is_r && !is_cpp)
    stop("'fun' must be an R function(n, shape) or an external pointer to a ",
         "C++ sampler.", call. = FALSE)
  if (is_r && !all(c("n", "shape") %in% names(formals(fun))))
    stop("The R sampler must have arguments 'n' and 'shape'.", call. = FALSE)
  if (!is.null(fixed)) {
    if (!is.numeric(fixed) || length(fixed) != 1 || !(fixed > 0))
      stop("'fixed' must be one positive number.", call. = FALSE)
  } else if (!(is.numeric(lower) && is.numeric(upper) && lower > 0 &&
               upper > lower && is.finite(upper))) {
    stop("Need 0 < lower < upper < Inf for the shape parameter.", call. = FALSE)
  }
  structure(list(fun = fun, code = if (is_r) .user_r_code else .user_cpp_code,
                 lower = lower, upper = upper, fixed = fixed, name = name,
                 type = if (is_r) "R" else "C++"),
            class = "funres_sampler")
}

#' @export
print.funres_sampler <- function(x, ...) {
  cat("<funres_sampler> ", x$name, " (", x$type, " sampler, mean 1); shape ",
      if (is.null(x$fixed)) paste0("in [", x$lower, ", ", x$upper, "]")
      else paste0("fixed at ", x$fixed), "\n", sep = "")
  invisible(x)
}

#' Specify the mechanistic functional response model
#'
#' Defines the stochastic foraging process from which the likelihood is built
#' by simulation (Okuyama 2012, 2026): the distribution of search times, the
#' distribution of handling times, whether attacked hosts remain in the arena,
#' and which parameters are fixed instead of estimated.
#'
#' **Search times** have mean \eqn{1/\lambda}, where \eqn{\lambda = a x^z} is
#' the encounter rate at host (or prey) density \eqn{x}. Their dispersion is
#' set by \eqn{k}, defined by \eqn{CV^2 = 1/k} for all built-in families
#' (larger \eqn{k}, more regular search). For the Gamma family \eqn{k} is the
#' shape parameter, as in Okuyama (2026); `"exponential"` is the Gamma family
#' with \eqn{k = 1}, the waiting-time distribution of the standard
#' consumer-resource model.
#'
#' **Handling times** have mean \eqn{h} and standard deviation \eqn{s} on the
#' natural scale (same units as \eqn{h}); `"constant"` fixes \eqn{s = 0}.
#'
#' **Process.** `"parasitoid"`: hosts are never removed; each encounter is
#' with a host chosen at random, so hosts can be attacked repeatedly
#' (superparasitism), and the response is the number of distinct hosts
#' attacked. `"predator"`: each encounter consumes one prey, which is removed,
#' so the encounter rate falls to \eqn{a (x - n)^z} after \eqn{n} captures.
#' The choice should follow the biology of the system rather than the taxon.
#' For a parasitoid, the predator process represents perfect and cost-free
#' discrimination of parasitized hosts: once attacked, a host is no longer
#' encountered, so it costs neither search nor handling time (Okuyama 2026,
#' Appendix S1). It does not represent host feeding, in which the host is
#' killed without being parasitized and would not be counted in the response.
#' Real behavior may lie between the two processes (partial or costly
#' discrimination), which neither process represents exactly. The process
#' can change the estimate of `z`, so when it is uncertain fit both and
#' compare them by AIC (they have the same number of parameters).
#'
#' **End of the trial.** A host or prey may still be in handling when the
#' trial ends. With `unfinished = "include"` (default) the attack counts as
#' soon as the host or prey is encountered, as in Okuyama (2026), who assumed
#' that parasitism was successful even if the trial ended during handling.
#' With `"exclude"` an attack counts only if its handling ends by `T`, so a
#' host or prey still being handled at the end is not counted (in the
#' parasitoid process, a host already attacked earlier stays counted). Choose
#' the convention that matches how the trials were scored; for example,
#' `"exclude"` when a host counts as parasitized only after oviposition is
#' complete.
#'
#' The distributions are assumptions of the analysis. No family is correct by
#' default; when the biology of the system suggests another distribution,
#' specify it here (built-in alternatives or [funres_sampler()]) and compare
#' the fits.
#'
#' @param search Distribution of search times: `"gamma"` (default),
#'   `"exponential"`, `"lognormal"`, `"weibull"` or a [funres_sampler()].
#' @param handling Distribution of handling times: `"lognormal"` (default),
#'   `"gamma"`, `"weibull"`, `"constant"` or a [funres_sampler()].
#' @param process `"parasitoid"` (default) or `"predator"`.
#' @param unfinished `"include"` (default) or `"exclude"`: whether a host or
#'   prey still being handled at the end of the trial counts in the response
#'   (see Details).
#' @param fix Named list of parameters held fixed, among `k`, `s` and `z`
#'   (for example `list(z = 1)` for a type II response).
#' @return An object of class `funres_model`.
#' @references
#' Okuyama, T. (2012). Biological Control 60: 103-107.
#' \doi{10.1016/j.biocontrol.2011.10.008}
#'
#' Okuyama, T. (2026). Journal of Applied Entomology. \doi{10.1111/jen.70148}
#' @seealso [fit_funres()], [simulate_response()], [funres_sampler()]
#' @examples
#' funres_model()                                   # default: Gamma search, lognormal handling
#' funres_model(search = "exponential", handling = "constant")
#' funres_model(search = "weibull", handling = "gamma", process = "predator")
#' funres_model(fix = list(z = 1))                  # type II response
#' funres_model(unfinished = "exclude")             # attacks count only if handling ends by T
#' @export
funres_model <- function(search = "gamma", handling = "lognormal",
                         process = c("parasitoid", "predator"), fix = list(),
                         unfinished = c("include", "exclude")) {
  process <- match.arg(process)
  unfinished <- match.arg(unfinished)
  if (!is.list(fix) || (length(fix) && is.null(names(fix))))
    stop("'fix' must be a named list, e.g. list(k = 1).", call. = FALSE)
  bad <- setdiff(names(fix), c("k", "s", "z"))
  if (length(bad))
    stop("Only k, s and z can be fixed (got: ", paste(bad, collapse = ", "),
         ").", call. = FALSE)
  for (nm in names(fix))
    if (!is.numeric(fix[[nm]]) || length(fix[[nm]]) != 1 || !is.finite(fix[[nm]]))
      stop("Fixed value of ", nm, " must be one finite number.", call. = FALSE)
  if (!is.null(fix$k) && !(fix$k > 0)) stop("k must be > 0.", call. = FALSE)
  if (!is.null(fix$s) && !(fix$s >= 0)) stop("s must be >= 0.", call. = FALSE)

  # ---- search ----
  if (inherits(search, "funres_sampler")) {
    s_spec <- list(name = search$name, code = search$code, fun = search$fun,
                   user = TRUE, lower = search$lower, upper = search$upper)
    if (!is.null(search$fixed)) fix$k <- search$fixed
  } else {
    search <- match.arg(search, c("gamma", "exponential", "lognormal", "weibull"))
    if (search == "exponential") {
      if (!is.null(fix$k) && fix$k != 1)
        stop("The exponential distribution has k = 1.", call. = FALSE)
      fix$k <- 1
    }
    s_spec <- list(name = search, user = FALSE, fun = NULL,
                   code = .search_codes[[if (search == "exponential") "gamma" else search]])
  }

  # ---- handling ----
  if (inherits(handling, "funres_sampler")) {
    h_spec <- list(name = handling$name, code = handling$code, fun = handling$fun,
                   user = TRUE, lower = handling$lower, upper = handling$upper)
    if (!is.null(handling$fixed)) fix$s <- handling$fixed
    if (!is.null(fix$s) && !(fix$s > 0))
      stop("The shape of a user handling sampler must be > 0.", call. = FALSE)
  } else {
    handling <- match.arg(handling, c("lognormal", "gamma", "weibull", "constant"))
    if (handling == "constant") {
      if (!is.null(fix$s) && fix$s != 0)
        stop("Constant handling time means s = 0.", call. = FALSE)
      fix$s <- 0
    }
    h_spec <- list(name = handling, user = FALSE, fun = NULL,
                   code = .handling_codes[[handling]])
  }

  structure(list(search = s_spec, handling = h_spec, process = process,
                 unfinished = unfinished,
                 fix = fix[c("k", "s", "z")[c("k", "s", "z") %in% names(fix)]]),
            class = "funres_model")
}

#' @export
print.funres_model <- function(x, ...) {
  fx <- x$fix
  free <- setdiff(c("a", "h", "z", "k", "s"), names(fx))
  cat("Mechanistic functional response model (simulation-based likelihood)\n")
  cat("  Process:        ", x$process,
      if (x$process == "parasitoid") " (hosts not removed; distinct hosts counted)"
      else " (consumed prey removed)", "\n", sep = "")
  cat("  Search times:   ", x$search$name, ", mean 1/(a x^z)",
      if (!x$search$user) ", CV^2 = 1/k" else ", shape k", "\n", sep = "")
  cat("  Handling times: ", x$handling$name, ", mean h",
      if (!x$handling$user) ", SD s" else ", shape s", "\n", sep = "")
  cat("  End of trial:   ", if (.excludes_unfinished(x))
        "a host or prey still in handling at T is not counted"
      else "a host or prey still in handling at T is counted", "\n", sep = "")
  cat("  Free parameters: ", paste(free, collapse = ", "),
      if (length(fx)) paste0("; fixed: ",
                             paste(names(fx), unlist(fx), sep = " = ", collapse = ", ")),
      "\n", sep = "")
  invisible(x)
}

# Default model (Okuyama 2026)
.default_model <- function() funres_model()

# TRUE when attacks whose handling ends after T are not counted (models created
# before this option existed have no 'unfinished' element: they include them)
.excludes_unfinished <- function(model) identical(model$unfinished, "exclude")

# Number of free parameters (z counted unless fixed)
.n_free <- function(model) 5L - length(model$fix)

# Simulated responses for one density with natural parameters
.simulate_x <- function(model, x, T, a, h, z, k, s, n_sim, cut = TRUE) {
  sim_trials_cpp(as.integer(x), a, h, z, k, s, T, as.integer(n_sim),
                 model$search$code, model$handling$code,
                 model$process == "predator", cut,
                 model$search$fun, model$handling$fun, FALSE,
                 .excludes_unfinished(model))
}

#' Simulate functional response trials
#'
#' Simulates the number of hosts parasitized (or prey consumed) in trials of
#' duration `T` at each host density, under a [funres_model()] and parameter
#' values. Useful to explore a model, to check a fit (compare the simulated
#' distributions with the data) and to design experiments.
#'
#' @param model A [funres_model()]; the default has Gamma search times and lognormal handling times.
#' @param dens Host (or prey) densities (positive integers).
#' @param T Trial duration, in the time unit of `h` and `s` (and the inverse of
#'   the unit of the encounter rate).
#' @param par Named numeric vector with `a`, `h`, `z`, `k` and `s` (fixed
#'   parameters of `model` may be omitted). For user samplers `k` and `s` are
#'   their shape parameters.
#' @param n_sim Number of simulated trials per density.
#' @return Data frame with columns `dens`, `trial` and `response`.
#' @examples
#' set.seed(1)
#' sims <- simulate_response(dens = c(2, 10, 50), T = 24,
#'                           par = c(a = 0.012, h = 0.63, z = 2.5, k = 0.9, s = 0.54),
#'                           n_sim = 200)
#' tapply(sims$response, sims$dens, mean)
#' @export
simulate_response <- function(model = funres_model(), dens, T, par, n_sim = 1000) {
  if (!inherits(model, "funres_model")) stop("'model' must be a funres_model().", call. = FALSE)
  p <- .complete_par(par, model)
  dens <- as.integer(dens)
  if (anyNA(dens) || any(dens < 1)) stop("Densities must be integers >= 1.", call. = FALSE)
  if (!(is.numeric(T) && length(T) == 1 && T > 0)) stop("T must be > 0.", call. = FALSE)
  out <- lapply(dens, function(x) {
    data.frame(dens = x, trial = seq_len(n_sim),
               response = .simulate_x(model, x, T, p[["a"]], p[["h"]], p[["z"]],
                                      p[["k"]], p[["s"]], n_sim))
  })
  do.call(rbind, out)
}

# Fill fixed parameters into a named vector c(a, h, z, k, s)
.complete_par <- function(par, model) {
  par <- unlist(par)
  for (nm in names(model$fix)) par[[nm]] <- model$fix[[nm]]
  miss <- setdiff(c("a", "h", "z", "k", "s"), names(par))
  if (length(miss)) stop("Missing parameter(s): ", paste(miss, collapse = ", "), call. = FALSE)
  par[c("a", "h", "z", "k", "s")]
}

# Simulated negative log-likelihood for natural parameters par = c(a, h, z, k, s).
# Densities are visited in order of appearance, as in Okuyama (2026). Returns
# 1e10 outside the valid domain so that optimizers discard the point.
.nll_sim <- function(par, model, data_spp, T, n_sim, cut = TRUE, detail = FALSE) {
  a <- par[["a"]]; h <- par[["h"]]; z <- par[["z"]]; k <- par[["k"]]; s <- par[["s"]]
  if (!all(is.finite(c(a, h, z, k, s))) || a <= 0 || h <= 0 || k <= 0 || s < 0)
    return(1e10)
  if (model$handling$user && !(s > 0)) return(1e10)
  minp <- .Machine$double.xmin
  nll <- 0
  n_floor <- 0L
  for (x in unique(data_spp$dens)) {
    if (x < 1) next
    y_obs <- data_spp$par[data_spp$dens == x]
    if (anyNA(y_obs) || any(y_obs < 0 | y_obs > x)) return(1e10)
    sim <- .simulate_x(model, x, T, a, h, z, k, s, n_sim, cut)
    frec <- tabulate(sim + 1L, nbins = x + 1L) / n_sim
    pr <- frec[y_obs + 1L]
    n_floor <- n_floor + sum(pr == 0)
    nll <- nll - sum(log(pmax(pr, minp)))
  }
  # n_floor: observed outcomes never produced by the simulation; each adds
  # -log(double.xmin) = 708.4 to the NLL
  if (detail) list(nll = nll, n_floor = n_floor) else nll
}

#' Check a user sampler
#'
#' Draws variates from a [funres_sampler()] with the simulation engine and
#' reports their mean (it should be close to 1), coefficient of variation and
#' the time per draw, at a few values of the shape parameter.
#'
#' @param sampler A [funres_sampler()].
#' @param shape Values of the shape parameter to test.
#' @param n Number of draws per value.
#' @return Data frame with one row per shape value.
#' @examples
#' r_gamma1 <- function(n, shape) rgamma(n, shape = shape, rate = shape)
#' test_sampler(funres_sampler(r_gamma1, 0.01, 100), shape = c(0.5, 2))
#' @export
test_sampler <- function(sampler, shape = NULL, n = 1e5) {
  if (!inherits(sampler, "funres_sampler")) stop("Need a funres_sampler().", call. = FALSE)
  if (is.null(shape))
    shape <- if (!is.null(sampler$fixed)) sampler$fixed
             else signif(exp(seq(log(sampler$lower), log(sampler$upper), length.out = 3)), 3)
  do.call(rbind, lapply(shape, function(sh) {
    tm <- system.time(x <- draw_times_cpp(as.integer(n), TRUE, sampler$code, 1, sh,
                                          sampler$fun))[["elapsed"]]
    data.frame(shape = sh, mean = mean(x), cv = stats::sd(x) / mean(x),
               min = min(x), us_per_draw = 1e6 * tm / n,
               mean_ok = abs(mean(x) - 1) < 5 * stats::sd(x) / sqrt(n) + 0.01)
  }))
}
