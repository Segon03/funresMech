#' funresMech: Mechanistic Functional Response Analysis with Simulation-Based Likelihood
#'
#' Fits mechanistic functional response models by simulation-based maximum
#' likelihood, the approach proposed by Okuyama (2012) and developed for
#' parasitoids by Okuyama (2026). The probability distribution of the number
#' of hosts parasitized (or prey consumed) in each trial is obtained by
#' simulating the foraging process, so the likelihood follows from explicit
#' assumptions about search, encounter and handling rather than from a
#' conventional count distribution.
#'
#' Main functions:
#' * [funres_model()]: the foraging process to simulate (distributions of
#'   search and handling times, parasitoid or predator process, fixed
#'   parameters); [funres_sampler()] adds distributions written by the user.
#' * [fit_funres()] with [funres_control()]: estimation, likelihood profile
#'   of the density-scaling exponent `z` and diagnostics.
#' * [simulate_response()]: simulated trials under a model.
#' * [run_app()]: interactive application.
#'
#' The simulation engine is written in C++ with 'Rcpp'.
#'
#' @references
#' Okuyama, T. (2012). A likelihood approach for functional response models.
#' *Biological Control*, 60(2), 103-107. \doi{10.1016/j.biocontrol.2011.10.008}
#'
#' Okuyama, T. (2026). Parametric assumptions in parasitoid functional
#' response analysis. *Journal of Applied Entomology*.
#' \doi{10.1111/jen.70148}
#'
#' @seealso [fit_funres()], [run_app()]
#' @useDynLib funresMech, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
