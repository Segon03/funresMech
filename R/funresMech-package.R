#' funresMech: Mechanistic Functional Response Analysis
#'
#' Implements the mechanistic, stochastic functional response model of
#' Okuyama (2012, 2026) for host-parasitoid systems. The probability
#' distribution of the number of parasitised hosts is obtained by simulating
#' the search-encounter-handling process (Gamma search times, Lognormal
#' handling times), and the parameters are estimated by maximum likelihood
#' with a simulated likelihood. The simulation engine is written in C++ with
#' 'Rcpp'.
#'
#' The package is used through the interactive application started with
#' [run_app()], which fits the model, computes the likelihood profile of the
#' density-scaling exponent `z` (95% interval by the likelihood-ratio
#' criterion), compares the model with free `z` against `z = 1` by AIC,
#' reports diagnostics and screens atypical trials.
#'
#' @references
#' Okuyama, T. (2012). A likelihood approach for functional response models.
#' *Biological Control*, 60(2), 103-107. \doi{10.1016/j.biocontrol.2011.10.008}
#'
#' Okuyama, T. (2026). Parametric assumptions in parasitoid functional
#' response analysis. *Journal of Applied Entomology*.
#' \doi{10.1111/jen.70148}
#'
#' @seealso [run_app()]
#' @useDynLib funresMech, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
