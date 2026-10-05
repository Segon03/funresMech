############################################################
# Log-scale parameterisation used by the optimiser
############################################################

# DEoptim mutates and crosses candidates by ADDING differences between
# members of the population, on the scale in which the bounds are given.
# On a linear scale [0, 5], values such as a ~ 1e-5 (needed at high z) or
# k ~ 0.01 occupy a tiny fraction of the range and are hardly ever
# explored. The search is therefore carried out on
#
#   theta1 = log(lambda_ref),  lambda_ref = a * H_ref^z
#            (encounter rate at the reference density H_ref; with z fixed
#             this is a one-to-one reparameterisation of 'a' that makes
#             the useful search range independent of z)
#   theta2 = log(h)
#   theta3 = log(k)
#   theta4 = u,  s = u^2   (as in Okuyama 2026; kept linear so that
#             s = 0, i.e. fixed handling time, remains reachable)
#
# H_ref is the geometric mean of the distinct host densities.

#' Bounds of the search on the theta scale
#'
#' @param k_max upper bound of the Gamma shape k (5 as in the validation
#'   against Okuyama 2026).
#' @noRd
.log_bounds <- function(k_max = 5) {
  list(
    lower = c(log_lambda_ref = log(1e-3), log_h = log(1e-3), log_k = log(1e-5), u_s = 0),
    upper = c(log_lambda_ref = log(1e3),  log_h = log(5),    log_k = log(k_max), u_s = 5)
  )
}

#' Reference density: geometric mean of the distinct host densities
#' @noRd
.h_ref <- function(data_spp) exp(mean(log(unique(data_spp$dens))))

#' theta -> natural parameters c(a, h, k, s)
#' @noRd
.theta_to_nat <- function(theta, z, H_ref) {
  c(a = exp(theta[[1]]) / H_ref^z, h = exp(theta[[2]]),
    k = exp(theta[[3]]), s = theta[[4]]^2)
}

#' Natural parameters c(a, h, k, s) -> theta
#' @noRd
.nat_to_theta <- function(par_nat, z, H_ref) {
  c(log_lambda_ref = log(par_nat[["a"]]) + z * log(H_ref),
    log_h = log(par_nat[["h"]]), log_k = log(par_nat[["k"]]),
    u_s = sqrt(par_nat[["s"]]))
}
