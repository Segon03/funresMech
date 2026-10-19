############################################################
# Parameterization used by the optimizer
############################################################

# Time is expressed in units of the trial duration T inside the optimizer,
# so that the search bounds do not depend on the time unit chosen by the
# user (minutes, hours, days):
#
#   lambda*_ref = a * H_ref^z * T  expected encounters per trial at the
#                                   reference density H_ref (geometric mean
#                                   of the distinct densities), ignoring
#                                   handling; one-to-one with 'a' for given z
#   h*          = h / T            mean handling time as a fraction of T
#   k                              dispersion of search times (CV^2 = 1/k);
#                                   dimensionless
#   cv_h        = s / h            coefficient of variation of handling time;
#                                   dimensionless
#   z                              density-scaling exponent
#
# DEoptim mutates candidates by adding differences between population members
# on the scale of the bounds, so the search is carried out on
#   theta = (log lambda*_ref, log h*, log k, sqrt(cv_h), log z)
# sqrt(cv_h) is kept linear so that cv_h = 0 (fixed handling time) remains
# reachable. For user samplers the shape parameters k and s are searched on
# the log scale within the bounds of funres_sampler(). Fixed parameters and z
# (when it is profiled) are not part of theta.

#' Default search bounds, on the dimensionless scale described above
#' @noRd
.default_bounds <- function() {
  list(lambda = c(1e-3, 1e4),   # expected encounters per trial at H_ref
       h      = c(1e-6, 1),     # h / T
       k      = c(1e-5, 100),   # CV^2 of search time = 1 / k
       cv_h   = c(0, 10),       # s / h
       z      = c(0.1, 20))
}

#' Search bounds for one analysis
#'
#' User bounds are given in natural units (h in the time unit of T, the
#' encounter rate per unit time at H_ref) and converted to the dimensionless
#' scale. Unknown names are an error.
#' @noRd
.make_bounds <- function(user = list(), T) {
  b <- .default_bounds()
  bad <- setdiff(names(user), c("lambda", "h", "k", "cv_h", "z"))
  if (length(bad))
    stop("Unknown bound(s): ", paste(bad, collapse = ", "),
         ". Use lambda, h, k, cv_h, z.", call. = FALSE)
  conv <- list(lambda = function(v) v * T, h = function(v) v / T,
               k = identity, cv_h = identity, z = identity)
  for (nm in names(user)) {
    v <- as.numeric(user[[nm]])
    if (length(v) != 2 || anyNA(v) || v[1] >= v[2] || v[1] < 0 ||
        (nm != "cv_h" && v[1] <= 0))
      stop("Bound '", nm, "' must be c(lower, upper) with 0 < lower < upper",
           if (nm == "cv_h") " (lower may be 0)" else "", ".", call. = FALSE)
    b[[nm]] <- conv[[nm]](v)
  }
  b
}

#' Reference density: geometric mean of the distinct host densities
#' @noRd
.h_ref <- function(data_spp) exp(mean(log(unique(data_spp$dens))))

#' Table of the parameters searched by the optimizer
#'
#' @param model funres_model()
#' @param bounds dimensionless bounds from .make_bounds()
#' @param z_free TRUE when z is estimated jointly (not profiled)
#' @return data.frame with name, lower, upper (theta scale) and natural_lower
#'   (TRUE when the lower bound is a natural boundary of the parameter space,
#'   so that ending there is not an optimization problem).
#' @noRd
.par_table <- function(model, bounds, z_free = FALSE) {
  rows <- list(
    data.frame(name = "lambda", lower = log(bounds$lambda[1]),
               upper = log(bounds$lambda[2]), natural_lower = FALSE),
    data.frame(name = "h", lower = log(bounds$h[1]), upper = log(bounds$h[2]),
               natural_lower = FALSE))
  if (is.null(model$fix$k)) {
    kb <- if (model$search$user) c(model$search$lower, model$search$upper) else bounds$k
    rows <- c(rows, list(data.frame(name = "k", lower = log(kb[1]), upper = log(kb[2]),
                                    natural_lower = FALSE)))
  }
  if (is.null(model$fix$s)) {
    if (model$handling$user) {
      rows <- c(rows, list(data.frame(name = "s", lower = log(model$handling$lower),
                                      upper = log(model$handling$upper),
                                      natural_lower = FALSE)))
    } else {
      rows <- c(rows, list(data.frame(name = "s", lower = sqrt(bounds$cv_h[1]),
                                      upper = sqrt(bounds$cv_h[2]),
                                      natural_lower = bounds$cv_h[1] == 0)))
    }
  }
  if (z_free && is.null(model$fix$z))
    rows <- c(rows, list(data.frame(name = "z", lower = log(bounds$z[1]),
                                    upper = log(bounds$z[2]), natural_lower = FALSE)))
  out <- do.call(rbind, rows)
  rownames(out) <- out$name
  out
}

#' theta -> natural parameters c(a, h, z, k, s)
#' @param z value of z when it is not in theta (profiled or fixed)
#' @noRd
.theta_to_nat <- function(theta, ptab, model, z, H_ref, T) {
  names(theta) <- ptab$name
  if ("z" %in% ptab$name) z <- exp(theta[["z"]])
  if (!is.null(model$fix$z)) z <- model$fix$z
  h <- exp(theta[["h"]]) * T
  a <- exp(theta[["lambda"]]) / (T * H_ref^z)
  k <- if (!is.null(model$fix$k)) model$fix$k else exp(theta[["k"]])
  s <- if (!is.null(model$fix$s)) {
    model$fix$s
  } else if (model$handling$user) {
    exp(theta[["s"]])
  } else {
    theta[["s"]]^2 * h
  }
  c(a = a, h = h, z = z, k = k, s = s)
}

#' natural parameters c(a, h, z, k, s) -> theta (inverse of .theta_to_nat)
#' @noRd
.nat_to_theta <- function(par, ptab, model, H_ref, T) {
  z <- par[["z"]]
  th <- c(lambda = log(par[["a"]] * H_ref^z * T), h = log(par[["h"]] / T))
  if ("k" %in% ptab$name) th[["k"]] <- log(par[["k"]])
  if ("s" %in% ptab$name)
    th[["s"]] <- if (model$handling$user) log(par[["s"]]) else sqrt(par[["s"]] / par[["h"]])
  if ("z" %in% ptab$name) th[["z"]] <- log(z)
  th[ptab$name]
}
