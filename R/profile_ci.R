############################################################
# Likelihood-ratio interval for z from a profile
############################################################

#' Threshold of the 95\% likelihood-ratio interval
#'
#' The criterion is 2 * (NLL(z) - NLL_min) <= qchisq(0.95, 1), i.e.
#' NLL(z) <= NLL_min + qchisq(0.95, 1) / 2 (= +1.92). Versions <= 1.0.4
#' omitted the "/ 2", which gives an interval of about 99.5\%.
#' @importFrom stats qchisq
#' @noRd
.umbral_ic95 <- function() qchisq(0.95, 1) / 2

#' 95\% interval of z by linear interpolation of the profile
#'
#' @param z_grid evaluated values of z (sorted).
#' @param nll_vals NLL minimised at each value of z_grid.
#' @param interpolate `"linear"` (between neighbouring points) or `"spline"`
#'   (shape-preserving cubic Hermite interpolation, PCHIP: it passes through
#'   the points and never goes beyond the values of two neighbouring points,
#'   so it creates no minimum or maximum that the points do not have; the
#'   limits are then searched on a grid of step 0.01).
#' @return list(z_hat, z_low, z_high, low_censored, high_censored,
#'   threshold, min_nll). A side is "censored" when the profile never
#'   crosses the threshold on that side within the grid: the real limit is
#'   then <= min(z_grid) or >= max(z_grid).
#' @noRd
.profile_ci <- function(z_grid, nll_vals, interpolate = c("linear", "spline")) {
  stopifnot(length(z_grid) == length(nll_vals), !is.unsorted(z_grid))
  interpolate <- match.arg(interpolate)
  if (interpolate == "spline" && length(z_grid) >= 3) {
    f <- .pchip(z_grid, nll_vals)
    zz <- sort(unique(c(z_grid, seq(min(z_grid), max(z_grid), by = 0.01))))
    ci <- .profile_ci(zz, f(zz), "linear")
    ci$z_hat <- z_grid[which.min(nll_vals)]
    ci$min_nll <- min(nll_vals)
    ci$threshold <- ci$min_nll + .umbral_ic95()
    return(ci)
  }
  min_nll   <- min(nll_vals)
  idx_min   <- which.min(nll_vals)
  threshold <- min_nll + .umbral_ic95()

  z_low <- NA_real_
  for (i in seq(idx_min, 1, by = -1)) {
    if (nll_vals[i] >= threshold) {
      z1 <- z_grid[i];     n1 <- nll_vals[i]
      z2 <- z_grid[i + 1]; n2 <- nll_vals[i + 1]
      z_low <- z1 + (z2 - z1) * (threshold - n1) / (n2 - n1)
      break
    }
  }

  z_high <- NA_real_
  for (i in seq(idx_min, length(z_grid), by = 1)) {
    if (nll_vals[i] >= threshold) {
      z1 <- z_grid[i - 1]; n1 <- nll_vals[i - 1]
      z2 <- z_grid[i];     n2 <- nll_vals[i]
      z_high <- z1 + (z2 - z1) * (threshold - n1) / (n2 - n1)
      break
    }
  }

  list(z_hat = z_grid[idx_min], z_low = z_low, z_high = z_high,
       low_censored  = is.na(z_low),
       high_censored = is.na(z_high),
       threshold = threshold, min_nll = min_nll)
}

#' NLL of the profile at z0 (exact if evaluated, otherwise linear
#' interpolation between neighbouring points)
#' @noRd
.nll_at_z <- function(z_grid, nll_vals, z0) {
  if (z0 < min(z_grid) || z0 > max(z_grid)) return(NA_real_)
  idx <- which.min(abs(z_grid - z0))
  if (abs(z_grid[idx] - z0) < 1e-8) return(nll_vals[idx])
  idx_low  <- max(which(z_grid <= z0))
  idx_high <- min(which(z_grid >= z0))
  z_low <- z_grid[idx_low];   n_low  <- nll_vals[idx_low]
  z_high <- z_grid[idx_high]; n_high <- nll_vals[idx_high]
  n_low + (n_high - n_low) * (z0 - z_low) / (z_high - z_low)
}

#' Shape-preserving piecewise cubic Hermite interpolation (PCHIP)
#'
#' Fritsch and Butland (1984): the slope at a point is the harmonic mean of the
#' two neighbouring secants if they have the same sign and 0 at a local
#' extremum, so the interpolant never leaves the range of two neighbouring
#' values. (`stats::splinefun(method = "monoH.FC")` is not used because it
#' only guarantees this for monotone data and can dip below the minimum of a
#' profile.)
#' @param x,y points, `x` increasing.
#' @return A function of the abscissae.
#' @noRd
.pchip <- function(x, y) {
  n <- length(x)
  h <- diff(x); d <- diff(y) / h
  m <- numeric(n)
  m[1] <- d[1]; m[n] <- d[n - 1]
  if (n > 2) for (k in 2:(n - 1)) {
    if (d[k - 1] * d[k] > 0) {
      w1 <- 2 * h[k] + h[k - 1]; w2 <- h[k] + 2 * h[k - 1]
      m[k] <- (w1 + w2) / (w1 / d[k - 1] + w2 / d[k])
    }
  }
  function(xx) {
    i <- pmin(pmax(findInterval(xx, x), 1L), n - 1L)
    tt <- (xx - x[i]) / h[i]
    (2 * tt^3 - 3 * tt^2 + 1) * y[i] + (tt^3 - 2 * tt^2 + tt) * h[i] * m[i] +
      (-2 * tt^3 + 3 * tt^2) * y[i + 1] + (tt^3 - tt^2) * h[i] * m[i + 1]
  }
}
