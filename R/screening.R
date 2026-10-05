############################################################
# Screening of atypical / influential trials (before fitting)
############################################################

# Trials that are very unlikely under the rest of the data can mask
# themselves: the fit absorbs them by shrinking k (Gamma shape of the
# search time), which in turn leaves z without an upper limit. The
# screening below is fast (seconds) and NEVER removes data: it flags
# trials to be checked and to be used in a sensitivity analysis (refit
# without them).
#
# It fits a Beta-Binomial whose mean is the closed-form functional
# response of Okuyama and, for each trial, computes the predictive
# probability of its outcome under the model fitted WITHOUT that trial
# (leave-one-out), with a Bonferroni threshold alpha / n.

#' @importFrom stats optim
#' @noRd
NULL

.mu_okuyama <- function(x, a, h, z, T = 1) {
  f <- a * x^z / (1 + a * h * x^z)
  x * (1 - exp(-f * T / x))
}

# Beta-Binomial log-pmf (alpha, beta parameterisation), no extra package
.dbbinom_log <- function(y, n, alpha, beta) {
  lchoose(n, y) + lbeta(y + alpha, n - y + beta) - lbeta(alpha, beta)
}

.bb_shapes <- function(pr, rho) {
  pr  <- pmin(pmax(pr, 1e-9), 1 - 1e-9)
  rho <- max(rho, 1e-6)
  list(alpha = pr * (1 - rho) / rho, beta = (1 - pr) * (1 - rho) / rho)
}

.negloglik_bb <- function(par, data_spp, T = 1) {
  a <- par[1]; h <- par[2]; z <- par[3]; rho <- par[4]
  if (a <= 0 || h <= 0 || z <= 0) return(1e10)
  pr <- .mu_okuyama(data_spp$dens, a, h, z, T) / data_spp$dens
  sh <- .bb_shapes(pr, rho)
  ll <- .dbbinom_log(data_spp$par, data_spp$dens, sh$alpha, sh$beta)
  if (any(!is.finite(ll))) return(1e10)
  -sum(ll)
}

# Multi-start L-BFGS-B (a single start can end in ABNORMAL_TERMINATION
# near rho -> 0); falls back to Nelder-Mead if none converges.
.fit_bb <- function(data_spp, T = 1) {
  lower <- c(1e-4, 1e-4, 0.1, 1e-6)
  upper <- c(5, 5, 5, 0.999)
  starts <- list(
    c(a = 0.1,  h = 0.05, z = 1,   rho = 0.01),
    c(a = 0.5,  h = 0.02, z = 1.5, rho = 0.05),
    c(a = 1.0,  h = 0.01, z = 2,   rho = 0.10),
    c(a = 0.3,  h = 0.10, z = 1,   rho = 0.20),
    c(a = 0.05, h = 0.03, z = 0.8, rho = 0.001)
  )
  intentos <- lapply(starts, function(st) {
    tryCatch(optim(par = st, fn = .negloglik_bb, data_spp = data_spp, T = T,
                   method = "L-BFGS-B", lower = lower, upper = upper),
             error = function(e) NULL)
  })
  ok <- Filter(function(f) !is.null(f) && f$convergence == 0, intentos)
  if (length(ok) > 0) return(ok[[which.min(vapply(ok, function(f) f$value, numeric(1)))]])

  validos <- Filter(Negate(is.null), intentos)
  if (length(validos) == 0) return(NULL)
  mejor <- validos[[which.min(vapply(validos, function(f) f$value, numeric(1)))]]
  fit_nm <- optim(par = mejor$par, fn = .negloglik_bb, data_spp = data_spp, T = T,
                  method = "Nelder-Mead", control = list(maxit = 5000))
  tryCatch(optim(par = pmin(pmax(fit_nm$par, lower), upper), fn = .negloglik_bb,
                 data_spp = data_spp, T = T, method = "L-BFGS-B",
                 lower = lower, upper = upper),
           error = function(e) fit_nm)
}

# screen_outliers
#
# @param data_spp data.frame with columns dens and par.
# @param T_exp duration of the trials.
# @param alpha global level; a trial is flagged when its two-sided
#   leave-one-out probability is < alpha / n.
# @return data.frame sorted by p_loo with columns row, dens, par,
#   expected_loo, p_loo, rho_without (over-dispersion without the trial),
#   rho_change (relative change), flagged (outlier) and
#   influences_dispersion (rho falls by more than half without it).
screen_outliers <- function(data_spp, T_exp, alpha = 0.05) {
  .check_data_spp(data_spp)
  n <- nrow(data_spp)
  full <- .fit_bb(data_spp, T_exp)
  if (is.null(full)) stop("The Beta-Binomial screening model could not be fitted.",
                          call. = FALSE)
  rho_all <- full$par[["rho"]]

  filas <- lapply(seq_len(n), function(i) {
    f <- .fit_bb(data_spp[-i, , drop = FALSE], T_exp)
    if (is.null(f)) return(NULL)
    p <- f$par
    N <- data_spp$dens[i]; y <- data_spp$par[i]
    sh  <- .bb_shapes(.mu_okuyama(N, p[["a"]], p[["h"]], p[["z"]], T_exp) / N, p[["rho"]])
    pmf <- exp(.dbbinom_log(0:N, N, sh$alpha, sh$beta))
    p_low  <- sum(pmf[seq_len(y + 1)])
    p_high <- sum(pmf[(y + 1):(N + 1)])
    data.frame(row = i, dens = N, par = y,
               expected_loo = sum((0:N) * pmf),
               p_loo = min(1, 2 * min(p_low, p_high)),
               rho_without = p[["rho"]],
               rho_change = p[["rho"]] / rho_all - 1)
  })
  r <- do.call(rbind, filas)
  r$flagged <- r$p_loo < alpha / n
  r$influences_dispersion <- r$rho_change < -0.5
  r <- r[order(r$p_loo), ]
  rownames(r) <- NULL
  attr(r, "rho_all")     <- rho_all
  attr(r, "p_threshold") <- alpha / n
  r
}
