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
# It fits a Beta-Binomial whose mean is the closed-form expected response of
# the mechanistic model (Okuyama 2026) and, for each trial, computes the
# predictive probability of its outcome under the model fitted WITHOUT that trial
# (leave-one-out), with a Bonferroni threshold alpha / n.

#' @importFrom stats optim
#' @noRd
NULL

.mu_parasitoid <- function(x, a, h, z, T = 1) {
  f <- a * x^z / (1 + a * h * x^z)
  x * (1 - exp(-f * T / x))
}

# Beta-Binomial log-pmf (alpha, beta parameterization), no extra package
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
  pr <- .mu_parasitoid(data_spp$dens, a, h, z, T) / data_spp$dens
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
  attempts <- lapply(starts, function(st) {
    tryCatch(optim(par = st, fn = .negloglik_bb, data_spp = data_spp, T = T,
                   method = "L-BFGS-B", lower = lower, upper = upper),
             error = function(e) NULL)
  })
  ok <- Filter(function(f) !is.null(f) && f$convergence == 0, attempts)
  if (length(ok) > 0) return(ok[[which.min(vapply(ok, function(f) f$value, numeric(1)))]])

  valid <- Filter(Negate(is.null), attempts)
  if (length(valid) == 0) return(NULL)
  best <- valid[[which.min(vapply(valid, function(f) f$value, numeric(1)))]]
  fit_nm <- optim(par = best$par, fn = .negloglik_bb, data_spp = data_spp, T = T,
                  method = "Nelder-Mead", control = list(maxit = 5000))
  tryCatch(optim(par = pmin(pmax(fit_nm$par, lower), upper), fn = .negloglik_bb,
                 data_spp = data_spp, T = T, method = "L-BFGS-B",
                 lower = lower, upper = upper),
           error = function(e) fit_nm)
}

#' Screen the data for atypical trials
#'
#' Flags trials whose outcome is very unlikely given the rest of the data,
#' before the mechanistic model is fitted. Such trials can dominate a
#' simulation-based fit: the model absorbs them by making search times
#' extremely variable (very small `k`), which in turn can leave the upper
#' limit of the interval of `z` open.
#'
#' A Beta-Binomial model whose mean is the closed-form expected response of
#' the parasitoid process (Okuyama 2026) is fitted to the data. For each
#' trial, the model is refitted without it and the two-sided predictive
#' probability of its outcome is computed (leave-one-out). A trial is flagged
#' when this probability is below `alpha / n` (Bonferroni). The screening
#' takes seconds. **It never removes data**: a flagged trial may be a
#' legitimate result of a highly variable process. Check it in the original
#' records and, if it is kept, report the analysis with and without it.
#'
#' @param data Data frame with one row per trial.
#' @param T Trial duration. Only the shape of the predicted distribution is
#'   used, so the result does not depend on the time unit.
#' @param alpha Global significance level of the screening.
#' @param dens,response Names of the columns with host density and number of
#'   hosts parasitized (prey consumed).
#' @return A data frame sorted from the least to the most probable trial,
#'   with columns `row` (row of `data`), `dens`, `par`, `expected_loo`
#'   (expected count without the trial), `p_loo` (leave-one-out
#'   probability), `rho_without` (Beta-Binomial overdispersion without the
#'   trial), `rho_change` (relative change of the overdispersion),
#'   `flagged` and `influences_dispersion` (the overdispersion falls by more
#'   than half without the trial). The overdispersion of the complete data
#'   and the threshold are stored in the attributes `rho_all` and
#'   `p_threshold`.
#' @seealso [fit_funres()]
#' @examples
#' d <- data.frame(dens = rep(c(5, 10, 20, 40), each = 6),
#'                 par  = c(2, 3, 2, 3, 2, 3, 5, 6, 5, 6, 5, 1,
#'                          9, 10, 9, 11, 10, 9, 17, 18, 16, 17, 18, 17))
#' scr <- screen_outliers(d, T = 24)
#' head(scr)
#' @export
screen_outliers <- function(data, T, alpha = 0.05, dens = "dens", response = "par") {
  if (!all(c(dens, response) %in% names(data)))
    stop("Columns '", dens, "' and '", response, "' not found in 'data'.", call. = FALSE)
  data_spp <- data.frame(dens = data[[dens]], par = data[[response]])
  .check_data_spp(data_spp)
  # The Beta-Binomial model is fitted with time measured in units of T
  # (a* = a T, h* = h / T), so that its fixed bounds and starting values do
  # not depend on the time unit of the experiment. The predicted
  # distribution, the only output used here, is the same in any unit.
  T_exp <- 1
  n <- nrow(data_spp)
  full <- .fit_bb(data_spp, T_exp)
  if (is.null(full)) stop("The Beta-Binomial screening model could not be fitted.",
                          call. = FALSE)
  rho_all <- full$par[["rho"]]

  rows <- lapply(seq_len(n), function(i) {
    f <- .fit_bb(data_spp[-i, , drop = FALSE], T_exp)
    if (is.null(f)) return(NULL)
    p <- f$par
    N <- data_spp$dens[i]; y <- data_spp$par[i]
    sh  <- .bb_shapes(.mu_parasitoid(N, p[["a"]], p[["h"]], p[["z"]], T_exp) / N, p[["rho"]])
    pmf <- exp(.dbbinom_log(0:N, N, sh$alpha, sh$beta))
    p_low  <- sum(pmf[seq_len(y + 1)])
    p_high <- sum(pmf[(y + 1):(N + 1)])
    data.frame(row = i, dens = N, par = y,
               expected_loo = sum((0:N) * pmf),
               p_loo = min(1, 2 * min(p_low, p_high)),
               rho_without = p[["rho"]],
               rho_change = p[["rho"]] / rho_all - 1)
  })
  r <- do.call(rbind, rows)
  r$flagged <- r$p_loo < alpha / n
  r$influences_dispersion <- r$rho_change < -0.5
  r <- r[order(r$p_loo), ]
  rownames(r) <- NULL
  attr(r, "rho_all")     <- rho_all
  attr(r, "p_threshold") <- alpha / n
  r
}
