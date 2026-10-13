############################################################
# Beta-Binomial reference model for the AIC comparison
############################################################

# Okuyama (2026) compares the simulation-based model with a Beta-Binomial
# model whose mean is the expected number of victims under the same
# functional response f(x) = a x^z / (1 + a h x^z):
#   parasitoid process: x {1 - exp(-f(x) T / x)}                (Eqs. 1-2)
#   predator process:   prey consumed from dH/dt = -f(H)        (Appendix S1,
#                       Script S2), solved for N_e = x - H(T)
# The model is fitted with time in units of T (a* = a T, h* = h / T), on a
# log/logit scale, by differential evolution followed by a local polish.
# The objective is deterministic.

# Expected prey consumed by one predator in time T (vectorised over x)
.mu_predator <- function(x, a, h, z, T = 1) {
  vapply(x, function(H0) {
    if (H0 <= 0) return(0)
    g <- function(Ne) {
      Ht <- H0 - Ne
      first <- if (abs(z - 1) < 1e-8) log(Ht / H0) / a
               else (Ht^(1 - z) - H0^(1 - z)) / (a * (1 - z))
      first - h * Ne + T
    }
    upper <- H0 * (1 - 1e-12)
    gu <- g(upper)
    if (is.nan(gu)) stop("predator mean not defined")
    if (gu >= 0) return(H0)                      # all prey consumed within T
    stats::uniroot(g, c(0, upper), f.lower = T, f.upper = gu, tol = 1e-10)$root
  }, numeric(1))
}

.mu_process <- function(process, x, a, h, z, T = 1) {
  if (process == "predator") .mu_predator(x, a, h, z, T) else .mu_okuyama(x, a, h, z, T)
}

# NLL of the Beta-Binomial model; theta = (log a*, log h*, [log z], logit rho)
.nll_bb_theta <- function(theta, data_spp, process, z_fixed = NULL) {
  a <- exp(theta[1]); h <- exp(theta[2])
  if (is.null(z_fixed)) { z <- exp(theta[3]); rho <- stats::plogis(theta[4]) }
  else { z <- z_fixed; rho <- stats::plogis(theta[3]) }
  ux <- unique(data_spp$dens)
  mu <- tryCatch(.mu_process(process, ux, a, h, z, 1), error = function(e) NULL)
  if (is.null(mu) || any(!is.finite(mu))) return(1e10)
  pr <- (mu / ux)[match(data_spp$dens, ux)]
  sh <- .bb_shapes(pr, rho)
  ll <- .dbbinom_log(data_spp$par, data_spp$dens, sh$alpha, sh$beta)
  if (any(!is.finite(ll))) return(1e10)
  -sum(ll)
}

# Fit of the reference Beta-Binomial model for a funres_model()
.fit_bb_model <- function(data_spp, model, itermax = 150, NP = 40) {
  z_fixed <- model$fix$z
  lower <- c(log(1e-8), log(1e-6), if (is.null(z_fixed)) log(0.1), -14)
  upper <- c(log(1e4),  log(1),    if (is.null(z_fixed)) log(20),   7)
  fn <- function(th) .nll_bb_theta(th, data_spp, model$process, z_fixed)
  de <- withCallingHandlers(
    DEoptim(fn, lower, upper, DEoptim.control(itermax = itermax, NP = NP, trace = FALSE)),
    warning = function(w) {
      if (grepl("ten times the length of the parameter vector", conditionMessage(w)))
        invokeRestart("muffleWarning")
    })
  pol <- tryCatch(stats::optim(de$optim$bestmem, fn, method = "Nelder-Mead",
                               control = list(maxit = 5000, reltol = 1e-12)),
                  error = function(e) NULL)
  best <- if (!is.null(pol) && pol$value < de$optim$bestval) pol$par else de$optim$bestmem
  val <- min(de$optim$bestval, if (!is.null(pol)) pol$value else Inf)
  best <- unname(best)
  par <- c(a_star = exp(best[1]), h_star = exp(best[2]),
           z = if (is.null(z_fixed)) exp(best[3]) else z_fixed,
           rho = stats::plogis(best[length(best)]))
  list(nll = val, par = par, n_par = length(lower))
}
