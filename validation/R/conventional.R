# ============================================================================
# conventional.R
#
# Conventional functional-response fits, shared by 22_conventional_methods.R
# (simulated data) and 25_real_data.R (published data): three mean functions
# (holling, rogers, eq2) and three error models (binomial and beta-binomial
# maximum likelihood, least squares), fitted with time in units of the trial
# duration T. See 22_conventional_methods.R for the definitions.
# ============================================================================

# ---- mean functions (time in units of T) ------------------------------------------
# Principal branch of the Lambert W function from the log of its argument
# (W(x) + log W(x) = log x), vectorized and safe for very large arguments.
lambertW_log <- function(lx) {
  w <- exp(pmin(lx, 0))                       # W(x) ~ x for small x (exact enough below exp(-30))
  big <- is.finite(lx) & lx > -30
  if (any(big)) {
    l <- lx[big]
    v <- ifelse(l < 1, exp(l), l - log(pmax(l, 1)))
    for (i in 1:100) {
      v_new <- v - (v + log(v) - l) / (1 + 1 / v)
      v_new <- ifelse(!is.finite(v_new) | v_new <= 0, v / 10, v_new)
      done <- max(abs(v_new - v) / v) < 1e-13
      v <- v_new
      if (done) break
    }
    w[big] <- v
  }
  w
}
mean_rogers <- function(N, a, h, z) {
  ae <- a * N^(z - 1)                                   # attack rate b N0^q
  lx <- log(ae * h * N) - ae * (1 - h * N)
  N - lambertW_log(lx) / (ae * h)
}
mean_eq2 <- function(N, a, h, z) funresMech:::.mu_parasitoid(N, a, h, z, 1)
mean_holling <- function(N, a, h, z) a * N^z / (1 + a * h * N^z)
MEANS <- list(holling = mean_holling, rogers = mean_rogers, eq2 = mean_eq2)

LO <- c(log(1e-8), log(1e-6), log(0.1), -14)
UP <- c(log(1e4), log(1), log(20), 7)

# Objective for theta = (log a*, log h*, [log z], [logit rho]); z fixed if z_fix given
objective <- function(theta, d, mean_fun, error, z_fix = NULL) {
  a <- exp(theta[1]); h <- exp(theta[2])
  z <- if (is.null(z_fix)) exp(theta[3]) else z_fix
  mu <- mean_fun(d$dens, a, h, z)
  if (any(!is.finite(mu))) return(1e10)
  mu <- pmin(pmax(mu, 0), d$dens)
  if (error == "ls") return(sum((d$par - mu)^2))
  p <- pmin(pmax(mu / d$dens, 1e-12), 1 - 1e-12)
  if (error == "binomial") return(-sum(stats::dbinom(d$par, d$dens, p, log = TRUE)))
  rho <- stats::plogis(theta[length(theta)])
  sh <- funresMech:::.bb_shapes(p, rho)
  ll <- funresMech:::.dbbinom_log(d$par, d$dens, sh$alpha, sh$beta)
  if (any(!is.finite(ll))) 1e10 else -sum(ll)
}

par_index <- function(error, z_free) c(1, 2, if (z_free) 3, if (error == "bb") 4)

fit_one <- function(d, mean_fun, error, z_fix = NULL) {
  idx <- par_index(error, is.null(z_fix))
  fn <- function(th) objective(th, d, mean_fun, error, z_fix)
  de <- suppressWarnings(DEoptim::DEoptim(fn, LO[idx], UP[idx],
                                          DEoptim::DEoptim.control(itermax = 150, NP = 40, trace = FALSE)))
  best <- de$optim$bestmem; val <- de$optim$bestval
  pol <- tryCatch(stats::optim(best, fn, method = "Nelder-Mead",
                               control = list(maxit = 5000, reltol = 1e-12)), error = function(e) NULL)
  if (!is.null(pol) && pol$value < val) { best <- pol$par; val <- pol$value }
  list(theta = unname(best), value = val, idx = idx)
}

# Value of the objective with z fixed, starting from a previous optimum
profile_value <- function(z, d, mean_fun, error, start) {
  fn <- function(th) objective(th, d, mean_fun, error, z)
  r <- stats::optim(start, fn, method = "Nelder-Mead", control = list(maxit = 3000, reltol = 1e-10))
  list(value = r$value, theta = r$par)
}

analyse <- function(d, mean_name, error, z_true = NA_real_, T = T_EXP) {
  mean_fun <- MEANS[[mean_name]]
  n <- nrow(d)
  f <- fit_one(d, mean_fun, error)
  th <- f$theta
  z_hat <- exp(th[3]); best <- f$value
  n_par <- length(f$idx)
  # statistic of a value of the objective relative to the best one
  stat <- function(v) if (error == "ls") (v - best) / (best / (n - n_par)) else 2 * (v - best)
  thr <- if (error == "ls") stats::qf(0.95, 1, n - n_par) else stats::qchisq(0.95, 1)
  # fit with z = 1 (test of z against 1)
  f1 <- fit_one(d, mean_fun, error, z_fix = 1)
  if (f1$value < best) f1$value <- best                     # z = 1 cannot beat the free fit
  p_z1 <- if (error == "ls") stats::pf(stat(f1$value), 1, n - n_par, lower.tail = FALSE)
          else stats::pchisq(stat(f1$value), 1, lower.tail = FALSE)
  # profile of z on a grid around the estimate, then the crossings by root finding
  start_other <- th[-3]
  prof_at <- function(z, start) profile_value(z, d, mean_fun, error, start)
  side <- function(direction) {
    zs <- exp(log(z_hat) + direction * seq(0.05, 6, by = 0.15))
    zs <- zs[zs >= 0.1 & zs <= 20]
    st <- start_other; prev_z <- z_hat; prev_s <- 0
    for (z in zs) {
      r <- prof_at(z, st); st <- r$theta
      if (r$value < best) { best <<- r$value }
      s <- stat(r$value)
      if (s >= thr) {
        g <- function(lz) stat(prof_at(exp(lz), st)$value) - thr
        root <- tryCatch(stats::uniroot(g, sort(c(log(prev_z), log(z))), tol = 1e-4)$root,
                         error = function(e) log(z))
        return(c(limit = exp(root), open = 0))
      }
      prev_z <- z; prev_s <- s
    }
    c(limit = NA_real_, open = 1)
  }
  lo <- side(-1); hi <- side(+1)
  a_hat <- exp(th[1]); h_hat <- exp(th[2])
  data.frame(mean = mean_name, error = error,
             z_hat = z_hat, a_hat = a_hat / T, h_hat = h_hat * T,
             rho_hat = if (error == "bb") stats::plogis(th[4]) else NA_real_,
             objective = f$value, n_par = n_par,
             nll = if (error == "ls") NA_real_ else f$value,
             z_low = unname(lo["limit"]), z_high = unname(hi["limit"]),
             low_open = unname(lo["open"]) == 1, high_open = unname(hi["open"]) == 1,
             p_z1 = p_z1,
             covers = if (is.na(z_true)) NA else
               (is.na(lo["limit"]) || lo["limit"] <= z_true) && (is.na(hi["limit"]) || z_true <= hi["limit"]),
             excludes_1 = (!is.na(lo["limit"]) && lo["limit"] > 1) || (!is.na(hi["limit"]) && hi["limit"] < 1),
             stringsAsFactors = FALSE)
}

