# Tests of the internal computational functions of the mechanistic model:
# the C++ simulation engine, the simulated likelihood, the search
# parameterization, the profile interval and the fit. They use funresMech:::
# to reach unexported functions.
#
# Rules: the random seed is set in the test (never inside the function under
# test), itermax, NP and n_sim are kept small to respect the CRAN time limits,
# and nothing is written to disk. The two most expensive tests (screening and
# a complete profile) are skipped on CRAN.

small_data <- data.frame(
  dens = rep(c(5, 10), each = 5),
  par  = c(1, 2, 1, 3, 2, 4, 5, 3, 6, 4)
)

# ---------------------------------------------------------------------------
# C++ engine versus the reference implementation in R
# ---------------------------------------------------------------------------

test_that("simulate_trial() returns a valid integer within range", {
  set.seed(42)
  x <- 10
  result <- funresMech:::simulate_trial(
    x = x, T = 1, a = 0.5, h = 0.1, z = 1, k = 1, s = 0
  )
  expect_type(result, "integer")
  expect_length(result, 1)
  expect_gte(result, 0)
  expect_lte(result, x)
})

test_that("simulate_trial() returns 0 at density 0", {
  set.seed(42)
  result <- funresMech:::simulate_trial(
    x = 0, T = 1, a = 0.5, h = 0.1, z = 1, k = 1, s = 0
  )
  expect_equal(result, 0)
})

test_that("simulate_trial() is stochastic: different seeds give different results", {
  args <- list(x = 20, T = 5, a = 0.3, h = 0.2, z = 1, k = 1, s = 0.3)
  set.seed(1)
  r1 <- do.call(funresMech:::simulate_trial, args)
  set.seed(2)
  r2 <- do.call(funresMech:::simulate_trial, args)
  set.seed(3)
  r3 <- do.call(funresMech:::simulate_trial, args)
  expect_true(length(unique(c(r1, r2, r3))) > 1)
})

test_that("the C++ engine reproduces simulate_trial() exactly with the same seed", {
  # Without the saturation stop both consume the same random numbers in the
  # same order, so the results must be identical (with s = 0 and with s > 0)
  cases <- list(
    list(x = 30, T = 2,  a = 0.5,   h = 0.2, z = 1.5, k = 0.5,  s = 0),
    list(x = 8,  T = 24, a = 0.01,  h = 0.6, z = 3,   k = 0.05, s = 0.8),
    list(x = 50, T = 24, a = 0.001, h = 0.5, z = 2.5, k = 0.04, s = 0.6)
  )
  for (cs in cases) {
    n <- 100L
    set.seed(7)
    engine <- funresMech:::sim_engine_cpp(
      as.integer(cs$x), cs$a, cs$h, cs$z, cs$k, cs$s, cs$T, n,
      stop_saturated = FALSE
    )
    set.seed(7)
    ref <- vapply(seq_len(n), function(i) {
      funresMech:::simulate_trial(cs$x, cs$T, cs$a, cs$h, cs$z, cs$k, cs$s)
    }, integer(1))
    expect_identical(as.integer(engine), as.integer(ref))
  }
})

test_that("stopping saturated trials does not change the simulated distribution", {
  # Parameters under which trials often saturate (all hosts parasitized)
  n <- 20000L
  set.seed(11)
  with_stop <- funresMech:::sim_engine_cpp(6L, 2, 0.05, 1, 1, 0.05, 5, n,
                                        stop_saturated = TRUE)
  set.seed(12)
  without_stop <- funresMech:::sim_engine_cpp(6L, 2, 0.05, 1, 1, 0.05, 5, n,
                                        stop_saturated = FALSE)
  tab <- rbind(tabulate(with_stop + 1L, nbins = 7), tabulate(without_stop + 1L, nbins = 7))
  tab <- tab[, colSums(tab) > 0, drop = FALSE]
  expect_gt(suppressWarnings(chisq.test(tab)$p.value), 0.001)
  expect_gt(mean(with_stop == 6L), 0.1)   # the scenario does saturate
})

test_that("stopping saturated trials avoids runaway loops (h ~ 0, k ~ 0)", {
  # Without the stop, 10 trials at H = 50 with h = k = 1e-6 need tens of
  # millions of loop iterations; with the stop, a few thousand
  loops <- funresMech:::count_loops_cpp(50L, 1, 1e-6, 1, 1e-6, 0, 24, 10L,
                                             stop_saturated = TRUE)
  expect_lt(loops, 1e5)
  set.seed(1)
  r <- funresMech:::sim_engine_cpp(50L, 1, 1e-6, 1, 1e-6, 0, 24, 10L, TRUE)
  expect_true(all(r == 50L))
})

test_that("the engine rejects invalid parameters with a clear error", {
  expect_error(funresMech:::sim_engine_cpp(0L, 1, 1, 1, 1, 0, 1, 10L), "H must be")
  expect_error(funresMech:::sim_engine_cpp(5L, -1, 1, 1, 1, 0, 1, 10L), "invalid parameters")
  expect_error(funresMech:::sim_engine_cpp(5L, 1, 0, 1, 1, 0.5, 1, 10L), "h must be")
})

# ---------------------------------------------------------------------------
# Simulated likelihood
# ---------------------------------------------------------------------------

test_that("simulate_distribution() returns valid probabilities", {
  set.seed(3)
  p <- funresMech:::simulate_distribution(x = 10, T = 1, a = 0.5, h = 0.1,
                                          z = 1, k = 1, s = 0.05, n_sim = 500)
  expect_length(p, 11)
  expect_named(p, as.character(0:10))
  expect_true(all(p > 0))
  expect_equal(sum(p), 1, tolerance = 1e-6)
})

test_that("negloglik_fixed_z() returns a finite numeric scalar", {
  set.seed(42)
  nll <- funresMech:::negloglik_fixed_z(
    par_vec = c(a = 0.3, h = 0.1, k = 1, s = 0.1),
    z_fixed = 1, data_spp = small_data, T = 1, n_sim = 50
  )
  expect_type(nll, "double")
  expect_length(nll, 1)
  expect_true(is.finite(nll))
})

test_that("negloglik_fixed_z() returns a large value outside the valid domain", {
  data_spp <- data.frame(dens = c(5, 10), par = c(1, 2))
  nll <- funresMech:::negloglik_fixed_z(
    par_vec = c(a = -1, h = 0.1, k = 1, s = 0.1),
    z_fixed = 1, data_spp = data_spp, T = 1, n_sim = 50
  )
  expect_true(is.infinite(nll) || nll >= 1e9)
  # more hosts parasitized than offered: impossible data
  nll2 <- funresMech:::negloglik_fixed_z(
    par_vec = c(a = 0.3, h = 0.1, k = 1, s = 0.1),
    z_fixed = 1, data_spp = data.frame(dens = 5, par = 7), T = 1, n_sim = 50
  )
  expect_gte(nll2, 1e9)
})

test_that("negloglik_fixed_z() responds to the parameters", {
  set.seed(42)
  nll_1 <- funresMech:::negloglik_fixed_z(
    par_vec = c(a = 0.2, h = 0.1, k = 1, s = 0.1),
    z_fixed = 1, data_spp = small_data, T = 1, n_sim = 50
  )
  nll_2 <- funresMech:::negloglik_fixed_z(
    par_vec = c(a = 1.5, h = 0.4, k = 3, s = 0.4),
    z_fixed = 2, data_spp = small_data, T = 1, n_sim = 50
  )
  expect_false(isTRUE(all.equal(nll_1, nll_2)))
})

test_that("negloglik_fixed_z() is reproducible with the seed set by the user", {
  f <- function() funresMech:::negloglik_fixed_z(
    c(a = 0.3, h = 0.1, k = 1, s = 0.1), z_fixed = 1,
    data_spp = small_data, T = 1, n_sim = 200)
  set.seed(99); n1 <- f()
  set.seed(99); n2 <- f()
  expect_identical(n1, n2)
})

# ---------------------------------------------------------------------------
# Dimensionless search parameterization
# ---------------------------------------------------------------------------

test_that("the dimensionless parameterization round-trips exactly", {
  H_ref <- funresMech:::.h_ref(data.frame(dens = c(2, 10, 50)))
  expect_equal(H_ref, exp(mean(log(c(2, 10, 50)))))
  b <- funresMech:::.make_bounds(list(), T = 24)
  for (m in list(funresMech::funres_model(),
                 funresMech::funres_model(search = "weibull", handling = "gamma"))) {
    for (zf in c(FALSE, TRUE)) {
      ptab <- funresMech:::.par_table(m, b, z_free = zf)
      for (z in c(0.5, 2.5, 4.5, 20)) {
        par <- c(a = 3.9e-5, h = 0.626, z = z, k = 0.026, s = 0.669)
        th <- funresMech:::.nat_to_theta(par, ptab, m, H_ref, 24)
        back <- funresMech:::.theta_to_nat(th, ptab, m, z, H_ref, 24)
        expect_equal(unname(back), unname(par), tolerance = 1e-12)
      }
    }
  }
})

test_that("the search bounds do not depend on the time unit", {
  # the same experiment in hours and in minutes gives the same theta
  H_ref <- 15.4
  m <- funresMech::funres_model()
  p_h   <- c(a = 0.012, h = 0.626, z = 2.5, k = 0.9, s = 0.54)          # hours
  p_min <- c(a = 0.012 / 60, h = 0.626 * 60, z = 2.5, k = 0.9, s = 0.54 * 60)
  ptab <- funresMech:::.par_table(m, funresMech:::.make_bounds(list(), 24))
  expect_equal(funresMech:::.nat_to_theta(p_h, ptab, m, H_ref, 24),
               funresMech:::.nat_to_theta(p_min, ptab, m, H_ref, 24 * 60))
  # user bounds in natural units
  b1 <- funresMech:::.make_bounds(list(h = c(0.01, 2)), T = 24)
  b2 <- funresMech:::.make_bounds(list(h = c(0.6, 120)), T = 24 * 60)
  expect_equal(b1$h, b2$h)
  expect_error(funresMech:::.make_bounds(list(foo = c(1, 2)), 1), "Unknown bound")
})

test_that("the default bounds cover the estimates for datasets D2 and D3", {
  b <- funresMech:::.make_bounds(list(), T = 24)
  m <- funresMech::funres_model()
  ptab <- funresMech:::.par_table(m, b)
  H_ref <- 15.4
  for (p in list(c(a = 0.001, h = 0.606, z = 3.1, k = 0.037, s = 0.682),
                 c(a = 0.012, h = 0.626, z = 2.5, k = 0.902, s = 0.542),
                 c(a = 1e-3 / H_ref^5, h = 0.6, z = 6, k = 0.05, s = 0.9))) {
    th <- funresMech:::.nat_to_theta(p, ptab, m, H_ref, 24)
    expect_true(all(th >= ptab$lower & th <= ptab$upper))
  }
})

test_that("fixed parameters are removed from the vector searched by the optimizer", {
  b <- funresMech:::.make_bounds(list(), T = 1)
  expect_equal(funresMech:::.par_table(funresMech::funres_model(), b)$name,
               c("lambda", "h", "k", "s"))
  expect_equal(funresMech:::.par_table(
    funresMech::funres_model(search = "exponential", handling = "constant"), b)$name,
    c("lambda", "h"))
  expect_equal(funresMech:::.par_table(funresMech::funres_model(), b, z_free = TRUE)$name,
               c("lambda", "h", "k", "s", "z"))
  expect_equal(funresMech:::.n_free(funresMech::funres_model(fix = list(z = 1))), 4L)
})

# ---------------------------------------------------------------------------
# Likelihood-ratio interval
# ---------------------------------------------------------------------------

test_that("the 95% interval threshold is qchisq(0.95, 1) / 2, not qchisq(0.95, 1)", {
  expect_equal(funresMech:::.ci95_threshold(), 1.920729, tolerance = 1e-6)
})

test_that(".profile_ci() recovers the exact interval of a quadratic profile", {
  # NLL(z) = 0.5 * ((z - 3) / 0.5)^2 crosses the threshold at 3 +/- 0.5 * 1.95996
  z <- seq(0.5, 6, by = 0.01)
  nll <- 100 + 0.5 * ((z - 3) / 0.5)^2
  ci <- funresMech:::.profile_ci(z, nll)
  expect_equal(ci$z_hat, 3)
  expect_equal(ci$z_low,  3 - 0.5 * qnorm(0.975), tolerance = 1e-3)
  expect_equal(ci$z_high, 3 + 0.5 * qnorm(0.975), tolerance = 1e-3)
  expect_false(ci$low_censored)
  expect_false(ci$high_censored)
})

test_that(".profile_ci() reports an open limit when the profile does not cross the threshold", {
  z <- seq(0.5, 6, by = 0.5)
  nll <- 100 + 0.5 * pmin(z - 2, 0)^2 * 4 + 0.02 * pmax(z - 2, 0)   # flat to the right
  ci <- funresMech:::.profile_ci(z, nll)
  expect_true(ci$high_censored)
  expect_true(is.na(ci$z_high))
  expect_false(ci$low_censored)
  expect_true(is.finite(ci$z_low))
})

test_that("an open upper limit is called 'not identified' only when the grid reached its cap", {
  expect_match(funresMech:::.high_open_text(5), "not reached within the grid")
  expect_match(funresMech:::.high_open_text(20), "not identified")
})

test_that(".refine_points() adds points only inside wide brackets", {
  # grid 1..5 in steps of 0.5 plus a jump to 8: the upper limit (5.12) falls in (5, 8)
  z <- c(1, 1.5, 2, 2.5, 3, 4, 5, 8)
  new_z <- funresMech:::.refine_points(z, z_low = 1.53, z_high = 5.12, step0 = 0.5)
  expect_equal(new_z, c(5.5, 6, 6.5, 7, 7.5))
  # open limits (NA): nothing to refine
  expect_length(funresMech:::.refine_points(z, NA_real_, NA_real_, 0.5), 0)
})

test_that(".refine_points() refines a narrow interval down to about 1/4 of its width", {
  z <- seq(0.5, 5, by = 0.5)
  # interval [1.53, 2.2]: width 0.67, step 0.15; the brackets (1.5, 2) and (2, 2.5) are subdivided
  new_z <- funresMech:::.refine_points(z, z_low = 1.53, z_high = 2.2, step0 = 0.5)
  expect_true(length(new_z) > 0)
  expect_true(all(new_z > 1.5 & new_z < 2.5))
  expect_true(max(diff(sort(c(z, new_z))[sort(c(z, new_z)) >= 1.5 & sort(c(z, new_z)) <= 2.5])) <= 0.31)
  # also around the minimum
  near <- funresMech:::.refine_points(z, z_low = 1.53, z_high = 2.2, step0 = 0.5, z_hat = 3)
  expect_true(any(near > 2.5 & near < 3.5))
  # never below 0.05, and with a cap on the number of new points
  expect_true(all(diff(sort(c(z, funresMech:::.refine_points(z, 1.99, 2.01, 0.5)))) >= 0.05 - 1e-9))
  expect_lte(length(funresMech:::.refine_points(z, 1.5, 4.5, 0.5, max_new = 3, z_hat = 3)), 3)
})

test_that(".pchip() interpolates the points and never dips below the minimum of a profile with a flat segment", {
  z <- c(0.75, 0.8, 0.85, 0.9)
  nll <- c(167.5, 123.83, 123.77, 158.4)              # the case in which monoH.FC dipped below the minimum
  f <- funresMech:::.pchip(z, nll)
  expect_equal(f(z), nll)
  zz <- seq(0.75, 0.9, by = 0.005)
  expect_gte(min(f(zz)), min(nll) - 1e-9)
  expect_lte(max(f(zz)), max(nll) + 1e-9)
  ci <- funresMech:::.profile_ci(c(0.5, z, 1), c(172, nll, 128), interpolate = "spline")
  expect_true(ci$z_low <= ci$z_hat && ci$z_hat <= ci$z_high)
})

test_that(".profile_ci() with spline interpolation stays within the points and matches the exact limits", {
  z <- seq(0.5, 5, by = 0.5)
  nll <- 100 + (z - 2.3)^2 * 1.3                      # quadratic: limits at 2.3 +/- sqrt(1.92 / 1.3)
  lin <- funresMech:::.profile_ci(z, nll)
  spl <- funresMech:::.profile_ci(z, nll, interpolate = "spline")
  exact <- 2.3 + c(-1, 1) * sqrt(funresMech:::.ci95_threshold() / 1.3)
  expect_lt(max(abs(c(spl$z_low, spl$z_high) - exact)), 0.12)   # same order as linear interpolation with a step of 0.5
  expect_equal(spl$z_hat, 2.5)
  expect_equal(spl$min_nll, min(nll))
  # a profile that does not cross stays open
  flat <- funresMech:::.profile_ci(z, rep(100, length(z)), interpolate = "spline")
  expect_true(flat$low_censored && flat$high_censored)
})

test_that(".ci_profile() leaves out the points at the probability floor", {
  prof <- data.frame(z = seq(1, 5, by = 1), nll = c(112, 101, 100, 101.5, 740),
                     n_floor = c(0L, 0L, 0L, 0L, 1L))
  ci <- funresMech:::.ci_profile(prof)
  # without the floored point the profile does not cross the threshold above (z = 5 is left out)
  expect_true(ci$high_censored)
  expect_equal(ci$z_hat, 3)
  # when every point is floored they cannot be left out, so all are used
  prof$n_floor[] <- 1L
  expect_false(is.null(funresMech:::.ci_profile(prof)))
})

test_that("floor_retry repeats the re-evaluations with more simulations when the floor is hit", {
  ctrl0 <- funres_control(n_sim = 2000, itermax = 50, n_reeval = 3, floor_retry = 0)
  expect_equal(ctrl0$floor_retry, 0L)
  expect_equal(funres_control()$floor_retry, 2L)
  calls <- integer(0)
  fake_nll <- function(par, model, data_spp, T, n_sim, cut = TRUE, detail = FALSE) {
    if (!detail) return(100)
    calls <<- c(calls, n_sim)
    # the floor is hit only with a small n_sim; four times more resolves it
    list(nll = 100, n_floor = if (n_sim < 4 * 300) 1L else 0L)
  }
  testthat::local_mocked_bindings(.nll_sim = fake_nll, .package = "funresMech")
  d <- data.frame(dens = rep(c(5, 10), each = 3), par = c(1, 2, 2, 3, 4, 4))
  ctrl <- funres_control(n_sim = 300, itermax = 3, NP = 8, n_reeval = 2, polish = 0)
  p <- funresMech:::.fit_point(1.5, d, 24, funres_model(), ctrl, funresMech:::.make_bounds(list(), 24),
                               funresMech:::.h_ref(d))
  expect_equal(p$n_floor, 0L)
  expect_equal(p$n_retry, 1L)
  expect_true(1200 %in% calls)
})

test_that(".nll_at_z() interpolates linearly and returns NA outside the grid", {
  z <- c(0.5, 1.5, 2.5)
  nll <- c(10, 6, 8)
  expect_equal(funresMech:::.nll_at_z(z, nll, 1.5), 6)
  expect_equal(funresMech:::.nll_at_z(z, nll, 1), 8)
  expect_true(is.na(funresMech:::.nll_at_z(z, nll, 3)))
})

# ---------------------------------------------------------------------------
# Fitting
# ---------------------------------------------------------------------------

test_that(".fit_point() returns the expected structure", {
  # small itermax and NP: CRAN limits the check time.
  set.seed(42)
  m <- funresMech::funres_model()
  ctrl <- funresMech::funres_control(n_sim = 100, itermax = 10,
                                     NP = 20, n_reeval = 2)
  fit <- suppressWarnings(funresMech:::.fit_point(
    z = 1, data_spp = small_data, T_exp = 1, model = m, ctrl = ctrl,
    bounds = funresMech:::.make_bounds(list(), 1),
    H_ref = funresMech:::.h_ref(small_data)))
  expect_type(fit, "list")
  expect_named(fit$par, c("a", "h", "z", "k", "s"))
  expect_true(all(fit$par[c("a", "h", "k")] > 0))
  expect_gte(fit$par[["s"]], 0)
  expect_true(is.finite(fit$nll_opt))
  expect_true(is.finite(fit$nll_reeval))
  expect_length(fit$at_bound, 4)
  expect_equal(fit$iter, 10)
})

test_that("fit_profile() returns the profile, z_hat, interval and AIC", {
  skip_on_cran()
  set.seed(42)
  fit <- suppressWarnings(funresMech:::fit_profile(
    data_spp = small_data, T_exp = 1, z_grid = c(0.5, 1, 2),
    n_sim = 100, itermax = 10, NP = 20, reltol = 1e-2, n_reeval = 2,
    extend_z = FALSE
  ))
  expect_named(fit, c("profile", "par", "nll", "aic", "ci", "extended", "refined", "notes"))
  # the profile NLL is the mean of the re-evaluations, not the optimizer minimum
  expect_true(all(c("nll", "nll_se", "nll_opt") %in% names(fit$profile)))
  expect_named(fit$par, c("a", "h", "z", "k", "s"))
  expect_equal(nrow(fit$profile), 3)
  expect_true(all(is.finite(fit$profile$nll)))
  # z_hat is the minimum of the profile and the reported NLL is its value
  expect_equal(fit$par[["z"]], fit$profile$z[which.min(fit$profile$nll)])
  expect_equal(fit$nll, min(fit$profile$nll))
  expect_equal(fit$aic$AIC_full, 2 * fit$nll + 10)
  expect_false(is.na(fit$aic$AIC_restricted))   # z = 1 is in the grid
})

test_that("fit_profile() adds z = 1 to the grid and validates the data", {
  skip_on_cran()
  set.seed(1)
  fit <- suppressWarnings(funresMech:::fit_profile(
    data_spp = small_data, T_exp = 1, z_grid = c(0.5, 1.5),
    n_sim = 100, itermax = 5, NP = 12, n_reeval = 1, extend_z = FALSE))
  expect_true(1 %in% fit$profile$z)
  expect_error(funresMech:::fit_profile(
    data_spp = data.frame(dens = 5, par = 9), T_exp = 1, z_grid = 1),
    "between 0 and the host density")
})

# ---------------------------------------------------------------------------
# Screening for atypical trials
# ---------------------------------------------------------------------------

test_that("screen_outliers() flags an atypical trial and does not flag consistent data", {
  skip_on_cran()
  # D2 (Kishani Farahani & Goldansaz 2013, FoRAGE): the only atypical trial
  # is 2 of 10 hosts parasitized (neighboring trials: 85-95%)
  d2 <- data.frame(
    dens = c(50, 50, 2, 2, 2, 4, 4, 6, 6, 6, 8, 8, 10, 10, 10, 15, 15, 15, 15,
             20, 20, 20, 20, 20, 25, 25, 25, 25, 25, 25, 30, 30, 30, 30, 30, 30,
             30, 30, 30, 30, 35, 35, 35, 35, 35, 35, 35, 35, 40, 40, 40, 40, 40,
             40, 40, 40, 40, 45, 45, 45, 45, 45, 45, 45, 45, 45, 50, 50, 50, 50,
             50, 50, 50),
    par  = c(31, 32, 0, 1, 2, 3, 4, 4, 5, 6, 7, 8, 2, 8, 10, 12, 13, 14, 15,
             15, 16, 17, 18, 19, 15, 18, 19, 20, 21, 22, 16, 17, 18, 19, 20, 21,
             22, 23, 24, 20, 21, 22, 23, 24, 25, 26, 27, 29, 21, 23, 24, 25, 26,
             27, 28, 30, 31, 20, 22, 23, 24, 25, 27, 28, 29, 30, 25, 26, 27, 28,
             29, 30, 33)
  )
  expect_equal(nrow(d2), 73)
  scr <- screen_outliers(d2, T = 24)
  expect_equal(sum(scr$flagged), 1)
  top <- scr[1, ]
  expect_equal(c(top$dens, top$par), c(10, 2))
  expect_true(top$influences_dispersion)

  # data generated without atypical trials
  set.seed(5)
  N <- rep(c(5, 10, 20, 40), each = 6)
  y <- rbinom(length(N), N, 0.5)
  scr0 <- screen_outliers(data.frame(dens = N, par = y), T = 24)
  expect_equal(sum(scr0$flagged), 0)
})

# ---------------------------------------------------------------------------
# Report template
# ---------------------------------------------------------------------------

test_that("the report template is installed and accessible", {
  path <- system.file("app", "report_template.Rmd", package = "funresMech")
  expect_true(nzchar(path))
  expect_true(file.exists(path))
})

test_that(".workers_ok() is TRUE under a sequential plan", {
  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)
  expect_true(funresMech:::.workers_ok())
})

test_that("the complete fit does not depend on the time unit (hours versus minutes)", {
  d <- data.frame(dens = rep(c(3, 6, 12, 24), each = 4),
                  par  = c(0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 6, 7, 6, 8, 10, 12))
  ctrl <- funres_control(n_sim = 300, itermax = 6, NP = 12, n_reeval = 2,
                         polish = 20, bb = FALSE, widen = 0)
  fit_in <- function(T) {
    set.seed(123)
    fit_funres(d, T = T, z_grid = c(1, 2, 3), control = ctrl, extend_z = FALSE, refine = FALSE)
  }
  h <- fit_in(24); m <- fit_in(24 * 60)
  expect_equal(h$profile$nll, m$profile$nll, tolerance = 1e-8)
  expect_equal(h$par[["z"]], m$par[["z"]])
  expect_equal(unname(m$par[["h"]]), unname(h$par[["h"]]) * 60, tolerance = 1e-6)
  expect_equal(unname(m$par[["a"]]), unname(h$par[["a"]]) / 60, tolerance = 1e-6)
  expect_equal(h$profile$k, m$profile$k, tolerance = 1e-8)
  expect_equal(unname(m$par[["s"]]), unname(h$par[["s"]]) * 60, tolerance = 1e-6)
})
