# Model specification, simulation and fitting interface (API of 1.2.0)

datos <- data.frame(dens = rep(c(2, 5, 10, 20), each = 4),
                    par  = c(0, 1, 1, 2, 2, 3, 3, 4, 5, 6, 6, 7, 8, 9, 10, 12))

test_that("funres_model() builds the model and validates the input", {
  m <- funres_model()
  expect_s3_class(m, "funres_model")
  expect_equal(m$search$code, 0L)
  expect_equal(m$handling$code, 0L)
  expect_length(m$fix, 0)
  expect_equal(funres_model(search = "exponential")$fix$k, 1)
  expect_equal(funres_model(handling = "constant")$fix$s, 0)
  expect_equal(funres_model(fix = list(z = 1))$fix$z, 1)
  expect_error(funres_model(search = "beta"))
  expect_error(funres_model(fix = list(a = 1)), "Only k, s and z")
  expect_error(funres_model(search = "exponential", fix = list(k = 2)), "k = 1")
  expect_error(funres_model(handling = "constant", fix = list(s = 1)), "s = 0")
  expect_output(print(funres_model(process = "predator")), "consumed prey removed")
})

test_that("default-model likelihood equals negloglik_fixed_z() of 1.1.0", {
  d3 <- read.csv(test_path("fixtures", "D3.csv"))
  p <- c(a = 0.012, h = 0.626, z = 2.5, k = 0.902, s = 0.542)
  set.seed(88)
  new <- funresMech:::.nll_sim(p, funres_model(), d3, 24, 3000, cut = FALSE)
  expect_equal(sprintf("%.8f", new), "132.04994217")
})

test_that("simulate_response() returns counts within [0, dens] for every family", {
  par <- c(a = 0.05, h = 0.3, z = 1.5, k = 0.8, s = 0.2)
  for (m in list(funres_model(), funres_model(search = "lognormal", handling = "gamma"),
                 funres_model(search = "weibull", handling = "weibull"),
                 funres_model(search = "exponential", handling = "constant"),
                 funres_model(process = "predator"))) {
    set.seed(1)
    s <- simulate_response(m, dens = c(3, 12), T = 24, par = par, n_sim = 200)
    expect_equal(nrow(s), 400)
    expect_true(all(s$response >= 0 & s$response <= s$dens))
  }
  expect_error(simulate_response(dens = 5, T = 1, par = c(a = 1)), "Missing parameter")
})

test_that("funres_sampler() and test_sampler() check user samplers", {
  r_gamma1 <- function(n, shape) rgamma(n, shape = shape, rate = shape)
  smp <- funres_sampler(r_gamma1, lower = 0.01, upper = 100, name = "gamma (R)")
  expect_s3_class(smp, "funres_sampler")
  set.seed(2)
  chk <- test_sampler(smp, shape = c(0.5, 2), n = 20000)
  expect_true(all(chk$mean_ok))
  expect_error(funres_sampler(function(x) x), "arguments 'n' and 'shape'")
  expect_error(funres_sampler(r_gamma1, lower = 0), "0 < lower")
  m <- funres_model(search = smp, handling = funres_sampler(r_gamma1, 0.1, 50))
  expect_true(m$search$user && m$handling$user)
})

test_that("fit_funres() fits every family, fixed parameters and both processes", {
  skip_on_cran()
  ctrl <- funres_control("exploratory", n_sim = 150, itermax = 6, NP = 12, n_reeval = 2)
  models <- list(funres_model(),
                 funres_model(search = "weibull", handling = "gamma"),
                 funres_model(search = "exponential", handling = "constant"),
                 funres_model(process = "predator"),
                 funres_model(fix = list(z = 1)),
                 funres_model(search = funres_sampler(function(n, shape)
                   rgamma(n, shape = shape, rate = shape), 0.01, 100)))
  for (m in models) {
    set.seed(3)
    f <- suppressWarnings(fit_funres(datos, T = 24, model = m, z_grid = c(1, 2),
                                     control = ctrl, extend_z = FALSE,
                                     refine = FALSE, time_unit = "h"))
    expect_s3_class(f, "funres_fit")
    expect_named(coef(f), c("a", "h", "z", "k", "s"))
    expect_true(all(is.finite(f$profile$nll)))
    expect_equal(f$n_par, funresMech:::.n_free(m))
    for (nm in names(m$fix)) expect_equal(coef(f)[[nm]], m$fix[[nm]])
    expect_true(is.finite(f$aic$AIC_BB))
    expect_equal(f$bb$n_par, if (is.null(m$fix$z)) 4 else 3)
    expect_output(print(f), "Estimates")
  }
})

test_that("fit_funres() is reproducible with set.seed() and validates input", {
  skip_on_cran()
  ctrl <- funres_control("exploratory", n_sim = 150, itermax = 5, NP = 12, n_reeval = 2)
  f <- function() suppressWarnings(fit_funres(datos, 24, z_grid = c(1, 2),
                                              control = ctrl, extend_z = FALSE))
  set.seed(4); f1 <- f()
  set.seed(4); f2 <- f()
  expect_identical(f1$profile, f2$profile)
  expect_error(fit_funres(datos, T = -1), "T must be")
  expect_error(fit_funres(data.frame(x = 1, y = 2), T = 1), "not found")
  expect_error(fit_funres(data.frame(dens = 5, par = 2.5), T = 1), "whole numbers")
  expect_error(funres_control(n_sim = 10), "n_sim")
})

test_that("screening of atypical trials does not depend on the time unit", {
  skip_on_cran()
  d2 <- read.csv(test_path("fixtures", "D2.csv"))
  a <- funresMech:::screen_outliers(d2, 24)
  b <- funresMech:::screen_outliers(d2, 24 * 60)
  expect_equal(a$p_loo, b$p_loo)
  expect_equal(sum(a$flagged), 1)
})

test_that("a search bound hit at the optimum is widened automatically", {
  # Deterministic: the point fit is replaced by a stub whose optimum lies on
  # the upper bound of the encounter rate while that bound is below 5
  # encounters per trial (dimensionless scale), and inside it otherwise. The
  # real optimiser is not used, so the test does not depend on the random
  # path of DEoptim (which differs between platforms).
  fake_point <- function(z, data_spp, T_exp, model, ctrl, bounds, H_ref,
                         z_free = FALSE, initialpop = NULL) {
    on_bound <- bounds$lambda[2] < 5
    lam <- if (on_bound) bounds$lambda[2] else 3
    nm <- c("lambda", "h", "k", "s")
    list(z = z, par = c(a = lam / (T_exp * H_ref^z), h = 0.5, z = z, k = 1, s = 0.2),
         theta = NULL, nll_opt = 100 + z, nll_reeval = 100.1 + z + if (on_bound) 5 else 0,
         nll_reeval_se = 0.05, n_floor = 0L,
         at_bound = stats::setNames(c(on_bound, FALSE, FALSE, FALSE), nm),
         bound_side = stats::setNames(c(if (on_bound) "upper" else "", "", "", ""), nm),
         iter = 10, nfeval = 100, late_gain = 0, spread = stats::setNames(rep(0.1, 4), nm),
         pop = NULL)
  }
  testthat::local_mocked_bindings(.fit_point = fake_point)
  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)

  base <- list(datos, 24, z_grid = c(1, 2), extend_z = FALSE, refine = FALSE)
  # upper bound of lambda: 0.05 per h -> 1.2 encounters per trial (hit);
  # after one widening x10 -> 12 (not hit)
  ctrl <- funres_control("exploratory", bounds = list(lambda = c(1e-4, 0.05)),
                         widen = 1, bb = FALSE)
  f <- do.call(fit_funres, c(base, list(control = ctrl)))
  expect_true("lambda(upper)" %in% f$widened)
  expect_equal(f$bounds$lambda[2], 0.05 * 24 * 10)
  expect_false(any(f$profile$at_bound))
  expect_true(any(grepl("widened automatically", f$notes)))

  ctrl0 <- funres_control("exploratory", bounds = list(lambda = c(1e-4, 0.05)),
                          widen = 0, bb = FALSE)
  f0 <- do.call(fit_funres, c(base, list(control = ctrl0)))
  expect_length(f0$widened, 0)
  expect_true(all(f0$profile$at_bound))
  expect_true(any(grepl("only a limit", f0$notes)))
  r0 <- reliability(f0)
  expect_equal(r0$status[r0$check == "Estimates away from search bounds"], "fail")
})

test_that("bound helpers widen and parse correctly", {
  b <- funresMech:::.make_bounds(list(), 1)
  hits <- funresMech:::.parse_hits("lambda(upper),k(lower),s(upper)")
  expect_equal(vapply(hits, `[[`, "", "name"), c("lambda", "k", "s"))
  w <- funresMech:::.widen_bounds(b, hits)
  expect_equal(w$lambda, b$lambda * c(1, 10))
  expect_equal(w$k, b$k * c(0.1, 1))
  expect_equal(w$cv_h, b$cv_h * c(1, 10))
  expect_length(funresMech:::.parse_hits(""), 0)
})

test_that("observed outcomes never simulated are counted and flagged", {
  d <- data.frame(dens = c(10, 10), par = c(0, 10))
  m <- funres_model()
  set.seed(6)
  r <- funresMech:::.nll_sim(c(a = 1e-4, h = 0.1, z = 1, k = 50, s = 0), m, d, 1, 200,
                             detail = TRUE)
  expect_equal(r$n_floor, 1L)                  # 10 of 10 is never simulated
  expect_gt(r$nll, 700)
})

test_that("reliability() and verify_fit() report the checks with a status", {
  skip_on_cran()
  ctrl <- funres_control("exploratory", n_sim = 200, itermax = 8, NP = 12, n_reeval = 3)
  set.seed(7)
  f <- suppressWarnings(fit_funres(datos, 24, z_grid = c(1, 2), control = ctrl,
                                   extend_z = FALSE, refine = FALSE))
  r <- reliability(f)
  expect_s3_class(r, "funres_reliability")
  expect_true(all(r$status %in% c("ok", "warning", "fail")))
  expect_true("Settings suitable for reporting" %in% r$check)
  expect_true(attr(r, "overall") %in% c("warning", "fail"))   # exploratory
  set.seed(8)
  v <- suppressWarnings(verify_fit(f, n_runs = 1, effort = 1.5))
  rv <- reliability(v)
  expect_true(all(c("Independent runs agree (NLL)", "More effort does not improve the optimum",
                    "Wider bounds do not improve the optimum") %in% rv$check))
  expect_length(v$verification_runs, 3)
  expect_output(print(rv), "Reliability of the fit")
  expect_output(print(v), "Reliability:")
})

test_that("predator mean equals the random predator equation and the ODE", {
  a <- 0.3; h <- 0.02; x <- c(5, 20, 100)
  ne <- funresMech:::.mu_predator(x, a, h, 1, 1)
  expect_equal(ne, x * (1 - exp(a * (h * ne - 1))), tolerance = 1e-8)   # Rogers
  ode <- function(H0, a, h, z, T, n = 1e5) {
    H <- H0; dt <- T / n
    for (i in seq_len(n)) H <- H - dt * a * H^z / (1 + a * h * H^z)
    H0 - H
  }
  expect_equal(funresMech:::.mu_predator(20, 0.01, 0.05, 2, 1), ode(20, 0.01, 0.05, 2, 1),
               tolerance = 1e-4)
})

test_that("Beta-Binomial reference reproduces the NLL implied by Okuyama (2026, Table S1)", {
  skip_on_cran()
  # AIC(BB) = AIC(mech) + Delta: D3 272.890 + 0.62; D2 336.589 + 10.3 (4 parameters)
  for (cs in list(list(f = "D3.csv", nll = (272.890 + 0.62 - 8) / 2),
                  list(f = "D2.csv", nll = (336.589 + 10.3 - 8) / 2))) {
    d <- read.csv(test_path("fixtures", cs$f))
    set.seed(1)
    b <- funresMech:::.fit_bb_model(d, funres_model())
    expect_lt(abs(b$nll - cs$nll), 0.15)
  }
})

test_that("joint estimation of z with an adaptive profile for the interval", {
  skip_on_cran()
  ctrl <- funres_control("exploratory", n_sim = 200, itermax = 8, NP = 15,
                         n_reeval = 2, bb = FALSE)
  set.seed(9)
  f <- suppressWarnings(fit_funres(datos, 24, control = ctrl, z_method = "joint",
                                   z_step = 0.5, refine = FALSE))
  expect_equal(f$z_method, "joint")
  expect_equal(nrow(f$joint), 1)
  expect_true(f$joint$z >= 0.1 && f$joint$z <= 20)
  expect_true(1 %in% round(f$profile$z, 6))            # z = 1 evaluated for the AIC
  expect_true(f$joint$z %in% f$profile$z)
  expect_equal(coef(f)[["z"]], f$profile$z[which.min(f$profile$nll)])
  expect_false(is.na(f$aic$AIC_restricted))
  r <- reliability(f)
  expect_equal(r$status[r$check == "Estimation of z by the validated method"], "warning")
})

test_that("defaults chosen in the phase 4 comparison: DEoptim + common random numbers + polish", {
  ctl <- funres_control()
  expect_equal(ctl$method, "DEoptim")
  expect_true(ctl$crn)
  expect_equal(ctl$polish, 200L)
})

test_that("alternative optimisers, common random numbers and polish run", {
  skip_on_cran()
  base <- list(datos, 24, z_grid = 1.5, extend_z = FALSE, refine = FALSE)
  for (meth in c("CMAES", "CRS2")) {
    skip_if_not_installed(if (meth == "CMAES") "cmaes" else "nloptr")
    ctrl <- funres_control("exploratory", n_sim = 200, itermax = 8, NP = 12,
                           n_reeval = 2, method = meth, bb = FALSE)
    set.seed(10)
    f <- suppressWarnings(do.call(fit_funres, c(base, list(control = ctrl))))
    expect_true(is.finite(f$nll))
  }
  ctrl <- funres_control("exploratory", n_sim = 200, itermax = 6, NP = 12,
                         n_reeval = 2, crn = TRUE, polish = 30, bb = FALSE)
  set.seed(11)
  f <- suppressWarnings(do.call(fit_funres, c(base, list(control = ctrl))))
  expect_true(is.finite(f$nll))
  # with common random numbers the reproducibility with set.seed() is kept
  set.seed(12); f1 <- suppressWarnings(do.call(fit_funres, c(base, list(control = ctrl))))
  set.seed(12); f2 <- suppressWarnings(do.call(fit_funres, c(base, list(control = ctrl))))
  expect_identical(f1$profile, f2$profile)
})

test_that("a jagged profile (extra local minima) fails the regularity check", {
  skip_on_cran()
  ctrl <- funres_control("exploratory", n_sim = 150, itermax = 4, NP = 12, n_reeval = 2, bb = FALSE)
  set.seed(1)
  f <- suppressWarnings(fit_funres(datos, 24, z_grid = c(1, 2), control = ctrl,
                                   extend_z = FALSE, refine = FALSE))
  f$profile <- f$profile[c(1, 2, 1, 2, 1), ]
  f$profile$z <- 1:5
  f$profile$nll <- c(20, 10, 18, 5, 15)          # minima at z = 2 and z = 4
  f$profile$nll_se <- 0.1
  r <- reliability(f)
  expect_equal(r$status[r$check == "Profile of z regular (one minimum)"], "fail")
  expect_match(r$detail[r$check == "Profile of z regular (one minimum)"], "z = 2")
})

test_that("an open lower limit of z extends the grid downwards", {
  # Deterministic stub of the point fit: the profile is a parabola in z.
  make_point <- function(curv, z0) {
    function(z, data_spp, T_exp, model, ctrl, bounds, H_ref,
             z_free = FALSE, initialpop = NULL) {
      nm <- c("lambda", "h", "k", "s")
      v <- 100 + curv * (z - z0)^2
      list(z = z, par = c(a = 1 / (T_exp * H_ref^z), h = 0.5, z = z, k = 1, s = 0.2),
           theta = NULL, nll_opt = v, nll_reeval = v, nll_reeval_se = 0.01,
           n_floor = 0L, at_bound = stats::setNames(rep(FALSE, 4), nm),
           bound_side = stats::setNames(rep("", 4), nm), iter = 10, nfeval = 100,
           late_gain = 0, spread = stats::setNames(rep(0.1, 4), nm), pop = NULL)
    }
  }
  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)
  ctrl <- funres_control("exploratory", widen = 0, bb = FALSE)
  base <- list(datos, 24, z_grid = c(0.5, 1, 1.5), refine = FALSE, control = ctrl)

  # minimum at the lowest grid value; crosses the threshold between 0.1 and 0.25
  testthat::local_mocked_bindings(.fit_point = make_point(20, 0.5))
  f <- do.call(fit_funres, base)
  expect_true(f$extended)
  expect_true(all(c(0.1, 0.25) %in% f$profile$z))
  expect_false(f$ci$low_censored)
  expect_gt(f$ci$z_low, 0.1); expect_lt(f$ci$z_low, 0.25)
  expect_false(any(f$profile$z > 1.5))           # upper limit was closed

  f0 <- do.call(fit_funres, c(base, list(extend_z = FALSE)))
  expect_true(f0$ci$low_censored)
  expect_true(any(grepl("not reached within the grid \\(z <= 0.5\\)", f0$notes)))

  # flat towards z = 0: still open after the extension, not identified
  testthat::local_mocked_bindings(.fit_point = make_point(0.5, 0))
  f1 <- do.call(fit_funres, base)
  expect_true(f1$ci$low_censored)
  expect_equal(min(f1$profile$z), 0.1)
  expect_true(any(grepl("not identified \\(profile flat down to z = 0.1\\)", f1$notes)))
  expect_match(funresMech:::.low_open_text(0.5), "extend the z grid downwards")
})
