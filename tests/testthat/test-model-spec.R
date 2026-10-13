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
