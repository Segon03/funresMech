# Convention for a host or prey still being handled at the end of the trial
# (funres_model(unfinished = ...)). The default ("include") is checked bit by
# bit against version 1.1.0 in test-regression.R and test-engine-families.R.

sim <- function(...) funresMech:::sim_trials_cpp(...)

test_that("funres_model() stores and prints the convention", {
  expect_equal(funres_model()$unfinished, "include")
  expect_equal(funres_model(unfinished = "exclude")$unfinished, "exclude")
  expect_error(funres_model(unfinished = "maybe"))
  expect_output(print(funres_model()), "still in handling at T is counted")
  expect_output(print(funres_model(unfinished = "exclude")), "is not counted")
  old <- funres_model(); old$unfinished <- NULL         # object saved before the option existed
  expect_false(funresMech:::.excludes_unfinished(old))
})

test_that("exclude reproduces the R reference implementation exactly", {
  for (cs in list(list(H = 5, a = 0.05, h = 2, z = 1, k = 0.7, s = 1.2, T = 24),
                  list(H = 30, a = 0.002, h = 0.6, z = 1.8, k = 3, s = 0, T = 24))) {
    set.seed(11)
    ref <- vapply(1:300, function(i) funresMech:::simulate_trial(cs$H, cs$T, cs$a, cs$h, cs$z,
                                                                  cs$k, cs$s, exclude = TRUE), 0L)
    set.seed(11)
    eng <- sim(cs$H, cs$a, cs$h, cs$z, cs$k, cs$s, cs$T, 300L, 0L, 0L, FALSE, FALSE, NULL, NULL,
               FALSE, TRUE)
    expect_identical(as.integer(eng), as.integer(ref))
  }
})

test_that("exclude counts at most one host less than include, trial by trial", {
  # Without the saturation cut-off both conventions use the same random
  # numbers in each trial, so they differ only by the last encounter.
  for (fam in list(c(0L, 0L), c(2L, 1L), c(1L, 3L))) {
    set.seed(3)
    inc <- sim(8L, 0.03, 3, 1.5, 0.8, 1, 24, 2000L, fam[1], fam[2], FALSE, FALSE, NULL, NULL, TRUE, FALSE)
    set.seed(3)
    exc <- sim(8L, 0.03, 3, 1.5, 0.8, 1, 24, 2000L, fam[1], fam[2], FALSE, FALSE, NULL, NULL, TRUE, TRUE)
    d <- inc - exc
    expect_true(all(d %in% c(0L, 1L)), info = paste(fam, collapse = "-"))
    expect_true(any(d == 1L), info = paste(fam, collapse = "-"))
  }
})

test_that("exact case: handling longer than the trial", {
  # With h = 2T and fixed handling, a consumer can attack at most once and
  # never finishes handling: exclude always gives 0 and include gives 1 with
  # probability P(first search time < T) = 1 - exp(-a x T) for exponential search.
  H <- 10L; a <- 0.01; T <- 24
  set.seed(5)
  inc <- sim(H, a, 2 * T, 1, 1, 0, T, 20000L, 0L, 3L, FALSE, TRUE, NULL, NULL, FALSE, FALSE)
  exc <- sim(H, a, 2 * T, 1, 1, 0, T, 20000L, 0L, 3L, FALSE, TRUE, NULL, NULL, FALSE, TRUE)
  expect_true(all(exc == 0L))
  expect_true(all(inc %in% c(0L, 1L)))
  expect_equal(mean(inc), 1 - exp(-a * H * T), tolerance = 0.02)
  for (pred in c(FALSE, TRUE)) {
    exc_p <- sim(H, a, 2 * T, 1, 1, 0, T, 2000L, 0L, 3L, pred, TRUE, NULL, NULL, FALSE, TRUE)
    expect_true(all(exc_p == 0L))
  }
})

test_that("the predator process and the likelihood accept the option", {
  set.seed(7)
  m_inc <- funres_model(process = "predator")
  m_exc <- funres_model(process = "predator", unfinished = "exclude")
  p <- c(a = 0.02, h = 2, z = 1, k = 1, s = 1)
  r_inc <- simulate_response(m_inc, dens = c(4, 16), T = 24, par = p, n_sim = 4000)
  r_exc <- simulate_response(m_exc, dens = c(4, 16), T = 24, par = p, n_sim = 4000)
  expect_true(all(r_exc$response <= r_exc$dens))
  expect_lt(mean(r_exc$response), mean(r_inc$response))
  d <- data.frame(dens = c(4, 4, 16, 16), par = c(2, 3, 5, 7))
  nll <- funresMech:::.nll_sim(p, m_exc, d, 24, 2000)
  expect_true(is.finite(nll) && nll > 0)
})
