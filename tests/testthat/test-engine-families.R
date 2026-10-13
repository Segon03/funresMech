# Generalized engine: distribution families, user samplers and the predator
# process. The default model is covered bit by bit in test-regression.R.

sim <- function(...) funresMech:::sim_trials_cpp(...)
draw <- function(...) funresMech:::draw_times_cpp(...)

fam_search   <- c(gamma = 0L, lognormal = 1L, weibull = 2L)
fam_handling <- c(lognormal = 0L, gamma = 1L, weibull = 2L, constant = 3L)

test_that("general loop reproduces the fast path of the default model exactly", {
  golden <- readRDS(test_path("fixtures", "golden_1.1.0.rds"))
  for (cs in golden$engine) {
    for (cut in c(TRUE, FALSE)) {
      set.seed(cs$seed)
      gen <- sim(cs$H, cs$a, cs$h, cs$z, cs$k, cs$s, cs$T, cs$n, 0L, 0L,
                 FALSE, cut, NULL, NULL, TRUE)
      expect_identical(as.integer(gen), as.integer(if (cut) cs$out_cut else cs$out_nocut))
    }
  }
})

test_that("search-time families have mean 1/lambda and CV^2 = 1/k", {
  n <- 2e5
  for (fam in names(fam_search)) {
    for (k in c(0.3, 1, 4)) {
      set.seed(1)
      x <- draw(n, TRUE, fam_search[[fam]], 2.5, k)
      expect_true(all(x > 0), info = fam)
      expect_equal(mean(x), 2.5, tolerance = 0.03, info = paste(fam, k))
      expect_equal(var(x) / mean(x)^2, 1 / k, tolerance = 0.06, info = paste(fam, k))
    }
  }
})

test_that("handling-time families have mean h and SD s (natural scale)", {
  n <- 2e5
  for (fam in setdiff(names(fam_handling), "constant")) {
    for (s in c(0.1, 0.6, 1.5)) {
      set.seed(2)
      x <- draw(n, FALSE, fam_handling[[fam]], 0.6, s)
      expect_true(all(x > 0), info = fam)
      expect_equal(mean(x), 0.6, tolerance = 0.03, info = paste(fam, s))
      expect_equal(sd(x), s, tolerance = 0.08, info = paste(fam, s))
    }
  }
  expect_true(all(draw(100, FALSE, fam_handling[["constant"]], 0.6, 0.4) == 0.6))
  expect_true(all(draw(100, FALSE, fam_handling[["lognormal"]], 0.6, 0) == 0.6))
})

test_that("Gamma search times follow the Gamma law; Weibull with k = 1 is exponential", {
  set.seed(3)
  x <- draw(5000, TRUE, 0L, 2, 0.5)
  expect_gt(suppressWarnings(ks.test(x, "pgamma", shape = 0.5, rate = 0.5 / 2)$p.value), 1e-3)
  set.seed(4)
  w <- draw(5000, TRUE, 2L, 2, 1)
  expect_gt(suppressWarnings(ks.test(w, "pexp", rate = 1 / 2)$p.value), 1e-3)
})

test_that("exact law: with k = 1, h = 0 and z = 1 both processes give Binomial(H, 1 - exp(-aT))", {
  # Parasitoid: each host receives Poisson(a T) attacks independently.
  # Predator: each prey is consumed at rate a independently (pure death).
  H <- 12L; a <- 0.08; T <- 10; n <- 40000L
  p <- 1 - exp(-a * T)
  expected <- dbinom(0:H, H, p)
  for (pred in c(FALSE, TRUE)) {
    set.seed(5)
    y <- sim(H, a, 0, 1, 1, 0, T, n, 0L, 3L, pred, TRUE)
    obs <- tabulate(y + 1L, nbins = H + 1L)
    keep <- expected * n >= 5
    chi <- sum((obs[keep] - n * expected[keep])^2 / (n * expected[keep]))
    expect_gt(pchisq(chi, df = sum(keep) - 1, lower.tail = FALSE), 1e-3)
  }
})

test_that("predator process: responses within [0, H] and lower than parasitoid when hosts are revisited", {
  set.seed(6)
  yp <- sim(10L, 0.5, 0.05, 1, 1, 0.02, 5, 5000L, 0L, 0L, TRUE, TRUE)
  set.seed(6)
  yq <- sim(10L, 0.5, 0.05, 1, 1, 0.02, 5, 5000L, 0L, 0L, FALSE, TRUE)
  expect_true(all(yp >= 0 & yp <= 10))
  # with removal every encounter is a new victim, so more victims on average
  expect_gt(mean(yp), mean(yq))
})

test_that("user R sampler (mean-1 variates) matches the built-in family in law", {
  r_gamma1 <- function(n, shape) rgamma(n, shape = shape, rate = shape)
  set.seed(7)
  user <- sim(20L, 0.012, 0.626, 2.5, 0.9, 0.542, 24, 20000L, 9L, 0L,
              FALSE, TRUE, r_gamma1, NULL)
  set.seed(8)
  built <- sim(20L, 0.012, 0.626, 2.5, 0.9, 0.542, 24, 20000L, 0L, 0L,
               FALSE, TRUE, NULL, NULL)
  tab <- rbind(tabulate(user + 1L, 21), tabulate(built + 1L, 21))
  tab <- tab[, colSums(tab) > 0, drop = FALSE]
  expect_gt(suppressWarnings(chisq.test(tab)$p.value), 1e-3)

  # reproducible with set.seed() and no repeated random numbers between blocks
  set.seed(9); x1 <- draw(600, TRUE, 9L, 1, 2, r_gamma1)
  set.seed(9); x2 <- draw(600, TRUE, 9L, 1, 2, r_gamma1)
  expect_identical(x1, x2)
  expect_false(isTRUE(all.equal(x1[1:256], x1[257:512])))
})

test_that("user C++ sampler through an external pointer", {
  ptr <- funresMech:::example_cpp_sampler()
  set.seed(10)
  x <- draw(20000, TRUE, 10L, 3, 0.7, ptr)
  expect_equal(mean(x), 3, tolerance = 0.03)
  expect_equal(var(x) / mean(x)^2, 1 / 0.7, tolerance = 0.06)
  set.seed(11)
  y <- sim(15L, 0.05, 0.3, 1.5, 0.7, 0.2, 24, 2000L, 10L, 10L, FALSE, TRUE, ptr, ptr)
  expect_true(all(y >= 0 & y <= 15))
})

test_that("invalid user samplers fail with a clear error", {
  bad_len <- function(n, shape) rep(1, n - 1)
  bad_neg <- function(n, shape) rep(-1, n)
  expect_error(draw(10, TRUE, 9L, 1, 1, bad_len), "exactly n values")
  expect_error(draw(10, TRUE, 9L, 1, 1, bad_neg), "finite and >= 0")
  expect_error(sim(5L, 1, 0.1, 1, 1, 0.1, 1, 10L, 7L, 0L), "unknown search family")
  expect_error(sim(5L, 1, 0.1, 1, 1, 0.1, 1, 10L, 0L, 7L), "unknown handling family")
})
