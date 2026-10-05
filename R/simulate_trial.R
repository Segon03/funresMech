############################################################
# Mechanistic model functions
############################################################

#' @importFrom stats rgamma runif rlnorm
#' @noRd
NULL

# simulate_trial
#
# Pure-R REFERENCE implementation of one foraging trial (Okuyama 2012,
# 2026). The fitting code does not use it: it runs on the C++ engine
# `motor_okuyama_cpp()` (src/motor.cpp). It is kept so that the tests can
# check, with the same seed, that the engine reproduces it exactly.
#
# s is the standard deviation of the handling time on the NATURAL scale
# (same units as h), as in Okuyama (2026); s = 0 gives a fixed handling
# time. The Lognormal has mean h and sd s, i.e.
# sdlog = sqrt(log(1 + (s / h)^2)).

simulate_trial <- function(x, T, a, h, z, k, s) {
  if (x < 1) return(0L)

  tocado <- logical(x)
  rate <- k * a * (x^z)

  if (s > 0) {
    sdlog   <- sqrt(log(1 + (s / h)^2))
    meanlog <- log(h) - 0.5 * sdlog^2
  }

  t <- 0
  while (t < T) {
    t <- t + rgamma(1, shape = k, scale = 1 / rate)
    if (t > T) break

    host_id <- floor(runif(1, 0, x)) + 1
    tocado[host_id] <- TRUE

    th <- if (s > 0) rlnorm(1, meanlog = meanlog, sdlog = sdlog) else h
    t <- t + th
  }

  sum(tocado)
}
