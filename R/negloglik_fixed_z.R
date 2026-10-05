############################################################
# Mechanistic model functions
############################################################

# negloglik_fixed_z
#
# Simulated negative log-likelihood (Okuyama 2026) for a fixed z.
# par_vec = c(a, h, k, s) on the natural scale (s = sd of the handling
# time, same units as h). Returns 1e10 for parameters outside the valid
# domain so that DEoptim discards them without failing.

negloglik_fixed_z <- function(par_vec, z_fixed, data_spp, T, n_sim = 1500,
                              corte_saturacion = TRUE) {

  a <- par_vec[[1]]
  h <- par_vec[[2]]
  k <- par_vec[[3]]
  s <- par_vec[[4]]

  if (!is.finite(a) || !is.finite(h) || !is.finite(k) || !is.finite(s) ||
      a <= 0 || h <= 0 || k <= 0 || s < 0) return(1e10)

  minp <- .Machine$double.xmin
  nll  <- 0

  # Densities are visited in order of appearance, as in Okuyama (2026);
  # only the order in which the random stream is consumed depends on it.
  for (x in unique(data_spp$dens)) {
    if (x < 1) next
    y_obs <- data_spp$par[data_spp$dens == x]
    if (anyNA(y_obs) || any(y_obs < 0 | y_obs > x)) return(1e10)

    sim  <- motor_okuyama_cpp(as.integer(x), a, h, z_fixed, k, s, T,
                              as.integer(n_sim), corte_saturacion)
    frec <- tabulate(sim + 1L, nbins = x + 1L) / n_sim
    nll  <- nll - sum(log(pmax(frec[y_obs + 1L], minp)))
  }

  nll
}
