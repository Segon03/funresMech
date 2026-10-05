############################################################
# Mechanistic model functions
############################################################

# simulate_distribution
#
# Probability of 0, 1, ..., x parasitised hosts at density x, estimated
# from n_sim simulated trials with the C++ engine. Probabilities that are
# zero in the simulation are replaced by .Machine$double.xmin (as in
# Okuyama 2026) so that the log-likelihood stays defined.

simulate_distribution <- function(x, T, a, h, z, k, s, n_sim = 1500,
                                  corte_saturacion = TRUE) {
  counts <- motor_okuyama_cpp(as.integer(x), a, h, z, k, s, T,
                              as.integer(n_sim), corte_saturacion)
  probs <- tabulate(counts + 1L, nbins = x + 1L) / n_sim
  names(probs) <- 0:x
  probs[probs == 0] <- .Machine$double.xmin
  probs
}
