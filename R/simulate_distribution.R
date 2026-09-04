############################################################
# Mechanistic model functions
############################################################

#simulate_distribution

simulate_distribution <- function(x, T, a, h, z, k, s, n_sim = 1500) {
  counts <- integer(n_sim)
  for (i in seq_len(n_sim)) {
    counts[i] <- simulate_trial(x, T, a, h, z, k, s)
  }
  
  tab <- table(counts)
  probs <- rep(0, x + 1)
  names(probs) <- 0:x
  probs[names(tab)] <- as.numeric(tab) / n_sim
  
  eps <- .Machine$double.xmin
  probs[probs == 0] <- eps
  
  probs
}