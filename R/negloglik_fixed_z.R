############################################################
# Mechanistic model functions
############################################################

# negloglik_fixed_z

negloglik_fixed_z <- function(par_vec, z_fixed, data_spp, T, n_sim = 1500) {
  
  a <- par_vec[1]
  h <- par_vec[2]
  k <- par_vec[3]
  s <- par_vec[4]
  
  if (a <= 0 || h <= 0 || k <= 0 || s < 0) return(1e10)
  
  dens_levels <- sort(unique(data_spp$dens))
  nll <- 0
  
  for (x in dens_levels) {
    subdat <- data_spp[data_spp$dens == x, ]
    y_obs <- subdat$par
    
    probs <- simulate_distribution(x, T, a, h, z_fixed, k, s, n_sim)
    
    for (y in y_obs) {
      p_y <- probs[as.character(y)]
      nll <- nll - log(p_y)
    }
  }
  
  nll
}