############################################################
# Mechanistic model functions
############################################################

# fit_full

fit_full <- function(data_spp, T_exp, itermax, NP, reltol, n_sim_profile) {
  
  lower <- c(a = 0.001, h = 0.001, z = 0.5, k = 0.5, s = 0.001)
  upper <- c(a = 2.0,   h = 0.5,   z = 3.0, k = 5.0, s = 0.5)
  
  set.seed(123)
  
  res <- DEoptim(
    fn = function(par) {
      a <- par[1]; h <- par[2]; z <- par[3]; k <- par[4]; s <- par[5]
      negloglik_fixed_z(c(a,h,k,s), z_fixed = z,
                        data_spp = data_spp, T = T_exp,
                        n_sim = n_sim_profile)
    },
    lower = lower,
    upper = upper,
    control = DEoptim.control(itermax = itermax,
                              NP = NP,
                              reltol = reltol,
                              trace = TRUE)
  )
  
  best <- res$optim$bestmem
  names(best) <- c("a","h","z","k","s")
  nll <- res$optim$bestval
  
  list(par = best, nll = nll)
}