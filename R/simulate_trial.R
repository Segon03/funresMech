############################################################
# Mechanistic model functions
############################################################

#simulate_trial

simulate_trial <- function(x, T, a, h, z, k, s) {
  t <- 0
  Nx <- rep(0, x)
  lambda <- a * (x^z)
  
  while (t < T) {
    ts <- rgamma(1, shape = k, rate = lambda * k)
    t <- t + ts
    if (t >= T) break
    
    host_id <- sample.int(x, 1)
    Nx[host_id] <- Nx[host_id] + 1
    
    if (s > 0) {
      sdlog <- s
      meanlog <- log(h) - 0.5 * sdlog^2
      th <- rlnorm(1, meanlog = meanlog, sdlog = sdlog)
    } else {
      th <- h
    }
    
    t <- t + th
  }
  
  sum(Nx >= 1)
}