# ============================================================================
# 02_engine_tests.R
#
# Tests of the simulation engine of funresMech 1.2.0 (S1 Text, section 2).
#
#   1. Anchor: with seed 88, the NLL of dataset D3 at z = 2.5 and the
#      parameters (a, h, k, s) = (0.012, 0.626, 0.902, 0.542), with 3,000
#      simulated trials per density and no saturation stop, equals the value
#      given by the supplementary code of the method (Okuyama 2026):
#      132.04994217.
#   2. The saturation stop does not change the simulated distribution:
#      chi-squared tests of the simulated frequencies with and without the
#      stop (200,000 simulated trials per case, three parameter sets).
#   3. The saturation stop removes the very long run times that occur when
#      h and k approach zero.
#   4. Moments of the built-in distributions of search and handling times.
#
# Run from the root of the package repository, with the package installed:
#   source("validation/R/02_engine_tests.R")
# Takes a few minutes. Output: validation/tables/engine_tests.csv
# ============================================================================

library(funresMech)
stopifnot(packageVersion("funresMech") >= "1.1.0.9000")
dir.create(file.path("validation", "tables"), showWarnings = FALSE, recursive = TRUE)

d3 <- read.csv(file.path("validation", "data", "D3.csv"))
d2 <- read.csv(file.path("validation", "data", "D2.csv"))
T_exp <- 24
results <- list()
record <- function(id, test, pass, detail) {
  results[[length(results) + 1]] <<- data.frame(id = id, test = test,
    result = if (isTRUE(pass)) "pass" else "fail", detail = detail, stringsAsFactors = FALSE)
  cat(sprintf("[%s] %s: %s\n    %s\n", id, if (isTRUE(pass)) "pass" else "FAIL", test, detail))
}

# ---- 1. Anchor ----------------------------------------------------------------
set.seed(88)
nll <- funresMech:::negloglik_fixed_z(c(0.012, 0.626, 0.902, 0.542), z_fixed = 2.5,
                                      data_spp = d3, T = T_exp, n_sim = 3000,
                                      stop_saturated = FALSE)
record(1, "NLL of D3 equals the value of the supplementary code of the method",
       sprintf("%.8f", nll) == "132.04994217", sprintf("NLL = %.8f (reference 132.04994217)", nll))

# ---- 2. Saturation stop: same distribution --------------------------------------
sat_cases <- list(list(H = 2,  a = 0.012, h = 0.626, z = 2.5, k = 0.902, s = 0.542),
                  list(H = 10, a = 0.3,   h = 0.05,  z = 1.2, k = 0.2,   s = 0.2),
                  list(H = 6,  a = 0.02,  h = 0.3,   z = 3.0, k = 0.05,  s = 0.6))
n_sat <- 200000
p_sat <- vapply(sat_cases, function(p) {
  set.seed(31); no_stop <- funresMech:::sim_engine_cpp(p$H, p$a, p$h, p$z, p$k, p$s, T_exp, n_sat, FALSE)
  set.seed(32); stop    <- funresMech:::sim_engine_cpp(p$H, p$a, p$h, p$z, p$k, p$s, T_exp, n_sat, TRUE)
  tab <- rbind(tabulate(no_stop + 1L, p$H + 1L), tabulate(stop + 1L, p$H + 1L))
  tab <- tab[, colSums(tab) > 0, drop = FALSE]
  if (ncol(tab) < 2) return(1)
  suppressWarnings(stats::chisq.test(tab)$p.value)
}, numeric(1))
record(2, "The saturation stop does not change the simulated distribution", all(p_sat > 0.001),
       sprintf("chi-squared p-values (200,000 trials each): %s", paste(sprintf("%.2f", p_sat), collapse = ", ")))

# ---- 3. Saturation stop: run time ------------------------------------------------
trap <- list(a = 1, h = 1e-6, z = 3, k = 1e-6, s = 0)
set.seed(4)
loops_stop <- funresMech:::count_loops_cpp(50, trap$a, trap$h, trap$z, trap$k, trap$s, T_exp, 10, TRUE)
t_stop <- system.time(funresMech:::negloglik_fixed_z(c(trap$a, trap$h, trap$k, trap$s), trap$z,
                                                     d2, T_exp, 10000))[["elapsed"]]
record(3, "With the saturation stop, h and k close to zero do not cause long run times",
       loops_stop < 1e4 && t_stop < 60,
       sprintf("loops in 10 trials with 50 hosts: %.0f; complete NLL of D2 (10,000 trials): %.2f s",
               loops_stop, t_stop))

# ---- 4. Moments of the built-in distributions ------------------------------------
# Search times: mean 1/lambda (here 2.5) and squared coefficient of variation 1/k.
# Handling times: mean h (here 0.6) and standard deviation s on the natural scale.
draw <- funresMech:::draw_times_cpp
n_m <- 200000
mom <- list()
for (fam in c(gamma = 0L, lognormal = 1L, weibull = 2L)) for (k in c(0.3, 1, 4)) {
  set.seed(1)
  x <- draw(n_m, TRUE, fam, 2.5, k)
  mom[[sprintf("search %s, k = %g", names(which(c(gamma = 0L, lognormal = 1L, weibull = 2L) == fam)), k)]] <-
    c(mean_ratio = mean(x) / 2.5, dispersion_ratio = (stats::var(x) / mean(x)^2) * k)
}
for (fam in c(lognormal = 0L, gamma = 1L, weibull = 2L)) for (s in c(0.1, 0.6, 1.5)) {
  set.seed(2)
  x <- draw(n_m, FALSE, fam, 0.6, s)
  mom[[sprintf("handling %s, s = %g", names(which(c(lognormal = 0L, gamma = 1L, weibull = 2L) == fam)), s)]] <-
    c(mean_ratio = mean(x) / 0.6, dispersion_ratio = stats::sd(x) / s)
}
mom <- do.call(rbind, mom)
record(4, "Built-in distributions have the specified mean and dispersion",
       all(abs(mom[, 1] - 1) < 0.03) && all(abs(mom[, 2] - 1) < 0.08),
       sprintf("ratio to the specified value: mean %.3f to %.3f; dispersion %.3f to %.3f",
               min(mom[, 1]), max(mom[, 1]), min(mom[, 2]), max(mom[, 2])))

# ---- Summary ---------------------------------------------------------------------
out <- do.call(rbind, results)
out$package_version <- as.character(packageVersion("funresMech"))
out$R_version <- R.version.string
write.csv(out, file.path("validation", "tables", "engine_tests.csv"), row.names = FALSE)
print(out[, c("id", "result", "test")], row.names = FALSE)
