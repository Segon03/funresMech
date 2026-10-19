# ---------------------------------------------------------------------------
# Generates golden_1.1.0.rds: reference outputs of the simulation engine and
# of the simulated likelihood of funresMech 1.1.0 (as released on CRAN,
# GitHub commit 52a779d). It is NOT run by testthat; it was run once with
# funresMech 1.1.0 installed, and test-regression.R checks that later
# versions reproduce these numbers EXACTLY for the default model
# (Gamma search times, Lognormal handling times).
#
# Run from this folder:  Rscript make_golden.R
# ---------------------------------------------------------------------------
stopifnot(as.character(utils::packageVersion("funresMech")) == "1.1.0")
# The internal name of the engine function changed after version 1.1.0;
# look it up so that this script still runs with 1.1.0 installed.
ns <- asNamespace("funresMech")
engine_fn <- get(grep("^(sim_engine|motor_[a-z]+)_cpp$", ls(ns, all.names = TRUE), value = TRUE)[1], ns)
nllfz <- funresMech:::negloglik_fixed_z

d2 <- read.csv("D2.csv"); d3 <- read.csv("D3.csv")
T_exp <- 24

engine_cases <- list(
  list(H = 30L, a = 0.5,   h = 0.2,   z = 1.5, k = 0.5,  s = 0,    T = 2),
  list(H = 8L,  a = 0.01,  h = 0.6,   z = 3,   k = 0.05, s = 0.8,  T = 24),
  list(H = 50L, a = 0.001, h = 0.5,   z = 2.5, k = 0.04, s = 0.6,  T = 24),
  list(H = 6L,  a = 2,     h = 0.05,  z = 1,   k = 1,    s = 0.05, T = 5),
  list(H = 20L, a = 0.012, h = 0.626, z = 2.5, k = 0.902, s = 0.542, T = 24)
)
engine <- lapply(seq_along(engine_cases), function(i) {
  cs <- engine_cases[[i]]
  out <- lapply(c(TRUE, FALSE), function(cut) {
    set.seed(1000 + i)
    engine_fn(cs$H, cs$a, cs$h, cs$z, cs$k, cs$s, cs$T, 500L, cut)
  })
  c(cs, list(seed = 1000 + i, n = 500L, out_cut = out[[1]], out_nocut = out[[2]]))
})

nll_cases <- list(
  list(data = "D3", par = c(0.012, 0.626, 0.902, 0.542), z = 2.5, n_sim = 3000L, cut = FALSE, seed = 88),
  list(data = "D3", par = c(0.012, 0.626, 0.902, 0.542), z = 2.5, n_sim = 3000L, cut = TRUE,  seed = 88),
  list(data = "D2", par = c(0.001, 0.606, 0.037, 0.682), z = 3.1, n_sim = 2000L, cut = TRUE,  seed = 21),
  list(data = "D2", par = c(0.02,  0.5,   1.2,   0.3),   z = 1.5, n_sim = 2000L, cut = TRUE,  seed = 22),
  list(data = "D3", par = c(0.3,   0.05,  1,     0),     z = 1,   n_sim = 2000L, cut = TRUE,  seed = 23)
)
nll <- lapply(nll_cases, function(cs) {
  d <- if (cs$data == "D2") d2 else d3
  set.seed(cs$seed)
  cs$nll <- nllfz(cs$par, z_fixed = cs$z, data_spp = d, T = T_exp,
                  n_sim = cs$n_sim, stop_saturated = cs$cut)
  cs
})

golden <- list(version = "1.1.0", T_exp = T_exp, engine = engine, nll = nll,
               R = R.version.string, created = format(Sys.time(), "%Y-%m-%d %H:%M"))
saveRDS(golden, "golden_1.1.0.rds", version = 2)
print(vapply(nll, function(x) sprintf("%s z=%.1f cut=%s: %.8f", x$data, x$z, x$cut, x$nll), ""))
