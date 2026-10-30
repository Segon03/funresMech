# ============================================================================
# 22_conventional_methods.R
#
# Conventional analyses of the data simulated from the mechanistic process
# (parasitoid process), to compare them with the simulation-based method on
# the same datasets. For each of the 432 simulated datasets of the simulation
# study (core, dense and extended-range designs; validation/data/
# simulated_datasets_parasitoid.csv), nine conventional fits are made:
#
#   mean function  x  error model
#   holling        x  binomial ML, beta-binomial ML, least squares
#   rogers         x  binomial ML, beta-binomial ML, least squares
#   eq2            x  binomial ML, beta-binomial ML, least squares
#
#   holling Generalized Holling disc equation, as in frair's flexp and
#           hollingsII (prey
#           replaced, so the density stays constant):
#             Ne = a N0^z T / (1 + a h N0^z), capped at N0.
#   rogers  Generalized Rogers random predator equation, as in frair's
#           flexpnr and rogersII (prey not replaced; attack rate b N0^q with z = q + 1),
#           the usual choice in functional-response studies:
#             Ne = N0 {1 - exp[a N0^(z-1) (h Ne - T)]}, solved with the
#             Lambert W function.
#   eq2     The mean used by Okuyama (2026) for parasitoids, whose hosts are
#           not removed:  mu = x {1 - exp[-f(x) T / x]},
#             f(x) = a x^z / (1 + a h x^z).
#
# Every fit is by differential evolution followed by a Nelder-Mead polish,
# with time in units of T and the same bounds as the package (a* in
# [1e-8, 1e4], h* = h / T in [1e-6, 1], z in [0.1, 20], beta-binomial rho in
# (0, 1)). The 95% interval of z is the profile interval: likelihood ratio
# (chi-square, 1 df) for the binomial and beta-binomial models, and the
# extra-sum-of-squares F test for least squares. z is tested against 1 with
# the same statistics.
#
# The results are merged with the simulation-based fits of the same datasets
# (validation/tables/simulation_study_fits.csv and extended_range_fits.csv).
#
# The first two are the flexible models of frair with and without
# replacement, the usual choices in functional-response studies; frair has
# no model for a parasitoid that leaves hosts in the arena and is counted by
# distinct hosts, which is what eq2 (Rogers' random parasite equation for
# z = 1) represents.
#
# Usage, from the root of the package repository:
#   QUICK <- TRUE;  source("validation/R/22_conventional_methods.R")   # 6 datasets
#   rm(list = ls()); source("validation/R/22_conventional_methods.R")    # all: about 2 h on one core, minutes on 20 or more
# Results: validation/tables/conventional_fits.csv and conventional_summary.csv.
# ============================================================================

library(funresMech)
source(file.path("validation", "R", "helpers.R"))
if (!exists("QUICK")) QUICK <- FALSE
WORKERS <- if (QUICK) 2L else max(1L, parallel::detectCores() - 1L)
DIR_TAB <- file.path("validation", "tables")
SUFFIX <- if (QUICK) "_quick" else ""

# Datasets: the parasitoid simulation study (core, dense and extended range) and,
# when they have been run, the studies of 23_new_scenarios.R (parasitoid with
# z = 0.5, predator process).
STUDIES <- c("parasitoid_z05", "predator")
D <- utils::read.csv(file.path("validation", "data", "simulated_datasets_parasitoid.csv"), stringsAsFactors = FALSE)
for (st in STUDIES) {
  f <- file.path("validation", "data", paste0("simulated_datasets_", st, if (QUICK) "_quick", ".csv"))
  if (file.exists(f)) D <- rbind(D, utils::read.csv(f, stringsAsFactors = FALSE))
}
keys <- unique(D[, c("source", "scenario", "dataset")])
if (QUICK) keys <- do.call(rbind, lapply(split(keys, keys$source), function(k) k[unique(c(1, nrow(k))), ]))
rownames(keys) <- NULL
cat("Datasets by source:\n"); print(table(keys$source))

# ---- mean functions (time in units of T) ------------------------------------------
# Principal branch of the Lambert W function from the log of its argument
# (W(x) + log W(x) = log x), vectorized and safe for very large arguments.
lambertW_log <- function(lx) {
  w <- exp(pmin(lx, 0))                       # W(x) ~ x for small x (exact enough below exp(-30))
  big <- is.finite(lx) & lx > -30
  if (any(big)) {
    l <- lx[big]
    v <- ifelse(l < 1, exp(l), l - log(pmax(l, 1)))
    for (i in 1:100) {
      v_new <- v - (v + log(v) - l) / (1 + 1 / v)
      v_new <- ifelse(!is.finite(v_new) | v_new <= 0, v / 10, v_new)
      done <- max(abs(v_new - v) / v) < 1e-13
      v <- v_new
      if (done) break
    }
    w[big] <- v
  }
  w
}
mean_rogers <- function(N, a, h, z) {
  ae <- a * N^(z - 1)                                   # attack rate b N0^q
  lx <- log(ae * h * N) - ae * (1 - h * N)
  N - lambertW_log(lx) / (ae * h)
}
mean_eq2 <- function(N, a, h, z) funresMech:::.mu_parasitoid(N, a, h, z, 1)
mean_holling <- function(N, a, h, z) a * N^z / (1 + a * h * N^z)
MEANS <- list(holling = mean_holling, rogers = mean_rogers, eq2 = mean_eq2)

LO <- c(log(1e-8), log(1e-6), log(0.1), -14)
UP <- c(log(1e4), log(1), log(20), 7)

# Objective for theta = (log a*, log h*, [log z], [logit rho]); z fixed if z_fix given
objective <- function(theta, d, mean_fun, error, z_fix = NULL) {
  a <- exp(theta[1]); h <- exp(theta[2])
  z <- if (is.null(z_fix)) exp(theta[3]) else z_fix
  mu <- mean_fun(d$dens, a, h, z)
  if (any(!is.finite(mu))) return(1e10)
  mu <- pmin(pmax(mu, 0), d$dens)
  if (error == "ls") return(sum((d$par - mu)^2))
  p <- pmin(pmax(mu / d$dens, 1e-12), 1 - 1e-12)
  if (error == "binomial") return(-sum(stats::dbinom(d$par, d$dens, p, log = TRUE)))
  rho <- stats::plogis(theta[length(theta)])
  sh <- funresMech:::.bb_shapes(p, rho)
  ll <- funresMech:::.dbbinom_log(d$par, d$dens, sh$alpha, sh$beta)
  if (any(!is.finite(ll))) 1e10 else -sum(ll)
}

par_index <- function(error, z_free) c(1, 2, if (z_free) 3, if (error == "bb") 4)

fit_one <- function(d, mean_fun, error, z_fix = NULL) {
  idx <- par_index(error, is.null(z_fix))
  fn <- function(th) objective(th, d, mean_fun, error, z_fix)
  de <- suppressWarnings(DEoptim::DEoptim(fn, LO[idx], UP[idx],
                                          DEoptim::DEoptim.control(itermax = 150, NP = 40, trace = FALSE)))
  best <- de$optim$bestmem; val <- de$optim$bestval
  pol <- tryCatch(stats::optim(best, fn, method = "Nelder-Mead",
                               control = list(maxit = 5000, reltol = 1e-12)), error = function(e) NULL)
  if (!is.null(pol) && pol$value < val) { best <- pol$par; val <- pol$value }
  list(theta = unname(best), value = val, idx = idx)
}

# Value of the objective with z fixed, starting from a previous optimum
profile_value <- function(z, d, mean_fun, error, start) {
  fn <- function(th) objective(th, d, mean_fun, error, z)
  r <- stats::optim(start, fn, method = "Nelder-Mead", control = list(maxit = 3000, reltol = 1e-10))
  list(value = r$value, theta = r$par)
}

analyse <- function(d, mean_name, error, z_true) {
  mean_fun <- MEANS[[mean_name]]
  n <- nrow(d)
  f <- fit_one(d, mean_fun, error)
  th <- f$theta
  z_hat <- exp(th[3]); best <- f$value
  n_par <- length(f$idx)
  # statistic of a value of the objective relative to the best one
  stat <- function(v) if (error == "ls") (v - best) / (best / (n - n_par)) else 2 * (v - best)
  thr <- if (error == "ls") stats::qf(0.95, 1, n - n_par) else stats::qchisq(0.95, 1)
  # fit with z = 1 (test of z against 1)
  f1 <- fit_one(d, mean_fun, error, z_fix = 1)
  if (f1$value < best) f1$value <- best                     # z = 1 cannot beat the free fit
  p_z1 <- if (error == "ls") stats::pf(stat(f1$value), 1, n - n_par, lower.tail = FALSE)
          else stats::pchisq(stat(f1$value), 1, lower.tail = FALSE)
  # profile of z on a grid around the estimate, then the crossings by root finding
  start_other <- th[-3]
  prof_at <- function(z, start) profile_value(z, d, mean_fun, error, start)
  side <- function(direction) {
    zs <- exp(log(z_hat) + direction * seq(0.05, 6, by = 0.15))
    zs <- zs[zs >= 0.1 & zs <= 20]
    st <- start_other; prev_z <- z_hat; prev_s <- 0
    for (z in zs) {
      r <- prof_at(z, st); st <- r$theta
      if (r$value < best) { best <<- r$value }
      s <- stat(r$value)
      if (s >= thr) {
        g <- function(lz) stat(prof_at(exp(lz), st)$value) - thr
        root <- tryCatch(stats::uniroot(g, sort(c(log(prev_z), log(z))), tol = 1e-4)$root,
                         error = function(e) log(z))
        return(c(limit = exp(root), open = 0))
      }
      prev_z <- z; prev_s <- s
    }
    c(limit = NA_real_, open = 1)
  }
  lo <- side(-1); hi <- side(+1)
  a_hat <- exp(th[1]); h_hat <- exp(th[2])
  data.frame(mean = mean_name, error = error,
             z_hat = z_hat, a_hat = a_hat / T_EXP, h_hat = h_hat * T_EXP,
             rho_hat = if (error == "bb") stats::plogis(th[4]) else NA_real_,
             objective = f$value, n_par = n_par,
             nll = if (error == "ls") NA_real_ else f$value,
             z_low = unname(lo["limit"]), z_high = unname(hi["limit"]),
             low_open = unname(lo["open"]) == 1, high_open = unname(hi["open"]) == 1,
             p_z1 = p_z1,
             covers = (is.na(lo["limit"]) || lo["limit"] <= z_true) && (is.na(hi["limit"]) || z_true <= hi["limit"]),
             excludes_1 = (!is.na(lo["limit"]) && lo["limit"] > 1) || (!is.na(hi["limit"]) && hi["limit"] < 1),
             stringsAsFactors = FALSE)
}

scenario_of <- function(id) {
  p <- regmatches(id, regexec("z([0-9.]+)_k([0-9.]+)_sh([0-9.]+)_hT([0-9.]+)_E([0-9.]+)_([a-z]+)$", id))[[1]]
  e <- scenario(z = as.numeric(p[2]), k = as.numeric(p[3]), s_over_h = as.numeric(p[4]),
                h_over_T = as.numeric(p[5]), E = as.numeric(p[6]), design = p[7],
                generating_model = sub("_z[0-9.]+_.*$", "", id))
  e$id <- id
  e
}

run_key <- function(i, means = names(MEANS)) {
  library(funresMech)
  key <- keys[i, ]
  d <- D[D$scenario == key$scenario & D$dataset == key$dataset & D$source == key$source, c("dens", "par")]
  e <- scenario_of(key$scenario)
  p0 <- true_parameters(e)
  set.seed(1000 + i)
  out <- do.call(rbind, lapply(means, function(m) do.call(rbind, lapply(c("binomial", "bb", "ls"),
    function(err) analyse(d, m, err, e$z)))))
  cbind(key[rep(1, nrow(out)), ], out,
        process = if (e$generating_model == "predator") "predator" else "parasitoid", z = e$z, a_true = p0[["a"]], h_true = p0[["h"]], k = e$k,
        s_over_h = e$s_over_h, h_over_T = e$h_over_T, E = e$E, design = e$design, n_trials = nrow(d),
        row.names = NULL)
}

# One file per dataset, so an interrupted run resumes where it stopped
DIR_RES <- file.path("validation", "results", if (QUICK) "conventional_quick" else "conventional")
dir.create(DIR_RES, recursive = TRUE, showWarnings = FALSE)
res_file <- function(i) file.path(DIR_RES, sprintf("%s__%s__d%d.rds", keys$source[i], keys$scenario[i], keys$dataset[i]))
# datasets with no file, or with a file that lacks some mean function
missing_means <- function(i) {
  f <- res_file(i)
  if (!file.exists(f)) return(names(MEANS))
  setdiff(names(MEANS), unique(readRDS(f)$mean))
}
todo <- lapply(seq_len(nrow(keys)), missing_means)
pending <- which(lengths(todo) > 0)
cat(sprintf("Datasets: %d in total, %d pending; %d workers.\n", nrow(keys), length(pending), WORKERS))
if (length(pending)) {
  cl <- parallel::makeCluster(WORKERS)
  parallel::clusterExport(cl, c("keys", "D", "MEANS", "LO", "UP", "objective", "par_index", "fit_one",
                                "profile_value", "analyse", "scenario_of", "lambertW_log", "mean_rogers",
                                "mean_eq2", "mean_holling", "scenario", "BASE_SCENARIO", "true_parameters",
                                "DESIGNS", "T_EXP", "run_key", "res_file", "DIR_RES", "todo"))
  t0 <- Sys.time()
  invisible(parallel::parLapplyLB(cl, pending, function(i) {
    new <- run_key(i, todo[[i]])
    old <- if (file.exists(res_file(i))) readRDS(res_file(i))
    saveRDS(rbind(old, new), res_file(i)); TRUE }))
  parallel::stopCluster(cl)
  cat(sprintf("Conventional fits of %d datasets in %.1f min\n", length(pending),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
R <- do.call(rbind, lapply(vapply(seq_len(nrow(keys)), res_file, ""), function(f) if (file.exists(f)) readRDS(f)))

# ---- merge with the simulation-based fits of the same datasets ----------------------
M1 <- utils::read.csv(file.path(DIR_TAB, "simulation_study_fits.csv"), stringsAsFactors = FALSE)
M1 <- M1[M1$arm %in% c("core_exploratory", "dense_design"), ]
M2 <- utils::read.csv(file.path(DIR_TAB, "extended_range_fits.csv"), stringsAsFactors = FALSE)
cols <- c("scenario", "dataset", "z_hat", "a_hat", "h_hat", "k_hat", "z_low", "z_high", "low_open",
          "high_open", "covers", "excludes_1", "nll_precise_hat")
M1$source <- ifelse(grepl("_dense$", M1$scenario), "dense", "core")
M2$source <- "extended_range"
M <- rbind(M1[, c("source", cols)], M2[, c("source", cols)])
for (st in STUDIES) {
  f <- file.path(DIR_TAB, paste0(st, "_fits", if (QUICK) "_quick", ".csv"))
  if (file.exists(f)) {
    M3 <- utils::read.csv(f, stringsAsFactors = FALSE); M3$source <- st
    M <- rbind(M, M3[, c("source", cols)])
  }
}
M <- merge(M, unique(R[, c("source", "scenario", "dataset", "process", "z", "a_true", "h_true", "k", "s_over_h",
                           "h_over_T", "E", "design", "n_trials")]))
M$mean <- "mechanistic"; M$error <- "simulation"; M$rho_hat <- NA_real_
M$objective <- M$nll_precise_hat; M$nll <- M$nll_precise_hat; M$n_par <- 5L; M$p_z1 <- NA_real_
M$nll_precise_hat <- NULL
A <- rbind(R[, names(R)], M[, names(R)])
A$method <- ifelse(A$mean == "mechanistic", "mechanistic", paste(A$mean, A$error, sep = "_"))
A$z_error <- A$z_hat - A$z
A$h_ratio <- A$h_hat / A$h_true
A$a_ratio <- A$a_hat / A$a_true
A$open <- A$low_open | A$high_open
A$AIC <- ifelse(is.na(A$nll), NA_real_, 2 * A$nll + 2 * A$n_par)
# conclusion about the type: z > 1 (type III), z < 1, or z = 1 not rejected (type II)
A$type_hat <- ifelse(A$excludes_1, ifelse(A$z_hat > 1, "above 1", "below 1"), "1")
A$type_true <- ifelse(A$z > 1, "above 1", ifelse(A$z < 1, "below 1", "1"))
A$correct_type <- A$type_hat == A$type_true
utils::write.csv(A, file.path(DIR_TAB, paste0("conventional_fits", SUFFIX, ".csv")), row.names = FALSE)

# ---- summaries ----------------------------------------------------------------------
summ <- function(g) {
  w <- wilson(g$covers)
  data.frame(n = nrow(g), coverage = w[["p"]], coverage_low = w[["lower"]], coverage_high = w[["upper"]],
             type_I_error = if (any(g$z == 1)) mean(g$excludes_1[g$z == 1]) else NA_real_,
             power = if (any(g$z > 1)) mean(g$excludes_1[g$z > 1]) else NA_real_,
             detects_z_below_1 = if (any(g$z < 1)) mean(g$excludes_1[g$z < 1]) else NA_real_,
             correct_type = mean(g$correct_type),
             median_z_error = stats::median(g$z_error), rmse_log_z = sqrt(mean(log(g$z_hat / g$z)^2)),
             median_h_ratio = stats::median(g$h_ratio), open_interval = mean(g$open))
}
blk <- ifelse(A$s_over_h != 0.5, "s/h = 1.5", ifelse(A$h_over_T != 0.025, "h/T = 0.1",
              ifelse(A$E != 20, "E = 5", paste("base,", A$design))))
groups <- list(all = A$process, k = paste0(A$process, ", k = ", A$k), z = paste0(A$process, ", z = ", A$z),
               block = paste0(A$process, ", ", blk))
S <- do.call(rbind, lapply(names(groups), function(gn) {
  sp <- split(A, list(A$method, groups[[gn]]), drop = TRUE)
  do.call(rbind, lapply(sp, function(g) cbind(grouping = gn, level = groups[[gn]][match(rownames(g)[1], rownames(A))],
                                               method = g$method[1], summ(g))))
}))
rownames(S) <- NULL
utils::write.csv(S, file.path(DIR_TAB, paste0("conventional_summary", SUFFIX, ".csv")), row.names = FALSE)
print(S[S$grouping == "all", c("level", "method", "n", "coverage", "type_I_error", "power", "median_z_error",
                                "rmse_log_z", "median_h_ratio", "open_interval")], digits = 3)
