############################################################
# Profile-likelihood fit of the mechanistic model
############################################################

#' @importFrom DEoptim DEoptim DEoptim.control
#' @importFrom future future value
#' @importFrom stats sd
#' @noRd
NULL

# fit_z_fixed
#
# One point of the profile: with z fixed, estimates (a, h, k, s) with
# DEoptim searching on the log scale (see R/reparam.R), and re-evaluates
# the NLL at the optimum with fresh random numbers. The NLL returned by
# DEoptim (`bestval`) is the MINIMUM of many noisy evaluations and is
# therefore biased downwards; the mean of `n_reeval` fresh evaluations
# (`nll_reeval`) and its standard error are reported alongside it.

fit_z_fixed <- function(z, data_spp, T_exp, n_sim = 3000, itermax = 50,
                        NP = 40, reltol = 1e-2, n_reeval = 3L,
                        bounds = .log_bounds()) {

  H_ref <- .h_ref(data_spp)

  fn <- function(theta) {
    p <- .theta_to_nat(theta, z, H_ref)
    negloglik_fixed_z(p, z_fixed = z, data_spp = data_spp, T = T_exp,
                      n_sim = n_sim)
  }

  res <- DEoptim(
    fn      = fn,
    lower   = bounds$lower,
    upper   = bounds$upper,
    control = DEoptim.control(itermax = itermax, NP = NP, reltol = reltol,
                              trace = FALSE)
  )

  theta   <- res$optim$bestmem
  par_nat <- .theta_to_nat(theta, z, H_ref)

  reevals <- vapply(seq_len(n_reeval), function(i) fn(theta), numeric(1))

  rango   <- bounds$upper - bounds$lower
  at_bound <- (theta <= bounds$lower + 1e-3 * rango) |
              (theta >= bounds$upper - 1e-3 * rango)
  names(at_bound) <- c("lambda_ref", "h", "k", "s")

  list(z = z, par = par_nat, nll = res$optim$bestval,
       nll_reeval = mean(reevals),
       nll_reeval_se = if (n_reeval > 1) sd(reevals) / sqrt(n_reeval) else NA_real_,
       at_bound = at_bound)
}

# TRUE when the current future plan can run this package's code: either the
# plan is sequential (same R session) or the workers have the INSTALLED
# package with the same version and the same code for the model fit.
# Workers are new R processes that load the installed package, so an older
# installed version (or devtools::load_all() with a different build) would
# fail or silently run different code.
.workers_ok <- function() {
  if (inherits(future::plan(), "sequential")) return(TRUE)
  mine <- paste(deparse(body(fit_z_fixed)), collapse = "\n")
  theirs <- tryCatch(
    value(future({
      f <- tryCatch(utils::getFromNamespace("fit_z_fixed", "funresMech"),
                    error = function(e) NULL)
      if (is.null(f)) NA_character_ else paste(deparse(body(f)), collapse = "\n")
    })),
    error = function(e) NA_character_)
  isTRUE(identical(mine, theirs))
}

# Evaluates several values of z, in parallel when the user's future plan
# allows it (the plan is set by the caller; see server.R). The random
# streams are derived from the caller's RNG state, so set.seed() makes
# the whole profile reproducible. If the workers cannot run this version
# of the package, it falls back to a sequential plan with a message.
# On the workers the function is taken from the installed namespace, because
# future does not export non-exported functions of a package.
.eval_profile_points <- function(zs, ...) {
  if (!.workers_ok()) {
    message("The parallel workers do not have this version of funresMech ",
            "installed; evaluating the profile sequentially. Install the ",
            "package to use parallel workers.")
    old_plan <- future::plan(future::sequential)
    on.exit(future::plan(old_plan), add = TRUE)
  }
  fs <- lapply(zs, function(z)
    future(utils::getFromNamespace("fit_z_fixed", "funresMech")(z, ...),
           seed = TRUE))
  lapply(fs, value)
}

.points_to_df <- function(pts) {
  do.call(rbind, lapply(pts, function(p) {
    data.frame(z = p$z, nll = p$nll,
               nll_reeval = p$nll_reeval, nll_reeval_se = p$nll_reeval_se,
               a = p$par[["a"]], h = p$par[["h"]], k = p$par[["k"]],
               s = p$par[["s"]],
               at_bound = any(p$at_bound),
               bound_names = paste(names(p$at_bound)[p$at_bound], collapse = ","),
               stringsAsFactors = FALSE)
  }))
}

#' Check the data of one species
#' @noRd
.check_data_spp <- function(data_spp) {
  if (!all(c("dens", "par") %in% names(data_spp)))
    stop("data_spp needs columns 'dens' and 'par'.", call. = FALSE)
  if (anyNA(data_spp$dens) || anyNA(data_spp$par))
    stop("Missing values in host density or parasitism.", call. = FALSE)
  if (any(data_spp$dens < 1))
    stop("Host densities must be >= 1.", call. = FALSE)
  if (any(data_spp$par < 0 | data_spp$par > data_spp$dens))
    stop("Parasitism must be between 0 and the host density.", call. = FALSE)
  invisible(TRUE)
}


# Text for an upper limit of z that the profile does not cross. It is only
# called "not identified" when the grid already reached the largest value
# tried by the automatic extension (the profile is flat up to there);
# otherwise the limit may simply lie beyond a short grid.
.high_open_text <- function(z_max, z_cap = 20) {
  if (z_max >= z_cap) {
    paste0("upper limit of z not identified (profile flat up to z = ", z_max, ")")
  } else {
    paste0("upper limit of z not reached within the grid (z >= ", z_max,
           "); extend the z grid")
  }
}

# Extra values of z inside the bracket where the profile crosses the 95%
# threshold, when that bracket is wider than the base grid step (e.g. after
# the automatic extension to z = 8, 10, ...). Linear interpolation across a
# wide bracket is crude because the profile is not linear; a few extra points
# make the limit as precise as the rest of the grid.
.refine_points <- function(z, z_low, z_high, step0, max_new = 12) {
  bracket_points <- function(zc) {
    if (is.na(zc) || !any(z < zc) || !any(z >= zc)) return(numeric(0))
    za <- max(z[z < zc]); zb <- min(z[z >= zc])
    if (zb - za <= 1.5 * step0) return(numeric(0))
    cand <- seq(za + step0, zb - 0.5 * step0, by = step0)
    cand[cand > za + 0.5 * step0]
  }
  new <- sort(unique(round(c(bracket_points(z_low), bracket_points(z_high)), 10)))
  new <- new[!new %in% round(z, 10)]
  utils::head(new, max_new)
}

# fit_profile
#
# Full analysis of one species:
#  * profile of the NLL over a grid of z (always including z = 1 when the
#    grid contains it), each point fitted with fit_z_fixed();
#  * z_hat = minimum of the profile (Okuyama 2026, Sec. 2.3), with the
#    parameters of the fit at that z;
#  * 95% likelihood-ratio interval (threshold +1.92), reported as
#    censored when the profile does not cross it (z is then not
#    identified on that side). If the upper limit is open and
#    `extend_z` is TRUE, the grid is extended to the values in `z_extend`;
#    when a limit falls in a bracket wider than the grid step, extra points
#    are evaluated inside it (`refine`);
#  * AIC of the model with z free (5 parameters) against z = 1 (4);
#  * diagnostic notes (small k, parameters on a bound, noisy NLL, open
#    interval).

fit_profile <- function(data_spp, T_exp, z_grid, n_sim = 3000, itermax = 50,
                        NP = 40, reltol = 1e-2, n_reeval = 3L,
                        extend_z = TRUE, z_extend = c(8, 10, 15, 20),
                        refine = TRUE) {

  .check_data_spp(data_spp)

  z_grid <- sort(unique(round(z_grid, 10)))
  step0  <- if (length(z_grid) > 1) stats::median(diff(z_grid)) else 1
  if (min(z_grid) <= 1 && max(z_grid) >= 1) z_grid <- sort(unique(c(z_grid, 1)))

  args <- list(data_spp = data_spp, T_exp = T_exp, n_sim = n_sim,
               itermax = itermax, NP = NP, reltol = reltol,
               n_reeval = n_reeval)

  pts <- do.call(.eval_profile_points, c(list(z_grid), args))
  prof <- .points_to_df(pts)
  ci <- .profile_ci(prof$z, prof$nll)

  extended <- FALSE
  if (isTRUE(extend_z) && ci$high_censored) {
    z_new <- z_extend[z_extend > max(prof$z)]
    if (length(z_new) > 0) {
      pts_new <- do.call(.eval_profile_points, c(list(z_new), args))
      prof <- rbind(prof, .points_to_df(pts_new))
      prof <- prof[order(prof$z), ]
      rownames(prof) <- NULL
      ci <- .profile_ci(prof$z, prof$nll)
      extended <- TRUE
    }
  }

  refined <- FALSE
  if (isTRUE(refine)) {
    z_new <- .refine_points(prof$z, ci$z_low, ci$z_high, step0)
    if (length(z_new) > 0) {
      pts_new <- do.call(.eval_profile_points, c(list(z_new), args))
      prof <- rbind(prof, .points_to_df(pts_new))
      prof <- prof[order(prof$z), ]
      rownames(prof) <- NULL
      ci <- .profile_ci(prof$z, prof$nll)
      refined <- TRUE
    }
  }

  i_best <- which.min(prof$nll)
  par <- c(a = prof$a[i_best], h = prof$h[i_best], z = prof$z[i_best],
           k = prof$k[i_best], s = prof$s[i_best])

  nll_1 <- .nll_at_z(prof$z, prof$nll, 1)
  aic_full <- 2 * prof$nll[i_best] + 2 * 5
  aic_restricted <- if (is.na(nll_1)) NA_real_ else 2 * nll_1 + 2 * 4

  notes <- character(0)
  if (prof$k[i_best] < 0.1)
    notes <- c(notes, paste0(
      "k is very small (", signif(prof$k[i_best], 2), "): search times are ",
      "extremely variable and z is weakly identified from above (ridge ",
      "between z and k)."))
  if (ci$high_censored)
    notes <- c(notes, paste0(
      "The profile does not cross the 95% threshold above z_hat: ",
      .high_open_text(max(prof$z), max(z_extend)), "."))
  if (ci$low_censored)
    notes <- c(notes, paste0(
      "The profile does not cross the 95% threshold below z_hat: the lower ",
      "limit of z is <= ", min(prof$z), "."))
  if (prof$at_bound[i_best])
    notes <- c(notes, paste0(
      "At z_hat the optimiser ended on a bound (",
      prof$bound_names[i_best], "); the estimate of that parameter is only ",
      "a limit."))
  gap <- prof$nll_reeval[i_best] - prof$nll[i_best]
  if (is.finite(gap) && gap > 2)
    notes <- c(notes, paste0(
      "The NLL re-evaluated at the optimum exceeds the optimiser's by ",
      round(gap, 2), " units: Monte Carlo noise is large; increase n_sim."))

  list(profile = prof, par = par, nll = prof$nll[i_best],
       aic = data.frame(AIC_full = aic_full, AIC_restricted = aic_restricted,
                        delta_AIC = aic_restricted - aic_full),
       ci = ci, extended = extended, refined = refined, notes = notes)
}
