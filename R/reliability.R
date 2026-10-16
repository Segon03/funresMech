############################################################
# Reliability checks of a fit
############################################################

# Thresholds of the checks (NLL units unless stated). They are provisional
# and will be calibrated with the simulation study (validation/, phase 5).
.rel_thresholds <- function() {
  list(se_warn = 0.25,          # Monte Carlo SE of the NLL at the optimum
       gap_warn = 2,            # re-evaluated minus optimiser NLL
       late_warn = 1,           # improvement in the last quarter of the run
       rerun_ok = 0.5,          # |NLL difference| between independent runs
       rerun_fail = 1.92,       #   (1.92 = half-width of the 95% threshold)
       better_fail = 1,         # NLL improvement with more effort / wider bounds
       par_warn = 0.25)         # relative change of a or h between runs
}

.rel_row <- function(check, status, value, detail, action) {
  data.frame(check = check, status = status, value = value, detail = detail,
             action = action, stringsAsFactors = FALSE)
}

#' Reliability checks of a fit
#'
#' Summarises, as a table with a status for each check (`"ok"`, `"warning"`
#' or `"fail"`) and the recommended action, the diagnostics available in a
#' fitted model. [verify_fit()] adds checks that require extra optimisations
#' (independent runs, more effort, wider bounds).
#'
#' A simulated likelihood is noisy and the optimiser can stop at a suboptimal
#' point. A fit should be reported only when no check fails, and warnings
#' should be understood and, when possible, removed by acting on them.
#'
#' @param fit A `funres_fit` from [fit_funres()].
#' @return Data frame with columns `check`, `status`, `value`, `detail` and
#'   `action`; attribute `overall` holds the worst status.
#' @seealso [verify_fit()]
#' @export
reliability <- function(fit) {
  if (!inherits(fit, "funres_fit")) stop("Need a funres_fit.", call. = FALSE)
  th <- .rel_thresholds()
  prof <- fit$profile
  i <- which.min(prof$nll)
  rows <- list()

  nf <- prof$n_floor[i]
  nf_other <- if (nrow(prof) > 1) sum(prof$n_floor[-i] > 0, na.rm = TRUE) else 0L
  rows$floor <- .rel_row(
    "Observed outcomes produced by the simulation",
    if (is.finite(nf) && nf > 0) "fail" else if (nf_other > 0) "warning" else "ok",
    if (is.finite(nf)) nf + nf_other else NA_real_,
    if (is.finite(nf) && nf > 0) paste(nf, "observed outcome(s) never simulated at the optimum")
    else if (nf_other > 0)
      paste(nf_other, "profile point(s) had observed outcomes never simulated, even after repeating",
            "the re-evaluation with more trials; they are left out of the interval of z")
    else "all observed outcomes have simulated probability > 0",
    if ((is.finite(nf) && nf > 0) || nf_other > 0)
      "Increase n_sim; check atypical trials; check that the model can produce these outcomes." else "")

  rows$bound <- .rel_row(
    "Estimates away from search bounds",
    if (prof$at_bound[i]) "fail" else "ok", as.numeric(prof$at_bound[i]),
    if (prof$at_bound[i]) paste("on bound:", prof$bound_names[i])
    else if (length(fit$widened)) paste("bounds widened automatically:", paste(fit$widened, collapse = ", "))
    else "no parameter on a bound",
    if (prof$at_bound[i]) "Widen the bound with funres_control(bounds = ...) and refit; if it persists, the parameter is not identified by the data." else "")

  se <- prof$nll_se[i]
  rows$se <- .rel_row(
    "Monte Carlo noise of the NLL",
    if (is.finite(se) && se > th$se_warn) "warning" else "ok", se,
    sprintf("SE of the NLL at the optimum = %.3f (threshold %.2f)", se, th$se_warn),
    if (is.finite(se) && se > th$se_warn) "Increase n_sim (noise falls with 1/sqrt(n_sim))." else "")

  gap <- prof$nll[i] - prof$nll_opt[i]
  rows$gap <- .rel_row(
    "Optimiser value vs re-evaluated NLL",
    if (is.finite(gap) && gap > th$gap_warn) "warning" else "ok", gap,
    sprintf("re-evaluated minus optimiser NLL = %.2f", gap),
    if (is.finite(gap) && gap > th$gap_warn) "Increase n_sim: the optimiser is exploiting Monte Carlo noise." else "")

  lg <- prof$late_gain[i]
  rows$late <- .rel_row(
    "Optimiser settled before the end of the run",
    if (is.finite(lg) && lg > th$late_warn) "warning" else "ok", lg,
    sprintf("improvement in the last quarter of the run = %.2f NLL units", lg),
    if (is.finite(lg) && lg > th$late_warn) "Increase itermax (and NP)." else "")

  if (nrow(prof) >= 3) {
    # Interior points lower than both neighbours (beyond Monte Carlo noise),
    # other than the global minimum: a smooth profile has a single minimum,
    # so extra local minima mean that the optimiser failed at the
    # neighbouring values of z.
    n <- nrow(prof)
    tol <- pmax(1, 2 * ifelse(is.finite(prof$nll_se), prof$nll_se, 0))
    loc <- vapply(2:(n - 1), function(j)
      prof$nll[j] < prof$nll[j - 1] - tol[j] && prof$nll[j] < prof$nll[j + 1] - tol[j], logical(1))
    idx <- (2:(n - 1))[loc]
    idx <- setdiff(idx, which.min(prof$nll))
    rows$smooth <- .rel_row(
      "Profile of z regular (one minimum)",
      if (length(idx)) "fail" else "ok", length(idx),
      if (length(idx)) paste0("extra local minima at z = ", paste(signif(prof$z[idx], 3), collapse = ", "))
      else "no extra local minima",
      if (length(idx)) "The optimiser failed at some values of z (the profile is jagged): increase itermax, NP and n_sim, and refit; the interval of z is not reliable." else "")
  }
  if (!is.null(fit$ci)) {
    open <- c(if (fit$ci$low_censored) "lower", if (fit$ci$high_censored) "upper")
    rows$ci <- .rel_row(
      "Interval of z closed",
      if (length(open)) "warning" else "ok", length(open),
      if (length(open)) paste("open", paste(open, collapse = " and "), "limit") else "both limits found",
      if (length(open)) "Extend the z grid; an open limit after extension means the data do not bound z on that side (often with very small k)." else "")
  }
  if (!fit$model$search$user && is.null(fit$model$fix$k)) {
    k <- fit$par[["k"]]
    rows$k <- .rel_row(
      "Search dispersion k not extreme",
      if (k < 0.1) "warning" else "ok", k,
      sprintf("k = %.3g", k),
      if (k < 0.1) "Very irregular search (ridge between z, a and k): check the atypical-trial screening; consider fixing k (e.g. exponential search) and compare by AIC." else "")
  }
  if (identical(fit$z_method, "joint"))
    rows$joint <- .rel_row("Estimation of z by the validated method", "warning", NA_real_,
                           "joint estimation of z (experimental)",
                           "In validation the joint method under-covered the true z (67% vs 96% for the profile); confirm z and its interval with z_method = 'profile'.")
  if (fit$control$preset == "exploratory")
    rows$preset <- .rel_row("Settings suitable for reporting", "warning", NA_real_,
                            "exploratory preset",
                            "Refit with funres_control('standard') or 'final' before reporting.")
  if (!is.null(fit$verification)) rows <- c(rows, list(fit$verification))

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  lev <- c("ok", "warning", "fail")
  attr(out, "overall") <- lev[max(match(out$status, lev))]
  class(out) <- c("funres_reliability", class(out))
  out
}

#' @export
print.funres_reliability <- function(x, ...) {
  sym <- c(ok = "[ok]  ", warning = "[warn]", fail = "[FAIL]")
  cat("Reliability of the fit: ", toupper(attr(x, "overall")), "\n", sep = "")
  for (j in seq_len(nrow(x))) {
    cat(sym[[x$status[j]]], " ", x$check[j], ": ", x$detail[j], "\n", sep = "")
    if (nzchar(x$action[j])) cat("        -> ", x$action[j], "\n", sep = "")
  }
  invisible(x)
}

#' Verify a fit with additional optimisations
#'
#' Refits the model at the estimate of `z` (i) in `n_runs` independent runs
#' with the same settings, (ii) with more effort (`effort` times more
#' simulated trials, generations and population) and (iii) with search bounds
#' `bounds_factor` times wider on each side. A reliable optimum is reproduced
#' by independent runs, is not improved by more effort and does not move when
#' the bounds are widened. The results are added to [reliability()].
#'
#' Cost: about `n_runs + 2` optimisations at one value of z, the second of
#' them `effort^2` times more expensive (more simulations per evaluation and
#' more generations).
#'
#' @param fit A `funres_fit`.
#' @param n_runs Number of independent runs with the same settings.
#' @param effort Multiplier of `n_sim` and `itermax` in the high-effort run.
#' @param bounds_factor Factor by which each search bound is widened.
#' @return The fit with a `verification` table added (see [reliability()]) and
#'   the refits in `verification_runs`.
#' @export
verify_fit <- function(fit, n_runs = 2, effort = 2, bounds_factor = 10) {
  if (!inherits(fit, "funres_fit")) stop("Need a funres_fit.", call. = FALSE)
  th <- .rel_thresholds()
  z_hat <- fit$par[["z"]]
  base <- list(data_spp = fit$data, T_exp = fit$T, model = fit$model,
               H_ref = fit$H_ref)
  ctrl_hi <- fit$control
  ctrl_hi$n_sim <- as.integer(ctrl_hi$n_sim * effort)
  ctrl_hi$itermax <- as.integer(ctrl_hi$itermax * effort)
  ctrl_hi$steptol <- ctrl_hi$itermax
  wide <- fit$bounds
  for (nm in c("lambda", "h", "k", "z")) wide[[nm]] <- wide[[nm]] * c(1 / bounds_factor, bounds_factor)
  wide$cv_h[2] <- wide$cv_h[2] * bounds_factor

  jobs <- c(rep(list(list(ctrl = fit$control, bounds = fit$bounds)), n_runs),
            list(list(ctrl = ctrl_hi, bounds = fit$bounds)),
            list(list(ctrl = fit$control, bounds = wide)))
  if (.uses_cpp_sampler(fit$model) && !inherits(future::plan(), "sequential")) {
    old_plan <- future::plan(future::sequential)
    on.exit(future::plan(old_plan), add = TRUE)
  }
  fs <- lapply(jobs, function(jb)
    future(do.call(utils::getFromNamespace(".fit_point", "funresMech"),
                   c(list(z = z_hat, ctrl = jb$ctrl, bounds = jb$bounds), base)),
           seed = TRUE))
  res <- lapply(fs, value)

  i <- which.min(fit$profile$nll)
  nll0 <- fit$profile$nll[i]
  runs <- res[seq_len(n_runs)]
  nll_runs <- c(nll0, vapply(runs, function(r) r$nll_reeval, 0))
  spread <- max(nll_runs) - min(nll_runs)
  par_runs <- rbind(fit$par, do.call(rbind, lapply(runs, function(r) r$par)))
  rel_par <- max(apply(log(par_runs[, c("a", "h"), drop = FALSE]), 2,
                       function(v) exp(max(v) - min(v)) - 1))
  hi <- res[[n_runs + 1]]; wd <- res[[n_runs + 2]]
  gain_hi <- nll0 - hi$nll_reeval
  gain_wd <- nll0 - wd$nll_reeval

  rows <- list(
    .rel_row("Independent runs agree (NLL)",
             if (spread > th$rerun_fail) "fail" else if (spread > th$rerun_ok) "warning" else "ok",
             spread, sprintf("range of NLL over %d runs = %.2f", n_runs + 1, spread),
             if (spread > th$rerun_ok) "Increase itermax/NP and n_sim; report only when independent runs agree." else ""),
    .rel_row("Independent runs agree (a, h)",
             if (rel_par > th$par_warn) "warning" else "ok", rel_par,
             sprintf("largest relative difference of a or h between runs = %.0f%%", 100 * rel_par),
             if (rel_par > th$par_warn) "Parameters weakly identified or optimiser unstable: more effort; interpret a and h with caution." else ""),
    .rel_row("More effort does not improve the optimum",
             if (gain_hi > th$better_fail) "fail" else "ok", gain_hi,
             sprintf("NLL improvement with %gx effort = %.2f", effort, gain_hi),
             if (gain_hi > th$better_fail) "The settings were insufficient: refit the whole profile with more effort (preset 'final')." else ""),
    .rel_row("Wider bounds do not improve the optimum",
             if (gain_wd > th$better_fail) "fail" else "ok", gain_wd,
             sprintf("NLL improvement with bounds %gx wider = %.2f", bounds_factor, gain_wd),
             if (gain_wd > th$better_fail) "Refit with wider bounds (funres_control(bounds = ...))." else ""))
  fit$verification <- do.call(rbind, rows)
  fit$verification_runs <- res
  fit
}
