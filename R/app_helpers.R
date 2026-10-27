############################################################
# Helpers of the shiny application (pure functions, testable without shiny)
############################################################

#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_point geom_errorbar
#'   geom_hline geom_vline geom_histogram geom_boxplot geom_violin annotate
#'   facet_wrap labs theme_minimal theme element_text element_blank
#'   element_line scale_colour_manual scale_shape_manual scale_x_continuous
#'   scale_fill_manual coord_cartesian margin
#' @importFrom rlang .data
#' @importFrom stats quantile
#' @noRd
NULL

# Palette shared by the interface and the figures
.app_col <- list(ink = "#1C2B3A", slate = "#4A5A6A", accent = "#1F5F8B",
                 accent_soft = "#C9DBE8", grid = "#E3E8ED",
                 ok = "#2E6F57", warning = "#9A6B12", fail = "#9E3B3B",
                 data = "#1C2B3A", sim = "#1F5F8B")

# ggplot theme of the application (also used for downloaded figures)
.app_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(text = element_text(colour = .app_col$ink),
          plot.title = element_text(size = base_size + 1, face = "bold"),
          plot.subtitle = element_text(size = base_size - 1, colour = .app_col$slate),
          axis.title = element_text(size = base_size - 1, colour = .app_col$slate),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = .app_col$grid, linewidth = 0.4),
          legend.position = "bottom",
          plot.margin = margin(8, 12, 8, 8))
}

# ---------------------------------------------------------------------------
# Model and control from the inputs of the interface (a plain list, so that
# the same code serves the server, the exported R script and the tests)
# ---------------------------------------------------------------------------

# R sampler from code typed by the user. The code must evaluate to a
# function(n, shape). It runs in the user's own R session.
.sampler_from_code <- function(code, lower, upper, name) {
  f <- tryCatch(eval(parse(text = code), envir = new.env(parent = globalenv())),
                error = function(e) stop("The sampler code could not be evaluated: ",
                                         conditionMessage(e), call. = FALSE))
  if (!is.function(f)) stop("The sampler code must evaluate to function(n, shape).",
                            call. = FALSE)
  funres_sampler(f, lower = lower, upper = upper, name = name)
}

.app_model <- function(o) {
  search <- if (identical(o$search, "custom"))
    .sampler_from_code(o$search_code, o$search_lower, o$search_upper, "custom (R)")
  else o$search
  handling <- if (identical(o$handling, "custom"))
    .sampler_from_code(o$handling_code, o$handling_lower, o$handling_upper, "custom (R)")
  else o$handling
  fix <- list()
  if (isTRUE(o$fix_k) && !identical(o$search, "exponential")) fix$k <- o$k_value
  if (isTRUE(o$fix_s) && !identical(o$handling, "constant")) fix$s <- o$s_value
  if (identical(o$z_mode, "fixed")) fix$z <- o$z_value
  funres_model(search = search, handling = handling, process = o$process, fix = fix,
               unfinished = if (is.null(o$unfinished)) "include" else o$unfinished)
}

.app_bounds <- function(o) {
  if (!isTRUE(o$custom_bounds)) return(list())
  list(lambda = c(o$b_lambda_lo, o$b_lambda_hi), h = c(o$b_h_lo, o$b_h_hi),
       k = c(o$b_k_lo, o$b_k_hi), cv_h = c(o$b_cv_lo, o$b_cv_hi),
       z = c(o$b_z_lo, o$b_z_hi))
}

.app_control <- function(o) {
  funres_control(o$preset, n_sim = o$n_sim, itermax = o$itermax, NP = o$NP,
                 n_reeval = o$n_reeval, method = o$method,
                 bounds = .app_bounds(o), crn = isTRUE(o$crn),
                 polish = if (isTRUE(o$polish)) 200L else 0L,
                 widen = if (isTRUE(o$widen)) 2L else 0L)
}

.app_z_grid <- function(o) {
  g <- seq(o$z_min, o$z_max, by = o$z_step)
  if (length(g) < 2) stop("The z grid needs at least two values.", call. = FALSE)
  g
}

# Default bounds in the units of the user (for display)
.bounds_in_units <- function(T) {
  b <- .default_bounds()
  data.frame(
    parameter = c("Encounter rate at the reference density",
                  "Mean handling time h", "Search dispersion k",
                  "Handling CV (s / h)", "Exponent z (joint estimation)"),
    lower = c(b$lambda[1] / T, b$h[1] * T, b$k[1], b$cv_h[1], b$z[1]),
    upper = c(b$lambda[2] / T, b$h[2] * T, b$k[2], b$cv_h[2], b$z[2]),
    stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Run-time estimate: times two evaluations of the likelihood at a typical
# point and multiplies by the evaluations of the whole analysis.
# ---------------------------------------------------------------------------
.estimate_runtime <- function(data_spp, T, model, control, n_points, workers = 1) {
  H_ref <- .h_ref(data_spp)
  par <- c(a = 5 / (T * H_ref^1.5), h = 0.02 * T, z = 1.5, k = 1, s = 0.01 * T)
  for (nm in names(model$fix)) par[[nm]] <- model$fix[[nm]]
  if (model$handling$user && !(par[["s"]] > 0)) par[["s"]] <- 1
  t1 <- system.time(for (i in 1:2)
    .nll_sim(par, model, data_spp, T, min(control$n_sim, 2000L)))[["elapsed"]] / 2
  t_eval <- t1 * control$n_sim / min(control$n_sim, 2000L)
  evals <- control$itermax * control$NP + control$n_reeval
  # + re-evaluations, start of parallel workers and Beta-Binomial reference
  t_eval * evals * ceiling(n_points / max(1, workers)) + 10 + if (workers > 1) 5 else 0
}

.format_duration <- function(sec) {
  if (!is.finite(sec)) return("unknown")
  if (sec < 90) return(paste0(round(sec), " s"))
  if (sec < 5400) return(paste0(round(sec / 60), " min"))
  paste0(round(sec / 3600, 1), " h")
}

# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------
# singular of the usual time units ("per hour", not "per hours")
.unit1 <- function(u) {
  map <- c(minutes = "minute", hours = "hour", days = "day", seconds = "second", weeks = "week")
  if (u %in% names(map)) map[[u]] else u
}

.estimates_table <- function(fit) {
  u <- fit$time_unit
  p <- fit$par
  fixed <- names(fit$model$fix)
  data.frame(
    Parameter = c("a", "h", "z", "k", "s"),
    Meaning = c(paste0("Encounter-rate coefficient (per ", .unit1(u), " per ", .a_density_unit(fit), ")"),
                paste0("Mean handling time (", u, ")"),
                "Density-scaling exponent",
                if (fit$model$search$user) "Shape of search sampler" else "Search dispersion (CV^2 = 1/k)",
                if (fit$model$handling$user) "Shape of handling sampler" else paste0("SD of handling time (", u, ")")),
    Estimate = vapply(p, function(v) format(signif(v, 4)), ""),
    Status = ifelse(names(p) %in% fixed, "fixed", "estimated"),
    stringsAsFactors = FALSE)
}

.interval_text <- function(fit) {
  if (is.null(fit$ci)) return("z fixed")
  f <- function(v, side) if (is.na(v)) "open" else format(round(v, 2), nsmall = 2)
  paste0("[", f(fit$ci$z_low), ", ", f(fit$ci$z_high), "]")
}

.aic_table <- function(fit) {
  a <- fit$aic
  rows <- data.frame(
    Model = c("Mechanistic (this fit)", "Mechanistic with z = 1",
              "Beta-Binomial, same expected value"),
    `Free parameters` = c(fit$n_par, fit$n_par - 1,
                          if (is.null(fit$bb)) NA else fit$bb$n_par),
    AIC = c(a$AIC_full, a$AIC_restricted, a$AIC_BB),
    check.names = FALSE, stringsAsFactors = FALSE)
  rows$`Delta AIC` <- rows$AIC - a$AIC_full
  rows[is.finite(rows$AIC), , drop = FALSE]
}

# One row of the model comparison table
.comparison_row <- function(fit, label, species) {
  m <- fit$model
  rel <- reliability(fit)
  data.frame(
    Label = label, Species = species, Process = m$process,
    `End of trial` = if (.excludes_unfinished(m)) "exclude" else "include",
    Search = m$search$name, Handling = m$handling$name,
    Fixed = if (length(m$fix)) paste(names(m$fix), unlist(m$fix), sep = " = ", collapse = ", ") else "none",
    `Free parameters` = fit$n_par, NLL = round(fit$nll, 3), AIC = round(fit$aic$AIC_full, 2),
    z = signif(fit$par[["z"]], 3), `z interval` = .interval_text(fit),
    Reliability = attr(rel, "overall"),
    check.names = FALSE, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Figures (ggplot; the interface converts them with plotly)
# ---------------------------------------------------------------------------
.plot_profile <- function(fit) {
  prof <- fit$profile
  if (nrow(prof) < 2) {
    return(ggplot() + annotate("text", x = 0, y = 0, label = "z is fixed: no profile",
                               colour = .app_col$slate) + theme_void_like())
  }
  m <- min(prof$nll)
  prof$dnll <- prof$nll - m
  prof$lo <- prof$dnll - 2 * ifelse(is.finite(prof$nll_se), prof$nll_se, 0)
  prof$hi <- prof$dnll + 2 * ifelse(is.finite(prof$nll_se), prof$nll_se, 0)
  prof$point <- ifelse(prof$at_bound, "estimate on a search bound", "profile point")
  thr <- .ci95_threshold()
  g <- ggplot(prof, aes(x = .data$z, y = .data$dnll)) +
    geom_hline(yintercept = thr, linetype = "dashed", colour = .app_col$slate, linewidth = 0.4) +
    geom_errorbar(aes(ymin = .data$lo, ymax = .data$hi), width = 0, colour = .app_col$accent_soft,
                  linewidth = 0.8) +
    geom_line(colour = .app_col$accent, linewidth = 0.7) +
    geom_point(aes(shape = .data$point, colour = .data$point), size = 2.2) +
    scale_colour_manual(values = c("profile point" = .app_col$accent,
                                   "estimate on a search bound" = .app_col$fail), name = NULL) +
    scale_shape_manual(values = c("profile point" = 16, "estimate on a search bound" = 4), name = NULL) +
    geom_vline(xintercept = fit$par[["z"]], colour = .app_col$ink, linewidth = 0.4) +
    labs(x = "Density-scaling exponent z", y = "NLL minus minimum",
         subtitle = "Bars: +/- 2 Monte Carlo SE. Dashed line: 95% threshold (1.92).") +
    coord_cartesian(ylim = c(-0.5, max(thr * 3, min(max(prof$dnll) + 1, 15)))) +
    .app_theme()
  if (!any(prof$at_bound)) g <- g + theme(legend.position = "none")
  if (!is.null(fit$ci)) {
    for (v in c(fit$ci$z_low, fit$ci$z_high))
      if (!is.na(v)) g <- g + geom_vline(xintercept = v, colour = .app_col$slate,
                                         linetype = "dotted", linewidth = 0.4)
  }
  g
}

theme_void_like <- function() {
  theme_minimal() + theme(axis.text = element_blank(), axis.title = element_blank(),
                          panel.grid = element_blank())
}

# Simulated trials at the fitted parameters for the observed densities
.simulate_fit <- function(fit, n_sim = 2000) {
  simulate_response(fit$model, dens = sort(unique(fit$data$dens)), T = fit$T,
                    par = fit$par, n_sim = n_sim)
}

.plot_response <- function(fit, sims) {
  qs <- do.call(rbind, lapply(split(sims$response, sims$dens), function(v)
    data.frame(mean = mean(v), lo = stats::quantile(v, 0.025, names = FALSE),
               hi = stats::quantile(v, 0.975, names = FALSE))))
  qs$dens <- as.numeric(names(split(sims$response, sims$dens)))
  ggplot() +
    geom_ribbon(data = qs, aes(x = .data$dens, ymin = .data$lo, ymax = .data$hi),
                fill = .app_col$accent_soft, alpha = 0.7) +
    geom_line(data = qs, aes(x = .data$dens, y = .data$mean), colour = .app_col$accent,
              linewidth = 0.8) +
    geom_point(data = fit$data, aes(x = .data$dens, y = .data$par), colour = .app_col$ink,
               size = 1.8, alpha = 0.75) +
    labs(x = "Host (prey) density", y = if (fit$model$process == "predator") "Prey consumed" else "Hosts parasitized",
         subtitle = "Line: simulated mean. Band: central 95% of simulated trials. Points: observed trials.") +
    .app_theme()
}

.plot_distribution <- function(fit, sims, type = c("histogram", "boxplot", "violin")) {
  type <- match.arg(type)
  ylab <- if (fit$model$process == "predator") "Prey consumed" else "Hosts parasitized"
  if (type == "histogram") {
    sp <- split(sims$response, sims$dens)
    pr <- do.call(rbind, Map(function(v, x) {
      tb <- table(factor(v, levels = 0:max(v, 0)))
      data.frame(dens = as.numeric(x), response = as.integer(names(tb)),
                 prob = as.numeric(tb) / length(v))
    }, sp, names(sp)))
    pr$panel <- factor(paste("density", pr$dens), levels = paste("density", sort(unique(pr$dens))))
    obs <- fit$data
    obs$panel <- factor(paste("density", obs$dens), levels = levels(pr$panel))
    # one point per trial: trials with the same count are stacked below the axis
    obs$stack <- stats::ave(obs$par, obs$dens, obs$par, FUN = seq_along)
    top <- tapply(pr$prob, pr$dens, max)
    obs$y <- -obs$stack * 0.08 * top[as.character(obs$dens)]
    return(ggplot(pr, aes(x = .data$response, y = .data$prob)) +
      ggplot2::geom_col(fill = .app_col$accent_soft, width = 0.9) +
      geom_hline(yintercept = 0, colour = .app_col$slate, linewidth = 0.3) +
      geom_point(data = obs, aes(x = .data$par, y = .data$y), colour = .app_col$ink, size = 1.2,
                 inherit.aes = FALSE) +
      facet_wrap(~ panel, scales = "free") +
      labs(x = ylab, y = "Simulated probability",
           subtitle = paste("Bars: simulated distribution at the estimates.",
                            "Points below the axis: observed trials, one point per trial.")) +
      .app_theme(9))
  }
  sims$dens_f <- factor(sims$dens, levels = sort(unique(sims$dens)))
  obs <- fit$data
  obs$dens_f <- factor(obs$dens, levels = levels(sims$dens_f))
  geom_main <- if (type == "boxplot")
    geom_boxplot(fill = .app_col$accent_soft, colour = .app_col$accent, outlier.size = 0.6)
  else geom_violin(fill = .app_col$accent_soft, colour = .app_col$accent, scale = "width")
  ggplot(sims, aes(x = .data$dens_f, y = .data$response)) + geom_main +
    geom_point(data = obs, aes(x = .data$dens_f, y = .data$par), colour = .app_col$ink,
               size = 1.6, alpha = 0.8, position = ggplot2::position_jitter(width = 0.12, height = 0)) +
    labs(x = "Host (prey) density", y = ylab,
         subtitle = "Simulated distribution at the estimates and observed trials (points).") +
    .app_theme()
}


# ---------------------------------------------------------------------------
# R script that reproduces the analysis of the interface
# ---------------------------------------------------------------------------
.r_script <- function(o, file_name, species_col, dens_col, resp_col, species,
                      excluded_rows = integer(0)) {
  dq <- function(x) deparse(x, width.cutoff = 500L)
  model_lines <- function(kind) {
    if (identical(o[[kind]], "custom"))
      c(sprintf("%s_fun <- %s", kind, paste(trimws(o[[paste0(kind, "_code")]]), collapse = "\n")),
        sprintf("%s <- funres_sampler(%s_fun, lower = %s, upper = %s)", kind, kind,
                o[[paste0(kind, "_lower")]], o[[paste0(kind, "_upper")]]))
    else sprintf("%s <- %s", kind, dq(o[[kind]]))
  }
  fix <- c(if (isTRUE(o$fix_k) && !identical(o$search, "exponential")) paste0("k = ", o$k_value),
           if (isTRUE(o$fix_s) && !identical(o$handling, "constant")) paste0("s = ", o$s_value),
           if (identical(o$z_mode, "fixed")) paste0("z = ", o$z_value))
  b <- .app_bounds(o)
  has_sp <- nzchar(species_col)
  c("# Analysis generated by the funresMech application",
    paste0("# funresMech ", utils::packageVersion("funresMech"), ", ",
           format(Sys.time(), "%Y-%m-%d %H:%M")),
    "# Same model, settings, seed and order of groups as in the application.",
    "library(funresMech)",
    "",
    sprintf("dat <- read.csv(%s, check.names = FALSE)   # adjust the path to your file", dq(file_name)),
    "dat$row <- seq_len(nrow(dat))",
    if (length(excluded_rows))
      c("# trials excluded as atypical (rows of the file)",
        sprintf("dat <- dat[!dat$row %%in%% %s, ]", dq(as.integer(excluded_rows)))),
    "",
    model_lines("search"), model_lines("handling"),
    sprintf("model <- funres_model(search = search, handling = handling, process = %s%s%s)",
            dq(o$process), if (length(fix)) paste0(", fix = list(", paste(fix, collapse = ", "), ")") else "",
            if (identical(o$unfinished, "exclude")) ', unfinished = "exclude"' else ""),
    sprintf("control <- funres_control(%s, n_sim = %d, itermax = %d, NP = %d, n_reeval = %d, method = %s%s%s%s%s)",
            dq(o$preset), as.integer(o$n_sim), as.integer(o$itermax), as.integer(o$NP),
            as.integer(o$n_reeval), dq(o$method),
            if (length(b)) paste0(", bounds = ", dq(b)) else "",
            if (isTRUE(o$crn)) "" else ", crn = FALSE",
            if (isTRUE(o$polish)) "" else ", polish = 0",
            if (!isTRUE(o$widen)) ", widen = 0" else ""),
    "",
    if (isTRUE(o$parallel) && o$workers > 1)
      sprintf("future::plan(future::multisession, workers = %d)", as.integer(o$workers)),
    sprintf("set.seed(%d)", as.integer(o$seed)),
    "fits <- list()",
    sprintf("for (g in %s) {", dq(if (has_sp) species else "all trials")),
    if (has_sp) sprintf("  d <- dat[dat[[%s]] == g, ]", dq(species_col)) else "  d <- dat",
    sprintf("  fits[[g]] <- fit_funres(d, T = %s, model = model, control = control,", o$T_exp),
    sprintf("                          z_grid = seq(%s, %s, by = %s), z_method = %s, z_step = %s,",
            o$z_min, o$z_max, o$z_step, dq(if (identical(o$z_mode, "joint")) "joint" else "profile"),
            o$z_step),
    sprintf("                          extend_z = %s, dens = %s, response = %s, time_unit = %s)",
            if (isTRUE(o$extend_z)) "TRUE" else "FALSE", dq(dens_col), dq(resp_col), dq(o$time_unit)),
    "}",
    "",
    "lapply(fits, print)",
    "fits <- lapply(fits, verify_fit)   # independent runs, more effort, wider bounds",
    "lapply(fits, reliability)")
}
