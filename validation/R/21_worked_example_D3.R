# ============================================================================
# 21_worked_example_D3.R
#
# Complete analysis of a published dataset, step by step, as shown in the
# worked example of the article: dataset D3 (third-instar hosts of
# Kishani Farahani and Goldansaz 2013, taken from the FoRAGE database; 64
# trials of 24 h). Each step saves what it prints, so the inputs, commands,
# outputs and figures of the example can be reproduced exactly.
#
#   1. data                validation/worked_example/01_data.txt
#   2. screening           validation/worked_example/02_screening.txt
#   3. model               validation/worked_example/03_model.txt
#   4. fit                 validation/worked_example/04_fit.txt
#   5. reliability checks  validation/worked_example/05_reliability.txt
#   6. verification        validation/worked_example/06_verification.txt
#   7. figures             validation/worked_example/fig_profile.png,
#                          fig_response.png, fig_distribution.png
#   8. end-of-trial convention (attacks whose handling ends after T not
#      counted), compared by AIC
#                          validation/worked_example/08_unfinished.txt
#
# Usage, from the root of the package repository with funresMech installed:
#   QUICK <- TRUE;  source("validation/R/21_worked_example_D3.R")  # check, a few minutes
#   rm(list = ls()); source("validation/R/21_worked_example_D3.R")   # complete run
# The complete run uses the standard preset and takes a few hours on a
# multicore computer (the profile of z is evaluated in parallel).
# ============================================================================

library(funresMech)
stopifnot(packageVersion("funresMech") >= "1.1.0.9000")

if (!exists("QUICK")) QUICK <- FALSE
WORKERS <- if (QUICK) 2L else max(1L, parallel::detectCores() - 1L)
DIR_OUT <- file.path("validation", if (QUICK) "worked_example_quick" else "worked_example")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

control <- if (QUICK) funres_control(n_reeval = 3, n_sim = 300, itermax = 8, NP = 16) else
  funres_control()
z_grid <- if (QUICK) c(0.5, 1, 2, 3) else seq(0.5, 5, by = 0.25)   # the default grid of fit_funres()

save_output <- function(file, expr) {
  out <- utils::capture.output(expr)
  writeLines(out, file.path(DIR_OUT, file))
  cat(out, sep = "\n")
  invisible(out)
}
t0 <- Sys.time()

# ---- 1. data -------------------------------------------------------------------
# One row per trial: host density (dens) and hosts parasitized (par).
d3 <- read.csv(file.path("validation", "data", "D3.csv"))
save_output("01_data.txt", {
  cat("Trials:", nrow(d3), " Densities:", length(unique(d3$dens)), "\n\n")
  print(head(d3))
  cat("\nHosts parasitized by density (mean, range):\n")
  print(do.call(rbind, lapply(split(d3$par, d3$dens), function(v)
    data.frame(trials = length(v), mean = round(mean(v), 2), min = min(v), max = max(v)))))
})

# ---- 2. screening of atypical trials -------------------------------------------
save_output("02_screening.txt", {
  scr <- screen_outliers(d3, T = 24)
  cat("Trials flagged as atypical:", sum(scr$flagged), "\n")
  if (any(scr$flagged)) print(scr[scr$flagged, ])
})

# ---- 3. model ------------------------------------------------------------------
model <- funres_model()          # Gamma search, lognormal handling, parasitoid process
save_output("03_model.txt", { print(model); cat("\n"); print(control) })

# ---- 4. fit --------------------------------------------------------------------
future::plan(future::multisession, workers = WORKERS)
set.seed(2026)
fit <- fit_funres(d3, T = 24, time_unit = "h", model = model, z_grid = z_grid, control = control)
save_output("04_fit.txt", print(fit))

# ---- 5. reliability checks (cheap) ---------------------------------------------
save_output("05_reliability.txt", print(reliability(fit)))

# ---- 6. verification: independent runs, more effort, wider bounds --------------
fit <- verify_fit(fit)
save_output("06_verification.txt", { print(fit); cat("\n"); print(reliability(fit)) })
saveRDS(fit, file.path(DIR_OUT, "fit_D3.rds"))

# ---- 7. figures ----------------------------------------------------------------
set.seed(1)
for (ty in c("profile", "response", "distribution"))
  ggplot2::ggsave(file.path(DIR_OUT, paste0("fig_", ty, ".png")), plot(fit, type = ty),
                  width = if (ty == "distribution") 9 else 6, height = if (ty == "distribution") 7 else 4.5,
                  dpi = 300, bg = "white")

# ---- 8. end-of-trial convention ------------------------------------------------
set.seed(2026)
fit_exc <- fit_funres(d3, T = 24, time_unit = "h", model = funres_model(unfinished = "exclude"),
                      z_grid = z_grid, control = control)
save_output("08_unfinished.txt", {
  print(fit_exc)
  cat(sprintf("\nAIC, attacks counted at the encounter (include): %.2f\n", fit$aic$AIC_full))
  cat(sprintf("AIC, attacks counted when handling ends by T (exclude): %.2f\n", fit_exc$aic$AIC_full))
  cat(sprintf("Difference (exclude - include): %.2f\n", fit_exc$aic$AIC_full - fit$aic$AIC_full))
})
saveRDS(fit_exc, file.path(DIR_OUT, "fit_D3_exclude.rds"))
future::plan(future::sequential)

cat(sprintf("\nDone in %.1f hours. Outputs in %s\n",
            as.numeric(difftime(Sys.time(), t0, units = "hours")), DIR_OUT))
