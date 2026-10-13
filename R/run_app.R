#' Run the funresMech Shiny App
#'
#' Launches the interactive application for mechanistic functional response
#' analysis by simulation-based likelihood.
#'
#' @details
#' The application covers the whole analysis: loading data with one row per
#' trial (and an optional species or group column) and screening atypical
#' trials; specifying the foraging model (parasitoid or predator process,
#' distributions of search and handling times, including samplers written in
#' R, fixed parameters, profile, joint or fixed estimation of `z`, optimiser,
#' Monte Carlo effort and search bounds); fitting with an estimate of the run
#' time; checking and verifying the reliability of each fit; results
#' (estimates, interval of `z`, AIC against `z = 1` and against a
#' Beta-Binomial model, likelihood profile, simulated response and
#' distributions); comparing alternative models by AIC; and exporting tables,
#' fitted objects, an R script that repeats the analysis and an HTML report.
#'
#' Results are reproducible: the seed is set in the application. Samplers
#' written in C++ are available from R code ([funres_sampler()]), not from the
#' application.
#'
#' @references
#' Okuyama, T. (2012). A likelihood approach for functional response models.
#' *Biological Control*, 60(2), 103-107. \doi{10.1016/j.biocontrol.2011.10.008}
#'
#' Okuyama, T. (2026). Parametric assumptions in parasitoid functional
#' response analysis. *Journal of Applied Entomology*.
#' \doi{10.1111/jen.70148}
#'
#' @param ... Additional arguments passed to [shiny::runApp()]
#'            (e.g., `port`, `host`, `launch.browser`).
#'
#' @return This function launches a Shiny app and does not return a value.
#'
#' @examples
#' if (interactive()) {
#'   run_app()
#' }
#'
#' @export
#' @import shiny
run_app <- function(...) {
  # Cargar los componentes UI y server
  app <- shiny::shinyApp(ui = ui(), server = server)
  shiny::runApp(app, ...)
}
