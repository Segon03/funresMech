#' Run the funresMech Shiny App
#'
#' Launches the interactive application for mechanistic functional response
#' analysis, based on the Okuyama model.
#'
#' @details
#' The application fits the model to a table with one row per trial (species,
#' host density and number of parasitised hosts) and offers the following
#' tabs: fitted parameters, likelihood profile of the density-scaling exponent
#' `z` with its 95% interval, stochastic curves and distributions, diagnostics
#' (small `k`, parameters on a bound, open intervals), screening of atypical
#' trials, and an HTML report. The parameter `s` is the standard deviation of
#' the handling time on the natural scale (same units as `h`).
#'
#' The profile is computed in parallel with the 'future' package, so results
#' are reproducible with [set.seed()] before launching the app.
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
  app <- shiny::shinyApp(ui = ui, server = server)
  shiny::runApp(app, ...)
}
