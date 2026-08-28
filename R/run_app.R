#' Run the funresMech Shiny App
#'
#' Launches the interactive application for mechanistic functional response
#' analysis, based on the Okuyama model.
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
