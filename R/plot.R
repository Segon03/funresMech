#' Plot a fitted functional response
#'
#' Draws the likelihood profile of `z`, the fitted functional response, or the
#' distribution of the counts predicted at each density, with the observed
#' trials. These are the same figures as in the interface ([run_app()]).
#'
#' * `"profile"`: negative log-likelihood (NLL) at each value of `z`, minus its
#'   minimum, with bars of plus or minus two Monte Carlo standard errors. The
#'   95% interval of `z` contains the values below the dashed line (1.92).
#' * `"response"`: mean response and central 95% of the trials simulated at the
#'   estimates, with the observed trials.
#' * `"distribution"`: simulated probability of each count at each density,
#'   with the observed counts. It shows whether the model reproduces the
#'   spread of the data, not only its mean.
#'
#' @param x A fit returned by [fit_funres()].
#' @param type `"profile"` (default), `"response"` or `"distribution"`.
#' @param n_sim Number of simulated trials per density for `"response"` and
#'   `"distribution"`.
#' @param ... Not used.
#' @return A ggplot object.
#' @seealso [fit_funres()], [simulate_response()]
#' @examples
#' \dontrun{
#' fit <- fit_funres(d, T = 24)
#' plot(fit)                       # profile of z
#' plot(fit, type = "response")
#' plot(fit, type = "distribution")
#' }
#' @export
plot.funres_fit <- function(x, type = c("profile", "response", "distribution"),
                            n_sim = 2000, ...) {
  type <- match.arg(type)
  switch(type,
         profile = .plot_profile(x),
         response = .plot_response(x, .simulate_fit(x, n_sim)),
         distribution = .plot_distribution(x, .simulate_fit(x, n_sim), "histogram"))
}
