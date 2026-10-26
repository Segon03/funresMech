############################################################
# User interface of the funresMech application (bslib, Bootstrap 5)
############################################################

#' @importFrom shiny actionButton checkboxInput conditionalPanel div downloadButton
#'   fileInput helpText htmlOutput numericInput p radioButtons selectInput
#'   selectizeInput tableOutput tags textAreaInput textInput uiOutput HTML
#'   verbatimTextOutput tagList span strong em
#' @importFrom plotly plotlyOutput
#' @noRd
NULL

# ---- theme ------------------------------------------------------------------
.app_bs_theme <- function() {
  th <- bslib::bs_theme(
    version = 5,
    bg = "#FFFFFF", fg = "#1C2B3A",
    primary = "#1F5F8B", secondary = "#4A5A6A", success = "#2E6F57",
    warning = "#9A6B12", danger = "#9E3B3B", info = "#1F5F8B",
    base_font = bslib::font_collection("Segoe UI", "Roboto", "Helvetica Neue",
                                       "Arial", "sans-serif"),
    code_font = bslib::font_collection("Consolas", "Menlo", "DejaVu Sans Mono",
                                       "monospace"),
    "font-size-base" = "0.875rem",
    "h1-font-size" = "1.3rem", "h2-font-size" = "1.15rem",
    "h3-font-size" = "1.05rem", "h4-font-size" = "0.975rem",
    "h5-font-size" = "0.925rem", "h6-font-size" = "0.875rem",
    "headings-font-weight" = 600,
    "border-radius" = "0.3rem",
    "navbar-padding-y" = "0.35rem",
    "card-cap-bg" = "#F4F6F8",
    "input-font-size" = "0.85rem",
    "btn-font-size" = "0.85rem"
  )
  bslib::bs_add_rules(th, "
    body { background: #F4F6F8; }
    .navbar { background-color: #1C2B3A !important; }
    .navbar .navbar-brand { font-weight: 600; letter-spacing: 0.01em; }
    .navbar .brand-sub { font-weight: 400; color: #AFC0CF; font-size: 0.8rem; margin-left: 0.6rem; }
    .navbar .nav-link { color: #D5DFE8 !important; font-size: 0.85rem; }
    .navbar .nav-link.active { color: #FFFFFF !important; border-bottom: 2px solid #8FB5D3; }
    .bslib-sidebar-layout > .sidebar { background: #FFFFFF; border-right: 1px solid #E3E8ED; }
    .card { border-color: #E3E8ED; box-shadow: none; }
    .card-header { font-weight: 600; font-size: 0.9rem; }
    .accordion-button { font-size: 0.875rem; font-weight: 600; padding: 0.55rem 0.8rem; }
    .accordion-body { padding: 0.6rem 0.8rem; }
    .form-label, label { font-size: 0.82rem; color: #33475B; margin-bottom: 0.15rem; }
    .form-group { margin-bottom: 0.6rem; }
    .help-block, .form-text { font-size: 0.78rem; color: #5B6B7A; }
    table { font-variant-numeric: tabular-nums; }
    .sidebar .table td:not(:first-child) { white-space: nowrap; }
    .table { font-size: 0.82rem; }
    .table th { color: #33475B; font-weight: 600; background: #F4F6F8; }
    textarea.code { font-family: Consolas, Menlo, 'DejaVu Sans Mono', monospace; font-size: 0.8rem; }
    .qmark { display: inline-block; width: 1.05em; height: 1.05em; line-height: 1.05em;
             border-radius: 50%; border: 1px solid #8A99A8; color: #5B6B7A; text-align: center;
             font-size: 0.7rem; margin-left: 0.3em; cursor: help; }
    .status-strip { display: flex; flex-wrap: wrap; gap: 0.5rem; }
    .status-item { flex: 1 1 15rem; background: #FFFFFF; border: 1px solid #E3E8ED;
                   border-left: 5px solid #8A99A8; border-radius: 0.3rem; padding: 0.5rem 0.75rem; }
    .status-item .sp { font-weight: 600; }
    .status-item .meta { font-size: 0.8rem; color: #4A5A6A; }
    .status-ok { border-left-color: #2E6F57; }
    .status-warning { border-left-color: #9A6B12; }
    .status-fail { border-left-color: #9E3B3B; }
    .pill { display: inline-block; padding: 0.05rem 0.45rem; border-radius: 0.6rem;
            font-size: 0.72rem; font-weight: 600; color: #FFFFFF; }
    .pill-ok { background: #2E6F57; } .pill-warning { background: #9A6B12; }
    .pill-fail { background: #9E3B3B; } .pill-none { background: #8A99A8; }
    .notice { border-left: 4px solid #9A6B12; background: #FBF6EC; padding: 0.5rem 0.75rem;
              font-size: 0.82rem; border-radius: 0.25rem; margin-bottom: 0.6rem; }
    .muted { color: #5B6B7A; font-size: 0.82rem; }
    .run-btn { width: 100%; margin-top: 0.4rem; }
    pre { font-size: 0.78rem; }
  ")
}

# Label with a help tooltip
.lab <- function(text, help) {
  span(text, bslib::tooltip(span(class = "qmark", "?"), help, placement = "right"))
}

.species_select <- function(id) {
  selectInput(id, "Species or group", choices = character(0))
}

# Color of the navigation bar: bslib >= 0.9.0 takes it through
# navbar_options(); older versions through the bg / inverse arguments.
.navbar_args <- function() {
  if (utils::packageVersion("bslib") >= "0.9.0") {
    list(navbar_options = bslib::navbar_options(bg = "#1C2B3A", theme = "dark"))
  } else {
    list(bg = "#1C2B3A", inverse = TRUE)
  }
}

ui <- function() {
  do.call(bslib::page_navbar, c(.navbar_args(), list(
    title = span("funresMech", span(class = "brand-sub",
                                     "Mechanistic functional response, simulation-based likelihood")),
    id = "nav",
    theme = .app_bs_theme(),
    window_title = "funresMech",
    fillable = FALSE,

    # ---- Data ---------------------------------------------------------------
    bslib::nav_panel(
      "Data",
      bslib::layout_sidebar(
        fillable = FALSE,
        sidebar = bslib::sidebar(
          width = 320, title = "Data file",
          fileInput("dataset", "CSV file, one row per trial", accept = ".csv"),
          uiOutput("column_select"),
          uiOutput("species_choice"),
          tags$hr(),
          strong("Atypical trials"),
          p(class = "muted", "Flags trials that are very unlikely under the rest of the data (leave-one-out Beta-Binomial, Bonferroni threshold). Nothing is removed unless you choose so in Model and fit."),
          actionButton("screen", "Screen trials", class = "btn-outline-primary btn-sm")
        ),
        bslib::card(bslib::card_header("Summary by species or group"),
                    tableOutput("data_summary")),
        bslib::card(bslib::card_header("Atypical trials"),
                    uiOutput("screen_ui")),
        bslib::card(bslib::card_header("First rows of the file"),
                    tableOutput("data_preview"))
      )
    ),

    # ---- Model and fit --------------------------------------------------------
    bslib::nav_panel(
      "Model and fit",
      bslib::layout_sidebar(
        fillable = FALSE,
        sidebar = bslib::sidebar(
          width = 380, title = "Model and settings",
          bslib::accordion(
            multiple = TRUE, open = c("Experiment", "Foraging model"),
            bslib::accordion_panel(
              "Experiment",
              selectInput("time_unit", "Time unit", c("minutes", "hours", "days", "other"),
                          selected = "hours"),
              conditionalPanel("input.time_unit == 'other'",
                               textInput("time_unit_other", "Name of the time unit", "time units")),
              numericInput("T_exp", .lab("Trial duration T", "Duration of each trial, in the time unit above. Handling times and their SD are reported in the same unit."),
                           value = 24, min = 0.001),
              checkboxInput("exclude_flagged", .lab("Exclude trials flagged as atypical",
                                                     "Refit without the trials flagged in Data. Use it as a sensitivity analysis and compare with the fit that keeps them."),
                            value = FALSE)
            ),
            bslib::accordion_panel(
              "Foraging model",
              radioButtons("process", .lab("Process", "Parasitoid: hosts stay in the arena and can be attacked again; the response is the number of distinct hosts attacked. Predator: each encounter consumes one prey, which is removed. Choose by the biology: a parasitoid that never revisits parasitized hosts behaves as the predator process."),
                           c("Parasitoid (no removal)" = "parasitoid",
                             "Predator (prey removed)" = "predator")),
              radioButtons("unfinished", .lab("Host or prey still in handling at T", "Include: an attack counts as soon as the host or prey is encountered, even if its handling ends after T (as in Okuyama 2026). Exclude: it counts only if handling ends by T. Choose the convention used to score the trials."),
                           c("Counted (include)" = "include",
                             "Not counted (exclude)" = "exclude")),
              selectInput("search", .lab("Search times", "Distribution of the time to the next encounter, with mean 1/(a x^z). Its dispersion k is defined by CV^2 = 1/k. Exponential = Gamma with k = 1."),
                          c("Gamma" = "gamma", "Exponential (k = 1)" = "exponential",
                            "Lognormal" = "lognormal", "Weibull" = "weibull",
                            "Custom (R code)" = "custom")),
              conditionalPanel(
                "input.search == 'custom'",
                textAreaInput("search_code", .lab("Sampler of search times", "R code that evaluates to function(n, shape) returning n non-negative values with mean 1. The shape parameter is estimated within the bounds below and reported as k."),
                              value = "function(n, shape) rgamma(n, shape = shape, rate = shape)",
                              rows = 3, resize = "vertical"),
                div(class = "row g-2",
                    div(class = "col", numericInput("search_lower", "Shape, lower", 0.01, min = 1e-8)),
                    div(class = "col", numericInput("search_upper", "Shape, upper", 100, min = 1e-8))),
                actionButton("test_search", "Test sampler", class = "btn-outline-secondary btn-sm"),
                tableOutput("test_search_out")
              ),
              selectInput("handling", .lab("Handling times", "Distribution of the handling time, with mean h and standard deviation s (same unit as T). Constant = s fixed at 0."),
                          c("Lognormal" = "lognormal", "Gamma" = "gamma", "Weibull" = "weibull",
                            "Constant (s = 0)" = "constant", "Custom (R code)" = "custom")),
              conditionalPanel(
                "input.handling == 'custom'",
                textAreaInput("handling_code", .lab("Sampler of handling times", "R code that evaluates to function(n, shape) returning n non-negative values with mean 1; they are multiplied by h. The shape parameter is estimated within the bounds below and reported as s."),
                              value = "function(n, shape) rgamma(n, shape = shape, rate = shape)",
                              rows = 3, resize = "vertical"),
                div(class = "row g-2",
                    div(class = "col", numericInput("handling_lower", "Shape, lower", 0.1, min = 1e-8)),
                    div(class = "col", numericInput("handling_upper", "Shape, upper", 1000, min = 1e-8))),
                actionButton("test_handling", "Test sampler", class = "btn-outline-secondary btn-sm"),
                tableOutput("test_handling_out")
              ),
              conditionalPanel(
                "input.search != 'exponential'",
                checkboxInput("fix_k", .lab("Fix k", "Hold the search dispersion (or the shape of a custom search sampler) at a known value instead of estimating it. Fewer free parameters make the optimization more stable."), FALSE),
                conditionalPanel("input.fix_k", numericInput("k_value", "k", 1, min = 1e-8))
              ),
              conditionalPanel(
                "input.handling != 'constant'",
                checkboxInput("fix_s", .lab("Fix s", "Hold the SD of handling time (or the shape of a custom handling sampler) at a known value."), FALSE),
                conditionalPanel("input.fix_s", numericInput("s_value", "s", 0, min = 0))
              )
            ),
            bslib::accordion_panel(
              "Exponent z",
              radioButtons("z_mode", .lab("Estimation of z", "Profile: z fixed at each grid value and the others optimized; the estimate is the minimum. Joint (experimental): z estimated with the other parameters, then profiled only around the estimate; in validation its intervals under-covered the true z (67% vs 96% for the profile). Fixed: z held at a value (e.g. 1 for a type II response)."),
                           c("Profile over a grid" = "profile", "Joint, then profile (experimental)" = "joint",
                             "Fixed" = "fixed")),
              conditionalPanel(
                "input.z_mode == 'profile'",
                div(class = "row g-2",
                    div(class = "col", numericInput("z_min", "From", 0.5, step = 0.1)),
                    div(class = "col", numericInput("z_max", "To", 5, step = 0.1)),
                    div(class = "col", numericInput("z_step", "Step", 0.25, min = 0.01, step = 0.05))),
                checkboxInput("extend_z", .lab("Extend the grid if a limit is open", "Evaluates z = 8, 10, 15 and 20 when the profile does not cross the 95% threshold above the estimate, and z = 0.25 and 0.1 when it does not cross it below."), TRUE)
              ),
              conditionalPanel("input.z_mode == 'joint'",
                               numericInput("z_step_joint", "Initial step of the profile", 0.25, min = 0.01)),
              conditionalPanel("input.z_mode == 'fixed'",
                               numericInput("z_value", "z", 1, min = 0.01))
            ),
            bslib::accordion_panel(
              "Optimization",
              selectInput("preset", .lab("Preset", "Exploratory: a quick first look, not for reporting. Standard: default. Final: the settings used in published analyses with this method (10000 simulations, 200 generations). Values below can be edited."),
                          c("Exploratory" = "exploratory", "Standard" = "standard", "Final" = "final"),
                          selected = "standard"),
              div(class = "row g-2",
                  div(class = "col-6", numericInput("n_sim", .lab("n_sim", "Simulated trials per density for each evaluation of the likelihood. More = less Monte Carlo noise."), 5000, min = 100, step = 500)),
                  div(class = "col-6", numericInput("n_reeval", .lab("Re-evaluations", "Fresh evaluations of the likelihood at each optimum; their mean is used for the profile and the interval."), 5, min = 1))),
              div(class = "row g-2",
                  div(class = "col-6", numericInput("itermax", .lab("Generations", "itermax: iterations of the optimizer."), 100, min = 5)),
                  div(class = "col-6", numericInput("NP", .lab("Population", "NP: candidate solutions per generation."), 40, min = 4))),
              selectInput("method", .lab("Method", "Global optimizer, with the same budget of evaluations (generations x population)."),
                          c("Differential evolution (DEoptim)" = "DEoptim",
                            "CMA-ES (package cmaes)" = "CMAES",
                            "CRS2 (package nloptr)" = "CRS2")),
              checkboxInput("widen", .lab("Widen bounds hit at the optimum", "If the optimum ends on a search bound, the bound is widened tenfold and the affected points refitted (up to twice)."), TRUE),
              checkboxInput("crn", .lab("Common random numbers", "The same random seed in every evaluation during the search, so the objective is deterministic (but stepwise). Re-evaluations use fresh numbers. Default: in validation it made the optimizer far more reliable."), TRUE),
              checkboxInput("polish", .lab("Polish the optimum (Nelder-Mead)", "200 extra evaluations of a local search from the best point. Default, together with common random numbers."), TRUE),
              numericInput("seed", .lab("Random seed", "Same seed and settings give the same results."), 1, min = 1, step = 1),
              bslib::input_switch("parallel", "Parallel computation", TRUE),
              conditionalPanel("input.parallel",
                               numericInput("workers", "Cores", max(1L, parallel::detectCores() - 1L), min = 1,
                                            max = parallel::detectCores()))
            ),
            bslib::accordion_panel(
              "Search bounds",
              p(class = "muted", "Defaults scale with T, so they do not depend on the time unit. Shown in your units:"),
              tableOutput("default_bounds"),
              checkboxInput("custom_bounds", "Set the bounds myself", FALSE),
              conditionalPanel(
                "input.custom_bounds",
                .bound_inputs("b_lambda", "Encounter rate at the reference density (per time unit)", 0.0001, 100),
                .bound_inputs("b_h", "Mean handling time h (time units)", 1e-5, 24),
                .bound_inputs("b_k", "Search dispersion k", 1e-5, 100),
                .bound_inputs("b_cv", "Handling CV s/h", 0, 10),
                .bound_inputs("b_z", "z (joint estimation)", 0.1, 20)
              )
            )
          ),
          actionButton("estimate", "Estimate run time", class = "btn-outline-secondary btn-sm run-btn"),
          actionButton("run", "Fit the model", class = "btn-primary run-btn")
        ),
        bslib::card(bslib::card_header("Specification"),
                    uiOutput("spec_ui")),
        bslib::card(bslib::card_header("Fits"),
                    uiOutput("status_strip"),
                    uiOutput("fit_messages"))
      )
    ),

    # ---- Reliability ----------------------------------------------------------
    bslib::nav_panel(
      "Reliability",
      bslib::layout_sidebar(
        fillable = FALSE,
        sidebar = bslib::sidebar(
          width = 320, title = "Verification",
          .species_select("rel_species"),
          p(class = "muted", "Refits at the estimate of z with independent runs, with more effort and with wider bounds. A reliable optimum is reproduced, not improved by more effort and not moved by wider bounds."),
          numericInput("v_runs", "Independent runs", 2, min = 1, max = 10),
          numericInput("v_effort", "Effort multiplier", 2, min = 1.2, step = 0.5),
          numericInput("v_bounds", "Bounds widened by", 10, min = 2),
          checkboxInput("v_all", "Verify all species or groups", FALSE),
          actionButton("verify", "Run verification", class = "btn-primary run-btn")
        ),
        bslib::card(bslib::card_header("Status of each fit"), uiOutput("status_strip2")),
        bslib::card(bslib::card_header("Checks and actions"), uiOutput("reliability_ui")),
        bslib::card(bslib::card_header("Likelihood profile of z"),
                    plotlyOutput("profile_plot_rel", height = "360px"))
      )
    ),

    # ---- Results --------------------------------------------------------------
    bslib::nav_panel(
      "Results",
      bslib::layout_sidebar(
        fillable = FALSE,
        sidebar = bslib::sidebar(
          width = 300, title = "Display",
          .species_select("res_species"),
          numericInput("plot_nsim", .lab("Simulated trials for figures", "Trials simulated at the estimates for each density, used only for the figures."), 2000, min = 200, step = 500),
          selectInput("dist_type", "Distribution view", c("Histogram" = "histogram", "Box plot" = "boxplot",
                                                          "Violin" = "violin")),
          tags$hr(),
          selectInput("fig_format", "Figure format", c("PDF" = "pdf", "PNG (300 dpi)" = "png")),
          downloadButton("dl_fig_profile", "Profile", class = "btn-outline-secondary btn-sm"),
          downloadButton("dl_fig_response", "Response", class = "btn-outline-secondary btn-sm"),
          downloadButton("dl_fig_dist", "Distributions", class = "btn-outline-secondary btn-sm")
        ),
        uiOutput("prelim_notice"),
        bslib::layout_columns(
          col_widths = c(6, 6),
          bslib::card(bslib::card_header("Estimates"), uiOutput("estimates_ui")),
          bslib::card(bslib::card_header("Model comparison by AIC"), uiOutput("aic_ui"))
        ),
        bslib::layout_columns(
          col_widths = c(6, 6),
          bslib::card(bslib::card_header("Likelihood profile of z"),
                      plotlyOutput("profile_plot", height = "340px")),
          bslib::card(bslib::card_header("Functional response"),
                      plotlyOutput("response_plot", height = "340px"))
        ),
        bslib::card(bslib::card_header("Distribution of the response at each density"),
                    plotlyOutput("dist_plot", height = "420px"))
      )
    ),

    # ---- Model comparison ---------------------------------------------------------
    bslib::nav_panel(
      "Model comparison",
      bslib::layout_sidebar(
        fillable = FALSE,
        sidebar = bslib::sidebar(
          width = 300, title = "Add fits",
          .species_select("cmp_species"),
          textInput("cmp_label", "Label", "Model 1"),
          actionButton("cmp_add", "Add current fit", class = "btn-primary btn-sm"),
          actionButton("cmp_clear", "Clear table", class = "btn-outline-secondary btn-sm"),
          tags$hr(),
          p(class = "muted", "Fit the same data with other distributions, processes or fixed parameters, add each fit here and compare them by AIC. Compare only fits of the same trials (same species and the same exclusions) and only fits whose reliability is not 'fail'.")
        ),
        bslib::card(bslib::card_header("Fits compared"), uiOutput("cmp_ui"))
      )
    ),

    # ---- Export ---------------------------------------------------------------
    bslib::nav_panel(
      "Export",
      bslib::layout_columns(
        col_widths = c(6, 6),
        bslib::card(
          bslib::card_header("Tables"),
          downloadButton("dl_estimates", "Estimates (CSV)", class = "btn-outline-primary btn-sm"),
          downloadButton("dl_profiles", "Likelihood profiles (CSV)", class = "btn-outline-primary btn-sm"),
          downloadButton("dl_reliability", "Reliability checks (CSV)", class = "btn-outline-primary btn-sm"),
          downloadButton("dl_comparison", "Model comparison (CSV)", class = "btn-outline-primary btn-sm"),
          downloadButton("dl_fits", "Fitted objects (RDS)", class = "btn-outline-primary btn-sm")
        ),
        bslib::card(
          bslib::card_header("Reproducibility and report"),
          p(class = "muted", "The R script repeats the analysis with the same model, settings and seed outside the application."),
          downloadButton("dl_script", "R script", class = "btn-outline-primary btn-sm"),
          checkboxInput("report_data", "Include the data in the report", FALSE),
          downloadButton("dl_report", "HTML report", class = "btn-primary btn-sm")
        )
      )
    ),

    # ---- Help -----------------------------------------------------------------
    bslib::nav_panel("Help", .help_panel()),

    bslib::nav_spacer(),
    bslib::nav_item(span(class = "small", style = "color:#AFC0CF;",
                         paste("version", utils::packageVersion("funresMech"))))
  )))
}

.bound_inputs <- function(id, label, lo, hi) {
  tagList(
    tags$div(class = "form-label", label),
    div(class = "row g-2",
        div(class = "col", numericInput(paste0(id, "_lo"), "lower", lo)),
        div(class = "col", numericInput(paste0(id, "_hi"), "upper", hi)))
  )
}

.help_panel <- function() {
  bslib::layout_columns(
    col_widths = c(5, 7),
    bslib::card(
      bslib::card_header("Workflow"),
      tags$ol(
        tags$li("Data: load a CSV with one row per trial (host or prey density and the number parasitized or consumed; optionally a species or group column). Screen for atypical trials."),
        tags$li("Model and fit: set T and its unit, choose the process and the distributions of search and handling times (built-in or your own), decide which parameters to fix, how to estimate z, and the optimization settings. Estimate the run time, then fit."),
        tags$li("Reliability: read the checks of each fit and run the verification. Act on every warning or failure before using the results."),
        tags$li("Results: estimates, interval of z, AIC against z = 1 and against a Beta-Binomial model, profile and simulated response."),
        tags$li("Model comparison: fit alternative specifications of the foraging process to the same data and compare them by AIC."),
        tags$li("Export: tables, fitted objects, an R script that repeats the analysis and an HTML report.")
      ),
      tags$h6("The model"),
      p("The likelihood of each observed count is the frequency of that count in trials simulated from the foraging process (Okuyama 2012, 2026): a consumer alternates search, with mean time 1/(a x^z) at density x, and handling, with mean time h, until the trial duration T is reached. z = 1 gives a type II response and z > 1 a type III response. The distributions of search and handling times are assumptions; no family is correct by default, and the application lets you choose them, write your own and compare the fits."),
      tags$h6("Defaults are a starting point"),
      p("A simulated likelihood is noisy and the optimizer can stop at a suboptimal point. The exploratory preset is for a first look only. Report results only from fits whose checks pass, preferably verified, and state the settings used.")
    ),
    bslib::card(
      bslib::card_header("Reading the reliability checks"),
      tags$table(class = "table table-sm",
        tags$thead(tags$tr(tags$th("Symptom"), tags$th("Likely cause"), tags$th("Action"))),
        tags$tbody(
          tags$tr(tags$td("Observed outcomes never simulated"), tags$td("Parameters far from the data, too few simulations, or an outcome the model cannot produce"), tags$td("Increase n_sim; check atypical trials; consider another process or distribution")),
          tags$tr(tags$td("Estimate on a search bound"), tags$td("Bound too narrow, or parameter not identified"), tags$td("Widen the bound and refit; if it persists, fix the parameter or interpret it as a limit")),
          tags$tr(tags$td("Large Monte Carlo SE or optimizer gap"), tags$td("Too few simulated trials"), tags$td("Increase n_sim")),
          tags$tr(tags$td("Still improving at the end of the run"), tags$td("Too few generations"), tags$td("Increase generations and population")),
          tags$tr(tags$td("Independent runs disagree"), tags$td("Optimizer stops at different points"), tags$td("More effort; do not report until runs agree")),
          tags$tr(tags$td("More effort or wider bounds improve the optimum"), tags$td("Settings insufficient"), tags$td("Refit the whole profile with the final preset or wider bounds")),
          tags$tr(tags$td("Profile of z with several minima"), tags$td("The optimizer failed at some values of z"), tags$td("Increase generations, population and n_sim; do not use the interval until the profile is regular")),
          tags$tr(tags$td("Open interval of z"), tags$td("The data do not bound z on that side (often with very small k)"), tags$td("Extend the grid; check atypical trials; report the limit as open")),
          tags$tr(tags$td("Very small k"), tags$td("Very irregular search, ridge between z, a and k"), tags$td("Check atypical trials; compare with k fixed (exponential search) by AIC"))
        )
      ),
      tags$h6("References"),
      p(class = "muted", "Okuyama T (2012) A likelihood approach for functional response models. Biological Control 60: 103-107. Okuyama T (2026) Parametric assumptions in parasitoid functional response analysis. Journal of Applied Entomology, doi:10.1111/jen.70148.")
    )
  )
}
