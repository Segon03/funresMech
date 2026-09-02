############################################################
# UI
############################################################

#' UI Definition for funresMech App
#'
#' @return A Shiny UI object
#' @importFrom shiny fluidPage tags HTML sidebarLayout sidebarPanel mainPanel tabsetPanel tabPanel actionButton numericInput fileInput checkboxInput checkboxGroupInput downloadButton br hr h2 h4 h5 p conditionalPanel
#' @importFrom shinythemes shinytheme
#' @importFrom shinyBS bsTooltip
#' @importFrom plotly plotlyOutput
#' @noRd
NULL

# CSS definition for the UI (compact styles)
compact_css <- "
  body, label, input, button, .selectize-input, .selectize-dropdown, .form-control {
    font-size: 13px !important;
  }
  h4, h5 {
    font-size: 15px !important;
    font-weight: 600;
  }
  .btn {
    padding: 4px 10px !important;
    font-size: 13px !important;
  }
  .shiny-output-error-validation {
    color: #d9534f;
    font-size: 13px;
  }
"

ui <- fluidPage(
  theme = shinytheme("flatly"),
  tags$head(tags$style(HTML(compact_css))),

  titlePanel("Okuyama Mechanistic Model for Functional Response"),

  sidebarLayout(
    sidebarPanel(
      fileInput("dataset", "Upload CSV dataset", accept = ".csv"),

      uiOutput("column_select"),

      numericInput("n_species", "Number of species to analyze",
                   value = 1, min = 1),

      numericInput("T_exp", "Experiment duration (hours)",
                   value = 1, min = 0.1, step = 0.1),

      actionButton("run", "Run analysis", class = "btn-primary btn-sm"),

      hr(),

      tags$h4("Advanced settings"),
      actionButton("toggle_advanced", "Show / Hide advanced settings",
                   class = "btn-info btn-sm"),

      conditionalPanel(
        condition = "input.toggle_advanced % 2 == 1",
        hr(),
        tags$h5("Optimization (DEoptim)"),
        numericInput("itermax", "Max iterations (itermax)",
                     value = 50, min = 5, max = 200),
        bsTooltip("itermax",
                  "Maximum number of DEoptim iterations. Higher values improve convergence but increase runtime.",
                  "right", options = list(container = "body")),

        numericInput("NP", "Population size (NP)",
                     value = 40, min = 10, max = 200),
        bsTooltip("NP",
                  "Number of candidate solutions in the DEoptim population. Larger NP explores more but is slower.",
                  "right", options = list(container = "body")),

        numericInput("reltol", "Relative tolerance (reltol)",
                     value = 1e-2, min = 1e-5, max = 1e-1),
        bsTooltip("reltol",
                  "Stopping criterion for DEoptim. Smaller values require more precision and more time.",
                  "right", options = list(container = "body")),

        hr(),
        tags$h5("Stochastic simulations"),
        numericInput("n_sim", "Simulations per density (n_sim)",
                     value = 3000, min = 500, max = 10000),
        bsTooltip("n_sim",
                  "Number of stochastic simulations per density. Higher values give smoother estimates but are slower.",
                  "right", options = list(container = "body")),

        hr(),
        tags$h5("Likelihood profile for z"),
        numericInput("z_min", "Minimum z",
                     value = 0.5, min = 0.1, max = 5),
        numericInput("z_max", "Maximum z",
                     value = 1.5, min = 0.2, max = 5),
        numericInput("z_step", "Step size for z",
                     value = 0.05, min = 0.01, max = 0.5),
        bsTooltip("z_min",
                  "Lower bound of z for the likelihood profile.",
                  "right", options = list(container = "body")),
        bsTooltip("z_max",
                  "Upper bound of z for the likelihood profile.",
                  "right", options = list(container = "body")),
        bsTooltip("z_step",
                  "Resolution of the z grid. Smaller steps give finer profiles but are slower.",
                  "right", options = list(container = "body")),

        hr(),
        tags$h5("Parallelization"),
        checkboxInput("use_parallel", "Use parallel computation (future::multisession)",
                      value = TRUE),
        numericInput("cores", "Number of CPU cores",
                     value = max(1, parallel::detectCores() - 1),
                     min = 1, max = parallel::detectCores()),
        bsTooltip("use_parallel",
                  "Enable parallel computation using future::multisession. Recommended for large datasets.",
                  "right", options = list(container = "body")),
        bsTooltip("cores",
                  "Number of CPU cores to use for parallel computation.",
                  "right", options = list(container = "body"))
      ),

      hr(),

      h4("Downloads"),
      downloadButton("download_profile", "Download likelihood profiles"),
      downloadButton("download_params", "Download fitted parameters"),

      hr(),
      h4("Report"),
      checkboxInput("include_dataset", "Include dataset in report", value = FALSE),
      checkboxGroupInput("report_elements", "Include in report:",
                         choices = list(
                           "Fitted parameters"   = "params",
                           "Likelihood profile"  = "profile",
                           "Stochastic curve"    = "curve",
                           "Histograms"          = "hist",
                           "Kernel density"      = "kernel",
                           "Boxplots"            = "box",
                           "Violins"             = "violin",
                           "Fan plot"            = "fan"
                         ),
                         selected = c("params", "profile", "curve")),
      downloadButton("download_report_html", "Generate HTML report")

    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Fitted parameters",
                 tableOutput("params_table")),

        tabPanel("Likelihood profile",
                 plotlyOutput("profile_plot"),
                 br(),
                 downloadButton("download_profile_jpg", "Download JPG"),
                 downloadButton("download_profile_pdf", "Download PDF")),

        tabPanel("Stochastic curve",
                 plotlyOutput("curve_plot"),
                 br(),
                 downloadButton("download_curve_jpg", "Download JPG"),
                 downloadButton("download_curve_pdf", "Download PDF")),

        tabPanel("Histograms",
                 plotlyOutput("hist_plot"),
                 br(),
                 downloadButton("download_hist_jpg", "Download JPG"),
                 downloadButton("download_hist_pdf", "Download PDF")),

        tabPanel("Kernel density",
                 plotlyOutput("kernel_plot"),
                 br(),
                 downloadButton("download_kernel_jpg", "Download JPG"),
                 downloadButton("download_kernel_pdf", "Download PDF")),

        tabPanel("Boxplots",
                 plotlyOutput("box_plot"),
                 br(),
                 downloadButton("download_box_jpg", "Download JPG"),
                 downloadButton("download_box_pdf", "Download PDF")),

        tabPanel("Violins",
                 plotlyOutput("violin_plot"),
                 br(),
                 downloadButton("download_violin_jpg", "Download JPG"),
                 downloadButton("download_violin_pdf", "Download PDF")),

        tabPanel("Fan plot",
                 plotlyOutput("fan_plot"),
                 br(),
                 downloadButton("download_fan_jpg", "Download JPG"),
                 downloadButton("download_fan_pdf", "Download PDF")),

        tabPanel("Model description",
                 h4("Okuyama mechanistic model"),
                 p("This app implements a mechanistic, stochastic functional response model inspired by Okuyama."),
                 p("The model describes a parasitoid (or predator) that repeatedly searches for hosts, encounters them, and handles them over a fixed experimental duration."),
                 tags$ul(
                   tags$li("Search time follows a Gamma distribution, controlled by parameter k and encounter rate lambda."),
                   tags$li("Encounter rate scales with host density x as lambda = a * x^z, where a is the baseline search efficiency and z controls how search scales with density."),
                   tags$li("Handling time follows a Lognormal distribution with mean h and variability s, or is constant when s = 0."),
                   tags$li("The parasitoid alternates between search and handling until the total time T is exhausted."),
                   tags$li("The number of hosts attacked at least once is recorded as the functional response.")
                 ),
                 h5("Parameters"),
                 tags$ul(
                   tags$li("a: baseline search rate; higher values mean more frequent encounters."),
                   tags$li("h: mean handling time per host; higher values mean slower processing and stronger saturation."),
                   tags$li("z: density-scaling exponent; z = 1 gives Type II-like behavior, z > 1 can produce Type III-like responses, z < 1 yields sublinear responses."),
                   tags$li("k: shape parameter of the Gamma search-time distribution; controls variability in search intervals."),
                   tags$li("s: standard deviation (log-scale) of handling time; controls variability in handling.")
                 ),
                 h5("Emergent functional responses"),
                 p("Because the model is generative, the functional response curve is not imposed by a fixed formula. It emerges from the simulated search-encounter-handling process."),
                 tags$ul(
                   tags$li("Type I, II, and III-like responses can all emerge depending on parameter values."),
                   tags$li("Sublinear, superlinear, saturating, and sigmoidal shapes are possible."),
                   tags$li("The model is suitable for both parasitoids and predators with discrete attack events.")
                 )
        ),

        tabPanel("Help / Workflow guide",
                 h4("Typical workflow"),
                 tags$ol(
                   tags$li("Upload a CSV dataset containing at least: species, host density, and parasitism columns."),
                   tags$li("Select the appropriate columns for species, host density, and parasitism in the sidebar."),
                   tags$li("Choose how many species to analyze and set the experiment duration (hours)."),
                   tags$li("Optionally open Advanced settings to tune DEoptim and stochastic simulation parameters."),
                   tags$li("Click 'Run analysis' to fit the mechanistic model and generate all plots."),
                   tags$li("Inspect fitted parameters, likelihood profiles, and stochastic curves in the main tabs."),
                   tags$li("Use the Downloads section to export fitted parameters, likelihood profiles, and reports.")
                 ),
                 h5("Recommended settings"),
                 tags$ul(
                   tags$li("Start with default Advanced settings (itermax = 50, NP = 40, reltol = 1e-2, n_sim = 3000)."),
                   tags$li("Increase itermax and NP if convergence seems poor or profiles are noisy."),
                   tags$li("Increase n_sim for smoother stochastic curves and distributions, at the cost of runtime."),
                   tags$li("Use parallel computation with multiple cores for large datasets or many species.")
                 ),
                 h5("Interpreting results"),
                 tags$ul(
                   tags$li("Fitted parameters summarize search efficiency, handling time, density scaling, and variability."),
                   tags$li("Likelihood profiles for z help assess whether the data support Type II-like (z ~ 1), Type III-like (z > 1), or sublinear (z < 1) behavior."),
                   tags$li("Stochastic curves show the mean mechanistic response with uncertainty bands."),
                   tags$li("Histograms, kernel densities, boxplots, violins, and fan plots reveal the full distribution of simulated parasitism across densities.")
                 )
        ),

        tabPanel("Report generator",
                 h4("Configure report contents"),
                 p("Use the checkboxes in the sidebar to select which elements to include in the report."),
                 p("Generate an HTML report using the button in the sidebar."),
                 p("The report will summarize the analysis for the current dataset and selected species, including fitted parameters, profiles, curves, and distributions, depending on your choices.")
        )
      )
    )
  )
)
