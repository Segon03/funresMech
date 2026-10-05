#' Server Logic for funresMech App
#'
#' @importFrom shiny reactiveValues observeEvent renderTable renderUI req
#'   showNotification withProgress incProgress downloadHandler renderPlot
#'   tagList selectInput validate need
#' @importFrom plotly renderPlotly ggplotly
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_point
#'   geom_histogram geom_density geom_boxplot geom_violin geom_smooth
#'   geom_jitter scale_x_continuous facet_grid vars theme_bw labs ggsave
#' @importFrom dplyr filter rename
#' @importFrom future plan multisession sequential
#' @importFrom parallel detectCores
#' @importFrom rmarkdown render
#' @importFrom stats quantile setNames
#' @importFrom magrittr %>%
#' @importFrom utils read.csv write.csv
#' @importFrom rlang sym .data
#' @noRd
NULL


# ---------------------------------------------------------------------------
# Helpers internos
# ---------------------------------------------------------------------------

#' Niveles de especie derivados de los datos
#' @noRd
.sp_levels <- function(x) {
  sort(unique(as.character(x)))
}

#' Desplazamientos horizontales genericos para separar especies en el eje x
#'
#' Reemplaza los offsets fijos c(Gp = -0.22, As = 0, Bs = 0.22), que solo
#' funcionaban con un conjunto de especies concreto.
#' @noRd
.sp_offsets <- function(sp_levels, span = 0.44) {
  n <- length(sp_levels)
  if (n <= 1L) {
    return(setNames(rep(0, max(n, 1L)), sp_levels))
  }
  setNames(seq(-span / 2, span / 2, length.out = n), sp_levels)
}

#' Niveles de densidad ordenados numericamente cuando es posible
#' @noRd
.dens_levels <- function(x) {
  num <- suppressWarnings(as.numeric(as.character(x)))
  if (anyNA(num)) sort(unique(as.character(x))) else sort(unique(num))
}


# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  rv <- reactiveValues(
    dat = NULL,
    params_all = NULL,
    profiles_all = NULL,
    curve_all = NULL,
    hist_all = NULL,
    aic_table = NULL,
    z_ci = NULL,
    diag = NULL,
    screen = NULL
  )

  output$column_select <- renderUI({
    req(input$dataset)
    dat <- read.csv(input$dataset$datapath)
    validate(need(ncol(dat) >= 3,
                  "El archivo debe tener al menos tres columnas."))
    tagList(
      selectInput("species_col", "Species column", choices = colnames(dat)),
      selectInput("dens_col",    "Density column", choices = colnames(dat)),
      selectInput("par_col",     "Parasitism column", choices = colnames(dat))
    )
  })

  observeEvent(input$run, {

    req(input$dataset, input$species_col, input$dens_col, input$par_col)

    dat <- read.csv(input$dataset$datapath)
    rv$dat <- dat

    species_list <- unique(dat[[input$species_col]])

    # Evita NAs si n_species supera el numero de especies disponibles
    n_sel <- min(as.integer(input$n_species), length(species_list))
    if (n_sel < 1L) {
      showNotification("No species found in the selected column.",
                       type = "error")
      return(invisible(NULL))
    }
    selected_species <- species_list[seq_len(n_sel)]

    all_params     <- list()
    all_profiles   <- list()
    all_curve_data <- list()
    all_hist_data  <- list()
    all_aic        <- list()
    all_z_ci       <- list()
    all_diag       <- list()
    all_screen     <- list()

    # Restaurar el plan de 'future' del usuario al salir, pase lo que pase
    old_plan <- plan()
    on.exit(plan(old_plan), add = TRUE)

    if (isTRUE(input$use_parallel)) {
      n_cores <- min(as.integer(input$cores),
                     max(1L, detectCores() - 1L))
      plan(multisession, workers = n_cores)
    } else {
      plan(sequential)
    }

    # An error in the analysis must not close the app: it is reported to the
    # user and the previous results are kept.
    ok <- tryCatch({
    withProgress(message = "Running mechanistic analysis...", value = 0, {

      for (sp in selected_species) {

        incProgress(1 / length(selected_species), detail = paste("Species:", sp))

        dat_sp <- dat %>%
          filter(.data[[input$species_col]] == sp) %>%
          rename(
            dens = !!sym(input$dens_col),
            par  = !!sym(input$par_col)
          )

        chk <- tryCatch(.check_data_spp(dat_sp), error = function(e) e)
        if (inherits(chk, "error")) {
          showNotification(paste0("Species ", sp, ": ", conditionMessage(chk)),
                           type = "error")
          next
        }

        # Screening of atypical trials BEFORE the fit (never removes data)
        if (isTRUE(input$run_screening)) {
          scr <- tryCatch(screen_outliers(dat_sp, input$T_exp),
                          error = function(e) NULL)
          if (!is.null(scr)) {
            scr$species <- sp
            all_screen[[as.character(sp)]] <- scr
            n_flag <- sum(scr$flagged)
            if (n_flag > 0) {
              showNotification(
                paste0("Species ", sp, ": ", n_flag, " atypical trial(s) ",
                       "flagged (see 'Data screening'). A single atypical ",
                       "trial can leave z without an upper limit; consider a ",
                       "sensitivity analysis without it."),
                type = "warning", duration = 15)
            }
          }
        }

        z_grid <- seq(input$z_min, input$z_max, by = input$z_step)

        fit <- fit_profile(
          data_spp = dat_sp,
          T_exp    = input$T_exp,
          z_grid   = z_grid,
          n_sim    = input$n_sim,
          itermax  = input$itermax,
          NP       = input$NP,
          reltol   = input$reltol,
          extend_z = isTRUE(input$extend_z)
        )

        a <- unname(fit$par["a"])
        h <- unname(fit$par["h"])
        z <- unname(fit$par["z"])
        k <- unname(fit$par["k"])
        s <- unname(fit$par["s"])

        profile <- fit$profile
        profile$species <- sp
        all_profiles[[as.character(sp)]] <- profile

        all_aic[[as.character(sp)]] <- data.frame(
          species          = sp,
          AIC_full         = fit$aic$AIC_full,
          AIC_restricted   = fit$aic$AIC_restricted,
          delta_AIC        = fit$aic$delta_AIC,
          stringsAsFactors = FALSE
        )

        ci <- fit$ci
        ci_note <- c(
          if (ci$low_censored)  paste0("lower limit <= ", min(profile$z)),
          if (ci$high_censored) .high_open_text(max(profile$z))
        )
        all_z_ci[[as.character(sp)]] <- data.frame(
          species    = sp,
          z_estimate = z,
          z_low      = ci$z_low,
          z_high     = ci$z_high,
          ci_note    = if (length(ci_note)) paste(ci_note, collapse = "; ") else "",
          stringsAsFactors = FALSE
        )

        all_diag[[as.character(sp)]] <- data.frame(
          species = sp,
          k_hat   = k,
          notes   = if (length(fit$notes)) paste(fit$notes, collapse = " | ") else "No warnings",
          stringsAsFactors = FALSE
        )

        params_df <- as.data.frame(t(fit$par))
        params_df$species <- sp
        params_df$nll     <- fit$nll
        all_params[[as.character(sp)]] <- params_df

        densities  <- sort(unique(dat_sp$dens))
        curve_data <- data.frame()
        hist_data  <- data.frame()

        for (x in densities) {
          sims <- motor_okuyama_cpp(as.integer(x), a, h, z, k, s,
                                    input$T_exp, as.integer(input$n_sim))

          curve_data <- rbind(
            curve_data,
            data.frame(
              dens    = x,
              mean    = mean(sims),
              lower   = unname(quantile(sims, 0.025)),
              upper   = unname(quantile(sims, 0.975)),
              species = sp,
              stringsAsFactors = FALSE
            )
          )

          hist_data <- rbind(
            hist_data,
            data.frame(
              dens    = x,
              paras   = sims,
              species = sp,
              stringsAsFactors = FALSE
            )
          )
        }

        all_curve_data[[as.character(sp)]] <- curve_data
        all_hist_data[[as.character(sp)]]  <- hist_data
      }
    })
    TRUE
    }, error = function(e) {
      showNotification(paste0("The analysis failed: ", conditionMessage(e)),
                       type = "error", duration = NULL)
      FALSE
    })
    if (!isTRUE(ok)) return(invisible(NULL))

    rv$params_all   <- do.call(rbind, all_params)
    rv$profiles_all <- do.call(rbind, all_profiles)
    rv$curve_all    <- do.call(rbind, all_curve_data)
    rv$hist_all     <- do.call(rbind, all_hist_data)
    rv$aic_table    <- do.call(rbind, all_aic)
    rv$z_ci         <- do.call(rbind, all_z_ci)
    rv$diag         <- do.call(rbind, all_diag)
    rv$screen       <- do.call(rbind, all_screen)
  })


  # -------------------------------------------------------------------------
  # Graficos
  # -------------------------------------------------------------------------

  plot_profile_static <- function() {
    req(rv$profiles_all)
    ggplot(rv$profiles_all,
           aes(x = .data$z, y = .data$nll, color = .data$species)) +
      geom_line(linewidth = 1) +
      theme_bw() +
      labs(x = "z", y = "Negative log-likelihood",
           title = "Likelihood profile for z")
  }

  plot_profile_plotly <- function() plot_profile_static()

  plot_curve_static <- function() {
    req(rv$curve_all, rv$dat)
    dat_long <- rv$dat %>%
      filter(.data[[input$species_col]] %in% unique(rv$curve_all$species)) %>%
      rename(
        dens    = !!sym(input$dens_col),
        par     = !!sym(input$par_col),
        species = !!sym(input$species_col)
      )

    ggplot(rv$curve_all,
           aes(x = .data$dens, y = .data$mean, color = .data$species)) +
      geom_line(linewidth = 1.2) +
      geom_ribbon(aes(ymin = .data$lower, ymax = .data$upper,
                      fill = .data$species),
                  alpha = 0.2, color = NA) +
      geom_point(data = dat_long,
                 aes(x = .data$dens, y = .data$par),
                 color = "black", size = 2, inherit.aes = FALSE) +
      theme_bw() +
      labs(x = "Host density", y = "Parasitized hosts",
           title = "Stochastic mechanistic curve")
  }

  plot_curve_plotly <- function() plot_curve_static()

  plot_hist_static <- function() {
    req(rv$hist_all)
    ggplot(rv$hist_all, aes(x = .data$paras)) +
      geom_histogram(binwidth = 1, fill = "skyblue", color = "white") +
      facet_grid(rows = vars(.data$species), cols = vars(.data$dens)) +
      theme_bw() +
      labs(x = "Parasitized hosts", y = "Count",
           title = "Histograms of simulated parasitism")
  }

  plot_hist_plotly <- function() plot_hist_static()

  plot_kernel_static <- function() {
    req(rv$hist_all)
    ggplot(rv$hist_all,
           aes(x = .data$paras, color = .data$species, fill = .data$species)) +
      geom_density(alpha = 0.4) +
      facet_grid(rows = vars(.data$species), cols = vars(.data$dens)) +
      theme_bw() +
      labs(x = "Parasitized hosts", y = "Density",
           title = "Kernel density of simulated parasitism")
  }

  plot_kernel_plotly <- function() plot_kernel_static()

  plot_box_static <- function() {
    req(rv$hist_all)
    df <- rv$hist_all
    df$dens    <- factor(df$dens,    levels = .dens_levels(df$dens))
    df$species <- factor(df$species, levels = .sp_levels(df$species))
    df$dens_species <- interaction(df$dens, df$species, sep = ".", drop = TRUE)

    ggplot(df, aes(x = .data$dens_species, y = .data$paras,
                   fill = .data$species)) +
      geom_boxplot() +
      theme_bw() +
      labs(title = "Box plots of simulated parasitism",
           x = "Density x Species", y = "Parasitized hosts")
  }

  plot_box_plotly <- function() plot_box_static()

  plot_violin_static <- function() {
    req(rv$hist_all)
    df <- rv$hist_all
    df$dens    <- factor(df$dens,    levels = .dens_levels(df$dens))
    df$species <- factor(df$species, levels = .sp_levels(df$species))
    df$dens_species <- interaction(df$dens, df$species, sep = "-", drop = TRUE)

    ggplot(df, aes(x = .data$dens_species, y = .data$paras,
                   fill = .data$species)) +
      geom_violin(alpha = 0.6, trim = FALSE) +
      geom_jitter(aes(color = .data$species),
                  width = 0.15, height = 0, alpha = 0.4, size = 1.5) +
      theme_bw() +
      labs(title = "Violin plots of simulated parasitism",
           x = "Density x Species", y = "Parasitized hosts")
  }

  plot_violin_plotly <- function() {
    req(rv$hist_all)
    df <- rv$hist_all

    dens_num <- suppressWarnings(as.numeric(as.character(df$dens)))
    if (anyNA(dens_num)) return(plot_violin_static())

    sp_levels <- .sp_levels(df$species)
    offsets   <- .sp_offsets(sp_levels)

    df$dens        <- dens_num
    df$species     <- factor(df$species, levels = sp_levels)
    df$dens_offset <- df$dens + offsets[as.character(df$species)]
    df$group_id    <- interaction(df$dens, df$species, drop = TRUE)

    ggplot(df, aes(x = .data$dens_offset, y = .data$paras,
                   fill = .data$species, group = .data$group_id)) +
      geom_violin(alpha = 0.6, trim = FALSE, width = 0.35) +
      geom_jitter(aes(color = .data$species),
                  width = 0.05, height = 0, alpha = 0.5, size = 1.5,
                  show.legend = FALSE) +
      scale_x_continuous(breaks = sort(unique(df$dens)),
                         labels = sort(unique(df$dens))) +
      theme_bw() +
      labs(title = "Violin plots of simulated parasitism",
           x = "Host density", y = "Parasitized hosts")
  }

  plot_fan_static <- function() {
    req(rv$hist_all, rv$dat)

    df <- rv$hist_all
    sp_levels <- .sp_levels(df$species)
    offsets   <- .sp_offsets(sp_levels, span = 0.4)

    df$species     <- factor(df$species, levels = sp_levels)
    df$dens_num    <- suppressWarnings(as.numeric(as.character(df$dens)))
    df$dens_offset <- df$dens_num + offsets[as.character(df$species)]

    dat_long <- rv$dat %>%
      filter(.data[[input$species_col]] %in% sp_levels) %>%
      rename(
        dens    = !!sym(input$dens_col),
        par     = !!sym(input$par_col),
        species = !!sym(input$species_col)
      )

    dat_long$species     <- factor(dat_long$species, levels = sp_levels)
    dat_long$dens_num    <- suppressWarnings(as.numeric(as.character(dat_long$dens)))
    dat_long$dens_offset <- dat_long$dens_num +
      offsets[as.character(dat_long$species)]

    ggplot() +
      geom_point(data = df,
                 aes(x = .data$dens_offset, y = .data$paras,
                     color = .data$species),
                 alpha = 0.3, size = 1.5) +
      geom_smooth(data = df,
                  aes(x = .data$dens_num, y = .data$paras,
                      color = .data$species),
                  method = "loess", formula = y ~ x, se = FALSE,
                  linewidth = 1.2) +
      geom_point(data = dat_long,
                 aes(x = .data$dens_offset, y = .data$par,
                     color = .data$species),
                 alpha = 0.8, size = 2.5) +
      theme_bw() +
      labs(x = "Host density", y = "Parasitized hosts",
           title = "Fan plot of simulated parasitism") +
      scale_x_continuous(breaks = sort(unique(df$dens_num)))
  }

  plot_fan_plotly <- function() {
    req(rv$hist_all, rv$dat)
    dat_long <- rv$dat %>%
      filter(.data[[input$species_col]] %in% unique(rv$hist_all$species)) %>%
      rename(
        dens    = !!sym(input$dens_col),
        par     = !!sym(input$par_col),
        species = !!sym(input$species_col)
      )

    ggplot(rv$hist_all,
           aes(x = .data$dens, y = .data$paras, color = .data$species)) +
      geom_point(alpha = 0.2) +
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                  linewidth = 1) +
      geom_point(data = dat_long,
                 aes(x = .data$dens, y = .data$par, color = .data$species),
                 alpha = 0.8, size = 2.5) +
      theme_bw() +
      labs(x = "Host density", y = "Parasitized hosts",
           title = "Fan plot of simulated parasitism")
  }


  # -------------------------------------------------------------------------
  # Outputs
  # -------------------------------------------------------------------------

  output$params_table <- renderTable({
    req(rv$params_all)
    rv$params_all
  }, rownames = FALSE)

  output$diag_table <- renderTable({
    req(rv$diag)
    rv$diag
  }, rownames = FALSE)

  output$screen_table <- renderTable({
    validate(need(!is.null(rv$screen),
                  "Screening was not run (enable it in the sidebar) or found nothing to report."))
    scr <- rv$screen
    scr <- scr[scr$flagged | scr$influences_dispersion, , drop = FALSE]
    validate(need(nrow(scr) > 0, "No atypical or influential trials were flagged."))
    scr[, c("species", "row", "dens", "par", "expected_loo", "p_loo",
            "rho_change", "flagged")]
  }, rownames = FALSE)

  output$profile_plot <- renderPlotly({ ggplotly(plot_profile_plotly()) })
  output$curve_plot   <- renderPlotly({ ggplotly(plot_curve_plotly())   })
  output$hist_plot    <- renderPlotly({ ggplotly(plot_hist_plotly())    })
  output$kernel_plot  <- renderPlotly({ ggplotly(plot_kernel_plotly())  })
  output$box_plot     <- renderPlotly({ ggplotly(plot_box_plotly())     })
  output$violin_plot  <- renderPlotly({ ggplotly(plot_violin_plotly())  })
  output$fan_plot     <- renderPlotly({ ggplotly(plot_fan_plotly())     })


  # -------------------------------------------------------------------------
  # Descargas de datos (CSV)
  #
  # Botones del sidebar (no confundir con las descargas JPG/PDF de cada
  # pestaña, que exportan el grafico). 'file' es una ruta temporal
  # administrada por shiny.
  # -------------------------------------------------------------------------

  output$download_profile <- downloadHandler(
    filename = "likelihood_profiles.csv",
    content = function(file) {
      req(rv$profiles_all)
      write.csv(rv$profiles_all, file, row.names = FALSE)
    }
  )

  output$download_params <- downloadHandler(
    filename = "fitted_parameters.csv",
    content = function(file) {
      req(rv$params_all)
      write.csv(rv$params_all, file, row.names = FALSE)
    }
  )


  # -------------------------------------------------------------------------
  # Descargas de graficos
  #
  # 'file' es siempre una ruta temporal administrada por shiny, de modo que
  # no se escribe en el espacio del usuario ni en el directorio del paquete.
  # -------------------------------------------------------------------------

  output$download_profile_jpg <- downloadHandler(
    filename = "likelihood_profile.jpg",
    content = function(file) {
      ggsave(file, plot = plot_profile_static(), device = "jpeg",
             width = 7, height = 5, dpi = 300)
    }
  )

  output$download_profile_pdf <- downloadHandler(
    filename = "likelihood_profile.pdf",
    content = function(file) {
      ggsave(file, plot = plot_profile_static(), device = "pdf",
             width = 7, height = 5)
    }
  )

  output$download_curve_jpg <- downloadHandler(
    filename = "stochastic_curve.jpg",
    content = function(file) {
      ggsave(file, plot = plot_curve_static(), device = "jpeg",
             width = 7, height = 5, dpi = 300)
    }
  )

  output$download_curve_pdf <- downloadHandler(
    filename = "stochastic_curve.pdf",
    content = function(file) {
      ggsave(file, plot = plot_curve_static(), device = "pdf",
             width = 7, height = 5)
    }
  )

  output$download_hist_jpg <- downloadHandler(
    filename = "histograms.jpg",
    content = function(file) {
      ggsave(file, plot = plot_hist_static(), device = "jpeg",
             width = 7, height = 5, dpi = 300)
    }
  )

  output$download_hist_pdf <- downloadHandler(
    filename = "histograms.pdf",
    content = function(file) {
      ggsave(file, plot = plot_hist_static(), device = "pdf",
             width = 7, height = 5)
    }
  )

  output$download_kernel_jpg <- downloadHandler(
    filename = "kernel_density.jpg",
    content = function(file) {
      ggsave(file, plot = plot_kernel_static(), device = "jpeg",
             width = 7, height = 5, dpi = 300)
    }
  )

  output$download_kernel_pdf <- downloadHandler(
    filename = "kernel_density.pdf",
    content = function(file) {
      ggsave(file, plot = plot_kernel_static(), device = "pdf",
             width = 7, height = 5)
    }
  )

  output$download_box_jpg <- downloadHandler(
    filename = "boxplots.jpg",
    content = function(file) {
      ggsave(file, plot = plot_box_static(), device = "jpeg",
             width = 7, height = 5, dpi = 300)
    }
  )

  output$download_box_pdf <- downloadHandler(
    filename = "boxplots.pdf",
    content = function(file) {
      ggsave(file, plot = plot_box_static(), device = "pdf",
             width = 7, height = 5)
    }
  )

  output$download_violin_jpg <- downloadHandler(
    filename = "violins.jpg",
    content = function(file) {
      ggsave(file, plot = plot_violin_static(), device = "jpeg",
             width = 7, height = 5, dpi = 300)
    }
  )

  output$download_violin_pdf <- downloadHandler(
    filename = "violins.pdf",
    content = function(file) {
      ggsave(file, plot = plot_violin_static(), device = "pdf",
             width = 7, height = 5)
    }
  )

  output$download_fan_jpg <- downloadHandler(
    filename = "fan_plot.jpg",
    content = function(file) {
      ggsave(file, plot = plot_fan_static(), device = "jpeg",
             width = 7, height = 5, dpi = 300)
    }
  )

  output$download_fan_pdf <- downloadHandler(
    filename = "fan_plot.pdf",
    content = function(file) {
      ggsave(file, plot = plot_fan_static(), device = "pdf",
             width = 7, height = 5)
    }
  )


  # -------------------------------------------------------------------------
  # Descarga del reporte
  #
  # La plantilla se copia a un directorio temporal antes de renderizar:
  # rmarkdown::render() escribe sus archivos intermedios junto al archivo de
  # entrada, que de otro modo seria el directorio del paquete instalado.
  # -------------------------------------------------------------------------

  output$download_report_html <- downloadHandler(
    filename = function() "Okuyama_report.html",
    content = function(file) {

      if (!requireNamespace("kableExtra", quietly = TRUE)) {
        showNotification(
          "El paquete 'kableExtra' es necesario para generar el reporte.",
          type = "error"
        )
        return(invisible(NULL))
      }

      rmd_src <- system.file("app", "report_template.Rmd",
                             package = "funresMech")
      if (!nzchar(rmd_src)) {
        stop("No se encontro 'report_template.Rmd' en el paquete instalado.")
      }

      tmp_dir <- tempfile("funresMech_report_")
      dir.create(tmp_dir)
      on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

      rmd_tmp <- file.path(tmp_dir, "report_template.Rmd")
      file.copy(rmd_src, rmd_tmp, overwrite = TRUE)

      dataset_renamed <- rv$dat %>%
        rename(
          species = !!sym(input$species_col),
          dens    = !!sym(input$dens_col),
          par     = !!sym(input$par_col)
        )

      params_list <- list(
        DATASET         = dataset_renamed,
        PARAMS_ALL      = rv$params_all,
        PROFILES_ALL    = rv$profiles_all,
        CURVE_ALL       = rv$curve_all,
        HIST_ALL        = rv$hist_all,
        AIC_TABLE       = rv$aic_table,
        Z_CI            = rv$z_ci,
        DIAG            = rv$diag,
        SCREEN          = rv$screen,
        REPORT_ELEMENTS = input$report_elements,
        INCLUDE_DATASET = input$include_dataset,
        T_EXP           = input$T_exp,
        ITERMAX         = input$itermax,
        NP              = input$NP,
        RELTOL          = input$reltol,
        N_SIM           = input$n_sim,
        Z_MIN           = input$z_min,
        Z_MAX           = input$z_max,
        Z_STEP          = input$z_step,
        CORES           = input$cores
      )

      render(
        input             = rmd_tmp,
        output_file       = file,
        intermediates_dir = tmp_dir,
        knit_root_dir     = tmp_dir,
        params            = params_list,
        envir             = new.env(parent = globalenv()),
        quiet             = TRUE
      )
    }
  )

}  # end server
