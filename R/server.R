#' Server Logic for funresMech App
#'
#' @importFrom shiny reactiveValues observeEvent renderTable renderUI req showNotification withProgress incProgress downloadHandler renderPlot
#' @importFrom plotly renderPlotly ggplotly
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_point geom_histogram geom_density geom_boxplot geom_violin geom_smooth geom_jitter scale_x_continuous facet_grid theme_bw labs ggsave position_dodge position_jitterdodge
#' @importFrom dplyr filter rename group_by summarise mutate
#' @importFrom tidyr pivot_longer
#' @importFrom DEoptim DEoptim DEoptim.control
#' @importFrom future plan multisession sequential
#' @importFrom rmarkdown render
#' @importFrom stats quantile qchisq rgamma rlnorm
#' @importFrom magrittr %>%
#' @importFrom utils read.csv write.csv
#' @importFrom rlang sym
#' @noRd
NULL

server <- function(input, output, session) {

  rv <- reactiveValues(
    dat = NULL,
    params_all = NULL,
    profiles_all = NULL,
    curve_all = NULL,
    hist_all = NULL,
    aic_table = NULL,
    z_ci = NULL
  )

  output$column_select <- renderUI({
    req(input$dataset)
    dat <- read.csv(input$dataset$datapath)
    tagList(
      selectInput("species_col", "Species column", choices = colnames(dat)),
      selectInput("dens_col", "Density column", choices = colnames(dat)),
      selectInput("par_col", "Parasitism column", choices = colnames(dat))
    )
  })

  observeEvent(input$run, {

    req(input$dataset, input$species_col, input$dens_col, input$par_col)

    dat <- read.csv(input$dataset$datapath)
    rv$dat <- dat

    species_list <- unique(dat[[input$species_col]])
    selected_species <- species_list[1:input$n_species]

    all_params <- list()
    all_profiles <- list()
    all_curve_data <- list()
    all_hist_data <- list()
    all_aic <- list()
    all_z_ci <- list()

    if (isTRUE(input$use_parallel)) {
      plan(multisession, workers = input$cores)
    } else {
      plan(sequential)
    }

    withProgress(message = "Running mechanistic analysis...", value = 0, {

      for (sp in selected_species) {

        incProgress(1 / length(selected_species), detail = paste("Species:", sp))

        dat_sp <- dat %>%
          filter(.data[[input$species_col]] == sp) %>%
          rename(
            dens = !!input$dens_col,
            par  = !!input$par_col
          )

        full_fit <- fit_full(
          data_spp       = dat_sp,
          T_exp          = input$T_exp,
          itermax        = input$itermax,
          NP             = input$NP,
          reltol         = input$reltol,
          n_sim_profile  = input$n_sim
        )

        a <- full_fit$par["a"]
        h <- full_fit$par["h"]
        z <- full_fit$par["z"]
        k <- full_fit$par["k"]
        s <- full_fit$par["s"]
        nll_full <- full_fit$nll

        z_grid <- seq(input$z_min, input$z_max, by = input$z_step)
        profile <- data.frame(z = z_grid, nll = NA)

        for (i in seq_along(z_grid)) {
          z_val <- z_grid[i]

          lower <- c(a = 0.001, h = 0.001, k = 0.5, s = 0.001)
          upper <- c(a = 2.0,   h = 0.5,   k = 5.0, s = 0.5)

          res <- DEoptim(
            fn = function(par) {
              negloglik_fixed_z(
                par_vec = par,
                z_fixed = z_val,
                data_spp = dat_sp,
                T = input$T_exp,
                n_sim = input$n_sim
              )
            },
            lower = lower,
            upper = upper,
            control = DEoptim.control(
              itermax = input$itermax,
              NP = input$NP,
              reltol = input$reltol,
              trace = FALSE
            )
          )

          profile$nll[i] <- res$optim$bestval
        }

        profile$species <- sp
        all_profiles[[sp]] <- profile

        AIC_full <- 2 * nll_full + 2 * 5
        z_vals <- profile$z
        nll_vals <- profile$nll

        if (min(z_vals) <= 1 && max(z_vals) >= 1) {
          idx <- which.min(abs(z_vals - 1))
          if (abs(z_vals[idx] - 1) < 0.01) {
            nll_1 <- nll_vals[idx]
          } else {
            idx_low <- max(which(z_vals <= 1))
            idx_high <- min(which(z_vals >= 1))
            if (idx_low > 0 && idx_high <= length(z_vals)) {
              z_low <- z_vals[idx_low]
              nll_low <- nll_vals[idx_low]
              z_high <- z_vals[idx_high]
              nll_high <- nll_vals[idx_high]
              nll_1 <- nll_low + (nll_high - nll_low) * (1 - z_low) / (z_high - z_low)
            } else {
              nll_1 <- NA
            }
          }
        } else {
          nll_1 <- NA
        }

        if (!is.na(nll_1)) {
          AIC_restricted <- 2 * nll_1 + 2 * 4
          delta_AIC <- AIC_restricted - AIC_full
        } else {
          AIC_restricted <- NA
          delta_AIC <- NA
        }

        aic_row <- data.frame(
          species = sp,
          AIC_full = AIC_full,
          AIC_restricted = AIC_restricted,
          delta_AIC = delta_AIC,
          stringsAsFactors = FALSE
        )
        all_aic[[sp]] <- aic_row

        min_nll <- min(profile$nll)
        threshold <- min_nll + qchisq(0.95, 1)
        idx_min <- which.min(nll_vals)

        z_low <- NA
        for (i in seq(idx_min, 1, by = -1)) {
          if (nll_vals[i] >= threshold) {
            if (i < length(z_vals)) {
              z1 <- z_vals[i]
              n1 <- nll_vals[i]
              z2 <- z_vals[i + 1]
              n2 <- nll_vals[i + 1]
              z_low <- z1 + (z2 - z1) * (threshold - n1) / (n2 - n1)
            } else {
              z_low <- z_vals[i]
            }
            break
          }
        }

        z_high <- NA
        for (i in seq(idx_min, length(z_vals), by = 1)) {
          if (nll_vals[i] >= threshold) {
            if (i > 1) {
              z1 <- z_vals[i - 1]
              n1 <- nll_vals[i - 1]
              z2 <- z_vals[i]
              n2 <- nll_vals[i]
              z_high <- z1 + (z2 - z1) * (threshold - n1) / (n2 - n1)
            } else {
              z_high <- z_vals[i]
            }
            break
          }
        }

        z_ci_row <- data.frame(
          species = sp,
          z_estimate = z,
          z_low = z_low,
          z_high = z_high,
          stringsAsFactors = FALSE
        )
        all_z_ci[[sp]] <- z_ci_row

        params_df <- as.data.frame(t(full_fit$par))
        params_df$species <- sp
        params_df$nll <- nll_full
        all_params[[sp]] <- params_df

        densities <- sort(unique(dat_sp$dens))
        curve_data <- data.frame()
        hist_data  <- data.frame()

        simulate_okuyama <- function(x, T, a, h, z, k, s, n_sim = input$n_sim) {
          counts <- integer(n_sim)
          for (i in seq_len(n_sim)) {
            counts[i] <- simulate_trial(x, T, a, h, z, k, s)
          }
          counts
        }

        for (x in densities) {
          sims <- simulate_okuyama(x, input$T_exp, a, h, z, k, s)

          curve_data <- rbind(
            curve_data,
            data.frame(
              dens  = x,
              mean  = mean(sims),
              lower = quantile(sims, 0.025),
              upper = quantile(sims, 0.975),
              species = sp
            )
          )

          hist_data <- rbind(
            hist_data,
            data.frame(
              dens = x,
              paras = sims,
              species = sp
            )
          )
        }

        all_curve_data[[sp]] <- curve_data
        all_hist_data[[sp]]  <- hist_data
      }
    })

    rv$params_all   <- do.call(rbind, all_params)
    rv$profiles_all <- do.call(rbind, all_profiles)
    rv$curve_all    <- do.call(rbind, all_curve_data)
    rv$hist_all     <- do.call(rbind, all_hist_data)
    rv$aic_table    <- do.call(rbind, all_aic)
    rv$z_ci         <- do.call(rbind, all_z_ci)

    plan(sequential)
  })


  plot_profile_static <- function() {
    req(rv$profiles_all)
    ggplot(rv$profiles_all, aes(x = z, y = nll, color = species)) +
      geom_line(linewidth = 1) +
      theme_bw() +
      labs(x = "z", y = "Negative log-likelihood", title = "Likelihood profile for z")
  }

  plot_profile_plotly <- function() {
    plot_profile_static()
  }

  plot_curve_static <- function() {
    req(rv$curve_all, rv$dat)
    dat_long <- rv$dat %>%
      filter(.data[[input$species_col]] %in% unique(rv$curve_all$species)) %>%
      rename(
        dens = !!input$dens_col,
        par  = !!input$par_col,
        species = !!input$species_col
      )

    ggplot(rv$curve_all, aes(x = dens, y = mean, color = species)) +
      geom_line(linewidth = 1.2) +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = species), alpha = 0.2, color = NA) +
      geom_point(data = dat_long, aes(x = dens, y = par), color = "black", size = 2) +
      theme_bw() +
      labs(x = "Host density", y = "Parasitized hosts", title = "Stochastic mechanistic curve")
  }

  plot_curve_plotly <- function() {
    plot_curve_static()
  }

  plot_hist_static <- function() {
    req(rv$hist_all)
    ggplot(rv$hist_all, aes(x = paras)) +
      geom_histogram(binwidth = 1, fill = "skyblue", color = "white") +
      facet_grid(species ~ dens) +
      theme_bw() +
      labs(x = "Parasitized hosts", y = "Count", title = "Histograms of simulated parasitism")
  }

  plot_hist_plotly <- function() {
    plot_hist_static()
  }

  plot_kernel_static <- function() {
    req(rv$hist_all)
    ggplot(rv$hist_all, aes(x = paras, color = species, fill = species)) +
      geom_density(alpha = 0.4) +
      facet_grid(species ~ dens) +
      theme_bw() +
      labs(x = "Parasitized hosts", y = "Density", title = "Kernel density of simulated parasitism")
  }

  plot_kernel_plotly <- function() {
    plot_kernel_static()
  }

  plot_box_static <- function() {
    req(rv$hist_all)
    df <- rv$hist_all
    df$dens <- factor(df$dens, levels = c(1, 5, 10))
    df$species <- factor(df$species, levels = c("Gp", "As", "Bs"))
    df$dens_species <- interaction(df$dens, df$species, sep = ".")

    ggplot(df, aes(x = dens_species, y = paras, fill = species)) +
      geom_boxplot() +
      theme_bw() +
      labs(title = "Box plots of simulated parasitism", x = "Density x Species", y = "Parasitized hosts")
  }

  plot_box_plotly <- function() {
    plot_box_static()
  }

  plot_violin_static <- function() {
    req(rv$hist_all)
    df <- rv$hist_all
    df$dens <- factor(df$dens, levels = c(1, 5, 10))
    df$species <- factor(df$species, levels = c("Gp", "As", "Bs"))
    df$dens_species <- interaction(df$dens, df$species, sep = "-")

    ggplot(df, aes(x = dens_species, y = paras, fill = species)) +
      geom_violin(alpha = 0.6, trim = FALSE) +
      geom_jitter(aes(color = species), width = 0.15, height = 0, alpha = 0.4, size = 1.5) +
      theme_bw() +
      labs(title = "Violin plots of simulated parasitism", x = "Density x Species", y = "Parasitized hosts")
  }

  plot_violin_plotly <- function() {
    req(rv$hist_all)
    df <- rv$hist_all

    df$dens <- as.numeric(as.character(df$dens))
    df$species <- factor(df$species, levels = c("Gp", "As", "Bs"))

    offsets <- c(Gp = -0.22, As = 0.00, Bs = 0.22)
    df$dens_offset <- df$dens + offsets[as.character(df$species)]

    df$group_id <- interaction(df$dens, df$species, drop = TRUE)

    ggplot(df, aes(x = df$dens_offset, y = paras, fill = species, group = df$group_id)) +
      geom_violin(alpha = 0.6, trim = FALSE, width = 0.35) +
      geom_jitter(aes(x = df$dens_offset, color = species), width = 0.05, height = 0, alpha = 0.5, size = 1.5, show.legend = FALSE) +
      scale_x_continuous(breaks = sort(unique(df$dens)), labels = sort(unique(df$dens))) +
      theme_bw() +
      labs(title = "Violin plots of simulated parasitism", x = "Host density", y = "Parasitized hosts")
  }

  plot_fan_static <- function() {
    req(rv$hist_all, rv$dat)

    df <- rv$hist_all
    df$species <- factor(df$species, levels = c("Gp", "As", "Bs"))
    df$dens_offset <- df$dens +
      ifelse(df$species == "Gp", -0.2,
             ifelse(df$species == "As", 0.0,
                    ifelse(df$species == "Bs", 0.2, 0)))

    dat_long <- rv$dat %>%
      filter(.data[[input$species_col]] %in% unique(df$species)) %>%
      rename(
        dens = !!input$dens_col,
        par  = !!input$par_col,
        species = !!input$species_col
      )

    dat_long$species <- factor(dat_long$species, levels = c("Gp", "As", "Bs"))
    dat_long$dens_offset <- dat_long$dens +
      ifelse(dat_long$species == "Gp", -0.2,
             ifelse(dat_long$species == "As", 0.0,
                    ifelse(dat_long$species == "Bs", 0.2, 0)))

    ggplot() +
      geom_point(data = df, aes(x = dens_offset, y = paras, color = species), alpha = 0.3, size = 1.5) +
      geom_smooth(data = df, aes(x = dens, y = paras, color = species), method = "loess", se = FALSE, linewidth = 1.2) +
      geom_point(data = dat_long, aes(x = dens_offset, y = par, color = species), alpha = 0.8, size = 2.5) +
      theme_bw() +
      labs(x = "Host density", y = "Parasitized hosts", title = "Fan plot of simulated parasitism") +
      scale_x_continuous(breaks = sort(unique(df$dens)))
  }

  plot_fan_plotly <- function() {
    req(rv$hist_all, rv$dat)
    dat_long <- rv$dat %>%
      filter(.data[[input$species_col]] %in% unique(rv$hist_all$species)) %>%
      rename(
        dens = !!input$dens_col,
        par  = !!input$par_col,
        species = !!input$species_col
      )

    ggplot(rv$hist_all, aes(x = dens, y = paras, color = species)) +
      geom_point(alpha = 0.2) +
      geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
      geom_point(data = dat_long, aes(x = dens, y = par, color = species), alpha = 0.8, size = 2.5) +
      theme_bw() +
      labs(x = "Host density", y = "Parasitized hosts", title = "Fan plot of simulated parasitism")
  }

  output$params_table <- renderTable({
    req(rv$params_all)
    rv$params_all
  }, rownames = FALSE)

  output$profile_plot <- renderPlotly({
    ggplotly(plot_profile_plotly())
  })

  output$curve_plot <- renderPlotly({
    ggplotly(plot_curve_plotly())
  })

  output$hist_plot <- renderPlotly({
    ggplotly(plot_hist_plotly())
  })

  output$kernel_plot <- renderPlotly({
    ggplotly(plot_kernel_plotly())
  })

  output$box_plot <- renderPlotly({
    ggplotly(plot_box_plotly())
  })

  output$violin_plot <- renderPlotly({
    ggplotly(plot_violin_plotly())
  })

  output$fan_plot <- renderPlotly({
    ggplotly(plot_fan_plotly())
  })

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
             width = 7, height = 5, useDingbats = FALSE)
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
             width = 7, height = 5, useDingbats = FALSE)
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
             width = 7, height = 5, useDingbats = FALSE)
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
             width = 7, height = 5, useDingbats = FALSE)
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
             width = 7, height = 5, useDingbats = FALSE)
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
             width = 7, height = 5, useDingbats = FALSE)
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
             width = 7, height = 5, useDingbats = FALSE)
    }
  )

  output$download_report_html <- downloadHandler(
    filename = function() "Okuyama_report.html",
    content = function(file) {

      dataset_renamed <- rv$dat %>%
        rename(
          species = !!rlang::sym(input$species_col),
          dens    = !!rlang::sym(input$dens_col),
          par     = !!rlang::sym(input$par_col)
        )

      params_list <- list(
        DATASET         = dataset_renamed,
        PARAMS_ALL      = rv$params_all,
        PROFILES_ALL    = rv$profiles_all,
        CURVE_ALL       = rv$curve_all,
        HIST_ALL        = rv$hist_all,
        AIC_TABLE       = rv$aic_table,
        Z_CI            = rv$z_ci,
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

      rmd_path <- system.file("app/report_template.Rmd", package = "funresMech")
      if (rmd_path == "") {
        rmd_path <- "inst/app/report_template.Rmd"
      }

      rmarkdown::render(
        input       = rmd_path,
        output_file = file,
        params      = params_list,
        envir       = new.env()
      )
    }
  )

}  # end server
