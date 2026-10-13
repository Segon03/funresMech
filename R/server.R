############################################################
# Server of the funresMech application
############################################################

#' @importFrom shiny reactiveValues reactive observeEvent observe renderTable renderUI
#'   req showNotification withProgress incProgress setProgress downloadHandler
#'   validate need updateSelectInput updateNumericInput updateTextInput isolate
#'   bindCache
#' @importFrom plotly renderPlotly ggplotly config
#' @importFrom future plan multisession sequential
#' @importFrom rmarkdown render
#' @importFrom utils read.csv write.csv head
#' @noRd
NULL

.pill <- function(status) {
  status <- if (is.null(status) || is.na(status)) "none" else status
  lab <- c(ok = "OK", warning = "Warning", fail = "Fail", none = "-")[[status]]
  span(class = paste0("pill pill-", status), lab)
}

# HTML table from a data frame (columns may contain HTML tags)
.html_table <- function(df, html_cols = character(0)) {
  tags$table(class = "table table-sm table-hover",
    tags$thead(tags$tr(lapply(names(df), tags$th))),
    tags$tbody(lapply(seq_len(nrow(df)), function(i)
      tags$tr(lapply(names(df), function(nm) {
        v <- df[[nm]][[i]]
        tags$td(if (nm %in% html_cols) v else as.character(v))
      })))))
}

.reliability_html <- function(rel) {
  rows <- lapply(seq_len(nrow(rel)), function(i) tags$tr(
    tags$td(.pill(rel$status[i])), tags$td(rel$check[i]), tags$td(rel$detail[i]),
    tags$td(class = "muted", rel$action[i])))
  tags$table(class = "table table-sm",
    tags$thead(tags$tr(tags$th("Status"), tags$th("Check"), tags$th("Result"), tags$th("Action"))),
    tags$tbody(rows))
}

.fig_save <- function(file, plot, format) {
  if (format == "png") ggplot2::ggsave(file, plot, device = "png", width = 7, height = 4.5, dpi = 300)
  else ggplot2::ggsave(file, plot, device = "pdf", width = 7, height = 4.5)
}

.plotly <- function(g) {
  config(ggplotly(g), displaylogo = FALSE,
         modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d"))
}

server <- function(input, output, session) {

  sim_cache <- new.env(parent = emptyenv())   # simulations for the figures (not reactive)
  test_out <- reactiveValues(search = NULL, handling = NULL)
  rv <- reactiveValues(fits = list(), screen = NULL, run_opts = NULL, run_species = NULL,
                       run_excluded = integer(0),
                       comparison = NULL, messages = character(0), runtime = NULL)

  # ---- data ------------------------------------------------------------------
  raw <- reactive({
    req(input$dataset)
    d <- tryCatch(read.csv(input$dataset$datapath, check.names = FALSE),
                  error = function(e) NULL)
    validate(need(!is.null(d), "The file could not be read as CSV."),
             need(ncol(d) >= 2, "The file needs at least two columns: density and response."))
    d
  })

  output$column_select <- renderUI({
    d <- raw()
    nm <- names(d)
    guess <- function(pat, default) {
      i <- grep(pat, nm, ignore.case = TRUE)
      if (length(i)) nm[i[1]] else default
    }
    tagList(
      selectInput("dens_col", "Host or prey density", nm, selected = guess("dens|host|prey|^x$", nm[1])),
      selectInput("par_col", "Hosts parasitized or prey consumed", nm,
                  selected = guess("par|eat|cons|attack|^y$", nm[min(2, length(nm))])),
      selectInput("species_col", "Species or group (optional)", c("(none)" = "", nm),
                  selected = guess("sp|species|group|treat", ""))
    )
  })

  mapped <- reactive({
    d <- raw()
    req(input$dens_col, input$par_col)
    sp <- if (!is.null(input$species_col) && nzchar(input$species_col))
      as.character(d[[input$species_col]]) else rep("all trials", nrow(d))
    data.frame(species = sp, dens = d[[input$dens_col]], par = d[[input$par_col]],
               row = seq_len(nrow(d)), stringsAsFactors = FALSE)
  })

  output$species_choice <- renderUI({
    sp <- unique(mapped()$species)
    selectizeInput("species", "Species or groups to analyse", sp, selected = sp, multiple = TRUE)
  })

  output$data_summary <- renderTable({
    m <- mapped()
    do.call(rbind, lapply(split(m, m$species), function(d) {
      chk <- tryCatch(.check_data_spp(d), error = function(e) conditionMessage(e))
      data.frame(Group = d$species[1], Trials = nrow(d), Densities = length(unique(d$dens)),
                 `Density range` = paste(range(d$dens), collapse = " to "),
                 `Mean response` = round(mean(d$par), 2),
                 Check = if (isTRUE(chk)) "valid" else chk, check.names = FALSE)
    }))
  }, spacing = "s")

  output$data_preview <- renderTable(head(raw(), 8), spacing = "s")

  observeEvent(input$screen, {
    m <- mapped(); req(input$species)
    out <- list()
    withProgress(message = "Screening trials", value = 0, {
      for (sp in input$species) {
        d <- m[m$species == sp, ]
        r <- tryCatch(screen_outliers(d[, c("dens", "par")], 1), error = function(e) NULL)
        if (!is.null(r)) {
          r$species <- sp
          r$file_row <- d$row[r$row]
          out[[sp]] <- r
        } else showNotification(paste0(sp, ": the screening model could not be fitted."), type = "warning")
        incProgress(1 / length(input$species))
      }
    })
    rv$screen <- do.call(rbind, out)
  })

  output$screen_ui <- renderUI({
    if (is.null(rv$screen)) return(p(class = "muted", "Not run yet. Use Screen trials in the sidebar."))
    s <- rv$screen[rv$screen$flagged | rv$screen$influences_dispersion, , drop = FALSE]
    if (!nrow(s)) return(p("No trial was flagged."))
    tagList(
      p(class = "muted", "Flagged: very unlikely given the other trials (two-sided leave-one-out probability below the Bonferroni threshold). Influential: overdispersion falls by more than half without the trial. Check these trials in the original records; to refit without the flagged ones, tick Exclude trials flagged as atypical in Model and fit."),
      .html_table(data.frame(Group = s$species, `File row` = s$file_row, Density = s$dens,
                             Response = s$par, `Expected (without it)` = round(s$expected_loo, 1),
                             `LOO probability` = signif(s$p_loo, 2),
                             Flagged = ifelse(s$flagged, "yes", "no"),
                             Influential = ifelse(s$influences_dispersion, "yes", "no"),
                             check.names = FALSE)))
  })

  # ---- options of the model ---------------------------------------------------
  time_unit <- reactive(if (identical(input$time_unit, "other")) input$time_unit_other else input$time_unit)

  opts <- reactive({
    list(search = input$search, handling = input$handling, process = input$process,
         search_code = input$search_code, search_lower = input$search_lower, search_upper = input$search_upper,
         handling_code = input$handling_code, handling_lower = input$handling_lower,
         handling_upper = input$handling_upper,
         fix_k = input$fix_k, k_value = input$k_value, fix_s = input$fix_s, s_value = input$s_value,
         z_mode = input$z_mode, z_value = input$z_value, z_min = input$z_min, z_max = input$z_max,
         z_step = if (identical(input$z_mode, "joint")) input$z_step_joint else input$z_step,
         extend_z = input$extend_z, preset = input$preset, n_sim = input$n_sim,
         itermax = input$itermax, NP = input$NP, n_reeval = input$n_reeval, method = input$method,
         crn = input$crn, polish = input$polish, widen = input$widen, seed = input$seed,
         parallel = input$parallel, workers = input$workers,
         custom_bounds = input$custom_bounds,
         b_lambda_lo = input$b_lambda_lo, b_lambda_hi = input$b_lambda_hi,
         b_h_lo = input$b_h_lo, b_h_hi = input$b_h_hi, b_k_lo = input$b_k_lo, b_k_hi = input$b_k_hi,
         b_cv_lo = input$b_cv_lo, b_cv_hi = input$b_cv_hi, b_z_lo = input$b_z_lo, b_z_hi = input$b_z_hi,
         T_exp = input$T_exp, time_unit = time_unit())
  })

  # Presets fill in the effort settings only when the preset really changes,
  # so that values edited by the user are never overwritten silently.
  last_preset <- shiny::reactiveVal("standard")
  observeEvent(input$preset, {
    if (identical(input$preset, last_preset())) return()
    last_preset(input$preset)
    p <- funres_control(input$preset)
    updateNumericInput(session, "n_sim", value = p$n_sim)
    updateNumericInput(session, "itermax", value = p$itermax)
    updateNumericInput(session, "NP", value = p$NP)
    updateNumericInput(session, "n_reeval", value = p$n_reeval)
  }, ignoreInit = TRUE)

  model <- reactive(tryCatch(.app_model(opts()), error = function(e) e))
  control <- reactive(tryCatch(.app_control(opts()), error = function(e) e))

  output$default_bounds <- renderTable({
    req(input$T_exp > 0)
    b <- .bounds_in_units(input$T_exp)
    fmt <- function(v) formatC(v, digits = 3, format = "g")
    data.frame(Parameter = b$parameter, Lower = fmt(b$lower), Upper = fmt(b$upper))
  }, spacing = "xs", align = "lrr")

  output$spec_ui <- renderUI({
    m <- model(); ctl <- control()
    if (inherits(m, "error")) return(div(class = "notice", conditionMessage(m)))
    if (inherits(ctl, "error")) return(div(class = "notice", conditionMessage(ctl)))
    free <- setdiff(c("a", "h", "z", "k", "s"), names(m$fix))
    tagList(
      tags$table(class = "table table-sm",
        tags$tbody(
          tags$tr(tags$th("Process"), tags$td(if (m$process == "parasitoid") "Parasitoid: hosts not removed, distinct hosts counted" else "Predator: consumed prey removed")),
          tags$tr(tags$th("Search times"), tags$td(paste0(m$search$name, ", mean 1/(a x^z)"))),
          tags$tr(tags$th("Handling times"), tags$td(paste0(m$handling$name, ", mean h"))),
          tags$tr(tags$th("Free parameters"), tags$td(paste(free, collapse = ", "))),
          tags$tr(tags$th("Fixed"), tags$td(if (length(m$fix)) paste(names(m$fix), unlist(m$fix), sep = " = ", collapse = ", ") else "none")),
          tags$tr(tags$th("z"), tags$td(switch(input$z_mode, profile = sprintf("profile from %s to %s in steps of %s", input$z_min, input$z_max, input$z_step),
                                                 joint = "joint estimate, then profile around it", fixed = paste("fixed at", input$z_value)))),
          tags$tr(tags$th("Optimisation"), tags$td(sprintf("%s; n_sim %d, %d generations, population %d, %d re-evaluations",
                                                         ctl$method, ctl$n_sim, ctl$itermax, ctl$NP, ctl$n_reeval))),
          if (!is.null(rv$runtime)) tags$tr(tags$th("Estimated run time"), tags$td(rv$runtime))
        )),
      if (ctl$preset == "exploratory")
        div(class = "notice", "Exploratory settings give a quick first look. Do not report results from them; refit with the standard or final preset.")
    )
  })

  # data of one species, with the exclusions chosen
  species_data <- function(sp) {
    d <- mapped(); d <- d[d$species == sp, ]
    if (isTRUE(input$exclude_flagged) && !is.null(rv$screen)) {
      fl <- rv$screen$file_row[rv$screen$species == sp & rv$screen$flagged]
      d <- d[!d$row %in% fl, ]
    }
    d
  }

  n_points <- function() switch(input$z_mode, fixed = 1,
                                joint = 12, length(seq(input$z_min, input$z_max, by = input$z_step)) + 1)

  observeEvent(input$estimate, {
    req(input$species)
    m <- model(); ctl <- control()
    if (inherits(m, "error") || inherits(ctl, "error")) return()
    d <- species_data(input$species[1])
    sec <- tryCatch(.estimate_runtime(d[, c("dens", "par")], input$T_exp, m, ctl, n_points(),
                                      if (isTRUE(input$parallel)) input$workers else 1),
                    error = function(e) NA_real_)
    rv$runtime <- paste0("about ", .format_duration(sec * length(input$species)),
                         " for ", length(input$species), " group(s), plus the Beta-Binomial reference")
  })

  # ---- samplers -------------------------------------------------------------
  test_sampler_out <- function(kind) {
    o <- opts()
    s <- tryCatch(.sampler_from_code(o[[paste0(kind, "_code")]], o[[paste0(kind, "_lower")]],
                                     o[[paste0(kind, "_upper")]], kind), error = function(e) e)
    if (inherits(s, "error")) return(data.frame(Problem = conditionMessage(s)))
    r <- tryCatch(test_sampler(s, n = 20000), error = function(e) e)
    if (inherits(r, "error")) return(data.frame(Problem = conditionMessage(r)))
    data.frame(Shape = signif(r$shape, 3), Mean = round(r$mean, 3), CV = round(r$cv, 3),
               `us/draw` = round(r$us_per_draw, 2), `Mean = 1` = ifelse(r$mean_ok, "yes", "NO"),
               check.names = FALSE)
  }
  observeEvent(input$test_search, test_out$search <- test_sampler_out("search"))
  observeEvent(input$test_handling, test_out$handling <- test_sampler_out("handling"))
  output$test_search_out <- renderTable(test_out$search, spacing = "xs")
  output$test_handling_out <- renderTable(test_out$handling, spacing = "xs")

  # ---- fit -------------------------------------------------------------------
  observeEvent(input$run, {
    req(input$species)
    m <- model(); ctl <- control()
    if (inherits(m, "error")) { showNotification(conditionMessage(m), type = "error"); return() }
    if (inherits(ctl, "error")) { showNotification(conditionMessage(ctl), type = "error"); return() }
    z_grid <- if (input$z_mode == "profile") tryCatch(.app_z_grid(opts()), error = function(e) e) else 1
    if (inherits(z_grid, "error")) { showNotification(conditionMessage(z_grid), type = "error"); return() }

    old_plan <- plan()
    on.exit(plan(old_plan), add = TRUE)
    if (isTRUE(input$parallel) && input$workers > 1) plan(multisession, workers = input$workers)
    else plan(sequential)

    set.seed(input$seed)
    fits <- list(); msgs <- character(0)
    total <- length(input$species) * n_points()
    withProgress(message = "Fitting", value = 0, {
      for (sp in input$species) {
        d <- species_data(sp)
        prog <- function(z) incProgress(1 / total, detail = sprintf("%s: z = %s", sp,
                                                                    if (is.na(z)) "joint" else format(z)))
        f <- tryCatch(
          fit_funres(d[, c("dens", "par")], input$T_exp, model = m, z_grid = z_grid,
                     control = ctl, extend_z = isTRUE(input$extend_z),
                     z_method = if (input$z_mode == "joint") "joint" else "profile",
                     z_step = input$z_step_joint, time_unit = time_unit(), progress = prog),
          error = function(e) e)
        if (inherits(f, "error")) {
          msgs <- c(msgs, paste0(sp, ": the fit failed (", conditionMessage(f), ")."))
        } else {
          f$species <- sp
          f$excluded_rows <- setdiff(mapped()$row[mapped()$species == sp], d$row)
          fits[[sp]] <- f
        }
      }
    })
    if (length(fits)) {
      rv$fits <- fits
      rv$run_opts <- opts()            # settings actually used (script, report)
      rv$run_species <- names(fits)
      rv$run_excluded <- unlist(lapply(fits, function(f) f$excluded_rows), use.names = FALSE)
      ch <- names(fits)
      for (id in c("rel_species", "res_species", "cmp_species"))
        updateSelectInput(session, id, choices = ch, selected = ch[1])
    }
    rv$messages <- msgs
    for (mm in msgs) showNotification(mm, type = "error", duration = NULL)
  })

  status_strip <- function() {
    if (!length(rv$fits)) return(p(class = "muted", "No fits yet. Set the model in the sidebar and use Fit the model."))
    div(class = "status-strip", lapply(names(rv$fits), function(sp) {
      f <- rv$fits[[sp]]
      st <- attr(reliability(f), "overall")
      div(class = paste0("status-item status-", st),
          div(span(class = "sp", sp), " ", .pill(st)),
          div(class = "meta", sprintf("z = %s, interval %s", signif(f$par[["z"]], 3), .interval_text(f))),
          div(class = "meta", sprintf("%d trials%s; %s", nrow(f$data),
                                      if (length(f$excluded_rows)) paste0(" (", length(f$excluded_rows), " excluded)") else "",
                                      if (is.null(f$verification)) "not verified" else "verified")))
    }))
  }
  output$status_strip <- renderUI(status_strip())
  output$status_strip2 <- renderUI(status_strip())
  output$fit_messages <- renderUI(if (length(rv$messages)) div(class = "notice", lapply(rv$messages, p)))

  fit_of <- function(id) {
    req(length(rv$fits), input[[id]])
    rv$fits[[input[[id]]]]
  }

  # ---- reliability -------------------------------------------------------------
  output$reliability_ui <- renderUI({
    f <- fit_of("rel_species")
    tagList(.reliability_html(reliability(f)),
            if (is.null(f$verification)) p(class = "muted", "Verification not run: the checks above use only the information of the fit. Run the verification before reporting."))
  })
  output$profile_plot_rel <- renderPlotly(.plotly(.plot_profile(fit_of("rel_species"))))

  observeEvent(input$verify, {
    req(length(rv$fits))
    sps <- if (isTRUE(input$v_all)) names(rv$fits) else input$rel_species
    old_plan <- plan()
    on.exit(plan(old_plan), add = TRUE)
    if (isTRUE(input$parallel) && input$workers > 1) plan(multisession, workers = input$workers)
    withProgress(message = "Verifying", value = 0, {
      for (sp in sps) {
        incProgress(0, detail = sp)
        v <- tryCatch(verify_fit(rv$fits[[sp]], n_runs = input$v_runs, effort = input$v_effort,
                                 bounds_factor = input$v_bounds), error = function(e) e)
        if (inherits(v, "error")) showNotification(paste0(sp, ": ", conditionMessage(v)), type = "error")
        else rv$fits[[sp]] <- v
        incProgress(1 / length(sps))
      }
    })
  })

  # ---- results -------------------------------------------------------------------
  # cached by group, number of trials and the fit itself (a new fit or a
  # verification gives a new key)
  sims_of <- function(sp) {
    f <- rv$fits[[sp]]
    key <- paste(sp, input$plot_nsim, paste(signif(f$par, 8), collapse = "_"), f$model$process)
    if (is.null(sim_cache[[key]])) sim_cache[[key]] <- .simulate_fit(f, input$plot_nsim)
    sim_cache[[key]]
  }

  output$prelim_notice <- renderUI({
    f <- fit_of("res_species")
    st <- attr(reliability(f), "overall")
    if (st == "ok" && !is.null(f$verification)) return(NULL)
    div(class = "notice",
        if (st == "fail") "Some reliability checks failed for this fit. Treat these results as unreliable and act on the checks in Reliability."
        else if (is.null(f$verification)) "Preliminary: the fit has not been verified. See Reliability before reporting."
        else "Some reliability checks gave warnings. See Reliability.")
  })

  output$estimates_ui <- renderUI({
    f <- fit_of("res_species")
    tagList(
      .html_table(.estimates_table(f)),
      p(strong("95% interval of z: "), .interval_text(f)),
      p(class = "muted", sprintf("NLL %.3f (Monte Carlo SE %.3f); %d free parameters; %d trials; T = %s %s; settings: %s.",
                                 f$nll, f$nll_se, f$n_par, nrow(f$data), f$T, f$time_unit, f$control$preset))
    )
  })

  output$aic_ui <- renderUI({
    f <- fit_of("res_species")
    a <- .aic_table(f)
    a$AIC <- round(a$AIC, 2); a$`Delta AIC` <- round(a$`Delta AIC`, 2)
    tagList(.html_table(a),
            p(class = "muted", "Delta AIC is relative to this fit: positive values favour this fit. The Beta-Binomial model has the same expected value as the mechanistic model (Okuyama 2026)."))
  })

  output$profile_plot <- renderPlotly(.plotly(.plot_profile(fit_of("res_species"))))
  output$response_plot <- renderPlotly({
    f <- fit_of("res_species")
    .plotly(.plot_response(f, sims_of(input$res_species)))
  })
  output$dist_plot <- renderPlotly({
    f <- fit_of("res_species")
    .plotly(.plot_distribution(f, sims_of(input$res_species), input$dist_type))
  })

  fig_handler <- function(name, make) downloadHandler(
    filename = function() paste0(name, "_", input$res_species, ".", input$fig_format),
    content = function(file) .fig_save(file, make(), input$fig_format))
  output$dl_fig_profile <- fig_handler("profile", function() .plot_profile(fit_of("res_species")))
  output$dl_fig_response <- fig_handler("response", function()
    .plot_response(fit_of("res_species"), sims_of(input$res_species)))
  output$dl_fig_dist <- fig_handler("distribution", function()
    .plot_distribution(fit_of("res_species"), sims_of(input$res_species), input$dist_type))

  # ---- model comparison ----------------------------------------------------------
  observeEvent(input$cmp_add, {
    f <- fit_of("cmp_species")
    row <- .comparison_row(f, input$cmp_label, input$cmp_species)
    row$Trials <- nrow(f$data)
    rv$comparison <- rbind(rv$comparison, row)
    n <- nrow(rv$comparison)
    updateTextInput(session, "cmp_label", value = paste("Model", n + 1))
  })
  observeEvent(input$cmp_clear, rv$comparison <- NULL)

  output$cmp_ui <- renderUI({
    cmp <- rv$comparison
    if (is.null(cmp)) return(p(class = "muted", "No fits added yet."))
    out <- do.call(rbind, lapply(split(cmp, cmp$Species), function(d) {
      d$`Delta AIC` <- round(d$AIC - min(d$AIC), 2)
      d
    }))
    mixed <- any(tapply(out$Trials, out$Species, function(v) length(unique(v)) > 1))
    out$Reliability <- lapply(out$Reliability, .pill)
    tagList(
      if (mixed) div(class = "notice", "Some fits of the same group used different trials (exclusions); their AIC values are not comparable."),
      .html_table(out, html_cols = "Reliability"),
      p(class = "muted", "Delta AIC within each group: 0 for the lowest AIC. Differences below about 2 do not distinguish the models.")
    )
  })

  # ---- export ----------------------------------------------------------------------
  need_fits <- function() validate(need(length(rv$fits) > 0, "Fit the model first."))
  csv_handler <- function(name, make) downloadHandler(
    filename = function() paste0(name, ".csv"),
    content = function(file) { need_fits(); write.csv(make(), file, row.names = FALSE) })

  output$dl_estimates <- csv_handler("funresMech_estimates", function()
    do.call(rbind, lapply(names(rv$fits), function(sp) {
      f <- rv$fits[[sp]]
      data.frame(species = sp, t(f$par), z_low = if (is.null(f$ci)) NA else f$ci$z_low,
                 z_high = if (is.null(f$ci)) NA else f$ci$z_high, nll = f$nll, nll_se = f$nll_se,
                 n_par = f$n_par, f$aic, trials = nrow(f$data), T = f$T, time_unit = f$time_unit,
                 reliability = attr(reliability(f), "overall"),
                 verified = !is.null(f$verification))
    })))
  output$dl_profiles <- csv_handler("funresMech_profiles", function()
    do.call(rbind, lapply(names(rv$fits), function(sp) cbind(species = sp, rv$fits[[sp]]$profile))))
  output$dl_reliability <- csv_handler("funresMech_reliability", function()
    do.call(rbind, lapply(names(rv$fits), function(sp) cbind(species = sp, as.data.frame(reliability(rv$fits[[sp]]))))))
  output$dl_comparison <- csv_handler("funresMech_comparison", function() rv$comparison)

  output$dl_fits <- downloadHandler(
    filename = function() "funresMech_fits.rds",
    content = function(file) { need_fits(); saveRDS(rv$fits, file) })

  # script of the last fit (or of the current settings before any fit)
  script_text <- function() {
    used <- !is.null(rv$run_opts)
    .r_script(if (used) rv$run_opts else opts(),
              if (is.null(input$dataset)) "data.csv" else input$dataset$name,
              if (is.null(input$species_col)) "" else input$species_col,
              input$dens_col, input$par_col,
              if (used) rv$run_species else input$species,
              if (used) rv$run_excluded else integer(0))
  }
  output$dl_script <- downloadHandler(
    filename = function() "funresMech_analysis.R",
    content = function(file) writeLines(script_text(), file))

  output$dl_report <- downloadHandler(
    filename = function() "funresMech_report.html",
    content = function(file) {
      need_fits()
      src <- system.file("app", "report_template.Rmd", package = "funresMech")
      tmp <- tempfile("funresMech_report_"); dir.create(tmp)
      on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
      file.copy(src, file.path(tmp, "report.Rmd"))
      sims <- lapply(names(rv$fits), sims_of); names(sims) <- names(rv$fits)
      withProgress(message = "Writing the report", value = 0.5,
        render(file.path(tmp, "report.Rmd"), output_file = file, intermediates_dir = tmp,
               knit_root_dir = tmp, envir = new.env(parent = globalenv()), quiet = TRUE,
               params = list(fits = rv$fits, sims = sims, comparison = rv$comparison,
                             screen = rv$screen, include_data = isTRUE(input$report_data),
                             script = script_text())))
    })
}
