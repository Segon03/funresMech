# Shiny application: server logic (shiny::testServer) and helpers.
# The whole workflow is run with tiny settings, sequentially.

app_inputs <- list(
  time_unit = "hours", T_exp = 24, exclude_flagged = FALSE, process = "parasitoid",
  search = "gamma", handling = "lognormal", fix_k = FALSE, fix_s = FALSE,
  z_mode = "profile", z_min = 1, z_max = 2, z_step = 1, extend_z = FALSE,
  z_step_joint = 0.5, z_value = 1, n_sim = 150, itermax = 4,
  NP = 12, n_reeval = 2, method = "DEoptim", crn = FALSE, polish = FALSE, widen = FALSE,
  seed = 1, parallel = FALSE, workers = 1, custom_bounds = FALSE,
  search_code = "function(n, shape) rgamma(n, shape = shape, rate = shape)",
  search_lower = 0.01, search_upper = 100, handling_code = "", handling_lower = 0.1,
  handling_upper = 100, k_value = 1, s_value = 0,
  plot_nsim = 200, dist_type = "histogram", report_data = FALSE,
  v_runs = 1, v_effort = 1.5, v_bounds = 10, v_all = FALSE, cmp_label = "Model 1")

app_data <- function() {
  d2 <- read.csv(test_path("fixtures", "D2.csv"))
  d3 <- read.csv(test_path("fixtures", "D3.csv"))
  d <- rbind(cbind(group = "A", d2[seq(1, 73, by = 4), ]), cbind(group = "B", d3[seq(1, 64, by = 4), ]))
  names(d) <- c("group", "hosts", "parasitized")
  f <- tempfile(fileext = ".csv")
  write.csv(d, f, row.names = FALSE)
  f
}

test_that("the user interface builds", {
  ui <- funresMech:::ui()
  html <- as.character(ui)
  expect_true(grepl("Model and fit", html))
  expect_true(grepl("Reliability", html))
  expect_false(grepl("Okuyama model", html, ignore.case = TRUE))
})

test_that("the server runs the whole workflow", {
  skip_on_cran()
  f <- app_data()
  shiny::testServer(funresMech:::server, {
    session$setInputs(dataset = data.frame(name = "data.csv", datapath = f, stringsAsFactors = FALSE))
    session$setInputs(dens_col = "hosts", par_col = "parasitized", species_col = "group",
                      species = c("A", "B"))
    expect_equal(sort(unique(mapped()$species)), c("A", "B"))
    do.call(session$setInputs, app_inputs)
    expect_s3_class(model(), "funres_model")
    expect_s3_class(control(), "funres_control")

    session$setInputs(screen = 1)
    expect_true(is.data.frame(rv$screen))

    session$setInputs(estimate = 1)
    expect_match(rv$runtime, "about")

    suppressWarnings(session$setInputs(run = 1))
    expect_equal(names(rv$fits), c("A", "B"))
    expect_length(rv$messages, 0)
    expect_equal(rv$run_opts$n_sim, 150)

    session$setInputs(rel_species = "A", res_species = "A", cmp_species = "A")
    expect_match(paste(output$reliability_ui$html), "Monte Carlo noise")
    expect_true(nchar(paste(output$estimates_ui$html)) > 0)
    expect_true(nchar(paste(output$aic_ui$html)) > 0)

    suppressWarnings(session$setInputs(verify = 1))
    expect_false(is.null(rv$fits$A$verification))

    session$setInputs(cmp_add = 1)
    expect_equal(nrow(rv$comparison), 1)

    scr <- script_text()
    expect_true(any(grepl("n_sim = 150", scr)))
    expect_silent(parse(text = scr))
  })
})

test_that("custom samplers from code and the R script of a custom model", {
  smp <- funresMech:::.sampler_from_code("function(n, shape) rgamma(n, shape, shape)", 0.1, 10, "x")
  expect_s3_class(smp, "funres_sampler")
  expect_error(funresMech:::.sampler_from_code("1 + ", 0.1, 10, "x"), "could not be evaluated")
  expect_error(funresMech:::.sampler_from_code("42", 0.1, 10, "x"), "must evaluate")
  o <- app_inputs; o$search <- "custom"; o$fix_k <- TRUE; o$k_value <- 2; o$custom_bounds <- TRUE
  o[c("b_lambda_lo", "b_lambda_hi", "b_h_lo", "b_h_hi", "b_k_lo", "b_k_hi", "b_cv_lo", "b_cv_hi",
      "b_z_lo", "b_z_hi")] <- list(1e-4, 100, 1e-5, 24, 1e-5, 100, 0, 10, 0.1, 20)
  m <- funresMech:::.app_model(o)
  expect_true(m$search$user)
  expect_equal(m$fix$k, 2)
  scr <- funresMech:::.r_script(o, "x.csv", "group", "hosts", "parasitized", c("A", "B"), c(3L, 7L))
  expect_silent(parse(text = scr))
  expect_true(any(grepl("funres_sampler", scr)))
  expect_true(any(grepl("bounds = list", scr)))
  expect_true(any(grepl("c(3L, 7L)", scr, fixed = TRUE)))
})

test_that("figures and tables of the application build for a fit", {
  skip_on_cran()
  d <- read.csv(test_path("fixtures", "D3.csv"))
  ctrl <- funres_control(n_sim = 150, itermax = 4, NP = 12, n_reeval = 2, bb = FALSE)
  set.seed(1)
  f <- suppressWarnings(fit_funres(d, 24, z_grid = c(1, 2, 3), control = ctrl, extend_z = FALSE,
                                   refine = FALSE, time_unit = "hours"))
  sims <- funresMech:::.simulate_fit(f, 100)
  for (g in list(funresMech:::.plot_profile(f), funresMech:::.plot_response(f, sims),
                 funresMech:::.plot_distribution(f, sims, "histogram"),
                 funresMech:::.plot_distribution(f, sims, "boxplot"),
                 funresMech:::.plot_distribution(f, sims, "violin")))
    expect_silent(ggplot2::ggplot_build(g))
  expect_equal(nrow(funresMech:::.estimates_table(f)), 5)
  expect_match(funresMech:::.estimates_table(f)$Meaning[1], "per hour per host")
  expect_output(print(f), "per host")
  # the distribution figure shows one point per trial, also for repeated counts
  f2 <- f; f2$data <- rbind(f$data, f$data[1:2, ])
  pts <- ggplot2::ggplot_build(plot(f2, type = "distribution", n_sim = 100))$data[[3]]
  expect_equal(nrow(pts), nrow(f2$data))
  expect_true(all(pts$y < 0))
  expect_true(nrow(funresMech:::.aic_table(f)) >= 1)
  expect_equal(ncol(funresMech:::.comparison_row(f, "m", "D3")), 13)
  expect_equal(funresMech:::.format_duration(30), "30 s")
  expect_equal(funresMech:::.format_duration(600), "10 min")
})

