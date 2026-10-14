# ============================================================================
# 12_validacion_ampliada.R
#
# Fase 5 del plan de respuesta a Okuyama: validación ampliada de funresMech
# 1.2.0 con parámetros conocidos. Brazos:
#
#   A. Simulación amplia, modelo correcto (búsqueda Gamma, manipulación
#      Lognormal, parasitoide). Núcleo factorial z x k x diseño y variaciones
#      de a un factor (s/h, h/T, nivel de encuentros). z verdadero fuera de la
#      grilla (salvo z = 1) para que la cobertura no dependa de que la grilla
#      contenga el valor verdadero (salvedad de la Fase 4). Cada ajuste:
#      preset exploratory + reliability() + verify_fit() + NLL precisa en el
#      estimado y en los parámetros verdaderos.
#   B. Desvíos del modelo: datos generados con búsqueda Weibull o Lognormal,
#      manipulación Gamma, proceso depredador y heterogeneidad individual en a.
#      Se ajusta el modelo por defecto y, cuando existe, el modelo correcto.
#   C. Presupuesto: exploratory (1 corrida), mejor de 3 corridas exploratory
#      (por punto del perfil, con reevaluación independiente) y standard.
#   D. 1.1.0 frente a 1.2.0 con los mismos datos y la misma grilla:
#      D-sim (datos simulados, = subconjunto de A) y D-real (D2 y D3, 3
#      semillas por versión, para ver la estabilidad entre corridas).
#      La 1.1.0 corre en procesos de R separados (callr) con su propia
#      biblioteca, instalada desde la carpeta funresMech_1.1.0_OR.
#   E. Calibración de los chequeos de fiabilidad: se arma al final con los
#      valores de cada chequeo (A, con verify_fit) frente a los resultados
#      reales (cobertura de z y déficit de optimización).
#   F. Reanálisis de D2 y D3 con el preset final (grilla 0.5-5, paso 0.25,
#      extensión y refinamiento), D2 con y sin el ensayo atípico (fila 13 del
#      archivo: dens 10, par 2), con verify_fit(), AIC frente a z = 1 y frente
#      al beta-binomial. Para el material suplementario.
#
# Uso (PC de 23-24 núcleos, R en Windows, desde funresMech_dev):
#   # reiniciar R (Ctrl+Shift+F10)
#   devtools::install()
#   rm(list = ls())
#   RAPIDO <- TRUE;  source("validation/R/12_validacion_ampliada.R")   # prueba (~15 min)
#   rm(list = ls())
#   source("validation/R/12_validacion_ampliada.R")                   # corrida larga
#
# Reanudable: cada trabajo se guarda en validation/results/fase5/. Si se
# interrumpe (o se apaga la PC), volver a correr el script hace solo lo que
# falta. Duración estimada con 22 núcleos: ~2 días (ver tabla al iniciar).
# Resultados: validation/tablas/fase5_*.csv (enviar a Claude).
# ============================================================================

library(funresMech)
library(future)
stopifnot(packageVersion("funresMech") >= "1.2.0")
for (p in c("callr")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

if (!exists("RAPIDO")) RAPIDO <- FALSE
N_WORKERS <- if (RAPIDO) 2L else max(1L, parallel::detectCores() - 1L)
DIR_RES <- file.path("validation", "results", if (RAPIDO) "fase5_rapido" else "fase5")
dir.create(DIR_RES, recursive = TRUE, showWarnings = FALSE)
dir.create("validation/tablas", showWarnings = FALSE)
SUFIJO <- if (RAPIDO) "_rapido" else ""

# ---- versión 1.1.0 en una biblioteca aparte ----------------------------------
if (!exists("RUTA_110")) RUTA_110 <- file.path("..", "funresMech_1.1.0_OR")
LIB_110 <- file.path(tools::R_user_dir("funresMech_validacion", "data"), "lib_1.1.0")
instalar_110 <- function() {
  ok <- file.exists(file.path(LIB_110, "funresMech", "DESCRIPTION")) &&
    utils::packageDescription("funresMech", lib.loc = LIB_110)$Version == "1.1.0"
  if (ok) return(invisible(TRUE))
  if (!file.exists(file.path(RUTA_110, "DESCRIPTION")))
    stop("No encuentro la versión 1.1.0 en ", normalizePath(RUTA_110, mustWork = FALSE),
         ". Defina RUTA_110 <- '<carpeta de funresMech_1.1.0_OR>' antes de source().")
  dsc <- read.dcf(file.path(RUTA_110, "DESCRIPTION"))
  dep <- unlist(strsplit(paste(dsc[, intersect(c("Imports", "LinkingTo", "Depends"),
                                               colnames(dsc))], collapse = ","), ","))
  dep <- trimws(sub("\\(.*", "", dep)); dep <- setdiff(dep[nzchar(dep)], c("R", "methods",
                                                       "stats", "utils", "parallel", "tools"))
  falta <- dep[!vapply(dep, requireNamespace, logical(1), quietly = TRUE)]
  if (length(falta)) install.packages(falta)
  dir.create(LIB_110, recursive = TRUE, showWarnings = FALSE)
  install.packages(normalizePath(RUTA_110), repos = NULL, type = "source", lib = LIB_110)
  stopifnot(utils::packageDescription("funresMech", lib.loc = LIB_110)$Version == "1.1.0")
  invisible(TRUE)
}
instalar_110()
cat("funresMech 1.1.0 disponible en", LIB_110, "\n")

# ---- diseño ------------------------------------------------------------------
T_EXP <- 24
DISENOS <- list(grande  = rep(c(2, 4, 8, 12, 16, 24, 32, 48), each = 8),   # 64 ensayos
                pequeno = rep(c(2, 6, 12, 24, 48), each = 6))              # 30 ensayos
GRILLA <- if (RAPIDO) c(0.5, 1, 2, 3) else seq(0.5, 5, by = 0.5)
GRILLA_REAL <- if (RAPIDO) c(0.5, 1, 2, 3) else seq(0.5, 5, by = 0.25)
N_DATOS <- if (RAPIDO) 1L else 8L          # conjuntos por escenario (A y B)
N_DATOS_CD <- if (RAPIDO) 1L else 4L       # de esos, cuántos entran en C y D
N_SEM_REAL <- if (RAPIDO) 1L else 3L       # semillas por versión en D-real
N_PREC <- if (RAPIDO) 2000L else 20000L    # simulaciones de la NLL precisa (x5)
N_INTERNOS_F <- 2L   # cada ajuste de F (preset final) usa 2 procesos para los puntos del perfil

ctrl_de <- function(preset, ...) {
  if (RAPIDO) {
    p <- switch(preset, exploratory = list(n_sim = 200, itermax = 5),
                standard = list(n_sim = 300, itermax = 8), final = list(n_sim = 400, itermax = 10))
    funres_control(preset, n_sim = p$n_sim, itermax = p$itermax, NP = 12, n_reeval = 2, ...)
  } else funres_control(preset, ...)
}

base <- list(z = 1, k = 1, sh = 0.5, hT = 0.025, E = 20, diseno = "grande",
             gen = "correcto", sigma_a = 0)
esc <- function(...) { e <- modifyList(base, list(...)); as.data.frame(e, stringsAsFactors = FALSE) }

# A: núcleo factorial + variaciones de a un factor (diseño grande)
Z_A <- c(1, 1.75, 2.75); K_A <- c(0.2, 1, 5)
ESC_A <- rbind(
  do.call(rbind, apply(expand.grid(z = Z_A, k = K_A, diseno = names(DISENOS),
                                   stringsAsFactors = FALSE), 1, function(r)
    esc(z = as.numeric(r[["z"]]), k = as.numeric(r[["k"]]), diseno = r[["diseno"]], bloque = "nucleo"))),
  do.call(rbind, lapply(list(list(sh = 1.5, bloque = "s/h = 1.5"),
                             list(hT = 0.1, bloque = "h/T = 0.1"),
                             list(E = 5, bloque = "E = 5")), function(v)
    do.call(rbind, apply(expand.grid(z = Z_A, k = K_A), 1, function(r)
      do.call(esc, c(list(z = r[["z"]], k = r[["k"]]), v)))))))
# B: desvíos del modelo (k = 0.3 para que las familias de búsqueda difieran)
GEN_B <- c("correcto", "busqueda_weibull", "busqueda_lognormal", "manip_gamma",
           "depredador", "heterogeneidad_a")
ESC_B <- do.call(rbind, lapply(GEN_B, function(g) do.call(rbind, lapply(c(1, 2.25), function(z)
  esc(z = z, k = 0.3, sh = if (g == "manip_gamma") 1 else 0.5, gen = g,
      sigma_a = if (g == "heterogeneidad_a") 0.5 else 0, bloque = "B")))))
if (RAPIDO) {
  ESC_A <- ESC_A[c(2, 14), ]
  ESC_B <- ESC_B[ESC_B$gen %in% c("busqueda_weibull", "depredador") & ESC_B$z == 2.25, ]
}
nombre_esc <- function(e) sprintf("%s_z%g_k%g_sh%g_hT%g_E%g_%s", e$gen, e$z, e$k, e$sh, e$hT,
                                  e$E, e$diseno)
ESC_A$id <- nombre_esc(ESC_A); ESC_B$id <- nombre_esc(ESC_B)
stopifnot(!anyDuplicated(ESC_A$id), !anyDuplicated(ESC_B$id))

modelo_gen <- function(e) switch(e$gen,
  busqueda_weibull   = funres_model(search = "weibull"),
  busqueda_lognormal = funres_model(search = "lognormal"),
  manip_gamma        = funres_model(handling = "gamma"),
  depredador         = funres_model(process = "predator"),
  funres_model())
verdad <- function(e) {
  dens <- DISENOS[[e$diseno]]
  H_ref <- exp(mean(log(unique(dens))))
  c(a = e$E / (T_EXP * H_ref^e$z), h = e$hT * T_EXP, z = e$z, k = e$k, s = e$sh * e$hT * T_EXP)
}
datos_sim <- function(e, j) {
  dens <- DISENOS[[e$diseno]]
  p <- verdad(e)
  set.seed(sum(utf8ToInt(e$id)) * 1000 + j)
  if (e$sigma_a > 0) {   # cada individuo con su propia a (media conservada)
    par <- vapply(dens, function(x) {
      pi <- p; pi[["a"]] <- p[["a"]] * exp(stats::rnorm(1, -e$sigma_a^2 / 2, e$sigma_a))
      simulate_response(funres_model(), x, T_EXP, pi, n_sim = 1)$response
    }, numeric(1))
    return(data.frame(dens = dens, par = par))
  }
  out <- do.call(rbind, lapply(unique(dens), function(x)
    simulate_response(modelo_gen(e), x, T_EXP, p, n_sim = sum(dens == x))))
  data.frame(dens = out$dens, par = out$response)
}
leer_real <- function(nombre) read.csv(file.path("tests", "testthat", "fixtures", paste0(nombre, ".csv")))

# ---- utilidades de los trabajos ----------------------------------------------
nll_precisa <- function(par, model, d) {
  mean(vapply(1:5, function(i) funresMech:::.nll_sim(par, model, d, T_EXP, N_PREC), 0))
}
aligerar <- function(f) { f$verification_runs <- NULL; f$joint <- NULL; f }
fila_fit <- function(f, z0 = NA_real_) {
  ci <- f$ci; p <- f$par
  data.frame(z_hat = p[["z"]], a_hat = p[["a"]], h_hat = p[["h"]], k_hat = p[["k"]],
             s_hat = p[["s"]], z_low = ci$z_low, z_high = ci$z_high,
             abierto_bajo = ci$low_censored, abierto_alto = ci$high_censored,
             cubre = (is.na(ci$z_low) || ci$z_low <= z0) && (is.na(ci$z_high) || z0 <= ci$z_high),
             excluye_1 = !is.na(ci$z_low) && ci$z_low > 1,
             nll = f$nll, nll_se = f$nll_se, dAIC_z1 = f$aic$delta_AIC,
             dAIC_BB = f$aic$delta_AIC_BB, n_puntos = nrow(f$profile),
             extendido = f$extended, cotas_ampliadas = paste(f$widened, collapse = ","),
             fiabilidad = attr(reliability(f), "overall"))
}
rel_larga <- function(f) {
  r <- as.data.frame(reliability(f))
  data.frame(check = r$check, status = r$status, value = suppressWarnings(as.numeric(r$value)))
}
ajustar_110 <- function(d, grilla, semilla) {
  callr::r(function(lib, d, grilla, semilla, T_EXP, rapido) {
    library(funresMech, lib.loc = lib)
    stopifnot(as.character(utils::packageVersion("funresMech", lib.loc = lib)) == "1.1.0")
    set.seed(semilla)
    t0 <- proc.time()[["elapsed"]]
    r <- if (rapido) suppressWarnings(funresMech:::fit_profile(d, T_EXP, grilla, n_sim = 200,
                                                                itermax = 5, NP = 12, n_reeval = 2))
         else suppressWarnings(funresMech:::fit_profile(d, T_EXP, grilla))   # valores por defecto de 1.1.0
    r$segundos <- proc.time()[["elapsed"]] - t0
    r
  }, args = list(lib = LIB_110, d = d, grilla = grilla, semilla = semilla, T_EXP = T_EXP,
              rapido = RAPIDO),
  libpath = c(LIB_110, .libPaths()))
}

# ---- lista de trabajos -------------------------------------------------------
# costo: unidades de "punto exploratory de 8 densidades" (~200 s en la PC de
# 23 núcleos), solo para ordenar (los largos primero) y estimar la duración.
tr <- list()
add <- function(brazo, id, j = NA_integer_, variante = "", costo = 1)
  tr[[length(tr) + 1]] <<- data.frame(brazo = brazo, id = id, dato = j, variante = variante,
                                      costo = costo, stringsAsFactors = FALSE)
pts <- length(GRILLA) + 3
for (nm in c("D2", "D2_sin_atipico", "D3"))
  add("F", nm, costo = 1.6 * 20 * (length(GRILLA_REAL) + 4 + 7) / N_INTERNOS_F)
for (v in c("1.1.0", "1.2.0")) for (nm in c("D2", "D3")) for (s in seq_len(N_SEM_REAL))
  add("Dreal", nm, s, v, costo = 1.6 * (length(GRILLA_REAL) + 3) * if (v == "1.1.0") 1.5 else 1)
for (i in seq_len(nrow(ESC_A))) for (j in seq_len(N_DATOS)) {
  e <- ESC_A[i, ]; f_dis <- if (e$diseno == "pequeno") 5 / 8 else 1
  add("A", e$id, j, costo = f_dis * (pts + 7.5))
  if (e$bloque == "nucleo" && e$diseno == "grande" && j <= N_DATOS_CD) {
    add("C", e$id, j, "standard", costo = 5 * pts)
    add("C", e$id, j, "exploratory_2", costo = pts)
    add("C", e$id, j, "exploratory_3", costo = pts)
    add("Dsim", e$id, j, "1.1.0", costo = 1.5 * pts)
  }
}
for (i in seq_len(nrow(ESC_B))) for (j in seq_len(N_DATOS)) {
  e <- ESC_B[i, ]
  add("B", e$id, j, "ajuste_por_defecto", costo = pts)
  if (!e$gen %in% c("correcto", "heterogeneidad_a")) add("B", e$id, j, "ajuste_correcto", costo = pts)
}
trabajos <- do.call(rbind, tr)
trabajos$archivo <- file.path(DIR_RES, gsub("[^A-Za-z0-9_.-]", "_", sprintf(
  "%s__%s__d%s__%s.rds", trabajos$brazo, trabajos$id, trabajos$dato, trabajos$variante)))
trabajos <- trabajos[order(-trabajos$costo), ]
pend <- trabajos[!file.exists(trabajos$archivo), ]

seg_pt <- 200
cat(sprintf("Trabajos: %d en total, %d pendientes; %d núcleos.\n", nrow(trabajos), nrow(pend), N_WORKERS))
if (nrow(pend))
  print(aggregate(cbind(trabajos = 1, horas = costo * seg_pt / 3600) ~ brazo, pend, sum), digits = 3)
if (nrow(pend)) {
  costo_nucleo <- ifelse(pend$brazo == "F", pend$costo * N_INTERNOS_F, pend$costo)
  cat(sprintf("Duración estimada: %.0f h (el trabajo más largo: %.0f h). Orientativa: puede variar +-30%%.\n",
              max(sum(costo_nucleo) * seg_pt / 3600 / N_WORKERS, max(pend$costo) * seg_pt / 3600),
              max(pend$costo) * seg_pt / 3600))
}

correr <- function(t) {
  library(funresMech)
  future::plan(future::sequential)
  t0 <- proc.time()[["elapsed"]]
  semilla <- sum(utf8ToInt(paste(t$brazo, t$id, t$dato, t$variante))) * 7 + 13
  out <- switch(t$brazo,
  A = , C = , Dsim = , B = {
    e <- if (t$brazo == "B") ESC_B[ESC_B$id == t$id, ] else ESC_A[ESC_A$id == t$id, ]
    d <- datos_sim(e, t$dato); p0 <- verdad(e)
    if (t$brazo == "Dsim") {
      r <- ajustar_110(d, GRILLA, semilla)
      m <- funres_model()
      res <- data.frame(z_hat = r$par[["z"]], a_hat = r$par[["a"]], h_hat = r$par[["h"]],
                        k_hat = r$par[["k"]], s_hat = r$par[["s"]], z_low = r$ci$z_low,
                        z_high = r$ci$z_high, abierto_bajo = r$ci$low_censored,
                        abierto_alto = r$ci$high_censored,
                        cubre = (is.na(r$ci$z_low) || r$ci$z_low <= e$z) &&
                          (is.na(r$ci$z_high) || e$z <= r$ci$z_high),
                        excluye_1 = !is.na(r$ci$z_low) && r$ci$z_low > 1,
                        nll = r$nll, nll_se = NA, dAIC_z1 = r$aic$delta_AIC, dAIC_BB = NA,
                        n_puntos = nrow(r$profile), extendido = r$extended, cotas_ampliadas = "",
                        fiabilidad = NA)
      res$nll_prec_hat <- nll_precisa(r$par, m, d)
      res$nll_prec_true <- nll_precisa(p0, m, d)
      list(res = res, perfil = r$profile)
    } else {
      m <- if (t$brazo == "B" && t$variante == "ajuste_correcto") modelo_gen(e) else funres_model()
      preset <- if (t$variante == "standard") "standard" else "exploratory"
      set.seed(semilla)
      f <- suppressWarnings(fit_funres(d, T_EXP, m, z_grid = GRILLA,
                                       control = ctrl_de(preset, bb = t$brazo %in% c("A", "B"))))
      if (t$brazo == "A") f <- verify_fit(f)
      res <- fila_fit(f, e$z)
      if (t$brazo %in% c("A", "C")) {
        res$nll_prec_hat <- nll_precisa(f$par, m, d)
        res$nll_prec_true <- nll_precisa(p0, m, d)
      }
      list(res = res, rel = rel_larga(f), fit = aligerar(f))
    }
  },
  Dreal = {
    d <- leer_real(t$id)
    if (t$variante == "1.1.0") {
      r <- ajustar_110(d, GRILLA_REAL, semilla)
      list(res = data.frame(z_hat = r$par[["z"]], a_hat = r$par[["a"]], h_hat = r$par[["h"]],
                            k_hat = r$par[["k"]], s_hat = r$par[["s"]], z_low = r$ci$z_low,
                            z_high = r$ci$z_high, nll = r$nll,
                            nll_prec_hat = nll_precisa(r$par, funres_model(), d)),
           perfil = r$profile)
    } else {
      set.seed(semilla)
      f <- suppressWarnings(fit_funres(d, T_EXP, z_grid = GRILLA_REAL,
                                       control = ctrl_de("exploratory", bb = FALSE)))
      res <- fila_fit(f)[, c("z_hat", "a_hat", "h_hat", "k_hat", "s_hat", "z_low", "z_high", "nll")]
      res$nll_prec_hat <- nll_precisa(f$par, funres_model(), d)
      list(res = res, perfil = f$profile, rel = rel_larga(f))
    }
  },
  F = {
    options(parallelly.maxWorkers.localhost = Inf)
    future::plan(future::multisession, workers = N_INTERNOS_F)
    on.exit(future::plan(future::sequential), add = TRUE)
    nm <- sub("_sin_atipico", "", t$id)
    d <- leer_real(nm)
    scr <- funresMech:::screen_outliers(d, T_EXP)
    quitadas <- integer(0)
    if (grepl("sin_atipico", t$id)) {
      quitadas <- scr$row[scr$flagged]
      stopifnot(identical(quitadas, 13L), d$dens[13] == 10, d$par[13] == 2)
      d <- d[-quitadas, ]
    }
    arch_fit <- sub("\\.rds$", "_ajuste.rds", t$archivo)
    if (file.exists(arch_fit)) f <- readRDS(arch_fit) else {
      set.seed(semilla)
      f <- fit_funres(d, T_EXP, z_grid = GRILLA_REAL, control = ctrl_de("final"))
      saveRDS(f, arch_fit)                     # punto de control antes de verify_fit()
    }
    f <- verify_fit(f)
    res <- cbind(data.frame(datos = t$id, ensayos = nrow(d),
                            filas_quitadas = paste(quitadas, collapse = ",")),
                 fila_fit(f), f$aic[, c("AIC_full", "AIC_restricted", "AIC_BB")])
    list(res = res, rel = rel_larga(f), fit = aligerar(f), screening = scr)
  })
  out$segundos <- proc.time()[["elapsed"]] - t0
  out$trabajo <- t[, c("brazo", "id", "dato", "variante")]
  saveRDS(out, t$archivo)
  TRUE
}

if (nrow(pend)) {
  # los trabajos F usan N_INTERNOS_F procesos propios: el resto del grupo se achica
  n_F <- sum(pend$brazo == "F")
  plan(multisession, workers = max(1L, N_WORKERS - n_F * (N_INTERNOS_F - 1L)))
  t_ini <- Sys.time()
  glob <- c("correr", "pend", "ESC_A", "ESC_B", "DISENOS", "T_EXP", "GRILLA", "GRILLA_REAL",
            "N_PREC", "N_INTERNOS_F", "RAPIDO", "LIB_110", "ctrl_de", "modelo_gen", "verdad", "datos_sim",
            "leer_real", "nll_precisa", "aligerar", "fila_fit", "rel_larga", "ajustar_110")
  fs <- vector("list", nrow(pend))
  for (i in seq_len(nrow(pend))) {
    fs[[i]] <- future(correr(pend[i, ]), seed = TRUE,
                      globals = c(mget(glob), list(i = i)))
    cat(sprintf("[%s] iniciado %d/%d: %s %s %s %s\n", format(Sys.time(), "%d %H:%M"), i,
                nrow(pend), pend$brazo[i], pend$id[i], pend$dato[i], pend$variante[i]))
  }
  errores <- 0
  for (i in seq_along(fs)) {
    tryCatch(value(fs[[i]]), error = function(e) {
      errores <<- errores + 1
      msg <- sprintf("%s | %s %s %s %s | %s", format(Sys.time()), pend$brazo[i], pend$id[i],
                     pend$dato[i], pend$variante[i], conditionMessage(e))
      message("Falló: ", msg)
      cat(msg, "\n", file = file.path(DIR_RES, "errores.log"), append = TRUE)
    })
  }
  plan(sequential)
  cat("Tiempo total:", format(Sys.time() - t_ini), "; errores:", errores, "\n")
}

# ---- tablas ------------------------------------------------------------------
cargar <- function(brazo) {
  a <- trabajos$archivo[trabajos$brazo == brazo & file.exists(trabajos$archivo)]
  lapply(a, readRDS)
}
fila <- function(o, e = NULL) {
  r <- cbind(o$trabajo, o$res, segundos = o$segundos)
  if (!is.null(e)) r <- cbind(r, e[, c("z", "k", "sh", "hT", "E", "diseno", "gen", "bloque")])
  r
}
escribir <- function(x, nombre) {
  if (!is.null(x) && nrow(x)) write.csv(x, sprintf("validation/tablas/fase5_%s%s.csv", nombre, SUFIJO),
                                        row.names = FALSE)
}
resumir <- function(x, por) {
  x$err_z <- x$z_hat - x$z
  x$abs_err_z <- abs(x$err_z)
  x$abierto <- x$abierto_bajo | x$abierto_alto
  x$z1_detectado <- x$excluye_1
  cols <- intersect(c("cubre", "abierto", "abs_err_z", "err_z", "z1_detectado", "deficit",
                      "fallo_opt", "segundos"), names(x))
  aggregate(x[, cols], x[, por, drop = FALSE], function(v) mean(as.numeric(v), na.rm = TRUE))
}
OA <- cargar("A")
A <- if (length(OA)) do.call(rbind, lapply(OA, function(o) fila(o, ESC_A[ESC_A$id == o$trabajo$id, ])))
if (!is.null(A)) {
  A$deficit <- A$nll_prec_hat - A$nll_prec_true     # > 0: el óptimo no llegó ni a la verdad
  A$fallo_opt <- A$deficit > 1
  A$h_rel <- A$h_hat / (A$hT * T_EXP)
  escribir(A, "A_simulacion")
  RA <- resumir(A, c("bloque", "diseno", "z", "k")); escribir(RA, "A_resumen")
  cat("\nA. Cobertura y error de z por bloque:\n")
  print(resumir(A, "bloque"), digits = 3)
}

OB <- cargar("B")
B <- if (length(OB)) do.call(rbind, lapply(OB, function(o) fila(o, ESC_B[ESC_B$id == o$trabajo$id, ])))
if (!is.null(B)) {
  escribir(B, "B_misespecificacion")
  RB <- resumir(B, c("gen", "variante", "z")); escribir(RB, "B_resumen")
  cat("\nB. Desvíos del modelo:\n"); print(RB, digits = 3)
}

# C: presupuesto. El brazo "exploratory" es el ajuste de A con el mismo conjunto.
OC <- cargar("C")
if (length(OC) && !is.null(A)) {
  C <- do.call(rbind, lapply(OC, function(o) fila(o, ESC_A[ESC_A$id == o$trabajo$id, ])))
  C$deficit <- C$nll_prec_hat - C$nll_prec_true
  C1 <- A[A$id %in% C$id & A$dato <= N_DATOS_CD, ]
  C1$variante <- "exploratory"
  comunes <- intersect(names(C), names(C1))
  CC <- rbind(C1[, comunes], C[, comunes])
  # mejor de 3: por punto del perfil, el mejor de las 3 corridas exploratory,
  # reevaluado con simulaciones nuevas (sin sesgo de selección)
  m3 <- list()
  for (k in unique(paste(C$id, C$dato))) {
    oa <- OA[vapply(OA, function(o) paste(o$trabajo$id, o$trabajo$dato) == k, logical(1))]
    oc <- OC[vapply(OC, function(o) paste(o$trabajo$id, o$trabajo$dato) == k &&
                      grepl("exploratory", o$trabajo$variante), logical(1))]
    fits <- lapply(c(oa, oc), `[[`, "fit")
    if (length(fits) < 3) next
    e <- ESC_A[ESC_A$id == oa[[1]]$trabajo$id, ]
    d <- fits[[1]]$data
    P <- do.call(rbind, lapply(fits, `[[`, "profile"))
    P <- do.call(rbind, lapply(split(P, round(P$z, 8)), function(q) q[which.min(q$nll), ]))
    ctl <- fits[[1]]$control
    set.seed(sum(utf8ToInt(k)))
    P$nll <- vapply(seq_len(nrow(P)), function(i) mean(vapply(seq_len(ctl$n_reeval), function(r)
      funresMech:::.nll_sim(c(a = P$a[i], h = P$h[i], z = P$z[i], k = P$k[i], s = P$s[i]),
                            funres_model(), d, T_EXP, ctl$n_sim), 0)), 0)
    P <- P[order(P$z), ]
    ci <- funresMech:::.profile_ci(P$z, P$nll)
    ib <- which.min(P$nll)
    par <- c(a = P$a[ib], h = P$h[ib], z = P$z[ib], k = P$k[ib], s = P$s[ib])
    m3[[k]] <- data.frame(id = e$id, dato = oa[[1]]$trabajo$dato, variante = "mejor_de_3",
                          z_hat = par[["z"]], z_low = ci$z_low, z_high = ci$z_high,
                          abierto_bajo = ci$low_censored, abierto_alto = ci$high_censored,
                          cubre = (is.na(ci$z_low) || ci$z_low <= e$z) && (is.na(ci$z_high) || e$z <= ci$z_high),
                          excluye_1 = !is.na(ci$z_low) && ci$z_low > 1,
                          nll_prec_hat = nll_precisa(par, funres_model(), d),
                          nll_prec_true = oa[[1]]$res$nll_prec_true,
                          segundos = sum(vapply(c(oa, oc), function(o) o$segundos, 0)))
  }
  if (length(m3)) {
    M3 <- do.call(rbind, m3)
    M3 <- merge(M3, ESC_A[, c("id", "z", "k", "sh", "hT", "E", "diseno", "gen", "bloque")], by = "id")
    M3$deficit <- M3$nll_prec_hat - M3$nll_prec_true
    comunes <- intersect(names(CC), names(M3))
    CC <- rbind(CC[, comunes], M3[, comunes])
  }
  # éxito: NLL precisa a menos de 0.5 (y de 1) de la mejor de todas las variantes
  mejor <- tapply(CC$nll_prec_hat, paste(CC$id, CC$dato), min)
  CC$brecha <- CC$nll_prec_hat - mejor[paste(CC$id, CC$dato)]
  CC$exito_05 <- CC$brecha < 0.5; CC$exito_1 <- CC$brecha < 1
  escribir(CC, "C_presupuesto")
  CC$abierto <- CC$abierto_bajo | CC$abierto_alto
  RC <- aggregate(cbind(exito_05, exito_1, brecha, cubre, abierto, deficit, segundos) ~ variante,
                  CC, function(v) mean(as.numeric(v)))
  escribir(RC, "C_resumen")
  cat("\nC. Éxito según el presupuesto:\n"); print(RC, digits = 3)
}

# D: 1.1.0 frente a 1.2.0
OD <- cargar("Dsim")
if (length(OD) && !is.null(A)) {
  D110 <- do.call(rbind, lapply(OD, function(o) fila(o, ESC_A[ESC_A$id == o$trabajo$id, ])))
  D110$version <- "1.1.0"
  D120 <- A[paste(A$id, A$dato) %in% paste(D110$id, D110$dato), ]
  D120$version <- "1.2.0"
  comunes <- intersect(names(D110), names(D120))
  DS <- rbind(D110[, comunes], D120[, comunes])
  DS$deficit <- DS$nll_prec_hat - DS$nll_prec_true
  mejor <- tapply(DS$nll_prec_hat, paste(DS$id, DS$dato), min)
  DS$brecha <- DS$nll_prec_hat - mejor[paste(DS$id, DS$dato)]
  DS$abierto <- DS$abierto_bajo | DS$abierto_alto
  DS$abs_err_z <- abs(DS$z_hat - DS$z)
  escribir(DS, "D_110_vs_120_simulados")
  RD <- aggregate(cbind(cubre, abierto, abs_err_z, excluye_1, deficit, brecha, segundos) ~ version + z,
                  DS, function(v) mean(as.numeric(v)))
  escribir(RD, "D_resumen")
  cat("\nD. 1.1.0 frente a 1.2.0 (datos simulados):\n"); print(RD, digits = 3)
}
ODR <- cargar("Dreal")
if (length(ODR)) {
  DR <- do.call(rbind, lapply(ODR, function(o) cbind(o$trabajo, o$res, segundos = o$segundos)))
  escribir(DR, "D_real")
  cat("\nD. D2 y D3 con cada versión (una fila por semilla):\n")
  print(DR[order(DR$id, DR$variante), c("id", "variante", "dato", "z_hat", "z_low", "z_high",
                                        "nll", "nll_prec_hat")], digits = 4)
}

# E: calibración de los chequeos de fiabilidad (datos de A, con verify_fit)
if (length(OA)) {
  E <- do.call(rbind, lapply(OA, function(o) {
    if (is.null(o$rel)) return(NULL)
    data.frame(id = o$trabajo$id, dato = o$trabajo$dato, o$rel,
               no_cubre = !o$res$cubre, fallo_opt = (o$res$nll_prec_hat - o$res$nll_prec_true) > 1,
               row.names = NULL)
  }))
  escribir(E, "E_chequeos")
  tasa <- function(x, y) if (length(x) && any(!is.na(x))) mean(x, na.rm = TRUE) else NA_real_
  EC <- do.call(rbind, lapply(split(E, E$check), function(q) {
    alerta <- q$status != "ok"
    data.frame(check = q$check[1], n = nrow(q), alertas = mean(alerta),
               sens_no_cubre = tasa(alerta[q$no_cubre]), espec_no_cubre = tasa(!alerta[!q$no_cubre]),
               sens_fallo_opt = tasa(alerta[q$fallo_opt]), espec_fallo_opt = tasa(!alerta[!q$fallo_opt]))
  }))
  # global (cualquier falla / cualquier advertencia o falla), sin el chequeo del preset
  G <- do.call(rbind, lapply(split(E, paste(E$id, E$dato)), function(q) {
    q <- q[!grepl("Settings suitable", q$check), ]
    data.frame(falla = any(q$status == "fail"), alerta = any(q$status != "ok"),
               no_cubre = q$no_cubre[1], fallo_opt = q$fallo_opt[1])
  }))
  EC <- rbind(EC, do.call(rbind, lapply(c("falla", "alerta"), function(v) data.frame(
    check = paste0("GLOBAL: cualquier ", v), n = nrow(G), alertas = mean(G[[v]]),
    sens_no_cubre = tasa(G[[v]][G$no_cubre]), espec_no_cubre = tasa(!G[[v]][!G$no_cubre]),
    sens_fallo_opt = tasa(G[[v]][G$fallo_opt]), espec_fallo_opt = tasa(!G[[v]][!G$fallo_opt])))))
  rownames(EC) <- NULL
  escribir(EC, "E_calibracion")
  cat("\nE. Chequeos de fiabilidad (sensibilidad y especificidad):\n"); print(EC, digits = 2)
}

# F: D2 y D3 con el preset final
OF <- cargar("F")
if (length(OF)) {
  FF <- do.call(rbind, lapply(OF, function(o) cbind(o$res, segundos = o$segundos)))
  escribir(FF, "F_D2_D3_final")
  escribir(do.call(rbind, lapply(OF, function(o) data.frame(datos = o$trabajo$id, o$rel,
                                                            row.names = NULL))), "F_chequeos")
  escribir(do.call(rbind, lapply(OF, function(o) data.frame(datos = o$trabajo$id, o$fit$profile,
                                                            row.names = NULL))), "F_perfiles")
  cat("\nF. D2 y D3 (preset final):\n")
  print(FF[, c("datos", "ensayos", "z_hat", "z_low", "z_high", "dAIC_z1", "dAIC_BB", "fiabilidad")],
        digits = 4)
}
cat("\nListo. Enviar a Claude los archivos validation/tablas/fase5_*", SUFIJO, ".csv\n", sep = "")
