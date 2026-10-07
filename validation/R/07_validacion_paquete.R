# ============================================================
# R/07_validacion_paquete.R
#
# VALIDACION DEL PAQUETE funresMech 1.1.0 (no del codigo del proyecto)
# con los datos reales D2 y D3 (Okuyama 2026, Tabla S1).
#
# Hasta 04k se valido el motor del proyecto (fun/motor_rcpp.R +
# fun/motor_okuyama.cpp). Este script ejecuta el CODIGO REAL DEL PAQUETE
# INSTALADO (funresMech:::negloglik_fixed_z, funresMech:::screen_outliers,
# funresMech:::fit_profile) con la misma configuracion que 04h y lo
# compara con lo ya validado y con Okuyama. Es lo mismo que hace la app
# Shiny al presionar "Run", pero sin interfaz, con semilla fija y con
# resultados guardados en CSV.
#
# FASES
#   A. Equivalencia de la NLL: mismos parametros y misma semilla ->
#      nll_motor() del proyecto vs funresMech:::negloglik_fixed_z().
#      Deben coincidir exactamente (diferencia 0).
#   B. Cribado de atipicos del paquete (screen_outliers): en D2 debe
#      marcar solo el ensayo 2/10; en D3 ninguno.
#   C. Perfil completo con funresMech:::fit_profile() en paralelo
#      (grilla 0.5-6.0 paso 0.1; el paquete extiende y refina solo).
#   D. Comparacion contra el proyecto (04h/04i) y contra Okuyama.
#
# USO (desde la carpeta del proyecto, la del .Rproj)
#   0) Instalar UNA vez el paquete (en una sesion de R nueva, sin otros
#      procesos usando el paquete). Desde GitHub (rama main = 1.1.0; queda
#      registrado el commit instalado, que este script guarda en los
#      resultados) o, cuando CRAN lo publique, desde CRAN:
#        install.packages("remotes")
#        remotes::install_github("Segon03/funresMech")
#        # alternativa: install.packages("funresMech")
#   1) Prueba de humo (minutos): Rscript R/07_validacion_paquete.R prueba
#   2) Corrida final, dos procesos a la vez para aprovechar los 23 nucleos
#      (cada uno con ~11 workers; el perfil de D2 tiene etapas con pocos
#      puntos y D3 ocupa los nucleos libres):
#        Rscript R/07_validacion_paquete.R D2 11
#        Rscript R/07_validacion_paquete.R D3 11
#      o, desde RStudio, source() del script (usa D2 y D3 en serie con
#      todos los nucleos menos uno).
#   3) Solo comparar con resultados ya guardados:
#        Rscript R/07_validacion_paquete.R comparar
#
# ARGUMENTOS (en cualquier orden): D2 / D3 (datasets; por defecto ambos),
#   un numero (workers), prueba (corrida de humo), comparar (solo fase D).
#
# NOTAS
#  * fit_profile no guarda puntos intermedios: si se corta, esa corrida se
#    repite. Evitar que la PC entre en suspension (powercfg /change
#    standby-timeout-ac 0) y no cerrar la sesion de R. Si ya existe
#    results/paquete_1.1.0/fit_<lab>.rds, ese dataset se salta.
#  * Semilla fija por dataset (SEEDS): el perfil es reproducible.
#  * reltol: se usa el default de DEoptim (sqrt(eps), el que uso 04h al
#    no pasar reltol). El default del paquete/app es 1e-2; con
#    steptol = itermax = 200 la diferencia practica es nula.
#  * n_reeval = 5 (como 04h; el default del paquete es 3). Solo afecta la
#    precision de nll_reeval, no el ajuste.
# ============================================================

suppressWarnings(suppressMessages({
  library(future)
}))

# ---- Configuracion ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
LABS        <- c("D2", "D3")
N_WORKERS   <- max(1, parallel::detectCores() - 1)
MODO_PRUEBA <- FALSE
SOLO_COMPARAR <- FALSE
if (length(args)) {
  labs_arg <- toupper(args[toupper(args) %in% c("D2", "D3")])
  if (length(labs_arg)) LABS <- unique(labs_arg)
  num <- suppressWarnings(as.integer(args[grepl("^[0-9]+$", args)]))
  if (length(num) && !is.na(num[1])) N_WORKERS <- max(1L, num[1])
  MODO_PRUEBA   <- any(tolower(args) == "prueba")
  SOLO_COMPARAR <- any(tolower(args) == "comparar")
}

T_EXP    <- 24
N_SIM    <- 10000                    # Okuyama (2026), Script S3
ITERMAX  <- 200                      # Okuyama (2026), Script S3
NP       <- 40
RELTOL   <- sqrt(.Machine$double.eps)  # default de DEoptim (lo que uso 04h)
N_REEVAL <- 5L
Z_GRID   <- round(seq(0.5, 6.0, by = 0.1), 2)
SEEDS    <- c(D2 = 20261006L, D3 = 20261007L)
Z_EQ     <- c(2.5, 4.0)              # puntos de la fase A
SEEDS_EQ <- c(1L, 2L, 3L)
N_SIM_EQ <- 10000
TOL_NLL  <- 0.62                     # tolerancia de ruido Monte Carlo (auditoria 3.3)
TOL_Z    <- 0.3                      # tolerancia para z_hat e IC (3 pasos de grilla)
TOL_REL  <- 0.05                     # tolerancia relativa para h

if (MODO_PRUEBA) {
  N_SIM <- 500; ITERMAX <- 3; N_REEVAL <- 2L
  Z_GRID <- c(1.5, 2.5, 3.5); N_SIM_EQ <- 500
}
SUF     <- if (MODO_PRUEBA) "_PRUEBA" else ""
dir_res <- file.path("results", paste0("paquete_1.1.0", SUF))
dir.create(dir_res, showWarnings = FALSE, recursive = TRUE)
dir.create("tablas", showWarnings = FALSE)

# ---- Paquete instalado ------------------------------------------------------------
# No depende de ninguna carpeta del paquete en esta PC: usa el paquete INSTALADO
# y registra de donde viene (CRAN o GitHub + commit) para dejarlo en los resultados.
cmd_instalar <- 'remotes::install_github("Segon03/funresMech")'
if (!requireNamespace("funresMech", quietly = TRUE))
  stop("funresMech no esta instalado. En una sesion de R nueva ejecutar:\n  ", cmd_instalar, call. = FALSE)
v_inst <- as.character(utils::packageVersion("funresMech"))
if (numeric_version(v_inst) < "1.1.0")
  stop("Se necesita funresMech >= 1.1.0 (instalado: ", v_inst, "). Reinstalar:\n  ",
       cmd_instalar, call. = FALSE)
desc <- utils::packageDescription("funresMech")
origen <- if (!is.null(desc$RemoteSha)) {
  sprintf("GitHub %s/%s @ %s", desc$RemoteUsername, desc$RemoteRepo, substr(desc$RemoteSha, 1, 7))
} else if (!is.null(desc$Repository)) desc$Repository else "instalacion local"
cat(sprintf("funresMech %s instalado (%s). Datasets: %s. Workers: %d.%s\n",
            v_inst, origen, paste(LABS, collapse = "+"), N_WORKERS,
            if (MODO_PRUEBA) "  [MODO PRUEBA]" else ""))
if (!file.exists(file.path("fun", "motor_rcpp.R")) || !dir.exists("data_clean"))
  stop("Correr desde la carpeta del proyecto (la que contiene fun/, data_clean/ y tablas/). ",
       "Carpeta actual: ", getwd(), call. = FALSE)

# ---- Datos ---------------------------------------------------------------------------
datasets <- setNames(lapply(c("D2", "D3"), function(lab) {
  f <- file.path("data_clean", paste0(lab, ".rds"))
  if (!file.exists(f)) stop("No se encontro ", f, call. = FALSE)
  as.data.frame(readRDS(f))
}), c("D2", "D3"))

leer <- function(f) read.csv(file.path("tablas", f), stringsAsFactors = FALSE)

# ============================================================
# FASE A - Equivalencia de la NLL (proyecto vs paquete)
# ============================================================
fase_A <- function(lab) {
  source(file.path("fun", "motor_rcpp.R"))
  cargar_motor_rcpp()
  datos <- datasets[[lab]]
  pp <- leer("ensayo2_rcpp_perfiles.csv")
  pp <- pp[pp$escala == "log" & pp$dataset == lab, ]
  filas <- list()
  for (z in Z_EQ) {
    fila <- pp[abs(pp$z - z) < 1e-8, ][1, ]
    if (nrow(fila) == 0 || is.na(fila$a)) next
    par <- c(fila$a, fila$h, fila$k, fila$s)
    for (sd in SEEDS_EQ) {
      set.seed(sd); n_proj <- nll_motor(par, z, datos, T_EXP, N_SIM_EQ)
      set.seed(sd); n_pkg  <- funresMech:::negloglik_fixed_z(par, z, datos, T_EXP, N_SIM_EQ)
      filas[[length(filas) + 1]] <- data.frame(
        dataset = lab, z = z, semilla = sd, n_sim = N_SIM_EQ,
        a = par[1], h = par[2], k = par[3], s = par[4],
        nll_proyecto = n_proj, nll_paquete = n_pkg,
        diferencia = n_pkg - n_proj,
        veredicto = if (isTRUE(abs(n_pkg - n_proj) < 1e-9)) "IDENTICO" else "DIFIERE")
    }
  }
  eq <- do.call(rbind, filas)
  write.csv(eq, file.path("tablas", paste0("paquete_equivalencia_nll_", lab, SUF, ".csv")),
            row.names = FALSE)
  cat(sprintf("\n[A] Equivalencia NLL %s: %d/%d casos identicos (max |dif| = %.3g)\n",
              lab, sum(eq$veredicto == "IDENTICO"), nrow(eq), max(abs(eq$diferencia))))
  eq
}

# ============================================================
# FASE B - Cribado de atipicos del paquete
# ============================================================
fase_B <- function(lab) {
  datos <- datasets[[lab]]
  t0 <- Sys.time()
  scr <- funresMech:::screen_outliers(datos, T_EXP)
  scr$dataset <- lab
  write.csv(scr, file.path("tablas", paste0("paquete_cribado_", lab, SUF, ".csv")),
            row.names = FALSE)
  marc <- sort(scr$row[scr$flagged])
  ref <- leer("diagnostico_cribado_D2_D3.csv")
  ref <- sort(ref$fila[ref$dataset == lab & ref$atipico])
  cat(sprintf("[B] Cribado %s: marcados = {%s}; proyecto = {%s} -> %s (%.0f s)\n",
              lab, paste(marc, collapse = ","), paste(ref, collapse = ","),
              if (identical(as.integer(marc), as.integer(ref))) "COINCIDE" else "DIFIERE",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  invisible(scr)
}

# ============================================================
# FASE C - Perfil completo con funresMech:::fit_profile()
# ============================================================
fase_C <- function(lab) {
  f_rds <- file.path(dir_res, paste0("fit_", lab, ".rds"))
  if (file.exists(f_rds)) {
    cat(sprintf("[C] %s: ya existe %s (se salta; borrarlo para repetir)\n", lab, f_rds))
    return(invisible(readRDS(f_rds)))
  }
  datos <- datasets[[lab]]

  # Estimacion de tiempo con el motor del paquete
  t_eval <- system.time(for (i in 1:3)
    funresMech:::negloglik_fixed_z(c(0.02, 0.5, 1, 0.5), 2, datos, T_EXP, N_SIM))[["elapsed"]] / 3
  evals <- NP * (ITERMAX + 1) + N_REEVAL
  n_pts <- length(Z_GRID)
  horas <- ceiling(n_pts / N_WORKERS) * evals * t_eval / 3600
  cat(sprintf(paste0("\n[C] %s: 1 evaluacion = %.2f s; %d puntos x %d evaluaciones; ",
                     "%d workers -> ~%.1f h solo la grilla base (mas extension/refinamiento si el IC queda abierto)\n"),
              lab, t_eval, n_pts, evals, N_WORKERS, horas))

  plan(multisession, workers = N_WORKERS)
  on.exit(plan(sequential), add = TRUE)
  if (!isTRUE(funresMech:::.workers_ok()))
    stop("Los workers no tienen la misma version instalada del paquete; reinstalar funresMech.",
         call. = FALSE)

  set.seed(SEEDS[[lab]])
  t0 <- Sys.time()
  res <- funresMech:::fit_profile(
    data_spp = datos, T_exp = T_EXP, z_grid = Z_GRID,
    n_sim = N_SIM, itermax = ITERMAX, NP = NP, reltol = RELTOL,
    n_reeval = N_REEVAL, extend_z = TRUE, refine = TRUE)
  res$meta <- list(dataset = lab, semilla = SEEDS[[lab]], n_sim = N_SIM, itermax = ITERMAX,
                   NP = NP, reltol = RELTOL, n_reeval = N_REEVAL, z_grid = Z_GRID,
                   workers = N_WORKERS, version_paquete = v_inst,
                   horas = as.numeric(difftime(Sys.time(), t0, units = "hours")),
                   origen_paquete = origen,
                   fecha = format(Sys.time(), "%Y-%m-%d %H:%M"))
  saveRDS(res, f_rds)
  write.csv(res$profile, file.path("tablas", paste0("paquete_perfil_", lab, SUF, ".csv")),
            row.names = FALSE)
  cat(sprintf("[C] %s listo en %.2f h: z_hat = %.2f, IC95 = [%s, %s], extendido = %s, refinado = %s\n",
              lab, res$meta$horas, res$ci$z_hat,
              format(round(res$ci$z_low, 2)),
              if (res$ci$high_censored) paste0(">= ", max(res$profile$z)) else format(round(res$ci$z_high, 2)),
              res$extended, res$refined))
  if (length(res$notes)) cat("    Notas del paquete:\n", paste0("     - ", res$notes, collapse = "\n"), "\n")
  invisible(res)
}

# ============================================================
# FASE D - Comparacion con el proyecto y con Okuyama
# ============================================================
fase_D <- function(lab) {
  f_rds <- file.path(dir_res, paste0("fit_", lab, ".rds"))
  if (!file.exists(f_rds)) { cat(sprintf("[D] %s: no hay %s\n", lab, f_rds)); return(invisible(NULL)) }
  res <- readRDS(f_rds)
  prof <- res$profile
  ci <- res$ci

  rp <- leer("ensayo2_rcpp_resultado.csv")
  rp <- rp[rp$dataset == lab & rp$escala == "log" & rp$tipo_nll == "nll", ][1, ]
  vo <- leer("ensayo2_rcpp_vs_okuyama.csv")
  vo <- vo[vo$dataset == lab & vo$escala == "log" & vo$tipo_nll == "nll", ][1, ]
  pp <- leer("ensayo2_rcpp_perfiles.csv")
  pp <- pp[pp$dataset == lab & pp$escala == "log", ]

  # Perfil punto a punto (z comun)
  m <- merge(data.frame(z = round(prof$z, 2), nll_paquete = prof$nll,
                        nll_reeval_paquete = prof$nll_reeval),
             data.frame(z = round(pp$z, 2), nll_proyecto = pp$nll,
                        nll_reeval_proyecto = pp$nll_reeval), by = "z")
  m$dif_nll <- m$nll_paquete - m$nll_proyecto
  m$dif_nll_reeval <- m$nll_reeval_paquete - m$nll_reeval_proyecto
  write.csv(m, file.path("tablas", paste0("paquete_perfil_comparacion_", lab, SUF, ".csv")),
            row.names = FALSE)

  # NLL del perfil del paquete en el z_hat del proyecto
  nll_en_zproj <- {
    i <- which(abs(prof$z - rp$z_hat) < 1e-8)
    if (length(i)) prof$nll[i[1]] else NA_real_
  }
  d_en_zproj <- nll_en_zproj - ci$min_nll

  ic_low_proj  <- rp$ic_interp_low
  ic_high_proj <- rp$ic_interp_high
  cmp_ic <- function(pkg, proj, tol) {
    if (is.na(pkg) && is.na(proj)) return("OK (ambos abiertos)")
    if (is.na(pkg) || is.na(proj)) return("REVISAR (uno abierto)")
    if (abs(pkg - proj) <= tol) "OK" else "REVISAR"
  }
  par <- res$par
  aic_bb <- rp$aic_bb
  chk <- function(x) if (isTRUE(x)) "OK" else "REVISAR"

  out <- data.frame(
    dataset = lab,
    z_hat_paquete = ci$z_hat, z_hat_proyecto = rp$z_hat, z_hat_okuyama = vo$z_hat_okuyama,
    z_hat_ok = chk(abs(ci$z_hat - rp$z_hat) <= TOL_Z || (is.finite(d_en_zproj) && d_en_zproj <= TOL_NLL)),
    nll_min_paquete = ci$min_nll, nll_min_proyecto = rp$nll_min,
    nll_paquete_en_zhat_proyecto_menos_min = d_en_zproj,
    ic_low_paquete = ci$z_low, ic_low_proyecto = ic_low_proj,
    ic_low_ok = cmp_ic(ci$z_low, ic_low_proj, TOL_Z),
    ic_high_paquete = ci$z_high, ic_high_proyecto = ic_high_proj,
    ic_high_censurado_paquete = ci$high_censored,
    ic_high_ok = cmp_ic(ci$z_high, ic_high_proj, TOL_Z),
    ic_okuyama = vo$ic_okuyama,
    h_paquete = par[["h"]], h_proyecto = rp$h,
    h_ok = chk(abs(par[["h"]] / rp$h - 1) <= TOL_REL),
    k_paquete = par[["k"]], k_proyecto = rp$k,
    a_paquete = par[["a"]], a_proyecto = rp$a,
    s_paquete = par[["s"]], s_proyecto = rp$s,
    aic_paquete = res$aic$AIC_full, aic_proyecto = rp$aic, aic_okuyama = vo$aic_okuyama,
    delta_aic_bb_menos_sim_paquete = aic_bb - res$aic$AIC_full,
    delta_aic_bb_menos_sim_proyecto = rp$delta_aic_bb_menos_sim,
    delta_aic_bb_menos_sim_okuyama = vo$delta_aic_bb_sim_okuyama,
    perfil_dif_nll_media = mean(m$dif_nll), perfil_dif_nll_sd = sd(m$dif_nll),
    perfil_dif_nll_maxabs = max(abs(m$dif_nll)),
    perfil_pct_dentro_1 = mean(abs(m$dif_nll) <= 1) * 100,
    n_puntos_comunes = nrow(m), n_puntos_paquete = nrow(prof),
    extendido = res$extended, refinado = res$refined,
    horas = res$meta$horas, version_paquete = res$meta$version_paquete,
    origen_paquete = if (is.null(res$meta$origen_paquete)) NA_character_ else res$meta$origen_paquete,
    stringsAsFactors = FALSE)
  write.csv(out, file.path("tablas", paste0("paquete_vs_proyecto_", lab, SUF, ".csv")),
            row.names = FALSE)

  cat(sprintf(paste0(
    "\n[D] %s  (paquete | proyecto | Okuyama)\n",
    "    z_hat : %.2f | %.2f | %.2f   -> %s\n",
    "    IC low: %s | %s   -> %s\n",
    "    IC high: %s | %s   -> %s   (Okuyama %s)\n",
    "    h     : %.3f | %.3f   -> %s\n",
    "    AIC   : %.3f | %.3f | %.3f\n",
    "    dAIC (BB - sim): %.2f | %.2f | %.2f\n",
    "    Perfil punto a punto (%d puntos): dif NLL media %.2f, sd %.2f, max %.2f; %.0f%% dentro de +-1\n"),
    lab, ci$z_hat, rp$z_hat, vo$z_hat_okuyama, out$z_hat_ok,
    format(round(ci$z_low, 2)), format(ic_low_proj), out$ic_low_ok,
    if (is.na(ci$z_high)) "abierto" else format(round(ci$z_high, 2)),
    if (is.na(ic_high_proj)) "abierto" else format(round(ic_high_proj, 2)), out$ic_high_ok, vo$ic_okuyama,
    par[["h"]], rp$h, out$h_ok,
    out$aic_paquete, out$aic_proyecto, out$aic_okuyama,
    out$delta_aic_bb_menos_sim_paquete, out$delta_aic_bb_menos_sim_proyecto,
    out$delta_aic_bb_menos_sim_okuyama,
    nrow(m), out$perfil_dif_nll_media, out$perfil_dif_nll_sd, out$perfil_dif_nll_maxabs,
    out$perfil_pct_dentro_1))
  invisible(out)
}

# ============================================================
# Ejecucion
# ============================================================
cat("\n== Validacion del paquete funresMech", v_inst, "==\n")
if (!SOLO_COMPARAR) {
  writeLines(capture.output(sessionInfo()),
             file.path(dir_res, paste0("sessionInfo_", paste(LABS, collapse = "_"), ".txt")))
  for (lab in LABS) {
    fase_A(lab)
    fase_B(lab)
  }
  for (lab in LABS) fase_C(lab)
}
for (lab in LABS) fase_D(lab)
cat("\nTerminado. Resultados en tablas/paquete_*", SUF, ".csv y ", dir_res, "\n", sep = "")
