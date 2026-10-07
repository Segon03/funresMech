# ============================================================
# R/04d_ensayo2_corr_OKU.R  (continua la Seccion 4 del plan)
#
# Reproduce el ajuste de D2 y D3 usando LOS MISMOS PARAMETROS DE
# OPTIMIZACION QUE OKUYAMA (2026), tal como aparecen en su propio
# material suplementario (papers/Okuyama_Supp_Mat/, Script S3 y
# S4: jen70148-sup-0006-datas3, jen70148-sup-0007-datas4):
#
#     T       = 24     (horas reales de exposicion del protocolo
#                        de Kishani Farahani & Goldansaz 2013;
#                        antes T_exp=1 en 04_ensayo2_datos_reales.R)
#     itermax = 200    (Okuyama: DEoptim.control(itermax = 200);
#                        antes 50)
#     NP      = 40     (sin cambios - ya coincidia con Okuyama,
#                        que no fija NP y usa el default 10*n_par)
#     reltol  : NO se fija (Okuyama tampoco lo hace en su Script
#               S3) -> se usa el default de DEoptim (~1.49e-8),
#               mucho mas estricto que el 1e-2 usado antes (ese
#               1e-2 permitia que DEoptim cortara "por las dudas"
#               apenas la poblacion dejaba de mejorar >1%, con
#               una verosimilitud ya de por si ruidosa por ser
#               simulada).
#
#     n_sim   = 3000   AJUSTADO A PROPOSITO, NO es el 10000 de
#               Okuyama (ver benchmark real corrido en la maquina
#               del usuario: n_sim=10000+itermax=200 daba ~114 h
#               de pared con 23 workers, ~5 dias). Se decidio
#               mantener n_sim=3000 (el mismo que ya se usaba en
#               el ajuste principal del Ensayo 2 original) porque
#               el diagnostico del bug de D2 (ver mas abajo) senalo
#               a itermax=50+reltol=1e-2+cotas h<=0.5 como la causa
#               real de la no convergencia - NO al n_sim. Mantener
#               itermax=200 (la correccion que si importa) mientras
#               se baja n_sim de 10000 a 3000 reduce el tiempo
#               estimado de ~114 h a ~34 h (~1.4 dias) de pared.
#               Si mas adelante se quiere maxima fidelidad a
#               Okuyama, subir N_SIM de nuevo a 10000 mas abajo.
#
# POR QUE HACE FALTA UN SCRIPT NUEVO (no alcanza con cambiar
# constantes en 04_ensayo2_datos_reales.R):
#
# 1) funresMech:::fit_full() (el que arma el ajuste PRINCIPAL con
#    z libre) trae las cotas de busqueda de DEoptim CODEADAS
#    ADENTRO del paquete: h en [0.001, 0.5]. Bajo T=24, Okuyama
#    encontro h=0.606 (D2) y h=0.626 (D3) - AMBOS por encima del
#    limite superior 0.5 que fit_full() nunca deja cruzar. Por
#    eso este script NO llama a fit_full(): arma su propia
#    llamada a DEoptim (mismo patron que fit_full() usa por
#    dentro - ver R/fit_full.R del paquete, confirmado leyendo el
#    codigo fuente en github.com/Segon03/funresMech) pero con
#    cotas mas anchas (0-5 para a,h,k, 0-3 para s), en el mismo
#    espiritu que las que uso Okuyama en su Script S3 (0 a 5 para
#    los 4 parametros).
#
# 2) Se verifico (ver chat) que el ajuste PRINCIPAL de D2 con la
#    configuracion anterior (T=1, itermax=50, reltol=1e-2,
#    n_sim=3000) quedo con NLL=183.667 - PEOR que el propio punto
#    de SU MISMO perfil de z en z=1.10 (NLL=175.07). Eso es
#    logicamente imposible si DEoptim hubiera convergido al
#    optimo real: el perfil en z=1.10 optimiza (a,h,k,s) con z
#    FIJO, un subespacio restringido del mismo espacio de 5
#    parametros que el ajuste principal explora libre - un
#    minimo libre NUNCA puede ser peor que uno restringido. Por
#    eso, ademas de las 4 correcciones de arriba, el ajuste
#    PRINCIPAL (el mas sensible y el que se cita en el paper) se
#    corre por TRIPLICADO (semillas independientes) y se toma el
#    de menor NLL.
#
# COSTO COMPUTACIONAL - LEER ANTES DE DEJARLO CORRIENDO:
# Con n_sim=3000 e itermax=200 (config. actual de este script) el
# benchmark real en la maquina de 24 nucleos dio ~34 h de pared
# (23 workers) - bajado a proposito desde las ~114 h que daba con
# n_sim=10000 (ver nota de arriba). Este script corre primero un
# benchmark de UNA sola evaluacion de verosimilitud para
# reconfirmar el tiempo total ANTES de lanzar todas las tareas -
# revisar ese numero en la consola/log antes de dejarlo
# desatendido muchas horas. Si el tiempo estimado sigue siendo
# excesivo, la palanca mas directa es bajar N_SIM aun mas
# (variable marcada abajo) - no hace falta tocar nada mas del
# script.
#
# Checkpointing por archivo (si se corta, correr de nuevo retoma
# solo lo que falta) + paralelizacion con future_lapply desde el
# arranque (evita el error de 04b v1, que corria todo secuencial
# en un solo proceso).
#
# IMPORTANTE: correr desde la carpeta del proyecto (usar el
# .Rproj), con 00_benchmark.R y 01_extract_forage_d2_d3.R ya
# corridos antes (se reutiliza data_clean/D2.rds y D3.rds - la
# extraccion de datos no depende de T, asi que no hace falta
# repetirla).
# ============================================================

library(funresMech)
library(future)
library(future.apply)
source(file.path("fun", "betabinom_model.R"))

dir.create(file.path("results", "ensayo2_corr_OKU"), showWarnings = FALSE, recursive = TRUE)
dir.create("tablas", showWarnings = FALSE)

# ---- Configuracion: EXACTA a Okuyama (2026), Script S3/S4 --------------
T_exp   <- 24     # horas reales de exposicion (antes T_exp=1)
itermax <- 200    # Okuyama: DEoptim.control(itermax = 200) (antes 50)
NP      <- 40     # sin cambios
N_SIM   <- 3000   # Ajustado a proposito (Okuyama usa 10000; ver nota
                   # arriba sobre el costo real medido: ~114h con 10000
                   # vs. ~34h con 3000). Subir a 10000 para maxima
                   # fidelidad si el tiempo lo permite mas adelante.
# reltol: no se fija a proposito (ver nota arriba) -> default de DEoptim.

n_restarts_principal <- 3   # solo para el ajuste principal (z libre),
                             # que fue el que quedo mal convergido antes

# Red de seguridad GENERICA (no depende de que el diagnostico de h de
# arriba sea la unica causa posible): si UNA sola evaluacion de la
# verosimilitud tarda mas de este limite, se aborta y se le devuelve a
# DEoptim una penalizacion enorme en vez de dejar que el worker quede
# atascado ahi por horas. El benchmark tipico (parametros razonables)
# dio ~7-8 s por evaluacion a n_sim=3000; 180 s es ~23x ese valor, un
# margen amplio para variacion normal sin dejar pasar un cuelgue real.
LIMITE_SEG_POR_EVAL <- 180

evaluar_con_limite <- function(expr_fn, limite_seg = LIMITE_SEG_POR_EVAL) {
  setTimeLimit(elapsed = limite_seg, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE))
  tryCatch(expr_fn(), error = function(e) 1e10)
}

# Cotas de busqueda: mismo espiritu que Okuyama Script S3 (0 a 5 para
# a, h, k, s). k=0 exacto rompe negloglik_fixed_z ("if (k<=0) return
# 1e10"), por eso el limite inferior es un epsilon chico, no 0 literal.
#
# OJO CON h: el piso de h NO puede ser un epsilon tan chico como el de
# a/k. En simulate_trial() (motor R puro de funresMech) cada vuelta del
# while(t<T) suma ~1/(a*x^z) (tiempo de busqueda) + h (tiempo de
# manejo). Si a*x^z es grande (va a pasar seguido, ya que DEoptim
# explora a en todo [0,5]) el tiempo de busqueda se vuelve
# insignificante y el numero de vueltas del loop queda gobernado por
# T/h. Con T=24 y h cerca de 0, eso puede disparar a decenas de
# millones de iteraciones por UN solo ensayo simulado (y hacen falta
# n_sim=3000 de esos por evaluacion) - en R puro (no C++, a diferencia
# del simulador de Okuyama) eso puede tardar horas para una sola
# evaluacion, y alcanza con que unas pocas de las 40 combinaciones
# iniciales de cada generacion caigan ahi para que un worker entero
# quede atascado. Por eso h NO baja de 0.01 (peor caso: T/h = 2400
# iteraciones, sin problema) aunque a y k si puedan usar un piso casi
# cero.
lower_ok <- c(a = 1e-6, h = 0.01, k = 1e-6, s = 0)
upper_ok <- c(a = 5,    h = 5,    k = 5,    s = 3)
z_lower  <- 0.01
z_upper  <- 6.0    # cubre holgadamente el rango ya usado en el perfil

# Grillas de perfil: se reutilizan las que ya se verificaron
# suficientes en el ensayo anterior (D3 necesito extenderse hasta
# 5.9 para que el IC cerrara del lado superior - ver 04b/04c).
z_grid <- list(
  D2 = seq(0.5, 3.0, by = 0.2),
  D3 = seq(0.5, 5.9, by = 0.2)
)

labs <- c("D2", "D3")
datasets <- setNames(lapply(labs, function(lab) {
  data_path <- file.path("data_clean", paste0(lab, ".rds"))
  if (!file.exists(data_path)) {
    stop("No se encontro ", data_path, ". Correr 01_extract_forage_d2_d3.R primero.")
  }
  readRDS(data_path)
}), labs)

# ---- Benchmark rapido: 1 evaluacion de negloglik_fixed_z a n_sim=10000 --
cat(sprintf("Benchmark: 1 evaluacion de la verosimilitud a n_sim = %d ...\n", N_SIM))
t0 <- proc.time()
invisible(funresMech:::negloglik_fixed_z(
  c(a = 0.02, h = 0.4, k = 1, s = 0.3), z_fixed = 1.5,
  data_spp = datasets[["D2"]], T = T_exp, n_sim = N_SIM
))
t_eval <- (proc.time() - t0)[3]

n_tareas_perfil    <- sum(lengths(z_grid))
n_tareas_principal <- length(labs) * n_restarts_principal
n_eval_deoptim     <- NP * (itermax + 1)  # aprox. evaluaciones por corrida DEoptim
horas_secuencial   <- (n_tareas_perfil + n_tareas_principal) * n_eval_deoptim * t_eval / 3600

cat(sprintf(
  paste0("1 evaluacion tardo %.2f s.\n",
         "Con NP=%d, itermax=%d (~%d evaluaciones por corrida de DEoptim),\n",
         "%d tareas de perfil + %d de ajuste principal (%d datasets x %d intentos):\n",
         "ESTIMADO SECUENCIAL (sin paralelizar): %.1f horas-nucleo.\n",
         "Este script SI paraleliza (future_lapply) - dividir aprox. por la\n",
         "cantidad de workers de abajo para una estimacion de tiempo de pared.\n"),
  t_eval, NP, itermax, n_eval_deoptim,
  n_tareas_perfil, n_tareas_principal, length(labs), n_restarts_principal,
  horas_secuencial
))

# ---- Lista de tareas ----------------------------------------------------
tareas <- list()
for (lab in labs) {
  for (i in seq_len(n_restarts_principal)) {
    tareas[[length(tareas) + 1]] <- list(type = "principal", lab = lab, intento = i)
  }
  for (z in z_grid[[lab]]) {
    tareas[[length(tareas) + 1]] <- list(type = "perfil", lab = lab, z = z)
  }
}
cat(sprintf("Ensayo 2 (corr_OKU): %d tareas a correr.\n", length(tareas)))

n_cores   <- parallel::detectCores()
n_workers <- max(1, n_cores - 1)
cat(sprintf("Paralelizando en %d workers (de %d cores detectados).\n", n_workers, n_cores))
plan(multisession, workers = n_workers)

correr_tarea <- function(tarea) {
  library(funresMech)
  lab <- tarea$lab
  data_spp <- datasets[[lab]]

  if (tarea$type == "principal") {
    out_file <- file.path("results", "ensayo2_corr_OKU",
                           sprintf("principal_%s_corrOKU_intento%d.rds", lab, tarea$intento))
    if (file.exists(out_file)) {
      return(sprintf("%s intento %d: ya existe, se salta", lab, tarea$intento))
    }

    res <- DEoptim::DEoptim(
      fn = function(par) {
        a <- par[1]; h <- par[2]; z <- par[3]; k <- par[4]; s <- par[5]
        evaluar_con_limite(function() {
          funresMech:::negloglik_fixed_z(c(a, h, k, s), z_fixed = z,
                                          data_spp = data_spp, T = T_exp,
                                          n_sim = N_SIM)
        })
      },
      lower = c(lower_ok["a"], lower_ok["h"], z_lower, lower_ok["k"], lower_ok["s"]),
      upper = c(upper_ok["a"], upper_ok["h"], z_upper, upper_ok["k"], upper_ok["s"]),
      control = DEoptim::DEoptim.control(itermax = itermax, NP = NP, trace = FALSE)
    )
    par <- res$optim$bestmem
    names(par) <- c("a", "h", "z", "k", "s")
    saveRDS(list(par = par, nll = res$optim$bestval), out_file)
    return(sprintf("%s intento %d: nll=%.3f z=%.3f h=%.3f",
                    lab, tarea$intento, res$optim$bestval, par["z"], par["h"]))
  }

  # tarea$type == "perfil"
  out_file <- file.path("results", "ensayo2_corr_OKU",
                         sprintf("perfil_%s_corrOKU_z%.2f.rds", lab, tarea$z))
  if (file.exists(out_file)) {
    return(sprintf("%s z=%.2f: ya existe, se salta", lab, tarea$z))
  }
  res <- DEoptim::DEoptim(
    fn = function(par) evaluar_con_limite(function() {
      funresMech:::negloglik_fixed_z(par, z_fixed = tarea$z,
                                      data_spp = data_spp, T = T_exp,
                                      n_sim = N_SIM)
    }),
    lower = lower_ok, upper = upper_ok,
    control = DEoptim::DEoptim.control(itermax = itermax, NP = NP, trace = FALSE)
  )
  saveRDS(list(z = tarea$z, nll = res$optim$bestval), out_file)
  sprintf("%s z=%.2f: nll=%.3f", lab, tarea$z, res$optim$bestval)
}

mensajes <- future_lapply(tareas, correr_tarea, future.seed = TRUE)
invisible(lapply(mensajes, cat, "\n"))

# ---- Beta-binomial bajo T=24 (rapido, no hace falta paralelizar) -------
# fun/betabinom_model.R ya acepta T como argumento - se reutiliza tal
# cual, solo cambia el T que se le pasa. Esto da a/h directamente
# comparables (misma escala) con los a/h que publica Okuyama.
for (lab in labs) {
  out_file <- file.path("results", "ensayo2_corr_OKU", sprintf("bb_%s_corrOKU.rds", lab))
  if (!file.exists(out_file)) {
    fit_bb <- fit_betabinom(datasets[[lab]], T = T_exp)
    saveRDS(fit_bb, out_file)
    cat(sprintf("%s beta-binomial (T=24): nll=%.3f convergence=%s\n",
                lab, fit_bb$value, fit_bb$convergence))
  } else {
    cat(sprintf("%s beta-binomial (T=24): ya existe, se salta\n", lab))
  }
}

cat("\nListo. Correr 04e_ensayo2_corr_OKU_reensamblar.R para armar la tabla final",
    "y la comparacion directa contra la Tabla S1 de Okuyama (2026).\n")
