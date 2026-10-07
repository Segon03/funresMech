# ============================================================
# R/03_ensayo1_resumen.R
#
# Lee todos los .rds de results/ensayo1/ (generados por
# 02_ensayo1_simulacion.R) y calcula, por nivel de rho:
#   - sesgo medio de z (z_hat - 1)
#   - cobertura del IC95% (proporcion de replicas cuyo IC contiene z=1)
#   - medias de a_hat, h_hat, k_hat, s_hat
#
# CORRECCION (25-sep-2026): el IC95% ya NO se lee de r$ci_low /
# r$ci_high (calculados en 02 con el umbral sin "/2" de la app
# v1.0.4, ver fun/profile_utils.R), sino que se RECALCULA aca a
# partir del perfil guardado (r$z_grid, r$nll_local) con el umbral
# correcto NLL_min + qchisq(.95,1)/2. No hace falta volver a correr
# 02_ensayo1_simulacion.R. Los IC viejos quedan en las columnas
# ci_low_app_v104 / ci_high_app_v104 para comparar.
#
# Ademas, un limite que el perfil no alcanza a cruzar dentro de su
# ventana local se marca como censurado (antes quedaba NA y la
# replica se descartaba en silencio de la cobertura con na.rm=TRUE).
# Para la cobertura, un limite censurado del lado izquierdo cuenta
# como "<= min(grilla)" y del lado derecho como ">= max(grilla)".
#
# AVISO: el nll de los ajustes con rho > 0 es enorme (miles), senal
# de que muchas observaciones tienen probabilidad simulada 0 (piso
# .Machine$double.xmin -> +708 por observacion). Ver el chat del
# 25-sep-2026: los resultados del Ensayo 1 con rho > 0 no son
# interpretables tal como estan y el ensayo se va a redisenar.
# ============================================================

library(dplyr)
source(file.path("fun", "profile_utils.R"))

archivos <- list.files(file.path("results", "ensayo1"), pattern = "\\.rds$",
                        full.names = TRUE)
if (length(archivos) == 0) {
  stop("No hay resultados en results/ensayo1/. Correr 02_ensayo1_simulacion.R primero.")
}

filas <- lapply(archivos, function(f) {
  r <- readRDS(f)
  tiene_perfil <- length(r$z_grid) > 1 && all(is.finite(r$nll_local))
  if (tiene_perfil) {
    ci <- ci_from_profile_okuyama(r$z_grid, r$nll_local)   # umbral correcto (+1.92)
    ci_low  <- if (ci$low_censurado)  min(r$z_grid) else ci$z_low
    ci_high <- if (ci$high_censurado) max(r$z_grid) else ci$z_high
    low_c <- ci$low_censurado; high_c <- ci$high_censurado
    contiene_1 <- (low_c || ci_low <= 1) && (high_c || ci_high >= 1)
  } else {
    ci_low <- ci_high <- NA_real_; low_c <- high_c <- NA; contiene_1 <- NA
  }
  data.frame(
    rho     = r$rho,
    rep_id  = r$rep_id,
    a_hat   = unname(r$fit$par["a"]),
    h_hat   = unname(r$fit$par["h"]),
    z_hat   = unname(r$fit$par["z"]),
    k_hat   = unname(r$fit$par["k"]),
    s_hat   = unname(r$fit$par["s"]),
    nll_fit = r$fit$nll,
    ci_low  = ci_low,
    ci_high = ci_high,
    ci_low_censurado  = low_c,
    ci_high_censurado = high_c,
    ci_contiene_1     = contiene_1,
    ci_low_app_v104  = r$ci_low,
    ci_high_app_v104 = r$ci_high,
    elapsed_sec = r$elapsed_sec
  )
})
res <- bind_rows(filas)

cat(sprintf("Total de ajustes leidos: %d (esperados: 60)\n", nrow(res)))
if (nrow(res) < 60) {
  cat("AVISO: todavia no estan los 60 ajustes. Este resumen es parcial;\n",
      "correr 02_ensayo1_simulacion.R de nuevo para completar (retoma solo).\n", sep = "")
}

resumen <- res %>%
  group_by(rho) %>%
  summarise(
    n                = n(),
    sesgo_z          = mean(z_hat - 1),
    de_z             = sd(z_hat),
    n_con_ic         = sum(!is.na(ci_contiene_1)),
    cobertura_z      = mean(ci_contiene_1, na.rm = TRUE),
    n_ic_censurado   = sum(ci_low_censurado | ci_high_censurado, na.rm = TRUE),
    mediana_nll_fit  = median(nll_fit),
    media_a          = mean(a_hat),
    media_h          = mean(h_hat),
    media_k          = mean(k_hat),
    media_s          = mean(s_hat),
    tiempo_medio_min = mean(elapsed_sec) / 60,
    .groups = "drop"
  )

print(resumen)

dir.create("tablas", showWarnings = FALSE)
write.csv(res,     file.path("tablas", "ensayo1_ajustes_individuales.csv"), row.names = FALSE)
write.csv(resumen, file.path("tablas", "ensayo1_resumen_por_rho.csv"), row.names = FALSE)

cat("\nGuardado en tablas/ensayo1_resumen_por_rho.csv y tablas/ensayo1_ajustes_individuales.csv\n")
