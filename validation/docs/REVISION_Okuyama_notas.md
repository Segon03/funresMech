# Revisión en respuesta a Okuyama (2026-10-13): notas de trabajo

Plan completo: `PLAN_respuesta_Okuyama.md` (carpeta `Okuyama_rev` y proyecto).

## Flujo de trabajo

* Línea base congelada: funresMech 1.1.0 tal como se publicó (CRAN 2026-10-06,
  GitHub commit 52a779d). Copia intacta en `Okuyama_rev/funresMech_1.1.0_OR`.
* Copia de trabajo con historial git: `Okuyama_rev/funresMech_dev`,
  rama `revision-okuyama`, etiqueta `baseline-1.1.0` en el commit base.
* El código se edita y se prueba en el entorno de Claude (R 4.3.3, Linux) y
  cada paso se aplica a `funresMech_dev` como commit. Las corridas largas
  (validación, banco de optimizadores) se ejecutan en la PC del autor (Windows,
  23 núcleos). Nada se sube a GitHub ni a CRAN hasta la Fase 9.

## Puerta de regresión (Fase 0.3)

`tests/testthat/test-regression.R` + `tests/testthat/fixtures/golden_1.1.0.rds`.

* Ancla contra el código de Okuyama (2026): D3, z = 2.5, par = (0.012, 0.626,
  0.902, 0.542), 3000 ensayos simulados, semilla 88, sin corte por saturación
  → NLL = 132.04994217 (reproducido con el paquete 1.1.0 instalado).
* Conteos simulados exactos del motor (5 casos × con/sin corte) y NLL exactas
  (5 casos D2/D3) de la versión 1.1.0.
* El modelo por defecto (búsqueda Gamma, manipulación Lognormal) debe
  reproducir estos valores bit a bit en todas las versiones siguientes.
* El valor 7944.851530 (D2) de la validación compara con la implementación en
  R puro de 1.0.4, cuya `s` era el desvío en escala log; no es reproducible con
  la definición actual de `s` y por eso no forma parte de la puerta.

## Inconsistencias registradas (a resolver en las fases indicadas)

| ID | Inconsistencia | Fase |
|---|---|---|
| H3 | El IC de z y z_hat usan `nll` (mínimo ruidoso del optimizador), no `nll_reeval`. | 3.1 |
| H4 | `n_reeval` = 3 por defecto en el código; el S1 Text dice que el NLL se reevaluó cinco veces. | 3.1 |
| H6 | El paquete compara AIC de z libre vs z = 1; el manuscrito dice que compara con beta-binomial. | 3.11 |
| H8 | Grilla de z por defecto en la app: 0.5–5 paso 0.25; manuscrito: 0.5–6.0 paso 0.1. | 3.8 / 8 |
| H1–H2 | Cotas absolutas (h ≤ 5, λ_ref en [1e-3, 1e3]) y T rotulado en horas: dependencia de la unidad de tiempo. | 2 |
| H14 | `reltol` de DEoptim no tenía efecto en 1.1.0: DEoptim lo evalúa cada `steptol` generaciones y `steptol` vale `itermax` por defecto, así que nunca había parada temprana. | resuelto en `funres_control()` (reltol y steptol explícitos; por defecto sin parada temprana) |
| H15 | Las cotas del modelo beta-binomial del screening (a, h <= 5) también son absolutas y dependen de la unidad de T. | 2 |

## Decisiones del autor (2026-10-13)

* D1: se exporta una API pública además de `run_app()`.
* D2: distribuciones en tres niveles: familias en C++, muestreadores en R y
  muestreadores en C++ (puntero externo).
* D3: modo depredador: verificar primero. Verificado: el Apéndice S1 y el
  Script S2 de Okuyama (2026) definen el proceso (presa consumida retirada,
  dH/dt = -P a H^z / (1 + a h H^z)); el propio apéndice indica que la elección
  depende del sistema (un parasitoide que no revisita hospedadores
  parasitados equivale al proceso depredador). Se implementa como
  `process = "predator"`.
* D4: versión 1.2.0.
* D5 (optimizadores), decidido tras la Fase 4: por defecto DEoptim con números
  aleatorios comunes y pulido Nelder-Mead (200 evaluaciones); z por perfil;
  la estimación conjunta queda como opción experimental (subcobertura 67%).
  Detalle en `FASE4_resultados.md`.
* D6 (presupuesto de cómputo): PC de 23-24 núcleos, hasta 4 días.

## Fase 1 (estado)

* Motor generalizado (`src/motor.cpp`): familias de búsqueda Gamma, Lognormal,
  Weibull, R y C++ del usuario (media 1/lambda, CV^2 = 1/k); manipulación
  Lognormal, Gamma, Weibull, constante, R y C++ (media h, desvío s); proceso
  parasitoide o depredador. Ruta rápida del modelo por defecto idéntica bit
  a bit a 1.1.0.
* API: `funres_model()`, `funres_sampler()`, `test_sampler()`,
  `simulate_response()`, `funres_control()`, `fit_funres()`.
* Parametrización adimensional (tiempo en unidades de T) y cotas por defecto
  independientes de la unidad (adelanto de la Fase 2).
* Perfil e IC sobre la media de las reevaluaciones (`nll`), con su error
  estándar (`nll_se`); el mínimo del optimizador queda como `nll_opt`
  (adelanto de la Fase 3.1).

## Fases 2 a 4 (estado al 2026-10-13)

* Fase 2: screening beta-binomial en unidades de T (invariante a la unidad);
  ampliación automática de cotas tocadas en el óptimo (x10, hasta 2 veces);
  DESCRIPTION y ayuda sin "Okuyama model"; versión 1.2.0.
* Fase 3: `reliability()` (semáforo y acción por chequeo), `verify_fit()`
  (corridas independientes, más esfuerzo, cotas más anchas), diagnóstico de
  piso de probabilidad (`n_floor`), AIC frente a beta-binomial con la misma
  media (parasitoide: Ecs. 1-2; depredador: ecuación del depredador aleatorio,
  Script S2). El beta-binomial reproduce la NLL implícita en la Tabla S1 de
  Okuyama (D3 132.68 vs 132.76; D2 169.46 vs 169.44).
* Fase 4 (código): `z_method = "joint"` con perfil adaptativo; métodos
  `DEoptim`, `CMAES`, `CRS2`; números aleatorios comunes (`crn`) y pulido
  Nelder-Mead (`polish`). Falta la comparación (script en la PC del autor).
* Umbrales de los chequeos: provisorios (`.rel_thresholds()`), a calibrar con
  la validación ampliada (Fase 5).

## Próximo paso en la PC del autor

1. `devtools::install()` y `devtools::test()` en `funresMech_dev` (Windows).
2. `source("validation/R/10_piloto_costos.R")` y enviar
   `validation/tablas/piloto_costos.csv`.

## Fase 6: interfaz Shiny (2026-10-13)

* Reescrita con bslib (Bootstrap 5); se quitan shinythemes, shinyBS, dplyr,
  magrittr y kableExtra. Tema sobrio (azul pizarra #1C2B3A, acento #1F5F8B,
  estados verde/ocre/rojo apagados), fuente del sistema a 14 px.
* Pestañas: Data (mapeo de columnas, grupos, screening con fila del archivo),
  Model and fit (unidad de tiempo, proceso, familias de búsqueda y
  manipulación, muestreadores propios en R con prueba, k/s/z fijos, z por
  perfil/conjunto/fijo, presets, método, crn, pulido, semilla, paralelo,
  cotas por defecto en unidades del usuario o propias, estimación del tiempo),
  Reliability (franja de estado, chequeos con acción, verificación), Results
  (estimaciones con unidades, IC de z, AIC vs z = 1 y vs beta-binomial,
  perfil con error MC, respuesta simulada, distribuciones; figuras PDF/PNG),
  Model comparison (tabla de AIC entre especificaciones, aviso si difieren los
  ensayos), Export (CSV, RDS, script R reproducible, reporte HTML), Help.
* El script R y el reporte usan la configuración con que se hizo el ajuste.
* Chequeo nuevo: perfil de z irregular (mínimos locales adicionales) = falla.
* Pruebas: tests/testthat/test-app.R (testServer del flujo completo).
* R CMD check (Linux): sin errores; la única advertencia es de la
  configuración regional del contenedor.

## Piloto de costos (PC de 8 núcleos, 2026-10-13)

* Una evaluación de la NLL cuesta ~0.09 s por cada 1000 simulaciones por
  densidad (D3: 0.19 s con 2000; 0.92 s con 10000); lineal en n_sim.
* Un punto del perfil: exploratory 6.5 min, standard ~31 min, final ~2 h
  (un núcleo). CMA-ES y CRS2 cuestan lo mismo que DEoptim con igual
  presupuesto de evaluaciones (346-413 s).
* Decisión (D6): Fases 4 y 5 en la PC de 23-24 núcleos, hasta 4 días de
  cómputo. Fase 4 (11_comparacion_optimizadores.R): ~7 h con 22 núcleos.
