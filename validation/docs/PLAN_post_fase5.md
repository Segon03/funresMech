# Plan de acción posterior a la Fase 5 (funresMech 1.2.0)

Fecha: 2026-10-16 · Reemplaza el orden de trabajo del plan del 2026-10-13 a partir de ahora.
Carpeta: `Okuyama_rev/funresMech_dev`, rama `revision-okuyama`. Nada se sube a GitHub ni a CRAN hasta el Paso 6.

---

## Estado al 2026-10-16

| Fase del plan original | Estado |
|---|---|
| 0. Línea base y puerta de regresión | Hecha. El modelo por defecto reproduce la 1.1.0 bit a bit (D3, NLL 132.04994217). |
| 1. Motor con distribuciones elegibles | Hecha. Familias en C++, muestreadores propios en R y en C++, parámetros fijos y modo depredador. |
| 2. Unidades, cotas y nombres | Hecha. Tiempo en unidades de T, cotas por defecto independientes de la unidad, ampliación automática de cotas y título sin "Okuyama model". **Falta** un test de invarianza de la unidad del ajuste completo (solo existe para el screening). |
| 3. Diagnósticos de fiabilidad | Hecha en código. `reliability()` y `verify_fit()`, con acción por chequeo. **Umbrales provisorios**: se calibran en el Paso 3. |
| 4. Optimizadores y z conjunta | Hecha. Por defecto, DEoptim con números aleatorios comunes y pulido; la z conjunta queda como opción experimental. No se evaluó ABC (ver D7). |
| 5. Validación ampliada | **Corrida terminada** el 2026-10-16: 679 de 679 trabajos, sin errores, 14 tablas. Falta verificar e interpretar. |
| 6. Interfaz Shiny | Hecha. Solo con `testServer`, sin `shinytest2`. |
| 7. Documentación | Pendiente. NEWS sin la 1.2.0, README viejo, sin viñetas. |
| 8. Manuscrito | Pendiente. |
| 9. Liberación y carta a Okuyama | Pendiente. |

---

## Paso 1. Verificar que los resultados de la Fase 5 son correctos

Antes de interpretar, se comprueba que las tablas dicen lo que deben decir. Cada chequeo deja un resultado (pasa o no pasa) en `FASE5_resultados.md`.

### 1.1 Integridad

- 679 archivos `.rds` y 14 tablas. Ya comprobado.
- Conteos esperados por brazo: A = 45 × 8, B = 160, C = 36 × 5, D = 36 × 2 más 12, E = 360 × 13 y F = 3. Ya comprobado.
- Explicar cada NA de las tablas: columnas vacías, límites abiertos y columnas que no aplican. Ya hecho en forma preliminar; queda dejarlo documentado.

### 1.2 Coherencia interna

- Los mismos conjuntos de datos se usan en A, C y D. Hay que comprobar que la fila de la 1.2.0 en D y la variante "exploratory" en C coinciden exactamente con la fila de A (z_hat, NLL y `nll_prec_true`).
- La NLL precisa en la verdad (`nll_prec_true`) debe coincidir entre las variantes de C, salvo ruido de Monte Carlo. Con 5 × 20 000 simulaciones se espera un desvío de ~0.1–0.2.
- El déficit es la NLL precisa del estimado menos la NLL precisa en la verdad. Debería ser ≤ 0 en la mayoría de los ajustes, porque el máximo de verosimilitud no puede ser peor que la verdad. Hay que contar la proporción de déficit mayor que 1, que define la falla de optimización.
- Cada intervalo debe contener a z_hat y respetar z_low ≤ z_hat ≤ z_high. Los intervalos abiertos tienen que estar marcados como tales.
- En los perfiles de F, el mínimo del perfil tiene que coincidir con z_hat. Los puntos agregados por la extensión y el refinamiento deben ser coherentes con la grilla.
- `D2_sin_atipico` tiene 72 ensayos y se quitó solo la fila 13.

### 1.3 Plausibilidad frente a lo conocido

- Con el modelo correcto y datos simulados, la cobertura del IC de z debería rondar el 95%. Con 360 ajustes, el error estándar de esa proporción es de ~1.2%. Se compara con la Fase 4 (26 de 27).
- En D real, la 1.2.0 debe dar resultados más estables entre semillas que la 1.1.0.
- D3 con el preset final debe ser compatible con la Tabla S1 de Okuyama: z alrededor de 2.5 y límite inferior cerca de 1.4–1.6.
- El beta-binomial de referencia debe reproducir la NLL implícita en la Tabla S1. Ya verificado en la Fase 3.
- Los tiempos por trabajo deben ser coherentes con el piloto de costos. Sirve para la sección de costo de cómputo del manuscrito.

### 1.4 Si algo no pasa

Se identifica la causa (código, script o datos) antes de seguir. Solo se repiten los trabajos afectados, borrando sus `.rds`: el script retoma solo lo que falta.

**Entregable:** sección "Verificación" de `validation/docs/FASE5_resultados.md`.

---

## Paso 2. Interpretar los resultados por brazo

Las métricas y criterios se fijan aquí, **antes** de mirar los números en detalle, para no elegirlos después.

| Brazo | Pregunta | Métricas | Criterio de lectura |
|---|---|---|---|
| A. Simulación amplia | ¿El método recupera los parámetros y da intervalos de z válidos en un rango amplio de parámetros y diseños? | Cobertura del IC de z; IC abiertos; sesgo y error de z, a y h; error de tipo I (z = 1 declarado mayor que 1, por IC y por ΔAIC > 2); potencia para z = 1.75 y 2.75; falla de optimización (déficit > 1) | Cobertura ≥ 90% en general. Se identifican los escenarios con baja cobertura y su causa (k chico, diseño pequeño, pocos encuentros, saturación). |
| B. Desvíos del modelo | ¿Qué pasa si las distribuciones o el proceso no son los supuestos? ¿Ajustar la familia correcta lo corrige? | Sesgo de z y cobertura, ajuste por defecto frente a ajuste correcto; ΔAIC entre ambos (¿el AIC elige el correcto?) | Se documenta para cada desvío si el error es grave y si la familia correcta lo corrige. La heterogeneidad en a no tiene modelo correcto y se declara como limitación. |
| C. Presupuesto | ¿Qué tan fiables son los valores por defecto? ¿Conviene el mejor de 3 o el preset standard? | Éxito (dentro de 0.5 y de 1 NLL de la mejor solución), déficit, cobertura y tiempo para exploratory, mejor de 3 y standard | Da la recomendación de presupuesto para el manuscrito y la viñeta. Decide D9 y D10. |
| D. 1.1.0 frente a 1.2.0 | ¿La 1.2.0 es mejor que la versión publicada? | Cobertura, IC abiertos, error de z, déficit y brecha; en D2 y D3, variabilidad entre semillas | Si mejora, se dice cuánto. Si no mejora en algo, se informa igual. |
| E. Calibración de chequeos | ¿Los chequeos detectan los ajustes que fallan? | Sensibilidad y especificidad de cada chequeo y del semáforo global; curvas sobre los umbrales de las variables continuas (EE de Monte Carlo, brecha optimizador–reevaluación, mejora tardía, dispersión entre corridas, mejora con más esfuerzo o cotas más anchas) | Se eligen umbrales que maximicen la sensibilidad sin una tasa de falsas alarmas inaceptable. Se quitan o reformulan los chequeos que no discriminan. |
| F. D2 y D3 con preset final | Resultados de los datos reales para el material suplementario | z_hat, IC, ΔAIC frente a z = 1 y frente al beta-binomial, chequeos y verificación; D2 con y sin la fila 13 | Se interpreta como coherencia del método, no como prueba de exactitud (punto C8 de Okuyama). |

**Entregables:**
- `FASE5_resultados.md`, con tablas resumen y conclusiones por brazo.
- Figuras de validación en borrador: cobertura por escenario, sesgo de z, éxito por presupuesto, 1.1.0 frente a 1.2.0 y curvas de los chequeos.
- Script reproducible `validation/R/13_analisis_fase5.R`, que genera todas las tablas y figuras desde los `.rds`.

---

## Paso 3. Decisiones derivadas de los resultados (cambios en el código)

1. **Umbrales de `.rel_thresholds()`**: reemplazar los provisorios por los calibrados en E y documentar de dónde sale cada uno. Ajustar los tests si cambian.
2. **Presupuesto por defecto** (D9 y D10, según C):
   - si el mejor de 3 mejora claramente, agregar una opción tipo `n_starts` en `funres_control()` y en la interfaz;
   - si el preset standard es necesario, decir explícitamente que exploratory solo sirve para una primera mirada (ya está así en la app; reforzarlo en la ayuda).
3. **Test de invarianza de la unidad** del ajuste completo: mismos datos con T en horas y en minutos, con la misma semilla. Debe dar las mismas estimaciones convertidas.
4. **Congelar el código** con una etiqueta `v1.2.0-rc1`. Cualquier cambio posterior que toque la estimación obliga a repetir lo afectado.

Si el Paso 3 cambia la estimación (por ejemplo, agregar `n_starts` como opción por defecto), se evalúa si hace falta una corrida corta de confirmación (D8).

---

## Paso 4. ¿Respondemos a Okuyama? Matriz de cumplimiento

Esta matriz se completa con la evidencia del Paso 2 y es la base de la carta y de la respuesta punto por punto.

| C | Comentario | Paquete | Interfaz | Evidencia | Manuscrito | Estado |
|---|---|---|---|---|---|---|
| C1 | Sin "Okuyama model" en el título | Hecho (DESCRIPTION, ayudas, reporte) | Hecho | — | Pendiente (título y texto) | Falta el manuscrito |
| C2 | Distribuciones elegibles y propias | Hecho (3 niveles, parámetros fijos, depredador) | Hecho (familias, editor de muestreadores con prueba) | Brazo B | Pendiente | Falta el manuscrito |
| C3 | Guía de estabilidad, convergencia y fallos | Hecho (`reliability()` y `verify_fit()` con acción) | Hecho (pestaña Reliability) | Brazo E | Pendiente (tabla síntoma, causa y acción) | Falta calibrar, la viñeta y el manuscrito |
| C4 | DEoptim frente a alternativas | Hecho (DEoptim, CMA-ES, CRS2, números aleatorios comunes y pulido) | Hecho (selector de método) | Fase 4 | Pendiente | **Falta ABC (D7)** |
| C5 | Cotas y unidad de T | Hecho (adimensional, ampliación automática, editables) | Hecho | Fase 3, chequeo de cotas en E, test del Paso 3 | Pendiente | Falta el test de invarianza |
| C6 | Corridas independientes y sensibilidad al esfuerzo; qué hacer si no coinciden | Hecho (`verify_fit()`) | Hecho (botón de verificación) | Brazos C y E | Pendiente | Falta calibrar y el manuscrito |
| C7 | z conjunta y perfil solo para el IC | Hecho (`z_method = "joint"`) | Hecho | Fase 4: cobertura 67% contra 96% del perfil | Pendiente (justificar por qué el perfil queda por defecto) | Respondido con evidencia |
| C8 | Datos reales a material suplementario; validación más amplia | — | — | Brazos A, B y F | Pendiente (D2 y D3 a material suplementario, validación al texto principal) | **Revisar alcance (D8)** |
| C9 | Advertir que los valores por defecto pueden fallar | Hecho (presets, notas, semáforo) | Hecho (aviso permanente en exploratory) | Brazo C | Pendiente | Falta el manuscrito |
| C10 | Coautoría solo si todo está atendido | — | — | Esta matriz | Carta | Paso 6 |

**Brechas conocidas frente al plan original (Fase 5.2 y 5.3) y frente a C8:**
- **Diseños:** se evaluaron 2 diseños (8 densidades × 8 réplicas y 5 × 6). No se variaron T, el número de densidades entre 4 y 13 ni las réplicas entre 5 y 20.
- **Desvíos no evaluados:** agotamiento de huevos o mortalidad, aprendizaje o tasa variable, selección no aleatoria de hospedadores y tratamiento del final del ensayo.
- **Método de muestreo:** se usó un diseño factorial con variaciones de a un factor, no un hipercubo latino. Es defendible porque se lee más fácil, pero conviene declararlo.

Estas brechas no impiden responder, pero hay que decidir si se cubren con una corrida adicional (D8) o si se declaran como alcance y limitaciones.

---

## Paso 5. Documentación y calidad del paquete (Fase 7)

1. `NEWS.md` de la 1.2.0, con la sección "Changes that alter results": cotas, IC sobre la NLL reevaluada, valores por defecto del optimizador y extensión de la grilla hacia abajo.
2. `README.md` reescrito.
3. Viñetas:
   - (a) inicio rápido;
   - (b) "How reliable is my fit?", con la tabla síntoma, causa y acción y los umbrales calibrados;
   - (c) distribuciones propias en R y en C++;
   - (d) unidades y cotas;
   - (e) resumen de la validación, con enlace a `validation/`.
4. `R CMD check --as-cran` limpio en Linux (acá) y en Windows (tu PC). Tests largos salteados en CRAN.
5. Opcional: `shinytest2` para un flujo completo de la interfaz.

Se puede avanzar en paralelo con el Paso 2, porque no depende de los números finales salvo en la viñeta (b).

---

## Paso 6. Manuscrito (Fase 8) y cierre (Fase 9)

1. Reescribir el manuscrito PLOS con control de cambios, inglés americano y las convenciones del proyecto:
   - nuevo título;
   - Design and Implementation con el motor generalizado;
   - Results centrado en la validación (A a E);
   - guía práctica;
   - D2 y D3 en el material suplementario (F);
   - Limitations actualizada.
2. Material suplementario: S1 Text (método y cotas), S2 File (repositorio y scripts), tablas y figuras de validación y reanálisis de D2 y D3.
3. Control del límite de palabras de PLOS Computational Biology.
4. Respuesta punto por punto a Okuyama, a partir de la matriz del Paso 4, y carta. Se envía **antes** de enviar al journal, para que decida sobre la coautoría.
5. Con tu aprobación:
   - subir la rama a GitHub;
   - publicar la release y la versión nueva de Zenodo;
   - enviar a CRAN;
   - citar el DOI final en el manuscrito.

---

## Decisiones pendientes

| ID | Decisión | Cuándo | Recomendación |
|---|---|---|---|
| D7 | ABC (pedido explícito de Okuyama en C4): ¿evaluarlo o justificar que no? | Después del Paso 2 | Justificarlo en el texto (ABC estima una posterior, no un máximo de verosimilitud, así que no es comparable con igual presupuesto), con una demostración pequeña en un subconjunto si el cómputo lo permite. |
| D8 | ¿Corrida adicional (Fase 5b) para diseños y desvíos no cubiertos? | Después del Paso 2 | Decidir según A y B. Si la cobertura es buena y estable entre los 2 diseños, alcanza con una corrida chica que varíe T y el número de densidades (más barata que la Fase 5). |
| D9 | ¿Agregar "mejor de N" como opción del paquete? | Después del brazo C | Según la mejora observada. |
| D10 | ¿Cambiar el preset por defecto o reforzar el aviso? | Después del brazo C | Según el éxito de exploratory frente a standard. |

## Orden y responsables

1. **Paso 1 y Paso 2** (yo, en este entorno, con las tablas ya copiadas): verificación, interpretación, script 13 y `FASE5_resultados.md`. Al terminar, te paso el resumen y las decisiones D7 a D10.
2. **Paso 3** (yo): cambios en el código. Vos corrés `devtools::install()` y `devtools::test()` en Windows.
3. **D8**: si hace falta la Fase 5b, la corrés en tu PC.
4. **Paso 5** (yo, en paralelo).
5. **Paso 6** (los dos): el manuscrito, por tandas y con tu revisión; la liberación, solo con tu aprobación.
