# Origen y verificación de D2 y D3 (Kishani Farahani & Goldansaz 2013)

Fecha: 2026-10-19. Fuentes revisadas: el artículo (Eur. J. Entomol. 2013, PDF en `referencias_biblio`), la base FoRAGE V5 (`FoRAGE_db_V5_Dec_20_2024_original_curves.csv` y `..._sources_and_meta.csv`) y `funresMech_validacion/data_clean/forage_D2_D3_raw.csv`.

## Correspondencia

| Conjunto | FoRAGE | Figura del artículo | Estadio del hospedador | Puntos |
|---|---|---|---|---|
| (no usado) | 56 | Fig. 1, L1 | 1.er estadio | 66 |
| **D2** | **57** | **Fig. 1, L2** | **2.º estadio** | **73** |
| **D3** | **58** | **Fig. 1, L3** | **3.er estadio** | **64** |

FoRAGE registra "Data origin = Point": los datos se digitalizaron punto por punto de la Figura 1.

## El ensayo atípico de D2 (2 parasitados de 10)

**Verificado:** el punto está en la Figura 1 (panel L2) del artículo, a densidad 10 y con ~2 hospedadores parasitados, muy por debajo del resto. No es un error de digitalización.

## Hallazgo importante: los datos digitalizados no son todos los ensayos

- **Diseño del artículo:** 13 densidades (2–50), con 25 réplicas para las densidades 2–20 y 20 réplicas para 25–50. Son **295 ensayos por estadio**.
- **D2 y D3 tienen 73 y 64 puntos.** En cada densidad aparecen casi solo los **valores distintos**: a densidad 2 hay 3 puntos (0, 1, 2) en vez de 25 ensayos, y a densidad 10 hay 3 puntos (2, 8, 10) en lugar de 25. En el gráfico de dispersión los ensayos con el mismo número de parasitados se superponen, y al digitalizarlos quedó un solo punto.
- **Consecuencias:**
  - Se pierden las frecuencias. Un valor raro (como el 2 de 10) pesa en la verosimilitud igual que el valor más frecuente.
  - La dispersión se infla, lo que se ve como k chico en D2.
  - Hay menos información de la que tuvo el experimento.
- **Okuyama (2026) usó los mismos datos de FoRAGE**, así que la limitación es compartida. Refuerza su pedido de que D2 y D3 se presenten como ilustración (material suplementario) y no como prueba de exactitud.

## Cómo se midió el parasitismo

- Una hembra naïve, apareada y de 1 día, por tubo de vidrio (7 cm de diámetro, 20 cm de alto), durante 24 h, a 25 °C. Después se retiró la avispa.
- **Parasitismo = avispas adultas emergidas por tubo.** No hubo disección. Se cuentan los hospedadores que produjeron un parasitoide, no las oviposiciones.
- **No se informa superparasitismo, discriminación de hospedadores ni host-feeding.** Con conteo por emergencia el superparasitismo no se puede detectar: de un hospedador superparasitado de un parasitoide solitario emerge un solo adulto.
- **La mortalidad de hospedadores o de parasitoides antes de la emergencia reduce el conteo.** Es una explicación posible y no comprobable del 2 de 10.
- *A. myeloenta* es sinovigénica. La fecundidad diaria medida en otro ensayo (5.0–5.5 adultos por día, con 30 hospedadores diarios) es mucho menor que lo observado en la respuesta funcional (25–33 parasitados en 24 h a densidad 50). El artículo no lo discute. No permite evaluar la limitación por huevos.

## El análisis original

Los autores concluyeron una respuesta de **tipo II** para los tres estadios (regresión logística de Juliano; ecuación de Rogers con Th ≈ 0.56–0.58 h). Con el modelo mecanístico, el proceso parasitoide da z > 1 (tipo III) y el proceso depredador da z ≈ 1 (tipo II), con mejor AIC para el parasitoide (ver `PROCESO_parasitoide_vs_depredador.md`). El tiempo de manipulación es coherente entre todos los análisis (0.56–0.62 h en el artículo; ~0.6 h en funresMech).

## Qué decir en el manuscrito sobre D2 y D3

1. Son datos digitalizados (FoRAGE, Figura 1 de Kishani Farahani & Goldansaz 2013). Los ensayos con valores repetidos se superponen en la figura, así que el conjunto tiene menos puntos que el experimento y no conserva las frecuencias.
2. El ensayo de 2 sobre 10 existe en la figura original.
3. El parasitismo se midió por emergencia de adultos, así que el superparasitismo no es observable.
4. Por todo esto, D2 y D3 ilustran el uso del paquete y su coherencia con la implementación de Okuyama. No sirven para estimar la biología de *A. myeloenta*.
