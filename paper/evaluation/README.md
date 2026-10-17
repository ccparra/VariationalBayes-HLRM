# Suite de evaluación del paper

Aquí vivirán los perfiles `smoke`, `pilot` y `full` y los experimentos E1--E5.
Este directorio sólo puede contener:

- configuraciones y expansión de tareas;
- escenarios específicos S1--S7c y preprocessing específico de Farms;
- métricas que requieren verdad simulada o un método de referencia;
- scheduler, manifiestos, checkpoints y agregación;
- construcción automática de tablas y figuras.

No se permiten copias de kernels, likelihoods, predicción, relabeling o
diagnósticos del paquete. Mientras no se cierre la puerta P de
`paper/EVALUATION_PLAN.md`, sólo se construirá infraestructura `smoke`; no se
ejecutarán resultados confirmatorios. SVI para LRM no forma parte de esa puerta.

Los presets que reproducen el PDF —prior principal, geometría de coeficientes,
S1--S7c, tamaños, semillas y grids— pertenecen aquí aunque utilicen
constructores genéricos del paquete. Las métricas contra verdad simulada o
contra MCMC también pertenecen aquí. En cambio, scoring predictivo sobre datos,
PPC configurable, log-likelihood pointwise y diagnósticos de un ajuste son
capacidades generales del paquete.

`R/scenarios.R` contiene la geometría, la prior y el catálogo declarativo E1:
`make_review_component_coefficients()`, `review_chlrm_prior()`,
`e1_scenario_catalog()` y `scenario_e1()`. Estas funciones sólo configuran
`conditional_dgp()`, `simulation_design()`, `simulation_sample()` y
`chlrm_prior()`; no duplican simulación ni inferencia del paquete. También
declara protocolos versionados E2--E5, separando el objetivo confirmatorio de
los tamaños reducidos usados por `smoke` y `pilot`.

El runner instala el paquete en una biblioteca temporal y ejecuta las pruebas
contra esa instalación. Al correr un archivo de prueba de manera aislada,
`pkgload::load_all(export_all = FALSE)` sigue disponible como modo de desarrollo
sin exponer símbolos privados.

`R/sbc.R` contiene rangos SBC y la reducción de CHLRM a cantidades invariantes
a permutaciones. Los coeficientes/varianzas/pesos crudos de componente y las
asignaciones no se usan como estimandos SBC sin relabeling explícito. El comando
`Rscript paper/evaluation/run-validation.R` ejecuta únicamente pruebas rápidas
de infraestructura y guarda el resumen en `paper/evaluation/output/`; no corre
réplicas confirmatorias. `R/metrics.R` implementa las métricas que sólo tienen
sentido al conocer la verdad o una referencia MCMC: ARI, CPE, entropía,
recuperación, cobertura, fidelidad posterior y regret. `R/aggregation.R` lee
artefactos terminados y no vuelve a ajustar modelos.

`R/tasks.R` expande E1 en 6 tareas `smoke`, 135 `pilot` o 5400 `full`. Sólo el
perfil smoke puede ejecutarse por ahora. `Rscript paper/evaluation/run-smoke.R`
instala el paquete, guarda un RDS atómico por tarea y reutiliza resultados
válidos al reanudar. `--fresh` regenera las seis tareas. El manifiesto separa
tiempo total de tarea, tiempo de ajuste, estimación de memoria y tamaño real del
RDS; `resource-summary.csv` los resume por método. `metrics-raw.csv` contiene
ARI, CPE, entropía y desempeño retenido por tarea; `metrics-summary.csv`
presenta medias y MCSE por escenario/método. Sus artefactos inspeccionables
viven en `output/smoke/`. Por usar sólo cuatro u ocho draws/iteraciones, esos
valores verifican integración y schemas, no son resultados científicos.

`parameter-recovery-raw.csv` conserva cada coordenada de coeficientes,
varianzas y pesos después de alinear componentes con la verdad;
`parameter-recovery-summary.csv` agrega por escenario/método/coordenada sesgo,
RMSE, cobertura, ancho, dispersión posterior y MCSE sin descartar el nivel
pointwise.

`R/provenance.R` genera `paper-run-provenance-1`. Cada runner escribe
`provenance.csv`/`provenance.rds` y añade a sus manifiestos huellas separadas de
configuración, código de cálculo, código completo de la suite, run y salida,
además de versión del paquete, commit y estado del árbol. Un RDS sólo se
reutiliza cuando su `run_fingerprint` coincide: cambios estadísticos invalidan
el caché, mientras cambios exclusivos de reporting regeneran salidas sin
repetir ajustes.

Las figuras de integración se construyen sin cargar el paquete ni ajustar:

```sh
Rscript paper/evaluation/run-smoke-figures.R
```

El comando genera seis PNG y `figures/figure-manifest.csv` a partir de los CSV
existentes. Todos se rotulan como `smoke` y no deben citarse como evidencia del
paper.

`Rscript paper/evaluation/run-e2-smoke.R` verifica de extremo a extremo la
selección de K bajo separación moderada y débil: CAVI/SVI sobre una grilla
reducida y MCMC sobredimensionado. Escribe `manifest.csv`, `fit-results.csv` y
`selection-results.csv` en `output/smoke/E2/`; `--fresh` fuerza regeneración.
Sus K, tamaños y presupuestos reducidos son deliberadamente no
confirmatorios.

Los runners restantes cubren las otras preguntas del review:

```sh
Rscript paper/evaluation/run-e3-smoke.R
Rscript paper/evaluation/run-e4-smoke.R
Rscript paper/evaluation/run-e5-smoke.R
```

E3 guarda la referencia CAVI y los deltas SVI; E4 conserva la comparación
pareada sobre prefixes anidados y la pendiente log--log; E5 ejecuta WAIC
marginal por finca, CV de cinco folds dentro de finca con MLPD y MSE,
estabilidad entre inicios y una corrida MCMC sobredimensionada. Todos aceptan
`--fresh`; sin esa opción reutilizan artefactos compatibles.

Antes de ejecutar E3 en tamaños científicos, el preflight materializa los
controles exactos y estima recursos sin ajustar ningún modelo:

```sh
Rscript paper/evaluation/run-e3-preflight.R --profile=pilot --stage=screen
Rscript paper/evaluation/run-e3-preflight.R --profile=pilot --stage=factorial
```

`screen` contiene una referencia CAVI, las ocho esquinas del factorial SVI y
su centro (10 tareas); `factorial` contiene una referencia CAVI y las 108
combinaciones SVI del piloto. Los artefactos de `output/preflight/E3/` incluyen
plan, presupuesto, configuración candidata, procedencia y una puerta de
ejecución explícitamente bloqueada. La estimación dimensional no se presenta
como RSS: primero se deben medir tiempo, RSS y disco con `screen`, revisar los
fallos y congelar la huella aceptada. El preflight nunca ejecuta CAVI/SVI.
La referencia CAVI usa cinco inicios y el mejor ELBO. Cada semilla SVI de E3 es
en cambio una corrida independiente con un inicio: así el diseño conserva
exactamente 36 × 20 = 720 corridas y puede estimar la proporción que no converge;
añadir cinco `starts` internos produciría erróneamente 3600 ajustes y ocultaría
fallos al seleccionar sólo el mejor.

Después de revisar el CSV, se debe crear el commit, volver a generar el
preflight desde ese árbol limpio y copiar deliberadamente la huella mostrada
en `provenance.csv` para congelarlo:

```sh
Rscript paper/evaluation/freeze-e3-preflight.R \
  --stage=screen \
  --accept-configuration-fingerprint=md5:HUELLA_REVISADA
Rscript paper/evaluation/run-e3-pilot.R --execute --stage=screen \
  --accept-configuration-fingerprint=md5:HUELLA_REVISADA
Rscript paper/evaluation/review-e3-pilot-screen.R
Rscript paper/evaluation/run-e3-pilot-figures.R
```

El runner admite sólo las etapas piloto `screen` y `factorial`, una tarea por
proceso y de forma secuencial; `full` permanece bloqueado. Registra tiempo de
pared, RSS máximo del sistema operativo, logs y RDS reanudables. Si una tarea
falla, se detiene y conserva las anteriores. La revisión del `screen` escribe
`screen-review.csv`, `screen-gate.csv`, `screen-ranking.csv` y el RDS que
vincula la evidencia aceptada con la congelación factorial. El gate exige dos
veces el RSS máximo observado dentro del presupuesto: no interpreta ese pase
de recursos como evidencia de calidad estadística.

El reporte gráfico genera tres PNG desde CSV/RDS terminados, sin volver a
ajustar modelos: comparación de ELBO/MLPD/ARI/CPE, perfil de tiempo/RSS/disco y
trayectorias del ELBO suavizado. Su manifiesto marca explícitamente
`scientific_result = FALSE` porque el `screen` es preliminar.

Para el factorial, primero se regenera su preflight desde el commit limpio. La
congelación y la ejecución requieren aceptar tanto su propia huella como la
huella de salida del `screen` aprobada:

```sh
Rscript paper/evaluation/run-e3-preflight.R --profile=pilot --stage=factorial
Rscript paper/evaluation/freeze-e3-preflight.R \
  --stage=factorial \
  --accept-configuration-fingerprint=md5:HUELLA_FACTORIAL \
  --accept-screen-output-fingerprint=md5:HUELLA_SCREEN
Rscript paper/evaluation/run-e3-pilot.R --execute --stage=factorial \
  --accept-configuration-fingerprint=md5:HUELLA_FACTORIAL \
  --accept-screen-output-fingerprint=md5:HUELLA_SCREEN
Rscript paper/evaluation/review-e3-pilot-factorial.R
Rscript paper/evaluation/run-e3-factorial-figures.R
```

El factorial piloto contiene 109 tareas. La configuración actual fuerza un
solo thread BLAS/OpenMP/MKL y no inicia procesos concurrentes. Su revisión
escribe agregados de las tres semillas para las 36 celdas, MCSE, dispersión,
frecuencia de ARI igual a cero, ranking exploratorio y auditoría de recursos.
Las figuras muestran mapas de MLPD/estabilidad, ARI/colapso y distribuciones de
tiempo/RSS. Ni el ranking ni las figuras desbloquean automáticamente `full`.

La suite puede escribir archivos y decidir su retención. El paquete sólo
devuelve objetos serializables o escribe un checkpoint cuando se invoca una
operación explícita para ello.
