# Roadmap maestro

## Dos productos, una sola dependencia

El repositorio produce dos entregables distintos:

1. un paquete R reutilizable para LRM, HLRM y CHLRM;
2. una suite científica que ejecuta la revisión computacional del paper.

La dependencia es unidireccional:

```text
paper/evaluation  ──usa──>  paquete instalado
       │                         │
       └── configura y reporta   └── calcula y devuelve objetos
```

El paquete nunca importa código de `paper/`, `thesis/` ni conoce identificadores
como E1--E5 o S1--S7c. La suite no copia kernels ni implementa variantes de
MCMC, CAVI, SVI, likelihood, predicción o relabeling.

## Autoridad documental

| Documento | Autoridad | No debe contener |
|---|---|---|
| [`API_CONTRACT.md`](API_CONTRACT.md) | interfaz pública y significado de resultados | escenarios, réplicas o rutas del paper |
| [`PACKAGE_ROADMAP.md`](PACKAGE_ROADMAP.md) | implementación, pruebas y publicación del paquete | tablas, figuras o configuraciones E1--E5 |
| [`paper/EVALUATION_PLAN.md`](paper/EVALUATION_PLAN.md) | requisitos del PDF, experimentos y criterios de aceptación | kernels o APIs alternativas |
| [`DEVELOPMENT_DECISIONS.md`](DEVELOPMENT_DECISIONS.md) | registro histórico y justificación | tareas vigentes |
| `README.md` | orientación y uso básico | sustituir alguno de los planes anteriores |

Si dos documentos discrepan, manda el documento especializado de esta tabla.
El registro histórico no puede reabrir una decisión vigente.

## Regla para asignar responsabilidades

Una capacidad pertenece al paquete cuando puede definirse sin saber qué
experimento, manuscrito o dataset motivó su uso y resulta útil para un usuario
externo. Pertenece a la evaluación cuando fija una decisión científica del
protocolo o combina muchos ajustes para responder una pregunta empírica.

Prueba práctica:

- si la función necesita conocer `E1`, `S4`, `Farms`, 200 réplicas, 5000 draws,
  una tabla o una ruta de salida, pertenece a `paper/evaluation/`;
- si recibe datos y parámetros arbitrarios y devuelve un objeto estadístico con
  significado estable, puede pertenecer al paquete;
- si sólo implementa una ecuación interna, permanece privada aunque sea
  reutilizada por varias clases.

## Frontera de responsabilidades

| Capacidad | Paquete R | Evaluación del paper |
|---|---|---|
| Validación de datos y parámetros | genérica | valida además la configuración del protocolo |
| Simulación | DGP condicional, jerárquico y prior-predictivo parametrizables | geometría de coeficientes, recetas S1--S7c, datasets anidados y semillas |
| Priors | constructores explícitos por modelo, presets generales y validación | valores exactos principal/sensibilidad del paper |
| Inferencia | MCMC, CAVI, SVI y sus contratos | métodos, cadenas, inicios y controles exigidos |
| Relabeling y draws | algoritmos generales | referencia y número de draws de cada experimento |
| Predicción y likelihood | primitivas con objetivo/unidad explícitos | folds, targets y tamaños concretos |
| Diagnósticos por ajuste | convergencia y costo de un ajuste | reglas de aceptación y resumen entre réplicas |
| Métricas observadas generales | scoring y PPC configurables | sesgo/RMSE/cobertura contra verdad, CPE y regret |
| Persistencia | objetos serializables y checkpoint explícito | rutas, manifiesto, reintentos y procedencia |
| Recursos | límites y estimadores por llamada | paralelismo y presupuesto global del estudio |
| Comunicación | documentación de métodos y objetos | tablas, figuras, respuesta a revisores y LaTeX |

La API no tiene que producir directamente las tablas del artículo. Sí debe
exponer suficiente información correcta y auditable para que la suite pueda
producirlas sin acceder al estado privado de las clases.

## Brecha de frontera que debe cerrarse antes de congelar la API

La frontera ya fue cerrada: los antiguos helpers del protocolo se retiraron de
`R/` porque mezclaban mecanismos reutilizables con recetas editoriales. La
separación vigente es:

- el paquete conservará `conditional_dgp()`, `hierarchical_dgp()` y
  `prior_predictive_dgp()`. Todos generarán diseño/respuesta según una
  especificación explícita y devolverán la verdad realizada;
- los constructores de prior serán explícitos y uniformes por modelo. Nunca
  tomarán automáticamente la verdad del simulador como información previa;
- `paper/evaluation/` construirá la geometría trigonométrica con separación
  `d`, fijará los hiperparámetros y declarará S1--S7c;
- la predictiva posterior seguirá siendo una operación sobre un ajuste y no se
  confundirá con los generadores Monte Carlo;
- los consumidores y tests usan el contrato nuevo; no quedaron wrappers
  históricos permanentes.

El escalado no es responsabilidad de los kernels, pero sí una capacidad
reutilizable del paquete mediante `data_scaling()` y un mapa explícito y
reversible. `fit_model()` puede orquestarlo cuando el usuario entrega esa
especificación. La suite decidirá cuándo usarlo: escalará Farms sin fuga entre
folds y no reescalará E1--E4, porque su escala forma parte del mecanismo
generador.

DIC fue retirado del paquete y de la evaluación final. WAIC permanece con
unidad/objetivo explícitos, acompañado por PSIS-LOO o CV/LOGO y MLPD cuando
corresponda; no se presentará como comparable entre escalas, unidades o
condicionamientos distintos.

## Secuencia vinculante

1. Ejecutar I0--I6 del orden atómico definido en
   [`PACKAGE_ROADMAP.md`](PACKAGE_ROADMAP.md): retirar DIC; construir la
   fachada tipada; unificar modelo/prior/likelihood/DGP; añadir escalado;
   completar evaluación predictiva; extender SVI a LRM; y cerrar la fachada
   idiomática instalable.
2. Desarrollar en paralelo sólo las recetas, métricas y orquestación propias de
   `paper/evaluation/`, sin copiar kernels ni depender de símbolos privados.
3. Instalar el paquete, aprobar `R CMD check` y verificar la Puerta P completa.
4. Ejecutar `smoke` y comprobar schemas, reanudación y procedencia.
5. Ejecutar `pilot` para fijar tiempo, RAM, almacenamiento y configuración SVI.
6. Congelar versión del paquete, configuración y semillas.
7. Ejecutar E3, E1, E2, E4 y E5 en perfil `full`.
8. Generar artefactos y actualizar el manuscrito desde resultados trazables.

No se ejecutará el perfil confirmatorio mientras exista una brecha que cambie
el estimando, el algoritmo o el contrato consumido por el paper.

## Estado sincronizado

### Logrado y probado

- P1 está cerrado para los motores existentes: relabeling CHLRM, MCMC
  multichain, diagnósticos, CAVI/SVI multistart y procedencia de ejecuciones.
- La semántica de resultados de P2 está cerrada: estado, RNG, validación,
  álgebra SPD, predicción, predictiva posterior y nombres aplicables por modelo
  sin campos ficticios `NULL`. La construcción pública tipada de
  modelo/prior/inicialización/control, `fit_model()` y `convergence-1` ya
  funcionan para los nueve motores.
- `posterior_draws()` expone todos los parámetros aplicables con una orientación
  uniforme; conserva cadena/iteración en MCMC y distingue muestras de (q) en
  CAVI/SVI. `posterior::as_draws()` produce `draws_df` y conserva además la
  procedencia MCMC original cuando el submuestreo deja iteraciones no contiguas.
- Existe un contrato único de log-likelihood con unidad, objetivo,
  condicionamiento y variables integradas explícitos.
- WAIC marginal de HLRM/CHLRM, predicción, PPC streaming configurable T1--T5,
  MCSE dependiente y reutilización de `y_rep` ya tienen infraestructura y
  pruebas.
- DIC, `include_dic`, sus puntos posteriores y campos de resultado fueron
  retirados; todos los consumidores activos usan `predictive-criteria-2`.
- Los CSV/RDS exploratorios generados antes de I0 no se reescriben manualmente:
  se consideran artefactos obsoletos y se regenerarán con la API congelada para
  no presentar resultados mezclados de dos contratos.
- `simulate_data()` implementa DGP condicional, jerárquico y prior-predictivo,
  diseño/tamaño tipados, train/test, procedencia y `tidy_truth()` en
  `simulation-result-2`; los simuladores productivos anteriores fueron retirados.
- La última validación del paquete ejecutó 150 bloques, sin fallos, errores,
  warnings ni skips. La validación rápida de infraestructura del paper ejecutó
  otros 35 bloques sin incidencias; ninguno es un experimento confirmatorio.

### Decidido y especificado, pero todavía no implementado

- ejecutar el factorial piloto E3: el `screen` de 10 tareas ya terminó y su
  gate reproducible aprobó completitud, procedencia y margen de dos veces el
  RSS máximo (496 MiB); proyecta 109 tareas secuenciales en 17.3 minutos y
  143.9 MiB. El nuevo runner exige congelación y aceptación explícita de las
  huellas del factorial y del `screen`;
- decidir la interfaz opcional fórmula/datos; los métodos base, accesores
  jerárquicos y CHLRM, salidas tidy, extractores orientados para
  `loo`/`bayesplot` y la conversión al paquete `posterior` ya están
  implementados;
- añadir vignettes y profiling empírico de RSS/tiempo antes del piloto.

### Pendiente para cerrar el paquete

- P3 está cerrado; quedan como tareas experimentales —no de la API— la
  agregación Monte Carlo entre réplicas y las configuraciones concretas;
- convertir la documentación inicial a roxygen, añadir vignettes y completar
  interoperabilidad; el paquete ya instala y `R CMD check --no-manual` está OK;
- medir RSS/tiempo en escenarios piloto y congelar límites seguros; el
  estimador previo, el presupuesto uniforme por ajuste y el registro smoke ya
  están implementados;
- decidir si la comodidad de `formula + data` compensa ampliar la superficie;
  la interfaz matricial y `grouped_regression_data()` ya son completas.

### Pendiente para la evaluación del paper

- E1--E5 ya tienen smoke que instala y consume el paquete, reanuda por
  artefacto y produce CSV inspeccionables;
- las configuraciones E1--E5 y las métricas puras principales ya están
  declaradas; el alineamiento paramétrico con verdad ya alimenta E1. Faltan la
  agregación Monte Carlo completa por parámetro, figuras y tablas finales;
- conservar `smoke` como prueba de integración, ejecutar luego `pilot`,
  congelar versión/configuración/semillas y sólo entonces ejecutar `full`.

### Siguiente incremento

I0--I5 están cerrados. I6 ya cubre la fachada tipada, la instalación, los
métodos S3, los extractores interoperables, `CAVI` como identificador único y
la migración de los consumidores activos. La frontera se verifica mediante
`dev/check-consumer-boundary.R`. Para cerrar I6 faltan documentación de usuario
en forma de vignettes y una decisión explícita sobre `formula + data`.

I2/P2.2 está cerrado. Los tres DGP, `simulation-result-2`, la migración de
consumidores, la verdad compacta opt-in, el presupuesto de memoria, la primitiva
observacional compartida y el smoke SBC están implementados. En CHLRM el smoke
usa cantidades invariantes a etiquetas y excluye parámetros crudos no
identificables.

I3/P2.3 está cerrado. `fit_model(..., scaling = data_scaling(...))` aprende el
mapa sólo con entrenamiento y aplica la misma transformación a datos, prior e
inicialización en MCMC/CAVI/SVI. Parámetros, `posterior_draws()`, resúmenes,
ajuste dentro de muestra, predicciones, PPC y log-likelihood regresan juntos a
la escala de entrada; WAIC se recalcula después del jacobiano por unidad.

I4/P3 está cerrado. Existen scoring retenido ELPD/MLPD, métricas predictivas de
test e intervalos, planes de CV por observación/dentro de grupo/por grupo/LOGO
y un runner secuencial con límite de ajustes; el escalado se aprende
independientemente en cada fold. I5 añadió SVI-LRM por observaciones, con
actualizaciones naturales y equivalencia CAVI a batch completo. La suite
completa valida 150 bloques sin
fallos, errores, advertencias ni skips.
`psis_loo()` usa la implementación estándar de `loo`, sólo acepta MCMC y exige
convergencia multichain por defecto; CAVI/SVI usan CV refitted. Siguen métricas
agregadas contra verdad y comparaciones entre métodos en la suite del paper,
no en la API.
En paralelo, el siguiente corte científico es congelar y ejecutar el factorial
piloto E3, ya que el `screen` y su revisión de recursos están completos; sólo
después se congelará E3 `full` y se ejecutarán los demás experimentos
confirmatorios. La lista completa I0--I7 y sus
criterios de cierre vive sólo en `PACKAGE_ROADMAP.md` para no mantener dos
planes operativos divergentes.
