# Código y experimentos de la tesis final

Todo se ejecuta desde la raíz del repositorio. `runs/` contiene los análisis
centrales vigentes y consume exclusivamente la API instalada de
`VariationalBayesHLRM`. Cada ejecución escribe resultados regenerables en su
propio `output/`. Las validaciones metodológicas nuevas viven en
`paper/evaluation/`; no se mantienen copias históricas rotas dentro de la tesis
final.

## Pipeline central

```sh
Rscript thesis/final/code/run_all.R
```

Antes de una corrida estándar, la ruta completa puede verificarse en segundos:

```sh
VB_HLRM_THESIS_PROFILE=smoke Rscript thesis/final/code/run_all.R
VB_HLRM_THESIS_PROFILE=smoke Rscript thesis/final/code/runs/04-svi-step-tuning/run.R
VB_HLRM_THESIS_PROFILE=smoke Rscript thesis/final/code/runs/05-select-k-chlrm/run.R
```

Los artefactos incluyen la columna o campo `profile`; un resultado `smoke` es
una prueba de integración y no evidencia científica para la tesis o el paper.

El pipeline corto ejecuta:

- `runs/01-farms-hlrm/`: ajuste jerárquico de fincas;
- `runs/02-farms-chlrm/`: mezcla y coclustering de fincas;
- `runs/03-simulation/`: recuperación con verdad conocida.

Los análisis costosos se ejecutan por separado:

```sh
Rscript thesis/final/code/runs/04-svi-step-tuning/run.R
Rscript thesis/final/code/runs/05-select-k-chlrm/run.R
```

Estos evalúan el esquema Robbins-Monro de SVI y el número de componentes del
CHLRM mediante multistart y validación dejando fincas completas fuera.

Los runs 01--03 usan priors fijas explícitas y separan prior de inicialización;
no reconstruyen priors empirical-Bayes con OLS. Los runs 04--05 son análisis
costosos de tasa SVI y selección de K y deben ejecutarse por separado. Todos
usan CAVI como nombre del método determinista, controles tipados, límites de
memoria y resultados versionados.

Para ejecuciones grandes se recomienda consultar primero
`estimate_fit_resources()` y comenzar con el perfil `smoke` del paper. La
matriz de coclustering se solicita sólo cuando el reporte la necesita.

Las decisiones estadísticas y su conexión con la implementación están en
[`MODEL_COHERENCE.md`](MODEL_COHERENCE.md).
