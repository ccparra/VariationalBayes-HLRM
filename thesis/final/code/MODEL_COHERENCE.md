# Coherencia estadística y decisiones de implementación

## HLRM

El modelo implementado es

$$
y_j\mid\beta_j,\sigma_j^2\sim N(X_j\beta_j,\sigma_j^2I),\qquad
\beta_j\mid\beta,\Sigma\sim N(\beta,\Sigma),
$$

con $\sigma_j^2\mid\xi^2\sim IG(\nu_0/2,\nu_0\xi^2/2)$, una prior normal para $\beta$, inverse-Wishart para $\Sigma$ y gamma para $\xi^2$.

Las actualizaciones MCMC y VI siguen esa parametrización shape-rate. En particular,

$$
\xi^2\mid-\sim Gamma\left(a_0+\frac{m\nu_0}{2},
b_0+\frac{\nu_0}{2}\sum_j\sigma_j^{-2}\right).
$$

Para VI se usa

$$
E\|y_j-X_j\beta_j\|^2=
\|y_j-X_jE\beta_j\|^2+tr(X_j'X_j Var(\beta_j)),
$$

sin volver a sumar el término cuadrático de la media.

## CHLRM de K fijo

La implementación adopta explícitamente una mezcla finita de $K$ fijo:

$$
z_j\sim Cat(\omega),\quad
\beta_k\mid\beta,\Sigma\sim N(\beta,\Sigma),\quad
y_j\mid z_j=k\sim N(X_j\beta_k,\sigma_k^2I).
$$

Por tanto, `Betaj` y `sigma2j` son vistas derivadas (`Betak[z,]` y `sigma2k[z]`), no parámetros con una distribución condicional adicional. Los componentes vacíos se actualizan desde su prior y siguen contando en las condicionales globales de `Beta`, `Sigma` y `xi2`. Mezclar aquí reglas de un modelo de componentes ocupados cambiaría silenciosamente la posterior.

## Convención de notación y nombres

La API pública usa nombres descriptivos en inglés y `snake_case`. Los kernels
mantienen símbolos cercanos a las ecuaciones sólo dentro de su estado
computacional; no se exponen como contrato de usuario. La correspondencia
canónica es:

| Símbolo | Significado | Ranura pública |
|---|---|---|
| $\beta$ | coeficientes poblacionales | `population_coefficients` |
| $\Sigma$ | covarianza poblacional | `population_covariance` |
| $\xi^2$ | escala jerárquica de varianzas | `variance_scale` |
| $\beta_j$ | coeficientes del grupo observado | `group_coefficients` |
| $\sigma_j^2$ | varianza residual del grupo | `group_variances` |
| $\beta_k$ | coeficientes del componente | `component_coefficients` |
| $\sigma_k^2$ | varianza residual del componente | `component_variances` |
| $\omega$ | pesos de mezcla | `mixture_weights` |
| $z_j$ | asignación de grupo | `group_allocations` |
| $\rho_{jk}$ | probabilidad/responsabilidad de asignación | `allocation_probabilities` |

En VI/SVI, `mean`, `covariance`, `shape`, `scale`, `rate`, `df` y
`concentration` identifican explícitamente la parametrización del factor. En
MCMC la misma ranura contiene `distribution = "posterior draws"` y `draws`.
ELBO, log-likelihood y log-joint no son parámetros: se devuelven dentro de
`objective`. Esta separación evita comparar por accidente un hiperparámetro
variacional con una muestra posterior.

## Identificabilidad y label switching

Para cualquier permutación $\pi$, transformar simultáneamente
$(\omega_k,\beta_k,\sigma_k^2,z_j)$ a
$(\omega_{\pi(k)},\beta_{\pi(k)},\sigma_{\pi(k)}^2,\pi(z_j))$
deja la distribución observada sin cambios. Los componentes no son identificables por su número de etiqueta.

La inferencia MCMC ofrece tres capas de resumen:

1. `incidence_matrix`: probabilidad posterior de que cada par de grupos comparta componente; es invariante a etiquetas.
2. `dahl`: una partición posterior representativa obtenida por mínima distancia cuadrática a esa matriz.
3. Relabeling húngaro conjunto respecto a una muestra de referencia para resumir parámetros por componente. La permutación se aplica a todos los parámetros y asignaciones relacionados.

El relabeling facilita interpretación, pero no crea identificabilidad matemática. Si dos componentes se solapan fuertemente, sus parámetros seguirán débilmente identificados. Deben reportarse también co-clustering, incertidumbre de asignación, tamaños efectivos y sensibilidad a $K$ y a la prior.

## Alcance de la validación

Las pruebas automáticas comprueban identidades algebraicas y ejecuciones
cortas. `thesis/final/code/runs/03-simulation/run.R` verifica recuperación con
verdad conocida y `paper/evaluation/` concentra ARI, CPE, entropía,
recuperación alineada, predicción retenida y estrés de label switching para
E1--E5. Estas verificaciones detectan regresiones de código; los perfiles
`smoke` no sustituyen cadenas largas, R-hat/ESS, análisis de sensibilidad ni
validación predictiva con réplicas suficientes para resultados definitivos.

## Escala de los datos

El ajuste puede trabajar en escala original o usar explícitamente
`fit_model(..., scaling = data_scaling(...))`. En este segundo caso se aprende
el mapa sólo con los datos de entrenamiento, se transforman de forma coherente
datos, prior e inicialización, y parámetros, draws, predicciones y densidades
se devuelven en la escala de entrada. En validación cruzada cada fold aprende
su propio mapa, evitando fuga de información. No existe estandarización
silenciosa dentro de los kernels; la decisión y el mapa quedan documentados en
el resultado.

## Parámetros globales

El Gibbs de HLRM y CHLRM conserva las muestras de `Beta`, `Sigma` y `xi2`, además de los parámetros locales o de componente. Esto permite evaluar recuperación y convergencia de toda la jerarquía. VI y SVI reportan los momentos de sus factores globales normal, inverse-Wishart y gamma.

## WAIC y unidad predictiva

Se usa la definición de varianza de WAIC:

$$
\mathrm{lppd}_i=\log\frac1S\sum_s p(y_i\mid\theta_s),\qquad
p_{\mathrm{WAIC},i}=\operatorname{Var}_s\{\log p(y_i\mid\theta_s)\}.
$$

La implementación usa la penalización por varianza y no selecciona punto a punto entre definiciones alternativas. HLRM distingue dos agregaciones condicionales:

- `waic_observation_conditional`: unidad observación, condicional en los parámetros del grupo;
- `waic_group_conditional`: unidad grupo, sumando la log-verosimilitud dentro del grupo antes de calcular WAIC.

La segunda respeta la dependencia dentro del grupo, pero sigue siendo condicional en los efectos locales observados. Predecir grupos completamente nuevos exige validación dejando grupos fuera y marginalizando sus parámetros; no debe confundirse con ninguno de estos dos cálculos.

CHLRM ya no usa esas dos rutas. `predictive_criteria()` devuelve
`criteria$waic_group_marginal`, con una fila por finca e integrando
una sola etiqueta compartida:

$$
\ell_j^{(s)}=\log\sum_k\omega_k^{(s)}\prod_i
N(y_{ij}\mid x_{ij}^{\top}\beta_k^{(s)},\sigma_k^{2(s)}).
$$

DIC fue retirado porque depende de un punto posterior y de la
parametrización; en CHLRM ese punto exige además resolver label switching en
una posterior multimodal. No se conserva como campo vacío ni como ruta de
compatibilidad. WAIC tampoco demuestra convergencia ni recuperación de K.
Los MSE/R² y PPC son condicionales a las fincas observadas y llevan otro
estimando. Por eso no comparten un contenedor GOF: `evaluate_in_sample_fit()`
devuelve `in-sample-fit-1`, `posterior_predictive_check()` devuelve
`posterior-predictive-check-2` y `predictive_criteria()` devuelve
`predictive-criteria-2`. Los nombres son uniformes, pero sólo aparecen las
cantidades definidas para el modelo y objetivo. Esta uniformidad no afirma que
WAIC condicional por observación y WAIC marginal por finca estimen lo mismo.

### Marginalización HLRM para una finca nueva

La prior implementada hace independientes beta_j y la varianza residual v_j
condicionalmente a los parámetros globales. Tras integrar beta_j:

$$
y_j\mid v_j,\beta,\Sigma\sim
N\left(X_j\beta,\;X_j\Sigma X_j^\top+v_j I\right),\qquad
v_j\mid\xi^2\sim IG\left(\nu_0/2,\nu_0\xi^2/2\right).
$$

Por tanto la likelihood marginal usada en `waic_group_marginal` es

$$
p(y_j\mid\beta,\Sigma,\xi^2)=\int_0^\infty
N\left(y_j\mid X_j\beta,X_j\Sigma X_j^\top+vI\right)
IG\left(v\mid\nu_0/2,\nu_0\xi^2/2\right)\,dv.
$$

No es, en general, una t multivariante: Sigma no está multiplicada por v en
esta jerarquía. La implementación cambia a precisión 1/v, utiliza cuadratura
Gamma y verifica el acuerdo entre órdenes crecientes. Conserva diagnóstico
numérico y falla si no converge o supera el presupuesto. El límite Sigma=0
tiene solución analítica, usada en las pruebas, junto con integración adaptativa
independiente para covarianzas no nulas.

`group_marginal=TRUE` activa esta evaluación con draws globales MCMC o q. LRM
puede agrupar su likelihood con `groups`, sin efectos locales que integrar. Comparar
WAIC de las tres familias requiere las mismas fincas, respuestas y escala;
no certifica convergencia ni exactitud de la aproximación variacional.

### Distribución predictiva y PPC

La posterior predictiva integra incertidumbre de parámetros y ruido residual.
En MCMC se aproxima con draws de la cadena; VI/SVI usan q en lugar de la
posterior exacta. Los PPC generan datasets replicados y comparan discrepancias;
no son por sí mismos toda la distribución ni p-valores frecuentistas calibrados.
Las tres familias exponen dos objetos distintos después de `fit()`. `predict()`
resume draws de la media condicional y produce intervalos creíbles puntuales;
`posterior_predict()` muestrea además el ruido residual y produce intervalos
predictivos puntuales y, opcionalmente, `y_rep`. Ninguno es una banda simultánea.

En HLRM una finca conocida reutiliza sus parámetros locales posteriores. Una
finca nueva genera, para cada draw global, una sola pareja
$(\beta_j,\sigma_j^2)$ y la comparte entre sus observaciones. En CHLRM una
finca conocida reutiliza su etiqueta posterior por draw; una nueva genera una
sola etiqueta desde $\omega^{(s)}$. Promediar predicciones después de este
muestreo integra la incertidumbre de asignación; clasificar con `z_pred` y
condicionar en esa clase no es equivalente.

`pp_check_from_prediction()` puede evaluar observaciones contra `y_rep` ya
almacenado. El PPC integrado usa la misma primitiva observacional en modo
streaming para limitar memoria. WAIC no se obtiene de `y_rep`: requiere
la matriz de log-verosimilitud con la unidad predictiva apropiada.

## Simulación generativa

La API común distingue `conditional_dgp()`, que fija parámetros locales o de
componente, de `hierarchical_dgp()`, que genera esos parámetros desde la
jerarquía implementada:

$$
\beta_j\sim N(\beta,\Sigma),\quad
\sigma_j^2=1/u_j,\quad
u_j\sim Gamma(\nu_0/2,\nu_0\xi^2/2),
$$

y el análogo por componentes de CHLRM. El segundo mecanismo es el apropiado para
evaluar cobertura y recuperación de parámetros poblacionales. El diseño se
declara por separado con `simulation_design()`; no hay covariables ni
desplazamientos de grupo ocultos. `prior_predictive_dgp()` genera además todos
los niveles desde la misma prior que recibe el ajuste.

## Diagnósticos y restricciones

`fit_model()` con `mcmc_control(chains >= 2)` calcula R-hat rank-normalizado,
bulk ESS y tail ESS mediante `posterior`. Un único ajuste no puede producir un
R-hat válido. En CHLRM aprende primero una referencia común sobre todas las
cadenas y relabela conjuntamente parámetros y asignaciones; calcular R-hat
sobre parámetros cuyos labels todavía conmutan carece de interpretación. Los
diagnósticos se recalculan después de cada bloque y las cadenas se extienden
desde su estado Gibbs/RNG exacto, sin repetir warm-up ni descartar draws.
`save_mcmc_checkpoint()` conserva un estado plano del `vb_model_fit` completo y
`resume_mcmc()` devuelve un nuevo ajuste; ambas rutas están probadas contra la
continuación exacta de cada cadena.

CAVI/SVI usa el mismo monitor de ELBO en LRM, HLRM y CHLRM: cambio relativo,
paciencia explícita y estado de terminación. `fit_model()` separa la
especificación inicial de la secuencia aleatoria, conserva todos los intentos
multistart y selecciona por ELBO únicamente entre ajustes convergentes y
finitos, salvo una anulación exploratoria explícita. En CHLRM los resultados se
alinean a una única referencia antes de comparar parámetros; el ELBO mismo es
invariante a permutaciones y no necesita relabeling.

Los hiperparámetros OLS siguen siendo empirical Bayes porque reutilizan los datos. Además, `nu0 = 2` induce una inverse-gamma con forma uno y media no finita. Estas decisiones permanecen por compatibilidad con los experimentos existentes, pero requieren análisis de sensibilidad y no deben denominarse automáticamente "no informativas".

Las escalas empirical-Bayes de CHLRM no se multiplican por el número total de observaciones. Bajo `Beta ~ N(Mu0, Lambda0)`, `Lambda0` tiene unidades de covarianza de coeficientes; bajo `Sigma ~ IW(n0, S0)`, `S0` tiene esas mismas unidades. Multiplicar ambas matrices por `n` haría que el prior se volviera más difuso al aumentar la muestra y no se desprende de la jerarquía declarada. HLRM y CHLRM usan la misma convención de escala.

### Contratos matriciales de la implementación

Las matrices que representan covarianzas, precisiones o escalas
inverse-Wishart deben ser simétricas definidas positivas. La implementación
calcula inversas, sistemas lineales y log-determinantes mediante Cholesky. Una
matriz indefinida o materialmente asimétrica se interpreta como violación del
modelo y detiene el ajuste antes de propagar valores inválidos.

La antigua operación de reemplazar siempre una matriz por
$(A+A^\top)/2$ no era un validador: podía ocultar un error en una actualización.
Ahora esa corrección sólo se permite cuando

$$
\max_{r,s}|A_{rs}-A_{sr}| \leq
\varepsilon\left(1+\max_{r,s}|A_{rs}|\right),
$$

con una tolerancia numérica explícita. Después se comprueba definitud positiva
antes de factorizar. Para $\Sigma\sim IW_p(n_0,S_0)$ se valida
$n_0>p-1$, condición de propiedad de la distribución; la existencia de
momentos adicionales debe comprobarse por separado cuando se pretendan
interpretar o reportar.
