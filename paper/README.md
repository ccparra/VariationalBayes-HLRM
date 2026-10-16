# Paper y revisión computacional

Este directorio contiene el segundo producto del repositorio: la evaluación
científica y los artefactos editoriales del paper. No contiene implementaciones
alternativas de LRM, HLRM, CHLRM, MCMC, VI o SVI.

- [`EVALUATION_PLAN.md`](EVALUATION_PLAN.md): matriz de trazabilidad frente al
  PDF de revisión.
- `evaluation/`: configuración, ejecución E1--E5, agregación y reportes.
- `manuscript/`: fuente del artículo y respuesta a revisores.
- `references/`: PDFs usados para fundamentar la revisión.

La evaluación debe instalar el paquete de la raíz y consumir únicamente su API
pública. Los resultados pesados y temporales no se versionan.
