# Local build

The thesis uses `EncabezadoTesisMSc.tex` as its root document.

Requirements:

- Tectonic 0.17 or newer.
- `make` (optional).

Compile from this directory with:

```sh
make
```

or directly with:

```sh
tectonic -X compile EncabezadoTesisMSc.tex --outdir ../document
```

The generated PDF is written to `thesis/final/document/`, keeping the source
directory free of compilation artifacts. Tectonic downloads its LaTeX bundle on the
first build and caches it for later offline builds.
