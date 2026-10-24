# Add expression-program scores with a sparse native fast path

The default Seurat contract is retained. For sparse assay layers,
complete feature lists, and \`k = FALSE\`, scop reproduces Seurat's
binning and control sampling and fuses all score calculations in one
sparse compiled pass. \`search\` and extra arguments are inert once
every requested feature is present. Other branches delegate to Seurat
unchanged.

## Usage

``` r
AddModuleScore(
  object,
  features,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = NULL,
  name = "Cluster",
  seed = 1,
  search = FALSE,
  slot = "data",
  ...
)
```

## Arguments

- object:

  A Seurat object.

- features:

  A list of feature vectors.

- pool, nbin, ctrl, k, assay, name, seed, search, slot:

  Seurat-compatible module-score parameters.

- ...:

  Additional arguments passed to Seurat on fallback.

## Value

The input Seurat object with module-score metadata columns.
