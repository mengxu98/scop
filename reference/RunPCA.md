# Run principal component analysis

Run principal component analysis

## Usage

``` r
RunPCA(object, ...)
```

## Arguments

- object:

  Object containing expression data.

- ...:

  Passed to methods. The \`"cpp"\` backend of the \`matrix\` and assay
  methods accepts \`cores\` to bound its BLAS calls; \`NULL\` keeps the
  process default. The \`Seurat\` method also accepts \`features\`,
  \`layer\`, \`reduction.name\` and the remaining \`Seurat::RunPCA()\`
  arguments.

## Value

PCA results or an object containing PCA results.
