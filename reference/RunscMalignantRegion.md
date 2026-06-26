# Run scMalignantFinder malignant spatial region identification

Use `scMalignantFinder` spatial utilities to score malignant signatures
and infer malignant spatial regions. For Seurat input, provide
`spatial.cols` when spatial-neighborhood refinement is requested and the
converted AnnData object does not already contain spatial coordinates.

## Usage

``` r
RunscMalignantRegion(
  srt = NULL,
  adata = NULL,
  h5ad = NULL,
  assay = "RNA",
  layer = "data",
  cells = NULL,
  signature_gmt,
  features = NULL,
  nclus = 3,
  define_feature = "Malignant_up",
  spatial_nn = TRUE,
  spatial.cols = NULL,
  spatial_key = "spatial",
  image = FALSE,
  norm_type = FALSE,
  prefix = "scMalignantFinder_",
  return_seurat = !is.null(srt),
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- adata:

  Optional Python AnnData object.

- h5ad:

  Optional path to an `.h5ad` file.

- assay:

  Assay used when `srt` is supplied. Default is `"RNA"`.

- layer:

  Layer used when `srt` is supplied. Default is `"data"`.

- cells:

  Optional cells to run. If supplied with `srt`, results are appended to
  these cells and other cells receive `NA`.

- signature_gmt:

  Path to the malignant signature `.gmt` file, such as
  `sc_malignant_deg.gmt` from the scMalignantFinder resources.

- features:

  Features in `adata.obs` used for region clustering. If `NULL`, uses
  `malignancy_probability` and `Malignant_up`, plus `image_score` when
  `image = TRUE`.

- nclus:

  Number of clusters used by the spatial region model.

- define_feature:

  Feature used to select the malignant region cluster.

- spatial_nn:

  Whether to refine labels by spatial neighbors.

- spatial.cols:

  Optional two metadata columns used as spatial coordinates for Seurat
  input.

- spatial_key:

  Key written to `adata.obsm` for spatial coordinates.

- image:

  Whether to call `scMalignantFinder.spatial.image_cal`. Default is
  `FALSE` because Seurat-to-AnnData conversion does not generally carry
  histology images.

- norm_type:

  Passed to `scMalignantFinder`. Use `TRUE` for raw counts that should
  be library-size normalized; use `FALSE` for already normalized input.
  Default is `FALSE`.

- prefix:

  Optional prefix for output metadata columns. Default preserves the
  original `scMalignantFinder` column names.

- return_seurat:

  Whether to return a Seurat object when `srt` is supplied. If `FALSE`,
  returns a data frame of predictions.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A Seurat object with malignant region metadata, or a data frame when
`return_seurat = FALSE`.

## Examples

``` r
if (FALSE) {
srt <- RunscMalignantRegion(
  srt,
  signature_gmt = "path/to/sc_malignant_deg.gmt",
  spatial.cols = c("x", "y")
)
}
```
