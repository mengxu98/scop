# Convert Seurat to an internal Giotto workflow object

Create a \`giotto2\` object from a Seurat object. The converter is
SCT-aware: raw counts remain the default Giotto input, while SCT
normalized values are optionally added as an extra Giotto expression
layer. The input Seurat object is not modified.

## Usage

``` r
SeuratToScopGiotto(
  srt,
  assay = NULL,
  layer = "counts",
  sct.assay = "SCT",
  use_sct = c("auto", "none", "normalized"),
  image = NULL,
  coord.cols = c("x", "y"),
  features = NULL,
  conversion_params = list(),
  use_official = TRUE,
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- layer:

  Assay layer used as the expression matrix.

- sct.assay:

  Name of the SCT assay.

- use_sct:

  How to handle SCT data. \`"auto"\` keeps counts as the main Giotto
  expression and records SCT availability. \`"none"\` ignores SCT.
  \`"normalized"\` adds SCT normalized values as an additional
  expression layer.

- image:

  Name of the Seurat spatial image used by the spatial workflow. If
  `NULL`, the first image is used when present.

- coord.cols:

  Metadata coordinate columns used by the spatial workflow when no image
  is available.

- features:

  Features used for PCA and clustering. If `NULL`, current variable
  features are used, falling back to all assay features.

- conversion_params:

  Additional parameters passed to `Giotto::createGiottoObject()`.

- use_official:

  Whether to try \`Giotto::seuratToGiottoV5()\` before falling back to
  the scop-controlled converter.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed for reproducibility. Default is `11`.

## Value

A \`giotto2\` object.

## Examples

``` r
data(visium_human_pancreas_sub)
spatial <- subset(
  visium_human_pancreas_sub,
  cells = colnames(visium_human_pancreas_sub)[1:80],
  features = rownames(visium_human_pancreas_sub)[1:300]
)
#> Warning: Not validating Centroids objects
#> Warning: Not validating Centroids objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating FOV objects
#> Warning: Not validating Seurat objects
g <- structure(
  list(
    giotto = list(
      umap = cbind(
        UMAP_1 = as.numeric(scale(spatial$x)),
        UMAP_2 = as.numeric(scale(spatial$y))
      )
    ),
    source = list(
      cells = colnames(spatial),
      features = rownames(spatial),
      coordinates = data.frame(
        cell_ID = colnames(spatial),
        sdimx = spatial$x,
        sdimy = spatial$y
      )
    ),
    results = list(
      cluster = list(
        table = data.frame(
          cluster = paste0("cluster_", (seq_len(ncol(spatial)) - 1) %% 3 + 1),
          row.names = colnames(spatial)
        )
      ),
      spatial_network = list(
        table = data.frame(
          from = colnames(spatial)[1:8],
          to = colnames(spatial)[2:9]
        )
      )
    ),
    active = "cluster"
  ),
  class = c("giotto2", "list")
)

GiottoPlot(g, plot_type = "cluster")

GiottoPlot(g, plot_type = "network")


if (
  requireNamespace("Giotto", quietly = TRUE) &&
    identical(Sys.getenv("SCOP_RUN_SPATIAL_BACKEND_EXAMPLES"), "true")
) {
g <- SeuratToScopGiotto(
  spatial,
  assay = "Spatial",
  layer = "counts",
  coord.cols = c("x", "y"),
  verbose = FALSE
)
}
```
