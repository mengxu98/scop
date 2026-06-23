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
if (FALSE) { # \dontrun{
# Convert Seurat/SCT data into a standalone Giotto workflow object.
g <- SeuratToScopGiotto(
  spatial_seurat,
  assay = "RNA",
  layer = "counts",
  use_sct = "normalized"
)

# Run the basic Giotto pipeline and plot directly from the Giotto object.
g <- RunGiottoWorkflow(g, steps = "basic")
GiottoPlot(g, plot_type = "cluster")
GiottoPlot(g, plot_type = "network")
GiottoPlot(g, plot_type = "dim")

# Add selected Giotto analyses without writing to the Seurat object.
g <- GiottoSpatialGenes(g, top_n = 50)
g <- GiottoCellProximity(g, group.by = "celltype", number_of_simulations = 100)
GiottoPlot(g, plot_type = "spatial_genes", top_n = 20)
GiottoPlot(g, plot_type = "cell_proximity")

# Export a result back to Seurat only when explicitly requested.
spatial_seurat <- AddGiottoToSeurat(
  spatial_seurat,
  g,
  result = "cluster",
  name = "Giotto_cluster"
)
} # }
```
