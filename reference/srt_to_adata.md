# Convert a Seurat object to an AnnData object

This function takes a Seurat object and converts it to an anndata object
using the reticulate package.

## Usage

``` r
srt_to_adata(
  srt,
  features = NULL,
  assay_x = "RNA",
  layer_x = "counts",
  assay_y = c("spliced", "unspliced"),
  layer_y = "counts",
  convert_tools = FALSE,
  convert_misc = FALSE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- features:

  Optional vector of features to include in the anndata object. Default
  is all features in `assay_x`.

- assay_x:

  Assay to convert as the main data matrix in the anndata object.
  Default is `"RNA"`.

- layer_x:

  Layer name for assay_x in the Seurat object. Default is `"counts"`.

- assay_y:

  Assays to convert as layers in the anndata object. Default is
  `c("spliced", "unspliced")`.

- layer_y:

  Layer names for the assay_y in the Seurat object. Default is
  `"counts"`.

- convert_tools:

  Whether to convert the tool-specific data. Default is `FALSE`.

- convert_misc:

  Whether to convert the miscellaneous data. Default is `FALSE`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A `anndata` object.

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
adata <- srt_to_adata(pancreas_sub)
adata

# Or save as a h5ad/loom file
adata$write_h5ad(
  "pancreas_sub.h5ad"
)
adata$write_loom(
  "pancreas_sub.loom",
  write_obsm_varm = TRUE
)
} # }
```
