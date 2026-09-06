# Convert a Seurat object to an AnnData object

Convert a Seurat object to an AnnData object

## Usage

``` r
srt_to_adata(
  srt,
  features = NULL,
  assay_x = "RNA",
  layer_x = "counts",
  assay_y = c("spliced", "unspliced"),
  layer_y = "counts",
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  convert_tools = FALSE,
  convert_misc = FALSE,
  prepare_env = TRUE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A `Seurat` object.

- features:

  Optional vector of features to include in the anndata object. Default
  is all features in `assay_x`.

- assay_x:

  Assay to convert as the main data matrix in the anndata object.

- layer_x:

  Layer name for assay_x in the Seurat object.

- assay_y:

  Assays to convert as layers in the anndata object.

- layer_y:

  Layer names for the assay_y in the Seurat object.

- reductions:

  Character vector specifying which Seurat reductions to convert into
  `obsm`. Default is `NULL`, which converts all available reductions.

- graphs:

  Character vector specifying which Seurat graphs to convert into
  `obsp`. Default is `NULL`, which converts all available graphs.

- neighbors:

  Character vector specifying which Seurat neighbor objects to convert
  into `obsp`. Default is `NULL`, which converts all available neighbor
  objects.

- convert_tools:

  Whether to convert the tool-specific data.

- convert_misc:

  Whether to convert the miscellaneous data.

- prepare_env:

  Whether to prepare and validate the Scanpy Python environment before
  conversion. Wrapper functions that already prepared the environment
  should pass `FALSE` to keep one reticulate binding.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A `anndata` object.

## See also

[adata_to_srt](https://mengxu98.github.io/scop/reference/adata_to_srt.md)

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
