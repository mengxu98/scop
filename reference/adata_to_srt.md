# Convert an anndata object to a seurat object using reticulate

Convert an anndata object to a seurat object using reticulate

## Usage

``` r
adata_to_srt(adata, verbose = TRUE)
```

## Arguments

- adata:

  A connected python anndata object.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
adata <- srt_to_adata(pancreas_sub)
adata <- RunPAGA(
  adata = adata,
  group_by = "SubCellType",
  linear_reduction = "X_pca",
  nonlinear_reduction = "X_umap"
)
srt <- adata_to_srt(adata)
srt

# Or convert a h5ad file to Seurat object
sc <- reticulate::import("scanpy")
adata <- sc$read_h5ad("pancreas.h5ad")
srt <- adata_to_srt(adata)
srt
} # }
```
