# Convert an anndata object to a seurat object

Convert an anndata object to a seurat object

## Usage

``` r
adata_to_srt(adata, verbose = TRUE)
```

## Arguments

- adata:

  An AnnData object. Can be a Python AnnData object (from
  `scanpy`/``` reticulate``), an  ```AnnDataR6`object from the`anndata`package, or an`InMemoryAnnData`object from the`anndataR\`
  package.

- verbose:

  Whether to print the message. Default is `TRUE`.

## See also

[srt_to_adata](https://mengxu98.github.io/scop/reference/srt_to_adata.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
adata <- srt_to_adata(pancreas_sub)
adata <- RunPAGA(
  adata = adata,
  group.by = "SubCellType",
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
