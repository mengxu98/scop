# Run PAGA analysis

PAGA is a graph-based method used to infer cellular trajectories. This
function runs the PAGA analysis on a Seurat object.

## Usage

``` r
RunPAGA(
  srt = NULL,
  adata = NULL,
  assay_x = "RNA",
  layer_x = "counts",
  assay_y = c("spliced", "unspliced"),
  layer_y = "counts",
  group_by = NULL,
  linear_reduction = NULL,
  nonlinear_reduction = NULL,
  basis = NULL,
  n_pcs = 30,
  n_neighbors = 30,
  use_rna_velocity = FALSE,
  vkey = "stochastic",
  embedded_with_PAGA = FALSE,
  paga_layout = "fr",
  threshold = 0.1,
  point_size = 20,
  infer_pseudotime = FALSE,
  root_group = NULL,
  root_cell = NULL,
  n_dcs = 10,
  n_branchings = 0,
  min_group_size = 0.01,
  palette = "Paired",
  palcolor = NULL,
  legend.position = "on data",
  cores = 1,
  show_plot = FALSE,
  save_plot = FALSE,
  plot_format = c("pdf", "png", "svg"),
  plot_dpi = 300,
  plot_prefix = "paga",
  dirpath = "./paga",
  return_seurat = !is.null(srt),
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object. Default is `NULL`. If provided, `adata` will be
  ignored.

- adata:

  An anndata object. Default is `NULL`.

- assay_x:

  Assay to convert in the anndata object.

- layer_x:

  Layer name for `assay_x` in the Seurat object.

- assay_y:

  Assay to convert in the anndata object.

- layer_y:

  Layer names for the `assay_y` in the Seurat object.

- group_by:

  Variable to use for grouping cells in the Seurat object.

- linear_reduction:

  Linear reduction method to use, e.g., `"PCA"`.

- nonlinear_reduction:

  Non-linear reduction method to use, e.g., `"UMAP"`.

- basis:

  The basis to use for reduction, e.g., `"UMAP"`.

- n_pcs:

  Number of principal components to use for linear reduction. Default is
  `30`.

- n_neighbors:

  Number of neighbors to use for constructing the KNN graph. Default is
  `30`.

- use_rna_velocity:

  Whether to use RNA velocity for PAGA analysis. Default is `FALSE`.

- vkey:

  The name of the RNA velocity data to use if `use_rna_velocity` is
  `TRUE`. Default is `"stochastic"`.

- embedded_with_PAGA:

  Whether to embed data using PAGA layout. Default is `FALSE`.

- paga_layout:

  The layout for plotting PAGA graph. See
  [layout](https://scanpy.readthedocs.io/en/stable/tutorials/plotting/advanced.html#paga)
  param in `scanpy.pl.paga` function.

- threshold:

  The threshold for plotting PAGA graph. Edges for weights below this
  threshold will not be drawn.

- point_size:

  The point size for plotting.

- infer_pseudotime:

  Whether to infer pseudotime.

- root_group:

  The group to use as the root for pseudotime inference.

- root_cell:

  The cell to use as the root for pseudotime inference.

- n_dcs:

  The number of diffusion components to use for pseudotime inference.

- n_branchings:

  Number of branchings to detect.

- min_group_size:

  The minimum size of a group (as a fraction of the total number of
  cells) to consider it as a potential branching point.

- palette:

  The palette to use for coloring cells.

- palcolor:

  A vector of colors to use as the palette.

- legend.position:

  Position of legend in plots. Can be `"on data"`, `"right margin"`,
  `"bottom right"`, etc. Default is `"on data"`.

- cores:

  The number of cores to use for parallelization with
  [foreach::foreach](https://rdrr.io/pkg/foreach/man/foreach.html).
  Default is `1`.

- show_plot:

  Whether to show the plot. Default is `FALSE`.

- save_plot:

  Whether to save plots to files. Default is `FALSE`.

- plot_format:

  Format for saved plots: `"png"` (default), `"pdf"`, or `"svg"`.

- plot_dpi:

  Resolution (DPI) for saved plots. Default is `300`.

- plot_prefix:

  Prefix for saved plot filenames. Default is "cellrank".

- dirpath:

  The directory to save the plots. Default is `"./cellrank"`.

- return_seurat:

  Whether to return a Seurat object instead of an anndata object.
  Default is `TRUE`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## See also

[srt_to_adata](https://mengxu98.github.io/scop/reference/srt_to_adata.md),
[PAGAPlot](https://mengxu98.github.io/scop/reference/PAGAPlot.md),
[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md),
[RunSCVELO](https://mengxu98.github.io/scop/reference/RunSCVELO.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
pancreas_sub <- RunPAGA(
  pancreas_sub,
  assay_x = "RNA",
  group_by = "SubCellType",
  linear_reduction = "PCA",
  nonlinear_reduction = "UMAP"
)
CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "draw_graph_fr"
)

PAGAPlot(pancreas_sub, reduction = "UMAP")

CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP",
  paga = pancreas_sub@misc$paga
)

pancreas_sub <- RunPAGA(
  pancreas_sub,
  group_by = "SubCellType",
  linear_reduction = "PCA",
  nonlinear_reduction = "UMAP",
  embedded_with_PAGA = TRUE,
  infer_pseudotime = TRUE,
  root_group = "Ductal"
)

FeatureDimPlot(
  pancreas_sub,
  features = "dpt_pseudotime",
  reduction = "PAGAUMAP2D"
)

PAGAPlot(pancreas_sub, reduction = "PAGAUMAP2D")

CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "PAGAUMAP2D",
  paga = pancreas_sub@misc$paga
)
} # }
```
