# Run Palantir analysis

Run Palantir analysis

## Usage

``` r
RunPalantir(
  srt = NULL,
  assay_x = "RNA",
  layer_x = "counts",
  assay_y = c("spliced", "unspliced"),
  layer_y = "counts",
  adata = NULL,
  group_by = NULL,
  linear_reduction = NULL,
  nonlinear_reduction = NULL,
  basis = NULL,
  n_pcs = 30,
  n_neighbors = 30,
  dm_n_components = 10,
  dm_alpha = 0,
  dm_n_eigs = NULL,
  early_group = NULL,
  early_cell = NULL,
  terminal_cells = NULL,
  terminal_groups = NULL,
  num_waypoints = 1200,
  scale_components = TRUE,
  use_early_cell_as_start = TRUE,
  adjust_early_cell = FALSE,
  adjust_terminal_cells = FALSE,
  max_iterations = 25,
  cores = 1,
  point_size = 20,
  palette = "Paired",
  palcolor = NULL,
  legend.position = "on data",
  show_plot = FALSE,
  save_plot = FALSE,
  plot_format = c("pdf", "png", "svg"),
  plot_dpi = 300,
  plot_prefix = "palantir",
  dirpath = "./",
  return_seurat = !is.null(srt),
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object. Default is `NULL`. If provided, `adata` will be
  ignored.

- assay_x:

  Assay to convert in the anndata object.

- layer_x:

  Layer name for `assay_x` in the Seurat object.

- assay_y:

  Assay to convert in the anndata object.

- layer_y:

  Layer names for the `assay_y` in the Seurat object.

- adata:

  An anndata object. Default is `NULL`.

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

- dm_n_components:

  The number of diffusion components to calculate.

- dm_alpha:

  Normalization parameter for the diffusion operator.

- dm_n_eigs:

  Number of eigen vectors to use.

- early_group:

  Name of the group to start Palantir analysis from.

- early_cell:

  Name of the cell to start Palantir analysis from.

- terminal_cells:

  Character vector specifying terminal cells for Palantir analysis.

- terminal_groups:

  Character vector specifying terminal groups for Palantir analysis.

- num_waypoints:

  Number of waypoints to be included.

- scale_components:

  Should the cell fate probabilities be scaled for each component
  independently?

- use_early_cell_as_start:

  Should the starting cell for each terminal group be set as early_cell?

- adjust_early_cell:

  Whether to adjust the early cell to the cell with the minimum
  pseudotime value.

- adjust_terminal_cells:

  Whether to adjust the terminal cells to the cells with the maximum
  pseudotime value for each terminal group.

- max_iterations:

  Maximum number of iterations for pseudotime convergence.

- cores:

  The number of cores to use for parallelization with
  [foreach::foreach](https://rdrr.io/pkg/foreach/man/foreach.html).
  Default is `1`.

- point_size:

  The point size for plotting.

- palette:

  Color palette name. Available palettes can be found in
  [thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html).
  Default is `"Paired"`.

- palcolor:

  Custom colors used to create a color palette. Default is `NULL`.

- legend.position:

  Position of legend in plots. Can be `"on data"`, `"right margin"`,
  `"bottom right"`, etc. Default is `"on data"`.

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

[srt_to_adata](https://mengxu98.github.io/scop/reference/srt_to_adata.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
pancreas_sub <- RunPalantir(
  pancreas_sub,
  group_by = "SubCellType",
  linear_reduction = "PCA",
  nonlinear_reduction = "UMAP",
  early_group = "Ductal",
  terminal_groups = c("Alpha", "Beta", "Delta", "Epsilon")
)

FeatureDimPlot(
  pancreas_sub,
  c("palantir_pseudotime", "palantir_diff_potential")
)

FeatureDimPlot(
  pancreas_sub,
  paste0(
    c("Alpha", "Beta", "Delta", "Epsilon"),
    "_diff_potential"
  )
)
} # }
```
