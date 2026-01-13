# Run WOT analysis

Run WOT analysis

## Usage

``` r
RunWOT(
  srt = NULL,
  assay_x = "RNA",
  layer_x = "counts",
  assay_y = c("spliced", "unspliced"),
  layer_y = "counts",
  adata = NULL,
  group_by = NULL,
  time_field = "Time",
  growth_iters = 3L,
  tmap_out = "tmaps/tmap_out",
  time_from = NULL,
  time_to = NULL,
  get_coupling = FALSE,
  recalculate = FALSE,
  palette = "Paired",
  palcolor = NULL,
  show_plot = FALSE,
  save_plot = FALSE,
  plot_format = c("pdf", "png", "svg"),
  plot_dpi = 300,
  plot_prefix = "wot",
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

- time_field:

  A character string specifying the column name in `adata.obs` or
  `srt@meta.data` that contains the time information.

- growth_iters:

  A number of growth iterations to perform during the OT Model
  computation. Default is `3`.

- tmap_out:

  A character string specifying the path to store the computed transport
  maps.

- time_from:

  The starting time point for trajectory and fate analysis.

- time_to:

  The ending time point for trajectory and fate analysis. If not
  provided, only trajectory and fate analysis for the specified
  `time_from` will be performed.

- get_coupling:

  Whether to compute and store the coupling matrix between the specified
  `time_from` and `time_to`. Default is `FALSE`.

- recalculate:

  Whether to recalculate the transport maps even if they already exist
  at the specified `tmap_out` location. Default is `FALSE`.

- palette:

  Color palette name. Available palettes can be found in
  [thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html).
  Default is `"Paired"`.

- palcolor:

  Custom colors used to create a color palette. Default is `NULL`.

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

## References

[Geoffrey et al. (2019)
Cell](https://doi.org/10.1016/j.cell.2019.01.006),
[GitHub](https://github.com/broadinstitute/wot)

## See also

[srt_to_adata](https://mengxu98.github.io/scop/reference/srt_to_adata.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)

print(range(pancreas_sub$Lineage1, na.rm = TRUE))

pancreas_sub <- RunWOT(
  pancreas_sub,
  group_by = "SubCellType",
  time_field = "Lineage1",
  time_from = min(pancreas_sub$Lineage1, na.rm = TRUE),
  time_to = max(pancreas_sub$Lineage1, na.rm = TRUE),
  get_coupling = TRUE,
  tmap_out = "tmaps/lineage_tmap"
)

pancreas_sub$Custom_Time <- sample(
  1:10,
  ncol(pancreas_sub),
  replace = TRUE
)
pancreas_sub <- RunWOT(
  pancreas_sub,
  group_by = "CellType",
  time_field = "Custom_Time",
  time_from = 1,
  time_to = 10,
  tmap_out = "tmaps/custom_tmap"
)
} # }
```
