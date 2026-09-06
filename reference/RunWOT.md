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
  group.by = NULL,
  time_field = "Time",
  growth_iters = 3L,
  tmap_out = "tmaps/tmap_out",
  time_from = NULL,
  time_to = NULL,
  get_coupling = FALSE,
  recalculate = FALSE,
  palette = "Chinese",
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

  A Seurat object. If provided, `adata` will be ignored.

- assay_x:

  Assay to convert as the main data matrix in the anndata object.

- layer_x:

  Layer name for assay_x in the Seurat object.

- assay_y:

  Assays to convert as layers in the anndata object.

- layer_y:

  Layer names for the assay_y in the Seurat object.

- adata:

  An anndata object.

- group.by:

  Metadata column(s) used to color cells.

- time_field:

  Column name in `adata.obs` or `srt@meta.data` that contains the time
  information.

- growth_iters:

  A number of growth iterations to perform during the OT Model
  computation.

- tmap_out:

  Path to store the computed transport maps.

- time_from:

  The starting time point for trajectory and fate analysis.

- time_to:

  The ending time point for trajectory and fate analysis. If not
  provided, only trajectory and fate analysis for the specified
  `time_from` will be performed.

- get_coupling:

  Whether to compute and store the coupling matrix between the specified
  `time_from` and `time_to`.

- recalculate:

  Whether to recalculate the transport maps even if they already exist
  at the specified `tmap_out` location.

- palette, palcolor:

  Palette name
  ([thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html))
  or custom colors.

- show_plot:

  Whether to show the plot.

- save_plot:

  Whether to save plots to files.

- plot_format:

  Format for saved plots: `"png"` (default), `"pdf"`, or `"svg"`.

- plot_dpi:

  Resolution (DPI) for saved plots.

- plot_prefix:

  Prefix for saved plot filenames.

- dirpath:

  The directory to save the plots.

- return_seurat:

  Whether to return a Seurat object instead of an anndata object.

- verbose:

  Whether to print the message. Default is `TRUE`.

## References

[Geoffrey et al. (2019)
Cell](https://doi.org/10.1016/j.cell.2019.01.006),
[GitHub](https://github.com/broadinstitute/wot)

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)

print(range(pancreas_sub$Lineage1, na.rm = TRUE))

pancreas_sub <- RunWOT(
  pancreas_sub,
  group.by = "SubCellType",
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
  group.by = "CellType",
  time_field = "Custom_Time",
  time_from = 1,
  time_to = 10,
  tmap_out = "tmaps/custom_tmap"
)
} # }
```
