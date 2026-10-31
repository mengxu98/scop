# Run SCExplorer

Run SCExplorer

## Usage

``` r
RunSCExplorer(
  base_dir = "SCExplorer",
  data_file = "data.hdf5",
  meta_file = "meta.hdf5",
  title = "SCExplorer",
  initial_dataset = NULL,
  initial_reduction = NULL,
  initial_group = NULL,
  initial_feature = NULL,
  initial_assay = NULL,
  initial_slot = NULL,
  initial_label = FALSE,
  initial_cell_palette = "Chinese",
  initial_feature_palette = "Spectral",
  initial_theme = "theme_scop",
  initial_size = 4,
  initial_ncol = 3,
  initial_arrange = NULL,
  initial_raster = NULL,
  create_script = TRUE,
  style_script = TRUE,
  overwrite = TRUE,
  return_app = TRUE,
  verbose = TRUE
)
```

## Arguments

- base_dir:

  The base directory of the SCExplorer app.

- data_file:

  HDF5 file that stores data matrices for each dataset.

- meta_file:

  HDF5 file that stores metadata for each dataset.

- title:

  The title of the SCExplorer app.

- initial_dataset:

  The initial dataset to be loaded into the app.

- initial_reduction:

  The initial dimensional reduction method to be loaded into the app.

- initial_group:

  The initial variable to group cells in the app.

- initial_feature:

  The initial feature to be loaded into the app.

- initial_assay:

  The initial assay to be loaded into the app.

- initial_slot:

  The initial layer to be loaded into the app.

- initial_label:

  Whether to add labels in the initial plot.

- initial_cell_palette:

  The initial color palette for cells.

- initial_feature_palette:

  The initial color palette for features.

- initial_theme:

  The initial theme for plots.

- initial_size:

  The initial size of plots.

- initial_ncol:

  The initial number of columns for arranging plots.

- initial_arrange:

  Whether to use "Row" as the initial arrangement.

- initial_raster:

  Whether to perform rasterization in the initial plot. By default, it
  is set to automatic, meaning it will be `TRUE` if the number of cells
  in the initial datasets exceeds 100,000.

- create_script:

  Whether to create the SCExplorer app script.

- style_script:

  Whether to style the SCExplorer app script.

- overwrite:

  Whether to overwrite existing data in the data file.

- return_app:

  Whether to return the SCExplorer app.

- verbose:

  Whether to print the message. Default is `TRUE`.

## See also

[CreateDataFile](https://mengxu98.github.io/scop/reference/CreateDataFile.md),
[CreateMetaFile](https://mengxu98.github.io/scop/reference/CreateMetaFile.md),
[PrepareSCExplorer](https://mengxu98.github.io/scop/reference/PrepareSCExplorer.md),
[FetchH5](https://mengxu98.github.io/scop/reference/FetchH5.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-10-31 18:38:33] Start standard processing workflow...
#> ℹ [2026-10-31 18:38:33] Checking a list of <Seurat>...
#> ! [2026-10-31 18:38:33] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-10-31 18:38:33] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-10-31 18:38:33] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-10-31 18:38:33] Use the separate HVF from `srt_list`
#> ℹ [2026-10-31 18:38:33] Number of available HVF: 2000
#> ℹ [2026-10-31 18:38:33] Finished check
#> ℹ [2026-10-31 18:38:33] Perform `ScaleData()`
#> ℹ [2026-10-31 18:38:33] Perform pca linear dimension reduction
#> ℹ [2026-10-31 18:38:34] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-10-31 18:38:34] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-10-31 18:38:34] Reorder clusters...
#> ℹ [2026-10-31 18:38:34] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-31 18:38:34] Perform umap nonlinear dimension reduction
#> ✔ [2026-10-31 18:38:44] Standard processing workflow completed
data(panc8_sub)
panc8_sub <- RunIntegration(
  panc8_sub,
  batch = "tech",
  integration_methods = "Harmony"
)
#> ◌ [2026-10-31 18:38:44] Run integration workflow...
#> ℹ [2026-10-31 18:38:44] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-10-31 18:38:45] Checking a list of <Seurat>...
#> ! [2026-10-31 18:38:45] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-31 18:38:45] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/5 of `srt_list`...
#> ℹ [2026-10-31 18:38:45] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ! [2026-10-31 18:38:45] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-31 18:38:45] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 2/5 of `srt_list`...
#> ℹ [2026-10-31 18:38:45] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ! [2026-10-31 18:38:45] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-31 18:38:45] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 3/5 of `srt_list`...
#> ℹ [2026-10-31 18:38:45] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ! [2026-10-31 18:38:45] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-31 18:38:45] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 4/5 of `srt_list`...
#> ℹ [2026-10-31 18:38:45] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ! [2026-10-31 18:38:46] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2026-10-31 18:38:46] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 5/5 of `srt_list`...
#> ℹ [2026-10-31 18:38:46] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-10-31 18:38:46] Use the separate HVF from `srt_list`
#> ℹ [2026-10-31 18:38:46] Number of available HVF: 2000
#> ℹ [2026-10-31 18:38:46] Finished check
#> ℹ [2026-10-31 18:38:48] Perform `Seurat::ScaleData()`
#> ℹ [2026-10-31 18:38:49] Perform linear dimension reduction("pca")
#> ℹ [2026-10-31 18:38:50] Perform Harmony integration
#> ℹ [2026-10-31 18:38:50] Using "Harmonypca" (1:20) as input
#> ℹ [2026-10-31 18:38:50] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-10-31 18:38:50] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-10-31 18:38:50] Reorder clusters...
#> ℹ [2026-10-31 18:38:50] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-31 18:38:50] Perform umap nonlinear dimension reduction using Harmony (1:20)
#> ℹ [2026-10-31 18:38:58] Perform umap nonlinear dimension reduction using Harmony (1:20)
#> ℹ [2026-10-31 18:39:05] Perform umap nonlinear dimension reduction using Harmonypca (1:20)
#> ✔ [2026-10-31 18:39:13] Harmony integration completed
panc8_sub <- RunStandardWorkflow(panc8_sub)
#> ℹ [2026-10-31 18:39:13] Start standard processing workflow...
#> ℹ [2026-10-31 18:39:13] Checking a list of <Seurat>...
#> ℹ [2026-10-31 18:39:14] Data 1/1 of the `srt_list` has been log-normalized
#> ℹ [2026-10-31 18:39:14] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-10-31 18:39:14] Use the separate HVF from `srt_list`
#> ℹ [2026-10-31 18:39:14] Number of available HVF: 2000
#> ℹ [2026-10-31 18:39:14] Finished check
#> ℹ [2026-10-31 18:39:14] Perform `ScaleData()`
#> ℹ [2026-10-31 18:39:14] Perform pca linear dimension reduction
#> ℹ [2026-10-31 18:39:15] Use stored estimated dimensions 1:26 for Standardpca
#> ℹ [2026-10-31 18:39:15] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-10-31 18:39:15] Reorder clusters...
#> ℹ [2026-10-31 18:39:15] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-31 18:39:15] Perform umap nonlinear dimension reduction
#> ✔ [2026-10-31 18:39:25] Standard processing workflow completed

PrepareSCExplorer(
  list(
    mouse_pancreas = pancreas_sub,
    human_pancreas = panc8_sub
  ),
  base_dir = paste0(tempdir(), "/SCExplorer")
)
#> ℹ [2026-10-31 18:39:25] Create SCExplorer base directory: /tmp/RtmpxjI8xq/SCExplorer
#> ℹ [2026-10-31 18:39:25] Prepare data for object: "mouse_pancreas"
#> ℹ [2026-10-31 18:39:25] Write the expression matrix to: /tmp/RtmpxjI8xq/SCExplorer/data.hdf5
#> ℹ [2026-10-31 18:39:27] Write the meta information to: /tmp/RtmpxjI8xq/SCExplorer/meta.hdf5
#> ℹ [2026-10-31 18:39:28] Prepare data for object: "human_pancreas"
#> ℹ [2026-10-31 18:39:28] Write the expression matrix to: /tmp/RtmpxjI8xq/SCExplorer/data.hdf5
#> ℹ [2026-10-31 18:39:32] Write the meta information to: /tmp/RtmpxjI8xq/SCExplorer/meta.hdf5

# Create the app.R script
app <- RunSCExplorer(
  base_dir = paste0(tempdir(), "/SCExplorer"),
  initial_dataset = "mouse_pancreas",
  initial_group = "CellType",
  initial_feature = "Ncoa2"
)
#> ℹ [2026-10-31 18:39:32] Create the SCExplorer app script: /tmp/RtmpxjI8xq/SCExplorer/app.R
#> ℹ [2026-10-31 18:39:32] Styling the script...
#> Loading required package: shiny
#> ✔ [2026-10-31 18:39:42] rhdf5, HDF5Array, shiny, ggplot2, ragg, htmlwidgets, plotly, bslib, promises, and thisplot installed successfully
#> 
#> Attaching package: ‘bslib’
#> The following object is masked from ‘package:utils’:
#> 
#>     page
#> 
#> Attaching package: ‘rlang’
#> The following object is masked from ‘package:Biobase’:
#> 
#>     exprs
#>           ⬢          .        ⬡             ⬢     .
#>              __  __    _              __       __
#>             / /_/ /_  (_)_____ ____  / /____  / /_
#>            / __/ __ ./ // ___// __ ./ // __ ./ __/
#>           / /_/ / / / /(__  )/ /_/ / // /_/ / /_
#>           .__/_/ /_/_//____// .___/_/ .____/.__/
#>                            /_/
#>       ⬡               ⬢      .        ⬡          ⬢
#> ------------------------------------------------------------
#> Version: 0.4.6 (2026-10-28 update)
#> Website: https://mengxu98.github.io/thisplot/
#> 
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(thisplot))
#>   or options(log_message.verbose = FALSE)
#> ------------------------------------------------------------
# Check files
list.files(paste0(tempdir(), "/SCExplorer"))
#> [1] "app.R"     "data.hdf5" "meta.hdf5"

# Run shiny app
# shiny::runApp(app)
# Note: If scop installed in the isolated environment using renv,
# add `renv::activate(project = "path/to/scop_env")` to the app.R script.


# You can deploy the app on the self-hosted shiny server
# (https://www.rstudio.com/products/shiny/shiny-server/).
# Or deploy the app on the website
# (https://www.shinyapps.io) for free:

# step1: install "rsconnect" package and authorize your account
# install.packages("rsconnect")
# library(rsconnect)
# setAccountInfo(
#   name = "<NAME>",
#   token = "<TOKEN>",
#   secret = "<SECRET>"
# )

### step2: deploy the app
# deployApp("./SCExplorer")
```
