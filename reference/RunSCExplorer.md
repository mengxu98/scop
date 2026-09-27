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
#> ℹ [2026-09-27 22:30:20] Start standard processing workflow...
#> ℹ [2026-09-27 22:30:20] Checking a list of <Seurat>...
#> ! [2026-09-27 22:30:20] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-27 22:30:20] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-27 22:30:20] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-27 22:30:20] Use the separate HVF from `srt_list`
#> ℹ [2026-09-27 22:30:20] Number of available HVF: 2000
#> ℹ [2026-09-27 22:30:20] Finished check
#> ℹ [2026-09-27 22:30:20] Perform `ScaleData()`
#> ℹ [2026-09-27 22:30:20] Perform pca linear dimension reduction
#> ℹ [2026-09-27 22:30:20] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-27 22:30:20] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-27 22:30:20] Reorder clusters...
#> ℹ [2026-09-27 22:30:20] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-27 22:30:20] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-27 22:30:28] Standard processing workflow completed
data(panc8_sub)
panc8_sub <- RunIntegration(
  panc8_sub,
  batch = "tech",
  integration_methods = "Harmony"
)
#> ◌ [2026-09-27 22:30:28] Run integration workflow...
#> ℹ [2026-09-27 22:30:29] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-09-27 22:30:29] Checking a list of <Seurat>...
#> ! [2026-09-27 22:30:29] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2026-09-27 22:30:29] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/5 of `srt_list`...
#> ℹ [2026-09-27 22:30:29] Perform `FindVariableFeatures()` on 1/5 of `srt_list`...
#> ! [2026-09-27 22:30:29] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2026-09-27 22:30:29] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 2/5 of `srt_list`...
#> ℹ [2026-09-27 22:30:29] Perform `FindVariableFeatures()` on 2/5 of `srt_list`...
#> ! [2026-09-27 22:30:29] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2026-09-27 22:30:29] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 3/5 of `srt_list`...
#> ℹ [2026-09-27 22:30:29] Perform `FindVariableFeatures()` on 3/5 of `srt_list`...
#> ! [2026-09-27 22:30:29] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2026-09-27 22:30:29] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 4/5 of `srt_list`...
#> ℹ [2026-09-27 22:30:30] Perform `FindVariableFeatures()` on 4/5 of `srt_list`...
#> ! [2026-09-27 22:30:30] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2026-09-27 22:30:30] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 5/5 of `srt_list`...
#> ℹ [2026-09-27 22:30:30] Perform `FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-09-27 22:30:30] Use the separate HVF from `srt_list`
#> ℹ [2026-09-27 22:30:30] Number of available HVF: 2000
#> ℹ [2026-09-27 22:30:30] Finished check
#> ℹ [2026-09-27 22:30:31] Perform `Seurat::ScaleData()`
#> ℹ [2026-09-27 22:30:31] Perform linear dimension reduction("pca")
#> ℹ [2026-09-27 22:30:32] Perform Harmony integration
#> ℹ [2026-09-27 22:30:32] Using "Harmonypca" (1:20) as input
#> ℹ [2026-09-27 22:30:32] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-09-27 22:30:32] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-09-27 22:30:32] Reorder clusters...
#> ℹ [2026-09-27 22:30:32] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-27 22:30:32] Perform umap nonlinear dimension reduction using Harmony (1:20)
#> ℹ [2026-09-27 22:30:39] Perform umap nonlinear dimension reduction using Harmony (1:20)
#> ℹ [2026-09-27 22:30:45] Perform umap nonlinear dimension reduction using Harmonypca (1:20)
#> ✔ [2026-09-27 22:30:51] Harmony integration completed
panc8_sub <- RunStandardWorkflow(panc8_sub)
#> ℹ [2026-09-27 22:30:51] Start standard processing workflow...
#> ℹ [2026-09-27 22:30:51] Checking a list of <Seurat>...
#> ℹ [2026-09-27 22:30:52] Data 1/1 of the `srt_list` has been log-normalized
#> ℹ [2026-09-27 22:30:52] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-27 22:30:52] Use the separate HVF from `srt_list`
#> ℹ [2026-09-27 22:30:52] Number of available HVF: 2000
#> ℹ [2026-09-27 22:30:52] Finished check
#> ℹ [2026-09-27 22:30:52] Perform `ScaleData()`
#> ℹ [2026-09-27 22:30:52] Perform pca linear dimension reduction
#> ℹ [2026-09-27 22:30:52] Use stored estimated dimensions 1:26 for Standardpca
#> ℹ [2026-09-27 22:30:53] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-27 22:30:53] Reorder clusters...
#> ℹ [2026-09-27 22:30:53] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-27 22:30:53] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-27 22:31:01] Standard processing workflow completed

PrepareSCExplorer(
  list(
    mouse_pancreas = pancreas_sub,
    human_pancreas = panc8_sub
  ),
  base_dir = paste0(tempdir(), "/SCExplorer")
)
#> ℹ [2026-09-27 22:31:01] Create SCExplorer base directory: /tmp/RtmpNIRZ7c/SCExplorer
#> ℹ [2026-09-27 22:31:01] Prepare data for object: "mouse_pancreas"
#> ℹ [2026-09-27 22:31:01] Write the expression matrix to: /tmp/RtmpNIRZ7c/SCExplorer/data.hdf5
#> ℹ [2026-09-27 22:31:03] Write the meta information to: /tmp/RtmpNIRZ7c/SCExplorer/meta.hdf5
#> ℹ [2026-09-27 22:31:03] Prepare data for object: "human_pancreas"
#> ℹ [2026-09-27 22:31:03] Write the expression matrix to: /tmp/RtmpNIRZ7c/SCExplorer/data.hdf5
#> ℹ [2026-09-27 22:31:05] Write the meta information to: /tmp/RtmpNIRZ7c/SCExplorer/meta.hdf5

# Create the app.R script
app <- RunSCExplorer(
  base_dir = paste0(tempdir(), "/SCExplorer"),
  initial_dataset = "mouse_pancreas",
  initial_group = "CellType",
  initial_feature = "Ncoa2"
)
#> ℹ [2026-09-27 22:31:07] Create the SCExplorer app script: /tmp/RtmpNIRZ7c/SCExplorer/app.R
#> ℹ [2026-09-27 22:31:07] Styling the script...
#> Loading required package: shiny
#> ✔ [2026-09-27 22:31:11] rhdf5, HDF5Array, shiny, ggplot2, ragg, htmlwidgets, plotly, bslib, promises, and thisplot installed successfully
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
#> Version: 0.4.6 (2026-09-26 update)
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
