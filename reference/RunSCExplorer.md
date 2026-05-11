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
  return_app = TRUE
)
```

## Arguments

- base_dir:

  The base directory of the SCExplorer app. Default is `"SCExplorer"`.

- data_file:

  The name of the HDF5 file that stores data matrices for each dataset.
  Default is `"data.hdf5"`.

- meta_file:

  The name of the HDF5 file that stores metadata for each dataset.
  Default is `"meta.hdf5"`.

- title:

  The title of the SCExplorer app. Default is `"SCExplorer"`.

- initial_dataset:

  The initial dataset to be loaded into the app. Default is `NULL`.

- initial_reduction:

  The initial dimensional reduction method to be loaded into the app.
  Default is `NULL`.

- initial_group:

  The initial variable to group cells in the app. Default is `NULL`.

- initial_feature:

  The initial feature to be loaded into the app. Default is `NULL`.

- initial_assay:

  The initial assay to be loaded into the app. Default is `NULL`.

- initial_slot:

  The initial layer to be loaded into the app. Default is `NULL`.

- initial_label:

  Whether to add labels in the initial plot. Default is `FALSE`.

- initial_cell_palette:

  The initial color palette for cells. Default is `"Chinese"`.

- initial_feature_palette:

  The initial color palette for features. Default is `"Spectral"`.

- initial_theme:

  The initial theme for plots. Default is `"theme_scop"`.

- initial_size:

  The initial size of plots. Default is `4`.

- initial_ncol:

  The initial number of columns for arranging plots. Default is `3`.

- initial_arrange:

  Whether to use "Row" as the initial arrangement. Default is `TRUE`.

- initial_raster:

  Whether to perform rasterization in the initial plot. By default, it
  is set to automatic, meaning it will be `TRUE` if the number of cells
  in the initial datasets exceeds 100,000. Default is `NULL`.

- create_script:

  Whether to create the SCExplorer app script. Default is `TRUE`.

- style_script:

  Whether to style the SCExplorer app script. Default is `TRUE`.

- overwrite:

  Whether to overwrite existing data in the data file. Default is
  `TRUE`.

- return_app:

  Whether to return the SCExplorer app. Default is `TRUE`.

## See also

[CreateDataFile](https://mengxu98.github.io/scop/reference/CreateDataFile.md),
[CreateMetaFile](https://mengxu98.github.io/scop/reference/CreateMetaFile.md),
[PrepareSCExplorer](https://mengxu98.github.io/scop/reference/PrepareSCExplorer.md),
[FetchH5](https://mengxu98.github.io/scop/reference/FetchH5.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2026-05-11 16:18:55] Start standard processing workflow...
#> ℹ [2026-05-11 16:18:56] Checking a list of <Seurat>...
#> ! [2026-05-11 16:18:56] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-11 16:18:56] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-11 16:18:57] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-11 16:18:58] Use the separate HVF from `srt_list`
#> ℹ [2026-05-11 16:18:58] Number of available HVF: 2000
#> ℹ [2026-05-11 16:18:58] Finished check
#> ℹ [2026-05-11 16:18:58] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-11 16:18:58] Perform pca linear dimension reduction
#> ℹ [2026-05-11 16:18:59] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-05-11 16:18:59] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-11 16:18:59] Reorder clusters...
#> ℹ [2026-05-11 16:19:00] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-11 16:19:00] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-11 16:19:00] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-05-11 16:19:05] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-05-11 16:19:09] Standard processing workflow completed
data(panc8_sub)
panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "Harmony"
)
#> ◌ [2026-05-11 16:19:09] Run integration workflow...
#> ℹ [2026-05-11 16:19:10] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-05-11 16:19:10] Checking a list of <Seurat>...
#> ! [2026-05-11 16:19:11] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2026-05-11 16:19:11] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/5 of `srt_list`...
#> ℹ [2026-05-11 16:19:12] Perform `Seurat::FindVariableFeatures()` on 1/5 of `srt_list`...
#> ! [2026-05-11 16:19:12] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2026-05-11 16:19:12] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 2/5 of `srt_list`...
#> ℹ [2026-05-11 16:19:14] Perform `Seurat::FindVariableFeatures()` on 2/5 of `srt_list`...
#> ! [2026-05-11 16:19:14] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2026-05-11 16:19:14] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 3/5 of `srt_list`...
#> ℹ [2026-05-11 16:19:16] Perform `Seurat::FindVariableFeatures()` on 3/5 of `srt_list`...
#> ! [2026-05-11 16:19:16] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2026-05-11 16:19:16] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 4/5 of `srt_list`...
#> ℹ [2026-05-11 16:19:18] Perform `Seurat::FindVariableFeatures()` on 4/5 of `srt_list`...
#> ! [2026-05-11 16:19:18] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2026-05-11 16:19:18] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 5/5 of `srt_list`...
#> ℹ [2026-05-11 16:19:20] Perform `Seurat::FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-05-11 16:19:20] Use the separate HVF from `srt_list`
#> ℹ [2026-05-11 16:19:20] Number of available HVF: 2000
#> ℹ [2026-05-11 16:19:21] Finished check
#> Warning: Layer ‘scale.data’ is empty
#> ℹ [2026-05-11 16:19:23] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-11 16:19:24] Perform linear dimension reduction("pca")
#> ℹ [2026-05-11 16:19:24] Perform Harmony integration
#> ℹ [2026-05-11 16:19:24] Using "Harmonypca" (1:20) as input
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 't': ‘Z_corr’ is not a valid field or method name for reference class “Rcpp_harmony”
panc8_sub <- standard_scop(panc8_sub)
#> ℹ [2026-05-11 16:19:25] Start standard processing workflow...
#> ℹ [2026-05-11 16:19:25] Checking a list of <Seurat>...
#> ! [2026-05-11 16:19:25] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-11 16:19:25] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-11 16:19:27] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-11 16:19:28] Use the separate HVF from `srt_list`
#> ℹ [2026-05-11 16:19:28] Number of available HVF: 2000
#> ℹ [2026-05-11 16:19:28] Finished check
#> ℹ [2026-05-11 16:19:28] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-11 16:19:29] Perform pca linear dimension reduction
#> ℹ [2026-05-11 16:19:29] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-05-11 16:19:30] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-11 16:19:30] Reorder clusters...
#> ℹ [2026-05-11 16:19:30] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-11 16:19:30] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-11 16:19:30] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-05-11 16:19:36] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-05-11 16:19:41] Standard processing workflow completed

PrepareSCExplorer(
  list(
    mouse_pancreas = pancreas_sub,
    human_pancreas = panc8_sub
  ),
  base_dir = "./SCExplorer"
)
#> ℹ [2026-05-11 16:19:41] Prepare data for object: "mouse_pancreas"
#> ℹ [2026-05-11 16:19:41] Write the expression matrix to: /home/runner/work/scop/scop/docs/reference/SCExplorer/data.hdf5
#> ℹ [2026-05-11 16:19:43] Write the meta information to: /home/runner/work/scop/scop/docs/reference/SCExplorer/meta.hdf5
#> ℹ [2026-05-11 16:19:43] Prepare data for object: "human_pancreas"
#> ℹ [2026-05-11 16:19:43] Write the expression matrix to: /home/runner/work/scop/scop/docs/reference/SCExplorer/data.hdf5
#> ℹ [2026-05-11 16:19:48] Write the meta information to: /home/runner/work/scop/scop/docs/reference/SCExplorer/meta.hdf5

# Create the app.R script
app <- RunSCExplorer(
  base_dir = "./SCExplorer",
  initial_dataset = "mouse_pancreas",
  initial_group = "CellType",
  initial_feature = "Ncoa2"
)
#> ℹ [2026-05-11 16:19:48] Create the SCExplorer app script: ./SCExplorer/app.R
#> ℹ [2026-05-11 16:19:48] Styling the script...
#> Loading required package: shiny
#> ✔ [2026-05-11 16:20:04] rhdf5, HDF5Array, shiny, ggplot2, ragg, htmlwidgets, plotly, bslib, promises, and thisplot installed successfully
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
#> Version: 0.3.8 (2026-04-23 update)
#> Website: https://mengxu98.github.io/thisplot/,
#> https://github.com/mengxu98/thisplot
#> 
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(thisplot))
#>   or options(log_message.verbose = FALSE)
#> ------------------------------------------------------------
# Check files
list.files("./SCExplorer")
#> [1] "app.R"     "data.hdf5" "meta.hdf5"

# Run shiny app
if (interactive()) {
  check_r("shiny")
  get_namespace_fun("shiny", "runApp")(app)
}
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
