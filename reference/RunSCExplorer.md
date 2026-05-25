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
#> ℹ [2026-05-25 05:38:36] Start standard processing workflow...
#> ℹ [2026-05-25 05:38:36] Checking a list of <Seurat>...
#> ! [2026-05-25 05:38:37] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-05-25 05:38:37] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-05-25 05:38:39] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-25 05:38:39] Use the separate HVF from `srt_list`
#> ℹ [2026-05-25 05:38:39] Number of available HVF: 2000
#> ℹ [2026-05-25 05:38:39] Finished check
#> ℹ [2026-05-25 05:38:39] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-25 05:38:40] Perform pca linear dimension reduction
#> ℹ [2026-05-25 05:38:40] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-05-25 05:38:41] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-25 05:38:41] Reorder clusters...
#> ℹ [2026-05-25 05:38:41] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-25 05:38:41] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-25 05:38:41] Perform umap nonlinear dimension reduction using Standardpca (1:23)
#> ℹ [2026-05-25 05:38:47] Perform umap nonlinear dimension reduction using Standardpca (1:23)
#> ✔ [2026-05-25 05:38:53] Standard processing workflow completed
data(panc8_sub)
panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "Harmony"
)
#> ◌ [2026-05-25 05:38:53] Run integration workflow...
#> ℹ [2026-05-25 05:38:54] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-05-25 05:38:54] Checking a list of <Seurat>...
#> ! [2026-05-25 05:38:55] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2026-05-25 05:38:55] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/5 of `srt_list`...
#> ℹ [2026-05-25 05:38:57] Perform `Seurat::FindVariableFeatures()` on 1/5 of `srt_list`...
#> ! [2026-05-25 05:38:57] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2026-05-25 05:38:57] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 2/5 of `srt_list`...
#> ℹ [2026-05-25 05:38:59] Perform `Seurat::FindVariableFeatures()` on 2/5 of `srt_list`...
#> ! [2026-05-25 05:38:59] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2026-05-25 05:38:59] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 3/5 of `srt_list`...
#> ℹ [2026-05-25 05:39:01] Perform `Seurat::FindVariableFeatures()` on 3/5 of `srt_list`...
#> ! [2026-05-25 05:39:02] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2026-05-25 05:39:02] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 4/5 of `srt_list`...
#> ℹ [2026-05-25 05:39:04] Perform `Seurat::FindVariableFeatures()` on 4/5 of `srt_list`...
#> ! [2026-05-25 05:39:04] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2026-05-25 05:39:04] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 5/5 of `srt_list`...
#> ℹ [2026-05-25 05:39:06] Perform `Seurat::FindVariableFeatures()` on 5/5 of `srt_list`...
#> ℹ [2026-05-25 05:39:06] Use the separate HVF from `srt_list`
#> ℹ [2026-05-25 05:39:07] Number of available HVF: 2000
#> ℹ [2026-05-25 05:39:07] Finished check
#> Warning: Layer ‘scale.data’ is empty
#> ℹ [2026-05-25 05:39:10] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-25 05:39:10] Perform linear dimension reduction("pca")
#> ℹ [2026-05-25 05:39:11] Perform Harmony integration
#> ℹ [2026-05-25 05:39:11] Using "Harmonypca" (1:20) as input
#> ℹ [2026-05-25 05:39:14] Adjust neighbor k from 20 to 20 for small-sample clustering
#> ℹ [2026-05-25 05:39:14] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-05-25 05:39:14] Reorder clusters...
#> ℹ [2026-05-25 05:39:15] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-25 05:39:15] Perform umap nonlinear dimension reduction using Harmony (1:20)
#> ℹ [2026-05-25 05:39:22] Perform umap nonlinear dimension reduction using Harmony (1:20)
#> ℹ [2026-05-25 05:39:28] Perform umap nonlinear dimension reduction using Harmonypca (1:20)
#> ✔ [2026-05-25 05:39:36] Harmony integration completed
panc8_sub <- standard_scop(panc8_sub)
#> ℹ [2026-05-25 05:39:36] Start standard processing workflow...
#> ℹ [2026-05-25 05:39:36] Checking a list of <Seurat>...
#> ℹ [2026-05-25 05:39:37] Data 1/1 of the `srt_list` has been log-normalized
#> ℹ [2026-05-25 05:39:37] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-05-25 05:39:38] Use the separate HVF from `srt_list`
#> ℹ [2026-05-25 05:39:38] Number of available HVF: 2000
#> ℹ [2026-05-25 05:39:38] Finished check
#> ℹ [2026-05-25 05:39:38] Perform `Seurat::ScaleData()`
#> ℹ [2026-05-25 05:39:39] Perform pca linear dimension reduction
#> ℹ [2026-05-25 05:39:39] Use stored estimated dimensions 1:27 for Standardpca
#> ℹ [2026-05-25 05:39:40] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-05-25 05:39:40] Reorder clusters...
#> ℹ [2026-05-25 05:39:40] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-05-25 05:39:40] Perform umap nonlinear dimension reduction
#> ℹ [2026-05-25 05:39:40] Perform umap nonlinear dimension reduction using Standardpca (1:27)
#> ℹ [2026-05-25 05:39:47] Perform umap nonlinear dimension reduction using Standardpca (1:27)
#> ✔ [2026-05-25 05:39:53] Standard processing workflow completed

PrepareSCExplorer(
  list(
    mouse_pancreas = pancreas_sub,
    human_pancreas = panc8_sub
  ),
  base_dir = "./SCExplorer"
)
#> ℹ [2026-05-25 05:39:53] Prepare data for object: "mouse_pancreas"
#> ℹ [2026-05-25 05:39:53] Write the expression matrix to: /home/runner/work/scop/scop/docs/reference/SCExplorer/data.hdf5
#> ℹ [2026-05-25 05:39:55] Write the meta information to: /home/runner/work/scop/scop/docs/reference/SCExplorer/meta.hdf5
#> ℹ [2026-05-25 05:39:56] Prepare data for object: "human_pancreas"
#> ℹ [2026-05-25 05:39:56] Write the expression matrix to: /home/runner/work/scop/scop/docs/reference/SCExplorer/data.hdf5
#> ℹ [2026-05-25 05:40:00] Write the meta information to: /home/runner/work/scop/scop/docs/reference/SCExplorer/meta.hdf5

# Create the app.R script
app <- RunSCExplorer(
  base_dir = "./SCExplorer",
  initial_dataset = "mouse_pancreas",
  initial_group = "CellType",
  initial_feature = "Ncoa2"
)
#> ℹ [2026-05-25 05:40:35] Create the SCExplorer app script: ./SCExplorer/app.R
#> ℹ [2026-05-25 05:40:35] Styling the script...
#> Loading required package: shiny
#> ✔ [2026-05-25 05:40:45] rhdf5, HDF5Array, shiny, ggplot2, ragg, htmlwidgets, plotly, bslib, promises, and thisplot installed successfully
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
#> Error in value[[3L]](cond): Package ‘thisutils’ version 0.4.5 cannot be unloaded:
#>  Error in unloadNamespace(package) : namespace ‘thisutils’ is imported by ‘thisplot’, ‘scop’ so cannot be unloaded
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
