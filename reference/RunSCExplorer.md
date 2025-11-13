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
  initial_cell_palette = "Paired",
  initial_feature_palette = "Spectral",
  initial_theme = "theme_scop",
  initial_size = 4,
  initial_ncol = 3,
  initial_arrange = NULL,
  initial_raster = NULL,
  create_script = TRUE,
  style_script = requireNamespace("styler", quietly = TRUE),
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

  The initial color palette for cells. Default is `"Paired"`.

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

  Whether to overwrite existing files. Default is `TRUE`.

- return_app:

  Whether to return the SCExplorer app. Default is `TRUE`.

## See also

[CreateDataFile](https://mengxu98.github.io/scop/reference/CreateDataFile.md),
[CreateMetaFile](https://mengxu98.github.io/scop/reference/CreateMetaFile.md),
[PrepareSCExplorer](https://mengxu98.github.io/scop/reference/PrepareSCExplorer.md),
[FetchH5](https://mengxu98.github.io/scop/reference/FetchH5.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
data(panc8_sub)
panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "Seurat"
)

PrepareSCExplorer(
  list(
    mouse_pancreas = pancreas_sub,
    human_pancreas = panc8_sub
  ),
  base_dir = "./SCExplorer"
)

# Create the app.R script
app <- RunSCExplorer(
  base_dir = "./SCExplorer",
  initial_dataset = "mouse_pancreas",
  initial_group = "CellType",
  initial_feature = "Ncoa2"
)
# Check files
list.files("./SCExplorer")

# Run shiny app
if (interactive()) {
  shiny::runApp(app)
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
} # }
```
