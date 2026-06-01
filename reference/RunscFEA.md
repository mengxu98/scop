# Run scFEA flux estimation for a Seurat object

`RunscFEA()` calls the bundled Python scFEA backend through reticulate.
The scFEA GNN architecture, loss function, and training loop are kept in
the backend. M168 human/mouse model files are downloaded from the
`mengxu98/datasets` repository and cached outside the package; the R
side prepares a Seurat expression matrix and stores flux / balance
outputs back into Seurat.

## Usage

``` r
RunscFEA(
  srt,
  assay = NULL,
  layer = "data",
  species = c("human", "mouse"),
  n_epoch = 100,
  sc_imputation = FALSE,
  assay_flux = "scFEAflux",
  assay_balance = "scFEAbalance",
  store_metadata = FALSE,
  data_dir = NULL,
  seed = 16,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  Assay to use as expression matrix. Default is `DefaultAssay(srt)`.

- layer:

  Assay layer to use. Default is `"data"`.

- species:

  One of `"human"` or `"mouse"`, selecting the M168 scFEA files.

- n_epoch:

  Number of scFEA training epochs.

- sc_imputation:

  Whether to run MAGIC imputation inside the scFEA backend.

- assay_flux:

  Name of the assay storing module flux scores.

- assay_balance:

  Name of the assay storing metabolite balance scores.

- store_metadata:

  Whether to also append flux and balance values to `srt@meta.data`.

- data_dir:

  Optional directory containing scFEA M168 CSV resources. If `NULL`,
  files are downloaded from `mengxu98/datasets` and cached with
  `tools::R_user_dir("scop", "data")`.

- seed:

  Random seed passed to R and the Python scFEA backend.

- verbose:

  Whether to print progress messages.

## Value

A Seurat object with `assay_flux`, `assay_balance`, and
`srt@tools[["scFEA"]]`.

## Details

scFEA is licensed for academic, non-commercial use. See the bundled
`inst/python/scfea/LICENSE` file for details. Data resources are
downloaded from <https://github.com/mengxu98/datasets/tree/main/scFEA>.
