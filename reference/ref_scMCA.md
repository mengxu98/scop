# Reference datasets for cell type annotation in single-cell RNA data

Reference datasets for cell type annotation in single-cell RNA data

## Source

[scMCA](https://github.com/ggjlab/scMCA)

## Examples

``` r
if (interactive()) {
  library(Seurat)
  check_r(c("ggjlab/scMCA"))
  ref_scMCA <- NormalizeData(get("ref.expr", envir = asNamespace("scMCA")))
  Encoding(colnames(ref_scMCA)) <- "latin1"
  colnames(ref_scMCA) <- iconv(colnames(ref_scMCA), "latin1", "UTF-8")
  # get_namespace_fun("usethis", "use_data")(ref_scMCA, compress = "xz")
}
```
