# Reference datasets for cell type annotation in single-cell RNA data

Reference datasets for cell type annotation in single-cell RNA data

## Usage

``` r
ref_scHCL

ref_scMCA

ref_scZCL
```

## Source

[scHCL](https://github.com/ggjlab/scHCL),
[scMCA](https://github.com/ggjlab/scMCA),
[scZCL](https://github.com/ggjlab/scZCL)

## Examples

``` r
if (FALSE) { # \dontrun{
if (interactive()) {
  library(Seurat)
  check_r(c("ggjlab/scZCL", "ggjlab/scHCL", "ggjlab/scMCA"))
  ref_scHCL <- NormalizeData(scHCL::ref.expr)
  ref_scMCA <- NormalizeData(scMCA::ref.expr)
  ref_scZCL <- NormalizeData(scZCL::ref.expr)
  Encoding(colnames(ref_scHCL)) <- "latin1"
  colnames(ref_scHCL) <- iconv(colnames(ref_scHCL), "latin1", "UTF-8")
  Encoding(colnames(ref_scMCA)) <- "latin1"
  colnames(ref_scMCA) <- iconv(colnames(ref_scMCA), "latin1", "UTF-8")
  Encoding(colnames(ref_scZCL)) <- "latin1"
  colnames(ref_scZCL) <- iconv(colnames(ref_scZCL), "latin1", "UTF-8")
  # usethis::use_data(ref_scHCL, compress = "xz")
  # usethis::use_data(ref_scMCA, compress = "xz")
  # usethis::use_data(ref_scZCL, compress = "xz")
}
} # }
```
