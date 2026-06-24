# Get expression data from `Assay5` or Seurat object

A re-implementation of the
[SeuratObject::GetAssayData](https://satijalab.github.io/seurat-object/reference/AssayData.html)
function to compatible with Assay5 objects.

## Usage

``` r
GetAssayData5(object, ...)

# S3 method for class 'Seurat'
GetAssayData5(object, layer = "counts", assay = NULL, ...)

# S3 method for class 'Assay5'
GetAssayData5(object, layer = "counts", ...)

# S3 method for class 'Assay'
GetAssayData5(object, layer = "counts", ...)
```

## Arguments

- object:

  An object

- ...:

  Additional arguments passed to
  [SeuratObject::GetAssayData](https://satijalab.github.io/seurat-object/reference/AssayData.html).

- layer:

  Name of layer to get or set

- assay:

  Specific assay to get data from or set data for; defaults to the
  [default
  assay](https://satijalab.github.io/seurat-object/reference/DefaultAssay.html)

## Value

A matrix or data frame containing the assay data.

## See also

[SeuratObject::GetAssayData](https://satijalab.github.io/seurat-object/reference/AssayData.html)

## Examples

``` r
data(pancreas_sub)
GetAssayData5(
  pancreas_sub,
  layer = "counts",
  assay = "RNA"
)[1:5, 1:5]
#> 5 x 5 sparse Matrix of class "dgCMatrix"
#>               AAACCTGAGCCTTGAT AAACCTGGTAAGTGGC AAACGGGAGATATGGT
#> Xkr4                         .                .                .
#> Mrpl15                       .                .                .
#> Npbwr1                       .                .                .
#> 4732440D04Rik                .                .                .
#> Gm26901                      1                .                .
#>               AAACGGGCAAAGAATC AAACGGGGTACAGTTC
#> Xkr4                         .                .
#> Mrpl15                       .                .
#> Npbwr1                       .                .
#> 4732440D04Rik                .                .
#> Gm26901                      .                .

data(panc8_sub)
GetAssayData5(
  panc8_sub,
  layer = "counts",
  assay = "RNA"
)[1:5, 1:5]
#> 5 x 5 sparse Matrix of class "dgCMatrix"
#>         D101_13 D101_22 D101_43  D101_50   D102_7
#> A1CF   3.017717       .       . .        .       
#> A4GALT .              .       . .        .       
#> AAAS   .              .       . .        .       
#> AACS   .              .       . 1.001958 2.007853
#> AADAC  .              .       . .        .       
```
