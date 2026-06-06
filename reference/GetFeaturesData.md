# Get features data

Get the data from the `Assay`, `Assay5` or `Seurat` object.

## Usage

``` r
GetFeaturesData(object, ...)

# S3 method for class 'Seurat'
GetFeaturesData(object, assay = NULL, ...)

# S3 method for class 'Assay'
GetFeaturesData(object, ...)

# S3 method for class 'Assay5'
GetFeaturesData(object, ...)
```

## Arguments

- object:

  A `Assay`, `Assay5` or `Seurat` object.

- ...:

  Additional arguments passed to the method.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

## Value

A data frame containing the features data.

## Examples

``` r
data(pancreas_sub)
features <- GetFeaturesData(pancreas_sub)
head(features)
```
