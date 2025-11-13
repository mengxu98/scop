# Add features data

Add features data to the `Assay`, `Assay5` or `Seurat` object.

## Usage

``` r
AddFeaturesData(object, ...)

# S3 method for class 'Seurat'
AddFeaturesData(object, features, assay = NULL, ...)

# S3 method for class 'Assay'
AddFeaturesData(object, features, ...)

# S3 method for class 'Assay5'
AddFeaturesData(object, features, ...)
```

## Arguments

- object:

  A `Assay`, `Assay5` or `Seurat` object.

- ...:

  Additional arguments passed to the method.

- features:

  Features data to add.

- assay:

  Assay name to use. Default is `NULL`.

## Value

A `Assay`, `Assay5` or `Seurat` object.

## Examples

``` r
data(pancreas_sub)
features <- GetFeaturesData(pancreas_sub)
pancreas_sub <- AddFeaturesData(pancreas_sub, features)
```
