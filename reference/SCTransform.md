# Apply SCTransform normalization

Apply SCTransform normalization

## Usage

``` r
SCTransform(object, ...)
```

## Arguments

- object:

  Object containing count data.

- ...:

  Passed to methods.

## Value

An object with SCTransform results.

## Details

The validated sparse \`vst.flavor = "v2"\` path uses native
corrected-count and residual kernels. Reference models, specified
residual features, memory-conserving mode, unsupported regression
designs, custom VST arguments, and other flavors transparently delegate
to \`Seurat::SCTransform\`.
