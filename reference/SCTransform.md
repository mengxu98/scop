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
residual features, memory-conserving or deferred-residual mode,
unsupported regression designs, custom VST arguments, and other flavors
transparently delegate to \`Seurat::SCTransform\`. Legacy Seurat v5
assays with \`do.correct.umi = FALSE\` also delegate to preserve their
output semantics.
