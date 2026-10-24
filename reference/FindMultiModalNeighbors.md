# Seurat weighted-nearest-neighbor construction

This wrapper preserves all Seurat multimodal-neighbor branches. The
Seurat implementation already uses optimized nearest-neighbor and SNN
kernels, so scop delegates rather than changing graph semantics.

## Usage

``` r
FindMultiModalNeighbors(object, ...)
```

## Arguments

- object:

  A Seurat object.

- ...:

  Arguments passed unchanged to \[Seurat::FindMultiModalNeighbors()\].

## Value

A Seurat object returned by \[Seurat::FindMultiModalNeighbors()\].
