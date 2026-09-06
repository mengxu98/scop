# Seurat clustering with transparent backend delegation

Seurat already performs Louvain and SLM clustering in compiled code and
delegates Leiden clustering to its selected compiled backend. scop
therefore preserves the complete Seurat contract for this step.

## Usage

``` r
FindClusters(object, ...)
```

## Arguments

- object:

  A Seurat object or graph accepted by Seurat.

- ...:

  Arguments passed unchanged to \[Seurat::FindClusters()\].

## Value

The value returned by \[Seurat::FindClusters()\].
