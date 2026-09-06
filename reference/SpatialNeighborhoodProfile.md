# Profile cell neighborhoods at several distances

Count labelled neighbors of selected cells without removing the
surrounding tissue. A native spatial index directly accumulates counts
at all distances; no cell-cell edge table is stored. This is an observed
composition summary, not a colocalization test, an optimal-radius
selector, or a communication model.

## Usage

``` r
SpatialNeighborhoodProfile(
  object,
  group.by,
  radii,
  cells = NULL,
  sample.by = NULL,
  image = NULL,
  coord.cols = c("col", "row"),
  cumulative = FALSE,
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object. The object is not modified.

- group.by:

  Metadata column containing nonmissing cell or spot labels.

- radii:

  Positive, strictly increasing distances in raw coordinate units. These
  are not necessarily micrometers. For spot-based assays, counts refer
  to labelled spots, not inferred numbers of cells.

- cells:

  Unique target cell IDs, in the desired output order. NULL uses all
  cells in the resolved image/context. Neighbors always come from the
  complete resolved context, not just these targets.

- sample.by:

  Metadata column identifying independent samples. NULL treats the
  resolved context as one sample; supply this for pooled metadata-only
  coordinates to prevent cross-sample neighbors.

- image:

  Spatial image name. Required when multiple images are present.

- coord.cols:

  Metadata coordinate columns when no image is available.

- cumulative:

  FALSE returns shells \[0, r1\], (r1, r2\], ...; TRUE returns
  cumulative neighborhoods \[0, r\] at each radius.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A data.frame with cell_id, sample, group, lower, radius, count, total
and fraction columns. The denominator includes all labelled neighbors in
that shell/radius. No neighbors gives count = 0 and fraction = NA. Self
is excluded by cell ID, but distinct cells at identical coordinates are
kept. Unselected cells have no output rows. The source and parameters
attributes describe coordinates and the calculation. There is no
tissue-edge correction or dedicated plot; the table can be used with
ordinary R plotting functions.

## Examples

``` r
# Constructed coordinates and labels, not a biological example.
counts <- matrix(1L, 2, 4, dimnames = list(c("g1", "g2"), letters[1:4]))
object <- SeuratObject::CreateSeuratObject(counts)
object$col <- c(0, 1, 3, 5)
object$row <- c(0, 0, 0, 0)
object$celltype <- c("T", "M", "M", "T")
profile <- SpatialNeighborhoodProfile(
  object, "celltype", radii = c(1, 3), cells = "a", verbose = FALSE
)
profile
#>   cell_id sample group lower radius count total fraction
#> 1       a    all     M     0      1     1     1        1
#> 2       a    all     T     0      1     0     1        0
#> 3       a    all     M     1      3     1     1        1
#> 4       a    all     T     1      3     0     1        0
```
