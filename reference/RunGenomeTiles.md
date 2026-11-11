# Count fragments in fixed-width genome tiles

Construct explicit non-overlapping genome tiles and delegate validated
sparse fragment counting to
[`RecountATAC()`](https://mengxu98.github.io/scop/reference/RecountATAC.md).
No peaks are needed for defining the tile universe; the source
ChromatinAssay supplies verified fragment mappings and reference
provenance.

## Usage

``` r
RunGenomeTiles(
  object,
  assay = NULL,
  seqlengths = NULL,
  chromosomes = NULL,
  bin.size = 5000L,
  cells = NULL,
  new.assay = NULL,
  process_n = 2000L,
  max.tiles = 2000000L,
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object with a validated ChromatinAssay and fragments.

- assay:

  Source chromatin assay; NULL uses the active assay.

- seqlengths:

  Optional named positive integer chromosome lengths. NULL uses declared
  assay lengths, which must all be known for selected contigs. Supplied
  lengths must agree with known assay lengths. Names are not used to
  guess assembly or chromosome style.

- chromosomes:

  Optional unique contig names in desired output order.

- bin.size:

  Positive integer tile width. The last tile per chromosome is truncated
  at that chromosome's end, never dropped or extended.

- cells:

  Selected source-assay cells; see
  [`RecountATAC()`](https://mengxu98.github.io/scop/reference/RecountATAC.md).

- new.assay:

  Optional unused assay name; NULL returns counts/ranges/receipt.

- process_n:

  FeatureMatrix chunk size; see
  [`RecountATAC()`](https://mengxu98.github.io/scop/reference/RecountATAC.md).

- max.tiles:

  Positive integer upper bound on tile count before allocation.

- verbose:

  Whether to print progress messages.

## Value

The
[`RecountATAC()`](https://mengxu98.github.io/scop/reference/RecountATAC.md)
return value with tile parameters in its receipt. Counts are overlapping
fragment records, not insertion endpoints or unique fragments across
bins. A fragment crossing a boundary can count in two bins.

## See also

[`RecountATAC()`](https://mengxu98.github.io/scop/reference/RecountATAC.md)
