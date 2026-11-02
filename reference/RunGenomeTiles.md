# Count fragments in fixed-width genome tiles

Generate non-overlapping genome tiles and delegate sparse fragment
counting to RecountATAC. No peak universe is required to define tiles;
the source assay provides fragment mappings and reference provenance.

## Usage

``` r
RunGenomeTiles(object, assay = NULL, seqlengths = NULL, chromosomes = NULL,
  bin.size = 5000L, cells = NULL, new.assay = NULL, process_n = 2000L,
  max.tiles = 2000000L, verbose = TRUE)
```

## Arguments

- object:

  Seurat object with validated ChromatinAssay and fragments.

- assay:

  Source assay; NULL selects active assay.

- seqlengths:

  Named positive integer reference lengths, or NULL to use known assay
  lengths. Supplied lengths must agree with known source lengths.

- chromosomes:

  Optional unique contigs in output order.

- bin.size:

  Positive integer tile width; last tile per chromosome is truncated.

- cells:

  Source-assay cells; see RecountATAC.

- new.assay:

  Unused assay name for writeback; NULL returns a list.

- process_n:

  Native FeatureMatrix chunk size.

- max.tiles:

  Positive upper bound on tile count before allocation.

- verbose:

  Whether to print progress messages.

## Value

RecountATAC's return value with tile parameters in its receipt. Counts
are overlapping fragment records, not insertion endpoints or unique
fragments across bins. A fragment crossing a boundary can count in both
bins.

## See also

[`RecountATAC`](https://mengxu98.github.io/scop/reference/RecountATAC.md)
