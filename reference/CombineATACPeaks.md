# Combine explicit peak sets into one common interval universe

Apply native
[`GenomicRanges::reduce()`](https://rdrr.io/pkg/IRanges/man/inter-range-methods.html)
within declared groups, optionally require support from several input
peak sets within each group, then take the union of retained group
intervals. No count matrices are combined.

## Usage

``` r
CombineATACPeaks(
  object,
  groups = NULL,
  min.support = 1L,
  min.width = 1L,
  max.width = Inf,
  genome = NULL
)
```

## Arguments

- object:

  A uniquely named list of `GRanges`. Each name identifies one
  peak-calling unit, not necessarily an independent biological
  replicate.

- groups:

  Optional named character vector assigning every peak-set name to one
  support group, for example a condition. With `NULL`, all peak sets
  form one group. Support filtering is done within groups before their
  union.

- min.support:

  Minimum number of distinct input peak sets supporting a reduced
  interval within a group. The default, one, does not require shared
  peaks across samples or conditions.

- min.width, max.width:

  Inclusive width limits applied to the final union.

- genome:

  Optional verified reference label. Known labels must agree; input sets
  without a reference label require this explicit declaration. Absent
  labels are not inferred from chromosome names. Known chromosome
  lengths must be compatible; the reduced union is checked again after
  it inherits known lengths from other inputs, before width filtering.

## Value

A list with `ranges` in native sorted order, named by canonical peak
IDs, and a `receipt` recording input sets, support groups, filtering and
backend versions. Recount fragments over these ranges before comparing
counts from previously different peak universes.

## Details

Reduction uses `ignore.strand = TRUE` and `min.gapwidth = 1L`, so
adjacent intervals are joined. Support counts distinct input sets
contributing to a connected reduced interval; it does not mean every
base has that support. Group-specific peaks can survive without support
in other conditions. Final metadata records contributing sets and
retained support groups, not caller p-values or a newly calculated
consensus significance score.

## See also

[`RecountATAC()`](https://mengxu98.github.io/scop/reference/RecountATAC.md)

## Examples

``` r
first <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 201), width = 100))
second <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(51, 401), width = 100))
common <- CombineATACPeaks(object = list(first = first, second = second), genome = "toy")
common$ranges
#> GRanges object with 3 ranges and 2 metadata columns:
#>                seqnames    ranges strand |      input_sets  support_groups
#>                   <Rle> <IRanges>  <Rle> | <CharacterList> <CharacterList>
#>     chr1-1-150     chr1     1-150      * |    first,second             all
#>   chr1-201-300     chr1   201-300      * |           first             all
#>   chr1-401-500     chr1   401-500      * |          second             all
#>   -------
#>   seqinfo: 1 sequence from toy genome; no seqlengths
common$receipt$input_sets
#>     name group n_peaks
#> 1  first   all       2
#> 2 second   all       2
```
