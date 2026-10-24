# Access spatial analysis results

Read the cluster assignments, parameters, and summary for a spatial
method from a Seurat object's stored tool result. If detailed results
were not stored, cluster assignments are read from the corresponding
metadata column. BANKSY retains a lightweight index for custom result
columns and exact sample membership. Sample retrieval requires this
index or stored per-sample results; sample identity is never inferred
from cluster labels.

Results describe cells present in `object`, in its current cell order,
rather than the original fit's full cell set. BANKSY domain counts and
PRECAST cell, domain, and sample counts are recomputed from the returned
assignments. Unfitted PRECAST cells are excluded. Parameters and the
native bundles in `object@tools` still describe the original fit; other
methods' stored summaries are omitted when their cell set changes.
PRECAST can be selected by its tool name or by `"PRECAST"` when that
method identifies a unique stored bundle.

New BANKSY fits and detailed SmoothClust, CHOIR, RareQ, and PRECAST
results retain a private metadata identity column, named by the result
index. Preserve this column to retrieve results after
[`subset()`](https://rdrr.io/r/base/subset.html) or `RenameCells()`;
names are mapped by saved identities, never by order or sample prefixes.
Missing, duplicate, or unknown identities cause an informative error.
SmoothClust and CHOIR's explicit `cell` columns also reflect the current
IDs. Legacy saved assignments support exact-ID subsets. Contradictory
metadata causes an error, but renames cannot always be detected without
provenance; rerun the method to establish identity tracking. Legacy
PRECAST and metadata-only results use their current metadata
assignments.

Storage-disabled methods other than BANKSY retain their existing storage
behavior. Custom metadata columns require a stored locator; the accessor
does not guess a column from its values or prefix.

## Usage

``` r
GetSpatialResult(object, method, sample = NULL)
```

## Arguments

- object:

  A `Seurat` object returned by a spatial method.

- method:

  Name of the method or tool entry, such as `"BANKSY"`.

- sample:

  Optional BANKSY sample label, or a sample label for another method
  that stored per-sample results.

## Value

A list with `clusters`, `parameters`, and `summary`. Values that were
not stored by the method are returned as `NULL`; BANKSY and PRECAST
summaries are recomputed from the returned cluster assignments,
including empty results.

## Examples

``` r
data(visium_human_pancreas_sub)
if (FALSE) { # \dontrun{
spatial <- RunBANKSY(
  visium_human_pancreas_sub,
  layer = "counts",
  features = rownames(visium_human_pancreas_sub)[1:200],
  verbose = FALSE
)
GetSpatialResult(spatial, "BANKSY")
} # }
```
