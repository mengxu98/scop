# Deprecated doublet-calling entry points

\`db_scDblFinder()\`, \`db_scds()\`, \`db_Scrublet()\`, and
\`db_DoubletDetection()\` were renamed to \[RunscDblFinder()\],
\[Runscds()\], \[RunScrublet()\], and \[RunDoubletDetection()\]. The
compatibility entry points remain available with a warning and will be
removed in version 1.0.0.

## Usage

``` r
db_scDblFinder(...)

db_scds(...)

db_Scrublet(...)

db_DoubletDetection(...)
```

## Arguments

- ...:

  Arguments forwarded unchanged to the replacement function.

## Value

The value returned by the replacement function.
