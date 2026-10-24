# Extract SecAct pattern-associated secreted proteins

Run `SecAct.signaling.pattern.gene` on a SpaCET object after
[`RunSecActSignalingPattern()`](https://mengxu98.github.io/scop/reference/RunSecActSignalingPattern.md).

## Usage

``` r
RunSecActPatternGenes(SpaCET_obj, n, verbose = TRUE)
```

## Arguments

- SpaCET_obj:

  A SpaCET object.

- n:

  Pattern index.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A matrix of pattern-associated secreted proteins.
