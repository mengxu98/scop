# Run SecAct spatial signaling pattern analysis

Run `SecAct.signaling.pattern` on a SpaCET object with existing SecAct
activity results.

## Usage

``` r
RunSecActSignalingPattern(
  SpaCET_obj,
  scale.factor = 1e+05,
  radius = 200,
  k,
  verbose = TRUE
)
```

## Arguments

- SpaCET_obj:

  A SpaCET object.

- scale.factor:

  Spot-level scale factor passed to `SecAct.activity.inference.ST`.

- radius:

  Radius cutoff.

- k:

  Number of NMF patterns, or candidate pattern numbers.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A SpaCET object with SecAct signaling pattern results.
