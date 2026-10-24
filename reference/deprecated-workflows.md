# Deprecated workflow entry points

\`standard_scop()\` and \`integration_scop()\` were renamed to
\[RunStandardWorkflow()\] and \[RunIntegration()\]. The compatibility
entry points remain available with a warning in releases before 1.0.0
and will be removed in version 1.0.0.

## Usage

``` r
standard_scop(...)

integration_scop(...)
```

## Arguments

- ...:

  Arguments forwarded unchanged to the replacement function.

## Value

The value returned by the replacement function.
