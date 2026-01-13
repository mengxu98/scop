# Detect outliers using MAD (Median Absolute Deviation)

Detect outliers using MAD (Median Absolute Deviation)

## Usage

``` r
is_outlier(
  x,
  nmads = 2.5,
  constant = 1.4826,
  type = c("both", "lower", "higher")
)
```

## Arguments

- x:

  Numeric vector.

- nmads:

  Number of median absolute deviations (MADs) from the median to define
  the boundaries for outliers. Default is `2.5`.

- constant:

  Constant factor to convert the MAD to a standard deviation. Default is
  `1.4826`, which is consistent with the MAD of a normal distribution.

- type:

  Type of outliers to detect. Available options are `"both"`, `"lower"`,
  or `"higher"`. If `type` is `"both"`, it detects both lower and higher
  outliers. If `type` is `"lower"`, it detects only lower outliers. If
  `type` is `"higher"`, it detects only higher outliers.

## Value

Numeric vector of indices indicating the positions of outliers in `x`.

## Examples

``` r
x <- c(1, 2, 3, 4, 5, 100)
is_outlier(x) # returns 6
#> [1] 6

x <- c(3, 4, 5, NA, 6, 7)
is_outlier(x, nmads = 1.5, type = "lower") # returns 4
#> [1] 4

x <- c(10, 20, NA, 15, 35)
is_outlier(x, nmads = 2, type = "higher") # returns 3, 5
#> [1] 3 5
```
