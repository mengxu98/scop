# scop logo

The scop logo, using ASCII or Unicode characters Use
[cli::ansi_strip](https://cli.r-lib.org/reference/ansi_strip.html) to
get rid of the colors.

## Usage

``` r
scop_logo(unicode = cli::is_utf8_output())

# S3 method for class 'scop_logo'
print(x, ...)
```

## Arguments

- unicode:

  Whether to use Unicode symbols on UTF-8 platforms. Default is
  [cli::is_utf8_output](https://cli.r-lib.org/reference/is_utf8_output.html).

- x:

  Input logo object.

- ...:

  Additional arguments passed to \[cat()\].

## References

<https://github.com/tidyverse/tidyverse/blob/main/R/logo.R>

## Examples

``` r
scop_logo()
#>           ⬢          .        ⬡             ⬢     .
#>                      _____ _________  ____
#>                     / ___// ___/ __ ./ __ .
#>                    (__  )/ /__/ /_/ / /_/ /
#>                   /____/ .___/.____/ .___/
#>                                   /_/
#>       ⬢               .      ⬡        .          ⬢
```
