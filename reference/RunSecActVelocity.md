# Run SecAct spatial signaling velocity

Run SecAct signaling velocity for spot-level ST or
single-cell-resolution ST SpaCET objects.

## Usage

``` r
RunSecActVelocity(
  SpaCET_obj,
  mode = c("spotST", "scST"),
  scale.factor = 1e+05,
  gene = NULL,
  signalMode = "receiving",
  radius = NULL,
  contourMap = FALSE,
  contourBins = 11,
  animated = FALSE,
  sender = NULL,
  secretedProtein = NULL,
  receiver = NULL,
  cellType_meta = NULL,
  CustomizedAreaCoordinates = NULL,
  show.coordinates = TRUE,
  colors = NULL,
  pointSize = 1,
  pointAlpha = 1,
  legend.position = "right",
  legend.size = 1,
  arrow.color = "#ff0099",
  arrow.width = 1,
  arrow.size = 0.3,
  verbose = TRUE
)
```

## Arguments

- SpaCET_obj:

  A SpaCET object.

- mode:

  `spotST` calls `SecAct.signaling.velocity.spotST`; `scST` calls
  `SecAct.signaling.velocity.scST`.

- scale.factor:

  Spot-level scale factor passed to `SecAct.activity.inference.ST`.

- gene:

  Secreted protein gene used by spot-level ST velocity.

- signalMode:

  `receiving` or `sending` for spot-level ST.

- radius:

  Spatial radius. Defaults to `200` for `spotST` and `20` for `scST`
  when `NULL`.

- contourMap, contourBins, animated:

  Plot options for spot-level ST.

- sender, secretedProtein, receiver, cellType_meta:

  Parameters for scST velocity.

- CustomizedAreaCoordinates, show.coordinates, colors, pointSize,
  pointAlpha, legend.position, legend.size, arrow.color, arrow.width,
  arrow.size:

  Plot options for scST velocity.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A ggplot object.
