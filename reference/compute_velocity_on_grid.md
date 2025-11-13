# Compute velocity on grid

Compute velocity on grid

## Usage

``` r
compute_velocity_on_grid(
  x_emb,
  v_emb,
  density = 1,
  smooth = 0.5,
  n_neighbors = ceiling(n_obs/50),
  min_mass = 1,
  scale = 1,
  adjust_for_stream = FALSE,
  cutoff_perc = 5
)
```

## Arguments

- x_emb:

  A matrix of dimension n_obs x n_dim specifying the embedding
  coordinates of the cells.

- v_emb:

  A matrix of dimension n_obs x n_dim specifying the velocity vectors of
  the cells.

- density:

  A numeric value specifying the density of the grid points along each
  dimension. Default is `1`.

- smooth:

  A numeric value specifying the smoothing factor for the velocity
  vectors. Default is `0.5`.

- n_neighbors:

  A numeric value specifying the number of nearest neighbors for each
  grid point. Default is `ceiling(n_obs / 50)`.

- min_mass:

  A numeric value specifying the minimum mass required for a grid point
  to be considered. Default is `1`.

- scale:

  A numeric value specifying the scaling factor for the velocity
  vectors. Default is `1`.

- adjust_for_stream:

  Whether to adjust the velocity vectors for streamlines. Default is
  `FALSE`.

- cutoff_perc:

  A numeric value specifying the percentile cutoff for removing
  low-density grid points. Default is `5`.

## References

<https://github.com/theislab/scvelo/blob/master/scvelo/plotting/velocity_embedding_grid.py>
