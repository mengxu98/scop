# Run LargeVis (Dimensionality Reduction with a LargeVis-like method)

Run LargeVis (Dimensionality Reduction with a LargeVis-like method)

## Usage

``` r
RunLargeVis(object, ...)

# S3 method for class 'Seurat'
RunLargeVis(
  object,
  reduction = "pca",
  dims = NULL,
  features = NULL,
  assay = NULL,
  layer = "data",
  perplexity = 50,
  n_neighbors = perplexity * 3,
  n_components = 2,
  metric = "euclidean",
  n_epochs = -1,
  learning_rate = 1,
  scale = "maxabs",
  init = "lvrandom",
  init_sdev = NULL,
  repulsion_strength = 7,
  negative_sample_rate = 5,
  nn_method = NULL,
  n_trees = 50,
  search_k = 2 * n_neighbors * n_trees,
  n_threads = NULL,
  n_sgd_threads = 0,
  grain_size = 1,
  kernel = "gauss",
  pca = NULL,
  pca_center = TRUE,
  pcg_rand = TRUE,
  fast_sgd = FALSE,
  batch = FALSE,
  opt_args = NULL,
  epoch_callback = NULL,
  pca_method = NULL,
  reduction.name = "largevis",
  reduction.key = "LargeVis_",
  verbose = TRUE,
  seed.use = 11L,
  ...
)

# Default S3 method
RunLargeVis(
  object,
  assay = NULL,
  perplexity = 50,
  n_neighbors = perplexity * 3,
  n_components = 2,
  metric = "euclidean",
  n_epochs = -1,
  learning_rate = 1,
  scale = "maxabs",
  init = "lvrandom",
  init_sdev = NULL,
  repulsion_strength = 7,
  negative_sample_rate = 5,
  nn_method = NULL,
  n_trees = 50,
  search_k = 2 * n_neighbors * n_trees,
  n_threads = NULL,
  n_sgd_threads = 0,
  grain_size = 1,
  kernel = "gauss",
  pca = NULL,
  pca_center = TRUE,
  pcg_rand = TRUE,
  fast_sgd = FALSE,
  batch = FALSE,
  opt_args = NULL,
  epoch_callback = NULL,
  pca_method = NULL,
  reduction.key = "LargeVis_",
  verbose = TRUE,
  seed.use = 11L,
  ...
)
```

## Arguments

- object:

  An object. This can be a Seurat object or a matrix-like object.

- ...:

  Additional arguments to be passed to
  [uwot::lvish](https://jlmelville.github.io/uwot/reference/lvish.html).

- reduction:

  The reduction to be used. Default is `"pca"`.

- dims:

  The dimensions to be used. Default is `NULL`.

- features:

  The features to be used. Default is `NULL`.

- assay:

  The assay to be used. Default is `NULL`.

- layer:

  The layer to be used. Default is `"data"`.

- perplexity:

  Controls the size of the local neighborhood used for manifold
  approximation. This is the analogous to `n_neighbors` in
  [`umap`](https://jlmelville.github.io/uwot/reference/umap.html).
  Change this, rather than `n_neighbors`.

- n_neighbors:

  The number of neighbors to use when calculating the `perplexity`.
  Usually set to three times the value of the `perplexity`. Must be at
  least as large as `perplexity`.

- n_components:

  The number of LargeVis components. Default is `2`.

- metric:

  Type of distance metric to use to find nearest neighbors. For
  `nn_method = "annoy"` this can be one of:

  - `"euclidean"` (the default)

  - `"cosine"`

  - `"manhattan"`

  - `"hamming"`

  - `"correlation"` (a distance based on the Pearson correlation)

  - `"categorical"` (see below)

  For `nn_method = "hnsw"` this can be one of:

  - `"euclidean"`

  - `"cosine"`

  - `"correlation"`

  If [rnndescent](https://cran.r-project.org/package=rnndescent) is
  installed and `nn_method = "nndescent"` is specified then many more
  metrics are avaiable, including:

  - `"braycurtis"`

  - `"canberra"`

  - `"chebyshev"`

  - `"dice"`

  - `"hamming"`

  - `"hellinger"`

  - `"jaccard"`

  - `"jensenshannon"`

  - `"kulsinski"`

  - `"rogerstanimoto"`

  - `"russellrao"`

  - `"sokalmichener"`

  - `"sokalsneath"`

  - `"spearmanr"`

  - `"symmetrickl"`

  - `"tsss"`

  - `"yule"`

  For more details see the package documentation of `rnndescent`. For
  `nn_method = "fnn"`, the distance metric is always "euclidean".

  If `X` is a data frame or matrix, then multiple metrics can be
  specified, by passing a list to this argument, where the name of each
  item in the list is one of the metric names above. The value of each
  list item should be a vector giving the names or integer ids of the
  columns to be included in a calculation, e.g.
  `metric = list(euclidean = 1:4, manhattan = 5:10)`.

  Each metric calculation results in a separate fuzzy simplicial set,
  which are intersected together to produce the final set. Metric names
  can be repeated. Because non-numeric columns are removed from the data
  frame, it is safer to use column names than integer ids.

  Factor columns can also be used by specifying the metric name
  `"categorical"`. Factor columns are treated different from numeric
  columns and although multiple factor columns can be specified in a
  vector, each factor column specified is processed individually. If you
  specify a non-factor column, it will be coerced to a factor.

  For a given data block, you may override the `pca` and `pca_center`
  arguments for that block, by providing a list with one unnamed item
  containing the column names or ids, and then any of the `pca` or
  `pca_center` overrides as named items, e.g.
  `metric = list(euclidean = 1:4, manhattan = list(5:10, pca_center = FALSE))`.
  This exists to allow mixed binary and real-valued data to be included
  and to have PCA applied to both, but with centering applied only to
  the real-valued data (it is typical not to apply centering to binary
  data before PCA is applied).

- n_epochs:

  Number of epochs to use during the optimization of the embedded
  coordinates. The default is calculate the number of epochs dynamically
  based on dataset size, to give the same number of edge samples as the
  LargeVis defaults. This is usually substantially larger than the UMAP
  defaults. If `n_epochs = 0`, then coordinates determined by `"init"`
  will be returned.

- learning_rate:

  Initial learning rate used in optimization of the coordinates.

- scale:

  Scaling to apply to `X` if it is a data frame or matrix:

  - `"none"` or `FALSE` or `NULL` No scaling.

  - `"Z"` or `"scale"` or `TRUE` Scale each column to zero mean and
    variance 1.

  - `"maxabs"` Center each column to mean 0, then divide each element by
    the maximum absolute value over the entire matrix.

  - `"range"` Range scale the entire matrix, so the smallest element is
    0 and the largest is 1.

  - `"colrange"` Scale each column in the range (0,1).

  For lvish, the default is `"maxabs"`, for consistency with LargeVis.

- init:

  Type of initialization for the coordinates. Options are:

  - `"spectral"` Spectral embedding using the normalized Laplacian of
    the fuzzy 1-skeleton, with Gaussian noise added.

  - `"normlaplacian"`. Spectral embedding using the normalized Laplacian
    of the fuzzy 1-skeleton, without noise.

  - `"random"`. Coordinates assigned using a uniform random distribution
    between -10 and 10.

  - `"lvrandom"`. Coordinates assigned using a Gaussian distribution
    with standard deviation 1e-4, as used in LargeVis (Tang et
    al., 2016) and t-SNE.

  - `"laplacian"`. Spectral embedding using the Laplacian Eigenmap
    (Belkin and Niyogi, 2002).

  - `"pca"`. The first two principal components from PCA of `X` if `X`
    is a data frame, and from a 2-dimensional classical MDS if `X` is of
    class `"dist"`.

  - `"spca"`. Like `"pca"`, but each dimension is then scaled so the
    standard deviation is 1e-4, to give a distribution similar to that
    used in t-SNE and LargeVis. This is an alias for
    `init = "pca", init_sdev = 1e-4`.

  - `"agspectral"` An "approximate global" modification of `"spectral"`
    which all edges in the graph to a value of 1, and then sets a random
    number of edges (`negative_sample_rate` edges per vertex) to 0.1, to
    approximate the effect of non-local affinities.

  - A matrix of initial coordinates.

  For spectral initializations, (`"spectral"`, `"normlaplacian"`,
  `"laplacian"`, `"agspectral"`), if more than one connected component
  is identified, no spectral initialization is attempted. Instead a
  PCA-based initialization is attempted. If `verbose = TRUE` the number
  of connected components are logged to the console. The existence of
  multiple connected components implies that a global view of the data
  cannot be attained with this initialization. Increasing the value of
  `n_neighbors` may help.

- init_sdev:

  If non-`NULL`, scales each dimension of the initialized coordinates
  (including any user-supplied matrix) to this standard deviation. By
  default no scaling is carried out, except when `init = "spca"`, in
  which case the value is `0.0001`. Scaling the input may help if the
  unscaled versions result in initial coordinates with large inter-point
  distances or outliers. This usually results in small gradients during
  optimization and very little progress being made to the layout.
  Shrinking the initial embedding by rescaling can help under these
  circumstances. Scaling the result of `init = "pca"` is usually
  recommended and `init = "spca"` as an alias for
  `init = "pca", init_sdev = 1e-4` but for the spectral initializations
  the scaled versions usually aren't necessary unless you are using a
  large value of `n_neighbors` (e.g. `n_neighbors = 150` or higher). For
  compatibility with recent versions of the Python UMAP package, if you
  are using `init = "spectral"`, then you should also set
  `init_sdev = "range"`, which will range scale each of the columns
  containing the initial data between 0-10. This is not set by default
  to maintain backwards compatibility with previous versions of uwot.

- repulsion_strength:

  Weighting applied to negative samples in low dimensional embedding
  optimization. Values higher than one will result in greater weight
  being given to negative samples.

- negative_sample_rate:

  The number of negative edge/1-simplex samples to use per positive
  edge/1-simplex sample in optimizing the low dimensional embedding.

- nn_method:

  Method for finding nearest neighbors. Options are:

  - `"fnn"`. Use exact nearest neighbors via the
    [FNN](https://cran.r-project.org/package=FNN) package.

  - `"annoy"` Use approximate nearest neighbors via the
    [RcppAnnoy](https://cran.r-project.org/package=RcppAnnoy) package.

  - `"hnsw"` Use approximate nearest neighbors with the Hierarchical
    Navigable Small World (HNSW) method (Malkov and Yashunin, 2018) via
    the [RcppHNSW](https://cran.r-project.org/package=RcppHNSW) package.
    `RcppHNSW` is not a dependency of this package: this option is only
    available if you have installed `RcppHNSW` yourself. Also, HNSW only
    supports the following arguments for `metric`: `"euclidean"`,
    `"cosine"` and `"correlation"`.

  - `"nndescent"` Use approximate nearest neighbors with the Nearest
    Neighbor Descent method (Dong et al., 2011) via the
    [rnndescent](https://cran.r-project.org/package=rnndescent) package.
    `rnndescent` is not a dependency of this package: this option is
    only available if you have installed `rnndescent` yourself.

  By default, if `X` has less than 4,096 vertices, the exact nearest
  neighbors are found. Otherwise, approximate nearest neighbors are
  used. You may also pass precalculated nearest neighbor data to this
  argument. It must be a list consisting of two elements:

  - `"idx"`. A `n_vertices x n_neighbors` matrix containing the integer
    indexes of the nearest neighbors in `X`. Each vertex is considered
    to be its own nearest neighbor, i.e. `idx[, 1] == 1:n_vertices`.

  - `"dist"`. A `n_vertices x n_neighbors` matrix containing the
    distances of the nearest neighbors.

  Multiple nearest neighbor data (e.g. from two different precomputed
  metrics) can be passed by passing a list containing the nearest
  neighbor data lists as items. The `n_neighbors` parameter is ignored
  when using precomputed nearest neighbor data.

- n_trees:

  Number of trees to build when constructing the nearest neighbor index.
  The more trees specified, the larger the index, but the better the
  results. With `search_k`, determines the accuracy of the Annoy nearest
  neighbor search. Only used if the `nn_method` is `"annoy"`. Sensible
  values are between `10` to `100`.

- search_k:

  Number of nodes to search during the neighbor retrieval. The larger k,
  the more the accurate results, but the longer the search takes. With
  `n_trees`, determines the accuracy of the Annoy nearest neighbor
  search. Only used if the `nn_method` is `"annoy"`.

- n_threads:

  Number of threads to use (except during stochastic gradient descent).
  Default is half the number of concurrent threads supported by the
  system. For nearest neighbor search, only applies if
  `nn_method = "annoy"`. If `n_threads > 1`, then the Annoy index will
  be temporarily written to disk in the location determined by
  [`tempfile`](https://rdrr.io/r/base/tempfile.html).

- n_sgd_threads:

  Number of threads to use during stochastic gradient descent. If set to
  \> 1, then be aware that if `batch = FALSE`, results will *not* be
  reproducible, even if `set.seed` is called with a fixed seed before
  running. Set to `"auto"` to use the same value as `n_threads`.

- grain_size:

  The minimum amount of work to do on each thread. If this value is set
  high enough, then less than `n_threads` or `n_sgd_threads` will be
  used for processing, which might give a performance improvement if the
  overhead of thread management and context switching was outweighing
  the improvement due to concurrent processing. This should be left at
  default (`1`) and work will be spread evenly over all the threads
  specified.

- kernel:

  Type of kernel function to create input probabilities. Can be one of
  `"gauss"` (the default) or `"knn"`. `"gauss"` uses the usual Gaussian
  weighted similarities. `"knn"` assigns equal probabilities to every
  edge in the nearest neighbor graph, and zero otherwise, using
  `perplexity` nearest neighbors. The `n_neighbors` parameter is ignored
  in this case.

- pca:

  If set to a positive integer value, reduce data to this number of
  columns using PCA. Doesn't applied if the distance `metric` is
  `"hamming"`, or the dimensions of the data is larger than the number
  specified (i.e. number of rows and columns must be larger than the
  value of this parameter). If you have \> 100 columns in a data frame
  or matrix, reducing the number of columns in this way may
  substantially increase the performance of the nearest neighbor search
  at the cost of a potential decrease in accuracy. In many t-SNE
  applications, a value of 50 is recommended, although there's no
  guarantee that this is appropriate for all settings.

- pca_center:

  If `TRUE`, center the columns of `X` before carrying out PCA. For
  binary data, it's recommended to set this to `FALSE`.

- pcg_rand:

  If `TRUE`, use the PCG random number generator (O'Neill, 2014) during
  optimization. Otherwise, use the faster (but probably less
  statistically good) Tausworthe "taus88" generator. The default is
  `TRUE`. This parameter has been superseded by `rng_type` – if both are
  set, `rng_type` takes precedence.

- fast_sgd:

  If `TRUE`, then the following combination of parameters is set:
  `pcg_rand = TRUE` and `n_sgd_threads = "auto"`. The default is
  `FALSE`. Setting this to `TRUE` will speed up the stochastic
  optimization phase, but give a potentially less accurate embedding,
  and which will not be exactly reproducible even with a fixed seed. For
  visualization, `fast_sgd = TRUE` will give perfectly good results. For
  more generic dimensionality reduction, it's safer to leave
  `fast_sgd = FALSE`. If `fast_sgd = TRUE`, then user-supplied values of
  `pcg_rand` and `n_sgd_threads`, are ignored.

- batch:

  If `TRUE`, then embedding coordinates are updated at the end of each
  epoch rather than during the epoch. In batch mode, results are
  reproducible with a fixed random seed even with `n_sgd_threads > 1`,
  at the cost of a slightly higher memory use. You may also have to
  modify `learning_rate` and increase `n_epochs`, so whether this
  provides a speed increase over the single-threaded optimization is
  likely to be dataset and hardware-dependent.

- opt_args:

  A list of optimizer parameters, used when `batch = TRUE`. The default
  optimization method used is Adam (Kingma and Ba, 2014).

  - `method` The optimization method to use. Either `"adam"` or `"sgd"`
    (stochastic gradient descent). Default: `"adam"`.

  - `beta1` (Adam only). The weighting parameter for the exponential
    moving average of the first moment estimator. Effectively the
    momentum parameter. Should be a floating point value between 0
    and 1. Higher values can smooth oscillatory updates in
    poorly-conditioned situations and may allow for a larger
    `learning_rate` to be specified, but too high can cause divergence.
    Default: `0.5`.

  - `beta2` (Adam only). The weighting parameter for the exponential
    moving average of the uncentered second moment estimator. Should be
    a floating point value between 0 and 1. Controls the degree of
    adaptivity in the step-size. Higher values put more weight on
    previous time steps. Default: `0.9`.

  - `eps` (Adam only). Intended to be a small value to prevent division
    by zero, but in practice can also affect convergence due to its
    interaction with `beta2`. Higher values reduce the effect of the
    step-size adaptivity and bring the behavior closer to stochastic
    gradient descent with momentum. Typical values are between 1e-8 and
    1e-3. Default: `1e-7`.

  - `alpha` The initial learning rate. Default: the value of the
    `learning_rate` parameter.

- epoch_callback:

  A function which will be invoked at the end of every epoch. Its
  signature should be: `(epoch, n_epochs, coords)`, where:

  - `epoch` The current epoch number (between `1` and `n_epochs`).

  - `n_epochs` Number of epochs to use during the optimization of the
    embedded coordinates.

  - `coords` The embedded coordinates as of the end of the current
    epoch, as a matrix with dimensions (N, `n_components`).

- pca_method:

  Method to carry out any PCA dimensionality reduction when the pca
  parameter is specified. Allowed values are: `"irlba"`, `"rsvd"`,
  `"bigstatsr"`, `"svd"`, `"auto"`. Uses `"irlba"`, unless more than 50
  case `"svd"` is used.

- reduction.name:

  The name of the reduction to be stored in the Seurat object. Default
  is `"largevis"`.

- reduction.key:

  The prefix for the column names of the LargeVis embeddings. Default is
  `"LargeVis_"`.

- verbose:

  If `TRUE`, log details to the console.

- seed.use:

  The random seed to be used. Default is `11`.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2025-11-13 12:25:48] Start standard scop workflow...
#> ℹ [2025-11-13 12:25:49] Checking a list of <Seurat> object...
#> ! [2025-11-13 12:25:49] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2025-11-13 12:25:49] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-13 12:25:51] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-13 12:25:52] Use the separate HVF from srt_list
#> ℹ [2025-11-13 12:25:52] Number of available HVF: 2000
#> ℹ [2025-11-13 12:25:52] Finished check
#> ℹ [2025-11-13 12:25:52] Perform `Seurat::ScaleData()`
#> ℹ [2025-11-13 12:25:53] Perform pca linear dimension reduction
#> StandardPC_ 1 
#> Positive:  Aplp1, Cpe, Gnas, Fam183b, Map1b, Hmgn3, Pcsk1n, Chga, Tuba1a, Bex2 
#>     Syt13, Isl1, 1700086L19Rik, Pax6, Chgb, Scgn, Rbp4, Scg3, Gch1, Camk2n1 
#>     Cryba2, Pcsk2, Pyy, Tspan7, Mafb, Hist3h2ba, Dbpht2, Abcc8, Rap1b, Slc38a5 
#> Negative:  Spp1, Anxa2, Sparc, Dbi, 1700011H14Rik, Wfdc2, Gsta3, Adamts1, Clu, Mgst1 
#>     Bicc1, Ldha, Vim, Cldn3, Cyr61, Rps2, Mt1, Ptn, Phgdh, Nudt19 
#>     Smtnl2, Smco4, Habp2, Mt2, Col18a1, Rpl12, Galk1, Cldn10, Acot1, Ccnd1 
#> StandardPC_ 2 
#> Positive:  Rbp4, Tagln2, Tuba1b, Fkbp2, Pyy, Pcsk2, Iapp, Tmem27, Meis2, Tubb4b 
#>     Pcsk1n, Dbpht2, Rap1b, Dynll1, Tubb2a, Sdf2l1, Scgn, 1700086L19Rik, Scg2, Abcc8 
#>     Atp1b1, Hspa5, Fam183b, Papss2, Slc38a5, Scg3, Mageh1, Tspan7, Ppp1r1a, Ociad2 
#> Negative:  Neurog3, Btbd17, Gadd45a, Ppp1r14a, Neurod2, Sox4, Smarcd2, Mdk, Pax4, Btg2 
#>     Sult2b1, Hes6, Grasp, Igfbpl1, Gpx2, Cbfa2t3, Foxa3, Shf, Mfng, Tmsb4x 
#>     Amotl2, Gdpd1, Cdc14b, Epb42, Rcor2, Cotl1, Upk3bl, Rbfox3, Cldn6, Cer1 
#> StandardPC_ 3 
#> Positive:  Nusap1, Top2a, Birc5, Aurkb, Cdca8, Pbk, Mki67, Tpx2, Plk1, Ccnb1 
#>     2810417H13Rik, Incenp, Cenpf, Ccna2, Prc1, Racgap1, Cdk1, Aurka, Cdca3, Hmmr 
#>     Spc24, Kif23, Sgol1, Cenpe, Cdc20, Hist1h1b, Cdca2, Mxd3, Kif22, Ska1 
#> Negative:  Anxa5, Pdzk1ip1, Acot1, Tpm1, Anxa2, Dcdc2a, Capg, Sparc, Ttr, Pamr1 
#>     Clu, Cxcl12, Ndrg2, Hnf1aos1, Gas6, Gsta3, Krt18, Ces1d, Atp1b1, Muc1 
#>     Hhex, Acadm, Spp1, Enpp2, Bcl2l14, Sat1, Smtnl2, 1700011H14Rik, Tgm2, Fam159a 
#> StandardPC_ 4 
#> Positive:  Glud1, Tm4sf4, Akr1c19, Cldn4, Runx1t1, Fev, Pou3f4, Gm43861, Pgrmc1, Arx 
#>     Cd200, Lrpprc, Hmgn3, Ppp1r14c, Pam, Etv1, Tsc22d1, Slc25a5, Akap17b, Pgf 
#>     Fam43a, Emb, Jun, Krt8, Dnajc12, Mid1ip1, Ids, Rgs17, Uchl1, Alcam 
#> Negative:  Ins2, Ins1, Ppp1r1a, Nnat, Calr, Sytl4, Sdf2l1, Iapp, Pdia6, Mapt 
#>     G6pc2, C2cd4b, Npy, Gng12, P2ry1, Ero1lb, Adra2a, Papss2, Arhgap36, Fam151a 
#>     Dlk1, Creld2, Gip, Tmem215, Gm27033, Cntfr, Prss53, C2cd4a, Lyve1, Ociad2 
#> StandardPC_ 5 
#> Positive:  Pdx1, Nkx6-1, Npepl1, Cldn4, Cryba2, Fev, Jun, Chgb, Gng12, Adra2a 
#>     Mnx1, Sytl4, Pdk3, Gm27033, Nnat, Chga, Ins2, 1110012L19Rik, Enho, Krt7 
#>     Mlxipl, Tmsb10, Flrt1, Pax4, Tubb3, Prrg2, Gars, Frzb, BC023829, Gm2694 
#> Negative:  Irx2, Irx1, Gcg, Ctxn2, Tmem27, Ctsz, Tmsb15l, Nap1l5, Pou6f2, Gria2 
#>     Ghrl, Peg10, Smarca1, Arx, Lrpap1, Rgs4, Ttr, Gast, Tmsb15b2, Serpina1b 
#>     Slc16a10, Wnk3, Ly6e, Auts2, Sct, Arg1, Dusp10, Sphkap, Dock11, Edn3 
#> ℹ [2025-11-13 12:25:54] Perform `Seurat::FindClusters()` with louvain and `cluster_resolution` = 0.6
#> ℹ [2025-11-13 12:25:54] Reorder clusters...
#> ℹ [2025-11-13 12:25:54] Perform umap nonlinear dimension reduction
#> ℹ [2025-11-13 12:25:54] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-13 12:25:54] UMAP will return its model
#> ℹ [2025-11-13 12:25:59] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-13 12:25:59] UMAP will return its model
#> ✔ [2025-11-13 12:26:04] Run scop standard workflow done
pancreas_sub <- RunLargeVis(
  object = pancreas_sub,
  features = SeuratObject::VariableFeatures(pancreas_sub)
)
#> 12:26:04 Read 1000 rows and found 2000 numeric columns
#> 12:26:04 Normalizing by max-abs
#> 12:26:04 Using FNN for neighbor search, n_neighbors = 150
#> 12:26:07 Commencing calibration for perplexity = 50 using 2 threads
#> 12:26:10 Initializing from random Gaussian with sd = 1e-4
#> 12:26:10 Commencing optimization for 254033 epochs, with 194344 positive edges
#> 12:26:10 Using rng type: pcg
#> 12:28:36 Optimization finished
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "largevis"
)
```
