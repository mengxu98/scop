# Run UMAP (Uniform Manifold Approximation and Projection)

Run UMAP (Uniform Manifold Approximation and Projection)

## Usage

``` r
RunUMAP2(object, ...)

# S3 method for class 'Seurat'
RunUMAP2(
  object,
  reduction = "pca",
  dims = NULL,
  features = NULL,
  neighbor = NULL,
  graph = NULL,
  assay = NULL,
  layer = "data",
  umap.method = "uwot",
  reduction.model = NULL,
  n_threads = NULL,
  return.model = FALSE,
  n.neighbors = 30L,
  n.components = 2L,
  metric = "cosine",
  n.epochs = 200L,
  spread = 1,
  min.dist = 0.3,
  set.op.mix.ratio = 1,
  local.connectivity = 1L,
  negative.sample.rate = 5L,
  a = NULL,
  b = NULL,
  learning.rate = 1,
  repulsion.strength = 1,
  reduction.name = "umap",
  reduction.key = "UMAP_",
  verbose = TRUE,
  seed.use = 11,
  ...
)

# Default S3 method
RunUMAP2(
  object,
  assay = NULL,
  umap.method = "uwot",
  reduction.model = NULL,
  n_threads = NULL,
  return.model = FALSE,
  n.neighbors = 30L,
  n.components = 2L,
  metric = "cosine",
  n.epochs = 200L,
  spread = 1,
  min.dist = 0.3,
  set.op.mix.ratio = 1,
  local.connectivity = 1L,
  negative.sample.rate = 5L,
  a = NULL,
  b = NULL,
  learning.rate = 1,
  repulsion.strength = 1,
  reduction.key = "UMAP_",
  verbose = TRUE,
  seed.use = 11L,
  ...
)
```

## Arguments

- object:

  An object. This can be a Seurat object, a matrix-like object, a
  Neighbor object, or a Graph object.

- ...:

  Additional arguments to be passed to UMAP.

- reduction:

  The reduction to be used. Default is `"pca"`.

- dims:

  The dimensions to be used. Default is `NULL`.

- features:

  The features to be used. Default is `NULL`.

- neighbor:

  The name of the Neighbor object to be used. Default is `NULL`.

- graph:

  The name of the Graph object to be used. Default is `NULL`.

- assay:

  The assay to be used. Default is `NULL`.

- layer:

  The layer to be used. Default is `"data"`.

- umap.method:

  The UMAP method to be used. Options are `"naive"` and `"uwot"`.
  Default is `"uwot"`.

- reduction.model:

  A DimReduc object containing a pre-trained UMAP model. Default is
  `NULL`.

- n_threads:

  Num of threads used.

- return.model:

  Whether to return the UMAP model. Default is `FALSE`.

- n.neighbors:

  A number of nearest neighbors to be used. Default is `30`.

- n.components:

  A number of UMAP components. Default is `2`.

- metric:

  The metric or a function to be used for distance calculations. When
  using a string, available metrics are: `euclidean`, `manhattan`. Other
  available generalized metrics are: cosine, pearson, pearson2. Note the
  triangle inequality may not be satisfied by some generalized metrics,
  hence knn search may not be optimal. When using metric.function as a
  function, the signature must be function(matrix, origin, target) and
  should compute a distance between the origin column and the target
  columns. Default is `"cosine"`.

- n.epochs:

  A number of iterations performed during layout optimization for UMAP.
  Default is `200`.

- spread:

  The spread parameter for UMAP, used during automatic estimation of a/b
  parameters. Default is `1`.

- min.dist:

  The minimum distance between UMAP embeddings, determines how close
  points appear in the final layout. Default is `0.3`.

- set.op.mix.ratio:

  Interpolate between (fuzzy) union and intersection as the set
  operation used to combine local fuzzy simplicial sets to obtain a
  global fuzzy simplicial sets. Both fuzzy set operations use the
  product t-norm. The value of this parameter should be between `0.0`
  and `1.0`; a value of `1.0` will use a pure fuzzy union, while `0.0`
  will use a pure fuzzy intersection.

- local.connectivity:

  The local connectivity, used during construction of fuzzy simplicial
  set. Default is `1`.

- negative.sample.rate:

  The negative sample rate for UMAP optimization. Determines how many
  non-neighbor points are used per point and per iteration during layout
  optimization. Default is `5`.

- a:

  The parameter a for UMAP optimization. Contributes to gradient
  calculations during layout optimization. When left at NA, a suitable
  value will be estimated automatically. Default is `NULL`.

- b:

  The parameter b for UMAP optimization. Details see parameter `a`.

- learning.rate:

  The initial value of "learning rate" of layout optimization. Default
  is `1`.

- repulsion.strength:

  A numeric value determines, together with alpha, the learning rate of
  layout optimization. Default is `1`.

- reduction.name:

  The name of the reduction to be stored in the Seurat object. Default
  is `"umap"`.

- reduction.key:

  The prefix for the column names of the UMAP embeddings. Default is
  `"UMAP_"`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed.use:

  The random seed to be used. Default is `11`.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2025-11-13 12:42:56] Start standard scop workflow...
#> ℹ [2025-11-13 12:42:57] Checking a list of <Seurat> object...
#> ! [2025-11-13 12:42:57] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2025-11-13 12:42:57] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-13 12:42:59] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-13 12:43:00] Use the separate HVF from srt_list
#> ℹ [2025-11-13 12:43:00] Number of available HVF: 2000
#> ℹ [2025-11-13 12:43:00] Finished check
#> ℹ [2025-11-13 12:43:00] Perform `Seurat::ScaleData()`
#> ℹ [2025-11-13 12:43:00] Perform pca linear dimension reduction
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
#> ℹ [2025-11-13 12:43:01] Perform `Seurat::FindClusters()` with louvain and `cluster_resolution` = 0.6
#> ℹ [2025-11-13 12:43:02] Reorder clusters...
#> ℹ [2025-11-13 12:43:02] Perform umap nonlinear dimension reduction
#> ℹ [2025-11-13 12:43:02] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-13 12:43:02] UMAP will return its model
#> ℹ [2025-11-13 12:43:06] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-13 12:43:06] UMAP will return its model
#> ✔ [2025-11-13 12:43:11] Run scop standard workflow done
pancreas_sub <- RunUMAP2(
  object = pancreas_sub,
  features = SeuratObject::VariableFeatures(pancreas_sub)
)
CellDimPlot(
  pancreas_sub,
  group.by = "CellType",
  reduction = "umap"
)
```
