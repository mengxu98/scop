# Standard single-cell and spatial workflow

Normalize, find variable features, scale, reduce dimensions, and cluster
a `Seurat` object. `workflow = "spatial"` adds spot QC, spatial variable
features, optional BayesSpace clustering, and optional deconvolution. It
is not a multi-slice integration orchestrator.

Objects that also contain a `ChromatinAssay` are preprocessed
sequentially (default assay, then chromatin).

## Usage

``` r
RunStandardWorkflow(
  srt,
  prefix = "Standard",
  workflow = c("single_cell", "spatial"),
  assay = NULL,
  image = NULL,
  coord.cols = c("x", "y"),
  do_spot_qc = TRUE,
  spot_qc_params = list(),
  do_spatial_variable_features = TRUE,
  spatial_variable_features_params = list(),
  do_spatial_cluster = FALSE,
  spatial_cluster_method = "BayesSpace",
  spatial_q = NULL,
  bayesspace_params = list(),
  reference = NULL,
  reference_label = NULL,
  reference_assay = NULL,
  do_deconvolution = !is.null(reference),
  deconvolution_method = "RCTD",
  deconvolution_params = list(),
  do_normalization = NULL,
  normalization_method = "LogNormalize",
  do_HVF_finding = TRUE,
  HVF_method = "vst",
  nHVF = 2000,
  HVF = NULL,
  do_scaling = TRUE,
  vars_to_regress = NULL,
  regression_model = "linear",
  linear_reduction = "pca",
  linear_reduction_dims = 50,
  linear_reduction_dims_use = NULL,
  linear_reduction_params = list(),
  force_linear_reduction = FALSE,
  nonlinear_reduction = "umap",
  nonlinear_reduction_dims = 2,
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  cores = 1L,
  verbose = TRUE,
  seed = 11,
  ...
)
```

## Arguments

- srt:

  A `Seurat` object.

- prefix:

  Prefix for intermediate object names.

- workflow:

  `"single_cell"` or `"spatial"` (basic single-image Visium-style).

- assay:

  Assay to use. `NULL` uses the default assay.

- image:

  Spatial image name. Required when multiple images are present; a
  single image is selected automatically when `NULL`.

- coord.cols:

  Metadata coordinate columns used when no image is available.

- do_spot_qc, spot_qc_params:

  Run
  [`RunSpotQC()`](https://mengxu98.github.io/scop/reference/RunSpotQC.md)
  and extra arguments.

- do_spatial_variable_features, spatial_variable_features_params:

  Run
  [`RunSpatialVariableFeatures()`](https://mengxu98.github.io/scop/reference/RunSpatialVariableFeatures.md).
  The workflow defaults `set_variable_features = FALSE` so expression
  HVFs are kept.

- do_spatial_cluster, spatial_cluster_method, spatial_q,
  bayesspace_params:

  Spatial clustering. Only `"BayesSpace"` is supported.
  `spatial_q = NULL` uses the number of ordinary spot clusters.

- reference, reference_label, reference_assay:

  Optional single-cell reference for deconvolution.

- do_deconvolution, deconvolution_method, deconvolution_params:

  Deconvolution (`"RCTD"`, `"SPOTlight"`, or `"Cell2location"`). `NULL`
  runs when `reference` or cell2location signatures are provided.

- do_normalization, normalization_method:

  Normalize if the assay has no scaled data (`NULL`). One of
  `"LogNormalize"`, `"SCT"`, `"TFIDF"`, `"scran"`. `"SCT"` switches
  downstream steps to the `"SCT"` assay.

- do_HVF_finding, HVF_method, nHVF, HVF:

  Highly variable features. `HVF_method` is `"vst"`, `"mvp"`, `"disp"`,
  or `"scran"`.

- do_scaling:

  Force scaling via ScaleData.

- vars_to_regress:

  Regressors used in scaling. `NULL` uses none.

- regression_model:

  `"linear"`, `"poisson"`, or `"negativebinomial"`.

- linear_reduction, linear_reduction_dims, linear_reduction_dims_use,
  linear_reduction_params, force_linear_reduction:

  Linear reduction (`"pca"`, `"svd"`, `"ica"`, `"nmf"`, `"mds"`,
  `"glmpca"`). `linear_reduction_dims_use = NULL` uses estimated
  dimensions, else the first 50.

- nonlinear_reduction, nonlinear_reduction_dims,
  nonlinear_reduction_params, force_nonlinear_reduction:

  Nonlinear reduction (`"umap"`, `"umap-naive"`, `"tsne"`, `"dm"`,
  `"phate"`, `"pacmap"`, `"trimap"`, `"largevis"`, `"fr"`).

- neighbor_metric, neighbor_k:

  Neighbor graph (`"euclidean"`, `"cosine"`, `"manhattan"`,
  `"hamming"`).

- cluster_algorithm, cluster_resolution:

  Clustering (`"louvain"`, `"slm"`, `"leiden"`). Larger
  `cluster_resolution` yields fewer clusters.

- cores:

  Number of CPU cores.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed.

- ...:

  Additional arguments passed to reduction methods.

## Value

A `Seurat` object. Spatial workflows store per-stage execution and
result-storage state in `srt@tools[["run_standard_spatial_workflow"]]`,
including the effective `set_variable_features` value, `result_stored`,
and `result_location`. A failed stage signals an error with the same
table in attribute `standard_spatial_stages`.

## Examples

``` r
library(Matrix)
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-09-06 22:35:31] Start standard processing workflow...
#> ℹ [2026-09-06 22:35:32] Checking a list of <Seurat>...
#> ! [2026-09-06 22:35:32] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 22:35:32] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:35:32] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:35:32] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:35:32] Number of available HVF: 2000
#> ℹ [2026-09-06 22:35:32] Finished check
#> ℹ [2026-09-06 22:35:32] Perform `ScaleData()`
#> ℹ [2026-09-06 22:35:32] Perform pca linear dimension reduction
#> ℹ [2026-09-06 22:35:33] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 22:35:33] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 22:35:33] Reorder clusters...
#> ℹ [2026-09-06 22:35:33] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:35:33] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 22:35:41] Standard processing workflow completed
CellDimPlot(
  pancreas_sub,
  group.by = "SubCellType"
)


# Use a combination of different linear
# or nonlinear dimension reduction methods
linear_reductions <- c(
  "pca", "nmf", "mds"
)
pancreas_sub <- RunStandardWorkflow(
  pancreas_sub,
  linear_reduction = linear_reductions,
  nonlinear_reduction = "umap"
)
#> ℹ [2026-09-06 22:35:41] Start standard processing workflow...
#> ℹ [2026-09-06 22:35:41] Checking a list of <Seurat>...
#> ℹ [2026-09-06 22:35:42] Data 1/1 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:35:42] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:35:42] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:35:42] Number of available HVF: 2000
#> ℹ [2026-09-06 22:35:42] Finished check
#> ℹ [2026-09-06 22:35:42] Perform `ScaleData()`
#> ℹ [2026-09-06 22:35:42] Perform pca linear dimension reduction
#> ℹ [2026-09-06 22:35:42] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 22:35:43] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 22:35:43] Reorder clusters...
#> ℹ [2026-09-06 22:35:43] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:35:43] Perform umap nonlinear dimension reduction
#> ℹ [2026-09-06 22:35:51] Perform nmf linear dimension reduction
#> ℹ [2026-09-06 22:35:51] Running NMF...
#> ℹ StandardBE_ 1 
#> ℹ Positive:  Ccnd1, Spp1, Eno1, Rps2, Mif, Ldha, Pebp1, Gapdh, Aldoa, Prdx1 
#> ℹ      Hspe1, Ybx1, Rplp1, Npm1, Mdk, Rpl12, Ptma, Cldn10, Clu, Mgst1 
#> ℹ      Tkt, Cd24a, Slc25a5, Rpl22l1, Krt8, Tpi1, Dbi, Ranbp1, Sox9, Rps12 
#> ℹ Negative:  Tmem108, Poc1a, Epn3, Wipi1, Tmcc3, Fgf12, Plekho1, Tecpr2, Zbtb4, Gm10941 
#> ℹ      Trf, Man1c1, Hmgcs1, Nipal1, Alpl, Larp1b, Kcnip3, Tnr, Gm15915, Cntln 
#> ℹ      Cbfa2t2, Sh2d4a, Bbc3, Megf6, Naaladl2, Fam46d, Hist2h2ac, Slc52a3, Prim2, 1810041L15Rik 
#> ℹ StandardBE_ 2 
#> ℹ Positive:  Spp1, Atp1b1, Sparc, Id2, Vim, Anxa2, Dbi, Mgst1, Gsta3, Mt1 
#> ℹ      Pdzk1ip1, Bicc1, 1700011H14Rik, Acot1, Nudt19, Adamts1, Rpl22l1, Sox9, Clu, Rps2 
#> ℹ      Jun, Ppp1r1b, Rplp1, S100a10, Cldn3, Krt8, Rps12, Lurap1l, Krt18, Cldn6 
#> ℹ Negative:  Aacs, Tmem108, Poc1a, Epn3, B830012L14Rik, Tmcc3, Rfc1, Wsb1, Plekho1, Ppp2r2b 
#> ℹ      Tecpr2, Zbtb4, Haus8, Gm5420, Man1c1, Hmgcs1, Nipal1, Jam3, Tcerg1, Pgap1 
#> ℹ      Alpl, Larp1b, Kcnip3, Tnr, Ptbp3, Gm15915, Cntln, Rbp2, Uchl5, Ahctf1 
#> ℹ StandardBE_ 3 
#> ℹ Positive:  Clu, Mt1, Ckb, Spp1, Mt2, Sparc, Gsta3, Tpi1, Aldoa, Cldn3 
#> ℹ      Cdkn1c, Krt18, 1700011H14Rik, Mgst1, Serpinh1, Prdx1, Tle6, Anxa2, Mif, Rbbp7 
#> ℹ      Muc1, Paics, Ambp, Ldha, Rps2, Myl12a, Litaf, Gapdh, Ptma, Cat 
#> ℹ Negative:  1110002L01Rik, Rpa3, Elovl6, Aacs, Tmem108, Poc1a, Nop56, B830012L14Rik, Tmcc3, Trib1 
#> ℹ      Fgf12, Plekho1, Slc20a1, Ppp2r2b, Zbtb4, Haus8, Gm10941, Trf, Pde4b, Gm5420 
#> ℹ      Man1c1, Nipal1, Jam3, Tcerg1, Pgap1, Snrpa1, Alpl, Larp1b, Kcnip3, Tnr 
#> ℹ StandardBE_ 4 
#> ℹ Positive:  Tpm1, Anxa5, Spp1, Myl12a, Vim, Krt19, Krt18, Ccnd2, Krt8, Sparc 
#> ℹ      Tnfrsf12a, Csrp1, S100a10, Anxa2, Cyr61, Clu, Cldn7, Cdkn1c, Tagln2, Nudt19 
#> ℹ      Myl12b, Jun, Tmsb10, Cldn6, Hsp90aa1, Myl9, Cd24a, Tsc22d1, 1700011H14Rik, Dbi 
#> ℹ Negative:  Elovl6, Aacs, Tmem108, Poc1a, B830012L14Rik, Tmcc3, Rfc1, Plekho1, Slc20a1, Tecpr2 
#> ℹ      Zbtb4, Gm10941, Tenm4, Man1c1, Nipal1, Jam3, Pgap1, Alpl, Kcnip3, Tnr 
#> ℹ      Ptbp3, Gm15915, Cntln, Ocln, Blvrb, Fras1, Rbp2, Cbfa2t2, Ptgfrn, Ugt8a 
#> ℹ StandardBE_ 5 
#> ℹ Positive:  2810417H13Rik, Hmgb2, Rrm2, Dut, Pcna, Tipin, Lig1, Tuba1b, Mcm5, H2afz 
#> ℹ      Tk1, Gmnn, Dek, Tyms, Mcm3, Rfc2, Ran, Orc6, Tubb5, Ranbp1 
#> ℹ      Mcm6, Uhrf1, Srsf2, Rfc3, Mcm7, Gins2, Rpa2, Siva1, Dnajc9, Smc2 
#> ℹ Negative:  1110002L01Rik, Aacs, Wipi1, B830012L14Rik, Tmcc3, Trib1, Fgf12, Plekho1, Ppp2r2b, Lama1 
#> ℹ      Tecpr2, Zbtb4, P4ha3, Tenm4, Trf, Gm5420, Man1c1, Nipal1, Jam3, Alpl 
#> ℹ      Mgst2, Kcnip3, Tnr, Gm15915, Cbfa2t2, Sh2d4a, Bbc3, Surf4, Ano6, Megf6 
#> ✔ [2026-09-06 22:36:10] NMF compute completed
#> ℹ [2026-09-06 22:36:10] Use stored estimated dimensions 1:50 for Standardnmf
#> ℹ [2026-09-06 22:36:10] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 22:36:11] Reorder clusters...
#> ℹ [2026-09-06 22:36:11] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:36:11] Perform umap nonlinear dimension reduction
#> ℹ [2026-09-06 22:36:18] Perform mds linear dimension reduction
#> Error in RunDimsEstimate(srt = srt, reduction = paste0(prefix, linear_reduction),     reduction_method = linear_reduction, use_stored = FALSE,     verbose = FALSE): No valid dimensions can be estimated for Standardmds
plist1 <- lapply(
  linear_reductions, function(lr) {
    CellDimPlot(
      pancreas_sub,
      group.by = "SubCellType",
      reduction = paste0(
        "Standard", lr, "UMAP2D"
      ),
      xlab = "", ylab = "",
      title = paste0(lr, "_umap"),
      legend.position = "none",
      theme_use = "theme_blank"
    )
  }
)
patchwork::wrap_plots(plist1)


nonlinear_reductions <- c(
  "umap", "tsne", "fr"
)
pancreas_sub <- RunStandardWorkflow(
  pancreas_sub,
  linear_reduction = "pca",
  nonlinear_reduction = nonlinear_reductions
)
#> ℹ [2026-09-06 22:36:20] Start standard processing workflow...
#> ℹ [2026-09-06 22:36:20] Checking a list of <Seurat>...
#> ℹ [2026-09-06 22:36:20] Data 1/1 of the `srt_list` has been log-normalized
#> ℹ [2026-09-06 22:36:20] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:36:20] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:36:20] Number of available HVF: 2000
#> ℹ [2026-09-06 22:36:20] Finished check
#> ℹ [2026-09-06 22:36:20] Perform `ScaleData()`
#> ℹ [2026-09-06 22:36:20] Perform pca linear dimension reduction
#> ℹ [2026-09-06 22:36:21] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-09-06 22:36:21] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 22:36:21] Reorder clusters...
#> ℹ [2026-09-06 22:36:21] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:36:21] Perform umap nonlinear dimension reduction
#> ℹ [2026-09-06 22:36:29] Perform tsne nonlinear dimension reduction
#> ℹ [2026-09-06 22:36:29] Perform tsne nonlinear dimension reduction using Standardpca (1:23)
#> ℹ [2026-09-06 22:36:31] Perform fr nonlinear dimension reduction
#> ℹ [2026-09-06 22:36:31] Perform fr nonlinear dimension reduction using Standardpca_SNN
#> ✔ [2026-09-06 22:36:33] Standard processing workflow completed
plist2 <- lapply(
  nonlinear_reductions, function(nr) {
    CellDimPlot(
      pancreas_sub,
      group.by = "SubCellType",
      reduction = paste0(
        "Standardpca", nr, "2D"
      ),
      xlab = "", ylab = "",
      title = paste0("pca_", nr),
      legend.position = "none",
      theme_use = "theme_blank"
    )
  }
)
patchwork::wrap_plots(plist2)


data(visium_human_pancreas_sub)
spatial <- RunStandardWorkflow(
  visium_human_pancreas_sub,
  workflow = "spatial",
  assay = "Spatial",
  do_spatial_cluster = FALSE,
  spatial_cluster_method = "BayesSpace",
  do_deconvolution = FALSE,
  deconvolution_method = "RCTD",
  linear_reduction_dims = 10,
  linear_reduction_dims_use = 1:5,
  nonlinear_reduction_dims = 2,
  spatial_variable_features_params = list(nfeatures = 50)
)
#> ℹ [2026-09-06 22:36:33] Start standard spot-level spatial workflow...
#> ◌ [2026-09-06 22:36:33] Running spot-level quality control
#> ✔ [2026-09-06 22:36:34] Spot QC completed: 1986 evaluated, 1907 Pass, 79 Fail
#> ℹ   Scope assay "Spatial", layer "counts"
#> ℹ   Saved metadata column `SpotQC`
#> ℹ   Plot returned object `SpatialSpotPlot(<returned_object>, group.by = "SpotQC")`
#> ℹ [2026-09-06 22:36:34] Start standard processing workflow...
#> ℹ [2026-09-06 22:36:34] Checking a list of <Seurat>...
#> ! [2026-09-06 22:36:34] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 22:36:34] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:36:34] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:36:34] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:36:34] Number of available HVF: 2000
#> ℹ [2026-09-06 22:36:34] Finished check
#> ℹ [2026-09-06 22:36:34] Perform `ScaleData()`
#> ℹ [2026-09-06 22:36:34] Perform pca linear dimension reduction
#> ℹ [2026-09-06 22:36:35] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 22:36:35] Reorder clusters...
#> ℹ [2026-09-06 22:36:35] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:36:35] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 22:36:44] Standard processing workflow completed
#> ◌ [2026-09-06 22:36:44] Running spatial variable feature detection
#> ✔ [2026-09-06 22:36:45] Spatial variable features completed: 2000 tested, 50 ranked in the top set
#> ℹ   Scope method "moran" ("cpp" backend); assay "Spatial", layer "data"; 1986 spots; coordinates "raw"
#> ℹ   Saved full result in returned object tool bundle `SpatialVariableFeatures`
#> ℹ   Plot returned object `SpatialVariableFeaturePlot(<returned_object>, plot_type = "combined", assay = "Spatial", image = "slice1", coord.cols = c("x", "y"))`
#> ✔ [2026-09-06 22:36:45] Standard spot-level spatial workflow completed
SpatialSpotPlot(spatial, group.by = "SpotQC")

SpatialSpotPlot(
  spatial,
  features = spatial@tools[["SpatialVariableFeatures"]]$summary$top_features[1:2]
)


spatial_bayes <- RunStandardWorkflow(
  visium_human_pancreas_sub,
  workflow = "spatial",
  assay = "Spatial",
  do_spatial_cluster = TRUE,
  spatial_cluster_method = "BayesSpace",
  spatial_q = 3,
  do_deconvolution = FALSE,
  deconvolution_method = "RCTD",
  bayesspace_params = list(
    n.PCs = 5,
    n.HVGs = 200,
    store_sce = FALSE,
    spatial_cluster_params = list(
      nrep = 200,
      burn.in = 50,
      thin = 10,
      save.chain = FALSE
    )
  )
)
#> ℹ [2026-09-06 22:36:45] Start standard spot-level spatial workflow...
#> ◌ [2026-09-06 22:36:45] Running spot-level quality control
#> ✔ [2026-09-06 22:36:45] Spot QC completed: 1986 evaluated, 1907 Pass, 79 Fail
#> ℹ   Scope assay "Spatial", layer "counts"
#> ℹ   Saved metadata column `SpotQC`
#> ℹ   Plot returned object `SpatialSpotPlot(<returned_object>, group.by = "SpotQC")`
#> ℹ [2026-09-06 22:36:45] Start standard processing workflow...
#> ℹ [2026-09-06 22:36:45] Checking a list of <Seurat>...
#> ! [2026-09-06 22:36:45] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-09-06 22:36:45] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:36:46] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-09-06 22:36:46] Use the separate HVF from `srt_list`
#> ℹ [2026-09-06 22:36:46] Number of available HVF: 2000
#> ℹ [2026-09-06 22:36:46] Finished check
#> ℹ [2026-09-06 22:36:46] Perform `ScaleData()`
#> ℹ [2026-09-06 22:36:46] Perform pca linear dimension reduction
#> ℹ [2026-09-06 22:36:46] Use stored estimated dimensions 1:30 for Standardpca
#> ℹ [2026-09-06 22:36:47] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-09-06 22:36:47] Reorder clusters...
#> ℹ [2026-09-06 22:36:48] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-09-06 22:36:48] Perform umap nonlinear dimension reduction
#> ✔ [2026-09-06 22:36:57] Standard processing workflow completed
#> ◌ [2026-09-06 22:36:57] Running spatial variable feature detection
#> ✔ [2026-09-06 22:36:57] Spatial variable features completed: 2000 tested, 2000 ranked in the top set
#> ℹ   Scope method "moran" ("cpp" backend); assay "Spatial", layer "data"; 1986 spots; coordinates "raw"
#> ℹ   Saved full result in returned object tool bundle `SpatialVariableFeatures`
#> ℹ   Plot returned object `SpatialVariableFeaturePlot(<returned_object>, plot_type = "combined", assay = "Spatial", image = "slice1", coord.cols = c("x", "y"))`
#> ℹ [2026-09-06 22:36:57] Convert <Seurat> to <SingleCellExperiment> for BayesSpace
#> ℹ [2026-09-06 22:36:57] Run BayesSpace spatial clustering with `q = 3`
#> Neighbors were identified for 1974 out of 1986 spots.
#> Fitting model...
#> Calculating labels using iterations 51 through 200.
#> ✔ [2026-09-06 22:37:03] BayesSpace completed: 1986 spots, 3 domains, domain size 129-1310
#> ℹ   Scope assay "Spatial", image "slice1", raw coordinates
#> ℹ   Saved metadata column `BayesSpace_cluster` and returned object tool bundle `BayesSpace`
#> ℹ   Plot returned object `SpatialSpotPlot(<returned_object>, group.by = "BayesSpace_cluster", image = "slice1")`
#> ✔ [2026-09-06 22:37:03] Standard spot-level spatial workflow completed
SpatialSpotPlot(spatial_bayes, group.by = "BayesSpace_cluster")
```
