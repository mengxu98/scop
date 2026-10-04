# Standard single-cell and spatial workflow

Normalize, find variable features, scale, reduce dimensions, and cluster
a `Seurat` object. `workflow = "spatial"` adds spot QC, spatial variable
features, optional spatial clustering, and optional deconvolution. It is
not a multi-slice integration orchestrator.

Objects that also contain a `ChromatinAssay` are preprocessed
sequentially (default assay, then chromatin).

## Usage

``` r
RunStandardWorkflow(
  object,
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
  spatial_cluster_params = list(),
  do_spatial_qc = FALSE,
  spatial_qc_params = list(),
  spanorm_params = list(),
  spatial_data_type = "auto",
  ...,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- prefix:

  Prefix for intermediate object names.

- workflow:

  `"single_cell"` or `"spatial"` (one explicitly selected spatial
  context).

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

  Spatial clustering with `"BayesSpace"`, `"BANKSY"`, or
  `"SmoothClust"`. `spatial_q = NULL` uses ordinary clusters for
  BayesSpace/SmoothClust. BANKSY uses its resolution parameter and
  rejects a supplied `spatial_q`.

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

  Number of CPU cores. `NULL` lets each backend pick its own default:
  the process OpenMP team (`OMP_NUM_THREADS`) for C++ kernels and one
  worker for R/Python parallel backends.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed.

- spatial_cluster_params:

  Named method-specific arguments. For BayesSpace, use either this list
  or the retained `bayesspace_params`, not both.

- do_spatial_qc, spatial_qc_params:

  Run spatially aware
  [`RunSpotSweeper()`](https://mengxu98.github.io/scop/reference/RunSpotSweeper.md)
  before preprocessing. Filtering is disabled; inspect its QC labels
  first.

- spanorm_params:

  Arguments for
  [`RunSpaNorm()`](https://mengxu98.github.io/scop/reference/RunSpaNorm.md)
  when `normalization_method = "SpaNorm"`. Creates a separate normalized
  assay; count-based backends continue to use the original assay.

- spatial_data_type:

  Observation type: `"auto"`, `"spot"`, `"bin"`, or `"cell"`. Auto uses
  import provenance or image evidence, otherwise unknown. Cell
  observations cannot use BayesSpace or spot deconvolution in this
  workflow.

- ...:

  Additional arguments passed to reduction methods.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

A `Seurat` object. Spatial workflows store per-stage execution and
result-storage state in `srt@tools[["run_standard_spatial_workflow"]]`,
including the effective `set_variable_features` value, `result_stored`,
and `result_location`. A failed stage signals an error with the same
table in attribute `standard_spatial_stages`.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#> ℹ [2026-10-04 22:54:44] Start standard processing workflow...
#> ℹ [2026-10-04 22:54:44] Checking a list of <Seurat>...
#> ! [2026-10-04 22:54:44] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-10-04 22:54:44] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-10-04 22:54:45] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-10-04 22:54:45] Use the separate HVF from `srt_list`
#> ℹ [2026-10-04 22:54:45] Number of available HVF: 2000
#> ℹ [2026-10-04 22:54:45] Finished check
#> ℹ [2026-10-04 22:54:45] Perform `ScaleData()`
#> ℹ [2026-10-04 22:54:45] Perform pca linear dimension reduction
#> ℹ [2026-10-04 22:54:45] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-10-04 22:54:45] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-10-04 22:54:46] Reorder clusters...
#> ℹ [2026-10-04 22:54:46] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-04 22:54:46] Perform umap nonlinear dimension reduction
#> ✔ [2026-10-04 22:54:55] Standard processing workflow completed
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
#> ℹ [2026-10-04 22:54:55] Start standard processing workflow...
#> ℹ [2026-10-04 22:54:55] Checking a list of <Seurat>...
#> ℹ [2026-10-04 22:54:55] Data 1/1 of the `srt_list` has been log-normalized
#> ℹ [2026-10-04 22:54:55] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-10-04 22:54:55] Use the separate HVF from `srt_list`
#> ℹ [2026-10-04 22:54:55] Number of available HVF: 2000
#> ℹ [2026-10-04 22:54:55] Finished check
#> ℹ [2026-10-04 22:54:55] Perform `ScaleData()`
#> ℹ [2026-10-04 22:54:56] Perform pca linear dimension reduction
#> ℹ [2026-10-04 22:54:56] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-10-04 22:54:56] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-10-04 22:54:56] Reorder clusters...
#> ℹ [2026-10-04 22:54:56] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-04 22:54:56] Perform umap nonlinear dimension reduction
#> ℹ [2026-10-04 22:55:05] Perform nmf linear dimension reduction
#> ℹ [2026-10-04 22:55:05] Running NMF...
#> ℹ StandardBE_ 1 
#> ℹ Positive:  Ccnd1, Spp1, Mdk, Rps2, Ldha, Pebp1, Cd24a, Dlk1, Krt8, Mgst1 
#> ℹ      Clu, Gapdh, Eno1, Prdx1, Cldn10, Mif, Cldn7, Npm1, Dbi, Vim 
#> ℹ      Sox9, Rpl12, Aldh1b1, Rplp1, Wfdc2, Krt18, Tkt, Aldoa, Hspe1, Ptma 
#> ℹ Negative:  Tmem108, Poc1a, Epn3, Wipi1, Tmcc3, Nhsl1, Fgf12, Tecpr2, Zbtb4, Plekho1 
#> ℹ      Gm10941, Trf, Man1c1, Hmgcs1, Nipal1, Jam3, Pgap1, Alpl, Tnr, Kcnip3 
#> ℹ      Gm15915, Rbp2, Cbfa2t2, Sh2d4a, Bbc3, Megf6, Naaladl2, Fam46d, Hist2h2ac, Tox2 
#> ℹ StandardBE_ 2 
#> ℹ Positive:  Spp1, Gsta3, Sparc, Vim, Atp1b1, Mt1, Dbi, Anxa2, Rps2, Id2 
#> ℹ      Rpl22l1, Rplp1, Mgst1, Clu, Sox9, Cldn6, Mdk, Pdzk1ip1, Bicc1, 1700011H14Rik 
#> ℹ      Rps12, S100a10, Cldn3, Rpl36a, Ppp1r1b, Adamts1, Serpinh1, Mt2, Ifitm2, Rpl39 
#> ℹ Negative:  Rpa3, Aacs, Tmem108, Poc1a, Epn3, Wipi1, B830012L14Rik, Tmcc3, Wsb1, Tecpr2 
#> ℹ      Zbtb4, Plekho1, Ppp2r2b, Haus8, Trf, Gm5420, Man1c1, Hmgcs1, Nipal1, Jam3 
#> ℹ      Tcerg1, Pgap1, Snrpa1, Alpl, Larp1b, Tnr, Kcnip3, Lsm12, Ptbp3, Gm15915 
#> ℹ StandardBE_ 3 
#> ℹ Positive:  Cck, Mdk, Gadd45a, Neurog3, Selm, Sox4, Btbd17, Tmsb4x, Btg2, Cldn6 
#> ℹ      Cotl1, Ptma, Jun, Ppp1r14a, Rps2, Ifitm2, Neurod2, Igfbpl1, Gnas, Krt7 
#> ℹ      Nkx6-1, Aplp1, Ppp3ca, Lrpap1, Rplp1, Hn1, Rps12, Mfng, BC023829, Smarcd2 
#> ℹ Negative:  Elovl6, Tmem108, Poc1a, Epn3, Nop56, Wipi1, B830012L14Rik, Rrp15, Rfc1, Fgf12 
#> ℹ      Lama1, Slc20a1, Tecpr2, Zbtb4, Ppp2r2b, Eif1ax, Fam162a, P4ha3, Gm10941, Tenm4 
#> ℹ      Pde4b, Gm5420, Man1c1, Hmgcs1, Pgap1, Mgst2, Larp1b, Tnr, Kcnip3, Lsm12 
#> ℹ StandardBE_ 4 
#> ℹ Positive:  Spp1, Cyr61, Krt18, Tpm1, Krt8, Myl12a, Vim, Jun, Anxa5, Tnfrsf12a 
#> ℹ      Csrp1, Sparc, Cldn7, Nudt19, Anxa2, Clu, Myl9, Atp1b1, Cldn3, Tagln2 
#> ℹ      S100a10, 1700011H14Rik, Cd24a, Rps2, Dbi, Id2, Lurap1l, Rplp1, Myl12b, Klf6 
#> ℹ Negative:  Rpa3, Elovl6, Aacs, Tmem108, Poc1a, Tmcc3, Rfc1, Lama1, Slc20a1, Tecpr2 
#> ℹ      Plekho1, Ppp2r2b, Gm10941, Tenm4, Pde4b, Man1c1, Nipal1, Jam3, Pgap1, Alpl 
#> ℹ      Mgst2, Tnr, Kcnip3, Ptbp3, Gm15915, Cntln, Ocln, Fras1, Rbp2, Cbfa2t2 
#> ℹ StandardBE_ 5 
#> ℹ Positive:  2810417H13Rik, Rrm2, Hmgb2, Dut, Pcna, Lig1, H2afz, Tipin, Tuba1b, Tk1 
#> ℹ      Mcm5, Dek, Tyms, Gmnn, Ran, Tubb5, Rfc2, Srsf2, Ranbp1, Orc6 
#> ℹ      Mcm3, Uhrf1, Gins2, Dnajc9, Mcm6, Siva1, Rfc3, Mcm7, Rpa2, Ptma 
#> ℹ Negative:  1110002L01Rik, Aacs, Wipi1, B830012L14Rik, Tmcc3, Trib1, Fgf12, Lama1, Plekho1, Ppp2r2b 
#> ℹ      Tenm4, Trf, Gm5420, Man1c1, Jam3, Mgst2, Tnr, Kcnip3, Gm15915, Cbfa2t2 
#> ℹ      Sh2d4a, Bbc3, Fkbp9, Ano6, Megf6, Prkcb, Fam46d, Tox2, Slc52a3, Ankrd2 
#> ✔ [2026-10-04 22:55:12] NMF compute completed
#> ℹ [2026-10-04 22:55:12] Use stored estimated dimensions 1:50 for Standardnmf
#> ℹ [2026-10-04 22:55:13] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-10-04 22:55:13] Reorder clusters...
#> ℹ [2026-10-04 22:55:13] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-04 22:55:13] Perform umap nonlinear dimension reduction
#> ℹ [2026-10-04 22:55:21] Perform mds linear dimension reduction
#> Error in RunDimsEstimate(object = srt, reduction = paste0(prefix, linear_reduction),     reduction_method = linear_reduction, use_stored = FALSE,     verbose = FALSE): No valid dimensions can be estimated for Standardmds
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
#> ℹ [2026-10-04 22:55:23] Start standard processing workflow...
#> ℹ [2026-10-04 22:55:23] Checking a list of <Seurat>...
#> ℹ [2026-10-04 22:55:23] Data 1/1 of the `srt_list` has been log-normalized
#> ℹ [2026-10-04 22:55:23] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-10-04 22:55:23] Use the separate HVF from `srt_list`
#> ℹ [2026-10-04 22:55:23] Number of available HVF: 2000
#> ℹ [2026-10-04 22:55:23] Finished check
#> ℹ [2026-10-04 22:55:23] Perform `ScaleData()`
#> ℹ [2026-10-04 22:55:23] Perform pca linear dimension reduction
#> ℹ [2026-10-04 22:55:24] Use stored estimated dimensions 1:23 for Standardpca
#> ℹ [2026-10-04 22:55:24] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-10-04 22:55:24] Reorder clusters...
#> ℹ [2026-10-04 22:55:24] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-04 22:55:24] Perform umap nonlinear dimension reduction
#> ℹ [2026-10-04 22:55:33] Perform tsne nonlinear dimension reduction
#> ℹ [2026-10-04 22:55:33] Perform tsne nonlinear dimension reduction using Standardpca (1:23)
#> ℹ [2026-10-04 22:55:35] Perform fr nonlinear dimension reduction
#> ℹ [2026-10-04 22:55:35] Perform fr nonlinear dimension reduction using Standardpca_SNN
#> ✔ [2026-10-04 22:55:36] Standard processing workflow completed
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
#> ℹ [2026-10-04 22:55:36] Start standard spot-level spatial workflow...
#> ◌ [2026-10-04 22:55:36] Running spot-level quality control
#> ✔ [2026-10-04 22:55:36] Spot QC completed: 1986 evaluated, 1907 Pass, 79 Fail
#> ℹ   Scope assay "Spatial", layer "counts"
#> ℹ   Saved metadata column `SpotQC`
#> ℹ   Plot returned object `SpatialSpotPlot(<returned_object>, group.by = "SpotQC")`
#> ℹ [2026-10-04 22:55:36] Start standard processing workflow...
#> ℹ [2026-10-04 22:55:36] Checking a list of <Seurat>...
#> ! [2026-10-04 22:55:36] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-10-04 22:55:36] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-10-04 22:55:37] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-10-04 22:55:37] Use the separate HVF from `srt_list`
#> ℹ [2026-10-04 22:55:37] Number of available HVF: 2000
#> ℹ [2026-10-04 22:55:37] Finished check
#> ℹ [2026-10-04 22:55:37] Perform `ScaleData()`
#> ℹ [2026-10-04 22:55:37] Perform pca linear dimension reduction
#> ℹ [2026-10-04 22:55:37] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-10-04 22:55:38] Reorder clusters...
#> ℹ [2026-10-04 22:55:38] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-04 22:55:38] Perform umap nonlinear dimension reduction
#> ✔ [2026-10-04 22:55:47] Standard processing workflow completed
#> ◌ [2026-10-04 22:55:47] Running spatial variable feature detection
#> ✔ [2026-10-04 22:55:48] Spatial variable features completed: 2000 tested, 50 ranked in the top set
#> ℹ   Scope method "moran" ("cpp" backend); assay "Spatial", layer "data"; 1986 spots; coordinates "raw"
#> ℹ   Saved full result in returned object tool bundle `SpatialVariableFeatures`
#> ℹ   Plot returned object `SpatialVariableFeaturePlot(<returned_object>, plot_type = "combined", assay = "Spatial", image = "slice1", coord.cols = c("x", "y"))`
#> ✔ [2026-10-04 22:55:48] Standard spot-level spatial workflow completed
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
#> ℹ [2026-10-04 22:55:48] Start standard spot-level spatial workflow...
#> ◌ [2026-10-04 22:55:48] Running spot-level quality control
#> ✔ [2026-10-04 22:55:48] Spot QC completed: 1986 evaluated, 1907 Pass, 79 Fail
#> ℹ   Scope assay "Spatial", layer "counts"
#> ℹ   Saved metadata column `SpotQC`
#> ℹ   Plot returned object `SpatialSpotPlot(<returned_object>, group.by = "SpotQC")`
#> ℹ [2026-10-04 22:55:48] Start standard processing workflow...
#> ℹ [2026-10-04 22:55:48] Checking a list of <Seurat>...
#> ! [2026-10-04 22:55:48] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-10-04 22:55:48] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-10-04 22:55:49] Perform `FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-10-04 22:55:49] Use the separate HVF from `srt_list`
#> ℹ [2026-10-04 22:55:49] Number of available HVF: 2000
#> ℹ [2026-10-04 22:55:49] Finished check
#> ℹ [2026-10-04 22:55:49] Perform `ScaleData()`
#> ℹ [2026-10-04 22:55:49] Perform pca linear dimension reduction
#> ℹ [2026-10-04 22:55:50] Use stored estimated dimensions 1:30 for Standardpca
#> ℹ [2026-10-04 22:55:50] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-10-04 22:55:50] Reorder clusters...
#> ℹ [2026-10-04 22:55:50] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-10-04 22:55:50] Perform umap nonlinear dimension reduction
#> ✔ [2026-10-04 22:56:00] Standard processing workflow completed
#> ◌ [2026-10-04 22:56:00] Running spatial variable feature detection
#> ✔ [2026-10-04 22:56:00] Spatial variable features completed: 2000 tested, 2000 ranked in the top set
#> ℹ   Scope method "moran" ("cpp" backend); assay "Spatial", layer "data"; 1986 spots; coordinates "raw"
#> ℹ   Saved full result in returned object tool bundle `SpatialVariableFeatures`
#> ℹ   Plot returned object `SpatialVariableFeaturePlot(<returned_object>, plot_type = "combined", assay = "Spatial", image = "slice1", coord.cols = c("x", "y"))`
#> ℹ [2026-10-04 22:56:01] Convert <Seurat> to <SingleCellExperiment> for BayesSpace
#> ℹ [2026-10-04 22:56:01] Run BayesSpace spatial clustering with `q = 3`
#> Neighbors were identified for 1974 out of 1986 spots.
#> Fitting model...
#> Calculating labels using iterations 51 through 200.
#> ✔ [2026-10-04 22:56:06] BayesSpace completed: 1986 spots, 3 domains, domain size 131-1308
#> ℹ   Scope assay "Spatial", image "slice1", raw coordinates
#> ℹ   Saved metadata column `BayesSpace_cluster` and returned object tool bundle `BayesSpace`
#> ℹ   Plot returned object `SpatialSpotPlot(<returned_object>, group.by = "BayesSpace_cluster", image = "slice1")`
#> ✔ [2026-10-04 22:56:06] Standard spot-level spatial workflow completed
SpatialSpotPlot(spatial_bayes, group.by = "BayesSpace_cluster")
```
