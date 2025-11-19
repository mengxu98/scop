# Run KNN prediction

This function performs KNN prediction to annotate cell types based on
reference scRNA-seq or bulk RNA-seq data.

## Usage

``` r
RunKNNPredict(
  srt_query,
  srt_ref = NULL,
  bulk_ref = NULL,
  query_group = NULL,
  ref_group = NULL,
  query_assay = NULL,
  ref_assay = NULL,
  query_reduction = NULL,
  ref_reduction = NULL,
  query_dims = 1:30,
  ref_dims = 1:30,
  query_collapsing = !is.null(query_group),
  ref_collapsing = TRUE,
  return_full_distance_matrix = FALSE,
  features = NULL,
  features_type = c("HVF", "DE"),
  feature_source = "both",
  nfeatures = 2000,
  DEtest_param = list(max.cells.per.ident = 200, test.use = "wilcox"),
  DE_threshold = "p_val_adj < 0.05",
  nn_method = NULL,
  distance_metric = "cosine",
  k = 30,
  filter_lowfreq = 0,
  prefix = "KNNPredict"
)
```

## Arguments

- srt_query:

  An object of class Seurat to be annotated with cell types.

- srt_ref:

  An object of class Seurat storing the reference cells.

- bulk_ref:

  A cell atlas matrix, where cell types are represented by columns and
  genes are represented by rows, for example, scop::ref_scHCL. Either
  \`srt_ref\` or \`bulk_ref\` must be provided.

- query_group:

  A character vector specifying the column name in the \`srt_query\`
  metadata that represents the cell grouping.

- ref_group:

  A character vector specifying the column name in the \`srt_ref\`
  metadata that represents the cell grouping.

- query_assay:

  A character vector specifying the assay to be used for the query data.
  Defaults to the default assay of the \`srt_query\` object.

- ref_assay:

  A character vector specifying the assay to be used for the reference
  data. Defaults to the default assay of the \`srt_ref\` object.

- query_reduction:

  A character vector specifying the dimensionality reduction method used
  for the query data. If NULL, the function will use the default
  reduction method specified in the \`srt_query\` object.

- ref_reduction:

  A character vector specifying the dimensionality reduction method used
  for the reference data. If NULL, the function will use the default
  reduction method specified in the \`srt_ref\` object.

- query_dims:

  A numeric vector specifying the dimensions to be used for the query
  data. Defaults to the first 30 dimensions.

- ref_dims:

  A numeric vector specifying the dimensions to be used for the
  reference data. Defaults to the first 30 dimensions.

- query_collapsing:

  A boolean value indicating whether the query data should be collapsed
  to group-level average expression values. If TRUE, the function will
  calculate the average expression values for each group in the query
  data and the annotation will be performed separately for each group.
  Otherwise it will use the raw expression values for each cell.

- ref_collapsing:

  A boolean value indicating whether the reference data should be
  collapsed to group-level average expression values. If TRUE, the
  function will calculate the average expression values for each group
  in the reference data and the annotation will be performed separately
  for each group. Otherwise it will use the raw expression values for
  each cell.

- return_full_distance_matrix:

  A boolean value indicating whether the full distance matrix should be
  returned. If TRUE, the function will return the distance matrix used
  for the KNN prediction, otherwise it will only return the annotated
  cell types.

- features:

  A character vector specifying the features (genes) to be used for the
  KNN prediction. If NULL, all the features in the query and reference
  data will be used.

- features_type:

  A character vector specifying the type of features to be used for the
  KNN prediction. Must be one of "HVF" (highly variable features) or
  "DE" (differentially expressed features). Defaults to "HVF".

- feature_source:

  A character vector specifying the source of the features to be used
  for the KNN prediction. Must be one of "both", "query", or "ref".
  Defaults to "both".

- nfeatures:

  An integer specifying the maximum number of features to be used for
  the KNN prediction. Defaults to 2000.

- DEtest_param:

  A list of parameters to be passed to the differential expression test
  function if \`features_type\` is set to "DE". Defaults to
  \`list(max.cells.per.ident = 200, test.use = "wilcox")\`.

- DE_threshold:

  Threshold used to filter the DE features. Default is \`"p_val \<
  0.05"\`. If using "roc" test, `DE_threshold` should be needs to be
  reassigned. e.g. "power \> 0.5".

- nn_method:

  A character vector specifying the method to be used for finding
  nearest neighbors. Must be one of "raw", "rann", or "annoy". Defaults
  to "raw".

- distance_metric:

  A character vector specifying the distance metric to be used for
  calculating similarity between cells. Must be one of "cosine",
  "euclidean", "manhattan", or "hamming". Defaults to "cosine".

- k:

  A number of nearest neighbors to be considered for the KNN prediction.
  Defaults to 30.

- filter_lowfreq:

  An integer specifying the threshold for filtering low-frequency cell
  types from the predicted results. Cell types with a frequency lower
  than \`filter_lowfreq\` will be labelled as "unreliable". Defaults to
  0, which means no filtering will be performed.

- prefix:

  A character vector specifying the prefix to be added to the resulting
  annotations. Defaults to "KNNPredict".

## Examples

``` r
# Annotate cells using bulk RNA-seq data
data(pancreas_sub)
data(ref_scMCA)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2025-11-19 14:54:59] Start standard scop workflow...
#> ℹ [2025-11-19 14:55:00] Checking a list of <Seurat> object...
#> ! [2025-11-19 14:55:00] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2025-11-19 14:55:00] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-19 14:55:02] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-19 14:55:03] Use the separate HVF from srt_list
#> ℹ [2025-11-19 14:55:03] Number of available HVF: 2000
#> ℹ [2025-11-19 14:55:03] Finished check
#> ℹ [2025-11-19 14:55:03] Perform `Seurat::ScaleData()`
#> ℹ [2025-11-19 14:55:03] Perform pca linear dimension reduction
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
#> ℹ [2025-11-19 14:55:04] Perform `Seurat::FindClusters()` with louvain and `cluster_resolution` = 0.6
#> ℹ [2025-11-19 14:55:04] Reorder clusters...
#> ℹ [2025-11-19 14:55:05] Perform umap nonlinear dimension reduction
#> ℹ [2025-11-19 14:55:05] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-19 14:55:05] UMAP will return its model
#> ℹ [2025-11-19 14:55:09] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-19 14:55:09] UMAP will return its model
#> ✔ [2025-11-19 14:55:14] Run scop standard workflow done

# Set the number of threads for RcppParallel
# details see: ?RcppParallel::setThreadOptions
# if (requireNamespace("RcppParallel", quietly = TRUE)) {
#   RcppParallel::setThreadOptions()
# }
pancreas_sub <- RunKNNPredict(
  srt_query = pancreas_sub,
  bulk_ref = ref_scMCA
)
#> ℹ [2025-11-19 14:55:14] Use [1] 549 features to calculate distance.
#> ℹ [2025-11-19 14:55:14] Detected query data type: "log_normalized_counts"
#> ℹ [2025-11-19 14:55:14] Detected reference data type: "log_normalized_counts"
#> ℹ [2025-11-19 14:55:14] Calculate similarity...
#> ℹ [2025-11-19 14:55:14] Use raw method to find neighbors
#> ℹ [2025-11-19 14:55:14] Predict cell type...
CellDimPlot(
  pancreas_sub,
  group.by = "KNNPredict_classification",
  label = TRUE
)


# Removal of low credible cell types from the predicted results
pancreas_sub <- RunKNNPredict(
  srt_query = pancreas_sub,
  bulk_ref = ref_scMCA,
  filter_lowfreq = 30
)
#> ℹ [2025-11-19 14:55:15] Use [1] 549 features to calculate distance.
#> ℹ [2025-11-19 14:55:15] Detected query data type: "log_normalized_counts"
#> ℹ [2025-11-19 14:55:15] Detected reference data type: "log_normalized_counts"
#> ℹ [2025-11-19 14:55:15] Calculate similarity...
#> ℹ [2025-11-19 14:55:15] Use raw method to find neighbors
#> ℹ [2025-11-19 14:55:15] Predict cell type...
CellDimPlot(
  pancreas_sub,
  group.by = "KNNPredict_classification",
  label = TRUE
)


# Annotate clusters using bulk RNA-seq data
pancreas_sub <- RunKNNPredict(
  srt_query = pancreas_sub,
  query_group = "SubCellType",
  bulk_ref = ref_scMCA
)
#> ℹ [2025-11-19 14:55:15] Use [1] 549 features to calculate distance.
#> ℹ [2025-11-19 14:55:16] Detected query data type: "log_normalized_counts"
#> ℹ [2025-11-19 14:55:16] Detected reference data type: "log_normalized_counts"
#> ℹ [2025-11-19 14:55:16] Calculate similarity...
#> ℹ [2025-11-19 14:55:16] Use raw method to find neighbors
#> ℹ [2025-11-19 14:55:16] Predict cell type...
CellDimPlot(
  pancreas_sub,
  group.by = "KNNPredict_classification",
  label = TRUE
)


# Annotate using single cell RNA-seq data
data(panc8_sub)
# Simply convert genes from human to mouse and preprocess the data
genenames <- make.unique(
  thisutils::capitalize(
    rownames(panc8_sub),
    force_tolower = TRUE
  )
)
names(genenames) <- rownames(panc8_sub)
panc8_sub <- RenameFeatures(
  panc8_sub,
  newnames = genenames
)
#> ℹ [2025-11-19 14:55:16] Rename features for the assay: RNA
panc8_sub <- CheckDataMerge(
  panc8_sub,
  batch = "tech"
)[["srt_merge"]]
#> ℹ [2025-11-19 14:55:16] Spliting `srt_merge` into `srt_list` by column "tech"...
#> ℹ [2025-11-19 14:55:17] Checking a list of <Seurat> object...
#> ! [2025-11-19 14:55:17] Data 1/5 of the `srt_list` is "unknown"
#> ℹ [2025-11-19 14:55:17] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/5 of the `srt_list`...
#> ℹ [2025-11-19 14:55:19] Perform `Seurat::FindVariableFeatures()` on the data 1/5 of the `srt_list`...
#> ! [2025-11-19 14:55:19] Data 2/5 of the `srt_list` is "unknown"
#> ℹ [2025-11-19 14:55:19] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 2/5 of the `srt_list`...
#> ℹ [2025-11-19 14:55:21] Perform `Seurat::FindVariableFeatures()` on the data 2/5 of the `srt_list`...
#> ! [2025-11-19 14:55:21] Data 3/5 of the `srt_list` is "unknown"
#> ℹ [2025-11-19 14:55:21] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 3/5 of the `srt_list`...
#> ℹ [2025-11-19 14:55:23] Perform `Seurat::FindVariableFeatures()` on the data 3/5 of the `srt_list`...
#> ! [2025-11-19 14:55:23] Data 4/5 of the `srt_list` is "unknown"
#> ℹ [2025-11-19 14:55:23] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 4/5 of the `srt_list`...
#> ℹ [2025-11-19 14:55:25] Perform `Seurat::FindVariableFeatures()` on the data 4/5 of the `srt_list`...
#> ! [2025-11-19 14:55:25] Data 5/5 of the `srt_list` is "unknown"
#> ℹ [2025-11-19 14:55:25] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 5/5 of the `srt_list`...
#> ℹ [2025-11-19 14:55:27] Perform `Seurat::FindVariableFeatures()` on the data 5/5 of the `srt_list`...
#> ℹ [2025-11-19 14:55:28] Use the separate HVF from srt_list
#> ℹ [2025-11-19 14:55:28] Number of available HVF: 2000
#> ℹ [2025-11-19 14:55:28] Finished check
panc8_sub <- SeuratObject::JoinLayers(panc8_sub)
pancreas_sub <- RunKNNPredict(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  ref_group = "celltype"
)
#> ℹ [2025-11-19 14:55:32] Use the HVF to calculate distance metric
#> ℹ [2025-11-19 14:55:32] Use [1] 632 features to calculate distance.
#> ℹ [2025-11-19 14:55:33] Detected query data type: "log_normalized_counts"
#> ℹ [2025-11-19 14:55:33] Detected reference data type: "log_normalized_counts"
#> ℹ [2025-11-19 14:55:34] Calculate similarity...
#> ℹ [2025-11-19 14:55:34] Use raw method to find neighbors
#> ℹ [2025-11-19 14:55:34] Predict cell type...
CellDimPlot(
  pancreas_sub,
  group.by = "KNNPredict_classification",
  label = TRUE
)

FeatureDimPlot(
  pancreas_sub,
  features = "KNNPredict_simil"
)


pancreas_sub <- RunKNNPredict(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  ref_group = "celltype",
  ref_collapsing = FALSE
)
#> ℹ [2025-11-19 14:55:34] Use the HVF to calculate distance metric
#> ℹ [2025-11-19 14:55:34] Use [1] 632 features to calculate distance.
#> ℹ [2025-11-19 14:55:34] Detected query data type: "log_normalized_counts"
#> ℹ [2025-11-19 14:55:35] Detected reference data type: "log_normalized_counts"
#> ℹ [2025-11-19 14:55:35] Calculate similarity...
#> ℹ [2025-11-19 14:55:35] Use raw method to find neighbors
#> ℹ [2025-11-19 14:55:35] Predict cell type...
CellDimPlot(
  pancreas_sub,
  group.by = "KNNPredict_classification",
  label = TRUE
)

FeatureDimPlot(
  pancreas_sub,
  features = "KNNPredict_prob"
)


pancreas_sub <- RunKNNPredict(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  query_group = "SubCellType",
  ref_group = "celltype"
)
#> ℹ [2025-11-19 14:55:36] Use the HVF to calculate distance metric
#> ℹ [2025-11-19 14:55:36] Use [1] 632 features to calculate distance.
#> ℹ [2025-11-19 14:55:36] Detected query data type: "log_normalized_counts"
#> ℹ [2025-11-19 14:55:36] Detected reference data type: "log_normalized_counts"
#> ℹ [2025-11-19 14:55:36] Calculate similarity...
#> ℹ [2025-11-19 14:55:36] Use raw method to find neighbors
#> ℹ [2025-11-19 14:55:36] Predict cell type...
CellDimPlot(
  pancreas_sub,
  group.by = "KNNPredict_classification",
  label = TRUE
)

FeatureDimPlot(
  pancreas_sub,
  features = "KNNPredict_simil"
)


# Annotate with DE gene instead of HVF
pancreas_sub <- RunKNNPredict(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  ref_group = "celltype",
  features_type = "DE",
  feature_source = "ref",
  DEtest_param = list(cores = 2)
)
#> ✔ [2025-11-19 14:55:36] immunogenomics/presto installed successfully
#> ℹ [2025-11-19 14:55:37] Data type is log-normalized
#> ℹ [2025-11-19 14:55:37] Start differential expression test
#> ℹ [2025-11-19 14:55:37] Find all markers(wilcox) among [1] 13 groups...
#> ℹ [2025-11-19 14:55:37] Using 2 cores
#> ⠙ [2025-11-19 14:55:37] Running [7/13] Processing: delta, acinar, alpha, activa…
#> ✔ [2025-11-19 14:55:37] Completed 13 tasks in 3.9s
#> 
#> ℹ [2025-11-19 14:55:37] Building results
#> ✔ [2025-11-19 14:55:41] Differential expression test completed
#> ℹ [2025-11-19 14:55:41] Use the DE features from AllMarkers_wilcox to calculate distance metric.
#> ℹ [2025-11-19 14:55:41] DE features number of the ref data: [1] 1998
#> ℹ [2025-11-19 14:55:41] Use [1] 1998 features to calculate distance.
#> ℹ [2025-11-19 14:55:42] Detected query data type: "log_normalized_counts"
#> ℹ [2025-11-19 14:55:42] Detected reference data type: "log_normalized_counts"
#> ℹ [2025-11-19 14:55:42] Calculate similarity...
#> ℹ [2025-11-19 14:55:42] Use raw method to find neighbors
#> ℹ [2025-11-19 14:55:42] Predict cell type...

CellDimPlot(
  pancreas_sub,
  group.by = "KNNPredict_classification",
  label = TRUE
)


FeatureDimPlot(
  pancreas_sub,
  features = "KNNPredict_simil"
)


pancreas_sub <- RunKNNPredict(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  query_group = "SubCellType",
  ref_group = "celltype",
  features_type = "DE",
  feature_source = "both",
  DEtest_param = list(cores = 2)
)
#> ✔ [2025-11-19 14:55:43] immunogenomics/presto installed successfully
#> ℹ [2025-11-19 14:55:43] Data type is log-normalized
#> ℹ [2025-11-19 14:55:43] Start differential expression test
#> ℹ [2025-11-19 14:55:43] Find all markers(wilcox) among [1] 8 groups...
#> ℹ [2025-11-19 14:55:43] Using 2 cores
#> ⠙ [2025-11-19 14:55:43] Running [4/8] Processing: Ductal, Beta, Pre-endocrine, …
#> ✔ [2025-11-19 14:55:43] Completed 8 tasks in 1.5s
#> 
#> ℹ [2025-11-19 14:55:43] Building results
#> ✔ [2025-11-19 14:55:45] Differential expression test completed
#> ℹ [2025-11-19 14:55:45] Use the DE features from AllMarkers_wilcox to calculate distance metric.
#> ℹ [2025-11-19 14:55:45] DE features number of the query data: [1] 1998
#> ✔ [2025-11-19 14:55:45] immunogenomics/presto installed successfully
#> ℹ [2025-11-19 14:55:46] Data type is log-normalized
#> ℹ [2025-11-19 14:55:46] Start differential expression test
#> ℹ [2025-11-19 14:55:46] Find all markers(wilcox) among [1] 13 groups...
#> ℹ [2025-11-19 14:55:46] Using 2 cores
#> ⠙ [2025-11-19 14:55:46] Running [7/13] Processing: delta, acinar, alpha, activa…
#> ✔ [2025-11-19 14:55:46] Completed 13 tasks in 3.6s
#> 
#> ℹ [2025-11-19 14:55:46] Building results
#> ✔ [2025-11-19 14:55:49] Differential expression test completed
#> ℹ [2025-11-19 14:55:49] Use the DE features from AllMarkers_wilcox to calculate distance metric.
#> ℹ [2025-11-19 14:55:49] DE features number of the ref data: [1] 352
#> ℹ [2025-11-19 14:55:49] Use [1] 102 features to calculate distance.
#> ℹ [2025-11-19 14:55:50] Detected query data type: "log_normalized_counts"
#> ℹ [2025-11-19 14:55:50] Detected reference data type: "log_normalized_counts"
#> ℹ [2025-11-19 14:55:50] Calculate similarity...
#> ℹ [2025-11-19 14:55:50] Use raw method to find neighbors
#> ℹ [2025-11-19 14:55:50] Predict cell type...

CellDimPlot(
  pancreas_sub,
  group.by = "KNNPredict_classification",
  label = TRUE
)


FeatureDimPlot(
  pancreas_sub,
  features = "KNNPredict_simil"
)
```
