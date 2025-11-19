# Single-cell reference mapping with Symphony method

Single-cell reference mapping with Symphony method

## Usage

``` r
RunSymphonyMap(
  srt_query,
  srt_ref,
  query_assay = NULL,
  ref_assay = srt_ref[[ref_pca]]@assay.used,
  ref_pca = NULL,
  ref_harmony = NULL,
  ref_umap = NULL,
  ref_group = NULL,
  projection_method = c("model", "knn"),
  nn_method = NULL,
  k = 30,
  distance_metric = "cosine",
  vote_fun = "mean"
)
```

## Arguments

- srt_query:

  An object of class Seurat storing the query cells.

- srt_ref:

  An object of class Seurat storing the reference cells.

- query_assay:

  A character string specifying the assay name for the query cells. If
  not provided, the default assay for the query object will be used.

- ref_assay:

  A character string specifying the assay name for the reference cells.
  If not provided, the default assay for the reference object will be
  used.

- ref_pca:

  The PCA reduction in the reference object to use for calculating the
  distance metric.

- ref_harmony:

  The Harmony reduction in the reference object to use for calculating
  the distance metric.

- ref_umap:

  A character string specifying the name of the UMAP reduction in the
  reference object. If not provided, the first UMAP reduction found in
  the reference object will be used.

- ref_group:

  A character string specifying a metadata column name in the reference
  object to use for grouping.

- projection_method:

  A character string specifying the projection method to use. Options
  are "model" and "knn". If "model" is selected, the function will try
  to use a pre-trained UMAP model in the reference object for
  projection. If "knn" is selected, the function will directly find the
  nearest neighbors using the distance metric.

- nn_method:

  A character string specifying the nearest neighbor search method to
  use. Options are "raw", "annoy", and "rann". If "raw" is selected, the
  function will use the brute-force method to find the nearest
  neighbors. If "annoy" is selected, the function will use the Annoy
  library for approximate nearest neighbor search. If "rann" is
  selected, the function will use the RANN library for approximate
  nearest neighbor search. If not provided, the function will choose the
  search method based on the size of the query and reference datasets.

- k:

  A number of nearest neighbors to find for each cell in the query
  object.

- distance_metric:

  A character string specifying the distance metric to use for
  calculating the pairwise distances between cells. Options include:
  "pearson", "spearman", "cosine", "correlation", "jaccard", "ejaccard",
  "dice", "edice", "hamman", "simple matching", and "faith". Additional
  distance metrics can also be used, such as "euclidean", "manhattan",
  "hamming", etc.

- vote_fun:

  A character string specifying the function to be used for aggregating
  the nearest neighbors in the reference object. Options are "mean",
  "median", "sum", "min", "max", "sd", "var", etc. If not provided, the
  default is "mean".

## Examples

``` r
data(panc8_sub)
panc8_sub <- standard_scop(panc8_sub)
#> ℹ [2025-11-19 15:11:25] Start standard scop workflow...
#> ℹ [2025-11-19 15:11:26] Checking a list of <Seurat> object...
#> ! [2025-11-19 15:11:26] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2025-11-19 15:11:26] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-19 15:11:28] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-19 15:11:29] Use the separate HVF from srt_list
#> ℹ [2025-11-19 15:11:29] Number of available HVF: 2000
#> ℹ [2025-11-19 15:11:29] Finished check
#> ℹ [2025-11-19 15:11:29] Perform `Seurat::ScaleData()`
#> ℹ [2025-11-19 15:11:30] Perform pca linear dimension reduction
#> StandardPC_ 1 
#> Positive:  CHGA, PCSK1N, G6PC2, PCSK1, IAPP, ARFGEF3, CRYBA2, PRUNE2, CDKN1C, SORL1 
#>     EDN3, CADM1, FXYD2, ELMO1, HADH, PAPPA2, GRIA3, RBP4, DLK1, ANXA6 
#>     HMGN2, GNAZ, AMPD2, IGF2, ROBO2, DNAJA4, PDK4, SEPT3, CD99L2, SYT17 
#> Negative:  IFITM3, ZFP36L1, SOX4, ANXA4, KRT7, TPM1, PMEPA1, SERPING1, TM4SF1, CD44 
#>     CDC42EP1, TMSB10, NFIB, SAT1, SDC4, SPTBN1, LCN2, KRT18, PDZK1IP1, MSN 
#>     SMAD3, CLDN10, CFTR, NOTCH2, KRT19, CTSH, SERPINA5, FLRT2, C3, EPS8 
#> StandardPC_ 2 
#> Positive:  SPARC, COL4A1, COL15A1, COL1A2, COL3A1, PXDN, PDGFRB, COL5A1, BGN, COL5A2 
#>     COL1A1, LAMA4, TIMP3, COL6A2, IGFBP4, AEBP1, SFRP2, THBS2, FBN1, COL6A1 
#>     CDH11, VCAN, SERPINE1, WNT5A, FN1, TPM2, FMOD, MMP2, SNAI1, DCN 
#> Negative:  KRT8, SPINK1, PRSS1, ELF3, GATM, MUC1, KRT18, CPA2, CTRB1, SDC4 
#>     PRSS3, CLDN4, LCN2, ANPEP, CPA1, PDZK1IP1, PLA2G1B, CTRC, CPB1, PNLIP 
#>     KLK1, CELA2A, CELA3A, KRT7, GSTA1, CD44, PNLIPRP1, PNLIPRP2, CELA3B, GSTA2 
#> StandardPC_ 3 
#> Positive:  FTO, SORL1, TBC1D24, CASR, PCYOX1, UTRN, ADH5, ENPP5, RNF14, PHKB 
#>     MAP1A, C2CD5, TTC17, RAB22A, PRR14L, AP3B1, MTR, HERC1, EXPH5, SMCHD1 
#>     ROBO1, ABHD10, PRUNE2, SPEN, BTBD3, IBTK, ARFGEF2, TSC1, PARP4, RMND5A 
#> Negative:  HSPB1, CELA3A, CELA3B, CLPS, CTRB1, SYCN, CELA2A, EIF4A1, VIM, PNLIPRP1 
#>     PLA2G1B, KLK1, CPA1, CTRC, DDIT4, PLTP, BGN, DYNLL2, ANGPTL4, COL6A2 
#>     IFITM1, IGFBP4, IGFBP2, TMSB10, PRSS1, CTRL, PDGFRB, CPA2, PRSS3, PXDN 
#> StandardPC_ 4 
#> Positive:  CPA2, PNLIP, PRSS1, CTRC, CPA1, CPB1, PLA2G1B, PNLIPRP2, PRSS3, BCAT1 
#>     CEL, KLK1, CELA2A, CTRB1, PNLIPRP1, SPINK1, GSTA2, MGST1, CELA3A, LDHB 
#>     ALB, CTRL, CELA3B, CLPS, ALDOB, REG3G, FAM129A, GSTA1, SYCN, CBS 
#> Negative:  CFTR, MMP7, KRT19, SERPINA5, TINAGL1, AQP1, SPP1, SERPING1, PMEPA1, KRT23 
#>     ALDH1A3, TSPAN8, PROM1, IGFBP7, VCAM1, LGALS4, ONECUT2, TRPV6, CCL2, ANXA3 
#>     TNFAIP2, CTSH, SDC1, SLC3A1, CLDN10, ANXA9, CCND1, KRT80, VNN1, PDGFD 
#> StandardPC_ 5 
#> Positive:  COL5A1, COL1A2, COL1A1, SFRP2, COL5A2, COL3A1, VCAN, FN1, PDGFRB, THBS2 
#>     FMOD, BGN, ANTXR1, MXRA8, COL6A1, AEBP1, TPM2, CDH11, DCN, ISLR 
#>     TGFB3, COL6A2, LTBP2, DDR2, EDNRA, ANO1, LTBP1, GFPT2, WNT5A, HEYL 
#> Negative:  CD93, PLVAP, PODXL, ACVRL1, ESAM, S1PR1, CXCR4, ECSCR, DYSF, CALCRL 
#>     ADGRF5, STC1, CD34, AFAP1L1, IFI27, SH3BP5, ACKR3, ANGPT2, DLL4, MMRN2 
#>     MCAM, PNP, IL3RA, SPARCL1, TCF4, FAM198B, RAPGEF5, ARHGAP31, P2RY6, F2RL3 
#> ℹ [2025-11-19 15:11:31] Perform `Seurat::FindClusters()` with louvain and `cluster_resolution` = 0.6
#> ℹ [2025-11-19 15:11:31] Reorder clusters...
#> ℹ [2025-11-19 15:11:31] Perform umap nonlinear dimension reduction
#> ℹ [2025-11-19 15:11:31] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-19 15:11:31] UMAP will return its model
#> ℹ [2025-11-19 15:11:36] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-19 15:11:36] UMAP will return its model
#> ✔ [2025-11-19 15:11:41] Run scop standard workflow done
srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
srt_ref <- integration_scop(
  srt_ref,
  batch = "tech",
  integration_method = "Harmony"
)
#> ◌ [2025-11-19 15:11:42] Run Harmony integration...
#> ℹ [2025-11-19 15:11:42] Spliting `srt_merge` into `srt_list` by column "tech"...
#> ℹ [2025-11-19 15:11:43] Checking a list of <Seurat> object...
#> ℹ [2025-11-19 15:11:43] Data 1/4 of the `srt_list` has been log-normalized
#> ℹ [2025-11-19 15:11:43] Perform `Seurat::FindVariableFeatures()` on the data 1/4 of the `srt_list`...
#> ℹ [2025-11-19 15:11:43] Data 2/4 of the `srt_list` has been log-normalized
#> ℹ [2025-11-19 15:11:43] Perform `Seurat::FindVariableFeatures()` on the data 2/4 of the `srt_list`...
#> ℹ [2025-11-19 15:11:44] Data 3/4 of the `srt_list` has been log-normalized
#> ℹ [2025-11-19 15:11:44] Perform `Seurat::FindVariableFeatures()` on the data 3/4 of the `srt_list`...
#> ℹ [2025-11-19 15:11:45] Data 4/4 of the `srt_list` has been log-normalized
#> ℹ [2025-11-19 15:11:45] Perform `Seurat::FindVariableFeatures()` on the data 4/4 of the `srt_list`...
#> ℹ [2025-11-19 15:11:45] Use the separate HVF from srt_list
#> ℹ [2025-11-19 15:11:45] Number of available HVF: 2000
#> ℹ [2025-11-19 15:11:46] Finished check
#> ℹ [2025-11-19 15:11:48] Perform `Seurat::ScaleData()`
#> ℹ [2025-11-19 15:11:48] Perform linear dimension reduction("pca")
#> ℹ [2025-11-19 15:11:49] Perform Harmony integration
#> ℹ [2025-11-19 15:11:49] Harmony integration using "Harmonypca", dims:1-10 as input
#> ✔ [2025-11-19 15:11:50] harmony installed successfully
#> ℹ [2025-11-19 15:11:50] Perform Seurat::FindClusters (louvain)
#> ℹ [2025-11-19 15:11:50] Reorder clusters...
#> ℹ [2025-11-19 15:11:51] Perform nonlinear dimension reduction (umap)
#> ℹ [2025-11-19 15:11:51] Non-linear dimensionality reduction (umap) using (Harmony) dims (1-10) as input
#> ℹ [2025-11-19 15:11:56] Non-linear dimensionality reduction (umap) using (Harmony) dims (1-10) as input
#> ✔ [2025-11-19 15:12:02] Run Harmony integration done
CellDimPlot(srt_ref, group.by = c("celltype", "tech"))


# Projection
srt_query <- RunSymphonyMap(
  srt_query = srt_query,
  srt_ref = srt_ref,
  ref_pca = "Harmonypca",
  ref_harmony = "Harmony",
  ref_umap = "HarmonyUMAP2D"
)
#> ◌ [2025-11-19 15:12:03] Installing: symphony...
#>  
#> → Will install 2 packages.
#> → All 2 packages (0 B) are cached.
#> + RANN       2.6.2 
#> + symphony   0.1.2 [bld][cmp] (GitHub: 91e2ca1)
#> ✔ All system requirements are already installed.
#>   
#> ℹ No downloads are needed, 2 pkgs are cached
#> ✔ Got RANN 2.6.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (43.84 kB)
#> ✔ Got symphony 0.1.2 (source) (4.19 MB)
#> ℹ Installing system requirements
#> ℹ Executing `sudo sh -c apt-get -y update`
#> Get:1 file:/etc/apt/apt-mirrors.txt Mirrorlist [144 B]
#> Hit:2 http://azure.archive.ubuntu.com/ubuntu noble InRelease
#> Hit:3 http://azure.archive.ubuntu.com/ubuntu noble-updates InRelease
#> Hit:6 https://packages.microsoft.com/repos/azure-cli noble InRelease
#> Hit:4 http://azure.archive.ubuntu.com/ubuntu noble-backports InRelease
#> Hit:5 http://azure.archive.ubuntu.com/ubuntu noble-security InRelease
#> Hit:7 https://packages.microsoft.com/ubuntu/24.04/prod noble InRelease
#> Reading package lists...
#> ℹ Executing `sudo sh -c apt-get -y install libicu-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> libicu-dev is already the newest version (74.2-1ubuntu3.1).
#> 0 upgraded, 0 newly installed, 0 to remove and 31 not upgraded.
#> ✔ Installed RANN 2.6.2  (1s)
#> ℹ Packaging symphony 0.1.2
#> ✔ Packaged symphony 0.1.2 (676ms)
#> ℹ Building symphony 0.1.2
#> ✔ Built symphony 0.1.2 (19.6s)
#> ✔ Installed symphony 0.1.2 (github::immunogenomics/symphony@91e2ca1) (1s)
#> ✔ 1 pkg + 44 deps: kept 43, added 2, dld 2 (NA B) [32.5s]
#> ✔ [2025-11-19 15:12:35] immunogenomics/symphony installed successfully
#> ℹ [2025-11-19 15:12:35] Data type is log-normalized
#> ℹ [2025-11-19 15:12:35] Detected srt_query data type: "log_normalized_counts"
#> ℹ [2025-11-19 15:12:36] Data type is log-normalized
#> ℹ [2025-11-19 15:12:36] Detected srt_ref data type: "log_normalized_counts"
#> ℹ [2025-11-19 15:12:36] Build reference
#> ℹ [2025-11-19 15:12:36] Saved embeddings
#> ℹ [2025-11-19 15:12:36] Saved soft cluster assignments
#> ℹ [2025-11-19 15:12:36] Saved variable gene information for 2000 genes
#> ℹ [2025-11-19 15:12:36] Saved PCA loadings
#> ℹ [2025-11-19 15:12:36] Saved metadata
#> ℹ [2025-11-19 15:12:36] Calculate final L2 normalized reference centroids (Y_cos)
#> ℹ [2025-11-19 15:12:36] Calculate reference compression terms (Nr and C)
#> ℹ [2025-11-19 15:12:36] Run mapQuery
#> ℹ [2025-11-19 15:12:36] Scaling and synchronizing query gene expression
#> ℹ [2025-11-19 15:12:36] Found 2000 reference variable genes in query dataset
#> ℹ [2025-11-19 15:12:36] Project query cells using reference gene loadings
#> ℹ [2025-11-19 15:12:36] Clustering query cells to reference centroids
#> ℹ [2025-11-19 15:12:36] Correcting query batch effects
#> ℹ [2025-11-19 15:12:36] Run UMAP projection
#> ℹ [2025-11-19 15:12:36] Use the reduction to calculate distance metric
#> ℹ [2025-11-19 15:12:36] Use raw method to find neighbors
#> ℹ [2025-11-19 15:12:37] Running UMAP projection
#> ℹ [2025-11-19 15:12:37] Run SymphonyMap finished
ProjectionPlot(
  srt_query = srt_query,
  srt_ref = srt_ref,
  query_group = "celltype",
  ref_group = "celltype"
)
#> Scale for x is already present.
#> Adding another scale for x, which will replace the existing scale.
#> Scale for y is already present.
#> Adding another scale for y, which will replace the existing scale.
```
