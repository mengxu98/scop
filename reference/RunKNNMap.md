# Single-cell reference mapping with KNN method

This function performs single-cell reference mapping using the K-nearest
neighbor (KNN) method. It takes two single-cell datasets as input:
srt_query and srt_ref. The function maps cells from the srt_query
dataset to the srt_ref dataset based on their similarity or distance.

## Usage

``` r
RunKNNMap(
  srt_query,
  srt_ref,
  query_assay = NULL,
  ref_assay = NULL,
  ref_umap = NULL,
  ref_group = NULL,
  features = NULL,
  nfeatures = 2000,
  query_reduction = NULL,
  ref_reduction = NULL,
  query_dims = 1:30,
  ref_dims = 1:30,
  projection_method = c("model", "knn"),
  nn_method = NULL,
  k = 30,
  distance_metric = "cosine",
  vote_fun = "mean"
)
```

## Arguments

- srt_query:

  An object of class Seurat to be annotated with cell types.

- srt_ref:

  A Seurat object or count matrix representing the reference object. If
  provided, the similarities will be calculated between cells from the
  query and reference objects. If not provided, the similarities will be
  calculated within the query object.

- query_assay:

  The assay to use for the query object. If not provided, the default
  assay of the query object will be used.

- ref_assay:

  The assay to use for the reference object. If not provided, the
  default assay of the reference object will be used.

- ref_umap:

  A character string specifying the name of the UMAP reduction in the
  reference object. If not provided, the first UMAP reduction found in
  the reference object will be used.

- ref_group:

  The grouping variable in the reference object. This variable will be
  used to group cells in the heatmap columns. If not provided, all cells
  will be treated as one group.

- features:

  A vector of feature names to include in the heatmap. If not provided,
  highly variable features (HVF) will be used.

- nfeatures:

  The maximum number of features to include in the heatmap. Default is
  `2000`.

- query_reduction:

  A character string specifying the name of a dimensionality reduction
  in the query object to use for calculating the distance metric.

- ref_reduction:

  A character string specifying the name of a dimensionality reduction
  in the reference object to use for calculating the distance metric.

- query_dims:

  A numeric vector specifying the dimension indices from the query
  reduction to be used for calculating the distance metric.

- ref_dims:

  A numeric vector specifying the dimension indices from the reference
  reduction to be used for calculating the distance metric.

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

  The distance metric to use for calculating the pairwise distances
  between cells. Options include: "pearson", "spearman", "cosine",
  "correlation", "jaccard", "ejaccard", "dice", "edice", "hamman",
  "simple matching", and "faith". Additional distance metrics can also
  be used, such as "euclidean", "manhattan", "hamming", etc.

- vote_fun:

  A character string specifying the function to be used for aggregating
  the nearest neighbors in the reference object. Options are "mean",
  "median", "sum", "min", "max", "sd", "var", etc. If not provided, the
  default is "mean".

## Value

A Seurat object with the projection results stored in the
"ref.embeddings" reduction. If `ref_group` is provided, the function
will also add a new metadata column called "predicted_ref_group" to the
query object.

## See also

[RunKNNPredict](https://mengxu98.github.io/scop/reference/RunKNNPredict.md),
[RunSingleR](https://mengxu98.github.io/scop/reference/RunSingleR.md),
[CellCorHeatmap](https://mengxu98.github.io/scop/reference/CellCorHeatmap.md)

## Examples

``` r
data(panc8_sub)
panc8_sub <- standard_scop(panc8_sub)
#> ℹ [2026-04-06 11:25:56] Start standard processing workflow...
#> ℹ [2026-04-06 11:25:57] Checking a list of <Seurat>...
#> ! [2026-04-06 11:25:57] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-06 11:25:57] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-06 11:25:59] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-06 11:26:00] Use the separate HVF from `srt_list`
#> ℹ [2026-04-06 11:26:00] Number of available HVF: 2000
#> ℹ [2026-04-06 11:26:00] Finished check
#> ℹ [2026-04-06 11:26:00] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-06 11:26:01] Perform pca linear dimension reduction
#> ℹ [2026-04-06 11:26:02] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-04-06 11:26:02] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-06 11:26:02] Reorder clusters...
#> ℹ [2026-04-06 11:26:03] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-06 11:26:03] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-06 11:26:03] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-04-06 11:26:08] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-04-06 11:26:13] Standard processing workflow completed
srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
srt_ref <- integration_scop(
  srt_ref,
  batch = "tech",
  integration_method = "Uncorrected"
)
#> ◌ [2026-04-06 11:26:13] Run integration workflow...
#> ℹ [2026-04-06 11:26:13] Split `srt_merge` into `srt_list` by "tech"
#> ℹ [2026-04-06 11:26:14] Checking a list of <Seurat>...
#> ℹ [2026-04-06 11:26:14] Data 1/4 of the `srt_list` has been log-normalized
#> ℹ [2026-04-06 11:26:14] Perform `Seurat::FindVariableFeatures()` on 1/4 of `srt_list`...
#> ℹ [2026-04-06 11:26:15] Data 2/4 of the `srt_list` has been log-normalized
#> ℹ [2026-04-06 11:26:15] Perform `Seurat::FindVariableFeatures()` on 2/4 of `srt_list`...
#> ℹ [2026-04-06 11:26:15] Data 3/4 of the `srt_list` has been log-normalized
#> ℹ [2026-04-06 11:26:15] Perform `Seurat::FindVariableFeatures()` on 3/4 of `srt_list`...
#> ℹ [2026-04-06 11:26:16] Data 4/4 of the `srt_list` has been log-normalized
#> ℹ [2026-04-06 11:26:16] Perform `Seurat::FindVariableFeatures()` on 4/4 of `srt_list`...
#> ℹ [2026-04-06 11:26:16] Use the separate HVF from `srt_list`
#> ℹ [2026-04-06 11:26:17] Number of available HVF: 2000
#> ℹ [2026-04-06 11:26:17] Finished check
#> ℹ [2026-04-06 11:26:18] Perform Uncorrected integration
#> ℹ [2026-04-06 11:26:19] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-06 11:26:19] Perform "pca" linear dimension reduction
#> ℹ [2026-04-06 11:26:21] Perform `Seurat::FindClusters()` with "louvain"
#> ℹ [2026-04-06 11:26:21] Reorder clusters...
#> ℹ [2026-04-06 11:26:23] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-06 11:26:23] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:20)
#> ℹ [2026-04-06 11:26:28] Perform umap nonlinear dimension reduction using Uncorrectedpca (1:20)
#> ✔ [2026-04-06 11:26:34] Uncorrected integration completed
CellDimPlot(
  srt_ref,
  group.by = c("celltype", "tech")
)


# Set the number of threads for RcppParallel
# details see: ?RcppParallel::setThreadOptions
# if (requireNamespace("RcppParallel", quietly = TRUE)) {
#   RcppParallel::setThreadOptions()
# }
# Projection
srt_query <- RunKNNMap(
  srt_query = srt_query,
  srt_ref = srt_ref,
  ref_umap = "UncorrectedUMAP2D"
)
#> ℹ [2026-04-06 11:26:35] Use the features to calculate distance metric
#> ℹ [2026-04-06 11:26:35] Data type is log-normalized
#> Warning: The following features were labelled as variable in 'var.features' but had no corresponding rank in `var.features.rank` and will therefore be ignored: 'A4GALT', 'AADAC', 'AARS2', 'ABCA1', 'ABCC3', 'ABCC4', 'ABCC8', 'ABCC9', 'ABCG2', 'ABHD15', 'ABI3', 'ABL2', 'ABTB2', 'ACACB', 'ACBD7', 'ACE', 'ACHE', 'ACSL3', 'ACSL4', 'ACTG1', 'ACTN1', 'ACTN4', 'ADAM12', 'ADAM9', 'ADAMTS1', 'ADAMTS2', 'ADAMTSL2', 'ADAT1', 'ADCK3', 'ADCY1', 'ADCY3', 'ADCY5', 'ADD3', 'ADM', 'ADRA2A', 'ADRA2C', 'ADRB1', 'ADRBK2', 'AFAP1L2', 'AGR2', 'AGRN', 'AHR', 'AIFM3', 'AIM1L', 'AIM1', 'AJAP1', 'AK5', 'AKAP7', 'ALDH1A1', 'ALDH1A2', 'ALDH2', 'ALDH3A2', 'ALDH3B1', 'ALDH4A1', 'ALOX5AP', 'ALPL', 'AMOTL1', 'AMOTL2', 'ANGPTL1', 'ANO6', 'ANO9', 'ANTXR2', 'AOC2', 'AP1S2', 'AP1S3', 'APOBEC2', 'APOD', 'AQP3', 'AQP4', 'ARHGAP26', 'ARHGAP29', 'ARHGAP6', 'ARHGEF17', 'ARHGEF2', 'ARHGEF3', 'ARHGEF40', 'ARHGEF6', 'ARL14', 'ARL4D', 'ARL6', 'ARMC9', 'ARNTL2', 'ARSJ', 'ARX', 'ASAP1', 'ASCL1', 'ASNS', 'ASPH', 'ASRGL1', 'ASS1', 'ASTN2', 'ASXL3', 'ATCAY', 'ATF3', 'ATF5', 'ATHL1', 'ATP10B', 'ATP11A', 'ATP1A1', 'ATP1B2', 'ATP2A3', 'ATP2B4', 'ATP8B1', 'AURKB', 'B3GNT2', 'B4GALT1', 'B4GALT5', 'BACE2', 'BACH2', 'BAIAP2L1', 'BAIAP3', 'BAMBI', 'BARX2', 'BASP1', 'BATF', 'BAZ1A', 'BBC3', 'BCHE', 'BCL2L15', 'BCL3', 'BCL6B', 'BCL6', 'BDH2', 'BHLHE41', 'BHMT2', 'BIK', 'BIRC5', 'BMF', 'BMP2', 'BMPR2', 'BNC2', 'BRIP1', 'BUB1', 'C1QL1', 'C1QTNF1', 'C1RL', 'C2CD4A', 'C2CD4B', 'CABP4', 'CACNA2D1', 'CADM2', 'CADPS2', 'CALCA', 'CALML4', 'CALU', 'CALY', 'CAMK1D', 'CAMKK1', 'CAMKK2', 'CAPN5', 'CAPN6', 'CARHSP1', 'CARTPT', 'CASP6', 'CAV2', 'CBR1', 'CCDC15', 'CCDC71L', 'CCL28', 'CCND2', 'CCNF', 'CCNG2', 'CD200', 'CD276', 'CD47', 'CD68', 'CD82', 'CDCA3', 'CDH17', 'CDH1', 'CDH23', 'CDH2', 'CDH3', 'CDHR2', 'CDK17', 'CDKN1A', 'CDKN2A', 'CDKN2B', 'CDR2', 'CDT1', 'CDX2', 'CEACAM1', 'CEBPB', 'CELF2', 'CENPF', 'CENPW', 'CERCAM', 'CFI', 'CFLAR', 'CGN', 'CHAC1', 'CHGB', 'CHN2', 'CHRM3', 'CHST1', 'CHSY1', 'CITED2', 'CITED4', 'CLCF1', 'CLCN4', 'CLDN11', 'CLDN19', 'CLDN3', 'CLDN5', 'CLDN6', 'CLDN7', 'CLEC11A', 'CLIC5', 'CLIC6', 'CLMN', 'CLMP', 'CLRN3', 'CLSPN', 'CMPK2', 'CMTM3', 'CMTM7', 'CNIH2', 'CNN2', 'CNN3', 'CNNM1', 'CNTN1', 'COCH', 'COL13A1', 'COL14A1', 'COL16A1', 'COL18A1', 'COL5A3', 'COL9A2', 'COMMD2', 'COMP', 'CORO1C', 'CORO2A', 'COX7A1', 'CPA3', 'CPLX2', 'CPM', 'CRABP2', 'CRADD', 'CREB3L1', 'CREB5', 'CREG2', 'CRIP1', 'CRLF1', 'CRTAC1', 'CRTAP', 'CTHRC1', 'CTNNA2', 'CTSB', 'CTSC', 'CTSS', 'CX3CL1', 'CXADR', 'CXCL12', 'CXCL16', 'CYB5A', 'CYBRD1', 'CYP20A1', 'CYP27B1', 'CYR61', 'CYTH3', 'DAB2IP', 'DACH1', 'DACT1', 'DAG1', 'DBN1', 'DCAF11', 'DCBLD2', 'DCUN1D2', 'DDAH2', 'DDC', 'DDIT3', 'DDX25', 'DDX51', 'DGKA', 'DHCR7', 'DHDH', 'DHODH', 'DHRS2', 'DIXDC1', 'DKK1', 'DMBT1', 'DMC1', 'DNAJB9', 'DNAJC22', 'DNER', 'DOCK4', 'DOCK5', 'DOK5', 'DPEP1', 'DPP9', 'DRAM1', 'DSC2', 'DSG2', 'DTNA', 'DUSP10', 'DUSP14', 'DUSP1', 'DUSP23', 'DUSP2', 'DUSP5', 'DZANK1', 'EBP', 'ECE1', 'ECEL1', 'EDN2', 'EFCAB7', 'EFNB1', 'EFNB3', 'EGFR', 'EGF', 'EGR1', 'EGR4', 'EHD1', 'EHD4', 'EHF', 'EIF4EBP1', 'ELAVL4', 'EMILIN1', 'EMP2', 'ENAH', 'ENG', 'ENPP2', 'EPAS1', 'EPDR1', 'EPHA2', 'EPHB2', 'EPHX1', 'EPS8L1', 'ERBB2', 'ERBB3', 'ERO1LB', 'ERO1L', 'ERRFI1', 'ESRRG', 'ETS2', 'ETV1', 'EVA1B', 'EVC2', 'EXT1', 'F10', 'F11R', 'F2RL1', 'F5', 'FA2H', 'FABP5', 'FADS1', 'FADS2', 'FAM101B', 'FAM105A', 'FAM107B', 'FAM129B', 'FAM134B', 'FAM155A', 'FAM159B', 'FAM161A', 'FAM162A', 'FAM163A', 'FAM20C', 'FAM217B', 'FAM227A', 'FAM43A', 'FAM46B', 'FAM46C', 'FAM73A', 'FAM84A', 'FAM84B', 'FBLIM1', 'FBLN1', 'FBLN5', 'FBP1', 'FBXL18', 'FBXO2', 'FBXO6', 'FCER1G', 'FDFT1', 'FERMT2', 'FEV', 'FEZ1', 'FFAR2', 'FFAR4', 'FGD2', 'FGD4', 'FGF7', 'FGFBP1', 'FGFR2', 'FGFR3', 'FGL2', 'FHL2', 'FKBP10', 'FKBP11', 'FKBP5', 'FLNA', 'FNDC4', 'FOLR1', 'FOSB', 'FOSL2', 'FOXC1', 'FOXM1', 'FOXQ1', 'FRAS1', 'FRMD4B', 'FSCN1', 'FXN', 'FXYD3', 'FXYD5', 'FXYD6', 'FZD5', 'G0S2', 'GAD2', 'GADD45A', 'GADD45B', 'GADD45G', 'GALE', 'GALM', 'GALNT2', 'GALNT3', 'GAL', 'GAS6', 'GAS7', 'GATA2', 'GCK', 'GCLC', 'GCNT1', 'GC', 'GDPD1', 'GFOD2', 'GFRA1', 'GGT1', 'GGT5', 'GHR', 'GINS3', 'GINS4', 'GJA1', 'GJA4', 'GJB1', 'GJB3', 'GJC1', 'GK5', 'GK', 'GLDC', 'GLI2', 'GLIPR2', 'GLT1D1', 'GLUL', 'GMDS', 'GMEB1', 'GMFG', 'GMNN', 'GNAL', 'GNPNAT1', 'GNPTAB', 'GOLM1', 'GOT1', 'GPC1', 'GPC6', 'GPR155', 'GPR160', 'GPR161', 'GPR4', 'GPRC5A', 'GPRC5B', 'GPSM1', 'GPT2', 'GPX3', 'GREB1', 'GSN', 'GSTM3', 'GSTM5', 'GSTP1', 'GTF2A1', 'GTPBP10', 'GUCA1B', 'GUCA2A', 'GUCY1A3', 'GULP1', 'GYLTL1B', 'H19', 'HABP2', 'HAP1', 'HBEGF', 'HCK', 'HDHD3', 'HEATR5A', 'HEBP1', 'HEG1', 'HEPH', 'HERPUD1', 'HES1', 'HIC1', 'HIF1A', 'HILPDA', 'HIP1', 'HIRA', 'HIST1H1C', 'HIST1H2BK', 'HIST2H2BE', 'HK1', 'HK2', 'HKDC1', 'HMGB3', 'HMGCS1', 'HMGN5', 'HMOX1', 'HN1', 'HNF1B', 'HNF4A', 'HNMT', 'HOMER2', 'HOPX', 'HOXB2', 'HOXB4', 'HPCAL1', 'HPGD', 'HPN', 'HS3ST1', 'HS3ST3A1', 'HS6ST2', 'HSD11B2', 'HSD17B11', 'HSPA1A', 'HSPA1B', 'HSPA5', 'HSPB6', 'HSPB8', 'HTRA1', 'HUNK', 'HYAL3', 'IBA57', 'ID1', 'ID4', 'IDH1', 'IDH2', 'IER2', 'IER3', 'IFI30', 'IFIT1', 'IFITM2', 'IFNGR1', 'IGFBP6', 'IGFBPL1', 'IGSF3', 'IKBIP', 'IKBKE', 'IKZF3', 'IL11', 'IL15RA', 'IL1RN', 'IL22RA1', 'IL33', 'IMPA2', 'INHBB', 'INPP4B', 'INSIG1', 'INSM1', 'IQGAP2', 'IRAK2', 'IRF1', 'IRS2', 'IRX1', 'IRX2', 'ISG20', 'ISL1', 'ITGA11', 'ITGA1', 'ITGA3', 'ITGA6', 'ITGAV', 'ITGB1', 'ITGB4', 'ITGB8', 'ITIH5', 'ITPR2', 'IVNS1ABP', 'IYD', 'JAM3', 'JARID2', 'JDP2', 'JMJD6', 'JPX', 'JUNB', 'JUN', 'JUP', 'KCNA5', 'KCNE3', 'KCNG1', 'KCNH6', 'KCNJ13', 'KCNJ5', 'KCNJ6', 'KCNJ8', 'KCNK16', 'KCNK5', 'KCNN2', 'KCNN3', 'KCNQ1OT1', 'KCNQ1', 'KCTD12', 'KDELC2', 'KDELR3', 'KIF11', 'KIF12', 'KIF20B', 'KIF21B', 'KIF5C', 'KIT', 'KLB', 'KLF2', 'KLF4', 'KLF5', 'KLF6', 'KLF7', 'KLHL28', 'KLHL31', 'KLHL3', 'KLK10', 'KLK11', 'KREMEN1', 'KRT15', 'KRT17', 'KRTCAP3', 'KYNU', 'L1TD1', 'L2HGDH', 'LAD1', 'LAIR1', 'LAMA1', 'LAMA3', 'LAMB1', 'LAMB2', 'LAMB3', 'LARP6', 'LASP1', 'LBH', 'LCP1', 'LDB2', 'LDLR', 'LEFTY1', 'LFNG', 'LGALS3', 'LGALS9', 'LHFPL2', 'LIFR', 'LIF', 'LIMA1', 'LIMK2', 'LIMS1', 'LIPC', 'LIPH', 'LITAF', 'LMO1', 'LMO2', 'LMO4', 'LMO7', 'LNX2', 'LOXL1', 'LOXL2', 'LPAR1', 'LRP1', 'LRP8', 'LRRC27', 'LRRC2', 'LRRC57', 'LRRC58', 'LRRIQ1', 'LRRK2', 'LRRTM2', 'LSR', 'LSS', 'LTBP4', 'LTB', 'LTF', 'LTV1', 'LURAP1L', 'LXN', 'LY6E', 'LY6H', 'LY96', 'LYN', 'LYPD1', 'LYRM7', 'LZTS1', 'MAB21L3', 'MACC1', 'MAFB', 'MAFF', 'MAL2', 'MALL', 'MAOA', 'MAOB', 'MAP3K9', 'MAPT', 'MARCKSL1', 'MAT1A', 'MCC', 'MCL1', 'MCM2', 'MCM6', 'MCOLN2', 'MDFIC', 'MDFI', 'MDGA1', 'MDK', 'MDM1', 'ME1', 'MEDAG', 'MEF2C', 'MEIS2', 'METRNL', 'METTL21A', 'MEX3A', 'MFAP2', 'MFAP4', 'MFGE8', 'MGAT4B', 'MGLL', 'MGP', 'MICAL2', 'MID1', 'MIR143HG', 'MKNK1', 'MLLT11', 'MLXIPL', 'MMACHC', 'MMP14', 'MNS1', 'MOB3B', 'MOSPD1', 'MOXD1', 'MPPED2', 'MPST', 'MPZL2', 'MRC2', 'MS4A2', 'MSMO1', 'MSRB1', 'MSRB3', 'MSX1', 'MTHFD2', 'MTUS1', 'MVP', 'MYADM', 'MYH9', 'MYL12A', 'MYL12B', 'MYL9', 'MYLK', 'MYO1C', 'MYO5B', 'MYOF', 'MYRF', 'MYZAP', 'N4BP3', 'NAA40', 'NAP1L2', 'NCAM1', 'NCF2', 'NCMAP', 'NCS1', 'NDRG1', 'NDRG2', 'NDUFA4L2', 'NDUFA4', 'NEDD9', 'NEK6', 'NEK8', 'NEURL3', 'NEUROD1', 'NEXN', 'NFATC4', 'NFKBIE', 'NFKBIZ', 'NHSL1', 'NICN1', 'NID1', 'NID2', 'NIPSNAP3B', 'NKD1', 'NLN', 'NLRC3', 'NMB', 'NOL9', 'NOTCH1', 'NOTCH3', 'NOTCH4', 'NPHS1', 'NPR1', 'NPTX2', 'NPY1R', 'NQO1', 'NR0B1', 'NR0B2', 'NR1D1', 'NR2F2', 'NR4A2', 'NR4A3', 'NR5A2', 'NR6A1', 'NRARP', 'NREP', 'NRP2', 'NSG1', 'NT5E', 'NTM', 'NTN4', 'NUAK2', 'NUGGC', 'NUSAP1', 'NYAP2', 'OAF', 'OAS3', 'OCLN', 'ODF2L', 'OGDH', 'OLFM2', 'OLFML1', 'OLFML2B', 'ORAI2', 'ORC4', 'OSGIN1', 'OSMR', 'OVOL1', 'OXR1', 'P2RX1', 'P4HA1', 'P4HA3', 'PACSIN2', 'PAG1', 'PAH', 'PALB2', 'PAMR1', 'PANK1', 'PAQR5', 'PARD6G', 'PARM1', 'PASK', 'PCDH10', 'PCDH17', 'PCDH18', 'PCDH7', 'PCK2', 'PCP4', 'PCSK2', 'PCSK5', 'PDCD4', 'PDDC1', 'PDE1A', 'PDE3B', 'PDE4A', 'PDE4C', 'PDE7B', 'PDE8B', 'PDGFA', 'PDGFB', 'PDGFC', 'PDIA4', 'PDLIM1', 'PDLIM2', 'PDLIM3', 'PDLIM4', 'PDP1', 'PDP2', 'PDX1', 'PEAK1', 'PECAM1', 'PEG10', 'PEMT', 'PGC', 'PGM1', 'PGM2L1', 'PHGDH', 'PHGR1', 'PHLDA1', 'PHLDA2', 'PHLDA3', 'PID1', 'PIEZO1', 'PIGN', 'PIK3IP1', 'PIK3R1', 'PIM1', 'PIM2', 'PIM3', 'PIWIL2', 'PKIB', 'PLA2G15', 'PLA2G4A', 'PLA2G4C', 'PLAT', 'PLBD1', 'PLCD1', 'PLCD3', 'P
#> ℹ [2026-04-06 11:26:36] Data type is log-normalized
#> Warning: The following features were labelled as variable in 'var.features' but had no corresponding rank in `var.features.rank` and will therefore be ignored: 'A4GALT', 'AADAC', 'AARS2', 'ABCA1', 'ABCC3', 'ABCC4', 'ABCC8', 'ABCC9', 'ABCG2', 'ABHD15', 'ABI3', 'ABL2', 'ABTB2', 'ACACB', 'ACBD7', 'ACE', 'ACHE', 'ACSL3', 'ACSL4', 'ACTG1', 'ACTN1', 'ACTN4', 'ADAM12', 'ADAM9', 'ADAMTS1', 'ADAMTS2', 'ADAMTSL2', 'ADAT1', 'ADCK3', 'ADCY1', 'ADCY3', 'ADCY5', 'ADD3', 'ADM', 'ADRA2A', 'ADRA2C', 'ADRB1', 'ADRBK2', 'AFAP1L2', 'AGR2', 'AGRN', 'AHR', 'AIFM3', 'AIM1L', 'AIM1', 'AJAP1', 'AK5', 'AKAP7', 'ALDH1A1', 'ALDH1A2', 'ALDH2', 'ALDH3A2', 'ALDH3B1', 'ALDH4A1', 'ALOX5AP', 'ALPL', 'AMOTL1', 'AMOTL2', 'ANGPTL1', 'ANO6', 'ANO9', 'ANTXR2', 'AOC2', 'AP1S2', 'AP1S3', 'APOBEC2', 'APOD', 'AQP3', 'AQP4', 'ARHGAP26', 'ARHGAP29', 'ARHGAP6', 'ARHGEF17', 'ARHGEF2', 'ARHGEF3', 'ARHGEF40', 'ARHGEF6', 'ARL14', 'ARL4D', 'ARL6', 'ARMC9', 'ARNTL2', 'ARSJ', 'ARX', 'ASAP1', 'ASCL1', 'ASNS', 'ASPH', 'ASRGL1', 'ASS1', 'ASTN2', 'ASXL3', 'ATCAY', 'ATF3', 'ATF5', 'ATHL1', 'ATP10B', 'ATP11A', 'ATP1A1', 'ATP1B2', 'ATP2A3', 'ATP2B4', 'ATP8B1', 'AURKB', 'B3GNT2', 'B4GALT1', 'B4GALT5', 'BACE2', 'BACH2', 'BAIAP2L1', 'BAIAP3', 'BAMBI', 'BARX2', 'BASP1', 'BATF', 'BAZ1A', 'BBC3', 'BCHE', 'BCL2L15', 'BCL3', 'BCL6B', 'BCL6', 'BDH2', 'BHLHE41', 'BHMT2', 'BIK', 'BIRC5', 'BMF', 'BMP2', 'BMPR2', 'BNC2', 'BRIP1', 'BUB1', 'C1QL1', 'C1QTNF1', 'C1RL', 'C2CD4A', 'C2CD4B', 'CABP4', 'CACNA2D1', 'CADM2', 'CADPS2', 'CALCA', 'CALML4', 'CALU', 'CALY', 'CAMK1D', 'CAMKK1', 'CAMKK2', 'CAPN5', 'CAPN6', 'CARHSP1', 'CARTPT', 'CASP6', 'CAV2', 'CBR1', 'CCDC15', 'CCDC71L', 'CCL28', 'CCND2', 'CCNF', 'CCNG2', 'CD200', 'CD276', 'CD47', 'CD68', 'CD82', 'CDCA3', 'CDH17', 'CDH1', 'CDH23', 'CDH2', 'CDH3', 'CDHR2', 'CDK17', 'CDKN1A', 'CDKN2A', 'CDKN2B', 'CDR2', 'CDT1', 'CDX2', 'CEACAM1', 'CEBPB', 'CELF2', 'CENPF', 'CENPW', 'CERCAM', 'CFI', 'CFLAR', 'CGN', 'CHAC1', 'CHGB', 'CHN2', 'CHRM3', 'CHST1', 'CHSY1', 'CITED2', 'CITED4', 'CLCF1', 'CLCN4', 'CLDN11', 'CLDN19', 'CLDN3', 'CLDN5', 'CLDN6', 'CLDN7', 'CLEC11A', 'CLIC5', 'CLIC6', 'CLMN', 'CLMP', 'CLRN3', 'CLSPN', 'CMPK2', 'CMTM3', 'CMTM7', 'CNIH2', 'CNN2', 'CNN3', 'CNNM1', 'CNTN1', 'COCH', 'COL13A1', 'COL14A1', 'COL16A1', 'COL18A1', 'COL5A3', 'COL9A2', 'COMMD2', 'COMP', 'CORO1C', 'CORO2A', 'COX7A1', 'CPA3', 'CPLX2', 'CPM', 'CRABP2', 'CRADD', 'CREB3L1', 'CREB5', 'CREG2', 'CRIP1', 'CRLF1', 'CRTAC1', 'CRTAP', 'CTHRC1', 'CTNNA2', 'CTSB', 'CTSC', 'CTSS', 'CX3CL1', 'CXADR', 'CXCL12', 'CXCL16', 'CYB5A', 'CYBRD1', 'CYP20A1', 'CYP27B1', 'CYR61', 'CYTH3', 'DAB2IP', 'DACH1', 'DACT1', 'DAG1', 'DBN1', 'DCAF11', 'DCBLD2', 'DCUN1D2', 'DDAH2', 'DDC', 'DDIT3', 'DDX25', 'DDX51', 'DGKA', 'DHCR7', 'DHDH', 'DHODH', 'DHRS2', 'DIXDC1', 'DKK1', 'DMBT1', 'DMC1', 'DNAJB9', 'DNAJC22', 'DNER', 'DOCK4', 'DOCK5', 'DOK5', 'DPEP1', 'DPP9', 'DRAM1', 'DSC2', 'DSG2', 'DTNA', 'DUSP10', 'DUSP14', 'DUSP1', 'DUSP23', 'DUSP2', 'DUSP5', 'DZANK1', 'EBP', 'ECE1', 'ECEL1', 'EDN2', 'EFCAB7', 'EFNB1', 'EFNB3', 'EGFR', 'EGF', 'EGR1', 'EGR4', 'EHD1', 'EHD4', 'EHF', 'EIF4EBP1', 'ELAVL4', 'EMILIN1', 'EMP2', 'ENAH', 'ENG', 'ENPP2', 'EPAS1', 'EPDR1', 'EPHA2', 'EPHB2', 'EPHX1', 'EPS8L1', 'ERBB2', 'ERBB3', 'ERO1LB', 'ERO1L', 'ERRFI1', 'ESRRG', 'ETS2', 'ETV1', 'EVA1B', 'EVC2', 'EXT1', 'F10', 'F11R', 'F2RL1', 'F5', 'FA2H', 'FABP5', 'FADS1', 'FADS2', 'FAM101B', 'FAM105A', 'FAM107B', 'FAM129B', 'FAM134B', 'FAM155A', 'FAM159B', 'FAM161A', 'FAM162A', 'FAM163A', 'FAM20C', 'FAM217B', 'FAM227A', 'FAM43A', 'FAM46B', 'FAM46C', 'FAM73A', 'FAM84A', 'FAM84B', 'FBLIM1', 'FBLN1', 'FBLN5', 'FBP1', 'FBXL18', 'FBXO2', 'FBXO6', 'FCER1G', 'FDFT1', 'FERMT2', 'FEV', 'FEZ1', 'FFAR2', 'FFAR4', 'FGD2', 'FGD4', 'FGF7', 'FGFBP1', 'FGFR2', 'FGFR3', 'FGL2', 'FHL2', 'FKBP10', 'FKBP11', 'FKBP5', 'FLNA', 'FNDC4', 'FOLR1', 'FOSB', 'FOSL2', 'FOXC1', 'FOXM1', 'FOXQ1', 'FRAS1', 'FRMD4B', 'FSCN1', 'FXN', 'FXYD3', 'FXYD5', 'FXYD6', 'FZD5', 'G0S2', 'GAD2', 'GADD45A', 'GADD45B', 'GADD45G', 'GALE', 'GALM', 'GALNT2', 'GALNT3', 'GAL', 'GAS6', 'GAS7', 'GATA2', 'GCK', 'GCLC', 'GCNT1', 'GC', 'GDPD1', 'GFOD2', 'GFRA1', 'GGT1', 'GGT5', 'GHR', 'GINS3', 'GINS4', 'GJA1', 'GJA4', 'GJB1', 'GJB3', 'GJC1', 'GK5', 'GK', 'GLDC', 'GLI2', 'GLIPR2', 'GLT1D1', 'GLUL', 'GMDS', 'GMEB1', 'GMFG', 'GMNN', 'GNAL', 'GNPNAT1', 'GNPTAB', 'GOLM1', 'GOT1', 'GPC1', 'GPC6', 'GPR155', 'GPR160', 'GPR161', 'GPR4', 'GPRC5A', 'GPRC5B', 'GPSM1', 'GPT2', 'GPX3', 'GREB1', 'GSN', 'GSTM3', 'GSTM5', 'GSTP1', 'GTF2A1', 'GTPBP10', 'GUCA1B', 'GUCA2A', 'GUCY1A3', 'GULP1', 'GYLTL1B', 'H19', 'HABP2', 'HAP1', 'HBEGF', 'HCK', 'HDHD3', 'HEATR5A', 'HEBP1', 'HEG1', 'HEPH', 'HERPUD1', 'HES1', 'HIC1', 'HIF1A', 'HILPDA', 'HIP1', 'HIRA', 'HIST1H1C', 'HIST1H2BK', 'HIST2H2BE', 'HK1', 'HK2', 'HKDC1', 'HMGB3', 'HMGCS1', 'HMGN5', 'HMOX1', 'HN1', 'HNF1B', 'HNF4A', 'HNMT', 'HOMER2', 'HOPX', 'HOXB2', 'HOXB4', 'HPCAL1', 'HPGD', 'HPN', 'HS3ST1', 'HS3ST3A1', 'HS6ST2', 'HSD11B2', 'HSD17B11', 'HSPA1A', 'HSPA1B', 'HSPA5', 'HSPB6', 'HSPB8', 'HTRA1', 'HUNK', 'HYAL3', 'IBA57', 'ID1', 'ID4', 'IDH1', 'IDH2', 'IER2', 'IER3', 'IFI30', 'IFIT1', 'IFITM2', 'IFNGR1', 'IGFBP6', 'IGFBPL1', 'IGSF3', 'IKBIP', 'IKBKE', 'IKZF3', 'IL11', 'IL15RA', 'IL1RN', 'IL22RA1', 'IL33', 'IMPA2', 'INHBB', 'INPP4B', 'INSIG1', 'INSM1', 'IQGAP2', 'IRAK2', 'IRF1', 'IRS2', 'IRX1', 'IRX2', 'ISG20', 'ISL1', 'ITGA11', 'ITGA1', 'ITGA3', 'ITGA6', 'ITGAV', 'ITGB1', 'ITGB4', 'ITGB8', 'ITIH5', 'ITPR2', 'IVNS1ABP', 'IYD', 'JAM3', 'JARID2', 'JDP2', 'JMJD6', 'JPX', 'JUNB', 'JUN', 'JUP', 'KCNA5', 'KCNE3', 'KCNG1', 'KCNH6', 'KCNJ13', 'KCNJ5', 'KCNJ6', 'KCNJ8', 'KCNK16', 'KCNK5', 'KCNN2', 'KCNN3', 'KCNQ1OT1', 'KCNQ1', 'KCTD12', 'KDELC2', 'KDELR3', 'KIF11', 'KIF12', 'KIF20B', 'KIF21B', 'KIF5C', 'KIT', 'KLB', 'KLF2', 'KLF4', 'KLF5', 'KLF6', 'KLF7', 'KLHL28', 'KLHL31', 'KLHL3', 'KLK10', 'KLK11', 'KREMEN1', 'KRT15', 'KRT17', 'KRTCAP3', 'KYNU', 'L1TD1', 'L2HGDH', 'LAD1', 'LAIR1', 'LAMA1', 'LAMA3', 'LAMB1', 'LAMB2', 'LAMB3', 'LARP6', 'LASP1', 'LBH', 'LCP1', 'LDB2', 'LDLR', 'LEFTY1', 'LFNG', 'LGALS3', 'LGALS9', 'LHFPL2', 'LIFR', 'LIF', 'LIMA1', 'LIMK2', 'LIMS1', 'LIPC', 'LIPH', 'LITAF', 'LMO1', 'LMO2', 'LMO4', 'LMO7', 'LNX2', 'LOXL1', 'LOXL2', 'LPAR1', 'LRP1', 'LRP8', 'LRRC27', 'LRRC2', 'LRRC57', 'LRRC58', 'LRRIQ1', 'LRRK2', 'LRRTM2', 'LSR', 'LSS', 'LTBP4', 'LTB', 'LTF', 'LTV1', 'LURAP1L', 'LXN', 'LY6E', 'LY6H', 'LY96', 'LYN', 'LYPD1', 'LYRM7', 'LZTS1', 'MAB21L3', 'MACC1', 'MAFB', 'MAFF', 'MAL2', 'MALL', 'MAOA', 'MAOB', 'MAP3K9', 'MAPT', 'MARCKSL1', 'MAT1A', 'MCC', 'MCL1', 'MCM2', 'MCM6', 'MCOLN2', 'MDFIC', 'MDFI', 'MDGA1', 'MDK', 'MDM1', 'ME1', 'MEDAG', 'MEF2C', 'MEIS2', 'METRNL', 'METTL21A', 'MEX3A', 'MFAP2', 'MFAP4', 'MFGE8', 'MGAT4B', 'MGLL', 'MGP', 'MICAL2', 'MID1', 'MIR143HG', 'MKNK1', 'MLLT11', 'MLXIPL', 'MMACHC', 'MMP14', 'MNS1', 'MOB3B', 'MOSPD1', 'MOXD1', 'MPPED2', 'MPST', 'MPZL2', 'MRC2', 'MS4A2', 'MSMO1', 'MSRB1', 'MSRB3', 'MSX1', 'MTHFD2', 'MTUS1', 'MVP', 'MYADM', 'MYH9', 'MYL12A', 'MYL12B', 'MYL9', 'MYLK', 'MYO1C', 'MYO5B', 'MYOF', 'MYRF', 'MYZAP', 'N4BP3', 'NAA40', 'NAP1L2', 'NCAM1', 'NCF2', 'NCMAP', 'NCS1', 'NDRG1', 'NDRG2', 'NDUFA4L2', 'NDUFA4', 'NEDD9', 'NEK6', 'NEK8', 'NEURL3', 'NEUROD1', 'NEXN', 'NFATC4', 'NFKBIE', 'NFKBIZ', 'NHSL1', 'NICN1', 'NID1', 'NID2', 'NIPSNAP3B', 'NKD1', 'NLN', 'NLRC3', 'NMB', 'NOL9', 'NOTCH1', 'NOTCH3', 'NOTCH4', 'NPHS1', 'NPR1', 'NPTX2', 'NPY1R', 'NQO1', 'NR0B1', 'NR0B2', 'NR1D1', 'NR2F2', 'NR4A2', 'NR4A3', 'NR5A2', 'NR6A1', 'NRARP', 'NREP', 'NRP2', 'NSG1', 'NT5E', 'NTM', 'NTN4', 'NUAK2', 'NUGGC', 'NUSAP1', 'NYAP2', 'OAF', 'OAS3', 'OCLN', 'ODF2L', 'OGDH', 'OLFM2', 'OLFML1', 'OLFML2B', 'ORAI2', 'ORC4', 'OSGIN1', 'OSMR', 'OVOL1', 'OXR1', 'P2RX1', 'P4HA1', 'P4HA3', 'PACSIN2', 'PAG1', 'PAH', 'PALB2', 'PAMR1', 'PANK1', 'PAQR5', 'PARD6G', 'PARM1', 'PASK', 'PCDH10', 'PCDH17', 'PCDH18', 'PCDH7', 'PCK2', 'PCP4', 'PCSK2', 'PCSK5', 'PDCD4', 'PDDC1', 'PDE1A', 'PDE3B', 'PDE4A', 'PDE4C', 'PDE7B', 'PDE8B', 'PDGFA', 'PDGFB', 'PDGFC', 'PDIA4', 'PDLIM1', 'PDLIM2', 'PDLIM3', 'PDLIM4', 'PDP1', 'PDP2', 'PDX1', 'PEAK1', 'PECAM1', 'PEG10', 'PEMT', 'PGC', 'PGM1', 'PGM2L1', 'PHGDH', 'PHGR1', 'PHLDA1', 'PHLDA2', 'PHLDA3', 'PID1', 'PIEZO1', 'PIGN', 'PIK3IP1', 'PIK3R1', 'PIM1', 'PIM2', 'PIM3', 'PIWIL2', 'PKIB', 'PLA2G15', 'PLA2G4A', 'PLA2G4C', 'PLAT', 'PLBD1', 'PLCD1', 'PLCD3', 'P
#> ℹ [2026-04-06 11:26:36] Use 2000 features to calculate distance
#> Warning: The following features were labelled as variable in 'var.features' but had no corresponding rank in `var.features.rank` and will therefore be ignored: 'A4GALT', 'AADAC', 'AARS2', 'ABCA1', 'ABCC3', 'ABCC4', 'ABCC8', 'ABCC9', 'ABCG2', 'ABHD15', 'ABI3', 'ABL2', 'ABTB2', 'ACACB', 'ACBD7', 'ACE', 'ACHE', 'ACSL3', 'ACSL4', 'ACTG1', 'ACTN1', 'ACTN4', 'ADAM12', 'ADAM9', 'ADAMTS1', 'ADAMTS2', 'ADAMTSL2', 'ADAT1', 'ADCK3', 'ADCY1', 'ADCY3', 'ADCY5', 'ADD3', 'ADM', 'ADRA2A', 'ADRA2C', 'ADRB1', 'ADRBK2', 'AFAP1L2', 'AGR2', 'AGRN', 'AHR', 'AIFM3', 'AIM1L', 'AIM1', 'AJAP1', 'AK5', 'AKAP7', 'ALDH1A1', 'ALDH1A2', 'ALDH2', 'ALDH3A2', 'ALDH3B1', 'ALDH4A1', 'ALOX5AP', 'ALPL', 'AMOTL1', 'AMOTL2', 'ANGPTL1', 'ANO6', 'ANO9', 'ANTXR2', 'AOC2', 'AP1S2', 'AP1S3', 'APOBEC2', 'APOD', 'AQP3', 'AQP4', 'ARHGAP26', 'ARHGAP29', 'ARHGAP6', 'ARHGEF17', 'ARHGEF2', 'ARHGEF3', 'ARHGEF40', 'ARHGEF6', 'ARL14', 'ARL4D', 'ARL6', 'ARMC9', 'ARNTL2', 'ARSJ', 'ARX', 'ASAP1', 'ASCL1', 'ASNS', 'ASPH', 'ASRGL1', 'ASS1', 'ASTN2', 'ASXL3', 'ATCAY', 'ATF3', 'ATF5', 'ATHL1', 'ATP10B', 'ATP11A', 'ATP1A1', 'ATP1B2', 'ATP2A3', 'ATP2B4', 'ATP8B1', 'AURKB', 'B3GNT2', 'B4GALT1', 'B4GALT5', 'BACE2', 'BACH2', 'BAIAP2L1', 'BAIAP3', 'BAMBI', 'BARX2', 'BASP1', 'BATF', 'BAZ1A', 'BBC3', 'BCHE', 'BCL2L15', 'BCL3', 'BCL6B', 'BCL6', 'BDH2', 'BHLHE41', 'BHMT2', 'BIK', 'BIRC5', 'BMF', 'BMP2', 'BMPR2', 'BNC2', 'BRIP1', 'BUB1', 'C1QL1', 'C1QTNF1', 'C1RL', 'C2CD4A', 'C2CD4B', 'CABP4', 'CACNA2D1', 'CADM2', 'CADPS2', 'CALCA', 'CALML4', 'CALU', 'CALY', 'CAMK1D', 'CAMKK1', 'CAMKK2', 'CAPN5', 'CAPN6', 'CARHSP1', 'CARTPT', 'CASP6', 'CAV2', 'CBR1', 'CCDC15', 'CCDC71L', 'CCL28', 'CCND2', 'CCNF', 'CCNG2', 'CD200', 'CD276', 'CD47', 'CD68', 'CD82', 'CDCA3', 'CDH17', 'CDH1', 'CDH23', 'CDH2', 'CDH3', 'CDHR2', 'CDK17', 'CDKN1A', 'CDKN2A', 'CDKN2B', 'CDR2', 'CDT1', 'CDX2', 'CEACAM1', 'CEBPB', 'CELF2', 'CENPF', 'CENPW', 'CERCAM', 'CFI', 'CFLAR', 'CGN', 'CHAC1', 'CHGB', 'CHN2', 'CHRM3', 'CHST1', 'CHSY1', 'CITED2', 'CITED4', 'CLCF1', 'CLCN4', 'CLDN11', 'CLDN19', 'CLDN3', 'CLDN5', 'CLDN6', 'CLDN7', 'CLEC11A', 'CLIC5', 'CLIC6', 'CLMN', 'CLMP', 'CLRN3', 'CLSPN', 'CMPK2', 'CMTM3', 'CMTM7', 'CNIH2', 'CNN2', 'CNN3', 'CNNM1', 'CNTN1', 'COCH', 'COL13A1', 'COL14A1', 'COL16A1', 'COL18A1', 'COL5A3', 'COL9A2', 'COMMD2', 'COMP', 'CORO1C', 'CORO2A', 'COX7A1', 'CPA3', 'CPLX2', 'CPM', 'CRABP2', 'CRADD', 'CREB3L1', 'CREB5', 'CREG2', 'CRIP1', 'CRLF1', 'CRTAC1', 'CRTAP', 'CTHRC1', 'CTNNA2', 'CTSB', 'CTSC', 'CTSS', 'CX3CL1', 'CXADR', 'CXCL12', 'CXCL16', 'CYB5A', 'CYBRD1', 'CYP20A1', 'CYP27B1', 'CYR61', 'CYTH3', 'DAB2IP', 'DACH1', 'DACT1', 'DAG1', 'DBN1', 'DCAF11', 'DCBLD2', 'DCUN1D2', 'DDAH2', 'DDC', 'DDIT3', 'DDX25', 'DDX51', 'DGKA', 'DHCR7', 'DHDH', 'DHODH', 'DHRS2', 'DIXDC1', 'DKK1', 'DMBT1', 'DMC1', 'DNAJB9', 'DNAJC22', 'DNER', 'DOCK4', 'DOCK5', 'DOK5', 'DPEP1', 'DPP9', 'DRAM1', 'DSC2', 'DSG2', 'DTNA', 'DUSP10', 'DUSP14', 'DUSP1', 'DUSP23', 'DUSP2', 'DUSP5', 'DZANK1', 'EBP', 'ECE1', 'ECEL1', 'EDN2', 'EFCAB7', 'EFNB1', 'EFNB3', 'EGFR', 'EGF', 'EGR1', 'EGR4', 'EHD1', 'EHD4', 'EHF', 'EIF4EBP1', 'ELAVL4', 'EMILIN1', 'EMP2', 'ENAH', 'ENG', 'ENPP2', 'EPAS1', 'EPDR1', 'EPHA2', 'EPHB2', 'EPHX1', 'EPS8L1', 'ERBB2', 'ERBB3', 'ERO1LB', 'ERO1L', 'ERRFI1', 'ESRRG', 'ETS2', 'ETV1', 'EVA1B', 'EVC2', 'EXT1', 'F10', 'F11R', 'F2RL1', 'F5', 'FA2H', 'FABP5', 'FADS1', 'FADS2', 'FAM101B', 'FAM105A', 'FAM107B', 'FAM129B', 'FAM134B', 'FAM155A', 'FAM159B', 'FAM161A', 'FAM162A', 'FAM163A', 'FAM20C', 'FAM217B', 'FAM227A', 'FAM43A', 'FAM46B', 'FAM46C', 'FAM73A', 'FAM84A', 'FAM84B', 'FBLIM1', 'FBLN1', 'FBLN5', 'FBP1', 'FBXL18', 'FBXO2', 'FBXO6', 'FCER1G', 'FDFT1', 'FERMT2', 'FEV', 'FEZ1', 'FFAR2', 'FFAR4', 'FGD2', 'FGD4', 'FGF7', 'FGFBP1', 'FGFR2', 'FGFR3', 'FGL2', 'FHL2', 'FKBP10', 'FKBP11', 'FKBP5', 'FLNA', 'FNDC4', 'FOLR1', 'FOSB', 'FOSL2', 'FOXC1', 'FOXM1', 'FOXQ1', 'FRAS1', 'FRMD4B', 'FSCN1', 'FXN', 'FXYD3', 'FXYD5', 'FXYD6', 'FZD5', 'G0S2', 'GAD2', 'GADD45A', 'GADD45B', 'GADD45G', 'GALE', 'GALM', 'GALNT2', 'GALNT3', 'GAL', 'GAS6', 'GAS7', 'GATA2', 'GCK', 'GCLC', 'GCNT1', 'GC', 'GDPD1', 'GFOD2', 'GFRA1', 'GGT1', 'GGT5', 'GHR', 'GINS3', 'GINS4', 'GJA1', 'GJA4', 'GJB1', 'GJB3', 'GJC1', 'GK5', 'GK', 'GLDC', 'GLI2', 'GLIPR2', 'GLT1D1', 'GLUL', 'GMDS', 'GMEB1', 'GMFG', 'GMNN', 'GNAL', 'GNPNAT1', 'GNPTAB', 'GOLM1', 'GOT1', 'GPC1', 'GPC6', 'GPR155', 'GPR160', 'GPR161', 'GPR4', 'GPRC5A', 'GPRC5B', 'GPSM1', 'GPT2', 'GPX3', 'GREB1', 'GSN', 'GSTM3', 'GSTM5', 'GSTP1', 'GTF2A1', 'GTPBP10', 'GUCA1B', 'GUCA2A', 'GUCY1A3', 'GULP1', 'GYLTL1B', 'H19', 'HABP2', 'HAP1', 'HBEGF', 'HCK', 'HDHD3', 'HEATR5A', 'HEBP1', 'HEG1', 'HEPH', 'HERPUD1', 'HES1', 'HIC1', 'HIF1A', 'HILPDA', 'HIP1', 'HIRA', 'HIST1H1C', 'HIST1H2BK', 'HIST2H2BE', 'HK1', 'HK2', 'HKDC1', 'HMGB3', 'HMGCS1', 'HMGN5', 'HMOX1', 'HN1', 'HNF1B', 'HNF4A', 'HNMT', 'HOMER2', 'HOPX', 'HOXB2', 'HOXB4', 'HPCAL1', 'HPGD', 'HPN', 'HS3ST1', 'HS3ST3A1', 'HS6ST2', 'HSD11B2', 'HSD17B11', 'HSPA1A', 'HSPA1B', 'HSPA5', 'HSPB6', 'HSPB8', 'HTRA1', 'HUNK', 'HYAL3', 'IBA57', 'ID1', 'ID4', 'IDH1', 'IDH2', 'IER2', 'IER3', 'IFI30', 'IFIT1', 'IFITM2', 'IFNGR1', 'IGFBP6', 'IGFBPL1', 'IGSF3', 'IKBIP', 'IKBKE', 'IKZF3', 'IL11', 'IL15RA', 'IL1RN', 'IL22RA1', 'IL33', 'IMPA2', 'INHBB', 'INPP4B', 'INSIG1', 'INSM1', 'IQGAP2', 'IRAK2', 'IRF1', 'IRS2', 'IRX1', 'IRX2', 'ISG20', 'ISL1', 'ITGA11', 'ITGA1', 'ITGA3', 'ITGA6', 'ITGAV', 'ITGB1', 'ITGB4', 'ITGB8', 'ITIH5', 'ITPR2', 'IVNS1ABP', 'IYD', 'JAM3', 'JARID2', 'JDP2', 'JMJD6', 'JPX', 'JUNB', 'JUN', 'JUP', 'KCNA5', 'KCNE3', 'KCNG1', 'KCNH6', 'KCNJ13', 'KCNJ5', 'KCNJ6', 'KCNJ8', 'KCNK16', 'KCNK5', 'KCNN2', 'KCNN3', 'KCNQ1OT1', 'KCNQ1', 'KCTD12', 'KDELC2', 'KDELR3', 'KIF11', 'KIF12', 'KIF20B', 'KIF21B', 'KIF5C', 'KIT', 'KLB', 'KLF2', 'KLF4', 'KLF5', 'KLF6', 'KLF7', 'KLHL28', 'KLHL31', 'KLHL3', 'KLK10', 'KLK11', 'KREMEN1', 'KRT15', 'KRT17', 'KRTCAP3', 'KYNU', 'L1TD1', 'L2HGDH', 'LAD1', 'LAIR1', 'LAMA1', 'LAMA3', 'LAMB1', 'LAMB2', 'LAMB3', 'LARP6', 'LASP1', 'LBH', 'LCP1', 'LDB2', 'LDLR', 'LEFTY1', 'LFNG', 'LGALS3', 'LGALS9', 'LHFPL2', 'LIFR', 'LIF', 'LIMA1', 'LIMK2', 'LIMS1', 'LIPC', 'LIPH', 'LITAF', 'LMO1', 'LMO2', 'LMO4', 'LMO7', 'LNX2', 'LOXL1', 'LOXL2', 'LPAR1', 'LRP1', 'LRP8', 'LRRC27', 'LRRC2', 'LRRC57', 'LRRC58', 'LRRIQ1', 'LRRK2', 'LRRTM2', 'LSR', 'LSS', 'LTBP4', 'LTB', 'LTF', 'LTV1', 'LURAP1L', 'LXN', 'LY6E', 'LY6H', 'LY96', 'LYN', 'LYPD1', 'LYRM7', 'LZTS1', 'MAB21L3', 'MACC1', 'MAFB', 'MAFF', 'MAL2', 'MALL', 'MAOA', 'MAOB', 'MAP3K9', 'MAPT', 'MARCKSL1', 'MAT1A', 'MCC', 'MCL1', 'MCM2', 'MCM6', 'MCOLN2', 'MDFIC', 'MDFI', 'MDGA1', 'MDK', 'MDM1', 'ME1', 'MEDAG', 'MEF2C', 'MEIS2', 'METRNL', 'METTL21A', 'MEX3A', 'MFAP2', 'MFAP4', 'MFGE8', 'MGAT4B', 'MGLL', 'MGP', 'MICAL2', 'MID1', 'MIR143HG', 'MKNK1', 'MLLT11', 'MLXIPL', 'MMACHC', 'MMP14', 'MNS1', 'MOB3B', 'MOSPD1', 'MOXD1', 'MPPED2', 'MPST', 'MPZL2', 'MRC2', 'MS4A2', 'MSMO1', 'MSRB1', 'MSRB3', 'MSX1', 'MTHFD2', 'MTUS1', 'MVP', 'MYADM', 'MYH9', 'MYL12A', 'MYL12B', 'MYL9', 'MYLK', 'MYO1C', 'MYO5B', 'MYOF', 'MYRF', 'MYZAP', 'N4BP3', 'NAA40', 'NAP1L2', 'NCAM1', 'NCF2', 'NCMAP', 'NCS1', 'NDRG1', 'NDRG2', 'NDUFA4L2', 'NDUFA4', 'NEDD9', 'NEK6', 'NEK8', 'NEURL3', 'NEUROD1', 'NEXN', 'NFATC4', 'NFKBIE', 'NFKBIZ', 'NHSL1', 'NICN1', 'NID1', 'NID2', 'NIPSNAP3B', 'NKD1', 'NLN', 'NLRC3', 'NMB', 'NOL9', 'NOTCH1', 'NOTCH3', 'NOTCH4', 'NPHS1', 'NPR1', 'NPTX2', 'NPY1R', 'NQO1', 'NR0B1', 'NR0B2', 'NR1D1', 'NR2F2', 'NR4A2', 'NR4A3', 'NR5A2', 'NR6A1', 'NRARP', 'NREP', 'NRP2', 'NSG1', 'NT5E', 'NTM', 'NTN4', 'NUAK2', 'NUGGC', 'NUSAP1', 'NYAP2', 'OAF', 'OAS3', 'OCLN', 'ODF2L', 'OGDH', 'OLFM2', 'OLFML1', 'OLFML2B', 'ORAI2', 'ORC4', 'OSGIN1', 'OSMR', 'OVOL1', 'OXR1', 'P2RX1', 'P4HA1', 'P4HA3', 'PACSIN2', 'PAG1', 'PAH', 'PALB2', 'PAMR1', 'PANK1', 'PAQR5', 'PARD6G', 'PARM1', 'PASK', 'PCDH10', 'PCDH17', 'PCDH18', 'PCDH7', 'PCK2', 'PCP4', 'PCSK2', 'PCSK5', 'PDCD4', 'PDDC1', 'PDE1A', 'PDE3B', 'PDE4A', 'PDE4C', 'PDE7B', 'PDE8B', 'PDGFA', 'PDGFB', 'PDGFC', 'PDIA4', 'PDLIM1', 'PDLIM2', 'PDLIM3', 'PDLIM4', 'PDP1', 'PDP2', 'PDX1', 'PEAK1', 'PECAM1', 'PEG10', 'PEMT', 'PGC', 'PGM1', 'PGM2L1', 'PHGDH', 'PHGR1', 'PHLDA1', 'PHLDA2', 'PHLDA3', 'PID1', 'PIEZO1', 'PIGN', 'PIK3IP1', 'PIK3R1', 'PIM1', 'PIM2', 'PIM3', 'PIWIL2', 'PKIB', 'PLA2G15', 'PLA2G4A', 'PLA2G4C', 'PLAT', 'PLBD1', 'PLCD1', 'PLCD3', 'P
#> ℹ [2026-04-06 11:26:36] Use raw method to find neighbors
#> ℹ [2026-04-06 11:26:36] Running UMAP projection
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
