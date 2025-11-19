# Projection Plot

This function generates a projection plot, which can be used to compare
two groups of cells in a dimensionality reduction space.

## Usage

``` r
ProjectionPlot(
  srt_query,
  srt_ref,
  query_group = NULL,
  ref_group = NULL,
  query_reduction = "ref.embeddings",
  ref_reduction = srt_query[[query_reduction]]@misc[["reduction.model"]] %||% NULL,
  query_param = list(palette = "Set1", cells.highlight = TRUE),
  ref_param = list(palette = "Paired"),
  xlim = NULL,
  ylim = NULL,
  pt.size = 0.8,
  stroke.highlight = 0.5
)
```

## Arguments

- srt_query:

  An object of class Seurat storing the query group cells.

- srt_ref:

  An object of class Seurat storing the reference group cells.

- query_group:

  The grouping variable for the query group cells.

- ref_group:

  The grouping variable for the reference group cells.

- query_reduction:

  The name of the reduction in the query group cells.

- ref_reduction:

  The name of the reduction in the reference group cells.

- query_param:

  A list of parameters for customizing the query group plot. Available
  parameters: palette (color palette for groups) and cells.highlight
  (whether to highlight cells).

- ref_param:

  A list of parameters for customizing the reference group plot.
  Available parameters: palette (color palette for groups) and
  cells.highlight (whether to highlight cells).

- xlim:

  The x-axis limits for the plot. If not provided, the limits will be
  calculated based on the data.

- ylim:

  The y-axis limits for the plot. If not provided, the limits will be
  calculated based on the data.

- pt.size:

  The size of the points in the plot.

- stroke.highlight:

  The size of the stroke highlight for cells.

## Examples

``` r
data(panc8_sub)
panc8_sub <- standard_scop(panc8_sub)
#> ℹ [2025-11-19 14:31:23] Start standard scop workflow...
#> ℹ [2025-11-19 14:31:24] Checking a list of <Seurat> object...
#> ! [2025-11-19 14:31:24] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2025-11-19 14:31:24] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-19 14:31:26] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-19 14:31:26] Use the separate HVF from srt_list
#> ℹ [2025-11-19 14:31:27] Number of available HVF: 2000
#> ℹ [2025-11-19 14:31:27] Finished check
#> ℹ [2025-11-19 14:31:27] Perform `Seurat::ScaleData()`
#> ℹ [2025-11-19 14:31:27] Perform pca linear dimension reduction
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
#> ℹ [2025-11-19 14:31:29] Perform `Seurat::FindClusters()` with louvain and `cluster_resolution` = 0.6
#> ℹ [2025-11-19 14:31:29] Reorder clusters...
#> ℹ [2025-11-19 14:31:29] Perform umap nonlinear dimension reduction
#> ℹ [2025-11-19 14:31:29] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-19 14:31:29] UMAP will return its model
#> ℹ [2025-11-19 14:31:33] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-19 14:31:33] UMAP will return its model
#> ✔ [2025-11-19 14:31:37] Run scop standard workflow done
srt_ref <- panc8_sub[, panc8_sub$tech != "fluidigmc1"]
srt_query <- panc8_sub[, panc8_sub$tech == "fluidigmc1"]
srt_ref <- integration_scop(
  srt_ref,
  batch = "tech",
  integration_method = "Uncorrected"
)
#> ◌ [2025-11-19 14:31:38] Run Uncorrected integration...
#> ℹ [2025-11-19 14:31:38] Spliting `srt_merge` into `srt_list` by column "tech"...
#> ℹ [2025-11-19 14:31:39] Checking a list of <Seurat> object...
#> ℹ [2025-11-19 14:31:39] Data 1/4 of the `srt_list` has been log-normalized
#> ℹ [2025-11-19 14:31:39] Perform `Seurat::FindVariableFeatures()` on the data 1/4 of the `srt_list`...
#> ℹ [2025-11-19 14:31:39] Data 2/4 of the `srt_list` has been log-normalized
#> ℹ [2025-11-19 14:31:39] Perform `Seurat::FindVariableFeatures()` on the data 2/4 of the `srt_list`...
#> ℹ [2025-11-19 14:31:40] Data 3/4 of the `srt_list` has been log-normalized
#> ℹ [2025-11-19 14:31:40] Perform `Seurat::FindVariableFeatures()` on the data 3/4 of the `srt_list`...
#> ℹ [2025-11-19 14:31:40] Data 4/4 of the `srt_list` has been log-normalized
#> ℹ [2025-11-19 14:31:40] Perform `Seurat::FindVariableFeatures()` on the data 4/4 of the `srt_list`...
#> ℹ [2025-11-19 14:31:41] Use the separate HVF from srt_list
#> ℹ [2025-11-19 14:31:41] Number of available HVF: 2000
#> ℹ [2025-11-19 14:31:42] Finished check
#> ℹ [2025-11-19 14:31:43] Perform Uncorrected integration
#> ℹ [2025-11-19 14:31:44] Perform `Seurat::ScaleData()`
#> ℹ [2025-11-19 14:31:44] Perform linear dimension reduction("pca")
#> ℹ [2025-11-19 14:31:45] Perform Seurat::FindClusters ("louvain")
#> ℹ [2025-11-19 14:31:46] Reorder clusters...
#> ℹ [2025-11-19 14:31:46] Perform nonlinear dimension reduction ("umap")
#> ℹ [2025-11-19 14:31:46] Non-linear dimensionality reduction (umap) using (Uncorrectedpca) dims (1-10) as input
#> ℹ [2025-11-19 14:31:50] Non-linear dimensionality reduction (umap) using (Uncorrectedpca) dims (1-10) as input
#> ✔ [2025-11-19 14:31:55] Run Uncorrected integration done
CellDimPlot(
  srt_ref,
  group.by = c("celltype", "tech")
)


# Projection
srt_query <- RunKNNMap(
  srt_query = srt_query,
  srt_ref = srt_ref,
  ref_umap = "UncorrectedUMAP2D"
)
#> ℹ [2025-11-19 14:31:56] Use the features to calculate distance metric
#> ℹ [2025-11-19 14:31:56] Data type is log-normalized
#> ℹ [2025-11-19 14:31:57] Data type is log-normalized
#> ℹ [2025-11-19 14:31:57] Use 2000 features to calculate distance
#> ℹ [2025-11-19 14:31:58] Use raw method to find neighbors
#> ℹ [2025-11-19 14:31:58] Running UMAP projection
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
