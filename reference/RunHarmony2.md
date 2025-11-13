# Run Harmony algorithm

This is a modified version of
[harmony::RunHarmony](https://rdrr.io/pkg/harmony/man/RunHarmony.html)
specifically designed for compatibility with
[RunSymphonyMap](https://mengxu98.github.io/scop/reference/RunSymphonyMap.md).

## Usage

``` r
RunHarmony2(object, ...)

# S3 method for class 'Seurat'
RunHarmony2(
  object,
  group.by.vars,
  assay = NULL,
  reduction = "pca",
  dims.use = 1:30,
  project.dim = TRUE,
  reduction.name = "Harmony",
  reduction.key = "Harmony_",
  verbose = TRUE,
  seed.use = 11L,
  ...
)
```

## Arguments

- object:

  A Seurat object.

- ...:

  Additional arguments to be passed to
  [harmony::RunHarmony](https://rdrr.io/pkg/harmony/man/RunHarmony.html).

- group.by.vars:

  The batch variable name.

- assay:

  The assay to be used. Default is `NULL`.

- reduction:

  The reduction to be used. Default is `"pca"`.

- dims.use:

  The dimensions to be used. Default is `1:30`.

- project.dim:

  Whether to project dimension reduction loadings. Default is `TRUE`.

- reduction.name:

  The name of the reduction to be stored in the Seurat object. Default
  is `"Harmony"`.

- reduction.key:

  The prefix for the column names of the Harmony embeddings. Default is
  `"Harmony_"`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed.use:

  The random seed to be used. Default is `11`.

## Examples

``` r
data(panc8_sub)
panc8_sub <- standard_scop(panc8_sub)
#> ℹ [2025-11-13 12:23:35] Start standard scop workflow...
#> ℹ [2025-11-13 12:23:35] Checking a list of <Seurat> object...
#> ! [2025-11-13 12:23:35] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2025-11-13 12:23:35] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-13 12:23:38] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-13 12:23:38] Use the separate HVF from srt_list
#> ℹ [2025-11-13 12:23:38] Number of available HVF: 2000
#> ℹ [2025-11-13 12:23:39] Finished check
#> ℹ [2025-11-13 12:23:39] Perform `Seurat::ScaleData()`
#> ℹ [2025-11-13 12:23:39] Perform pca linear dimension reduction
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
#> ℹ [2025-11-13 12:23:41] Perform `Seurat::FindClusters()` with louvain and `cluster_resolution` = 0.6
#> ℹ [2025-11-13 12:23:41] Reorder clusters...
#> ℹ [2025-11-13 12:23:41] Perform umap nonlinear dimension reduction
#> ℹ [2025-11-13 12:23:41] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-13 12:23:41] UMAP will return its model
#> ℹ [2025-11-13 12:23:46] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-13 12:23:46] UMAP will return its model
#> ✔ [2025-11-13 12:23:51] Run scop standard workflow done
panc8_sub <- RunHarmony2(
  panc8_sub,
  group.by.vars = "tech",
  reduction = "pca"
)
#> ◌ [2025-11-13 12:23:51] Installing: harmony...
#>  
#> → Will install 2 packages.
#> → All 2 packages (0 B) are cached.
#> + RhpcBLASctl   0.23-42 
#> + harmony       1.2.4   
#>   
#> ℹ No downloads are needed, 2 pkgs are cached
#> ✔ Installed RhpcBLASctl 0.23-42  (1s)
#> ✔ Installed harmony 1.2.4  (1.1s)
#> ✔ 1 pkg + 29 deps: kept 28, added 2 [2.6s]
#> ✔ [2025-11-13 12:23:54] harmony installed successfully
#> Transposing data matrix
#> Initializing state using k-means centroids initialization
#> Harmony 1/10
#> Harmony 2/10
#> Harmony 3/10
#> Harmony 4/10
#> Harmony 5/10
#> Harmony 6/10
#> Harmony 7/10
#> Harmony converged after 7 iterations

CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype"),
  reduction = "pca"
)


CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype"),
  reduction = "Harmony"
)


panc8_sub <- standard_scop(
  panc8_sub,
  prefix = "Harmony",
  linear_reduction = "Harmony"
)
#> ℹ [2025-11-13 12:23:55] Start standard scop workflow...
#> ℹ [2025-11-13 12:23:55] Checking a list of <Seurat> object...
#> ℹ [2025-11-13 12:23:56] Data 1/1 of the `srt_list` has been log-normalized
#> ℹ [2025-11-13 12:23:56] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-13 12:23:57] Use the separate HVF from srt_list
#> ℹ [2025-11-13 12:23:57] Number of available HVF: 2000
#> ℹ [2025-11-13 12:23:57] Finished check
#> ℹ [2025-11-13 12:23:57] Perform `Seurat::ScaleData()`
#> ℹ [2025-11-13 12:23:58] Perform Harmony linear dimension reduction
#> ℹ [2025-11-13 12:23:58] `linear_reduction` Harmony is already existed. Skip calculation
#> ℹ [2025-11-13 12:23:58] Perform `Seurat::FindClusters()` with louvain and `cluster_resolution` = 0.6
#> ℹ [2025-11-13 12:23:58] Reorder clusters...
#> ℹ [2025-11-13 12:23:58] Perform umap nonlinear dimension reduction
#> ℹ [2025-11-13 12:23:58] Non-linear dimensionality reduction (umap) using (HarmonyHarmony) dims (1-30) as input
#> ℹ [2025-11-13 12:23:58] UMAP will return its model
#> ℹ [2025-11-13 12:24:04] Non-linear dimensionality reduction (umap) using (HarmonyHarmony) dims (1-30) as input
#> ℹ [2025-11-13 12:24:04] UMAP will return its model
#> ✔ [2025-11-13 12:24:09] Run scop standard workflow done

CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype"),
  reduction = "StandardpcaUMAP2D"
)


CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype"),
  reduction = "HarmonyUMAP2D"
)
```
