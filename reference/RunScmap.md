# Annotate single cells using scmap.

Annotate single cells using scmap.

## Usage

``` r
RunScmap(
  srt_query,
  srt_ref,
  ref_group = NULL,
  query_assay = "RNA",
  ref_assay = "RNA",
  method = "scmapCluster",
  nfeatures = 500,
  threshold = 0.5,
  k = 10
)
```

## Arguments

- srt_query:

  An object of class Seurat to be annotated with cell types.

- srt_ref:

  An object of class Seurat storing the reference cells.

- ref_group:

  A character vector specifying the column name in the \`srt_ref\`
  metadata that represents the cell grouping.

- query_assay:

  A character vector specifying the assay to be used for the query data.
  Defaults to the default assay of the \`srt_query\` object.

- ref_assay:

  A character vector specifying the assay to be used for the reference
  data. Defaults to the default assay of the \`srt_ref\` object.

- method:

  The method to be used for scmap analysis. Can be any of "scmapCluster"
  or "scmapCell". The default value is "scmapCluster".

- nfeatures:

  The number of top features to be selected. The default value is 500.

- threshold:

  The threshold value on similarity to determine if a cell is assigned
  to a cluster. This should be a value between 0 and 1. The default
  value is 0.5.

- k:

  Number of clusters per group for k-means clustering when method is
  "scmapCell".

## Examples

``` r
data(panc8_sub)
panc8_sub <- standard_scop(panc8_sub)
#> ℹ [2025-11-19 15:00:31] Start standard scop workflow...
#> ℹ [2025-11-19 15:00:32] Checking a list of <Seurat> object...
#> ! [2025-11-19 15:00:32] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2025-11-19 15:00:32] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-19 15:00:34] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-19 15:00:35] Use the separate HVF from srt_list
#> ℹ [2025-11-19 15:00:35] Number of available HVF: 2000
#> ℹ [2025-11-19 15:00:35] Finished check
#> ℹ [2025-11-19 15:00:35] Perform `Seurat::ScaleData()`
#> ℹ [2025-11-19 15:00:36] Perform pca linear dimension reduction
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
#> ℹ [2025-11-19 15:00:37] Perform `Seurat::FindClusters()` with louvain and `cluster_resolution` = 0.6
#> ℹ [2025-11-19 15:00:37] Reorder clusters...
#> ℹ [2025-11-19 15:00:37] Perform umap nonlinear dimension reduction
#> ℹ [2025-11-19 15:00:37] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-19 15:00:37] UMAP will return its model
#> ℹ [2025-11-19 15:00:42] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-19 15:00:42] UMAP will return its model
#> ✔ [2025-11-19 15:00:47] Run scop standard workflow done

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
#> ℹ [2025-11-19 15:00:47] Rename features for the assay: RNA
panc8_sub <- CheckDataMerge(
  panc8_sub,
  batch = "tech"
)[["srt_merge"]]
#> ℹ [2025-11-19 15:00:48] Spliting `srt_merge` into `srt_list` by column "tech"...
#> ℹ [2025-11-19 15:00:49] Checking a list of <Seurat> object...
#> ℹ [2025-11-19 15:00:49] Data 1/5 of the `srt_list` has been log-normalized
#> ℹ [2025-11-19 15:00:49] Perform `Seurat::FindVariableFeatures()` on the data 1/5 of the `srt_list`...
#> ℹ [2025-11-19 15:00:49] Data 2/5 of the `srt_list` has been log-normalized
#> ℹ [2025-11-19 15:00:49] Perform `Seurat::FindVariableFeatures()` on the data 2/5 of the `srt_list`...
#> ℹ [2025-11-19 15:00:50] Data 3/5 of the `srt_list` has been log-normalized
#> ℹ [2025-11-19 15:00:50] Perform `Seurat::FindVariableFeatures()` on the data 3/5 of the `srt_list`...
#> ℹ [2025-11-19 15:00:50] Data 4/5 of the `srt_list` has been log-normalized
#> ℹ [2025-11-19 15:00:50] Perform `Seurat::FindVariableFeatures()` on the data 4/5 of the `srt_list`...
#> ℹ [2025-11-19 15:00:51] Data 5/5 of the `srt_list` has been log-normalized
#> ℹ [2025-11-19 15:00:51] Perform `Seurat::FindVariableFeatures()` on the data 5/5 of the `srt_list`...
#> ℹ [2025-11-19 15:00:51] Use the separate HVF from srt_list
#> ℹ [2025-11-19 15:00:52] Number of available HVF: 2000
#> ℹ [2025-11-19 15:00:52] Finished check

data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2025-11-19 15:00:54] Start standard scop workflow...
#> ℹ [2025-11-19 15:00:55] Checking a list of <Seurat> object...
#> ! [2025-11-19 15:00:55] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2025-11-19 15:00:55] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-19 15:00:57] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-19 15:00:58] Use the separate HVF from srt_list
#> ℹ [2025-11-19 15:00:58] Number of available HVF: 2000
#> ℹ [2025-11-19 15:00:58] Finished check
#> ℹ [2025-11-19 15:00:58] Perform `Seurat::ScaleData()`
#> ℹ [2025-11-19 15:00:59] Perform pca linear dimension reduction
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
#> ℹ [2025-11-19 15:01:00] Perform `Seurat::FindClusters()` with louvain and `cluster_resolution` = 0.6
#> ℹ [2025-11-19 15:01:00] Reorder clusters...
#> ℹ [2025-11-19 15:01:00] Perform umap nonlinear dimension reduction
#> ℹ [2025-11-19 15:01:00] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-19 15:01:00] UMAP will return its model
#> ℹ [2025-11-19 15:01:05] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-19 15:01:05] UMAP will return its model
#> ✔ [2025-11-19 15:01:09] Run scop standard workflow done
pancreas_sub <- RunScmap(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  ref_group = "celltype",
  method = "scmapCluster"
)
#> ◌ [2025-11-19 15:01:09] Installing: scmap...
#>  
#> → Will install 3 packages.
#> → All 3 packages (0 B) are cached.
#> + googleVis      0.7.3   
#> + randomForest   4.7-1.2 
#> + scmap          1.32.0  [bld][cmp]
#> ✔ All system requirements are already installed.
#>   
#> ℹ No downloads are needed, 3 pkgs are cached
#> ✔ Got randomForest 4.7-1.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (218.82 kB)
#> ✔ Got scmap 1.32.0 (source) (2.34 MB)
#> ✔ Got googleVis 0.7.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (498.39 kB)
#> ℹ Installing system requirements
#> ℹ Executing `sudo sh -c apt-get -y update`
#> Get:1 file:/etc/apt/apt-mirrors.txt Mirrorlist [144 B]
#> Hit:6 https://packages.microsoft.com/repos/azure-cli noble InRelease
#> Hit:2 http://azure.archive.ubuntu.com/ubuntu noble InRelease
#> Hit:7 https://packages.microsoft.com/ubuntu/24.04/prod noble InRelease
#> Hit:3 http://azure.archive.ubuntu.com/ubuntu noble-updates InRelease
#> Hit:4 http://azure.archive.ubuntu.com/ubuntu noble-backports InRelease
#> Hit:5 http://azure.archive.ubuntu.com/ubuntu noble-security InRelease
#> Reading package lists...
#> ℹ Executing `sudo sh -c apt-get -y install libicu-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> libicu-dev is already the newest version (74.2-1ubuntu3.1).
#> 0 upgraded, 0 newly installed, 0 to remove and 31 not upgraded.
#> ✔ Installed googleVis 0.7.3  (33ms)
#> ✔ Installed randomForest 4.7-1.2  (48ms)
#> ℹ Building scmap 1.32.0
#> ✔ Built scmap 1.32.0 (22.2s)
#> ✔ Installed scmap 1.32.0  (1.1s)
#> ✔ 1 pkg + 54 deps: kept 52, added 3, dld 3 (3.06 MB) [27.8s]
#> ✔ [2025-11-19 15:01:37] scmap installed successfully
#> ℹ [2025-11-19 15:01:37] Data type is log-normalized
#> ℹ [2025-11-19 15:01:37] Detected srt_query data type: log_normalized_counts
#> ℹ [2025-11-19 15:01:39] Data type is log-normalized
#> ℹ [2025-11-19 15:01:39] Detected srt_ref data type: log_normalized_counts
#> ℹ [2025-11-19 15:01:41] Perform selectFeatures
#> ℹ [2025-11-19 15:01:41] Perform indexCluster
#> ℹ [2025-11-19 15:01:42] Perform scmapCluster
CellDimPlot(
  pancreas_sub,
  group.by = "scmap_annotation"
)


pancreas_sub <- RunScmap(
  srt_query = pancreas_sub,
  srt_ref = panc8_sub,
  ref_group = "celltype",
  method = "scmapCell"
)
#> ✔ [2025-11-19 15:01:42] scmap installed successfully
#> ℹ [2025-11-19 15:01:43] Data type is log-normalized
#> ℹ [2025-11-19 15:01:43] Detected srt_query data type: log_normalized_counts
#> ℹ [2025-11-19 15:01:44] Data type is log-normalized
#> ℹ [2025-11-19 15:01:44] Detected srt_ref data type: log_normalized_counts
#> ℹ [2025-11-19 15:01:46] Perform selectFeatures
#> ℹ [2025-11-19 15:01:47] Perform indexCell
#> ℹ [2025-11-19 15:01:47] Perform scmapCell
#> ℹ [2025-11-19 15:01:48] Perform scmapCell2Cluster
CellDimPlot(
  pancreas_sub,
  group.by = "scmap_annotation"
)
```
