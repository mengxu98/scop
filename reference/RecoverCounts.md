# Attempt to recover raw counts from the normalized matrix

Attempt to recover raw counts from the normalized matrix

## Usage

``` r
RecoverCounts(
  srt,
  assay = NULL,
  trans = c("expm1", "exp", "none"),
  min_count = c(1, 2, 3),
  tolerance = 0.1,
  sf = NULL,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  Name of assay to recover counts.

- trans:

  The transformation function to applied when data is presumed to be
  log-normalized.

- min_count:

  Minimum UMI count of genes.

- tolerance:

  When recovering the raw counts, the nCount of each cell is
  theoretically calculated as an integer. However, due to decimal point
  preservation during normalization, the calculated nCount is usually a
  floating point number close to the integer. The tolerance is its
  difference from the integer. Default is `0.1`

- sf:

  Set the scaling factor manually.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
#> ℹ [2025-11-13 12:03:05] Start standard scop workflow...
#> ℹ [2025-11-13 12:03:06] Checking a list of <Seurat> object...
#> ! [2025-11-13 12:03:06] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2025-11-13 12:03:06] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-13 12:03:08] Perform `Seurat::FindVariableFeatures()` on the data 1/1 of the `srt_list`...
#> ℹ [2025-11-13 12:03:09] Use the separate HVF from srt_list
#> ℹ [2025-11-13 12:03:09] Number of available HVF: 2000
#> ℹ [2025-11-13 12:03:09] Finished check
#> ℹ [2025-11-13 12:03:09] Perform `Seurat::ScaleData()`
#> ℹ [2025-11-13 12:03:10] Perform pca linear dimension reduction
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
#> ℹ [2025-11-13 12:03:10] Perform `Seurat::FindClusters()` with louvain and `cluster_resolution` = 0.6
#> ℹ [2025-11-13 12:03:11] Reorder clusters...
#> ℹ [2025-11-13 12:03:11] Perform umap nonlinear dimension reduction
#> ℹ [2025-11-13 12:03:11] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-13 12:03:11] UMAP will return its model
#> ℹ [2025-11-13 12:03:14] Non-linear dimensionality reduction (umap) using (Standardpca) dims (1-50) as input
#> ℹ [2025-11-13 12:03:15] UMAP will return its model
#> ✔ [2025-11-13 12:03:18] Run scop standard workflow done
raw_counts <- GetAssayData5(
  pancreas_sub,
  assay = "RNA",
  layer = "counts"
)

# Normalized the data
pancreas_sub <- Seurat::NormalizeData(pancreas_sub)
#> Normalizing layer: counts

# Now replace counts with the log-normalized data matrix
data <- GetAssayData5(
  pancreas_sub,
  assay = "RNA",
  layer = "data"
)
new_pancreas_sub <- SeuratObject::SetAssayData(
  object = pancreas_sub,
  layer = "counts",
  new.data = data,
  assay = "RNA"
)
# Recover the counts and compare with the raw counts matrix
pancreas_sub <- RecoverCounts(new_pancreas_sub)
#> ℹ [2025-11-13 12:03:20] Data type is log-normalized
#> ℹ [2025-11-13 12:03:20] The data is presumed to be log-normalized
#> ℹ [2025-11-13 12:03:20] Perform "expm1" on the raw data
new_counts <- GetAssayData5(
  pancreas_sub,
  assay = "RNA",
  layer = "counts"
)
identical(raw_counts, new_counts)
#> [1] TRUE
```
