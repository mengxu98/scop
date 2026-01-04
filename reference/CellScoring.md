# Cell scoring

This function performs cell scoring on a Seurat object. It calculates
scores for a given set of features and adds the scores as metadata to
the Seurat object.

## Usage

``` r
CellScoring(
  srt,
  features = NULL,
  layer = "data",
  assay = NULL,
  split.by = NULL,
  IDtype = "symbol",
  species = "Homo_sapiens",
  db = "GO_BP",
  termnames = NULL,
  db_update = FALSE,
  db_version = "latest",
  convert_species = TRUE,
  Ensembl_version = NULL,
  mirror = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  method = "Seurat",
  classification = TRUE,
  name = "",
  new_assay = FALSE,
  seed = 11,
  cores = 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- features:

  A named list of feature lists for scoring. If NULL, `db` will be used
  to create features sets.

- layer:

  The layer of the Seurat object to use for scoring. Default is
  `"data"`.

- assay:

  The assay of the Seurat object to use for scoring. Default is `NULL`,
  in which case the default assay of the object is used.

- split.by:

  A cell metadata variable used for splitting the Seurat object into
  subsets and performing scoring on each subset. Default is `NULL`.

- IDtype:

  A character vector specifying the type of gene IDs in the `srt` object
  or `geneID` argument. This argument is used to convert the gene IDs to
  a different type if `IDtype` is different from `result_IDtype`.

- species:

  A character vector specifying the species for which the analysis is
  performed.

- db:

  A character vector specifying the name of the database to be used for
  enrichment analysis.

- termnames:

  A vector of term names to be used from the database. Default is
  `NULL`, in which case all features from the database are used.

- db_update:

  Whether the gene annotation databases should be forcefully updated. If
  set to FALSE, the function will attempt to load the cached databases
  instead. Default is `FALSE`.

- db_version:

  A character vector specifying the version of the database to be used.
  This argument is ignored if `db_update` is `TRUE`. Default is
  `"latest"`.

- convert_species:

  Whether to use a species-converted database when the annotation is
  missing for the specified species. Default is `TRUE`.

- Ensembl_version:

  Ensembl database version. If NULL, use the current release version.

- mirror:

  Specify an Ensembl mirror to connect to. The valid options here are
  `"www"`, `"uswest"`, `"useast"`, `"asia"`.

- minGSSize:

  The minimum size of a gene set to be considered in the enrichment
  analysis.

- maxGSSize:

  The maximum size of a gene set to be considered in the enrichment
  analysis.

- method:

  The method to use for scoring. Can be "Seurat", "AUCell", or "UCell".
  Default is `"Seurat"`.

- classification:

  Whether to perform classification based on the scores. Default is
  `TRUE`.

- name:

  The name of the assay to store the scores in. Only used if new_assay
  is TRUE. Default is `""`.

- new_assay:

  Whether to create a new assay for storing the scores. Default is
  `FALSE`.

- seed:

  The random seed for reproducibility. Default is `11`.

- cores:

  The number of cores to use for parallelization with
  [foreach::foreach](https://rdrr.io/pkg/foreach/man/foreach.html).
  Default is `1`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments to be passed to the scoring methods.

## See also

[PrepareDB](https://mengxu98.github.io/scop/reference/PrepareDB.md),
[ListDB](https://mengxu98.github.io/scop/reference/ListDB.md),
[RunDynamicFeatures](https://mengxu98.github.io/scop/reference/RunDynamicFeatures.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
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
features_all <- rownames(pancreas_sub)
pancreas_sub <- CellScoring(
  pancreas_sub,
  features = list(
    A = features_all[1:100],
    B = features_all[101:200]
  ),
  method = "Seurat",
  name = "test"
)
#> ⠙ [2026-01-04 08:37:06] Running [1/2] Processing: 1  ETA:  0s
#> ✔ [2026-01-04 08:37:06] Completed 2 tasks in 132ms
#> 
CellDimPlot(pancreas_sub, "test_classification")
#> Warning: No shared levels found between `names(values)` of the manual scale and the
#> data's fill values.


FeatureDimPlot(pancreas_sub, "test_A")


if (FALSE) { # \dontrun{
data(panc8_sub)
panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "Seurat"
)
CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype")
)

panc8_sub <- CellScoring(
  panc8_sub,
  layer = "data",
  assay = "RNA",
  db = "GO_BP",
  species = "Homo_sapiens",
  minGSSize = 10,
  maxGSSize = 100,
  method = "Seurat",
  name = "GO",
  new_assay = TRUE
)

panc8_sub <- integration_scop(
  panc8_sub,
  assay = "GO",
  batch = "tech",
  integration_method = "Seurat"
)
CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype")
)

pancreas_sub <- CellScoring(
  pancreas_sub,
  layer = "data",
  assay = "RNA",
  db = "GO_BP",
  species = "Mus_musculus",
  termnames = panc8_sub[["GO"]]@meta.features[, "termnames"],
  method = "Seurat",
  name = "GO",
  new_assay = TRUE
)
pancreas_sub <- standard_scop(
  pancreas_sub,
  assay = "GO"
)
CellDimPlot(pancreas_sub, "SubCellType")

pancreas_sub[["tech"]] <- "Mouse"
panc_merge <- integration_scop(
  srt_list = list(panc8_sub, pancreas_sub),
  assay = "GO",
  batch = "tech", integration_method = "Seurat"
)
CellDimPlot(
  srt = panc_merge,
  group.by = c("tech", "celltype", "SubCellType", "Phase")
)

genenames <- make.unique(
  thisutils::capitalize(
    rownames(panc8_sub[["RNA"]])
  ),
  force_tolower = TRUE
)
names(genenames) <- rownames(panc8_sub)
panc8_sub <- RenameFeatures(
  panc8_sub,
  newnames = genenames,
  assay = "RNA"
)
head(rownames(panc8_sub))
panc_merge <- integration_scop(
  srt_list = list(panc8_sub, pancreas_sub),
  assay = "RNA",
  batch = "tech", integration_method = "Seurat"
)
CellDimPlot(
  srt = panc_merge,
  group.by = c("tech", "celltype", "SubCellType", "Phase")
)
} # }
```
