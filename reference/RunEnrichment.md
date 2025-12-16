# Perform the enrichment analysis (over-representation) on the genes

Perform the enrichment analysis (over-representation) on the genes

## Usage

``` r
RunEnrichment(
  srt = NULL,
  group_by = NULL,
  test.use = "wilcox",
  DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05",
  geneID = NULL,
  geneID_groups = NULL,
  geneID_exclude = NULL,
  IDtype = "symbol",
  result_IDtype = "symbol",
  species = "Homo_sapiens",
  db = "GO_BP",
  db_update = FALSE,
  db_version = "latest",
  db_combine = FALSE,
  convert_species = TRUE,
  Ensembl_version = NULL,
  mirror = NULL,
  TERM2GENE = NULL,
  TERM2NAME = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  unlimited_db = c("Chromosome", "GeneType", "TF", "Enzyme", "CSPA"),
  GO_simplify = FALSE,
  GO_simplify_cutoff = "p.adjust < 0.05",
  simplify_method = "Wang",
  simplify_similarityCutoff = 0.7,
  cores = 1,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object containing the results of differential expression
  analysis (RunDEtest). If specified, the genes and groups will be
  extracted from the Seurat object automatically. If not specified, the
  `geneID` and `geneID_groups` arguments must be provided.

- group_by:

  A character vector specifying the grouping variable in the Seurat
  object. This argument is only used if `srt` is specified.

- test.use:

  A character vector specifying the test to be used in differential
  expression analysis. This argument is only used if `srt` is specified.

- DE_threshold:

  A character vector specifying the filter condition for differential
  expression analysis. This argument is only used if `srt` is specified.

- geneID:

  A character vector specifying the gene IDs.

- geneID_groups:

  A factor vector specifying the group labels for each gene.

- geneID_exclude:

  A character vector specifying the gene IDs to be excluded from the
  analysis.

- IDtype:

  A character vector specifying the type of gene IDs in the `srt` object
  or `geneID` argument. This argument is used to convert the gene IDs to
  a different type if `IDtype` is different from `result_IDtype`.

- result_IDtype:

  A character vector specifying the desired type of gene ID to be used
  in the output. This argument is used to convert the gene IDs from
  `IDtype` to `result_IDtype`.

- species:

  A character vector specifying the species for which the analysis is
  performed.

- db:

  A character vector specifying the name of the database to be used for
  enrichment analysis.

- db_update:

  Whether the gene annotation databases should be forcefully updated. If
  set to FALSE, the function will attempt to load the cached databases
  instead. Default is `FALSE`.

- db_version:

  A character vector specifying the version of the database to be used.
  This argument is ignored if `db_update` is `TRUE`. Default is
  `"latest"`.

- db_combine:

  Whether to combine multiple databases into one. If TRUE, all database
  specified by `db` will be combined as one named "Combined".

- convert_species:

  Whether to use a species-converted database when the annotation is
  missing for the specified species. Default is `TRUE`.

- Ensembl_version:

  Ensembl database version. If NULL, use the current release version.

- mirror:

  Specify an Ensembl mirror to connect to. The valid options here are
  `"www"`, `"uswest"`, `"useast"`, `"asia"`.

- TERM2GENE:

  A data frame specifying the gene-term mapping for a custom database.
  The first column should contain the term IDs, and the second column
  should contain the gene IDs.

- TERM2NAME:

  A data frame specifying the term-name mapping for a custom database.
  The first column should contain the term IDs, and the second column
  should contain the corresponding term names.

- minGSSize:

  The minimum size of a gene set to be considered in the enrichment
  analysis.

- maxGSSize:

  The maximum size of a gene set to be considered in the enrichment
  analysis.

- unlimited_db:

  A character vector specifying the names of databases that do not have
  size restrictions.

- GO_simplify:

  Whether to simplify the GO terms. If `TRUE`, additional results with
  simplified GO terms will be returned.

- GO_simplify_cutoff:

  A character vector specifying the filter condition for simplification
  of GO terms. This argument is only used if `GO_simplify` is `TRUE`.

- simplify_method:

  A character vector specifying the method to be used for simplification
  of GO terms. This argument is only used if `GO_simplify` is `TRUE`.

- simplify_similarityCutoff:

  The similarity cutoff for simplification of GO terms. This argument is
  only used if `GO_simplify` is `TRUE`.

- cores:

  The number of cores to use for parallelization with
  [foreach::foreach](https://rdrr.io/pkg/foreach/man/foreach.html).
  Default is `1`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

If input is a Seurat object, returns the modified Seurat object with the
enrichment result stored in the tools slot.

If input is a geneID vector with or without geneID_groups, return the
enrichment result directly.

Enrichment result is a list with the following component:

- `enrichment`: A data.frame containing all enrichment results.

- `results`: A list of `enrichResult` objects from the DOSE package.

- `geneMap`: A data.frame containing the ID mapping table for input gene
  IDs.

- `input`: A data.frame containing the input gene IDs and gene ID
  groups.

- `DE_threshold`: A specific threshold for differential expression
  analysis (only returned if input is a Seurat object).

## See also

[PrepareDB](https://mengxu98.github.io/scop/reference/PrepareDB.md),
[ListDB](https://mengxu98.github.io/scop/reference/ListDB.md),
[EnrichmentPlot](https://mengxu98.github.io/scop/reference/EnrichmentPlot.md),
[RunGSEA](https://mengxu98.github.io/scop/reference/RunGSEA.md),
[GSEAPlot](https://mengxu98.github.io/scop/reference/GSEAPlot.md)

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
pancreas_sub <- RunDEtest(
  pancreas_sub,
  group_by = "CellType"
)
#> Error: The `slot` argument of `Assays()` was deprecated in SeuratObject 5.0.0
#> and is now defunct.
#> ℹ Please use `LayerData()` instead.
pancreas_sub <- RunEnrichment(
  pancreas_sub,
  group_by = "CellType",
  DE_threshold = "p_val_adj < 0.05",
  db = "GO_BP",
  species = "Mus_musculus"
)
#>  
#> → Will install 29 packages.
#> → Will download 3 CRAN packages (1.93 MB), cached: 26 (0 B).
#> + DOSE                4.4.0  [bld]
#> + GO.db               3.22.0 [bld]
#> + GOSemSim            2.36.0 [bld][cmp]
#> + R.methodsS3         1.8.2  
#> + R.oo                1.27.1 
#> + R.utils             2.13.0 
#> + ape                 5.8-1  [bld][cmp][dl] (1.42 MB)
#> + clusterProfiler     4.18.3 [bld]
#> + enrichplot          1.30.4 [bld]
#> + fastmatch           1.1-6  
#> + fgsea               1.36.0 [bld][cmp]
#> + fontBitstreamVera   0.1.1  
#> + fontLiberation      0.1.0  
#> + fontquiver          0.2.1  
#> + gdtools             0.4.4  [bld][cmp][dl] (73.90 kB) + ✔ libcairo2-dev, ✔ libfontconfig1-dev, ✔ libfreetype6-dev
#> + ggforce             0.5.0  
#> + ggiraph             0.9.2  [bld][cmp][dl] (429.12 kB) + ✔ libpng-dev
#> + ggnewscale          0.5.2  
#> + ggtangle            0.0.9  
#> + ggtree              4.0.1  [bld]
#> + gson                0.1.0  
#> + org.Hs.eg.db        3.22.0 [bld]
#> + quarto              1.5.1  
#> + qvalue              2.42.0 [bld]
#> + scatterpie          0.2.6  
#> + tidydr              0.0.6  
#> + tidytree            0.4.6  
#> + treeio              1.34.0 [bld]
#> + tweenr              2.0.3  
#> ✔ All system requirements are already installed.
#>   
#> ℹ Getting 3 pkgs (1.93 MB), 26 cached
#> ✔ Got ggtree 4.0.1 (source) (370.24 kB)
#> ✔ Got GOSemSim 2.36.0 (source) (610.99 kB)
#> ✔ Got treeio 1.34.0 (source) (701.64 kB)
#> ✔ Got fastmatch 1.1-6 (x86_64-pc-linux-gnu-ubuntu-24.04) (35.95 kB)
#> ✔ Got tweenr 2.0.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (461.68 kB)
#> ✔ Got ggnewscale 0.5.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (351.15 kB)
#> ✔ Got fontBitstreamVera 0.1.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (699.47 kB)
#> ✔ Got quarto 1.5.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (544.37 kB)
#> ✔ Got R.oo 1.27.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (996.09 kB)
#> ✔ Got ggforce 0.5.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.95 MB)
#> ✔ Got qvalue 2.42.0 (source) (2.77 MB)
#> ✔ Got fontquiver 0.2.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.28 MB)
#> ✔ Got scatterpie 0.2.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (150.90 kB)
#> ✔ Got ggtangle 0.0.9 (x86_64-pc-linux-gnu-ubuntu-24.04) (257.09 kB)
#> ✔ Got R.methodsS3 1.8.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (82.67 kB)
#> ✔ Got ggiraph 0.9.2 (source) (430.83 kB)
#> ✔ Got R.utils 2.13.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.45 MB)
#> ✔ Got gson 0.1.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (232.54 kB)
#> ✔ Got gdtools 0.4.4 (source) (73.99 kB)
#> ✔ Got tidytree 0.4.6 (x86_64-pc-linux-gnu-ubuntu-24.04) (343.02 kB)
#> ✔ Got ape 5.8-1 (source) (1.43 MB)
#> ✔ Got fontLiberation 0.1.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (4.54 MB)
#> ✔ Got DOSE 4.4.0 (source) (5.75 MB)
#> ✔ Got fgsea 1.36.0 (source) (6.17 MB)
#> ✔ Got GO.db 3.22.0 (source) (25.25 MB)
#> ✔ Got org.Hs.eg.db 3.22.0 (source) (104.97 MB)
#> ℹ Installing system requirements
#> ℹ Executing `sudo sh -c apt-get -y update`
#> Get:1 file:/etc/apt/apt-mirrors.txt Mirrorlist [144 B]
#> Hit:2 http://azure.archive.ubuntu.com/ubuntu noble InRelease
#> Hit:6 https://packages.microsoft.com/repos/azure-cli noble InRelease
#> Hit:3 http://azure.archive.ubuntu.com/ubuntu noble-updates InRelease
#> Hit:4 http://azure.archive.ubuntu.com/ubuntu noble-backports InRelease
#> Hit:5 http://azure.archive.ubuntu.com/ubuntu noble-security InRelease
#> Hit:7 https://packages.microsoft.com/ubuntu/24.04/prod noble InRelease
#> Reading package lists...
#> ℹ Executing `sudo sh -c apt-get -y install libcairo2-dev libfontconfig1-dev libfreetype6-dev libpng-dev libx11-dev libcurl4-openssl-dev libssl-dev make libglpk-dev libxml2-dev pandoc libicu-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> libcairo2-dev is already the newest version (1.18.0-3build1).
#> libfontconfig1-dev is already the newest version (2.15.0-1.1ubuntu2).
#> libfreetype-dev is already the newest version (2.13.2+dfsg-1build3).
#> libpng-dev is already the newest version (1.6.43-5ubuntu0.1).
#> libx11-dev is already the newest version (2:1.8.7-1build1).
#> libcurl4-openssl-dev is already the newest version (8.5.0-2ubuntu10.6).
#> libssl-dev is already the newest version (3.0.13-0ubuntu3.6).
#> make is already the newest version (4.3-4.1build2).
#> libglpk-dev is already the newest version (5.0-1build2).
#> libxml2-dev is already the newest version (2.9.14+dfsg-1.3ubuntu3.6).
#> pandoc is already the newest version (3.1.3+ds-2).
#> libicu-dev is already the newest version (74.2-1ubuntu3.1).
#> 0 upgraded, 0 newly installed, 0 to remove and 49 not upgraded.
#> ℹ Building ape 5.8-1
#> ℹ Building qvalue 2.42.0
#> ℹ Building GO.db 3.22.0
#> ℹ Building org.Hs.eg.db 3.22.0
#> ✔ Built qvalue 2.42.0 (4.3s)
#> ✔ Installed fastmatch 1.1-6  (25ms)
#> ℹ Building fgsea 1.36.0
#> ✔ Built ape 5.8-1 (34.7s)
#> ✔ Installed ape 5.8-1  (1.1s)
#> ✔ Installed fontBitstreamVera 0.1.1  (37ms)
#> ✔ Installed fontLiberation 0.1.0  (1.1s)
#> ✔ Installed fontquiver 0.2.1  (1.1s)
#> ℹ Building gdtools 0.4.4
#> ✔ Built fgsea 1.36.0 (36.3s)
#> ✔ Installed ggforce 0.5.0  (56ms)
#> ✔ Installed ggnewscale 0.5.2  (32ms)
#> ✔ Installed ggtangle 0.0.9  (32ms)
#> ✔ Installed gson 0.1.0  (80ms)
#> ✔ Installed quarto 1.5.1  (47ms)
#> ✔ Installed R.methodsS3 1.8.2  (28ms)
#> ✔ Installed R.oo 1.27.1  (45ms)
#> ✔ Installed R.utils 2.13.0  (47ms)
#> ✔ Installed scatterpie 0.2.6  (32ms)
#> ✔ Installed tidydr 0.0.6  (31ms)
#> ✔ Installed tidytree 0.4.6  (101ms)
#> ℹ Building treeio 1.34.0
#> ✔ Built treeio 1.34.0 (6s)
#> ✔ Installed tweenr 2.0.3  (36ms)
#> ✔ Installed fgsea 1.36.0  (169ms)
#> ✔ Installed qvalue 2.42.0  (50ms)
#> ✔ Installed treeio 1.34.0  (63ms)
#> ✔ Built GO.db 3.22.0 (48.7s)
#> ✔ Installed GO.db 3.22.0  (534ms)
#> ℹ Building GOSemSim 2.36.0
#> ✔ Built gdtools 0.4.4 (15.4s)
#> ✔ Installed gdtools 0.4.4  (40ms)
#> ℹ Building ggiraph 0.9.2
#> ✔ Built GOSemSim 2.36.0 (15.8s)
#> ✔ Installed GOSemSim 2.36.0  (53ms)
#> ℹ Building DOSE 4.4.0
#> ✔ Built DOSE 4.4.0 (14.1s)
#> ✔ Installed DOSE 4.4.0  (59ms)
#> ✔ Built ggiraph 0.9.2 (38.9s)
#> ✔ Installed ggiraph 0.9.2  (83ms)
#> ℹ Building ggtree 4.0.1
#> ✔ Built ggtree 4.0.1 (6.5s)
#> ✔ Installed ggtree 4.0.1  (1.1s)
#> ℹ Building enrichplot 1.30.4
#> ✔ Built enrichplot 1.30.4 (11.3s)
#> ✔ Installed enrichplot 1.30.4  (1s)
#> ℹ Building clusterProfiler 4.18.3
#> ✔ Built clusterProfiler 4.18.3 (12s)
#> ✔ Installed clusterProfiler 4.18.3  (31ms)
#> ✔ Built org.Hs.eg.db 3.22.0 (4m 38s)
#> ✔ Installed org.Hs.eg.db 3.22.0  (2.3s)
#> ✔ 1 pkg + 155 deps: kept 124, added 29, dld 26 (162.90 MB) [5m 1.4s]
#> Error in RunEnrichment(pancreas_sub, group_by = "CellType", DE_threshold = "p_val_adj < 0.05",     db = "GO_BP", species = "Mus_musculus"): Cannot find the DEtest result for the group 'CellType'. You may perform
#> RunDEtest first.
EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP",
  group_by = "CellType",
  plot_type = "comparison"
)
#> Error in EnrichmentPlot(pancreas_sub, db = "GO_BP", group_by = "CellType",     plot_type = "comparison"): No enrichment result found. You may perform RunEnrichment first

if (FALSE) { # \dontrun{
pancreas_sub <- RunEnrichment(
  pancreas_sub,
  group_by = "CellType",
  DE_threshold = "p_val_adj < 0.05",
  db = c("MSigDB", "MSigDB_MH"),
  species = "Mus_musculus"
)
EnrichmentPlot(
  pancreas_sub,
  db = "MSigDB",
  group_by = "CellType",
  plot_type = "comparison"
)
EnrichmentPlot(
  pancreas_sub,
  db = "MSigDB_MH",
  group_by = "CellType",
  plot_type = "comparison"
)

# Remove redundant GO terms
pancreas_sub <- RunEnrichment(
  pancreas_sub,
  group_by = "CellType",
  db = "GO_BP",
  GO_simplify = TRUE,
  species = "Mus_musculus"
)
EnrichmentPlot(
  pancreas_sub,
  db = "GO_BP_sim",
  group_by = "CellType",
  plot_type = "comparison"
)

# Or use "geneID" and "geneID_groups" as input to run enrichment
de_df <- dplyr::filter(
  pancreas_sub@tools$DEtest_CellType$AllMarkers_wilcox,
  p_val_adj < 0.05
)
enrich_out <- RunEnrichment(
  geneID = de_df[["gene"]],
  geneID_groups = de_df[["group1"]],
  db = "GO_BP",
  species = "Mus_musculus"
)
EnrichmentPlot(
  res = enrich_out,
  db = "GO_BP",
  plot_type = "comparison"
)

# Use a combined database
pancreas_sub <- RunEnrichment(
  pancreas_sub,
  group_by = "CellType",
  db = c(
    "KEGG", "WikiPathway", "Reactome", "PFAM", "MP"
  ),
  db_combine = TRUE,
  species = "Mus_musculus"
)
EnrichmentPlot(
  pancreas_sub,
  db = "Combined",
  group_by = "CellType",
  plot_type = "comparison"
)
} # }
```
