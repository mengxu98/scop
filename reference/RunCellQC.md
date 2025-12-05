# Run cell-level quality control for single cell RNA-seq data.

This function handles multiple quality control methods for single-cell
RNA-seq data.

## Usage

``` r
RunCellQC(
  srt,
  assay = "RNA",
  split.by = NULL,
  return_filtered = FALSE,
  qc_metrics = c("doublets", "outlier", "umi", "gene", "mito", "ribo", "ribo_mito_ratio",
    "species"),
  db_method = "scDblFinder",
  db_rate = NULL,
  db_coefficient = 0.01,
  outlier_threshold = c("log10_nCount:lower:2.5", "log10_nCount:higher:5",
    "log10_nFeature:lower:2.5", "log10_nFeature:higher:5", "featurecount_dist:lower:2.5"),
  outlier_n = 1,
  UMI_threshold = 3000,
  gene_threshold = 1000,
  mito_threshold = 20,
  mito_pattern = c("MT-", "Mt-", "mt-"),
  mito_gene = NULL,
  ribo_threshold = 50,
  ribo_pattern = c("RP[SL]\\d+\\w{0,1}\\d*$", "Rp[sl]\\d+\\w{0,1}\\d*$",
    "rp[sl]\\d+\\w{0,1}\\d*$"),
  ribo_gene = NULL,
  ribo_mito_ratio_range = c(1, Inf),
  species = NULL,
  species_gene_prefix = NULL,
  species_percent = 95,
  seed = 11
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  The name of the assay to be used for doublet-calling. Default is
  `"RNA"`.

- split.by:

  Name of the sample variable to split the Seurat object. Default is
  `NULL`.

- return_filtered:

  Logical indicating whether to return a cell-filtered Seurat object.
  Default is `FALSE`.

- qc_metrics:

  A character vector specifying the quality control metrics to be
  applied. Default is
  `c("doublets", "outlier", "umi", "gene", "mito", "ribo", "ribo_mito_ratio", "species")`.

- db_method:

  Method used for doublet-calling. Can be one of `"scDblFinder"`,
  `"Scrublet"`, `"DoubletDetection"`, `"scds_cxds"`, `"scds_bcds"`,
  `"scds_hybrid"`.

- db_rate:

  The expected doublet rate. Default is calculated as
  `ncol(srt) / 1000 * 0.01`.

- db_coefficient:

  The coefficient used to calculate the doublet rate. Default is `0.01`.
  Doublet rate is calculated as `ncol(srt) / 1000 * db_coefficient`.

- outlier_threshold:

  A character vector specifying the outlier threshold. Default is
  `c("log10_nCount:lower:2.5", "log10_nCount:higher:5", "log10_nFeature:lower:2.5", "log10_nFeature:higher:5", "featurecount_dist:lower:2.5")`.
  See
  [scuttle::isOutlier](https://rdrr.io/pkg/scuttle/man/isOutlier.html).

- outlier_n:

  Minimum number of outlier metrics that meet the conditions for
  determining outlier cells. Default is `1`.

- UMI_threshold:

  UMI number threshold. Cells that exceed this threshold will be
  considered as kept. Default is `3000`.

- gene_threshold:

  Gene number threshold. Cells that exceed this threshold will be
  considered as kept. Default is `1000`.

- mito_threshold:

  Percentage of UMI counts of mitochondrial genes. Cells that exceed
  this threshold will be considered as discarded. Default is `20`.

- mito_pattern:

  Regex patterns to match the mitochondrial genes. Default is
  `c("MT-", "Mt-", "mt-")`.

- mito_gene:

  A defined mitochondrial genes. If features provided, will ignore the
  `mito_pattern` matching. Default is `NULL`.

- ribo_threshold:

  Percentage of UMI counts of ribosomal genes. Cells that exceed this
  threshold will be considered as discarded. Default is `50`.

- ribo_pattern:

  Regex patterns to match the ribosomal genes. Default is
  `c("RP[SL]\\d+\\w{0,1}\\d*$", "Rp[sl]\\d+\\w{0,1}\\d*$", "rp[sl]\\d+\\w{0,1}\\d*$")`.

- ribo_gene:

  A defined ribosomal genes. If features provided, will ignore the
  `ribo_pattern` matching. Default is `NULL`.

- ribo_mito_ratio_range:

  A numeric vector specifying the range of ribosomal/mitochondrial gene
  expression ratios for ribo_mito_ratio outlier cells. Default is
  `c(1, Inf)`.

- species:

  Species used as the suffix of the QC metrics. The first is the species
  of interest. Default is `NULL`.

- species_gene_prefix:

  Species gene prefix used to calculate QC metrics for each species.
  Default is `NULL`.

- species_percent:

  Percentage of UMI counts of the first species. Cells that exceed this
  threshold will be considered as kept. Default is `95`.

- seed:

  Set a random seed. Default is `11`.

## Value

Returns Seurat object with the QC results stored in the meta.data layer.

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
pancreas_sub <- RunCellQC(pancreas_sub)
#>  
#> → Will install 57 packages.
#> → Will download 1 CRAN package (14.98 kB), cached: 56 (0 B).
#> + AnnotationFilter    1.34.0    [bld]
#> + BiocIO              1.20.0    [bld]
#> + BiocSingular        1.26.1    [bld][cmp]
#> + Cairo               1.7-0      + ✔ libcairo2-dev
#> + ClusterR            1.3.5     
#> + ExperimentHub       3.0.0     [bld]
#> + GenomeInfoDb        1.46.1    [bld]
#> + GenomicAlignments   1.46.0    [bld][cmp]
#> + GenomicFeatures     1.62.0    [bld]
#> + HDF5Array           1.38.0    [bld] + ✔ make
#> + ProtGenerics        1.42.0    [bld]
#> + RCurl               1.98-1.17  + ✔ make, ✔ libcurl4-openssl-dev
#> + RcppML              0.3.7     
#> + Rhdf5lib            1.32.0    [bld][cmp] + ✔ make
#> + Rhtslib             3.6.0     [bld][cmp] + ✔ libbz2-dev, ✔ libcurl4-openssl-dev, ✔ liblzma-dev
#> + Rsamtools           2.26.0    [bld][cmp] + ✔ make
#> + ScaledMatrix        1.18.0    [bld]
#> + UCSC.utils          1.6.0     [bld]
#> + V8                  8.0.1      + ✖ libnode-dev
#> + XML                 3.99-0.20  + ✔ libxml2-dev
#> + alabaster.base      1.10.0    [bld][cmp] + ✔ make
#> + alabaster.matrix    1.10.0    [bld][cmp]
#> + alabaster.ranges    1.10.0    [bld]
#> + alabaster.sce       1.10.0    [bld]
#> + alabaster.schemas   1.10.0    [bld]
#> + alabaster.se        1.10.0    [bld]
#> + beachmat            2.26.0    [bld][cmp]
#> + beeswarm            0.4.0     
#> + benchmarkme         1.0.8     
#> + benchmarkmeData     1.0.4     
#> + bluster             1.20.0    [bld][cmp]
#> + cigarillo           1.0.0     [bld][cmp]
#> + dqrng               0.4.1     
#> + edgeR               4.8.0     [bld][cmp]
#> + ensembldb           2.34.0    [bld]
#> + ggbeeswarm          0.7.3     
#> + ggrastr             1.0.2     
#> + gmp                 0.7-5      + ✖ libgmp3-dev
#> + gypsum              1.6.0     [bld]
#> + h5mread             1.2.1     [bld][cmp]
#> + jsonvalidate        1.5.0     
#> + mbkmeans            1.26.0    [bld][cmp]
#> + metapod             1.18.0    [bld][cmp]
#> + pheatmap            1.0.13    
#> + restfulr            0.0.16    [bld][cmp][dl] (14.98 kB)
#> + rhdf5               2.54.0    [bld][cmp] + ✔ make
#> + rhdf5filters        1.22.0    [bld][cmp] + ✔ make
#> + rsvd                1.0.5     
#> + rtracklayer         1.70.0    [bld][cmp]
#> + scDblFinder         1.24.0    [bld]
#> + scRNAseq            2.24.0    [bld]
#> + scater              1.38.0    [bld]
#> + scran               1.38.0    [bld][cmp]
#> + scuttle             1.20.0    [bld][cmp]
#> + vipor               0.4.7     
#> + viridis             0.6.5     
#> + xgboost             3.1.2.1    + ✔ make
#> → Will install 2 system packages:
#> + libgmp3-dev  - gmp
#> + libnode-dev  - V8 
#> ℹ Getting 1 pkg (14.98 kB), 56 cached
#> ✔ Got alabaster.sce 1.10.0 (source) (226.14 kB)
#> ✔ Got alabaster.ranges 1.10.0 (source) (232.01 kB)
#> ✔ Got alabaster.se 1.10.0 (source) (231.15 kB)
#> ✔ Got alabaster.matrix 1.10.0 (source) (282.83 kB)
#> ✔ Got alabaster.schemas 1.10.0 (source) (247.50 kB)
#> ✔ Got restfulr 0.0.16 (source) (15.03 kB)
#> ✔ Got dqrng 0.4.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (215.52 kB)
#> ✔ Got beeswarm 0.4.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (78.13 kB)
#> ✔ Got alabaster.base 1.10.0 (source) (342.94 kB)
#> ✔ Got Cairo 1.7-0 (x86_64-pc-linux-gnu-ubuntu-24.04) (93.28 kB)
#> ✔ Got UCSC.utils 1.6.0 (source) (236.81 kB)
#> ✔ Got ScaledMatrix 1.18.0 (source) (314.08 kB)
#> ✔ Got gypsum 1.6.0 (source) (284.67 kB)
#> ✔ Got benchmarkme 1.0.8 (x86_64-pc-linux-gnu-ubuntu-24.04) (122.90 kB)
#> ✔ Got benchmarkmeData 1.0.4 (x86_64-pc-linux-gnu-ubuntu-24.04) (282.00 kB)
#> ✔ Got ClusterR 1.3.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.20 MB)
#> ✔ Got RCurl 1.98-1.17 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.07 MB)
#> ✔ Got scuttle 1.20.0 (source) (1.03 MB)
#> ✔ Got Rsamtools 2.26.0 (source) (1.92 MB)
#> ✔ Got GenomicAlignments 1.46.0 (source) (2.26 MB)
#> ✔ Got edgeR 4.8.0 (source) (3.07 MB)
#> ✔ Got pheatmap 1.0.13 (x86_64-pc-linux-gnu-ubuntu-24.04) (78.32 kB)
#> ✔ Got jsonvalidate 1.5.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (177.69 kB)
#> ✔ Got AnnotationFilter 1.34.0 (source) (325.42 kB)
#> ✔ Got gmp 0.7-5 (x86_64-pc-linux-gnu-ubuntu-24.04) (330.19 kB)
#> ✔ Got scater 1.38.0 (source) (4.57 MB)
#> ✔ Got BiocSingular 1.26.1 (source) (618.28 kB)
#> ✔ Got scran 1.38.0 (source) (1.83 MB)
#> ✔ Got ExperimentHub 3.0.0 (source) (503.09 kB)
#> ✔ Got ProtGenerics 1.42.0 (source) (11.78 kB)
#> ✔ Got mbkmeans 1.26.0 (source) (300.82 kB)
#> ✔ Got GenomicFeatures 1.62.0 (source) (574.45 kB)
#> ✔ Got RcppML 0.3.7 (x86_64-pc-linux-gnu-ubuntu-24.04) (206.79 kB)
#> ✔ Got ggbeeswarm 0.7.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.77 MB)
#> ✔ Got rtracklayer 1.70.0 (source) (4.10 MB)
#> ✔ Got beachmat 2.26.0 (source) (383.08 kB)
#> ✔ Got cigarillo 1.0.0 (source) (258.57 kB)
#> ✔ Got rhdf5filters 1.22.0 (source) (1.19 MB)
#> ✔ Got metapod 1.18.0 (source) (333.89 kB)
#> ✔ Got ggrastr 1.0.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.68 MB)
#> ✔ Got BiocIO 1.20.0 (source) (232.37 kB)
#> ✔ Got scDblFinder 1.24.0 (source) (2.25 MB)
#> ✔ Got scRNAseq 2.24.0 (source) (311.65 kB)
#> ✔ Got HDF5Array 1.38.0 (source) (8.36 MB)
#> ✔ Got rsvd 1.0.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (3.59 MB)
#> ✔ Got h5mread 1.2.1 (source) (2.48 MB)
#> ✔ Got XML 3.99-0.20 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.83 MB)
#> ✔ Got GenomeInfoDb 1.46.1 (source) (3.65 MB)
#> ✔ Got viridis 0.6.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (3.01 MB)
#> ✔ Got vipor 0.4.7 (x86_64-pc-linux-gnu-ubuntu-24.04) (4.58 MB)
#> ✔ Got rhdf5 2.54.0 (source) (1.31 MB)
#> ✔ Got ensembldb 2.34.0 (source) (3.54 MB)
#> ✔ Got Rhtslib 3.6.0 (source) (5.18 MB)
#> ✔ Got V8 8.0.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (12.16 MB)
#> ✔ Got bluster 1.20.0 (source) (3.26 MB)
#> ✔ Got Rhdf5lib 1.32.0 (source) (12.07 MB)
#> ℹ Installing system requirements
#> ℹ Executing `sudo sh -c apt-get -y update`
#> Get:1 file:/etc/apt/apt-mirrors.txt Mirrorlist [144 B]
#> Hit:2 http://azure.archive.ubuntu.com/ubuntu noble InRelease
#> Hit:3 http://azure.archive.ubuntu.com/ubuntu noble-updates InRelease
#> Hit:4 http://azure.archive.ubuntu.com/ubuntu noble-backports InRelease
#> Hit:5 http://azure.archive.ubuntu.com/ubuntu noble-security InRelease
#> Hit:6 https://packages.microsoft.com/repos/azure-cli noble InRelease
#> Hit:7 https://packages.microsoft.com/ubuntu/24.04/prod noble InRelease
#> Reading package lists...
#> ℹ Executing `sudo sh -c apt-get -y install libcairo2-dev libgmp3-dev make libcurl4-openssl-dev libnode-dev libxml2-dev libbz2-dev liblzma-dev perl pandoc libssl-dev libglpk-dev libpng-dev libfreetype6-dev libjpeg-dev libtiff-dev libwebp-dev libicu-dev libfontconfig1-dev libfribidi-dev libharfbuzz-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> libcairo2-dev is already the newest version (1.18.0-3build1).
#> make is already the newest version (4.3-4.1build2).
#> libcurl4-openssl-dev is already the newest version (8.5.0-2ubuntu10.6).
#> libxml2-dev is already the newest version (2.9.14+dfsg-1.3ubuntu3.6).
#> libbz2-dev is already the newest version (1.0.8-5.1build0.1).
#> liblzma-dev is already the newest version (5.6.1+really5.4.5-1ubuntu0.2).
#> perl is already the newest version (5.38.2-3.2ubuntu0.2).
#> pandoc is already the newest version (3.1.3+ds-2).
#> libssl-dev is already the newest version (3.0.13-0ubuntu3.6).
#> libglpk-dev is already the newest version (5.0-1build2).
#> libpng-dev is already the newest version (1.6.43-5build1).
#> libfreetype-dev is already the newest version (2.13.2+dfsg-1build3).
#> libjpeg-dev is already the newest version (8c-2ubuntu11).
#> libtiff-dev is already the newest version (4.5.1+git230720-4ubuntu2.4).
#> libwebp-dev is already the newest version (1.3.2-0.4build3).
#> libicu-dev is already the newest version (74.2-1ubuntu3.1).
#> libfontconfig1-dev is already the newest version (2.15.0-1.1ubuntu2).
#> libfribidi-dev is already the newest version (1.0.13-3build1).
#> libharfbuzz-dev is already the newest version (8.3.0-2build2).
#> The following additional packages will be installed:
#> libnode109 libuv1-dev node-acorn node-busboy node-cjs-module-lexer
#>   node-undici node-xtend nodejs nodejs-doc
#> Suggested packages:
#>   npm
#> The following NEW packages will be installed:
#> libgmp3-dev libnode-dev libnode109 libuv1-dev node-acorn node-busboy
#>   node-cjs-module-lexer node-undici node-xtend nodejs nodejs-doc
#> 0 upgraded, 11 newly installed, 0 to remove and 49 not upgraded.
#> Need to get 16.6 MB of archives.
#> After this operation, 74.2 MB of additional disk space will be used.
#> Get:1 file:/etc/apt/apt-mirrors.txt Mirrorlist [144 B]
#> Get:2 http://azure.archive.ubuntu.com/ubuntu noble-updates/main amd64 libgmp3-dev amd64 2:6.3.0+dfsg-2ubuntu6.1 [2310 B]
#> Get:3 http://azure.archive.ubuntu.com/ubuntu noble/main amd64 libuv1-dev amd64 1.48.0-1.1build1 [136 kB]
#> Get:4 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 node-xtend all 4.0.2-3 [3902 B]
#> Get:5 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 nodejs amd64 18.19.1+dfsg-6ubuntu5 [306 kB]
#> Get:6 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 node-acorn all 8.8.1+ds+~cs25.17.7-2 [115 kB]
#> Get:7 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 node-cjs-module-lexer all 1.2.3+dfsg-1 [32.1 kB]
#> Get:8 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 node-busboy all 1.6.0+~cs2.6.0-2 [17.3 kB]
#> Get:9 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 node-undici all 5.26.3+dfsg1+~cs23.10.12-2 [325 kB]
#> Get:10 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 libnode109 amd64 18.19.1+dfsg-6ubuntu5 [11.6 MB]
#> Get:11 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 libnode-dev amd64 18.19.1+dfsg-6ubuntu5 [501 kB]
#> Get:12 http://azure.archive.ubuntu.com/ubuntu noble/universe amd64 nodejs-doc all 18.19.1+dfsg-6ubuntu5 [3552 kB]
#> Fetched 16.6 MB in 0s (111 MB/s)
#> Selecting previously unselected package libgmp3-dev:amd64.
#> (Reading database ...
#> (Reading database ... 5%(Reading database ... 10%(Reading database ... 15%(Reading database ... 20%(Reading database ... 25%(Reading database ... 30%(Reading database ... 35%(Reading database ... 40%(Reading database ... 45%(Reading database ... 50%(Reading database ... 55%
#> (Reading database ... 60%
#> (Reading database ... 65%
#> (Reading database ... 70%
#> (Reading database ... 75%
#> (Reading database ... 80%
#> (Reading database ... 85%
#> (Reading database ... 90%
#> (Reading database ... 95%
#> (Reading database ... 100%(Reading database ... 255135 files and directories currently installed.)
#> Preparing to unpack .../00-libgmp3-dev_2%3a6.3.0+dfsg-2ubuntu6.1_amd64.deb ...
#> Unpacking libgmp3-dev:amd64 (2:6.3.0+dfsg-2ubuntu6.1) ...
#> Selecting previously unselected package libuv1-dev:amd64.
#> Preparing to unpack .../01-libuv1-dev_1.48.0-1.1build1_amd64.deb ...
#> Unpacking libuv1-dev:amd64 (1.48.0-1.1build1) ...
#> Selecting previously unselected package node-xtend.
#> Preparing to unpack .../02-node-xtend_4.0.2-3_all.deb ...
#> Unpacking node-xtend (4.0.2-3) ...
#> Selecting previously unselected package nodejs.
#> Preparing to unpack .../03-nodejs_18.19.1+dfsg-6ubuntu5_amd64.deb ...
#> Unpacking nodejs (18.19.1+dfsg-6ubuntu5) ...
#> Selecting previously unselected package node-acorn.
#> Preparing to unpack .../04-node-acorn_8.8.1+ds+~cs25.17.7-2_all.deb ...
#> Unpacking node-acorn (8.8.1+ds+~cs25.17.7-2) ...
#> Selecting previously unselected package node-cjs-module-lexer.
#> Preparing to unpack .../05-node-cjs-module-lexer_1.2.3+dfsg-1_all.deb ...
#> Unpacking node-cjs-module-lexer (1.2.3+dfsg-1) ...
#> Selecting previously unselected package node-busboy.
#> Preparing to unpack .../06-node-busboy_1.6.0+~cs2.6.0-2_all.deb ...
#> Unpacking node-busboy (1.6.0+~cs2.6.0-2) ...
#> Selecting previously unselected package node-undici.
#> Preparing to unpack .../07-node-undici_5.26.3+dfsg1+~cs23.10.12-2_all.deb ...
#> Unpacking node-undici (5.26.3+dfsg1+~cs23.10.12-2) ...
#> Selecting previously unselected package libnode109:amd64.
#> Preparing to unpack .../08-libnode109_18.19.1+dfsg-6ubuntu5_amd64.deb ...
#> Unpacking libnode109:amd64 (18.19.1+dfsg-6ubuntu5) ...
#> Selecting previously unselected package libnode-dev.
#> Preparing to unpack .../09-libnode-dev_18.19.1+dfsg-6ubuntu5_amd64.deb ...
#> Unpacking libnode-dev (18.19.1+dfsg-6ubuntu5) ...
#> Selecting previously unselected package nodejs-doc.
#> Preparing to unpack .../10-nodejs-doc_18.19.1+dfsg-6ubuntu5_all.deb ...
#> Unpacking nodejs-doc (18.19.1+dfsg-6ubuntu5) ...
#> Setting up libuv1-dev:amd64 (1.48.0-1.1build1) ...
#> Setting up node-cjs-module-lexer (1.2.3+dfsg-1) ...
#> Setting up nodejs-doc (18.19.1+dfsg-6ubuntu5) ...
#> Setting up libgmp3-dev:amd64 (2:6.3.0+dfsg-2ubuntu6.1) ...
#> Setting up node-xtend (4.0.2-3) ...
#> Setting up node-busboy (1.6.0+~cs2.6.0-2) ...
#> Setting up node-undici (5.26.3+dfsg1+~cs23.10.12-2) ...
#> Setting up node-acorn (8.8.1+ds+~cs25.17.7-2) ...
#> Setting up libnode109:amd64 (18.19.1+dfsg-6ubuntu5) ...
#> Setting up nodejs (18.19.1+dfsg-6ubuntu5) ...
#> update-alternatives: using /usr/bin/nodejs to provide /usr/bin/js (js) in auto mode
#> Setting up libnode-dev (18.19.1+dfsg-6ubuntu5) ...
#> Processing triggers for libc-bin (2.39-0ubuntu8.6) ...
#> Processing triggers for man-db (2.12.0-4build2) ...
#> Not building database; man-db/auto-update is not 'true'.
#> Running kernel seems to be up-to-date.
#> 
#> Restarting services...
#> Service restarts being deferred:
#>  systemctl restart networkd-dispatcher.service
#> No containers need to be restarted.
#> 
#> No user sessions are running outdated binaries.
#> 
#> No VM guests are running outdated hypervisor (qemu) binaries on this host.
#> ℹ Building alabaster.schemas 1.10.0
#> ℹ Building AnnotationFilter 1.34.0
#> ℹ Building beachmat 2.26.0
#> ℹ Building BiocIO 1.20.0
#> ✔ Built alabaster.schemas 1.10.0 (1.3s)
#> ℹ Building bluster 1.20.0
#> ✔ Built BiocIO 1.20.0 (4.6s)
#> ℹ Building cigarillo 1.0.0
#> ✔ Built AnnotationFilter 1.34.0 (9.1s)
#> ℹ Building edgeR 4.8.0
#> ✔ Built cigarillo 1.0.0 (9.7s)
#> ℹ Building ExperimentHub 3.0.0
#> ✔ Built bluster 1.20.0 (22.2s)
#> ℹ Building gypsum 1.6.0
#> ✔ Built edgeR 4.8.0 (14.9s)
#> ℹ Building metapod 1.18.0
#> ✔ Built ExperimentHub 3.0.0 (11.1s)
#> ℹ Building ProtGenerics 1.42.0
#> ✔ Built gypsum 1.6.0 (3.7s)
#> ℹ Building Rhdf5lib 1.32.0
#> ✔ Built ProtGenerics 1.42.0 (4.3s)
#> ℹ Building Rhtslib 3.6.0
#> ✔ Built metapod 1.18.0 (24s)
#> ℹ Building ScaledMatrix 1.18.0
#> ✔ Built ScaledMatrix 1.18.0 (16.6s)
#> ℹ Building UCSC.utils 1.6.0
#> ✔ Built UCSC.utils 1.6.0 (3.8s)
#> ✔ Installed beeswarm 0.4.0  (31ms)
#> ✔ Installed benchmarkme 1.0.8  (35ms)
#> ✔ Installed benchmarkmeData 1.0.4  (35ms)
#> ✔ Installed Cairo 1.7-0  (33ms)
#> ✔ Installed ClusterR 1.3.5  (59ms)
#> ✔ Installed dqrng 0.4.1  (44ms)
#> ✔ Installed ggbeeswarm 0.7.3  (52ms)
#> ✔ Installed ggrastr 1.0.2  (68ms)
#> ✔ Installed gmp 0.7-5  (38ms)
#> ✔ Installed jsonvalidate 1.5.0  (1.1s)
#> ✔ Installed pheatmap 1.0.13  (32ms)
#> ✔ Installed RcppML 0.3.7  (37ms)
#> ✔ Installed RCurl 1.98-1.17  (55ms)
#> ✔ Installed rsvd 1.0.5  (59ms)
#> ✔ Installed V8 8.0.1  (379ms)
#> ✔ Installed vipor 0.4.7  (129ms)
#> ✔ Installed viridis 0.6.5  (63ms)
#> ✔ Installed xgboost 3.1.2.1  (122ms)
#> ✔ Installed XML 3.99-0.20  (73ms)
#> ℹ Building restfulr 0.0.16
#> ✔ Built restfulr 0.0.16 (7.9s)
#> ✔ Installed restfulr 0.0.16  (54ms)
#> ✔ Installed alabaster.schemas 1.10.0  (43ms)
#> ✔ Installed AnnotationFilter 1.34.0  (47ms)
#> ✔ Installed BiocIO 1.20.0  (44ms)
#> ✔ Installed bluster 1.20.0  (109ms)
#> ✔ Installed cigarillo 1.0.0  (45ms)
#> ✔ Installed edgeR 4.8.0  (87ms)
#> ✔ Installed ExperimentHub 3.0.0  (43ms)
#> ✔ Installed gypsum 1.6.0  (1s)
#> ✔ Installed metapod 1.18.0  (83ms)
#> ✔ Installed ProtGenerics 1.42.0  (33ms)
#> ✔ Installed ScaledMatrix 1.18.0  (44ms)
#> ✔ Installed UCSC.utils 1.6.0  (37ms)
#> ℹ Building GenomeInfoDb 1.46.1
#> ✔ Built GenomeInfoDb 1.46.1 (9.4s)
#> ✔ Installed GenomeInfoDb 1.46.1  (98ms)
#> ✔ Built Rhtslib 3.6.0 (1m 31.1s)
#> ✔ Installed Rhtslib 3.6.0  (209ms)
#> ℹ Building Rsamtools 2.26.0
#> ✔ Built beachmat 2.26.0 (2m 3s)
#> ✔ Installed beachmat 2.26.0  (256ms)
#> ℹ Building BiocSingular 1.26.1
#> ℹ Building scuttle 1.20.0
#> ✔ Built BiocSingular 1.26.1 (26.6s)
#> ✔ Installed BiocSingular 1.26.1  (70ms)
#> ✔ Built Rsamtools 2.26.0 (41.1s)
#> ✔ Installed Rsamtools 2.26.0  (126ms)
#> ℹ Building GenomicAlignments 1.46.0
#> ✔ Built scuttle 1.20.0 (58.6s)
#> ✔ Installed scuttle 1.20.0  (122ms)
#> ℹ Building scater 1.38.0
#> ℹ Building scran 1.38.0
#> ✔ Built GenomicAlignments 1.46.0 (27.3s)
#> ✔ Installed GenomicAlignments 1.46.0  (73ms)
#> ℹ Building rtracklayer 1.70.0
#> ✔ Built scater 1.38.0 (24.9s)
#> ✔ Installed scater 1.38.0  (105ms)
#> ✔ Built rtracklayer 1.70.0 (44.9s)
#> ✔ Installed rtracklayer 1.70.0  (1.1s)
#> ℹ Building GenomicFeatures 1.62.0
#> ✔ Built scran 1.38.0 (1m 6.6s)
#> ✔ Installed scran 1.38.0  (120ms)
#> ℹ Building scDblFinder 1.24.0
#> ✔ Built GenomicFeatures 1.62.0 (22.9s)
#> ✔ Installed GenomicFeatures 1.62.0  (37ms)
#> ℹ Building ensembldb 2.34.0
#> ✔ Built Rhdf5lib 1.32.0 (4m 1.4s)
#> ✔ Installed Rhdf5lib 1.32.0  (1.4s)
#> ℹ Building mbkmeans 1.26.0
#> ℹ Building rhdf5filters 1.22.0
#> ✔ Built scDblFinder 1.24.0 (32.2s)
#> ✔ Installed scDblFinder 1.24.0  (1.1s)
#> ✔ Built rhdf5filters 1.22.0 (14.2s)
#> ✔ Installed rhdf5filters 1.22.0  (43ms)
#> ℹ Building rhdf5 2.54.0
#> ✔ Built ensembldb 2.34.0 (31.1s)
#> ✔ Installed ensembldb 2.34.0  (79ms)
#> ✔ Built rhdf5 2.54.0 (12.1s)
#> ✔ Installed rhdf5 2.54.0  (144ms)
#> ℹ Building alabaster.base 1.10.0
#> ℹ Building h5mread 1.2.1
#> ✔ Built mbkmeans 1.26.0 (39.9s)
#> ✔ Installed mbkmeans 1.26.0  (75ms)
#> ✔ Built h5mread 1.2.1 (15.9s)
#> ✔ Installed h5mread 1.2.1  (1.1s)
#> ℹ Building HDF5Array 1.38.0
#> ✔ Built HDF5Array 1.38.0 (12s)
#> ✔ Installed HDF5Array 1.38.0  (92ms)
#> ✔ Built alabaster.base 1.10.0 (1m 5.2s)
#> ✔ Installed alabaster.base 1.10.0  (1.2s)
#> ℹ Building alabaster.matrix 1.10.0
#> ℹ Building alabaster.ranges 1.10.0
#> ✔ Built alabaster.ranges 1.10.0 (5.6s)
#> ✔ Installed alabaster.ranges 1.10.0  (35ms)
#> ✔ Built alabaster.matrix 1.10.0 (15.6s)
#> ✔ Installed alabaster.matrix 1.10.0  (33ms)
#> ℹ Building alabaster.se 1.10.0
#> ✔ Built alabaster.se 1.10.0 (12.9s)
#> ✔ Installed alabaster.se 1.10.0  (23ms)
#> ℹ Building alabaster.sce 1.10.0
#> ✔ Built alabaster.sce 1.10.0 (13.6s)
#> ✔ Installed alabaster.sce 1.10.0  (37ms)
#> ℹ Building scRNAseq 2.24.0
#> ✔ Built scRNAseq 2.24.0 (21.3s)
#> ✔ Installed scRNAseq 2.24.0  (26ms)
#> ✔ 1 pkg + 202 deps: kept 145, added 57, dld 56 (101.78 MB) [7m 25.9s]
CellStatPlot(
  pancreas_sub,
  stat.by = c(
    "db_qc", "outlier_qc",
    "umi_qc", "gene_qc",
    "mito_qc", "ribo_qc",
    "ribo_mito_ratio_qc", "species_qc"
  ),
  plot_type = "upset",
  stat_level = "Fail"
)
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?

table(pancreas_sub$CellQC)
#> 
#> Pass Fail 
#>  955   45 

data(ifnb_sub)
ifnb_sub <- RunCellQC(
  srt = ifnb_sub,
  split.by = "stim",
  UMI_threshold = 1000,
  gene_threshold = 550
)
CellStatPlot(
  srt = ifnb_sub,
  stat.by = c(
    "db_qc", "outlier_qc",
    "umi_qc", "gene_qc",
    "mito_qc", "ribo_qc",
    "ribo_mito_ratio_qc", "species_qc"
  ),
  plot_type = "upset",
  stat_level = "Fail"
)


table(ifnb_sub$CellQC)
#> 
#> Pass Fail 
#> 1382  618 
```
