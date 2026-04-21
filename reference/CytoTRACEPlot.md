# Plot CytoTRACE 2 Results

Plot CytoTRACE 2 Results

## Usage

``` r
CytoTRACEPlot(
  srt,
  reduction = NULL,
  group.by = NULL,
  combine = TRUE,
  nrow = NULL,
  ncol = NULL,
  byrow = TRUE,
  pt.size = NULL,
  pt.alpha = 1,
  palette = "Chinese",
  palcolor = NULL,
  theme_use = "theme_scop",
  theme_args = list(),
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- reduction:

  Which dimensionality reduction to use. If not specified, will use the
  reduction returned by
  [DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- group.by:

  Name of one or more meta.data columns to group (color) cells by.

- combine:

  Combine plots into a single `patchwork` object. If `FALSE`, return a
  list of ggplot objects.

- nrow:

  Number of rows in the combined plot. Default is `NULL`, which means
  determined automatically based on the number of plots.

- ncol:

  Number of columns in the combined plot. Default is `NULL`, which means
  determined automatically based on the number of plots.

- byrow:

  Whether to arrange the plots by row in the combined plot. Default is
  `TRUE`.

- pt.size:

  The size of the points in the plot.

- pt.alpha:

  The transparency of the data points. Default is `1`.

- palette:

  Color palette name. Available palettes can be found in
  [thisplot::show_palettes](https://mengxu98.github.io/thisplot/reference/show_palettes.html).
  Default is `"Chinese"`.

- palcolor:

  Custom colors used to create a color palette. Default is `NULL`.

- theme_use:

  Theme used. Can be a character string or a theme function. Default is
  `"theme_scop"`.

- theme_args:

  Other arguments passed to the `theme_use`. Default is
  [`list()`](https://rdrr.io/r/base/list.html).

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments to be passed to
  [CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md)
  and
  [FeatureDimPlot](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md).

## Value

If `combine = TRUE`, returns a `patchwork` object combining all plots.
If `combine = FALSE`, returns a named list of ggplot objects:

- `Score`: UMAP plot colored by score computed by CytoTRACE2;

- `Potency`: UMAP plot colored by potency category computed by
  CytoTRACE2;

- `Relative`: UMAP plot colored by relative score computed by
  CytoTRACE2;

- `Phenotype`: UMAP plot colored by phenotype (if `group.by` is
  provided);

- `Boxplot`: Boxplot of score computed by CytoTRACE2 corresponding to
  phenotype (if `group.by` is provided).

## See also

[RunCytoTRACE](https://mengxu98.github.io/scop/reference/RunCytoTRACE.md),
[CellDimPlot](https://mengxu98.github.io/scop/reference/CellDimPlot.md),
[FeatureDimPlot](https://mengxu98.github.io/scop/reference/FeatureDimPlot.md)

## Examples

``` r
if (thisutils::check_ci_env()) {
  data(pancreas_sub)
  pancreas_sub <- standard_scop(pancreas_sub)
  pancreas_sub <- RunCytoTRACE(
    pancreas_sub,
    species = "Mus_musculus"
  )

  CytoTRACEPlot(
    pancreas_sub,
    group.by = "CellType"
  )

  plots <- CytoTRACEPlot(
    pancreas_sub,
    group.by = "CellType",
    combine = FALSE
  )
  plots$Boxplot
}
#> ℹ [2026-04-21 06:32:14] Start standard processing workflow...
#> ℹ [2026-04-21 06:32:15] Checking a list of <Seurat>...
#> ! [2026-04-21 06:32:15] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-21 06:32:15] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-21 06:32:17] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-21 06:32:17] Use the separate HVF from `srt_list`
#> ℹ [2026-04-21 06:32:17] Number of available HVF: 2000
#> ℹ [2026-04-21 06:32:17] Finished check
#> ℹ [2026-04-21 06:32:18] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-21 06:32:18] Perform pca linear dimension reduction
#> ℹ [2026-04-21 06:32:19] Use stored estimated dimensions 1:20 for Standardpca
#> ℹ [2026-04-21 06:32:20] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-21 06:32:20] Reorder clusters...
#> ℹ [2026-04-21 06:32:20] Skip `log1p()` because `layer = data` is not "counts"
#> ℹ [2026-04-21 06:32:20] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-21 06:32:20] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ℹ [2026-04-21 06:32:23] Perform umap nonlinear dimension reduction using Standardpca (1:20)
#> ✔ [2026-04-21 06:32:26] Standard processing workflow completed
#> ◌ [2026-04-21 06:32:26] Running CytoTRACE2
#> ℹ [2026-04-21 06:32:26] Package CytoTRACE2 is not installed. Installing from GitHub...
#>  
#> → Will install 62 packages.
#> → Will download 1 CRAN package (2.11 MB), cached: 61 (0 B).
#> + CytoTRACE2         1.1.0     [bld][cmp] (GitHub: 4a398a4)
#> + HiClimR            2.2.1      + ✔ libnetcdf-dev
#> + RANN               2.6.2     
#> + ROCR               1.0-12    
#> + RcppEigen          0.3.4.0.2 
#> + RcppHNSW           0.6.0     
#> + RcppParallel       5.1.11-2   + ✔ make
#> + RcppProgress       0.4.2     
#> + RcppTOML           0.2.3     
#> + Rfast              2.1.5.2   
#> + Rtsne              0.17      
#> + Seurat             5.4.0     [bld][cmp][dl] (2.11 MB)
#> + SeuratObject       5.4.0     
#> + bitops             1.0-9     
#> + caTools            1.18.3    
#> + cowplot            1.2.0     
#> + crosstalk          1.2.2     
#> + data.table         1.18.2.1  
#> + deldir             2.0-4     
#> + dotCall64          1.2       
#> + fastDummies        1.7.5     
#> + fitdistrplus       1.2-6     
#> + future             1.70.0    
#> + future.apply       1.20.2    
#> + ggrepel            0.9.8     
#> + globals            0.19.1    
#> + goftest            1.2-3     
#> + gplots             3.3.0     
#> + gridExtra          2.3       
#> + gtools             3.9.5     
#> + here               1.0.2     
#> + httpuv             1.6.17     + ✔ make, ✔ zlib1g-dev
#> + ica                1.0-3     
#> + later              1.4.8     
#> + lazyeval           0.2.3     
#> + listenv            0.10.1    
#> + lmtest             0.9-40    
#> + miniUI             0.1.2     
#> + ncdf4              1.24       + ✔ libnetcdf-dev
#> + otel               0.2.0     
#> + parallelly         1.47.0    
#> + pbapply            1.7-4     
#> + plotly             4.12.0    
#> + polyclip           1.10-7    
#> + progressr          0.19.0    
#> + promises           1.5.0     
#> + reticulate         1.46.0     + ✔ python3
#> + scattermore        1.2       
#> + sctransform        0.4.3     
#> + shiny              1.13.0    
#> + sourcetools        0.1.7-2   
#> + spam               2.11-3    
#> + spatstat.data      3.1-9     
#> + spatstat.explore   3.8-0     
#> + spatstat.geom      3.7-3     
#> + spatstat.random    3.4-5     
#> + spatstat.sparse    3.1-0     
#> + spatstat.univar    3.1-7     
#> + spatstat.utils     3.2-2     
#> + tensor             1.5.1     
#> + zigg               0.0.2     
#> + zoo                1.8-15    
#> ✔ All system requirements are already installed.
#>   
#> ℹ Getting 1 pkg (2.11 MB), 61 cached
#> ✔ Cached copy of Seurat 5.4.0 (source) is the latest build
#> ✔ Got zigg 0.0.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (25.57 kB)
#> ✔ Got ncdf4 1.24 (x86_64-pc-linux-gnu-ubuntu-24.04) (281.21 kB)
#> ✔ Got HiClimR 2.2.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (575.27 kB)
#> ✔ Got Rfast 2.1.5.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (3.00 MB)
#> ✔ Got CytoTRACE2 1.1.0 (source) (182.95 MB)
#> ℹ Installing system requirements
#> ℹ Executing `sudo sh -c apt-get -y update`
#> Get:1 file:/etc/apt/apt-mirrors.txt Mirrorlist [144 B]
#> Hit:2 http://azure.archive.ubuntu.com/ubuntu noble InRelease
#> Hit:6 https://packages.microsoft.com/repos/azure-cli noble InRelease
#> Hit:3 http://azure.archive.ubuntu.com/ubuntu noble-updates InRelease
#> Hit:4 http://azure.archive.ubuntu.com/ubuntu noble-backports InRelease
#> Hit:5 http://azure.archive.ubuntu.com/ubuntu noble-security InRelease
#> Hit:7 https://packages.microsoft.com/repos/edge stable InRelease
#> Hit:8 https://packages.microsoft.com/ubuntu/24.04/prod noble InRelease
#> Hit:9 https://dl.google.com/linux/chrome-stable/deb stable InRelease
#> Reading package lists...
#> ℹ Executing `sudo sh -c apt-get -y install libnetcdf-dev make zlib1g-dev python3 libcurl4-openssl-dev libssl-dev cmake libuv1-dev libglpk-dev libxml2-dev pandoc libpng-dev libicu-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> libnetcdf-dev is already the newest version (1:4.9.2-5ubuntu4).
#> libnetcdf-dev set to manually installed.
#> make is already the newest version (4.3-4.1build2).
#> zlib1g-dev is already the newest version (1:1.3.dfsg-3.1ubuntu2.1).
#> python3 is already the newest version (3.12.3-0ubuntu2.1).
#> libcurl4-openssl-dev is already the newest version (8.5.0-2ubuntu10.8).
#> libssl-dev is already the newest version (3.0.13-0ubuntu3.9).
#> cmake is already the newest version (3.28.3-1build7).
#> libuv1-dev is already the newest version (1.48.0-1.1build1).
#> libglpk-dev is already the newest version (5.0-1build2).
#> libxml2-dev is already the newest version (2.9.14+dfsg-1.3ubuntu3.7).
#> pandoc is already the newest version (3.1.3+ds-2).
#> libpng-dev is already the newest version (1.6.43-5ubuntu0.5).
#> libicu-dev is already the newest version (74.2-1ubuntu3.1).
#> 0 upgraded, 0 newly installed, 0 to remove and 19 not upgraded.
#> ✔ Installed bitops 1.0-9  (59ms)
#> ✔ Installed caTools 1.18.3  (77ms)
#> ✔ Installed cowplot 1.2.0  (119ms)
#> ✔ Installed crosstalk 1.2.2  (176ms)
#> ✔ Installed data.table 1.18.2.1  (142ms)
#> ✔ Installed deldir 2.0-4  (73ms)
#> ✔ Installed dotCall64 1.2  (65ms)
#> ✔ Installed fastDummies 1.7.5  (66ms)
#> ✔ Installed fitdistrplus 1.2-6  (65ms)
#> ✔ Installed future 1.70.0  (68ms)
#> ✔ Installed future.apply 1.20.2  (108ms)
#> ✔ Installed ggrepel 0.9.8  (71ms)
#> ✔ Installed globals 0.19.1  (67ms)
#> ✔ Installed goftest 1.2-3  (66ms)
#> ✔ Installed gplots 3.3.0  (1.1s)
#> ✔ Installed gridExtra 2.3  (1.1s)
#> ✔ Installed gtools 3.9.5  (100ms)
#> ✔ Installed here 1.0.2  (68ms)
#> ✔ Installed HiClimR 2.2.1  (71ms)
#> ✔ Installed httpuv 1.6.17  (73ms)
#> ✔ Installed ica 1.0-3  (65ms)
#> ✔ Installed later 1.4.8  (62ms)
#> ✔ Installed lazyeval 0.2.3  (105ms)
#> ✔ Installed listenv 0.10.1  (104ms)
#> ✔ Installed lmtest 0.9-40  (66ms)
#> ✔ Installed miniUI 0.1.2  (64ms)
#> ✔ Installed ncdf4 1.24  (62ms)
#> ✔ Installed otel 0.2.0  (62ms)
#> ✔ Installed parallelly 1.47.0  (64ms)
#> ✔ Installed pbapply 1.7-4  (102ms)
#> ✔ Installed plotly 4.12.0  (104ms)
#> ✔ Installed polyclip 1.10-7  (66ms)
#> ✔ Installed progressr 0.19.0  (66ms)
#> ✔ Installed promises 1.5.0  (68ms)
#> ✔ Installed RANN 2.6.2  (65ms)
#> ✔ Installed RcppHNSW 0.6.0  (1s)
#> ✔ Installed RcppEigen 0.3.4.0.2  (1.1s)
#> ✔ Installed RcppParallel 5.1.11-2  (72ms)
#> ✔ Installed RcppProgress 0.4.2  (69ms)
#> ✔ Installed RcppTOML 0.2.3  (66ms)
#> ✔ Installed reticulate 1.46.0  (67ms)
#> ✔ Installed Rfast 2.1.5.2  (67ms)
#> ✔ Installed ROCR 1.0-12  (106ms)
#> ✔ Installed Rtsne 0.17  (68ms)
#> ✔ Installed scattermore 1.2  (64ms)
#> ✔ Installed sctransform 0.4.3  (67ms)
#> ✔ Installed SeuratObject 5.4.0  (1.1s)
#> ✔ Installed Seurat 5.4.0  (1.2s)
#> ✔ Installed shiny 1.13.0  (162ms)
#> ✔ Installed sourcetools 0.1.7-2  (72ms)
#> ✔ Installed spam 2.11-3  (67ms)
#> ✔ Installed spatstat.data 3.1-9  (76ms)
#> ✔ Installed spatstat.explore 3.8-0  (75ms)
#> ✔ Installed spatstat.geom 3.7-3  (69ms)
#> ✔ Installed spatstat.random 3.4-5  (65ms)
#> ✔ Installed spatstat.sparse 3.1-0  (108ms)
#> ✔ Installed spatstat.univar 3.1-7  (68ms)
#> ✔ Installed spatstat.utils 3.2-2  (66ms)
#> ✔ Installed tensor 1.5.1  (70ms)
#> ✔ Installed zigg 0.0.2  (63ms)
#> ✔ Installed zoo 1.8-15  (47ms)
#> ℹ Packaging CytoTRACE2 1.1.0
#> ✔ Packaged CytoTRACE2 1.1.0 (4.8s)
#> ℹ Building CytoTRACE2 1.1.0
#> ✔ Built CytoTRACE2 1.1.0 (9.5s)
#> ✔ Installed CytoTRACE2 1.1.0 (github::digitalcytometry/cytotrace2@4a398a4) (1.1s)
#> ✔ 1 pkg + 145 deps: kept 82, added 62, dld 5 (NA B) [53.3s]
#> Warning: replacing previous import ‘data.table::first’ by ‘dplyr::first’ when loading ‘CytoTRACE2’
#> Warning: replacing previous import ‘data.table::between’ by ‘dplyr::between’ when loading ‘CytoTRACE2’
#> Warning: replacing previous import ‘data.table::last’ by ‘dplyr::last’ when loading ‘CytoTRACE2’
#> cytotrace2: Started loading data
#> Dataset contains 15998 genes and 1000 cells.
#> The number of cells in your dataset is less than 1000. Fast mode has been disabled.
#> The passed subsample size is greater than the number of cells in dataset.
#> Now setting subsample size to 1000
#> cytotrace2: Running on 1 subsample(s) approximately of length 1000
#> cytotrace2: Started running on subsample(s). This will take a few minutes.
#> cytotrace2: Started preprocessing.
#> 12486 input genes mapped to model genes.
#> cytotrace2: Started prediction.
#> This section will run using  1 / 4 core(s).
#> cytotrace2: Started postprocessing.
#> cytotrace2: Running with slow mode (subsamples are processed sequentially)
#> Number of cores for KNN: 1
#> cytotrace2: Finished
#> ✔ [2026-04-21 06:34:41] CytoTRACE2 computed successfully
```
