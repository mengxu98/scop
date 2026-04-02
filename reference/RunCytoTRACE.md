# Run CytoTRACE 2

Run CytoTRACE 2

## Usage

``` r
RunCytoTRACE(object, ...)

# S3 method for class 'Seurat'
RunCytoTRACE(
  object,
  assay = NULL,
  layer = c("counts", "data"),
  species = c("Homo_sapiens", "Mus_musculus"),
  batch_size = 10000,
  smooth_batch_size = 1000,
  parallelize_models = TRUE,
  parallelize_smoothing = TRUE,
  cores = 1,
  seed = 11,
  verbose = TRUE,
  ...
)

# Default S3 method
RunCytoTRACE(
  object,
  species = c("Homo_sapiens", "Mus_musculus"),
  batch_size = 10000,
  smooth_batch_size = 1000,
  parallelize_models = TRUE,
  parallelize_smoothing = TRUE,
  cores = 1,
  seed = 11,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  An object. This can be a Seurat object or a matrix-like object (genes
  as rows, cells as columns).

- ...:

  Additional arguments to be passed to
  [`CytoTRACE2::cytotrace2`](https://rdrr.io/pkg/CytoTRACE2/man/cytotrace2.html).

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used.

- layer:

  Which layer to use. Default is `"counts"`.

- species:

  The species of the input data. Currently supported values are
  `"human"` and `"mouse"`. Default is `"human"`.

- batch_size:

  The number of cells to process at once, including subsampling for KNN
  smoothing. No subsampling if `NULL`. Default is `10000` (recommended
  for input data size \> 10K cells).

- smooth_batch_size:

  The number of cells to subsample further within the batch_size for the
  smoothing by diffusion step of the pipeline. No subsampling if `NULL`.
  Default is `1000` (recommended for input data size \> 1K cells).

- parallelize_models:

  Whether to run the prediction function on models in parallel on
  multiple threads. Default is `TRUE`.

- parallelize_smoothing:

  Whether to run the smoothing function on subsamples in parallel on
  multiple threads. Default is `TRUE`.

- cores:

  The number of cores to use for parallelization with
  [foreach::foreach](https://rdrr.io/pkg/foreach/man/foreach.html).
  Default is `1`.

- seed:

  Random seed for reproducibility. Default is `11`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

When the input is a Seurat object, the function returns a Seurat object
with the following metadata columns added:

- `CytoTRACE2_Score`: The final predicted cellular potency score (0-1)

- `CytoTRACE2_Potency`: The final predicted cellular potency category
  (Differentiated, Unipotent, Oligopotent, Multipotent, Pluripotent,
  Totipotent)

- `CytoTRACE2_Relative`: The predicted relative order (normalized to
  0-1)

- `preKNN_CytoTRACE2_Score`: The potency score before KNN smoothing

- `preKNN_CytoTRACE2_Potency`: The potency category before KNN smoothing

When the input is a matrix or data.frame, the function returns a
data.frame with the same columns as above, with cell IDs as row names.

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
}
#> ℹ [2026-04-02 16:33:33] Start standard processing workflow...
#> ℹ [2026-04-02 16:33:34] Checking a list of <Seurat>...
#> ! [2026-04-02 16:33:34] Data 1/1 of the `srt_list` is "unknown"
#> ℹ [2026-04-02 16:33:34] Perform `NormalizeData()` with `normalization.method = 'LogNormalize'` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 16:33:36] Perform `Seurat::FindVariableFeatures()` on 1/1 of `srt_list`...
#> ℹ [2026-04-02 16:33:36] Use the separate HVF from `srt_list`
#> ℹ [2026-04-02 16:33:36] Number of available HVF: 2000
#> ℹ [2026-04-02 16:33:36] Finished check
#> ℹ [2026-04-02 16:33:37] Perform `Seurat::ScaleData()`
#> ℹ [2026-04-02 16:33:37] Perform pca linear dimension reduction
#> ℹ [2026-04-02 16:33:40] Use stored estimated dimensions 1:50 for Standardpca
#> ℹ [2026-04-02 16:33:41] Perform `Seurat::FindClusters()` with `cluster_algorithm = 'louvain'` and `cluster_resolution = 0.6`
#> ℹ [2026-04-02 16:33:41] Reorder clusters...
#> ℹ [2026-04-02 16:33:41] Skip `log1p()` because `layer = data` is not "counts"
#> ! [2026-04-02 16:33:41] <packageNotFoundError in loadNamespace(x): there is no package called ‘proxyC’>
#> ! [2026-04-02 16:33:41] Error when performing `Seurat::FindClusters()`. Skip it
#> ℹ [2026-04-02 16:33:41] Perform umap nonlinear dimension reduction
#> ℹ [2026-04-02 16:33:41] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ℹ [2026-04-02 16:33:44] Perform umap nonlinear dimension reduction using Standardpca (1:50)
#> ✔ [2026-04-02 16:33:47] Standard processing workflow completed
#> ◌ [2026-04-02 16:33:47] Running CytoTRACE2
#> ℹ [2026-04-02 16:33:47] Package CytoTRACE2 is not installed. Installing from GitHub...
#>  
#> → Will install 59 packages.
#> → All 59 packages (0 B) are cached.
#> + CytoTRACE2         1.1.0    [bld][cmp] (GitHub: 4a398a4)
#> + FNN                1.1.4.1  
#> + HiClimR            2.2.1     + ✔ libnetcdf-dev
#> + RANN               2.6.2    
#> + ROCR               1.0-12   
#> + RSpectra           0.16-2   
#> + RcppAnnoy          0.0.23   
#> + RcppHNSW           0.6.0    
#> + RcppParallel       5.1.11-2  + ✔ make
#> + RcppTOML           0.2.3    
#> + Rfast              2.1.5.2  
#> + Rtsne              0.17     
#> + Seurat             5.4.0    
#> + SeuratObject       5.3.0    
#> + bitops             1.0-9    
#> + caTools            1.18.3   
#> + cowplot            1.2.0    
#> + deldir             2.0-4    
#> + doParallel         1.0.17   
#> + dotCall64          1.2      
#> + fastDummies        1.7.5    
#> + fitdistrplus       1.2-6    
#> + foreach            1.5.2    
#> + future.apply       1.20.2   
#> + ggrepel            0.9.8    
#> + ggridges           0.5.7    
#> + goftest            1.2-3    
#> + gplots             3.3.0    
#> + gridExtra          2.3      
#> + gtools             3.9.5    
#> + here               1.0.2    
#> + ica                1.0-3    
#> + iterators          1.0.14   
#> + lmtest             0.9-40   
#> + miniUI             0.1.2    
#> + ncdf4              1.24      + ✔ libnetcdf-dev
#> + patchwork          1.3.2    
#> + pbapply            1.7-4    
#> + plyr               1.8.9    
#> + polyclip           1.10-7   
#> + progressr          0.19.0   
#> + reshape2           1.4.5    
#> + reticulate         1.45.0    + ✔ python3
#> + rprojroot          2.1.1    
#> + scattermore        1.2      
#> + sctransform        0.4.3    
#> + sp                 2.2-1    
#> + spam               2.11-3   
#> + spatstat.data      3.1-9    
#> + spatstat.explore   3.8-0    
#> + spatstat.geom      3.7-3    
#> + spatstat.random    3.4-5    
#> + spatstat.sparse    3.1-0    
#> + spatstat.univar    3.1-7    
#> + spatstat.utils     3.2-2    
#> + tensor             1.5.1    
#> + uwot               0.2.4    
#> + zigg               0.0.2    
#> + zoo                1.8-15   
#> ✔ All system requirements are already installed.
#>   
#> ℹ No downloads are needed, 59 pkgs are cached
#> ✔ Got ncdf4 1.24 (x86_64-pc-linux-gnu-ubuntu-24.04) (281.21 kB)
#> ✔ Got zigg 0.0.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (25.57 kB)
#> ✔ Got HiClimR 2.2.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (575.27 kB)
#> ✔ Got Rfast 2.1.5.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.99 MB)
#> ✔ Got CytoTRACE2 1.1.0 (source) (182.95 MB)
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
#> ℹ Executing `sudo sh -c apt-get -y install libnetcdf-dev make python3 libcurl4-openssl-dev libssl-dev cmake libuv1-dev zlib1g-dev libglpk-dev libxml2-dev pandoc libpng-dev libicu-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> libnetcdf-dev is already the newest version (1:4.9.2-5ubuntu4).
#> libnetcdf-dev set to manually installed.
#> make is already the newest version (4.3-4.1build2).
#> python3 is already the newest version (3.12.3-0ubuntu2.1).
#> libcurl4-openssl-dev is already the newest version (8.5.0-2ubuntu10.8).
#> libssl-dev is already the newest version (3.0.13-0ubuntu3.7).
#> cmake is already the newest version (3.28.3-1build7).
#> libuv1-dev is already the newest version (1.48.0-1.1build1).
#> zlib1g-dev is already the newest version (1:1.3.dfsg-3.1ubuntu2.1).
#> libglpk-dev is already the newest version (5.0-1build2).
#> libxml2-dev is already the newest version (2.9.14+dfsg-1.3ubuntu3.7).
#> pandoc is already the newest version (3.1.3+ds-2).
#> libpng-dev is already the newest version (1.6.43-5ubuntu0.5).
#> libicu-dev is already the newest version (74.2-1ubuntu3.1).
#> 0 upgraded, 0 newly installed, 0 to remove and 82 not upgraded.
#> ✔ Installed bitops 1.0-9  (47ms)
#> ✔ Installed caTools 1.18.3  (68ms)
#> ✔ Installed cowplot 1.2.0  (111ms)
#> ✔ Installed deldir 2.0-4  (132ms)
#> ✔ Installed doParallel 1.0.17  (91ms)
#> ✔ Installed dotCall64 1.2  (96ms)
#> ✔ Installed fastDummies 1.7.5  (63ms)
#> ✔ Installed fitdistrplus 1.2-6  (68ms)
#> ✔ Installed FNN 1.1.4.1  (62ms)
#> ✔ Installed foreach 1.5.2  (60ms)
#> ✔ Installed future.apply 1.20.2  (62ms)
#> ✔ Installed ggrepel 0.9.8  (103ms)
#> ✔ Installed ggridges 0.5.7  (103ms)
#> ✔ Installed goftest 1.2-3  (61ms)
#> ✔ Installed gridExtra 2.3  (1s)
#> ✔ Installed gplots 3.3.0  (1.1s)
#> ✔ Installed gtools 3.9.5  (64ms)
#> ✔ Installed here 1.0.2  (62ms)
#> ✔ Installed HiClimR 2.2.1  (64ms)
#> ✔ Installed ica 1.0-3  (62ms)
#> ✔ Installed iterators 1.0.14  (94ms)
#> ✔ Installed lmtest 0.9-40  (61ms)
#> ✔ Installed miniUI 0.1.2  (63ms)
#> ✔ Installed ncdf4 1.24  (59ms)
#> ✔ Installed patchwork 1.3.2  (61ms)
#> ✔ Installed pbapply 1.7-4  (61ms)
#> ✔ Installed plyr 1.8.9  (98ms)
#> ✔ Installed polyclip 1.10-7  (97ms)
#> ✔ Installed progressr 0.19.0  (62ms)
#> ✔ Installed RANN 2.6.2  (62ms)
#> ✔ Installed RcppAnnoy 0.0.23  (62ms)
#> ✔ Installed RcppHNSW 0.6.0  (62ms)
#> ✔ Installed RcppParallel 5.1.11-2  (74ms)
#> ✔ Installed RcppTOML 0.2.3  (100ms)
#> ✔ Installed reshape2 1.4.5  (65ms)
#> ✔ Installed reticulate 1.45.0  (68ms)
#> ✔ Installed Rfast 2.1.5.2  (68ms)
#> ✔ Installed ROCR 1.0-12  (65ms)
#> ✔ Installed rprojroot 2.1.1  (64ms)
#> ✔ Installed RSpectra 0.16-2  (101ms)
#> ✔ Installed Rtsne 0.17  (97ms)
#> ✔ Installed scattermore 1.2  (59ms)
#> ✔ Installed sctransform 0.4.3  (62ms)
#> ✔ Installed Seurat 5.4.0  (71ms)
#> ✔ Installed SeuratObject 5.3.0  (73ms)
#> ✔ Installed sp 2.2-1  (1.1s)
#> ✔ Installed spam 2.11-3  (1.1s)
#> ✔ Installed spatstat.data 3.1-9  (72ms)
#> ✔ Installed spatstat.explore 3.8-0  (125ms)
#> ✔ Installed spatstat.geom 3.7-3  (275ms)
#> ✔ Installed spatstat.random 3.4-5  (256ms)
#> ✔ Installed spatstat.sparse 3.1-0  (63ms)
#> ✔ Installed spatstat.univar 3.1-7  (63ms)
#> ✔ Installed spatstat.utils 3.2-2  (63ms)
#> ✔ Installed tensor 1.5.1  (97ms)
#> ✔ Installed uwot 0.2.4  (66ms)
#> ✔ Installed zigg 0.0.2  (62ms)
#> ✔ Installed zoo 1.8-15  (48ms)
#> ℹ Packaging CytoTRACE2 1.1.0
#> ✔ Packaged CytoTRACE2 1.1.0 (4.7s)
#> ℹ Building CytoTRACE2 1.1.0
#> ✔ Built CytoTRACE2 1.1.0 (9.7s)
#> ✔ Installed CytoTRACE2 1.1.0 (github::digitalcytometry/cytotrace2@4a398a4) (1.1s)
#> ✔ 1 pkg + 143 deps: kept 83, added 59, dld 5 (NA B) [51.4s]
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
#> ✔ [2026-04-02 16:35:59] CytoTRACE2 computed successfully
```
