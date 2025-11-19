# Check and install R packages

Check and install R packages

## Usage

``` r
check_r(packages, lib = .libPaths()[1], force = FALSE, verbose = TRUE)
```

## Arguments

- packages:

  Package to be installed. Package source can be CRAN, Bioconductor or
  Github. By default, the package name is extracted according to the
  `packages` parameter.

- lib:

  The location of the library directories where to install the packages.

- force:

  Whether to force the installation of packages. Default is `FALSE`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Examples

``` r
check_r(c("Seurat", "reticulate"))
#> ◌ [2025-11-19 15:14:03] Installing: Seurat...
#>  
#> → Will install 26 packages.
#> → All 26 packages (0 B) are cached.
#> + ROCR               1.0-11 
#> + Seurat             5.3.1  
#> + SeuratObject       5.2.0  
#> + caTools            1.18.3 
#> + deldir             2.0-4  
#> + dotCall64          1.2    
#> + fastDummies        1.7.5  
#> + fitdistrplus       1.2-4  
#> + goftest            1.2-3  
#> + gplots             3.2.0  
#> + gtools             3.9.5  
#> + ica                1.0-3  
#> + leidenbase         0.1.35 
#> + miniUI             0.1.2  
#> + progressr          0.18.0 
#> + scattermore        1.2    
#> + sctransform        0.4.2  
#> + spam               2.11-1 
#> + spatstat.data      3.1-9  
#> + spatstat.explore   3.5-3  
#> + spatstat.geom      3.6-0  
#> + spatstat.random    3.4-2  
#> + spatstat.sparse    3.1-0  
#> + spatstat.univar    3.1-5  
#> + spatstat.utils     3.2-0  
#> + tensor             1.5.1  
#> ✔ All system requirements are already installed.
#>   
#> ℹ No downloads are needed, 26 pkgs are cached
#> ✔ Got dotCall64 1.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (31.26 kB)
#> ✔ Got goftest 1.2-3 (x86_64-pc-linux-gnu-ubuntu-24.04) (58.64 kB)
#> ✔ Got caTools 1.18.3 (x86_64-pc-linux-gnu-ubuntu-24.04) (205.11 kB)
#> ✔ Got deldir 2.0-4 (x86_64-pc-linux-gnu-ubuntu-24.04) (278.91 kB)
#> ✔ Got sctransform 0.4.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (565.94 kB)
#> ✔ Got gtools 3.9.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (361.11 kB)
#> ✔ Got progressr 0.18.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (412.00 kB)
#> ✔ Got spatstat.sparse 3.1-0 (x86_64-pc-linux-gnu-ubuntu-24.04) (213.31 kB)
#> ✔ Got fitdistrplus 1.2-4 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.96 MB)
#> ✔ Got fastDummies 1.7.5 (x86_64-pc-linux-gnu-ubuntu-24.04) (41.76 kB)
#> ✔ Got ica 1.0-3 (x86_64-pc-linux-gnu-ubuntu-24.04) (85.54 kB)
#> ✔ Got tensor 1.5.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (13.96 kB)
#> ✔ Got gplots 3.2.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (496.91 kB)
#> ✔ Got spatstat.random 3.4-2 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.22 MB)
#> ✔ Got ROCR 1.0-11 (x86_64-pc-linux-gnu-ubuntu-24.04) (465.58 kB)
#> ✔ Got scattermore 1.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (339.10 kB)
#> ✔ Got spatstat.utils 3.2-0 (x86_64-pc-linux-gnu-ubuntu-24.04) (401.56 kB)
#> ✔ Got SeuratObject 5.2.0 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.08 MB)
#> ✔ Got leidenbase 0.1.35 (x86_64-pc-linux-gnu-ubuntu-24.04) (1.33 MB)
#> ✔ Got miniUI 0.1.2 (x86_64-pc-linux-gnu-ubuntu-24.04) (35.20 kB)
#> ✔ Got spam 2.11-1 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.75 MB)
#> ✔ Got Seurat 5.3.1 (x86_64-pc-linux-gnu-ubuntu-24.04) (2.70 MB)
#> ✔ Got spatstat.explore 3.5-3 (x86_64-pc-linux-gnu-ubuntu-24.04) (3.51 MB)
#> ✔ Got spatstat.geom 3.6-0 (x86_64-pc-linux-gnu-ubuntu-24.04) (4.18 MB)
#> ✔ Got spatstat.data 3.1-9 (x86_64-pc-linux-gnu-ubuntu-24.04) (4.17 MB)
#> ℹ Installing system requirements
#> ℹ Executing `sudo sh -c apt-get -y update`
#> Get:1 file:/etc/apt/apt-mirrors.txt Mirrorlist [144 B]
#> Hit:2 http://azure.archive.ubuntu.com/ubuntu noble InRelease
#> Hit:6 https://packages.microsoft.com/repos/azure-cli noble InRelease
#> Hit:7 https://packages.microsoft.com/ubuntu/24.04/prod noble InRelease
#> Hit:3 http://azure.archive.ubuntu.com/ubuntu noble-updates InRelease
#> Hit:4 http://azure.archive.ubuntu.com/ubuntu noble-backports InRelease
#> Hit:5 http://azure.archive.ubuntu.com/ubuntu noble-security InRelease
#> Reading package lists...
#> ℹ Executing `sudo sh -c apt-get -y install libcurl4-openssl-dev libssl-dev make zlib1g-dev libglpk-dev libxml2-dev pandoc libpng-dev python3 libicu-dev`
#> Reading package lists...
#> Building dependency tree...
#> Reading state information...
#> libcurl4-openssl-dev is already the newest version (8.5.0-2ubuntu10.6).
#> libssl-dev is already the newest version (3.0.13-0ubuntu3.6).
#> make is already the newest version (4.3-4.1build2).
#> zlib1g-dev is already the newest version (1:1.3.dfsg-3.1ubuntu2.1).
#> libglpk-dev is already the newest version (5.0-1build2).
#> libxml2-dev is already the newest version (2.9.14+dfsg-1.3ubuntu3.6).
#> pandoc is already the newest version (3.1.3+ds-2).
#> libpng-dev is already the newest version (1.6.43-5build1).
#> python3 is already the newest version (3.12.3-0ubuntu2.1).
#> libicu-dev is already the newest version (74.2-1ubuntu3.1).
#> 0 upgraded, 0 newly installed, 0 to remove and 31 not upgraded.
#> ✔ Installed caTools 1.18.3  (58ms)
#> ✔ Installed deldir 2.0-4  (73ms)
#> ✔ Installed dotCall64 1.2  (97ms)
#> ✔ Installed fastDummies 1.7.5  (160ms)
#> ✔ Installed fitdistrplus 1.2-4  (66ms)
#> ✔ Installed goftest 1.2-3  (64ms)
#> ✔ Installed gplots 3.2.0  (65ms)
#> ✔ Installed gtools 3.9.5  (68ms)
#> ✔ Installed ica 1.0-3  (65ms)
#> ✔ Installed leidenbase 0.1.35  (62ms)
#> ✔ Installed miniUI 0.1.2  (148ms)
#> ✔ Installed progressr 0.18.0  (64ms)
#> ✔ Installed ROCR 1.0-11  (67ms)
#> ✔ Installed scattermore 1.2  (65ms)
#> ✔ Installed sctransform 0.4.2  (66ms)
#> ✔ Installed Seurat 5.3.1  (72ms)
#> ✔ Installed SeuratObject 5.2.0  (73ms)
#> ✔ Installed spam 2.11-1  (112ms)
#> ✔ Installed spatstat.data 3.1-9  (81ms)
#> ✔ Installed spatstat.explore 3.5-3  (75ms)
#> ✔ Installed spatstat.geom 3.6-0  (69ms)
#> ✔ Installed spatstat.random 3.4-2  (67ms)
#> ✔ Installed spatstat.sparse 3.1-0  (67ms)
#> ✔ Installed spatstat.univar 3.1-5  (64ms)
#> ✔ Installed spatstat.utils 3.2-0  (63ms)
#> ✔ Installed tensor 1.5.1  (88ms)
#> ✔ 1 pkg + 135 deps: kept 109, added 26, dld 25 (28.90 MB) [9.6s]
#> ✔ [2025-11-19 15:14:13] Seurat and reticulate installed successfully
```
