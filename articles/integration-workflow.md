# Data integration workflow

## What you will do

- **Question:** does correction reduce technology effects while
  retaining cell types?
- **Input:** All 1,600 cells and 12,932 supplied genes in `panc8_sub`;
  `tech` identifies five technologies and `celltype` supplies evaluation
  labels.
- **Output:** Uncorrected and Harmony reductions on the same cells,
  marker checks and a reloadable integrated object.
- **How to run:** start a fresh R session and run the chunks in order.
  Printed tables and figures below are generated from that run;
  code-only alternatives are identified explicitly.

Use integration when technical structure prevents comparison across
samples or platforms, and retain an unintegrated baseline. If every
condition occurs in only one batch, batch and biology are confounded: no
attractive embedding can resolve that experimental-design problem.

## Prepare the required Harmony source

SCOP’s
[`RunHarmony2()`](https://mengxu98.github.io/scop/reference/RunHarmony2.md)
checks the official GitHub package origin. A CRAN installation with the
same package name does not satisfy that check and can trigger an
installation attempt in the middle of an analysis. Install the official
GitHub source before starting this page, then restart R:

``` r

remotes::install_github(
  "immunogenomics/harmony",
  upgrade = "never"
)
```

The tutorial’s dependency preflight is read-only. It fails with a clear
setup requirement instead of allowing an implicit installation during
fitting.

## Prepare a Merged Object

The example uses the bundled `panc8_sub` object, a subset of human
pancreas datasets. The `tech` column records the batch or technology
label. The `celltype` column is used as known biology for diagnostics,
not as a target to force during integration.

``` r

library(scop)
set.seed(11)
harmony_source <- utils::packageDescription("harmony")
stopifnot(
  identical(harmony_source$RemoteUsername, "immunogenomics"),
  identical(harmony_source$RemoteRepo, "harmony")
)
data.frame(
  package = "harmony", version = harmony_source$Version,
  source = paste(harmony_source$RemoteUsername, harmony_source$RemoteRepo, sep = "/"),
  revision = harmony_source$RemoteSha
)
#>   package version                 source
#> 1 harmony   2.1.0 immunogenomics/harmony
#>                                   revision
#> 1 7e21bd6629b628099735395f57a47d11f873b2ee

data(panc8_sub)
input_cells <- colnames(panc8_sub)
stopifnot(!anyNA(panc8_sub$tech), !anyNA(panc8_sub$celltype))
data.frame(genes = nrow(panc8_sub), cells = ncol(panc8_sub))
#>   genes cells
#> 1 12932  1600
table(panc8_sub$tech)
#> 
#>     celseq    celseq2 fluidigmc1     indrop  smartseq2 
#>        200        200        200        800        200
table(panc8_sub$celltype)
#> 
#>             acinar activated-stellate              alpha               beta 
#>                164                 50                483                428 
#>              delta             ductal        endothelial            epsilon 
#>                111                237                 26                  3 
#>              gamma         macrophage               mast quiescent-stellate 
#>                 67                  9                  7                 14 
#>            schwann 
#>                  1
```

Always inspect the unintegrated embedding or PCA first. This gives you a
baseline for deciding whether integration removes technical structure
without erasing cell-type structure.

``` r

panc8_sub <- RunStandardWorkflow(panc8_sub,
  assay = "RNA", seed = 11, cores = 1,
  verbose = FALSE
)

CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype"),
  reduction = "UMAP",
  theme_use = "theme_blank"
)
```

![](integration-workflow_files/figure-html/unintegrated-check-1.png)

## Run an Integration Method

Harmony is a common first pass when the object is already a merged
Seurat object with a batch column.

``` r

panc8_harmony <- RunIntegration(
  srt_merge = panc8_sub,
  batch = "tech",
  integration_methods = "Harmony",
  nonlinear_reduction_dims = 2,
  seed = 11, verbose = FALSE
)
#> Warning: Different features in new layer data than already exists for
#> scale.data

stopifnot(identical(colnames(panc8_harmony), input_cells))
SeuratObject::Reductions(panc8_harmony)
#> [1] "Standardpca"       "StandardpcaUMAP2D" "StandardUMAP2D"   
#> [4] "Harmonypca"        "Harmony"           "HarmonyUMAP2D"    
#> [7] "HarmonypcaUMAP2D"
```

Use the method-specific integrated UMAP for plotting. The exact
reduction name depends on the selected method; for Harmony, the examples
use `HarmonyUMAP2D`. When comparing methods, report the reduction name
explicitly. A plot labeled only as UMAP is ambiguous once an object
contains uncorrected and integrated embeddings.

``` r

CellDimPlot(
  panc8_harmony,
  group.by = c("tech", "celltype"),
  reduction = "HarmonyUMAP2D",
  theme_use = "theme_blank"
)
```

![](integration-workflow_files/figure-html/plot-harmony-1.png)

## Choose a Method

Use method choice to match the data and dependency budget:

- `Harmony`, `fastMNN`, `MNN`, `ComBat`, and `Seurat` methods are
  practical R-side choices for many scRNA-seq batch problems.
- `scVI`, `BBKNN`, `Scanorama`, `GLUE`, and `MultiMAP` require
  Python-side dependencies but can be useful for large datasets or
  multimodal settings.
- `WNN` and cross-modality methods are better suited to multimodal RNA +
  ATAC or paired-modality objects than to ordinary batch correction.
- `Uncorrected` is useful as a baseline and should be kept in reports
  when showing that integration changed the object.

The following alternatives are recipes, not additional measured results
in this page. To compare executed methods, use the [integration
benchmark](https://mengxu98.github.io/scop/articles/integration-benchmark-workflow.md).
Other methods follow the same entry-point pattern, but may require
additional runtime dependencies or matrix backends:

``` r

panc8_rpca <- RunIntegration(
  srt_merge = panc8_sub,
  batch = "tech",
  integration_methods = "RPCA"
)

panc8_fastmnn <- RunIntegration(
  srt_merge = panc8_sub,
  batch = "tech",
  integration_methods = "fastMNN"
)
```

## Inspect Integration Quality

Do not judge integration from batch mixing alone. Also check known
labels, marker genes, and whether rare cell types remain separated.

``` r

CellDimPlot(
  panc8_harmony,
  group.by = "tech",
  reduction = "HarmonyUMAP2D",
  theme_use = "theme_blank"
)
```

![](integration-workflow_files/figure-html/inspect-quality-1.png)

``` r


CellDimPlot(
  panc8_harmony,
  group.by = "celltype",
  reduction = "HarmonyUMAP2D",
  label = TRUE,
  label_repel = TRUE, label.size = 3,
  theme_use = "theme_blank"
)
```

![](integration-workflow_files/figure-html/inspect-quality-2.png)

``` r


FeatureDimPlot(
  panc8_harmony,
  features = c("GCG", "SST", "KRT19"),
  reduction = "HarmonyUMAP2D"
)
```

![](integration-workflow_files/figure-html/inspect-quality-3.png)

When comparing methods, keep the input object fixed and compare the same
metadata labels, marker genes, and diagnostic plots across methods. The
unintegrated baseline is useful because it shows what the integration
changed. Good integration mixes batches within the same biological state
while keeping distinct cell types and marker programs separated.
Overcorrection is visible when expected marker-defined groups collapse
together.

INS is absent from the supplied feature-filtered object, so this marker
check uses measured GCG, SST and KRT19. Missing input genes must not be
interpreted as measured zero expression.

## Save, reload, and check the hand-off

The code below saves to a temporary directory. Choose a persistent
output folder when saving your analysis. Save the object **and**
plain-text metadata; a figure alone cannot reproduce the analysis. No
cells are sampled for this export.

``` r

output_dir <- file.path(tempdir(), "integration-workflow")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
result_file <- file.path(output_dir, "analysis.rds")
saveRDS(panc8_harmony, result_file)
utils::write.csv(panc8_harmony@meta.data, file.path(output_dir, "cell-metadata.csv"))
reloaded <- readRDS(result_file)
stopifnot(
  identical(colnames(reloaded), colnames(panc8_harmony)),
  identical(reloaded@meta.data, panc8_harmony@meta.data),
  identical(
    SeuratObject::Reductions(reloaded),
    SeuratObject::Reductions(panc8_harmony)
  )
)
data.frame(
  file = basename(list.files(output_dir, full.names = TRUE)),
  bytes = file.info(list.files(output_dir, full.names = TRUE))$size
)
#>                file    bytes
#> 1      analysis.rds 48008781
#> 2 cell-metadata.csv   186170
```

## Use your own data and recover common errors

- Replace the `data(...)` input with `readRDS("your-object.rds")`, or
  create a Seurat object from a gene-by-cell raw-count matrix. Check
  that cell IDs are unique and metadata row names equal the expression
  column names.
- Set the assay and grouping columns explicitly. Species, gene
  identifiers and reference labels must match your experiment; do not
  copy pancreas labels to unrelated tissue.
- If a reduction or metadata column is missing, inspect
  `SeuratObject::Reductions(object)` and `colnames(object@meta.data)`,
  then run its producer chunk. A plot does not create the result it
  needs.
- If a dependency is unavailable, install it in the analysis
  environment, restart R and rerun from the input chunk. A skipped or
  failed method is not evidence that it worked. Preserve warnings and
  the exact error in your report.

## Analysis provenance

The receipt records the R seed. Python methods retain their documented
backend seeds and parameters; setting the R seed alone does not control
every backend. Package versions, BLAS libraries and backend changes can
alter numerical results. Keep this receipt with the exported object.

``` r

source_revision <- Sys.getenv("GITHUB_SHA")
if (!nzchar(source_revision) && nzchar(Sys.which("git"))) {
  source_revision <- tryCatch(system2("git", c("rev-parse", "HEAD"),
    stdout = TRUE,
    stderr = FALSE
  )[1], error = function(e) NA_character_)
}
article_source <- knitr::current_input(dir = TRUE)
build_receipt <- data.frame(
  article = "integration-workflow",
  built_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  scop = as.character(utils::packageVersion("scop")),
  R = as.character(getRversion()),
  source_revision = source_revision,
  source_has_local_changes = if (nzchar(Sys.which("git"))) {
    length(system2("git", c("status", "--porcelain"), stdout = TRUE, stderr = FALSE)) > 0
  } else {
    NA
  },
  article_md5 = unname(tools::md5sum(article_source)),
  R_seed = 11,
  cells_retained = ncol(panc8_harmony),
  genes_in_active_assay = nrow(panc8_harmony),
  stringsAsFactors = FALSE
)
build_receipt
#>                article               built_utc  scop     R
#> 1 integration-workflow 2026-11-11 12:37:25 UTC 0.9.2 4.6.2
#>                            source_revision source_has_local_changes
#> 1 fce961e0e303ffc7bcb7aa68b6f72914fc1474e6                    FALSE
#>                        article_md5 R_seed cells_retained genes_in_active_assay
#> 1 258fa719de827f5739c0606e0dfd1e7c     11           1600                 12932
utils::write.csv(build_receipt, file.path(output_dir, "build-receipt.csv"),
  row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session-info.txt"))
```

Show the complete R session used for this page

``` r

sessionInfo()
#> R version 4.6.2 (2026-10-29)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] scop_0.9.2
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_1.2.1    dplyr_1.2.1         farver_2.1.2       
#>  [4] S7_0.2.2            fastmap_1.2.0       RANN_2.6.3         
#>  [7] digest_0.6.39       dotCall64_1.2       lifecycle_1.0.5    
#> [10] cluster_2.1.8.3     SeuratObject_5.4.0  magrittr_2.0.5     
#> [13] compiler_4.6.2      rlang_1.3.0         sass_0.4.10        
#> [16] tools_4.6.2         igraph_2.3.4        yaml_2.3.12        
#> [19] knitr_1.52          labeling_0.4.3      htmlwidgets_1.6.4  
#> [22] sp_2.2-4            plyr_1.8.9          RColorBrewer_1.1-3 
#> [25] KernSmooth_2.23-27  Rtsne_0.17          withr_3.0.3        
#> [28] desc_1.4.3          grid_4.6.2          RcppHNSW_0.7.0     
#> [31] fastDummies_1.7.6   future_1.82.0       progressr_1.1.0    
#> [34] ggplot2_4.0.3       globals_0.20.0      scales_1.4.0       
#> [37] MASS_7.3-66         ggridges_0.5.7      dichromat_2.0-1    
#> [40] cli_3.6.6           rmarkdown_2.32      ragg_1.5.2         
#> [43] thisutils_0.5.3     generics_0.1.4      Seurat_6.0.0       
#> [46] otel_0.3.0          ica_1.0-3           future.apply_1.21.0
#> [49] RSpectra_0.16-2     sctransform_0.4.3   httr_1.4.9         
#> [52] reshape2_1.4.5      pbapply_1.7-5       cachem_1.1.0       
#> [55] stringr_1.6.0       parallel_4.6.2      proxyC_0.6.0       
#> [58] matrixStats_1.5.0   vctrs_0.7.3         Matrix_1.7-6       
#> [61] jsonlite_2.0.0      patchwork_1.3.2     thisplot_0.4.7     
#> [64] ggrepel_0.9.8       irlba_2.4.1         scattermore_1.2    
#> [67] listenv_1.1.0       systemfonts_1.3.2   harmony_2.1.0      
#> [70] jquerylib_0.1.4     pls_2.9-0           glue_1.8.1         
#> [73] parallelly_1.50.0   spam_2.11-4         pkgdown_2.2.1      
#> [76] codetools_0.2-20    uwot_0.2.5          cowplot_1.2.0      
#> [79] stringi_1.8.9       RcppAnnoy_0.0.23    gtable_0.3.6       
#> [82] tibble_3.3.1        pillar_1.11.1       htmltools_0.5.9    
#> [85] R6_2.6.1            textshaping_1.0.5   evaluate_1.0.5     
#> [88] lattice_0.23-1      png_0.1-9           RhpcBLASctl_0.23-42
#> [91] bslib_0.12.0        Rcpp_1.1.2          gridExtra_2.3.1    
#> [94] xfun_0.62           fs_2.1.0            pkgconfig_2.0.3
```

## What to Report

A compact integration report should include:

- the batch column and biological annotation columns used for
  evaluation;
- the integration method and major parameters;
- the integrated reduction names used for plots;
- before and after plots colored by batch and biological labels;
- marker-gene checks for expected cell types;
- any method-specific dependencies or environment setup needed to
  reproduce the run.

Integration is a preprocessing decision, not a biological result by
itself. Use it only when the integrated representation improves
downstream interpretation without hiding real biological differences.
