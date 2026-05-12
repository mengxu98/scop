# Convert homologous gene symbols in expression objects

Convert feature names between species with
[GeneConvert](https://mengxu98.github.io/scop/reference/GeneConvert.md)
and collapse duplicated target homologs by summing expression values.
The Seurat method rebuilds the selected assay from the converted counts
matrix and keeps cell metadata and spatial images when present.

## Usage

``` r
ConvertHomologs(
  object,
  species_from,
  species_to,
  geneID_from_IDtype = "symbol",
  geneID_to_IDtype = "symbol",
  assay = NULL,
  layer = "counts",
  multi_mapping = c("first"),
  keep_unmapped = FALSE,
  collapse_fun = c("sum"),
  Ensembl_version = NULL,
  biomart = NULL,
  mirror = NULL,
  max_tries = 5,
  verbose = TRUE
)

# S3 method for class 'Seurat'
ConvertHomologs(
  object,
  species_from,
  species_to,
  geneID_from_IDtype = "symbol",
  geneID_to_IDtype = "symbol",
  assay = NULL,
  layer = "counts",
  multi_mapping = c("first"),
  keep_unmapped = FALSE,
  collapse_fun = c("sum"),
  Ensembl_version = NULL,
  biomart = NULL,
  mirror = NULL,
  max_tries = 5,
  verbose = TRUE
)

# S3 method for class 'matrix'
ConvertHomologs(
  object,
  species_from,
  species_to,
  geneID_from_IDtype = "symbol",
  geneID_to_IDtype = "symbol",
  assay = NULL,
  layer = "counts",
  multi_mapping = c("first"),
  keep_unmapped = FALSE,
  collapse_fun = c("sum"),
  Ensembl_version = NULL,
  biomart = NULL,
  mirror = NULL,
  max_tries = 5,
  verbose = TRUE
)

# S3 method for class 'Matrix'
ConvertHomologs(
  object,
  species_from,
  species_to,
  geneID_from_IDtype = "symbol",
  geneID_to_IDtype = "symbol",
  assay = NULL,
  layer = "counts",
  multi_mapping = c("first"),
  keep_unmapped = FALSE,
  collapse_fun = c("sum"),
  Ensembl_version = NULL,
  biomart = NULL,
  mirror = NULL,
  max_tries = 5,
  verbose = TRUE
)

# Default S3 method
ConvertHomologs(
  object,
  species_from,
  species_to,
  geneID_from_IDtype = "symbol",
  geneID_to_IDtype = "symbol",
  assay = NULL,
  layer = "counts",
  multi_mapping = c("first"),
  keep_unmapped = FALSE,
  collapse_fun = c("sum"),
  Ensembl_version = NULL,
  biomart = NULL,
  mirror = NULL,
  max_tries = 5,
  verbose = TRUE
)
```

## Arguments

- object:

  A `Seurat` object or a gene-by-cell matrix.

- species_from:

  Latin names for animals of the input geneID. e.g. `"Homo_sapiens"`,
  `"Mus_musculus"`.

- species_to:

  Latin names for animals of the output geneID. e.g. `"Homo_sapiens"`,
  `"Mus_musculus"`.

- geneID_from_IDtype:

  Gene ID type of the input `geneID`. e.g. `"symbol"`, `"ensembl_id"`,
  `"entrez_id"`

- geneID_to_IDtype:

  Gene ID type(s) to convert to. e.g. `"symbol"`, `"ensembl_id"`,
  `"entrez_id"`.

- assay:

  Assay to convert when `object` is a `Seurat` object. If `NULL`, the
  default assay is used.

- layer:

  Assay layer used for conversion. Default `"counts"`.

- multi_mapping:

  How to handle source genes mapped to multiple target homologs.
  `"first"` keeps the first target homolog for each source gene.

- keep_unmapped:

  Whether to keep unmapped source genes with their original names.

- collapse_fun:

  Function used to collapse duplicated target homologs. Currently only
  `"sum"` is supported.

- Ensembl_version:

  An integer specifying the Ensembl version. Default is `NULL`. If
  `NULL`, the latest version will be used.

- biomart:

  The name of the BioMart database that you want to connect to. Possible
  options include `"ensembl"`, `"protists_mart"`, `"fungi_mart"`, and
  `"plants_mart"`.

- mirror:

  Specify an Ensembl mirror to connect to. The valid options here are
  `"www"`, `"uswest"`, `"useast"`, `"asia"`.

- max_tries:

  The maximum number of attempts to connect with the BioMart service.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A converted object of the same high-level type as `object`. The mapping
table is stored in `@tools$ConvertHomologs` for Seurat objects and in
the `"ConvertHomologs"` attribute for matrix inputs.

## See also

[AnnotateFeatures](https://mengxu98.github.io/scop/reference/AnnotateFeatures.md),
ConvertHomologs\]

## Examples

``` r
data(pancreas_sub)
pancreas_human <- ConvertHomologs(
  pancreas_sub,
  species_from = "Mus_musculus",
  species_to = "Homo_sapiens"
)
#> ℹ [2026-05-12 03:41:36] Connect to the Ensembl archives...
#> ! [2026-05-12 03:41:47] <error/httr2_failure>
#> !                       Error in `req_perform()`:
#> !                       ! Failed to perform HTTP request.
#> !                       Caused by error in `curl::curl_fetch_memory()`:
#> !                       ! Timeout was reached [www.ensembl.org]:
#> !                       Send failure: Broken pipe
#> !                       ---
#> !                       Backtrace:
#> !                            ▆
#> !                         1. ├─base::tryCatch(...)
#> !                         2. │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#> !                         3. │   ├─base (local) tryCatchOne(...)
#> !                         4. │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#> !                         5. │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#> !                         6. │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#> !                         7. │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#> !                         8. ├─base::withCallingHandlers(...)
#> !                         9. ├─base::saveRDS(...)
#> !                        10. ├─base::do.call(...)
#> !                        11. ├─base (local) `<fn>`(...)
#> !                        12. └─global `<fn>`(...)
#> !                        13.   └─pkgdown::build_site(...)
#> !                        14.     └─pkgdown:::build_site_local(...)
#> !                        15.       └─pkgdown::build_reference(...)
#> !                        16.         ├─pkgdown:::unwrap_purrr_error(...)
#> !                        17.         │ └─base::withCallingHandlers(...)
#> !                        18.         └─purrr::map(...)
#> !                        19.           └─purrr:::map_("list", .x, .f, ..., .progress = .progress)
#> !                        20.             ├─purrr:::with_indexed_errors(...)
#> !                        21.             │ └─base::withCallingHandlers(...)
#> !                        22.             ├─purrr:::call_with_cleanup(...)
#> !                        23.             └─pkgdown (local) .f(.x[[i]], ...)
#> !                        24.               ├─base::withCallingHandlers(...)
#> !                        25.               └─pkgdown:::data_reference_topic(...)
#> !                        26.                 └─pkgdown:::run_examples(...)
#> !                        27.                   └─pkgdown:::highlight_examples(code, topic, env = env)
#> !                        28.                     └─downlit::evaluate_and_highlight(...)
#> !                        29.                       └─evaluate::evaluate(code, child_env(env), new_device = TRUE, output_handler = output_handler)
#> !                        30.                         ├─base::withRestarts(...)
#> !                        31.                         │ └─base (local) withRestartList(expr, restarts)
#> !                        32.                         │   ├─base (local) withOneRestart(withRestartList(expr, restarts[-nr]), restarts[[nr]])
#> !                        33.                         │   │ └─base (local) doWithOneRestart(return(expr), restart)
#> !                        34.                         │   └─base (local) withRestartList(expr, restarts[-nr])
#> !                        35.                         │     └─base (local) withOneRestart(expr, restarts[[1L]])
#> !                        36.                         │       └─base (local) doWithOneRestart(return(expr), restart)
#> !                        37.                         ├─evaluate:::with_handlers(...)
#> !                        38.                         │ ├─base::eval(call)
#> !                        39.                         │ │ └─base::eval(call)
#> !                        40.                         │ └─base::withCallingHandlers(...)
#> !                        41.                         ├─base::withVisible(eval(expr, envir))
#> !                        42.                         └─base::eval(expr, envir)
#> !                        43.                           └─base::eval(expr, envir)
#> !                        44.                             └─scop::ConvertHomologs(...)
#> !                        45.                               └─scop:::ConvertHomologs.Seurat(...)
#> !                        46.                                 └─scop:::ConvertHomologs.default(...)
#> !                        47.                                   └─scop::GeneConvert(...)
#> !                        48.                                     ├─thisutils::try_get(...)
#> !                        49.                                     │ ├─base::tryCatch(...)
#> !                        50.                                     │ │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#> !                        51.                                     │ │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#> !                        52.                                     │ │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#> !                        53.                                     │ └─base::eval.parent(substitute(expr))
#> !                        54.                                     │   └─base::eval(expr, p)
#> !                        55.                                     │     └─base::eval(expr, p)
#> !                        56.                                     └─biomaRt::listEnsemblArchives()
#> !                        57.                                       └─biomaRt:::.listEnsemblArchives(http_config = list())
#> !                        58.                                         └─biomaRt:::.checkArchiveList(http_config)
#> !                        59.                                           └─biomaRt:::.getArchiveList(http_config = http_config)
#> !                        60.                                             └─httr2::req_perform(html_request)
#> ! [2026-05-12 03:41:47] Get errors when connecting with EnsemblArchives...
#> ! [2026-05-12 03:41:48] Retrying...
#> ! [2026-05-12 03:41:58] <error/httr2_failure>
#> !                       Error in `req_perform()`:
#> !                       ! Failed to perform HTTP request.
#> !                       Caused by error in `curl::curl_fetch_memory()`:
#> !                       ! Timeout was reached [www.ensembl.org]:
#> !                       Operation timed out after 10003 milliseconds with 0 bytes received
#> !                       ---
#> !                       Backtrace:
#> !                            ▆
#> !                         1. ├─base::tryCatch(...)
#> !                         2. │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#> !                         3. │   ├─base (local) tryCatchOne(...)
#> !                         4. │   │ └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#> !                         5. │   └─base (local) tryCatchList(expr, names[-nh], parentenv, handlers[-nh])
#> !                         6. │     └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#> !                         7. │       └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#> !                         8. ├─base::withCallingHandlers(...)
#> !                         9. ├─base::saveRDS(...)
#> !                        10. ├─base::do.call(...)
#> !                        11. ├─base (local) `<fn>`(...)
#> !                        12. └─global `<fn>`(...)
#> !                        13.   └─pkgdown::build_site(...)
#> !                        14.     └─pkgdown:::build_site_local(...)
#> !                        15.       └─pkgdown::build_reference(...)
#> !                        16.         ├─pkgdown:::unwrap_purrr_error(...)
#> !                        17.         │ └─base::withCallingHandlers(...)
#> !                        18.         └─purrr::map(...)
#> !                        19.           └─purrr:::map_("list", .x, .f, ..., .progress = .progress)
#> !                        20.             ├─purrr:::with_indexed_errors(...)
#> !                        21.             │ └─base::withCallingHandlers(...)
#> !                        22.             ├─purrr:::call_with_cleanup(...)
#> !                        23.             └─pkgdown (local) .f(.x[[i]], ...)
#> !                        24.               ├─base::withCallingHandlers(...)
#> !                        25.               └─pkgdown:::data_reference_topic(...)
#> !                        26.                 └─pkgdown:::run_examples(...)
#> !                        27.                   └─pkgdown:::highlight_examples(code, topic, env = env)
#> !                        28.                     └─downlit::evaluate_and_highlight(...)
#> !                        29.                       └─evaluate::evaluate(code, child_env(env), new_device = TRUE, output_handler = output_handler)
#> !                        30.                         ├─base::withRestarts(...)
#> !                        31.                         │ └─base (local) withRestartList(expr, restarts)
#> !                        32.                         │   ├─base (local) withOneRestart(withRestartList(expr, restarts[-nr]), restarts[[nr]])
#> !                        33.                         │   │ └─base (local) doWithOneRestart(return(expr), restart)
#> !                        34.                         │   └─base (local) withRestartList(expr, restarts[-nr])
#> !                        35.                         │     └─base (local) withOneRestart(expr, restarts[[1L]])
#> !                        36.                         │       └─base (local) doWithOneRestart(return(expr), restart)
#> !                        37.                         ├─evaluate:::with_handlers(...)
#> !                        38.                         │ ├─base::eval(call)
#> !                        39.                         │ │ └─base::eval(call)
#> !                        40.                         │ └─base::withCallingHandlers(...)
#> !                        41.                         ├─base::withVisible(eval(expr, envir))
#> !                        42.                         └─base::eval(expr, envir)
#> !                        43.                           └─base::eval(expr, envir)
#> !                        44.                             └─scop::ConvertHomologs(...)
#> !                        45.                               └─scop:::ConvertHomologs.Seurat(...)
#> !                        46.                                 └─scop:::ConvertHomologs.default(...)
#> !                        47.                                   └─scop::GeneConvert(...)
#> !                        48.                                     ├─thisutils::try_get(...)
#> !                        49.                                     │ ├─base::tryCatch(...)
#> !                        50.                                     │ │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#> !                        51.                                     │ │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#> !                        52.                                     │ │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#> !                        53.                                     │ └─base::eval.parent(substitute(expr))
#> !                        54.                                     │   └─base::eval(expr, p)
#> !                        55.                                     │     └─base::eval(expr, p)
#> !                        56.                                     └─biomaRt::listEnsemblArchives()
#> !                        57.                                       └─biomaRt:::.listEnsemblArchives(http_config = list())
#> !                        58.                                         └─biomaRt:::.checkArchiveList(http_config)
#> !                        59.                                           └─biomaRt:::.getArchiveList(http_config = http_config)
#> !                        60.                                             └─httr2::req_perform(html_request)
#> ! [2026-05-12 03:41:58] Get errors when connecting with EnsemblArchives...
#> ! [2026-05-12 03:41:59] Retrying...
#> ℹ [2026-05-12 03:42:07] Using the 115 version of ensembl database...
#> ℹ [2026-05-12 03:42:07] Downloading the ensembl database from https://sep2025.archive.ensembl.org...
#> ℹ [2026-05-12 03:42:08] Searching the dataset mmusculus ...
#> ℹ [2026-05-12 03:42:09] Connecting to the dataset mmusculus_gene_ensembl ...
#> ℹ [2026-05-12 03:42:10] Converting the geneIDs...
#> ℹ [2026-05-12 03:42:20] 14809 genes mapped with "ensembl_symbol"
#> ℹ [2026-05-12 03:42:23] 13 genes mapped with "entrez_symbol"
#> ℹ [2026-05-12 03:42:27] 14 genes mapped with "uniprot_symbol"
#> ℹ [2026-05-12 03:42:30] ==============================
#> ℹ                       14836 genes mapped
#> ℹ                       1162 genes unmapped
#> ℹ                       ==============================
#> ℹ [2026-05-12 03:42:35] Converted 13052 source genes to 12903 target homologs
rownames(pancreas_human)[1:5]
#> [1] "C2orf68"  "PGCKA1"   "C11orf58" "C3orf80"  "C8orf33" 
```
