# Prefetch cell cycle genes

Based on the human cell cycle genes, the cell cycle genes of the
corresponding species were captured by homologous gene conversion.

## Usage

``` r
CycGenePrefetch(
  species = "Homo_sapiens",
  Ensembl_version = NULL,
  mirror = NULL,
  max_tries = 5,
  use_cached_gene = TRUE,
  verbose = TRUE
)
```

## Arguments

- species:

  Latin names for animals, i.e., `"Homo_sapiens"`, `"Mus_musculus"`

- Ensembl_version:

  An integer specifying the Ensembl version. Default is `NULL`. If
  `NULL`, the latest version will be used.

- mirror:

  Specify an Ensembl mirror to connect to. The valid options here are
  `"www"`, `"uswest"`, `"useast"`, `"asia"`.

- max_tries:

  The maximum number of attempts to connect with the BioMart service.

- use_cached_gene:

  Whether to use previously cached cell cycle gene conversion results
  for the species. Default is `TRUE`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A list of S-phase and G2M-phase genes.

## See also

[GeneConvert](https://mengxu98.github.io/scop/reference/GeneConvert.md)

## Examples

``` r
ccgenes <- CycGenePrefetch("Homo_sapiens")
#> ℹ [2026-02-27 15:23:04] Prefetching cell cycle genes for "Homo_sapiens" ...
#> ✔ [2026-02-27 15:23:04] Cell cycle gene prefetching completed "Homo_sapiens"
str(ccgenes)
#> List of 3
#>  $ res: NULL
#>  $ S  : chr [1:43] "MCM5" "PCNA" "TYMS" "FEN1" ...
#>  $ G2M: chr [1:54] "HMGB2" "CDK1" "NUSAP1" "UBE2C" ...

ccgenes <- CycGenePrefetch("Mus_musculus")
#> ℹ [2026-02-27 15:23:04] Prefetching cell cycle genes for "Mus_musculus" ...
#> ℹ [2026-02-27 15:23:24] Connect to the Ensembl archives...
#> ! [2026-02-27 15:23:35] <error/httr2_failure>
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
#> !                        44.                             └─scop::CycGenePrefetch("Mus_musculus")
#> !                        45.                               └─scop::GeneConvert(...)
#> !                        46.                                 ├─thisutils::try_get(...)
#> !                        47.                                 │ ├─base::tryCatch(...)
#> !                        48.                                 │ │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#> !                        49.                                 │ │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#> !                        50.                                 │ │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#> !                        51.                                 │ └─base::eval.parent(substitute(expr))
#> !                        52.                                 │   └─base::eval(expr, p)
#> !                        53.                                 │     └─base::eval(expr, p)
#> !                        54.                                 └─biomaRt::listEnsemblArchives()
#> !                        55.                                   └─biomaRt:::.listEnsemblArchives(http_config = list())
#> !                        56.                                     └─biomaRt:::.checkArchiveList(http_config)
#> !                        57.                                       └─biomaRt:::.getArchiveList(http_config = http_config)
#> !                        58.                                         └─httr2::req_perform(html_request)
#> ! [2026-02-27 15:23:35] Get errors when connecting with EnsemblArchives...
#> ! [2026-02-27 15:23:36] Retrying...
#> ! [2026-02-27 15:23:46] <error/httr2_failure>
#> !                       Error in `req_perform()`:
#> !                       ! Failed to perform HTTP request.
#> !                       Caused by error in `curl::curl_fetch_memory()`:
#> !                       ! Timeout was reached [www.ensembl.org]:
#> !                       Operation timed out after 10002 milliseconds with 0 bytes received
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
#> !                        44.                             └─scop::CycGenePrefetch("Mus_musculus")
#> !                        45.                               └─scop::GeneConvert(...)
#> !                        46.                                 ├─thisutils::try_get(...)
#> !                        47.                                 │ ├─base::tryCatch(...)
#> !                        48.                                 │ │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#> !                        49.                                 │ │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#> !                        50.                                 │ │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#> !                        51.                                 │ └─base::eval.parent(substitute(expr))
#> !                        52.                                 │   └─base::eval(expr, p)
#> !                        53.                                 │     └─base::eval(expr, p)
#> !                        54.                                 └─biomaRt::listEnsemblArchives()
#> !                        55.                                   └─biomaRt:::.listEnsemblArchives(http_config = list())
#> !                        56.                                     └─biomaRt:::.checkArchiveList(http_config)
#> !                        57.                                       └─biomaRt:::.getArchiveList(http_config = http_config)
#> !                        58.                                         └─httr2::req_perform(html_request)
#> ! [2026-02-27 15:23:46] Get errors when connecting with EnsemblArchives...
#> ! [2026-02-27 15:23:47] Retrying...
#> ! [2026-02-27 15:23:57] <error/httr2_failure>
#> !                       Error in `req_perform()`:
#> !                       ! Failed to perform HTTP request.
#> !                       Caused by error in `curl::curl_fetch_memory()`:
#> !                       ! Timeout was reached [www.ensembl.org]:
#> !                       Operation timed out after 10002 milliseconds with 0 bytes received
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
#> !                        44.                             └─scop::CycGenePrefetch("Mus_musculus")
#> !                        45.                               └─scop::GeneConvert(...)
#> !                        46.                                 ├─thisutils::try_get(...)
#> !                        47.                                 │ ├─base::tryCatch(...)
#> !                        48.                                 │ │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#> !                        49.                                 │ │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#> !                        50.                                 │ │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#> !                        51.                                 │ └─base::eval.parent(substitute(expr))
#> !                        52.                                 │   └─base::eval(expr, p)
#> !                        53.                                 │     └─base::eval(expr, p)
#> !                        54.                                 └─biomaRt::listEnsemblArchives()
#> !                        55.                                   └─biomaRt:::.listEnsemblArchives(http_config = list())
#> !                        56.                                     └─biomaRt:::.checkArchiveList(http_config)
#> !                        57.                                       └─biomaRt:::.getArchiveList(http_config = http_config)
#> !                        58.                                         └─httr2::req_perform(html_request)
#> ! [2026-02-27 15:23:57] Get errors when connecting with EnsemblArchives...
#> ! [2026-02-27 15:23:58] Retrying...
#> ! [2026-02-27 15:24:09] <error/httr2_failure>
#> !                       Error in `req_perform()`:
#> !                       ! Failed to perform HTTP request.
#> !                       Caused by error in `curl::curl_fetch_memory()`:
#> !                       ! Timeout was reached [www.ensembl.org]:
#> !                       Operation timed out after 10002 milliseconds with 0 bytes received
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
#> !                        44.                             └─scop::CycGenePrefetch("Mus_musculus")
#> !                        45.                               └─scop::GeneConvert(...)
#> !                        46.                                 ├─thisutils::try_get(...)
#> !                        47.                                 │ ├─base::tryCatch(...)
#> !                        48.                                 │ │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#> !                        49.                                 │ │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#> !                        50.                                 │ │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#> !                        51.                                 │ └─base::eval.parent(substitute(expr))
#> !                        52.                                 │   └─base::eval(expr, p)
#> !                        53.                                 │     └─base::eval(expr, p)
#> !                        54.                                 └─biomaRt::listEnsemblArchives()
#> !                        55.                                   └─biomaRt:::.listEnsemblArchives(http_config = list())
#> !                        56.                                     └─biomaRt:::.checkArchiveList(http_config)
#> !                        57.                                       └─biomaRt:::.getArchiveList(http_config = http_config)
#> !                        58.                                         └─httr2::req_perform(html_request)
#> ! [2026-02-27 15:24:09] Get errors when connecting with EnsemblArchives...
#> ! [2026-02-27 15:24:10] Retrying...
#> ! [2026-02-27 15:24:20] <error/httr2_failure>
#> !                       Error in `req_perform()`:
#> !                       ! Failed to perform HTTP request.
#> !                       Caused by error in `curl::curl_fetch_memory()`:
#> !                       ! Timeout was reached [www.ensembl.org]:
#> !                       Operation timed out after 10002 milliseconds with 0 bytes received
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
#> !                        44.                             └─scop::CycGenePrefetch("Mus_musculus")
#> !                        45.                               └─scop::GeneConvert(...)
#> !                        46.                                 ├─thisutils::try_get(...)
#> !                        47.                                 │ ├─base::tryCatch(...)
#> !                        48.                                 │ │ └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#> !                        49.                                 │ │   └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#> !                        50.                                 │ │     └─base (local) doTryCatch(return(expr), name, parentenv, handler)
#> !                        51.                                 │ └─base::eval.parent(substitute(expr))
#> !                        52.                                 │   └─base::eval(expr, p)
#> !                        53.                                 │     └─base::eval(expr, p)
#> !                        54.                                 └─biomaRt::listEnsemblArchives()
#> !                        55.                                   └─biomaRt:::.listEnsemblArchives(http_config = list())
#> !                        56.                                     └─biomaRt:::.checkArchiveList(http_config)
#> !                        57.                                       └─biomaRt:::.getArchiveList(http_config = http_config)
#> !                        58.                                         └─httr2::req_perform(html_request)
#> ! [2026-02-27 15:24:20] Get errors when connecting with EnsemblArchives...
#> Error in try_get(expr = {    biomaRt::listEnsemblArchives()}, max_tries = max_tries, error_message = "Get errors when connecting with EnsemblArchives..."): <error/httr2_failure> Error in `req_perform()`: ! Failed to perform HTTP
#> request. Caused by error in `curl::curl_fetch_memory()`: ! Timeout was reached
#> [www.ensembl.org]: Operation timed out after 10002 milliseconds with 0 bytes
#> received --- Backtrace:  ▆  1. ├─base::tryCatch(...)  2. │ └─base (local)
#> tryCatchList(expr, classes, parentenv, handlers)  3. │ ├─base (local)
#> tryCatchOne(...)  4. │ │ └─base (local) doTryCatch(return(expr), name,
#> parentenv, handler)  5. │ └─base (local) tryCatchList(expr, names[-nh],
#> parentenv, handlers[-nh])  6. │ └─base (local) tryCatchOne(expr, names,
#> parentenv, handlers[[1L]])  7. │ └─base (local) doTryCatch(return(expr), name,
#> parentenv, handler)  8. ├─base::withCallingHandlers(...)  9.
#> ├─base::saveRDS(...)  10. ├─base::do.call(...)  11. ├─base (local) `<fn>`(...)
#> 12. └─global `<fn>`(...)  13.  └─pkgdown::build_site(...)  14.
#> └─pkgdown:::build_site_local(...)  15.  └─pkgdown::build_reference(...)  16.
#> ├─pkgdown:::unwrap_purrr_error(...)  17.  │ └─base::withCallingHandlers(...)
#> 18.  └─purrr::map(...)  19.  └─purrr:::map_("list", .x, .f, ..., .progress =
#> .progress)  20.  ├─purrr:::with_indexed_errors(...)  21.  │
#> └─base::withCallingHandlers(...)  22.  ├─purrr:::call_with_cleanup(...)  23.
#> └─pkgdown (local) .f(.x[[i]], ...)  24.  ├─base::withCallingHandlers(...)  25.
#> └─pkgdown:::data_reference_topic(...)  26.  └─pkgdown:::run_examples(...)  27.
#> └─pkgdown:::highlight_examples(code, topic, env = env)  28.
#> └─downlit::evaluate_and_highlight(...)  29.  └─evaluate::evaluate(code,
#> child_env(env), new_device = TRUE, output_handler = output_handler)  30.
#> ├─base::withRestarts(...)  31.  │ └─base (local) withRestartList(expr,
#> restarts)  32.  │ ├─base (local) withOneRestart(withRestartList(expr,
#> restarts[-nr]), restarts[[nr]])  33.  │ │ └─base (local)
#> doWithOneRestart(return(expr), restart)  34.  │ └─base (local)
#> withRestartList(expr, restarts[-nr])  35.  │ └─base (local) withOneRestart(expr,
#> restarts[[1L]])  36.  │ └─base (local) doWithOneRestart(return(expr), restart)
#> 37.  ├─evaluate:::with_handlers(...)  38.  │ ├─base::eval(call)  39.  │ │
#> └─base::eval(call)  40.  │ └─base::withCallingHandlers(...)  41.
#> ├─base::withVisible(eval(expr, envir))  42.  └─base::eval(expr, envir)  43.
#> └─base::eval(expr, envir)  44.  └─scop::CycGenePrefetch("Mus_musculus")  45.
#> └─scop::GeneConvert(...)  46.  ├─thisutils::try_get(...)  47.  │
#> ├─base::tryCatch(...)  48.  │ │ └─base (local) tryCatchList(expr, classes,
#> parentenv, handlers)  49.  │ │ └─base (local) tryCatchOne(expr, names,
#> parentenv, handlers[[1L]])  50.  │ │ └─base (local) doTryCatch(return(expr),
#> name, parentenv, handler)  51.  │ └─base::eval.parent(substitute(expr))  52.  │
#> └─base::eval(expr, p)  53.  │ └─base::eval(expr, p)  54.
#> └─biomaRt::listEnsemblArchives()  55.
#> └─biomaRt:::.listEnsemblArchives(http_config = list())  56.
#> └─biomaRt:::.checkArchiveList(http_config)  57.
#> └─biomaRt:::.getArchiveList(http_config = http_config)  58.
#> └─httr2::req_perform(html_request)
str(ccgenes)
#> List of 3
#>  $ res: NULL
#>  $ S  : chr [1:43] "MCM5" "PCNA" "TYMS" "FEN1" ...
#>  $ G2M: chr [1:54] "HMGB2" "CDK1" "NUSAP1" "UBE2C" ...
```
