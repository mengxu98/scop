# Run CytoSPACE spatial assignment

Assign reference single cells to spatial transcriptomics spots using a
native R/C++ implementation of the default CytoSPACE spot-level
workflow.

## Usage

``` r
RunCytoSPACE(
  srt,
  reference,
  reference_label,
  assay = NULL,
  reference_assay = NULL,
  layer = "counts",
  reference_layer = "counts",
  features = NULL,
  cell_fractions = NULL,
  n_cells_per_spot = NULL,
  mean_cell_numbers = 5,
  scRNA_max_transcripts_per_cell = 1500,
  sampling_method = "duplicates",
  seed = 1,
  prefix = "CytoSPACE",
  store_results = TRUE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- reference:

  Reference `Seurat` object containing annotated single cells.

- reference_label:

  Metadata column in `reference` with cell type labels.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- reference_assay:

  Assay used in `reference`.

- layer, reference_layer:

  Assay layers used for spatial and reference expression.

- features:

  Features used for assignment. If `NULL`, shared features are used.

- cell_fractions:

  Optional cell-type fractions. Provide a named numeric vector, one-row
  matrix/data.frame, or a spot-by-cell-type matrix/data.frame.
  Spot-level rows are aggregated to the global composition used by the
  default CytoSPACE assignment workflow.

- n_cells_per_spot:

  Optional number of cells assigned to each spatial spot. If `NULL`,
  counts are estimated from spatial RNA reads with `mean_cell_numbers`.

- mean_cell_numbers:

  Mean number of cells per spot. Default `5`, matching the CytoSPACE
  Visium default.

- scRNA_max_transcripts_per_cell:

  Maximum reference transcripts per cell before assignment. Default
  `1500`, matching CytoSPACE.

- sampling_method:

  Sampling method. Only `"duplicates"` is supported in the package
  runtime.

- seed:

  Random seed used for deterministic reference downsampling and
  duplicate sampling.

- prefix:

  Prefix for metadata columns.

- store_results:

  Whether to store detailed assignment results in `srt@tools`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A `Seurat` object with CytoSPACE metadata columns and detailed results
stored in `srt@tools[["CytoSPACE"]]`.

## Examples

``` r
data(visium_human_pancreas_sub)
data(pancreas_sub)
pancreas_human <- ConvertHomologs(
  pancreas_sub,
  species_from = "Mus_musculus",
  species_to = "Homo_sapiens",
  verbose = FALSE
)
#> ! [2026-05-23 14:17:04] <error/httr2_http_405>
#> !                       Error in `httr2::req_perform()`:
#> !                       ! HTTP 405 Method Not Allowed.
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
#> !                        56.                                     └─biomaRt::getBM(...)
#> !                        57.                                       └─biomaRt:::.submitQueryXML(...)
#> !                        58.                                         └─httr2::req_perform(req)
#> ! [2026-05-23 14:17:04] Get errors when retrieving information from the BioMart database
#> ! [2026-05-23 14:17:05] Retrying...
#> ! [2026-05-23 14:17:06] <error/httr2_http_405>
#> !                       Error in `httr2::req_perform()`:
#> !                       ! HTTP 405 Method Not Allowed.
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
#> !                        56.                                     └─biomaRt::getBM(...)
#> !                        57.                                       └─biomaRt:::.submitQueryXML(...)
#> !                        58.                                         └─httr2::req_perform(req)
#> ! [2026-05-23 14:17:06] Get errors when retrieving information from the BioMart database
#> ! [2026-05-23 14:17:07] Retrying...
#> ! [2026-05-23 14:17:08] <error/httr2_http_405>
#> !                       Error in `httr2::req_perform()`:
#> !                       ! HTTP 405 Method Not Allowed.
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
#> !                        56.                                     └─biomaRt::getBM(...)
#> !                        57.                                       └─biomaRt:::.submitQueryXML(...)
#> !                        58.                                         └─httr2::req_perform(req)
#> ! [2026-05-23 14:17:08] Get errors when retrieving information from the BioMart database
#> ! [2026-05-23 14:17:09] Retrying...
#> ! [2026-05-23 14:17:10] <error/httr2_http_405>
#> !                       Error in `httr2::req_perform()`:
#> !                       ! HTTP 405 Method Not Allowed.
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
#> !                        56.                                     └─biomaRt::getBM(...)
#> !                        57.                                       └─biomaRt:::.submitQueryXML(...)
#> !                        58.                                         └─httr2::req_perform(req)
#> ! [2026-05-23 14:17:10] Get errors when retrieving information from the BioMart database
#> ! [2026-05-23 14:17:11] Retrying...
#> ! [2026-05-23 14:17:12] <error/httr2_http_405>
#> !                       Error in `httr2::req_perform()`:
#> !                       ! HTTP 405 Method Not Allowed.
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
#> !                        56.                                     └─biomaRt::getBM(...)
#> !                        57.                                       └─biomaRt:::.submitQueryXML(...)
#> !                        58.                                         └─httr2::req_perform(req)
#> ! [2026-05-23 14:17:12] Get errors when retrieving information from the BioMart database
#> Error in try_get(expr = {    biomaRt::getBM(mart = mart1, attributes = c(from_attr, "ensembl_gene_id"),         filters = from_attr, values = list(geneID))}, max_tries = max_tries, error_message = "Get errors when retrieving information from the BioMart database"): <error/httr2_http_405> Error in `httr2::req_perform()`: ! HTTP 405
#> Method Not Allowed. --- Backtrace:  ▆  1. ├─base::tryCatch(...)  2. │ └─base
#> (local) tryCatchList(expr, classes, parentenv, handlers)  3. │ ├─base (local)
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
#> └─base::eval(expr, envir)  44.  └─scop::ConvertHomologs(...)  45.
#> └─scop:::ConvertHomologs.Seurat(...)  46.
#> └─scop:::ConvertHomologs.default(...)  47.  └─scop::GeneConvert(...)  48.
#> ├─thisutils::try_get(...)  49.  │ ├─base::tryCatch(...)  50.  │ │ └─base
#> (local) tryCatchList(expr, classes, parentenv, handlers)  51.  │ │ └─base
#> (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])  52.  │ │ └─base
#> (local) doTryCatch(return(expr), name, parentenv, handler)  53.  │
#> └─base::eval.parent(substitute(expr))  54.  │ └─base::eval(expr, p)  55.  │
#> └─base::eval(expr, p)  56.  └─biomaRt::getBM(...)  57.
#> └─biomaRt:::.submitQueryXML(...)  58.  └─httr2::req_perform(req)
spatial <- RunCytoSPACE(
  visium_human_pancreas_sub,
  reference = pancreas_human,
  reference_label = "CellType",
  mean_cell_numbers = 1,
  verbose = FALSE
)
#> Error: object 'pancreas_human' not found

SpatialSpotPlot(
  visium_human_pancreas_sub,
  group.by = "coda_label",
  theme_use = "theme_scop"
)


SpatialSpotPlot(
  spatial,
  group.by = "CytoSPACE_dominant_type",
  theme_use = "theme_scop"
)
#> Error: object 'spatial' not found
```
