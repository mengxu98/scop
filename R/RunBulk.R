#' @title Unified bulk analysis entrypoint
#'
#' @description
#' `RunBulk()` provides one bulk-strategy entrypoint for Seurat-centered
#' workflows. Methods are selected with one character vector, so users do not
#' need to choose a separate module argument.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams standard_scop
#' @param srt A `Seurat` object used for pseudobulk mode.
#' @param bulk_se A `SummarizedExperiment` object used for true bulk mode.
#' @param method Character vector specifying methods to run.
#' Supported values are `"de_limma_voom"`, `"de_edgeR_qlf"`, `"de_DESeq2"`,
#' `"de_dream"`, `"deconv_MuSiC"`, `"deconv_BisqueRNA"`,
#' `"deconv_BayesPrism"`, and `"csde_TOAST"`.
#' @param sample.by Metadata column containing sample IDs in pseudobulk mode.
#' @param condition.by Metadata column containing condition labels.
#' Required when `.arg method` includes `"de"` or `"csde"`.
#' Optional for deconvolution-only runs.
#' @param group.by Optional metadata column used to stratify DE by subgroup.
#' @param ref_srt Optional `Seurat` reference for deconvolution/CSDE.
#' If omitted in pseudobulk mode, `srt` is used as reference.
#' @param celltype.by Metadata column in `ref_srt` defining reference cell types.
#' Required when deconvolution or CSDE is requested.
#' @param bulk_assay Assay name in `bulk_se` used as counts matrix.
#' @param ref_assay Assay name in `ref_srt` for reference profiles.
#' @param ref_layer Layer name in `ref_srt` for reference counts.
#' @param condition1 First condition label for pairwise comparison.
#' @param condition2 Second condition label for pairwise comparison.
#' @param de_markers_type DE comparison mode.
#' One of `"single"` (one pairwise comparison),
#' `"all"` (each condition vs the rest),
#' or `"paired"` (all pairwise condition comparisons).
#' @param markers_type Alias of `.arg de_markers_type` for consistency with [RunDEtest].
#' If provided, it overrides `.arg de_markers_type`.
#' @param run_enrichment Whether to run [RunEnrichment] from bulk DE results.
#' @param run_gsea Whether to run [RunGSEA] from bulk DE results.
#' @param enrichment_args Named list of additional parameters for [RunEnrichment].
#' @param gsea_args Named list of additional parameters for [RunGSEA].
#' @param de_args Named list of additional parameters for DE method functions.
#' @param deconv_args Named list of additional parameters for deconvolution method functions.
#' Current deconvolution wrappers use `backend = "internal"` by default and
#' record the computational engine in the result bundle.
#' @param csde_args Named list of additional parameters for CSDE method functions.
#' Current `csde_TOAST` uses `backend = "limma_interaction"` by default and
#' records the computational engine in the result bundle.
#' @param method_args Named list of method-specific parameters. Names should be
#' RunBulk method names, such as `"de_edgeR_qlf"` or `"csde_TOAST"`.
#' Method-specific values override `.arg de_args`, `.arg deconv_args`,
#' `.arg csde_args`, and `.arg ...`.
#' @param ... Additional parameters forwarded to method functions.
#'
#' @return
#' Returns a `Seurat` object in pseudobulk mode with results in
#' `srt@tools[["Bulk"]]`. If `.arg run_enrichment` or `.arg run_gsea` is
#' `TRUE`, successful pathway results are also stored in the standard
#' `srt@tools[["Enrichment_Bulk_wilcox"]]` and
#' `srt@tools[["GSEA_Bulk_wilcox"]]` slots used by [EnrichmentPlot] and
#' [GSEAPlot]. True bulk mode returns a `SummarizedExperiment` object
#' with results in `metadata(bulk_se)[["Bulk"]]`.
#'
#' @examples
#' # ----------------------------
#' # Example 1: Seurat -> pseudobulk (DE) using scop built-in dataset
#' # ----------------------------
#' if (requireNamespace("edgeR", quietly = TRUE)) {
#'   data("pancreas_sub", package = "scop")
#'   srt <- pancreas_sub
#'
#'   # Build sample/condition labels for pseudobulk aggregation.
#'   srt$sample_id <- paste0(
#'     "sample_",
#'     (seq_len(ncol(srt)) - 1L) %% 8L + 1L
#'   )
#'   srt$condition <- ifelse(
#'     srt$sample_id %in% paste0("sample_", 1:4),
#'     "A",
#'     "B"
#'   )
#'
#'   srt <- RunBulk(
#'     srt = srt,
#'     method = "de_edgeR_qlf",
#'     sample.by = "sample_id",
#'     condition.by = "condition",
#'     verbose = FALSE
#'   )
#'   srt@tools$Bulk$status$de$status
#' }
#'
#' # ----------------------------
#' # Example 2: pure bulk (SummarizedExperiment -> DE)
#' # ----------------------------
#' if (requireNamespace("edgeR", quietly = TRUE) &&
#'   requireNamespace("SummarizedExperiment", quietly = TRUE) &&
#'   requireNamespace("S4Vectors", quietly = TRUE)) {
#'   data("pancreas_sub", package = "scop")
#'   srt <- pancreas_sub
#'   srt$sample_id <- paste0(
#'     "sample_",
#'     (seq_len(ncol(srt)) - 1L) %% 8L + 1L
#'   )
#'   srt$condition <- ifelse(
#'     srt$sample_id %in% paste0("sample_", 1:4),
#'     "A",
#'     "B"
#'   )
#'
#'   counts <- GetAssayData5(srt, layer = "counts")
#'   sid <- as.character(srt$sample_id)
#'   sid_levels <- unique(sid)
#'   bulk_counts <- do.call(cbind, lapply(sid_levels, function(x) {
#'     Matrix::rowSums(counts[, sid == x, drop = FALSE])
#'   }))
#'   colnames(bulk_counts) <- sid_levels
#'   rownames(bulk_counts) <- rownames(counts)
#'
#'   sample_condition <- tapply(
#'     as.character(srt$condition),
#'     sid,
#'     function(x) unique(x)[1]
#'   )
#'   bulk_se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(counts = as.matrix(bulk_counts)),
#'     colData = S4Vectors::DataFrame(
#'       condition = as.character(sample_condition[colnames(bulk_counts)]),
#'       row.names = colnames(bulk_counts)
#'     )
#'   )
#'
#'   bulk_se <- RunBulk(
#'     bulk_se = bulk_se,
#'     method = "de_edgeR_qlf",
#'     condition.by = "condition",
#'     verbose = FALSE
#'   )
#'   S4Vectors::metadata(bulk_se)$Bulk$status$de$status
#' }
#'
#' # ----------------------------
#' # Example 3: deconvolution + CSDE in one call
#' # ----------------------------
#' if (requireNamespace("limma", quietly = TRUE)) {
#'   data("pancreas_sub", package = "scop")
#'   srt <- pancreas_sub
#'   srt$sample_id <- paste0(
#'     "sample_",
#'     (seq_len(ncol(srt)) - 1L) %% 8L + 1L
#'   )
#'   srt$condition <- ifelse(
#'     srt$sample_id %in% paste0("sample_", 1:4),
#'     "A",
#'     "B"
#'   )
#'   if ("CellType" %in% colnames(srt@meta.data)) {
#'     srt$celltype <- as.character(srt$CellType)
#'   } else {
#'     srt$celltype <- as.character(SeuratObject::Idents(srt))
#'   }
#'
#'   srt <- RunBulk(
#'     srt = srt,
#'     method = c("deconv_MuSiC", "csde_TOAST"),
#'     sample.by = "sample_id",
#'     condition.by = "condition",
#'     condition1 = "A",
#'     condition2 = "B",
#'     celltype.by = "celltype",
#'     verbose = FALSE
#'   )
#'   srt@tools$Bulk$status$deconv$status
#'   srt@tools$Bulk$status$csde$status
#' }
#'
#' @export
RunBulk <- function(
  srt = NULL,
  bulk_se = NULL,
  method = "de_edgeR_qlf",
  sample.by = NULL,
  condition.by = NULL,
  group.by = NULL,
  ref_srt = NULL,
  celltype.by = NULL,
  assay = NULL,
  layer = "counts",
  bulk_assay = "counts",
  ref_assay = NULL,
  ref_layer = "counts",
  condition1 = NULL,
  condition2 = NULL,
  de_markers_type = "single",
  markers_type = NULL,
  run_enrichment = FALSE,
  run_gsea = FALSE,
  enrichment_args = list(),
  gsea_args = list(),
  de_args = list(),
  deconv_args = list(),
  csde_args = list(),
  method_args = list(),
  verbose = TRUE,
  ...
) {
  mode <- .bulk_detect_mode(srt = srt, bulk_se = bulk_se)
  method_plan <- .bulk_validate_method_spec(method = method)
  method_modules <- vapply(method_plan, `[[`, character(1), "module")
  needs_condition <- any(method_modules %in% c("de", "csde"))
  if (isTRUE(needs_condition) && (is.null(condition.by) || !nzchar(condition.by))) {
    .bulk_abort(
      "{.arg condition.by} is required when {.arg method} includes {.val de} or {.val csde}."
    )
  }
  if (!is.null(markers_type)) {
    de_markers_type <- markers_type
  }
  de_markers_type <- .bulk_normalize_de_markers_type(de_markers_type)
  if (!is.list(enrichment_args) || !is.list(gsea_args)) {
    .bulk_abort("{.arg enrichment_args} and {.arg gsea_args} must be named lists.")
  }
  extra_args <- list(...)
  method_args <- .bulk_normalize_method_args(method_args)

  ctx <- .bulk_build_context(
    mode = mode,
    srt = srt,
    bulk_se = bulk_se,
    sample.by = sample.by,
    condition.by = condition.by,
    group.by = group.by,
    assay = assay,
    layer = layer,
    bulk_assay = bulk_assay
  )

  needs_pairwise_condition <- isTRUE("csde" %in% method_modules) ||
    isTRUE("de" %in% method_modules && identical(de_markers_type, "single"))

  if (isTRUE("de" %in% method_modules) &&
    !identical(de_markers_type, "single") &&
    (!is.null(condition1) || !is.null(condition2))) {
    log_message(
      "{.arg condition1} and {.arg condition2} are ignored when {.arg de_markers_type} is {.val all} or {.val paired}.",
      message_type = "warning",
      verbose = verbose
    )
  }

  if (isTRUE(needs_pairwise_condition) && (is.null(condition1) || is.null(condition2))) {
    unique_condition <- unique(as.character(ctx$condition_global))
    if (length(unique_condition) > 2) {
      .bulk_abort(
        "More than two conditions were detected. Please provide {.arg condition1} and {.arg condition2}."
      )
    }
  }

  needs_reference <- any(method_modules %in% c("deconv", "csde"))
  ref_use <- ref_srt
  if (needs_reference && is.null(ref_use) && identical(mode, "pseudobulk")) {
    ref_use <- srt
  }
  if (needs_reference && is.null(ref_use)) {
    .bulk_abort(
      "{.arg ref_srt} is required when {.arg method} includes a deconvolution or CSDE method in true bulk mode."
    )
  }
  if (needs_reference && is.null(celltype.by)) {
    .bulk_abort(
      "{.arg celltype.by} is required when {.arg method} includes {.val deconv} or {.val csde}."
    )
  }
  if (needs_reference) {
    ctx$reference <- .bulk_build_reference_profiles(
      ref_srt = ref_use,
      celltype.by = celltype.by,
      assay = ref_assay %||% assay,
      layer = ref_layer
    )
  }

  bulk_store <- list(
    input = list(
      mode = mode,
      sample.by = sample.by,
      condition.by = condition.by,
      condition.by_used = ctx$condition_col,
      group.by = group.by,
      assay = ctx$assay,
      layer = ctx$layer,
      bulk_assay = ctx$bulk_assay,
      celltype.by = celltype.by
    ),
    de = list(
      active_method = NULL,
      methods = list(),
      results = .bulk_empty_de_result(),
      parameters = list(),
      details = list()
    ),
    deconv = list(
      active_method = NULL,
      methods = list(),
      results = .bulk_empty_deconv_result(),
      parameters = list(),
      details = list()
    ),
    csde = list(
      active_method = NULL,
      methods = list(),
      results = .bulk_empty_csde_result(),
      parameters = list(),
      details = list()
    ),
    enrichment = list(status = "skipped", reason = "Not requested", results = NULL),
    gsea = list(status = "skipped", reason = "Not requested", results = NULL),
    status = list(),
    params = list(
      method = vapply(method_plan, `[[`, character(1), "method"),
      condition1 = condition1,
      condition2 = condition2,
      de_markers_type = de_markers_type,
      run_enrichment = run_enrichment,
      run_gsea = run_gsea,
      enrichment_args = enrichment_args,
      gsea_args = gsea_args
    )
  )

  bulk_store$results <- list(
    de = bulk_store$de$results,
    deconv = bulk_store$deconv$results,
    csde = bulk_store$csde$results
  )
  bulk_store$active_method <- list(
    de = NULL,
    deconv = NULL,
    csde = NULL
  )
  bulk_store$methods <- list(
    de = list(),
    deconv = list(),
    csde = list()
  )
  bulk_store$status$methods <- list()

  module_order <- c("de", "deconv", "csde")
  method_plan <- method_plan[order(match(method_modules, module_order), seq_along(method_plan))]
  for (task in method_plan) {
    module_name <- task$module
    method_name <- task$method
    if (identical(module_name, "de")) {
      bundle <- .RunBulk_dispatch_de(
        ctx = ctx,
        method_name = method_name,
        condition1 = condition1,
        condition2 = condition2,
        de_markers_type = de_markers_type,
        de_args = .bulk_collect_method_args(
          method_name = method_name,
          global_args = extra_args,
          module_args = de_args,
          method_args = method_args
        ),
        verbose = verbose
      )
      bundle$results <- .bulk_coerce_de_schema(bundle$results)
    } else if (identical(module_name, "deconv")) {
      bundle <- .RunBulk_dispatch_deconv(
        ctx = ctx,
        method_name = method_name,
        deconv_args = .bulk_collect_method_args(
          method_name = method_name,
          global_args = extra_args,
          module_args = deconv_args,
          method_args = method_args
        ),
        verbose = verbose
      )
      bundle$results <- .bulk_coerce_deconv_schema(bundle$results)
    } else {
      if (!"deconv" %in% method_modules) {
        .bulk_abort(
          "{.val csde_TOAST} requires a deconvolution method in the same {.fn RunBulk} call."
        )
      }
      deconv_bundle <- bulk_store$deconv$methods[[bulk_store$deconv$active_method]]
      if (is.null(deconv_bundle)) {
        .bulk_abort(
          "{.val csde_TOAST} requires successful deconvolution results in the same {.fn RunBulk} call."
        )
      }
      bundle <- .RunBulk_dispatch_csde(
        ctx = ctx,
        method_name = method_name,
        deconv_bundle = deconv_bundle,
        condition1 = condition1,
        condition2 = condition2,
        csde_args = .bulk_collect_method_args(
          method_name = method_name,
          global_args = extra_args,
          module_args = csde_args,
          method_args = method_args
        ),
        verbose = verbose
      )
      bundle$results <- .bulk_coerce_csde_schema(bundle$results)
    }

    bulk_store[[module_name]]$active_method <- method_name
    bulk_store[[module_name]]$methods[[method_name]] <- bundle
    bulk_store[[module_name]]$results <- bundle$results
    bulk_store[[module_name]]$parameters <- bundle$parameters %||% list()
    bulk_store[[module_name]]$details <- bundle$details %||% list()
    bulk_store$active_method[[module_name]] <- method_name
    bulk_store$methods[[module_name]][[method_name]] <- bundle
    bulk_store$results[[module_name]] <- bundle$results
    bulk_store$status[[module_name]] <- list(
      method = method_name,
      status = bundle$status %||% "failed",
      reason = bundle$reason %||% NULL
    )
    bulk_store$status$methods[[method_name]] <- bulk_store$status[[module_name]]
  }

  if (isTRUE(run_enrichment)) {
    bulk_store$enrichment <- .RunBulk_dispatch_enrichment(
      de_results = bulk_store$results$de,
      enrichment_args = enrichment_args,
      verbose = verbose
    )
  }
  if (isTRUE(run_gsea)) {
    bulk_store$gsea <- .RunBulk_dispatch_gsea(
      de_results = bulk_store$results$de,
      gsea_args = gsea_args,
      verbose = verbose
    )
  }

  if (identical(mode, "pseudobulk")) {
    srt <- .bulk_store_downstream_tool_slots(
      srt = srt,
      bulk_store = bulk_store
    )
    srt@tools[["Bulk"]] <- bulk_store
    return(srt)
  }

  metadata_bulk <- S4Vectors::metadata(bulk_se)
  metadata_bulk[["Bulk"]] <- bulk_store
  S4Vectors::metadata(bulk_se) <- metadata_bulk
  bulk_se
}

.RunBulk_dispatch_de <- function(
  ctx,
  method_name,
  condition1 = NULL,
  condition2 = NULL,
  de_markers_type = "single",
  de_args = list(),
  verbose = TRUE
) {
  method_map <- list(
    de_limma_voom = RunBulk_de_limma_voom,
    de_edgeR_qlf = RunBulk_de_edgeR,
    de_DESeq2 = RunBulk_de_deseq2,
    de_dream = RunBulk_de_dream
  )
  method_fun <- method_map[[method_name]]
  if (is.null(method_fun)) {
    .bulk_abort(
      paste0("Unsupported DE method: ", method_name)
    )
  }

  per_group <- list()
  results_list <- list()
  fail_messages <- character(0)

  for (grp in names(ctx$counts_by_group)) {
    condition_vec <- ctx$condition_by_group[[grp]]
    comparison_specs <- tryCatch(
      .bulk_prepare_de_comparisons(
        condition = condition_vec,
        condition1 = condition1,
        condition2 = condition2,
        markers_type = de_markers_type
      ),
      error = function(e) e
    )
    if (inherits(comparison_specs, "error")) {
      fail_messages <- c(
        fail_messages,
        paste0(grp, ": ", comparison_specs$message)
      )
      next
    }

    group_bundles <- list()
    for (spec in comparison_specs) {
      args <- utils::modifyList(
        list(
          count_matrix = ctx$counts_by_group[[grp]],
          condition = spec$condition,
          sample_data = ctx$sample_meta_by_group[[grp]],
          condition1 = spec$condition1,
          condition2 = spec$condition2,
          verbose = verbose
        ),
        de_args
      )
      bundle <- tryCatch(
        invoke_fun(
          method_fun,
          args[names(args) %in% names(formals(method_fun))]
        ),
        error = function(e) {
          .bulk_method_failed(reason = e$message)
        }
      )
      if (!is.list(bundle) || is.null(bundle$status)) {
        bundle <- .bulk_method_failed(reason = "Method did not return a valid bundle.")
      }
      group_bundles[[spec$label]] <- bundle

      if (identical(bundle$status, "success")) {
        df <- as.data.frame(bundle$results)
        if (nrow(df) > 0) {
          if (!"gene" %in% colnames(df)) {
            df$gene <- rownames(df)
          }
          df$group1 <- grp
          df$group2 <- spec$label
          df$comparison <- .bulk_compose_comparison_label(grp, spec$label)
          df$method <- method_name
          results_list[[paste0(grp, "::", spec$label)]] <- .bulk_coerce_de_schema(df)
        }
      } else {
        fail_messages <- c(
          fail_messages,
          paste0(grp, "[", spec$label, "]: ", bundle$reason %||% "failed")
        )
      }
    }

    per_group[[grp]] <- list(
      markers_type = de_markers_type,
      comparisons = group_bundles
    )
  }

  results <- if (length(results_list) == 0) {
    .bulk_empty_de_result()
  } else {
    out <- do.call(rbind, results_list)
    rownames(out) <- NULL
    out
  }

  if (nrow(results) > 0) {
    .bulk_method_success(
      results = results,
      details = list(groups = per_group),
      parameters = list(
        method = method_name,
        de_markers_type = de_markers_type
      )
    )
  } else {
    .bulk_method_failed(
      reason = paste(fail_messages, collapse = "; "),
      details = list(groups = per_group),
      parameters = list(
        method = method_name,
        de_markers_type = de_markers_type
      ),
      results = .bulk_empty_de_result()
    )
  }
}

.bulk_prepare_de_comparisons <- function(
  condition,
  condition1 = NULL,
  condition2 = NULL,
  markers_type = "single"
) {
  markers_type <- .bulk_normalize_de_markers_type(markers_type)
  cond <- as.character(condition)
  cond_use <- cond[!is.na(cond)]
  levels_use <- unique(cond_use)

  if (identical(markers_type, "single")) {
    pair <- .bulk_resolve_condition_pair(
      condition = cond,
      condition1 = condition1,
      condition2 = condition2,
      strict_two_levels = TRUE
    )
    return(list(list(
      condition = cond,
      condition1 = pair$condition1,
      condition2 = pair$condition2,
      label = paste0(pair$condition2, "_vs_", pair$condition1)
    )))
  }

  if (length(levels_use) < 2) {
    .bulk_abort("At least two condition levels are required for DE comparison.")
  }

  if (identical(markers_type, "all")) {
    return(lapply(levels_use, function(lv) {
      cond_binary <- ifelse(!is.na(cond) & cond == lv, lv, "others")
      list(
        condition = cond_binary,
        condition1 = "others",
        condition2 = lv,
        label = paste0(lv, "_vs_others")
      )
    }))
  }

  pairs <- utils::combn(levels_use, 2, simplify = FALSE)
  lapply(pairs, function(pr) {
    list(
      condition = cond,
      condition1 = pr[[1]],
      condition2 = pr[[2]],
      label = paste0(pr[[2]], "_vs_", pr[[1]])
    )
  })
}

.bulk_normalize_de_markers_type <- function(markers_type) {
  if (!is.character(markers_type) || length(markers_type) != 1) {
    .bulk_abort("{.arg de_markers_type} must be a single string.")
  }
  key <- gsub("[^a-z]", "", tolower(markers_type))
  map <- c(
    single = "single",
    all = "all",
    onevsrest = "all",
    onevsall = "all",
    paired = "paired",
    pairwise = "paired"
  )
  out <- unname(map[[key]])
  if (is.null(out) || is.na(out)) {
    .bulk_abort(
      "Unsupported {.arg de_markers_type}. Use one of {.val single}, {.val all}, or {.val paired}."
    )
  }
  out
}

.RunBulk_dispatch_deconv <- function(
  ctx,
  method_name,
  deconv_args = list(),
  verbose = TRUE
) {
  method_map <- list(
    deconv_MuSiC = RunBulk_deconv_music,
    deconv_BisqueRNA = RunBulk_deconv_bisque,
    deconv_BayesPrism = RunBulk_deconv_bayesprism
  )
  method_fun <- method_map[[method_name]]
  if (is.null(method_fun)) {
    .bulk_abort(
      paste0("Unsupported deconvolution method: ", method_name)
    )
  }

  args <- utils::modifyList(
    list(
      count_matrix = ctx$counts_global,
      reference_matrix = ctx$reference$profiles,
      verbose = verbose
    ),
    deconv_args
  )
  bundle <- tryCatch(
    invoke_fun(
      method_fun,
      args[names(args) %in% names(formals(method_fun))]
    ),
    error = function(e) {
      .bulk_method_failed(reason = e$message)
    }
  )

  if (!is.list(bundle) || is.null(bundle$status)) {
    bundle <- .bulk_method_failed(reason = "Method did not return a valid bundle.")
  }
  bundle$parameters <- utils::modifyList(
    list(method = method_name),
    bundle$parameters %||% list()
  )
  bundle$results <- .bulk_coerce_deconv_schema(bundle$results)
  if (nrow(bundle$results) > 0) {
    bundle$results$method <- method_name
  }
  bundle
}

.RunBulk_dispatch_csde <- function(
  ctx,
  method_name,
  deconv_bundle,
  condition1 = NULL,
  condition2 = NULL,
  csde_args = list(),
  verbose = TRUE
) {
  method_map <- list(
    csde_TOAST = RunBulk_csde_toast
  )
  method_fun <- method_map[[method_name]]
  if (is.null(method_fun)) {
    .bulk_abort(
      paste0("Unsupported CSDE method: ", method_name)
    )
  }

  prop_matrix <- deconv_bundle$details$proportion_matrix %||%
    .bulk_deconv_results_to_matrix(deconv_bundle$results)
  if (is.null(prop_matrix) || nrow(prop_matrix) == 0 || ncol(prop_matrix) == 0) {
    return(.bulk_method_failed(
      reason = "CSDE requires successful deconvolution results.",
      parameters = list(method = method_name),
      results = .bulk_empty_csde_result()
    ))
  }

  args <- utils::modifyList(
    list(
      count_matrix = ctx$counts_global,
      condition = ctx$condition_global,
      proportions = prop_matrix,
      condition1 = condition1,
      condition2 = condition2,
      verbose = verbose
    ),
    csde_args
  )
  bundle <- tryCatch(
    invoke_fun(
      method_fun,
      args[names(args) %in% names(formals(method_fun))]
    ),
    error = function(e) {
      .bulk_method_failed(reason = e$message)
    }
  )
  if (!is.list(bundle) || is.null(bundle$status)) {
    bundle <- .bulk_method_failed(reason = "Method did not return a valid bundle.")
  }
  bundle$parameters <- utils::modifyList(
    list(method = method_name),
    bundle$parameters %||% list()
  )
  bundle$results <- .bulk_coerce_csde_schema(bundle$results)
  if (nrow(bundle$results) > 0) {
    bundle$results$method <- method_name
  }
  bundle
}

.RunBulk_dispatch_enrichment <- function(
  de_results,
  enrichment_args = list(),
  verbose = TRUE
) {
  DE_threshold <- enrichment_args$DE_threshold %||%
    "avg_log2FC > 0 & p_val_adj < 0.05"
  de_results <- tryCatch(
    .bulk_filter_de_results(
      de_results = de_results,
      DE_threshold = DE_threshold
    ),
    error = function(e) e
  )
  if (inherits(de_results, "error")) {
    return(list(
      status = "failed",
      reason = de_results$message,
      results = NULL
    ))
  }
  de_use <- .bulk_prepare_de_for_pathway(
    de_results = de_results,
    require_score = FALSE
  )
  if (is.null(de_use) || nrow(de_use) == 0) {
    return(list(
      status = "skipped",
      reason = "No DE results available for enrichment.",
      results = NULL
    ))
  }

  run_args <- utils::modifyList(
    enrichment_args,
    list(
      srt = NULL,
      geneID = de_use$gene,
      geneID_groups = de_use$comparison,
      verbose = verbose
    )
  )
  run_args <- run_args[names(run_args) %in% names(formals(RunEnrichment))]
  res <- tryCatch(
    invoke_fun(RunEnrichment, run_args),
    error = function(e) {
      e
    }
  )
  if (inherits(res, "error")) {
    return(list(
      status = "failed",
      reason = res$message,
      results = NULL
    ))
  }
  res[["DE_threshold"]] <- DE_threshold
  list(status = "success", reason = NULL, results = res)
}

.bulk_store_downstream_tool_slots <- function(
  srt,
  bulk_store,
  test.use = "wilcox"
) {
  if (!inherits(srt, "Seurat")) {
    return(srt)
  }
  if (
    is.list(bulk_store$enrichment) &&
      identical(bulk_store$enrichment$status, "success") &&
      !is.null(bulk_store$enrichment$results)
  ) {
    srt@tools[[paste("Enrichment", "Bulk", test.use, sep = "_")]] <-
      bulk_store$enrichment$results
  }
  if (
    is.list(bulk_store$gsea) &&
      identical(bulk_store$gsea$status, "success") &&
      !is.null(bulk_store$gsea$results)
  ) {
    srt@tools[[paste("GSEA", "Bulk", test.use, sep = "_")]] <-
      bulk_store$gsea$results
  }
  srt
}

.RunBulk_dispatch_gsea <- function(
  de_results,
  gsea_args = list(),
  verbose = TRUE
) {
  DE_threshold <- gsea_args$DE_threshold %||% "p_val_adj < 0.05"
  de_results <- tryCatch(
    .bulk_filter_de_results(
      de_results = de_results,
      DE_threshold = DE_threshold
    ),
    error = function(e) e
  )
  if (inherits(de_results, "error")) {
    return(list(
      status = "failed",
      reason = de_results$message,
      results = NULL
    ))
  }
  de_use <- .bulk_prepare_de_for_pathway(
    de_results = de_results,
    require_score = TRUE
  )
  if (is.null(de_use) || nrow(de_use) == 0) {
    return(list(
      status = "skipped",
      reason = "No DE results available for GSEA.",
      results = NULL
    ))
  }

  run_args <- utils::modifyList(
    gsea_args,
    list(
      srt = NULL,
      geneID = de_use$gene,
      geneScore = de_use$avg_log2FC,
      geneID_groups = de_use$comparison,
      verbose = verbose
    )
  )
  run_args <- run_args[names(run_args) %in% names(formals(RunGSEA))]
  res <- tryCatch(
    invoke_fun(RunGSEA, run_args),
    error = function(e) {
      e
    }
  )
  if (inherits(res, "error")) {
    return(list(
      status = "failed",
      reason = res$message,
      results = NULL
    ))
  }
  res[["DE_threshold"]] <- DE_threshold
  list(status = "success", reason = NULL, results = res)
}

.bulk_filter_de_results <- function(
  de_results,
  DE_threshold = NULL
) {
  if (is.null(de_results) || nrow(de_results) == 0) {
    return(de_results)
  }
  if (is.null(DE_threshold)) {
    return(de_results)
  }
  DE_threshold <- as.character(DE_threshold)
  if (length(DE_threshold) != 1 || !nzchar(DE_threshold)) {
    .bulk_abort("{.arg DE_threshold} must be a single non-empty expression.")
  }
  df <- as.data.frame(de_results, stringsAsFactors = FALSE)
  keep <- with(df, eval(rlang::parse_expr(DE_threshold)))
  if (length(keep) == 1) {
    keep <- rep(isTRUE(keep), nrow(df))
  }
  if (length(keep) != nrow(df)) {
    .bulk_abort("{.arg DE_threshold} must evaluate to one value per DE row.")
  }
  keep <- !is.na(keep) & as.logical(keep)
  df[keep, , drop = FALSE]
}

.bulk_prepare_de_for_pathway <- function(
  de_results,
  require_score = FALSE
) {
  if (is.null(de_results) || nrow(de_results) == 0) {
    return(NULL)
  }
  df <- as.data.frame(de_results, stringsAsFactors = FALSE)
  if (!"gene" %in% colnames(df)) {
    return(NULL)
  }

  df$gene <- as.character(df$gene)
  df <- df[!is.na(df$gene) & nzchar(df$gene), , drop = FALSE]
  if (nrow(df) == 0) {
    return(NULL)
  }

  df$comparison <- .bulk_make_comparison_label(df)
  keep_cols <- c("gene", "comparison")
  if (isTRUE(require_score)) {
    if (!"avg_log2FC" %in% colnames(df)) {
      return(NULL)
    }
    df$avg_log2FC <- suppressWarnings(as.numeric(df$avg_log2FC))
    df <- df[!is.na(df$avg_log2FC), , drop = FALSE]
    if (nrow(df) == 0) {
      return(NULL)
    }
    # Keep one ranking score per gene in each comparison.
    ord <- order(df$comparison, df$gene, -abs(df$avg_log2FC))
    df <- df[ord, , drop = FALSE]
    df <- df[!duplicated(paste(df$comparison, df$gene, sep = "::")), , drop = FALSE]
    keep_cols <- c(keep_cols, "avg_log2FC")
  } else {
    df <- df[!duplicated(paste(df$comparison, df$gene, sep = "::")), , drop = FALSE]
  }
  rownames(df) <- NULL
  df[, keep_cols, drop = FALSE]
}

.bulk_get_de_results <- function(srt) {
  if (!inherits(srt, "Seurat") || is.null(srt@tools[["Bulk"]])) {
    return(NULL)
  }
  bulk <- srt@tools[["Bulk"]]
  de_res <- bulk[["results"]][["de"]] %||% bulk[["de"]][["results"]]
  if (is.null(de_res) || nrow(de_res) == 0) {
    return(NULL)
  }
  .bulk_coerce_de_schema(de_res)
}

.bulk_prepare_de_for_downstream <- function(
  srt,
  DE_threshold,
  require_score = FALSE
) {
  de_results <- .bulk_get_de_results(srt)
  if (is.null(de_results) || nrow(de_results) == 0) {
    return(NULL)
  }
  de_results <- .bulk_filter_de_results(
    de_results = de_results,
    DE_threshold = DE_threshold
  )
  .bulk_prepare_de_for_pathway(
    de_results = de_results,
    require_score = require_score
  )
}

.bulk_make_comparison_label <- function(df) {
  if ("comparison" %in% colnames(df)) {
    out <- as.character(df$comparison)
    out[is.na(out) | !nzchar(out)] <- "All"
    return(out)
  }
  if (!"group1" %in% colnames(df) && !"group2" %in% colnames(df)) {
    return(rep("All", nrow(df)))
  }
  g1 <- if ("group1" %in% colnames(df)) as.character(df$group1) else rep("", nrow(df))
  g2 <- if ("group2" %in% colnames(df)) as.character(df$group2) else rep("", nrow(df))

  has_g2 <- !is.na(g2) & nzchar(g2)
  out <- ifelse(has_g2, g2, g1)
  has_multi_group1 <- length(unique(g1[!is.na(g1) & nzchar(g1)])) > 1
  if (has_multi_group1) {
    out <- ifelse(
      has_g2 & !is.na(g1) & nzchar(g1),
      paste0(g1, "::", g2),
      out
    )
  }
  out[is.na(out) | !nzchar(out)] <- "All"
  out
}

.bulk_compose_comparison_label <- function(group, comparison) {
  group <- as.character(group)
  comparison <- as.character(comparison)
  has_group <- !is.na(group) && nzchar(group) && !identical(group, "All")
  has_comparison <- !is.na(comparison) && nzchar(comparison)
  if (has_group && has_comparison) {
    return(paste0(group, "::", comparison))
  }
  if (has_comparison) {
    return(comparison)
  }
  if (has_group) {
    return(group)
  }
  "All"
}

.bulk_detect_mode <- function(srt = NULL, bulk_se = NULL) {
  has_srt <- !is.null(srt)
  has_bulk <- !is.null(bulk_se)
  if (isTRUE(has_srt) == isTRUE(has_bulk)) {
    .bulk_abort(
      "Provide exactly one of {.arg srt} or {.arg bulk_se}."
    )
  }
  if (has_srt) {
    if (!inherits(srt, "Seurat")) {
      .bulk_abort("{.arg srt} must be a {.cls Seurat} object.")
    }
    return("pseudobulk")
  }
  if (!inherits(bulk_se, "SummarizedExperiment")) {
    .bulk_abort("{.arg bulk_se} must be a {.cls SummarizedExperiment} object.")
  }
  "pure_bulk"
}

.bulk_validate_method_spec <- function(method) {
  if (is.null(method) || length(method) == 0) {
    .bulk_abort("{.arg method} must be a non-empty character vector.")
  }
  if (is.list(method)) {
    .bulk_abort(
      "{.arg method} must be a character vector such as {.val de_edgeR_qlf}, not a named module list."
    )
  }
  if (!is.character(method)) {
    .bulk_abort("{.arg method} must be a character vector.")
  }
  method <- method[!is.na(method) & nzchar(method)]
  if (length(method) == 0) {
    .bulk_abort("{.arg method} must contain at least one method name.")
  }

  normalized <- lapply(method, .bulk_canonical_method)
  canonical_names <- vapply(normalized, `[[`, character(1), "method")
  dup <- duplicated(canonical_names)
  if (any(dup)) {
    .bulk_abort(
      paste0(
        "Duplicated RunBulk method(s): ",
        paste(unique(canonical_names[dup]), collapse = ", ")
      )
    )
  }
  names(normalized) <- canonical_names
  normalized
}

.bulk_normalize_method_args <- function(method_args) {
  if (is.null(method_args)) {
    return(list())
  }
  if (!is.list(method_args)) {
    .bulk_abort("{.arg method_args} must be a named list.")
  }
  if (length(method_args) == 0) {
    return(list())
  }
  if (is.null(names(method_args)) || any(!nzchar(names(method_args)))) {
    .bulk_abort("{.arg method_args} must be named by RunBulk method.")
  }

  out <- list()
  for (nm in names(method_args)) {
    canonical <- .bulk_canonical_method(nm)$method
    val <- method_args[[nm]]
    if (!is.list(val)) {
      .bulk_abort(
        paste0("method_args[['", nm, "']] must be a list.")
      )
    }
    out[[canonical]] <- utils::modifyList(out[[canonical]] %||% list(), val)
  }
  out
}

.bulk_collect_method_args <- function(
  method_name,
  global_args = list(),
  module_args = list(),
  method_args = list()
) {
  if (is.null(global_args)) {
    global_args <- list()
  }
  if (is.null(module_args)) {
    module_args <- list()
  }
  specific_args <- method_args[[method_name]] %||% list()
  if (!is.list(global_args) || !is.list(module_args) || !is.list(specific_args)) {
    .bulk_abort("RunBulk tuning arguments must be lists.")
  }
  utils::modifyList(
    utils::modifyList(global_args, module_args),
    specific_args
  )
}

.bulk_check_r_packages <- function(packages, method_name) {
  err <- tryCatch(
    {
      check_r(packages, verbose = FALSE)
      NULL
    },
    error = function(e) e
  )
  if (inherits(err, "error")) {
    return(.bulk_method_failed(
      reason = paste0(
        "Package(s) required for ",
        method_name,
        " are missing or failed to install: ",
        err$message
      )
    ))
  }
  NULL
}

.bulk_canonical_method <- function(method_name) {
  key <- gsub("[^a-z0-9]", "", tolower(method_name))
  method_map <- list(
    delimmavoom = list(module = "de", method = "de_limma_voom"),
    limmavoom = list(module = "de", method = "de_limma_voom"),
    deedgerqlf = list(module = "de", method = "de_edgeR_qlf"),
    edgerqlf = list(module = "de", method = "de_edgeR_qlf"),
    edger = list(module = "de", method = "de_edgeR_qlf"),
    dedeseq2 = list(module = "de", method = "de_DESeq2"),
    deseq2 = list(module = "de", method = "de_DESeq2"),
    dedream = list(module = "de", method = "de_dream"),
    dream = list(module = "de", method = "de_dream"),
    deconvmusic = list(module = "deconv", method = "deconv_MuSiC"),
    music = list(module = "deconv", method = "deconv_MuSiC"),
    deconvbisquerna = list(module = "deconv", method = "deconv_BisqueRNA"),
    bisquerna = list(module = "deconv", method = "deconv_BisqueRNA"),
    deconvbayesprism = list(module = "deconv", method = "deconv_BayesPrism"),
    bayesprism = list(module = "deconv", method = "deconv_BayesPrism"),
    csdetoast = list(module = "csde", method = "csde_TOAST"),
    toast = list(module = "csde", method = "csde_TOAST")
  )

  out <- method_map[[key]]
  if (is.null(out)) {
    .bulk_abort(
      paste0(
        "Unsupported RunBulk method: ",
        method_name,
        ". Supported methods are de_limma_voom, de_edgeR_qlf, de_DESeq2, de_dream, ",
        "deconv_MuSiC, deconv_BisqueRNA, deconv_BayesPrism, and csde_TOAST."
      )
    )
  }
  out
}

.bulk_build_context <- function(
  mode,
  srt = NULL,
  bulk_se = NULL,
  sample.by = NULL,
  condition.by = NULL,
  group.by = NULL,
  assay = NULL,
  layer = "counts",
  bulk_assay = "counts"
) {
  if (identical(mode, "pseudobulk")) {
    return(.bulk_build_context_from_seurat(
      srt = srt,
      sample.by = sample.by,
      condition.by = condition.by,
      group.by = group.by,
      assay = assay,
      layer = layer
    ))
  }
  .bulk_build_context_from_bulk_se(
    bulk_se = bulk_se,
    condition.by = condition.by,
    group.by = group.by,
    bulk_assay = bulk_assay
  )
}

.bulk_build_context_from_seurat <- function(
  srt,
  sample.by,
  condition.by = NULL,
  group.by = NULL,
  assay = NULL,
  layer = "counts"
) {
  if (is.null(sample.by) || !sample.by %in% colnames(srt@meta.data)) {
    .bulk_abort(
      "{.arg sample.by} must be a valid metadata column in pseudobulk mode."
    )
  }
  if (!is.null(condition.by) && !condition.by %in% colnames(srt@meta.data)) {
    .bulk_abort(
      "{.arg condition.by} must be a valid metadata column in pseudobulk mode."
    )
  }
  if (!is.null(group.by) && !group.by %in% colnames(srt@meta.data)) {
    .bulk_abort(
      "{.arg group.by} must be a valid metadata column in pseudobulk mode."
    )
  }

  assay_use <- assay %||% SeuratObject::DefaultAssay(srt)
  counts <- GetAssayData5(
    object = srt,
    assay = assay_use,
    layer = layer
  )
  meta <- srt@meta.data
  sample_vec <- as.character(meta[[sample.by]])
  condition_col_use <- condition.by %||% ".bulk_condition"
  condition_vec <- if (is.null(condition.by)) {
    rep("All", length(sample_vec))
  } else {
    as.character(meta[[condition.by]])
  }
  valid <- !is.na(sample_vec) & !is.na(condition_vec)
  if (sum(valid) == 0) {
    .bulk_abort("No valid cells remained after filtering NA sample/condition labels.")
  }
  cells_use <- rownames(meta)[valid]
  counts <- counts[, cells_use, drop = FALSE]
  meta <- meta[cells_use, , drop = FALSE]
  sample_vec <- as.character(meta[[sample.by]])
  condition_vec <- if (is.null(condition.by)) {
    rep("All", length(sample_vec))
  } else {
    as.character(meta[[condition.by]])
  }

  pair_df <- unique(data.frame(
    sample = sample_vec,
    condition = condition_vec,
    stringsAsFactors = FALSE
  ))
  sample_dup <- pair_df$sample[duplicated(pair_df$sample)]
  if (length(sample_dup) > 0) {
    .bulk_abort(
      paste0(
        "Each sample must map to exactly one condition. Conflicting sample(s): ",
        paste(unique(sample_dup), collapse = ", ")
      )
    )
  }
  condition_map <- stats::setNames(pair_df$condition, pair_df$sample)

  counts_global <- .bulk_aggregate_counts(counts = counts, groups = sample_vec)
  condition_global <- condition_map[colnames(counts_global)]
  names(condition_global) <- colnames(counts_global)

  sample_meta_global <- .bulk_build_sample_metadata(
    meta = meta,
    sample_vec = sample_vec,
    sample_levels = colnames(counts_global),
    condition = condition_global,
    sample.by = sample.by,
    condition_col = condition_col_use
  )

  counts_by_group <- list()
  condition_by_group <- list()
  sample_meta_by_group <- list()
  if (is.null(group.by)) {
    counts_by_group[["All"]] <- counts_global
    condition_by_group[["All"]] <- condition_global
    sample_meta_by_group[["All"]] <- sample_meta_global
  } else {
    group_vec <- as.character(meta[[group.by]])
    group_levels <- unique(group_vec[!is.na(group_vec)])
    for (grp in group_levels) {
      idx <- which(group_vec %in% grp)
      if (length(idx) == 0) {
        next
      }
      counts_grp <- .bulk_aggregate_counts(
        counts = counts[, idx, drop = FALSE],
        groups = sample_vec[idx]
      )
      if (ncol(counts_grp) == 0) {
        next
      }
      cond_grp <- condition_map[colnames(counts_grp)]
      names(cond_grp) <- colnames(counts_grp)
      sample_meta_grp <- sample_meta_global[colnames(counts_grp), , drop = FALSE]
      counts_by_group[[grp]] <- counts_grp
      condition_by_group[[grp]] <- cond_grp
      sample_meta_by_group[[grp]] <- sample_meta_grp
    }
    if (length(counts_by_group) == 0) {
      .bulk_abort("No valid subgroup was found in {.arg group.by} for pseudobulk aggregation.")
    }
  }

  list(
    input_mode = "pseudobulk",
    counts_global = counts_global,
    condition_global = condition_global,
    counts_by_group = counts_by_group,
    condition_by_group = condition_by_group,
    sample_meta_global = sample_meta_global,
    sample_meta_by_group = sample_meta_by_group,
    group_levels = names(counts_by_group),
    assay = assay_use,
    layer = layer,
    bulk_assay = NULL,
    condition_col = condition_col_use
  )
}

.bulk_build_context_from_bulk_se <- function(
  bulk_se,
  condition.by = NULL,
  group.by = NULL,
  bulk_assay = "counts"
) {
  assay_names <- SummarizedExperiment::assayNames(bulk_se)
  if (!bulk_assay %in% assay_names) {
    .bulk_abort(
      paste0(
        "bulk_assay '",
        bulk_assay,
        "' not found in bulk_se. Available assays: ",
        paste(assay_names, collapse = ", ")
      )
    )
  }

  counts_global <- SummarizedExperiment::assay(bulk_se, bulk_assay)
  sample_meta_global <- as.data.frame(SummarizedExperiment::colData(bulk_se))
  condition_col_use <- condition.by %||% ".bulk_condition"
  if (!is.null(condition.by) && !condition.by %in% colnames(sample_meta_global)) {
    .bulk_abort(
      "{.arg condition.by} must be present in {.arg colData(bulk_se)}."
    )
  }
  if (is.null(rownames(sample_meta_global)) || any(rownames(sample_meta_global) == "")) {
    rownames(sample_meta_global) <- colnames(counts_global)
  }
  if (!all(colnames(counts_global) %in% rownames(sample_meta_global))) {
    .bulk_abort("Column names of bulk counts must align with row names of colData.")
  }
  sample_meta_global <- sample_meta_global[colnames(counts_global), , drop = FALSE]
  condition_global <- if (is.null(condition.by)) {
    rep("All", ncol(counts_global))
  } else {
    as.character(sample_meta_global[[condition.by]])
  }
  names(condition_global) <- colnames(counts_global)
  sample_meta_global[[condition_col_use]] <- condition_global

  counts_by_group <- list()
  condition_by_group <- list()
  sample_meta_by_group <- list()
  if (is.null(group.by)) {
    counts_by_group[["All"]] <- counts_global
    condition_by_group[["All"]] <- condition_global
    sample_meta_by_group[["All"]] <- sample_meta_global
  } else {
    if (!group.by %in% colnames(sample_meta_global)) {
      .bulk_abort(
        "{.arg group.by} must be present in {.arg colData(bulk_se)}."
      )
    }
    group_vec <- as.character(sample_meta_global[[group.by]])
    group_levels <- unique(group_vec[!is.na(group_vec)])
    for (grp in group_levels) {
      sample_idx <- which(group_vec %in% grp)
      if (length(sample_idx) == 0) {
        next
      }
      sample_names <- rownames(sample_meta_global)[sample_idx]
      counts_grp <- counts_global[, sample_names, drop = FALSE]
      cond_grp <- condition_global[sample_names]
      counts_by_group[[grp]] <- counts_grp
      condition_by_group[[grp]] <- cond_grp
      sample_meta_by_group[[grp]] <- sample_meta_global[sample_names, , drop = FALSE]
    }
    if (length(counts_by_group) == 0) {
      .bulk_abort("No valid subgroup was found in {.arg group.by} for bulk data.")
    }
  }

  list(
    input_mode = "pure_bulk",
    counts_global = counts_global,
    condition_global = condition_global,
    counts_by_group = counts_by_group,
    condition_by_group = condition_by_group,
    sample_meta_global = sample_meta_global,
    sample_meta_by_group = sample_meta_by_group,
    group_levels = names(counts_by_group),
    assay = NULL,
    layer = NULL,
    bulk_assay = bulk_assay,
    condition_col = condition_col_use
  )
}

.bulk_build_reference_profiles <- function(
  ref_srt,
  celltype.by,
  assay = NULL,
  layer = "counts"
) {
  if (!inherits(ref_srt, "Seurat")) {
    .bulk_abort("{.arg ref_srt} must be a {.cls Seurat} object.")
  }
  if (!celltype.by %in% colnames(ref_srt@meta.data)) {
    .bulk_abort(
      "{.arg celltype.by} must be present in {.arg ref_srt@meta.data}."
    )
  }

  assay_use <- assay %||% SeuratObject::DefaultAssay(ref_srt)
  ref_counts <- GetAssayData5(
    object = ref_srt,
    assay = assay_use,
    layer = layer
  )
  cell_type <- as.character(ref_srt@meta.data[[celltype.by]])
  valid <- !is.na(cell_type) & nzchar(cell_type)
  if (sum(valid) == 0) {
    .bulk_abort("No valid cell types were found for reference profile construction.")
  }
  ref_counts <- ref_counts[, valid, drop = FALSE]
  cell_type <- cell_type[valid]

  ref_profile <- .bulk_aggregate_counts(
    counts = ref_counts,
    groups = cell_type
  )
  ref_profile <- as.matrix(ref_profile)
  libsize <- Matrix::colSums(ref_profile)
  ref_profile <- t(t(ref_profile) / pmax(libsize, 1)) * 1e6

  list(
    profiles = ref_profile,
    cell_types = colnames(ref_profile),
    n_cells = table(cell_type),
    assay = assay_use,
    layer = layer
  )
}

.bulk_aggregate_counts <- function(counts, groups) {
  groups <- as.character(groups)
  valid <- !is.na(groups) & nzchar(groups)
  if (sum(valid) == 0) {
    return(matrix(
      0,
      nrow = nrow(counts),
      ncol = 0,
      dimnames = list(rownames(counts), NULL)
    ))
  }

  counts <- counts[, valid, drop = FALSE]
  groups <- groups[valid]
  index_list <- split(seq_along(groups), groups)
  agg <- lapply(index_list, function(i) {
    Matrix::rowSums(counts[, i, drop = FALSE])
  })
  out <- do.call(cbind, agg)
  if (is.null(dim(out))) {
    out <- matrix(out, ncol = 1)
  }
  rownames(out) <- rownames(counts)
  colnames(out) <- names(index_list)
  out
}

.bulk_build_sample_metadata <- function(
  meta,
  sample_vec,
  sample_levels,
  condition,
  sample.by,
  condition_col
) {
  constant_cols <- vapply(
    meta,
    .bulk_is_sample_constant,
    logical(1),
    sample_vec = sample_vec
  )
  sample_meta <- meta[
    match(sample_levels, sample_vec),
    names(constant_cols)[constant_cols],
    drop = FALSE
  ]
  rownames(sample_meta) <- sample_levels
  sample_meta[[sample.by]] <- sample_levels
  sample_meta[[condition_col]] <- as.character(condition[sample_levels])
  sample_meta
}

.bulk_is_sample_constant <- function(values, sample_vec) {
  tryCatch(
    {
      values_by_sample <- split(values, sample_vec)
      all(vapply(
        values_by_sample,
        function(x) {
          x <- x[!is.na(x)]
          length(unique(x)) <= 1
        },
        logical(1)
      ))
    },
    error = function(e) FALSE
  )
}

.bulk_resolve_condition_pair <- function(
  condition,
  condition1 = NULL,
  condition2 = NULL,
  strict_two_levels = TRUE
) {
  cond <- as.character(condition)
  cond <- cond[!is.na(cond)]
  levels_use <- unique(cond)
  if (length(levels_use) < 2) {
    .bulk_abort("At least two condition levels are required for comparison.")
  }

  if (is.null(condition1) || is.null(condition2)) {
    if (strict_two_levels && length(levels_use) != 2) {
      .bulk_abort(
        "Multiple condition levels were detected. Please set {.arg condition1} and {.arg condition2}."
      )
    }
    condition1 <- condition1 %||% levels_use[[1]]
    condition2 <- condition2 %||% levels_use[[2]]
  }

  if (identical(condition1, condition2)) {
    .bulk_abort("{.arg condition1} and {.arg condition2} must be different.")
  }
  if (!condition1 %in% levels_use || !condition2 %in% levels_use) {
    .bulk_abort(
      paste0(
        "condition1/condition2 must exist in condition levels: ",
        paste(levels_use, collapse = ", ")
      )
    )
  }
  list(condition1 = condition1, condition2 = condition2)
}

.bulk_fit_proportions <- function(
  count_matrix,
  reference_matrix,
  transform = c("log1p", "sqrt", "none"),
  min_overlap = 10
) {
  transform <- match.arg(transform)
  if (!is.numeric(min_overlap) || length(min_overlap) != 1 || is.na(min_overlap) || min_overlap < 1) {
    .bulk_abort("{.arg min_overlap} must be a positive number.")
  }
  genes <- intersect(rownames(count_matrix), rownames(reference_matrix))
  if (length(genes) < min_overlap) {
    .bulk_abort("Insufficient overlapping genes between bulk and reference matrices.")
  }

  bulk_use <- as.matrix(count_matrix[genes, , drop = FALSE])
  ref_use <- as.matrix(reference_matrix[genes, , drop = FALSE])
  if (identical(transform, "log1p")) {
    bulk_use <- log1p(bulk_use)
    ref_use <- log1p(ref_use)
  } else if (identical(transform, "sqrt")) {
    bulk_use <- sqrt(pmax(bulk_use, 0))
    ref_use <- sqrt(pmax(ref_use, 0))
  }

  prop <- matrix(
    0,
    nrow = ncol(bulk_use),
    ncol = ncol(ref_use),
    dimnames = list(colnames(bulk_use), colnames(ref_use))
  )
  x <- as.matrix(ref_use)
  for (i in seq_len(ncol(bulk_use))) {
    y <- as.numeric(bulk_use[, i])
    fit <- stats::lm.fit(x = x, y = y)
    coef_use <- fit$coefficients
    coef_use[is.na(coef_use)] <- 0
    coef_use[coef_use < 0] <- 0
    if (sum(coef_use) <= 0) {
      coef_use <- rep(1, length(coef_use))
    }
    prop[i, ] <- coef_use / sum(coef_use)
  }
  prop
}

.bulk_prop_matrix_to_long <- function(prop_matrix, method_name) {
  if (is.null(prop_matrix) || nrow(prop_matrix) == 0 || ncol(prop_matrix) == 0) {
    return(.bulk_empty_deconv_result())
  }
  df <- as.data.frame(as.table(prop_matrix), stringsAsFactors = FALSE)
  colnames(df) <- c("sample", "cell_type", "proportion")
  df$method <- method_name
  .bulk_coerce_deconv_schema(df)
}

.bulk_deconv_results_to_matrix <- function(deconv_results) {
  if (is.null(deconv_results) || nrow(deconv_results) == 0) {
    return(NULL)
  }
  df <- deconv_results[, c("sample", "cell_type", "proportion"), drop = FALSE]
  as.matrix(stats::xtabs(
    proportion ~ sample + cell_type,
    data = df
  ))
}

.bulk_coerce_de_schema <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(.bulk_empty_de_result())
  }
  if (!"gene" %in% colnames(df)) {
    df$gene <- rownames(df)
  }
  required <- c("group1", "group2", "avg_log2FC", "p_val", "p_val_adj", "method")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    for (nm in missing) {
      if (nm %in% c("avg_log2FC", "p_val", "p_val_adj")) {
        df[[nm]] <- NA_real_
      } else {
        df[[nm]] <- NA_character_
      }
    }
  }
  keep <- c("gene", "group1", "group2", "avg_log2FC", "p_val", "p_val_adj", "method")
  optional <- intersect(c("pct.1", "pct.2", "ave_expr", "comparison"), colnames(df))
  out <- df[, c(keep, optional), drop = FALSE]
  rownames(out) <- NULL
  out
}

.bulk_coerce_deconv_schema <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(.bulk_empty_deconv_result())
  }
  required <- c("sample", "cell_type", "proportion", "method")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    for (nm in missing) {
      if (identical(nm, "proportion")) {
        df[[nm]] <- NA_real_
      } else {
        df[[nm]] <- NA_character_
      }
    }
  }
  out <- df[, required, drop = FALSE]
  rownames(out) <- NULL
  out
}

.bulk_coerce_csde_schema <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(.bulk_empty_csde_result())
  }
  required <- c("gene", "cell_type", "group1", "group2", "effect", "p_val", "p_val_adj", "method")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    for (nm in missing) {
      if (nm %in% c("effect", "p_val", "p_val_adj")) {
        df[[nm]] <- NA_real_
      } else {
        df[[nm]] <- NA_character_
      }
    }
  }
  out <- df[, required, drop = FALSE]
  rownames(out) <- NULL
  out
}

.bulk_empty_de_result <- function() {
  data.frame(
    gene = character(),
    group1 = character(),
    group2 = character(),
    avg_log2FC = numeric(),
    p_val = numeric(),
    p_val_adj = numeric(),
    method = character(),
    stringsAsFactors = FALSE
  )
}

.bulk_empty_deconv_result <- function() {
  data.frame(
    sample = character(),
    cell_type = character(),
    proportion = numeric(),
    method = character(),
    stringsAsFactors = FALSE
  )
}

.bulk_empty_csde_result <- function() {
  data.frame(
    gene = character(),
    cell_type = character(),
    group1 = character(),
    group2 = character(),
    effect = numeric(),
    p_val = numeric(),
    p_val_adj = numeric(),
    method = character(),
    stringsAsFactors = FALSE
  )
}

.bulk_method_success <- function(
  results,
  details = list(),
  parameters = list(),
  reason = NULL
) {
  list(
    status = "success",
    reason = reason,
    results = results,
    details = details,
    parameters = parameters
  )
}

.bulk_method_failed <- function(
  reason,
  details = list(),
  parameters = list(),
  results = data.frame()
) {
  list(
    status = "failed",
    reason = reason,
    results = results,
    details = details,
    parameters = parameters
  )
}

.bulk_abort <- function(...) {
  log_message(..., message_type = "error")
}
