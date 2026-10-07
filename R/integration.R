#' @title Integration workflow
#'
#' @description
#' Integrate single-cell data with one or more methods. For `ChromatinAssay`,
#' the workflow uses TFIDF + SVD/LSI; `Uncorrected` is supported directly,
#' `Harmony5` is redirected to `Harmony`, and `Seurat`/`RPCA` are not supported.
#'
#' @md
#' @inheritParams CheckDataList
#' @inheritParams CheckDataMerge
#' @inheritParams RunStandardWorkflow
#' @inheritParams thisutils::log_message
#' @inheritParams scop-params
#' @param scale_within_batch Scale within each batch. Only used by
#' `"Uncorrected"`, `"Seurat"`, `"MNN"`, `"Harmony"`, `"BBKNN"`, `"CSS"`, `"ComBat"`.
#' @param integration_methods Method(s) to run. Multiple methods require
#' `append = TRUE` and are applied sequentially. For `ChromatinAssay`, prefer
#' `"Uncorrected"` or `"Harmony5"`.
#' @param integration_method Deprecated alias of `integration_methods`.
#' @param compute_lisi,lisi_label_colnames,lisi_reduction,lisi_dims,lisi_prefix,lisi_tool_name,lisi_perplexity,lisi_tol,lisi_max_iter,lisi_knn_algorithm,lisi_cores,lisi_max_dense_bytes
#' LISI scores. `lisi_label_colnames = NULL` uses `batch` when it is a single
#' metadata column; `lisi_reduction = NULL` uses [DefaultReduction()].
#' @param compute_metrics,metrics_batch_col,metrics_celltype_col,metrics_reduction,metrics_cluster_col,metrics_tool_name,metrics_k_graph
#' Integration summary metrics on the selected reduction.
#' @param append Append integrated results to `srt_merge`.
#' @param ... Passed to the integration method functions.
#'
#' @return A `Seurat` object. For `ChromatinAssay`, names follow the ATAC
#' convention (`*lsi`, `*UMAP2D`, cluster aliases, `ATAC_default_*`).
#'
#' @seealso
#' [Seurat_integrate],
#' [scVI_integrate],
#' [MultiMAP_integrate],
#' [GLUE_integrate],
#' [MNN_integrate],
#' [fastMNN_integrate],
#' [Harmony_integrate],
#' [Scanorama_integrate],
#' [BBKNN_integrate],
#' [CSS_integrate],
#' [Coralysis_integrate],
#' [LIGER_integrate],
#' [Conos_integrate],
#' [ComBat_integrate],
#' [RunIntegrationBenchmark],
#' [IntegrationBenchmarkPlot]
#'
#' @export
#' @examples
#' data(panc8_sub)
#' panc8_sub <- RunIntegration(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_methods = "Harmony",
#'   nHVF = 500,
#'   linear_reduction_dims = 20,
#'   linear_reduction_dims_use = 1:10,
#'   nonlinear_reduction_dims = 2,
#'   compute_lisi = TRUE,
#'   lisi_label_colnames = "tech",
#'   lisi_perplexity = 10
#' )
#' CellDimPlot(
#'   panc8_sub,
#'   group.by = c("tech", "celltype"),
#'   reduction = "HarmonyUMAP2D"
#' )
#'
#' IntegrationBenchmarkPlot(panc8_sub, plot_type = "box")
#'
#' panc8_sub <- RunIntegration(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_methods = "LIGER"
#' )
#' panc8_sub <- RunIntegration(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_methods = "Harmony",
#'   compute_lisi = TRUE,
#'   lisi_label_colnames = "tech"
#' )
#' IntegrationBenchmarkPlot(panc8_sub, plot_type = "box")
#'
#' data("pbmcmultiome_sub", package = "scop")
#' pbmcmultiome_sub$batch <- rep(c("batch1", "batch2"), length.out = ncol(pbmcmultiome_sub))
#' pbmcmultiome_sub <- RunIntegration(
#'   pbmcmultiome_sub,
#'   batch = "batch",
#'   assay = "peaks",
#'   integration_methods = "Harmony5",
#'   normalization_method = "TFIDF"
#' )
#'
#' integration_methods <- c(
#'   "Uncorrected", "Seurat", "CCA", "RPCA",
#'   "MNN", "fastMNN", "Harmony", "Harmony5",
#'   "Coralysis", "LIGER", "Conos", "ComBat"
#' )
#' p_list <- list()
#' for (method in integration_methods) {
#'   panc8_sub <- RunIntegration(
#'     panc8_sub,
#'     batch = "tech",
#'     integration_methods = method,
#'     linear_reduction_dims_use = 1:50,
#'     nonlinear_reduction = "umap"
#'   )
#'   p_list[[method]] <- CellDimPlot(
#'     panc8_sub,
#'     group.by = c("tech", "celltype"),
#'     reduction = paste0(method, "UMAP2D"),
#'     xlab = "", ylab = "",
#'     title = method,
#'     legend.position = "none",
#'     theme_use = "theme_blank"
#'   )
#' }
#'
#' # Python-backed methods prepare a scVI/scvi-tools environment and run model
#' # training, so keep them separate from ordinary example checks.
#' \dontrun{
#' if (reticulate::py_module_available("scvi")) {
#'   panc8_sub <- RunIntegration(
#'     panc8_sub,
#'     batch = "tech",
#'     integration_methods = "scVI",
#'     train_params = list(max_epochs = 2L),
#'     nonlinear_reduction = "umap"
#'   )
#'   panc8_sub <- RunIntegration(
#'     panc8_sub,
#'     batch = "tech",
#'     integration_methods = "scVI5",
#'     IntegrateLayers_params = list(max_epochs = 2L),
#'     nonlinear_reduction = "umap"
#'   )
#' }
#' }
#'
#' nonlinear_reductions <- c(
#'   "umap", "tsne", "dm", "phate",
#'   "pacmap", "trimap", "largevis", "fr"
#' )
#' panc8_sub <- RunIntegration(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_methods = "Seurat",
#'   linear_reduction_dims_use = 1:50,
#'   nonlinear_reduction = nonlinear_reductions
#' )
#' for (nr in nonlinear_reductions) {
#'   print(
#'     CellDimPlot(
#'       panc8_sub,
#'       group.by = c("tech", "celltype"),
#'       reduction = paste0("Seurat", nr, "2D"),
#'       xlab = "", ylab = "", title = nr,
#'       legend.position = "none", theme_use = "theme_blank"
#'     )
#'   )
#' }
RunIntegration <- function(
  srt_merge = NULL,
  batch,
  append = TRUE,
  srt_list = NULL,
  assay = NULL,
  integration_methods = c(
    "Uncorrected",
    "Seurat",
    "CCA",
    "RPCA",
    "scVI",
    "PeakVI",
    "PoissonVI",
    "WNN",
    "MultiMAP",
    "GLUE",
    "scVI5",
    "MNN",
    "fastMNN",
    "fastMNN5",
    "Harmony",
    "Harmony5",
    "Scanorama",
    "BBKNN",
    "CSS",
    "Coralysis",
    "LIGER",
    "Conos",
    "ComBat"
  ),
  compute_lisi = FALSE,
  lisi_label_colnames = NULL,
  lisi_reduction = NULL,
  lisi_dims = NULL,
  lisi_prefix = NULL,
  lisi_tool_name = NULL,
  lisi_perplexity = 30,
  lisi_tol = 1e-5,
  lisi_max_iter = 50,
  lisi_knn_algorithm = c("auto", "brute_force", "clustered"),
  lisi_cores = NULL,
  lisi_max_dense_bytes = Inf,
  compute_metrics = FALSE,
  metrics_batch_col = NULL,
  metrics_celltype_col = NULL,
  metrics_reduction = NULL,
  metrics_cluster_col = NULL,
  metrics_tool_name = NULL,
  metrics_k_graph = 15,
  do_normalization = NULL,
  normalization_method = "LogNormalize",
  do_HVF_finding = TRUE,
  HVF_source = "separate",
  HVF_method = "vst",
  nHVF = 2000,
  HVF_min_intersection = 1,
  HVF = NULL,
  do_scaling = TRUE,
  vars_to_regress = NULL,
  regression_model = "linear",
  scale_within_batch = FALSE,
  linear_reduction = "pca",
  linear_reduction_dims = 50,
  linear_reduction_dims_use = NULL,
  linear_reduction_params = list(),
  force_linear_reduction = FALSE,
  nonlinear_reduction = "umap",
  nonlinear_reduction_dims = c(2, 3),
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  seed = 11,
  verbose = TRUE,
  integration_method = NULL,
  ...
) {
  integration_methods_missing <- missing(integration_methods)
  integration_method_missing <- missing(integration_method)

  if (!isTRUE(integration_method_missing)) {
    if (!isTRUE(integration_methods_missing)) {
      log_message(
        "Supply only one of `integration_methods` and the deprecated `integration_method`.",
        message_type = "error"
      )
    }
    .Deprecated(
      new = "integration_methods",
      package = "scop",
      msg = paste0(
        "`integration_method` is deprecated; use `integration_methods` instead. ",
        "It will be removed in scop 1.0.0."
      )
    )
    integration_methods <- integration_method
  }

  log_message(
    "Run integration workflow...",
    message_type = "running",
    text_color = "blue",
    verbose = verbose
  )

  if (is.null(srt_list) && is.null(srt_merge)) {
    log_message(
      "{.arg srt_list} or {.arg srt_merge} must be provided",
      message_type = "error"
    )
  }

  args <- as.list(match.call())[-1]
  new_env <- new.env(parent = parent.frame())
  args <- lapply(args, function(x) eval(x, envir = new_env))

  formals <- mget(names(formals()))
  formals <- formals[names(formals) != "..."]
  args <- utils::modifyList(formals, args)

  integration_method_choices <- eval(
    base::formals(RunIntegration)[["integration_methods"]]
  )
  integration_methods <- if (
    isTRUE(integration_methods_missing) &&
      isTRUE(integration_method_missing)
  ) {
    integration_method_choices[[1]]
  } else {
    match.arg(
      args[["integration_methods"]],
      choices = integration_method_choices,
      several.ok = TRUE
    )
  }
  args[["integration_method"]] <- NULL
  args[["integration_methods"]] <- integration_methods

  if (length(integration_methods) > 1L) {
    if (!isTRUE(args[["append"]])) {
      log_message(
        "Multiple {.arg integration_methods} values require {.arg append = TRUE}",
        message_type = "error"
      )
    }
    return(run_integration_methods(args, integration_methods))
  }

  integration_method <- integration_methods[[1]]
  args[["integration_methods"]] <- NULL
  args[["integration_method"]] <- integration_method

  assay_requested <- args[["assay"]] %||% NULL
  assay_source <- args[["srt_merge"]] %||% NULL
  if (
    is.null(assay_source) &&
      !is.null(args[["srt_list"]]) &&
      length(args[["srt_list"]]) > 0
  ) {
    assay_source <- args[["srt_list"]][[1]]
  }
  if (
    identical(integration_method, "Seurat") &&
      inherits(assay_source, "Seurat")
  ) {
    assay_use <- assay_requested %||% SeuratObject::DefaultAssay(assay_source)
    if (inherits(assay_source[[assay_use]], "ChromatinAssay")) {
      log_message(
        "`integration_method = 'Seurat'` is not supported for `ChromatinAssay` in the current TFIDF/rlsi workflow. Please use `Uncorrected` or `Harmony5` (auto-switches to `Harmony`).",
        message_type = "error"
      )
    }
  }
  if (
    identical(integration_method, "RPCA") &&
      inherits(assay_source, "Seurat")
  ) {
    assay_use <- assay_requested %||% SeuratObject::DefaultAssay(assay_source)
    if (inherits(assay_source[[assay_use]], "ChromatinAssay")) {
      log_message(
        "`integration_method = 'RPCA'` is not supported for `ChromatinAssay` in the current implementation.",
        message_type = "error"
      )
    }
  }
  if (
    identical(integration_method, "PeakVI") &&
      inherits(assay_source, "Seurat")
  ) {
    assay_use <- assay_requested %||% SeuratObject::DefaultAssay(assay_source)
    if (!inherits(assay_source[[assay_use]], "ChromatinAssay")) {
      log_message(
        "`integration_method = 'PeakVI'` requires a `ChromatinAssay`.",
        message_type = "error"
      )
    }
  }
  if (
    identical(integration_method, "PoissonVI") &&
      inherits(assay_source, "Seurat")
  ) {
    assay_use <- assay_requested %||% SeuratObject::DefaultAssay(assay_source)
    if (!inherits(assay_source[[assay_use]], "ChromatinAssay")) {
      log_message(
        "`integration_method = 'PoissonVI'` requires a `ChromatinAssay`.",
        message_type = "error"
      )
    }
  }
  if (
    identical(integration_method, "WNN") &&
      inherits(assay_source, "Seurat")
  ) {
    assays_available <- SeuratObject::Assays(assay_source)
    chrom_assays <- assays_available[vapply(
      assays_available,
      function(x) inherits(assay_source[[x]], "ChromatinAssay"),
      logical(1)
    )]
    rna_assays <- setdiff(assays_available, chrom_assays)
    if (length(chrom_assays) == 0 || length(rna_assays) == 0) {
      log_message(
        "`integration_method = 'WNN'` requires both an RNA assay and a `ChromatinAssay` in the same Seurat object.",
        message_type = "error"
      )
    }
  }
  if (
    identical(integration_method, "MultiMAP") &&
      inherits(assay_source, "Seurat")
  ) {
    assays_available <- SeuratObject::Assays(assay_source)
    chrom_assays <- assays_available[vapply(
      assays_available,
      function(x) inherits(assay_source[[x]], "ChromatinAssay"),
      logical(1)
    )]
    rna_assays <- setdiff(assays_available, chrom_assays)
    if (length(chrom_assays) == 0 || length(rna_assays) == 0) {
      log_message(
        "`integration_method = 'MultiMAP'` requires both an RNA assay and a `ChromatinAssay` in the same Seurat object.",
        message_type = "error"
      )
    }
  }
  if (
    identical(integration_method, "GLUE") &&
      inherits(assay_source, "Seurat")
  ) {
    assays_available <- SeuratObject::Assays(assay_source)
    chrom_assays <- assays_available[vapply(
      assays_available,
      function(x) inherits(assay_source[[x]], "ChromatinAssay"),
      logical(1)
    )]
    rna_assays <- setdiff(assays_available, chrom_assays)
    if (length(chrom_assays) == 0 || length(rna_assays) == 0) {
      log_message(
        "`integration_method = 'GLUE'` requires both an RNA assay and a `ChromatinAssay` in the same Seurat object.",
        message_type = "error"
      )
    }
  }
  if (
    identical(integration_method, "Harmony5") &&
      inherits(assay_source, "Seurat")
  ) {
    assay_use <- assay_requested %||% SeuratObject::DefaultAssay(assay_source)
    if (inherits(assay_source[[assay_use]], "ChromatinAssay")) {
      log_message(
        "{.arg integration_method = 'Harmony5'} is not compatible with {.cls ChromatinAssay} in current Seurat v5 workflow. Automatically switch to {.val Harmony}",
        message_type = "warning",
        verbose = verbose
      )
      integration_method <- "Harmony"
      args[["integration_method"]] <- "Harmony"
    }
  }

  method_map <- list(
    Uncorrected = Uncorrected_integrate,
    Seurat = Seurat_integrate,
    CCA = CCA_integrate,
    RPCA = RPCA_integrate,
    scVI = scVI_integrate,
    PeakVI = scVI_integrate,
    PoissonVI = scVI_integrate,
    WNN = WNN_integrate,
    MultiMAP = MultiMAP_integrate,
    GLUE = GLUE_integrate,
    scVI5 = scVI5_integrate,
    MNN = MNN_integrate,
    fastMNN = fastMNN_integrate,
    fastMNN5 = fastMNN5_integrate,
    Harmony = Harmony_integrate,
    Harmony5 = Harmony5_integrate,
    Scanorama = Scanorama_integrate,
    BBKNN = BBKNN_integrate,
    CSS = CSS_integrate,
    Coralysis = Coralysis_integrate,
    LIGER = LIGER_integrate,
    Conos = Conos_integrate,
    ComBat = ComBat_integrate
  )

  assay_use <- args[["assay"]] %||%
    (if (!is.null(args[["srt_merge"]])) {
      SeuratObject::DefaultAssay(args[["srt_merge"]])
    } else if (!is.null(args[["srt_list"]]) && length(args[["srt_list"]]) > 0) {
      SeuratObject::DefaultAssay(args[["srt_list"]][[1]])
    })
  is_assay5 <- !is.null(assay_use) &&
    inherits(assay_source, "Seurat") &&
    inherits(Seurat::GetAssay(assay_source, assay = assay_use), "Assay5")
  v5_only_methods <- c("CCA", "RPCA", "fastMNN5", "Harmony5", "scVI5")
  if (integration_method %in% v5_only_methods && !isTRUE(is_assay5)) {
    log_message(
      "{.arg integration_method = '{integration_method}'} requires an {.cls Assay5} (Seurat v5) assay, but the assay {.val {assay_use}} is not Assay5. Please upgrade to Seurat v5 or use an alternative integration method.",
      message_type = "error"
    )
  }

  integrate_fun <- method_map[[integration_method]]
  append_requested <- isTRUE(args[["append"]])
  srt_merge_raw <- args[["srt_merge"]] %||% NULL
  if (
    append_requested &&
      !is.null(srt_merge_raw) &&
      inherits(srt_merge_raw, "Seurat") &&
      !integration_method %in% c("WNN", "MultiMAP", "GLUE")
  ) {
    canonical_linear_reductions <- c("pca", "svd", "ica", "nmf", "mds", "glmpca")
    linear_reduction_arg <- args[["linear_reduction"]] %||% "pca"
    can_drop_reductions <- all(linear_reduction_arg %in% canonical_linear_reductions)
    assay_keep <- args[["assay"]] %||% SeuratObject::DefaultAssay(srt_merge_raw)
    if (length(assay_keep) == 1L && assay_keep %in% SeuratObject::Assays(srt_merge_raw)) {
      SeuratObject::DefaultAssay(srt_merge_raw) <- assay_keep
      args[["srt_merge"]] <- Seurat::DietSeurat(
        object = srt_merge_raw,
        assays = assay_keep,
        dimreducs = if (isTRUE(can_drop_reductions)) NULL else SeuratObject::Reductions(srt_merge_raw),
        graphs = NULL,
        misc = FALSE
      )
      SeuratObject::DefaultAssay(args[["srt_merge"]]) <- assay_keep
    }
  }
  if (identical(integration_method, "PeakVI")) {
    args[["model"]] <- "PEAKVI"
  }
  if (identical(integration_method, "PoissonVI")) {
    args[["model"]] <- "POISSONVI"
  }
  if (
    "append" %in% names(args) && "append" %in% names(formals(integrate_fun))
  ) {
    args[["append"]] <- FALSE
  }
  srt_integrated <- invoke_fun(
    integrate_fun,
    args[names(args) %in% names(formals(integrate_fun))]
  )
  if (length(batch) == 1 && batch %in% colnames(srt_integrated@meta.data)) {
    srt_integrated@misc[["integration_batch"]] <- batch
  }
  if (
    inherits(
      srt_integrated[[SeuratObject::DefaultAssay(srt_integrated)]],
      "ChromatinAssay"
    )
  ) {
    srt_integrated <- standardize_atac(
      srt = srt_integrated,
      prefix = integration_method
    )
    reduction_linear_name <- tryCatch(
      DefaultReduction(
        srt_integrated,
        pattern = paste0("^", integration_method, "(lsi|svd|Harmony|Harmony5)$")
      ),
      error = function(...) character(0)
    )
    expected_umap <- paste0(integration_method, "UMAP2D")
    if (
      length(reduction_linear_name) == 1 &&
        nzchar(reduction_linear_name) &&
        !expected_umap %in% names(srt_integrated@reductions)
    ) {
      dims_fallback <- tryCatch(
        resolve_linear_dims_use(
          srt = srt_integrated,
          reduction = reduction_linear_name,
          linear_reduction_dims_use = linear_reduction_dims_use,
          normalization_method = normalization_method,
          reduction_method = "svd",
          verbose = FALSE
        ),
        error = function(...) {
          seq_len(min(
            10L,
            ncol(Seurat::Embeddings(
              srt_integrated,
              reduction = reduction_linear_name
            ))
          ))
        }
      )
      srt_integrated <- run_nonlinear_reduction(
        srt = srt_integrated,
        prefix = integration_method,
        reduction_use = reduction_linear_name,
        reduction_dims = dims_fallback,
        graph_use = NULL,
        nonlinear_reduction = "umap",
        nonlinear_reduction_dims = 2L,
        nonlinear_reduction_params = nonlinear_reduction_params,
        force_nonlinear_reduction = FALSE,
        seed = seed,
        verbose = verbose
      )
      srt_integrated <- standardize_atac(
        srt = srt_integrated,
        prefix = integration_method
      )
    }
  }

  integrated_default_reduction <- tryCatch(
    DefaultReduction(srt_integrated),
    error = function(...) character(0)
  )
  baseline_reduction <- character(0)
  pca_reduction <- tryCatch(
    DefaultReduction(srt_integrated, pattern = "pca"),
    error = function(...) character(0)
  )
  if (length(pca_reduction) == 1 && nzchar(pca_reduction)) {
    pca_dims_use <- tryCatch(
      resolve_linear_dims_use(
        srt = srt_integrated,
        reduction = pca_reduction,
        linear_reduction_dims_use = linear_reduction_dims_use,
        normalization_method = normalization_method,
        reduction_method = "pca",
        verbose = FALSE
      ),
      error = function(...) {
        seq_len(min(
          30L,
          ncol(Seurat::Embeddings(
            srt_integrated,
            reduction = pca_reduction
          ))
        ))
      }
    )
    baseline_nr <- nonlinear_reduction[[1]] %||% "umap"
    baseline_nr_dim <- 2L
    baseline_prefix <- pca_reduction
    srt_integrated <- RunDimsReduction(
      object = srt_integrated,
      prefix = baseline_prefix,
      reduction_use = pca_reduction,
      reduction_dims = pca_dims_use,
      nonlinear_reduction = baseline_nr,
      nonlinear_reduction_dims = baseline_nr_dim,
      nonlinear_reduction_params = nonlinear_reduction_params,
      force_nonlinear_reduction = FALSE,
      verbose = verbose,
      seed = seed
    )
    baseline_reduction <- paste0(
      baseline_prefix,
      toupper(gsub("-.*", "", baseline_nr)),
      baseline_nr_dim,
      "D"
    )
  }

  lisi_tool_name_use <- NULL
  lisi_prefix_map <- NULL
  if (isTRUE(compute_lisi)) {
    if (is.null(lisi_label_colnames)) {
      if (length(batch) == 1 && batch %in% colnames(srt_integrated@meta.data)) {
        lisi_label_colnames <- batch
      } else {
        log_message(
          "{.arg lisi_label_colnames} must be provided when {.arg batch} is not a single metadata column name",
          message_type = "error"
        )
      }
    }

    if (is.null(lisi_reduction)) {
      lisi_reductions <- unique(c(
        baseline_reduction,
        integrated_default_reduction
      ))
    } else {
      lisi_reductions <- unique(as.character(lisi_reduction))
    }
    lisi_reductions <- lisi_reductions[nzchar(lisi_reductions)]
    lisi_prefix_use <- lisi_prefix %||% lisi_reductions
    if (length(lisi_prefix_use) == 1 && length(lisi_reductions) > 1) {
      lisi_prefix_use <- rep(lisi_prefix_use, length(lisi_reductions))
    }
    lisi_tool_name_use <- lisi_tool_name %||%
      if (length(lisi_reductions) > 1) {
        "LISI"
      } else {
        paste0(lisi_prefix_use[[1]], "_LISI")
      }
    lisi_prefix_map <- stats::setNames(lisi_prefix_use, lisi_reductions)

    srt_integrated <- RunLISI(
      object = srt_integrated,
      reductions = lisi_reductions,
      dims = lisi_dims,
      label_colnames = lisi_label_colnames,
      prefix = lisi_prefix_use,
      tool_name = lisi_tool_name_use,
      perplexity = lisi_perplexity,
      tol = lisi_tol,
      max_iter = lisi_max_iter,
      knn_algorithm = lisi_knn_algorithm,
      cores = lisi_cores,
      max_dense_bytes = lisi_max_dense_bytes,
      verbose = verbose
    )
  }

  if (isTRUE(compute_metrics)) {
    metrics_batch_col <- metrics_batch_col %||%
      if (length(batch) == 1 && batch %in% colnames(srt_integrated@meta.data)) {
        batch
      } else {
        NULL
      }
    if (
      is.null(metrics_batch_col) &&
        is.null(metrics_celltype_col)
    ) {
      log_message(
        "At least one of {.arg metrics_batch_col} or {.arg metrics_celltype_col} must be available when {.arg compute_metrics = TRUE}",
        message_type = "error"
      )
    }
    if (
      !is.null(metrics_batch_col) &&
        !metrics_batch_col %in% colnames(srt_integrated@meta.data)
    ) {
      log_message(
        "{.arg metrics_batch_col} must be present in {.arg srt_integrated@meta.data}",
        message_type = "error"
      )
    }
    if (
      !is.null(metrics_celltype_col) &&
        !metrics_celltype_col %in% colnames(srt_integrated@meta.data)
    ) {
      log_message(
        "{.arg metrics_celltype_col} must be present in {.arg srt_integrated@meta.data}",
        message_type = "error"
      )
    }
    metrics_reduction_use <- metrics_reduction %||% integrated_default_reduction
    if (
      is.null(metrics_reduction_use) ||
        !metrics_reduction_use %in% SeuratObject::Reductions(srt_integrated)
    ) {
      log_message(
        "{.arg metrics_reduction} must refer to an existing reduction in the integrated object",
        message_type = "error"
      )
    }
    if (is.null(metrics_cluster_col)) {
      cluster_candidates <- unique(c(
        srt_integrated@misc[["ATAC_default_cluster_col"]] %||% NULL,
        paste0(integration_method, "clusters"),
        paste0(integration_method, linear_reduction, "clusters"),
        sub("UMAP2D$", "clusters", metrics_reduction_use),
        sub("UMAP3D$", "clusters", metrics_reduction_use),
        sub("TSNE2D$", "clusters", metrics_reduction_use),
        sub("DM2D$", "clusters", metrics_reduction_use),
        sub("PHATE2D$", "clusters", metrics_reduction_use),
        sub("PACMAP2D$", "clusters", metrics_reduction_use),
        sub("TRIMAP2D$", "clusters", metrics_reduction_use),
        sub("LARGEVIS2D$", "clusters", metrics_reduction_use),
        sub("FR2D$", "clusters", metrics_reduction_use)
      ))
      cluster_candidates <- cluster_candidates[
        cluster_candidates %in% colnames(srt_integrated@meta.data)
      ]
      metrics_cluster_col <- cluster_candidates[[1]] %||% NULL
    }
    metrics_tool_name <- metrics_tool_name %||%
      paste0(integration_method, "_metrics")
    metrics_lisi_prefix <- NULL
    if (
      !is.null(lisi_prefix_map) &&
        metrics_reduction_use %in% names(lisi_prefix_map)
    ) {
      metrics_lisi_prefix <- lisi_prefix_map[[metrics_reduction_use]]
    }
    metrics_summary <- collect_integration_metrics(
      srt = srt_integrated,
      reduction = metrics_reduction_use,
      batch_col = metrics_batch_col,
      celltype_col = metrics_celltype_col,
      cluster_col = metrics_cluster_col,
      lisi_tool_name = lisi_tool_name_use,
      lisi_prefix = metrics_lisi_prefix,
      k_graph = metrics_k_graph
    )
    srt_integrated@tools[[metrics_tool_name]] <- list(
      summary = metrics_summary,
      reduction = metrics_reduction_use,
      batch_col = metrics_batch_col,
      celltype_col = metrics_celltype_col,
      cluster_col = metrics_cluster_col,
      k_graph = metrics_k_graph,
      integration_method = integration_method,
      lisi_tool_name = lisi_tool_name_use,
      lisi_prefix = metrics_lisi_prefix
    )
  }

  if (isTRUE(append_requested) && !is.null(srt_merge_raw)) {
    assay_use <- assay %||% SeuratObject::DefaultAssay(srt_integrated)
    append_slots <- methods::slotNames(srt_integrated)
    if (inherits(srt_integrated[[assay_use]], "ChromatinAssay")) {
      append_pattern <- paste0(
        integration_method,
        "|pca|PCA|svd|SVD|lsi|LSI|UMAP2D|clusters|Default_reduction|LISI|integration_batch|ATAC_default_linear_reduction|ATAC_default_cluster_col"
      )
      append_slots <- intersect(
        append_slots,
        c("reductions", "meta.data", "misc", "tools")
      )
    } else {
      append_pattern <- paste0(
        assay_use,
        "|",
        integration_method,
        "|pca|PCA|svd|SVD|lsi|LSI|UMAP2D|clusters|Default_reduction|LISI|integration_batch|ATAC_default_linear_reduction|ATAC_default_cluster_col"
      )
    }
    srt_integrated <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_integrated,
      slots = append_slots,
      pattern = append_pattern,
      overwrite = TRUE,
      verbose = FALSE
    )
  }

  log_message(
    "{.pkg {integration_method}} integration completed",
    message_type = "success",
    text_color = "green",
    verbose = verbose
  )

  return(srt_integrated)
}


#' @title Uncorrected integration function
#'
#' @inheritParams RunIntegration
#'
#' @export
Uncorrected_integrate <- function(
  srt_merge = NULL,
  batch = NULL,
  append = TRUE,
  srt_list = NULL,
  assay = NULL,
  do_normalization = NULL,
  normalization_method = "LogNormalize",
  do_HVF_finding = TRUE,
  HVF_source = "separate",
  HVF_method = "vst",
  nHVF = 2000,
  HVF_min_intersection = 1,
  HVF = NULL,
  do_scaling = TRUE,
  vars_to_regress = NULL,
  regression_model = "linear",
  scale_within_batch = FALSE,
  linear_reduction = "pca",
  linear_reduction_dims = 50,
  linear_reduction_dims_use = NULL,
  linear_reduction_params = list(),
  force_linear_reduction = FALSE,
  nonlinear_reduction = "umap",
  nonlinear_reduction_dims = c(2, 3),
  nonlinear_reduction_params = list(),
  force_nonlinear_reduction = TRUE,
  neighbor_metric = "euclidean",
  neighbor_k = 20L,
  cluster_algorithm = "louvain",
  cluster_resolution = 0.6,
  verbose = TRUE,
  seed = 11
) {
  if (length(linear_reduction) > 1) {
    log_message(
      "Only the first method in the {.arg linear_reduction} will be used",
      message_type = "warning",
      verbose = verbose
    )
    linear_reduction <- linear_reduction[1]
  }
  reduc_test <- c("pca", "ica", "nmf", "mds", "glmpca")
  if (!is.null(srt_merge)) {
    reduc_test <- c(reduc_test, SeuratObject::Reductions(srt_merge))
  }
  if (any(!linear_reduction %in% reduc_test)) {
    log_message(
      "{.arg linear_reduction} must be one of {.val {reduc_test}}",
      message_type = "error"
    )
  }
  if (
    !is.null(linear_reduction_dims_use) &&
      max(linear_reduction_dims_use) > linear_reduction_dims
  ) {
    linear_reduction_dims <- max(linear_reduction_dims_use)
  }

  validate_nonlinear_reductions(nonlinear_reduction)
  cluster_algorithm_index <- resolve_cluster_algorithm_index(cluster_algorithm)

  set.seed(seed)
  validate_integration_input_cells(srt_list, srt_merge)
  srt_merge_raw <- srt_merge
  if (!is.null(srt_list)) {
    checked <- CheckDataList(
      srt_list = srt_list,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      verbose = verbose,
      seed = seed
    )
    srt_list <- checked[["srt_list"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
    srt_merge <- Reduce(merge, srt_list)
    SeuratObject::VariableFeatures(srt_merge) <- HVF
  }
  if (is.null(srt_list) && !is.null(srt_merge)) {
    checked <- CheckDataMerge(
      srt_merge = srt_merge,
      batch = batch,
      assay = assay,
      do_normalization = do_normalization,
      do_HVF_finding = do_HVF_finding,
      normalization_method = normalization_method,
      HVF_source = HVF_source,
      HVF_method = HVF_method,
      nHVF = nHVF,
      HVF_min_intersection = HVF_min_intersection,
      HVF = HVF,
      vars_to_regress = vars_to_regress,
      verbose = verbose,
      seed = seed
    )
    srt_merge <- checked[["srt_merge"]]
    HVF <- checked[["HVF"]]
    assay <- checked[["assay"]]
    type <- checked[["type"]]
  }

  if (normalization_method == "TFIDF") {
    log_message(
      "{.arg normalization_method} is {.val TFIDF}. Use {.pkg lsi} workflow..."
    )
    do_scaling <- FALSE
    linear_reduction <- "svd"
  }

  log_message(
    "Perform {.pkg Uncorrected} integration"
  )
  scale_features <- rownames(
    GetAssayData5(
      srt_merge,
      layer = "scale.data",
      assay = SeuratObject::DefaultAssay(srt_merge)
    )
  )
  if (
    isTRUE(do_scaling) || (is.null(do_scaling) && any(!HVF %in% scale_features))
  ) {
    if (normalization_method != "SCT") {
      log_message(
        "Perform {.fn Seurat::ScaleData}",
        verbose = verbose
      )
      assay_merge <- SeuratObject::DefaultAssay(srt_merge)
      if (inherits(srt_merge[[assay_merge]], "Assay5")) {
        srt_merge[[assay_merge]] <- SeuratObject::JoinLayers(
          srt_merge[[assay_merge]]
        )
      }
      srt_merge <- ScaleData(
        object = srt_merge,
        split.by = if (isTRUE(scale_within_batch)) batch else NULL,
        assay = SeuratObject::DefaultAssay(srt_merge),
        features = HVF,
        vars.to.regress = vars_to_regress,
        model.use = regression_model,
        verbose = FALSE
      )
    }
  }

  log_message(
    "Perform {.val {linear_reduction}} linear dimension reduction",
    verbose = verbose
  )
  srt_merge <- RunDimsReduction(
    srt_merge,
    prefix = "Uncorrected",
    features = HVF,
    assay = SeuratObject::DefaultAssay(srt_merge),
    linear_reduction = linear_reduction,
    linear_reduction_dims = linear_reduction_dims,
    linear_reduction_params = linear_reduction_params,
    force_linear_reduction = force_linear_reduction,
    verbose = verbose,
    seed = seed
  )
  if (is.null(linear_reduction_dims_use)) {
    linear_reduction_dims_use <- resolve_linear_dims_use(
      srt = srt_merge,
      reduction = paste0("Uncorrected", linear_reduction),
      normalization_method = normalization_method,
      reduction_method = linear_reduction
    )
  }

  srt_merge <- find_neighbors_and_clusters(
    srt = srt_merge,
    reduction = paste0("Uncorrected", linear_reduction),
    dims_use = linear_reduction_dims_use,
    graph_prefix = "Uncorrected_",
    graph_snn = "Uncorrected_SNN",
    cluster_colname = paste0("Uncorrected", linear_reduction, "clusters"),
    HVF = HVF,
    neighbor_metric = neighbor_metric,
    neighbor_k = neighbor_k,
    cluster_algorithm = cluster_algorithm,
    cluster_algorithm_index = cluster_algorithm_index,
    cluster_resolution = cluster_resolution,
    verbose = verbose
  )

  srt_merge <- run_nonlinear_reduction(
    srt = srt_merge,
    prefix = "Uncorrected",
    reduction_use = paste0("Uncorrected", linear_reduction),
    reduction_dims = linear_reduction_dims_use,
    graph_use = "Uncorrected_SNN",
    nonlinear_reduction = nonlinear_reduction,
    nonlinear_reduction_dims = nonlinear_reduction_dims,
    nonlinear_reduction_params = nonlinear_reduction_params,
    force_nonlinear_reduction = force_nonlinear_reduction,
    seed = seed,
    verbose = verbose
  )

  SeuratObject::DefaultAssay(srt_merge) <- assay
  SeuratObject::VariableFeatures(srt_merge) <- srt_merge@misc[[
    "Uncorrected_HVF"
  ]] <- HVF

  if (isTRUE(append) && !is.null(srt_merge_raw)) {
    srt_merge_raw <- srt_append(
      srt_raw = srt_merge_raw,
      srt_append = srt_merge,
      pattern = paste0(assay, "|Uncorrected|Default_reduction"),
      overwrite = TRUE,
      verbose = FALSE
    )
    return(srt_merge_raw)
  } else {
    return(srt_merge)
  }
}
