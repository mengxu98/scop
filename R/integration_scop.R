#' @title The integration workflow
#'
#' @description
#' Integrate single-cell RNA-seq data using various integration methods.
#' For `ChromatinAssay`, the current workflow uses `TFIDF + SVD/LSI` preprocessing.
#' In this setting, `Uncorrected` is supported directly and `Harmony5` will be
#' automatically redirected to the legacy `Harmony` workflow. `Seurat` and `RPCA`
#' are currently not supported for `ChromatinAssay`.
#'
#' @md
#' @inheritParams CheckDataList
#' @inheritParams CheckDataMerge
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param scale_within_batch Whether to scale data within each batch.
#' Only valid when the `integration_method` is one of `"Uncorrected"`,
#' `"Seurat"`, `"MNN"`, `"Harmony"`, `"BBKNN"`, `"CSS"`, `"ComBat"`.
#' @param integration_method A character vector specifying the integration method to use.
#' Supported methods are: `"Uncorrected"`, `"Seurat"`, `"CCA"`, `"RPCA"`, `"scVI"`,
#' `"PeakVI"`, `"PoissonVI"`, `"WNN"`, `"MultiMAP"`, `"GLUE"`, `"scVI5"`, `"MNN"`,
#' `"fastMNN"`, `"fastMNN5"`, `"Harmony"`, `"Harmony5"`,
#' `"Scanorama"`, `"BBKNN"`, `"CSS"`, `"Coralysis"`, `"LIGER"`, `"Conos"`, `"ComBat"`.
#' Default is `"Uncorrected"`. For `ChromatinAssay`, prefer `"Uncorrected"` or
#' `"Harmony5"`; the latter is automatically switched to `"Harmony"`.
#' @param compute_lisi Whether to compute LISI scores on the integrated result.
#' Default is `FALSE`.
#' @param lisi_label_colnames Character vector of metadata columns used to compute
#' LISI. If `NULL` and `compute_lisi = TRUE`, `batch` will be used when it is a
#' single metadata column name.
#' @param lisi_reduction Dimensional reduction used for LISI computation.
#' Default is `NULL`, which uses [DefaultReduction()] from the integrated object.
#' @param lisi_dims Dimensions used from `lisi_reduction`. Default is `NULL`,
#' which uses all available dimensions.
#' @param lisi_prefix Prefix used when storing LISI metadata columns.
#' Default is `NULL`, which uses `lisi_reduction`.
#' @param lisi_tool_name Name of the tool entry used to store LISI results.
#' Default is `NULL`, which uses `paste0(lisi_prefix, "_LISI")`.
#' @param lisi_perplexity Effective neighborhood size used by LISI.
#' Default is `30`.
#' @param lisi_tol Tolerance used in the LISI binary search. Default is `1e-5`.
#' @param lisi_max_iter Maximum iterations used in the LISI binary search.
#' Default is `50`.
#' @param compute_metrics Whether to compute integration summary metrics on the
#' selected reduction. Default is `FALSE`.
#' @param metrics_batch_col Metadata column used for batch-mixing metrics.
#' Default is `NULL`, which uses `batch` when it is a single metadata column name.
#' @param metrics_celltype_col Metadata column used for biological conservation
#' metrics. Default is `NULL`.
#' @param metrics_reduction Reduction used for integration metric computation.
#' Default is `NULL`, which uses [DefaultReduction()] from the integrated object.
#' @param metrics_cluster_col Metadata column used as cluster labels for
#' `celltype_NMI`, `celltype_ARI`, and `celltype_purity`. Default is `NULL`,
#' which resolves the integrated cluster column automatically when possible.
#' @param metrics_tool_name Name of the tool entry used to store integration
#' metrics. Default is `NULL`, which uses `paste0(integration_method, "_metrics")`.
#' @param metrics_k_graph Number of neighbors used for graph-connectivity
#' computation. Default is `15`.
#' @param append Whether the integrated data will be appended to the original Seurat object (`srt_merge`).
#' Default is `TRUE`.
#' @param ... Additional arguments to be passed to the integration method functions.
#'
#' @return A `Seurat` object.
#'
#' For `ChromatinAssay`, integrated outputs are additionally normalized to the
#' ATAC naming convention used in `scop`, including `*lsi`, `*UMAP2D`, cluster
#' aliases, and `ATAC_default_*` metadata when available.
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
#' [ComBat_integrate]
#'
#' @export
#' @examples
#' data(panc8_sub)
#' panc8_sub <- integration_scop(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_method = "Harmony",
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
#' LISIPlot(panc8_sub)
#'
#' panc8_sub <- integration_scop(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_method = "LIGER"
#' )
#' panc8_sub <- integration_scop(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_method = "Harmony",
#'   compute_lisi = TRUE,
#'   lisi_label_colnames = "tech"
#' )
#' LISIPlot(
#'   panc8_sub,
#'   features = c("HarmonypcaUMAP2D_tech_LISI", "HarmonyUMAP2D_tech_LISI")
#' )
#'
#' data("pbmcmultiome_sub", package = "scop")
#' pbmcmultiome_sub$batch <- rep(c("batch1", "batch2"), length.out = ncol(pbmcmultiome_sub))
#' pbmcmultiome_sub <- integration_scop(
#'   pbmcmultiome_sub,
#'   batch = "batch",
#'   assay = "peaks",
#'   integration_method = "Harmony5",
#'   normalization_method = "TFIDF"
#' )
#'
#' integration_methods <- c(
#'   "Uncorrected", "Seurat", "CCA", "RPCA",
#'   "MNN", "fastMNN", "fastMNN5", "Harmony", "Harmony5",
#'   "Scanorama", "BBKNN", "CSS", "Coralysis", "LIGER", "Conos", "ComBat"
#' )
#' p_list <- list()
#' for (method in integration_methods) {
#'   panc8_sub <- integration_scop(
#'     panc8_sub,
#'     batch = "tech",
#'     integration_method = method,
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
#' \dontrun{
#' # Python-backed methods prepare a scVI/scvi-tools environment and run model
#' # training, so keep them separate from ordinary example checks.
#' panc8_sub <- integration_scop(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_method = "scVI",
#'   train_params = list(max_epochs = 2L),
#'   nonlinear_reduction = "umap"
#' )
#' panc8_sub <- integration_scop(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_method = "scVI5",
#'   IntegrateLayers_params = list(max_epochs = 2L),
#'   nonlinear_reduction = "umap"
#' )
#' }
#'
#' nonlinear_reductions <- c(
#'   "umap", "tsne", "dm", "phate",
#'   "pacmap", "trimap", "largevis", "fr"
#' )
#' panc8_sub <- integration_scop(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_method = "Seurat",
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
integration_scop <- function(
  srt_merge = NULL,
  batch,
  append = TRUE,
  srt_list = NULL,
  assay = NULL,
  integration_method = c(
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
  ...
) {
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
  integration_method <- match.arg(integration_method)

  args <- as.list(match.call())[-1]
  new_env <- new.env(parent = parent.frame())
  args <- lapply(args, function(x) eval(x, envir = new_env))

  formals <- mget(names(formals()))
  formals <- formals[names(formals) != "..."]
  args <- utils::modifyList(formals, args)

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
      stop(
        "`integration_method = 'Seurat'` is not supported for `ChromatinAssay` in the current TFIDF/rlsi workflow. Please use `Uncorrected` or `Harmony5` (auto-switches to `Harmony`).",
        call. = FALSE
      )
    }
  }
  if (
    identical(integration_method, "RPCA") &&
      inherits(assay_source, "Seurat")
  ) {
    assay_use <- assay_requested %||% SeuratObject::DefaultAssay(assay_source)
    if (inherits(assay_source[[assay_use]], "ChromatinAssay")) {
      stop(
        "`integration_method = 'RPCA'` is not supported for `ChromatinAssay` in the current implementation.",
        call. = FALSE
      )
    }
  }
  if (
    identical(integration_method, "PeakVI") &&
      inherits(assay_source, "Seurat")
  ) {
    assay_use <- assay_requested %||% SeuratObject::DefaultAssay(assay_source)
    if (!inherits(assay_source[[assay_use]], "ChromatinAssay")) {
      stop(
        "`integration_method = 'PeakVI'` requires a `ChromatinAssay`.",
        call. = FALSE
      )
    }
  }
  if (
    identical(integration_method, "PoissonVI") &&
      inherits(assay_source, "Seurat")
  ) {
    assay_use <- assay_requested %||% SeuratObject::DefaultAssay(assay_source)
    if (!inherits(assay_source[[assay_use]], "ChromatinAssay")) {
      stop(
        "`integration_method = 'PoissonVI'` requires a `ChromatinAssay`.",
        call. = FALSE
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
      stop(
        "`integration_method = 'WNN'` requires both an RNA assay and a `ChromatinAssay` in the same Seurat object.",
        call. = FALSE
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
      stop(
        "`integration_method = 'MultiMAP'` requires both an RNA assay and a `ChromatinAssay` in the same Seurat object.",
        call. = FALSE
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
      stop(
        "`integration_method = 'GLUE'` requires both an RNA assay and a `ChromatinAssay` in the same Seurat object.",
        call. = FALSE
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
      srt = srt_integrated,
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
      srt = srt_integrated,
      reductions = lisi_reductions,
      dims = lisi_dims,
      label_colnames = lisi_label_colnames,
      prefix = lisi_prefix_use,
      tool_name = lisi_tool_name_use,
      perplexity = lisi_perplexity,
      tol = lisi_tol,
      max_iter = lisi_max_iter,
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

collect_integration_metrics <- function(
  srt,
  reduction,
  batch_col,
  celltype_col = NULL,
  cluster_col = NULL,
  lisi_tool_name = NULL,
  lisi_prefix = NULL,
  k_graph = 15
) {
  emb <- Seurat::Embeddings(srt, reduction = reduction)
  summary_list <- list()
  if (!is.null(batch_col) && batch_col %in% colnames(srt@meta.data)) {
    summary_list[["batch_ASW_mixing"]] <- metric_silhouette(
      embeddings = emb,
      labels = srt[[batch_col, drop = TRUE]],
      maximize = FALSE
    )
  }
  if (!is.null(celltype_col) && celltype_col %in% colnames(srt@meta.data)) {
    summary_list[["celltype_ASW"]] <- metric_silhouette(
      embeddings = emb,
      labels = srt[[celltype_col, drop = TRUE]],
      maximize = TRUE
    )
    summary_list[["celltype_graph_connectivity"]] <- metric_graph_connectivity(
      embeddings = emb,
      labels = srt[[celltype_col, drop = TRUE]],
      k = k_graph
    )
  }
  if (
    !is.null(cluster_col) &&
      cluster_col %in% colnames(srt@meta.data) &&
      !is.null(celltype_col) &&
      celltype_col %in% colnames(srt@meta.data)
  ) {
    metrics <- classification_metrics_compute(
      predicted = srt[[cluster_col, drop = TRUE]],
      truth = srt[[celltype_col, drop = TRUE]]
    )
    summary_list[["celltype_NMI"]] <- metrics[["nmi"]]
    summary_list[["celltype_ARI"]] <- metrics[["ari"]]
    summary_list[["celltype_purity"]] <- metrics[["purity"]]
  }
  if (!is.null(lisi_tool_name) && lisi_tool_name %in% names(srt@tools)) {
    lisi_res <- srt@tools[[lisi_tool_name]]
    if (!is.null(lisi_res$label_colnames) && !is.null(lisi_res$scores)) {
      for (label in lisi_res$label_colnames) {
        cols <- grep(
          paste0("_", label, "_LISI$"),
          colnames(lisi_res$scores),
          value = TRUE
        )
        if (!is.null(lisi_prefix)) {
          cols <- cols[grepl(
            paste0("(^|\\.)", make.names(lisi_prefix), "_"),
            cols
          )]
        }
        if (length(cols) > 0) {
          summary_list[[paste0(label, "_LISI_mean")]] <- mean(
            as.numeric(as.matrix(lisi_res$scores[, cols, drop = FALSE])),
            na.rm = TRUE
          )
        }
      }
    }
  }
  data.frame(
    metric = names(summary_list),
    value = unlist(summary_list),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

find_neighbors_and_clusters <- function(
  srt,
  reduction,
  dims_use,
  graph_prefix,
  graph_snn,
  cluster_colname,
  HVF,
  neighbor_metric,
  neighbor_k,
  cluster_algorithm,
  cluster_algorithm_index,
  cluster_resolution,
  run_find_neighbors = TRUE,
  verbose
) {
  neighbor_k_use <- min(as.integer(neighbor_k), max(1L, ncol(srt) - 1L))
  if (!identical(neighbor_k_use, as.integer(neighbor_k))) {
    log_message(
      "Adjust neighbor k from {.val {neighbor_k}} to {.val {neighbor_k_use}} for small-sample clustering",
      verbose = verbose
    )
  }
  srt <- tryCatch(
    {
      if (isTRUE(run_find_neighbors)) {
        srt <- Seurat::FindNeighbors(
          object = srt,
          reduction = reduction,
          dims = dims_use,
          annoy.metric = neighbor_metric,
          k.param = neighbor_k_use,
          graph.name = paste0(graph_prefix, c("KNN", "SNN")),
          verbose = FALSE
        )
      }

      log_message(
        "Perform {.fn Seurat::FindClusters} with {.val {cluster_algorithm}}",
        verbose = verbose
      )
      srt <- Seurat::FindClusters(
        object = srt,
        resolution = cluster_resolution,
        algorithm = cluster_algorithm_index,
        leiden_method = "igraph",
        graph.name = graph_snn,
        verbose = FALSE
      )
      log_message("Reorder clusters...")
      srt <- srt_reorder(
        srt,
        features = HVF,
        reorder_by = "seurat_clusters",
        layer = "data"
      )
      srt[["seurat_clusters"]] <- NULL
      srt[[cluster_colname]] <- SeuratObject::Idents(srt)
      srt
    },
    error = function(error) {
      err_msg <- conditionMessage(error)
      err_msg <- gsub("{", "{{", err_msg, fixed = TRUE)
      err_msg <- gsub("}", "}}", err_msg, fixed = TRUE)
      log_message(err_msg, message_type = "warning", verbose = verbose)
      log_message(
        "Error when performing {.fn Seurat::FindClusters}. Skip this step",
        message_type = "warning",
        verbose = verbose
      )
      srt
    }
  )

  return(srt)
}

resolve_linear_dims_use <- function(
  srt,
  reduction,
  linear_reduction_dims_use = NULL,
  normalization_method = "LogNormalize",
  reduction_method = NULL,
  verbose = FALSE
) {
  if (!is.null(linear_reduction_dims_use)) {
    return(linear_reduction_dims_use)
  }
  RunDimsEstimate(
    srt = srt,
    reduction = reduction,
    reduction_method = reduction_method,
    skip_first = normalization_method == "TFIDF",
    use_stored = TRUE,
    verbose = verbose
  )
}

run_nonlinear_reduction <- function(
  srt,
  prefix,
  reduction_use = NULL,
  reduction_dims = NULL,
  graph_use = NULL,
  neighbor_use = NULL,
  nonlinear_reduction,
  nonlinear_reduction_dims,
  nonlinear_reduction_params,
  force_nonlinear_reduction,
  seed,
  verbose
) {
  if (
    !is.null(reduction_use) && reduction_use %in% SeuratObject::Reductions(srt)
  ) {
    available_dims <- seq_len(
      ncol(Seurat::Embeddings(srt, reduction = reduction_use))
    )
    reduction_dims <- intersect(reduction_dims, available_dims)
    if (length(reduction_dims) == 0) {
      log_message(
        "No valid dimensions remain for {.arg reduction_use = '{reduction_use}'}",
        message_type = "warning",
        verbose = verbose
      )
      return(srt)
    }
  }
  srt <- tryCatch(
    {
      for (nr in nonlinear_reduction) {
        params_use <- nonlinear_reduction_params
        if (nr %in% c("fr")) {
          params_use[["n.neighbors"]] <- NULL
        }
        for (n in nonlinear_reduction_dims) {
          srt <- RunDimsReduction(
            srt,
            prefix = prefix,
            reduction_use = reduction_use,
            reduction_dims = reduction_dims,
            graph_use = graph_use,
            neighbor_use = neighbor_use,
            nonlinear_reduction = nr,
            nonlinear_reduction_dims = n,
            nonlinear_reduction_params = params_use,
            force_nonlinear_reduction = force_nonlinear_reduction,
            verbose = verbose,
            seed = seed
          )
        }
      }
      srt
    },
    error = function(error) {
      err_msg <- conditionMessage(error)
      err_msg <- gsub("{", "{{", err_msg, fixed = TRUE)
      err_msg <- gsub("}", "}}", err_msg, fixed = TRUE)
      log_message(err_msg, message_type = "warning", verbose = verbose)
      log_message(
        "Error when performing nonlinear dimension reduction. Skip this step",
        message_type = "warning",
        verbose = verbose
      )
      srt
    }
  )

  return(srt)
}
metric_graph_connectivity <- function(
  embeddings,
  labels,
  k = 15,
  backend = "r"
) {
  backend <- match.arg(backend, "r")

  embeddings <- as.matrix(embeddings)
  storage.mode(embeddings) <- "double"
  if (nrow(embeddings) != length(labels)) {
    log_message(
      "{.arg embeddings} rows must match {.arg labels} length",
      message_type = "error"
    )
  }

  labels <- as.factor(labels)
  keep <- !is.na(labels)
  embeddings <- embeddings[keep, , drop = FALSE]
  labels <- droplevels(labels[keep])
  if (nrow(embeddings) < 3 || nlevels(labels) < 1) {
    return(NA_real_)
  }

  k_use <- min(as.integer(k), nrow(embeddings) - 1L)
  if (k_use < 1) {
    return(NA_real_)
  }

  edges <- graph_conn_edges_r(
    embeddings = embeddings,
    k = k_use
  )

  graph_conn_score(
    edges = edges,
    labels = labels
  )
}

graph_conn_edges_from_index <- function(
  index,
  k,
  remove_self = TRUE
) {
  rows <- seq_len(nrow(index))
  nn <- lapply(rows, function(i) {
    idx <- as.integer(index[i, ])
    idx <- idx[!is.na(idx) & idx >= 1L]
    if (isTRUE(remove_self)) {
      idx <- idx[idx != i]
    }
    idx[seq_len(min(k, length(idx)))]
  })

  cbind(
    rep(rows, lengths(nn)),
    unlist(nn, use.names = FALSE)
  )
}

graph_conn_edges_r <- function(embeddings, k) {
  check_r("BiocNeighbors", verbose = FALSE)

  knn <- BiocNeighbors::findKNN(
    embeddings,
    k = k,
    BNPARAM = BiocNeighbors::KmknnParam(distance = "Euclidean"),
    num.threads = 1L
  )
  graph_conn_edges_from_index(
    index = knn$index,
    k = k,
    remove_self = FALSE
  )
}


graph_conn_score <- function(edges, labels) {
  if (is.null(edges) || nrow(edges) == 0L) {
    return(NA_real_)
  }

  graph <- igraph::graph_from_edgelist(edges, directed = FALSE)
  graph <- igraph::simplify(graph)
  per_label <- tapply(seq_along(labels), labels, function(idx) {
    if (length(idx) <= 1) {
      return(1)
    }
    comps <- igraph::components(igraph::induced_subgraph(graph, vids = idx))
    max(comps$csize) / length(idx)
  })
  mean(unlist(per_label), na.rm = TRUE)
}
metric_silhouette <- function(embeddings, labels, maximize = TRUE) {
  check_r("cluster", verbose = FALSE)
  labels <- as.factor(labels)
  keep <- !is.na(labels)
  embeddings <- embeddings[keep, , drop = FALSE]
  labels <- droplevels(labels[keep])
  if (nrow(embeddings) < 3 || nlevels(labels) < 2) {
    return(NA_real_)
  }
  sil <- cluster::silhouette(
    x = as.integer(labels),
    dist = stats::dist(embeddings)
  )
  score <- mean(sil[, "sil_width"], na.rm = TRUE)
  if (isTRUE(maximize)) {
    return(score)
  }
  1 - abs(score)
}
