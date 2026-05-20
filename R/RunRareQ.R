#' @title RareQ rare-cell population detection
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param reduction Reduction used to build nearest neighbors when the required
#' `{assay}.nn` neighbor slot is absent or `force_recalc = TRUE`. If `NULL`,
#' [DefaultReduction()] is used.
#' @param dims Dimensions from `reduction` used for nearest-neighbor search.
#' @param k.param Number of nearest neighbors to compute with
#' [Seurat::FindNeighbors()] when neighbor search is needed.
#' @param k Number of nearest neighbors used by RareQ to compute Q values.
#' @param Q_cut Q-value threshold passed to `RareQ::FindRare()`.
#' @param ratio Merge-ratio threshold passed to `RareQ::FindRare()`.
#' @param max_iter Maximum number of RareQ propagation iterations.
#' @param run_neighbors Whether to build the required Seurat neighbor slot if
#' it is missing.
#' @param force_recalc Whether to rebuild the Seurat neighbor slot before
#' running RareQ.
#' @param find_neighbors_params Additional named parameters passed to
#' [Seurat::FindNeighbors()] when neighbor search is run.
#' @param rare_threshold Cluster-size threshold used to mark rare clusters. A
#' value smaller than 1 is treated as a fraction of cells; a value of 1 or
#' larger is treated as a cell count. Set to `NULL` to skip rare flags.
#' @param prefix Prefix used for metadata columns.
#' @param cluster_colname,q_colname,size_colname,rare_colname Metadata column
#' names for RareQ clusters, Q values, cluster sizes, and rare-cluster flags.
#' @param tool_name Name of the `srt@tools` entry.
#'
#' @return A `Seurat` object with RareQ results in metadata and
#' `srt@tools[[tool_name]]`.
#' @export
#'
#' @references
#' Fa, B. et al. Cell neighborhood topology directs rare cell population
#' identification. \emph{Nature Communications} (2026).
#' \doi{10.1038/s41467-026-71180-x}
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(
#'   pancreas_sub,
#'   verbose = FALSE
#' )
#' pancreas_sub <- RunRareQ(
#'   pancreas_sub,
#'   dims = 1:20
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "RareQ_cluster"
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "RareQ_is_rare"
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   features = "RareQ_Q"
#' )
RunRareQ <- function(
  srt,
  assay = NULL,
  reduction = "pca",
  dims = 1:30,
  k.param = 20,
  k = 6,
  Q_cut = 0.6,
  ratio = 0.2,
  max_iter = 100,
  run_neighbors = TRUE,
  force_recalc = FALSE,
  find_neighbors_params = list(),
  rare_threshold = 0.01,
  prefix = "RareQ",
  cluster_colname = paste0(prefix, "_cluster"),
  q_colname = paste0(prefix, "_Q"),
  size_colname = paste0(prefix, "_cluster_size"),
  rare_colname = paste0(prefix, "_is_rare"),
  tool_name = "RareQ",
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 2) {
    log_message(
      "{.arg k} must be a single number >= 2",
      message_type = "error"
    )
  }
  if (
    !is.numeric(k.param) ||
      length(k.param) != 1L ||
      is.na(k.param) ||
      k.param < 2
  ) {
    log_message(
      "{.arg k.param} must be a single number >= 2",
      message_type = "error"
    )
  }
  if (
    !is.numeric(Q_cut) ||
      length(Q_cut) != 1L ||
      is.na(Q_cut) ||
      Q_cut < 0 ||
      Q_cut > 1
  ) {
    log_message(
      "{.arg Q_cut} must be a single number between 0 and 1",
      message_type = "error"
    )
  }
  if (
    !is.numeric(ratio) ||
      length(ratio) != 1L ||
      is.na(ratio) ||
      ratio < 0 ||
      ratio > 1
  ) {
    log_message(
      "{.arg ratio} must be a single number between 0 and 1",
      message_type = "error"
    )
  }
  if (
    !is.numeric(max_iter) ||
      length(max_iter) != 1L ||
      is.na(max_iter) ||
      max_iter < 1
  ) {
    log_message(
      "{.arg max_iter} must be a single number >= 1",
      message_type = "error"
    )
  }
  if (!is.numeric(dims) || length(dims) == 0L || any(!is.finite(dims))) {
    log_message(
      "{.arg dims} must be a numeric vector of reduction dimensions",
      message_type = "error"
    )
  }
  if (!is.null(rare_threshold)) {
    if (
      !is.numeric(rare_threshold) ||
        length(rare_threshold) != 1L ||
        is.na(rare_threshold) ||
        rare_threshold < 0
    ) {
      log_message(
        "{.arg rare_threshold} must be a single number >= 0 or NULL",
        message_type = "error"
      )
    }
  }
  if (
    !is.list(find_neighbors_params) ||
      (length(find_neighbors_params) > 0L &&
        (is.null(names(find_neighbors_params)) ||
          any(!nzchar(names(find_neighbors_params)))))
  ) {
    log_message(
      "{.arg find_neighbors_params} must be a named list",
      message_type = "error"
    )
  }

  check_r("fabotao/RareQ", verbose = FALSE)

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!assay %in% Seurat::Assays(srt)) {
    log_message(
      "{.arg assay} {.val {assay}} is not present in {.cls Seurat}",
      message_type = "error"
    )
  }
  neighbor_slot <- paste0(assay, ".nn")
  old_assay <- SeuratObject::DefaultAssay(srt)
  SeuratObject::DefaultAssay(srt) <- assay

  k <- as.integer(k)
  k.param <- as.integer(k.param)
  max_iter <- as.integer(max_iter)

  need_neighbors <- isTRUE(force_recalc) ||
    !neighbor_slot %in% names(srt@neighbors)
  if (!need_neighbors) {
    nn_idx <- srt@neighbors[[neighbor_slot]]@nn.idx
    if (is.null(nn_idx) || ncol(nn_idx) < k) {
      if (isTRUE(run_neighbors)) {
        need_neighbors <- TRUE
        k.param <- max(k.param, k)
      } else {
        log_message(
          "Neighbor slot {.val {neighbor_slot}} contains fewer than {.arg k} neighbors",
          message_type = "error"
        )
      }
    }
  }

  dims_use <- NULL
  if (isTRUE(need_neighbors)) {
    if (!isTRUE(run_neighbors)) {
      log_message(
        "Missing neighbor slot {.val {neighbor_slot}}. Set {.arg run_neighbors = TRUE} or run {.fn Seurat::FindNeighbors} first",
        message_type = "error"
      )
    }
    reduction <- if (is.null(reduction)) {
      DefaultReduction(srt)
    } else {
      DefaultReduction(srt, pattern = reduction)
    }
    emb <- SeuratObject::Embeddings(srt, reduction = reduction)
    dims_use <- dims[dims >= 1 & dims <= ncol(emb)]
    if (length(dims_use) == 0L) {
      log_message(
        "No valid {.arg dims} are available in reduction {.val {reduction}}",
        message_type = "error"
      )
    }
    find_neighbors_args <- utils::modifyList(
      list(
        object = srt,
        reduction = reduction,
        dims = dims_use,
        assay = assay,
        k.param = k.param,
        return.neighbor = TRUE,
        compute.SNN = FALSE,
        prune.SNN = 0,
        verbose = verbose
      ),
      find_neighbors_params
    )
    log_message(
      "Build {.pkg Seurat} nearest neighbors for {.pkg RareQ} using reduction {.val {reduction}}",
      verbose = verbose
    )
    srt <- do.call(Seurat::FindNeighbors, find_neighbors_args)
    if (!neighbor_slot %in% names(srt@neighbors)) {
      log_message(
        "{.pkg Seurat} did not create the expected neighbor slot {.val {neighbor_slot}}",
        message_type = "error"
      )
    }
  }

  nn_idx <- srt@neighbors[[neighbor_slot]]@nn.idx
  if (is.null(nn_idx) || ncol(nn_idx) < k) {
    log_message(
      "Neighbor slot {.val {neighbor_slot}} must contain at least {.arg k} neighbors",
      message_type = "error"
    )
  }

  log_message(
    "Run {.pkg RareQ} with {.arg k = {k}}, {.arg Q_cut = {Q_cut}}, and {.arg ratio = {ratio}}",
    verbose = verbose
  )
  q_values <- get_namespace_fun("RareQ", "ComputeQ")(
    sc_object = srt,
    assay = assay,
    k = k
  )
  clusters <- get_namespace_fun("RareQ", "FindRare")(
    sc_object = srt,
    assay = assay,
    k = k,
    Q_cut = Q_cut,
    ratio = ratio,
    max_iter = max_iter
  )
  clusters <- as.character(clusters)
  names(clusters) <- colnames(srt)
  q_values <- as.numeric(q_values)
  names(q_values) <- colnames(srt)

  cluster_sizes <- table(clusters)
  cell_cluster_size <- as.integer(cluster_sizes[clusters])
  metadata <- data.frame(
    RareQ_cluster = clusters,
    RareQ_Q = q_values,
    RareQ_cluster_size = cell_cluster_size,
    row.names = colnames(srt),
    stringsAsFactors = FALSE
  )
  colnames(metadata) <- c(cluster_colname, q_colname, size_colname)
  rare_cutoff <- NULL
  if (!is.null(rare_threshold)) {
    rare_cutoff <- if (rare_threshold < 1) {
      ncol(srt) * rare_threshold
    } else {
      rare_threshold
    }
    metadata[[rare_colname]] <- cell_cluster_size <= rare_cutoff
  }
  srt <- Seurat::AddMetaData(srt, metadata = metadata)

  cluster_summary <- data.frame(
    cluster = names(cluster_sizes),
    n_cells = as.integer(cluster_sizes),
    fraction = as.integer(cluster_sizes) / ncol(srt),
    stringsAsFactors = FALSE
  )
  cluster_summary <- cluster_summary[
    order(cluster_summary$n_cells),
    ,
    drop = FALSE
  ]
  srt@tools[[tool_name]] <- list(
    clusters = clusters,
    q_values = q_values,
    cluster_summary = cluster_summary,
    neighbor_slot = neighbor_slot,
    parameters = list(
      assay = assay,
      reduction = reduction,
      dims = dims_use %||% dims,
      k.param = k.param,
      k = k,
      Q_cut = Q_cut,
      ratio = ratio,
      max_iter = max_iter,
      run_neighbors = run_neighbors,
      force_recalc = force_recalc,
      rare_threshold = rare_threshold,
      rare_cutoff = rare_cutoff,
      prefix = prefix,
      cluster_colname = cluster_colname,
      q_colname = q_colname,
      size_colname = size_colname,
      rare_colname = rare_colname
    )
  )

  SeuratObject::DefaultAssay(srt) <- old_assay
  log_message(
    "{.pkg RareQ} clusters stored in metadata column {.val {cluster_colname}}",
    verbose = verbose
  )
  srt
}
