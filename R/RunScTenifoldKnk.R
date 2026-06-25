#' @title Run scTenifoldKnk in-silico knockout analysis
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param layer Assay layer used as the count matrix.
#' @param gKO Gene symbol or symbols to knock out. All genes must be present
#' after optional feature and QC filtering.
#' @param features Optional genes to retain before running network construction.
#' If supplied, `gKO` is always retained when present in the input assay.
#' @param qc Whether to apply scTenifoldKnk-style quality control.
#' @param qc_mt_threshold Maximum mitochondrial read fraction per cell.
#' @param qc_min_library_size Minimum library size per cell.
#' @param qc_min_cells Minimum number of expressing cells required per gene.
#' @param nc_lambda,nc_nNet,nc_nCells,nc_nComp,nc_scaleScores,nc_symmetric,nc_q
#' Network construction parameters forwarded to `scTenifoldNet::makeNetworks()`.
#' @param td_K,td_maxIter,td_maxError,td_nDecimal Tensor decomposition
#' parameters forwarded to `scTenifoldNet::tensorDecomposition()`.
#' @param ma_nDim Manifold-alignment dimension forwarded to
#' `scTenifoldNet::manifoldAlignment()`.
#' @param cores Number of cores used by native network-construction workers and
#' forwarded to downstream linear algebra where applicable.
#' @param backend `r` calls `scTenifoldKnk::scTenifoldKnk()` directly and is the
#' default high-consistency path. `cpp` follows the upstream
#' `scTenifoldNet`/`scTenifoldKnk` network construction, tensor decomposition,
#' manifold alignment, and differential-regulation steps while keeping input
#' handling and result storage inside `scop`.
#' @param store_networks Whether to keep WT/KO tensor networks in
#' `srt@tools`.
#' @param store_manifold Whether to keep manifold-alignment coordinates in
#' `srt@tools`.
#' @param tool_name Name of the `srt@tools` entry.
#'
#' @return A `Seurat` object with scTenifoldKnk results stored in
#' `srt@tools[[tool_name]]`.
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' gene_use <- "Pdx1"
#' counts <- GetAssayData5(
#'   pancreas_sub,
#'   assay = "RNA",
#'   layer = "counts"
#' )
#' detected <- names(
#'   sort(Matrix::rowSums(counts > 0),
#'     decreasing = TRUE
#'   )
#' )
#' features_use <- unique(c(gene_use, head(detected, 300)))
#'
#' pancreas_sub <- RunscTenifoldKnk(
#'   pancreas_sub,
#'   gKO = gene_use,
#'   features = features_use,
#'   qc = FALSE,
#'   nc_nNet = 3,
#'   nc_nCells = 200,
#'   td_maxIter = 200,
#'   store_networks = FALSE,
#'   store_manifold = TRUE
#' )
#'
#' dr <- pancreas_sub@tools$scTenifoldKnk$diffRegulation
#' head(dr)
#'
#' scTenifoldKnkPlot(pancreas_sub, plot_type = "effect")
RunscTenifoldKnk <- function(
  srt,
  gKO,
  assay = NULL,
  layer = "counts",
  features = NULL,
  qc = TRUE,
  qc_mt_threshold = 0.1,
  qc_min_library_size = 1000,
  qc_min_cells = 25,
  nc_lambda = 0,
  nc_nNet = 10,
  nc_nCells = 500,
  nc_nComp = 3,
  nc_scaleScores = TRUE,
  nc_symmetric = FALSE,
  nc_q = 0.9,
  td_K = 3,
  td_maxIter = 1000,
  td_maxError = 1e-05,
  td_nDecimal = 3,
  ma_nDim = 2,
  cores = 1,
  backend = c("r", "cpp"),
  store_networks = TRUE,
  store_manifold = TRUE,
  tool_name = "scTenifoldKnk",
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (missing(gKO) || length(gKO) == 0L || any(is.na(gKO))) {
    log_message(
      "{.arg gKO} must contain at least one gene symbol",
      message_type = "error"
    )
  }
  backend <- match.arg(backend)
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  gKO <- unique(as.character(gKO))

  numeric_params <- list(
    nc_lambda = c(value = nc_lambda, lower = 0, upper = 1),
    nc_q = c(value = nc_q, lower = 0, upper = 1),
    td_maxError = c(value = td_maxError, lower = 0, upper = 1)
  )
  for (param_name in names(numeric_params)) {
    param <- numeric_params[[param_name]]
    if (
      length(param[["value"]]) != 1L ||
        !is.numeric(param[["value"]]) ||
        is.na(param[["value"]]) ||
        param[["value"]] < param[["lower"]] ||
        param[["value"]] > param[["upper"]]
    ) {
      log_message(
        "{.arg {param_name}} must be a single number between {.val {param[['lower']]}} and {.val {param[['upper']]}}",
        message_type = "error"
      )
    }
  }
  integer_params <- list(
    nc_nNet = c(value = nc_nNet, lower = 1),
    nc_nCells = c(value = nc_nCells, lower = 1),
    nc_nComp = c(value = nc_nComp, lower = 3),
    td_K = c(value = td_K, lower = 1),
    td_maxIter = c(value = td_maxIter, lower = 1),
    td_nDecimal = c(value = td_nDecimal, lower = 0),
    ma_nDim = c(value = ma_nDim, lower = 1),
    cores = c(value = cores, lower = 1)
  )
  for (param_name in names(integer_params)) {
    param <- integer_params[[param_name]]
    if (
      length(param[["value"]]) != 1L ||
        !is.numeric(param[["value"]]) ||
        is.na(param[["value"]]) ||
        param[["value"]] < param[["lower"]] ||
        param[["value"]] != as.integer(param[["value"]])
    ) {
      log_message(
        "{.arg {param_name}} must be a single integer greater than or equal to {.val {param[['lower']]}}",
        message_type = "error"
      )
    }
  }

  count_matrix <- GetAssayData5(srt, assay = assay, layer = layer)
  if (is.null(rownames(count_matrix)) || is.null(colnames(count_matrix))) {
    log_message(
      "{.arg layer} must contain feature and cell names",
      message_type = "error"
    )
  }

  missing_ko <- setdiff(gKO, rownames(count_matrix))
  if (length(missing_ko) > 0L) {
    log_message(
      "{.arg gKO} genes are absent from the selected assay/layer: {.val {missing_ko}}",
      message_type = "error"
    )
  }

  if (!is.null(features)) {
    features <- unique(c(as.character(features), gKO))
    missing_features <- setdiff(features, rownames(count_matrix))
    features <- intersect(features, rownames(count_matrix))
    if (length(features) == 0L) {
      log_message(
        "No requested {.arg features} are present in the selected assay/layer",
        message_type = "error"
      )
    }
    if (length(missing_features) > 0L) {
      log_message(
        "Ignoring {.val {length(missing_features)}} requested features absent from the selected assay/layer",
        message_type = "warning",
        verbose = verbose
      )
    }
    count_matrix <- count_matrix[features, , drop = FALSE]
  }

  if (!inherits(count_matrix, "Matrix")) {
    count_matrix <- Matrix::Matrix(as.matrix(count_matrix), sparse = TRUE)
  }
  count_matrix <- methods::as(count_matrix, "dgCMatrix")
  qc_summary <- list(
    before = c(genes = nrow(count_matrix), cells = ncol(count_matrix)),
    after = c(genes = nrow(count_matrix), cells = ncol(count_matrix)),
    applied = FALSE
  )

  if (isTRUE(qc)) {
    before <- c(genes = nrow(count_matrix), cells = ncol(count_matrix))
    library_size <- Matrix::colSums(count_matrix)
    keep_library <- library_size >= qc_min_library_size
    count_matrix <- count_matrix[, keep_library, drop = FALSE]
    if (ncol(count_matrix) < 3L) {
      log_message(
        "Fewer than three cells remain after library-size filtering",
        message_type = "error"
      )
    }

    library_size <- Matrix::colSums(count_matrix)
    n_genes <- Matrix::colSums(count_matrix != 0)
    genes_lm <- stats::lm(n_genes ~ library_size)
    genes_pred <- as.data.frame(stats::predict(
      genes_lm,
      data.frame(library_size = library_size),
      interval = "prediction"
    ))

    mt_genes <- grep("^MT-", toupper(rownames(count_matrix)))
    if (length(mt_genes) > 0L) {
      mt_counts <- Matrix::colSums(count_matrix[mt_genes, , drop = FALSE])
      mt_proportion <- mt_counts / library_size
      mt_lm <- stats::lm(mt_counts ~ library_size)
      mt_pred <- as.data.frame(stats::predict(
        mt_lm,
        data.frame(library_size = library_size),
        interval = "prediction"
      ))
      selected_cells <- (mt_counts > mt_pred$lwr &
        mt_counts < mt_pred$upr &
        n_genes > genes_pred$lwr &
        n_genes < genes_pred$upr &
        mt_proportion <= qc_mt_threshold &
        library_size < 2 * mean(library_size))
    } else {
      selected_cells <- (n_genes > genes_pred$lwr &
        n_genes < genes_pred$upr &
        library_size < 2 * mean(library_size))
    }

    count_matrix <- count_matrix[, selected_cells, drop = FALSE]
    gene_keep <- Matrix::rowSums(count_matrix != 0) >= qc_min_cells
    count_matrix <- count_matrix[gene_keep, , drop = FALSE]

    missing_ko <- setdiff(gKO, rownames(count_matrix))
    if (length(missing_ko) > 0L) {
      log_message(
        "{.arg gKO} genes are absent after scTenifoldKnk QC filtering: {.val {missing_ko}}",
        message_type = "error"
      )
    }
    if (nrow(count_matrix) == 0L || ncol(count_matrix) == 0L) {
      log_message(
        "No genes or cells remain after scTenifoldKnk QC filtering",
        message_type = "error"
      )
    }

    count_matrix <- methods::as(count_matrix, "dgCMatrix")
    qc_summary <- list(
      before = before,
      after = c(genes = nrow(count_matrix), cells = ncol(count_matrix)),
      applied = TRUE,
      mt_threshold = qc_mt_threshold,
      min_library_size = qc_min_library_size,
      min_cells = qc_min_cells
    )
  }

  if (!all(gKO %in% rownames(count_matrix))) {
    log_message(
      "{.arg gKO} must be present in the retained count matrix",
      message_type = "error"
    )
  }
  if (nrow(count_matrix) <= nc_nComp) {
    log_message(
      "{.arg nc_nComp} must be lower than the number of retained genes",
      message_type = "error"
    )
  }
  if (length(gKO) >= nrow(count_matrix)) {
    log_message(
      "{.arg gKO} cannot cover all retained genes",
      message_type = "error"
    )
  }

  log_message(
    "Run {.pkg scTenifoldKnk} knockout for {.val {gKO}} using {.val {backend}} backend",
    verbose = verbose
  )

  result <- switch(
    backend,
    cpp = {
      check_r(c("scTenifoldNet", "MASS"), verbose = FALSE)

      log_message(
        "Construct scTenifoldNet network ensemble",
        verbose = verbose
      )
      geneList <- rownames(count_matrix)
      nGenes <- length(geneList)
      nCol <- ncol(count_matrix)
      if (nGenes <= 0L) {
        log_message(
          "Gene names are required for scTenifoldNet network construction",
          message_type = "error"
        )
      }
      if (!(nc_nComp > 1L && nc_nComp < nGenes)) {
        log_message(
          "{.arg nc_nComp} should be greater or equal than 2 and lower than the total number of genes",
          message_type = "error"
        )
      }
      wt_networks <- lapply(
        seq_len(nc_nNet),
        function(i) {
          Z <- sample(x = seq_len(nCol), size = nc_nCells, replace = TRUE)
          Z <- as.matrix(count_matrix[, Z])
          Z <- Z[rowSums(Z) > 0, , drop = FALSE]
          Z <- get_namespace_fun("scTenifoldNet", "pcNet")(
            X = Z,
            nComp = as.integer(nc_nComp),
            scaleScores = isTRUE(nc_scaleScores),
            symmetric = isTRUE(nc_symmetric),
            q = nc_q,
            verbose = FALSE,
            nCores = as.integer(cores)
          )
          if (nrow(Z) == nGenes && identical(rownames(Z), geneList)) {
            return(Z)
          }

          O <- matrix(data = 0, nrow = nGenes, ncol = nGenes)
          rownames(O) <- colnames(O) <- geneList
          O[rownames(Z), colnames(Z)] <- as.matrix(Z)
          methods::as(O, "dgCMatrix")
        }
      )

      RhpcBLASctl::omp_set_num_threads(as.integer(cores))
      RhpcBLASctl::blas_set_num_threads(as.integer(cores))

      log_message(
        "Denoise network ensemble with tensor decomposition",
        verbose = verbose
      )
      gene_list <- rownames(wt_networks[[1]])
      same_names <- !is.null(gene_list) &&
        !is.null(colnames(wt_networks[[1]])) &&
        identical(gene_list, colnames(wt_networks[[1]])) &&
        all(vapply(
          wt_networks,
          function(x) {
            identical(rownames(x), gene_list) &&
              identical(colnames(x), gene_list)
          },
          logical(1)
        ))

      log_message(
        if (!same_names) {
          "Network names are not identical across tensors; falling back to {.pkg scTenifoldNet} tensor decomposition"
        } else {
          "Use {.pkg scTenifoldNet} tensor decomposition"
        },
        message_type = if (!same_names) "warning" else NULL,
        verbose = verbose
      )
      wt <- get_namespace_fun("scTenifoldNet", "tensorDecomposition")(
        xList = wt_networks,
        nDecimal = as.integer(td_nDecimal),
        K = as.integer(td_K),
        maxError = td_maxError,
        maxIter = as.integer(td_maxIter)
      )

      wt <- as.matrix(wt$X)
      wt <- getFromNamespace("strictDirection", "scTenifoldKnk")(
        wt,
        lambda = nc_lambda
      )
      wt <- as.matrix(wt)
      diag(wt) <- 0
      wt <- t(wt)
      ko <- wt
      ko[gKO, ] <- 0
      ma <- get_namespace_fun("scTenifoldNet", "manifoldAlignment")(
        wt,
        ko,
        d = as.integer(ma_nDim),
        nCores = as.integer(cores)
      )
      dr <- getFromNamespace("dRegulation", "scTenifoldKnk")(ma, gKO)

      list(
        tensorNetworks = if (isTRUE(store_networks)) {
          list(
            WT = Matrix::Matrix(wt),
            KO = Matrix::Matrix(ko)
          )
        } else {
          NULL
        },
        manifoldAlignment = if (isTRUE(store_manifold)) ma else NULL,
        diffRegulation = dr
      )
    },
    r = {
      check_r(c("scTenifoldKnk", "scTenifoldNet"), verbose = FALSE)

      wt <- get_namespace_fun("scTenifoldNet", "makeNetworks")(
        X = count_matrix,
        q = nc_q,
        nNet = as.integer(nc_nNet),
        nCells = as.integer(nc_nCells),
        scaleScores = isTRUE(nc_scaleScores),
        symmetric = isTRUE(nc_symmetric),
        nComp = as.integer(nc_nComp),
        nCores = as.integer(cores)
      )
      wt <- get_namespace_fun("scTenifoldNet", "tensorDecomposition")(
        xList = wt,
        K = as.integer(td_K),
        maxError = td_maxError,
        maxIter = as.integer(td_maxIter),
        nDecimal = as.integer(td_nDecimal)
      )
      wt <- wt$X
      wt <- getFromNamespace("strictDirection", "scTenifoldKnk")(
        wt,
        lambda = nc_lambda
      )
      wt <- as.matrix(wt)
      diag(wt) <- 0
      wt <- t(wt)

      ko <- wt
      ko[gKO, ] <- 0
      ma <- get_namespace_fun("scTenifoldNet", "manifoldAlignment")(
        wt,
        ko,
        d = as.integer(ma_nDim),
        nCores = as.integer(cores)
      )
      dr <- getFromNamespace("dRegulation", "scTenifoldKnk")(ma, gKO)

      list(
        tensorNetworks = if (isTRUE(store_networks)) {
          list(
            WT = Matrix::Matrix(wt),
            KO = Matrix::Matrix(ko)
          )
        } else {
          NULL
        },
        manifoldAlignment = if (isTRUE(store_manifold)) ma else NULL,
        diffRegulation = dr
      )
    }
  )
  if (!is.list(result) || is.null(result$diffRegulation)) {
    log_message(
      "{.pkg scTenifoldKnk} backend {.val {backend}} did not return a valid {.field diffRegulation} table",
      message_type = "error"
    )
  }

  if (!isTRUE(store_networks)) {
    result$tensorNetworks <- NULL
  }
  if (!isTRUE(store_manifold)) {
    result$manifoldAlignment <- NULL
  }

  srt@tools[[tool_name]] <- list(
    diffRegulation = result$diffRegulation,
    result = result,
    qc_summary = qc_summary,
    parameters = list(
      assay = assay,
      layer = layer,
      features = features,
      gKO = gKO,
      qc = qc,
      qc_mt_threshold = qc_mt_threshold,
      qc_min_library_size = qc_min_library_size,
      qc_min_cells = qc_min_cells,
      nc_lambda = nc_lambda,
      nc_nNet = nc_nNet,
      nc_nCells = nc_nCells,
      nc_nComp = nc_nComp,
      nc_scaleScores = nc_scaleScores,
      nc_symmetric = nc_symmetric,
      nc_q = nc_q,
      td_K = td_K,
      td_maxIter = td_maxIter,
      td_maxError = td_maxError,
      td_nDecimal = td_nDecimal,
      ma_nDim = ma_nDim,
      cores = cores,
      backend = backend,
      store_networks = store_networks,
      store_manifold = store_manifold
    )
  )

  log_message(
    "{.pkg scTenifoldKnk} results stored in {.code srt@tools[[{tool_name}]]}",
    message_type = "success",
    verbose = verbose
  )
  srt
}
