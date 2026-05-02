#' @title Run scTenifoldKnk in-silico knockout analysis
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
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
#' @param backend `optimized` uses the local scop implementation with native
#' equivalent covariance-based network construction, direct sparse network
#' assembly, controlled per-gene eigensolver parallelism, and helpers for
#' tensor decomposition, manifold matrix construction, directionality, and
#' distance calculation.
#' `upstream` calls `scTenifoldKnk::scTenifoldKnk()` directly for comparison.
#' @param store_networks Whether to keep WT/KO tensor networks in
#' `srt@tools`.
#' @param store_manifold Whether to keep manifold-alignment coordinates in
#' `srt@tools`.
#' @param tool_name Name of the `srt@tools` entry.
#'
#' @return A `Seurat` object with scTenifoldKnk results stored in
#' `srt@tools[[tool_name]]`.
#' @export
RunScTenifoldKnk <- function(
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
  backend = c("optimized", "upstream"),
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

  sctenifold_validate_parameters(
    nc_lambda = nc_lambda,
    nc_nNet = nc_nNet,
    nc_nCells = nc_nCells,
    nc_nComp = nc_nComp,
    nc_q = nc_q,
    td_K = td_K,
    td_maxIter = td_maxIter,
    td_maxError = td_maxError,
    td_nDecimal = td_nDecimal,
    ma_nDim = ma_nDim,
    cores = cores
  )

  count_matrix <- sctenifold_get_count_matrix(
    srt = srt,
    assay = assay,
    layer = layer,
    features = features,
    gKO = gKO,
    verbose = verbose
  )
  qc_summary <- list(
    before = c(genes = nrow(count_matrix), cells = ncol(count_matrix)),
    after = c(genes = nrow(count_matrix), cells = ncol(count_matrix)),
    applied = FALSE
  )

  if (isTRUE(qc)) {
    qc_res <- sctenifold_qc_filter(
      count_matrix = count_matrix,
      gKO = gKO,
      mt_threshold = qc_mt_threshold,
      min_library_size = qc_min_library_size,
      min_cells = qc_min_cells
    )
    count_matrix <- qc_res$count_matrix
    qc_summary <- qc_res$summary
  }

  sctenifold_validate_matrix_for_run(
    count_matrix = count_matrix,
    gKO = gKO,
    nc_nCells = nc_nCells,
    nc_nComp = nc_nComp
  )

  log_message(
    "Run {.pkg scTenifoldKnk} knockout for {.val {paste(gKO, collapse = ', ')}} using {.val {backend}} backend",
    verbose = verbose
  )

  result <- switch(
    backend,
    optimized = sctenifold_run_optimized(
      count_matrix = count_matrix,
      gKO = gKO,
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
      verbose = verbose
    ),
    upstream = sctenifold_run_upstream(
      count_matrix = count_matrix,
      gKO = gKO,
      qc = FALSE,
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
      cores = cores
    )
  )

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

sctenifold_validate_parameters <- function(
  nc_lambda,
  nc_nNet,
  nc_nCells,
  nc_nComp,
  nc_q,
  td_K,
  td_maxIter,
  td_maxError,
  td_nDecimal,
  ma_nDim,
  cores
) {
  sctenifold_validate_scalar(nc_lambda, "nc_lambda", lower = 0, upper = 1)
  sctenifold_validate_scalar(nc_q, "nc_q", lower = 0, upper = 1)
  sctenifold_validate_scalar(td_maxError, "td_maxError", lower = 0, upper = 1)
  sctenifold_validate_integer(nc_nNet, "nc_nNet", lower = 1)
  sctenifold_validate_integer(nc_nCells, "nc_nCells", lower = 1)
  sctenifold_validate_integer(nc_nComp, "nc_nComp", lower = 3)
  sctenifold_validate_integer(td_K, "td_K", lower = 1)
  sctenifold_validate_integer(td_maxIter, "td_maxIter", lower = 1)
  sctenifold_validate_integer(td_nDecimal, "td_nDecimal", lower = 0)
  sctenifold_validate_integer(ma_nDim, "ma_nDim", lower = 1)
  sctenifold_validate_integer(cores, "cores", lower = 1)
}

sctenifold_validate_scalar <- function(x, name, lower = -Inf, upper = Inf) {
  if (length(x) != 1L || !is.numeric(x) || is.na(x) || x < lower || x > upper) {
    log_message(
      "{.arg {name}} must be a single number between {.val {lower}} and {.val {upper}}",
      message_type = "error"
    )
  }
}

sctenifold_validate_integer <- function(x, name, lower = -Inf, upper = Inf) {
  if (length(x) != 1L || !is.numeric(x) || is.na(x) || x < lower || x > upper || x != as.integer(x)) {
    log_message(
      "{.arg {name}} must be a single integer between {.val {lower}} and {.val {upper}}",
      message_type = "error"
    )
  }
}

sctenifold_get_count_matrix <- function(srt, assay, layer, features, gKO, verbose = TRUE) {
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
  methods::as(count_matrix, "dgCMatrix")
}

sctenifold_qc_filter <- function(
  count_matrix,
  gKO,
  mt_threshold = 0.1,
  min_library_size = 1000,
  min_cells = 25
) {
  before <- c(genes = nrow(count_matrix), cells = ncol(count_matrix))
  library_size <- Matrix::colSums(count_matrix)
  keep_library <- library_size >= min_library_size
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
    selected_cells <- (
      mt_counts > mt_pred$lwr &
        mt_counts < mt_pred$upr &
        n_genes > genes_pred$lwr &
        n_genes < genes_pred$upr &
        mt_proportion <= mt_threshold &
        library_size < 2 * mean(library_size)
    )
  } else {
    selected_cells <- (
      n_genes > genes_pred$lwr &
        n_genes < genes_pred$upr &
        library_size < 2 * mean(library_size)
    )
  }

  count_matrix <- count_matrix[, selected_cells, drop = FALSE]
  gene_keep <- Matrix::rowSums(count_matrix != 0) >= min_cells
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

  list(
    count_matrix = methods::as(count_matrix, "dgCMatrix"),
    summary = list(
      before = before,
      after = c(genes = nrow(count_matrix), cells = ncol(count_matrix)),
      applied = TRUE,
      mt_threshold = mt_threshold,
      min_library_size = min_library_size,
      min_cells = min_cells
    )
  )
}

sctenifold_validate_matrix_for_run <- function(count_matrix, gKO, nc_nCells, nc_nComp) {
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
}

sctenifold_run_optimized <- function(
  count_matrix,
  gKO,
  nc_lambda,
  nc_nNet,
  nc_nCells,
  nc_nComp,
  nc_scaleScores,
  nc_symmetric,
  nc_q,
  td_K,
  td_maxIter,
  td_maxError,
  td_nDecimal,
  ma_nDim,
  cores,
  verbose = TRUE
) {
  check_r("scTenifoldNet", verbose = FALSE)
  if (!requireNamespace("MASS", quietly = TRUE)) {
    log_message(
      "{.pkg MASS} is required for scTenifoldKnk differential regulation statistics",
      message_type = "error"
    )
  }
  sctenifold_check_cpp()

  log_message(
    "Construct scTenifoldNet network ensemble",
    verbose = verbose
  )
  wt <- sctenifold_make_networks(
    X = count_matrix,
    q = nc_q,
    nNet = as.integer(nc_nNet),
    nCells = as.integer(nc_nCells),
    scaleScores = isTRUE(nc_scaleScores),
    symmetric = isTRUE(nc_symmetric),
    nComp = as.integer(nc_nComp),
    nCores = as.integer(cores)
  )

  RhpcBLASctl::omp_set_num_threads(as.integer(cores))
  RhpcBLASctl::blas_set_num_threads(as.integer(cores))

  log_message(
    "Denoise network ensemble with tensor decomposition",
    verbose = verbose
  )
  wt <- sctenifold_tensor_decomposition(
    xList = wt,
    K = as.integer(td_K),
    maxError = td_maxError,
    maxIter = as.integer(td_maxIter),
    nDecimal = as.integer(td_nDecimal)
  )

  wt <- as.matrix(wt$X)
  wt <- sctenifold_strict_direction_cpp(
    x = wt,
    lambda = nc_lambda
  )
  diag(wt) <- 0
  wt <- t(wt)

  ko <- wt
  ko[gKO, ] <- 0

  log_message(
    "Align WT and KO network manifolds",
    verbose = verbose
  )
  ma <- sctenifold_manifold_alignment(
    X = wt,
    Y = ko,
    d = as.integer(ma_nDim),
    nCores = as.integer(cores)
  )

  dr <- sctenifold_diff_regulation(
    manifold_output = ma,
    gKO = gKO
  )

  list(
    tensorNetworks = list(
      WT = Matrix::Matrix(wt),
      KO = Matrix::Matrix(ko)
    ),
    manifoldAlignment = ma,
    diffRegulation = dr
  )
}

sctenifold_run_upstream <- function(
  count_matrix,
  gKO,
  qc,
  qc_mt_threshold,
  qc_min_library_size,
  qc_min_cells,
  nc_lambda,
  nc_nNet,
  nc_nCells,
  nc_nComp,
  nc_scaleScores,
  nc_symmetric,
  nc_q,
  td_K,
  td_maxIter,
  td_maxError,
  td_nDecimal,
  ma_nDim,
  cores
) {
  check_r("scTenifoldKnk", verbose = FALSE)
  scTenifoldKnk::scTenifoldKnk(
    countMatrix = count_matrix,
    qc = isTRUE(qc),
    gKO = gKO,
    qc_mtThreshold = qc_mt_threshold,
    qc_minLSize = qc_min_library_size,
    qc_minCells = as.integer(qc_min_cells),
    nc_lambda = nc_lambda,
    nc_nNet = as.integer(nc_nNet),
    nc_nCells = as.integer(nc_nCells),
    nc_nComp = as.integer(nc_nComp),
    nc_scaleScores = isTRUE(nc_scaleScores),
    nc_symmetric = isTRUE(nc_symmetric),
    nc_q = nc_q,
    td_K = as.integer(td_K),
    td_maxIter = as.integer(td_maxIter),
    td_maxError = td_maxError,
    td_nDecimal = as.integer(td_nDecimal),
    ma_nDim = as.integer(ma_nDim),
    nCores = as.integer(cores)
  )
}

sctenifold_make_networks <- function(
  X,
  nNet = 10,
  nCells = 500,
  nComp = 3,
  scaleScores = TRUE,
  symmetric = FALSE,
  q = 0.95,
  nCores = parallel::detectCores()
) {
  geneList <- rownames(X)
  nGenes <- length(geneList)
  nCol <- ncol(X)
  if (nGenes <= 0L) {
    log_message(
      "Gene names are required for scTenifoldNet network construction",
      message_type = "error"
    )
  }
  if (!(nComp > 1L && nComp < nGenes)) {
    log_message(
      "{.arg nComp} should be greater or equal than 2 and lower than the total number of genes",
      message_type = "error"
    )
  }

  lapply(
    seq_len(nNet),
    function(W) {
      Z <- sample(x = seq_len(nCol), size = nCells, replace = TRUE)
      Z <- as.matrix(X[, Z])
      Z <- Z[apply(Z, 1, sum) > 0, , drop = FALSE]
      if (nrow(Z) >= ncol(Z)) {
        Z <- sctenifold_pcnet_covariance(
          X = Z,
          nComp = nComp,
          scaleScores = scaleScores,
          symmetric = symmetric,
          q = q,
          nCores = nCores
        )
      } else {
        Z <- get_namespace_fun("scTenifoldNet", "pcNet")(
          X = Z,
          nComp = nComp,
          scaleScores = scaleScores,
          symmetric = symmetric,
          q = q,
          verbose = FALSE,
          nCores = nCores
        )
      }
      if (nrow(Z) == nGenes && identical(rownames(Z), geneList)) {
        return(Z)
      }

      O <- matrix(data = 0, nrow = nGenes, ncol = nGenes)
      rownames(O) <- colnames(O) <- geneList
      O[rownames(Z), colnames(Z)] <- as.matrix(Z)
      methods::as(O, "dgCMatrix")
    }
  )
}

sctenifold_pcnet_covariance <- function(
  X,
  nComp = 3,
  scaleScores = TRUE,
  symmetric = FALSE,
  q = 0,
  nCores = parallel::detectCores()
) {
  if (!all(Matrix::rowSums(X) > 0)) {
    log_message(
      "Quality control has not been applied over the matrix",
      message_type = "error"
    )
  }
  xClass <- class(X)[[1]]
  validClass <- xClass %in% c("matrix", "dgCMatrix")
  if (!validClass) {
    log_message(
      "Input should be a matrix with cells as columns and genes as rows",
      message_type = "error"
    )
  }
  if (nComp < 2L || nComp >= nrow(X)) {
    log_message(
      "{.arg nComp} should be greater or equal than 2 and lower than the total number of genes",
      message_type = "error"
    )
  }

  gNames <- rownames(X)
  RhpcBLASctl::omp_set_num_threads(1L)
  RhpcBLASctl::blas_set_num_threads(1L)
  A <- sctenifold_pcnet_covariance_sparse_cpp(
    x = as.matrix(X),
    n_comp = as.integer(nComp),
    scale_scores = isTRUE(scaleScores),
    symmetric = isTRUE(symmetric),
    q = q,
    n_threads = as.integer(nCores)
  )
  colnames(A) <- rownames(A) <- gNames
  A
}

sctenifold_check_cpp <- function() {
  ok <- exists("sctenifold_strict_direction_cpp", mode = "function") &&
    exists("sctenifold_pair_distances_cpp", mode = "function") &&
    exists("sctenifold_tensor_decomposition_cpp", mode = "function") &&
    exists("sctenifold_manifold_matrix_cpp", mode = "function") &&
    exists("sctenifold_pcnet_covariance_raw_cpp", mode = "function") &&
    exists("sctenifold_pcnet_covariance_sparse_cpp", mode = "function") &&
    isTRUE(is.loaded("_scop_sctenifold_strict_direction_cpp")) &&
    isTRUE(is.loaded("_scop_sctenifold_pair_distances_cpp")) &&
    isTRUE(is.loaded("_scop_sctenifold_tensor_decomposition_cpp")) &&
    isTRUE(is.loaded("_scop_sctenifold_manifold_matrix_cpp")) &&
    isTRUE(is.loaded("_scop_sctenifold_pcnet_covariance_raw_cpp")) &&
    isTRUE(is.loaded("_scop_sctenifold_pcnet_covariance_sparse_cpp"))
  if (!ok) {
    log_message(
      "{.arg backend = 'optimized'} requires the compiled {.pkg scop} shared library. Reinstall the package to build native code.",
      message_type = "error"
    )
  }
}

sctenifold_tensor_decomposition <- function(
  xList,
  nDecimal = 1,
  K = 5,
  maxError = 1e-05,
  maxIter = 1000
) {
  sctenifold_check_cpp()
  if (!is.list(xList) || length(xList) == 0L) {
    log_message(
      "{.arg xList} must contain at least one network",
      message_type = "error"
    )
  }

  gene_list <- rownames(xList[[1]])
  same_names <- !is.null(gene_list) &&
    !is.null(colnames(xList[[1]])) &&
    identical(gene_list, colnames(xList[[1]])) &&
    all(vapply(
      xList,
      function(x) {
        identical(rownames(x), gene_list) && identical(colnames(x), gene_list)
      },
      logical(1)
    ))

  if (!same_names) {
    log_message(
      "Network names are not identical across tensors; falling back to {.pkg scTenifoldNet} tensor decomposition",
      message_type = "warning"
    )
    return(scTenifoldNet::tensorDecomposition(
      xList = xList,
      nDecimal = nDecimal,
      K = K,
      maxError = maxError,
      maxIter = maxIter
    ))
  }

  n_genes <- length(gene_list)
  n_net <- length(xList)
  set.seed(1)
  init_u <- list(
    matrix(stats::rnorm(n_genes * K), nrow = n_genes, ncol = K),
    matrix(stats::rnorm(n_genes * K), nrow = n_genes, ncol = K),
    matrix(stats::rnorm(K), nrow = 1, ncol = K),
    matrix(stats::rnorm(n_net * K), nrow = n_net, ncol = K)
  )
  x_dense <- lapply(
    xList,
    function(x) {
      x <- as.matrix(x)
      storage.mode(x) <- "double"
      x
    }
  )

  out <- sctenifold_tensor_decomposition_cpp(
    x_list = x_dense,
    init_u = init_u,
    max_iter = as.integer(maxIter),
    tol = maxError
  )
  out <- round(out, digits = nDecimal)
  out <- Matrix::Matrix(out, sparse = TRUE)
  rownames(out) <- colnames(out) <- gene_list
  list(X = out)
}

sctenifold_manifold_alignment <- function(X, Y, d = 30, nCores = parallel::detectCores()) {
  sctenifold_check_cpp()
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]
  check_r("RSpectra", verbose = FALSE)
  W <- sctenifold_manifold_matrix_cpp(
    x = as.matrix(X),
    y = as.matrix(Y)
  )
  RhpcBLASctl::omp_set_num_threads(nCores)
  RhpcBLASctl::blas_set_num_threads(nCores)
  E <- suppressWarnings(RSpectra::eigs(W, d * 2, "SR"))
  E$values <- suppressWarnings(as.numeric(E$values))
  E$vectors <- suppressWarnings(apply(E$vectors, 2, as.numeric))
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- E$vectors[, newOrder]
  E$vectors <- E$vectors[, E$values > 1e-08]
  alignedNet <- E$vectors[, seq_len(d)]
  colnames(alignedNet) <- paste0("NLMA ", seq_len(d))
  rownames(alignedNet) <- c(paste0("X_", sharedGenes), paste0("Y_", sharedGenes))
  alignedNet
}

sctenifold_diff_regulation <- function(manifold_output, gKO) {
  sctenifold_check_cpp()
  aligned <- as.matrix(manifold_output)
  if (nrow(aligned) %% 2L != 0L || is.null(rownames(aligned))) {
    log_message(
      "Unexpected scTenifoldNet manifold-alignment output",
      message_type = "error"
    )
  }

  n_genes <- nrow(aligned) / 2L
  x_names <- rownames(aligned)[seq_len(n_genes)]
  y_names <- rownames(aligned)[seq.int(n_genes + 1L, nrow(aligned))]
  gene_list <- sub("^X_", "", x_names)
  expected_gene_list <- sub("^Y_", "", y_names)

  if (!all(grepl("^X_", x_names)) ||
      !all(grepl("^Y_", y_names)) ||
      !identical(gene_list, expected_gene_list)) {
    log_message(
      "Genes are not ordered as expected in scTenifoldNet manifold output",
      message_type = "error"
    )
  }

  distances <- sctenifold_pair_distances_cpp(aligned)
  names(distances) <- gene_list

  positive <- distances[is.finite(distances) & distances > 0]
  if (length(positive) > 1L && length(unique(positive)) > 1L) {
    lambda_values <- seq(-2, 2, length.out = 1000)
    lambda_values <- lambda_values[lambda_values != 0]
    bc <- try(
      MASS::boxcox(positive ~ 1, plot = FALSE, lambda = lambda_values),
      silent = TRUE
    )
    if (inherits(bc, "try-error")) {
      normalized_distances <- distances
    } else {
      bc_lambda <- bc$x[which.max(bc$y)]
      if (bc_lambda < 0) {
        normalized_distances <- 1 / (distances^bc_lambda)
      } else {
        normalized_distances <- distances^bc_lambda
      }
    }
  } else {
    normalized_distances <- distances
  }

  z_score <- as.numeric(scale(normalized_distances))
  d_out <- data.frame(
    gene = gene_list,
    distance = as.numeric(distances),
    Z = z_score,
    stringsAsFactors = FALSE
  )
  d_out <- d_out[order(d_out$distance, decreasing = TRUE), , drop = FALSE]
  background_idx <- setdiff(seq_len(nrow(d_out)), seq_len(length(gKO)))
  if (length(background_idx) == 0L) {
    log_message(
      "Not enough non-KO genes to compute scTenifoldKnk differential regulation statistics",
      message_type = "error"
    )
  }
  distance_mean <- mean(d_out$distance[background_idx]^2)
  if (!is.finite(distance_mean) || distance_mean <= 0) {
    log_message(
      "Unable to compute finite scTenifoldKnk background distance",
      message_type = "error"
    )
  }
  fc <- (d_out$distance^2) / distance_mean
  d_out$FC <- fc
  d_out$p.value <- stats::pchisq(q = fc, df = 1, lower.tail = FALSE)
  d_out$p.adj <- stats::p.adjust(d_out$p.value, method = "fdr")
  rownames(d_out) <- NULL
  d_out
}
