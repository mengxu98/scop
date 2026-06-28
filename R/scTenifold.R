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
    count_matrix <- Matrix::Matrix(
      if (is.data.frame(count_matrix)) as.matrix(count_matrix) else count_matrix,
      sparse = TRUE
    )
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
      for (pkg in c("scTenifoldNet", "RSpectra", "RhpcBLASctl", "MASS")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
          log_message(
            "Package {.pkg {pkg}} is required for {.arg backend = 'cpp'}",
            message_type = "error"
          )
        }
      }

      log_message(
        "Construct scTenifoldNet network ensemble",
        verbose = verbose
      )
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
      wt <- sctenifold_tensor_cpp(
        wt,
        k = as.integer(td_K),
        max_iter = as.integer(td_maxIter),
        tol = td_maxError,
        n_decimal = as.integer(td_nDecimal),
        verbose = verbose
      )
      wt <- sctenifold_strict_direction(
        as.matrix(wt),
        lambda = nc_lambda
      )
      diag(wt) <- 0
      wt <- t(wt)
      ko <- wt
      ko[gKO, ] <- 0
      ma <- sctenifold_manifold_cpp(
        wt,
        ko,
        d = as.integer(ma_nDim),
        cores = as.integer(cores)
      )
      dr <- sctenifold_dregulation(ma, gKO)

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

sctenifold_tensor_cpp <- function(
  x_list,
  k,
  max_iter,
  tol,
  n_decimal,
  verbose
) {
  gene_names <- rownames(x_list[[1]])
  n_genes <- length(gene_names)
  n_net <- length(x_list)
  set.seed(1)
  init_u <- list(
    matrix(stats::rnorm(n_genes * k), nrow = n_genes, ncol = k),
    matrix(stats::rnorm(n_genes * k), nrow = n_genes, ncol = k),
    matrix(stats::rnorm(k), nrow = 1, ncol = k),
    matrix(stats::rnorm(n_net * k), nrow = n_net, ncol = k)
  )
  tensor <- sctenifold_tensor_decomposition(
    lapply(x_list, as.matrix),
    init_u = init_u,
    max_iter = as.integer(max_iter),
    tol = tol,
    verbose = verbose
  )
  tensor <- round(tensor, as.integer(n_decimal))
  rownames(tensor) <- colnames(tensor) <- gene_names
  methods::as(tensor, "dgCMatrix")
}

sctenifold_manifold_cpp <- function(x, y, d, cores) {
  shared_genes <- intersect(rownames(x), rownames(y))
  x <- as.matrix(x[shared_genes, shared_genes, drop = FALSE])
  y <- as.matrix(y[shared_genes, shared_genes, drop = FALSE])
  w <- sctenifold_manifold_matrix(x, y)
  RhpcBLASctl::omp_set_num_threads(as.integer(cores))
  RhpcBLASctl::blas_set_num_threads(as.integer(cores))
  eig <- suppressWarnings(RSpectra::eigs(w, as.integer(d) * 2L, "SR"))
  eig$values <- suppressWarnings(as.numeric(eig$values))
  eig$vectors <- suppressWarnings(apply(eig$vectors, 2, as.numeric))
  new_order <- order(eig$values)
  eig$values <- eig$values[new_order]
  eig$vectors <- eig$vectors[, new_order, drop = FALSE]
  eig$vectors <- eig$vectors[, eig$values > 1e-08, drop = FALSE]
  aligned <- eig$vectors[, seq_len(as.integer(d)), drop = FALSE]
  colnames(aligned) <- paste0("NLMA ", seq_len(as.integer(d)))
  rownames(aligned) <- c(paste0("X_", shared_genes), paste0("Y_", shared_genes))
  aligned
}

sctenifold_dregulation <- function(manifold_output, gko) {
  gene_list <- rownames(manifold_output)
  gene_list <- gene_list[grepl("^X_", gene_list)]
  gene_list <- sub("^X_", "", gene_list)
  n_genes <- length(gene_list)
  expected <- rownames(manifold_output)
  expected <- expected[grepl("^Y_", expected)]
  expected <- sub("^Y_", "", expected)
  if (n_genes != nrow(manifold_output) / 2 || !all(expected == gene_list)) {
    log_message(
      "Genes are not ordered as expected in manifold alignment output",
      message_type = "error"
    )
  }
  d_metric <- sctenifold_pair_distances(manifold_output)
  lambda_values <- seq(-2, 2, length.out = 1000)
  lambda_values <- lambda_values[lambda_values != 0]
  bc <- try(
    MASS::boxcox(
      d_metric[d_metric > 0] ~ 1,
      plot = FALSE,
      lambda = lambda_values
    ),
    silent = TRUE
  )
  n_d <- if (inherits(bc, "try-error")) {
    d_metric
  } else {
    lambda <- bc$x[which.max(bc$y)]
    if (lambda < 0) 1 / (d_metric^lambda) else d_metric^lambda
  }
  z <- scale(n_d)
  d_out <- data.frame(gene = gene_list, distance = d_metric, Z = z)
  d_out <- d_out[order(d_out$distance, decreasing = TRUE), ]
  fc <- (d_out$distance^2) / mean((d_out$distance[-seq_len(length(gko))]^2))
  p_values <- stats::pchisq(q = fc, df = 1, lower.tail = FALSE)
  d_out$FC <- fc
  d_out$p.value <- p_values
  d_out$p.adj <- stats::p.adjust(p_values, method = "fdr")
  as.data.frame.array(d_out)
}

#' @title Run scTenifoldNet network comparison
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param object A `Seurat` object or a raw count matrix with genes in rows and
#' cells in columns.
#' @param y A second raw count matrix. Required when `object` is a matrix and
#' ignored when `object` is a `Seurat` object.
#' @param group.by Metadata column used to split a `Seurat` object into the two
#' conditions being compared.
#' @param condition1,condition2 Condition labels from `group.by`. If omitted,
#' the first two group levels are used.
#' @param assay,layer Assay and layer used as the count matrix when `object` is
#' a `Seurat` object.
#' @param features Optional genes to retain before running the comparison.
#' @param qc Whether to apply scTenifoldNet-style quality control.
#' @param qc_min_library_size,qc_remove_outlier_cells,qc_min_pct,qc_max_mt_ratio
#' Quality-control parameters forwarded to `scTenifoldNet::scQC()`.
#' @param nc_nNet,nc_nCells,nc_nComp,nc_symmetric,nc_scaleScores,nc_q Network
#' construction parameters forwarded to `scTenifoldNet::makeNetworks()`.
#' @param td_K,td_nDecimal,td_maxIter,td_maxError Tensor decomposition
#' parameters forwarded to `scTenifoldNet::tensorDecomposition()`.
#' @param ma_nDim Manifold-alignment dimension forwarded to
#' `scTenifoldNet::manifoldAlignment()`.
#' @param cores Number of cores forwarded to `scTenifoldNet::scTenifoldNet()`.
#' @param store_networks Whether to keep tensor networks in the stored result
#' when `object` is a `Seurat` object.
#' @param store_manifold Whether to keep manifold-alignment coordinates in the
#' stored result when `object` is a `Seurat` object.
#' @param tool_name Name of the `object@tools` entry when `object` is a `Seurat`
#' object.
#'
#' @return A scTenifoldNet result list for matrix input, or a `Seurat` object
#' with results stored in `object@tools[[tool_name]]`.
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' counts <- GetAssayData5(pancreas_sub, assay = "RNA", layer = "counts")
#' detected <- names(sort(Matrix::rowSums(counts > 0), decreasing = TRUE))
#' features_use <- head(detected, 300)
#'
#' pancreas_sub <- RunscTenifoldNet(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   condition1 = "Ductal",
#'   condition2 = "Endocrine",
#'   features = features_use,
#'   qc = FALSE,
#'   nc_nNet = 3,
#'   nc_nCells = 200,
#'   td_maxIter = 200,
#'   ma_nDim = 2,
#'   store_networks = FALSE,
#'   store_manifold = TRUE
#' )
#'
#' dr <- pancreas_sub@tools$scTenifoldNet$diffRegulation
#' head(dr)
#'
#' scTenifoldNetPlot(pancreas_sub, plot_type = "effect")
RunscTenifoldNet <- function(
  object,
  y = NULL,
  group.by = NULL,
  condition1 = NULL,
  condition2 = NULL,
  assay = NULL,
  layer = "counts",
  features = NULL,
  qc = TRUE,
  qc_min_library_size = 1000,
  qc_remove_outlier_cells = TRUE,
  qc_min_pct = 0.05,
  qc_max_mt_ratio = 0.1,
  nc_nNet = 10,
  nc_nCells = 500,
  nc_nComp = 3,
  nc_symmetric = FALSE,
  nc_scaleScores = TRUE,
  nc_q = 0.05,
  td_K = 3,
  td_nDecimal = 1,
  td_maxIter = 1000,
  td_maxError = 1e-05,
  ma_nDim = 30,
  cores = 1,
  store_networks = TRUE,
  store_manifold = TRUE,
  tool_name = "scTenifoldNet",
  verbose = TRUE
) {
  numeric_params <- list(
    nc_q = c(value = nc_q, lower = 0, upper = 1),
    td_maxError = c(value = td_maxError, lower = 0, upper = 1),
    qc_min_pct = c(value = qc_min_pct, lower = 0, upper = 1),
    qc_max_mt_ratio = c(value = qc_max_mt_ratio, lower = 0, upper = 1)
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
    nc_nComp = c(value = nc_nComp, lower = 2),
    td_K = c(value = td_K, lower = 1),
    td_nDecimal = c(value = td_nDecimal, lower = 0),
    td_maxIter = c(value = td_maxIter, lower = 1),
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

  is_seurat <- inherits(object, "Seurat")
  assay <- assay %||% if (is_seurat) SeuratObject::DefaultAssay(object) else NULL
  input_summary <- list(type = if (is_seurat) "Seurat" else class(object)[[1]])

  if (is_seurat) {
    if (is.null(group.by) || !group.by %in% colnames(object[[]])) {
      log_message(
        "{.arg group.by} must identify a metadata column when {.arg object} is a {.cls Seurat} object",
        message_type = "error"
      )
    }
    groups <- object[[group.by]][, 1]
    condition_levels <- levels(groups) %||% unique(as.character(groups))
    condition1 <- condition1 %||% condition_levels[[1]]
    condition2 <- condition2 %||% condition_levels[[2]]
    if (
      is.null(condition1) ||
        is.null(condition2) ||
        identical(condition1, condition2)
    ) {
      log_message(
        "{.arg condition1} and {.arg condition2} must identify two different groups",
        message_type = "error"
      )
    }
    cells_x <- rownames(object[[]])[as.character(groups) == condition1]
    cells_y <- rownames(object[[]])[as.character(groups) == condition2]
    if (length(cells_x) == 0L || length(cells_y) == 0L) {
      log_message(
        "Both selected conditions must contain cells",
        message_type = "error"
      )
    }
    count_matrix <- GetAssayData5(object, assay = assay, layer = layer)
    x <- count_matrix[, cells_x, drop = FALSE]
    y <- count_matrix[, cells_y, drop = FALSE]
    input_summary$group.by <- group.by
    input_summary$condition1 <- condition1
    input_summary$condition2 <- condition2
    input_summary$cells <- stats::setNames(
      c(length(cells_x), length(cells_y)),
      c(condition1, condition2)
    )
  } else {
    x <- object
    if (is.null(y)) {
      log_message(
        "{.arg y} is required when {.arg object} is a matrix",
        message_type = "error"
      )
    }
  }

  if (is.null(rownames(x)) || is.null(rownames(y))) {
    log_message(
      "Both input matrices must contain gene names as row names",
      message_type = "error"
    )
  }
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("X_cell_", seq_len(ncol(x)))
  }
  if (is.null(colnames(y))) {
    colnames(y) <- paste0("Y_cell_", seq_len(ncol(y)))
  }

  if (!is.null(features)) {
    features <- unique(as.character(features))
    shared_features <- Reduce(intersect, list(features, rownames(x), rownames(y)))
    if (length(shared_features) == 0L) {
      log_message(
        "No requested {.arg features} are present in both inputs",
        message_type = "error"
      )
    }
    missing_features <- setdiff(features, shared_features)
    if (length(missing_features) > 0L) {
      log_message(
        "Ignoring {.val {length(missing_features)}} requested features absent from at least one input",
        message_type = "warning",
        verbose = verbose
      )
    }
    x <- x[shared_features, , drop = FALSE]
    y <- y[shared_features, , drop = FALSE]
  }

  if (!inherits(x, "Matrix")) {
    x <- Matrix::Matrix(
      if (is.data.frame(x)) as.matrix(x) else x,
      sparse = TRUE
    )
  }
  if (!inherits(y, "Matrix")) {
    y <- Matrix::Matrix(
      if (is.data.frame(y)) as.matrix(y) else y,
      sparse = TRUE
    )
  }
  x <- methods::as(x, "dgCMatrix")
  y <- methods::as(y, "dgCMatrix")

  if (nrow(x) <= nc_nComp || nrow(y) <= nc_nComp) {
    log_message(
      "{.arg nc_nComp} must be lower than the number of genes in both inputs",
      message_type = "error"
    )
  }

  log_message(
    "Run {.pkg scTenifoldNet} comparison using upstream implementation",
    verbose = verbose
  )

  check_r("cailab-tamu/scTenifoldNet", verbose = FALSE)
  result <- get_namespace_fun("scTenifoldNet", "scTenifoldNet")(
    X = x,
    Y = y,
    qc = qc,
    qc_minLibSize = as.integer(qc_min_library_size),
    qc_removeOutlierCells = isTRUE(qc_remove_outlier_cells),
    qc_minPCT = qc_min_pct,
    qc_maxMTratio = qc_max_mt_ratio,
    nc_nNet = as.integer(nc_nNet),
    nc_nCells = as.integer(nc_nCells),
    nc_nComp = as.integer(nc_nComp),
    nc_symmetric = isTRUE(nc_symmetric),
    nc_scaleScores = isTRUE(nc_scaleScores),
    nc_q = nc_q,
    td_K = as.integer(td_K),
    td_nDecimal = as.integer(td_nDecimal),
    td_maxIter = as.integer(td_maxIter),
    td_maxError = td_maxError,
    ma_nDim = as.integer(ma_nDim),
    nCores = as.integer(cores)
  )

  if (!is.list(result) || is.null(result$diffRegulation)) {
    log_message(
      "{.pkg scTenifoldNet} did not return a valid {.field diffRegulation} table",
      message_type = "error"
    )
  }
  result$qc_summary <- result$qc_summary %||% list(applied = qc)

  parameters <- list(
    assay = assay,
    layer = layer,
    features = features,
    qc = qc,
    qc_min_library_size = qc_min_library_size,
    qc_remove_outlier_cells = qc_remove_outlier_cells,
    qc_min_pct = qc_min_pct,
    qc_max_mt_ratio = qc_max_mt_ratio,
    nc_nNet = nc_nNet,
    nc_nCells = nc_nCells,
    nc_nComp = nc_nComp,
    nc_symmetric = nc_symmetric,
    nc_scaleScores = nc_scaleScores,
    nc_q = nc_q,
    td_K = td_K,
    td_nDecimal = td_nDecimal,
    td_maxIter = td_maxIter,
    td_maxError = td_maxError,
    ma_nDim = ma_nDim,
    cores = cores
  )

  if (!is_seurat) {
    result$parameters <- parameters
    result$input_summary <- input_summary
    return(result)
  }

  stored_result <- result
  if (!isTRUE(store_networks)) {
    stored_result$tensorNetworks <- NULL
  }
  if (!isTRUE(store_manifold)) {
    stored_result$manifoldAlignment <- NULL
  }
  object@tools[[tool_name]] <- list(
    diffRegulation = result$diffRegulation,
    result = stored_result,
    qc_summary = result$qc_summary %||% list(applied = qc),
    input_summary = input_summary,
    parameters = c(
      parameters,
      list(
        store_networks = store_networks,
        store_manifold = store_manifold
      )
    )
  )

  log_message(
    "{.pkg scTenifoldNet} results stored in {.code object@tools[[{tool_name}]]}",
    message_type = "success",
    verbose = verbose
  )
  object
}
