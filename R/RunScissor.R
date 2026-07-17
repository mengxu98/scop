#' @title Run Scissor phenotype-associated cell selection
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param srt A `Seurat` object containing single-cell expression data.
#' @param bulk_dataset A bulk expression matrix-like object or a
#' `SummarizedExperiment`. Rows are genes and columns are bulk samples.
#' @param phenotype Phenotype annotation for bulk samples. For
#' `family = "binomial"`, character or factor values are converted to 0/1
#' according to `positive`. If `NULL` and `bulk_dataset` is a
#' `SummarizedExperiment`, `condition.by` is used.
#' @param condition.by Column in `colData(bulk_dataset)` used as phenotype
#' when `phenotype = NULL`.
#' @param positive Positive phenotype level for binomial Scissor. If `NULL`,
#' the second sorted level is used.
#' @param bulk_assay Assay name used when `bulk_dataset` is a
#' `SummarizedExperiment`.
#' @param layer Assay layer used from `srt`.
#' @param features Optional genes used before intersecting bulk and
#' single-cell features.
#' @param family Regression family passed to Scissor.
#' @param backend Scissor backend. `"r"` calls the upstream package path and
#' `"cpp"` uses the optimized `scop` path. Legacy aliases `"original"` and
#' `"scop"` are accepted as `"r"` and `"cpp"`, respectively.
#' @param alpha Scissor alpha search values. If `NULL`, Scissor's default grid
#' is used.
#' @param cutoff Maximum selected-cell fraction used by Scissor's alpha search.
#' @param tag Optional phenotype labels for printed/stored summaries.
#' @param graph Existing Seurat graph used by the `"cpp"` backend. If `NULL`,
#' an assay SNN graph or `standard_scop()` graph is reused when present,
#' otherwise a temporary Scissor-like graph is built.
#' @param dims PCA dimensions used when a temporary graph is built.
#' @param nfeatures Number of variable features used when a temporary graph is
#' built.
#' @param seed Random seed used by Scissor's alpha loop.
#' @param prefix Prefix for metadata column names.
#' @param tool_name Name of the `srt@tools` entry.
#' @param store_inputs Whether to store Scissor regression inputs in
#' `srt@tools`. Default is `FALSE` to avoid large objects.
#'
#' @return A `Seurat` object with Scissor status and coefficient columns in
#' metadata and a Scissor result bundle in `srt@tools[[tool_name]]`.
#' @export
#'
#' @seealso [ScissorPlot]
#'
#' @references
#' Sun, D. et al. Identifying phenotype-associated subpopulations by
#' integrating bulk and single-cell sequencing data. \emph{Nature
#' Biotechnology} (2021). \doi{10.1038/s41587-021-01091-3}
#'
#' @examples
#' data(panc8_sub)
#' data(islet_bulk)
#' panc8_sub <- standard_scop(panc8_sub, verbose = FALSE)
#' panc8_sub <- RunScissor(
#'   panc8_sub,
#'   bulk_dataset = islet_bulk,
#'   condition.by = "condition",
#'   positive = "bfa",
#'   family = "binomial",
#'   features = head(
#'     intersect(
#'       rownames(panc8_sub),
#'       rownames(SummarizedExperiment::assay(islet_bulk, "counts"))
#'     ), 1000
#'   ),
#'   alpha = 0.2,
#'   cutoff = 0.5
#' )
#'
#' ScissorPlot(
#'   panc8_sub,
#'   xlab = "UMAP_1",
#'   ylab = "UMAP_2"
#' )
RunScissor <- function(
  srt,
  bulk_dataset,
  phenotype = NULL,
  condition.by = NULL,
  positive = NULL,
  assay = NULL,
  bulk_assay = "counts",
  layer = "counts",
  features = NULL,
  family = c("gaussian", "binomial", "cox"),
  backend = c("cpp", "r"),
  alpha = NULL,
  cutoff = 0.2,
  tag = NULL,
  graph = NULL,
  dims = 1:10,
  nfeatures = 2000,
  seed = 123,
  prefix = "Scissor",
  tool_name = "Scissor",
  store_inputs = FALSE,
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (missing(bulk_dataset) || is.null(bulk_dataset)) {
    log_message(
      "{.arg bulk_dataset} must be provided",
      message_type = "error"
    )
  }
  family <- match.arg(family)
  backend <- backend[[1]]
  backend <- match.arg(backend, choices = c("cpp", "r", "scop", "original"))
  if (identical(backend, "scop")) {
    backend <- "cpp"
  } else if (identical(backend, "original")) {
    backend <- "r"
  }
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  alpha <- alpha %||%
    c(
      0.005,
      0.01,
      0.05,
      0.1,
      0.2,
      0.3,
      0.4,
      0.5,
      0.6,
      0.7,
      0.8,
      0.9
    )

  check_r("preprocessCore", verbose = FALSE)
  if (identical(backend, "r")) {
    check_r("sunduanchen/Scissor", verbose = FALSE)
  }

  bulk_matrix <- scissor_bulk_matrix(
    bulk_dataset = bulk_dataset,
    bulk_assay = bulk_assay
  )
  phenotype_info <- scissor_phenotype(
    bulk_dataset = bulk_dataset,
    phenotype = phenotype,
    condition.by = condition.by,
    positive = positive,
    family = family,
    tag = tag
  )
  phenotype <- phenotype_info$phenotype
  tag <- phenotype_info$tag

  if (length(phenotype) != ncol(bulk_matrix)) {
    log_message(
      "{.arg phenotype} length must match bulk sample number",
      message_type = "error"
    )
  }
  names(phenotype) <- colnames(bulk_matrix)

  sc_matrix <- GetAssayData5(
    object = srt,
    assay = assay,
    layer = layer
  )
  if (is.null(rownames(sc_matrix)) || is.null(colnames(sc_matrix))) {
    log_message(
      "Single-cell expression matrix must contain gene and cell names",
      message_type = "error"
    )
  }
  if (!is.null(features)) {
    features <- unique(as.character(features))
    keep_features <- intersect(
      features,
      intersect(rownames(bulk_matrix), rownames(sc_matrix))
    )
    if (length(keep_features) == 0L) {
      log_message(
        "No requested {.arg features} are shared by bulk and single-cell data",
        message_type = "error"
      )
    }
    bulk_matrix <- bulk_matrix[keep_features, , drop = FALSE]
    sc_matrix <- sc_matrix[keep_features, , drop = FALSE]
  }

  inputs <- scissor_build_inputs(
    srt = srt,
    bulk_matrix = bulk_matrix,
    sc_matrix = sc_matrix,
    assay = assay,
    features = features,
    graph = graph,
    dims = dims,
    nfeatures = nfeatures,
    verbose = verbose
  )
  quality_check <- stats::quantile(inputs$X)
  if (quality_check[[3]] < 0.01) {
    log_message(
      "The median correlation between single-cell and bulk samples is relatively low",
      message_type = "warning",
      verbose = verbose
    )
  }

  if (identical(backend, "r")) {
    scissor_fun <- get_namespace_fun("Scissor", "Scissor")
    load_file <- tempfile("scissor_inputs_", fileext = ".RData")
    on.exit(unlink(load_file), add = TRUE)
    X <- inputs$X
    Y <- phenotype
    network <- as.matrix(inputs$network)
    Expression_bulk <- inputs$Expression_bulk
    Expression_cell <- inputs$Expression_cell
    save(
      X,
      Y,
      network,
      Expression_bulk,
      Expression_cell,
      file = load_file
    )
    call_original <- function() {
      scissor_fun(
        bulk_dataset = matrix(0, nrow = 1, ncol = 1),
        sc_dataset = matrix(0, nrow = 1, ncol = 1),
        phenotype = phenotype,
        tag = tag,
        alpha = alpha,
        cutoff = cutoff,
        family = family,
        Load_file = load_file
      )
    }
    if (isTRUE(verbose)) {
      result <- call_original()
    } else {
      utils::capture.output(result <- call_original())
    }
    input_summary <- inputs$summary
    input_summary$backend <- "r"
  } else {
    result <- NULL
    fit_cv <- NULL
    for (i in seq_along(alpha)) {
      set.seed(seed)
      fit0 <- scissor_fit_apml1(
        x = inputs$X,
        y = phenotype,
        family = family,
        penalty = "Net",
        alpha = alpha[[i]],
        Omega = inputs$network,
        nlambda = 100,
        nfolds = min(10L, nrow(inputs$X)),
        inzero = identical(family, "cox")
      )
      fit1 <- scissor_fit_apml1(
        x = inputs$X,
        y = phenotype,
        family = family,
        penalty = "Net",
        alpha = alpha[[i]],
        Omega = inputs$network,
        lambda = fit0$lambda.min
      )
      if (identical(family, "binomial")) {
        coefs <- as.numeric(fit1$Beta[seq_len(ncol(inputs$X)) + 1L])
      } else {
        coefs <- as.numeric(fit1$Beta)
      }
      names(coefs) <- colnames(inputs$X)
      cell_pos <- names(coefs)[coefs > 0]
      cell_neg <- names(coefs)[coefs < 0]
      percentage <- (length(cell_pos) + length(cell_neg)) / ncol(inputs$X)

      log_message(
        "Scissor alpha {.val {alpha[[i]]}} selected {.val {length(cell_pos)}} positive and {.val {length(cell_neg)}} negative cells ({.val {round(percentage * 100, 3)}}%)",
        verbose = verbose
      )
      result <- list(
        para = list(
          alpha = alpha[[i]],
          lambda = fit0$lambda.min,
          family = family
        ),
        Coefs = coefs,
        Scissor_pos = cell_pos,
        Scissor_neg = cell_neg
      )
      fit_cv <- fit0
      if (percentage < cutoff) {
        break
      }
    }
    result$fit <- fit_cv$fit
    input_summary <- inputs$summary
  }

  coefs <- result$Coefs
  names(coefs) <- names(coefs) %||% colnames(srt)
  coefs_full <- stats::setNames(rep(0, ncol(srt)), colnames(srt))
  coefs_full[names(coefs)] <- coefs
  status <- stats::setNames(rep("Background", ncol(srt)), colnames(srt))
  status[intersect(result$Scissor_pos, colnames(srt))] <- "Scissor+"
  status[intersect(result$Scissor_neg, colnames(srt))] <- "Scissor-"
  status <- factor(status, levels = c("Scissor+", "Scissor-", "Background"))

  meta <- data.frame(
    coef = coefs_full,
    status = status,
    selected = status != "Background",
    row.names = colnames(srt)
  )
  colnames(meta) <- paste0(prefix, c("_coef", "_status", "_selected"))
  srt <- Seurat::AddMetaData(srt, metadata = meta)

  bundle <- list(
    result = result,
    parameters = list(
      backend = backend,
      family = family,
      alpha = alpha,
      cutoff = cutoff,
      tag = tag,
      assay = assay,
      layer = layer,
      bulk_assay = bulk_assay,
      condition.by = condition.by,
      positive = phenotype_info$positive,
      prefix = prefix
    ),
    input_summary = input_summary,
    quality_check = quality_check
  )
  if (isTRUE(store_inputs)) {
    bundle$inputs <- inputs
  }
  srt@tools[[tool_name]] <- bundle
  srt <- Seurat::LogSeuratCommand(srt)

  log_message(
    "{.pkg Scissor} stored {.val {length(result$Scissor_pos)}} Scissor+ and {.val {length(result$Scissor_neg)}} Scissor- cells",
    message_type = "success",
    verbose = verbose
  )
  srt
}

scissor_bulk_matrix <- function(bulk_dataset, bulk_assay = "counts") {
  if (inherits(bulk_dataset, "SummarizedExperiment")) {
    assay_names <- SummarizedExperiment::assayNames(bulk_dataset)
    if (!bulk_assay %in% assay_names) {
      log_message(
        "{.arg bulk_assay} {.val {bulk_assay}} is not present in {.arg bulk_dataset}",
        message_type = "error"
      )
    }
    bulk_matrix <- SummarizedExperiment::assay(bulk_dataset, bulk_assay)
  } else {
    bulk_matrix <- bulk_dataset
  }
  if (inherits(bulk_matrix, "data.frame")) {
    bulk_matrix <- as.matrix(bulk_matrix)
  }
  if (!inherits(bulk_matrix, c("matrix", "Matrix"))) {
    log_message(
      "{.arg bulk_dataset} must be matrix-like or a {.cls SummarizedExperiment}",
      message_type = "error"
    )
  }
  if (is.null(rownames(bulk_matrix)) || is.null(colnames(bulk_matrix))) {
    log_message(
      "{.arg bulk_dataset} must contain gene and sample names",
      message_type = "error"
    )
  }
  as.matrix(bulk_matrix)
}

scissor_phenotype <- function(
  bulk_dataset,
  phenotype,
  condition.by,
  positive,
  family,
  tag
) {
  if (is.null(phenotype)) {
    if (!inherits(bulk_dataset, "SummarizedExperiment")) {
      log_message(
        "{.arg phenotype} is required unless {.arg bulk_dataset} is a {.cls SummarizedExperiment}",
        message_type = "error"
      )
    }
    cdata <- as.data.frame(SummarizedExperiment::colData(bulk_dataset))
    if (is.null(condition.by) || !condition.by %in% colnames(cdata)) {
      log_message(
        "{.arg condition.by} must be a column in {.fn colData}",
        message_type = "error"
      )
    }
    phenotype <- cdata[[condition.by]]
  }

  if (identical(family, "cox")) {
    phenotype <- as.matrix(phenotype)
    if (ncol(phenotype) != 2L) {
      log_message(
        "{.arg phenotype} must have two columns for {.val cox}",
        message_type = "error"
      )
    }
    return(list(phenotype = phenotype, tag = tag, positive = positive))
  }

  if (identical(family, "binomial")) {
    level_values <- sort(unique(as.character(phenotype)))
    if (length(level_values) != 2L) {
      log_message(
        "{.arg phenotype} must contain exactly two levels for {.val binomial}",
        message_type = "error"
      )
    }
    positive <- positive %||% level_values[[2]]
    if (!positive %in% level_values) {
      log_message(
        "{.arg positive} must be one phenotype level",
        message_type = "error"
      )
    }
    negative <- setdiff(level_values, positive)
    levels_use <- c(negative, positive)
    phenotype <- as.integer(factor(
      as.character(phenotype),
      levels = levels_use
    )) -
      1L
    tag <- tag %||% levels_use
    return(list(phenotype = phenotype, tag = tag, positive = positive))
  }

  phenotype <- as.numeric(phenotype)
  if (anyNA(phenotype)) {
    log_message(
      "{.arg phenotype} must be numeric for {.val gaussian}",
      message_type = "error"
    )
  }
  tag <- tag %||% names(table(phenotype))
  list(phenotype = phenotype, tag = tag, positive = positive)
}

scissor_build_inputs <- function(
  srt,
  bulk_matrix,
  sc_matrix,
  assay,
  features,
  graph,
  dims,
  nfeatures,
  verbose
) {
  common <- intersect(rownames(bulk_matrix), rownames(sc_matrix))
  if (!is.null(features)) {
    common <- intersect(common, features)
  }
  if (length(common) == 0L) {
    log_message(
      "There are no common genes between single-cell and bulk data",
      message_type = "error"
    )
  }

  network_info <- scissor_network(
    srt = srt,
    sc_matrix = sc_matrix,
    assay = assay,
    graph = graph,
    dims = dims,
    nfeatures = nfeatures,
    verbose = verbose
  )

  dataset0 <- cbind(
    bulk_matrix[common, , drop = FALSE],
    as.matrix(sc_matrix[common, , drop = FALSE])
  )
  dataset1 <- get_namespace_fun(
    "preprocessCore",
    "normalize.quantiles"
  )(dataset0)
  rownames(dataset1) <- rownames(dataset0)
  colnames(dataset1) <- colnames(dataset0)

  bulk_n <- ncol(bulk_matrix)
  expression_bulk <- dataset1[, seq_len(bulk_n), drop = FALSE]
  expression_cell <- dataset1[, (bulk_n + 1L):ncol(dataset1), drop = FALSE]
  X <- stats::cor(expression_bulk, expression_cell)

  list(
    X = X,
    network = network_info$network,
    Expression_bulk = expression_bulk,
    Expression_cell = expression_cell,
    summary = list(
      backend = "cpp",
      genes = length(common),
      features = common,
      bulk_samples = ncol(bulk_matrix),
      cells = ncol(sc_matrix),
      reused_graph = network_info$reused_graph,
      graph = network_info$graph
    )
  )
}

scissor_network <- function(
  srt,
  sc_matrix,
  assay,
  graph,
  dims,
  nfeatures,
  verbose
) {
  graph_name <- graph %||% scissor_default_graph(srt, assay)
  if (!is.null(graph_name)) {
    network <- srt@graphs[[graph_name]]
    network <- network[colnames(sc_matrix), colnames(sc_matrix), drop = FALSE]
    network <- scissor_binary_network(network)
    return(list(network = network, reused_graph = TRUE, graph = graph_name))
  }

  log_message(
    "Build a temporary Scissor-style SNN graph",
    verbose = verbose
  )
  tmp <- Seurat::CreateSeuratObject(counts = sc_matrix)
  tmp <- NormalizeData(tmp, verbose = FALSE)
  tmp <- FindVariableFeatures(
    tmp,
    selection.method = "vst",
    nfeatures = nfeatures,
    verbose = FALSE
  )
  tmp <- ScaleData(tmp, verbose = FALSE)
  tmp <- RunPCA(
    tmp,
    features = SeuratObject::VariableFeatures(tmp),
    verbose = FALSE
  )
  dims_use <- dims[dims <= ncol(SeuratObject::Embeddings(tmp, "pca"))]
  if (length(dims_use) == 0L) {
    log_message(
      "No PCA dimensions are available for Scissor graph construction",
      message_type = "error"
    )
  }
  tmp <- FindNeighbors(tmp, dims = dims_use, verbose = FALSE)
  network <- scissor_binary_network(
    tmp@graphs[[paste0(SeuratObject::DefaultAssay(tmp), "_snn")]]
  )
  list(network = network, reused_graph = FALSE, graph = NULL)
}

scissor_binary_network <- function(network) {
  network <- methods::as(network, "dgCMatrix")
  Matrix::diag(network) <- 0
  network <- Matrix::drop0(network)
  network@x[] <- 1
  network
}

scissor_default_graph <- function(srt, assay) {
  if (length(srt@graphs) == 0L) {
    return(NULL)
  }
  graph_names <- names(srt@graphs)
  candidates <- c(
    paste0(assay, "_snn"),
    paste0(assay, "_SNN"),
    paste0(assay, "_nn"),
    paste0(assay, "_NN"),
    paste0("Standard", "pca", "_SNN"),
    paste0("Standard", "pca", "_KNN")
  )
  hit <- candidates[candidates %in% names(srt@graphs)]
  if (length(hit) > 0L) {
    return(hit[[1]])
  }
  snn <- grep("_snn$", graph_names, value = TRUE, ignore.case = TRUE)
  if (length(snn) > 0L) {
    return(snn[[1]])
  }
  nn <- grep("_knn$|_nn$", graph_names, value = TRUE, ignore.case = TRUE)
  if (length(nn) > 0L) {
    return(nn[[1]])
  }
  NULL
}

scissor_fit_apml1 <- function(
  x,
  y,
  family = c("gaussian", "binomial", "cox"),
  penalty = c("Lasso", "Enet", "Net"),
  Omega = NULL,
  alpha = 1.0,
  lambda = NULL,
  nlambda = 50,
  rlambda = NULL,
  wbeta = rep(1, ncol(x)),
  sgn = rep(1, ncol(x)),
  nfolds = 1,
  foldid = NULL,
  inzero = TRUE,
  isd = FALSE,
  iysd = FALSE,
  keep.beta = FALSE,
  ifast = TRUE,
  thresh = 1e-7,
  maxit = 1e+5
) {
  family <- match.arg(family)
  penalty <- match.arg(penalty)

  fit <- switch(family,
    "gaussian" = {
      scissor_validate_native_net(
        family = family,
        penalty = penalty,
        Omega = Omega,
        alpha = alpha,
        rlambda = rlambda,
        wbeta = wbeta,
        sgn = sgn,
        iysd = iysd,
        keep.beta = keep.beta
      )
      foldid <- scissor_resolve_foldid(foldid, nfolds, nrow(x))
      scissor_gaussian_net_fit_cpp(
        x = x,
        y = y,
        omega = Omega,
        alpha = alpha,
        lambda = lambda,
        nlambda = nlambda,
        foldid = foldid,
        inzero = inzero,
        isd = isd,
        thresh = thresh,
        maxit = maxit
      )
    },
    "binomial" = {
      scissor_validate_native_net(
        family = family,
        penalty = penalty,
        Omega = Omega,
        alpha = alpha,
        rlambda = rlambda,
        wbeta = wbeta,
        sgn = sgn,
        keep.beta = keep.beta
      )
      foldid <- scissor_resolve_foldid(foldid, nfolds, nrow(x))
      scissor_binomial_net_fit_cpp(
        x = x,
        y = y,
        omega = Omega,
        alpha = alpha,
        lambda = lambda,
        nlambda = nlambda,
        foldid = foldid,
        inzero = inzero,
        isd = isd,
        thresh = thresh,
        maxit = maxit,
        threshP = 1e-5
      )
    },
    "cox" = scissor_fit_cox(
      x = x,
      y = y,
      Omega = Omega,
      alpha = alpha,
      lambda = lambda,
      nlambda = nlambda,
      rlambda = rlambda,
      nfolds = nfolds,
      foldid = foldid,
      inzero = inzero,
      wbeta = abs(wbeta),
      sgn = sgn,
      isd = isd,
      keep.beta = keep.beta,
      ifast = ifast,
      thresh = thresh,
      maxit = maxit
    )
  )
  fit$family <- family

  class(fit) <- "scissor_fit_apml1"
  fit
}


scissor_validate_native_net <- function(
  family,
  penalty,
  Omega,
  alpha,
  rlambda = NULL,
  wbeta = NULL,
  sgn = NULL,
  iysd = FALSE,
  keep.beta = FALSE
) {
  if (!identical(penalty, "Net")) {
    log_message(
      "The native {.val {family}} Scissor wrapper supports {.val Net} penalty only",
      message_type = "error"
    )
  }
  if (is.null(Omega) || alpha == 1.0) {
    log_message(
      "The native {.val {family}} Scissor wrapper requires {.arg Omega} and {.arg alpha < 1}",
      message_type = "error"
    )
  }
  if (
    !is.null(rlambda) ||
      isTRUE(iysd) ||
      isTRUE(keep.beta) ||
      any(abs(wbeta) != 1) ||
      any(sgn != 1)
  ) {
    log_message(
      "The native {.val {family}} Scissor wrapper supports the RunScissor Net path only",
      message_type = "error"
    )
  }
}


scissor_resolve_foldid <- function(foldid, nfolds, n) {
  if (!is.null(foldid)) {
    return(foldid)
  }
  if (nfolds <= 1L) {
    return(NULL)
  }
  sample(rep(seq_len(nfolds), length.out = n))
}

scissor_fit_cox <- function(
  x,
  y,
  Omega = NULL,
  alpha = 1.0,
  lambda = NULL,
  nlambda = 100,
  rlambda = NULL,
  nfolds = 1,
  foldid = NULL,
  inzero = TRUE,
  wbeta = rep(1, ncol(x)),
  sgn = rep(1, ncol(x)),
  isd = FALSE,
  keep.beta = FALSE,
  ifast = TRUE,
  thresh = 1e-6,
  maxit = 1e+5
) {
  N0 <- nrow(x)
  p <- ncol(x)
  ifast <- as.integer(ifast)

  aPen <- ifelse(all(wbeta > 0), TRUE, FALSE)

  if (is.null(lambda)) {
    ilambda <- 1
    if (is.null(rlambda)) {
      rlambda <- ifelse(N0 > p, 0.0001, 0.01)
    }
    lambda <- (rlambda)^(c(0:(nlambda - 1)) / (nlambda - 1))
  } else {
    ilambda <- 0
    nlambda <- length(lambda)
  }

  if (is.null(Omega)) {
    penalty <- ifelse(alpha == 1, "Lasso", "Enet")
    adaptive <- ifelse(any(wbeta != 1), TRUE, FALSE)
  } else {
    penalty <- ifelse(alpha == 1, "Lasso", "Net")
    adaptive <- c(
      ifelse(any(wbeta != 1), TRUE, FALSE),
      ifelse(any(sgn != 1), TRUE, FALSE)
    )

    Omega <- rbind(0, cbind(0, Omega))
    sgn1 <- c(1, sgn)

    if (any(Matrix::diag(Omega) != 0)) {
      Matrix::diag(Omega) <- 0
    }
    if (any(Omega < 0)) {
      Omega <- abs(Omega)
    }

    if (inherits(Omega, "dgCMatrix")) {
      W <- OmegaSC(Omega, sgn1)
      W$loc <- W$loc + 1
    } else {
      W <- OmegaC(Omega, sgn1)
      W$loc <- W$loc + 1
    }
    rm(Omega)
  }

  prep0 <- scissor_prepare_cox(x, y)
  out <- switch(penalty,
    "Net" = NetCoxC(
      prep0$x,
      prep0$tevent,
      alpha,
      lambda,
      nlambda,
      ilambda,
      wbeta,
      W$Omega,
      W$loc,
      W$nadj,
      prep0$N,
      prep0$nevent,
      prep0$nevent1,
      prep0$loc1,
      prep0$n,
      p,
      N0,
      thresh,
      maxit,
      ifast
    ),
    EnetCoxC(
      prep0$x,
      prep0$tevent,
      alpha,
      lambda,
      nlambda,
      ilambda,
      wbeta,
      prep0$N,
      prep0$nevent,
      prep0$nevent1,
      prep0$loc1,
      prep0$n,
      p,
      N0,
      thresh,
      maxit,
      ifast
    )
  )
  nlambdai <- out$nlambda
  if (nlambdai == 0) {
    return(NULL)
  }
  lambdai <- lambda[1:nlambdai]

  temi <- out$Beta[, 1:nlambdai]
  temi[is.na(temi)] <- 0.0
  out$Beta <- Matrix::Matrix(temi, sparse = TRUE)
  out$BetaSTD <- Matrix::Matrix(out$BetaSTD[, 1:nlambdai], sparse = TRUE)
  out$nzero <- as.numeric(Matrix::colSums(out$Beta != 0))
  out$flag <- out$flag[1:nlambdai]

  if (nfolds == 1 & is.null(foldid)) {
    fit <- data.frame(lambda = lambdai, nzero = out$nzero)
    if (!isd) {
      return(list(
        Beta = out$Beta,
        fit = fit,
        penalty = penalty,
        adaptive = adaptive,
        flag = out$flag
      ))
    } else {
      return(list(
        Beta = out$BetaSTD,
        fit = fit,
        penalty = penalty,
        adaptive = adaptive,
        flag = out$flag
      ))
    }
  } else {
    if (is.null(foldid)) {
      foldid <- sample(rep(seq(nfolds), length = N0))
    } else {
      nfolds <- max(foldid)
    }
    tb <- table(foldid)
    N0i <- numeric(nfolds)
    for (i in 1:nfolds) {
      N0i[i] <- sum(tb[-i])
    }

    prepk <- list()
    for (i in 1:nfolds) {
      temid <- which(foldid != i)
      prepk[[i]] <- scissor_prepare_cox(x[temid, ], y[temid, ])
    }
    weighti <- as.vector(tapply(y[, "status"], foldid, sum))

    outi <- list()
    cvPL <- matrix(NA, nrow = nfolds, ncol = nlambdai)
    for (i in 1:nfolds) {
      outi[[i]] <- switch(penalty,
        "Net" = cvNetCoxC(
          prepk[[i]]$x,
          prepk[[i]]$tevent,
          alpha,
          lambdai,
          nlambdai,
          wbeta,
          W$Omega,
          W$loc,
          W$nadj,
          prepk[[i]]$N,
          prepk[[i]]$nevent,
          prepk[[i]]$nevent1,
          prepk[[i]]$loc1,
          prepk[[i]]$n,
          p,
          N0i[i],
          thresh,
          maxit,
          0,
          prep0$x,
          prep0$N,
          prep0$nevent,
          prep0$nevent1,
          prep0$loc1,
          prep0$n
        ),
        cvEnetCoxC(
          prepk[[i]]$x,
          prepk[[i]]$tevent,
          alpha,
          lambdai,
          nlambdai,
          wbeta,
          prepk[[i]]$N,
          prepk[[i]]$nevent,
          prepk[[i]]$nevent1,
          prepk[[i]]$loc1,
          prepk[[i]]$n,
          p,
          N0i[i],
          thresh,
          maxit,
          0,
          prep0$x,
          prep0$N,
          prep0$nevent,
          prep0$nevent1,
          prep0$loc1,
          prep0$n
        )
      )
      cvPL[i, 1:outi[[i]]$nlambda] <- outi[[i]]$lf[1:outi[[i]]$nlambda] -
        outi[[i]]$ll[1:outi[[i]]$nlambda]
    }

    temi <- colSums(is.infinite(cvPL))
    if (any(temi > 0)) {
      nlambdai <- max(min(which(temi > 0)) - 1, 1)
    }

    cvraw <- cvPL / weighti
    nfoldi <- colSums(!is.na(cvraw)) # rm(cvPL) #
    cvm <- apply(cvraw, 2, stats::weighted.mean, w = weighti, na.rm = TRUE)
    cvse <- sqrt(
      apply(
        sweep(cvraw, 2, cvm, "-")^2,
        2,
        stats::weighted.mean,
        w = weighti,
        na.rm = TRUE
      ) /
        (nfoldi - 1)
    )

    cvraw <- cvraw[1:nlambdai]
    cvm <- cvm[1:nlambdai]
    cvse <- cvse[1:nlambdai]

    indexi <- which.max(cvm)
    indexij <- which(cvm >= (cvm[indexi] - cvse[indexi]))[1]
    temi <- rep("", nlambdai)
    temi[indexi] <- "*"

    temCV <- data.frame(
      lambda = lambdai[1:nlambdai],
      cvm = cvm,
      cvse = cvse,
      nzero = out$nzero[1:nlambdai],
      index = temi[1:nlambdai],
      stringsAsFactors = FALSE
    )

    if (!inzero) {
      rm(outi)
      temCV$cvm <- -temCV$cvm
      if (!keep.beta) {
        return(list(
          Beta = out$Beta[, indexi],
          fit = temCV,
          lambda.min = lambdai[indexi],
          penalty = penalty,
          adaptive = adaptive,
          flag = out$flag[1:nlambdai]
        ))
      } else {
        return(list(
          Beta = out$Beta[, 1:nlambdai],
          fit = temCV,
          lambda.min = lambdai[indexi],
          penalty = penalty,
          adaptive = adaptive,
          flag = out$flag[1:nlambdai]
        ))
      }
    }

    il0 <- indexi
    cvm <- list()
    cv.max <- rep(NA, nlambdai)
    repeat {
      numi <- out$nzero[il0]
      Betai <- sapply(outi, function(x) {
        x$Beta[, il0]
      })
      Betai[is.na(Betai)] <- 0
      BetaSTDi <- sapply(outi, function(x) {
        x$BetaSTD[, il0]
      })

      Betao <- colSums(Betai != 0)
      numi2 <- min(max(Betao), numi)

      if (numi2 > 0) {
        cvPL <- matrix(NA, nrow = nfolds, ncol = numi2)
        i <- 1
        for (i in 1:nfolds) {
          Betaj <- Betai[, i]
          BetaSTDj <- BetaSTDi[, i]
          numj <- min(Betao[i], numi)
          if (numj == 0) {
            cvPL[i, ] <- cvTrimCoxC(
              c(0.0, 0.0),
              numj,
              numi2,
              c(0, 0),
              prep0$x,
              prep0$N,
              prep0$nevent,
              prep0$nevent1,
              prep0$loc1,
              prep0$n,
              prepk[[i]]$x,
              prepk[[i]]$N,
              prepk[[i]]$nevent,
              prepk[[i]]$nevent1,
              prepk[[i]]$loc1,
              prepk[[i]]$n,
              0,
              1
            )
          } else {
            BetaSTDjj <- BetaSTDj
            BetaSTDjj[wbeta == 0] <- max(abs(BetaSTDj)) + 1
            temo <- rank(-abs(BetaSTDjj), ties.method = "min")

            temo <- data.frame(temo[which(temo <= numj)], which(temo <= numj))
            temo <- temo[order(temo[, 1]), ]
            cvPL[i, ] <- cvTrimCoxC(
              Betaj[temo[, 2]],
              numj,
              numi2,
              temo[, 2] - 1,
              prep0$x,
              prep0$N,
              prep0$nevent,
              prep0$nevent1,
              prep0$loc1,
              prep0$n,
              prepk[[i]]$x,
              prepk[[i]]$N,
              prepk[[i]]$nevent,
              prepk[[i]]$nevent1,
              prepk[[i]]$loc1,
              prepk[[i]]$n,
              0,
              1
            )
          }
        }
      } else {
        cvPL <- matrix(NA, nrow = nfolds, ncol = 1)
        for (i in 1:nfolds) {
          cvPL[i, ] <- cvTrimCoxC(
            c(0.0, 0.0),
            0,
            0,
            c(0, 0),
            prep0$x,
            prep0$N,
            prep0$nevent,
            prep0$nevent1,
            prep0$loc1,
            prep0$n,
            prepk[[i]]$x,
            prepk[[i]]$N,
            prepk[[i]]$nevent,
            prepk[[i]]$nevent1,
            prepk[[i]]$loc1,
            prepk[[i]]$n,
            0,
            1
          )
        }
      }

      cvraw <- cvPL / weighti
      nfoldi <- colSums(!is.na(cvraw))
      rm(cvPL)
      cvm[[il0]] <- apply(cvraw, 2, stats::weighted.mean, w = weighti, na.rm = TRUE)
      temi <- cvm[[il0]]
      if (aPen) {
        cv.max[il0] <- max(temi)
      } else {
        cv.max[il0] <- ifelse(
          length(temi) > sum(wbeta == 0),
          max(temi[-c(1:sum(wbeta == 0))]),
          temi[sum(wbeta == 0)]
        )
      }

      il1 <- c(il0 - 1, il0 + 1)
      for (j in 1:2) {
        if (il1[j] >= 1 & il1[j] <= nlambdai) {
          if (is.na(cv.max[il1[j]])) {
            numi <- out$nzero[il1[j]]
            Betai <- sapply(outi, function(x) {
              x$Beta[, il1[j]]
            })
            Betai[is.na(Betai)] <- 0
            BetaSTDi <- sapply(outi, function(x) {
              x$BetaSTD[, il1[j]]
            })

            Betao <- colSums(Betai != 0)
            numi2 <- min(max(Betao), numi)

            if (numi2 > 0) {
              cvPL <- matrix(NA, nrow = nfolds, ncol = numi2)
              for (i in 1:nfolds) {
                Betaj <- Betai[, i]
                BetaSTDj <- BetaSTDi[, i]
                numj <- min(Betao[i], numi)
                if (numj == 0) {
                  cvPL[i, ] <- cvTrimCoxC(
                    c(0.0, 0.0),
                    numj,
                    numi2,
                    c(0, 0),
                    prep0$x,
                    prep0$N,
                    prep0$nevent,
                    prep0$nevent1,
                    prep0$loc1,
                    prep0$n,
                    prepk[[i]]$x,
                    prepk[[i]]$N,
                    prepk[[i]]$nevent,
                    prepk[[i]]$nevent1,
                    prepk[[i]]$loc1,
                    prepk[[i]]$n,
                    0,
                    1
                  )
                } else {
                  BetaSTDjj <- BetaSTDj
                  BetaSTDjj[wbeta == 0] <- max(abs(BetaSTDj)) + 1
                  temo <- rank(-abs(BetaSTDjj), ties.method = "min")

                  temo <- data.frame(
                    temo[which(temo <= numj)],
                    which(temo <= numj)
                  )
                  temo <- temo[order(temo[, 1]), ]
                  cvPL[i, ] <- cvTrimCoxC(
                    Betaj[temo[, 2]],
                    numj,
                    numi2,
                    temo[, 2] - 1,
                    prep0$x,
                    prep0$N,
                    prep0$nevent,
                    prep0$nevent1,
                    prep0$loc1,
                    prep0$n,
                    prepk[[i]]$x,
                    prepk[[i]]$N,
                    prepk[[i]]$nevent,
                    prepk[[i]]$nevent1,
                    prepk[[i]]$loc1,
                    prepk[[i]]$n,
                    0,
                    1
                  )
                }
              }
            } else {
              cvPL <- matrix(NA, nrow = nfolds, ncol = 1)
              for (i in 1:nfolds) {
                cvPL[i, ] <- cvTrimCoxC(
                  c(0.0, 0.0),
                  0,
                  0,
                  c(0, 0),
                  prep0$x,
                  prep0$N,
                  prep0$nevent,
                  prep0$nevent1,
                  prep0$loc1,
                  prep0$n,
                  prepk[[i]]$x,
                  prepk[[i]]$N,
                  prepk[[i]]$nevent,
                  prepk[[i]]$nevent1,
                  prepk[[i]]$loc1,
                  prepk[[i]]$n,
                  0,
                  1
                )
              }
            }

            cvraw <- cvPL / weighti
            nfoldi <- colSums(!is.na(cvraw))
            rm(cvPL)
            cvm[[il1[j]]] <- apply(
              cvraw,
              2,
              stats::weighted.mean,
              w = weighti,
              na.rm = TRUE
            )

            temi <- cvm[[il1[j]]]
            if (aPen) {
              cv.max[il1[j]] <- max(temi)
            } else {
              cv.max[il1[j]] <- ifelse(
                length(temi) > sum(wbeta == 0),
                max(temi[-c(1:sum(wbeta == 0))]),
                temi[sum(wbeta == 0)]
              )
            }
          }
        } else {
          break
        }
      }
      if (il0 == which.max(cv.max)) {
        break
      } else {
        il0 <- which.max(cv.max)
      }
    }
    index0 <- which.max(cv.max)

    Beta0 <- out$Beta[, index0]
    BetaSTD0 <- out$BetaSTD[, index0]

    temi <- cvm[[index0]]
    if (aPen) {
      cuti <- which.max(temi)
    } else {
      cuti <- ifelse(
        length(temi) > sum(wbeta == 0),
        which.max(temi[-c(1:sum(wbeta == 0))]) + sum(wbeta == 0),
        sum(wbeta == 0)
      )
    }
    Beta0j <- out$BetaSTD[, index0]
    Beta0j[which(wbeta == 0)] <- max(abs(Beta0j)) + 1

    Beta0[abs(Beta0j) <= sort(abs(Beta0j), TRUE)[cuti + 1]] <- 0
    BetaSTD0[abs(Beta0j) <= sort(abs(Beta0j), TRUE)[cuti + 1]] <- 0

    temCV0 <- data.frame(
      lambda = lambdai[index0],
      cvm = cv.max[index0],
      nzero = cuti
    )

    temCV$cvm <- -temCV$cvm
    temCV0$cvm <- -temCV0$cvm
    if (!keep.beta) {
      if (!isd) {
        return(list(
          Beta = out$Beta[, indexi],
          Beta0 = Beta0,
          fit = temCV,
          fit0 = temCV0,
          lambda.min = lambdai[indexi],
          lambda.opt = lambdai[index0],
          penalty = penalty,
          adaptive = adaptive,
          flag = out$flag[1:nlambdai]
        ))
      } else {
        return(list(
          Beta = out$BetaSTD[, indexi],
          Beta0 = BetaSTD0,
          fit = temCV,
          fit0 = temCV0,
          lambda.min = lambdai[indexi],
          lambda.opt = lambdai[index0],
          penalty = penalty,
          adaptive = adaptive,
          flag = out$flag[1:nlambdai]
        ))
      }
    } else {
      if (!isd) {
        return(list(
          Beta = out$Beta[, 1:nlambdai],
          Beta0 = Beta0,
          fit = temCV,
          fit0 = temCV0,
          lambda.min = lambdai[indexi],
          lambda.opt = lambdai[index0],
          penalty = penalty,
          adaptive = adaptive,
          flag = out$flag[1:nlambdai]
        ))
      } else {
        return(list(
          Beta = out$BetaSTD[, 1:nlambdai],
          Beta0 = BetaSTD0,
          fit = temCV,
          fit0 = temCV0,
          lambda.min = lambdai[indexi],
          lambda.opt = lambdai[index0],
          penalty = penalty,
          adaptive = adaptive,
          flag = out$flag[1:nlambdai]
        ))
      }
    }
  }
}

scissor_prepare_cox <- function(x, y) {
  N0 <- nrow(x)
  oi <- order(y[, "status"], decreasing = TRUE)
  x <- x[oi, ]
  y <- y[oi, ]
  oi <- order(y[, "time"])
  x <- x[oi, ]
  y <- y[oi, ]

  i1 <- which(y[, "status"] == 1)
  mi1 <- min(i1) - 1
  if (mi1 != 0) {
    x <- x[-c(1:mi1), ]
    y <- y[-c(1:mi1), ]
  }
  ty <- y[, "time"]
  tevent <- y[, "status"]
  N <- nrow(x)
  n1 <- sum(y[, "status"])

  dty <- duplicated(ty) # ties

  if (any(dty)) {
    tevent0 <- tevent
    tevent0[which(dty)] <- 0

    ievent <- cumsum(tevent0)
    loc1 <- which(tevent0 == 1)
    nevent <- table(ievent)
    n <- length(unique(ievent))
    nevent1 <- tapply(tevent == 1, ievent, sum)
  } else {
    ievent <- cumsum(tevent)
    loc1 <- which(tevent == 1)
    nevent <- table(ievent)
    n <- length(unique(ievent))
    nevent1 <- rep(1, n)
  }

  return(list(
    x = x,
    N0 = N0,
    tevent = tevent,
    N = N,
    nevent = nevent,
    nevent1 = nevent1,
    loc1 = loc1,
    n = n
  ))
}
