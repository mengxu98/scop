grn_matrix_from_object <- function(
  object,
  genes_in = c("rows", "columns")
) {
  genes_in <- match.arg(genes_in)
  if (
    !inherits(object, "Matrix") && !is.matrix(object) && !is.data.frame(object)
  ) {
    log_message(
      "{.arg object} must be a matrix-like expression object or a {.cls Seurat} object",
      message_type = "error"
    )
  }
  mat <- if (inherits(object, "Matrix")) object else as.matrix(object)
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    log_message(
      "{.arg object} must have row and column names for GRN inference",
      message_type = "error"
    )
  }
  if (identical(genes_in, "rows")) {
    mat <- if (inherits(mat, "Matrix")) Matrix::t(mat) else t(mat)
  }
  if (inherits(mat, "sparseMatrix") && !inherits(mat, "dgCMatrix")) {
    mat <- methods::as(mat, "dgCMatrix")
  }
  if (nrow(mat) < 2L || ncol(mat) < 2L) {
    log_message(
      "GRN inference requires at least two samples and two genes",
      message_type = "error"
    )
  }
  mat
}

normalize_grn_inputs <- function(
  grn_matrix,
  regulators = NULL,
  targets = NULL
) {
  genes <- colnames(grn_matrix)
  regulators <- unique(as.character(regulators %||% genes))
  targets <- unique(as.character(targets %||% genes))
  regulators <- intersect(regulators, genes)
  targets <- intersect(genes, targets)
  if (length(regulators) == 0L) {
    log_message(
      "No candidate regulators are present in the expression matrix",
      message_type = "error"
    )
  }
  if (length(targets) == 0L) {
    log_message(
      "No target genes are present in the expression matrix",
      message_type = "error"
    )
  }
  list(regulators = regulators, targets = targets)
}

write_grn_adjacency <- function(
  adjacency,
  output_file = NULL,
  force = FALSE
) {
  if (is.null(output_file)) {
    return(adjacency)
  }
  if (file.exists(output_file) && isFALSE(force)) {
    return(utils::read.delim(
      output_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    adjacency,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  adjacency
}

filter_grn_adjacency_file <- function(
  input_file,
  output_file,
  targets,
  verbose = TRUE
) {
  adjacency <- utils::read.delim(
    input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!"target" %in% colnames(adjacency)) {
    log_message(
      "GRN adjacency file must contain a {.field target} column",
      message_type = "error"
    )
  }
  adjacency <- adjacency[
    adjacency[["target"]] %in% targets,
    ,
    drop = FALSE
  ]
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    adjacency,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  log_message(
    "Filtered GRN adjacency to {.val {nrow(adjacency)}} edges",
    verbose = verbose
  )
  invisible(output_file)
}

grn_python_packages <- function(module = c("grnboost2", "regdiffusion")) {
  module <- match.arg(module)
  if (identical(module, "regdiffusion")) {
    return(c("regdiffusion"))
  }
  c(
    "pyscenic==0.12.1",
    "arboreto==0.1.6",
    "ctxcore==0.2.0",
    "numpy==1.23.5",
    "dask==2024.2.1",
    "distributed==2024.2.1",
    "pyarrow"
  )
}

cap_grn_edges_per_target <- function(adjacency, max_edges_per_target = Inf) {
  max_edges_per_target <- suppressWarnings(as.numeric(max_edges_per_target))
  adjacency <- adjacency[
    is.finite(adjacency[["importance"]]) & adjacency[["importance"]] > 0,
    ,
    drop = FALSE
  ]
  if (nrow(adjacency) == 0L) {
    return(adjacency)
  }
  adjacency <- adjacency[
    order(adjacency[["target"]], -adjacency[["importance"]], adjacency[["TF"]]),
    ,
    drop = FALSE
  ]
  if (
    length(max_edges_per_target) == 1L &&
      is.finite(max_edges_per_target) &&
      max_edges_per_target > 0
  ) {
    adjacency <- do.call(
      rbind,
      lapply(
        split(adjacency, adjacency[["target"]]),
        utils::head,
        as.integer(max_edges_per_target)
      )
    )
  }
  rownames(adjacency) <- NULL
  adjacency[
    order(-adjacency[["importance"]], adjacency[["TF"]], adjacency[["target"]]),
    ,
    drop = FALSE
  ]
}

fill_grnboost_edges <- function(
  adjacency,
  grn_matrix,
  regulators,
  targets,
  max_edges_per_target = 50,
  boost_weight = 0.5,
  correlation_weight = 0.8,
  covariance_weight = 0.4,
  correlation_method = c("pearson", "spearman")
) {
  correlation_method <- match.arg(correlation_method)
  boost_weight <- suppressWarnings(as.numeric(boost_weight))
  correlation_weight <- suppressWarnings(as.numeric(correlation_weight))
  covariance_weight <- suppressWarnings(as.numeric(covariance_weight))
  if (length(boost_weight) != 1L || is.na(boost_weight) || boost_weight < 0) {
    boost_weight <- 0
  }
  if (
    length(correlation_weight) != 1L ||
      is.na(correlation_weight) ||
      correlation_weight < 0
  ) {
    correlation_weight <- 0
  }
  if (
    length(covariance_weight) != 1L ||
      is.na(covariance_weight) ||
      covariance_weight < 0
  ) {
    covariance_weight <- 0
  }
  if (boost_weight <= 0 && correlation_weight <= 0 && covariance_weight <= 0) {
    return(adjacency)
  }
  expr <- log1p(as.matrix(grn_matrix))
  reg_mat <- expr[, regulators, drop = FALSE]
  target_mat <- expr[, targets, drop = FALSE]
  cor_mat <- abs(suppressWarnings(stats::cor(
    reg_mat,
    target_mat,
    method = correlation_method,
    use = "pairwise.complete.obs"
  )))
  cor_mat[!is.finite(cor_mat)] <- 0
  cov_mat <- abs(stats::cov(reg_mat, target_mat, use = "pairwise.complete.obs"))
  cov_mat[!is.finite(cov_mat)] <- 0
  filled <- data.frame(
    TF = rep(rownames(cor_mat), times = ncol(cor_mat)),
    target = rep(colnames(cor_mat), each = nrow(cor_mat)),
    importance = as.numeric(cor_mat),
    covariance = as.numeric(cov_mat),
    stringsAsFactors = FALSE
  )
  filled <- filled[
    filled[["TF"]] != filled[["target"]] & filled[["importance"]] > 0,
    ,
    drop = FALSE
  ]
  if (nrow(filled) == 0L) {
    return(adjacency)
  }
  filled[["edge"]] <- paste(filled[["TF"]], filled[["target"]], sep = "\r")
  adjacency[["edge"]] <- paste(
    adjacency[["TF"]],
    adjacency[["target"]],
    sep = "\r"
  )
  filled[["boost"]] <- 0
  edge_match <- match(adjacency[["edge"]], filled[["edge"]])
  keep <- !is.na(edge_match)
  filled[["boost"]][edge_match[keep]] <- adjacency[["importance"]][keep]
  boost_score <- filled[["boost"]]
  cor_score <- filled[["importance"]]
  if (max(boost_score) > 0) {
    boost_score <- boost_score / max(boost_score)
  }
  if (max(cor_score) > 0) {
    cor_score <- cor_score / max(cor_score)
  }
  cov_score <- filled[["covariance"]]
  if (max(cov_score) > 0) {
    cov_score <- cov_score / max(cov_score)
  }
  filled[["importance"]] <- boost_weight *
    boost_score +
    correlation_weight * cor_score +
    covariance_weight * cov_score
  filled <- filled[, c("TF", "target", "importance"), drop = FALSE]
  cap_grn_edges_per_target(
    filled,
    max_edges_per_target = max_edges_per_target
  )
}

grnboost <- function(
  grn_matrix,
  regulators,
  targets = NULL,
  max_edges_per_target = Inf,
  n_rounds = 5000,
  learning_rate = 0.01,
  max_depth = 3,
  max_features = 0.1,
  subsample = 0.9,
  early_stop_window_length = 25,
  exclude_self = TRUE,
  importance_norm_power = 0,
  correlation_fill = FALSE,
  boost_weight = 0.5,
  correlation_weight = 0.8,
  covariance_weight = 0.4,
  correlation_method = c("pearson", "spearman"),
  cores = 1,
  seed = 1234,
  output_file = NULL,
  force = FALSE,
  verbose = TRUE
) {
  if (!is.null(output_file) && file.exists(output_file) && isFALSE(force)) {
    log_message("Reusing existing GRNBoost2 output", verbose = verbose)
    return(utils::read.delim(
      output_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  correlation_method <- match.arg(correlation_method)
  cores <- suppressWarnings(as.integer(cores))
  if (length(cores) != 1L || is.na(cores) || cores < 1L) {
    cores <- 1L
  }
  inputs <- normalize_grn_inputs(
    grn_matrix,
    regulators = regulators,
    targets = targets
  )
  genes_keep <- union(inputs[["regulators"]], inputs[["targets"]])
  grn_matrix <- grn_matrix[, genes_keep, drop = FALSE]
  gene_names <- colnames(grn_matrix)
  max_edges_per_target_cpp <- suppressWarnings(as.numeric(max_edges_per_target))
  max_edges_per_target_cpp <- if (
    length(max_edges_per_target_cpp) != 1L ||
      !is.finite(max_edges_per_target_cpp) ||
      max_edges_per_target_cpp <= 0
  ) {
    0L
  } else {
    as.integer(max_edges_per_target_cpp)
  }
  regulator_idx <- as.integer(match(inputs[["regulators"]], gene_names))
  target_idx <- as.integer(match(inputs[["targets"]], gene_names))
  n_targets <- length(target_idx)
  if (n_targets == 0L) {
    log_message(
      "No target genes remain after filtering",
      message_type = "error"
    )
  }
  expr <- if (inherits(grn_matrix, "dgCMatrix")) {
    density <- Matrix::nnzero(grn_matrix) /
      (nrow(grn_matrix) * ncol(grn_matrix))
    if (
      max_edges_per_target_cpp == 0L || (is.finite(density) && density < 0.05)
    ) {
      log_message(
        "Detected sparse input, using sparse-native GRNBoost2 backend",
        verbose = verbose
      )
      grn_matrix
    } else {
      log_message(
        "Detected medium-density input, using dense GRNBoost2 backend",
        verbose = verbose
      )
      as.matrix(grn_matrix)
    }
  } else {
    as.matrix(grn_matrix)
  }
  is_sparse <- inherits(expr, "dgCMatrix")
  run_chunk <- function(target_idx_chunk) {
    if (is_sparse) {
      grnboost_tree_sparse(
        expr = expr,
        regulator_idx = regulator_idx,
        target_idx = as.integer(target_idx_chunk),
        n_rounds = as.integer(max(1L, n_rounds)),
        learning_rate = learning_rate,
        max_edges_per_target = max_edges_per_target_cpp,
        max_depth = as.integer(max(1L, max_depth)),
        max_features = max_features,
        subsample = subsample,
        early_stop_window_length = as.integer(max(
          0L,
          early_stop_window_length
        )),
        random_seed = as.integer(seed %||% 1234L),
        exclude_self = isTRUE(exclude_self)
      )
    } else {
      grnboost_tree(
        expr = expr,
        regulator_idx = regulator_idx,
        target_idx = as.integer(target_idx_chunk),
        n_rounds = as.integer(max(1L, n_rounds)),
        learning_rate = learning_rate,
        max_edges_per_target = max_edges_per_target_cpp,
        max_depth = as.integer(max(1L, max_depth)),
        max_features = max_features,
        subsample = subsample,
        early_stop_window_length = as.integer(max(
          0L,
          early_stop_window_length
        )),
        random_seed = as.integer(seed %||% 1234L),
        exclude_self = isTRUE(exclude_self)
      )
    }
  }
  edge_idx <- if (cores > 1L) {
    log_message(
      "Running C++ GRNBoost2 for {.val {n_targets}} target genes",
      verbose = verbose
    )
    if (is_sparse) {
      grnboost_tree_sparse_parallel(
        expr = expr,
        regulator_idx = regulator_idx,
        target_idx = target_idx,
        n_rounds = as.integer(max(1L, n_rounds)),
        learning_rate = learning_rate,
        max_edges_per_target = max_edges_per_target_cpp,
        max_depth = as.integer(max(1L, max_depth)),
        max_features = max_features,
        subsample = subsample,
        early_stop_window_length = as.integer(max(
          0L,
          early_stop_window_length
        )),
        random_seed = as.integer(seed %||% 1234L),
        exclude_self = isTRUE(exclude_self),
        cores = as.integer(cores)
      )
    } else {
      grnboost_tree_parallel(
        expr = expr,
        regulator_idx = regulator_idx,
        target_idx = target_idx,
        n_rounds = as.integer(max(1L, n_rounds)),
        learning_rate = learning_rate,
        max_edges_per_target = max_edges_per_target_cpp,
        max_depth = as.integer(max(1L, max_depth)),
        max_features = max_features,
        subsample = subsample,
        early_stop_window_length = as.integer(max(
          0L,
          early_stop_window_length
        )),
        random_seed = as.integer(seed %||% 1234L),
        exclude_self = isTRUE(exclude_self),
        cores = as.integer(cores)
      )
    }
  } else {
    run_chunk(target_idx)
  }
  if (nrow(edge_idx) == 0L) {
    log_message(
      "C++ GRNBoost2 inference returned no edges",
      message_type = "error"
    )
  }
  adjacency <- data.frame(
    TF = gene_names[edge_idx[["regulator"]]],
    target = gene_names[edge_idx[["target"]]],
    importance = edge_idx[["importance"]],
    stringsAsFactors = FALSE
  )
  importance_norm_power <- suppressWarnings(as.numeric(importance_norm_power))
  if (
    length(importance_norm_power) == 1L &&
      !is.na(importance_norm_power) &&
      importance_norm_power > 0
  ) {
    target_total <- stats::ave(
      abs(adjacency[["importance"]]),
      adjacency[["target"]],
      FUN = sum
    )
    target_total[!is.finite(target_total) | target_total <= 0] <- 1
    adjacency[["importance"]] <- adjacency[["importance"]] /
      (target_total^importance_norm_power)
  }
  if (isTRUE(correlation_fill)) {
    adjacency <- fill_grnboost_edges(
      adjacency = adjacency,
      grn_matrix = grn_matrix,
      regulators = inputs[["regulators"]],
      targets = inputs[["targets"]],
      max_edges_per_target = max_edges_per_target,
      boost_weight = boost_weight,
      correlation_weight = correlation_weight,
      covariance_weight = covariance_weight,
      correlation_method = correlation_method
    )
  }
  adjacency <- adjacency[
    order(-adjacency[["importance"]], adjacency[["TF"]], adjacency[["target"]]),
    ,
    drop = FALSE
  ]
  write_grn_adjacency(
    adjacency,
    output_file = output_file,
    force = force
  )
}

grnboost_python <- function(
  grn_matrix,
  regulators,
  targets = NULL,
  max_edges_per_target = Inf,
  n_rounds = 5000,
  learning_rate = 0.01,
  max_depth = 3,
  max_features = 0.1,
  subsample = 0.9,
  early_stop_window_length = 25,
  output_file = NULL,
  work_dir = tempdir(),
  prefix = "grnboost2",
  envname = NULL,
  conda = "auto",
  cores = 1,
  seed = 1234,
  force = FALSE,
  verbose = TRUE
) {
  envname <- envname %||% "scenic_env"
  output_file <- output_file %||%
    file.path(work_dir, paste0(prefix, "_adj.tsv"))
  if (file.exists(output_file) && isFALSE(force)) {
    log_message("Reusing existing Python GRNBoost2 output", verbose = verbose)
    return(utils::read.delim(
      output_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  inputs <- normalize_grn_inputs(
    grn_matrix,
    regulators = regulators,
    targets = targets
  )
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  expr_csv <- file.path(work_dir, paste0(prefix, "_expression_mtx.csv"))
  regulators_file <- file.path(work_dir, paste0(prefix, "_regulators.txt"))
  raw_file <- if (identical(inputs[["targets"]], colnames(grn_matrix))) {
    output_file
  } else {
    file.path(work_dir, paste0(prefix, "_adj_raw.tsv"))
  }
  genes_keep <- union(inputs[["regulators"]], inputs[["targets"]])
  grn_sub <- grn_matrix[, genes_keep, drop = FALSE]
  utils::write.csv(as.matrix(grn_sub), file = expr_csv, row.names = TRUE)
  writeLines(inputs[["regulators"]], regulators_file, useBytes = TRUE)
  check_python(
    packages = grn_python_packages("grnboost2"),
    envname = envname,
    conda = conda,
    verbose = verbose
  )
  conda_resolved <- resolve_conda(conda)
  python_path <- conda_python(conda = conda_resolved, envname = envname)
  assert_python_runtime_switchable(
    python_path,
    restart_hint = python_runtime_restart_hint(
      envname = envname,
      modules = "scenic"
    )
  )
  configure_python_runtime(python_path)
  functions <- reticulate::import_from_path(
    "functions",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )
  functions$RunSCENICGrn(
    expression_mtx = expr_csv,
    regulators = regulators_file,
    adj_output = raw_file,
    n_rounds = as.integer(n_rounds),
    learning_rate = as.numeric(learning_rate),
    max_depth = as.integer(max_depth),
    max_features = as.numeric(max_features),
    subsample = as.numeric(subsample),
    early_stop_window_length = as.integer(early_stop_window_length),
    cores = as.integer(cores),
    seed = as.integer(seed),
    force = TRUE,
    verbose = isTRUE(verbose)
  )
  if (!identical(raw_file, output_file)) {
    filter_grn_adjacency_file(
      raw_file,
      output_file,
      inputs[["targets"]],
      verbose = verbose
    )
  }
  adjacency <- utils::read.delim(
    output_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (is.finite(max_edges_per_target) && max_edges_per_target > 0) {
    adjacency <- cap_grn_edges_per_target(adjacency, max_edges_per_target)
  }
  adjacency
}

regdiffusion_python <- function(
  grn_matrix,
  regulators,
  targets = NULL,
  output_file = NULL,
  work_dir = tempdir(),
  prefix = "regdiffusion",
  envname = NULL,
  conda = "auto",
  force = FALSE,
  verbose = TRUE,
  ...
) {
  envname <- envname %||% "scenic_env"
  output_file <- output_file %||%
    file.path(work_dir, paste0(prefix, "_adj.tsv"))
  if (file.exists(output_file) && isFALSE(force)) {
    log_message(
      "Reusing existing Python RegDiffusion output",
      verbose = verbose
    )
    return(utils::read.delim(
      output_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  inputs <- normalize_grn_inputs(
    grn_matrix,
    regulators = regulators,
    targets = targets
  )
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  expr_csv <- file.path(work_dir, paste0(prefix, "_expression_mtx.csv"))
  regulators_file <- file.path(work_dir, paste0(prefix, "_regulators.txt"))
  raw_file <- if (identical(inputs[["targets"]], colnames(grn_matrix))) {
    output_file
  } else {
    file.path(work_dir, paste0(prefix, "_adj_raw.tsv"))
  }
  genes_keep <- union(inputs[["regulators"]], inputs[["targets"]])
  grn_sub <- grn_matrix[, genes_keep, drop = FALSE]
  utils::write.csv(as.matrix(grn_sub), file = expr_csv, row.names = TRUE)
  writeLines(inputs[["regulators"]], regulators_file, useBytes = TRUE)
  check_python(
    packages = grn_python_packages("regdiffusion"),
    envname = envname,
    conda = conda,
    verbose = verbose
  )
  conda_resolved <- resolve_conda(conda)
  python_path <- conda_python(conda = conda_resolved, envname = envname)
  assert_python_runtime_switchable(
    python_path,
    restart_hint = python_runtime_restart_hint(
      envname = envname,
      modules = "regdiffusion"
    )
  )
  configure_python_runtime(python_path)
  functions <- reticulate::import_from_path(
    "functions",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )
  functions$RunSCENICRegDiffusion(
    expression_mtx = expr_csv,
    regulators = regulators_file,
    adj_output = raw_file,
    force = TRUE,
    verbose = isTRUE(verbose),
    ...
  )
  if (!identical(raw_file, output_file)) {
    filter_grn_adjacency_file(
      raw_file,
      output_file,
      inputs[["targets"]],
      verbose = verbose
    )
  }
  utils::read.delim(output_file, stringsAsFactors = FALSE, check.names = FALSE)
}

genie3 <- function(
  grn_matrix,
  regulators,
  targets = NULL,
  max_edges_per_target = 50,
  output_file = NULL,
  cores = 1,
  force = FALSE,
  verbose = TRUE,
  ...
) {
  if (!is.null(output_file) && file.exists(output_file) && isFALSE(force)) {
    log_message("Reusing existing GENIE3 output", verbose = verbose)
    return(utils::read.delim(
      output_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  check_r("GENIE3", verbose = FALSE)
  inputs <- normalize_grn_inputs(
    grn_matrix,
    regulators = regulators,
    targets = targets
  )
  expr <- Matrix::t(grn_matrix)
  weight_matrix <- GENIE3::GENIE3(
    exprMatrix = as.matrix(expr),
    regulators = inputs[["regulators"]],
    targets = inputs[["targets"]],
    nCores = as.integer(max(1L, cores)),
    ...
  )
  adjacency <- GENIE3::getLinkList(weightMatrix = weight_matrix)
  colnames(adjacency)[1:3] <- c("TF", "target", "importance")
  adjacency <- adjacency[, c("TF", "target", "importance"), drop = FALSE]
  adjacency <- cap_grn_edges_per_target(
    adjacency,
    max_edges_per_target = max_edges_per_target
  )
  if (nrow(adjacency) == 0L) {
    log_message("GENIE3 returned no edges", message_type = "error")
  }
  write_grn_adjacency(
    adjacency,
    output_file = output_file,
    force = force
  )
}

#' @title Infer gene regulatory networks with a selected backend method
#'
#' @description
#' Unified GRN inference entry point for GRNBoost2, GENIE3, RegDiffusion, and
#' GNIPLR. It returns a standardized adjacency table with at least `TF`,
#' `target`, and `importance` columns.
#'
#' @param object A Seurat object or expression matrix.
#' @param assay Assay used when `object` is a Seurat object.
#' @param layer Assay layer used when `object` is a Seurat object.
#' @param regulators Candidate transcription factor genes.
#' @param targets Optional target genes. If `NULL`, all genes are considered.
#' @param genes_in Matrix orientation for matrix inputs. `"rows"` means genes x
#' cells; `"columns"` means cells x genes.
#' @param grn_method GRN inference method.
#' @param backend Runtime backend. `"cpp"` is available for GRNBoost2 and
#' GNIPLR; `"python"` is required for RegDiffusion.
#' @param output_file Optional path where the adjacency table is written.
#' @param work_dir Working directory used by Python backends.
#' @param prefix Prefix for temporary backend files.
#' @param max_edges_per_target Maximum incoming regulator edges retained per
#' target.
#' @param n_rounds Number of boosting rounds for GRNBoost2-like inference.
#' @param learning_rate GRNBoost2-like tree ensemble learning rate.
#' @param max_depth Maximum depth of each regression tree.
#' @param max_features Fraction of candidate regulators sampled at each split.
#' @param subsample Fraction of cells sampled for each boosting round.
#' @param early_stop_window_length Out-of-bag improvement window used for
#' GRNBoost2 early stopping.
#' @param exclude_self Whether GRNBoost2-like inference excludes a target gene
#' from its own regulator feature set.
#' @param correlation_threshold Relative correlation filter used by GNIPLR.
#' @param lasso_degree Polynomial degree used by GNIPLR.
#' @param lasso_alpha LASSO regularization strength used by GNIPLR.
#' @param max_lag Maximum lag used by GNIPLR.
#' @param envname Python environment used by Python backends. If `NULL`,
#' GNIPLR uses the default `"scop_env"` through [RunGNIPLR()], while pySCENIC
#' GRNBoost2 and RegDiffusion use the isolated `"scenic_env"` environment.
#' @param conda Conda-compatible executable used by Python backends.
#' @param cores Number of workers used by supported methods.
#' @param seed Random seed passed to supported backends.
#' @param force Whether to rebuild existing `output_file`.
#' @param verbose Whether to print progress messages.
#' @param ... Additional backend-specific arguments.
#'
#' @return A data frame with standardized GRN edges.
#'
#' @examples
#' data(pancreas_sub)
#' expr <- GetAssayData5(
#'   pancreas_sub,
#'   assay = SeuratObject::DefaultAssay(pancreas_sub),
#'   layer = "counts"
#' )
#' expr <- as.matrix(expr[, seq_len(8)])
#' expr <- expr[
#'   names(sort(apply(expr, 1, stats::var), decreasing = TRUE))[seq_len(5)],
#' ]
#'
#' gniplr_grn <- RunGRN(
#'   expr,
#'   genes_in = "rows",
#'   grn_method = "gniplr",
#'   backend = "cpp",
#'   correlation_threshold = 0,
#'   lasso_degree = 1,
#'   max_lag = 1,
#'   max_edges_per_target = 2,
#'   verbose = FALSE
#' )
#' @export
RunGRN <- function(object, ...) {
  UseMethod("RunGRN", object)
}

#' @rdname RunGRN
#' @export
RunGRN.Seurat <- function(
  object,
  assay = NULL,
  layer = "counts",
  regulators = NULL,
  targets = NULL,
  grn_method = c("grnboost2", "regdiffusion", "genie3", "gniplr"),
  backend = c("cpp", "python"),
  ...
) {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  expr <- GetAssayData5(object, assay = assay, layer = layer)
  RunGRN.default(
    expr,
    regulators = regulators,
    targets = targets,
    genes_in = "rows",
    grn_method = grn_method,
    backend = backend,
    ...
  )
}

#' @rdname RunGRN
#' @export
RunGRN.matrix <- function(object, ...) {
  RunGRN.default(object, ...)
}

#' @rdname RunGRN
#' @export
RunGRN.default <- function(
  object,
  regulators = NULL,
  targets = NULL,
  genes_in = c("rows", "columns"),
  grn_method = c("grnboost2", "regdiffusion", "genie3", "gniplr"),
  backend = c("cpp", "python"),
  output_file = NULL,
  work_dir = tempdir(),
  prefix = "grn",
  max_edges_per_target = Inf,
  n_rounds = 5000,
  learning_rate = 0.01,
  max_depth = 3,
  max_features = 0.1,
  subsample = 0.9,
  early_stop_window_length = 25,
  exclude_self = TRUE,
  correlation_threshold = 0.3,
  lasso_degree = 30,
  lasso_alpha = 0.1,
  max_lag = 3,
  envname = NULL,
  conda = "auto",
  cores = 1,
  seed = 1234,
  force = FALSE,
  verbose = TRUE,
  ...
) {
  grn_method <- match.arg(grn_method)
  backend <- match.arg(backend)
  grn_matrix <- grn_matrix_from_object(object, genes_in = genes_in)
  regulators <- regulators %||% colnames(grn_matrix)
  targets <- targets %||% colnames(grn_matrix)
  if (identical(grn_method, "grnboost2")) {
    grnboost_envname <- if (identical(backend, "python")) {
      envname %||% "scenic_env"
    } else {
      envname
    }
    return(RunGRNBoost2(
      grn_matrix,
      regulators = regulators,
      targets = targets,
      genes_in = "columns",
      backend = if (identical(backend, "python")) "python" else "cpp",
      max_edges_per_target = max_edges_per_target,
      n_rounds = n_rounds,
      learning_rate = learning_rate,
      max_depth = max_depth,
      max_features = max_features,
      subsample = subsample,
      early_stop_window_length = early_stop_window_length,
      exclude_self = exclude_self,
      output_file = output_file,
      work_dir = work_dir,
      prefix = prefix,
      envname = grnboost_envname,
      conda = conda,
      cores = cores,
      seed = seed,
      force = force,
      verbose = verbose,
      ...
    ))
  }
  if (identical(grn_method, "regdiffusion")) {
    if (!identical(backend, "python")) {
      log_message(
        "{.arg grn_method = 'regdiffusion'} is only available with {.arg backend = 'python'}.",
        message_type = "error"
      )
    }
    regdiffusion_envname <- envname %||% "scenic_env"
    return(regdiffusion_python(
      grn_matrix = grn_matrix,
      regulators = regulators,
      targets = targets %||% colnames(grn_matrix),
      output_file = output_file,
      work_dir = work_dir,
      prefix = prefix,
      envname = regdiffusion_envname,
      conda = conda,
      force = force,
      verbose = verbose,
      ...
    ))
  }
  if (identical(grn_method, "gniplr")) {
    return(RunGNIPLR(
      grn_matrix,
      targets = targets,
      genes_in = "columns",
      correlation_threshold = correlation_threshold,
      lasso_degree = lasso_degree,
      lasso_alpha = lasso_alpha,
      max_lag = max_lag,
      backend = backend,
      max_edges_per_target = max_edges_per_target,
      output_file = output_file,
      work_dir = work_dir,
      prefix = prefix,
      envname = envname,
      conda = conda,
      force = force,
      verbose = verbose,
      ...
    ))
  }
  RunGENIE3(
    grn_matrix,
    regulators = regulators,
    targets = targets,
    genes_in = "columns",
    max_edges_per_target = max_edges_per_target,
    output_file = output_file,
    cores = cores,
    force = force,
    verbose = verbose,
    ...
  )
}
