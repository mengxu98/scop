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
    adjacency[["target"]] %in% targets, ,
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

grn_python_packages <- function() {
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
    is.finite(adjacency[["importance"]]) & adjacency[["importance"]] > 0, ,
    drop = FALSE
  ]
  if (nrow(adjacency) == 0L) {
    return(adjacency)
  }
  adjacency <- adjacency[
    order(adjacency[["target"]], -adjacency[["importance"]], adjacency[["TF"]]), ,
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
    order(-adjacency[["importance"]], adjacency[["TF"]], adjacency[["target"]]), ,
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
    filled[["TF"]] != filled[["target"]] & filled[["importance"]] > 0, ,
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
    dense_max_entries <- getOption(
      "scop.grnboost_dense_max_entries",
      1e7
    )
    dense_max_entries <- suppressWarnings(as.numeric(dense_max_entries[[1L]]))
    if (!is.finite(dense_max_entries) || dense_max_entries < 0) {
      dense_max_entries <- 1e7
    }
    use_dense <- as.double(nrow(grn_matrix)) * as.double(ncol(grn_matrix)) <=
      dense_max_entries
    if (!isTRUE(use_dense)) {
      log_message(
        "Detected sparse input, using sparse-native GRNBoost2 backend",
        verbose = verbose
      )
      grn_matrix
    } else {
      log_message(
        "Using dense GRNBoost2 backend within the configured matrix-size limit",
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
    order(-adjacency[["importance"]], adjacency[["TF"]], adjacency[["target"]]), ,
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
  cores <- suppressWarnings(as.integer(cores))
  if (length(cores) != 1L || is.na(cores) || cores < 1L) {
    cores <- 1L
  }
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
  PrepareEnv(
    envname = envname,
    conda = conda,
    modules = "scenic",
    verbose = verbose
  )
  check_python(
    packages = grn_python_packages(),
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
  functions <- scop_python_import("functions", convert = TRUE)
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

#' @title Infer gene regulatory networks with a selected method
#'
#' @description
#' Run transcriptomic and multiomic network methods through one entry point.
#' Complete results retain native scores, direction, sign evidence, groups,
#' runtime identity, and method parameters.
#'
#' @param object A Seurat object or expression matrix.
#' @param assay Assay used when `object` is a Seurat object.
#' @param layer Assay layer used when `object` is a Seurat object.
#' @param regulators Candidate transcription factor genes.
#' @param targets Optional target genes. If `NULL`, all genes are considered.
#' @param genes_in Matrix orientation for matrix inputs. `"rows"` means genes x
#' cells; `"columns"` means cells x genes.
#' @param methods One registered GRN method name, matched without case sensitivity.
#' @param return_type Output to return: `"seurat"`, `"edges"`, or `"result"`.
#'   With `NULL`, Seurat inputs return the updated object and matrix inputs
#'   return an edge table.
#' @param backend Runtime backend. `"cpp"` is available for GRNBoost2 and
#' GNIPLR; `"python"` is required for RegDiffusion.
#' @param resources Named method-specific resources. For PIDC, SCSGL, SCODE,
#'   SINCERITIES, and the Python GNIPLR result adapter, set \code{source_dir} to a
#'   trusted local checkout of the linked multiCSN analysis code. The linked
#'   GNIPLR implementation uses its fixed LASSO settings.
#' @param output_file Optional adjacency file for `return_type = "edges"`.
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
#' GRNBoost2 uses `"scenic_env"` on the legacy path; RegDiffusion uses
#' its own `"regdiffusion_env"` environment.
#' @param conda Conda-compatible executable used by Python backends.
#' @param cores Number of workers used by supported methods.
#' @param seed Random seed passed to supported backends.
#' @param force Whether to rebuild existing `output_file`.
#' @param verbose Whether to print progress messages.
#' @param ... Arguments passed to other functions.
#' @return The Seurat object with the complete network in its tools, a
#'   `TF`, `target`, and `importance` edge table, or a complete `GRNResult`.
#' @references
#' \url{https://github.com/mengxu98/multiCSN_manuscript/tree/42ed2c50d11ba0255326fa32ce59b97c7bb4a8b8/code}
#'
#' \url{https://github.com/zyllluck/GNIPLR}
#'
#' \url{https://doi.org/10.1016/j.cels.2017.08.014}
#'
#' \url{https://github.com/Single-Cell-Graph-Learning/scSGL}
#'
#' \url{https://github.com/hmatsu1226/SCODE}
#'
#' \url{https://github.com/CABSEL/SINCERITIES}
#'
#' @export
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
#'   methods = "gniplr",
#'   backend = "cpp",
#'   correlation_threshold = 0,
#'   lasso_degree = 1,
#'   max_lag = 1,
#'   max_edges_per_target = 2,
#'   verbose = FALSE
#' )
RunGRN <- function(object, ...) {
  argument_names <- as.character(names(as.list(match.call(expand.dots = FALSE)$...)))
  allowed <- c("methods", if ("methods" %in% argument_names) "method")
  prefixes <- setdiff(argument_names[nzchar(argument_names)], allowed)
  if (any(startsWith("methods", prefixes))) {
    stop("Use 'methods' to select a GRN method.", call. = FALSE)
  }
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
  methods = "grnboost2",
  backend = c("cpp", "python"),
  ...,
  return_type = NULL, input_mode = NULL, atac_assay = NULL, atac_layer = NULL,
  group.by = NULL, groups = NULL, pseudotime = NULL, resources = list(),
  method_args = list(), tool_name = NULL, overwrite = FALSE
) {
  RunGRN.default(
    object,
    regulators = regulators,
    targets = targets,
    genes_in = "rows",
    methods = methods,
    backend = if (missing(backend)) NULL else backend,
    ...,
    assay = assay, layer = layer, return_type = return_type, input_mode = input_mode,
    atac_assay = atac_assay, atac_layer = atac_layer, group.by = group.by,
    groups = groups, pseudotime = pseudotime, resources = resources,
    method_args = method_args, tool_name = tool_name, overwrite = overwrite
  )
}

RunGRN.default <- function(
  object, regulators = NULL, targets = NULL,
  genes_in = c("rows", "columns"),
  methods = "grnboost2",
  backend = c("cpp", "python"), output_file = NULL, work_dir = tempdir(), prefix = "grn",
  max_edges_per_target = Inf, n_rounds = 5000, learning_rate = 0.01,
  max_depth = 3, max_features = 0.1, subsample = 0.9, early_stop_window_length = 25,
  exclude_self = TRUE, correlation_threshold = 0.3, lasso_degree = 30,
  lasso_alpha = 0.1, max_lag = 3, envname = NULL, conda = "auto", cores = 1,
  seed = 1234, force = FALSE, verbose = TRUE, ...,
  return_type = NULL, input_mode = NULL, assay = NULL, layer = "counts",
  atac_assay = NULL, atac_layer = NULL, group.by = NULL, groups = NULL,
  pseudotime = NULL, resources = list(), method_args = list(),
  tool_name = NULL, overwrite = FALSE
) {
  dots <- list(...)
  argument_names <- as.character(names(dots))
  prefixes <- argument_names[nzchar(argument_names)]
  if (any(startsWith("grn_method", prefixes))) {
    stop("Use 'methods' to select a GRN method.", call. = FALSE)
  }
  if (any(startsWith("return_type", prefixes) | startsWith("return_data", prefixes) |
      startsWith("return_seurat", prefixes))) {
    stop("Use 'return_type' to select the GRN output.", call. = FALSE)
  }
  if (length(methods) != 1L || is.na(methods)) stop("Choose exactly one GRN method.", call. = FALSE)
  is_seurat <- inherits(object, "Seurat")
  return_type <- return_type %||% if (is_seurat) "seurat" else "edges"
  if (!is.character(return_type) || length(return_type) != 1L ||
      !return_type %in% c("seurat", "edges", "result")) {
    stop("return_type must be one of 'seurat', 'edges', or 'result'.", call. = FALSE)
  }
  if (return_type == "seurat" && !is_seurat) stop("return_type = 'seurat' requires a Seurat object.", call. = FALSE)
  if (return_type != "edges" && !is.null(output_file)) stop("output_file requires return_type = 'edges'.", call. = FALSE)
  spec <- grn_method_spec(methods, if (missing(backend) || is.null(backend)) NULL else backend)
  if (spec$method == "pando") {
    restore_pando_session <- grn_pando_preserve_session()
    on.exit(restore_pando_session(), add = TRUE)
  }
  scientific <- list(max_edges_per_target = max_edges_per_target, n_rounds = n_rounds,
    learning_rate = learning_rate, max_depth = max_depth, max_features = max_features,
    subsample = subsample, early_stop_window_length = early_stop_window_length,
    exclude_self = exclude_self, correlation_threshold = correlation_threshold,
    lasso_degree = lasso_degree, lasso_alpha = lasso_alpha, max_lag = max_lag)
  assay <- if (is_seurat) assay %||% SeuratObject::DefaultAssay(object) else assay
  legacy <- spec$method %in% c("grnboost2", "genie3", "gniplr") &&
    return_type == "edges" && is.null(group.by) && is.null(groups) &&
    is.null(pseudotime) && is.null(input_mode) && !length(resources) && !length(method_args) && spec$backend %in% c("cpp", "python")
  if (legacy) {
    expression <- if (is_seurat) GetAssayData5(object, assay = assay, layer = layer) else object
    old_backend <- if (missing(backend) || is.null(backend)) "cpp" else backend
    return(do.call(grn_legacy_inference, c(list(object = expression, regulators = regulators,
      targets = targets, genes_in = genes_in, grn_method = spec$method, backend = old_backend,
      output_file = output_file, work_dir = work_dir, prefix = prefix), scientific,
      list(envname = envname, conda = conda, cores = cores, seed = seed, force = force, verbose = verbose), dots)))
  }
  if (!is.list(method_args) || !is.list(resources) || (length(method_args) && (is.null(names(method_args)) || any(!nzchar(names(method_args))) || anyDuplicated(names(method_args))))) stop("method_args and resources must be named lists.", call. = FALSE)
  explicit <- intersect(names(match.call(expand.dots = FALSE)), names(scientific))
  arguments <- c(scientific[explicit], dots, method_args)
  if (anyDuplicated(names(arguments))) stop("Method arguments repeat an explicit top-level argument.", call. = FALSE)
  runtime <- arguments$runtime %||% list()
  arguments$runtime <- NULL
  runtime_isolation <- FALSE
  if (spec$method == "pando") {
    if (!is.list(runtime) || (length(runtime) && (is.null(names(runtime)) ||
        anyDuplicated(names(runtime)) || any(!names(runtime) %in% c("rscript", "r_library", "timeout"))))) {
      stop("Pando runtime options must specify an R interpreter, library or timeout.", call. = FALSE)
    }
    if (is.null(resources$genome)) stop("Pando requires resources$genome matching the chromatin coordinates.", call. = FALSE)
    runtime_isolation <- any(c("rscript", "r_library", "timeout") %in% names(runtime)) ||
      "multiCSN" %in% loadedNamespaces()
  }
  unknown <- setdiff(names(arguments), spec$allowed_args)
  if (length(unknown)) stop("Unknown argument(s) for ", spec$method, ": ", paste(unknown, collapse = ", "), call. = FALSE)
  tool_name <- tool_name %||% paste0("GRN_", spec$method)
  if (return_type == "seurat") {
    if (!is.character(tool_name) || length(tool_name) != 1L || is.na(tool_name) || !nzchar(tool_name)) {
      stop("tool_name must be one non-empty name.", call. = FALSE)
    }
    if (tool_name %in% names(object@tools) && !isTRUE(overwrite)) {
      stop("The GRN tool_name already exists. Choose another name or use overwrite = TRUE.", call. = FALSE)
    }
  }
  input_mode <- input_mode %||% if (length(spec$input_modes) == 1L) spec$input_modes else "rna"
  if (length(input_mode) != 1L || !input_mode %in% spec$input_modes) stop("Unsupported input_mode for ", spec$method, ".", call. = FALSE)
  if (input_mode == "multiome" && !is_seurat) stop("Multiome inference requires a Seurat object.", call. = FALSE)
  if (!is.null(output_file) && file.exists(output_file) && !isTRUE(force)) stop("The output_file already exists. Choose a new path or use force = TRUE.", call. = FALSE)
  if (is_seurat) {
    if (!assay %in% SeuratObject::Assays(object)) stop("Unknown RNA assay.", call. = FALSE)
    rna_layers <- SeuratObject::Layers(object[[assay]], search = layer)
    if (length(rna_layers) != 1L || !identical(rna_layers, layer)) stop("Select one exact RNA layer. Use SeuratObject::JoinLayers() to combine split layers first.", call. = FALSE)
  }
  if (input_mode == "multiome") {
    candidates <- SeuratObject::Assays(object)[vapply(SeuratObject::Assays(object), function(name) inherits(object[[name]], "ChromatinAssay"), logical(1))]
    if (is.null(atac_assay)) {
      if (length(candidates) != 1L) stop("Specify atac_assay for multiome inference.", call. = FALSE)
      atac_assay <- candidates[[1L]]
    }
    if (!atac_assay %in% SeuratObject::Assays(object)) stop("Unknown ATAC assay.", call. = FALSE)
    atac_layer <- atac_layer %||% if (spec$method == "multicsn") "data" else "counts"
    atac_layers <- SeuratObject::Layers(object[[atac_assay]], search = atac_layer)
    if (length(atac_layers) != 1L || !identical(atac_layers, atac_layer)) stop("Select one exact ATAC layer. Use SeuratObject::JoinLayers() to combine split layers first.", call. = FALSE)
  }
  expression <- if (is_seurat) GetAssayData5(object, assay = assay, layer = layer) else object
  expression <- Matrix::t(grn_matrix_from_object(expression, genes_in = genes_in))
  cells <- colnames(expression)
  genes <- rownames(expression)
  if (anyNA(cells) || anyNA(genes) || any(!nzchar(cells)) || any(!nzchar(genes)) || anyDuplicated(cells) || anyDuplicated(genes)) stop("GRN input requires unique, complete gene and cell identifiers.", call. = FALSE)
  values <- if (inherits(expression, "sparseMatrix")) expression@x else as.vector(expression)
  if (!is.numeric(values) || any(!is.finite(values))) stop("GRN expression values must be finite numbers.", call. = FALSE)
  if (is.null(group.by)) {
    group_cells <- list(all = cells)
  } else {
    labels <- if (is_seurat && is.character(group.by) && length(group.by) == 1L) {
      if (!group.by %in% colnames(object[[]])) stop("Unknown group.by metadata column.", call. = FALSE)
      stats::setNames(as.character(object[[]][cells, group.by]), cells)
    } else group.by
    if (is.null(names(labels)) || anyDuplicated(names(labels)) || !all(cells %in% names(labels))) stop("Matrix group.by must be named by cell ID.", call. = FALSE)
    labels <- as.character(labels[cells])
    if (anyNA(labels) || any(!nzchar(labels))) stop("GRN group labels must be complete.", call. = FALSE)
    group_cells <- stats::setNames(lapply(unique(labels), function(label) cells[labels == label]), unique(labels))
  }
  if (!is.null(groups)) {
    if (!length(groups) || anyDuplicated(groups) || !all(groups %in% names(group_cells))) stop("groups must select existing unique group names.", call. = FALSE)
    group_cells <- group_cells[groups]
  }
  if (return_type == "edges" && length(group_cells) != 1L) stop("Use return_type = 'result' or 'seurat' for multiple groups.", call. = FALSE)
  if (spec$method %in% c("scenicplus", "linger") && length(group_cells) != 1L && !is.null(arguments$resume_dir)) {
    stop(if (spec$method == "scenicplus") "SCENIC+ " else "LINGER ",
      "resume_dir requires exactly one selected group.", call. = FALSE)
  }
  resume_linger <- spec$method == "linger" && !is.null(arguments$resume_dir)
  if (resume_linger) {
    resume_dir <- arguments$resume_dir
    if (!is.character(resume_dir) || length(resume_dir) != 1L || is.na(resume_dir) ||
        !dir.exists(resume_dir) || !file.exists(file.path(resume_dir, "prepare.completed.json")) ||
        dir.exists(file.path(resume_dir, "prepare.completed.json"))) {
      stop("LINGER resume_dir must contain a completed prepare checkpoint.", call. = FALSE)
    }
    arguments$resume_dir <- normalizePath(resume_dir, winslash = "/", mustWork = TRUE)
  }
  cells <- cells[cells %in% unlist(group_cells, use.names = FALSE)]
  expression <- expression[, cells, drop = FALSE]
  if (is.character(pseudotime) && length(pseudotime) == 1L && is_seurat) {
    if (!pseudotime %in% colnames(object[[]])) stop("Unknown pseudotime metadata column.", call. = FALSE)
    pseudotime <- stats::setNames(object[[]][cells, pseudotime], cells)
  }
  if (!is.null(pseudotime)) {
    if (is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime, ncol = 1L, dimnames = list(names(pseudotime), "trajectory"))
    pseudotime <- as.matrix(pseudotime)
    if (!is.numeric(pseudotime) || is.null(rownames(pseudotime)) || anyDuplicated(rownames(pseudotime)) || !all(cells %in% rownames(pseudotime))) stop("pseudotime must have numeric values and unique cell IDs.", call. = FALSE)
    pseudotime <- pseudotime[cells, , drop = FALSE]
    if (any(!is.finite(pseudotime) & !is.na(pseudotime)) || any(rowSums(is.finite(pseudotime)) == 0L)) stop("Each cell needs a finite pseudotime on at least one trajectory.", call. = FALSE)
  } else if (isTRUE(spec$time)) stop("This method requires named pseudotime values.", call. = FALSE)
  regulators <- unique(as.character(regulators %||% genes))
  select_pando_targets <- spec$method == "pando" && is.null(targets)
  targets <- unique(as.character(targets %||% genes))
  if (!length(intersect(regulators, genes)) || !length(intersect(targets, genes))) stop("No candidate regulators or targets occur in the input.", call. = FALSE)
  if (length(cores) != 1L || !is.finite(cores) || cores < 1 || cores != as.integer(cores)) stop("cores must be a positive integer.", call. = FALSE)
  if (length(seed) != 1L || !is.finite(seed)) stop("seed must be a finite number.", call. = FALSE)
  config <- list(method = spec$method, backend = spec$backend, input_mode = input_mode,
    method_args = arguments, resources = resources, runtime = runtime, runtime_isolation = runtime_isolation, assay = assay, layer = layer,
    atac_assay = atac_assay, atac_layer = atac_layer, group.by = group.by, groups = group_cells,
    cores = as.integer(cores), seed = as.integer(seed), work_dir = work_dir, envname = envname,
    conda = conda, force = force, verbose = verbose, runner = spec$runner, prepare = spec$prepare,
    implementation = paste(spec$method, spec$backend, sep = "/"), package_versions = c(scop = as.character(utils::packageVersion("scop"))))
  input <- list(object = if (is_seurat) object else NULL, expression = expression,
    cells = cells, groups = group_cells, group = NULL, pseudotime = pseudotime,
    regulators = regulators, targets = if (select_pando_targets) NULL else targets, input_mode = input_mode)
  if (spec$runner %in% c("grn_run_rna", "grn_run_time")) {
    config$external_sources <- grn_external_sources(spec$method, resources)
  }
  if (spec$method == "scenicplus") grn_scenicplus_check_resources(resources)
  if (spec$method == "linger") linger_reference <- grn_linger_check_resources(input, config)
  profile <- if (spec$method == "multicsn" && input_mode == "rna") "multicsn_rna" else spec$profile
  early_runtime <- spec$method %in% c("scenicplus", "linger")
  if (early_runtime) config$runtime <- grn_resolve_runtime(config, profile)
  config$resource_identity <- switch(spec$method,
    linger = grn_linger_resource_identity(input, config, linger_reference),
    scenicplus = grn_scenicplus_resource_identity(config),
    grn_resource_identity(if (is.null(config$external_sources)) resources else
      c(resources[setdiff(names(resources), "source_dir")], list(source_files = as.list(config$external_sources)))))
  if (!early_runtime) {
    config$runtime <- if (spec$method == "pando" && !runtime_isolation) grn_pando_current_runtime(config) else grn_resolve_runtime(config, profile)
  }
  compute <- function() {
    attempt <- tempfile(paste0(prefix, "-"), tmpdir = work_dir)
    dir.create(attempt, recursive = TRUE, showWarnings = FALSE)
    run_config <- config
    run_config$work_dir <- attempt
    prepared <- if (is.null(spec$prepare)) input else get(spec$prepare, mode = "function")(input, run_config)
    had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_rng) old_rng <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(if (had_rng) assign(".Random.seed", old_rng, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    fitted <- lapply(seq_along(group_cells), function(index) {
      group_input <- prepared
      group_input$group <- names(group_cells)[[index]]
      group_input$cells <- group_cells[[index]]
      group_input$expression <- prepared$expression[, group_input$cells, drop = FALSE]
      group_input$pseudotime <- if (is.null(pseudotime)) NULL else pseudotime[group_input$cells, , drop = FALSE]
      group_config <- run_config
      group_config$work_dir <- file.path(attempt, sprintf("group-%04d", index))
      dir.create(group_config$work_dir, recursive = TRUE, showWarnings = FALSE)
      set.seed(config$seed)
      value <- get(spec$runner, mode = "function")(group_input, group_config)
      value$chains <- value$chains %||% list(status = "not_supported")
      value$gene_sets <- value$gene_sets %||% list(status = "not_supported")
      if (isTRUE(spec$allow_cell_filter)) {
        if (anyDuplicated(value$cells) || !length(value$cells) || !all(value$cells %in% group_input$cells)) stop("The GRN backend returned invalid retained cells.", call. = FALSE)
        value$input_cells <- group_input$cells
      } else if (!identical(as.character(value$cells), as.character(group_input$cells))) stop("The GRN backend returned cells from a different group.", call. = FALSE)
      value
    })
    names(fitted) <- names(group_cells)
    artifacts <- unname(do.call(c, lapply(fitted, function(value) value$execution$artifacts %||% list())))
    result <- structure(list(schema_version = "1.0", method = spec$method,
      implementation = unique(vapply(fitted, `[[`, character(1), "implementation")),
      package_versions = lapply(fitted, `[[`, "package_versions"),
      parameters = list(method_args = arguments, backend = spec$backend, seed = config$seed, cores = config$cores),
      input_summary = list(cells = cells, features = genes, groups = group_cells),
      resources = config$resource_identity, execution = list(status = "completed", artifacts = artifacts, work_dir = attempt),
      groups = fitted), class = c("GRNResult", "list"))
    grn_validate_result(result)
    result
  }
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  if (resume_linger) config$force <- TRUE
  result <- grn_cached_run(input, config, compute)
  if (return_type == "result") return(result)
  if (return_type == "seurat") {
    object@tools[[tool_name]] <- result
    return(object)
  }
  edges <- GetGRN(result)
  adjacency <- data.frame(TF = edges$source, target = edges$target, importance = edges$raw_score, stringsAsFactors = FALSE)
  write_grn_adjacency(adjacency, output_file = output_file, force = force)
}

grn_run_scenic <- function(input, config) {
  if (is.null(input$object)) {
    object <- SeuratObject::CreateSeuratObject(counts = input$expression)
    assay <- SeuratObject::DefaultAssay(object)
    layer <- "counts"
  } else {
    object <- input$object[, input$cells]
    assay <- config$assay
    layer <- config$layer
  }
  native <- do.call(RunSCENIC, c(list(object = object, assay = assay, layer = layer,
    regulators = input$regulators, targets = input$targets,
    ranking_dbs = config$resources$ranking_dbs, motif_annotations = config$resources$motif_annotation,
    work_dir = config$work_dir, backend = config$backend, cores = config$cores,
    seed = config$seed, return_seurat = FALSE, envname = config$envname,
    conda = config$conda, verbose = config$verbose), config$method_args))
  list(cells = input$cells, native = native,
    layer_specs = list(tf_gene = list(status = "computed", table_path = "adjacency", source = "TF", target = "target",
      score = "importance", source_type = "TF", target_type = "gene", score_type = "feature_importance", directed = TRUE,
      sign_column = NULL, sign_basis = "unknown", rank_by = "score")),
    gene_sets = list(status = "computed", table_path = "regulon_list"),
    parameters = native$parameters, implementation = paste0("SCENIC/", config$backend),
    package_versions = c(scop = as.character(utils::packageVersion("scop"))), execution = list(status = "completed", artifacts = list()))
}

grn_run_legacy <- function(input, config) {
  args <- config$method_args
  if (config$method == "genie3") {
    native <- do.call(RunGENIE3, c(list(object = input$expression, regulators = input$regulators,
      targets = input$targets, genes_in = "rows", cores = config$cores, verbose = config$verbose), args))
  } else {
    native <- do.call(grn_legacy_inference, c(list(object = input$expression,
      regulators = input$regulators, targets = input$targets, genes_in = "rows", grn_method = config$method,
      backend = config$backend, cores = config$cores, seed = config$seed, verbose = config$verbose), args))
  }
  attr(native, "grn_matrix") <- NULL
  relevant <- switch(config$method,
    grnboost2 = c("max_edges_per_target", "n_rounds", "learning_rate", "max_depth", "max_features", "subsample", "early_stop_window_length", "exclude_self"),
    gniplr = c("max_edges_per_target", "correlation_threshold", "lasso_degree", "lasso_alpha", "max_lag"),
    genie3 = "max_edges_per_target")
  defaults <- lapply(formals(grn_legacy_inference)[relevant], eval, envir = baseenv())
  if (config$method == "genie3") defaults <- c(defaults, lapply(formals(GENIE3::GENIE3)[c("treeMethod", "K", "nTrees")], eval, envir = baseenv()))
  parameters <- utils::modifyList(defaults, args)
  parameters$cores <- config$cores
  parameters$seed <- config$seed
  spec <- list(status = "computed", table_path = character(), source = "TF", target = "target",
      score = "importance", source_type = "TF", target_type = "gene", score_type = "importance", directed = TRUE,
      sign_column = NULL, sign_basis = "unknown", rank_by = "score")
  if (config$method == "gniplr") {
    spec$rows <- which(native$TF %in% input$regulators & native$target %in% input$targets)
    parameters$candidate_filter <- "postfit_view"
  }
  list(cells = input$cells, native = native, layer_specs = list(tf_gene = spec),
    parameters = parameters, implementation = paste0("scop::", config$method, "/", config$backend),
    package_versions = c(scop = as.character(utils::packageVersion("scop"))), execution = list(status = "completed", artifacts = list()))
}

grn_legacy_inference <- function(
  object,
  regulators = NULL,
  targets = NULL,
  genes_in = c("rows", "columns"),
  grn_method = c("grnboost2", "genie3", "gniplr"),
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
