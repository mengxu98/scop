#' @title Infer gene regulatory networks with GRNBoost2
#'
#' @description
#' Run GRNBoost2-like regulatory network inference and return a
#' standardized adjacency table with columns `TF`, `target`, and `importance`.
#'
#' @param object A Seurat object or expression matrix.
#' @param assay Assay used when `object` is a Seurat object.
#' @param layer Assay layer used when `object` is a Seurat object.
#' @param regulators Candidate transcription factor genes.
#' @param targets Optional target genes. If `NULL`, all genes in the GRN matrix
#' are considered as candidate targets.
#' @param genes_in Matrix orientation for matrix inputs. `"rows"` means genes x
#' cells; `"columns"` means cells x genes.
#' @param backend Runtime backend. Supports `"cpp"` and `"python"`.
#' @param max_edges_per_target Maximum incoming regulator edges retained per
#' target. The default `Inf` keeps all positive-importance links, matching
#' arboreto GRNBoost2 output.
#' @param n_rounds Number of boosting rounds for GRNBoost2-like tree ensemble
#' inference. The default follows arboreto `SGBM_KWARGS`.
#' @param learning_rate GRNBoost2-like tree ensemble learning rate.
#' @param max_depth Maximum depth of each regression tree.
#' @param max_features Fraction of candidate regulators sampled at each tree
#' split.
#' @param subsample Fraction of cells sampled for each boosting round.
#' @param early_stop_window_length Out-of-bag improvement window used for
#' GRNBoost2 early stopping.
#' @param exclude_self Whether GRNBoost2-like inference excludes a target
#' gene from its own regulator feature set.
#' @param importance_norm_power Power used to normalize GRNBoost2-like
#' edge importance by the total importance for each target. Set to `0` to
#' disable normalization.
#' @param correlation_fill Whether GRNBoost2-like inference should fill
#' missing TF-target candidates with expression-correlation support before
#' applying `max_edges_per_target`.
#' @param correlation_weight Weight of the correlation-fill score relative to
#' boosting importance.
#' @param boost_weight Weight of the boosting score used when
#' correlation fill is enabled.
#' @param covariance_weight Weight of the covariance score used when
#' correlation fill is enabled.
#' @param correlation_method Correlation method used for GRNBoost2-like
#' candidate filling.
#' @param output_file Optional path where the adjacency table is written.
#' @param work_dir Working directory used by Python backends.
#' @param prefix Prefix for temporary Python backend files.
#' @param envname Python environment used by Python backends. If `NULL`,
#' `backend = "python"` uses the isolated `"scenic_env"` environment.
#' @param conda Conda-compatible executable used by Python backends.
#' @param cores Number of workers used by GRNBoost2.
#' @param seed Random seed passed to supported backends.
#' @param force Whether to rebuild existing `output_file`.
#' @param verbose Whether to print progress messages.
#' @param ... Arguments passed to other functions.
#' @return A data frame with columns `TF`, `target`, and `importance`.
#' @export
RunGRNBoost2 <- function(object, ...) {
  UseMethod("RunGRNBoost2", object)
}

#' @rdname RunGRNBoost2
#' @export
RunGRNBoost2.Seurat <- function(
  object,
  assay = NULL,
  layer = "counts",
  regulators = NULL,
  targets = NULL,
  backend = c("cpp", "python"),
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
  output_file = NULL,
  work_dir = tempdir(),
  prefix = "grnboost2",
  envname = NULL,
  conda = "auto",
  cores = 1,
  seed = 1234,
  force = FALSE,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  counts <- GetAssayData5(object, assay = assay, layer = layer)
  RunGRNBoost2.default(
    counts,
    regulators = regulators,
    targets = targets,
    genes_in = "rows",
    backend = backend,
    max_edges_per_target = max_edges_per_target,
    n_rounds = n_rounds,
    learning_rate = learning_rate,
    max_depth = max_depth,
    max_features = max_features,
    subsample = subsample,
    early_stop_window_length = early_stop_window_length,
    exclude_self = exclude_self,
    importance_norm_power = importance_norm_power,
    correlation_fill = correlation_fill,
    boost_weight = boost_weight,
    correlation_weight = correlation_weight,
    covariance_weight = covariance_weight,
    correlation_method = correlation_method,
    output_file = output_file,
    work_dir = work_dir,
    prefix = prefix,
    envname = envname,
    conda = conda,
    cores = cores,
    seed = seed,
    force = force,
    verbose = verbose,
    ...
  )
}

#' @rdname RunGRNBoost2
#' @export
RunGRNBoost2.default <- function(
  object,
  regulators = NULL,
  targets = NULL,
  genes_in = c("rows", "columns"),
  backend = c("cpp", "python"),
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
  seed = 1234,
  output_file = NULL,
  work_dir = tempdir(),
  prefix = "grnboost2",
  envname = NULL,
  conda = "auto",
  cores = 1,
  force = FALSE,
  verbose = TRUE,
  ...
) {
  backend <- match.arg(backend)
  grn_matrix <- grn_matrix_from_object(object, genes_in = genes_in)
  regulators <- regulators %||% colnames(grn_matrix)
  targets <- targets %||% colnames(grn_matrix)
  if (identical(backend, "python")) {
    envname <- envname %||% "scenic_env"
    return(grnboost_python(
      grn_matrix = grn_matrix,
      regulators = regulators,
      targets = targets,
      max_edges_per_target = max_edges_per_target,
      n_rounds = n_rounds,
      learning_rate = learning_rate,
      max_depth = max_depth,
      max_features = max_features,
      subsample = subsample,
      early_stop_window_length = early_stop_window_length,
      output_file = output_file,
      work_dir = work_dir,
      prefix = prefix,
      envname = envname,
      conda = conda,
      cores = cores,
      seed = seed,
      force = force,
      verbose = verbose
    ))
  }
  grnboost(
    grn_matrix = grn_matrix,
    regulators = regulators,
    targets = targets,
    max_edges_per_target = max_edges_per_target,
    n_rounds = n_rounds,
    learning_rate = learning_rate,
    max_depth = max_depth,
    max_features = max_features,
    subsample = subsample,
    early_stop_window_length = early_stop_window_length,
    exclude_self = exclude_self,
    importance_norm_power = importance_norm_power,
    correlation_fill = correlation_fill,
    boost_weight = boost_weight,
    correlation_weight = correlation_weight,
    covariance_weight = covariance_weight,
    correlation_method = correlation_method,
    cores = cores,
    seed = seed,
    output_file = output_file,
    force = force,
    verbose = verbose
  )
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
