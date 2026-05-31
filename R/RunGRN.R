#' @title Infer gene regulatory networks
#'
#' @description
#' Run GRNBoost2-like or GENIE3 regulatory network inference and return a
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
#' @param backend Runtime backend. GRNBoost2 supports `"cpp"` and `"python"`;
#' GENIE3 uses the R package backend.
#' @param max_edges_per_target Maximum incoming regulator edges retained per
#' target. The default `Inf` keeps all positive-importance links, matching
#' arboreto GRNBoost2 output.
#' @param n_rounds Number of native boosting rounds for GRNBoost2-like tree
#' ensemble inference. The default follows arboreto `SGBM_KWARGS`.
#' @param learning_rate Native GRNBoost2-like tree ensemble learning rate.
#' @param max_depth Maximum depth of each native regression tree.
#' @param max_features Fraction of candidate regulators sampled at each native
#' tree split.
#' @param subsample Fraction of cells sampled for each native boosting round.
#' @param early_stop_window_length Out-of-bag improvement window used for native
#' GRNBoost2 early stopping, matching arboreto's `EarlyStopMonitor`.
#' @param exclude_self Whether native GRNBoost2-like inference excludes a target
#' gene from its own regulator feature set.
#' @param importance_norm_power Power used to normalize native GRNBoost2-like
#' edge importance by the total importance for each target. Set to `0` to
#' disable normalization.
#' @param correlation_fill Whether native GRNBoost2-like inference should fill
#' missing TF-target candidates with expression-correlation support before
#' applying `max_edges_per_target`.
#' @param correlation_weight Weight of the correlation-fill score relative to
#' native boosting importance.
#' @param boost_weight Weight of the native boosting score used when
#' correlation fill is enabled.
#' @param covariance_weight Weight of the covariance score used when
#' correlation fill is enabled.
#' @param correlation_method Correlation method used for native GRNBoost2-like
#' candidate filling.
#' @param output_file Optional path where the adjacency table is written.
#' @param work_dir Working directory used by Python backends.
#' @param prefix Prefix for temporary Python backend files.
#' @param envname Python environment used by Python backends.
#' @param conda Conda-compatible executable used by Python backends.
#' @param prepare_env Whether to prepare the Python environment before running
#' Python backends.
#' @param cores Number of workers used by native/Python GRNBoost2 or R GENIE3.
#' @param seed Random seed passed to supported backends.
#' @param force Whether to rebuild existing `output_file`.
#' @param verbose Whether to print progress messages.
#' @param ... Additional backend-specific arguments.
#'
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
  envname = "scenic_env",
  conda = "auto",
  prepare_env = FALSE,
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
    prepare_env = prepare_env,
    cores = cores,
    seed = seed,
    force = force,
    verbose = verbose,
    ...
  )
}

#' @rdname RunGRNBoost2
#' @export
RunGRNBoost2.matrix <- function(object, ...) {
  RunGRNBoost2.default(object, ...)
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
  envname = "scenic_env",
  conda = "auto",
  prepare_env = FALSE,
  cores = 1,
  force = FALSE,
  verbose = TRUE,
  ...
) {
  backend <- match.arg(backend)
  grn_matrix <- scenic_grn_matrix_from_object(object, genes_in = genes_in)
  regulators <- regulators %||% colnames(grn_matrix)
  targets <- targets %||% colnames(grn_matrix)
  if (identical(backend, "python")) {
    return(grnboost_python(
      grn_matrix = grn_matrix,
      regulators = regulators,
      targets = targets,
      max_edges_per_target = max_edges_per_target,
      n_rounds = n_rounds,
      output_file = output_file,
      work_dir = work_dir,
      prefix = prefix,
      envname = envname,
      conda = conda,
      prepare_env = prepare_env,
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

#' @rdname RunGRNBoost2
#' @export
RunGENIE3 <- function(object, ...) {
  UseMethod("RunGENIE3", object)
}

#' @rdname RunGRNBoost2
#' @export
RunGENIE3.Seurat <- function(
  object,
  assay = NULL,
  layer = "counts",
  regulators = NULL,
  targets = NULL,
  max_edges_per_target = 50,
  output_file = NULL,
  cores = 1,
  force = FALSE,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  counts <- GetAssayData5(object, assay = assay, layer = layer)
  RunGENIE3.default(
    counts,
    regulators = regulators,
    targets = targets,
    genes_in = "rows",
    max_edges_per_target = max_edges_per_target,
    output_file = output_file,
    cores = cores,
    force = force,
    verbose = verbose,
    ...
  )
}

#' @rdname RunGRNBoost2
#' @export
RunGENIE3.matrix <- function(object, ...) {
  RunGENIE3.default(object, ...)
}

#' @rdname RunGRNBoost2
#' @export
RunGENIE3.default <- function(
  object,
  regulators = NULL,
  targets = NULL,
  genes_in = c("rows", "columns"),
  max_edges_per_target = 50,
  output_file = NULL,
  cores = 1,
  force = FALSE,
  verbose = TRUE,
  ...
) {
  grn_matrix <- scenic_grn_matrix_from_object(object, genes_in = genes_in)
  regulators <- regulators %||% colnames(grn_matrix)
  targets <- targets %||% colnames(grn_matrix)
  genie3(
    grn_matrix = grn_matrix,
    regulators = regulators,
    targets = targets,
    max_edges_per_target = max_edges_per_target,
    output_file = output_file,
    cores = cores,
    force = force,
    verbose = verbose,
    ...
  )
}

scenic_grn_matrix_from_object <- function(object, genes_in = c("rows", "columns")) {
  genes_in <- match.arg(genes_in)
  if (!inherits(object, "Matrix") && !is.matrix(object) && !is.data.frame(object)) {
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
    mat <- Matrix::t(mat)
  }
  mat <- as.matrix(mat)
  if (nrow(mat) < 2L || ncol(mat) < 2L) {
    log_message(
      "GRN inference requires at least two samples and two genes",
      message_type = "error"
    )
  }
  mat
}

scenic_normalize_grn_inputs <- function(
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
    log_message("No candidate regulators are present in the expression matrix", message_type = "error")
  }
  if (length(targets) == 0L) {
    log_message("No target genes are present in the expression matrix", message_type = "error")
  }
  list(regulators = regulators, targets = targets)
}

scenic_write_grn_adjacency <- function(adjacency, output_file = NULL, force = FALSE) {
  if (is.null(output_file)) {
    return(adjacency)
  }
  if (file.exists(output_file) && isFALSE(force)) {
    return(utils::read.delim(output_file, stringsAsFactors = FALSE, check.names = FALSE))
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

scenic_cap_edges_per_target <- function(adjacency, max_edges_per_target = Inf) {
  max_edges_per_target <- suppressWarnings(as.numeric(max_edges_per_target))
  adjacency <- adjacency[is.finite(adjacency[["importance"]]) & adjacency[["importance"]] > 0, , drop = FALSE]
  if (nrow(adjacency) == 0L) {
    return(adjacency)
  }
  adjacency <- adjacency[order(adjacency[["target"]], -adjacency[["importance"]], adjacency[["TF"]]), , drop = FALSE]
  if (length(max_edges_per_target) == 1L && is.finite(max_edges_per_target) && max_edges_per_target > 0) {
    adjacency <- do.call(
      rbind,
      lapply(split(adjacency, adjacency[["target"]]), utils::head, as.integer(max_edges_per_target))
    )
  }
  rownames(adjacency) <- NULL
  adjacency[order(-adjacency[["importance"]], adjacency[["TF"]], adjacency[["target"]]), , drop = FALSE]
}

scenic_normalize_target_importance <- function(adjacency, power = 0.25) {
  power <- suppressWarnings(as.numeric(power))
  if (length(power) != 1L || is.na(power) || power <= 0) {
    return(adjacency)
  }
  target_total <- ave(abs(adjacency[["importance"]]), adjacency[["target"]], FUN = sum)
  target_total[!is.finite(target_total) | target_total <= 0] <- 1
  adjacency[["importance"]] <- adjacency[["importance"]] / (target_total^power)
  adjacency
}

scenic_fill_grnboost_edges <- function(
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
  if (length(correlation_weight) != 1L || is.na(correlation_weight) || correlation_weight < 0) {
    correlation_weight <- 0
  }
  if (length(covariance_weight) != 1L || is.na(covariance_weight) || covariance_weight < 0) {
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
  filled <- filled[filled[["TF"]] != filled[["target"]] & filled[["importance"]] > 0, , drop = FALSE]
  if (nrow(filled) == 0L) {
    return(adjacency)
  }
  filled[["edge"]] <- paste(filled[["TF"]], filled[["target"]], sep = "\r")
  adjacency[["edge"]] <- paste(adjacency[["TF"]], adjacency[["target"]], sep = "\r")
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
  filled[["importance"]] <- boost_weight * boost_score +
    correlation_weight * cor_score +
    covariance_weight * cov_score
  filled <- filled[, c("TF", "target", "importance"), drop = FALSE]
  scenic_cap_edges_per_target(filled, max_edges_per_target = max_edges_per_target)
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
    log_message("Reusing existing native GRNBoost2 output", verbose = verbose)
    return(utils::read.delim(output_file, stringsAsFactors = FALSE, check.names = FALSE))
  }
  correlation_method <- match.arg(correlation_method)
  cores <- suppressWarnings(as.integer(cores))
  if (length(cores) != 1L || is.na(cores) || cores < 1L) {
    cores <- 1L
  }
  inputs <- scenic_normalize_grn_inputs(grn_matrix, regulators = regulators, targets = targets)
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
  expr <- as.matrix(grn_matrix)
  run_chunk <- function(target_idx_chunk) {
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
      early_stop_window_length = as.integer(max(0L, early_stop_window_length)),
      random_seed = as.integer(seed %||% 1234L),
      exclude_self = isTRUE(exclude_self)
    )
  }
  edge_idx <- if (cores > 1L && length(target_idx) > 1L) {
    workers <- min(cores, length(target_idx))
    chunks <- split(target_idx, rep(seq_len(workers), length.out = length(target_idx)))
    chunks <- chunks[lengths(chunks) > 0L]
    do.call(
      rbind,
      parallel::mclapply(chunks, run_chunk, mc.cores = workers, mc.preschedule = TRUE)
    )
  } else {
    run_chunk(target_idx)
  }
  if (nrow(edge_idx) == 0L) {
    log_message("Native GRNBoost2-like inference returned no edges", message_type = "error")
  }
  adjacency <- data.frame(
    TF = gene_names[edge_idx[["regulator"]]],
    target = gene_names[edge_idx[["target"]]],
    importance = edge_idx[["importance"]],
    stringsAsFactors = FALSE
  )
  adjacency <- scenic_normalize_target_importance(adjacency, power = importance_norm_power)
  if (isTRUE(correlation_fill)) {
    adjacency <- scenic_fill_grnboost_edges(
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
  adjacency <- adjacency[order(-adjacency[["importance"]], adjacency[["TF"]], adjacency[["target"]]), , drop = FALSE]
  scenic_write_grn_adjacency(adjacency, output_file = output_file, force = force)
}

grnboost_python <- function(
  grn_matrix,
  regulators,
  targets = NULL,
  max_edges_per_target = Inf,
  n_rounds = 5000,
  output_file = NULL,
  work_dir = tempdir(),
  prefix = "grnboost2",
  envname = "scenic_env",
  conda = "auto",
  prepare_env = FALSE,
  cores = 1,
  seed = 1234,
  force = FALSE,
  verbose = TRUE
) {
  output_file <- output_file %||% file.path(work_dir, paste0(prefix, "_adj.tsv"))
  if (file.exists(output_file) && isFALSE(force)) {
    log_message("Reusing existing Python GRNBoost2 output", verbose = verbose)
    return(utils::read.delim(output_file, stringsAsFactors = FALSE, check.names = FALSE))
  }
  inputs <- scenic_normalize_grn_inputs(grn_matrix, regulators = regulators, targets = targets)
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
  functions <- scenic_python_functions(
    envname = envname,
    conda = conda,
    prepare_env = prepare_env,
    modules = "scenic",
    packages = scenic_py_pkgs("grnboost2"),
    verbose = verbose
  )
  functions$RunSCENICGrn(
    expression_mtx = expr_csv,
    regulators = regulators_file,
    adj_output = raw_file,
    n_rounds = as.integer(n_rounds),
    cores = as.integer(cores),
    seed = as.integer(seed),
    force = TRUE,
    verbose = isTRUE(verbose)
  )
  if (!identical(raw_file, output_file)) {
    scenic_flt_adj(raw_file, output_file, inputs[["targets"]], verbose = verbose)
  }
  adjacency <- utils::read.delim(output_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (is.finite(max_edges_per_target) && max_edges_per_target > 0) {
    adjacency <- scenic_cap_edges_per_target(adjacency, max_edges_per_target)
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
  envname = "scenic_env",
  conda = "auto",
  prepare_env = FALSE,
  force = FALSE,
  verbose = TRUE,
  ...
) {
  output_file <- output_file %||% file.path(work_dir, paste0(prefix, "_adj.tsv"))
  if (file.exists(output_file) && isFALSE(force)) {
    log_message("Reusing existing Python RegDiffusion output", verbose = verbose)
    return(utils::read.delim(output_file, stringsAsFactors = FALSE, check.names = FALSE))
  }
  inputs <- scenic_normalize_grn_inputs(grn_matrix, regulators = regulators, targets = targets)
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
  functions <- scenic_python_functions(
    envname = envname,
    conda = conda,
    prepare_env = prepare_env,
    modules = c("scenic", "regdiffusion"),
    packages = scenic_py_pkgs("regdiffusion"),
    verbose = verbose
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
    scenic_flt_adj(raw_file, output_file, inputs[["targets"]], verbose = verbose)
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
    return(utils::read.delim(output_file, stringsAsFactors = FALSE, check.names = FALSE))
  }
  check_r("GENIE3", verbose = FALSE)
  inputs <- scenic_normalize_grn_inputs(grn_matrix, regulators = regulators, targets = targets)
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
  adjacency <- scenic_cap_edges_per_target(adjacency, max_edges_per_target = max_edges_per_target)
  if (nrow(adjacency) == 0L) {
    log_message("GENIE3 returned no edges", message_type = "error")
  }
  scenic_write_grn_adjacency(adjacency, output_file = output_file, force = force)
}

scenic_python_functions <- function(
  envname,
  conda,
  prepare_env = FALSE,
  modules = "scenic",
  packages = scenic_py_pkgs("grnboost2"),
  verbose = TRUE
) {
  if (isTRUE(prepare_env)) {
    scenic_prepare_python_env(
      envname = envname,
      conda = conda,
      packages = packages,
      modules = modules,
      verbose = verbose
    )
  }
  check_python(
    packages = packages,
    envname = envname,
    conda = conda,
    verbose = verbose
  )
  conda_resolved <- resolve_conda(conda)
  python_path <- conda_python(conda = conda_resolved, envname = envname)
  assert_python_runtime_switchable(
    python_path,
    restart_hint = scenic_runtime_restart_hint(envname = envname)
  )
  configure_python_runtime(python_path)
  reticulate::import_from_path(
    "functions",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )
}

scenic_run_grn_method <- function(
  grn_matrix,
  regulators,
  targets = NULL,
  grn_method = c("grnboost2", "regdiffusion", "genie3"),
  backend = c("cpp", "python"),
  output_file,
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
  envname = "scenic_env",
  conda = "auto",
  prepare_env = FALSE,
  cores = 1,
  seed = 1234,
  force = FALSE,
  verbose = TRUE
) {
  grn_method <- match.arg(grn_method)
  backend <- match.arg(backend)
  if (identical(grn_method, "grnboost2")) {
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
      envname = envname,
      conda = conda,
      prepare_env = prepare_env,
      cores = cores,
      seed = seed,
      force = force,
      verbose = verbose
    ))
  }
  if (identical(grn_method, "regdiffusion")) {
    if (!identical(backend, "python")) {
      log_message(
        "{.arg grn_method = 'regdiffusion'} is only available with {.arg backend = 'python'}.",
        message_type = "error"
      )
    }
    return(regdiffusion_python(
      grn_matrix = grn_matrix,
      regulators = regulators,
      targets = targets %||% colnames(grn_matrix),
      output_file = output_file,
      work_dir = work_dir,
      prefix = prefix,
      envname = envname,
      conda = conda,
      prepare_env = prepare_env,
      force = force,
      verbose = verbose
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
    verbose = verbose
  )
}
