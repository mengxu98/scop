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
#' @param envname Python environment used by Python backends. If `NULL`,
#' `backend = "python"` uses the isolated `"scenic_env"` environment.
#' @param conda Conda-compatible executable used by Python backends.
#' @param cores Number of workers used by native/Python GRNBoost2.
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
