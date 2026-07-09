#' @title Infer gene regulatory networks with GNIPLR
#'
#' @param object A Seurat object or expression matrix.
#' @param assay Assay used when `object` is a Seurat object.
#' @param layer Assay layer used when `object` is a Seurat object.
#' @param targets Optional target genes. If `NULL`, all genes are considered.
#' @param genes_in Matrix orientation for matrix inputs. `"rows"` means genes x
#' cells; `"columns"` means cells x genes.
#' @param correlation_threshold Relative correlation filter used before lagged
#' regression. A gene pair is tested when its absolute correlation is at least
#' this fraction of the maximum absolute correlation for the regulator.
#' @param lasso_degree Polynomial degree used by the LASSO projection step.
#' @param lasso_alpha LASSO regularization strength used by the projection step.
#' @param max_lag Maximum lag used by Granger-style lagged regression. Values
#' above `3` are capped at `3` to match the GNIPLR reference implementation.
#' @param backend Runtime backend. Supports `"cpp"` and `"python"`.
#' @param max_edges_per_target Maximum incoming regulator edges retained per
#' target. The default `Inf` keeps all positive-importance links.
#' @param output_file Optional path where the adjacency table is written.
#' @param work_dir Working directory used by the Python backend.
#' @param prefix Prefix for temporary files.
#' @param envname Python environment used by the Python backend.
#' @param conda Conda-compatible executable used by the Python backend.
#' @param force Whether to rebuild existing `output_file`.
#' @param verbose Whether to print progress messages.
#' @param ... Additional backend-specific arguments.
#'
#' @references
#' Zhang Y, Chang X, Liu X. Inference of gene regulatory networks using
#' pseudo-time series data. Bioinformatics. 2021;37(16):2423-2431.
#' doi:10.1093/bioinformatics/btab099.
#'
#' @return A data frame with columns `TF`, `target`, `importance`, and
#' `pvalue`. The original GNIPLR p-value matrix is stored in
#' `attr(result, "grn_matrix")` when the network is newly inferred.
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
#' grn <- RunGNIPLR(
#'   expr,
#'   genes_in = "rows",
#'   correlation_threshold = 0,
#'   lasso_degree = 1,
#'   max_lag = 1,
#'   max_edges_per_target = 2
#' )
#' head(grn)
#' attr(grn, "grn_matrix")[1:3, 1:3]
#' @export
RunGNIPLR <- function(object, ...) {
  UseMethod("RunGNIPLR", object)
}

#' @rdname RunGNIPLR
#' @export
RunGNIPLR.Seurat <- function(
  object,
  assay = NULL,
  layer = "counts",
  targets = NULL,
  correlation_threshold = 0.3,
  lasso_degree = 30,
  lasso_alpha = 0.1,
  max_lag = 3,
  backend = c("cpp", "python"),
  max_edges_per_target = Inf,
  output_file = NULL,
  work_dir = tempdir(),
  prefix = "gniplr",
  envname = "scop_env",
  conda = "auto",
  force = FALSE,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  expr <- GetAssayData5(object, assay = assay, layer = layer)
  RunGNIPLR.default(
    expr,
    targets = targets,
    genes_in = "rows",
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
  )
}

#' @rdname RunGNIPLR
#' @export
RunGNIPLR.matrix <- function(object, ...) {
  RunGNIPLR.default(object, ...)
}

#' @rdname RunGNIPLR
#' @export
RunGNIPLR.default <- function(
  object,
  targets = NULL,
  genes_in = c("rows", "columns"),
  correlation_threshold = 0.3,
  lasso_degree = 30,
  lasso_alpha = 0.1,
  max_lag = 3,
  backend = c("cpp", "python"),
  max_edges_per_target = Inf,
  output_file = NULL,
  work_dir = tempdir(),
  prefix = "gniplr",
  envname = "scop_env",
  conda = "auto",
  force = FALSE,
  verbose = TRUE,
  ...
) {
  backend <- match.arg(backend)
  genes_in <- match.arg(genes_in)
  if (!is.null(output_file) && file.exists(output_file) && isFALSE(force)) {
    return(utils::read.delim(
      output_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  expr <- if (inherits(object, "Matrix")) {
    as.matrix(object)
  } else {
    as.matrix(object)
  }
  if (identical(genes_in, "columns")) {
    expr <- t(expr)
  }
  if (is.null(rownames(expr)) || any(!nzchar(rownames(expr)))) {
    log_message(
      "{.arg object} must have gene names as row names",
      message_type = "error"
    )
  }
  if (is.null(colnames(expr)) || any(!nzchar(colnames(expr)))) {
    colnames(expr) <- paste0("cell", seq_len(ncol(expr)))
  }
  if (!is.numeric(expr)) {
    storage.mode(expr) <- "double"
  }
  expr[!is.finite(expr)] <- 0
  targets <- intersect(targets %||% rownames(expr), rownames(expr))
  if (length(targets) == 0L) {
    log_message(
      "No target genes are present in the expression matrix",
      message_type = "error"
    )
  }
  log_message(
    "Running {.pkg GNIPLR} with {.arg backend = {backend}} on {.val {nrow(expr)}} genes and {.val {ncol(expr)}} cells",
    verbose = verbose
  )
  if (identical(backend, "cpp")) {
    res <- gniplr_cpp(
      expression = expr,
      target_idx = as.integer(match(targets, rownames(expr))),
      correlation_threshold = as.numeric(correlation_threshold),
      lasso_degree = as.integer(lasso_degree),
      lasso_alpha = as.numeric(lasso_alpha),
      max_lag = as.integer(max_lag)
    )
    adjacency <- as.data.frame(res[["adjacency"]], stringsAsFactors = FALSE)
    if (nrow(adjacency) > 0L) {
      adjacency[["TF"]] <- rownames(expr)[adjacency[["regulator"]]]
      adjacency[["target"]] <- rownames(expr)[adjacency[["target"]]]
      adjacency <- adjacency[, c("TF", "target", "importance", "pvalue"), drop = FALSE]
    }
  } else {
    check_python(
      packages = c(
        "numpy>=0",
        "pandas>=0",
        "scipy>=0",
        "scikit-learn>=0",
        "statsmodels>=0"
      ),
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
        modules = "scanpy"
      )
    )
    configure_python_runtime(python_path)
    functions <- reticulate::import_from_path(
      "functions",
      path = system.file("python", package = "scop", mustWork = TRUE),
      convert = TRUE
    )
    res <- functions$RunGNIPLR(
      expression = expr,
      genes = rownames(expr),
      targets = targets,
      correlation_threshold = as.numeric(correlation_threshold),
      lasso_degree = as.integer(lasso_degree),
      lasso_alpha = as.numeric(lasso_alpha),
      max_lag = as.integer(max_lag),
      verbose = isTRUE(verbose),
      ...
    )
    adjacency <- as.data.frame(res[["adjacency"]], stringsAsFactors = FALSE)
  }
  if (nrow(adjacency) == 0L) {
    adjacency <- data.frame(
      TF = character(),
      target = character(),
      importance = numeric(),
      pvalue = numeric(),
      stringsAsFactors = FALSE
    )
  }
  adjacency <- adjacency[, c("TF", "target", "importance", "pvalue"), drop = FALSE]
  adjacency[["importance"]] <- suppressWarnings(as.numeric(adjacency[["importance"]]))
  adjacency[["pvalue"]] <- suppressWarnings(as.numeric(adjacency[["pvalue"]]))
  adjacency <- cap_grn_edges_per_target(
    adjacency,
    max_edges_per_target = max_edges_per_target
  )
  grn_matrix <- as.matrix(res[["grn_matrix"]])
  dimnames(grn_matrix) <- list(rownames(expr), rownames(expr))
  attr(adjacency, "grn_matrix") <- grn_matrix
  write_grn_adjacency(
    adjacency,
    output_file = output_file,
    force = force
  )
}
