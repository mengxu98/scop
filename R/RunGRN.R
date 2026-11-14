#' @title Infer gene regulatory networks
#'
#' @md
#' @inheritParams RunGRNBoost2
#' @inheritParams RunGNIPLR
#' @param targets Candidate target genes. `NULL` uses all genes for most methods;
#'   Pando selects variable genes.
#' @param methods One method name (case insensitive): `multiCSN`, `GENIE3`,
#'   `GRNBoost2`, `GNIPLR`, `PPCOR`, `PIDC`, `SCORPION`, `SCSGL`, `LEAP`,
#'   `SCODE`, `SINCERITIES`, `Pando`, `LINGER`, `SCENIC+`, `RegDiffusion`, or `SCENIC`.
#' @param backend Method-specific implementation. `NULL` uses the method default.
#' @param return_type `"seurat"`, `"edges"`, or `"result"`. The default is
#'   `"seurat"` for Seurat input and `"edges"` for matrix input.
#' @param resources Named backend resources. Set `source_dir` to the
#'   analysis checkout in References for PIDC, SCSGL, SCODE,
#'   and Python GNIPLR.
#' @param output_file Optional adjacency file for `return_type = "edges"`.
#' @param work_dir Directory for method outputs and logs.
#' @param prefix Prefix for run directories.
#' @param max_edges_per_target Maximum incoming edges per target for methods that support it.
#' @param envname Python environment. `NULL` uses the method default.
#' @param cores Number of workers for methods that support parallel execution.
#' @param force Overwrite an existing `output_file`.
#' @param input_mode `"rna"` or `"multiome"`. Multiome input requires a Seurat object.
#' @param atac_assay,atac_layer Chromatin assay and exact layer. Join split layers first.
#' @param group.by Metadata column or cell-named vector defining separate networks.
#' @param groups Group labels to fit, in the requested order.
#' @param pseudotime Metadata column, cell-named numeric vector, or cell-by-trajectory matrix.
#' @param method_args Named backend arguments. `runtime` accepts `python` and
#'   `timeout` for subprocesses. `resume_dir` resumes one LINGER or SCENIC+ group.
#' @param tool_name,overwrite Name and replacement policy for the Seurat tools entry.
#' @param ... Additional parameters to other functions.
#' @return A Seurat object with the network in `tools`, an edge table (`TF`,
#'   `target`, `importance`), or a `GRNResult`. Multiple groups require `"seurat"` or `"result"`.
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
#' @seealso [GetGRN()], [GRNPlot()]
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
  if (any(startsWith("methods", prefixes) | startsWith("grn_method", prefixes))) {
    stop("Use 'methods' to select a GRN method.", call. = FALSE)
  }
  output_args <- setdiff(prefixes, "return_type")
  if (any(startsWith("return_type", output_args) | startsWith("return_data", output_args) |
    startsWith("return_seurat", output_args))) {
    stop("Use 'return_type' to select the GRN output.", call. = FALSE)
  }
  if ("method" %in% argument_names) {
    args <- list(...)
    method <- args[names(args) %in% "method"]
    args$method <- NULL
    if (!is.null(args$method_args) && !is.list(args$method_args)) {
      stop("method_args must be a named list.", call. = FALSE)
    }
    args$method_args <- c(method, args$method_args)
    return(do.call(RunGRN, c(list(object = object), args)))
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
  ...
) {
  RunGRN.default(
    object,
    regulators = regulators,
    targets = targets,
    genes_in = "rows",
    methods = methods,
    backend = if (missing(backend)) NULL else backend,
    assay = assay,
    layer = layer,
    ...
  )
}

#' @rdname RunGRN
#' @export
RunGRN.default <- function(
  object,
  regulators = NULL,
  targets = NULL,
  genes_in = c("rows", "columns"),
  methods = "grnboost2",
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
  return_type = NULL,
  input_mode = NULL,
  assay = NULL,
  layer = "counts",
  atac_assay = NULL,
  atac_layer = NULL,
  group.by = NULL,
  groups = NULL,
  pseudotime = NULL,
  resources = list(),
  method_args = list(),
  tool_name = NULL,
  overwrite = FALSE,
  ...
) {
  dots <- list(...)
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
  scientific <- list(
    max_edges_per_target = max_edges_per_target, n_rounds = n_rounds,
    learning_rate = learning_rate, max_depth = max_depth, max_features = max_features,
    subsample = subsample, early_stop_window_length = early_stop_window_length,
    exclude_self = exclude_self, correlation_threshold = correlation_threshold,
    lasso_degree = lasso_degree, lasso_alpha = lasso_alpha, max_lag = max_lag
  )
  assay <- if (is_seurat) assay %||% SeuratObject::DefaultAssay(object) else assay
  if (!is.list(method_args) || !is.list(resources) || (length(method_args) && (is.null(names(method_args)) || any(!nzchar(names(method_args))) || anyDuplicated(names(method_args))))) stop("method_args and resources must be named lists.", call. = FALSE)
  explicit <- intersect(names(match.call(expand.dots = FALSE)), names(scientific))
  arguments <- c(scientific[explicit], dots, method_args)
  if (anyDuplicated(names(arguments))) stop("Method arguments repeat an explicit top-level argument.", call. = FALSE)
  runtime <- arguments$runtime %||% list()
  arguments$runtime <- NULL
  if (length(runtime) && (!is.list(runtime) || is.null(names(runtime)) ||
    anyDuplicated(names(runtime)) || any(!names(runtime) %in% c("python", "timeout")))) {
    stop("method_args$runtime accepts only python and timeout for Python subprocesses.", call. = FALSE)
  }
  python_methods <- c("grnboost2", "gniplr", "scsgl", "regdiffusion", "linger", "scenicplus")
  if (length(runtime) && (!spec$method %in% python_methods || spec$backend == "cpp")) {
    stop("method_args$runtime is only used by Python subprocesses.", call. = FALSE)
  }
  unknown <- setdiff(names(arguments), spec$args)
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
  input_mode <- input_mode %||% if (length(spec$modes) == 1L) spec$modes else "rna"
  if (length(input_mode) != 1L || !input_mode %in% spec$modes) stop("Unsupported input_mode for ", spec$method, ".", call. = FALSE)
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
    } else {
      group.by
    }
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
      "resume_dir requires exactly one selected group.",
      call. = FALSE
    )
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
  } else if (spec$method %in% c("leap", "scode", "sincerities")) stop("This method requires named pseudotime values.", call. = FALSE)
  if (spec$method != "scenic" || length(regulators)) {
    regulators <- intersect(as.character(regulators %||% genes), genes)
  } else {
    regulators <- NULL
  }
  select_pando_targets <- spec$method == "pando" && is.null(targets)
  targets <- intersect(as.character(targets %||% genes), genes)
  if ((!is.null(regulators) && !length(regulators)) || !length(targets)) stop("No candidate regulators or targets occur in the input.", call. = FALSE)
  edge_limit <- arguments$max_edges_per_target
  if (!is.null(edge_limit) && (!is.numeric(edge_limit) || length(edge_limit) != 1L ||
    is.na(edge_limit) || edge_limit < 1 || (is.finite(edge_limit) && edge_limit != floor(edge_limit)))) {
    stop("max_edges_per_target must be a positive integer or Inf.", call. = FALSE)
  }
  if (length(cores) != 1L || !is.finite(cores) || cores < 1 || cores != as.integer(cores)) stop("cores must be a positive integer.", call. = FALSE)
  if (length(seed) != 1L || !is.finite(seed)) stop("seed must be a finite number.", call. = FALSE)
  config <- list(
    method = spec$method, backend = spec$backend,
    method_args = arguments, resources = resources, runtime = runtime, assay = assay, layer = layer,
    atac_assay = atac_assay, atac_layer = atac_layer,
    cores = as.integer(cores), seed = as.integer(seed), work_dir = work_dir, envname = envname,
    conda = conda, verbose = verbose
  )
  input <- list(
    object = if (is_seurat) object else NULL, expression = expression,
    cells = cells, groups = group_cells, group = NULL, pseudotime = pseudotime,
    regulators = regulators, targets = if (select_pando_targets) NULL else targets, input_mode = input_mode
  )
  config$external_sources <- if (spec$backend != "cpp") grn_external_sources(spec$method, resources)
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  attempt <- tempfile(paste0(prefix, "-"), tmpdir = work_dir)
  dir.create(attempt, recursive = TRUE, showWarnings = FALSE)
  attempt <- normalizePath(attempt, winslash = "/", mustWork = TRUE)
  run_config <- config
  run_config$work_dir <- attempt
  prepared <- if (spec$method == "multicsn") {
    grn_prepare_multicsn(input, run_config)
  } else {
    input
  }
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
    value <- spec$runner(group_input, group_config)
    value$chains <- value$chains %||% list(status = "not_supported")
    value$gene_sets <- value$gene_sets %||% list(status = "not_supported")
    if (spec$method == "linger") {
      if (anyDuplicated(value$cells) || !length(value$cells) || !all(value$cells %in% group_input$cells)) stop("The GRN backend returned invalid retained cells.", call. = FALSE)
      value$input_cells <- group_input$cells
    } else if (!identical(as.character(value$cells), as.character(group_input$cells))) stop("The GRN backend returned cells from a different group.", call. = FALSE)
    value
  })
  names(fitted) <- names(group_cells)
  artifacts <- unname(do.call(c, lapply(fitted, function(value) value$execution$artifacts %||% list())))
  result <- structure(list(
    schema_version = "1.0", method = spec$method,
    implementation = unique(vapply(fitted, `[[`, character(1), "implementation")),
    package_versions = lapply(fitted, `[[`, "package_versions"),
    parameters = list(method_args = arguments, backend = spec$backend, seed = config$seed, cores = config$cores),
    input_summary = list(cells = cells, features = genes, groups = group_cells),
    resources = c(
      resources[setdiff(names(resources), "genome")],
      if (is.character(resources[["genome"]])) list(genome = resources[["genome"]])
    ),
    execution = list(status = "completed", artifacts = artifacts, work_dir = attempt),
    groups = fitted
  ), class = c("GRNResult", "list"))
  if (return_type == "result") {
    return(result)
  }
  if (return_type == "seurat") {
    object@tools[[tool_name]] <- result
    return(object)
  }
  native <- fitted[[1L]]$native
  if (identical(spec$runner, grn_run_direct)) {
    adjacency <- native
  } else {
    edges <- GetGRN(result)
    adjacency <- data.frame(TF = edges$source, target = edges$target, importance = edges$raw_score, stringsAsFactors = FALSE)
  }
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
  native <- do.call(RunSCENIC, c(list(
    object = object, assay = assay, layer = layer,
    regulators = input$regulators, targets = input$targets,
    ranking_dbs = config$resources$ranking_dbs, motif_annotations = config$resources$motif_annotation,
    work_dir = config$work_dir, backend = config$backend, cores = config$cores,
    seed = config$seed, return_seurat = FALSE, envname = config$envname,
    conda = config$conda, verbose = config$verbose
  ), config$method_args))
  list(
    cells = input$cells, native = native,
    layer_specs = list(tf_gene = list(
      status = "computed", table_path = "adjacency", source = "TF", target = "target",
      score = "importance", source_type = "TF", target_type = "gene", score_type = "feature_importance", directed = TRUE,
      sign_column = NULL, sign_basis = "unknown", rank_by = "score"
    )),
    gene_sets = list(status = "computed", table_path = "regulon_list"),
    parameters = native$parameters, implementation = paste0("SCENIC/", config$backend),
    package_versions = c(scop = as.character(utils::packageVersion("scop"))), execution = list(status = "completed", artifacts = list())
  )
}

grn_run_direct <- function(input, config) {
  args <- config$method_args
  method <- config$method
  common <- list(
    object = input$expression, targets = input$targets, genes_in = "rows",
    verbose = config$verbose
  )
  if (method == "gniplr") {
    limit <- args$max_edges_per_target %||% Inf
    args$max_edges_per_target <- NULL
    native <- do.call(RunGNIPLR, c(common, list(backend = "cpp"), args))
    network <- attr(native, "grn_matrix")
    native <- cap_grn_edges_per_target(
      native[native$TF %in% input$regulators, , drop = FALSE], limit
    )
    attr(native, "grn_matrix") <- network
    parameters <- utils::modifyList(list(
      correlation_threshold = 0.3, lasso_degree = 30, lasso_alpha = 0.1,
      max_lag = 3, max_edges_per_target = limit, candidate_filter = "postfit_view"
    ), args)
  } else {
    common$regulators <- input$regulators
    common$cores <- config$cores
    if (method == "grnboost2") {
      fit <- RunGRNBoost2
      common$backend <- "cpp"
      common$seed <- config$seed
    } else {
      fit <- RunGENIE3
    }
    native <- do.call(fit, c(common, args))
    defaults <- if (method == "grnboost2") {
      formals(RunGRNBoost2.default)[c("max_edges_per_target", "n_rounds", "learning_rate", "max_depth",
        "max_features", "subsample", "early_stop_window_length", "exclude_self", "importance_norm_power",
        "correlation_fill", "boost_weight", "correlation_weight", "covariance_weight", "correlation_method")]
    } else {
      c(formals(RunGENIE3.default)["max_edges_per_target"], formals(GENIE3::GENIE3)[c("treeMethod", "K", "nTrees")])
    }
    parameters <- utils::modifyList(lapply(defaults, eval, envir = baseenv()), args)
  }
  parameters$cores <- config$cores
  parameters$seed <- config$seed
  spec <- list(
    status = "computed", table_path = character(), source = "TF", target = "target",
    score = "importance", source_type = "TF", target_type = "gene", score_type = "importance", directed = TRUE,
    sign_column = NULL, sign_basis = "unknown", rank_by = "score"
  )
  list(
    cells = input$cells, native = native, layer_specs = list(tf_gene = spec),
    parameters = parameters, implementation = paste0("scop::", method, "/", config$backend),
    package_versions = c(scop = as.character(utils::packageVersion("scop"))), execution = list(status = "completed", artifacts = list())
  )
}

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

grn_method_spec <- function(method, backend = NULL) {
  method <- tolower(gsub("[^a-zA-Z0-9+]", "", method))
  if (identical(method, "scenic+")) {
    method <- "scenicplus"
  }
  methods <- c(
    "grnboost2", "regdiffusion", "genie3", "gniplr", "ppcor", "pidc",
    "scorpion", "scsgl", "leap", "scode", "sincerities", "multicsn",
    "pando", "linger", "scenicplus", "scenic"
  )
  if (!method %in% methods) {
    stop("Unknown GRN method. Choose one of: ", paste(methods, collapse = ", "), call. = FALSE)
  }
  boost_cpp_args <- c("exclude_self", "importance_norm_power", "correlation_fill", "boost_weight",
    "correlation_weight", "covariance_weight", "correlation_method")
  boost_args <- c("max_edges_per_target", "n_rounds", "learning_rate", "max_depth",
    "max_features", "subsample", "early_stop_window_length")
  spec <- switch(method,
    grnboost2 = list(backends = c("cpp", "python", "native", "official"), runner = grn_run_rna, args = boost_args),
    regdiffusion = list(backends = c("python", "native"), runner = grn_run_rna, args = c(
      "T", "start_noise", "end_noise", "time_dim", "celltype_dim", "hidden_dims", "init_coef", "lr_nn", "lr_adj",
      "weight_decay_nn", "weight_decay_adj", "sparse_loss_coef", "adj_dropout", "batch_size", "n_steps", "train_split",
      "train_split_seed", "device", "compile", "eval_on_n_steps", "gradient_accumulation", "use_amp",
      "memory_efficient", "top_gene_percentile", "workers"
    )),
    genie3 = list(backends = c("cpp", "r", "native"), runner = grn_run_direct, args = c("nTrees", "treeMethod", "K", "max_edges_per_target")),
    gniplr = list(backends = c("cpp", "python", "native", "official"), runner = grn_run_rna, args = c("lasso_degree", "lasso_alpha", "max_edges_per_target")),
    ppcor = list(backends = c("r", "native"), runner = grn_run_rna, args = "method"),
    pidc = list(backends = c("r", "native"), runner = grn_run_rna, args = "log_scale"),
    scorpion = list(backends = c("r", "native"), runner = grn_run_rna, args = c("n_pc", "gamma_value", "show_progress", "tf_motifs", "ppi_network", "pad_small_input")),
    scsgl = list(backends = c("python", "native"), runner = grn_run_rna, args = c("pos_density", "neg_density")),
    leap = list(backends = c("r", "native"), runner = grn_run_time, args = c("max_lag_prop", "MAC_cutoff")),
    scode = list(backends = c("r", "native"), runner = grn_run_time, args = c("trajectory", "latent_dim", "n_iter", "n_runs")),
    sincerities = list(backends = c("r", "native"), runner = grn_run_time, args = c("trajectory", "max_time_points", "method", "noDIAG", "SIGN", "normalize")),
    multicsn = list(backends = "r", runner = infer_multicsn, modes = c("rna", "multiome"), args = unique(unlist(grn_multicsn_arguments(), use.names = FALSE))),
    pando = list(backends = c("r", "native"), runner = infer_pando, modes = "multiome", args = names(grn_pando_defaults())),
    linger = list(backends = c("python", "native"), runner = infer_linger, modes = "multiome", args = c("sample.by", "label.by", "min_genes", "min_cells", "genome", "activation", "patch_missing_keys", "resume_dir")),
    scenicplus = list(backends = c("python", "official", "native"), runner = grn_run_scenicplus, modes = "multiome", args = c("gsea_n_perm", "min_target_genes", "data_preparation", "motif_enrichment", "inference", "resume_dir")),
    scenic = list(backends = c("python", "cpp"), runner = grn_run_scenic, args = setdiff(names(formals(RunSCENIC)), c("object", "srt", "assay", "layer", "regulators", "targets", "work_dir", "backend", "cores", "seed", "return_seurat", "envname", "conda", "verbose", "ranking_dbs", "motif_annotations")))
  )
  backend <- backend %||% spec$backends[[1L]]
  if (length(backend) != 1L || !backend %in% spec$backends) {
    stop("Unsupported backend for ", method, ". Choose: ", paste(spec$backends, collapse = ", "), call. = FALSE)
  }
  if (backend == "cpp" && method %in% c("grnboost2", "gniplr")) {
    spec$runner <- grn_run_direct
    spec$args <- if (method == "grnboost2") c(boost_args, boost_cpp_args) else
      c(spec$args, "correlation_threshold", "max_lag")
  }
  spec$method <- method
  spec$backend <- backend
  spec$modes <- spec$modes %||% "rna"
  spec
}

grn_external_sources <- function(method, resources) {
  files <- switch(method,
    gniplr = c(entry = "code/methods/GNIPLR.py"),
    pidc = c(entry = "code/vendor/PIDC/PIDC.R"),
    scode = c(entry = "code/vendor/SCODE/SCODE.R"),
    scsgl = c(
      entry = "code/vendor/SCSGL/run_scSGL.py",
      paste0("code/vendor/SCSGL/pysrc/", c(
        "__init__.py", "associations/__init__.py", "associations/correlation.py",
        "graphlearning/__init__.py", "graphlearning/signed.py", "graphlearning/unsigned.py",
        "utils/__init__.py", "utils/rowsum_matrix.py"
      ))
    ),
    NULL
  )
  if (is.null(files)) {
    return(NULL)
  }
  root <- resources$source_dir
  if (!is.character(root) || length(root) != 1L || is.na(root) || !dir.exists(root)) {
    stop(method, " requires resources$source_dir pointing to a multiCSN_manuscript checkout. See ?RunGRN.", call. = FALSE)
  }
  paths <- file.path(root, files)
  missing <- !file.exists(paths) | dir.exists(paths)
  if (any(missing)) {
    stop("Missing ", method, " upstream source files: ", paste(files[missing], collapse = ", "), call. = FALSE)
  }
  stats::setNames(normalizePath(paths, winslash = "/", mustWork = TRUE), names(files))
}

grn_rna_matrix_edges <- function(matrix) {
  index <- which(!is.na(matrix), arr.ind = TRUE)
  data.frame(
    regulator = rownames(matrix)[index[, 1L]],
    target = colnames(matrix)[index[, 2L]],
    weight = as.numeric(matrix[index]), stringsAsFactors = FALSE
  )
}

grn_run_rna <- function(input, config) {
  method <- config$method
  arguments <- config$method_args
  expression <- input$expression
  genes <- rownames(expression)
  regulators <- input$regulators
  targets <- input$targets
  keep <- union(regulators, targets)
  expression <- as.matrix(expression[genes %in% keep, , drop = FALSE])
  if (nrow(expression) < 2L) stop("RNA inference needs at least two candidate genes.", call. = FALSE)
  execution <- list(status = "completed", artifacts = list(), fit_genes = rownames(expression), candidate_policy = "fit_union_then_filter")
  if (method %in% c("ppcor", "pidc", "scorpion")) {
    fit <- switch(method,
      ppcor = infer_ppcor(expression, arguments),
      pidc = infer_pidc(expression, regulators, targets, arguments, config),
      scorpion = infer_scorpion(expression, arguments, config$cores)
    )
    edges <- fit$edges
    extra <- fit$extra %||% list()
    arguments <- fit$parameters
    implementation <- fit$implementation
    versions <- fit$package_versions
    edges <- edges[edges$regulator != edges$target & edges$regulator %in% regulators & edges$target %in% targets, , drop = FALSE]
    edges <- edges[order(if (method == "ppcor") -abs(edges$weight) else -edges$weight, method = "radix"), , drop = FALSE]
    rownames(edges) <- NULL
  } else {
    limit <- arguments$max_edges_per_target
    arguments$max_edges_per_target <- NULL
    python <- config$runtime$python %||% conda_python(
      envname = config$envname %||% env_module_environment(method),
      conda = config$conda %||% "auto"
    )
    dir.create(config$work_dir, recursive = TRUE, showWarnings = FALSE)
    expression_path <- file.path(config$work_dir, "expression.tsv")
    output_path <- file.path(config$work_dir, "edges.tsv")
    receipt_path <- file.path(config$work_dir, "native-receipt.json")
    request_path <- file.path(config$work_dir, "request.json")
    utils::write.table(data.frame(gene = rownames(expression), expression, check.names = FALSE), expression_path,
      sep = "\t", row.names = FALSE, quote = TRUE
    )
    request <- list(
      method = method, expression = normalizePath(expression_path), output = output_path,
      receipt = receipt_path, work_dir = config$work_dir,
      regulators = I(regulators), targets = I(targets), method_args = arguments, seed = config$seed, cores = config$cores, verbose = config$verbose
    )
    if (method %in% c("gniplr", "scsgl")) {
      request$upstream_source <- config$external_sources[["entry"]]
    }
    jsonlite::write_json(request, request_path, auto_unbox = TRUE, null = "null", pretty = TRUE)
    script <- system.file("grn", paste0(method, ".py"), package = "scop", mustWork = TRUE)
    process <- runner_run_process(python, c(script, normalizePath(request_path)),
      work_dir = config$work_dir, timeout = config$runtime$timeout %||% Inf,
      env = c(OMP_NUM_THREADS = as.character(config$cores), OPENBLAS_NUM_THREADS = "1", PYTHONHASHSEED = as.character(config$seed), PYTHONDONTWRITEBYTECODE = "1")
    )
    edges <- utils::read.delim(output_path,
      stringsAsFactors = FALSE, check.names = FALSE,
      colClasses = c(regulator = "character", target = "character"), na.strings = character()
    )
    receipt <- jsonlite::read_json(receipt_path, simplifyVector = TRUE)
    extra <- list()
    implementation <- receipt$implementation
    versions <- unlist(receipt$package_versions)
    arguments <- receipt$parameters
    if (!is.null(limit)) {
      arguments$max_edges_per_target <- limit
      if (is.finite(limit)) {
        names(edges) <- c("TF", "target", "importance")
        edges <- cap_grn_edges_per_target(edges, limit)
        names(edges) <- c("regulator", "target", "weight")
      }
    }
    execution <- list(
      status = "completed", artifacts = list(receipt = receipt_path, edges = output_path), process = process,
      fit_genes = rownames(expression), candidate_policy = "fit_union_then_filter"
    )
  }
  if (any(!is.finite(edges$weight))) stop("Native RNA inference returned non-finite edge scores.", call. = FALSE)
  signed <- method %in% c("ppcor", "scsgl")
  directed <- !method %in% c("ppcor", "pidc", "scsgl")
  specification <- list(
    status = "computed", table_path = "edges", source = "regulator", target = "target", score = "weight",
    source_type = if (directed) "TF" else "gene", target_type = "gene", directed = directed,
    score_type = switch(method,
      ppcor = "partial_correlation",
      pidc = "context_information",
      scorpion = "scorpion_score",
      scsgl = "signed_graph_weight",
      gniplr = "negative_log10_p",
      grnboost2 = "feature_importance",
      regdiffusion = "diffusion_importance"
    ),
    sign_column = if (signed) "weight" else NULL, sign_basis = if (signed) "association" else "unknown",
    rank_by = if (signed) "absolute_score" else "score"
  )
  list(
    cells = input$cells, native = c(list(edges = edges), extra), layer_specs = list(tf_gene = specification),
    chains = list(status = "not_supported"), gene_sets = list(status = "not_computed"),
    parameters = arguments, implementation = implementation, package_versions = versions, execution = execution
  )
}

grn_time_data <- function(input, config) {
  expression <- as.matrix(input$expression)
  time <- input$pseudotime
  genes <- rownames(expression)
  regulators <- input$regulators
  targets <- input$targets
  if (!is.null(regulators)) genes <- intersect(genes, regulators)
  if (!is.null(targets)) genes <- union(genes, intersect(rownames(expression), targets))
  expression <- expression[genes, , drop = FALSE]
  if (nrow(expression) < 2L || ncol(expression) < 3L) {
    stop("Time methods require at least two candidate genes and three cells.", call. = FALSE)
  }
  if (config$method != "leap") {
    trajectory <- config$method_args$trajectory
    if (is.null(trajectory) && ncol(time) != 1L) {
      stop("Select one trajectory with method_args$trajectory for SCODE or SINCERITIES.", call. = FALSE)
    }
    if (!is.null(trajectory)) {
      if (!is.character(trajectory) || length(trajectory) != 1L ||
        is.null(colnames(time)) || !trajectory %in% colnames(time)) {
        stop("method_args$trajectory must name a pseudotime column.", call. = FALSE)
      }
      time <- time[, trajectory, drop = FALSE]
    }
    if (any(!is.finite(time))) {
      stop("The selected trajectory must provide a finite time for every selected cell.", call. = FALSE)
    }
    if (length(unique(as.numeric(time))) < 2L) {
      stop("The selected trajectory needs at least two distinct times.", call. = FALSE)
    }
    order <- order(time[, 1L], colnames(expression), method = "radix")
    expression <- expression[, order, drop = FALSE]
    time <- time[order, , drop = FALSE]
  }
  list(
    expression = expression, pseudotime = time,
    regulators = regulators, targets = targets
  )
}

grn_run_time <- function(input, config) {
  data <- grn_time_data(input, config)
  fit <- switch(config$method,
    leap = infer_leap,
    scode = infer_scode,
    sincerities = infer_sincerities
  )
  inputs <- list(expression = data$expression, pseudotime = data$pseudotime)
  inputs <- c(inputs, switch(config$method,
    scode = list(source = config$external_sources[["entry"]], work_dir = config$work_dir, seed = config$seed),
    sincerities = list(cores = config$cores, verbose = config$verbose)
  ))
  arguments <- config$method_args
  arguments$trajectory <- NULL
  rng_kind <- RNGkind()
  had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_rng) {
    rng_state <- get(".Random.seed", envir = .GlobalEnv)
  }
  on.exit(
    {
      do.call(RNGkind, as.list(rng_kind))
      if (had_rng) {
        assign(".Random.seed", rng_state, envir = .GlobalEnv)
      } else {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  set.seed(config$seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  native <- do.call(fit, c(inputs, arguments))
  edges <- native$edges
  if (!is.null(data$regulators)) {
    edges <- edges[edges$regulator %in% data$regulators, , drop = FALSE]
  }
  if (!is.null(data$targets)) {
    edges <- edges[edges$target %in% data$targets, , drop = FALSE]
  }
  edges <- edges[order(-abs(edges$weight), edges$regulator, edges$target, method = "radix"), , drop = FALSE]
  rownames(edges) <- NULL
  native$edges <- edges
  sign_column <- if (config$method == "sincerities" && native$parameters$SIGN == 0L) NULL else "signed_weight"
  score_type <- switch(config$method,
    leap = "lagged_correlation",
    scode = "dynamic_coefficient",
    sincerities = "normalized_distribution_regression_weight"
  )
  if (config$method == "sincerities" && !isTRUE(native$parameters$normalize)) score_type <- "distribution_regression_weight"
  sign_basis <- switch(config$method,
    leap = "lagged_correlation",
    scode = "dynamic_coefficient",
    sincerities = if (is.null(sign_column)) "unknown" else "partial_spearman_correlation"
  )
  list(
    cells = input$cells,
    native = native,
    layer_specs = list(tf_gene = list(
      status = "computed", table_path = "edges", source = "regulator", target = "target",
      score = "weight", source_type = if (is.null(input$regulators)) "gene" else "TF", target_type = "gene", score_type = score_type,
      directed = TRUE, sign_column = sign_column, sign_basis = sign_basis, rank_by = "absolute_score"
    )),
    chains = list(status = "not_supported"),
    gene_sets = list(status = "not_supported"),
    parameters = native$parameters,
    implementation = native$implementation,
    package_versions = native$package_versions,
    execution = list(status = "completed", artifacts = list())
  )
}
