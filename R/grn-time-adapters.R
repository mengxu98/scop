grn_time_data <- function(input, config) {
  expression <- as.matrix(input$expression)
  time <- input$pseudotime
  if (!is.matrix(time) || !is.numeric(time) || is.null(rownames(time)) ||
      anyDuplicated(rownames(time)) || !setequal(rownames(time), colnames(expression))) {
    stop("Time methods require a numeric pseudotime matrix with unique cell IDs matching expression.", call. = FALSE)
  }
  time <- time[colnames(expression), , drop = FALSE]
  if (!ncol(time) || any(rowSums(is.finite(time)) == 0L) || any(is.infinite(time))) {
    stop("Every cell must have a finite pseudotime in at least one trajectory.", call. = FALSE)
  }
  if (anyDuplicated(rownames(expression)) || anyDuplicated(colnames(expression)) ||
      anyNA(expression) || any(!is.finite(expression))) {
    stop("Time methods require finite expression and unique gene and cell IDs.", call. = FALSE)
  }
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
  list(expression = expression, pseudotime = time,
       regulators = regulators, targets = targets)
}

grn_run_time <- function(input, config) {
  data <- grn_time_data(input, config)
  sources <- if (config$method %in% c("scode", "sincerities")) {
    config$external_sources %||% grn_external_sources(config$method, config$resources)
  }
  runtime <- if (is.null(config$runtime$kind)) grn_resolve_runtime(config, "r") else config$runtime
  config$runtime <- runtime
  dir.create(config$work_dir, recursive = TRUE, showWarnings = FALSE)
  request <- file.path(config$work_dir, "time-input.rds")
  output <- file.path(config$work_dir, "time-output.rds")
  script <- system.file("grn", "time", "run.R", package = "scop")
  if (!nzchar(script)) stop("The installed scop time-method runner is missing.", call. = FALSE)
  saveRDS(list(data = data, method = config$method, sources = sources,
               method_args = config$method_args, seed = config$seed,
               cores = config$cores, verbose = config$verbose), request)
  command <- runtime$rscript
  if (is.null(command)) command <- runtime$executable
  environment <- character()
  if (length(runtime$r_library)) environment <- c(R_LIBS_USER = paste(runtime$r_library, collapse = .Platform$path.sep))
  execution <- grn_run_process(command, c("--vanilla", script, request, output), config, env = environment)
  if (!file.exists(output)) stop("The time-method runner did not produce its result.", call. = FALSE)
  native <- readRDS(output)
  sign_column <- if (config$method == "sincerities" && native$parameters$SIGN == 0L) NULL else "signed_weight"
  score_type <- switch(config$method, leap = "lagged_correlation", scode = "dynamic_coefficient",
                       sincerities = "normalized_distribution_regression_weight")
  if (config$method == "sincerities" && !isTRUE(native$parameters$normalize)) score_type <- "distribution_regression_weight"
  sign_basis <- switch(config$method, leap = "lagged_correlation", scode = "dynamic_coefficient",
                       sincerities = if (is.null(sign_column)) "unknown" else "partial_spearman_correlation")
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
    execution = list(status = "completed", artifacts = list(input = request, output = output), process = execution)
  )
}
