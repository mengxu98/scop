grn_rna_method_args <- function(method) {
  switch(method,
    ppcor = c("method"),
    pidc = c("log_scale"),
    scorpion = c("n_pc", "gamma_value", "show_progress", "tf_motifs", "ppi_network", "pad_small_input"),
    scsgl = c("pos_density", "neg_density"),
    gniplr = c("lasso_degree", "lasso_alpha"),
    grnboost2 = c("early_stop_window_length"),
    regdiffusion = c("T", "start_noise", "end_noise", "time_dim", "celltype_dim", "hidden_dims", "init_coef", "lr_nn", "lr_adj", "weight_decay_nn", "weight_decay_adj", "sparse_loss_coef", "adj_dropout", "batch_size", "n_steps", "train_split", "train_split_seed", "device", "compile", "eval_on_n_steps", "gradient_accumulation", "use_amp", "memory_efficient", "top_gene_percentile", "workers"),
    stop("Unsupported RNA method: ", method, call. = FALSE)
  )
}

grn_rna_script <- function(...) {
  path <- system.file("grn", "rna", ..., package = "scop")
  if (!nzchar(path) || !file.exists(path)) {
    stop("The installed scop package is missing its RNA runner files.", call. = FALSE)
  }
  path
}

grn_rna_matrix_edges <- function(matrix, source_names = rownames(matrix), target_names = colnames(matrix)) {
  index <- which(!is.na(matrix), arr.ind = TRUE)
  data.frame(
    regulator = as.character(source_names[index[, 1L]]),
    target = as.character(target_names[index[, 2L]]),
    weight = as.numeric(matrix[index]), stringsAsFactors = FALSE
  )
}

grn_run_rna <- function(input, config) {
  method <- config$method
  arguments <- config$method_args
  unknown <- setdiff(names(arguments), grn_rna_method_args(method))
  if (length(unknown)) stop("Unknown ", method, " argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  expression <- input$expression
  genes <- rownames(expression)
  if (is.null(genes) || is.null(colnames(expression)) || anyDuplicated(genes) || anyDuplicated(colnames(expression))) {
    stop("RNA inference requires unique gene and cell identifiers.", call. = FALSE)
  }
  regulators <- intersect(as.character(input$regulators), genes)
  targets <- intersect(as.character(input$targets), genes)
  keep <- union(regulators, targets)
  expression <- as.matrix(expression[genes %in% keep, , drop = FALSE])
  if (!is.numeric(expression) || any(!is.finite(expression))) {
    stop("RNA expression must contain finite numeric values.", call. = FALSE)
  }
  if (nrow(expression) < 2L || !length(regulators) || !length(targets)) {
    stop("RNA inference needs at least two candidate genes and nonempty regulators and targets.", call. = FALSE)
  }
  if (method %in% c("gniplr", "pidc", "scsgl")) {
    config$external_sources <- config$external_sources %||% grn_external_sources(method, config$resources)
  }
  dir.create(config$work_dir, recursive = TRUE, showWarnings = FALSE)
  if (method %in% c("ppcor", "pidc", "scorpion") && length(config$runtime)) {
    runtime <- config$runtime
    rscript <- runtime$rscript %||% runtime$executable
    libraries <- runtime$r_library %||% .libPaths()
    alternate <- (!is.null(rscript) && !identical(normalizePath(rscript), normalizePath(file.path(R.home("bin"), "Rscript")))) ||
      !identical(normalizePath(libraries), normalizePath(.libPaths()))
    if (alternate) {
      if (is.null(runtime$executable)) runtime <- grn_resolve_runtime(config, method)
      request_path <- file.path(config$work_dir, "r-request.rds")
      output_path <- file.path(config$work_dir, "r-result.rds")
      payload <- list(input = input[c("expression", "regulators", "targets", "cells")], config = config, libraries = runtime$r_library)
      payload$input$expression <- expression
      payload$config$runtime <- list()
      saveRDS(payload, request_path)
      config$runtime <- runtime
      process <- grn_run_process(runtime$rscript, c("--vanilla", grn_rna_script("run.R"), normalizePath(request_path), normalizePath(output_path, mustWork = FALSE)), config)
      result <- readRDS(output_path)
      result$execution$process <- process
      return(result)
    }
  }
  extra <- list()
  execution <- list(status = "completed", artifacts = list(), fit_genes = rownames(expression), candidate_policy = "fit_union_then_filter")
  versions <- character()
  if (method %in% c("ppcor", "pidc", "scorpion")) {
    package <- switch(method, ppcor = "ppcor", pidc = "pbapply", scorpion = "SCORPION")
    if (!requireNamespace(package, quietly = TRUE)) stop("Install the ", package, " package before running ", method, ".", call. = FALSE)
    versions <- stats::setNames(as.character(utils::packageVersion(package)), package)
    if (identical(method, "ppcor")) {
      correlation <- arguments$method %||% "spearman"
      if (!correlation %in% c("spearman", "pearson", "kendall")) stop("Invalid ppcor method.", call. = FALSE)
      result <- ppcor::pcor(t(expression), method = correlation)
      edges <- grn_rna_matrix_edges(result$estimate)
      index <- cbind(match(edges$regulator, rownames(result$estimate)), match(edges$target, colnames(result$estimate)))
      edges$p_value <- as.numeric(result$p.value[index])
      edges$statistic <- as.numeric(result$statistic[index])
      extra <- list(sample_size = result$n, variables = result$gp, method = result$method)
      arguments <- list(method = correlation)
      implementation <- "ppcor::pcor"
    } else if (identical(method, "pidc")) {
      if (!requireNamespace("purrr", quietly = TRUE)) stop("Install purrr before running PIDC.", call. = FALSE)
      nonconstant <- apply(expression, 1L, stats::sd) > 0
      expression <- expression[nonconstant, , drop = FALSE]
      regulators <- intersect(regulators, rownames(expression))
      targets <- intersect(targets, rownames(expression))
      if (length(regulators) < 2L || !length(targets)) stop("PIDC needs at least two nonconstant regulators and one target.", call. = FALSE)
      environment <- new.env(parent = baseenv())
      environment$is <- methods::is
      sys.source(config$external_sources[["entry"]], envir = environment)
      arguments <- list(log_scale = arguments$log_scale %||% FALSE)
      if (!is.logical(arguments$log_scale) || length(arguments$log_scale) != 1L || is.na(arguments$log_scale)) stop("PIDC log_scale must be TRUE or FALSE.", call. = FALSE)
      result <- environment$PIDC(expression, regulators = regulators, targets = targets,
        logScale = arguments$log_scale, ncores = config$cores, diag = "zero", verbose = config$verbose)
      edges <- grn_rna_matrix_edges(result)
      implementation <- "multiCSN-manuscript::PIDC-R"
      versions <- c(versions, purrr = as.character(utils::packageVersion("purrr")))
    } else {
      original_genes <- rownames(expression)
      pad_small <- arguments$pad_small_input %||% TRUE
      if (nrow(expression) < 10L && isTRUE(pad_small)) {
        scale <- stats::sd(as.numeric(expression))
        if (!is.finite(scale) || scale <= 0) scale <- 1
        signal <- colMeans(expression)
        padding <- t(vapply(seq_len(10L - nrow(expression)), function(i) signal + i * seq_along(signal) * scale * 1e-9, numeric(ncol(expression))))
        rownames(padding) <- paste0(".SCORPION_PAD_", seq_len(nrow(padding)))
        colnames(padding) <- colnames(expression)
        expression <- rbind(expression, padding)
      }
      arguments <- utils::modifyList(list(n_pc = min(25L, nrow(expression) - 1L, ncol(expression) - 1L), gamma_value = 10, show_progress = FALSE, tf_motifs = NULL, ppi_network = NULL, pad_small_input = pad_small), arguments)
      if (!is.numeric(arguments$n_pc) || length(arguments$n_pc) != 1L || !is.finite(arguments$n_pc) || arguments$n_pc < 1L || arguments$n_pc != as.integer(arguments$n_pc) || arguments$n_pc >= min(dim(expression))) stop("SCORPION n_pc must be a positive integer smaller than both expression dimensions.", call. = FALSE)
      if (!is.numeric(arguments$gamma_value) || length(arguments$gamma_value) != 1L || !is.finite(arguments$gamma_value) || arguments$gamma_value <= 0) stop("SCORPION gamma_value must be positive and finite.", call. = FALSE)
      result <- SCORPION::scorpion(tfMotifs = arguments$tf_motifs, gexMatrix = expression,
        ppiNet = arguments$ppi_network, nCores = as.integer(config$cores),
        nPC = as.integer(arguments$n_pc), gammaValue = arguments$gamma_value,
        showProgress = arguments$show_progress)
      scores <- if (is.list(result) && !is.null(result$regNet)) result$regNet else result
      edges <- grn_rna_matrix_edges(as.matrix(scores))
      if (is.list(result)) extra <- result[setdiff(names(result), "regNet")]
      extra$input_preprocessing <- list(original_genes = original_genes, padded_genes = setdiff(rownames(expression), original_genes))
      implementation <- "SCORPION::scorpion"
    }
    edges <- edges[edges$regulator != edges$target & edges$regulator %in% regulators & edges$target %in% targets, , drop = FALSE]
    edges <- edges[order(if (method == "ppcor") -abs(edges$weight) else -edges$weight, method = "radix"), , drop = FALSE]
    rownames(edges) <- NULL
  } else {
    if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite before running a Python RNA backend.", call. = FALSE)
    if (is.null(config$runtime$executable)) {
      config$runtime <- grn_resolve_runtime(config, method)
    }
    expression_path <- file.path(config$work_dir, "expression.tsv")
    output_path <- file.path(config$work_dir, "edges.tsv")
    receipt_path <- file.path(config$work_dir, "native-receipt.json")
    request_path <- file.path(config$work_dir, "request.json")
    utils::write.table(data.frame(gene = rownames(expression), expression, check.names = FALSE), expression_path,
      sep = "\t", row.names = FALSE, quote = TRUE)
    request <- list(method = method, expression = normalizePath(expression_path), output = normalizePath(output_path, mustWork = FALSE),
      receipt = normalizePath(receipt_path, mustWork = FALSE), work_dir = normalizePath(config$work_dir),
      regulators = I(regulators), targets = I(targets), method_args = arguments, seed = config$seed, cores = config$cores, verbose = config$verbose)
    if (method %in% c("gniplr", "scsgl")) {
      request$upstream_source <- config$external_sources[["entry"]]
    }
    jsonlite::write_json(request, request_path, auto_unbox = TRUE, null = "null", pretty = TRUE)
    python <- config$runtime$python %||% config$runtime$executable
    process <- grn_run_process(python, c(grn_rna_script("run.py"), normalizePath(request_path)), config,
      env = c(OMP_NUM_THREADS = as.character(config$cores), OPENBLAS_NUM_THREADS = "1", PYTHONHASHSEED = as.character(config$seed), PYTHONDONTWRITEBYTECODE = "1"))
    edges <- utils::read.delim(output_path, stringsAsFactors = FALSE, check.names = FALSE,
      colClasses = c(regulator = "character", target = "character"), na.strings = character())
    receipt <- jsonlite::read_json(receipt_path, simplifyVector = TRUE)
    implementation <- receipt$implementation
    versions <- unlist(receipt$package_versions)
    arguments <- receipt$parameters
    execution <- list(status = "completed", artifacts = list(receipt = receipt_path, edges = output_path), process = process,
      fit_genes = rownames(expression), candidate_policy = "fit_union_then_filter")
  }
  if (any(!is.finite(edges$weight))) stop("Native RNA inference returned non-finite edge scores.", call. = FALSE)
  signed <- method %in% c("ppcor", "scsgl")
  directed <- !method %in% c("ppcor", "pidc", "scsgl")
  specification <- list(status = "computed", table_path = "edges", source = "regulator", target = "target", score = "weight",
    source_type = if (directed) "TF" else "gene", target_type = "gene", directed = directed,
    score_type = switch(method, ppcor = "partial_correlation", pidc = "context_information", scorpion = "scorpion_score", scsgl = "signed_graph_weight", gniplr = "negative_log10_p", grnboost2 = "feature_importance", regdiffusion = "diffusion_importance"),
    sign_column = if (signed) "weight" else NULL, sign_basis = if (signed) "association" else "unknown",
    rank_by = if (signed) "absolute_score" else "score")
  list(cells = input$cells, native = c(list(edges = edges), extra), layer_specs = list(tf_gene = specification),
    chains = list(status = "not_supported"), gene_sets = list(status = "not_computed"),
    parameters = arguments, implementation = implementation, package_versions = versions, execution = execution)
}
