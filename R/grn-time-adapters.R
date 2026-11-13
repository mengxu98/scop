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
  saveRDS(list(data = data, method = config$method, sources = sources,
               method_args = config$method_args, seed = config$seed, work_dir = normalizePath(config$work_dir, mustWork = TRUE),
               cores = config$cores, verbose = config$verbose), request)
  command <- runtime$rscript
  if (is.null(command)) command <- runtime$executable
  environment <- character()
  if (length(runtime$r_library)) environment <- c(R_LIBS_USER = paste(runtime$r_library, collapse = .Platform$path.sep))
  execution <- grn_run_process(command, grn_r_worker_args(grn_fit_time, request, output), config, env = environment)
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

grn_fit_time <- function(request) {
  work_dir <- request$work_dir
  data <- request$data
  regulators <- data$regulators
  targets <- data$targets
  expression <- data$expression
  time <- data$pseudotime
  method_args <- request$method_args
  defaults <- switch(request$method,
    leap = list(max_lag_prop = 1/3, MAC_cutoff = 0),
    scode = list(latent_dim = min(4L, nrow(expression)), n_iter = 100L, n_runs = 5L),
    sincerities = list(max_time_points = NA_integer_, method = 1L, noDIAG = 1L, SIGN = 0L, normalize = TRUE)
  )
  if (is.null(defaults)) stop("Unknown time method.")
  unknown <- setdiff(names(method_args), c(names(defaults), "trajectory"))
  if (length(unknown)) stop("Unknown time-method argument(s): ", paste(unknown, collapse = ", "))
  settings <- utils::modifyList(defaults, method_args[intersect(names(method_args), names(defaults))])
  set.seed(request$seed)
  native <- list()
  if (request$method %in% c("scode", "sincerities")) {
    upstream <- new.env(parent = globalenv())
    sys.source(request$sources[["entry"]], envir = upstream)
  }
  if (request$method == "leap") {
    if (!requireNamespace("LEAP", quietly = TRUE)) stop("Install LEAP to run this method.")
    if (length(settings$max_lag_prop) != 1L || !is.finite(settings$max_lag_prop) ||
        settings$max_lag_prop <= 0 || settings$max_lag_prop >= 1 ||
        length(settings$MAC_cutoff) != 1L || !is.finite(settings$MAC_cutoff) || settings$MAC_cutoff < 0) {
      stop("LEAP requires 0 < max_lag_prop < 1 and a nonnegative MAC_cutoff.")
    }
    branches <- lapply(seq_len(ncol(time)), function(branch) {
      cells <- which(is.finite(time[, branch]))
      cells[order(time[cells, branch], colnames(expression)[cells], method = "radix")]
    })
    branches <- unique(branches[lengths(branches) >= 3L])
    if (!length(branches)) stop("LEAP needs a trajectory with at least three finite-time cells.")
    native$branches <- lapply(branches, function(cells) {
      LEAP::MAC_counter(expression[, cells, drop = FALSE], max_lag_prop = settings$max_lag_prop,
                        MAC_cutoff = settings$MAC_cutoff, lag_matrix = FALSE, symmetric = FALSE)
    })
    native$ordered_cells <- lapply(branches, function(cells) colnames(expression)[cells])
    tables <- lapply(native$branches, function(result) {
      result <- as.matrix(result)
      if (!nrow(result)) return(data.frame(regulator = character(), target = character(), weight = numeric(), lag = numeric(), signed_weight = numeric()))
      data.frame(regulator = rownames(expression)[result[, "Row gene index"]],
                 target = rownames(expression)[result[, "Column gene index"]],
                 weight = as.numeric(result[, "Correlation"]), lag = as.numeric(result[, "Lag"]),
                 signed_weight = as.numeric(result[, "Correlation"]))
    })
    edges <- do.call(rbind, tables)
    edges <- edges[is.finite(edges$weight), , drop = FALSE]
    edges <- edges[order(-abs(edges$weight), method = "radix"), , drop = FALSE]
    edges <- edges[!duplicated(edges[c("regulator", "target")]), , drop = FALSE]
    settings$branches <- native$ordered_cells
    native$implementation <- "LEAP::MAC_counter"
    native$package_versions <- c(LEAP = as.character(utils::packageVersion("LEAP")))
  } else if (request$method == "scode") {
    for (field in c("latent_dim", "n_iter", "n_runs")) {
      value <- settings[[field]]
      if (length(value) != 1L || !is.finite(value) || value < 1 || value != as.integer(value)) stop("SCODE ", field, " must be a positive integer.")
    }
    if (settings$latent_dim > min(nrow(expression), ncol(expression))) stop("SCODE latent_dim exceeds the available dimensions.")
    times <- as.numeric(time[, 1L])
    normalized <- (times - min(times)) / diff(range(times))
    settings$time_origin <- min(times)
    settings$time_range <- diff(range(times))
    settings$normalized_time <- stats::setNames(normalized, colnames(expression))
    expression_file <- file.path(work_dir, "scode-expression.tsv")
    time_file <- file.path(work_dir, "scode-time.tsv")
    utils::write.table(expression, expression_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    utils::write.table(data.frame(cell = seq_len(ncol(expression)) - 1L, time = normalized), time_file,
                       sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    native$runs <- lapply(seq_len(settings$n_runs), function(run) {
      set.seed(request$seed + run - 1L)
      directory <- file.path(work_dir, paste0("scode-run-", run))
      upstream$run_scode(expression_file, time_file, directory, nrow(expression), settings$latent_dim, ncol(expression), settings$n_iter)
      matrix <- as.matrix(utils::read.table(file.path(directory, "A.txt"), sep = "\t", header = FALSE))
      if (!identical(dim(matrix), rep(nrow(expression), 2L)) || any(!is.finite(matrix))) stop("SCODE returned an invalid coefficient matrix.")
      dimnames(matrix) <- list(rownames(expression), rownames(expression))
      list(A = matrix, RSS = scan(file.path(directory, "RSS.txt"), quiet = TRUE), seed = request$seed + run - 1L)
    })
    native$coefficient_matrix <- Reduce(`+`, lapply(native$runs, `[[`, "A")) / settings$n_runs
    index <- which(row(native$coefficient_matrix) != col(native$coefficient_matrix), arr.ind = TRUE)
    edges <- data.frame(regulator = colnames(native$coefficient_matrix)[index[, 2L]],
                        target = rownames(native$coefficient_matrix)[index[, 1L]],
                        weight = native$coefficient_matrix[index], signed_weight = native$coefficient_matrix[index])
    native$implementation <- "SCODE::run_scode:multiCSN_manuscript"
    native$package_versions <- c(MASS = as.character(utils::packageVersion("MASS")))
  } else {
    for (field in c("method", "noDIAG", "SIGN")) {
      value <- settings[[field]]
      allowed <- if (field == "method") 1:3 else 0:1
      if (length(value) != 1L || !is.finite(value) || !value %in% allowed) stop("Invalid SINCERITIES ", field, ".")
    }
    times <- as.numeric(time[, 1L])
    if (!is.logical(settings$normalize) || length(settings$normalize) != 1L || is.na(settings$normalize)) stop("SINCERITIES normalize must be TRUE or FALSE.")
    cap <- settings$max_time_points
    if (length(cap) != 1L || (!is.na(cap) && (!is.finite(cap) || cap < 5L || cap != as.integer(cap)))) stop("max_time_points must be NA or an integer of at least five.")
    if (!is.na(cap) && length(unique(times)) > cap) {
      ranks <- rank(times, ties.method = "first")
      times <- floor((ranks - 1) * cap / length(ranks))
    }
    ordered_time <- sort(unique(times))
    if (length(ordered_time) < 5L) stop("SINCERITIES requires at least five distinct time points or bins.")
    keep <- apply(expression, 1L, stats::sd) > 0
    expression <- expression[keep, , drop = FALSE]
    if (nrow(expression) < 2L) stop("SINCERITIES needs at least two nonconstant genes.")
    data <- list(time = ordered_time, singleCELLdata = lapply(ordered_time, function(timepoint) t(expression[, times == timepoint, drop = FALSE])),
                 numGENES = nrow(expression), totDATA = t(expression))
    native$fit <- upstream$run_SINCERITIES_parallel(data, method = settings$method, noDIAG = settings$noDIAG,
                                                   SIGN = settings$SIGN, cores = request$cores, verbose = request$verbose)
    native$fit <- lapply(native$fit, unname)
    if (settings$SIGN == 1L) dimnames(native$fit$adj_matrix) <- rep(list(rownames(expression)), 2L)
    scores <- native$fit$adj_matrix
    dimnames(scores) <- list(rownames(expression), rownames(expression))
    divisor <- if (isTRUE(settings$normalize)) max(abs(scores)) else 1
    if (!is.finite(divisor) || divisor == 0) divisor <- 1
    native$binned_time <- stats::setNames(times, colnames(expression))
    index <- which(row(scores) != col(scores), arr.ind = TRUE)
    edges <- data.frame(regulator = rownames(scores)[index[, 1L]], target = colnames(scores)[index[, 2L]],
                        weight = abs(scores[index]) / divisor, signed_weight = scores[index] / divisor)
    settings$normalization_divisor <- divisor
    settings$time_points <- ordered_time
    native$implementation <- "SINCERITIES::run_SINCERITIES_parallel:multiCSN_manuscript"
    native$package_versions <- c(glmnet = as.character(utils::packageVersion("glmnet")), ppcor = as.character(utils::packageVersion("ppcor")))
  }
  if (!is.null(regulators)) edges <- edges[edges$regulator %in% regulators, , drop = FALSE]
  if (!is.null(targets)) edges <- edges[edges$target %in% targets, , drop = FALSE]
  edges <- edges[order(-abs(edges$weight), edges$regulator, edges$target, method = "radix"), , drop = FALSE]
  rownames(edges) <- NULL
  native$edges <- edges
  native$parameters <- settings
  native$genes <- rownames(expression)
  native$cells <- colnames(expression)
  native$source_sha256 <- c(
    grn_fit_time = digest::digest(body(sys.function()), algo = "sha256"),
    stats::setNames(vapply(request$sources, function(path) digest::digest(file = path, algo = "sha256"), character(1)), basename(as.character(request$sources)))
  )
  native
}
