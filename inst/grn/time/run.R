arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) stop("Expected an input RDS and an output RDS.")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1L]])
request <- readRDS(arguments[[1L]])
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
  expression_file <- file.path(dirname(arguments[[2L]]), "scode-expression.tsv")
  time_file <- file.path(dirname(arguments[[2L]]), "scode-time.tsv")
  utils::write.table(expression, expression_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  utils::write.table(data.frame(cell = seq_len(ncol(expression)) - 1L, time = normalized), time_file,
                     sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  native$runs <- lapply(seq_len(settings$n_runs), function(run) {
    set.seed(request$seed + run - 1L)
    directory <- file.path(dirname(arguments[[2L]]), paste0("scode-run-", run))
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
sources <- c(script, request$sources)
native$source_sha256 <- stats::setNames(vapply(sources, function(path) digest::digest(file = path, algo = "sha256"), character(1)), basename(sources))
saveRDS(native, arguments[[2L]])
