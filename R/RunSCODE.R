#' @rdname RunGRN-methods
#' @export
RunSCODE <- function(object, pseudotime, ...) {
  RunGRN(object, methods = "scode", pseudotime = pseudotime, ...)
}

infer_scode <- function(
  expression,
  pseudotime,
  source,
  work_dir,
  seed,
  latent_dim = min(4L, nrow(expression)),
  n_iter = 100L,
  n_runs = 5L
) {
  settings <- list(latent_dim = latent_dim, n_iter = n_iter, n_runs = n_runs)
  for (field in names(settings)) {
    value <- settings[[field]]
    if (length(value) != 1L || !is.finite(value) || value < 1 || value != as.integer(value)) {
      stop("SCODE ", field, " must be a positive integer.", call. = FALSE)
    }
  }
  if (latent_dim > min(nrow(expression), ncol(expression))) {
    stop("SCODE latent_dim exceeds the available dimensions.", call. = FALSE)
  }

  upstream <- new.env(parent = baseenv())
  upstream$library <- function(package) invisible(loadNamespace(as.character(substitute(package))))
  upstream$ginv <- MASS::ginv
  upstream$lm <- stats::lm
  upstream$runif <- stats::runif
  upstream$read.table <- utils::read.table
  upstream$write.table <- utils::write.table
  times <- as.numeric(pseudotime[, 1L])
  normalized <- (times - min(times)) / diff(range(times))
  settings$time_origin <- min(times)
  settings$time_range <- diff(range(times))
  settings$normalized_time <- stats::setNames(normalized, colnames(expression))
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  expression_file <- file.path(work_dir, "scode-expression.tsv")
  time_file <- file.path(work_dir, "scode-time.tsv")
  utils::write.table(
    expression, expression_file,
    sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )
  utils::write.table(
    data.frame(cell = seq_len(ncol(expression)) - 1L, time = normalized),
    time_file,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )

  runs <- lapply(seq_len(n_runs), function(run) {
    run_seed <- seed + run - 1L
    set.seed(run_seed)
    directory <- file.path(work_dir, paste0("scode-run-", run))
    arguments <- c(expression_file, time_file, directory, nrow(expression), latent_dim, ncol(expression), n_iter)
    upstream$commandArgs <- function(trailingOnly) arguments
    upstream$system <- function(command) dir.create(directory)
    sys.source(source, envir = upstream)
    coefficients <- as.matrix(utils::read.table(file.path(directory, "A.txt"), sep = "\t", header = FALSE))
    if (!identical(dim(coefficients), rep(nrow(expression), 2L)) || any(!is.finite(coefficients))) {
      stop("SCODE returned an invalid coefficient matrix.", call. = FALSE)
    }
    dimnames(coefficients) <- list(rownames(expression), rownames(expression))
    list(A = coefficients, RSS = scan(file.path(directory, "RSS.txt"), quiet = TRUE), seed = run_seed)
  })
  coefficients <- Reduce(`+`, lapply(runs, `[[`, "A")) / n_runs
  index <- which(row(coefficients) != col(coefficients), arr.ind = TRUE)
  edges <- data.frame(
    regulator = colnames(coefficients)[index[, 2L]],
    target = rownames(coefficients)[index[, 1L]],
    weight = coefficients[index],
    signed_weight = coefficients[index]
  )

  list(
    runs = runs,
    coefficient_matrix = coefficients,
    implementation = "SCODE::SCODE.R",
    package_versions = c(MASS = as.character(utils::packageVersion("MASS"))),
    edges = edges,
    parameters = settings,
    genes = rownames(expression),
    cells = colnames(expression)
  )
}
