run_sincerities <- function(
  expression,
  pseudotime,
  source,
  cores,
  verbose,
  max_time_points = NA_integer_,
  method = 1L,
  noDIAG = 1L,
  SIGN = 0L,
  normalize = TRUE
) {
  flags <- list(method = method, noDIAG = noDIAG, SIGN = SIGN)
  for (field in names(flags)) {
    value <- flags[[field]]
    allowed <- if (field == "method") 1:3 else 0:1
    if (length(value) != 1L || !is.finite(value) || !value %in% allowed) {
      stop("Invalid SINCERITIES ", field, ".", call. = FALSE)
    }
  }
  if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize)) {
    stop("SINCERITIES normalize must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(max_time_points) != 1L || (!is.na(max_time_points) &&
    (!is.finite(max_time_points) || max_time_points < 5L || max_time_points != as.integer(max_time_points)))) {
    stop("max_time_points must be NA or an integer of at least five.", call. = FALSE)
  }

  times <- as.numeric(pseudotime[, 1L])
  if (!is.na(max_time_points) && length(unique(times)) > max_time_points) {
    ranks <- rank(times, ties.method = "first")
    times <- floor((ranks - 1) * max_time_points / length(ranks))
  }
  ordered_time <- sort(unique(times))
  if (length(ordered_time) < 5L) {
    stop("SINCERITIES requires at least five distinct time points or bins.", call. = FALSE)
  }
  expression <- expression[apply(expression, 1L, stats::sd) > 0, , drop = FALSE]
  if (nrow(expression) < 2L) {
    stop("SINCERITIES needs at least two nonconstant genes.", call. = FALSE)
  }

  upstream <- new.env(parent = globalenv())
  sys.source(source, envir = upstream)
  data <- list(
    time = ordered_time,
    singleCELLdata = lapply(ordered_time, function(point) t(expression[, times == point, drop = FALSE])),
    numGENES = nrow(expression),
    totDATA = t(expression)
  )
  fit <- upstream$run_SINCERITIES_parallel(
    data,
    method = method, noDIAG = noDIAG, SIGN = SIGN,
    cores = cores, verbose = verbose
  )
  fit <- lapply(fit, unname)
  if (SIGN == 1L) {
    dimnames(fit$adj_matrix) <- rep(list(rownames(expression)), 2L)
  }
  scores <- fit$adj_matrix
  dimnames(scores) <- list(rownames(expression), rownames(expression))
  divisor <- if (isTRUE(normalize)) max(abs(scores)) else 1
  if (!is.finite(divisor) || divisor == 0) {
    divisor <- 1
  }
  index <- which(row(scores) != col(scores), arr.ind = TRUE)
  edges <- data.frame(
    regulator = rownames(scores)[index[, 1L]],
    target = colnames(scores)[index[, 2L]],
    weight = abs(scores[index]) / divisor,
    signed_weight = scores[index] / divisor
  )

  list(
    fit = fit,
    binned_time = stats::setNames(times, colnames(expression)),
    implementation = "SINCERITIES::run_SINCERITIES_parallel:multiCSN_manuscript",
    package_versions = c(
      glmnet = as.character(utils::packageVersion("glmnet")),
      ppcor = as.character(utils::packageVersion("ppcor"))
    ),
    edges = edges,
    parameters = list(
      max_time_points = max_time_points, method = method, noDIAG = noDIAG,
      SIGN = SIGN, normalize = normalize, normalization_divisor = divisor,
      time_points = ordered_time
    ),
    genes = rownames(expression),
    cells = colnames(expression)
  )
}
