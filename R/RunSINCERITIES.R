infer_sincerities <- function(
  expression,
  pseudotime,
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

  check_r(c("glmnet", "ppcor"), strict = TRUE, verbose = FALSE)
  fit <- sincerities_fit(expression, times, ordered_time, method, noDIAG, SIGN, cores, verbose)
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
    implementation = "scop::RunSINCERITIES",
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

sincerities_fit <- function(expression, times, time_points, method, noDIAG, SIGN, cores, verbose) {
  groups <- lapply(time_points, function(point) expression[, times == point, drop = FALSE])
  genes <- seq_len(nrow(expression))
  intervals <- seq_len(length(time_points) - 1L)
  map <- function(fun) {
    if (cores == 1L) {
      return(lapply(genes, fun))
    }
    thisutils::parallelize_fun(
      x = genes, fun = fun, cores = cores, throw_error = TRUE,
      clean_result = FALSE, verbose = verbose
    )
  }
  distances <- map(function(gene) {
    vapply(intervals, function(i) {
      stats::ks.test(groups[[i]][gene, ], groups[[i + 1L]][gene, ])$statistic
    }, numeric(1))
  })
  distances <- do.call(cbind, distances) / diff(time_points)
  predictors <- distances[-length(intervals), , drop = FALSE]
  responses <- distances[-1L, , drop = FALSE]
  variable <- apply(predictors, 2L, stats::sd) > 0
  folds <- seq_len(nrow(predictors))
  alphas <- switch(as.integer(method), 0, seq(0, 1, 0.1), 1)

  results <- map(function(gene) {
    y <- responses[, gene]
    beta <- numeric(length(genes))
    keep <- variable
    if (noDIAG == 1L) {
      keep[gene] <- FALSE
    }
    columns <- which(keep)
    if (!length(columns) || max(tabulate(match(y, unique(y)))) >= length(y) - 1L) {
      return(list(lambda = NA_real_, alpha = NA_real_, beta = beta))
    }
    if (length(columns) == 1L) {
      stop(
        "SINCERITIES needs at least two variable predictors per fitted target; only one remains for ",
        rownames(expression)[gene], ".",
        call. = FALSE
      )
    }
    x <- predictors[, columns, drop = FALSE]
    fits <- lapply(alphas, function(alpha) {
      glmnet::cv.glmnet(
        x, y, alpha = alpha, nfolds = length(folds), foldid = folds, keep = TRUE,
        lower.limits = 0, upper.limits = Inf, grouped = FALSE
      )
    })
    errors <- vapply(fits, function(fit) fit$cvm[match(fit$lambda.min, fit$lambda)], numeric(1))
    best <- which(is.finite(errors))
    if (!length(best)) {
      stop("SINCERITIES cross-validation produced no finite error for ", rownames(expression)[gene], ".", call. = FALSE)
    }
    best <- max(best[errors[best] == min(errors[best])])
    cv <- fits[[best]]
    coefficients <- stats::coef(cv, s = "lambda.min")
    if (length(coefficients@i) > 1L) {
      beta[columns[coefficients@i[-1L]]] <- coefficients@x[-1L]
    }
    list(lambda = cv$lambda.min, alpha = alphas[best], beta = beta)
  })

  adjacency <- do.call(cbind, lapply(results, `[[`, "beta"))
  if (SIGN == 1L) {
    adjacency <- adjacency * sign(ppcor::pcor(t(expression), method = "spearman")$estimate)
    dimnames(adjacency) <- rep(list(rownames(expression)), 2L)
  }
  list(
    DISTANCE_matrix = unname(distances), adj_matrix = adjacency,
    lambda = unname(vapply(results, `[[`, numeric(1), "lambda")),
    alpha = unname(vapply(results, `[[`, numeric(1), "alpha"))
  )
}
