sincerities_map <- function(x, fun, cores = 1L, verbose = FALSE) {
  if (cores <= 1L) return(lapply(x, fun))
  unname(thisutils::parallelize_fun(x = x, fun = fun, cores = cores, throw_error = TRUE, clean_result = FALSE, verbose = verbose))
}

sincerities_fit <- function(
  DATA,
  method = 1,
  noDIAG = 1,
  SIGN = 0,
  cores = 1L,
  verbose = FALSE
) {
  if (!method %in% c(1, 2, 3, 4)) {
    stop("Choose regularization regression strategy with 1,2,3,4.")
  }
  if (!noDIAG %in% c(0, 1)) {
    stop("noDIAG should be either 0 or 1")
  }
  if (!SIGN %in% c(0, 1)) {
    stop("SIGN should be either 0 or 1")
  }

  single_cell_data <- DATA$singleCELLdata
  time <- DATA$time
  numGENES <- DATA$numGENES
  num_time_points <- length(time)

  if (num_time_points < 5) {
    stop("** DATA with a number of time points < 5. Please run SINCERITIES_CROSS_VALIDATION function **")
  }

  distance_cols <- sincerities_map(
    x = seq_len(numGENES),
    cores = cores,
    verbose = verbose,
    fun = function(gi) {
      values <- numeric(num_time_points - 1)
      for (ti in seq_len(num_time_points - 1)) {
        data_ti <- t(single_cell_data[[ti]])
        data_ti_plus1 <- t(single_cell_data[[ti + 1]])
        p1 <- data_ti[gi, ]
        p2 <- data_ti_plus1[gi, ]
        values[ti] <- stats::ks.test(p1, p2)$statistic
      }
      values
    }
  )

  DISTANCE_matrix <- do.call(cbind, distance_cols)
  deltaT <- replicate(ncol(DISTANCE_matrix), time[2:length(time)] - time[1:(length(time) - 1)])
  DISTANCE_matrix <- DISTANCE_matrix / deltaT

  if (method == 1) {
    alphas <- 0
  } else if (method == 2) {
    alphas <- seq(0, 1, 0.1)
  } else if (method == 3) {
    alphas <- 1
  } else {
    stop("Manual alpha selection is not supported in benchmark mode.")
  }

  X_matrix <- DISTANCE_matrix[1:(num_time_points - 2), , drop = FALSE]
  predictor_keep_global <- apply(X_matrix, 2, function(values) {
    values <- suppressWarnings(as.numeric(values))
    keep <- is.finite(values)
    if (sum(keep) < 2) {
      return(FALSE)
    }
    stats::sd(values[keep]) > 0
  })
  nfold <- nrow(X_matrix)
  foldid <- seq_len(nfold)

  gene_results <- sincerities_map(
    x = seq_len(numGENES),
    cores = cores,
    verbose = verbose,
    fun = function(gi) {

      Y_vector <- DISTANCE_matrix[2:(num_time_points - 1), gi]
      full_beta <- numeric(ncol(X_matrix))
      y_keep <- is.finite(Y_vector)
      if (sum(y_keep) < 2 || stats::sd(Y_vector[y_keep]) == 0) {
        return(list(
          lambda = NA_real_,
          alpha = NA_real_,
          beta = full_beta
        ))
      }

      predictor_keep <- predictor_keep_global
      if (noDIAG == 1 && gi <= length(predictor_keep)) {
        predictor_keep[gi] <- FALSE
      }
      valid_cols <- which(predictor_keep)
      if (!length(valid_cols)) {
        return(list(
          lambda = NA_real_,
          alpha = NA_real_,
          beta = full_beta
        ))
      }

      X_fit <- X_matrix[, valid_cols, drop = FALSE]
      lambda <- numeric(length(alphas))
      cvERROR <- numeric(length(alphas))
      beta <- matrix(data = 0, nrow = ncol(X_matrix), ncol = length(alphas))

      for (test in seq_along(alphas)) {
        CV_results <- tryCatch(
          glmnet::cv.glmnet(
            X_fit,
            Y_vector,
            alpha = alphas[test],
            nfolds = nfold,
            foldid = foldid,
            keep = TRUE,
            lower.limits = 0,
            upper.limits = Inf,
            grouped = FALSE
          ),
          error = function(...) NULL
        )
        if (is.null(CV_results)) {
          lambda[test] <- NA_real_
          cvERROR[test] <- Inf
          next
        }
        lambda[test] <- CV_results$lambda.min
        cvERROR[test] <- CV_results$cvm[CV_results$lambda == CV_results$lambda.min][1]
        coef.CV_results <- stats::coef(CV_results, s = "lambda.min")
        if (length(coef.CV_results@i) > 1) {
          selected_idx <- coef.CV_results@i[-1]
          beta[valid_cols[selected_idx], test] <- coef.CV_results@x[-1]
        }
      }

      if (!any(is.finite(cvERROR))) {
        return(list(
          lambda = NA_real_,
          alpha = NA_real_,
          beta = full_beta
        ))
      }

      minIdx <- max(which(cvERROR == min(cvERROR)))
      list(
        lambda = lambda[minIdx],
        alpha = alphas[minIdx],
        beta = beta[, minIdx]
      )
    }
  )

  pred_lambda_min <- do.call(cbind, lapply(gene_results, `[[`, "beta"))

  if (SIGN == 1) {
    parcorr_matrix <- ppcor::pcor(DATA$totDATA, method = "spearman")$estimate
    pred_lambda_min <- pred_lambda_min * sign(parcorr_matrix)
  }

  list(
    DISTANCE_matrix = DISTANCE_matrix,
    adj_matrix = pred_lambda_min,
    lambda = vapply(gene_results, `[[`, numeric(1), "lambda"),
    alpha = vapply(gene_results, `[[`, numeric(1), "alpha")
  )
}
