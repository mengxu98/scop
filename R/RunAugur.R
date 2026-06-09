#' @title Prioritize perturbed cell types using Augur
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param celltype.by Metadata column or vector defining the cell types used by
#' Augur.
#' @param label.by Metadata column or vector defining the labels to predict,
#' such as condition, stimulation, sample type, or technology.
#' @param assay Assay used to extract the feature matrix. If `NULL`, the
#' default assay is used.
#' @param layer Assay layer used as the feature matrix.
#' @param backend Backend used to run Augur. `"cpp"` uses a parity-preserving
#' `scop` R/C++ backend that keeps Augur's feature selection, sampling, random
#' forest, and metric semantics without requiring the Augur package. `"r"`
#' calls the native `Augur::calculate_auc` implementation.
#' @param features Features used by Augur. If `NULL`, all features in `assay`
#' are used.
#' @param n_subsamples,subsample_size,folds,min_cells,var_quantile,feature_perc,select_var,augur_mode,classifier,rf_params,lr_params
#' Arguments passed to `Augur::calculate_auc`.
#' @param cores Number of cores used by Augur.
#' @param prefix Prefix for metadata columns written to `srt@meta.data`.
#' @param tool_name Name of the `srt@tools` entry used to store Augur results.
#' @param add_meta Whether to write `prefix_auc` and `prefix_rank` metadata
#' columns back to each cell by matching `celltype.by`.
#' @param ... Additional arguments passed to `Augur::calculate_auc`.
#'
#' @return A `Seurat` object with native Augur results stored in
#' `srt@tools[[tool_name]]`. When `add_meta = TRUE`, cell-level metadata columns
#' `prefix_auc` and `prefix_rank` are added for use with existing `scop`
#' plotting functions such as [FeatureDimPlot()].
#'
#' @export
#'
#' @references
#' Skinnider, M.A., Squair, J.W., Kathe, C., et al. (2021). Cell type
#' prioritization in single-cell data. \emph{Nature Biotechnology}, 39,
#' 30-34. \doi{10.1038/s41587-020-0605-1}
#'
#' @examples
#' data(panc8_sub)
#' panc8_sub <- subset(panc8_sub, subset = tech %in% c("celseq", "celseq2"))
#' panc8_sub <- standard_scop(panc8_sub, verbose = FALSE)
#' panc8_sub <- RunAugur(
#'   panc8_sub,
#'   celltype.by = "celltype",
#'   label.by = "tech",
#'   n_subsamples = 5,
#'   subsample_size = 20,
#'   min_cells = 20,
#'   cores = 1,
#'   verbose = FALSE,
#'   rf_params = list(
#'     trees = 20,
#'     mtry = 2,
#'     min_n = NULL,
#'     importance = "accuracy"
#'   )
#' )
#'
#' panc8_sub@tools$Augur$AUC
#' FeatureDimPlot(
#'   panc8_sub,
#'   features = "augur_auc",
#'   reduction = "StandardUMAP2D",
#'   bg_cutoff = -Inf
#' )
RunAugur <- function(
  srt,
  celltype.by,
  label.by,
  assay = NULL,
  layer = "counts",
  backend = c("cpp", "r"),
  features = NULL,
  n_subsamples = 50,
  subsample_size = 20,
  folds = 3,
  min_cells = NULL,
  var_quantile = 0.5,
  feature_perc = 0.5,
  cores = 1,
  select_var = TRUE,
  augur_mode = c("default", "velocity", "permute"),
  classifier = c("rf", "lr"),
  rf_params = list(
    trees = 100,
    mtry = 2,
    min_n = NULL,
    importance = "accuracy"
  ),
  lr_params = list(
    mixture = 1,
    penalty = "auto"
  ),
  prefix = "augur",
  tool_name = "Augur",
  add_meta = TRUE,
  verbose = TRUE,
  ...
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  backend <- match.arg(backend)
  augur_mode <- match.arg(augur_mode)
  classifier <- match.arg(classifier)
  if (identical(backend, "r")) {
    check_r("neurorestore/Augur", verbose = FALSE)
  }
  assay <- assay %||% SeuratObject::DefaultAssay(srt)

  if (
    !is.character(prefix) ||
      length(prefix) != 1L ||
      is.na(prefix) ||
      !nzchar(prefix)
  ) {
    log_message(
      "{.arg prefix} must be a non-empty string",
      message_type = "error"
    )
  }
  if (
    !is.character(tool_name) ||
      length(tool_name) != 1L ||
      is.na(tool_name) ||
      !nzchar(tool_name)
  ) {
    log_message(
      "{.arg tool_name} must be a non-empty string",
      message_type = "error"
    )
  }

  if (length(celltype.by) == 1L && celltype.by %in% colnames(srt@meta.data)) {
    cell_types <- srt[[celltype.by, drop = TRUE]]
    celltype_col <- celltype.by
  } else if (length(celltype.by) == ncol(srt)) {
    cell_types <- celltype.by
    celltype_col <- "cell_type"
  } else {
    log_message(
      "{.arg celltype.by} must be one metadata column or one value per cell",
      message_type = "error"
    )
  }

  if (length(label.by) == 1L && label.by %in% colnames(srt@meta.data)) {
    labels <- srt[[label.by, drop = TRUE]]
    label_col <- label.by
  } else if (length(label.by) == ncol(srt)) {
    labels <- label.by
    label_col <- "label"
  } else {
    log_message(
      "{.arg label.by} must be one metadata column or one value per cell",
      message_type = "error"
    )
  }

  cell_types <- as.character(cell_types)
  labels <- as.character(labels)
  keep_cells <- !is.na(cell_types) & !is.na(labels)
  if (!all(keep_cells)) {
    log_message(
      "Drop {.val {sum(!keep_cells)}} cells with missing {.arg celltype.by} or {.arg label.by}",
      verbose = verbose
    )
  }
  if (sum(keep_cells) == 0L) {
    log_message(
      "No cells remain after removing missing Augur labels",
      message_type = "error"
    )
  }
  if (length(unique(labels[keep_cells])) < 2L) {
    log_message(
      "{.arg label.by} must contain at least two non-missing labels",
      message_type = "error"
    )
  }
  if (length(unique(cell_types[keep_cells])) < 1L) {
    log_message(
      "{.arg celltype.by} must contain at least one non-missing cell type",
      message_type = "error"
    )
  }

  expr <- GetAssayData5(srt, assay = assay, layer = layer)
  features <- features %||% rownames(expr)
  features <- unique(features)
  features <- intersect(features, rownames(expr))
  if (length(features) == 0L) {
    log_message(
      "No requested {.arg features} are present in {.arg assay}",
      message_type = "error"
    )
  }
  expr <- expr[features, keep_cells, drop = FALSE]
  expr <- expr[Matrix::rowSums(expr != 0) > 0, , drop = FALSE]
  if (nrow(expr) == 0L) {
    log_message(
      "No non-zero features remain for Augur",
      message_type = "error"
    )
  }

  meta <- data.frame(
    cell_type = factor(cell_types[keep_cells]),
    label = factor(labels[keep_cells]),
    row.names = colnames(srt)[keep_cells]
  )

  log_message(
    "Run {.pkg Augur} on {.val {ncol(expr)}} cells, {.val {nrow(expr)}} features, and {.val {nlevels(meta$cell_type)}} cell types",
    verbose = verbose
  )

  if (identical(backend, "r")) {
    augur <- get_namespace_fun("Augur", "calculate_auc")(
      input = expr,
      meta = meta,
      label_col = "label",
      cell_type_col = "cell_type",
      n_subsamples = n_subsamples,
      subsample_size = subsample_size,
      folds = folds,
      min_cells = min_cells,
      var_quantile = var_quantile,
      feature_perc = feature_perc,
      n_threads = cores,
      show_progress = verbose,
      select_var = select_var,
      augur_mode = augur_mode,
      classifier = classifier,
      rf_params = rf_params,
      lr_params = lr_params,
      ...
    )
  } else {
    augur <- augur_cpp(
      expr = expr,
      meta = meta,
      label_col = "label",
      cell_type_col = "cell_type",
      n_subsamples = n_subsamples,
      subsample_size = subsample_size,
      folds = folds,
      min_cells = min_cells,
      var_quantile = var_quantile,
      feature_perc = feature_perc,
      cores = cores,
      verbose = verbose,
      select_var = select_var,
      augur_mode = augur_mode,
      classifier = classifier,
      rf_params = rf_params,
      lr_params = lr_params,
      ...
    )
  }

  srt@tools[[tool_name]] <- augur
  srt@tools[[tool_name]][["input"]] <- list(
    assay = assay,
    layer = layer,
    backend = backend,
    features = rownames(expr),
    cells = colnames(expr),
    celltype.by = celltype_col,
    label.by = label_col,
    prefix = prefix
  )

  if (isTRUE(add_meta)) {
    auc_table <- as.data.frame(augur$AUC)
    if (!all(c("cell_type", "auc") %in% colnames(auc_table))) {
      log_message(
        "{.pkg Augur} did not return expected {.field AUC} columns",
        message_type = "error"
      )
    }
    auc_values <- auc_table[["auc"]]
    names(auc_values) <- as.character(auc_table[["cell_type"]])
    rank_values <- rank(-auc_values, ties.method = "average")

    auc_meta <- rep(NA_real_, ncol(srt))
    rank_meta <- rep(NA_real_, ncol(srt))
    names(auc_meta) <- colnames(srt)
    names(rank_meta) <- colnames(srt)
    matched_auc <- auc_values[cell_types]
    matched_rank <- rank_values[cell_types]
    auc_meta[] <- unname(matched_auc)
    rank_meta[] <- unname(matched_rank)

    srt[[paste0(prefix, "_auc")]] <- auc_meta
    srt[[paste0(prefix, "_rank")]] <- rank_meta
  }

  return(srt)
}

augur_row_sds <- function(mat) {
  sds <- MatrixGenerics::rowSds(mat)
  sds[is.na(sds)] <- 0
  sds
}

augur_select_random <- function(mat, feature_perc = 0.5) {
  if (feature_perc < 1) {
    features <- rownames(mat)
    keep <- sample(features, floor(nrow(mat) * feature_perc))
    mat <- mat[keep, , drop = FALSE]
  }
  mat
}

augur_select_variance <- function(
  mat,
  var_quantile = 0.5,
  filter_negative_residuals = FALSE
) {
  sds <- augur_row_sds(mat)
  sds[is.na(sds)] <- 0
  mat <- mat[sds > 0, , drop = FALSE]
  if (nrow(mat) == 0L) {
    return(mat)
  }

  if (var_quantile < 1 || isTRUE(filter_negative_residuals)) {
    means <- Matrix::rowMeans(mat)
    sds <- sds[sds > 0]
    cvs <- means / sds
    lower <- stats::quantile(cvs, 0.01, na.rm = TRUE)
    upper <- stats::quantile(cvs, 0.99, na.rm = TRUE)
    keep <- dplyr::between(cvs, lower, upper)
    if (sum(keep) < 3L) {
      return(mat)
    }

    cv0 <- cvs[keep]
    mean0 <- means[keep]
    if (any(mean0 < 0)) {
      model <- stats::loess(cv0 ~ mean0)
    } else {
      fit1 <- stats::loess(cv0 ~ mean0)
      fit2 <- stats::loess(cv0 ~ log(mean0))
      check_r("lmtest", verbose = FALSE)
      cox <- get_namespace_fun("lmtest", "coxtest")(fit1, fit2)
      probs <- cox[["Pr(>|z|)"]]
      model <- if (probs[1] < probs[2]) fit1 else fit2
    }
    genes <- rownames(mat)[keep]
    residuals <- stats::setNames(model$residuals, genes)

    if (isTRUE(filter_negative_residuals)) {
      genes <- names(residuals)[residuals > 0]
    } else {
      genes <- names(residuals)[
        residuals > stats::quantile(residuals, var_quantile, na.rm = TRUE)
      ]
    }
    mat <- mat[genes, , drop = FALSE]
  }

  mat
}

augur_cpp <- function(
  expr,
  meta,
  label_col = "label",
  cell_type_col = "cell_type",
  n_subsamples = 50,
  subsample_size = 20,
  folds = 3,
  min_cells = NULL,
  var_quantile = 0.5,
  feature_perc = 0.5,
  cores = 1,
  verbose = TRUE,
  select_var = TRUE,
  augur_mode = c("default", "velocity", "permute"),
  classifier = c("rf", "lr"),
  rf_params = list(
    trees = 100,
    mtry = 2,
    min_n = NULL,
    importance = "accuracy"
  ),
  lr_params = list(
    mixture = 1,
    penalty = "auto"
  ),
  ...
) {
  unused <- list(...)
  if (length(unused) > 0L) {
    log_message(
      "{.arg ...} is not supported by {.arg backend = \"cpp\"}",
      message_type = "error"
    )
  }
  classifier <- match.arg(classifier)
  augur_mode <- match.arg(augur_mode)
  if (!identical(classifier, "rf")) {
    log_message(
      "{.arg backend = \"cpp\"} currently supports {.arg classifier = \"rf\"} only",
      message_type = "error"
    )
  }
  if (n_subsamples > 1 && subsample_size / folds < 2) {
    log_message(
      "{.arg subsample_size} / {.arg folds} must be greater than or equal to 2",
      message_type = "error"
    )
  }
  min_cells <- min_cells %||% subsample_size
  check_r(
    c(
      "rsample",
      "yardstick",
      "tibble",
      "purrr",
      "MatrixGenerics",
      "sparseMatrixStats",
      "randomForest",
      "lmtest"
    ),
    verbose = FALSE
  )

  meta <- droplevels(meta)
  labels <- meta[[label_col]]
  cell_types <- meta[[cell_type_col]]
  if (!is.numeric(labels)) {
    labels <- factor(as.character(labels))
    meta[[label_col]] <- labels
  }
  if (length(dim(expr)) != 2 || !all(dim(expr) > 0)) {
    log_message(
      "Expression matrix has at least one dimension of size zero",
      message_type = "error"
    )
  }
  if (nrow(meta) != ncol(expr)) {
    log_message(
      "Number of cells in metadata does not match number of cells in expression",
      message_type = "error"
    )
  }
  if (dplyr::n_distinct(labels) == 1L) {
    log_message(
      "Only one label provided: {.val {unique(labels)}}",
      message_type = "error"
    )
  }
  if (any(is.na(labels))) {
    log_message(
      "Labels contain {.val {sum(is.na(labels))}} missing values",
      message_type = "error"
    )
  }
  if (any(is.na(cell_types))) {
    log_message(
      "Cell types contain {.val {sum(is.na(cell_types))}} missing values",
      message_type = "error"
    )
  }
  if (is.numeric(labels)) {
    log_message(
      "{.arg backend = \"cpp\"} currently supports classification labels only",
      message_type = "error"
    )
  }
  mode <- "classification"
  multiclass <- dplyr::n_distinct(labels) > 2L
  if (augur_mode == "velocity") {
    feature_perc <- 1
    var_quantile <- 1
  } else if (augur_mode == "permute" && n_subsamples < 100) {
    n_subsamples <- 500
  }

  res <- thisutils::parallelize_fun(
    unique(cell_types),
    fun = function(cell_type) {
      y <- labels[cell_types == cell_type]
      if (min(table(y)) < min_cells) {
        warning(
          "skipping cell type ",
          cell_type,
          ": minimum number of cells (",
          min(table(y)),
          ") is less than ",
          min_cells
        )
        return(list())
      }
      X <- expr[, cell_types == cell_type, drop = FALSE]
      if (nrow(X) >= 1000 && isTRUE(select_var)) {
        X <- augur_select_variance(
          X,
          var_quantile,
          filter_negative_residuals = FALSE
        )
      }

      multi_metric <- yardstick::metric_set(
        yardstick::accuracy,
        yardstick::precision,
        yardstick::recall,
        yardstick::sens,
        yardstick::spec,
        yardstick::npv,
        yardstick::ppv,
        yardstick::roc_auc
      )
      prob_select <- 3
      estimator <- ifelse(multiclass, "macro", "binary")
      if (multiclass) {
        prob_select <- seq(3, 3 + dplyr::n_distinct(labels) - 1)
      }
      metric_fun <- function(x) {
        withCallingHandlers(
          multi_metric(
            x,
            truth = true,
            estimate = pred,
            prob_select,
            estimator = estimator
          ),
          warning = function(w) {
            if (grepl("external vector in selections", conditionMessage(w))) {
              invokeRestart("muffleWarning")
            }
          }
        )
      }
      n_iter <- ifelse(n_subsamples < 1, 1, n_subsamples)
      tmp_results <- vector("list", n_iter)
      tmp_importances <- vector("list", n_iter)
      for (subsample_idx in seq_len(n_iter)) {
        set.seed(subsample_idx)
        if (augur_mode == "permute") {
          y <- sample(y)
        }
        if (n_subsamples < 1) {
          if (nrow(X) >= 1000 && feature_perc < 1) {
            X0 <- augur_select_random(X, feature_perc)
          } else {
            X0 <- X
          }
          X0 <- Matrix::t(X0)
          X0 <- as.matrix(X0)
          X0 <- as.data.frame(X0)
          X0 <- tibble::repair_names(X0)
          X0[["label"]] <- y
        } else {
          subsample_idxs <- unlist(
            lapply(
              levels(y),
              function(level) sample(which(y == level), subsample_size)
            ),
            use.names = FALSE
          )
          y0 <- y[subsample_idxs]
          if (nrow(X) >= 1000 && feature_perc < 1) {
            X0 <- augur_select_random(X, feature_perc)
          } else {
            X0 <- X
          }
          if (inherits(X0, "dgCMatrix")) {
            X0 <- augur_subsample_cpp(X0, as.integer(subsample_idxs))
            X0 <- as.data.frame(X0)
          } else {
            X0 <- X0[, subsample_idxs, drop = FALSE]
            X0 <- Matrix::t(X0)
            keep <- MatrixGenerics::colVars(X0) > 0
            X0 <- X0[, keep, drop = FALSE]
            X0 <- as.matrix(X0)
            X0 <- as.data.frame(X0)
          }
          X0 <- tibble::repair_names(X0)
          X0[["label"]] <- y0
        }

        cv <- rsample::vfold_cv(X0, v = folds, strata = "label")
        folded <- dplyr::mutate(
          cv,
          fits = purrr::map(
            splits,
            function(split) {
              train <- rsample::analysis(split)
              target_indexes <- which(colnames(train) == "label")
              if (identical(mode, "classification")) {
                y_var <- as.factor(t(train[, target_indexes]))
              } else {
                y_var <- t(train[, target_indexes])
              }
              x_var <- train[, -target_indexes]
              original_seed <- .Random.seed
              set.seed(1)
              forest <- randomForest::randomForest(
                y = y_var,
                x = x_var,
                importance = TRUE,
                localImp = TRUE,
                ntree = rf_params$trees,
                mtry = rf_params$mtry,
                min_n = rf_params$min_n,
                type = mode
              )
              .Random.seed <<- original_seed
              forest
            }
          ),
          pred = purrr::map2(
            splits,
            fits,
            function(split, model) {
              test <- rsample::assessment(split)
              tbl <- tibble::tibble(
                true = test$label,
                pred = stats::predict(model, test),
                prob = stats::predict(model, test, type = "prob")
              )
              tbl <- cbind(tbl, tbl$prob)
              tbl[["prob"]] <- NULL
              prob_cols <- match(levels(test$label), colnames(tbl), nomatch = 0L)
              prob_cols <- prob_cols[prob_cols > 0L]
              colnames(tbl)[prob_cols] <- paste0(".pred_", colnames(tbl)[prob_cols])
              tbl
            }
          )
        )

        eval <- dplyr::mutate(
          folded,
          metrics = purrr::map(pred, metric_fun)
        )
        eval <- eval[["metrics"]]
        result <- purrr::map2_df(
          eval,
          row.names(folded),
          function(.x, .y) dplyr::mutate(.x, fold = .y)
        )
        names(result) <- gsub("\\.", "", names(result))
        result <- dplyr::mutate(
          result,
          cell_type = cell_type,
          subsample_idx = subsample_idx
        )
        result <- dplyr::select(
          result,
          cell_type,
          subsample_idx,
          fold,
          metric,
          estimator,
          estimate
        )

        importance <- purrr::map(
          dplyr::pull(folded, fits),
          function(model) {
            importance <- as.data.frame(model$importance)
            importance <- tibble::rownames_to_column(importance, "gene")
            impval_name <- if (identical(rf_params$importance, "accuracy")) {
              "MeanDecreaseAccuracy"
            } else {
              "MeanDecreaseGini"
            }
            dplyr::rename(
              importance,
              importance = dplyr::all_of(impval_name)
            )
          }
        )
        importance <- purrr::map2_df(
          importance,
          seq_along(importance),
          function(.x, .y) dplyr::mutate(.x, fold = .y)
        )
        importance <- dplyr::mutate(
          importance,
          cell_type = cell_type,
          subsample_idx = subsample_idx
        )
        importance <- dplyr::select(
          importance,
          cell_type,
          subsample_idx,
          fold,
          gene,
          importance
        )

        tmp_results[[subsample_idx]] <- result
        tmp_importances[[subsample_idx]] <- importance
      }
      list(
        results = dplyr::bind_rows(tmp_results),
        importances = dplyr::bind_rows(tmp_importances)
      )
    },
    cores = cores,
    verbose = verbose
  )

  valid <- vapply(
    res,
    function(x) is.list(x) && all(c("results", "importances") %in% names(x)),
    logical(1)
  )
  if (!any(valid)) {
    if (all(lengths(res) == 0)) {
      log_message(
        "No cell type had at least {.val {min_cells}} cells in all conditions",
        message_type = "error"
      )
    }
    log_message(
      "No Augur results were produced; check the cell-type task errors above",
      message_type = "error"
    )
  }
  res <- res[valid]
  feature_importances <- dplyr::bind_rows(purrr::map(res, "importances"))
  results <- dplyr::bind_rows(purrr::map(res, "results"))
  if (!all(c("metric", "estimate", "cell_type", "subsample_idx") %in% colnames(results))) {
    log_message(
      "Augur results did not include the expected metric columns",
      message_type = "error"
    )
  }
  AUCs <- results
  AUCs <- dplyr::filter(AUCs, metric == "roc_auc")
  AUCs <- dplyr::group_by(AUCs, cell_type, subsample_idx)
  AUCs <- dplyr::summarise(AUCs, estimate = mean(estimate), .groups = "drop")
  AUCs <- dplyr::group_by(AUCs, cell_type)
  AUCs <- dplyr::summarise(AUCs, auc = mean(estimate), .groups = "drop")
  AUCs <- dplyr::arrange(AUCs, dplyr::desc(auc))

  list(
    X = expr,
    y = labels,
    cell_types = cell_types,
    parameters = list(
      n_subsamples = n_subsamples,
      subsample_size = subsample_size,
      folds = folds,
      min_cells = min_cells,
      var_quantile = var_quantile,
      feature_perc = feature_perc,
      cores = cores,
      classifier = classifier,
      rf_params = rf_params
    ),
    results = results,
    feature_importance = feature_importances,
    AUC = AUCs
  )
}
