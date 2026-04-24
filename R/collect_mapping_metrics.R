collect_mapping_metrics <- function(
  srt,
  predicted_col,
  truth_col,
  probability_col = NULL,
  rare_threshold = 0.05
) {
  predicted <- srt[[predicted_col, drop = TRUE]]
  truth <- srt[[truth_col, drop = TRUE]]
  keep <- !is.na(predicted) & !is.na(truth)
  predicted <- predicted[keep]
  truth <- truth[keep]
  probability <- NULL
  if (!is.null(probability_col) && probability_col %in% colnames(srt@meta.data)) {
    probability <- srt[[probability_col, drop = TRUE]][keep]
  }
  per_class <- metric_class_table(predicted = predicted, truth = truth)
  summary <- data.frame(
    metric = c(
      "accuracy",
      "macro_f1",
      "purity",
      "nmi",
      "ari",
      "mean_confidence",
      "rare_recall"
    ),
    value = c(
      metric_accuracy(predicted = predicted, truth = truth),
      metric_macro_f1(predicted = predicted, truth = truth),
      metric_purity(predicted = predicted, truth = truth),
      metric_nmi(predicted = predicted, truth = truth),
      metric_ari(predicted = predicted, truth = truth),
      if (is.null(probability)) NA_real_ else mean(probability, na.rm = TRUE),
      metric_weighted_recall(
        predicted = predicted,
        truth = truth,
        rare_threshold = rare_threshold
      )
    ),
    stringsAsFactors = FALSE
  )
  list(
    summary = summary,
    per_class = per_class,
    confusion = as.data.frame.matrix(table(predicted, truth))
  )
}

majority_map <- function(predicted, truth) {
  predicted <- as.character(predicted)
  truth <- as.character(truth)
  keep <- !is.na(predicted) & !is.na(truth)
  predicted <- predicted[keep]
  truth <- truth[keep]
  if (length(predicted) == 0) {
    return(character(0))
  }
  tab <- table(predicted, truth)
  apply(tab, 1, function(x) colnames(tab)[which.max(x)])
}

apply_majority_map <- function(predicted, truth) {
  mapping <- majority_map(predicted = predicted, truth = truth)
  mapped <- unname(mapping[as.character(predicted)])
  mapped[is.na(mapped)] <- "unclassified"
  factor(mapped)
}
