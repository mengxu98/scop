metric_weighted_recall <- function(predicted, truth, rare_threshold = 0.05) {
  classification_metrics_compute(
    predicted = predicted,
    truth = truth,
    rare_threshold = rare_threshold
  )[["rare_recall"]]
}
