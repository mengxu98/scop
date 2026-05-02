metric_accuracy <- function(predicted, truth) {
  classification_metrics_compute(
    predicted = predicted,
    truth = truth
  )[["accuracy"]]
}
