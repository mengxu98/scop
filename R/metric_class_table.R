metric_class_table <- function(predicted, truth) {
  classification_metrics_compute(
    predicted = predicted,
    truth = truth
  )[["class_table"]]
}
