metric_macro_f1 <- function(predicted, truth) {
  classification_metrics_compute(
    predicted = predicted,
    truth = truth
  )[["macro_f1"]]
}
