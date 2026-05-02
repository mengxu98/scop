metric_purity <- function(predicted, truth) {
  classification_metrics_compute(
    predicted = predicted,
    truth = truth
  )[["purity"]]
}
