metric_ari <- function(predicted, truth) {
  classification_metrics_compute(
    predicted = predicted,
    truth = truth
  )[["ari"]]
}

choose2 <- function(x) {
  x * (x - 1) / 2
}
