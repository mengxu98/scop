metric_accuracy <- function(predicted, truth) {
  predicted <- as.character(predicted)
  truth <- as.character(truth)
  keep <- !is.na(predicted) & !is.na(truth)
  if (!any(keep)) {
    return(NA_real_)
  }
  mean(predicted[keep] == truth[keep])
}
