metric_purity <- function(predicted, truth) {
  predicted <- as.character(predicted)
  truth <- as.character(truth)
  keep <- !is.na(predicted) & !is.na(truth)
  predicted <- predicted[keep]
  truth <- truth[keep]
  if (length(predicted) == 0) {
    return(NA_real_)
  }
  tab <- table(predicted, truth)
  sum(apply(tab, 1, max)) / sum(tab)
}
