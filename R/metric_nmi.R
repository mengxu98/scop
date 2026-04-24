metric_nmi <- function(predicted, truth) {
  predicted <- as.character(predicted)
  truth <- as.character(truth)
  keep <- !is.na(predicted) & !is.na(truth)
  predicted <- predicted[keep]
  truth <- truth[keep]
  if (length(predicted) == 0) {
    return(NA_real_)
  }
  tab <- table(predicted, truth)
  n <- sum(tab)
  pi <- rowSums(tab) / n
  pj <- colSums(tab) / n
  pij <- tab / n
  non_zero <- pij > 0
  mi <- sum(pij[non_zero] * log(pij[non_zero] / (pi[row(tab)][non_zero] * pj[col(tab)][non_zero])))
  hx <- metric_entropy(rowSums(tab))
  hy <- metric_entropy(colSums(tab))
  if ((hx + hy) == 0) {
    return(NA_real_)
  }
  2 * mi / (hx + hy)
}
