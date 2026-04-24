metric_ari <- function(predicted, truth) {
  predicted <- as.character(predicted)
  truth <- as.character(truth)
  keep <- !is.na(predicted) & !is.na(truth)
  predicted <- predicted[keep]
  truth <- truth[keep]
  if (length(predicted) == 0) {
    return(NA_real_)
  }
  tab <- table(predicted, truth)
  a <- rowSums(tab)
  b <- colSums(tab)
  sum_nij <- sum(choose2(tab))
  sum_a <- sum(choose2(a))
  sum_b <- sum(choose2(b))
  n <- sum(tab)
  expected <- sum_a * sum_b / choose2(n)
  max_index <- (sum_a + sum_b) / 2
  if (max_index == expected) {
    return(NA_real_)
  }
  (sum_nij - expected) / (max_index - expected)
}

choose2 <- function(x) {
  x * (x - 1) / 2
}
