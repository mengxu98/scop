metric_silhouette <- function(embeddings, labels, maximize = TRUE) {
  labels <- as.factor(labels)
  keep <- !is.na(labels)
  embeddings <- embeddings[keep, , drop = FALSE]
  labels <- droplevels(labels[keep])
  if (nrow(embeddings) < 3 || nlevels(labels) < 2) {
    return(NA_real_)
  }
  sil <- cluster::silhouette(
    x = as.integer(labels),
    dist = stats::dist(embeddings)
  )
  score <- mean(sil[, "sil_width"], na.rm = TRUE)
  if (isTRUE(maximize)) {
    return(score)
  }
  1 - abs(score)
}
