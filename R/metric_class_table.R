metric_class_table <- function(predicted, truth) {
  predicted <- as.character(predicted)
  truth <- as.character(truth)
  keep <- !is.na(predicted) & !is.na(truth)
  predicted <- predicted[keep]
  truth <- truth[keep]
  classes <- sort(unique(c(predicted, truth)))
  res <- lapply(classes, function(cls) {
    tp <- sum(predicted == cls & truth == cls)
    fp <- sum(predicted == cls & truth != cls)
    fn <- sum(predicted != cls & truth == cls)
    precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
    recall <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
    f1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) {
      NA_real_
    } else {
      2 * precision * recall / (precision + recall)
    }
    support <- sum(truth == cls)
    data.frame(
      class = cls,
      precision = precision,
      recall = recall,
      f1 = f1,
      support = support,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, res)
}
