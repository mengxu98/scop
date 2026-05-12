classification_metrics_compute <- function(predicted, truth, rare_threshold = 0.05) {
  predicted <- as.character(predicted)
  truth <- as.character(truth)
  keep <- !is.na(predicted) & !is.na(truth)
  predicted <- predicted[keep]
  truth <- truth[keep]

  if (length(predicted) == 0L) {
    return(list(
      accuracy = NA_real_,
      macro_f1 = NA_real_,
      purity = NA_real_,
      nmi = NA_real_,
      ari = NA_real_,
      rare_recall = NA_real_,
      class_table = NULL
    ))
  }

  classes <- sort(unique(c(predicted, truth)))
  classification_metrics(
    predicted = predicted,
    truth = truth,
    classes = classes,
    rare_threshold = rare_threshold
  )
}
