metric_weighted_recall <- function(predicted, truth, rare_threshold = 0.05) {
  class_df <- metric_class_table(predicted = predicted, truth = truth)
  support <- class_df$support / sum(class_df$support)
  rare_df <- class_df[support <= rare_threshold, , drop = FALSE]
  if (nrow(rare_df) == 0) {
    return(NA_real_)
  }
  mean(rare_df$recall, na.rm = TRUE)
}
