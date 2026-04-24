metric_macro_f1 <- function(predicted, truth) {
  class_df <- metric_class_table(predicted = predicted, truth = truth)
  mean(class_df$f1, na.rm = TRUE)
}
