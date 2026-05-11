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
  if (classification_metrics_available()) {
    classification_metrics(
      predicted = predicted,
      truth = truth,
      classes = classes,
      rare_threshold = rare_threshold
    )
  } else {
    class_table <- classification_metrics_class_table_r(
      predicted = predicted,
      truth = truth,
      classes = classes
    )
    tab <- table(predicted, truth)
    n <- sum(tab)
    pi <- rowSums(tab) / n
    pj <- colSums(tab) / n
    pij <- tab / n
    non_zero <- pij > 0
    mi <- sum(pij[non_zero] * log(pij[non_zero] / (pi[row(tab)][non_zero] * pj[col(tab)][non_zero])))
    hx <- metric_entropy(rowSums(tab))
    hy <- metric_entropy(colSums(tab))
    a <- rowSums(tab)
    b <- colSums(tab)
    sum_nij <- sum(choose2(tab))
    sum_a <- sum(choose2(a))
    sum_b <- sum(choose2(b))
    expected <- sum_a * sum_b / choose2(n)
    max_index <- (sum_a + sum_b) / 2
    support <- class_table$support / sum(class_table$support)
    rare_df <- class_table[support <= rare_threshold, , drop = FALSE]
    list(
      accuracy = mean(predicted == truth),
      macro_f1 = mean(class_table$f1, na.rm = TRUE),
      purity = sum(apply(tab, 1, max)) / n,
      nmi = if ((hx + hy) == 0) NA_real_ else 2 * mi / (hx + hy),
      ari = if (max_index == expected) NA_real_ else (sum_nij - expected) / (max_index - expected),
      rare_recall = if (nrow(rare_df) == 0) NA_real_ else mean(rare_df$recall, na.rm = TRUE),
      class_table = class_table
    )
  }
}

classification_metrics_class_table_r <- function(predicted, truth, classes = sort(unique(c(predicted, truth)))) {
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

classification_metrics_available <- function() {
  exists("classification_metrics", mode = "function") &&
    isTRUE(is.loaded("_scop_classification_metrics"))
}
