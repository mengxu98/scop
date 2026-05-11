
#' @title Propeller differential abundance wrapper
#'
#' @description
#' Method-specific implementation used by [RunProportionTest] when
#' `proportion_method = "propeller"`.
#' This implementation works on sample-level proportions using a propeller-style
#' transformed test and stores standardized outputs for plotting.
#'
#' @md
#' @inheritParams RunProportionTest
#' @param n_bootstrap Number of bootstrap iterations for confidence intervals.
#'
#' @return A method result bundle used internally by [RunProportionTest].
#'
#' @export
RunPropeller <- function(
  srt,
  group.by,
  split.by,
  sample.by,
  comparison = NULL,
  n_bootstrap = 1000,
  seed = 11,
  verbose = TRUE
) {
  meta_data <- validate_proportion_inputs(
    srt = srt,
    group.by = group.by,
    split.by = split.by,
    sample.by = sample.by,
    require_sample = TRUE
  )

  comparisons_condition <- parse_proportion_comparisons(
    meta_data = meta_data,
    split.by = split.by,
    comparison = comparison,
    include_bidirectional = TRUE
  )

  engine <- "internal"

  results_list <- list()
  for (i in seq_len(nrow(comparisons_condition))) {
    cluster_1 <- comparisons_condition[i, 1]
    cluster_2 <- comparisons_condition[i, 2]
    comparison_name <- paste0(cluster_1, "_vs_", cluster_2)

    prop_res <- sample_level_proportion_test(
      meta_data = meta_data,
      group.by = group.by,
      split.by = split.by,
      sample.by = sample.by,
      cluster_1 = cluster_1,
      cluster_2 = cluster_2,
      n_bootstrap = n_bootstrap,
      transform = "logit",
      seed = seed + i,
      verbose = verbose
    )

    results_list[[comparison_name]] <- standardize_proportion_result(
      prop_res,
      cluster_1 = cluster_1,
      cluster_2 = cluster_2,
      comparison_name = comparison_name,
      method = "propeller"
    )
  }

  list(
    method = "propeller",
    results = results_list,
    result_levels = "group",
    details = list(engine = engine),
    parameters = list(
      sample.by = sample.by,
      n_bootstrap = n_bootstrap,
      transform = "logit",
      engine = engine
    )
  )
}
