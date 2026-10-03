#' @title Propeller differential abundance wrapper
#'
#' @description
#' Method-specific implementation used by [RunProportionTest] when
#' `proportion_method = "propeller"`.
#' Uses the optional speckle backend: `getTransformedProps()` computes the
#' transformed proportions and `propeller.ttest()` performs limma empirical
#' Bayes moderated tests. Each requested comparison is a two-group test;
#' multiple conditions produce pairwise results, not an omnibus ANOVA.
#' P-values, FDR, and `Tstatistic` come from speckle. `obs_log2FD` and bootstrap
#' intervals describe the ratio of mean untransformed sample proportions,
#' with a `1e-5` pseudocount; these intervals are not moderated-test intervals.
#'
#' @md
#' @inheritParams RunProportionTest
#' @param sample.by Metadata column identifying biological samples. Required
#' when calling `RunPropeller()` directly. Fully paired sample IDs across two
#' conditions add a sample blocking effect to the speckle design matrix;
#' partially paired comparisons are rejected. At least two biological samples
#' per condition are required.
#' @param n_bootstrap Number of sample bootstrap iterations for descriptive
#' confidence intervals. Paired samples are resampled together. Use zero to
#' omit bootstrap intervals.
#' @param transform Proportion transformation, `"logit"` (default) or `"asin"`.
#' Logit uses speckle's count pseudocount of 0.5 before computing proportions.
#' @param robust Whether speckle uses robust empirical Bayes variance shrinkage.
#' @param trend Whether speckle fits a mean-variance trend.
#'
#' @return A method result bundle used internally by [RunProportionTest].
#'
#' @export
RunPropeller <- function(
  object,
  group.by,
  split.by,
  sample.by,
  comparison = NULL,
  n_bootstrap = 1000,
  seed = 11,
  verbose = TRUE,
  srt = NULL,
  transform = c("logit", "asin"),
  robust = TRUE,
  trend = FALSE
) {
  srt <- resolve_deprecated_srt(object, srt, missing(object))
  transform <- match.arg(transform)
  if (length(n_bootstrap) != 1L || is.na(n_bootstrap) || !is.finite(n_bootstrap) ||
      n_bootstrap < 0 || n_bootstrap != as.integer(n_bootstrap)) {
    stop("n_bootstrap must be a non-negative integer.", call. = FALSE)
  }
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

  propeller_check_r()
  engine <- "speckle"

  results_list <- list()
  backend_results <- list()
  for (i in seq_len(nrow(comparisons_condition))) {
    cluster_1 <- comparisons_condition[i, 1]
    cluster_2 <- comparisons_condition[i, 2]
    comparison_name <- paste0(cluster_1, "_vs_", cluster_2)

    prop_res <- run_propeller_speckle(
      meta_data = meta_data,
      group.by = group.by,
      split.by = split.by,
      sample.by = sample.by,
      cluster_1 = cluster_1,
      cluster_2 = cluster_2,
      n_bootstrap = n_bootstrap,
      transform = transform,
      robust = robust,
      trend = trend,
      seed = seed + i,
      verbose = verbose
    )

    results_list[[comparison_name]] <- standardize_proportion_result(
      prop_res$result,
      cluster_1 = cluster_1,
      cluster_2 = cluster_2,
      comparison_name = comparison_name,
      method = "propeller"
    )
    backend_results[[comparison_name]] <- prop_res$backend
  }

  list(
    method = "propeller",
    results = results_list,
    result_levels = "group",
    details = list(
      engine = engine,
      speckle_version = as.character(utils::packageVersion("speckle")),
      backend_results = backend_results
    ),
    parameters = list(
      sample.by = sample.by,
      n_bootstrap = n_bootstrap,
      transform = transform,
      robust = robust,
      trend = trend,
      engine = engine
    )
  )
}

propeller_check_r <- function() {
  status <- check_r("speckle", verbose = FALSE)
  if (!isTRUE(status[["speckle"]])) {
    stop("RunPropeller requires the optional Bioconductor package 'speckle'.", call. = FALSE)
  }
  invisible(TRUE)
}

propeller_get_fun <- function(fun) {
  out <- get_namespace_fun("speckle", fun)
  if (!is.function(out)) {
    stop(sprintf("The installed speckle backend does not provide '%s'.", fun), call. = FALSE)
  }
  out
}

run_propeller_speckle <- function(
  meta_data, group.by, split.by, sample.by, cluster_1, cluster_2,
  n_bootstrap, transform, robust, trend, seed, verbose
) {
  dat <- meta_data[, c(group.by, split.by, sample.by), drop = FALSE]
  colnames(dat) <- c("clusters", "condition", "sample")
  dat <- dat[dat$condition %in% c(cluster_1, cluster_2), , drop = FALSE]
  dat[] <- lapply(dat, as.character)
  if (anyNA(dat$clusters) || any(!nzchar(dat$clusters))) {
    stop("Propeller cell type labels must be non-missing and non-empty.", call. = FALSE)
  }
  sample_info <- proportion_sample_condition_keys(dat, "sample", "condition")
  pairs <- sample_info$pairs
  paired <- proportion_pair_is_paired(pairs, cluster_1, cluster_2)
  if (any(table(factor(pairs$condition, levels = c(cluster_1, cluster_2))) < 2L)) {
    stop("Propeller requires at least two biological samples per condition.", call. = FALSE)
  }
  if (length(unique(dat$clusters)) < 2L) {
    stop("Propeller requires at least two cell types in the comparison.", call. = FALSE)
  }

  props <- propeller_get_fun("getTransformedProps")(
    clusters = factor(dat$clusters), sample = factor(sample_info$cell_keys),
    transform = transform
  )
  pairs <- pairs[match(colnames(props$Proportions), pairs$sample_key), , drop = FALSE]
  condition <- factor(pairs$condition, levels = c(cluster_1, cluster_2))
  design <- stats::model.matrix(~ 0 + condition)
  if (paired) {
    donor <- stats::model.matrix(~ factor(pairs$sample))[, -1, drop = FALSE]
    design <- cbind(design, donor)
  }
  rownames(design) <- pairs$sample_key
  contrasts <- c(1, -1, rep(0, ncol(design) - 2L))
  upstream <- propeller_get_fun("propeller.ttest")(
    prop.list = props, design = design, contrasts = contrasts,
    robust = robust, trend = trend, sort = FALSE
  )

  result <- upstream
  result$engine <- "speckle"
  result$clusters <- rownames(upstream)
  prop_mat <- props$Proportions[result$clusters, , drop = FALSE]
  group_1 <- which(pairs$condition == cluster_1)
  group_2 <- which(pairs$condition == cluster_2)
  # Align donor IDs before both effect calculation and paired bootstrap.
  if (paired) {
    group_2 <- group_2[match(pairs$sample[group_1], pairs$sample[group_2])]
  }
  pseudocount <- 1e-5
  result$obs_log2FD <- log2((rowMeans(prop_mat[, group_1, drop = FALSE]) + pseudocount) /
    (rowMeans(prop_mat[, group_2, drop = FALSE]) + pseudocount))
  result$boot_mean_log2FD <- result$boot_CI_2.5 <- result$boot_CI_97.5 <- NA_real_
  if (!is.null(seed)) set.seed(seed)
  if (n_bootstrap > 0) {
    for (i in seq_len(nrow(result))) {
      v1 <- as.numeric(prop_mat[i, group_1])
      v2 <- as.numeric(prop_mat[i, group_2])
      boot_result <- if (paired) {
        boot <- replicate(n_bootstrap, {
          idx <- sample.int(length(v1), length(v1), replace = TRUE)
          log2((mean(v1[idx]) + pseudocount) / (mean(v2[idx]) + pseudocount))
        })
        list(boot_mean_log2FD = mean(boot),
          boot_CI_2.5 = as.numeric(stats::quantile(boot, 0.025)),
          boot_CI_97.5 = as.numeric(stats::quantile(boot, 0.975)))
      } else {
        proportion_bootstrap_stats(v1, v2, n_bootstrap, pseudocount, verbose)
      }
      boot_columns <- c("boot_mean_log2FD", "boot_CI_2.5", "boot_CI_97.5")
      result[i, boot_columns] <- boot_result[boot_columns]
    }
  }
  list(result = result, backend = list(
    results = upstream, design = design, contrasts = contrasts,
    sample_data = pairs, counts = props$Counts, paired = paired
  ))
}
