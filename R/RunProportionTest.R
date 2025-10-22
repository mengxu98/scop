#' @title Proportion Test
#'
#' @description
#' [RunProportionTest] performs a Monte-carlo permutation test to quantify the cell proportion differences between each condition.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt A Seurat object containing count data and metadata.
#' @param group.by The name of the metadata column in `srt` containing cluster labels or cell type names.
#' @param split.by The name of the metadata column in `srt` that contains sample identifiers.
#' @param comparison Optional: specify comparisons to perform.
#' @param n_permutations Number of permutations for the test.
#' @param FDR_threshold FDR value cutoff for significance.
#' @param log2FD_threshold Absolute value of log2FD cutoff for significance.
#' @param include_all_cells Whether to include all cell types in the complete grid (default: FALSE).
#'
#' @export
#'
#' @seealso
#' [ProportionTestPlot]
#'
#' @references
#' [Paper: Miller SA](https://doi.org/10.1158/0008-5472.can-20-3562),
#' [scProportionTest](https://github.com/rpolicastro/scProportionTest)
#'
#' @examples
#' data(pancreas_sub)
#' # Default behavior: only include cell types present in comparison groups
#' pancreas_sub <- RunProportionTest(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   split.by = "Phase",
#'   comparison = list(c("G2M", "G1"))
#' )
#'
#' # Include all cell types from the dataset
#' pancreas_sub <- RunProportionTest(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   split.by = "Phase",
#'   comparison = list(c("G2M", "G1")),
#'   include_all_cells = TRUE
#' )
#'
#' ProportionTestPlot(
#'   pancreas_sub
#' )
RunProportionTest <- function(
    srt,
    group.by,
    split.by,
    comparison = NULL,
    n_permutations = 1000,
    FDR_threshold = 0.05,
    log2FD_threshold = log2(1.5),
    include_all_cells = FALSE,
    verbose = TRUE) {
  meta_data <- srt@meta.data
  if (!split.by %in% colnames(meta_data)) {
    log_message(
      "{.val {split.by}} does not exist in the seurat object meta data",
      message_type = "error"
    )
  }
  if (!group.by %in% colnames(meta_data)) {
    log_message(
      "{.val {group.by}} does not exist in the seurat object meta data",
      message_type = "error"
    )
  }

  conditions <- as.character(unique(meta_data[, split.by]))

  if (!is.null(comparison)) {
    if (is.list(comparison)) {
      if (all(sapply(comparison, function(x) is.character(x) && length(x) == 2))) {
        comparisons_condition <- do.call(rbind, comparison)
      } else {
        comparison <- unlist(comparison)
        comparisons_condition <- do.call(
          rbind, lapply(comparison, function(x) {
            strsplit(x, "_vs_")[[1]]
          })
        )
      }
    } else if (is.character(comparison)) {
      comparisons_condition <- do.call(
        rbind, lapply(comparison, function(x) {
          strsplit(x, "_vs_")[[1]]
        })
      )
    }
  } else {
    n <- length(conditions)
    r <- 2
    if (r > n) {
      comparisons_condition <- matrix(nrow = 0, ncol = r)
    } else {
      comparisons_condition <- expand.grid(
        rep(list(seq_len(n)), r)
      )
      comparisons_condition <- comparisons_condition[
        apply(comparisons_condition, 1, function(x) length(unique(x)) == r),
        , drop = FALSE
      ]
      comparisons_condition <- matrix(
        conditions[as.matrix(comparisons_condition)],
        ncol = r
      )
    }
  }

  log_message(
    "Start proportion test",
    verbose = verbose
  )

  results_list <- list()

  for (i in seq_len(nrow(comparisons_condition))) {
    for (direction in 1:2) {
      if (direction == 1) {
        cluster_1 <- comparisons_condition[i, 1]
        cluster_2 <- comparisons_condition[i, 2]
        comparison_name <- paste0(cluster_1, "_vs_", cluster_2)
      } else {
        cluster_1 <- comparisons_condition[i, 2]
        cluster_2 <- comparisons_condition[i, 1]
        comparison_name <- paste0(cluster_1, "_vs_", cluster_2)
      }

      log_message(
        "Running comparison: {.val {cluster_2}} vs {.val {cluster_1}}",
        verbose = verbose
      )

      test_result <- .permutation_test(
        srt = srt,
        group.by = group.by,
        split.by = split.by,
        cluster_1 = cluster_1,
        cluster_2 = cluster_2,
        n_permutations = n_permutations,
        include_all_cells = include_all_cells
      )

      results_list[[comparison_name]] <- test_result
      results_list[[comparison_name]]$group1 <- cluster_1
      results_list[[comparison_name]]$group2 <- cluster_2
      results_list[[comparison_name]]$comparison <- comparison_name
    }
  }

  srt@tools[["ProportionTest"]][["results"]] <- results_list
  srt@tools[["ProportionTest"]][["parameters"]] <- list(
    group.by = group.by,
    split.by = split.by,
    n_permutations = n_permutations,
    FDR_threshold = FDR_threshold,
    log2FD_threshold = log2FD_threshold,
    include_all_cells = include_all_cells
  )

  log_message(
    "Proportion test completed",
    message_type = "success",
    verbose = verbose
  )

  return(srt)
}

.permutation_test <- function(
    srt,
    group.by,
    split.by,
    cluster_1,
    cluster_2,
    n_permutations,
    include_all_cells = FALSE) {
  meta_data <- srt@meta.data
  meta_data <- meta_data[, c(split.by, group.by)]

  colnames(meta_data) <- c("samples", "clusters")

  meta_data$clusters <- as.character(meta_data$clusters)
  comparison_data <- meta_data[meta_data$samples %in% c(cluster_1, cluster_2), ]
  if (include_all_cells) {
    cluster_cases <- unique(meta_data$clusters)
  } else {
    cluster_cases <- unique(comparison_data$clusters)
  }

  comparison_data$count <- 1
  obs_diff <- stats::aggregate(
    count ~ samples + clusters,
    data = comparison_data,
    FUN = length,
    drop = FALSE
  )

  if (include_all_cells) {
    all_clusters <- unique(meta_data$clusters)
    complete_grid <- expand.grid(
      samples = c(cluster_1, cluster_2),
      clusters = all_clusters,
      stringsAsFactors = FALSE
    )
  } else {
    complete_grid <- expand.grid(
      samples = c(cluster_1, cluster_2),
      clusters = cluster_cases,
      stringsAsFactors = FALSE
    )
  }

  obs_diff <- merge(
    complete_grid, obs_diff,
    by = c("samples", "clusters"), all.x = TRUE
  )
  obs_diff$count[is.na(obs_diff$count)] <- 0

  sample_totals <- stats::aggregate(
    count ~ samples,
    data = obs_diff,
    FUN = sum
  )
  obs_diff <- merge(
    obs_diff, sample_totals,
    by = "samples",
    suffixes = c("", "_total")
  )
  obs_diff$fraction <- obs_diff$count / obs_diff$count_total

  obs_diff_wide <- stats::reshape(
    obs_diff[, c("clusters", "samples", "fraction")],
    idvar = "clusters",
    timevar = "samples",
    direction = "wide"
  )

  colnames(obs_diff_wide)[-1] <- gsub(
    "fraction.", "", colnames(obs_diff_wide)[-1]
  )

  obs_diff_wide$obs_log2FD <- log2(obs_diff_wide[[cluster_2]]) - log2(obs_diff_wide[[cluster_1]])

  perm_results <- matrix(NA, nrow(obs_diff_wide), n_permutations)
  rownames(perm_results) <- sort(cluster_cases)

  for (i in seq_len(n_permutations)) {
    permuted <- comparison_data
    permuted$samples <- sample(permuted$samples)

    permuted$count <- 1
    permuted_count <- stats::aggregate(
      count ~ samples + clusters,
      data = permuted,
      FUN = length,
      drop = FALSE
    )

    permuted_count <- merge(
      complete_grid, permuted_count,
      by = c("samples", "clusters"), all.x = TRUE
    )
    permuted_count$count[is.na(permuted_count$count)] <- 0

    sample_totals_perm <- stats::aggregate(
      count ~ samples,
      data = permuted_count, FUN = sum
    )
    permuted_count <- merge(
      permuted_count, sample_totals_perm,
      by = "samples", suffixes = c("", "_total")
    )
    permuted_count$fraction <- permuted_count$count / permuted_count$count_total

    permuted_wide <- stats::reshape(
      permuted_count[, c("clusters", "samples", "fraction")],
      idvar = "clusters",
      timevar = "samples",
      direction = "wide"
    )

    colnames(permuted_wide)[-1] <- gsub(
      "fraction.", "", colnames(permuted_wide)[-1]
    )

    permuted_wide$perm_log2FD <- log2(
      permuted_wide[[cluster_2]]
    ) - log2(permuted_wide[[cluster_1]])

    perm_results[, i] <- permuted_wide$perm_log2FD
  }

  increased <- rowSums(
    apply(perm_results, 2, function(x) obs_diff_wide$obs_log2FD <= x)
  )
  increased <- (increased + 1) / (n_permutations + 1)

  decreased <- rowSums(
    apply(perm_results, 2, function(x) obs_diff_wide$obs_log2FD >= x)
  )
  decreased <- (decreased + 1) / (n_permutations + 1)

  obs_diff_wide$pval <- ifelse(
    obs_diff_wide$obs_log2FD > 0, increased, decreased
  )
  obs_diff_wide$FDR <- stats::p.adjust(obs_diff_wide$pval, "fdr")

  boot_results <- matrix(NA, nrow(obs_diff_wide), n_permutations)
  rownames(boot_results) <- sort(cluster_cases)

  for (i in seq_len(n_permutations)) {
    booted <- comparison_data

    for (sample in unique(booted$samples)) {
      sample_idx <- booted$samples == sample
      booted$clusters[sample_idx] <- sample(
        booted$clusters[sample_idx],
        replace = TRUE
      )
    }

    booted$count <- 1
    booted_count <- stats::aggregate(
      count ~ samples + clusters,
      data = booted,
      FUN = length,
      drop = FALSE
    )

    booted_count <- merge(
      complete_grid, booted_count,
      by = c("samples", "clusters"), all.x = TRUE
    )
    booted_count$count[is.na(booted_count$count)] <- 0

    sample_totals_boot <- stats::aggregate(
      count ~ samples,
      data = booted_count, FUN = sum
    )
    booted_count <- merge(
      booted_count, sample_totals_boot,
      by = "samples", suffixes = c("", "_total")
    )
    booted_count$fraction <- booted_count$count / booted_count$count_total

    booted_wide <- stats::reshape(
      booted_count[, c("clusters", "samples", "fraction")],
      idvar = "clusters",
      timevar = "samples",
      direction = "wide"
    )

    colnames(booted_wide)[-1] <- gsub(
      "fraction.", "", colnames(booted_wide)[-1]
    )

    booted_wide$boot_log2FD <- log2(
      booted_wide[[cluster_2]]
    ) - log2(booted_wide[[cluster_1]])

    boot_results[, i] <- booted_wide$boot_log2FD
  }

  boot_results[!is.finite(boot_results)] <- NA
  boot_mean <- rowMeans(boot_results, na.rm = TRUE)
  boot_ci <- t(
    apply(
      boot_results, 1, function(x) {
        stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
      }
    )
  )
  boot_ci <- as.data.frame(boot_ci)
  colnames(boot_ci) <- c("boot_CI_2.5", "boot_CI_97.5")

  obs_diff_wide$boot_mean_log2FD <- boot_mean
  obs_diff_wide <- cbind(obs_diff_wide, boot_ci)

  return(obs_diff_wide)
}
