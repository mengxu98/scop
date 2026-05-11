#' @title Proportion Test
#'
#' @description
#' [RunProportionTest] performs differential abundance testing for cell proportions.
#' The function acts as a dispatcher and routes to one of the method-specific
#' implementations: permutation, milo, sccoda, or propeller.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams CellDimPlot
#' @param comparison Optional: specify comparisons to perform.
#' @param proportion_method Differential abundance method.
#' One of `"permutation"`, `"milo"`, `"sccoda"`, or `"propeller"`.
#' This argument is required.
#' Alias values such as `"permutation_test"` and `"perm"` are accepted and
#' normalized to `"permutation"`.
#' @param split.by Metadata column that identifies the condition groups to
#' compare. For sample-level methods, if `split.by` is omitted and `sample.by`
#' is provided, `sample.by` is treated as the condition column and virtual
#' samples are created within each condition.
#' @param sample.by Metadata column that identifies biological samples.
#' For `"milo"`, `"sccoda"`, and `"propeller"`, when `sample.by` is omitted
#' or identical to `split.by`, virtual samples are created within each
#' `split.by` group for convenience.
#' @param pseudo_sample_n Number of virtual samples per `split.by` group when
#' a sample-level method has no usable `sample.by`.
#' @param n_permutations Number of permutations for permutation-based test.
#' @param FDR_threshold FDR value cutoff for significance.
#' @param log2FD_threshold Absolute value of log2FD cutoff for significance.
#' @param include_all_cells Whether to include all cell types in the complete grid
#' for permutation mode.
#' @param seed Random seed.
#' @param ... Additional arguments passed to the selected method function.
#'
#' @export
#'
#' @seealso
#' [RunPermutation], [RunMilo], [RunscCODA], [RunPropeller],
#' [ProportionTestPlot]
#'
#' @references
#' [Miller et al. paper](https://doi.org/10.1158/0008-5472.can-20-3562),
#' [scProportionTest](https://github.com/rpolicastro/scProportionTest),
#' [miloR](https://bioconductor.org/packages/miloR),
#' [scCODA](https://github.com/theislab/scCODA),
#' [propeller/speckle](https://bioconductor.org/packages/speckle)
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- RunProportionTest(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   split.by = "Phase",
#'   proportion_method = "permutation",
#'   comparison = list(c("G2M", "G1"))
#' )
#'
#' ProportionTestPlot(
#'   pancreas_sub
#' )
RunProportionTest <- function(
  srt,
  group.by,
  split.by = NULL,
  comparison = NULL,
  proportion_method,
  sample.by = NULL,
  pseudo_sample_n = 3L,
  n_permutations = 1000,
  FDR_threshold = 0.05,
  log2FD_threshold = log2(1.5),
  include_all_cells = FALSE,
  seed = 11,
  verbose = TRUE,
  ...
) {
  proportion_method <- normalize_proportion_method(proportion_method)

  pseudo_sample_info <- list(enabled = FALSE)
  if (is.null(split.by)) {
    if (!is.null(sample.by) && nzchar(sample.by)) {
      split.by <- sample.by
      sample.by <- NULL
    } else {
      log_message(
        "{.arg split.by} must be provided unless {.arg sample.by} names the condition column",
        message_type = "error"
      )
    }
  }

  require_sample <- proportion_method %in% c("milo", "sccoda", "propeller")
  if (
    require_sample &&
      (is.null(sample.by) ||
        !nzchar(sample.by) ||
        identical(sample.by, split.by))
  ) {
    pseudo <- add_proportion_pseudo_samples(
      srt = srt,
      split.by = split.by,
      pseudo_sample_n = pseudo_sample_n,
      proportion_method = proportion_method,
      seed = seed,
      verbose = verbose
    )
    srt <- pseudo$srt
    sample.by <- pseudo$sample.by
    pseudo_sample_info <- pseudo$info
  }
  validate_proportion_inputs(
    srt = srt,
    group.by = group.by,
    split.by = split.by,
    sample.by = sample.by,
    require_sample = require_sample
  )

  method_map <- list(
    permutation = RunPermutation,
    milo = RunMilo,
    sccoda = RunscCODA,
    propeller = RunPropeller
  )
  method_fun <- method_map[[proportion_method]]
  if (is.null(method_fun)) {
    log_message(
      "Unsupported {.arg proportion_method}: {.val {proportion_method}}",
      message_type = "error"
    )
  }

  method_args <- utils::modifyList(
    list(
      srt = srt,
      group.by = group.by,
      split.by = split.by,
      comparison = comparison,
      sample.by = sample.by,
      pseudo_sample_n = pseudo_sample_n,
      n_permutations = n_permutations,
      include_all_cells = include_all_cells,
      seed = seed,
      verbose = verbose
    ),
    list(...)
  )

  log_message(
    "Start proportion test ({.val {proportion_method}})",
    verbose = verbose
  )

  method_bundle <- invoke_fun(
    method_fun,
    method_args[names(method_args) %in% names(formals(method_fun))]
  )

  if (!is.list(method_bundle) || is.null(method_bundle[["results"]])) {
    log_message(
      "Method {.val {proportion_method}} did not return valid results",
      message_type = "error"
    )
  }

  method_bundle$method <- proportion_method
  method_bundle$result_levels <- method_bundle$result_levels %||% "group"

  run_params <- list(
    group.by = group.by,
    split.by = split.by,
    sample.by = sample.by,
    pseudo_sample_n = pseudo_sample_n,
    pseudo_sample = pseudo_sample_info,
    comparison = comparison,
    n_permutations = n_permutations,
    FDR_threshold = FDR_threshold,
    log2FD_threshold = log2FD_threshold,
    include_all_cells = include_all_cells,
    seed = seed,
    proportion_method = proportion_method
  )

  method_bundle$parameters <- utils::modifyList(
    run_params,
    method_bundle$parameters %||% list()
  )

  if (isTRUE(pseudo_sample_info$enabled)) {
    srt@meta.data[[pseudo_sample_info$column]] <- NULL
  }

  if (is.null(srt@tools[["ProportionTest"]])) {
    srt@tools[["ProportionTest"]] <- list()
  }

  methods_store <- srt@tools[["ProportionTest"]][["methods"]] %||% list()
  methods_store[[proportion_method]] <- method_bundle

  srt@tools[["ProportionTest"]][["results"]] <- method_bundle[["results"]]
  srt@tools[["ProportionTest"]][["parameters"]] <- run_params
  srt@tools[["ProportionTest"]][["active_method"]] <- proportion_method
  srt@tools[["ProportionTest"]][["methods"]] <- methods_store

  log_message(
    "Proportion test completed ({.val {proportion_method}})",
    message_type = "success",
    verbose = verbose
  )

  return(srt)
}

#' @title Permutation-based proportion test
#'
#' @description
#' Method-specific implementation used by [RunProportionTest] when
#' `proportion_method = "permutation"`.
#' This method is a permutation-based statistical test for compositional shift
#' rather than a dedicated biological model.
#' Bootstrap confidence intervals are reported as uncertainty estimates for
#' observed log2 fold-differences.
#'
#' @md
#' @inheritParams RunProportionTest
#' @return A method result bundle used internally by [RunProportionTest].
#'
#' @export
RunPermutation <- function(
  srt,
  group.by,
  split.by,
  comparison = NULL,
  n_permutations = 1000,
  include_all_cells = FALSE,
  verbose = TRUE
) {
  meta_data <- validate_proportion_inputs(
    srt = srt,
    group.by = group.by,
    split.by = split.by,
    require_sample = FALSE
  )

  comparisons_condition <- parse_proportion_comparisons(
    meta_data = meta_data,
    split.by = split.by,
    comparison = comparison,
    include_bidirectional = TRUE
  )

  results_list <- list()
  for (i in seq_len(nrow(comparisons_condition))) {
    cluster_1 <- comparisons_condition[i, 1]
    cluster_2 <- comparisons_condition[i, 2]
    comparison_name <- paste0(cluster_1, "_vs_", cluster_2)

    log_message(
      "Running comparison: {.val {cluster_2}} vs {.val {cluster_1}}",
      verbose = verbose
    )

    test_result <- permutation_test(
      srt = srt,
      group.by = group.by,
      split.by = split.by,
      cluster_1 = cluster_1,
      cluster_2 = cluster_2,
      n_permutations = n_permutations,
      include_all_cells = include_all_cells,
      verbose = verbose
    )

    results_list[[comparison_name]] <- standardize_proportion_result(
      test_result,
      cluster_1 = cluster_1,
      cluster_2 = cluster_2,
      comparison_name = comparison_name,
      method = "permutation"
    )
  }

  list(
    method = "permutation",
    results = results_list,
    result_levels = "group",
    details = list(),
    parameters = list(
      n_permutations = n_permutations,
      include_all_cells = include_all_cells
    )
  )
}

normalize_proportion_method <- function(proportion_method) {
  if (missing(proportion_method) || is.null(proportion_method)) {
    log_message(
      "{.arg proportion_method} must be provided",
      message_type = "error"
    )
  }

  method <- tolower(as.character(proportion_method)[1])
  alias <- c(
    scproportiontest = "permutation",
    scproportion_method = "permutation",
    montecarlo = "permutation",
    perm = "permutation",
    permutation_test = "permutation",
    runpermutation = "permutation",
    run_permutation = "permutation",
    runmilo = "milo",
    run_milo = "milo",
    runsccoda = "sccoda",
    run_sccoda = "sccoda",
    runpropeller = "propeller",
    run_propeller = "propeller"
  )
  if (method %in% names(alias)) {
    method <- unname(alias[[method]])
  }

  allowed <- c("permutation", "milo", "sccoda", "propeller")
  if (!method %in% allowed) {
    log_message(
      "{.arg proportion_method} must be one of {.val {allowed}}",
      message_type = "error"
    )
  }

  method
}

validate_proportion_inputs <- function(
  srt,
  group.by,
  split.by,
  sample.by = NULL,
  require_sample = FALSE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }

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

  if (isTRUE(require_sample)) {
    if (is.null(sample.by) || !nzchar(sample.by)) {
      log_message(
        "{.arg sample.by} must be provided for this {.arg proportion_method}",
        message_type = "error"
      )
    }
    if (!sample.by %in% colnames(meta_data)) {
      log_message(
        "{.val {sample.by}} does not exist in the seurat object meta data",
        message_type = "error"
      )
    }
  }

  meta_data
}

add_proportion_pseudo_samples <- function(
  srt,
  split.by,
  pseudo_sample_n = 3L,
  proportion_method = NULL,
  seed = 11,
  verbose = TRUE
) {
  meta_data <- srt@meta.data
  if (!split.by %in% colnames(meta_data)) {
    log_message(
      "{.val {split.by}} does not exist in the seurat object meta data",
      message_type = "error"
    )
  }

  pseudo_sample_n <- as.integer(pseudo_sample_n[1])
  if (is.na(pseudo_sample_n) || pseudo_sample_n < 2L) {
    log_message(
      "{.arg pseudo_sample_n} must be an integer greater than or equal to 2",
      message_type = "error"
    )
  }

  split_values <- as.character(meta_data[[split.by]])
  valid_values <- unique(split_values[!is.na(split_values)])
  if (length(valid_values) < 2) {
    log_message(
      "Need at least two groups in {.arg split.by} to construct virtual samples",
      message_type = "error"
    )
  }

  tab <- table(split_values)
  small_groups <- names(tab)[tab < 2]
  if (length(small_groups) > 0) {
    log_message(
      "Cannot construct virtual samples because these {.arg split.by} groups contain fewer than 2 cells: {.val {small_groups}}",
      message_type = "error"
    )
  }

  column <- ".scop_pseudo_sample"
  idx <- 1L
  while (column %in% colnames(meta_data)) {
    idx <- idx + 1L
    column <- paste0(".scop_pseudo_sample", idx)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  pseudo_sample <- rep(NA_character_, nrow(meta_data))
  names(pseudo_sample) <- rownames(meta_data)
  for (condition in valid_values) {
    cells <- which(split_values == condition)
    cells <- sample(cells)
    n_rep <- min(pseudo_sample_n, length(cells))
    rep_id <- rep(seq_len(n_rep), length.out = length(cells))
    rep_id <- sample(rep_id, length(rep_id))
    prefix <- make.names(condition)
    pseudo_sample[cells] <- paste0(prefix, "_pseudo", rep_id)
  }

  srt@meta.data[[column]] <- pseudo_sample

  method_label <- proportion_method %||% "sample-level proportion testing"
  log_message(
    "Construct virtual {.arg sample.by} column {.val {column}} within {.arg split.by} groups for {.val {method_label}}; use biological sample metadata when available.",
    message_type = "warning",
    verbose = verbose
  )

  list(
    srt = srt,
    sample.by = column,
    info = list(
      enabled = TRUE,
      column = column,
      split.by = split.by,
      proportion_method = proportion_method,
      pseudo_sample_n = pseudo_sample_n
    )
  )
}

parse_proportion_comparisons <- function(
  meta_data,
  split.by,
  comparison = NULL,
  include_bidirectional = TRUE
) {
  coerce_comparison_pairs <- function(comparison) {
    parse_vs <- function(x) {
      parts <- strsplit(x, "_vs_", fixed = TRUE)[[1]]
      if (length(parts) != 2) {
        log_message(
          "Invalid comparison value: {.val {x}}. Expected format {.val A_vs_B}",
          message_type = "error"
        )
      }
      trimws(parts)
    }

    if (is.list(comparison)) {
      if (
        all(vapply(
          comparison,
          function(x) is.character(x) && length(x) == 2,
          logical(1)
        ))
      ) {
        out <- do.call(rbind, lapply(comparison, function(x) trimws(x)))
      } else {
        comp_vec <- unlist(comparison)
        out <- do.call(rbind, lapply(comp_vec, parse_vs))
      }
    } else if (is.character(comparison)) {
      if (
        length(comparison) == 2 && !any(grepl("_vs_", comparison, fixed = TRUE))
      ) {
        out <- matrix(trimws(comparison), nrow = 1)
      } else {
        out <- do.call(rbind, lapply(comparison, parse_vs))
      }
    } else {
      log_message(
        "{.arg comparison} must be NULL, a character vector, or a list",
        message_type = "error"
      )
    }

    if (ncol(out) != 2) {
      log_message(
        "{.arg comparison} must define two groups for each comparison",
        message_type = "error"
      )
    }

    out
  }

  conditions <- unique(as.character(meta_data[[split.by]]))
  conditions <- conditions[!is.na(conditions)]

  if (length(conditions) < 2) {
    log_message(
      "Need at least two groups in {.arg split.by} to run proportion test",
      message_type = "error"
    )
  }

  if (is.null(comparison)) {
    base_pairs <- t(utils::combn(conditions, 2))
  } else {
    base_pairs <- coerce_comparison_pairs(comparison)

    missing_groups <- setdiff(unique(as.vector(base_pairs)), conditions)
    if (length(missing_groups) > 0) {
      log_message(
        "Comparison group{?s} not found in {.arg split.by}: {.val {missing_groups}}",
        message_type = "error"
      )
    }
  }

  if (isTRUE(include_bidirectional)) {
    reverse_pairs <- base_pairs[, c(2, 1), drop = FALSE]
    pairs <- unique(rbind(base_pairs, reverse_pairs))
  } else {
    pairs <- unique(base_pairs)
  }

  colnames(pairs) <- c("group1", "group2")
  pairs
}

sample_level_proportion_test <- function(
  meta_data,
  group.by,
  split.by,
  sample.by,
  cluster_1,
  cluster_2,
  n_bootstrap = 1000,
  transform = c("raw", "logit", "asin"),
  pseudocount = 1e-5,
  seed = NULL,
  verbose = FALSE
) {
  transform <- match.arg(transform)

  dat <- meta_data[, c(group.by, split.by, sample.by), drop = FALSE]
  colnames(dat) <- c("clusters", "condition", "sample")
  dat$clusters <- as.character(dat$clusters)
  dat$condition <- as.character(dat$condition)
  dat$sample <- as.character(dat$sample)

  dat <- dat[dat$condition %in% c(cluster_1, cluster_2), , drop = FALSE]

  sample_condition <- stats::aggregate(
    condition ~ sample,
    data = dat,
    FUN = function(x) {
      ux <- unique(x)
      ux[1]
    }
  )
  rownames(sample_condition) <- sample_condition$sample

  sample_levels <- sample_condition$sample
  cluster_levels <- sort(unique(dat$clusters))

  count_tab <- stats::xtabs(~ sample + clusters, data = dat)
  count_mat <- matrix(
    0,
    nrow = length(sample_levels),
    ncol = length(cluster_levels),
    dimnames = list(sample_levels, cluster_levels)
  )
  count_mat[rownames(count_tab), colnames(count_tab)] <- as.matrix(count_tab)

  sample_totals <- rowSums(count_mat)
  sample_totals[sample_totals == 0] <- 1
  prop_mat <- count_mat / sample_totals

  samples_1 <- sample_condition$sample[sample_condition$condition == cluster_1]
  samples_2 <- sample_condition$sample[sample_condition$condition == cluster_2]

  trans_fun <- switch(
    transform,
    raw = function(x) x,
    logit = function(x) log((x + pseudocount) / (1 - x + pseudocount)),
    asin = function(x) asin(sqrt(x))
  )

  if (!is.null(seed)) {
    set.seed(seed)
  }

  out <- lapply(cluster_levels, function(ct) {
    v1 <- as.numeric(prop_mat[samples_1, ct, drop = TRUE])
    v2 <- as.numeric(prop_mat[samples_2, ct, drop = TRUE])

    m1 <- mean(v1, na.rm = TRUE)
    m2 <- mean(v2, na.rm = TRUE)
    obs_log2FD <- log2((m2 + pseudocount) / (m1 + pseudocount))

    t1 <- trans_fun(v1)
    t2 <- trans_fun(v2)

    pval <- tryCatch(
      {
        if (length(v1) == 0 || length(v2) == 0) {
          NA_real_
        } else if (length(unique(c(v1, v2))) <= 1) {
          1
        } else if (length(v1) >= 2 && length(v2) >= 2) {
          stats::t.test(t2, t1)$p.value
        } else {
          stats::wilcox.test(t2, t1, exact = FALSE)$p.value
        }
      },
      error = function(e) NA_real_
    )

    if (n_bootstrap > 0 && length(v1) > 0 && length(v2) > 0) {
      boot_result <- proportion_bootstrap_stats(
        v1 = v1,
        v2 = v2,
        n_bootstrap = n_bootstrap,
        pseudocount = pseudocount,
        verbose = verbose
      )
      boot_mean_log2FD <- boot_result[["boot_mean_log2FD"]]
      boot_CI_2.5 <- boot_result[["boot_CI_2.5"]]
      boot_CI_97.5 <- boot_result[["boot_CI_97.5"]]
    } else {
      boot_mean_log2FD <- NA_real_
      boot_CI_2.5 <- NA_real_
      boot_CI_97.5 <- NA_real_
    }

    data.frame(
      clusters = ct,
      obs_log2FD = obs_log2FD,
      pval = pval,
      boot_mean_log2FD = boot_mean_log2FD,
      boot_CI_2.5 = boot_CI_2.5,
      boot_CI_97.5 = boot_CI_97.5,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  out$FDR <- stats::p.adjust(out$pval, method = "fdr")
  out
}

standardize_proportion_result <- function(
  x,
  cluster_1,
  cluster_2,
  comparison_name,
  method
) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)

  cluster_col <- find_first_col(
    df,
    c(
      "clusters",
      "cluster",
      "cell_type",
      "CellType",
      "celltype",
      "feature",
      "name"
    )
  )
  if (is.null(cluster_col)) {
    if (
      !is.null(rownames(df)) &&
        !all(rownames(df) %in% c("", as.character(seq_len(nrow(df)))))
    ) {
      df$clusters <- rownames(df)
    } else {
      df$clusters <- paste0("cluster_", seq_len(nrow(df)))
    }
  } else if (cluster_col != "clusters") {
    df$clusters <- as.character(df[[cluster_col]])
  }

  obs_col <- find_first_col(
    df,
    c("obs_log2FD", "log2FD", "logFC", "coef", "effect")
  )
  if (is.null(obs_col)) {
    df$obs_log2FD <- NA_real_
  } else {
    df$obs_log2FD <- suppressWarnings(as.numeric(df[[obs_col]]))
  }

  pval_col <- find_first_col(
    df,
    c("pval", "p_val", "p.value", "PValue", "P.Value")
  )
  if (is.null(pval_col)) {
    df$pval <- NA_real_
  } else {
    df$pval <- suppressWarnings(as.numeric(df[[pval_col]]))
  }

  fdr_col <- find_first_col(
    df,
    c("FDR", "adj.P.Val", "SpatialFDR", "q_value", "qval")
  )
  if (is.null(fdr_col)) {
    df$FDR <- stats::p.adjust(df$pval, method = "fdr")
  } else {
    df$FDR <- suppressWarnings(as.numeric(df[[fdr_col]]))
  }

  boot_mean_col <- find_first_col(
    df,
    c("boot_mean_log2FD", "boot_mean", "mean_log2FD")
  )
  if (is.null(boot_mean_col)) {
    df$boot_mean_log2FD <- NA_real_
  } else {
    df$boot_mean_log2FD <- suppressWarnings(as.numeric(df[[boot_mean_col]]))
  }

  ci_low_col <- find_first_col(
    df,
    c("boot_CI_2.5", "CI_2.5", "ci_low", "lower", "hdi_2.5")
  )
  if (is.null(ci_low_col)) {
    df$boot_CI_2.5 <- NA_real_
  } else {
    df$boot_CI_2.5 <- suppressWarnings(as.numeric(df[[ci_low_col]]))
  }

  ci_high_col <- find_first_col(
    df,
    c("boot_CI_97.5", "CI_97.5", "ci_high", "upper", "hdi_97.5")
  )
  if (is.null(ci_high_col)) {
    df$boot_CI_97.5 <- NA_real_
  } else {
    df$boot_CI_97.5 <- suppressWarnings(as.numeric(df[[ci_high_col]]))
  }

  effect_col <- find_first_col(
    df,
    c("effect", "coef", "logFC", "log2FD", "obs_log2FD")
  )
  if (is.null(effect_col)) {
    df$effect <- df$obs_log2FD
  } else {
    df$effect <- suppressWarnings(as.numeric(df[[effect_col]]))
  }

  ip_col <- find_first_col(
    df,
    c("inclusion_prob", "inclusion_probability", "probability", "pip")
  )
  if (is.null(ip_col)) {
    df$inclusion_prob <- NA_real_
  } else {
    df$inclusion_prob <- suppressWarnings(as.numeric(df[[ip_col]]))
  }

  credible_col <- find_first_col(
    df,
    c("credible", "is_credible", "is_significant")
  )
  if (is.null(credible_col)) {
    df$credible <- !is.na(df$inclusion_prob) & df$inclusion_prob >= 0.95
  } else {
    df$credible <- as.logical(df[[credible_col]])
  }

  hdi_low_col <- find_first_col(
    df,
    c("hdi_2.5", "hdi_low", "hdi_lower", "ci_low", "boot_CI_2.5")
  )
  if (is.null(hdi_low_col)) {
    df$hdi_2.5 <- suppressWarnings(as.numeric(df$boot_CI_2.5))
  } else {
    df$hdi_2.5 <- suppressWarnings(as.numeric(df[[hdi_low_col]]))
  }

  hdi_high_col <- find_first_col(
    df,
    c("hdi_97.5", "hdi_high", "hdi_upper", "ci_high", "boot_CI_97.5")
  )
  if (is.null(hdi_high_col)) {
    df$hdi_97.5 <- suppressWarnings(as.numeric(df$boot_CI_97.5))
  } else {
    df$hdi_97.5 <- suppressWarnings(as.numeric(df[[hdi_high_col]]))
  }

  nhood_col <- find_first_col(
    df,
    c("neighborhood", "Nhood", "nhood", "nhood_index")
  )
  if (!is.null(nhood_col) && nhood_col != "neighborhood") {
    df$neighborhood <- as.character(df[[nhood_col]])
  }

  df$group1 <- cluster_1
  df$group2 <- cluster_2
  df$comparison <- comparison_name
  df$method <- method

  core_cols <- c(
    "clusters",
    "obs_log2FD",
    "pval",
    "FDR",
    "boot_mean_log2FD",
    "boot_CI_2.5",
    "boot_CI_97.5",
    "group1",
    "group2",
    "comparison",
    "method"
  )
  for (nm in core_cols) {
    if (!nm %in% colnames(df)) {
      df[[nm]] <- NA
    }
  }

  df <- df[,
    unique(c(core_cols, setdiff(colnames(df), core_cols))),
    drop = FALSE
  ]
  rownames(df) <- NULL
  df
}

find_first_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    return(NULL)
  }
  hit[1]
}

permutation_test <- function(
  srt,
  group.by,
  split.by,
  cluster_1,
  cluster_2,
  n_permutations,
  include_all_cells = FALSE,
  pseudocount = 1e-8,
  verbose = FALSE
) {
  meta_data <- srt@meta.data[, c(split.by, group.by), drop = FALSE]
  colnames(meta_data) <- c("samples", "clusters")
  meta_data[["clusters"]] <- as.character(meta_data[["clusters"]])
  comparison_data <- meta_data[
    meta_data[["samples"]] %in% c(cluster_1, cluster_2),
    ,
    drop = FALSE
  ]
  if (isTRUE(include_all_cells)) {
    cluster_cases <- unique(meta_data[["clusters"]])
  } else {
    cluster_cases <- unique(comparison_data[["clusters"]])
  }

  sample_ids <- ifelse(comparison_data[["samples"]] == cluster_1, 1L, 2L)
  cluster_ids <- match(comparison_data[["clusters"]], cluster_cases)
  res <- proportion_permutation(
    sample_ids = as.integer(sample_ids),
    cluster_ids = as.integer(cluster_ids),
    cluster_levels = as.character(cluster_cases),
    n_permutations = as.integer(n_permutations),
    pseudocount = pseudocount,
    verbose = verbose
  )
  colnames(res)[colnames(res) == "fraction_1"] <- cluster_1
  colnames(res)[colnames(res) == "fraction_2"] <- cluster_2
  res[["FDR"]] <- stats::p.adjust(res[["pval"]], "fdr")
  res[,
    c(
      "clusters",
      cluster_1,
      cluster_2,
      "obs_log2FD",
      "pval",
      "FDR",
      "boot_mean_log2FD",
      "boot_CI_2.5",
      "boot_CI_97.5"
    ),
    drop = FALSE
  ]
}
