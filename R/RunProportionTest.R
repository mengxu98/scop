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
#' @param sample.by Metadata column that identifies biological samples.
#' Required for `"milo"`, `"sccoda"`, and `"propeller"`.
#' @param n_permutations Number of permutations for permutation-based test.
#' Also used as default bootstrap iterations in method fallbacks.
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
#' [RunProportionTestPermutation], [RunProportionTestMilo],
#' [RunProportionTestScCODA], [RunProportionTestPropeller], [ProportionTestPlot]
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
    split.by,
    comparison = NULL,
    proportion_method,
    sample.by = NULL,
    n_permutations = 1000,
    FDR_threshold = 0.05,
    log2FD_threshold = log2(1.5),
    include_all_cells = FALSE,
    seed = 11,
    verbose = TRUE,
    ...) {
  proportion_method <- .normalize_proportion_method(proportion_method)

  require_sample <- proportion_method %in% c("milo", "sccoda", "propeller")
  .validate_proportion_inputs(
    srt = srt,
    group.by = group.by,
    split.by = split.by,
    sample.by = sample.by,
    require_sample = require_sample
  )

  method_map <- list(
    permutation = RunProportionTestPermutation,
    milo = RunProportionTestMilo,
    sccoda = RunProportionTestScCODA,
    propeller = RunProportionTestPropeller
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
#' @param backend Permutation backend. `"cpp"` is the default and uses a native
#' permutation/bootstrap loop. `"r"` uses the original R implementation. The C++
#' backend is statistically equivalent but does not reproduce the exact random
#' draws from R's `sample()`.
#'
#' @return A method result bundle used internally by [RunProportionTest].
#'
#' @export
RunProportionTestPermutation <- function(
    srt,
    group.by,
    split.by,
    comparison = NULL,
    n_permutations = 1000,
    include_all_cells = FALSE,
    backend = c("cpp", "r"),
    verbose = TRUE) {
  backend <- match.arg(backend)
  meta_data <- .validate_proportion_inputs(
    srt = srt,
    group.by = group.by,
    split.by = split.by,
    require_sample = FALSE
  )

  comparisons_condition <- .parse_proportion_comparisons(
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

    if (identical(backend, "cpp")) {
      test_result <- .permutation_test_cpp(
        srt = srt,
        group.by = group.by,
        split.by = split.by,
        cluster_1 = cluster_1,
        cluster_2 = cluster_2,
        n_permutations = n_permutations,
        include_all_cells = include_all_cells
      )
    } else {
      test_result <- .permutation_test(
        srt = srt,
        group.by = group.by,
        split.by = split.by,
        cluster_1 = cluster_1,
        cluster_2 = cluster_2,
        n_permutations = n_permutations,
        include_all_cells = include_all_cells
      )
    }

    results_list[[comparison_name]] <- .standardize_proportion_result(
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
      include_all_cells = include_all_cells,
      backend = backend
    )
  )
}

#' @title Milo differential abundance wrapper
#'
#' @description
#' Method-specific implementation used by [RunProportionTest] when
#' `proportion_method = "milo"`.
#' The function always returns a group-level summary and additionally stores a
#' neighborhood-level result list under `neighborhood_results`.
#'
#' @md
#' @inheritParams RunProportionTest
#' @param milo_k Number of nearest neighbors used for Milo graph building.
#' @param milo_d Number of dimensions used by Milo.
#' @param n_bootstrap Number of bootstrap iterations used by fallback summary.
#'
#' @return A method result bundle used internally by [RunProportionTest].
#'
#' @export
RunProportionTestMilo <- function(
    srt,
    group.by,
    split.by,
    sample.by,
    comparison = NULL,
    milo_k = 20L,
    milo_d = 30L,
    n_bootstrap = 500,
    seed = 11,
    verbose = TRUE) {
  meta_data <- .validate_proportion_inputs(
    srt = srt,
    group.by = group.by,
    split.by = split.by,
    sample.by = sample.by,
    require_sample = TRUE
  )

  comparisons_condition <- .parse_proportion_comparisons(
    meta_data = meta_data,
    split.by = split.by,
    comparison = comparison,
    include_bidirectional = TRUE
  )

  milo_neighborhood <- .run_milo_da_with_milor(
    srt = srt,
    group.by = group.by,
    split.by = split.by,
    sample.by = sample.by,
    comparisons_condition = comparisons_condition,
    milo_k = milo_k,
    milo_d = milo_d,
    verbose = verbose
  )
  milo_graph_data <- NULL
  if (!is.null(milo_neighborhood) && is.list(milo_neighborhood)) {
    milo_graph_data <- milo_neighborhood[["graph_data"]]
    milo_neighborhood <- milo_neighborhood[["results"]]
  }
  engine <- if (is.null(milo_neighborhood)) "fallback" else "miloR"

  results_list <- list()
  neighborhood_results <- list()

  for (i in seq_len(nrow(comparisons_condition))) {
    cluster_1 <- comparisons_condition[i, 1]
    cluster_2 <- comparisons_condition[i, 2]
    comparison_name <- paste0(cluster_1, "_vs_", cluster_2)

    group_res <- .sample_level_proportion_test(
      meta_data = meta_data,
      group.by = group.by,
      split.by = split.by,
      sample.by = sample.by,
      cluster_1 = cluster_1,
      cluster_2 = cluster_2,
      n_bootstrap = n_bootstrap,
      transform = "raw",
      seed = seed + i
    )
    group_res <- .standardize_proportion_result(
      group_res,
      cluster_1 = cluster_1,
      cluster_2 = cluster_2,
      comparison_name = comparison_name,
      method = "milo"
    )
    results_list[[comparison_name]] <- group_res

    nhood_res <- milo_neighborhood[[comparison_name]] %||%
      .milo_neighborhood_fallback(group_res)

    nhood_res <- .standardize_proportion_result(
      nhood_res,
      cluster_1 = cluster_1,
      cluster_2 = cluster_2,
      comparison_name = comparison_name,
      method = "milo"
    )

    if (!"neighborhood" %in% colnames(nhood_res)) {
      nhood_res$neighborhood <- paste0("nhood_", seq_len(nrow(nhood_res)))
    }
    neighborhood_results[[comparison_name]] <- nhood_res
  }

  list(
    method = "milo",
    results = results_list,
    neighborhood_results = neighborhood_results,
    result_levels = c("group", "neighborhood"),
    details = list(
      neighborhood_results = neighborhood_results,
      milo_graph_data = milo_graph_data,
      engine = engine
    ),
    parameters = list(
      sample.by = sample.by,
      milo_k = milo_k,
      milo_d = milo_d,
      n_bootstrap = n_bootstrap,
      engine = engine
    )
  )
}

#' @title scCODA differential abundance wrapper
#'
#' @description
#' Method-specific implementation used by [RunProportionTest] when
#' `proportion_method = "sccoda"`.
#' The function calls Python helper `ScCODA` in `inst/python/functions.py`
#' via `reticulate` and falls back to sample-level testing if Python execution
#' is unavailable.
#'
#' @md
#' @inheritParams RunProportionTest
#' @param reference_cell_type Optional reference cell type for scCODA.
#' @param credible_effect_threshold Inclusion probability threshold for
#' credible effects.
#' @param n_mcmc_samples Number of MCMC samples requested in scCODA.
#' @param n_bootstrap Number of bootstrap iterations used by fallback summary.
#'
#' @return A method result bundle used internally by [RunProportionTest].
#'
#' @export
RunProportionTestScCODA <- function(
    srt,
    group.by,
    split.by,
    sample.by,
    comparison = NULL,
    reference_cell_type = NULL,
    credible_effect_threshold = 0.95,
    n_mcmc_samples = 20000L,
    n_bootstrap = 500,
    seed = 11,
    verbose = TRUE) {
  meta_data <- .validate_proportion_inputs(
    srt = srt,
    group.by = group.by,
    split.by = split.by,
    sample.by = sample.by,
    require_sample = TRUE
  )

  comparisons_condition <- .parse_proportion_comparisons(
    meta_data = meta_data,
    split.by = split.by,
    comparison = comparison,
    include_bidirectional = TRUE
  )
  comparison_names <- apply(comparisons_condition, 1, function(x) {
    paste0(x[1], "_vs_", x[2])
  })

  composition_data <- .build_composition_input(
    meta_data = meta_data,
    group.by = group.by,
    split.by = split.by,
    sample.by = sample.by
  )

  py_output <- .run_sccoda_python(
    counts = composition_data$counts,
    metadata = composition_data$metadata,
    comparison_names = comparison_names,
    reference_cell_type = reference_cell_type,
    credible_effect_threshold = credible_effect_threshold,
    n_mcmc_samples = n_mcmc_samples,
    seed = seed,
    verbose = verbose
  )

  results_list <- list()
  for (i in seq_len(nrow(comparisons_condition))) {
    cluster_1 <- comparisons_condition[i, 1]
    cluster_2 <- comparisons_condition[i, 2]
    comparison_name <- paste0(cluster_1, "_vs_", cluster_2)

    sccoda_res <- .extract_sccoda_comparison(py_output, comparison_name)
    if (is.null(sccoda_res)) {
      sccoda_res <- .sample_level_proportion_test(
        meta_data = meta_data,
        group.by = group.by,
        split.by = split.by,
        sample.by = sample.by,
        cluster_1 = cluster_1,
        cluster_2 = cluster_2,
        n_bootstrap = n_bootstrap,
        transform = "logit",
        seed = seed + i
      )
      sccoda_res$inclusion_prob <- pmax(0, pmin(1, 1 - sccoda_res$pval))
      sccoda_res$credible <- sccoda_res$inclusion_prob >= credible_effect_threshold
    }

    std <- .standardize_proportion_result(
      sccoda_res,
      cluster_1 = cluster_1,
      cluster_2 = cluster_2,
      comparison_name = comparison_name,
      method = "sccoda"
    )

    if (!"inclusion_prob" %in% colnames(std)) {
      std$inclusion_prob <- pmax(0, pmin(1, 1 - std$pval))
    }
    if (!"credible" %in% colnames(std)) {
      std$credible <- std$inclusion_prob >= credible_effect_threshold
    }

    results_list[[comparison_name]] <- std
  }

  list(
    method = "sccoda",
    results = results_list,
    result_levels = "group",
    details = list(
      python = py_output,
      reference_cell_type = reference_cell_type,
      credible_effect_threshold = credible_effect_threshold
    ),
    parameters = list(
      sample.by = sample.by,
      reference_cell_type = reference_cell_type,
      credible_effect_threshold = credible_effect_threshold,
      n_mcmc_samples = n_mcmc_samples
    )
  )
}

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
RunProportionTestPropeller <- function(
    srt,
    group.by,
    split.by,
    sample.by,
    comparison = NULL,
    n_bootstrap = 1000,
    seed = 11,
    verbose = TRUE) {
  meta_data <- .validate_proportion_inputs(
    srt = srt,
    group.by = group.by,
    split.by = split.by,
    sample.by = sample.by,
    require_sample = TRUE
  )

  comparisons_condition <- .parse_proportion_comparisons(
    meta_data = meta_data,
    split.by = split.by,
    comparison = comparison,
    include_bidirectional = TRUE
  )

  engine <- if (requireNamespace("speckle", quietly = TRUE)) {
    "speckle-compatible"
  } else {
    "internal"
  }

  results_list <- list()
  for (i in seq_len(nrow(comparisons_condition))) {
    cluster_1 <- comparisons_condition[i, 1]
    cluster_2 <- comparisons_condition[i, 2]
    comparison_name <- paste0(cluster_1, "_vs_", cluster_2)

    prop_res <- .sample_level_proportion_test(
      meta_data = meta_data,
      group.by = group.by,
      split.by = split.by,
      sample.by = sample.by,
      cluster_1 = cluster_1,
      cluster_2 = cluster_2,
      n_bootstrap = n_bootstrap,
      transform = "logit",
      seed = seed + i
    )

    results_list[[comparison_name]] <- .standardize_proportion_result(
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

.normalize_proportion_method <- function(proportion_method) {
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
    permutation_test = "permutation"
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

.validate_proportion_inputs <- function(
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

.parse_proportion_comparisons <- function(
    meta_data,
    split.by,
    comparison = NULL,
    include_bidirectional = TRUE
) {
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
    base_pairs <- .coerce_comparison_pairs(comparison)

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

.coerce_comparison_pairs <- function(comparison) {
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
    if (all(vapply(comparison, function(x) is.character(x) && length(x) == 2, logical(1)))) {
      out <- do.call(rbind, lapply(comparison, function(x) trimws(x)))
    } else {
      comp_vec <- unlist(comparison)
      out <- do.call(rbind, lapply(comp_vec, parse_vs))
    }
  } else if (is.character(comparison)) {
    if (length(comparison) == 2 && !any(grepl("_vs_", comparison, fixed = TRUE))) {
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

.sample_level_proportion_test <- function(
    meta_data,
    group.by,
    split.by,
    sample.by,
    cluster_1,
    cluster_2,
    n_bootstrap = 1000,
    transform = c("raw", "logit", "asin"),
    pseudocount = 1e-5,
    seed = NULL
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

    # Bootstrap computation - use C++ when available
    if (n_bootstrap > 0 && length(v1) > 0 && length(v2) > 0) {
      if (run_proportion_bootstrap_stats_cpp_available()) {
        boot_result <- run_proportion_bootstrap_stats_cpp(
          v1 = v1,
          v2 = v2,
          n_bootstrap = n_bootstrap,
          pseudocount = pseudocount
        )
        boot_mean_log2FD <- boot_result[["boot_mean_log2FD"]]
        boot_CI_2.5 <- boot_result[["boot_CI_2.5"]]
        boot_CI_97.5 <- boot_result[["boot_CI_97.5"]]
      } else {
        boot <- rep(NA_real_, n_bootstrap)
        for (b in seq_len(n_bootstrap)) {
          b1 <- sample(v1, size = length(v1), replace = TRUE)
          b2 <- sample(v2, size = length(v2), replace = TRUE)
          boot[b] <- log2((mean(b2, na.rm = TRUE) + pseudocount) /
            (mean(b1, na.rm = TRUE) + pseudocount))
        }
        boot_mean_log2FD <- mean(boot, na.rm = TRUE)
        boot_CI_2.5 <- stats::quantile(boot, probs = 0.025, na.rm = TRUE, names = FALSE)
        boot_CI_97.5 <- stats::quantile(boot, probs = 0.975, na.rm = TRUE, names = FALSE)
      }
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

.standardize_proportion_result <- function(
    x,
    cluster_1,
    cluster_2,
    comparison_name,
    method
) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)

  cluster_col <- .find_first_col(
    df,
    c("clusters", "cluster", "cell_type", "CellType", "celltype", "feature", "name")
  )
  if (is.null(cluster_col)) {
    if (!is.null(rownames(df)) && !all(rownames(df) %in% c("", as.character(seq_len(nrow(df)))))) {
      df$clusters <- rownames(df)
    } else {
      df$clusters <- paste0("cluster_", seq_len(nrow(df)))
    }
  } else if (cluster_col != "clusters") {
    df$clusters <- as.character(df[[cluster_col]])
  }

  obs_col <- .find_first_col(df, c("obs_log2FD", "log2FD", "logFC", "coef", "effect"))
  if (is.null(obs_col)) {
    df$obs_log2FD <- NA_real_
  } else {
    df$obs_log2FD <- suppressWarnings(as.numeric(df[[obs_col]]))
  }

  pval_col <- .find_first_col(df, c("pval", "p_val", "p.value", "PValue", "P.Value"))
  if (is.null(pval_col)) {
    df$pval <- NA_real_
  } else {
    df$pval <- suppressWarnings(as.numeric(df[[pval_col]]))
  }

  fdr_col <- .find_first_col(df, c("FDR", "adj.P.Val", "SpatialFDR", "q_value", "qval"))
  if (is.null(fdr_col)) {
    df$FDR <- stats::p.adjust(df$pval, method = "fdr")
  } else {
    df$FDR <- suppressWarnings(as.numeric(df[[fdr_col]]))
  }

  boot_mean_col <- .find_first_col(df, c("boot_mean_log2FD", "boot_mean", "mean_log2FD"))
  if (is.null(boot_mean_col)) {
    df$boot_mean_log2FD <- NA_real_
  } else {
    df$boot_mean_log2FD <- suppressWarnings(as.numeric(df[[boot_mean_col]]))
  }

  ci_low_col <- .find_first_col(df, c("boot_CI_2.5", "CI_2.5", "ci_low", "lower", "hdi_2.5"))
  if (is.null(ci_low_col)) {
    df$boot_CI_2.5 <- NA_real_
  } else {
    df$boot_CI_2.5 <- suppressWarnings(as.numeric(df[[ci_low_col]]))
  }

  ci_high_col <- .find_first_col(df, c("boot_CI_97.5", "CI_97.5", "ci_high", "upper", "hdi_97.5"))
  if (is.null(ci_high_col)) {
    df$boot_CI_97.5 <- NA_real_
  } else {
    df$boot_CI_97.5 <- suppressWarnings(as.numeric(df[[ci_high_col]]))
  }

  effect_col <- .find_first_col(df, c("effect", "coef", "logFC", "log2FD", "obs_log2FD"))
  if (is.null(effect_col)) {
    df$effect <- df$obs_log2FD
  } else {
    df$effect <- suppressWarnings(as.numeric(df[[effect_col]]))
  }

  ip_col <- .find_first_col(
    df,
    c("inclusion_prob", "inclusion_probability", "probability", "pip")
  )
  if (is.null(ip_col)) {
    df$inclusion_prob <- NA_real_
  } else {
    df$inclusion_prob <- suppressWarnings(as.numeric(df[[ip_col]]))
  }

  credible_col <- .find_first_col(df, c("credible", "is_credible", "is_significant"))
  if (is.null(credible_col)) {
    df$credible <- !is.na(df$inclusion_prob) & df$inclusion_prob >= 0.95
  } else {
    df$credible <- as.logical(df[[credible_col]])
  }

  hdi_low_col <- .find_first_col(df, c("hdi_2.5", "hdi_low", "hdi_lower", "ci_low", "boot_CI_2.5"))
  if (is.null(hdi_low_col)) {
    df$hdi_2.5 <- suppressWarnings(as.numeric(df$boot_CI_2.5))
  } else {
    df$hdi_2.5 <- suppressWarnings(as.numeric(df[[hdi_low_col]]))
  }

  hdi_high_col <- .find_first_col(df, c("hdi_97.5", "hdi_high", "hdi_upper", "ci_high", "boot_CI_97.5"))
  if (is.null(hdi_high_col)) {
    df$hdi_97.5 <- suppressWarnings(as.numeric(df$boot_CI_97.5))
  } else {
    df$hdi_97.5 <- suppressWarnings(as.numeric(df[[hdi_high_col]]))
  }

  nhood_col <- .find_first_col(df, c("neighborhood", "Nhood", "nhood", "nhood_index"))
  if (!is.null(nhood_col) && nhood_col != "neighborhood") {
    df$neighborhood <- as.character(df[[nhood_col]])
  }

  df$group1 <- cluster_1
  df$group2 <- cluster_2
  df$comparison <- comparison_name
  df$method <- method

  core_cols <- c(
    "clusters", "obs_log2FD", "pval", "FDR",
    "boot_mean_log2FD", "boot_CI_2.5", "boot_CI_97.5",
    "group1", "group2", "comparison", "method"
  )
  for (nm in core_cols) {
    if (!nm %in% colnames(df)) {
      df[[nm]] <- NA
    }
  }

  df <- df[, unique(c(core_cols, setdiff(colnames(df), core_cols))), drop = FALSE]
  rownames(df) <- NULL
  df
}

.find_first_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    return(NULL)
  }
  hit[1]
}

.milo_neighborhood_fallback <- function(group_result) {
  out <- group_result
  out$neighborhood <- paste0("nhood_", make.names(out$clusters))
  out
}

.run_milo_da_with_milor <- function(
    srt,
    group.by,
    split.by,
    sample.by,
    comparisons_condition,
    milo_k = 20L,
    milo_d = 30L,
    verbose = TRUE
) {
  if (!requireNamespace("miloR", quietly = TRUE)) {
    return(NULL)
  }

  tryCatch(
    {
      sce <- Seurat::as.SingleCellExperiment(srt)
      reduced_name <- NULL

      if ("PCA" %in% SingleCellExperiment::reducedDimNames(sce)) {
        reduced_name <- "PCA"
      } else if ("pca" %in% names(srt@reductions)) {
        SingleCellExperiment::reducedDim(sce, "PCA") <- SeuratObject::Embeddings(srt, "pca")
        reduced_name <- "PCA"
      }

      if (is.null(reduced_name)) {
        log_message(
          "No PCA reduction available for {.pkg miloR}; use fallback Milo summary",
          message_type = "warning",
          verbose = verbose
        )
        return(NULL)
      }

      milo_obj <- miloR::Milo(sce)
      n_dim <- min(milo_d, ncol(SingleCellExperiment::reducedDim(sce, reduced_name)))
      milo_obj <- miloR::buildGraph(
        milo_obj,
        k = as.integer(milo_k),
        d = as.integer(n_dim),
        reduced.dim = reduced_name
      )
      milo_obj <- miloR::makeNhoods(
        milo_obj,
        k = as.integer(milo_k),
        d = as.integer(n_dim),
        refined = TRUE,
        reduced_dims = reduced_name
      )

      cdata <- as.data.frame(SummarizedExperiment::colData(milo_obj))
      cdata[[sample.by]] <- as.character(cdata[[sample.by]])
      cdata[[split.by]] <- as.character(cdata[[split.by]])
      cdata[[group.by]] <- as.character(cdata[[group.by]])

      milo_obj <- miloR::countCells(
        milo_obj,
        meta.data = cdata,
        sample = sample.by
      )

      design_df <- unique(cdata[, c(sample.by, split.by), drop = FALSE])
      rownames(design_df) <- design_df[[sample.by]]
      design_df[[split.by]] <- as.factor(design_df[[split.by]])

      output <- list()
      graph_data <- list()
      for (i in seq_len(nrow(comparisons_condition))) {
        cluster_1 <- comparisons_condition[i, 1]
        cluster_2 <- comparisons_condition[i, 2]
        comparison_name <- paste0(cluster_1, "_vs_", cluster_2)

        keep <- design_df[[split.by]] %in% c(cluster_1, cluster_2)
        dsub <- droplevels(design_df[keep, , drop = FALSE])
        if (length(unique(dsub[[split.by]])) < 2) {
          next
        }

        design_formula <- stats::as.formula(paste0("~", split.by))
        contrast <- paste0(split.by, cluster_2, "-", split.by, cluster_1)

        da <- miloR::testNhoods(
          milo_obj,
          design = design_formula,
          design.df = dsub,
          model.contrasts = contrast
        )

        da$neighborhood <- if ("Nhood" %in% colnames(da)) {
          paste0("nhood_", da$Nhood)
        } else {
          paste0("nhood_", seq_len(nrow(da)))
        }

        if (!"clusters" %in% colnames(da)) {
          da$clusters <- da$neighborhood
        }

        output[[comparison_name]] <- da

        node_df <- data.frame(
          neighborhood = as.character(da$neighborhood),
          clusters = as.character(da$clusters),
          stringsAsFactors = FALSE
        )
        if (nrow(node_df) > 0) {
          theta <- seq(0, 2 * pi, length.out = nrow(node_df) + 1)[seq_len(nrow(node_df))]
          node_df$x <- cos(theta)
          node_df$y <- sin(theta)

          edge_df <- data.frame(
            from = node_df$neighborhood,
            to = c(node_df$neighborhood[-1], node_df$neighborhood[1]),
            stringsAsFactors = FALSE
          )
        } else {
          node_df$x <- numeric(0)
          node_df$y <- numeric(0)
          edge_df <- data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE)
        }
        graph_data[[comparison_name]] <- list(nodes = node_df, edges = edge_df)
      }
      list(results = output, graph_data = graph_data)
    },
    error = function(e) {
      log_message(
        "{.pkg miloR} execution failed, fallback to internal Milo summary: {.val {e$message}}",
        message_type = "warning",
        verbose = verbose
      )
      NULL
    }
  )
}

.build_composition_input <- function(meta_data, group.by, split.by, sample.by) {
  dat <- meta_data[, c(group.by, split.by, sample.by), drop = FALSE]
  colnames(dat) <- c("cluster", "condition", "sample")
  dat$cluster <- as.character(dat$cluster)
  dat$condition <- as.character(dat$condition)
  dat$sample <- as.character(dat$sample)

  sample_meta <- stats::aggregate(
    condition ~ sample,
    data = dat,
    FUN = function(x) {
      ux <- unique(x)
      ux[1]
    }
  )
  sample_meta <- sample_meta[order(sample_meta$sample), , drop = FALSE]
  rownames(sample_meta) <- sample_meta$sample

  sample_levels <- sample_meta$sample
  cluster_levels <- sort(unique(dat$cluster))

  count_tab <- stats::xtabs(~ sample + cluster, data = dat)
  count_mat <- matrix(
    0,
    nrow = length(sample_levels),
    ncol = length(cluster_levels),
    dimnames = list(sample_levels, cluster_levels)
  )
  count_mat[rownames(count_tab), colnames(count_tab)] <- as.matrix(count_tab)

  list(
    counts = as.data.frame(count_mat, stringsAsFactors = FALSE),
    metadata = data.frame(
      sample = sample_levels,
      condition = sample_meta[sample_levels, "condition"],
      stringsAsFactors = FALSE
    )
  )
}

.run_sccoda_python <- function(
    counts,
    metadata,
    comparison_names,
    reference_cell_type = NULL,
    credible_effect_threshold = 0.95,
    n_mcmc_samples = 20000L,
    seed = 11,
    verbose = TRUE
) {
  tryCatch(
    {
      PrepareEnv()
      check_python("sccoda", verbose = FALSE)

      functions <- reticulate::import_from_path(
        "functions",
        path = system.file("python", package = "scop", mustWork = TRUE),
        convert = TRUE
      )

      functions$ScCODA(
        counts = counts,
        metadata = metadata,
        condition_key = "condition",
        sample_key = "sample",
        comparisons = as.list(comparison_names),
        reference_cell_type = reference_cell_type %||% "",
        credible_effect_threshold = as.double(credible_effect_threshold),
        random_seed = as.integer(seed),
        mcmc_samples = as.integer(n_mcmc_samples),
        verbose = verbose
      )
    },
    error = function(e) {
      log_message(
        "scCODA python execution unavailable, using fallback summary: {.val {e$message}}",
        message_type = "warning",
        verbose = verbose
      )
      NULL
    }
  )
}

.extract_sccoda_comparison <- function(py_output, comparison_name) {
  if (is.null(py_output) || is.null(py_output[["results"]])) {
    return(NULL)
  }

  comp <- py_output[["results"]][[comparison_name]]
  if (is.null(comp)) {
    return(NULL)
  }

  if (is.data.frame(comp)) {
    return(comp)
  }

  if (is.list(comp) && length(comp) > 0) {
    comp_df <- tryCatch(
      {
        as.data.frame(comp, stringsAsFactors = FALSE)
      },
      error = function(e) {
        tryCatch(
          {
            do.call(
              rbind,
              lapply(comp, function(x) {
                as.data.frame(x, stringsAsFactors = FALSE)
              })
            )
          },
          error = function(e2) NULL
        )
      }
    )
    return(comp_df)
  }

  NULL
}

.permutation_test <- function(
    srt,
    group.by,
    split.by,
    cluster_1,
    cluster_2,
    n_permutations,
    include_all_cells = FALSE,
    pseudocount = 1e-8) {
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

  obs_diff_wide$obs_log2FD <- log2(
    (obs_diff_wide[[cluster_2]] + pseudocount) /
      (obs_diff_wide[[cluster_1]] + pseudocount)
  )

  perm_results <- matrix(NA_real_, nrow(obs_diff_wide), n_permutations)
  rownames(perm_results) <- obs_diff_wide$clusters

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
      (permuted_wide[[cluster_2]] + pseudocount) /
        (permuted_wide[[cluster_1]] + pseudocount)
    )

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

  boot_results <- matrix(NA_real_, nrow(obs_diff_wide), n_permutations)
  rownames(boot_results) <- obs_diff_wide$clusters

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
      (booted_wide[[cluster_2]] + pseudocount) /
        (booted_wide[[cluster_1]] + pseudocount)
    )

    boot_results[, i] <- booted_wide$boot_log2FD
  }

  boot_results[!is.finite(boot_results)] <- NA
  obs_diff_wide$boot_mean_log2FD <- rowMeans(boot_results, na.rm = TRUE)

  boot_ci <- t(
    apply(
      boot_results, 1, function(x) {
        stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
      }
    )
  )
  boot_ci <- as.data.frame(boot_ci)
  colnames(boot_ci) <- c("boot_CI_2.5", "boot_CI_97.5")

  cbind(obs_diff_wide, boot_ci)
}

.permutation_test_cpp_available <- function() {
  exists("proportion_permutation_cpp", mode = "function") &&
    isTRUE(is.loaded("_scop_proportion_permutation_cpp"))
}

.permutation_test_cpp <- function(
    srt,
    group.by,
    split.by,
    cluster_1,
    cluster_2,
    n_permutations,
    include_all_cells = FALSE,
    pseudocount = 1e-8) {
  if (!.permutation_test_cpp_available()) {
    log_message(
      "{.arg backend = 'cpp'} requires the compiled {.pkg scop} shared library. Reinstall the package to build native code.",
      message_type = "error"
    )
  }

  meta_data <- srt@meta.data[, c(split.by, group.by), drop = FALSE]
  colnames(meta_data) <- c("samples", "clusters")
  meta_data[["clusters"]] <- as.character(meta_data[["clusters"]])
  comparison_data <- meta_data[
    meta_data[["samples"]] %in% c(cluster_1, cluster_2), ,
    drop = FALSE
  ]
  if (isTRUE(include_all_cells)) {
    cluster_cases <- unique(meta_data[["clusters"]])
  } else {
    cluster_cases <- unique(comparison_data[["clusters"]])
  }

  sample_ids <- ifelse(comparison_data[["samples"]] == cluster_1, 1L, 2L)
  cluster_ids <- match(comparison_data[["clusters"]], cluster_cases)
  res <- proportion_permutation_cpp(
    sample_ids = as.integer(sample_ids),
    cluster_ids = as.integer(cluster_ids),
    cluster_levels = as.character(cluster_cases),
    n_permutations = as.integer(n_permutations),
    pseudocount = pseudocount
  )
  colnames(res)[colnames(res) == "fraction_1"] <- cluster_1
  colnames(res)[colnames(res) == "fraction_2"] <- cluster_2
  res[["FDR"]] <- stats::p.adjust(res[["pval"]], "fdr")
  res[, c(
    "clusters", cluster_1, cluster_2, "obs_log2FD", "pval", "FDR",
    "boot_mean_log2FD", "boot_CI_2.5", "boot_CI_97.5"
  ), drop = FALSE]
}
