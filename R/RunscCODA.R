
#' @title scCODA differential abundance
#'
#' @md
#' @inheritParams RunProportionTest
#' @param reference_cell_type Optional reference cell type for scCODA.
#' @param credible_effect_threshold Inclusion probability threshold for
#' credible effects.
#' @param n_mcmc_samples Number of MCMC samples requested in scCODA.
#' @param envname Name of the conda-compatible environment used by scCODA.
#' Defaults to `"scop_sccoda_env"` to keep the TensorFlow/scCODA stack isolated
#' from the default Python environment.
#' @param conda The path or command name of a conda-compatible executable.
#'
#' @return A method result bundle used internally by [RunProportionTest].
#'
#' @export
RunscCODA <- function(
  srt,
  group.by,
  split.by,
  sample.by,
  comparison = NULL,
  reference_cell_type = NULL,
  credible_effect_threshold = 0.95,
  n_mcmc_samples = 20000L,
  envname = "scop_sccoda_env",
  conda = "auto",
  seed = 11,
  verbose = TRUE
) {
  check_r("reticulate", verbose = FALSE)
  PrepareEnv(envname = envname, conda = conda, modules = "sccoda")
  conda <- resolve_conda(conda)
  envname <- get_envname(envname)
  env_path <- conda_env_path(envname = envname, conda = conda)
  python_path <- if (!is.null(env_path)) {
    conda_env_python_path(env_path)
  } else {
    NULL
  }
  if (is.null(python_path) || !file.exists(python_path)) {
    log_message(
      "Unable to locate the Python executable for {.file {envname}}",
      message_type = "error"
    )
  }
  sccoda_preflight <- tryCatch(
    suppressWarnings(system2(
      python_path,
      c(
        "-X",
        "faulthandler",
        "-c",
        shQuote(
          "import sccoda.util.comp_ana; import sccoda.util.cell_composition_data"
        )
      ),
      stdout = TRUE,
      stderr = TRUE
    )),
    error = function(e) {
      structure(e$message, status = 1L)
    }
  )
  sccoda_status <- attr(sccoda_preflight, "status") %||% 0L
  if (!identical(as.integer(sccoda_status), 0L)) {
    detail <- paste(
      utils::tail(as.character(sccoda_preflight), 8),
      collapse = " | "
    )
    log_message(
      "{.pkg scCODA} import preflight failed in {.file {envname}}: {.val {detail}}",
      message_type = "error"
    )
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
  comparison_names <- apply(comparisons_condition, 1, function(x) {
    paste0(x[1], "_vs_", x[2])
  })

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

  functions <- reticulate::import_from_path(
    "functions",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )
  py_output <- functions$ScCODA(
    counts = as.data.frame(count_mat, stringsAsFactors = FALSE),
    metadata = data.frame(
      sample = sample_levels,
      condition = sample_meta[sample_levels, "condition"],
      stringsAsFactors = FALSE
    ),
    condition_key = "condition",
    sample_key = "sample",
    comparisons = as.list(comparison_names),
    reference_cell_type = reference_cell_type %||% "",
    credible_effect_threshold = as.double(credible_effect_threshold),
    random_seed = as.integer(seed),
    mcmc_samples = as.integer(n_mcmc_samples),
    verbose = verbose
  )

  results_list <- list()
  for (i in seq_len(nrow(comparisons_condition))) {
    cluster_1 <- comparisons_condition[i, 1]
    cluster_2 <- comparisons_condition[i, 2]
    comparison_name <- paste0(cluster_1, "_vs_", cluster_2)

    sccoda_res <- NULL
    if (!is.null(py_output) && !is.null(py_output[["results"]])) {
      comp <- py_output[["results"]][[comparison_name]]
      if (is.data.frame(comp)) {
        sccoda_res <- comp
      } else if (is.list(comp) && length(comp) > 0) {
        sccoda_res <- tryCatch(
          {
            if (all(vapply(comp, is.list, logical(1)))) {
              do.call(
                rbind,
                lapply(comp, function(x) {
                  as.data.frame(x, stringsAsFactors = FALSE)
                })
              )
            } else {
              as.data.frame(comp, stringsAsFactors = FALSE)
            }
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
      }
    }
    if (is.null(sccoda_res)) {
      log_message(
        "scCODA did not return results for comparison {.val {comparison_name}}",
        message_type = "error"
      )
    }

    std <- standardize_proportion_result(
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
