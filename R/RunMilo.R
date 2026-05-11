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
#' @param n_bootstrap Number of bootstrap iterations used by the group-level summary.
#' @return A method result bundle used internally by [RunProportionTest].
#'
#' @export
RunMilo <- function(
  srt,
  group.by,
  split.by,
  sample.by,
  comparison = NULL,
  milo_k = 20L,
  milo_d = 30L,
  n_bootstrap = 500,
  seed = 11,
  verbose = TRUE
) {
  check_r(
    c("miloR", "SingleCellExperiment", "SummarizedExperiment"),
    verbose = FALSE
  )

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

  milo_neighborhood <- run_milo_da(
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
  engine <- "miloR"

  results_list <- list()
  neighborhood_results <- list()

  for (i in seq_len(nrow(comparisons_condition))) {
    cluster_1 <- comparisons_condition[i, 1]
    cluster_2 <- comparisons_condition[i, 2]
    comparison_name <- paste0(cluster_1, "_vs_", cluster_2)

    group_res <- sample_level_proportion_test(
      meta_data = meta_data,
      group.by = group.by,
      split.by = split.by,
      sample.by = sample.by,
      cluster_1 = cluster_1,
      cluster_2 = cluster_2,
      n_bootstrap = n_bootstrap,
      transform = "raw",
      seed = seed + i,
      verbose = verbose
    )
    group_res <- standardize_proportion_result(
      group_res,
      cluster_1 = cluster_1,
      cluster_2 = cluster_2,
      comparison_name = comparison_name,
      method = "milo"
    )
    results_list[[comparison_name]] <- group_res

    nhood_res <- milo_neighborhood[[comparison_name]]
    if (is.null(nhood_res)) {
      log_message(
        "Milo did not return neighborhood results for comparison {.val {comparison_name}}",
        message_type = "error"
      )
    }

    nhood_res <- standardize_proportion_result(
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

run_milo_da <- function(
  srt,
  group.by,
  split.by,
  sample.by,
  comparisons_condition,
  milo_k = 20L,
  milo_d = 30L,
  verbose = TRUE
) {
  tryCatch(
    {
      sce <- Seurat::as.SingleCellExperiment(srt)
      reduced_name <- NULL

      if ("PCA" %in% SingleCellExperiment::reducedDimNames(sce)) {
        reduced_name <- "PCA"
      } else {
        pca_reduction <- tryCatch(
          DefaultReduction(srt, pattern = "pca"),
          error = function(e) NULL
        )
        if (
          !is.null(pca_reduction) && pca_reduction %in% names(srt@reductions)
        ) {
          pca_embedding <- SeuratObject::Embeddings(srt, pca_reduction)
          pca_embedding <- pca_embedding[colnames(sce), , drop = FALSE]
          SingleCellExperiment::reducedDim(sce, "PCA") <- pca_embedding
          reduced_name <- "PCA"
        }
      }

      if (is.null(reduced_name)) {
        log_message(
          "No PCA reduction available for {.pkg miloR}; run PCA before {.fn RunMilo}",
          message_type = "error"
        )
      }

      milo_obj <- miloR::Milo(sce)
      n_dim <- min(
        milo_d,
        ncol(SingleCellExperiment::reducedDim(sce, reduced_name))
      )
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

        design_formula <- stats::as.formula(paste0("~ 0 + ", split.by))
        design_cols <- colnames(stats::model.matrix(design_formula, dsub))
        contrast_cols <- stats::setNames(design_cols, levels(dsub[[split.by]]))
        contrast <- paste0(
          contrast_cols[[cluster_2]],
          "-",
          contrast_cols[[cluster_1]]
        )

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
          theta <- seq(0, 2 * pi, length.out = nrow(node_df) + 1)[seq_len(nrow(
            node_df
          ))]
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
          edge_df <- data.frame(
            from = character(0),
            to = character(0),
            stringsAsFactors = FALSE
          )
        }
        graph_data[[comparison_name]] <- list(nodes = node_df, edges = edge_df)
      }
      list(results = output, graph_data = graph_data)
    },
    error = function(e) {
      log_message(
        "{.pkg miloR} execution failed: {.val {e$message}}",
        message_type = "error"
      )
    }
  )
}
