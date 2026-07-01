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
#' @param reduction Dimensional reduction used for Milo graph construction.
#' If `NULL`, the default PCA-like reduction is used.
#' @param backend Backend used to compute Milo neighborhoods and tests.
#' `"r"` calls the upstream `miloR` implementation. `"cpp"` uses native
#' `scop` graph/neighborhood/counting kernels and keeps edgeR for the
#' neighborhood-level quasi-likelihood test.
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
  reduction = NULL,
  backend = c("r", "cpp"),
  n_bootstrap = 500,
  seed = 11,
  verbose = TRUE
) {
  backend <- match.arg(backend)
  if (identical(backend, "r")) {
    check_r(
      c("miloR", "SingleCellExperiment", "SummarizedExperiment"),
      verbose = FALSE
    )
  } else {
    check_r(c("edgeR", "limma", "BiocNeighbors"), verbose = FALSE)
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

  milo_neighborhood <- if (identical(backend, "r")) {
    run_milo_da_r(
      srt = srt,
      group.by = group.by,
      split.by = split.by,
      sample.by = sample.by,
      comparisons_condition = comparisons_condition,
      milo_k = milo_k,
      milo_d = milo_d,
      reduction = reduction,
      seed = seed,
      verbose = verbose
    )
  } else {
    run_milo_da_cpp(
      srt = srt,
      group.by = group.by,
      split.by = split.by,
      sample.by = sample.by,
      comparisons_condition = comparisons_condition,
      milo_k = milo_k,
      milo_d = milo_d,
      reduction = reduction,
      seed = seed,
      verbose = verbose
    )
  }
  milo_graph_data <- NULL
  if (!is.null(milo_neighborhood) && is.list(milo_neighborhood)) {
    milo_graph_data <- milo_neighborhood[["graph_data"]]
    milo_neighborhood <- milo_neighborhood[["results"]]
  }

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
      backend = backend
    ),
    parameters = list(
      sample.by = sample.by,
      milo_k = milo_k,
      milo_d = milo_d,
      reduction = reduction,
      n_bootstrap = n_bootstrap,
      backend = backend
    )
  )
}

run_milo_da_r <- function(
  srt,
  group.by,
  split.by,
  sample.by,
  comparisons_condition,
  milo_k = 20L,
  milo_d = 30L,
  reduction = NULL,
  seed = 11,
  verbose = TRUE
) {
  tryCatch(
    {
      sce <- Seurat::as.SingleCellExperiment(srt)
      reduced_name <- NULL

      if (!is.null(reduction) && reduction %in% names(srt@reductions)) {
        pca_embedding <- SeuratObject::Embeddings(srt, reduction)
        pca_embedding <- pca_embedding[colnames(sce), , drop = FALSE]
        SingleCellExperiment::reducedDim(sce, "PCA") <- pca_embedding
        reduced_name <- "PCA"
      } else if ("PCA" %in% SingleCellExperiment::reducedDimNames(sce)) {
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
      set.seed(seed)
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
      safe_split_by <- ".scop_milo_condition"
      condition_map <- setNames(
        make.names(unique(cdata[[split.by]]), unique = TRUE),
        unique(cdata[[split.by]])
      )
      cdata[[safe_split_by]] <- unname(condition_map[cdata[[split.by]]])

      milo_obj <- miloR::countCells(
        milo_obj,
        meta.data = cdata,
        sample = sample.by
      )

      nhood_metadata <- extract_milo_r_neighborhood_metadata(milo_obj)

      design_df <- unique(cdata[, c(sample.by, safe_split_by), drop = FALSE])
      rownames(design_df) <- design_df[[sample.by]]
      design_df[[safe_split_by]] <- as.factor(design_df[[safe_split_by]])

      output <- list()
      graph_data <- list()
      for (i in seq_len(nrow(comparisons_condition))) {
        cluster_1 <- comparisons_condition[i, 1]
        cluster_2 <- comparisons_condition[i, 2]
        comparison_name <- paste0(cluster_1, "_vs_", cluster_2)

        keep <- design_df[[safe_split_by]] %in% unname(condition_map[c(cluster_1, cluster_2)])
        dsub <- droplevels(design_df[keep, , drop = FALSE])
        if (length(unique(dsub[[safe_split_by]])) < 2) {
          next
        }

        design_formula <- stats::as.formula(paste0("~ 0 + ", safe_split_by))
        design_cols <- colnames(stats::model.matrix(design_formula, dsub))
        contrast_cols <- stats::setNames(design_cols, levels(dsub[[safe_split_by]]))
        contrast <- paste0(
          contrast_cols[[condition_map[[cluster_2]]]],
          "-",
          contrast_cols[[condition_map[[cluster_1]]]]
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
      graph_data[[".metadata"]] <- nhood_metadata
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

run_milo_da_cpp <- function(
  srt,
  group.by,
  split.by,
  sample.by,
  comparisons_condition,
  milo_k = 20L,
  milo_d = 30L,
  reduction = NULL,
  seed = 11,
  verbose = TRUE
) {
  tryCatch(
    {
      pca_reduction <- reduction %||% tryCatch(
        DefaultReduction(srt, pattern = "pca"),
        error = function(e) NULL
      )
      if (is.null(pca_reduction) || !pca_reduction %in% names(srt@reductions)) {
        log_message(
          "No PCA reduction available for {.arg backend = 'cpp'}; run PCA before {.fn RunMilo}",
          message_type = "error"
        )
      }

      coords <- SeuratObject::Embeddings(srt, pca_reduction)
      n_dim <- min(as.integer(milo_d), ncol(coords))
      coords <- as.matrix(coords[, seq_len(n_dim), drop = FALSE])
      storage.mode(coords) <- "double"

      n_cells <- nrow(coords)
      if (n_cells < 2L) {
        log_message(
          "{.fn RunMilo} requires at least two cells for {.arg backend = 'cpp'}",
          message_type = "error"
        )
      }
      knn_k <- max(1L, min(as.integer(milo_k), n_cells - 1L))
      log_message(
        "Computing {.pkg Milo} cpp KNN graph...",
        verbose = verbose
      )
      knn_raw <- BiocNeighbors::findKNN(
        coords,
        k = knn_k,
        get.index = TRUE,
        get.distance = TRUE,
        BNPARAM = BiocNeighbors::KmknnParam()
      )
      knn <- list(idx = knn_raw[["index"]], dist = knn_raw[["distance"]])

      set.seed(seed)
      n_seeds <- max(1L, floor(0.1 * n_cells))
      random_vertices <- sample.int(n_cells, n_seeds)
      refined <- milo_refined_vertices(
        random_vertices = random_vertices,
        coords = coords,
        k = knn_k
      )
      sampled_vertices <- unique(as.integer(refined))

      meta_data <- srt@meta.data
      meta_data[[sample.by]] <- as.character(meta_data[[sample.by]])
      meta_data[[split.by]] <- as.character(meta_data[[split.by]])
      meta_data[[group.by]] <- as.character(meta_data[[group.by]])
      safe_split_by <- ".scop_milo_condition"
      condition_map <- setNames(
        make.names(unique(meta_data[[split.by]]), unique = TRUE),
        unique(meta_data[[split.by]])
      )
      meta_data[[safe_split_by]] <- unname(condition_map[meta_data[[split.by]]])
      samples <- unique(meta_data[[sample.by]])

      nhood <- milo_nhood_counts_cpp(
        knn_idx = knn[["idx"]],
        sampled_vertices = sampled_vertices,
        sample_id = match(meta_data[[sample.by]], samples),
        n_samples = length(samples),
        k_dist = knn[["dist"]][, knn_k]
      )
      count_matrix <- nhood[["counts"]]
      colnames(count_matrix) <- samples
      rownames(count_matrix) <- seq_len(nrow(count_matrix))
      nhood_members <- lapply(nhood[["members"]], function(i) {
        rownames(coords)[as.integer(i)]
      })
      names(nhood_members) <- paste0("nhood_", seq_along(nhood_members))
      nhood_metadata <- list(
        sampled_vertices = as.integer(nhood[["sampled_vertices"]]),
        sampled_cells = rownames(coords)[as.integer(nhood[["sampled_vertices"]])],
        members = nhood_members,
        counts = count_matrix,
        sample_names = samples,
        cell_names = rownames(coords),
        k_distance = nhood[["k_distance"]],
        random_vertices = random_vertices
      )

      design_df <- unique(meta_data[, c(sample.by, safe_split_by), drop = FALSE])
      rownames(design_df) <- design_df[[sample.by]]
      design_df <- design_df[colnames(count_matrix), , drop = FALSE]
      design_df[[safe_split_by]] <- as.factor(design_df[[safe_split_by]])

      output <- list()
      graph_data <- list()
      for (i in seq_len(nrow(comparisons_condition))) {
        cluster_1 <- comparisons_condition[i, 1]
        cluster_2 <- comparisons_condition[i, 2]
        comparison_name <- paste0(cluster_1, "_vs_", cluster_2)

        keep <- design_df[[safe_split_by]] %in% unname(condition_map[c(cluster_1, cluster_2)])
        dsub <- droplevels(design_df[keep, , drop = FALSE])
        if (length(unique(dsub[[safe_split_by]])) < 2) {
          next
        }
        counts_sub <- count_matrix[, rownames(dsub), drop = FALSE]
        cell_sizes <- colSums(counts_sub)

        design_formula <- stats::as.formula(paste0("~ 0 + ", safe_split_by))
        x_model <- stats::model.matrix(design_formula, dsub)
        contrast_cols <- stats::setNames(colnames(x_model), levels(dsub[[safe_split_by]]))
        contrast <- paste0(
          contrast_cols[[condition_map[[cluster_2]]]],
          "-",
          contrast_cols[[condition_map[[cluster_1]]]]
        )

        dge <- edgeR::DGEList(counts = counts_sub, lib.size = cell_sizes)
        dge <- edgeR::calcNormFactors(dge, method = "TMM")
        dge <- edgeR::estimateDisp(dge, x_model)
        fit <- edgeR::glmQLFit(dge, x_model, robust = TRUE, legacy = TRUE)
        qlf <- edgeR::glmQLFTest(
          fit,
          contrast = limma::makeContrasts(contrasts = contrast, levels = x_model)
        )
        da <- as.data.frame(edgeR::topTags(qlf, sort.by = "none", n = Inf))
        da$Nhood <- seq_len(nrow(da))
        da$SpatialFDR <- milo_weighted_fdr_cpp(
          pvalues = da$PValue,
          weights = 1 / nhood[["k_distance"]]
        )
        da$neighborhood <- paste0("nhood_", da$Nhood)
        da$clusters <- da$neighborhood
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
          edge_df <- data.frame(
            from = character(0),
            to = character(0),
            stringsAsFactors = FALSE
          )
        }
        graph_data[[comparison_name]] <- list(nodes = node_df, edges = edge_df)
      }
      graph_data[[".metadata"]] <- nhood_metadata
      list(results = output, graph_data = graph_data)
    },
    error = function(e) {
      log_message(
        "{.pkg Milo} cpp backend execution failed: {.val {e$message}}",
        message_type = "error"
      )
    }
  )
}

extract_milo_r_neighborhood_metadata <- function(milo_obj) {
  cell_names <- colnames(milo_obj)
  nhood_counts <- tryCatch(as.matrix(miloR::nhoodCounts(milo_obj)), error = function(e) NULL)
  nhood_index <- tryCatch(miloR::nhoodIndex(milo_obj), error = function(e) NULL)
  nhood_matrix <- tryCatch(miloR::nhoods(milo_obj), error = function(e) NULL)

  members <- list()
  sampled_vertices <- integer(0)
  if (!is.null(nhood_matrix)) {
    if (nrow(nhood_matrix) == length(cell_names)) {
      members <- lapply(seq_len(ncol(nhood_matrix)), function(j) {
        cell_names[which(nhood_matrix[, j] != 0)]
      })
    } else if (ncol(nhood_matrix) == length(cell_names)) {
      members <- lapply(seq_len(nrow(nhood_matrix)), function(i) {
        cell_names[which(nhood_matrix[i, ] != 0)]
      })
    }
  }
  if (is.list(nhood_index) && length(nhood_index) > 0L) {
    sampled_vertices <- vapply(nhood_index, function(i) {
      idx <- as.integer(i)
      idx <- idx[!is.na(idx) & idx >= 1L]
      if (length(idx) == 0L) NA_integer_ else idx[1L]
    }, integer(1))
  }

  names(members) <- paste0("nhood_", seq_along(members))
  if (length(sampled_vertices) == 0L) {
    sampled_vertices <- rep(NA_integer_, length(members))
  }

  list(
    sampled_vertices = sampled_vertices,
    sampled_cells = ifelse(
      !is.na(sampled_vertices) & sampled_vertices >= 1L & sampled_vertices <= length(cell_names),
      cell_names[sampled_vertices],
      NA_character_
    ),
    members = members,
    counts = nhood_counts,
    sample_names = if (!is.null(nhood_counts)) colnames(nhood_counts) else character(0),
    cell_names = cell_names
  )
}

milo_refined_vertices <- function(random_vertices, coords, k) {
  vertex_knn <- BiocNeighbors::findKNN(
    coords,
    k = k,
    subset = as.integer(random_vertices),
    get.index = TRUE,
    get.distance = FALSE,
    BNPARAM = BiocNeighbors::KmknnParam()
  )
  nh_reduced_dims <- milo_neighborhood_medians_cpp(coords, vertex_knn[["index"]])
  colnames(nh_reduced_dims) <- colnames(coords)
  rownames(nh_reduced_dims) <- paste0("nh_", seq_len(nrow(nh_reduced_dims)))
  all_reduced_dims <- rbind(nh_reduced_dims, coords)
  nn_mat <- BiocNeighbors::findKNN(
    all_reduced_dims,
    k = nrow(nh_reduced_dims) + 1L,
    subset = seq_len(nrow(nh_reduced_dims)),
    get.index = TRUE,
    get.distance = FALSE,
    BNPARAM = BiocNeighbors::KmknnParam()
  )[["index"]]

  nh_ixs <- seq_len(nrow(nh_reduced_dims))
  i <- 1L
  sampled_vertices <- rep(0L, nrow(nn_mat))
  while (any(sampled_vertices <= max(nh_ixs))) {
    update_ix <- which(sampled_vertices <= max(nh_ixs))
    sampled_vertices[update_ix] <- nn_mat[update_ix, i]
    i <- i + 1L
    if (i > ncol(nn_mat)) {
      break
    }
  }
  sampled_vertices - max(nh_ixs)
}
