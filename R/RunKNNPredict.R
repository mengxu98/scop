#' @title Run KNN prediction
#'
#' @description
#' This function performs KNN prediction to annotate cell types based on reference scRNA-seq or bulk RNA-seq data.
#'
#' @param srt_query An object of class Seurat to be annotated with cell types.
#' @param srt_ref An object of class Seurat storing the reference cells.
#' @param bulk_ref A cell atlas matrix, where cell types are represented by columns and genes are represented by rows,
#' for example, scop::ref_scHCL. Either `srt_ref` or `bulk_ref` must be provided.
#' @param query_group A character vector specifying the column name in the `srt_query` metadata that represents the cell grouping.
#' @param ref_group A character vector specifying the column name in the `srt_ref` metadata that represents the cell grouping.
#' @param query_assay A character vector specifying the assay to be used for the query data.
#' Defaults to the default assay of the `srt_query` object.
#' @param ref_assay A character vector specifying the assay to be used for the reference data.
#' Defaults to the default assay of the `srt_ref` object.
#' @param query_reduction A character vector specifying the dimensionality reduction method used for the query data.
#' If NULL, the function will use the default reduction method specified in the `srt_query` object.
#' @param ref_reduction A character vector specifying the dimensionality reduction method used for the reference data.
#' If NULL, the function will use the default reduction method specified in the `srt_ref` object.
#' @param query_dims A numeric vector specifying the dimensions to be used for the query data.
#' Defaults to the first 30 dimensions.
#' @param ref_dims A numeric vector specifying the dimensions to be used for the reference data.
#' Defaults to the first 30 dimensions.
#' @param query_collapsing A boolean value indicating whether the query data should be collapsed to group-level average expression values. If TRUE, the function will calculate the average expression values for each group in the query data and the annotation will be performed separately for each group. Otherwise it will use the raw expression values for each cell.
#' @param ref_collapsing A boolean value indicating whether the reference data should be collapsed to group-level average expression values.
#' If TRUE, the function will calculate the average expression values for each group in the reference data and the annotation will be performed separately for each group.
#' Otherwise it will use the raw expression values for each cell.
#' @param return_full_distance_matrix A boolean value indicating whether the full distance matrix should be returned.
#' If TRUE, the function will return the distance matrix used for the KNN prediction, otherwise it will only return the annotated cell types.
#' @param features A character vector specifying the features (genes) to be used for the KNN prediction. If NULL, all the features in the query and reference data will be used.
#' @param features_type A character vector specifying the type of features to be used for the KNN prediction. Must be one of "HVF" (highly variable features) or "DE" (differentially expressed features). Defaults to "HVF".
#' @param feature_source A character vector specifying the source of the features to be used for the KNN prediction. Must be one of "both", "query", or "ref". Defaults to "both".
#' @param nfeatures An integer specifying the maximum number of features to be used for the KNN prediction. Defaults to 2000.
#' @param DEtest_param A list of parameters to be passed to the differential expression test function if `features_type` is set to "DE". Defaults to `list(max.cells.per.ident = 200, test.use = "wilcox")`.
#' @param DE_threshold Threshold used to filter the DE features.
#' Default is `"p_val < 0.05"`. If using "roc" test, \code{DE_threshold} should be needs to be reassigned. e.g. "power > 0.5".
#' @param nn_method A character vector specifying the method to be used for finding nearest neighbors. Must be one of "raw", "rann", or "annoy". Defaults to "raw".
#' @param distance_metric A character vector specifying the distance metric to be used for calculating similarity between cells. Must be one of "cosine", "euclidean", "manhattan", or "hamming". Defaults to "cosine".
#' @param k A number of nearest neighbors to be considered for the KNN prediction. Defaults to 30.
#' @param filter_lowfreq An integer specifying the threshold for filtering low-frequency cell types from the predicted results. Cell types with a frequency lower than `filter_lowfreq` will be labelled as "unreliable". Defaults to 0, which means no filtering will be performed.
#' @param prefix A character vector specifying the prefix to be added to the resulting annotations. Defaults to "KNNPredict".
#'
#' @export
#'
#' @examples
#' # Annotate cells using bulk RNA-seq data
#' data(pancreas_sub)
#' data(ref_scMCA)
#' pancreas_sub <- standard_scop(pancreas_sub)
#'
#' # Set the number of threads for RcppParallel
#' # details see: ?RcppParallel::setThreadOptions
#' # if (requireNamespace("RcppParallel", quietly = TRUE)) {
#' #   RcppParallel::setThreadOptions()
#' # }
#' pancreas_sub <- RunKNNPredict(
#'   srt_query = pancreas_sub,
#'   bulk_ref = ref_scMCA
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "KNNPredict_classification",
#'   label = TRUE
#' )
#'
#' # Removal of low credible cell types from the predicted results
#' pancreas_sub <- RunKNNPredict(
#'   srt_query = pancreas_sub,
#'   bulk_ref = ref_scMCA,
#'   filter_lowfreq = 30
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "KNNPredict_classification",
#'   label = TRUE
#' )
#'
#' # Annotate clusters using bulk RNA-seq data
#' pancreas_sub <- RunKNNPredict(
#'   srt_query = pancreas_sub,
#'   query_group = "SubCellType",
#'   bulk_ref = ref_scMCA
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "KNNPredict_classification",
#'   label = TRUE
#' )
#'
#' # Annotate using single cell RNA-seq data
#' data(panc8_sub)
#' # Simply convert genes from human to mouse and preprocess the data
#' genenames <- make.unique(
#'   thisutils::capitalize(
#'     rownames(panc8_sub),
#'     force_tolower = TRUE
#'   )
#' )
#' names(genenames) <- rownames(panc8_sub)
#' panc8_sub <- RenameFeatures(
#'   panc8_sub,
#'   newnames = genenames
#' )
#' panc8_sub <- CheckDataMerge(
#'   panc8_sub,
#'   batch = "tech"
#' )[["srt_merge"]]
#' panc8_sub <- SeuratObject::JoinLayers(panc8_sub)
#' pancreas_sub <- RunKNNPredict(
#'   srt_query = pancreas_sub,
#'   srt_ref = panc8_sub,
#'   ref_group = "celltype"
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "KNNPredict_classification",
#'   label = TRUE
#' )
#' FeatureDimPlot(
#'   pancreas_sub,
#'   features = "KNNPredict_simil"
#' )
#'
#' pancreas_sub <- RunKNNPredict(
#'   srt_query = pancreas_sub,
#'   srt_ref = panc8_sub,
#'   ref_group = "celltype",
#'   ref_collapsing = FALSE
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "KNNPredict_classification",
#'   label = TRUE
#' )
#' FeatureDimPlot(
#'   pancreas_sub,
#'   features = "KNNPredict_prob"
#' )
#'
#' pancreas_sub <- RunKNNPredict(
#'   srt_query = pancreas_sub,
#'   srt_ref = panc8_sub,
#'   query_group = "SubCellType",
#'   ref_group = "celltype"
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "KNNPredict_classification",
#'   label = TRUE
#' )
#' FeatureDimPlot(
#'   pancreas_sub,
#'   features = "KNNPredict_simil"
#' )
#'
#' # Annotate with DE gene instead of HVF
#' pancreas_sub <- RunKNNPredict(
#'   srt_query = pancreas_sub,
#'   srt_ref = panc8_sub,
#'   ref_group = "celltype",
#'   features_type = "DE",
#'   feature_source = "ref",
#'   DEtest_param = list(cores = 2)
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "KNNPredict_classification",
#'   label = TRUE
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   features = "KNNPredict_simil"
#' )
#'
#' pancreas_sub <- RunKNNPredict(
#'   srt_query = pancreas_sub,
#'   srt_ref = panc8_sub,
#'   query_group = "SubCellType",
#'   ref_group = "celltype",
#'   features_type = "DE",
#'   feature_source = "both",
#'   DEtest_param = list(cores = 2)
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "KNNPredict_classification",
#'   label = TRUE
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   features = "KNNPredict_simil"
#' )
RunKNNPredict <- function(
    srt_query,
    srt_ref = NULL,
    bulk_ref = NULL,
    query_group = NULL,
    ref_group = NULL,
    query_assay = NULL,
    ref_assay = NULL,
    query_reduction = NULL,
    ref_reduction = NULL,
    query_dims = 1:30,
    ref_dims = 1:30,
    query_collapsing = !is.null(query_group),
    ref_collapsing = TRUE,
    return_full_distance_matrix = FALSE,
    features = NULL,
    features_type = c("HVF", "DE"),
    feature_source = "both",
    nfeatures = 2000,
    DEtest_param = list(
      max.cells.per.ident = 200,
      test.use = "wilcox"
    ),
    DE_threshold = "p_val_adj < 0.05",
    nn_method = NULL,
    distance_metric = "cosine",
    k = 30,
    filter_lowfreq = 0,
    prefix = "KNNPredict") {
  query_assay <- query_assay %||% SeuratObject::DefaultAssay(srt_query)
  features_type <- match.arg(features_type)
  if (is.null(query_reduction) + is.null(ref_reduction) == 1) {
    log_message(
      "query_reduction and ref_reduction must be both provided",
      message_type = "error"
    )
  }
  if (is.null(query_reduction) + is.null(ref_reduction) == 2) {
    use_reduction <- FALSE
  } else {
    use_reduction <- TRUE
  }
  if (isFALSE(use_reduction)) {
    if (length(features) == 0) {
      if (features_type == "HVF" && feature_source %in% c("both", "query")) {
        if (length(SeuratObject::VariableFeatures(srt_query, assay = query_assay)) == 0) {
          log_message("Perform {.fn Seurat::FindVariableFeatures} on the query data...")
          srt_query <- Seurat::FindVariableFeatures(
            srt_query,
            nfeatures = nfeatures,
            assay = query_assay,
            verbose = FALSE
          )
        }
        features_query <- SeuratObject::VariableFeatures(srt_query, assay = query_assay)
      } else if (features_type == "DE" && feature_source %in% c("both", "query")) {
        if (is.null(query_group)) {
          log_message(
            "'query_group' must be provided when 'features_type' is 'DE' and 'feature_source' is 'both' or 'query'",
            message_type = "error"
          )
        } else {
          layer <- paste0("DEtest_", query_group)
          DEtest_param[["force"]] <- TRUE
          if (
            !layer %in% names(srt_query@tools) ||
              length(grep(
                pattern = "AllMarkers",
                names(srt_query@tools[[layer]])
              )) ==
                0
          ) {
            srt_query <- do.call(
              RunDEtest,
              c(
                list(srt = srt_query, group_by = query_group),
                DEtest_param
              )
            )
          }
          if ("test.use" %in% names(DEtest_param)) {
            test.use <- DEtest_param[["test.use"]]
          } else {
            test.use <- "wilcox"
          }
          index <- grep(
            pattern = paste0("AllMarkers_", test.use),
            names(srt_query@tools[[layer]])
          )[1]
          de <- names(srt_query@tools[[layer]])[index]
          log_message(
            "Use the DE features from ",
            de,
            " to calculate distance metric."
          )
          de_df <- srt_query@tools[[layer]][[de]]
          de_df <- de_df[
            with(de_df, eval(rlang::parse_expr(DE_threshold))), ,
            drop = FALSE
          ]
          rownames(de_df) <- seq_len(nrow(de_df))
          de_df <- de_df[
            order(de_df[["avg_log2FC"]], decreasing = TRUE), ,
            drop = FALSE
          ]
          de_top <- de_df[!duplicated(de_df[["gene"]]), , drop = FALSE]
          stat <- sort(table(de_top$group1))
          stat <- stat[stat > 0]
          mat <- matrix(FALSE, nrow = max(stat), ncol = length(stat))
          colnames(mat) <- names(stat)
          for (g in names(stat)) {
            mat[1:stat[g], g] <- TRUE
          }
          nfeatures <- sum(cumsum(Matrix::rowSums(mat)) <= nfeatures)
          if (test.use == "roc") {
            features_query <- unlist(by(
              de_top,
              list(de_top[["group1"]]),
              function(x) {
                x <- x[order(x[["power"]], decreasing = TRUE), , drop = FALSE]
                utils::head(x[["gene"]], nfeatures)
              }
            ))
          } else {
            features_query <- unlist(by(
              de_top,
              list(de_top[["group1"]]),
              function(x) {
                x <- x[order(x[["p_val"]], decreasing = FALSE), , drop = FALSE]
                utils::head(x[["gene"]], nfeatures)
              }
            ))
          }
          log_message(
            "DE features number of the query data: ",
            length(features_query)
          )
        }
      } else {
        features_query <- rownames(srt_query[[query_assay]])
      }
    }
  }

  if (!is.null(bulk_ref)) {
    if (length(features) == 0) {
      features <- features_query
    }
    features_common <- intersect(features, rownames(bulk_ref))
    log_message("Use ", length(features_common), " features to calculate distance.")
    ref <- Matrix::t(bulk_ref[features_common, , drop = FALSE])
  } else if (!is.null(srt_ref)) {
    ref_assay <- ref_assay %||% SeuratObject::DefaultAssay(srt_ref)
    if (!is.null(ref_group)) {
      if (length(ref_group) == ncol(srt_ref)) {
        srt_ref[["ref_group"]] <- ref_group
      } else if (length(ref_group) == 1) {
        if (!ref_group %in% colnames(srt_ref@meta.data)) {
          log_message(
            "ref_group must be one of the column names in the meta.data",
            message_type = "error"
          )
        } else {
          srt_ref[["ref_group"]] <- srt_ref[[ref_group]]
        }
      } else {
        log_message(
          "Length of ref_group must be one or length of srt_ref.",
          message_type = "error"
        )
      }
    } else {
      log_message(
        "ref_group must be provided.",
        message_type = "error"
      )
    }

    drop_cell <- colnames(srt_ref)[is.na(srt_ref[["ref_group", drop = TRUE]])]
    if (length(drop_cell) > 0) {
      log_message("Drop ", length(drop_cell), " cells with NA in the ref_group")
      srt_ref <- srt_ref[, setdiff(colnames(srt_ref), drop_cell)]
    }
    if (isTRUE(use_reduction)) {
      log_message("Use the reduction to calculate distance metric.")
      if (
        !is.null(query_dims) &&
          !is.null(ref_dims) &&
          length(query_dims) == length(ref_dims)
      ) {
        query <- Embeddings(srt_query, reduction = query_reduction)[
          ,
          query_dims
        ]
        ref <- Embeddings(srt_ref, reduction = ref_reduction)[, ref_dims]
      } else {
        log_message(
          "query_dims and ref_dims must be provided with the same length.",
          message_type = "error"
        )
      }
    } else {
      if (length(features) == 0) {
        if (features_type == "HVF" && feature_source %in% c("both", "ref")) {
          log_message("Use the HVF to calculate distance metric")
          if (length(SeuratObject::VariableFeatures(srt_ref, assay = ref_assay)) == 0) {
            srt_ref <- Seurat::FindVariableFeatures(
              srt_ref,
              nfeatures = nfeatures,
              assay = ref_assay
            )
          }
          features_ref <- SeuratObject::VariableFeatures(srt_ref, assay = ref_assay)
        } else if (
          features_type == "DE" && feature_source %in% c("both", "ref")
        ) {
          layer <- paste0("DEtest_", ref_group)
          DEtest_param[["force"]] <- TRUE
          if (
            !layer %in% names(srt_ref@tools) ||
              length(grep(
                pattern = "AllMarkers",
                names(srt_ref@tools[[layer]])
              )) ==
                0
          ) {
            srt_ref <- do.call(
              RunDEtest,
              c(list(srt = srt_ref, group_by = ref_group), DEtest_param)
            )
          }
          if ("test.use" %in% names(DEtest_param)) {
            test.use <- DEtest_param[["test.use"]]
          } else {
            test.use <- "wilcox"
          }
          index <- grep(
            pattern = paste0("AllMarkers_", test.use),
            names(srt_ref@tools[[layer]])
          )[1]
          de <- names(srt_ref@tools[[layer]])[index]
          log_message(
            "Use the DE features from ",
            de,
            " to calculate distance metric."
          )
          de_df <- srt_ref@tools[[layer]][[de]]
          de_df <- de_df[
            with(de_df, eval(rlang::parse_expr(DE_threshold))), ,
            drop = FALSE
          ]
          rownames(de_df) <- seq_len(nrow(de_df))
          de_df <- de_df[
            order(de_df[["avg_log2FC"]], decreasing = TRUE), ,
            drop = FALSE
          ]
          de_top <- de_df[!duplicated(de_df[["gene"]]), , drop = FALSE]
          stat <- sort(table(de_top$group1))
          stat <- stat[stat > 0]
          mat <- matrix(FALSE, nrow = max(stat), ncol = length(stat))
          colnames(mat) <- names(stat)
          for (g in names(stat)) {
            mat[1:stat[g], g] <- TRUE
          }
          nfeatures <- sum(cumsum(Matrix::rowSums(mat)) <= nfeatures)
          if (test.use == "roc") {
            features_ref <- unlist(by(
              de_top,
              list(de_top[["group1"]]),
              function(x) {
                x <- x[order(x[["power"]], decreasing = TRUE), , drop = FALSE]
                utils::head(x[["gene"]], nfeatures)
              }
            ))
          } else {
            features_ref <- unlist(by(
              de_top,
              list(de_top[["group1"]]),
              function(x) {
                x <- x[order(x[["p_val"]], decreasing = FALSE), , drop = FALSE]
                utils::head(x[["gene"]], nfeatures)
              }
            ))
          }
          log_message("DE features number of the ref data: ", length(features_ref))
        } else {
          features_ref <- rownames(srt_ref[[ref_assay]])
        }
        features <- intersect(features_ref, features_query)
      }
      features_common <- Reduce(
        intersect,
        list(
          features,
          rownames(srt_query[[query_assay]]),
          rownames(srt_ref[[ref_assay]])
        )
      )
      log_message(
        "Use ",
        length(features_common),
        " features to calculate distance."
      )
      if (isTRUE(ref_collapsing)) {
        ref <- Seurat::AverageExpression(
          # ref <- Seurat::PseudobulkExpression( # require run JoinLayers
          object = srt_ref,
          features = features_common,
          layer = "data",
          assays = ref_assay,
          group.by = "ref_group",
          verbose = FALSE
        )[[1]]
        ref <- Matrix::t(log1p(ref))
      } else {
        ref <- Matrix::t(
          GetAssayData5(
            srt_ref,
            layer = "data",
            assay = ref_assay
          )[features_common, ]
        )
      }
    }
  } else {
    log_message(
      "srt_ref or bulk_ref must be provided at least one",
      message_type = "error"
    )
  }

  if (!inherits(ref, "matrix")) {
    ref <- as_matrix(ref)
  }
  k <- min(c(k, nrow(ref)))

  if (!is.null(query_group)) {
    if (length(query_group) == ncol(srt_query)) {
      srt_query[["query_group"]] <- query_group
    } else if (length(query_group) == 1) {
      if (!query_group %in% colnames(srt_query@meta.data)) {
        log_message(
          "query_group must be one of the column names in the meta.data",
          message_type = "error"
        )
      } else {
        srt_query[["query_group"]] <- srt_query[[query_group]]
      }
    } else {
      log_message(
        "Length of query_group must be one or length of srt_query.",
        message_type = "error"
      )
    }
  }
  if (isFALSE(use_reduction)) {
    query_assay <- query_assay %||% SeuratObject::DefaultAssay(srt_query)
    if (isTRUE(query_collapsing)) {
      if (is.null(query_group)) {
        log_message(
          "{.arg query_group} must be provided when query_collapsing is TRUE",
          message_type = "error"
        )
      }
      query <- Seurat::AverageExpression(
        # query <- Seurat::PseudobulkExpression( # require run JoinLayers
        object = srt_query,
        features = colnames(ref),
        layer = "data",
        assays = query_assay,
        group.by = "query_group",
        verbose = FALSE
      )[[1]]
      query <- Matrix::t(log1p(query))
    } else {
      query <- Matrix::t(
        GetAssayData5(
          srt_query,
          layer = "data",
          assay = query_assay
        )[colnames(ref), , drop = FALSE]
      )
    }
  }

  if (isFALSE(use_reduction)) {
    status_dat <- CheckDataType(query, verbose = FALSE)
    log_message("Detected query data type: {.val {status_dat}}")
    status_ref <- CheckDataType(ref, verbose = FALSE)
    log_message("Detected reference data type: {.val {status_ref}}")
    if (status_ref != status_dat || any(status_dat == "unknown", status_ref == "unknown")) {
      log_message(
        "Data type is unknown or different between query and reference",
        message_type = "warning"
      )
    }
  }

  log_message("Calculate similarity...")

  if (is.null(nn_method)) {
    if (as.numeric(nrow(query)) * as.numeric(nrow(ref)) >= 1e8) {
      nn_method <- "annoy"
    } else {
      nn_method <- "raw"
    }
  }
  log_message("Use {.pkg {nn_method}} method to find neighbors")
  if (!nn_method %in% c("raw", "annoy", "rann")) {
    log_message("nn_method must be one of raw, rann and annoy",
      message_type = "error"
    )
  }
  if (
    nn_method == "annoy" &&
      !distance_metric %in% c("euclidean", "cosine", "manhattan", "hamming")
  ) {
    log_message(
      "distance_metric must be one of euclidean, cosine, manhattan, and hamming when nn_method='annoy'",
      message_type = "error"
    )
  }
  if (isTRUE(return_full_distance_matrix) && nn_method != "raw") {
    log_message(
      "Distance matrix will not be returned besause nn_method is not 'raw'",
      message_type = "warning"
    )
    return_full_distance_matrix <- FALSE
  }
  simil_methods <- c(
    "cosine",
    "pearson",
    "spearman",
    "correlation",
    "jaccard",
    "ejaccard",
    "dice",
    "edice",
    "hamman",
    "simple matching",
    "faith"
  )
  dist_methods <- c(
    "euclidean",
    "chisquared",
    "kullback",
    "manhattan",
    "maximum",
    "canberra",
    "minkowski",
    "hamming"
  )
  if (!(distance_metric %in% c(simil_methods, dist_methods))) {
    log_message(
      "{.val {distance_metric}} method is invalid",
      message_type = "error"
    )
  }

  if (nn_method %in% c("annoy", "rann")) {
    query.neighbor <- Seurat::FindNeighbors(
      query = query,
      object = ref,
      k.param = k,
      nn.method = nn_method,
      annoy.metric = distance_metric,
      return.neighbor = TRUE
    )
    match_k <- query.neighbor@nn.idx
    rownames(match_k) <- rownames(query)
    match_k_cell <- apply(match_k, c(1, 2), function(x) rownames(ref)[x])
    match_k_distance <- query.neighbor@nn.dist
    rownames(match_k_distance) <- rownames(query)
  } else {
    if (distance_metric %in% c(simil_methods, "pearson", "spearman")) {
      if (distance_metric %in% c("pearson", "spearman")) {
        if (distance_metric == "spearman") {
          ref <- Matrix::t(apply(ref, 1, rank))
          query <- Matrix::t(apply(query, 1, rank))
        }
        distance_metric <- "correlation"
      }
      d <- 1 -
        proxyC::simil(
          x = SeuratObject::as.sparse(ref),
          y = SeuratObject::as.sparse(query),
          method = distance_metric,
          use_nan = TRUE
        )
    } else if (distance_metric %in% dist_methods) {
      d <- proxyC::dist(
        x = SeuratObject::as.sparse(ref),
        y = SeuratObject::as.sparse(query),
        method = distance_metric,
        use_nan = TRUE
      )
    }
    if (k == 1) {
      match_k_cell <- as_matrix(
        apply(d, 2, function(x) {
          names(x)[order(x, decreasing = FALSE)[1]]
        })
      )
      match_k_distance <- as_matrix(
        apply(d, 2, function(x) x[order(x, decreasing = FALSE)[1]])
      )
    } else {
      match_k_cell <- Matrix::t(
        as_matrix(
          apply(d, 2, function(x) names(x)[order(x, decreasing = FALSE)[1:k]])
        )
      )
      match_k_distance <- Matrix::t(
        as_matrix(
          apply(d, 2, function(x) x[order(x, decreasing = FALSE)[1:k]])
        )
      )
    }
  }

  log_message("Predict cell type...")
  match_prob <- NULL
  if (!is.null(srt_ref) && (isFALSE(ref_collapsing) || isTRUE(use_reduction))) {
    level <- as.character(unique(srt_ref[["ref_group", drop = TRUE]]))
    if (k == 1) {
      match_best <- srt_ref[["ref_group", drop = TRUE]][match_k_cell[, 1]]
      names(match_best) <- names(match_k_cell[, 1])
    } else {
      rn <- rownames(match_k_cell)
      match_k_cell <- matrix(
        srt_ref[["ref_group", drop = TRUE]][match_k_cell],
        nrow = nrow(match_k_cell),
        ncol = ncol(match_k_cell)
      )
      rownames(match_k_cell) <- rn
      match_freq <- apply(match_k_cell, 1, table)
      if (!inherits(match_freq, "list")) {
        match_freq <- as.list(
          stats::setNames(
            object = rep(k, nrow(match_k_cell)), rn
          )
        )
        match_freq <- lapply(
          stats::setNames(
            names(match_freq),
            names(match_freq)
          ),
          function(x) stats::setNames(k, match_k_cell[x, 1])
        )
      }
      match_prob <- do.call(
        rbind,
        lapply(match_freq, function(x) {
          x[level[!level %in% names(x)]] <- 0
          x <- x / sum(x)
          return(x)
        })
      )
      match_prob <- as_matrix(match_prob)
      rownames(match_prob) <- names(match_freq)
      match_best <- apply(
        match_prob,
        1,
        function(x) names(x)[order(x, decreasing = TRUE)][1]
      )
    }
  } else {
    match_best <- match_k_cell[, 1]
  }

  result <- list(
    features = features,
    nn_method = nn_method,
    distance_metric = distance_metric,
    k = k,
    other_params = list(
      query_group = query_group,
      query_reduction = query_reduction,
      query_assay = query_assay,
      query_dims = query_dims,
      query_collapsing = query_collapsing,
      ref_group = ref_group,
      ref_reduction = ref_reduction,
      ref_assay = ref_assay,
      ref_dims = ref_dims,
      ref_collapsing = ref_collapsing
    )
  )
  result[["match_best"]] <- match_best
  if (!is.null(match_prob)) {
    result[["match_prob"]] <- match_prob
  }
  result[["match_k_cell"]] <- match_k_cell
  result[["match_k_distance"]] <- match_k_distance
  if (isTRUE(return_full_distance_matrix)) {
    result[["distance_matrix"]] <- d[seq_len(nrow(d)), , drop = FALSE]
  }

  srt_query@tools[[paste0(prefix, "_classification")]] <- result

  if (isTRUE(query_collapsing)) {
    query_index <- as.character(
      srt_query[["query_group", drop = TRUE]]
    )
    cell_type_to_id <- split(colnames(srt_query), query_index)
    classification <- unlist(
      lapply(
        names(match_best), function(ct) {
          rep(match_best[ct], length(cell_type_to_id[[ct]]))
        }
      )
    )
    names(classification) <- unlist(cell_type_to_id)
  } else {
    query_index <- colnames(srt_query)
    classification <- match_best[query_index]
    names(classification) <- query_index
  }

  srt_query[[paste0(prefix, "_classification")]] <- classification

  if (!is.null(match_prob)) {
    if (isTRUE(query_collapsing)) {
      prob <- unlist(lapply(names(match_best), function(ct) {
        rep(apply(match_prob, 1, max)[ct], length(cell_type_to_id[[ct]]))
      }))
      names(prob) <- unlist(cell_type_to_id)
    } else {
      prob <- apply(match_prob, 1, max)[query_index]
      names(prob) <- query_index
    }
    srt_query[[paste0(prefix, "_prob")]] <- prob
  } else {
    distance <- match_k_distance[, 1]
    if (distance_metric %in% c(simil_methods, "pearson", "spearman")) {
      if (isTRUE(query_collapsing)) {
        simil <- unlist(lapply(names(match_best), function(ct) {
          rep((1 - distance)[ct], length(cell_type_to_id[[ct]]))
        }))
        names(simil) <- unlist(cell_type_to_id)
      } else {
        simil <- (1 - distance)[query_index]
        names(simil) <- query_index
      }
      srt_query[[paste0(prefix, "_simil")]] <- simil
    } else {
      if (isTRUE(query_collapsing)) {
        dist <- unlist(lapply(names(match_best), function(ct) {
          rep(distance[ct], length(cell_type_to_id[[ct]]))
        }))
        names(dist) <- unlist(cell_type_to_id)
      } else {
        dist <- distance[query_index]
        names(dist) <- query_index
      }
      srt_query[[paste0(prefix, "_dist")]] <- dist
    }
  }

  if (is.numeric(filter_lowfreq) && filter_lowfreq > 0) {
    drop <- table(srt_query[[paste0(prefix, "_classification"), drop = TRUE]])
    drop <- names(drop)[drop <= filter_lowfreq]
    srt_query[[paste0(prefix, "_classification"), drop = TRUE]][
      srt_query[[paste0(prefix, "_classification"), drop = TRUE]] %in% drop
    ] <- "unreliable"
  }

  return(srt_query)
}
