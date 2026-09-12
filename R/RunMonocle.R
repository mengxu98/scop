#' @title Run Monocle2 analysis
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams RunStandardWorkflow
#' @inheritParams GroupHeatmap
#' @inheritParams CellDimPlot
#' @inheritParams scop-params
#' @param expressionFamily The distribution family to use for modeling gene expression.
#' @param features Features to use.
#' Defaults to NULL, in which case features were determined by `feature_type`.
#' @param feature_type The type of features to use in the analysis.
#' Possible values are "HVF" for highly variable features or "Disp" for features selected based on dispersion.
#' @param disp_filter Filter to use when `feature_type` is "Disp".
#' @param max_components The maximum number of dimensions to use for dimensionality reduction.
#' @param reduction_method The dimensionality reduction method to use.
#' Possible values are `"DDRTree"`, `"ICA"`, `"tSNE"`,
#' `"SimplePPT"`, `"L1-graph"`, `"SGL-tree"`.
#' @param norm_method The normalization method to use.
#' Possible values are `"log"` and `"none"`.
#' @param residualModelFormulaStr A model formula specifying the effects to subtract.
#' @param pseudo_expr Amount to increase expression values before dimensionality reduction.
#' @param root_state The state to use as the root of the trajectory.
#' If NULL, the trajectory is rooted at the first state after an initial
#' ordering.
#' @param backend Deprecated and ignored. Monocle now accelerates cell ordering
#' natively with C++, so the separate `"cpp"` backend has been removed.
#' @param n_neighbors Deprecated and ignored.
#' @param ddrtree_maxIter Optional maximum iteration count passed to
#' `DDRTree::DDRTree()` when `reduction_method = "DDRTree"`. Lower values can
#' speed up exploratory runs but may reduce agreement with the default Monocle2
#' trajectory.
#' @param ddrtree_ncenter Optional number of DDRTree centers. This can change
#' trajectory topology and is intended for advanced exploratory use.
#' @param ddrtree_tol Optional convergence tolerance passed to `DDRTree::DDRTree()`.
#' @param show_plot Whether to print diagnostic plots during the run.
#'
#' @export
#' @seealso
#' [RunSlingshot],
#' [CellDimPlot],
#' [FeatureDimPlot]
#'
#' @return
#' A Seurat object with the Monocle2 analysis results added to the @tools slot.
#'
#' @examples
#' if (interactive()) {
#'   data(pancreas_sub)
#'   pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#'   pancreas_sub <- RunMonocle2(
#'     pancreas_sub,
#'     group.by = "SubCellType"
#'   )
#'   names(pancreas_sub@tools$Monocle2)
#'   trajectory <- pancreas_sub@tools$Monocle2$trajectory
#'
#'   p1 <- CellDimPlot(
#'     pancreas_sub,
#'     group.by = "Monocle2_State",
#'     reduction = "DDRTree",
#'     label = TRUE,
#'     theme_use = "theme_blank"
#'   )
#'   p1
#'
#'   p1 + trajectory
#'
#'   FeatureDimPlot(
#'     pancreas_sub,
#'     features = "Monocle2_Pseudotime",
#'     reduction = "UMAP",
#'     theme_use = "theme_blank"
#'   )
#'
#'   pancreas_sub <- RunMonocle2(
#'     pancreas_sub,
#'     feature_type = "Disp",
#'     disp_filter = "mean_expression >= 0.01 & dispersion_empirical >= 1 * dispersion_fit"
#'   )
#'   trajectory <- pancreas_sub@tools$Monocle2$trajectory
#'   p2 <- CellDimPlot(
#'     pancreas_sub,
#'     group.by = "Monocle2_State",
#'     reduction = "DDRTree",
#'     label = TRUE,
#'     theme_use = "theme_blank"
#'   )
#'   p2
#'
#'   p2 + trajectory
#'
#'   FeatureDimPlot(
#'     pancreas_sub,
#'     features = "Monocle2_Pseudotime",
#'     reduction = "UMAP",
#'     theme_use = "theme_blank"
#'   )
#' }
RunMonocle2 <- function(
  object,
  assay = NULL,
  layer = "counts",
  group.by = NULL,
  expressionFamily = "negbinomial.size",
  features = NULL,
  feature_type = "HVF",
  disp_filter = "mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit",
  max_components = 2,
  reduction_method = "DDRTree",
  norm_method = "log",
  residualModelFormulaStr = NULL,
  pseudo_expr = 1,
  root_state = NULL,
  backend = NULL,
  n_neighbors = NULL,
  ddrtree_maxIter = NULL,
  ddrtree_ncenter = NULL,
  ddrtree_tol = NULL,
  show_plot = FALSE,
  xlab = NULL,
  ylab = NULL,
  seed = 11,
  verbose = TRUE,
  srt = NULL
) {
  srt <- resolve_deprecated_srt(object, srt, missing(object))
  if (!is.null(backend)) {
    log_message(
      "{.arg backend} is deprecated and ignored: {.pkg monocle} now
      accelerates cell ordering natively.",
      message_type = "warning"
    )
  }
  if (!is.null(n_neighbors)) {
    log_message(
      "{.arg n_neighbors} is deprecated and ignored.",
      message_type = "warning"
    )
  }
  log_message("Run {.pkg monocle2}...", verbose = verbose)
  set.seed(seed)
  check_r(
    c("mengxu98/monocle", "DDRTree", "BiocGenerics", "Biobase", "VGAM"),
    verbose = FALSE
  )

  if (!"package:DDRTree" %in% search()) {
    suppressMessages(
      suppressWarnings(
        attachNamespace("DDRTree")
      )
    )
  }
  if (!"package:Biobase" %in% search()) {
    suppressMessages(
      suppressWarnings(
        attachNamespace("Biobase")
      )
    )
  }

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  expr_matrix <- SeuratObject::as.sparse(
    GetAssayData5(
      srt,
      assay = assay,
      layer = layer
    )
  )
  f_data <- data.frame(
    gene_short_name = row.names(expr_matrix),
    row.names = row.names(expr_matrix)
  )
  pd <- methods::new("AnnotatedDataFrame", data = srt@meta.data)
  fd <- methods::new("AnnotatedDataFrame", data = f_data)
  cds <- get_namespace_fun("monocle", "newCellDataSet")(
    expr_matrix,
    phenoData = pd,
    featureData = fd,
    expressionFamily = do.call(
      get(
        expressionFamily,
        envir = getNamespace("VGAM")
      ),
      args = list()
    )
  )
  need_dispersion <- monocle2_needs_dispersion(
    features = features,
    feature_type = feature_type,
    show_plot = show_plot
  )
  uses_negbinomial <- any(c("negbinomial", "negbinomial.size") %in% expressionFamily)
  if (isTRUE(uses_negbinomial)) {
    cds <- get_namespace_fun("BiocGenerics", "estimateSizeFactors")(cds)
    if (isTRUE(need_dispersion)) {
      cds <- get_namespace_fun("BiocGenerics", "estimateDispersions")(cds)
    }
  }
  dispersion_table <- NULL
  if (is.null(features)) {
    if (feature_type == "HVF") {
      features <- SeuratObject::VariableFeatures(
        srt,
        assay = assay
      )
      if (length(features) == 0) {
        features <- SeuratObject::VariableFeatures(
          FindVariableFeatures(srt, assay = assay, verbose = FALSE),
          assay = assay
        )
      }
    }
    if (feature_type == "Disp") {
      dispersion_table <- get_namespace_fun("monocle", "dispersionTable")(cds)
      features <- subset(
        dispersion_table,
        eval(rlang::parse_expr(disp_filter))
      )$gene_id
    }
  }
  log_message(
    "{.val {length(features)}} features selected",
    verbose = verbose
  )
  cds <- get_namespace_fun("monocle", "setOrderingFilter")(cds, features)
  if (isTRUE(show_plot) && isTRUE(uses_negbinomial)) {
    p <- get_namespace_fun("monocle", "plot_ordering_genes")(cds)
    print(p)
  }

  reduce_args <- monocle2_reduce_dimension_args(
    cds = cds,
    max_components = max_components,
    reduction_method = reduction_method,
    norm_method = norm_method,
    residualModelFormulaStr = residualModelFormulaStr,
    pseudo_expr = pseudo_expr,
    ddrtree_maxIter = ddrtree_maxIter,
    ddrtree_ncenter = ddrtree_ncenter,
    ddrtree_tol = ddrtree_tol
  )
  cds <- do.call(get_namespace_fun("monocle", "reduceDimension"), reduce_args)
  cds <- get_namespace_fun("monocle", "orderCells")(cds)

  embeddings <- Matrix::t(cds@reducedDimS)
  colnames(embeddings) <- paste0(
    cds@dim_reduce_type, "_", seq_len(ncol(embeddings))
  )
  srt[[cds@dim_reduce_type]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings,
    key = paste0(cds@dim_reduce_type, "_"),
    assay = assay
  )
  srt[["Monocle2_State"]] <- cds[["State"]]
  if (cds@dim_reduce_type == "ICA") {
    reduced_dim_coords <- as.data.frame(Matrix::t(cds@reducedDimS))
  } else if (cds@dim_reduce_type %in% c("simplePPT", "DDRTree")) {
    reduced_dim_coords <- as.data.frame(Matrix::t(cds@reducedDimK))
  }
  edge_df <- igraph::as_data_frame(cds@minSpanningTree)
  edge_df[, c("x", "y")] <- reduced_dim_coords[edge_df[["from"]], 1:2]
  edge_df[, c("xend", "yend")] <- reduced_dim_coords[edge_df[["to"]], 1:2]
  trajectory <- ggplot2::geom_segment(
    data = edge_df,
    aes(x = x, y = y, xend = xend, yend = yend)
  )

  if (isTRUE(show_plot)) {
    p1 <- CellDimPlot(
      srt,
      group.by = "Monocle2_State",
      reduction = cds@dim_reduce_type,
      label = TRUE,
      force = TRUE,
      xlab = xlab,
      ylab = ylab
    ) +
      trajectory
    p3 <- FeatureDimPlot(
      srt,
      features = "Monocle2_Pseudotime",
      reduction = cds@dim_reduce_type,
      xlab = xlab,
      ylab = ylab
    ) +
      trajectory
    if (!is.null(group.by)) {
      p2 <- CellDimPlot(
        srt,
        group.by = group.by,
        reduction = cds@dim_reduce_type,
        label = TRUE,
        force = TRUE,
        xlab = xlab,
        ylab = ylab
      ) +
        trajectory
      print(p1 + p2 + p3)
    } else {
      print(p1 + p3)
    }
  }

  if (is.null(root_state)) {
    root_state <- sort(unique(cds[["State"]]))[1]
  }
  cds <- get_namespace_fun("monocle", "orderCells")(cds, root_state = root_state)
  srt[["Monocle2_Pseudotime"]] <- cds[["Pseudotime"]]
  srt@tools$Monocle2 <- list(
    cds = cds,
    features = features,
    trajectory = trajectory
  )

  log_message(
    "{.pkg monocle2} completed",
    message_type = "success",
    verbose = verbose
  )

  return(srt)
}

monocle2_needs_dispersion <- function(features, feature_type, show_plot) {
  isTRUE(show_plot) || (is.null(features) && identical(feature_type, "Disp"))
}

monocle2_reduce_dimension_args <- function(
  cds,
  max_components = 2,
  reduction_method = "DDRTree",
  norm_method = "log",
  residualModelFormulaStr = NULL,
  pseudo_expr = 1,
  ddrtree_maxIter = NULL,
  ddrtree_ncenter = NULL,
  ddrtree_tol = NULL
) {
  args <- list(
    cds = cds,
    max_components = max_components,
    reduction_method = reduction_method,
    norm_method = norm_method,
    residualModelFormulaStr = residualModelFormulaStr,
    pseudo_expr = pseudo_expr
  )
  if (identical(reduction_method, "DDRTree")) {
    if (!is.null(ddrtree_maxIter)) {
      args[["maxIter"]] <- as.integer(ddrtree_maxIter)
    }
    if (!is.null(ddrtree_ncenter)) {
      args[["ncenter"]] <- as.integer(ddrtree_ncenter)
    }
    if (!is.null(ddrtree_tol)) {
      args[["tol"]] <- as.numeric(ddrtree_tol)
    }
  }
  args
}

#' @title Run Monocle3 analysis
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams RunStandardWorkflow
#' @inheritParams FeatureDimPlot
#' @inheritParams GroupHeatmap
#' @inheritParams CellDimPlot
#' @param clusters The cluster variable in the Seurat object to use for analysis.
#' Defaults to NULL, in which case use Monocle clusters is used.
#' @param graph Graph slot in the Seurat object to use for analysis.
#' Defaults to NULL, in which case Monocle graph is used.
#' @param partition_qval The q-value threshold for partitioning cells.
#' Defaults to 0.05.
#' @param k The number of nearest neighbors to consider for clustering.
#' Defaults to 50.
#' @param cluster_method The clustering method to use.
#' Defaults to "louvain".
#' @param num_iter The number of iterations for clustering.
#' Defaults to 2.
#' @param resolution The resolution parameter for clustering.
#' Defaults to NULL.
#' @param use_partition Whether to use partitions to learn disjoint graph in each partition.
#' Defaults to NULL, in which case TRUE is used.
#' @param close_loop Whether to close loops in the graph.
#' Defaults to TRUE.
#' @param root_pr_nodes The root nodes to order cells.
#' If not specified, user will be prompted for input.
#' Defaults to NULL.
#' @param root_cells The root cells to order cells.
#' If not specified, user will be prompted for input.
#' Defaults to NULL.
#' @param show_plot Whether to print diagnostic plots during the run.
#'
#' @export
#'
#' @examples
#' if (interactive()) {
#'   data(pancreas_sub)
#'   pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#'   pancreas_sub <- RunMonocle3(
#'     pancreas_sub,
#'     reduction = "UMAP"
#'   )
#'
#'   pancreas_sub <- RunMonocle3(
#'     pancreas_sub,
#'     reduction = "UMAP",
#'     group.by = "CellType"
#'   )
#'   names(pancreas_sub@tools$Monocle3)
#'   trajectory <- pancreas_sub@tools$Monocle3$trajectory
#'   milestones <- pancreas_sub@tools$Monocle3$milestones
#'
#'   CellDimPlot(
#'     pancreas_sub,
#'     group.by = "Monocle3_partitions",
#'     reduction = "UMAP",
#'     label = TRUE,
#'     theme_use = "theme_blank"
#'   ) +
#'     trajectory +
#'     milestones
#'   CellDimPlot(
#'     pancreas_sub,
#'     group.by = "Monocle3_clusters",
#'     reduction = "UMAP",
#'     label = TRUE,
#'     theme_use = "theme_blank"
#'   ) +
#'     trajectory
#'   FeatureDimPlot(
#'     pancreas_sub,
#'     features = "Monocle3_Pseudotime",
#'     reduction = "UMAP",
#'     theme_use = "theme_blank"
#'   ) +
#'     trajectory
#'
#'   if (FALSE) {
#'     # Select the lineage using monocle3::choose_graph_segments
#'     cds <- pancreas_sub@tools$Monocle3$cds
#'     cds_sub <- thisutils::get_namespace_fun(
#'       "monocle3", "choose_graph_segments"
#'     )(
#'       cds,
#'       starting_pr_node = NULL,
#'       ending_pr_nodes = NULL
#'     )
#'     pancreas_sub$Lineages_1 <- NA
#'     pancreas_sub$Lineages_1[colnames(
#'       cds_sub
#'     )] <- pancreas_sub$Monocle3_Pseudotime[colnames(cds_sub)]
#'     CellDimPlot(
#'       pancreas_sub,
#'       group.by = "SubCellType",
#'       lineages = "Lineages_1",
#'       lineages_span = 0.1,
#'       theme_use = "theme_blank"
#'     )
#'   }
#'
#'   # Use Seurat clusters to infer the trajectories
#'   pancreas_sub <- RunStandardWorkflow(pancreas_sub)
#'   CellDimPlot(
#'     pancreas_sub,
#'     group.by = c("Standardclusters", "CellType"),
#'     label = TRUE,
#'     theme_use = "theme_blank"
#'   )
#'
#'   pancreas_sub <- RunMonocle3(
#'     pancreas_sub,
#'     clusters = "Standardclusters"
#'   )
#'
#'   trajectory <- pancreas_sub@tools$Monocle3$trajectory
#'   CellDimPlot(
#'     pancreas_sub,
#'     group.by = "Monocle3_partitions",
#'     reduction = "StandardUMAP2D",
#'     label = TRUE, theme_use = "theme_blank"
#'   ) + trajectory
#'
#'   CellDimPlot(
#'     pancreas_sub,
#'     group.by = "Monocle3_clusters",
#'     reduction = "StandardUMAP2D",
#'     label = TRUE, theme_use = "theme_blank"
#'   ) + trajectory
#'   FeatureDimPlot(
#'     pancreas_sub,
#'     features = "Monocle3_Pseudotime",
#'     reduction = "StandardUMAP2D",
#'     theme_use = "theme_blank"
#'   ) + trajectory
#'
#'   # Use custom graphs and cell clusters to infer
#'   # the partitions and trajectories, respectively
#'   pancreas_sub <- RunStandardWorkflow(
#'     pancreas_sub,
#'     cluster_resolution = 5
#'   )
#'   CellDimPlot(
#'     pancreas_sub,
#'     group.by = c("Standardclusters", "CellType"),
#'     label = TRUE
#'   )
#'   pancreas_sub <- RunMonocle3(
#'     pancreas_sub,
#'     clusters = "Standardclusters",
#'     graph = "Standardpca_SNN"
#'   )
#'   trajectory <- pancreas_sub@tools$Monocle3$trajectory
#'   CellDimPlot(
#'     pancreas_sub,
#'     group.by = "Monocle3_partitions",
#'     reduction = "StandardUMAP2D",
#'     label = TRUE, theme_use = "theme_blank"
#'   ) + trajectory
#'   CellDimPlot(
#'     pancreas_sub,
#'     group.by = "Monocle3_clusters",
#'     reduction = "StandardUMAP2D",
#'     label = TRUE, theme_use = "theme_blank"
#'   ) + trajectory
#'   FeatureDimPlot(
#'     pancreas_sub,
#'     features = "Monocle3_Pseudotime",
#'     reduction = "StandardUMAP2D",
#'     theme_use = "theme_blank"
#'   ) + trajectory
#' }
RunMonocle3 <- function(
  object,
  group.by = NULL,
  assay = NULL,
  layer = "counts",
  reduction = NULL,
  clusters = NULL,
  graph = NULL,
  partition_qval = 0.05,
  k = 50,
  cluster_method = "louvain",
  num_iter = 2,
  resolution = NULL,
  use_partition = NULL,
  close_loop = TRUE,
  root_pr_nodes = NULL,
  root_cells = NULL,
  show_plot = FALSE,
  xlab = NULL,
  ylab = NULL,
  seed = 11,
  verbose = TRUE,
  srt = NULL
) {
  srt <- resolve_deprecated_srt(object, srt, missing(object))
  log_message("Run {.pkg monocle3}...", verbose = verbose)
  set.seed(seed)
  check_r("cole-trapnell-lab/monocle3", verbose = FALSE)

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  expr_matrix <- SeuratObject::as.sparse(
    GetAssayData5(
      srt,
      assay = assay,
      layer = layer
    )
  )
  p_data <- srt@meta.data
  f_data <- data.frame(
    gene_short_name = row.names(expr_matrix),
    row.names = row.names(expr_matrix)
  )
  cds <- get_namespace_fun("monocle3", "new_cell_data_set")(
    expression_data = expr_matrix,
    cell_metadata = p_data,
    gene_metadata = f_data
  )
  if (!"Size_Factor" %in% colnames(cds@colData)) {
    size_factor <- paste0("nCount_", assay)
    if (size_factor %in% colnames(srt@meta.data)) {
      cds[["Size_Factor"]] <- cds[[size_factor, drop = TRUE]]
    }
  }
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  SingleCellExperiment::reducedDims(cds)[["UMAP"]] <- SeuratObject::Embeddings(
    srt[[reduction]]
  )
  loadings <- SeuratObject::Loadings(
    object = srt[[reduction]]
  )
  if (length(loadings) > 0) {
    methods::slot(
      object = cds, name = "reduce_dim_aux"
    )[["gene_loadings"]] <- loadings
  }
  stdev <- SeuratObject::Stdev(object = srt[[reduction]])
  if (length(stdev) > 0) {
    methods::slot(
      object = cds, name = "reduce_dim_aux"
    )[["prop_var_expl"]] <- stdev
  }

  if (!is.null(clusters)) {
    if (!is.null(graph)) {
      g <- igraph::graph_from_adjacency_matrix(
        adjmatrix = srt[[graph]],
        weighted = TRUE
      )
      cluster_result <- list(
        g = g,
        relations = NULL,
        distMatrix = "matrix",
        coord = NULL,
        edge_links = NULL,
        optim_res = list(
          membership = as.integer(as.factor(srt[[clusters, drop = TRUE]])),
          modularity = NA_real_
        )
      )
      if (length(unique(cluster_result$optim_res$membership)) > 1) {
        cluster_graph_res <- get_namespace_fun(
          "monocle3", "compute_partitions"
        )(
          cluster_result$g,
          cluster_result$optim_res,
          partition_qval
        )
        partitions <- igraph::components(
          cluster_graph_res$cluster_g
        )$membership[cluster_result$optim_res$membership]
        partitions <- as.factor(partitions)
      } else {
        partitions <- rep(1, ncol(srt))
      }
      names(partitions) <- colnames(cds)
      cds@clusters[["UMAP"]] <- list(
        cluster_result = cluster_result,
        partitions = partitions,
        clusters = as.factor(srt[[clusters, drop = TRUE]])
      )
      cds[["clusters"]] <- cds[[clusters]]
      add_citation <- get_namespace_fun(
        "monocle3", "add_citation"
      )
      cds <- add_citation(cds, "clusters")
      cds <- add_citation(cds, "partitions")
    } else {
      cds <- get_namespace_fun("monocle3", "cluster_cells")(
        cds,
        reduction_method = "UMAP",
        partition_qval = partition_qval,
        k = k,
        cluster_method = cluster_method,
        num_iter = num_iter,
        resolution = resolution
      )
      cds[["clusters"]] <- cds@clusters[["UMAP"]]$clusters <- as.factor(srt[[
        clusters,
        drop = TRUE
      ]])
    }
  } else {
    cds <- get_namespace_fun("monocle3", "cluster_cells")(
      cds,
      reduction_method = "UMAP",
      partition_qval = partition_qval,
      k = k,
      cluster_method = cluster_method,
      num_iter = num_iter,
      resolution = resolution
    )
    cds[["clusters"]] <- cds@clusters[["UMAP"]]$clusters
  }
  srt[["Monocle3_clusters"]] <- cds@clusters[["UMAP"]]$clusters
  srt[["Monocle3_partitions"]] <- cds@clusters[["UMAP"]]$partitions

  if (is.null(use_partition)) {
    use_partition <- TRUE
  }
  cds <- get_namespace_fun("monocle3", "learn_graph")(
    cds = cds,
    use_partition = use_partition,
    close_loop = close_loop
  )

  reduced_dim_coords <- Matrix::t(cds@principal_graph_aux[["UMAP"]]$dp_mst)
  edge_df <- igraph::as_data_frame(
    cds@principal_graph[["UMAP"]]
  )
  edge_df[, c("x", "y")] <- reduced_dim_coords[edge_df[["from"]], 1:2]
  edge_df[, c("xend", "yend")] <- reduced_dim_coords[edge_df[["to"]], 1:2]
  mst_branch_nodes <- get_namespace_fun(
    "monocle3", "branch_nodes"
  )(cds, "UMAP")
  mst_leaf_nodes <- get_namespace_fun(
    "monocle3", "leaf_nodes"
  )(cds, "UMAP")
  mst_root_nodes <- get_namespace_fun(
    "monocle3", "root_nodes"
  )(cds, "UMAP")
  pps <- c(mst_branch_nodes, mst_leaf_nodes, mst_root_nodes)
  point_df <- data.frame(
    nodes = names(pps),
    x = reduced_dim_coords[pps, 1],
    y = reduced_dim_coords[pps, 2]
  )
  point_df[, "is_branch"] <- names(pps) %in% names(mst_branch_nodes)
  trajectory <- list(
    geom_segment(data = edge_df, aes(x = x, y = y, xend = xend, yend = yend))
  )
  milestones <- list(
    geom_point(
      data = point_df[point_df[["is_branch"]] == FALSE, , drop = FALSE],
      aes(x = x, y = y),
      shape = 21,
      color = "white",
      fill = "black",
      size = 3,
      stroke = 1
    ),
    geom_point(
      data = point_df[point_df[["is_branch"]] == TRUE, , drop = FALSE],
      aes(x = x, y = y),
      shape = 21,
      color = "white",
      fill = "red",
      size = 3,
      stroke = 1
    ),
    ggnewscale::new_scale_color(),
    ggrepel::geom_text_repel(
      data = point_df,
      aes(x = x, y = y, label = nodes, color = is_branch),
      fontface = "bold",
      min.segment.length = 0,
      point.size = 3,
      max.overlaps = 100,
      bg.color = "white",
      bg.r = 0.1,
      size = 3.5
    ),
    scale_color_manual(
      values = stats::setNames(c("red", "black"), nm = c(TRUE, FALSE))
    )
  )

  p1 <- CellDimPlot(
    srt,
    group.by = "Monocle3_partitions",
    reduction = reduction,
    label = TRUE,
    force = TRUE,
    xlab = xlab,
    ylab = ylab
  ) +
    trajectory +
    milestones
  p2 <- CellDimPlot(
    srt,
    group.by = "Monocle3_clusters",
    reduction = reduction,
    label = TRUE,
    force = TRUE,
    xlab = xlab,
    ylab = ylab
  ) +
    trajectory
  p3 <- p2 +
    milestones
  if (!is.null(group.by)) {
    p4 <- CellDimPlot(
      srt,
      group.by = group.by,
      reduction = reduction,
      label = TRUE,
      force = TRUE,
      xlab = xlab,
      ylab = ylab
    ) +
      trajectory
    p5 <- p4 +
      milestones
  }
  if (isTRUE(show_plot)) {
    tryCatch(
      {
        if (!is.null(group.by)) {
          print(p1 + p3 + p5)
        } else {
          print(p1 + p3)
        }
      },
      error = function(e) {
        log_message(
          "Failed to print trajectory plots: {.val {conditionMessage(e)}}",
          message_type = "warning",
          verbose = verbose
        )
      }
    )
  }


  cds <- get_namespace_fun("monocle3", "order_cells")(
    cds,
    root_pr_nodes = root_pr_nodes,
    root_cells = root_cells
  )
  pseudotime <- cds@principal_graph_aux[["UMAP"]]$pseudotime
  pseudotime[is.infinite(pseudotime)] <- NA
  srt[["Monocle3_Pseudotime"]] <- pseudotime
  srt@tools$Monocle3 <- list(
    cds = cds,
    trajectory = trajectory,
    milestones = milestones
  )

  p6 <- FeatureDimPlot(
    srt,
    features = "Monocle3_Pseudotime",
    reduction = reduction,
    xlab = xlab,
    ylab = ylab
  ) +
    theme(legend.position = "none") +
    trajectory
  if (isTRUE(show_plot)) {
    tryCatch(
      {
        if (!is.null(group.by)) {
          print((p1 + p2) / (p4 + p6))
        } else {
          print(p1 + p2 + p6)
        }
      },
      error = function(e) {
        log_message(
          "Failed to print pseudotime plots: {.val {conditionMessage(e)}}",
          message_type = "warning",
          verbose = verbose
        )
      }
    )
  }

  log_message(
    "{.pkg monocle3} completed",
    message_type = "success",
    verbose = verbose
  )
  return(srt)
}
