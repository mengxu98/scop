#' Run Monocle2 analysis
#'
#' Runs the Monocle2 algorithm on a Seurat object.
#'
#' @param srt A Seurat object.
#' @param assay The name of the assay in the Seurat object to use for analysis. Defaults to NULL, in which case the default assay of the object is used.
#' @param layer The layer in the Seurat object to use for analysis. Default is "counts".
#' @param expressionFamily The distribution family to use for modeling gene expression. Default is "negbinomial.size".
#' @param features A vector of gene names or indices specifying the features to use in the analysis. Defaults to NULL, in which case features were determined by \code{feature_type}.
#' @param feature_type The type of features to use in the analysis. Possible values are "HVF" for highly variable features
#'   or "Disp" for features selected based on dispersion. Default is "HVF".
#' @param disp_filter A string specifying the filter to use when \code{feature_type} is "Disp". Default is
#'   "mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit".
#' @param max_components The maximum number of dimensions to use for dimensionality reduction. Default is 2.
#' @param reduction_method The dimensionality reduction method to use. Possible values are "DDRTree" and "UMAP". Default is "DDRTree".
#' @param norm_method The normalization method to use. Possible values are "log" and "none". Default is "log".
#' @param residualModelFormulaStr A model formula specifying the effects to subtract. Default is NULL.
#' @param pseudo_expr Amount to increase expression values before dimensionality reduction. Default is 1.
#' @param root_state The state to use as the root of the trajectory. If NULL, will prompt for user input.
#' @param seed An integer specifying the random seed to use. Default is 11.
#'
#' @examples
#' if (interactive()) {
#'   data("pancreas_sub")
#'   pancreas_sub <- RunMonocle2(srt = pancreas_sub)
#'   names(pancreas_sub@tools$Monocle2)
#'   trajectory <- pancreas_sub@tools$Monocle2$trajectory
#'
#'   CellDimPlot(pancreas_sub, group.by = "Monocle2_State", reduction = "DDRTree", label = TRUE, theme_use = "theme_blank") + trajectory
#'   CellDimPlot(pancreas_sub, group.by = "Monocle2_State", reduction = "UMAP", label = TRUE, theme_use = "theme_blank")
#'   FeatureDimPlot(pancreas_sub, features = "Monocle2_Pseudotime", reduction = "UMAP", theme_use = "theme_blank")
#'
#'   pancreas_sub <- RunMonocle2(
#'     srt = pancreas_sub,
#'     feature_type = "Disp", disp_filter = "mean_expression >= 0.01 & dispersion_empirical >= 1 * dispersion_fit"
#'   )
#'   trajectory <- pancreas_sub@tools$Monocle2$trajectory
#'   CellDimPlot(pancreas_sub, group.by = "Monocle2_State", reduction = "DDRTree", label = TRUE, theme_use = "theme_blank") + trajectory
#'   CellDimPlot(pancreas_sub, group.by = "Monocle2_State", reduction = "UMAP", label = TRUE, theme_use = "theme_blank")
#'   FeatureDimPlot(pancreas_sub, features = "Monocle2_Pseudotime", reduction = "UMAP", theme_use = "theme_blank")
#' }
#' @importFrom Seurat DefaultAssay CreateDimReducObject GetAssayData VariableFeatures FindVariableFeatures
#' @importFrom SeuratObject as.sparse
#' @importFrom igraph as_data_frame
#' @importFrom ggplot2 geom_segment
#' @export
RunMonocle2 <- function(
    srt,
    assay = NULL,
    layer = "counts",
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
    seed = 11) {
  set.seed(seed)
  check_r(c("monocle", "DDRTree", "BiocGenerics", "Biobase", "VGAM"))

  if (!"package:DDRTree" %in% search()) {
    attachNamespace("DDRTree")
  }

  assay <- assay %||% DefaultAssay(srt)
  expr_matrix <- SeuratObject::as.sparse(
    Seurat::GetAssayData(srt, assay = assay, layer = layer)
  )
  p_data <- srt@meta.data
  f_data <- data.frame(
    gene_short_name = row.names(expr_matrix),
    row.names = row.names(expr_matrix)
  )
  pd <- new("AnnotatedDataFrame", data = p_data)
  fd <- new("AnnotatedDataFrame", data = f_data)
  cds <- monocle::newCellDataSet(
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
  if (any(c("negbinomial", "negbinomial.size") %in% expressionFamily)) {
    cds <- BiocGenerics::estimateSizeFactors(cds)
    cds <- BiocGenerics::estimateDispersions(cds)
  }
  if (is.null(features)) {
    if (feature_type == "HVF") {
      features <- VariableFeatures(srt, assay = assay)
      if (length(features) == 0) {
        features <- VariableFeatures(
          FindVariableFeatures(srt, assay = assay),
          assay = assay
        )
      }
    }
    if (feature_type == "Disp") {
      features <- subset(
        monocle::dispersionTable(cds),
        eval(rlang::parse_expr(disp_filter))
      )$gene_id
    }
  }
  message("features number: ", length(features))
  cds <- monocle::setOrderingFilter(cds, features)
  p <- monocle::plot_ordering_genes(cds)
  print(p)

  cds <- monocle::reduceDimension(
    cds = cds,
    max_components = max_components,
    reduction_method = reduction_method,
    norm_method = norm_method,
    residualModelFormulaStr = residualModelFormulaStr,
    pseudo_expr = pseudo_expr
  )
  cds <- orderCells(cds)

  embeddings <- t(cds@reducedDimS)
  colnames(embeddings) <- paste0(cds@dim_reduce_type, "_", 1:ncol(embeddings))
  srt[[cds@dim_reduce_type]] <- CreateDimReducObject(
    embeddings = embeddings,
    key = paste0(cds@dim_reduce_type, "_"),
    assay = assay
  )
  srt[["Monocle2_State"]] <- cds[["State"]]
  if (cds@dim_reduce_type == "ICA") {
    reduced_dim_coords <- as.data.frame(t(cds@reducedDimS))
  } else if (cds@dim_reduce_type %in% c("simplePPT", "DDRTree")) {
    reduced_dim_coords <- as.data.frame(t(cds@reducedDimK))
  }
  edge_df <- as_data_frame(cds@minSpanningTree)
  edge_df[, c("x", "y")] <- reduced_dim_coords[edge_df[["from"]], 1:2]
  edge_df[, c("xend", "yend")] <- reduced_dim_coords[edge_df[["to"]], 1:2]
  trajectory <- geom_segment(
    data = edge_df,
    aes(x = x, y = y, xend = xend, yend = yend)
  )
  p <- CellDimPlot(
    srt,
    group.by = "Monocle2_State",
    reduction = reduction_method,
    label = TRUE,
    force = TRUE
  ) +
    trajectory
  print(p)
  if (is.null(root_state)) {
    root_state <- utils::select.list(
      sort(unique(cds[["State"]])),
      title = "Select the root state to order cells:"
    )
    if (root_state == "" || length(root_state) == 0) {
      root_state <- NULL
    }
  }
  cds <- orderCells(cds, root_state = root_state)
  srt[["Monocle2_State"]] <- cds[["State"]]
  srt[["Monocle2_Pseudotime"]] <- cds[["Pseudotime"]]
  srt@tools$Monocle2 <- list(
    cds = cds,
    features = features,
    trajectory = trajectory
  )

  p1 <- CellDimPlot(
    srt,
    group.by = "Monocle2_State",
    reduction = reduction_method,
    label = TRUE,
    force = TRUE
  ) +
    trajectory
  p2 <- FeatureDimPlot(
    srt,
    features = "Monocle2_Pseudotime",
    reduction = reduction_method
  ) +
    trajectory
  print(p1)
  Sys.sleep(1)
  print(p2)
  return(srt)
}

orderCells <- function(
    cds,
    root_state = NULL,
    num_paths = NULL,
    reverse = NULL) {
  if (class(cds)[1] != "CellDataSet") {
    stop("Error cds is not of type 'CellDataSet'")
  }
  if (is.null(cds@dim_reduce_type)) {
    stop(
      "Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function."
    )
  }
  if (any(c(length(cds@reducedDimS) == 0, length(cds@reducedDimK) == 0))) {
    stop(
      "Error: dimension reduction didn't prodvide correct results. Please check your reduceDimension() step and ensure correct dimension reduction are performed before calling this function."
    )
  }
  root_cell <- monocle:::select_root_cell(cds, root_state, reverse)
  cds@auxOrderingData <- new.env(hash = TRUE)
  if (cds@dim_reduce_type == "ICA") {
    if (is.null(num_paths)) {
      num_paths <- 1
    }
    adjusted_S <- t(cds@reducedDimS)
    dp <- Matrix::as.matrix(stats::dist(adjusted_S))
    cellPairwiseDistances(cds) <- Matrix::as.matrix(stats::dist(adjusted_S))
    gp <- igraph::graph.adjacency(dp, mode = "undirected", weighted = TRUE)
    dp_mst <- igraph::minimum.spanning.tree(gp)
    monocle::minSpanningTree(cds) <- dp_mst
    next_node <<- 0
    res <- monocle:::pq_helper(
      dp_mst,
      use_weights = FALSE,
      root_node = root_cell
    )
    cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
    order_list <- monocle:::extract_good_branched_ordering(
      res$subtree,
      res$root,
      monocle::cellPairwiseDistances(cds),
      num_paths,
      FALSE
    )
    cc_ordering <- order_list$ordering_df
    row.names(cc_ordering) <- cc_ordering$sample_name
    monocle::minSpanningTree(cds) <- igraph::as.undirected(
      order_list$cell_ordering_tree
    )
    cds[["Pseudotime"]] <- cc_ordering[
      row.names(Biobase::pData(cds)),
    ]$pseudo_time
    cds[["State"]] <- cc_ordering[row.names(Biobase::pData(cds)), ]$cell_state
    mst_branch_nodes <- igraph::V(
      monocle::minSpanningTree(cds)
    )[which(igraph::degree(monocle::minSpanningTree(cds)) > 2)]$name
    monocle::minSpanningTree(cds) <- dp_mst
    cds@auxOrderingData[[
      cds@dim_reduce_type
    ]]$cell_ordering_tree <- igraph::as.undirected(
      order_list$cell_ordering_tree
    )
  } else if (cds@dim_reduce_type == "DDRTree") {
    if (is.null(num_paths) == FALSE) {
      message(
        "Warning: num_paths only valid for method 'ICA' in reduceDimension()"
      )
    }
    cc_ordering <- extract_ddrtree_ordering(cds, root_cell)
    cds[["Pseudotime"]] <- cc_ordering[
      row.names(Biobase::pData(cds)),
    ]$pseudo_time
    K_old <- monocle::reducedDimK(cds)
    old_dp <- monocle::cellPairwiseDistances(cds)
    old_mst <- monocle::minSpanningTree(cds)
    old_A <- monocle::reducedDimA(cds)
    old_W <- monocle::reducedDimW(cds)
    cds <- project2MST(cds, monocle:::project_point_to_line_segment)
    monocle::minSpanningTree(cds) <- cds@auxOrderingData[[
      cds@dim_reduce_type
    ]]$pr_graph_cell_proj_tree
    root_cell_idx <- which(
      igraph::V(old_mst)$name == root_cell,
      arr.ind = TRUE
    )
    cells_mapped_to_graph_root <- which(
      cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex ==
        root_cell_idx
    )
    if (length(cells_mapped_to_graph_root) == 0) {
      cells_mapped_to_graph_root <- root_cell_idx
    }
    cells_mapped_to_graph_root <- igraph::V(
      monocle::minSpanningTree(cds)
    )[cells_mapped_to_graph_root]$name
    tip_leaves <- names(
      which(
        igraph::degree(
          monocle::minSpanningTree(cds)
        ) ==
          1
      )
    )
    root_cell <- cells_mapped_to_graph_root[
      cells_mapped_to_graph_root %in% tip_leaves
    ][1]
    if (is.na(root_cell)) {
      root_cell <- monocle:::select_root_cell(cds, root_state, reverse)
    }
    cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
    cc_ordering_new_pseudotime <- extract_ddrtree_ordering(cds, root_cell)
    cds[["Pseudotime"]] <- cc_ordering_new_pseudotime[
      row.names(Biobase::pData(cds)),
    ]$pseudo_time
    if (is.null(root_state) == TRUE) {
      closest_vertex <- cds@auxOrderingData[[
        "DDRTree"
      ]]$pr_graph_cell_proj_closest_vertex
      cds[["State"]] <- cc_ordering[closest_vertex[, 1], ]$cell_state
    }
    cds@reducedDimK <- K_old
    cds@cellPairwiseDistances <- old_dp
    cds@minSpanningTree <- old_mst
    cds@reducedDimA <- old_A
    cds@reducedDimW <- old_W
    mst_branch_nodes <- igraph::V(
      monocle::minSpanningTree(cds)
    )[which(igraph::degree(monocle::minSpanningTree(cds)) > 2)]$name
  } else if (cds@dim_reduce_type == "SimplePPT") {
    if (is.null(num_paths) == FALSE) {
      message(
        "Warning: num_paths only valid for method 'ICA' in reduceDimension()"
      )
    }
    cc_ordering <- extract_ddrtree_ordering(cds, root_cell)
    cds[["Pseudotime"]] <- cc_ordering[
      row.names(Biobase::pData(cds)),
    ]$pseudo_time
    cds[["State"]] <- cc_ordering[row.names(Biobase::pData(cds)), ]$cell_state
    mst_branch_nodes <- igraph::V(
      monocle::minSpanningTree(cds)
    )[which(igraph::degree(monocle::minSpanningTree(cds)) > 2)]$name
  }
  cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points <- mst_branch_nodes
  cds
}
project2MST <- function(cds, Projection_Method) {
  dp_mst <- monocle::minSpanningTree(cds)
  Z <- monocle::reducedDimS(cds)
  Y <- monocle::reducedDimK(cds)
  cds <- monocle:::findNearestPointOnMST(cds)
  closest_vertex <- cds@auxOrderingData[[
    "DDRTree"
  ]]$pr_graph_cell_proj_closest_vertex
  closest_vertex_names <- colnames(Y)[closest_vertex]
  closest_vertex_df <- Matrix::as.matrix(closest_vertex)
  row.names(closest_vertex_df) <- row.names(closest_vertex)
  tip_leaves <- names(which(igraph::degree(dp_mst) == 1))
  if (!is.function(Projection_Method)) {
    P <- Y[, closest_vertex]
  } else {
    P <- matrix(rep(0, length(Z)), nrow = nrow(Z))
    for (i in 1:length(closest_vertex)) {
      neighbors <- names(
        igraph::V(dp_mst)[suppressWarnings(nei(
          closest_vertex_names[i],
          mode = "all"
        ))]
      )
      projection <- NULL
      distance <- NULL
      Z_i <- Z[, i]
      for (neighbor in neighbors) {
        if (closest_vertex_names[i] %in% tip_leaves) {
          tmp <- monocle:::projPointOnLine(
            Z_i,
            Y[, c(closest_vertex_names[i], neighbor)]
          )
        } else {
          tmp <- Projection_Method(
            Z_i,
            Y[, c(closest_vertex_names[i], neighbor)]
          )
        }
        projection <- rbind(projection, tmp)
        distance <- c(distance, stats::dist(rbind(Z_i, tmp)))
      }
      if (!inherits(projection, "matrix")) {
        projection <- Matrix::as.matrix(projection)
      }
      P[, i] <- projection[which(distance == min(distance))[1], ]
    }
  }
  colnames(P) <- colnames(Z)
  dp <- Matrix::as.matrix(stats::dist(t(P)))
  min_dist <- min(dp[dp != 0])
  dp <- dp + min_dist
  diag(dp) <- 0
  monocle::cellPairwiseDistances(cds) <- dp
  gp <- igraph::graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- igraph::minimum.spanning.tree(gp)
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree <- dp_mst
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_dist <- P
  cds@auxOrderingData[[
    "DDRTree"
  ]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df
  cds
}

extract_ddrtree_ordering <- function(cds, root_cell, verbose = TRUE) {
  dp <- monocle::cellPairwiseDistances(cds)
  dp_mst <- monocle::minSpanningTree(cds)
  curr_state <- 1
  res <- list(subtree = dp_mst, root = root_cell)
  states <- rep(1, ncol(dp))
  names(states) <- igraph::V(dp_mst)$name
  pseudotimes <- rep(0, ncol(dp))
  names(pseudotimes) <- igraph::V(dp_mst)$name
  parents <- rep(NA, ncol(dp))
  names(parents) <- igraph::V(dp_mst)$name
  mst_traversal <- igraph::graph.dfs(
    dp_mst,
    root = root_cell,
    mode = "all",
    unreachable = FALSE,
    father = TRUE
  )
  mst_traversal$father <- as.numeric(mst_traversal$father)
  curr_state <- 1
  for (i in 1:length(mst_traversal$order)) {
    curr_node <- mst_traversal$order[i]
    curr_node_name <- igraph::V(dp_mst)[curr_node]$name
    if (is.na(mst_traversal$father[curr_node]) == FALSE) {
      parent_node <- mst_traversal$father[curr_node]
      parent_node_name <- igraph::V(dp_mst)[parent_node]$name
      parent_node_pseudotime <- pseudotimes[parent_node_name]
      parent_node_state <- states[parent_node_name]
      curr_node_pseudotime <- parent_node_pseudotime +
        dp[curr_node_name, parent_node_name]
      if (igraph::degree(dp_mst, v = parent_node_name) > 2) {
        curr_state <- curr_state + 1
      }
    } else {
      parent_node <- NA
      parent_node_name <- NA
      curr_node_pseudotime <- 0
    }
    curr_node_state <- curr_state
    pseudotimes[curr_node_name] <- curr_node_pseudotime
    states[curr_node_name] <- curr_node_state
    parents[curr_node_name] <- parent_node_name
  }
  ordering_df <- data.frame(
    sample_name = names(states),
    cell_state = factor(states),
    pseudo_time = as.vector(pseudotimes),
    parent = parents
  )
  row.names(ordering_df) <- ordering_df$sample_name
  return(ordering_df)
}

#' Run Monocle3 analysis
#'
#' Runs the Monocle3 algorithm on a Seurat object.
#'
#' @param srt A Seurat object.
#' @param assay The name of the assay in the Seurat object to use for analysis. Defaults to NULL, in which case the default assay of the object is used.
#' @param layer The layer in the Seurat object to use for analysis. Default is "counts".
#' @param reduction The reduction used. Defaults to NULL, in which case the default reduction of the Seurat object is used.
#' @param clusters The cluster variable in the Seurat object to use for analysis. Defaults to NULL, in which case use Monocle clusters is used.
#' @param graph The name of the graph slot in the Seurat object to use for analysis. Defaults to NULL, in which case Monocle graph is used.
#' @param partition_qval The q-value threshold for partitioning cells. Defaults to 0.05.
#' @param k The number of nearest neighbors to consider for clustering. Defaults to 50.
#' @param cluster_method The clustering method to use. Defaults to "louvain".
#' @param num_iter The number of iterations for clustering. Defaults to 2.
#' @param resolution The resolution parameter for clustering. Defaults to NULL.
#' @param use_partition Whether to use partitions to learn disjoint graph in each partition. If not specified, user will be prompted for input. Defaults to NULL.
#' @param close_loop Whether to close loops in the graph. Defaults to TRUE.
#' @param root_pr_nodes The root nodes to order cells. If not specified, user will be prompted for input. Defaults to NULL.
#' @param root_cells The root cells to order cells. If not specified, user will be prompted for input. Defaults to NULL.
#' @param seed The random seed to use for reproducibility. Defaults to 11.
#'
#' @examples
#' if (interactive()) {
#'   data("pancreas_sub")
#'   # Use Monocle clusters to infer the trajectories
#'   pancreas_sub <- RunMonocle3(
#'     srt = pancreas_sub,
#'     reduction = "UMAP"
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
#'   ## Select the lineage using monocle3::choose_graph_segments
#'   # cds <- pancreas_sub@tools$Monocle3$cds
#'   # cds_sub <- monocle3::choose_graph_segments(
#'   # cds,
#'   # starting_pr_node = NULL,
#'   # ending_pr_nodes = NULL
#'   # )
#'   # pancreas_sub$Lineages_1 <- NA
#'   # pancreas_sub$Lineages_1[colnames(cds_sub)] <- pancreas_sub$Monocle3_Pseudotime[colnames(cds_sub)]
#'   # CellDimPlot(
#'   # pancreas_sub,
#'   # group.by = "SubCellType",
#'   # lineages = "Lineages_1",
#'   # lineages_span = 0.1,
#'   # theme_use = "theme_blank"
#'   # )
#'
#'   # Use Seurat clusters to infer the trajectories
#'   pancreas_sub <- standard_scop(pancreas_sub)
#'   CellDimPlot(
#'     pancreas_sub,
#'     group.by = c("Standardclusters", "CellType"),
#'     label = TRUE,
#'     theme_use = "theme_blank"
#'   )
#'   pancreas_sub <- RunMonocle3(
#'     srt = pancreas_sub,
#'     clusters = "Standardclusters"
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
#'
#'   # Use custom graphs and cell clusters to infer the partitions and trajectories, respectively
#'   pancreas_sub <- standard_scop(
#'     pancreas_sub,
#'     cluster_resolution = 5
#'   )
#'   CellDimPlot(
#'     pancreas_sub,
#'     group.by = c("Standardclusters", "CellType"),
#'     label = TRUE
#'   )
#'   pancreas_sub <- RunMonocle3(
#'     srt = pancreas_sub,
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
#' @importFrom SeuratObject as.sparse Embeddings Loadings Stdev
#' @importFrom Seurat DefaultAssay GetAssayData
#' @importFrom igraph as_data_frame
#' @importFrom ggplot2 geom_segment
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggnewscale new_scale_color
#' @export
RunMonocle3 <- function(
    srt,
    assay = NULL,
    layer = "counts",
    reduction = DefaultReduction(srt),
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
    seed = 11) {
  set.seed(seed)
  if (
    !requireNamespace("monocle3", quietly = TRUE) ||
      utils::packageVersion("monocle3") < package_version("1.2.0")
  ) {
    check_r("cole-trapnell-lab/monocle3", force = TRUE)
  }
  assay <- assay %||% DefaultAssay(srt)
  expr_matrix <- SeuratObject::as.sparse(
    Seurat::GetAssayData(srt, assay = assay, layer = layer)
  )
  p_data <- srt@meta.data
  f_data <- data.frame(
    gene_short_name = row.names(expr_matrix),
    row.names = row.names(expr_matrix)
  )
  cds <- monocle3::new_cell_data_set(
    expression_data = expr_matrix,
    cell_metadata = p_data,
    gene_metadata = f_data
  )
  if (!"Size_Factor" %in% colnames(cds@colData)) {
    size.factor <- paste0("nCount_", assay)
    if (size.factor %in% colnames(srt@meta.data)) {
      cds[["Size_Factor"]] <- cds[[size.factor, drop = TRUE]]
    }
  }
  reduction <- reduction %||% DefaultReduction(srt)
  SingleCellExperiment::reducedDims(cds)[["UMAP"]] <- Embeddings(srt[[
    reduction
  ]])
  loadings <- Loadings(object = srt[[reduction]])
  if (length(loadings) > 0) {
    slot(object = cds, name = "reduce_dim_aux")[["gene_loadings"]] <- loadings
  }
  stdev <- Stdev(object = srt[[reduction]])
  if (length(stdev) > 0) {
    slot(object = cds, name = "reduce_dim_aux")[["prop_var_expl"]] <- stdev
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
        cluster_graph_res <- monocle3:::compute_partitions(
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
      cds <- monocle3:::add_citation(cds, "clusters")
      cds <- monocle3:::add_citation(cds, "partitions")
    } else {
      cds <- monocle3::cluster_cells(
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
    cds <- monocle3::cluster_cells(
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
  p1 <- CellDimPlot(
    srt,
    "Monocle3_partitions",
    reduction = reduction,
    label = FALSE,
    force = TRUE
  )
  p2 <- CellDimPlot(
    srt,
    "Monocle3_clusters",
    reduction = reduction,
    label = FALSE,
    force = TRUE
  )
  print(p1)
  Sys.sleep(1)
  print(p2)
  if (is.null(use_partition)) {
    use_partition <- utils::select.list(
      c(TRUE, FALSE),
      title = "Whether to use partitions to learn disjoint graph in each partition?"
    )
    if (use_partition == "" || length(use_partition) == 0) {
      use_partition <- TRUE
    }
  }
  cds <- monocle3::learn_graph(
    cds = cds,
    use_partition = use_partition,
    close_loop = close_loop
  )

  reduced_dim_coords <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst)
  edge_df <- as_data_frame(cds@principal_graph[["UMAP"]])
  edge_df[, c("x", "y")] <- reduced_dim_coords[edge_df[["from"]], 1:2]
  edge_df[, c("xend", "yend")] <- reduced_dim_coords[edge_df[["to"]], 1:2]
  mst_branch_nodes <- monocle3:::branch_nodes(cds, "UMAP")
  mst_leaf_nodes <- monocle3:::leaf_nodes(cds, "UMAP")
  mst_root_nodes <- monocle3:::root_nodes(cds, "UMAP")
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
    new_scale_color(),
    geom_text_repel(
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
  p <- CellDimPlot(
    srt,
    group.by = "Monocle3_partitions",
    reduction = reduction,
    label = FALSE,
    force = TRUE
  ) +
    trajectory +
    milestones
  print(p)

  if (is.null(root_pr_nodes) && is.null(root_cells)) {
    root_pr_nodes <- utils::select.list(
      names(pps),
      title = "Select the root nodes to order cells, or leave blank for interactive selection:",
      multiple = TRUE
    )
    if (root_pr_nodes == "" || length(root_pr_nodes) == 0) {
      root_pr_nodes <- NULL
    }
  }
  cds <- monocle3::order_cells(
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

  p1 <- CellDimPlot(
    srt,
    group.by = "Monocle3_partitions",
    reduction = reduction,
    label = FALSE,
    force = TRUE
  ) +
    trajectory
  p2 <- FeatureDimPlot(
    srt,
    features = "Monocle3_Pseudotime",
    reduction = reduction
  ) +
    theme(legend.position = "none") +
    trajectory
  print(p1)
  Sys.sleep(1)
  print(p2)
  return(srt)
}
