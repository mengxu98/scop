#' @title Run Monocle2 analysis
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams standard_scop
#' @inheritParams GroupHeatmap
#' @inheritParams CellDimPlot
#' @param expressionFamily The distribution family to use for modeling gene expression.
#' Default is `"negbinomial.size"`.
#' @param features A character vector of features to use.
#' Defaults to NULL, in which case features were determined by `feature_type`.
#' @param feature_type The type of features to use in the analysis.
#' Possible values are "HVF" for highly variable features or "Disp" for features selected based on dispersion.
#' Default is `"HVF"`.
#' @param disp_filter A string specifying the filter to use when `feature_type` is "Disp".
#' Default is `"mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit"`.
#' @param max_components The maximum number of dimensions to use for dimensionality reduction.
#' Default is `2`.
#' @param reduction_method The dimensionality reduction method to use.
#' Possible values are `"DDRTree"`, `"ICA"`, `"tSNE"`,
#' `"SimplePPT"`, `"L1-graph"`, `"SGL-tree"`.
#' Default is `"DDRTree"`.
#' @param norm_method The normalization method to use.
#' Possible values are `"log"` and `"none"`.
#' Default is `"log"`.
#' @param residualModelFormulaStr A model formula specifying the effects to subtract.
#' Default is NULL.
#' @param pseudo_expr Amount to increase expression values before dimensionality reduction.
#' Default is 1.
#' @param root_state The state to use as the root of the trajectory.
#' If NULL, the R backend prompts for user input, and the C++ backend prompts
#' in interactive sessions after initial ordering. In non-interactive C++ runs,
#' the first cell is used. For `backend = "cpp"`, `root_state` can also match a
#' C++ trajectory state id after initial ordering, or a `group.by` label when
#' `group.by` is provided.
#' @param backend Backend used to compute the trajectory. `"r"` keeps the
#' original Monocle2 workflow and remains the default. `"cpp"` keeps Monocle2
#' dimensional reduction, uses native C++ ordering for DDRTree, and falls back
#' to Monocle2 ordering for other reduction methods.
#' @param n_neighbors Deprecated compatibility parameter for the C++ backend.
#' The current C++ backend reuses Monocle2's learned minimum spanning tree and
#' ignores this value.
#' @param ddrtree_maxIter Optional maximum iteration count passed to
#' `DDRTree::DDRTree()` when `reduction_method = "DDRTree"`. Lower values can
#' speed up exploratory runs but may reduce agreement with the default Monocle2
#' trajectory.
#' @param ddrtree_ncenter Optional number of DDRTree centers. This can change
#' trajectory topology and is intended for advanced exploratory use.
#' @param ddrtree_tol Optional convergence tolerance passed to `DDRTree::DDRTree()`.
#' @param show_plot Whether to print diagnostic plots during the run.
#' Default is `FALSE`.
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
#'   pancreas_sub <- standard_scop(pancreas_sub)
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
  srt,
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
  backend = c("r", "cpp"),
  n_neighbors = 30,
  ddrtree_maxIter = NULL,
  ddrtree_ncenter = NULL,
  ddrtree_tol = NULL,
  show_plot = FALSE,
  xlab = NULL,
  ylab = NULL,
  seed = 11,
  verbose = TRUE
) {
  backend <- match.arg(backend)
  if (identical(backend, "cpp")) {
    return(run_monocle2_cpp(
      srt = srt,
      assay = assay,
      layer = layer,
      group.by = group.by,
      expressionFamily = expressionFamily,
      features = features,
      feature_type = feature_type,
      disp_filter = disp_filter,
      max_components = max_components,
      reduction_method = reduction_method,
      norm_method = norm_method,
      residualModelFormulaStr = residualModelFormulaStr,
      pseudo_expr = pseudo_expr,
      root_state = root_state,
      n_neighbors = n_neighbors,
      ddrtree_maxIter = ddrtree_maxIter,
      ddrtree_ncenter = ddrtree_ncenter,
      ddrtree_tol = ddrtree_tol,
      show_plot = show_plot,
      xlab = xlab,
      ylab = ylab,
      seed = seed,
      verbose = verbose
    ))
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
  uses_negbinomial <- any(c("negbinomial", "negbinomial.size") %in% expressionFamily)
  if (isTRUE(uses_negbinomial)) {
    cds <- get_namespace_fun("BiocGenerics", "estimateSizeFactors")(cds)
    cds <- get_namespace_fun("BiocGenerics", "estimateDispersions")(cds)
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
          Seurat::FindVariableFeatures(srt, assay = assay, verbose = FALSE),
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
      print(p1 + p2)
    } else {
      print(p1)
    }
  }

  if (is.null(root_state)) {
    if (interactive()) {
      root_state <- utils::select.list(
        sort(unique(cds[["State"]])),
        title = "Select the root state to order cells:"
      )
      if (root_state == "" || length(root_state) == 0) {
        root_state <- NULL
      }
    } else {
      root_state <- sort(unique(cds[["State"]]))[1]
    }
  }
  cds <- get_namespace_fun("monocle", "orderCells")(
    cds, root_state = root_state
  )
  srt[["Monocle2_State"]] <- cds[["State"]]
  srt[["Monocle2_Pseudotime"]] <- cds[["Pseudotime"]]
  srt@tools$Monocle2 <- list(
    cds = cds,
    features = features,
    trajectory = trajectory
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
  log_message(
    "{.pkg monocle2} completed",
    message_type = "success",
    verbose = verbose
  )

  return(srt)
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

run_monocle2_cpp <- function(
  srt,
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
  n_neighbors = 30,
  ddrtree_maxIter = NULL,
  ddrtree_ncenter = NULL,
  ddrtree_tol = NULL,
  show_plot = FALSE,
  xlab = NULL,
  ylab = NULL,
  seed = 11,
  verbose = TRUE
) {
  log_message(
    "Run {.pkg monocle2} hybrid cpp backend...",
    verbose = verbose
  )
  set.seed(seed)
  check_r(
    c("mengxu98/monocle", "DDRTree", "BiocGenerics", "Biobase", "VGAM"),
    verbose = FALSE
  )

  for (pkg in c("DDRTree", "Biobase")) {
    if (!paste0("package:", pkg) %in% search()) {
      suppressMessages(
        suppressWarnings(
          attachNamespace(pkg)
        )
      )
    }
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
  need_dispersion <- isTRUE(show_plot) || (is.null(features) && identical(feature_type, "Disp"))
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
          Seurat::FindVariableFeatures(srt, assay = assay, verbose = FALSE),
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

  embeddings <- Matrix::t(cds@reducedDimS)
  colnames(embeddings) <- paste0(
    cds@dim_reduce_type, "_", seq_len(ncol(embeddings))
  )
  srt[[cds@dim_reduce_type]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings,
    key = paste0(cds@dim_reduce_type, "_"),
    assay = assay
  )

  initialized_ordering <- monocle2_cpp_initialize_ordering(cds)
  cds <- initialized_ordering[["cds"]]
  initial_state <- Biobase::pData(cds)[["State"]]
  root_state_resolved <- monocle2_cpp_resolve_root_state(
    srt = srt,
    root_state = root_state,
    group.by = group.by,
    initial_state = initial_state
  )
  if (identical(cds@dim_reduce_type, "DDRTree")) {
    tree <- monocle2_cpp_order_cells(
      cds,
      root_state = root_state_resolved,
      initial_root_state = root_state_resolved,
      projection = initialized_ordering[["projection"]]
    )
    states <- factor(tree[["state"]])
    pseudotime <- as.numeric(tree[["pseudotime"]])
    engine <- "monocle_mst_cpp_order"
  } else {
    cds <- get_namespace_fun("monocle", "orderCells")(cds, root_state = root_state_resolved)
    states <- Biobase::pData(cds)[["State"]]
    pseudotime <- Biobase::pData(cds)[["Pseudotime"]]
    engine <- "monocle_orderCells_fallback"
  }
  srt[["Monocle2_State"]] <- states
  srt[["Monocle2_Pseudotime"]] <- pseudotime

  pdata <- Biobase::pData(cds)
  pdata[["State"]] <- states
  pdata[["Pseudotime"]] <- pseudotime
  Biobase::pData(cds) <- pdata

  if (cds@dim_reduce_type == "ICA") {
    reduced_dim_coords <- as.data.frame(Matrix::t(cds@reducedDimS))
  } else if (cds@dim_reduce_type %in% c("simplePPT", "DDRTree")) {
    reduced_dim_coords <- as.data.frame(Matrix::t(cds@reducedDimK))
  } else {
    reduced_dim_coords <- as.data.frame(Matrix::t(cds@reducedDimS))
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
      print(p1 + p2)
    } else {
      print(p1)
    }
  }

  srt@tools$Monocle2 <- list(
    cds = cds,
    features = features,
    trajectory = trajectory,
    graph = edge_df,
    dispersion_table = dispersion_table,
    backend = "cpp",
    parameters = list(
      backend = "cpp",
      engine = engine,
      reduction = cds@dim_reduce_type,
      n_neighbors = as.integer(n_neighbors),
      ddrtree_maxIter = ddrtree_maxIter,
      ddrtree_ncenter = ddrtree_ncenter,
      ddrtree_tol = ddrtree_tol,
      root_state = root_state_resolved,
      root_cell = names(which.min(pseudotime)),
      state_count = length(levels(states))
    )
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
  log_message(
    "{.pkg monocle2} hybrid cpp backend completed",
    message_type = "success",
    verbose = verbose
  )

  srt
}

monocle2_cpp_order_cells <- function(cds, root_state = NULL, initial_root_state = NULL, projection = NULL) {
  if (cds@dim_reduce_type != "DDRTree") {
    select_root_cell <- get_namespace_fun("monocle", "select_root_cell")
    root_cell <- select_root_cell(cds, root_state, reverse = NULL)
    return(monocle2_cpp_extract_ordering(cds, root_cell = root_cell))
  }

  old_mst <- get_namespace_fun("monocle", "minSpanningTree")(cds)
  if (is.null(projection)) {
    projection <- monocle2_cpp_project_cells_to_mst(cds)
  }
  root_cell <- monocle2_cpp_select_root_cell(
    cds = cds,
    root_state = root_state,
    projection = projection,
    initial_root_state = initial_root_state
  )

  root_cell_idx <- which(igraph::V(old_mst)$name == root_cell)
  cells_mapped_to_graph_root <- which(
    projection[["closest_vertex"]] == root_cell_idx
  )
  if (length(cells_mapped_to_graph_root) == 0L) {
    cells_mapped_to_graph_root <- root_cell_idx
  }
  cells_mapped_to_graph_root <- projection[["vertex_names"]][cells_mapped_to_graph_root]
  tip_leaves <- projection[["tip_leaves"]]
  projected_root_cell <- cells_mapped_to_graph_root[cells_mapped_to_graph_root %in% tip_leaves][1]
  if (is.na(projected_root_cell)) {
    projected_root_cell <- cells_mapped_to_graph_root[1]
  }

  tree <- monocle2_cpp_order_projection(projection, root_cell = projected_root_cell)
  tree[["state"]] <- Biobase::pData(cds)[["State"]]

  tree
}

monocle2_cpp_select_root_cell <- function(cds, root_state = NULL, projection = NULL, initial_root_state = NULL) {
  mst <- get_namespace_fun("monocle", "minSpanningTree")(cds)
  if (is.null(root_state)) {
    diameter <- igraph::get_diameter(mst)
    return(names(diameter[1]))
  }
  pdata <- Biobase::pData(cds)
  if (is.null(pdata[["State"]])) {
    stop("Error: State has not yet been set. Please call orderCells() without specifying root_state, then try this call again.")
  }
  root_cell_candidates <- which(as.character(pdata[["State"]]) == as.character(root_state))
  if (length(root_cell_candidates) == 0L) {
    stop(paste("Error: no cells for State =", root_state))
  }
  if (is.null(projection)) {
    projection <- monocle2_cpp_project_cells_to_mst(cds)
  }
  previous_root_cell <- cds@auxOrderingData[[cds@dim_reduce_type]][["root_cell"]]
  use_min_pseudotime <- !is.null(previous_root_cell) &&
    !is.na(match(previous_root_cell, rownames(pdata))) &&
    as.character(pdata[previous_root_cell, "State"]) == as.character(root_state)
  root_vertex_idx <- monocle2_select_root_by_state_cpp(
    coords = get_namespace_fun("monocle", "reducedDimS")(cds),
    candidate_cells = as.integer(root_cell_candidates),
    pseudotime = as.numeric(pdata[["Pseudotime"]]),
    closest_vertex = as.integer(projection[["closest_vertex"]]),
    use_min_pseudotime = isTRUE(use_min_pseudotime)
  )
  igraph::V(mst)[root_vertex_idx]$name
}

monocle2_cpp_initialize_ordering <- function(cds) {
  if (cds@dim_reduce_type != "DDRTree") {
    return(list(
      cds = get_namespace_fun("monocle", "orderCells")(cds),
      projection = NULL
    ))
  }

  select_root_cell <- get_namespace_fun("monocle", "select_root_cell")
  root_cell <- select_root_cell(cds, root_state = NULL, reverse = NULL)
  graph_ordering <- monocle2_cpp_extract_ordering(
    cds,
    root_cell = root_cell,
    reorder_to_cells = FALSE
  )

  old_mst <- get_namespace_fun("monocle", "minSpanningTree")(cds)
  projection <- monocle2_cpp_project_cells_to_mst(cds)

  root_cell_idx <- which(igraph::V(old_mst)$name == root_cell)
  cells_mapped_to_graph_root <- which(
    projection[["closest_vertex"]] == root_cell_idx
  )
  if (length(cells_mapped_to_graph_root) == 0L) {
    cells_mapped_to_graph_root <- root_cell_idx
  }
  cells_mapped_to_graph_root <- projection[["vertex_names"]][cells_mapped_to_graph_root]
  tip_leaves <- projection[["tip_leaves"]]
  projected_root_cell <- cells_mapped_to_graph_root[cells_mapped_to_graph_root %in% tip_leaves][1]
  if (is.na(projected_root_cell)) {
    projected_root_cell <- cells_mapped_to_graph_root[1]
  }

  projected_ordering <- monocle2_cpp_order_projection(projection, root_cell = projected_root_cell)
  closest_vertex <- projection[["closest_vertex"]]
  initial_state <- graph_ordering[["state"]][closest_vertex]

  pdata <- Biobase::pData(cds)
  pdata[["Pseudotime"]] <- as.numeric(projected_ordering[["pseudotime"]])
  pdata[["State"]] <- factor(initial_state)
  Biobase::pData(cds) <- pdata
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_tree <- projection[["graph"]]
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_dist <- projection[["projected"]]
  cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex <- matrix(
    closest_vertex,
    ncol = 1,
    dimnames = list(projection[["vertex_names"]], NULL)
  )
  cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- projected_root_cell

  list(cds = cds, projection = projection)
}

monocle2_cpp_project_cells_to_mst <- function(cds) {
  mst <- get_namespace_fun("monocle", "minSpanningTree")(cds)
  z <- get_namespace_fun("monocle", "reducedDimS")(cds)
  y <- get_namespace_fun("monocle", "reducedDimK")(cds)
  closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  if (is.null(closest_vertex)) {
    cds <- get_namespace_fun("monocle", "findNearestPointOnMST")(cds)
    closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
  }
  closest_vertex <- as.integer(closest_vertex[, 1])
  graph_edges <- igraph::as_edgelist(mst, names = FALSE)
  if (is.null(dim(graph_edges))) {
    graph_edges <- matrix(graph_edges, ncol = 2L)
  }
  projection <- monocle2_project_cells_to_mst_cpp(
    z = z,
    y = y,
    graph_edges = matrix(as.integer(graph_edges), ncol = 2L),
    closest_vertex = closest_vertex
  )
  vertex_names <- colnames(z)
  colnames(projection[["projected"]]) <- vertex_names
  rownames(projection[["projected"]]) <- rownames(z)
  projection[["vertex_names"]] <- vertex_names
  graph <- igraph::graph_from_edgelist(
    matrix(vertex_names[as.integer(projection[["edges"]])], ncol = 2L),
    directed = FALSE
  )
  igraph::E(graph)$weight <- projection[["weights"]]
  projection[["graph"]] <- graph
  projection[["tip_leaves"]] <- names(which(igraph::degree(graph) == 1))
  projection
}

monocle2_cpp_order_projection <- function(projection, root_cell) {
  root_idx <- match(root_cell, projection[["vertex_names"]])
  if (is.na(root_idx)) {
    root_idx <- 1L
  }
  tree <- monocle2_order_from_weighted_edges_cpp(
    n_cells = length(projection[["vertex_names"]]),
    edges = projection[["edges"]],
    weights = projection[["weights"]],
    root_cell = as.integer(root_idx)
  )
  names(tree[["pseudotime"]]) <- projection[["vertex_names"]]
  names(tree[["state"]]) <- projection[["vertex_names"]]
  tree
}

monocle2_cpp_extract_ordering <- function(cds, root_cell, reorder_to_cells = TRUE) {
  dp <- get_namespace_fun("monocle", "cellPairwiseDistances")(cds)
  mst <- get_namespace_fun("monocle", "minSpanningTree")(cds)
  edges <- igraph::as_edgelist(mst, names = FALSE)
  if (is.null(dim(edges))) {
    edges <- matrix(edges, ncol = 2L)
  }
  root_idx <- match(root_cell, igraph::V(mst)$name)
  if (is.na(root_idx)) {
    root_idx <- 1L
  }
  vertex_names <- igraph::V(mst)$name
  dp_mat <- if (is.matrix(dp)) dp else as.matrix(dp)
  if (!identical(rownames(dp_mat), vertex_names) || !identical(colnames(dp_mat), vertex_names)) {
    dp_mat <- dp_mat[vertex_names, vertex_names, drop = FALSE]
  }
  tree <- monocle2_order_from_mst_cpp(
    distances = dp_mat,
    edges = matrix(as.integer(edges), ncol = 2L),
    root_cell = as.integer(root_idx)
  )
  names(tree[["pseudotime"]]) <- vertex_names
  names(tree[["state"]]) <- vertex_names
  if (isTRUE(reorder_to_cells)) {
    tree[["pseudotime"]] <- tree[["pseudotime"]][rownames(Biobase::pData(cds))]
    tree[["state"]] <- tree[["state"]][rownames(Biobase::pData(cds))]
  }
  tree
}

monocle2_cpp_resolve_root_state <- function(srt, root_state = NULL, group.by = NULL, initial_state = NULL) {
  if (is.null(root_state) && !is.null(initial_state) && interactive()) {
    root_state <- utils::select.list(
      sort(unique(as.character(initial_state))),
      title = "Select the root state to order cells:"
    )
    if (root_state == "" || length(root_state) == 0L) {
      root_state <- NULL
    }
  }
  if (!is.null(root_state) && !is.null(group.by) && group.by %in% colnames(srt@meta.data)) {
    group_cells <- which(as.character(srt@meta.data[[group.by]]) %in% as.character(root_state))
    if (length(group_cells) > 0L && !is.null(initial_state)) {
      state_table <- sort(table(as.character(initial_state[group_cells])), decreasing = TRUE)
      if (length(state_table) > 0L) {
        return(names(state_table)[1L])
      }
    }
  }
  root_state
}

#' @title Run Monocle3 analysis
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams standard_scop
#' @inheritParams FeatureDimPlot
#' @inheritParams GroupHeatmap
#' @inheritParams CellDimPlot
#' @param clusters The cluster variable in the Seurat object to use for analysis.
#' Defaults to NULL, in which case use Monocle clusters is used.
#' @param graph The name of the graph slot in the Seurat object to use for analysis.
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
#' If not specified, user will be prompted for input.
#' Defaults to NULL.
#' @param close_loop Whether to close loops in the graph.
#' Defaults to TRUE.
#' @param root_pr_nodes The root nodes to order cells.
#' If not specified, user will be prompted for input.
#' Defaults to NULL.
#' @param root_cells The root cells to order cells.
#' If not specified, user will be prompted for input.
#' Defaults to NULL.
#' @param show_plot Whether to print diagnostic plots during the run.
#' Default is `FALSE`.
#'
#' @export
#'
#' @examples
#' if (interactive()) {
#'   data(pancreas_sub)
#'   pancreas_sub <- standard_scop(pancreas_sub)
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
#'   pancreas_sub <- standard_scop(pancreas_sub)
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
  srt,
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
  verbose = TRUE
) {
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
    if (interactive()) {
      use_partition <- utils::select.list(
        c(TRUE, FALSE),
        title = "Whether to use partitions to learn disjoint graph in each partition?"
      )
      if (use_partition == "" || length(use_partition) == 0) {
        use_partition <- TRUE
      }
    } else {
      use_partition <- TRUE
    }
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

  if (is.null(root_pr_nodes) && is.null(root_cells)) {
    if (interactive()) {
      root_pr_nodes <- utils::select.list(
        names(pps),
        title = "Select the root nodes to order cells, or leave blank for interactive selection:",
        multiple = TRUE
      )
      if (root_pr_nodes == "" || length(root_pr_nodes) == 0) {
        root_pr_nodes <- NULL
      }
    } else {
      root_pr_nodes <- names(pps)[1]
    }
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
