#' Cluster a spatial sketch and project labels to the full assay
#'
#' Use Seurat's sketch/PCA/clustering/ProjectData path. The selected assay must
#' match the complete selected image context. Full counts are preserved; scaled
#' expression is created only for the sketch. The result records sampled cells
#' separately from projected labels. These are expression clusters, not inferred
#' spatial domains, and sketch sampling is not spatially stratified.
#' A single observed sketch cluster projects as the same constant label; this
#' explicit one-class case does not imply evidence for multiple populations.
#'
#' @param object A Seurat object with one selected spatial assay/context.
#' @param assay,image,coord.cols See SpatialDataInfo.
#' @param ncells Maximum number of sampled observations (at least 3).
#' @param nfeatures Maximum number of variable genes.
#' @param npcs Number of PCA components; must be smaller than the sketch size
#'   and number of usable features.
#' @param method Seurat sketch method: LeverageScore or Uniform.
#' @param resolution Seurat clustering resolution.
#' @param prefix Prefix for the new assay, reductions, graphs, metadata and tools.
#'   Existing outputs with these names are rejected; choose a new prefix to rerun.
#' @param seed Random seed passed to sampling, PCA and clustering.
#' @param max_dense_gb Maximum estimated GB for the sketch scaled matrix alone.
#'   This is a preflight guard, not a bound on total process memory.
#' @param verbose Whether to print progress.
#' @return A Seurat object with prefix_cluster for sampled cells and
#'   prefix_projected for every cell in the selected full assay. The tools entry
#'   records selected cells, projected cells, source and exact output names.
#' @seealso ReadSpatialData, SpatialDataInfo
#' @export
RunSpatialSketch <- function(object, assay = NULL, image = NULL,
                              coord.cols = c("x", "y"), ncells = 50000L,
                              nfeatures = 2000L, npcs = 30L,
                              method = c("LeverageScore", "Uniform"),
                              resolution = 0.6, prefix = "SpatialSketch",
                              seed = 11L, max_dense_gb = 2, verbose = TRUE) {
  info <- SpatialDataInfo(object, assay, image, coord.cols)
  assay <- info$assay
  method <- match.arg(method)
  validate_scalar_string(prefix, "prefix")
  if (!grepl("^[A-Za-z][A-Za-z0-9]*$", prefix)) stop("prefix must be alphanumeric and start with a letter", call. = FALSE)
  ncells <- validate_scalar_integer(ncells, "ncells", min = 3L)
  nfeatures <- validate_scalar_integer(nfeatures, "nfeatures", min = 2L)
  npcs <- validate_scalar_integer(npcs, "npcs", min = 1L)
  if (length(max_dense_gb) != 1L || !is.numeric(max_dense_gb) || !is.finite(max_dense_gb) || max_dense_gb <= 0) {
    stop("max_dense_gb must be positive and finite", call. = FALSE)
  }
  if (length(resolution) != 1L || !is.numeric(resolution) || !is.finite(resolution) || resolution <= 0) {
    stop("resolution must be positive and finite", call. = FALSE)
  }
  if (!setequal(colnames(object[[assay]]), info$cells)) {
    stop("Sketch assay must match the selected image; subset to that context before running", call. = FALSE)
  }
  layers <- SeuratObject::Layers(object[[assay]])
  if (!"counts" %in% layers || any(grepl("^counts[.]", layers))) {
    stop("Sketch requires one explicit counts layer; resolve split layers before running", call. = FALSE)
  }
  sketch_assay <- paste0(prefix, "Assay")
  pca <- paste0(prefix, "PCA")
  full_pca <- paste0(prefix, "FullPCA")
  graph <- paste0(prefix, c("NN", "SNN"))
  cluster <- paste0(prefix, "_cluster")
  projected <- paste0(prefix, "_projected")
  score <- paste0(prefix, "_leverage")
  targets <- c(sketch_assay, pca, full_pca, graph, cluster, projected, paste0(projected, ".score"), score)
  if (any(targets %in% c(names(object), names(object[[]]))) || !is.null(object@tools[[prefix]])) {
    stop("Sketch outputs already exist; choose a new prefix", call. = FALSE)
  }
  ncells <- min(ncells, length(info$cells))
  nfeatures <- min(nfeatures, nrow(object[[assay]]))
  dense_gb <- 8 * as.double(ncells) * nfeatures / 1024^3
  if (dense_gb > max_dense_gb) stop("Sketch scaled matrix exceeds max_dense_gb; reduce ncells or nfeatures", call. = FALSE)
  if (npcs >= min(ncells, nfeatures)) stop("npcs must be smaller than sketch cells and features", call. = FALSE)
  default_assay <- SeuratObject::DefaultAssay(object)
  # Every output is local until full projections have been validated.
  work <- object
  work <- Seurat::NormalizeData(work, assay = assay, verbose = verbose)
  work <- Seurat::FindVariableFeatures(work, assay = assay, nfeatures = nfeatures, verbose = verbose)
  features <- SeuratObject::VariableFeatures(work[[assay]])
  if (length(features) <= npcs) stop("Too few variable features for npcs", call. = FALSE)
  work <- Seurat::SketchData(work, assay = assay, ncells = ncells,
    sketched.assay = sketch_assay, method = method, var.name = score,
    features = features, seed = seed, verbose = verbose)
  sampled <- colnames(work[[sketch_assay]])
  if (!length(sampled) || anyDuplicated(sampled) || !all(sampled %in% info$cells)) {
    stop("Seurat returned invalid sketch cell IDs", call. = FALSE)
  }
  work <- Seurat::ScaleData(work, assay = sketch_assay, features = features, verbose = verbose)
  work <- Seurat::RunPCA(work, assay = sketch_assay, features = features,
    npcs = npcs, reduction.name = pca, seed.use = seed, verbose = verbose)
  work <- Seurat::FindNeighbors(work, reduction = pca, dims = seq_len(npcs),
    k.param = min(20L, length(sampled) - 1L), graph.name = graph, verbose = verbose)
  work <- Seurat::FindClusters(work, graph.name = graph[2], resolution = resolution,
    cluster.name = cluster, random.seed = seed, verbose = verbose)
  sampled_labels <- as.character(work[[]][sampled, cluster, drop = TRUE])
  if (anyNA(sampled_labels) || any(!nzchar(sampled_labels))) stop("Sketch clusters are incomplete", call. = FALSE)
  single_class <- length(unique(sampled_labels)) == 1L
  work <- Seurat::ProjectData(work, assay = assay, sketched.assay = sketch_assay,
    sketched.reduction = pca, full.reduction = full_pca, dims = seq_len(npcs),
    refdata = if (single_class) NULL else stats::setNames(list(cluster), projected),
    k.weight = min(50L, length(sampled)),
    recompute.neighbors = TRUE, recompute.weights = TRUE, verbose = verbose)
  if (single_class) {
    work[[projected]] <- stats::setNames(rep(sampled_labels[1], length(info$cells)), info$cells)
  }
  labels <- work[[]][info$cells, projected, drop = TRUE]
  embedding <- SeuratObject::Embeddings(work[[full_pca]])
  if (length(labels) != length(info$cells) || anyNA(labels) ||
      !setequal(rownames(embedding), info$cells) || any(!is.finite(embedding))) stop("Incomplete full-data projection", call. = FALSE)
  # Seurat also writes generic clustering metadata and a projection cache.
  # Retain only this run's named outputs; preserve any caller-owned generic state.
  work@tools$TransferSketchLabels <- object@tools$TransferSketchLabels
  work[["seurat_clusters"]] <- object[[]][["seurat_clusters"]]
  SeuratObject::Idents(work) <- SeuratObject::Idents(object)
  work@tools[[prefix]] <- spatial_tag_coordinate_contract(list(
    method = "SeuratSketch", status = "completed", source = info,
    parameters = list(assay = assay, ncells = ncells, nfeatures = nfeatures,
      npcs = npcs, method = method, resolution = resolution, seed = seed,
      label_projection = if (single_class) "constant single observed cluster" else "Seurat weighted neighbors",
      estimated_sketch_dense_gb = dense_gb, backend = "Seurat",
      backend_version = as.character(utils::packageVersion("Seurat"))),
    sampled_cells = sampled, projected_cells = info$cells,
    outputs = list(assay = sketch_assay, pca = pca, full_pca = full_pca,
      cluster = cluster, projected = projected),
    summary = list(n_sampled = length(sampled), n_projected = length(labels))))
  SeuratObject::DefaultAssay(work) <- default_assay
  spatial_run_receipt(done = "Spatial sketch and full-data projection completed",
    scope = sprintf("%s sampled; %s projected observations", length(sampled), length(labels)),
    saved = paste0("tools[[\"", prefix, "\"]]"), verbose = verbose)
  work
}
