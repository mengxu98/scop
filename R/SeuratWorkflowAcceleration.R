#' Seurat clustering with transparent backend delegation
#'
#' Seurat already performs Louvain and SLM clustering in compiled code and
#' delegates Leiden clustering to its selected compiled backend. scop therefore
#' preserves the complete Seurat contract for this step.
#'
#' @param object A Seurat object or graph accepted by Seurat.
#' @param ... Arguments passed unchanged to [Seurat::FindClusters()].
#'
#' @return The value returned by [Seurat::FindClusters()].
#' @export
FindClusters <- function(object, ...) {
  Seurat::FindClusters(object = object, ...)
}

#' Seurat weighted-nearest-neighbor construction
#'
#' This wrapper preserves all Seurat multimodal-neighbor branches. The Seurat
#' implementation already uses optimized nearest-neighbor and SNN kernels, so
#' scop delegates rather than changing graph semantics.
#'
#' @param object A Seurat object.
#' @param ... Arguments passed unchanged to [Seurat::FindMultiModalNeighbors()].
#'
#' @return A Seurat object returned by [Seurat::FindMultiModalNeighbors()].
#' @export
FindMultiModalNeighbors <- function(object, ...) {
  Seurat::FindMultiModalNeighbors(object = object, ...)
}

module_score_context <- function(object, assay, slot) {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  if (!assay %in% names(object)) {
    return(NULL)
  }
  assay_object <- object[[assay]]
  if (inherits(assay_object, "Assay") && !inherits(assay_object, "StdAssay")) {
    if (!slot %in% methods::slotNames(assay_object)) {
      return(NULL)
    }
    data_use <- list(methods::slot(assay_object, slot))
  } else {
    layers <- SeuratObject::Layers(assay_object, search = slot)
    if (length(layers) == 0L) {
      return(NULL)
    }
    data_use <- lapply(layers, function(layer) {
      SeuratObject::LayerData(assay_object, layer = layer)
    })
  }
  if (!all(vapply(data_use, inherits, logical(1), what = "dgCMatrix"))) {
    return(NULL)
  }
  data_use <- lapply(data_use, function(data) {
    if (is.null(rownames(data))) {
      rownames(data) <- rownames(assay_object)
    }
    data
  })
  list(assay = assay, data = data_use)
}

module_score_supported <- function(
  context,
  features,
  pool,
  nbin,
  ctrl,
  k,
  search,
  dots
) {
  !is.null(context) &&
    is.list(features) &&
    length(features) > 0L &&
    all(vapply(features, is.character, logical(1))) &&
    all(vapply(features, length, integer(1)) > 0L) &&
    all(vapply(context$data, function(data) {
      all(vapply(features, function(x) all(x %in% rownames(data)), logical(1)))
    }, logical(1))) &&
    (is.null(pool) || (is.character(pool) && all(vapply(
      context$data,
      function(data) all(pool %in% rownames(data)),
      logical(1)
    )))) &&
    is.numeric(nbin) && length(nbin) == 1L && is.finite(nbin) &&
    nbin >= 1 && nbin == as.integer(nbin) &&
    is.numeric(ctrl) && length(ctrl) == 1L && is.finite(ctrl) && ctrl >= 1 &&
    identical(k, FALSE)
}

module_score_native <- function(
  object,
  context,
  features,
  pool,
  nbin,
  ctrl,
  name,
  seed
) {
  if (!is.null(seed)) {
    set.seed(seed = seed)
  }
  features <- lapply(features, unique)
  score_parts <- lapply(context$data, function(data) {
    pool_use <- pool %||% rownames(data)
    data_avg <- Matrix::rowMeans(data[pool_use, , drop = FALSE])
    data_avg <- data_avg[order(data_avg)]
    data_cut <- ggplot2::cut_number(
      data_avg + stats::rnorm(length(data_avg)) / 1e30,
      n = nbin,
      labels = FALSE,
      right = FALSE
    )
    names(data_cut) <- names(data_avg)

    control_sets <- vector("list", length(features))
    for (i in seq_along(features)) {
      for (feature in features[[i]]) {
        candidates <- data_cut[data_cut == data_cut[[feature]]]
        control_sets[[i]] <- c(
          control_sets[[i]],
          names(sample(candidates, size = ctrl, replace = FALSE))
        )
      }
      control_sets[[i]] <- unique(control_sets[[i]])
    }
    scores <- module_score_sparse(
      expr = data,
      feature_sets = lapply(features, match, table = rownames(data)),
      control_sets = lapply(control_sets, match, table = rownames(data))
    )
    rownames(scores) <- colnames(data)
    scores
  })
  scores <- do.call(rbind, score_parts)
  colnames(scores) <- paste0(name, seq_along(features))
  object[[colnames(scores)]] <- as.data.frame(scores)
  SeuratObject::CheckGC()
  object
}

#' Add expression-program scores with a sparse native fast path
#'
#' The default Seurat contract is retained. For sparse assay layers, complete
#' feature lists, and `k = FALSE`, scop reproduces Seurat's binning and control
#' sampling and fuses all score calculations in one sparse compiled pass.
#' `search` and extra arguments are inert once every requested feature is
#' present. Other branches delegate to Seurat unchanged.
#'
#' @param object A Seurat object.
#' @param features A list of feature vectors.
#' @param pool,nbin,ctrl,k,assay,name,seed,search,slot Seurat-compatible
#'   module-score parameters.
#' @param ... Additional arguments passed to Seurat on fallback.
#'
#' @return The input Seurat object with module-score metadata columns.
#' @export
AddModuleScore <- function(
  object,
  features,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = NULL,
  name = "Cluster",
  seed = 1,
  search = FALSE,
  slot = "data",
  ...
) {
  dots <- list(...)
  fallback <- function() {
    do.call(
      Seurat::AddModuleScore,
      c(
        list(
          object = object,
          features = features,
          pool = pool,
          nbin = nbin,
          ctrl = ctrl,
          k = k,
          assay = assay,
          name = name,
          seed = seed,
          search = search,
          slot = slot
        ),
        dots
      )
    )
  }
  if (!inherits(object, "Seurat")) {
    return(fallback())
  }
  context <- module_score_context(object, assay, slot)
  if (!module_score_supported(
    context = context,
    features = features,
    pool = pool,
    nbin = nbin,
    ctrl = ctrl,
    k = k,
    search = search,
    dots = dots
  )) {
    return(fallback())
  }
  module_score_native(
    object = object,
    context = context,
    features = features,
    pool = pool,
    nbin = as.integer(nbin),
    ctrl = as.integer(ctrl),
    name = name,
    seed = seed
  )
}

#' Cell-cycle scoring with native sparse module scoring
#'
#' This follows Seurat's cell-cycle classification contract while routing the
#' two module scores through [AddModuleScore()] when its native sparse contract
#' is met. Unsupported scoring branches transparently use Seurat.
#'
#' @param object A Seurat object.
#' @param s.features,g2m.features Feature vectors for S and G2M phases.
#' @param ctrl Number of control features; defaults to the smaller feature set.
#' @param set.ident Whether to set cell identities to the inferred phase.
#' @param ... Arguments passed to [AddModuleScore()].
#'
#' @return A Seurat object containing `S.Score`, `G2M.Score`, and `Phase`.
#' @export
CellCycleScoring <- function(
  object,
  s.features,
  g2m.features,
  ctrl = NULL,
  set.ident = FALSE,
  ...
) {
  if (is.null(ctrl)) {
    ctrl <- min(length(s.features), length(g2m.features))
  }
  if (any(grepl("Cell.Cycle", colnames(object[[]]), fixed = TRUE))) {
    return(Seurat::CellCycleScoring(
      object = object,
      s.features = s.features,
      g2m.features = g2m.features,
      ctrl = ctrl,
      set.ident = set.ident,
      ...
    ))
  }
  object <- AddModuleScore(
    object = object,
    features = list(S.Score = s.features, G2M.Score = g2m.features),
    name = "Cell.Cycle",
    ctrl = ctrl,
    ...
  )
  scores <- object[[c("Cell.Cycle1", "Cell.Cycle2")]]
  phase <- ifelse(
    scores[[1L]] < 0 & scores[[2L]] < 0,
    "G1",
    ifelse(
      scores[[1L]] == scores[[2L]],
      "Undecided",
      ifelse(scores[[1L]] > scores[[2L]], "S", "G2M")
    )
  )
  metadata <- data.frame(
    S.Score = scores[[1L]],
    G2M.Score = scores[[2L]],
    Phase = phase,
    row.names = rownames(scores)
  )
  object[[colnames(metadata)]] <- metadata
  if (set.ident) {
    object[["old.ident"]] <- SeuratObject::Idents(object)
    SeuratObject::Idents(object) <- "Phase"
  }
  object
}

#' Select integration features with exact vectorized ranking
#'
#' Uses a vectorized rank calculation when every object already has unique
#' variable features in its default assay. Other branches delegate to Seurat.
#'
#' @param object.list List of Seurat objects.
#' @param nfeatures Number of integration features to return.
#' @param assay Assay names, one per object.
#' @param verbose Show Seurat progress messages on delegated paths.
#' @param fvf.nfeatures Features requested if Seurat must calculate them.
#' @param ... Passed to Seurat on delegated paths.
#'
#' @return A character vector of integration features.
#' @export
SelectIntegrationFeatures <- function(
  object.list,
  nfeatures = 2000,
  assay = NULL,
  verbose = TRUE,
  fvf.nfeatures = 2000,
  ...
) {
  fallback <- function() {
    get("SelectIntegrationFeatures", asNamespace("Seurat"))(
      object.list = object.list,
      nfeatures = nfeatures,
      assay = assay,
      verbose = verbose,
      fvf.nfeatures = fvf.nfeatures,
      ...
    )
  }
  if (
    !length(object.list) ||
      !is.numeric(nfeatures) || length(nfeatures) != 1L ||
      !is.finite(nfeatures) || nfeatures < 1 ||
      nfeatures != as.integer(nfeatures) ||
      !requireNamespace("matrixStats", quietly = TRUE)
  ) {
    return(fallback())
  }

  defaults <- tryCatch(
    vapply(object.list, SeuratObject::DefaultAssay, character(1)),
    error = function(e) NULL
  )
  if (is.null(defaults)) {
    return(fallback())
  }
  if (is.null(assay)) {
    assay <- defaults
  } else if (
    !is.character(assay) || length(assay) != length(object.list) ||
      anyNA(assay) || !identical(as.character(assay), as.character(defaults))
  ) {
    return(fallback())
  }

  vf.list <- tryCatch(
    lapply(seq_along(object.list), function(i) {
      SeuratObject::VariableFeatures(object.list[[i]], assay = assay[[i]])
    }),
    error = function(e) NULL
  )
  if (
    is.null(vf.list) || any(lengths(vf.list) == 0L) ||
      any(vapply(vf.list, anyDuplicated, integer(1)) > 0L)
  ) {
    return(fallback())
  }

  frequency <- sort(
    table(unlist(vf.list, use.names = FALSE)),
    decreasing = TRUE
  )
  for (i in seq_along(object.list)) {
    genes <- tryCatch(rownames(object.list[[i]][[assay[[i]]]]), error = function(e) NULL)
    if (is.null(genes)) {
      return(fallback())
    }
    frequency <- frequency[names(frequency) %in% genes]
  }
  if (!length(frequency)) {
    return(fallback())
  }

  tie.value <- frequency[min(nfeatures, length(frequency))]
  above <- names(frequency[frequency > tie.value])
  tied <- names(frequency[frequency == tie.value])
  candidates <- c(above, tied)
  ranks <- vapply(
    vf.list,
    function(vf) match(candidates, vf),
    integer(length(candidates))
  )
  ranks <- matrix(ranks, nrow = length(candidates), ncol = length(vf.list))
  median.rank <- matrixStats::rowMedians(ranks, na.rm = TRUE)
  if (length(above)) {
    above <- above[order(median.rank[seq_along(above)])]
  }
  if (length(tied)) {
    tied.index <- length(above) + seq_along(tied)
    tied <- tied[order(median.rank[tied.index])]
  }
  c(above, utils::head(tied, nfeatures - length(above)))
}
