sct_fastrowscale <- function(mat_dense, do.scale, do.center, scale.max) {
  cm <- if (isTRUE(do.center)) {
    rowMeans(mat_dense)
  } else {
    rep(0, nrow(mat_dense))
  }
  mat_dense <- mat_dense - cm
  if (isTRUE(do.scale)) {
    n <- ncol(mat_dense)
    rv <- sqrt(pmax(rowSums(mat_dense^2) / pmax(n - 1, 1), 0))
    rv[rv == 0] <- 1
    mat_dense <- mat_dense / rv
    if (!is.na(scale.max)) {
      mat_dense[mat_dense > scale.max] <- scale.max
      mat_dense[mat_dense < -scale.max] <- -scale.max
    }
  }
  mat_dense
}

# Linear-model residualization + row scaling over feature x cell blocks,
# honoring split.by groups. Regression runs per group as in Seurat; with
# use.umi the residual space is mapped back through log1p after a per-gene
# shift, matching RegressOutMatrix.
sct_scale_general <- function(mat, # sparse or dense, features x all cells
                              latent_df, # data.frame aligned to colnames(mat) or NULL
                              split_levels, # factor over colnames(mat) or NULL
                              do.scale,
                              do.center,
                              scale.max,
                              use.umi = FALSE) {
  cells <- colnames(mat)
  groups <- if (is.null(split_levels)) {
    list(all = cells)
  } else {
    lapply(levels(droplevels(split_levels)), function(lvl) {
      cells[as.character(split_levels) == lvl]
    })
  }
  blocks <- lapply(groups, function(gcells) {
    sub <- mat[, gcells, drop = FALSE]
    y <- methods::as(sub, "matrix")
    storage.mode(y) <- "double"
    if (!is.null(latent_df)) {
      y <- sct_regress_out(y, latent_df[gcells, , drop = FALSE])
      if (isTRUE(use.umi)) {
        y <- log1p(y - apply(y, 1, min))
      }
    }
    sct_fastrowscale(y, do.scale, do.center, scale.max)
  })
  out <- do.call(cbind, blocks)[, cells, drop = FALSE]
  out
}

#' @export
ScaleData.Seurat <- function(
  object,
  features = NULL,
  assay = NULL,
  vars.to.regress = NULL,
  latent.data = NULL,
  split.by = NULL,
  model.use = "linear",
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  verbose = TRUE,
  ...
) {
  dots <- list(...)
  if (
    !identical(model.use, "linear") ||
      !all(vapply(
        list(use.umi, do.scale, do.center),
        function(x) is.logical(x) && length(x) == 1L && !is.na(x),
        logical(1)
      ))
  ) {
    return(do.call(
      utils::getFromNamespace("ScaleData.Seurat", "Seurat"),
      c(
        list(
          object = object,
          features = features,
          assay = assay,
          vars.to.regress = vars.to.regress,
          split.by = split.by,
          model.use = model.use,
          use.umi = use.umi,
          do.scale = do.scale,
          do.center = do.center,
          scale.max = scale.max,
          block.size = block.size,
          min.cells.to.block = min.cells.to.block,
          verbose = verbose
        ),
        dots
      )
    ))
  }
  assay_name <- if (is.null(assay)) {
    SeuratObject::DefaultAssay(object)
  } else {
    assay[1L]
  }
  assay_obj <- methods::slot(object, "assays")[[assay_name]]
  if (length(dots) > 0L || inherits(assay_obj, "SCTAssay")) {
    return(do.call(
      utils::getFromNamespace("ScaleData.Seurat", "Seurat"),
      c(
        list(
          object = object,
          features = features,
          assay = assay_name,
          vars.to.regress = vars.to.regress,
          split.by = split.by,
          model.use = model.use,
          use.umi = use.umi,
          do.scale = do.scale,
          do.center = do.center,
          scale.max = scale.max,
          block.size = block.size,
          min.cells.to.block = min.cells.to.block,
          verbose = verbose
        ),
        dots
      )
    ))
  }
  features <- if (is.null(features)) {
    SeuratObject::VariableFeatures(object)
  } else {
    features
  }
  if (length(features) == 0L) {
    features <- rownames(assay_obj)
  }
  all_genes <- rownames(assay_obj)
  features <- intersect(features, all_genes)
  features <- features[order(match(features, all_genes))]
  idx <- match(features, all_genes) - 1L

  # use.umi switches the source layer to counts, matching ScaleData.default.
  want_counts <- isTRUE(use.umi)
  if (inherits(assay_obj, "Assay") && !inherits(assay_obj, "StdAssay")) {
    src_layer <- if (want_counts) "counts" else "data"
    data_mat <- tryCatch(
      SeuratObject::GetAssayData(assay_obj, layer = src_layer),
      error = function(e) methods::slot(assay_obj, src_layer)
    )
    if (!inherits(data_mat, "dgCMatrix")) {
      data_mat <- tryCatch(
        methods::as(data_mat, "dgCMatrix"),
        error = function(e) NULL
      )
      if (is.null(data_mat)) {
        log_message(
          "ScaleData.Seurat requires the {src_layer} layer convertible to dgCMatrix.",
          message_type = "error"
        )
      }
    }
  } else {
    if (!inherits(assay_obj, "Assay5")) {
      log_message(
        "ScaleData.Seurat requires an Assay5 object for the layers path.",
        message_type = "error"
      )
    }
    prefix <- if (want_counts) "counts" else "data"
    if (length(SeuratObject::Layers(assay_obj, search = prefix)) == 0L) {
      log_message(
        "ScaleData.Seurat requires an Assay5 object with a {prefix} layer.",
        message_type = "error"
      )
    }
    src_layers <- SeuratObject::Layers(assay_obj, search = prefix)
    src_layers <- src_layers[grepl(paste0("^", prefix, "(\\.|$)"), src_layers)]
    if (length(src_layers) == 0L) {
      log_message(
        "ScaleData.Seurat requires an Assay5 object with a {prefix} layer.",
        message_type = "error"
      )
    }
    data_mats <- lapply(src_layers, function(src_layer) {
      SeuratObject::LayerData(
        assay_obj,
        layer = src_layer,
        features = features
      )
    })
    data_mat <- if (length(data_mats) == 1L) {
      data_mats[[1L]]
    } else {
      SeuratObject::StitchMatrix(
        x = data_mats[[1L]],
        y = data_mats[-1L],
        rowmap = methods::slot(assay_obj, "features")[
          features,
          src_layers,
          drop = FALSE
        ],
        colmap = methods::slot(assay_obj, "cells")[
          ,
          src_layers,
          drop = FALSE
        ]
      )
    }
    data_mat <- data_mat[features, colnames(object), drop = FALSE]
    idx <- seq_along(features) - 1L
  }

  sct_latent_df <- NULL
  regress_vars <- vars.to.regress
  if (length(regress_vars) > 0L) {
    meta <- methods::slot(object, "meta.data")
    found_meta <- intersect(regress_vars, colnames(meta))
    cell_order <- colnames(object)
    sct_latent_df <- if (length(found_meta)) {
      meta[match(cell_order, rownames(meta)), found_meta, drop = FALSE]
    } else {
      data.frame(row.names = cell_order)
    }
    # requested covariates that are assay features come from expression rows
    feature_vars <- setdiff(vars.to.regress, colnames(sct_latent_df))
    feature_vars <- intersect(feature_vars, rownames(assay_obj))
    if (length(feature_vars) > 0L) {
      expr_rows <- as.matrix(data_mat[
        match(feature_vars, rownames(data_mat)),
        cell_order,
        drop = FALSE
      ])
      df_feature <- as.data.frame(t(expr_rows))
      sct_latent_df <- cbind(sct_latent_df, df_feature)
    }
    if (!is.null(latent.data)) {
      extra <- as.data.frame(latent.data)
      missing_cells <- setdiff(cell_order, rownames(extra))
      if (length(missing_cells) > 0L) {
        log_message(
          "{.val {missing_cells}} are missing rows in {.val latent.data}.",
          message_type = "error"
        )
      }
      extra <- extra[cell_order, , drop = FALSE]
      extra <- extra[, !colnames(extra) %in% colnames(sct_latent_df), drop = FALSE]
      if (ncol(extra) > 0L) {
        sct_latent_df <- cbind(sct_latent_df, extra)
      }
    }
    missing_vars <- setdiff(regress_vars, colnames(sct_latent_df))
    if (length(missing_vars) == length(regress_vars) &&
      length(regress_vars) > 0L
    ) {
      log_message(
        "None of the requested variables to regress are present in the object.",
        message_type = "error"
      )
    }
    if (length(missing_vars) > 0L) {
      log_message(
        "Requested variables to regress not in object: {missing_vars}.",
        message_type = "warning",
        verbose = verbose
      )
    }
    if (any(!stats::complete.cases(sct_latent_df))) {
      log_message(
        "Regression variables contain missing values; ScaleData requires complete cases.",
        message_type = "error"
      )
    }
  }

  fast_path <-
    is.null(sct_latent_df) &&
      is.null(split.by) &&
      isFALSE(use.umi) &&
      isTRUE(do.scale) &&
      isTRUE(do.center)
  if (fast_path) {
    result <- scale_sparse_full(data_mat, idx, scale.max)
    dimnames(result) <- list(features, colnames(object))
  } else {
    sub <- data_mat[idx + 1L, , drop = FALSE]
    split_levels <- if (!is.null(split.by)) {
      meta <- methods::slot(object, "meta.data")
      factor(meta[match(colnames(object), rownames(meta)), split.by[1]])
    } else {
      NULL
    }
    general <- sct_scale_general(
      sub,
      sct_latent_df,
      split_levels,
      do.scale,
      do.center,
      scale.max,
      use.umi = want_counts
    )
    dimnames(general) <- list(features, colnames(object))
    result <- general
  }

  if (inherits(assay_obj, "Assay") && !inherits(assay_obj, "StdAssay")) {
    methods::slot(assay_obj, "scale.data") <- result
    methods::slot(object, "assays")[[assay_name]] <- assay_obj
    assay <- assay_name
    return(SeuratObject::LogSeuratCommand(object))
  }
  assay_obj@layers[["scale.data"]] <- result
  if (!"scale.data" %in% colnames(assay_obj@cells)) {
    cm <- assay_obj@cells
    methods::slot(assay_obj@cells, ".Data") <- cbind(
      cm,
      matrix(
        TRUE,
        nrow = nrow(cm),
        ncol = 1,
        dimnames = list(rownames(cm), "scale.data")
      )
    )
    fm <- assay_obj@features
    methods::slot(assay_obj@features, ".Data") <- cbind(
      fm,
      matrix(
        rownames(fm) %in% features,
        nrow = nrow(fm),
        ncol = 1,
        dimnames = list(rownames(fm), "scale.data")
      )
    )
  }
  methods::slot(object, "assays")[[assay_name]] <- assay_obj
  assay <- assay_name
  SeuratObject::LogSeuratCommand(object)
}

#' Scale expression data
#'
#' @param object Object containing expression data.
#' @param ... Passed to methods.
#'
#' @return The input object with scaled data.
#' @export
ScaleData <- function(object, ...) {
  UseMethod("ScaleData")
}

#' @export
ScaleData.default <- function(
  object,
  features = NULL,
  vars.to.regress = NULL,
  latent.data = NULL,
  split.by = NULL,
  model.use = "linear",
  use.umi = FALSE,
  do.scale = TRUE,
  do.center = TRUE,
  scale.max = 10,
  block.size = 1000,
  min.cells.to.block = 3000,
  verbose = TRUE,
  ...
) {
  fallback <- function() {
    utils::getFromNamespace("ScaleData.default", "Seurat")(
      object = object,
      features = features,
      vars.to.regress = vars.to.regress,
      latent.data = latent.data,
      split.by = split.by,
      model.use = model.use,
      use.umi = use.umi,
      do.scale = do.scale,
      do.center = do.center,
      scale.max = scale.max,
      block.size = block.size,
      min.cells.to.block = min.cells.to.block,
      verbose = verbose,
      ...
    )
  }
  if (
    !inherits(object, "dgCMatrix") ||
      !identical(model.use, "linear") ||
      !all(vapply(
        list(use.umi, do.scale, do.center),
        function(x) is.logical(x) && length(x) == 1L && !is.na(x),
        logical(1)
      )) ||
      length(list(...)) > 0L
  ) {
    return(fallback())
  }
  features <- features %||% rownames(object)
  features <- intersect(features, rownames(object))
  data <- object[features, , drop = FALSE]
  latent <- NULL
  if (!is.null(vars.to.regress)) {
    latent <- if (is.null(latent.data)) {
      data.frame(row.names = colnames(data))
    } else {
      as.data.frame(latent.data)[colnames(data), , drop = FALSE]
    }
    feature_vars <- intersect(vars.to.regress, rownames(data))
    if (length(feature_vars)) {
      latent <- cbind(latent, as.data.frame(t(data[feature_vars, , drop = FALSE])))
    }
    present <- intersect(vars.to.regress, colnames(latent))
    if (!length(present)) {
      return(fallback())
    }
    latent <- latent[, present, drop = FALSE]
    if (any(!stats::complete.cases(latent))) {
      return(fallback())
    }
  }
  split_levels <- if (is.null(split.by)) NULL else factor(split.by)
  if (!is.null(split_levels) && length(split_levels) != ncol(data)) {
    return(fallback())
  }
  if (
    is.null(latent) && is.null(split_levels) &&
      isTRUE(do.scale) && isTRUE(do.center)
  ) {
    out <- scale_sparse_full(data, seq_len(nrow(data)) - 1L, scale.max)
    dimnames(out) <- dimnames(data)
    return(out)
  }
  out <- sct_scale_general(
    data,
    latent,
    split_levels,
    do.scale,
    do.center,
    scale.max,
    use.umi = use.umi
  )
  out[is.na(out)] <- 0
  dimnames(out) <- dimnames(data)
  out
}
