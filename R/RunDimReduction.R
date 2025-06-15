#' Run dimensionality reduction
#'
#' @param srt A Seurat object.
#' @param prefix The prefix used to name the result.
#' @param features Use features expression data to run linear or nonlinear dimensionality reduction.
#' @param assay Specific assay to get data from.
#' @param layer Specific layer to get data from.
#' @param linear_reduction Method of linear dimensionality reduction. Options are "pca", "ica", "nmf", "mds", "glmpca".
#' @param linear_reduction_dims Total number of dimensions to compute and store for \code{linear_reduction}.
#' @param linear_reduction_params Other parameters passed to the \code{linear_reduction} method.
#' @param force_linear_reduction Whether force to do linear dimensionality reduction.
#' @param nonlinear_reduction Method of nonlinear dimensionality reduction. Options are "umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis"
#' @param nonlinear_reduction_dims Total number of dimensions to compute and store for \code{nonlinear_reduction}.
#' @param reduction_use Which dimensional reduction to use as input for \code{nonlinear_reduction}.
#' @param reduction_dims Which dimensions to use as input for \code{nonlinear_reduction}, used only if \code{features} is \code{NULL}.
#' @param neighbor_use Name of neighbor to use for the \code{nonlinear_reduction}.
#' @param graph_use Name of graph to use for the \code{nonlinear_reduction}.
#' @param nonlinear_reduction_params  Other parameters passed to the \code{nonlinear_reduction} method.
#' @param force_nonlinear_reduction Whether force to do nonlinear dimensionality reduction.
#' @param verbose Show messages.
#' @param seed Set a seed.
#'
#' @export
RunDimReduction <- function(
    srt,
    prefix = "",
    features = NULL,
    assay = NULL,
    layer = "data",
    linear_reduction = NULL,
    linear_reduction_dims = 50,
    linear_reduction_params = list(),
    force_linear_reduction = FALSE,
    nonlinear_reduction = NULL,
    nonlinear_reduction_dims = 2,
    reduction_use = NULL,
    reduction_dims = NULL,
    graph_use = NULL,
    neighbor_use = NULL,
    nonlinear_reduction_params = list(),
    force_nonlinear_reduction = TRUE,
    verbose = TRUE,
    seed = 11) {
  set.seed(seed)
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (inherits(srt[[assay]], "ChromatinAssay")) {
    type <- "Chromatin"
  } else {
    type <- "RNA"
  }
  linear_reduction_dims <- min(
    linear_reduction_dims,
    nrow(srt[[assay]]) - 1,
    ncol(srt[[assay]]) - 1,
    na.rm = TRUE
  )
  nonlinear_reduction_dims <- min(
    nonlinear_reduction_dims,
    nrow(srt[[assay]]) - 1,
    ncol(srt[[assay]]) - 1,
    na.rm = TRUE
  )
  if (!is.null(linear_reduction)) {
    if (
      any(
        !linear_reduction %in%
          c("pca", "svd", "ica", "nmf", "mds", "glmpca", SeuratObject::Reductions(srt))
      ) ||
        length(linear_reduction) > 1
    ) {
      stop(
        "'linear_reduction' must be one of 'pca','svd', 'ica', 'nmf', 'mds', 'glmpca'."
      )
    }
  }
  if (!is.null(nonlinear_reduction)) {
    if (
      any(
        !nonlinear_reduction %in%
          c(
            "umap",
            "umap-naive",
            "tsne",
            "dm",
            "phate",
            "pacmap",
            "trimap",
            "largevis",
            "fr",
            SeuratObject::Reductions(srt)
          )
      ) ||
        length(nonlinear_reduction) > 1
    ) {
      stop(
        "'nonlinear_reduction' must be one of 'umap', 'tsne', 'dm', 'phate', 'pacmap', 'trimap', 'largevis', 'fr'."
      )
    }
    if (
      is.null(features) &&
        is.null(reduction_use) &&
        is.null(neighbor_use) &&
        is.null(graph_use)
    ) {
      stop(
        "'features', 'reduction_use', 'neighbor_use', or 'graph_use' must be provided when running non-linear dimensionality reduction."
      )
    }
    if (nonlinear_reduction %in% c("fr")) {
      if (!is.null(graph_use)) {
        message(
          "Non-linear dimensionality reduction(",
          nonlinear_reduction,
          ") using Graphs(",
          graph_use,
          ") as input"
        )
      } else if (!is.null(neighbor_use)) {
        message(
          "Non-linear dimensionality reduction(",
          nonlinear_reduction,
          ") using Neighbors(",
          neighbor_use,
          ") as input"
        )
      } else if (!is.null(features)) {
        message(
          "Non-linear dimensionality reduction(",
          nonlinear_reduction,
          ") using Features(length:",
          length(features),
          ") as input"
        )
      } else if (!is.null(reduction_use)) {
        message(
          "Non-linear dimensionality reduction(",
          nonlinear_reduction,
          ") using Reduction(",
          reduction_use,
          ", dims:",
          min(reduction_dims),
          "-",
          max(reduction_dims),
          ") as input"
        )
      }
    } else {
      if (!is.null(features)) {
        message(
          "Non-linear dimensionality reduction(",
          nonlinear_reduction,
          ") using Features(length:",
          length(features),
          ") as input"
        )
      } else if (!is.null(reduction_use)) {
        message(
          "Non-linear dimensionality reduction(",
          nonlinear_reduction,
          ") using Reduction(",
          reduction_use,
          ", dims:",
          min(reduction_dims),
          "-",
          max(reduction_dims),
          ") as input"
        )
      } else if (!is.null(neighbor_use)) {
        message(
          "Non-linear dimensionality reduction(",
          nonlinear_reduction,
          ") using Neighbors(",
          neighbor_use,
          ") as input"
        )
      } else if (!is.null(graph_use)) {
        message(
          "Non-linear dimensionality reduction(",
          nonlinear_reduction,
          ") using Graphs(",
          graph_use,
          ") as input"
        )
      }
    }
  }
  if (!is.null(linear_reduction)) {
    if (!isTRUE(force_linear_reduction)) {
      if (linear_reduction %in% SeuratObject::Reductions(srt)) {
        if (srt[[linear_reduction]]@assay.used == assay) {
          message(
            "linear_reduction(",
            linear_reduction,
            ") is already existed. Skip calculation."
          )
          reduc <- srt[[linear_reduction]]
          SeuratObject::Key(reduc) <- paste0(prefix, linear_reduction, "_")
          srt[[paste0(prefix, linear_reduction)]] <- reduc
          srt@misc[["Default_reduction"]] <- paste0(prefix, linear_reduction)
          return(srt)
        } else {
          message(
            "assay.used is ",
            srt[[linear_reduction]]@assay.used,
            ", which is not the same as the ",
            assay,
            " specified. Recalculate the linear reduction(pca)"
          )
          linear_reduction <- "pca"
        }
      }
    }
    if (is.null(features) || length(features) == 0) {
      message("No features provided. Use variable features.")
      if (length(SeuratObject::DefaultAssay(srt, assay = assay)) == 0) {
        srt <- Seurat::FindVariableFeatures(srt, assay = assay, verbose = FALSE)
      }
      features <- SeuratObject::VariableFeatures(srt, assay = assay)
    }
    fun_use <- switch(linear_reduction,
      "pca" = "RunPCA",
      "svd" = "RunSVD",
      "ica" = "RunICA",
      "nmf" = "RunNMF",
      "mds" = "RunMDS",
      "glmpca" = "RunGLMPCA"
    )
    key_use <- switch(linear_reduction,
      "pca" = "PC_",
      "svd" = "LSI_",
      "ica" = "IC_",
      "nmf" = "BE_",
      "mds" = "MDS_",
      "glmpca" = "GLMPC_"
    )
    components_nm <- switch(linear_reduction,
      "pca" = "npcs",
      "svd" = "n",
      "ica" = "nics",
      "nmf" = "nbes",
      "mds" = "nmds",
      "glmpca" = "L"
    )
    params <- list(
      object = srt,
      assay = assay,
      layer = layer,
      features = features,
      components_nm = linear_reduction_dims,
      reduction.name = paste0(prefix, linear_reduction),
      reduction.key = paste0(prefix, key_use),
      verbose = verbose,
      seed.use = seed
    )
    if (fun_use %in% c("RunSVD", "RunICA")) {
      params <- params[!names(params) %in% "layer"]
    }
    if (fun_use == "RunGLMPCA") {
      params[["layer"]] <- "counts"
    }
    names(params)[names(params) == "components_nm"] <- components_nm
    for (nm in names(linear_reduction_params)) {
      params[[nm]] <- linear_reduction_params[[nm]]
    }
    srt <- invoke_fun(.fn = fun_use, .args = params)

    if (is.null(rownames(srt[[paste0(prefix, linear_reduction)]]))) {
      rownames(
        srt[[paste0(prefix, linear_reduction)]]@cell.embeddings
      ) <- colnames(srt)
    }
    if (linear_reduction == "pca") {
      pca.out <- srt[[paste0(prefix, linear_reduction)]]
      center <- rowMeans(
        GetAssayData5(
          object = srt,
          layer = "scale.data",
          assay = assay
        )[features, , drop = FALSE]
      )
      model <- list(
        sdev = pca.out@stdev,
        rotation = pca.out@feature.loadings,
        center = center,
        scale = FALSE,
        x = pca.out@cell.embeddings
      )
      class(model) <- "prcomp"
      srt@reductions[[paste0(prefix, linear_reduction)]]@misc[[
        "model"
      ]] <- model
    }
    if (linear_reduction %in% c("glmpca", "nmf")) {
      dims_estimate <- 1:linear_reduction_dims
    } else {
      dim_est <- tryCatch(
        expr = {
          min(
            intrinsicDimension::maxLikGlobalDimEst(
              data = Embeddings(
                srt,
                reduction = paste0(prefix, linear_reduction)
              ),
              k = 20
            )[["dim.est"]],
            ncol(Embeddings(srt, reduction = paste0(prefix, linear_reduction)))
          )
        },
        error = function(e) {
          message(
            "Can not estimate intrinsic dimensions with maxLikGlobalDimEst."
          )
          return(NA)
        }
      )
      if (!is.na(dim_est)) {
        dims_estimate <- seq_len(max(
          min(
            ncol(Embeddings(srt, reduction = paste0(prefix, linear_reduction))),
            10
          ),
          ceiling(dim_est)
        ))
      } else {
        dims_estimate <- seq_len(min(
          ncol(Embeddings(srt, reduction = paste0(prefix, linear_reduction))),
          30
        ))
      }
    }
    srt@reductions[[paste0(prefix, linear_reduction)]]@misc[[
      "dims_estimate"
    ]] <- dims_estimate
    srt@misc[["Default_reduction"]] <- paste0(prefix, linear_reduction)
  } else if (!is.null(nonlinear_reduction)) {
    if (!isTRUE(force_nonlinear_reduction)) {
      if (nonlinear_reduction %in% SeuratObject::Reductions(srt)) {
        if (srt[[nonlinear_reduction]]@assay.used == assay) {
          message(
            "nonlinear_reduction(",
            nonlinear_reduction,
            ") is already existed. Skip calculation."
          )
          reduc <- srt[[nonlinear_reduction]]
          SeuratObject::Key(reduc) <- paste0(prefix, nonlinear_reduction, "_")
          srt[[paste0(prefix, nonlinear_reduction)]] <- reduc
          srt@misc[["Default_reduction"]] <- paste0(prefix, nonlinear_reduction)
          return(srt)
        } else {
          message(
            "assay.used is ",
            srt[[nonlinear_reduction]]@assay.used,
            ", which is not the same as the ",
            assay,
            " specified. Recalculate the nonlinear reduction(umap)"
          )
          nonlinear_reduction <- "umap"
        }
      }
    }
    # if (!is.null(neighbor_use) && !nonlinear_reduction %in% c("umap", "umap-naive", "fr")) {
    #   stop("'neighbor_use' only support 'umap', 'umap-naive' or 'fr' method")
    # }
    # if (!is.null(graph_use) && !nonlinear_reduction %in% c("umap", "umap-naive", "fr")) {
    #   stop("'graph_use' only support 'umap', 'umap-naive' or 'fr' method")
    # }
    fun_use <- switch(nonlinear_reduction,
      "umap" = "RunUMAP2",
      "umap-naive" = "RunUMAP2",
      "tsne" = "RunTSNE",
      "dm" = "RunDM",
      "phate" = "RunPHATE",
      "pacmap" = "RunPaCMAP",
      "trimap" = "RunTriMap",
      "largevis" = "RunLargeVis",
      "fr" = "RunFR"
    )
    components_nm <- switch(nonlinear_reduction,
      "umap" = "n.components",
      "umap-naive" = "n.components",
      "tsne" = "dim.embed",
      "dm" = "ndcs",
      "phate" = "n_components",
      "pacmap" = "n_components",
      "trimap" = "n_components",
      "largevis" = "n_components",
      "fr" = "ndim"
    )
    other_params <- switch(nonlinear_reduction,
      "umap" = list(umap.method = "uwot", return.model = TRUE),
      "umap-naive" = list(umap.method = "naive", return.model = TRUE),
      "tsne" = list(
        tsne.method = "Rtsne",
        num_threads = 0,
        check_duplicates = FALSE
      ),
      "dm" = list(),
      "phate" = list(),
      "pacmap" = list(),
      "trimap" = list(),
      "largevis" = list(),
      "fr" = list()
    )
    nonlinear_reduction_sim <- toupper(gsub(
      pattern = "-.*",
      replacement = "",
      x = nonlinear_reduction
    ))
    params <- list(
      object = srt,
      assay = assay,
      layer = layer,
      components_nm = nonlinear_reduction_dims,
      features = features,
      reduction = reduction_use,
      dims = reduction_dims,
      reduction.name = paste0(
        prefix,
        nonlinear_reduction_sim,
        nonlinear_reduction_dims,
        "D"
      ),
      reduction.key = paste0(
        prefix,
        nonlinear_reduction_sim,
        nonlinear_reduction_dims,
        "D_"
      ),
      verbose = verbose,
      seed.use = seed
    )
    if (!is.null(neighbor_use)) {
      params[["neighbor"]] <- neighbor_use
    }
    if (!is.null(graph_use)) {
      params[["graph"]] <- graph_use
    }
    names(params)[names(params) == "components_nm"] <- components_nm
    for (nm in names(other_params)) {
      params[[nm]] <- other_params[[nm]]
    }
    for (nm in names(nonlinear_reduction_params)) {
      params[[nm]] <- nonlinear_reduction_params[[nm]]
    }
    srt <- invoke_fun(.fn = fun_use, .args = params)

    srt@reductions[[paste0(
      prefix,
      nonlinear_reduction_sim,
      nonlinear_reduction_dims,
      "D"
    )]]@misc[["reduction_dims"]] <- reduction_dims
    srt@reductions[[paste0(
      prefix,
      nonlinear_reduction_sim,
      nonlinear_reduction_dims,
      "D"
    )]]@misc[["reduction_use"]] <- reduction_use
    srt@misc[["Default_reduction"]] <- paste0(prefix, nonlinear_reduction_sim)
  }
  return(srt)
}
