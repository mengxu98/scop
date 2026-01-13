#' @title Run dimensionality reduction
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams GroupHeatmap
#' @param prefix The prefix used to name the result.
#' @param linear_reduction Method of linear dimensionality reduction.
#' Options are `"pca"`, `"ica"`, `"nmf"`, `"mds"`, `"glmpca"`.
#' @param linear_reduction_dims Total number of dimensions to compute and store for `linear_reduction`.
#' @param linear_reduction_params Other parameters passed to the `linear_reduction` method.
#' @param force_linear_reduction Whether force to do linear dimensionality reduction.
#' @param nonlinear_reduction Method of nonlinear dimensionality reduction.
#' Options are `"umap"`, `"umap-naive"`, `"tsne"`, `"dm"`, `"phate"`, `"pacmap"`, `"trimap"`, `"largevis"`.
#' @param nonlinear_reduction_dims Total number of dimensions to compute and store for `nonlinear_reduction`.
#' @param reduction_use Which dimensional reduction to use as input for `nonlinear_reduction`.
#' @param reduction_dims Which dimensions to use as input for `nonlinear_reduction`, used only if `features` is `NULL`.
#' @param neighbor_use Name of neighbor to use for the `nonlinear_reduction`.
#' @param graph_use Name of graph to use for the `nonlinear_reduction`.
#' @param nonlinear_reduction_params  Other parameters passed to the `nonlinear_reduction` method.
#' @param force_nonlinear_reduction Whether force to do nonlinear dimensionality reduction.
#'
#' @seealso
#' [DefaultReduction]
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

  linear_reductions <- c(
    "pca", "svd", "ica", "nmf", "mds", "glmpca"
  )
  reduction_exist <- SeuratObject::Reductions(srt)
  if (!is.null(linear_reduction)) {
    if (any(!linear_reduction %in% c(linear_reductions, reduction_exist)) || length(linear_reduction) > 1) {
      log_message(
        "{.arg linear_reduction} must be one of {.val {linear_reductions}}",
        message_type = "error"
      )
    }
  }
  nonlinear_reductions <- c(
    "umap", "umap-naive", "tsne", "dm", "phate", "pacmap", "trimap", "largevis", "fr"
  )
  if (!is.null(nonlinear_reduction)) {
    if (any(!nonlinear_reduction %in% c(nonlinear_reductions, reduction_exist)) || length(nonlinear_reduction) > 1) {
      log_message(
        "{.arg nonlinear_reduction} must be one of {.val {nonlinear_reductions}}",
        message_type = "error"
      )
    }

    if (is.null(features) && is.null(reduction_use) && is.null(neighbor_use) && is.null(graph_use)) {
      log_message(
        "'features', 'reduction_use', 'neighbor_use', or 'graph_use' must be provided when running non-linear dimensionality reduction",
        message_type = "error"
      )
    }

    if (nonlinear_reduction %in% c("fr")) {
      if (!is.null(graph_use)) {
        log_message(
          "Non-linear dimensionality reduction ({.pkg {nonlinear_reduction}}) using ({.pkg {graph_use}}) as input"
        )
      } else if (!is.null(neighbor_use)) {
        log_message(
          "Non-linear dimensionality reduction ({.pkg {nonlinear_reduction}}) using ({.pkg {neighbor_use}}) as input"
        )
      } else if (!is.null(features)) {
        log_message(
          "Non-linear dimensionality reduction ({.pkg {nonlinear_reduction}}) using ({.val {length(features)}} features) as input"
        )
      } else if (!is.null(reduction_use)) {
        log_message(
          "Non-linear dimensionality reduction ({.pkg {nonlinear_reduction}}) using ({.pkg {reduction_use}}) dims ({.val {min(reduction_dims)}}-{.val {max(reduction_dims)}}) as input"
        )
      }
    } else {
      if (!is.null(features)) {
        log_message(
          "Non-linear dimensionality reduction ({.pkg {nonlinear_reduction}}) using ({.val {length(features)}} features) as input"
        )
      } else if (!is.null(reduction_use)) {
        log_message(
          "Non-linear dimensionality reduction ({.pkg {nonlinear_reduction}}) using ({.pkg {reduction_use}}) dims ({.val {min(reduction_dims)}}-{.val {max(reduction_dims)}}) as input"
        )
      } else if (!is.null(neighbor_use)) {
        log_message(
          "Non-linear dimensionality reduction ({.pkg {nonlinear_reduction}}) using ({.pkg {neighbor_use}}) as input"
        )
      } else if (!is.null(graph_use)) {
        log_message(
          "Non-linear dimensionality reduction ({.pkg {nonlinear_reduction}}) using ({.pkg {graph_use}}) as input"
        )
      }
    }
  }

  if (!is.null(linear_reduction)) {
    if (isFALSE(force_linear_reduction)) {
      if (linear_reduction %in% reduction_exist) {
        if (srt[[linear_reduction]]@assay.used == assay) {
          log_message(
            "{.arg linear_reduction} {.pkg {linear_reduction}} is already existed. Skip calculation"
          )
          reduc <- srt[[linear_reduction]]
          SeuratObject::Key(reduc) <- paste0(prefix, linear_reduction, "_")
          srt[[paste0(prefix, linear_reduction)]] <- reduc
          srt@misc[["Default_reduction"]] <- paste0(prefix, linear_reduction)
          return(srt)
        } else {
          log_message(
            "assay.used is {.val {srt[[linear_reduction]]@assay.used}}, which is not the same as the {.val {assay}} specified. Recalculate the linear reduction (pca)"
          )
          linear_reduction <- "pca"
        }
      }
    }

    if (is.null(features) || length(features) == 0) {
      log_message("No features provided. Use variable features.")
      if (length(SeuratObject::DefaultAssay(srt)) == 0) {
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
    if (fun_use == "RunPCA") {
      params[["slot"]] <- params[["layer"]]
      params[["verbose"]] <- FALSE
      params <- params[!names(params) %in% "layer"]
    }
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
      pca_out <- srt[[paste0(prefix, linear_reduction)]]
      center <- rowMeans(
        GetAssayData5(
          object = srt,
          layer = "scale.data",
          assay = assay
        )[features, , drop = FALSE]
      )
      model <- list(
        sdev = pca_out@stdev,
        rotation = pca_out@feature.loadings,
        center = center,
        scale = FALSE,
        x = pca_out@cell.embeddings
      )
      class(model) <- "prcomp"
      srt@reductions[[paste0(prefix, linear_reduction)]]@misc[["model"]] <- model
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
          log_message(
            "Can not estimate intrinsic dimensions with {.pkg maxLikGlobalDimEst}",
            message_type = "warning"
          )
          NA
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
    srt@reductions[[paste0(prefix, linear_reduction)]]@misc[["dims_estimate"]] <- dims_estimate
    srt@misc[["Default_reduction"]] <- paste0(prefix, linear_reduction)
  } else if (!is.null(nonlinear_reduction)) {
    if (isFALSE(force_nonlinear_reduction)) {
      if (nonlinear_reduction %in% reduction_exist) {
        if (srt[[nonlinear_reduction]]@assay.used == assay) {
          log_message(
            "{.arg nonlinear_reduction} {.pkg {nonlinear_reduction}} is already existed. Skip calculation"
          )
          reduc <- srt[[nonlinear_reduction]]
          SeuratObject::Key(reduc) <- paste0(prefix, nonlinear_reduction, "_")
          srt[[paste0(prefix, nonlinear_reduction)]] <- reduc
          srt@misc[["Default_reduction"]] <- paste0(prefix, nonlinear_reduction)
          return(srt)
        } else {
          log_message(
            "assay.used is {.val {srt[[nonlinear_reduction]]@assay.used}}, which is not the same as the {.val {assay}} specified. Recalculate the nonlinear reduction (umap)"
          )
          nonlinear_reduction <- "umap"
        }
      }
    }

    fun_use <- switch(
      EXPR = nonlinear_reduction,
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
    components_nm <- switch(
      EXPR = nonlinear_reduction,
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
    other_params <- switch(
      EXPR = nonlinear_reduction,
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
    nonlinear_reduction_sim <- toupper(
      gsub(
        pattern = "-.*",
        replacement = "",
        x = nonlinear_reduction
      )
    )
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
      verbose = FALSE,
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
