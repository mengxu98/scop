#' Check and report the type of data
#'
#' This function checks the type of data and returns a string indicating the type of data. It checks for the presence of infinite values, negative values, and whether the values are floats or integers.
#'
#' @param srt An object of class 'Seurat'.
#' @param data The input data. If not provided, it will be extracted from the the 'srt' object.
#' @param layer The layer in the 'srt' object from which to extract the data. Default is "data".
#' @param assay The assay to extract the data from. If not provided, the default assay will be used.
#'
#' @return A string indicating the type of data. Possible values are: "raw_counts", "log_normalized_counts", "raw_normalized_counts", or "unknown".
#'
#' @export
check_data_type <- function(
    srt,
    data = NULL,
    layer = "data",
    assay = NULL) {
  if (is.null(data)) {
    assay <- assay %||% SeuratObject::DefaultAssay(srt)
    data <- Seurat::GetAssayData(
      srt,
      layer = layer,
      assay = assay
    )
  }
  isfinite <- all(is.finite(range(data, na.rm = TRUE)))
  if (inherits(data, "dgCMatrix")) {
    isfloat <- any(data@x %% 1 != 0, na.rm = TRUE)
  } else {
    isfloat <- any(
      data[, sample(seq_len(ncol(data)), min(ncol(data), 1000))] %% 1 != 0,
      na.rm = TRUE
    )
  }
  islog <- is.finite(expm1(x = max(data, na.rm = TRUE)))
  isnegative <- any(data < 0)

  if (!isTRUE(isfinite)) {
    warning("Infinite values detected!", immediate. = TRUE)
    return("unknown")
  } else if (isTRUE(isnegative)) {
    warning("Negative values detected!", immediate. = TRUE)
    return("unknown")
  } else {
    if (!isfloat) {
      return("raw_counts")
    } else if (isfloat && islog) {
      return("log_normalized_counts")
    } else if (isfloat && !islog) {
      if (isFALSE(isnegative)) {
        return("raw_normalized_counts")
      } else {
        return("unknown")
      }
    }
  }
}

#' Check and preprocess a list of seurat objects
#'
#' This function checks and preprocesses a list of seurat objects. It performs various checks on the input, including verification of input types, assay type consistency, feature name consistency, and batch column consistency. It also performs data normalization and variable feature finding based on the specified parameters. Finally, it prepares the data for integration analysis based on the highly variable features.
#'
#' @param srt_list A list of Seurat objects to be checked and preprocessed.
#' @param batch A character string specifying the batch variable name.
#' @param assay The name of the assay to be used for downstream analysis.
#' @param do_normalization A logical value indicating whether data normalization should be performed.
#' @param normalization_method The normalization method to be used. Possible values are "LogNormalize", "SCT", and "TFIDF". Default is "LogNormalize".
#' @param do_HVF_finding A logical value indicating whether highly variable feature (HVF) finding should be performed. Default is TRUE.
#' @param HVF_source The source of highly variable features. Possible values are "global" and "separate". Default is "separate".
#' @param HVF_method The method for selecting highly variable features. Default is "vst".
#' @param nHVF The number of highly variable features to select. Default is 2000.
#' @param HVF_min_intersection The feature needs to be present in batches for a minimum number of times in order to be considered as highly variable. The default value is 1.
#' @param HVF A vector of highly variable features. Default is NULL.
#' @param vars_to_regress A vector of variable names to include as additional regression variables. Default is NULL.
#' @param seed An integer specifying the random seed for reproducibility. Default is 11.
#'
#' @return A list containing the preprocessed seurat objects, the highly variable features, the assay name, and the type of assay (e.g., "RNA" or "Chromatin").
#'
#' @importFrom Signac RunTFIDF
#' @export
#'
check_srt_list <- function(
    srt_list,
    batch,
    assay = NULL,
    do_normalization = NULL,
    normalization_method = "LogNormalize",
    do_HVF_finding = TRUE,
    HVF_source = "separate",
    HVF_method = "vst",
    nHVF = 2000,
    HVF_min_intersection = 1,
    HVF = NULL,
    vars_to_regress = NULL,
    seed = 11) {
  cat(paste0("[", Sys.time(), "]", " Checking srt_list... ...\n"))
  set.seed(seed)

  if (
    !inherits(srt_list, "list") ||
      any(sapply(srt_list, function(x) !inherits(x, "Seurat")))
  ) {
    stop("'srt_list' is not a list of Seurat object.")
  }
  if (!normalization_method %in% c("LogNormalize", "SCT", "TFIDF")) {
    stop(
      "'normalization_method' must be one of: 'LogNormalize', 'SCT', 'TFIDF'"
    )
  }
  if (normalization_method %in% c("SCT")) {
    check_r("glmGamPoi")
  }
  if (!HVF_source %in% c("global", "separate")) {
    stop("'HVF_source' must be one of: 'global', 'separate'")
  }
  if (any(sapply(srt_list, ncol) < 2)) {
    stop(paste0(
      "Seurat objects in srt_list contain less than 2 cells. srt_list index: ",
      paste0(which(sapply(srt_list, ncol) < 2), collapse = ",")
    ))
  }

  if (is.null(assay)) {
    default_assay <- unique(sapply(srt_list, SeuratObject::DefaultAssay))
    if (length(default_assay) != 1) {
      stop(
        "The default assay name of the Seurat object in the srt_list is inconsistent."
      )
    } else {
      assay <- default_assay
    }
  }

  assay_type <- unique(
    sapply(
      srt_list, function(srt) {
        class(
          Seurat::GetAssay(
            srt,
            assay = assay
          )
        )
      }
    )
  )
  if (length(assay_type) != 1) {
    stop("The assay type of the Seurat object in the srt_list is inconsistent.")
  } else {
    if (assay_type == "Assay" || assay_type == "Assay5") {
      type <- "RNA"
    } else if (assay_type == "ChromatinAssay") {
      type <- "Chromatin"
    } else {
      type <- "Unknown"
    }
  }

  features_list <- lapply(
    srt_list, function(srt) {
      sort(
        rownames(
          Seurat::GetAssay(
            srt,
            assay = assay
          )
        )
      )
    }
  )
  if (length(unique(features_list)) != 1) {
    if (type == "Chromatin") {
      warning(
        "The peaks in assay ",
        assay,
        " is different between batches. Creating a common set of peaks and may take a long time..."
      )
      srt_merge <- Reduce(merge, srt_list)
      srt_list <- Seurat::SplitObject(
        object = srt_merge,
        split.by = batch
      )
    }
    cf <- Reduce(
      intersect,
      lapply(
        srt_list, function(srt) {
          rownames(
            Seurat::GetAssay(
              srt,
              assay = assay
            )
          )
        }
      )
    )
    warning(
      "'srt_list' have different feature names! Will subset the common features(",
      length(cf),
      ") for downstream analysis!",
      immediate. = TRUE
    )
    for (i in seq_along(srt_list)) {
      srt_list[[i]][[assay]] <- subset(srt_list[[i]][[assay]], features = cf)
    }
  }

  celllist <- unlist(lapply(srt_list, colnames))
  if (length(celllist) != length(unique(celllist))) {
    stop("'srt_list' have duplicated cell names!")
  }

  if (length(batch) != 1 && length(batch) != length(srt_list)) {
    stop(
      "'batch' must be a character to specify the batch column in the meta.data or a vector of the same length of the srt_list!"
    )
  }
  if (length(batch) == length(srt_list)) {
    srt_list_tmp <- list()
    for (bat in unique(batch)) {
      srt_list_tmp[[bat]] <- Reduce(merge, srt_list[batch == bat])
    }
    srt_list <- srt_list_tmp
  } else {
    if (
      !all(sapply(srt_list, function(x) {
        batch %in% colnames(x@meta.data)
      }))
    ) {
      stop(paste0(
        "batch column('",
        batch,
        "') was not found in one or more object of the srt_list!"
      ))
    }
    for (i in seq_along(srt_list)) {
      u <- unique(srt_list[[i]][[batch, drop = TRUE]])
      if (length(u) > 1) {
        x <- Seurat::SplitObject(srt_list[[i]], split.by = batch)
        srt_list[[i]] <- character(0)
        srt_list <- c(srt_list, x)
      }
    }
    srt_list <- srt_list[sapply(srt_list, length) > 0]
    srt_list_batch <- sapply(
      srt_list,
      function(x) unique(x[[batch, drop = TRUE]])
    )
    batch_to_merge <- names(which(table(srt_list_batch) > 1))
    if (length(batch_to_merge) > 0) {
      for (b in batch_to_merge) {
        index <- which(srt_list_batch == b)
        srt_list_tmp <- Reduce(merge, srt_list[index])
        for (i in index) {
          srt_list[[i]] <- character(0)
        }
        srt_list <- c(srt_list, srt_list_tmp)
      }
    }
    srt_list <- srt_list[sapply(srt_list, length) > 0]
  }

  for (i in seq_along(srt_list)) {
    if (!assay %in% SeuratObject::Assays(srt_list[[i]])) {
      stop(paste0("srt_list ", i, " does not contain '", assay, "' assay."))
    }
    SeuratObject::DefaultAssay(srt_list[[i]]) <- assay
    if (isTRUE(do_normalization)) {
      if (normalization_method == "LogNormalize") {
        cat(
          "Perform NormalizeData(LogNormalize) on the data ",
          i,
          "/",
          length(srt_list),
          " of the srt_list...\n",
          sep = ""
        )
        srt_list[[i]] <- NormalizeData(
          object = srt_list[[i]],
          assay = assay,
          normalization.method = "LogNormalize",
          verbose = FALSE
        )
      }
      if (normalization_method == "TFIDF") {
        cat(
          "Perform RunTFIDF on the data ",
          i,
          "/",
          length(srt_list),
          " of the srt_list...\n",
          sep = ""
        )
        srt_list[[i]] <- RunTFIDF(
          object = srt_list[[i]],
          assay = assay,
          verbose = FALSE
        )
      }
    } else if (is.null(do_normalization)) {
      status <- check_data_type(srt_list[[i]], layer = "data", assay = assay)
      if (status == "log_normalized_counts") {
        cat(
          "Data ",
          i,
          "/",
          length(srt_list),
          " of the srt_list has been log-normalized.\n",
          sep = ""
        )
      }
      if (status %in% c("raw_counts", "raw_normalized_counts")) {
        if (normalization_method == "LogNormalize") {
          cat(
            "Data ",
            i,
            "/",
            length(srt_list),
            " of the srt_list is ",
            status,
            ". Perform NormalizeData(LogNormalize) on the data ...\n",
            sep = ""
          )
          srt_list[[i]] <- Seurat::NormalizeData(
            object = srt_list[[i]],
            assay = assay,
            normalization.method = "LogNormalize",
            verbose = FALSE
          )
        }
        if (normalization_method == "TFIDF") {
          cat(
            "Data ",
            i,
            "/",
            length(srt_list),
            " of the srt_list is ",
            status,
            ". Perform RunTFIDF on the data ...\n",
            sep = ""
          )
          srt_list[[i]] <- RunTFIDF(
            object = srt_list[[i]],
            assay = assay,
            verbose = FALSE
          )
        }
      }
      if (status == "unknown") {
        warning(
          "Can not determine whether data ",
          i,
          " is log-normalized...",
          immediate. = TRUE
        )
      }
    }
    if (is.null(HVF)) {
      if (
        isTRUE(do_HVF_finding) ||
          is.null(do_HVF_finding) ||
          length(VariableFeatures(srt_list[[i]], assay = assay)) == 0
      ) {
        # if (type == "RNA") {
        cat(
          "Perform FindVariableFeatures on the data ",
          i,
          "/",
          length(srt_list),
          " of the srt_list...\n",
          sep = ""
        )
        srt_list[[i]] <- Seurat::FindVariableFeatures(
          srt_list[[i]],
          assay = assay,
          nfeatures = nHVF,
          selection.method = HVF_method,
          verbose = FALSE
        )
        # }
        # if (type == "Chromatin") {
        #   cat("Perform FindTopFeatures on the data ", i, "/", length(srt_list), " of the srt_list...\n", sep = "")
        #   srt_list[[i]] <- FindTopFeatures(srt_list[[i]], assay = assay, min.cutoff = HVF_min_cutoff, verbose = FALSE)
        # }
      }
    }

    if (normalization_method %in% c("SCT") && type == "RNA") {
      if (
        isTRUE(do_normalization) ||
          isTRUE(do_HVF_finding) ||
          !"SCT" %in% SeuratObject::Assays(srt_list[[i]])
      ) {
        cat("Perform SCTransform on the data", i, "of the srt_list...\n")
        srt_list[[i]] <- Seurat::SCTransform(
          object = srt_list[[i]],
          variable.features.n = nHVF,
          vars.to.regress = vars_to_regress,
          assay = assay,
          method = "glmGamPoi",
          new.assay.name = "SCT",
          verbose = FALSE
        )
      } else {
        SeuratObject::DefaultAssay(srt_list[[i]]) <- "SCT"
      }
      if (
        !"residual_variance" %in%
          colnames(
            GetFeaturesData(srt_list[[i]], assay = "SCT")
          )
      ) {
        if (length(srt_list[[i]]@assays$SCT@SCTModel.list) > 1) {
          index <- which(sapply(
            srt_list[[i]]@assays$SCT@SCTModel.list,
            function(x) nrow(x@cell.attributes) == ncol(srt_list[[i]])
          ))
        } else {
          index <- 1
        }
        model <- srt_list[[i]]@assays$SCT@SCTModel.list[[index]]
        feature.attr <- Seurat::SCTResults(
          object = model,
          slot = "feature.attributes"
        )
      } else {
        feature.attr <- GetFeaturesData(srt_list[[i]], assay = "SCT")
      }
      nfeatures <- min(nHVF, nrow(x = feature.attr))
      top.features <- rownames(x = feature.attr)[head(
        order(feature.attr$residual_variance, decreasing = TRUE),
        n = nfeatures
      )]
      VariableFeatures(
        srt_list[[i]],
        assay = SeuratObject::DefaultAssay(srt_list[[i]])
      ) <- top.features
      # srt_list[[i]]@assays$SCT@meta.features <- feature.attr
      srt_list[[i]] <- AddFeaturesData(
        srt_list[[i]],
        features = feature.attr,
        assay = "SCT"
      )
    }
  }

  if (is.null(HVF)) {
    if (HVF_source == "global") {
      cat("Use the global HVF from merged dataset...\n")
      srt_merge <- Reduce(merge, srt_list)
      # if (type == "RNA") {
      srt_merge <- Seurat::FindVariableFeatures(
        srt_merge,
        assay = SeuratObject::DefaultAssay(srt_merge),
        nfeatures = nHVF,
        selection.method = HVF_method,
        verbose = FALSE
      )
      # }
      # if (type == "Chromatin") {
      #   srt_merge <- FindTopFeatures(srt_merge, assay = DefaultAssay(srt_merge), min.cutoff = HVF_min_cutoff, verbose = FALSE)
      # }
      HVF <- SeuratObject::VariableFeatures(srt_merge)
    }
    if (HVF_source == "separate") {
      cat("Use the separate HVF from srt_list...\n")
      # if (type == "RNA") {
      HVF <- Seurat::SelectIntegrationFeatures(
        object.list = srt_list,
        nfeatures = nHVF,
        verbose = FALSE
      )
      HVF_sort <- sort(
        table(unlist(lapply(srt_list, SeuratObject::VariableFeatures))),
        decreasing = TRUE
      )
      HVF_filter <- HVF_sort[HVF_sort >= HVF_min_intersection]
      HVF <- intersect(HVF, names(HVF_filter))
      # }
      # if (type == "Chromatin") {
      #   nHVF <- min(sapply(srt_list, function(srt) length(VariableFeatures(srt))))
      #   HVF_sort <- sort(table(unlist(lapply(srt_list, VariableFeatures))), decreasing = TRUE)
      #   HVF_filter <- HVF_sort[HVF_sort >= HVF_min_intersection]
      #   HVF <- names(head(HVF_filter, nHVF))
      # }
      if (length(HVF) == 0) {
        stop("No HVF available.")
      }
    }
  } else {
    cf <- Reduce(
      intersect,
      lapply(srt_list, function(srt) {
        rownames(
          Seurat::GetAssayData(
            srt,
            layer = "counts",
            assay = SeuratObject::DefaultAssay(srt)
          )
        )
      })
    )
    HVF <- HVF[HVF %in% cf]
  }
  message("Number of available HVF: ", length(HVF))

  hvf_sum <- lapply(srt_list, function(srt) {
    colSums(
      Seurat::GetAssayData(
        srt,
        layer = "counts",
        assay = SeuratObject::DefaultAssay(srt)
      )[
        HVF, ,
        drop = FALSE
      ]
    )
  })
  cell_all <- unlist(unname(hvf_sum))
  cell_abnormal <- names(cell_all)[cell_all == 0]
  if (length(cell_abnormal) > 0) {
    warning(
      "Some cells do not express any of the highly variable features: ",
      paste(cell_abnormal, collapse = ","),
      immediate. = TRUE
    )
  }

  if (normalization_method == "SCT" && type == "RNA") {
    srt_list <- Seurat::PrepSCTIntegration(
      object.list = srt_list,
      anchor.features = HVF,
      assay = "SCT",
      verbose = FALSE
    )
  }
  cat(paste0("[", Sys.time(), "]", " Finished checking.\n"))

  return(
    list(
      srt_list = srt_list,
      HVF = HVF,
      assay = assay,
      type = type
    )
  )
}

#' Check and preprocess a merged seurat object
#'
#' This function checks and preprocesses a merged seurat object.
#'
#' @inheritParams check_srt_list
#' @inheritParams Integration_scop
#' @param srt_merge A merged Seurat object that includes the batch information.
#'
#' @seealso \link{check_srt_list}
#'
#' @export
check_srt_merge <- function(
    srt_merge,
    batch = NULL,
    assay = NULL,
    do_normalization = NULL,
    normalization_method = "LogNormalize",
    do_HVF_finding = TRUE,
    HVF_source = "separate",
    HVF_method = "vst",
    nHVF = 2000,
    HVF_min_intersection = 1,
    HVF = NULL,
    vars_to_regress = NULL,
    seed = 11) {
  if (!inherits(srt_merge, "Seurat")) {
    stop("'srt_merge' is not a Seurat object.")
  }
  if (length(batch) != 1) {
    stop(
      "'batch' must be provided to specify the batch column in the meta.data"
    )
  }
  if (!batch %in% colnames(srt_merge@meta.data)) {
    stop(paste0("No batch column('", batch, "') found in the meta.data"))
  }
  if (!is.factor(srt_merge[[batch, drop = TRUE]])) {
    srt_merge[[batch, drop = TRUE]] <- factor(
      srt_merge[[batch, drop = TRUE]],
      levels = unique(srt_merge[[batch, drop = TRUE]])
    )
  }
  assay <- assay %||% SeuratObject::DefaultAssay(srt_merge)
  srt_merge_raw <- srt_merge

  cat(paste0(
    "[",
    Sys.time(),
    "]",
    " Spliting srt_merge into srt_list by column ",
    batch,
    "... ...\n"
  ))
  srt_list <- Seurat::SplitObject(
    object = srt_merge_raw,
    split.by = batch
  )

  checked <- check_srt_list(
    srt_list = srt_list,
    batch = batch,
    assay = assay,
    do_normalization = do_normalization,
    do_HVF_finding = do_HVF_finding,
    normalization_method = normalization_method,
    HVF_source = HVF_source,
    HVF_method = HVF_method,
    nHVF = nHVF,
    HVF_min_intersection = HVF_min_intersection,
    HVF = HVF,
    vars_to_regress = vars_to_regress,
    seed = seed
  )
  srt_list <- checked[["srt_list"]]
  HVF <- checked[["HVF"]]
  assay <- checked[["assay"]]
  type <- checked[["type"]]
  srt_merge <- Reduce(merge, srt_list)

  srt_merge <- SrtAppend(
    srt_raw = srt_merge,
    srt_append = srt_merge_raw,
    pattern = "",
    slots = "reductions",
    overwrite = TRUE,
    verbose = FALSE
  )
  if (normalization_method == "SCT" && type == "RNA") {
    SeuratObject::DefaultAssay(srt_merge) <- "SCT"
  } else {
    SeuratObject::DefaultAssay(srt_merge) <- assay
  }
  SeuratObject::VariableFeatures(srt_merge) <- HVF

  return(
    list(
      srt_merge = srt_merge,
      srt_list = srt_list,
      HVF = HVF,
      assay = assay,
      type = type
    )
  )
}
