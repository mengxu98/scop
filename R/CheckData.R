#' @title Check and report the type of data in `Seurat` object
#'
#' @description
#' This function checks and returns a string indicating the type of data.
#' It checks for the presence of infinite values, negative values,
#' and whether the values are floats or integers.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param object A `Seurat` object or a matrix.
#'
#' @return
#' A string indicating the type of data.
#' Possible values are:
#' `"raw_counts"`, `"log_normalized_counts"`, `"raw_normalized_counts"`, or `"unknown"`.
#'
#' @export
CheckDataType <- function(object, ...) {
  UseMethod(generic = "CheckDataType", object = object)
}

#' @md
#' @param layer The layer in the `srt` object from which to extract the data.
#' Default is `"data"`.
#' @param assay The assay to extract the data from.
#' If not provided, the default assay will be used.
#'
#' @rdname CheckDataType
#' @method CheckDataType Seurat
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' CheckDataType(pancreas_sub)
CheckDataType.Seurat <- function(
    object,
    layer = "data",
    assay = NULL,
    verbose = TRUE,
    ...) {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  data <- GetAssayData5(
    object,
    layer = layer,
    assay = assay,
    ...
  )
  CheckDataType(data, verbose = verbose)
}

#' @rdname CheckDataType
#' @method CheckDataType default
#' @export
CheckDataType.default <- function(
    object,
    verbose = TRUE,
    ...) {
  isfinite <- all(is.finite(range(object, na.rm = TRUE)))
  if (inherits(object, "dgCMatrix")) {
    isfloat <- any(object@x %% 1 != 0, na.rm = TRUE)
  } else {
    isfloat <- any(
      object[, sample(seq_len(ncol(object)), min(ncol(object), 1000))] %% 1 != 0,
      na.rm = TRUE
    )
  }
  islog <- is.finite(expm1(x = max(object, na.rm = TRUE)))
  isnegative <- any(object < 0)

  if (isFALSE(isfinite)) {
    log_message(
      "Infinite values detected",
      message_type = "warning",
      verbose = verbose
    )
    return("unknown")
  } else if (isTRUE(isnegative)) {
    log_message(
      "Negative values detected",
      message_type = "warning",
      verbose = verbose
    )
    return("unknown")
  } else {
    if (!isfloat) {
      log_message(
        "Data type is raw counts",
        verbose = verbose
      )
      return("raw_counts")
    } else if (isfloat && islog) {
      log_message(
        "Data type is log-normalized",
        verbose = verbose
      )
      return("log_normalized_counts")
    } else if (isfloat && !islog) {
      if (isFALSE(isnegative)) {
        log_message(
          "Data type is normalized without log transformation",
          verbose = verbose
        )
        return("raw_normalized_counts")
      } else {
        log_message(
          "Can not determine whether data type is log-normalized",
          message_type = "warning",
          verbose = verbose
        )
        return("unknown")
      }
    }
  }
}

#' @title Check and preprocess a list of `Seurat` objects
#'
#' @description
#' This function checks and preprocesses a list of `Seurat` objects.
#' It performs various checks on the input, including verification of input types,
#' assay type consistency, feature name consistency, and batch column consistency.
#' It also performs data normalization and variable feature finding based on the specified parameters.
#' Finally, it prepares the data for integration analysis based on the highly variable features.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt_list A list of `Seurat` objects to be checked and preprocessed.
#' @param batch A character string specifying the batch variable name.
#' @param assay The name of the assay to be used for downstream analysis.
#' @param do_normalization Whether data normalization should be performed.
#' Default is `TRUE`.
#' @param normalization_method The normalization method to be used.
#' Possible values are `"LogNormalize"`, `"SCT"`, and `"TFIDF"`.
#' Default is `"LogNormalize"`.
#' @param do_HVF_finding Whether highly variable feature (HVF) finding should be performed.
#' Default is `TRUE`.
#' @param HVF_source The source of highly variable features.
#' Possible values are `"global"` and `"separate"`.
#' Default is `"separate"`.
#' @param HVF_method The method for selecting highly variable features.
#' Default is `"vst"`.
#' @param nHVF The number of highly variable features to select.
#' Default is `2000`.
#' @param HVF_min_intersection The feature needs to be present in batches for a minimum number of times in order to be considered as highly variable.
#' Default is `1`.
#' @param HVF A vector of highly variable features.
#' Default is `NULL`.
#' @param vars_to_regress A vector of variable names to include as additional regression variables.
#' Default is `NULL`.
#' @param seed An integer specifying the random seed for reproducibility.
#' Default is `11`.
#'
#' @return
#' A list containing the preprocessed `Seurat` objects,
#' the highly variable features, the assay name,
#' and the type of assay.
#'
#' @export
CheckDataList <- function(
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
    verbose = TRUE,
    seed = 11) {
  log_message(
    "Checking a list of {.cls Seurat} object...",
    verbose = verbose
  )
  set.seed(seed)

  if (!inherits(srt_list, "list") || any(sapply(srt_list, function(x) !inherits(x, "Seurat")))) {
    log_message(
      "{.arg srt_list} is not a list of {.cls Seurat} objects",
      message_type = "error"
    )
  }
  normalization_methods <- c("LogNormalize", "SCT", "TFIDF")
  if (!normalization_method %in% normalization_methods) {
    log_message(
      "{.arg normalization_method} must be one of: {.val {normalization_methods}}",
      message_type = "error"
    )
  }
  hvf_sources <- c("global", "separate")
  if (!HVF_source %in% hvf_sources) {
    log_message(
      "{.arg HVF_source} must be one of: {.val {hvf_sources}}",
      message_type = "error"
    )
  }
  which_less_2 <- which(sapply(srt_list, ncol) < 2)
  if (length(which_less_2) > 0) {
    log_message(
      "{.cls Seurat} objects in {.arg srt_list} contain less than 2 cells. {.arg srt_list} index: {.val {which_less_2}}",
      message_type = "error"
    )
  }

  if (is.null(assay)) {
    default_assay <- unique(sapply(srt_list, SeuratObject::DefaultAssay))
    if (length(default_assay) != 1) {
      log_message(
        "The default assay name of the Seurat object in the {.arg srt_list} is inconsistent",
        message_type = "error"
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
    log_message(
      "The assay type of the {.cls Seurat} object in the {.arg srt_list} is inconsistent",
      message_type = "error"
    )
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
      log_message(
        "The peaks in assay {.val {assay}} is different between batches. Creating a common set...",
        message_type = "warning",
        verbose = verbose
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
    log_message(
      "{.arg srt_list} have different feature names! Will subset the common features ({.val {length(cf)}}) for downstream analysis",
      message_type = "warning",
      verbose = verbose
    )
    for (i in seq_along(srt_list)) {
      srt_list[[i]][[assay]] <- subset(srt_list[[i]][[assay]], features = cf)
    }
  }

  celllist <- unlist(lapply(srt_list, colnames))
  if (length(celllist) != length(unique(celllist))) {
    log_message(
      "{.arg srt_list} have duplicated cell names",
      message_type = "error",
      verbose = verbose
    )
  }

  if (length(batch) != 1 && length(batch) != length(srt_list)) {
    log_message(
      "{.arg batch} must be a character to specify the batch column in the meta.data or a vector of the same length of the {.arg srt_list}",
      message_type = "error"
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
      log_message(
        "batch column ({.val {batch}}) was not found in one or more object of the {.arg srt_list}",
        message_type = "error"
      )
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
      log_message(
        "{.arg srt_list} {.val {i}} does not contain {.val {assay}} assay",
        message_type = "error"
      )
    }
    SeuratObject::DefaultAssay(srt_list[[i]]) <- assay
    status <- suppressWarnings(
      CheckDataType(
        srt_list[[i]],
        layer = "data",
        assay = assay,
        verbose = FALSE
      )
    )
    if (status == "log_normalized_counts") {
      log_message(
        "Data {.val {i}}/{.val {length(srt_list)}} of the {.arg srt_list} has been log-normalized",
        verbose = verbose
      )
    } else if (status %in% c("raw_counts", "raw_normalized_counts", "unknown")) {
      log_message(
        "Data {.val {i}}/{.val {length(srt_list)}} of the {.arg srt_list} is {.val {status}}",
        message_type = "warning",
        verbose = verbose
      )
      do_normalization <- TRUE
    }
    if (isTRUE(do_normalization)) {
      if (normalization_method == "LogNormalize") {
        log_message(
          "Perform {.fn NormalizeData} with {.arg normalization.method = 'LogNormalize'} on the data {.val {i}}/{.val {length(srt_list)}} of the {.arg srt_list}...",
          verbose = verbose
        )
        srt_list[[i]] <- NormalizeData(
          object = srt_list[[i]],
          assay = assay,
          normalization.method = "LogNormalize",
          verbose = FALSE
        )
        srt_list[[i]] <- ScaleData(
          object = srt_list[[i]],
          assay = assay,
          verbose = FALSE
        )
      }
      if (normalization_method == "TFIDF") {
        log_message(
          "Perform {.fn RunTFIDF} on the data {.val {i}}/{.val {length(srt_list)}} of the {.arg srt_list}...",
          verbose = verbose
        )
        srt_list[[i]] <- Signac::RunTFIDF(
          object = srt_list[[i]],
          assay = assay,
          verbose = FALSE
        )
      }
    }

    if (is.null(HVF)) {
      if (
        isTRUE(do_HVF_finding) ||
          is.null(do_HVF_finding) ||
          length(SeuratObject::VariableFeatures(srt_list[[i]], assay = assay)) == 0
      ) {
        if (type == "RNA") {
          log_message(
            "Perform {.fn Seurat::FindVariableFeatures} on the data {.val {i}}/{.val {length(srt_list)}} of the {.arg srt_list}...",
            verbose = verbose
          )
          srt_list[[i]] <- Seurat::FindVariableFeatures(
            srt_list[[i]],
            assay = assay,
            nfeatures = nHVF,
            selection.method = HVF_method,
            verbose = FALSE
          )
        }
        if (type == "Chromatin") {
          log_message(
            "Perform {.fn FindTopFeatures} on the data {.val {i}}/{.val {length(srt_list)}} of the {.arg srt_list}...",
            verbose = verbose
          )
          srt_list[[i]] <- Signac::FindTopFeatures(
            srt_list[[i]],
            assay = assay,
            min.cutoff = "q5",
            verbose = FALSE
          )
        }
      }
    }

    if (normalization_method == "SCT" && type == "RNA") {
      check_r("glmGamPoi", verbose = FALSE)
      if (
        isTRUE(do_normalization) ||
          isTRUE(do_HVF_finding) ||
          !"SCT" %in% SeuratObject::Assays(srt_list[[i]])
      ) {
        log_message(
          "Perform {.fn Seurat::SCTransform} on the data {.val {i}}/{.val {length(srt_list)}} of the {.arg srt_list}...",
          verbose = verbose
        )
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
      feature_attr_sct <- GetFeaturesData(srt_list[[i]], assay = "SCT")
      if (!"residual_variance" %in% colnames(feature_attr_sct)) {
        if (length(srt_list[[i]]@assays$SCT@SCTModel.list) > 1) {
          index <- which(
            sapply(
              srt_list[[i]]@assays$SCT@SCTModel.list,
              function(x) nrow(x@cell.attributes) == ncol(srt_list[[i]])
            )
          )
        } else {
          index <- 1
        }
        model <- srt_list[[i]]@assays$SCT@SCTModel.list[[index]]
        feature_attr_sct <- Seurat::SCTResults(
          object = model,
          slot = "feature.attributes"
        )
      }
      nfeatures <- min(nHVF, nrow(x = feature_attr_sct))
      top_features <- rownames(x = feature_attr_sct)[utils::head(
        order(feature_attr_sct$residual_variance, decreasing = TRUE),
        n = nfeatures
      )]
      SeuratObject::VariableFeatures(
        srt_list[[i]],
        assay = SeuratObject::DefaultAssay(srt_list[[i]])
      ) <- top_features
      srt_list[[i]] <- AddFeaturesData(
        srt_list[[i]],
        features = feature_attr_sct,
        assay = "SCT"
      )
    }
  }

  if (is.null(HVF)) {
    if (HVF_source == "global") {
      log_message(
        "Use the global HVF from merged dataset",
        verbose = verbose
      )
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
      log_message(
        "Use the separate HVF from srt_list",
        verbose = verbose
      )
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
      #   nHVF <- min(sapply(srt_list, function(srt) length(SeuratObject::VariableFeatures(srt))))
      #   HVF_sort <- sort(table(unlist(lapply(srt_list, VariableFeatures))), decreasing = TRUE)
      #   HVF_filter <- HVF_sort[HVF_sort >= HVF_min_intersection]
      #   HVF <- names(utils::head(HVF_filter, nHVF))
      # }
      if (length(HVF) == 0) {
        log_message(
          "No HVF available",
          message_type = "error",
          verbose = verbose
        )
      }
    }
  } else {
    cf <- Reduce(
      intersect,
      lapply(srt_list, function(srt) {
        rownames(
          GetAssayData5(
            srt,
            layer = "counts",
            assay = SeuratObject::DefaultAssay(srt)
          )
        )
      })
    )
    HVF <- HVF[HVF %in% cf]
  }
  log_message(
    "Number of available HVF: {.val {length(HVF)}}",
    verbose = verbose
  )

  hvf_sum <- lapply(srt_list, function(srt) {
    Matrix::colSums(
      GetAssayData5(
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
    log_message(
      "Some cells do not express any of the highly variable features: {.val {cell_abnormal}}",
      message_type = "warning",
      verbose = verbose
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
  log_message(
    "Finished check",
    verbose = verbose
  )

  return(
    list(
      srt_list = srt_list,
      HVF = HVF,
      assay = assay,
      type = type
    )
  )
}

#' @title Check and preprocess a merged seurat object
#'
#' @description
#' This function checks and preprocesses a merged seurat object.
#'
#' @inheritParams CheckDataList
#' @inheritParams integration_scop
#' @param srt_merge A merged `Seurat` object that includes the batch information.
#'
#' @seealso [CheckDataList]
#'
#' @export
CheckDataMerge <- function(
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
    verbose = TRUE,
    seed = 11) {
  if (!inherits(srt_merge, "Seurat")) {
    log_message(
      "{.arg srt_merge} is not a Seurat object",
      message_type = "error"
    )
  }
  if (length(batch) != 1) {
    log_message(
      "{.arg batch} must be provided to specify the batch column in meta.data",
      message_type = "error"
    )
  }
  if (!batch %in% colnames(srt_merge@meta.data)) {
    log_message(
      "No batch column: {.val {batch}} found in the meta.data",
      message_type = "error"
    )
  }
  if (!is.factor(srt_merge[[batch, drop = TRUE]])) {
    srt_merge[[batch, drop = TRUE]] <- factor(
      srt_merge[[batch, drop = TRUE]],
      levels = unique(srt_merge[[batch, drop = TRUE]])
    )
  }
  assay <- assay %||% SeuratObject::DefaultAssay(srt_merge)
  srt_merge_raw <- srt_merge

  log_message(
    "Spliting {.arg srt_merge} into {.arg srt_list} by column {.val {batch}}...",
    verbose = verbose
  )
  srt_list <- Seurat::SplitObject(
    object = srt_merge_raw,
    split.by = batch
  )

  checked <- CheckDataList(
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

  srt_merge <- srt_append(
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
