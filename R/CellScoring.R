#' CellScoring
#'
#' This function performs cell scoring on a Seurat object.
#' It calculates scores for a given set of features and adds the scores as metadata to the Seurat object.
#'
#' @md
#' @inheritParams RunEnrichment
#' @param srt A Seurat object
#' @param features A named list of feature lists for scoring. If NULLL, \code{db} will be used to create features sets.
#' @param layer The layer of the Seurat object to use for scoring. Defaults to "data".
#' @param assay The assay of the Seurat object to use for scoring. Defaults to NULL, in which case the default assay of the object is used.
#' @param split.by A cell metadata variable used for splitting the Seurat object into subsets and performing scoring on each subset. Defaults to NULL.
#' @param termnames A vector of term names to be used from the database. Defaults to NULL, in which case all features from the database are used.
#' @param method The method to use for scoring. Can be "Seurat", "AUCell", or "UCell". Defaults to "Seurat".
#' @param classification Whether to perform classification based on the scores. Defaults to TRUE.
#' @param name The name of the assay to store the scores in. Only used if new_assay is TRUE. Defaults to an empty string.
#' @param new_assay Whether to create a new assay for storing the scores. Defaults to FALSE.
#' @param BPPARAM The BiocParallel parameter object. Defaults to [BiocParallel::bpparam()].
#' @param seed The random seed for reproducibility. Defaults to 11.
#' @param ... Additional arguments to be passed to the scoring methods.
#'
#' @seealso \code{\link{PrepareDB}} \code{\link{ListDB}}
#'
#' @export
#'
#' @examples
#' data("pancreas_sub")
#' ccgenes <- CC_GenePrefetch("Mus_musculus")
#' pancreas_sub <- CellScoring(
#'   srt = pancreas_sub,
#'   features = list(S = ccgenes$S, G2M = ccgenes$G2M),
#'   method = "Seurat",
#'   name = "CC"
#' )
#' CellDimPlot(pancreas_sub, "CC_classification")
#' FeatureDimPlot(pancreas_sub, "CC_G2M")
#'
#' \dontrun{
#' data("panc8_sub")
#' panc8_sub <- integration_scop(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_method = "Seurat"
#' )
#' CellDimPlot(
#'   panc8_sub,
#'   group.by = c("tech", "celltype")
#'  )
#'
#' panc8_sub <- CellScoring(
#'   srt = panc8_sub,
#'   layer = "data",
#'   assay = "RNA",
#'   db = "GO_BP",
#'   species = "Homo_sapiens",
#'   minGSSize = 10,
#'   maxGSSize = 100,
#'   method = "Seurat",
#'   name = "GO",
#'   new_assay = TRUE
#' )
#'
#' panc8_sub <- integration_scop(
#'   panc8_sub,
#'   assay = "GO",
#'   batch = "tech",
#'   integration_method = "Seurat"
#' )
#' CellDimPlot(
#'   panc8_sub,
#'   group.by = c("tech", "celltype")
#'  )
#'
#' pancreas_sub <- CellScoring(
#'   srt = pancreas_sub,
#'   layer = "data",
#'   assay = "RNA",
#'   db = "GO_BP",
#'   species = "Mus_musculus",
#'   termnames = panc8_sub[["GO"]]@meta.features[, "termnames"],
#'   method = "Seurat",
#'   name = "GO",
#'   new_assay = TRUE
#' )
#' pancreas_sub <- standard_scop(
#'   pancreas_sub,
#'   assay = "GO"
#'  )
#' CellDimPlot(pancreas_sub, "SubCellType")
#'
#' pancreas_sub[["tech"]] <- "Mouse"
#' panc_merge <- integration_scop(
#'   srt_list = list(panc8_sub, pancreas_sub),
#'   assay = "GO",
#'   batch = "tech", integration_method = "Seurat"
#' )
#' CellDimPlot(
#'   srt = panc_merge,
#'   group.by = c("tech", "celltype", "SubCellType", "Phase")
#' )
#'
#' genenames <- make.unique(
#'   thisutils::capitalize(
#'     rownames(panc8_sub[["RNA"]])
#'   ),
#'   force_tolower = TRUE
#' )
#' names(genenames) <- rownames(panc8_sub)
#' panc8_sub <- RenameFeatures(
#'   srt = panc8_sub,
#'   newnames = genenames,
#'   assay = "RNA"
#' )
#' head(rownames(panc8_sub))
#' panc_merge <- integration_scop(
#'   srt_list = list(panc8_sub, pancreas_sub),
#'   assay = "RNA",
#'   batch = "tech", integration_method = "Seurat"
#' )
#' CellDimPlot(
#'   srt = panc_merge,
#'   group.by = c("tech", "celltype", "SubCellType", "Phase")
#' )
#' }
CellScoring <- function(
    srt,
    features = NULL,
    layer = "data",
    assay = NULL,
    split.by = NULL,
    IDtype = "symbol",
    species = "Homo_sapiens",
    db = "GO_BP",
    termnames = NULL,
    db_update = FALSE,
    db_version = "latest",
    convert_species = TRUE,
    Ensembl_version = 103,
    mirror = NULL,
    minGSSize = 10,
    maxGSSize = 500,
    method = "Seurat",
    classification = TRUE,
    name = "",
    new_assay = FALSE,
    BPPARAM = BiocParallel::bpparam(),
    seed = 11,
    ...) {
  set.seed(seed)
  bpprogressbar(BPPARAM) <- TRUE
  bpRNGseed(BPPARAM) <- seed

  if (!method %in% c("Seurat", "AUCell", "UCell")) {
    log_message(
      "method must be 'Seurat', 'AUCell'or 'UCell'.",
      message_type = "error"
    )
  }
  assay <- assay %||% DefaultAssay(srt)
  if (layer == "counts") {
    status <- check_data_type(srt, layer = "counts", assay = assay)
    if (status != "raw_counts") {
      log_message(
        "Data is not raw counts",
        message_type = "warning"
      )
    }
  }
  if (layer == "data") {
    status <- check_data_type(srt, layer = "data", assay = assay)
    if (status == "raw_counts") {
      log_message(
        "Data is raw counts. Perform NormalizeData(LogNormalize) on the data ..."
      )
      srt <- NormalizeData(
        object = srt,
        assay = assay,
        normalization.method = "LogNormalize",
        verbose = FALSE
      )
    }
    if (status == "raw_normalized_counts") {
      log_message(
        "Data is normalized without log transformation. Perform NormalizeData(LogNormalize) on the data..."
      )
      srt <- NormalizeData(
        object = srt,
        assay = assay,
        normalization.method = "LogNormalize",
        verbose = FALSE
      )
    }
    if (status == "unknown") {
      log_message(
        "Can not determine whether data ",
        i,
        " is log-normalized...\n",
        message_type = "warning"
      )
    }
  }
  if (name == "" && isTRUE(new_assay)) {
    log_message(
      "name must be specified when new_assay=TRUE",
      message_type = "error"
    )
  }
  if (is.null(features)) {
    for (single_db in db) {
      db_list <- PrepareDB(
        species = species,
        db = single_db,
        db_update = db_update,
        db_version = db_version,
        db_IDtypes = IDtype,
        convert_species = convert_species,
        Ensembl_version = Ensembl_version,
        mirror = mirror
      )
      TERM2GENE_tmp <- db_list[[species]][[single_db]][["TERM2GENE"]][, c(
        "Term",
        IDtype
      )]
      TERM2NAME_tmp <- db_list[[species]][[single_db]][["TERM2NAME"]]
      dup <- duplicated(TERM2GENE_tmp)
      na <- Matrix::rowSums(is.na(TERM2GENE_tmp)) > 0
      TERM2GENE_tmp <- TERM2GENE_tmp[!(dup | na), , drop = FALSE]
      TERM2NAME_tmp <- TERM2NAME_tmp[
        TERM2NAME_tmp[, "Term"] %in% TERM2GENE_tmp[, "Term"],
        ,
        drop = FALSE
      ]

      TERM2GENE_tmp <- unique(TERM2GENE_tmp)
      TERM2NAME_tmp <- unique(TERM2NAME_tmp)
      rownames(TERM2NAME_tmp) <- TERM2NAME_tmp[, "Term"]
      features_tmp <- split(
        TERM2GENE_tmp[, IDtype],
        TERM2NAME_tmp[TERM2GENE_tmp[, "Term"], "Name"]
      )
      if (is.null(termnames)) {
        GSSize <- sapply(features_tmp, length)
        features_tmp <- features_tmp[GSSize >= minGSSize & GSSize <= maxGSSize]
      } else {
        if (length(intersect(termnames, names(features_tmp)) > 0)) {
          features_tmp <- features_tmp[intersect(
            termnames,
            names(features_tmp)
          )]
        } else {
          log_message(
            "None of termnames found in the db: ", single_db,
            message_type = "error"
          )
        }
      }
      features <- c(features, features_tmp)
    }
  }

  if (!is.list(features) || length(names(features)) == 0) {
    log_message(
      "'features' must be a named list",
      message_type = "error"
    )
  }
  expressed <- names(
    which(
      Matrix::rowSums(
        GetAssayData5(
          srt,
          layer = layer,
          assay = assay
        ) >
          0
      ) >
        0
    )
  )
  features <- lapply(
    stats::setNames(names(features), names(features)),
    function(x) {
      features[[x]][features[[x]] %in% expressed]
    }
  )
  filtered_none <- names(which(sapply(features, length) == 0))
  if (length(filtered_none) > 0) {
    log_message(
      paste0(
        "The following list of features were filtered because none of features were found in the srt assay:\n",
        paste0(filtered_none, collapse = ", ")
      ),
      message_type = "warning"
    )
  }
  features <- features[!names(features) %in% filtered_none]
  features_raw <- features
  names(features) <- make.names(names(features))
  log_message(
    "Number of feature lists to be scored: ", length(features)
  )

  time_start <- Sys.time()
  log_message("Start CellScoring")
  log_message("Workers: ", bpworkers(BPPARAM))

  if (!is.null(split.by)) {
    split_list <- Seurat::SplitObject(srt, split.by = split.by)
  } else {
    split_list <- list(srt)
  }
  scores_list <- list()
  features_nm_list <- list()
  for (i in seq_along(split_list)) {
    srt_sp <- split_list[[i]]
    if (method == "Seurat") {
      srt_tmp <- AddModuleScore2(
        srt_sp,
        features = features,
        name = name,
        layer = layer,
        assay = assay,
        BPPARAM = BPPARAM,
        ...
      )
      if (name != "") {
        scores <- srt_tmp[[paste0(name, seq_along(features))]]
      } else {
        scores <- srt_tmp[[paste0("X", seq_along(features))]]
      }
      features_nm <- features_raw
      colnames(scores) <- make.names(paste(name, names(features_nm), sep = "_"))
    } else if (method == "UCell") {
      check_r("UCell")
      srt_tmp <- UCell::AddModuleScore_UCell(
        srt_sp,
        features = features,
        name = name,
        slot = layer,
        assay = assay,
        BPPARAM = BPPARAM,
        ...
      )
      filtered <- names(features)[
        !paste0(names(features), name) %in% colnames(srt_tmp@meta.data)
      ]
      if (length(filtered) > 0) {
        log_message(
          paste0(
            "The following list of features were filtered when scoring:\n",
            paste0(filtered, collapse = ", ")
          ),
          message_type = "warning"
        )
      }
      features_keep <- features[!names(features) %in% filtered]
      features_nm <- features_raw[!names(features) %in% filtered]
      scores <- srt_tmp[[paste0(names(features_keep), name)]]
      colnames(scores) <- make.names(paste(name, names(features_nm), sep = "_"))
    } else if (method == "AUCell") {
      check_r("AUCell")
      CellRank <- AUCell::AUCell_buildRankings(
        Matrix::as.matrix(
          GetAssayData5(
            srt_sp,
            layer = layer,
            assay = assay
          )
        ),
        BPPARAM = BPPARAM,
        plotStats = FALSE
      )
      cells_AUC <- AUCell::AUCell_calcAUC(
        geneSets = features,
        rankings = CellRank,
        ...
      )
      filtered <- names(features)[
        !names(features) %in% rownames(AUCell::getAUC(cells_AUC))
      ]
      if (length(filtered) > 0) {
        log_message(
          "The following list of features were filtered when scoring:",
          paste0(filtered, collapse = ", "),
          message_type = "warning"
        )
      }
      features_keep <- features[!names(features) %in% filtered]
      features_nm <- features_raw[!names(features) %in% filtered]
      scores <- as.data.frame(
        Matrix::t(
          AUCell::getAUC(cells_AUC)
        )
      )[, names(features_keep)]
      colnames(scores) <- make.names(
        paste(name, names(features_nm), sep = "_")
      )
    }
    features_nm_list[[i]] <- stats::setNames(
      object = names(features_nm),
      nm = colnames(scores)
    )
    scores_list[[i]] <- scores
  }
  features_used <- Reduce(intersect, lapply(scores_list, colnames))
  features_nm_used <- Reduce(intersect, features_nm_list)
  scores_mat <- do.call(
    rbind,
    lapply(
      scores_list,
      function(x) {
        x[intersect(rownames(x), colnames(srt)), features_used, drop = FALSE]
      }
    )
  )

  if (isTRUE(new_assay)) {
    srt[[name]] <- Seurat::CreateAssayObject(
      counts = Matrix::t(
        Matrix::as.matrix(
          scores_mat[colnames(srt), , drop = FALSE]
        )
      )
    )
    srt[[name]] <- Seurat::AddMetaData(
      object = srt[[name]],
      metadata = data.frame(
        termnames = features_nm_used[colnames(scores_mat)]
      )
    )
  } else {
    srt <- Seurat::AddMetaData(
      object = srt,
      metadata = scores_mat
    )
  }

  if (isTRUE(classification)) {
    assignments <- apply(
      scores_mat,
      MARGIN = 1,
      FUN = function(x) {
        if (all(x < 0)) {
          return(NA)
        } else {
          if (length(which(x == max(x))) > 1) {
            return(NA)
          } else {
            return(names(features)[which(x == max(x))])
          }
        }
      }
    )
    srt[[paste0(name, "_classification")]] <- assignments[rownames(scores_mat)]
  }

  time_end <- Sys.time()
  log_message("CellScoring done")
  log_message(
    "Elapsed time: ",
    format(
      round(difftime(time_end, time_start), 2),
      format = "%Y-%m-%d %H:%M:%S"
    )
  )

  return(srt)
}

AddModuleScore2 <- function(
    object,
    layer = "data",
    features,
    pool = NULL,
    nbin = 24,
    ctrl = 100,
    k = FALSE,
    assay = NULL,
    name = "Cluster",
    seed = 1,
    search = FALSE,
    BPPARAM = BiocParallel::bpparam(),
    ...) {
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  assay.old <- SeuratObject::DefaultAssay(object = object)
  assay <- assay %||% assay.old
  SeuratObject::DefaultAssay(object = object) <- assay
  assay.data <- GetAssayData5(object = object, layer = layer)
  features.old <- features
  if (k) {
    .NotYetUsed(arg = "k")
    features <- list()
    for (i in as.numeric(
      x = names(x = table(object@kmeans.obj[[1]]$cluster))
    )) {
      features[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == i))
    }
    cluster.length <- length(x = features)
  } else {
    if (is.null(x = features)) {
      log_message(
        "Missing input feature list",
        message_type = "error"
      )
    }
    features <- lapply(X = features, FUN = function(x) {
      missing.features <- setdiff(x = x, y = rownames(x = object))
      if (length(x = missing.features) > 0) {
        log_message(
          "The following features are not present in the object: ",
          paste(missing.features, collapse = ", "),
          ifelse(
            test = search,
            yes = ", attempting to find updated synonyms",
            no = ", not searching for symbol synonyms"
          ),
          message_type = "warning"
        )
        if (search) {
          tryCatch(
            expr = {
              updated.features <- Seurat::UpdateSymbolList(
                symbols = missing.features,
                ...
              )
              names(x = updated.features) <- missing.features
              for (miss in names(x = updated.features)) {
                index <- which(x == miss)
                x[index] <- updated.features[miss]
              }
            },
            error = function(...) {
              log_message(
                "Could not reach HGNC's gene names database",
                message_type = "warning"
              )
            }
          )
          missing.features <- setdiff(x = x, y = rownames(x = object))
          if (length(x = missing.features) > 0) {
            log_message(
              "The following features are still not present in the object: ",
              paste(missing.features, collapse = ", "),
              message_type = "warning"
            )
          }
        }
      }
      return(intersect(x = x, y = rownames(x = object)))
    })
    cluster.length <- length(x = features)
  }
  if (!all(check_length(values = features))) {
    log_message(paste(
      "Could not find enough features in the object from the following feature lists:",
      paste(names(x = which(x = !check_length(values = features)))),
      "Attempting to match case..."
    ),
    message_type = "warning"
  )
    features <- lapply(
      X = features.old,
      FUN = Seurat::CaseMatch,
      match = rownames(x = object)
    )
  }
  if (!all(check_length(values = features))) {
    log_message(paste(
      "The following feature lists do not have enough features present in the object:",
      paste(names(x = which(x = !check_length(values = features)))),
      "exiting..."
    ),
    message_type = "error"
  )
  }
  pool <- pool %||% rownames(x = object)
  data_avg <- Matrix::rowMeans(
    x = assay.data[pool, , drop = FALSE]
  )
  data_avg <- data_avg[order(data_avg)]
  data.cut <- cut_number(
    x = data_avg + stats::rnorm(n = length(data_avg)) / 1e+30,
    n = nbin,
    labels = FALSE,
    right = FALSE
  )
  names(x = data.cut) <- names(x = data_avg)

  scores <- bplapply(
    1:cluster.length,
    function(i) {
      features.use <- features[[i]]
      ctrl.use <- unlist(
        lapply(
          1:length(features.use),
          function(j) {
            data.cut[which(data.cut == data.cut[features.use[j]])]
          }
        )
      )
      ctrl.use <- names(
        sample(
          ctrl.use,
          size = min(ctrl * length(features.use), length(ctrl.use)),
          replace = FALSE
        )
      )
      ctrl.scores_i <- Matrix::colMeans(
        x = assay.data[ctrl.use, , drop = FALSE]
      )
      features.scores_i <- Matrix::colMeans(
        x = assay.data[features.use, , drop = FALSE]
      )
      return(list(ctrl.scores_i, features.scores_i))
    },
    BPPARAM = BPPARAM
  )
  ctrl.scores <- do.call(rbind, lapply(scores, function(x) x[[1]]))
  features.scores <- do.call(rbind, lapply(scores, function(x) x[[2]]))

  features.scores.use <- features.scores - ctrl.scores
  rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  rownames(x = features.scores.use) <- colnames(x = object)
  object[[colnames(x = features.scores.use)]] <- features.scores.use
  SeuratObject::CheckGC()
  SeuratObject::DefaultAssay(object = object) <- assay.old

  return(object)
}

check_length <- function(
    values,
    cutoff = 0) {
  return(
    vapply(
      X = values,
      FUN = function(x) {
        return(length(x = x) > cutoff)
      },
      FUN.VALUE = logical(1)
    )
  )
}
