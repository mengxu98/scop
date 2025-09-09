#' @title Cell scoring
#'
#' @description
#' This function performs cell scoring on a Seurat object.
#' It calculates scores for a given set of features and adds the scores as metadata to the Seurat object.
#'
#' @md
#' @inheritParams RunEnrichment
#' @param srt A Seurat object.
#' @param features A named list of feature lists for scoring. If NULL, \code{db} will be used to create features sets.
#' @param layer The layer of the Seurat object to use for scoring. Defaults to "data".
#' @param assay The assay of the Seurat object to use for scoring. Defaults to NULL, in which case the default assay of the object is used.
#' @param split.by A cell metadata variable used for splitting the Seurat object into subsets and performing scoring on each subset. Defaults to NULL.
#' @param termnames A vector of term names to be used from the database. Defaults to NULL, in which case all features from the database are used.
#' @param method The method to use for scoring. Can be "Seurat", "AUCell", or "UCell". Defaults to "Seurat".
#' @param classification Whether to perform classification based on the scores. Defaults to TRUE.
#' @param name The name of the assay to store the scores in. Only used if new_assay is TRUE. Defaults to an empty string.
#' @param new_assay Whether to create a new assay for storing the scores. Defaults to FALSE.
#' @param seed The random seed for reproducibility. Defaults to 11.
#' @inheritParams thisutils::parallelize_fun
#' @param ... Additional arguments to be passed to the scoring methods.
#'
#' @seealso
#' [PrepareDB], [ListDB], [RunDynamicFeatures]
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' features_all <- rownames(pancreas_sub)
#' pancreas_sub <- CellScoring(
#'   pancreas_sub,
#'   features = list(
#'     A = features_all[1:100],
#'     B = features_all[101:200]
#'   ),
#'   method = "Seurat",
#'   name = "test"
#' )
#' CellDimPlot(pancreas_sub, "test_classification")
#'
#' FeatureDimPlot(pancreas_sub, "test_A")
#'
#' \dontrun{
#' data(panc8_sub)
#' panc8_sub <- integration_scop(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_method = "Seurat"
#' )
#' CellDimPlot(
#'   panc8_sub,
#'   group.by = c("tech", "celltype")
#' )
#'
#' panc8_sub <- CellScoring(
#'   panc8_sub,
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
#' )
#'
#' pancreas_sub <- CellScoring(
#'   pancreas_sub,
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
#' )
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
#'   panc8_sub,
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
    Ensembl_version = NULL,
    mirror = NULL,
    minGSSize = 10,
    maxGSSize = 500,
    method = "Seurat",
    classification = TRUE,
    name = "",
    new_assay = FALSE,
    seed = 11,
    cores = 1,
    verbose = TRUE,
    ...) {
  log_message(
    "Start cell scoring",
    verbose = verbose
  )
  set.seed(seed)

  score_methods <- c("Seurat", "AUCell", "UCell")
  if (!method %in% score_methods) {
    log_message(
      "{.arg method} must be one of {.val {score_methods}}",
      message_type = "error",
      verbose = verbose
    )
  }
  assay <- assay %||% DefaultAssay(srt)
  if (layer == "counts") {
    status <- CheckDataType(srt, layer = "counts", assay = assay)
    if (status != "raw_counts") {
      log_message(
        "Data is not raw counts",
        message_type = "warning",
        verbose = verbose
      )
    }
  }
  if (layer == "data") {
    status <- CheckDataType(srt, layer = "data", assay = assay)
    if (status == "raw_counts") {
      log_message(
        "Perform {.fn NormalizeData} with {.arg normalization.method = 'LogNormalize'} on the data"
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
        "Perform {.fn NormalizeData} with {.arg normalization.method = 'LogNormalize'} on the data"
      )
      srt <- NormalizeData(
        object = srt,
        assay = assay,
        normalization.method = "LogNormalize",
        verbose = FALSE
      )
    }
  }
  if (name == "" && isTRUE(new_assay)) {
    log_message(
      "{.arg name} must be specified when {.arg new_assay = TRUE}",
      message_type = "error",
      verbose = verbose
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
      db_data <- db_list[[species]][[single_db]]
      term2gene_tmp <- db_data[["TERM2GENE"]][, c("Term", IDtype)]
      term2name_tmp <- db_data[["TERM2NAME"]]
      dup <- duplicated(term2gene_tmp)
      na <- Matrix::rowSums(is.na(term2gene_tmp)) > 0
      term2gene_tmp <- term2gene_tmp[!(dup | na), , drop = FALSE]
      term2name_tmp <- term2name_tmp[
        term2name_tmp[, "Term"] %in% term2gene_tmp[, "Term"], ,
        drop = FALSE
      ]

      term2gene_tmp <- unique(term2gene_tmp)
      term2name_tmp <- unique(term2name_tmp)
      rownames(term2name_tmp) <- term2name_tmp[, "Term"]
      features_tmp <- split(
        term2gene_tmp[, IDtype],
        term2name_tmp[term2gene_tmp[, "Term"], "Name"]
      )
      if (is.null(termnames)) {
        gssize <- sapply(features_tmp, length)
        features_tmp <- features_tmp[gssize >= minGSSize & gssize <= maxGSSize]
      } else {
        if (length(intersect(termnames, names(features_tmp)) > 0)) {
          features_tmp <- features_tmp[intersect(
            termnames,
            names(features_tmp)
          )]
        } else {
          log_message(
            "None of termnames found in the db: {.val {single_db}}",
            message_type = "error",
            verbose = verbose
          )
        }
      }
      features <- c(features, features_tmp)
    }
  }

  if (!is.list(features) || length(names(features)) == 0) {
    log_message(
      "{.arg features} must be a named list",
      message_type = "error",
      verbose = verbose
    )
  }
  expr_data <- GetAssayData5(
    srt,
    layer = layer,
    assay = assay,
    verbose = FALSE
  )
  features_expressed <- names(
    which(
      Matrix::rowSums(expr_data > 0) > 0
    )
  )
  features <- lapply(
    stats::setNames(names(features), names(features)),
    function(x) {
      features[[x]][features[[x]] %in% features_expressed]
    }
  )
  features_none <- names(which(sapply(features, length) == 0))
  if (length(features_none) > 0) {
    log_message(
      "The following features were filtered because not found in the srt assay: {.val {features_none}}",
      message_type = "warning",
      verbose = verbose
    )
  }
  features <- features[!names(features) %in% features_none]
  features_raw <- features
  names(features) <- make.names(names(features))
  log_message(
    "Number of feature lists to be scored: {.val {length(features)}}",
    verbose = verbose
  )

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
        cores = cores,
        verbose = verbose,
        ...
      )
      if (name != "") {
        scores <- srt_tmp[[paste0(name, seq_along(features))]]
      } else {
        scores <- srt_tmp[[paste0("X", seq_along(features))]]
      }
      features_nm <- features_raw
    } else if (method == "UCell") {
      check_r("UCell")
      srt_tmp <- UCell::AddModuleScore_UCell(
        srt_sp,
        features = features,
        name = name,
        slot = layer,
        assay = assay,
        ...
      )
      filtered <- names(features)[
        !paste0(names(features), name) %in% colnames(srt_tmp@meta.data)
      ]
      if (length(filtered) > 0) {
        log_message(
          "The following features were filtered when scoring: {.val {filtered}}",
          message_type = "warning",
          verbose = verbose
        )
      }
      features_keep <- features[!names(features) %in% filtered]
      features_nm <- features_raw[!names(features) %in% filtered]
      scores <- srt_tmp[[paste0(names(features_keep), name)]]
    } else if (method == "AUCell") {
      check_r("AUCell")
      cell_rank <- AUCell::AUCell_buildRankings(
        as_matrix(
          GetAssayData5(
            srt_sp,
            layer = layer,
            assay = assay,
            verbose = FALSE
          )
        ),
        plotStats = FALSE
      )
      cells_auc <- AUCell::AUCell_calcAUC(
        geneSets = features,
        rankings = cell_rank,
        ...
      )
      filtered <- names(features)[
        !names(features) %in% rownames(AUCell::getAUC(cells_auc))
      ]
      if (length(filtered) > 0) {
        log_message(
          "The following features were filtered when scoring: {.val {filtered}}",
          message_type = "warning",
          verbose = verbose
        )
      }
      features_keep <- features[!names(features) %in% filtered]
      features_nm <- features_raw[!names(features) %in% filtered]
      scores <- as.data.frame(
        Matrix::t(
          AUCell::getAUC(cells_auc)
        )
      )[, names(features_keep)]
    }
    colnames(scores) <- make.names(
      paste(name, names(features_nm), sep = "_")
    )
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
        as_matrix(
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

  log_message(
    "Cell scoring completed",
    message_type = "success",
    verbose = verbose
  )

  return(srt)
}

AddModuleScore2 <- function(
    object,
    features,
    assay = NULL,
    layer = "data",
    pool = NULL,
    nbin = 24,
    ctrl = 100,
    name = "Cluster",
    seed = 1,
    search = FALSE,
    cores = 1,
    verbose = TRUE,
    ...) {
  set.seed(seed = seed)
  assay_raw <- SeuratObject::DefaultAssay(object = object)
  assay <- assay %||% assay_raw
  SeuratObject::DefaultAssay(object = object) <- assay
  expr_data <- GetAssayData5(
    object = object,
    layer = layer,
    verbose = FALSE
  )
  features_raw <- features

  features <- lapply(
    X = features, FUN = function(x) {
      features_missing <- setdiff(x, y = rownames(object))
      if (length(features_missing) > 0) {
        log_message(
          "The following features were not found in the object: {.val {features_missing}}",
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
              features_updated <- Seurat::UpdateSymbolList(
                symbols = features_missing,
                ...
              )
              names(features_updated) <- features_missing
              for (miss in names(features_updated)) {
                index <- which(x == miss)
                x[index] <- features_updated[miss]
              }
            },
            error = function(...) {
              log_message(
                "Could not reach HGNC's gene names database",
                message_type = "warning"
              )
            }
          )
          features_missing <- setdiff(x, y = rownames(object))
          if (length(features_missing) > 0) {
            log_message(
              "The following features were not found in the object: {.val {features_missing}}",
              message_type = "warning"
            )
          }
        }
      }
      intersect(x, y = rownames(object))
    }
  )

  cluster_length <- length(features)
  if (!all(check_length(values = features))) {
    log_message(
      "Could not find enough features in the object from the following feature lists: {.val {names(which(!check_length(values = features)))}}\n",
      "Attempting to match case...",
      message_type = "warning"
    )
    features <- lapply(
      X = features_raw,
      FUN = Seurat::CaseMatch,
      match = rownames(object)
    )
  }
  if (!all(check_length(values = features))) {
    log_message(
      "The following feature lists do not have enough features present in the object: {.val {names(which(!check_length(values = features)))}}",
      message_type = "error"
    )
  }
  pool <- pool %||% rownames(object)
  data_avg <- Matrix::rowMeans(
    expr_data[pool, , drop = FALSE]
  )
  data_avg <- data_avg[order(data_avg)]
  data_cut <- ggplot2::cut_number(
    data_avg + stats::rnorm(n = length(data_avg)) / 1e+30,
    n = nbin,
    labels = FALSE,
    right = FALSE
  )
  names(data_cut) <- names(data_avg)

  scores <- parallelize_fun(
    seq_len(cluster_length),
    function(i) {
      features_use <- features[[i]]
      ctrl_use <- unlist(
        lapply(
          seq_len(length(features_use)),
          function(j) {
            data_cut[which(data_cut == data_cut[features_use[j]])]
          }
        )
      )
      ctrl_use <- names(
        sample(
          ctrl_use,
          size = min(ctrl * length(features_use), length(ctrl_use)),
          replace = FALSE
        )
      )
      ctrl_scores_i <- Matrix::colMeans(
        expr_data[ctrl_use, , drop = FALSE]
      )
      features_scores_i <- Matrix::colMeans(
        expr_data[features_use, , drop = FALSE]
      )
      list(ctrl_scores_i, features_scores_i)
    },
    cores = cores,
    verbose = verbose
  )
  ctrl_scores <- do.call(rbind, lapply(scores, function(x) x[[1]]))
  features_scores <- do.call(rbind, lapply(scores, function(x) x[[2]]))

  features_scores_use <- features_scores - ctrl_scores
  rownames(features_scores_use) <- paste0(name, seq_len(cluster_length))
  features_scores_use <- as.data.frame(t(features_scores_use))
  rownames(features_scores_use) <- colnames(object)
  object[[colnames(features_scores_use)]] <- features_scores_use
  SeuratObject::CheckGC()
  SeuratObject::DefaultAssay(object = object) <- assay_raw

  return(object)
}

check_length <- function(
    values,
    cutoff = 0) {
  vapply(
    X = values,
    FUN = function(x) {
      length(x) > cutoff
    },
    FUN.VALUE = logical(1)
  )
}
