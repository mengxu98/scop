#' CellScoring
#'
#' This function performs cell scoring on a Seurat object.
#' It calculates scores for a given set of features and adds the scores as metadata to the Seurat object.
#'
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
#' @param BPPARAM The BiocParallel parameter object. Defaults to BiocParallel::bpparam().
#' @param seed The random seed for reproducibility. Defaults to 11.
#' @param ... Additional arguments to be passed to the scoring methods.
#'
#' @seealso \code{\link{PrepareDB}} \code{\link{ListDB}}
#'
#' @examples
#' data("pancreas_sub")
#' ccgenes <- CC_GenePrefetch("Mus_musculus")
#' pancreas_sub <- CellScoring(
#'   srt = pancreas_sub,
#'   features = list(S = ccgenes$S, G2M = ccgenes$G2M),
#'   method = "Seurat", name = "CC"
#' )
#' CellDimPlot(pancreas_sub, "CC_classification")
#' FeatureDimPlot(pancreas_sub, "CC_G2M")
#'
#' \dontrun{
#' data("panc8_sub")
#' panc8_sub <- Integration_scop(panc8_sub,
#'   batch = "tech", integration_method = "Seurat"
#' )
#' CellDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' panc8_sub <- CellScoring(
#'   srt = panc8_sub, layer = "data", assay = "RNA",
#'   db = "GO_BP", species = "Homo_sapiens",
#'   minGSSize = 10, maxGSSize = 100,
#'   method = "Seurat", name = "GO", new_assay = TRUE
#' )
#' panc8_sub <- Integration_scop(panc8_sub,
#'   assay = "GO",
#'   batch = "tech", integration_method = "Seurat"
#' )
#' CellDimPlot(panc8_sub, group.by = c("tech", "celltype"))
#'
#' pancreas_sub <- CellScoring(
#'   srt = pancreas_sub, layer = "data", assay = "RNA",
#'   db = "GO_BP", species = "Mus_musculus",
#'   termnames = panc8_sub[["GO"]]@meta.features[, "termnames"],
#'   method = "Seurat", name = "GO", new_assay = TRUE
#' )
#' pancreas_sub <- Standard_scop(pancreas_sub, assay = "GO")
#' CellDimPlot(pancreas_sub, "SubCellType")
#'
#' pancreas_sub[["tech"]] <- "Mouse"
#' panc_merge <- Integration_scop(
#'   srtList = list(panc8_sub, pancreas_sub),
#'   assay = "GO",
#'   batch = "tech", integration_method = "Seurat"
#' )
#' CellDimPlot(
#'   srt = panc_merge,
#'   group.by = c("tech", "celltype", "SubCellType", "Phase")
#' )
#'
#' genenames <- make.unique(
#'   capitalize(rownames(panc8_sub[["RNA"]]),
#'   force_tolower = TRUE
#' )
#' panc8_sub <- RenameFeatures(
#'   srt = panc8_sub,
#'   newnames = genenames,
#'   assay = "RNA"
#' )
#' head(rownames(panc8_sub))
#' panc_merge <- Integration_scop(
#'   srtList = list(panc8_sub, pancreas_sub),
#'   assay = "RNA",
#'   batch = "tech", integration_method = "Seurat"
#' )
#' CellDimPlot(
#'   srt = panc_merge,
#'   group.by = c("tech", "celltype", "SubCellType", "Phase")
#' )
#' }
#'
#' @importFrom BiocParallel bpprogressbar<- bpRNGseed<- bpworkers
#' @importFrom Seurat AddModuleScore AddMetaData
#' @export
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
    stop("method must be 'Seurat', 'AUCell'or 'UCell'.")
  }
  assay <- assay %||% DefaultAssay(srt)
  if (layer == "counts") {
    status <- check_DataType(srt, layer = "counts", assay = assay)
    if (status != "raw_counts") {
      warning("Data is not raw counts", immediate. = TRUE)
    }
  }
  if (layer == "data") {
    status <- check_DataType(srt, layer = "data", assay = assay)
    if (status == "raw_counts") {
      message(
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
      message(
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
      warning(
        "Can not determine whether data ",
        i,
        " is log-normalized...\n",
        immediate. = TRUE
      )
    }
  }
  if (name == "" && isTRUE(new_assay)) {
    stop("name must be specified when new_assay=TRUE")
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
      na <- rowSums(is.na(TERM2GENE_tmp)) > 0
      TERM2GENE_tmp <- TERM2GENE_tmp[!(dup | na), , drop = FALSE]
      TERM2NAME_tmp <- TERM2NAME_tmp[
        TERM2NAME_tmp[, "Term"] %in% TERM2GENE_tmp[, "Term"], ,
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
          stop("None of termnames found in the db: ", single_db)
        }
      }
      features <- c(features, features_tmp)
    }
  }

  if (!is.list(features) || length(names(features)) == 0) {
    stop("'features' must be a named list")
  }
  expressed <- names(
    which(
      rowSums(
        Seurat::GetAssayData(
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
    setNames(names(features), names(features)),
    function(x) {
      features[[x]][features[[x]] %in% expressed]
    }
  )
  filtered_none <- names(which(sapply(features, length) == 0))
  if (length(filtered_none) > 0) {
    warning(
      "The following list of features were filtered because none of features were found in the srt assay:\n",
      paste0(filtered_none, collapse = ", "),
      immediate. = TRUE
    )
  }
  features <- features[!names(features) %in% filtered_none]
  features_raw <- features
  names(features) <- make.names(names(features))
  message("Number of feature lists to be scored: ", length(features))

  time_start <- Sys.time()
  message(paste0("[", time_start, "] ", "Start CellScoring"))
  message("Workers: ", bpworkers(BPPARAM))

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
      check_R("UCell")
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
        warning(
          "The following list of features were filtered when scoring:\n",
          paste0(filtered, collapse = ", "),
          immediate. = TRUE
        )
      }
      features_keep <- features[!names(features) %in% filtered]
      features_nm <- features_raw[!names(features) %in% filtered]
      scores <- srt_tmp[[paste0(names(features_keep), name)]]
      colnames(scores) <- make.names(paste(name, names(features_nm), sep = "_"))
    } else if (method == "AUCell") {
      check_R("AUCell")
      CellRank <- AUCell::AUCell_buildRankings(
        Matrix::as.matrix(
          Seurat::GetAssayData(
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
        warning(
          "The following list of features were filtered when scoring:",
          paste0(filtered, collapse = ", "),
          immediate. = TRUE
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
    features_nm_list[[i]] <- setNames(
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
  message(paste0("[", time_end, "] ", "CellScoring done"))
  message(
    "Elapsed time:",
    format(
      round(difftime(time_end, time_start), 2),
      format = "%Y-%m-%d %H:%M:%S"
    )
  )

  return(srt)
}
