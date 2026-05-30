#' @title Cell scoring
#'
#' @description
#' This function performs cell scoring on a Seurat object.
#' It calculates scores for a given set of features and stores them in `meta.data`
#' and/or a score assay, depending on `new_assay` and `store_metadata`.
#'
#' @md
#' @inheritParams thisutils::parallelize_fun
#' @inheritParams CellDimPlot
#' @inheritParams RunEnrichment
#' @inheritParams standard_scop
#' @inheritParams FeatureDimPlot
#' @inheritParams RunGSVA
#' @param features A named list of feature lists for scoring.
#' If `NULL`, `db` will be used to create features sets.
#' @param termnames A vector of term names to be used from the database. Default is `NULL`, in which case all features from the database are used.
#' @param method The method to use for scoring. Can be `"Seurat"`, `"AUCell"`,
#' `"UCell"`, `"GSVA"`, `"ssGSEA"`, `"zscore"`, `"PLAGE"`, or `"VISION"`.
#' Multiple methods can be supplied at once; in that case each method is run
#' separately and stored with a method suffix such as `"GO_AUCell"` or
#' `"GO_GSVA"`. Default is `"Seurat"`.
#' @param backend Scoring backend. `"cpp"` is the default for supported methods.
#' `"r"` uses the original package implementation. `"cpp"` currently supports
#' `method = "Seurat"`, `method = "AUCell"`, `method = "GSVA"`,
#' `method = "ssGSEA"`, `method = "zscore"`, and `method = "PLAGE"`.
#' `method = "UCell"` and `method = "VISION"` fall back to `"r"` when
#' `backend` is not explicitly set.
#' @param cpp_strategy C++ AUCell ranking strategy. `"sparse"` ranks non-zero
#' genes and approximates zero ties, `"topk"` ranks only genes that can contribute
#' to AUCell AUC, and `"full"` ranks all genes.
#' @param classification Whether to perform classification based on the scores. Default is `TRUE`.
#' @param name The name of the assay to store the scores in. Only used if new_assay is TRUE. Default is `""`.
#' @param new_assay Whether to create a new assay for storing the scores. Default is `FALSE`.
#' @param store_metadata Whether to also store score columns in `meta.data`.
#' When `NULL`, manual `features = list(...)` input is stored in `meta.data` by default,
#' while database-derived results stay assay-only when `new_assay = TRUE`.
#' @param ... Additional arguments to be passed to the scoring methods.
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
#'   method = "AUCell",
#'   name = "test"
#' )
#' CellDimPlot(pancreas_sub, "test_classification")
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   features = "test_A"
#' )
#'
#' pancreas_sub <- CellScoring(
#'   pancreas_sub,
#'   features = list(A = features_all[1:100]),
#'   method = c("AUCell", "GSVA")
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   stat.by = c("AUCell_A", "GSVA_A"),
#'   group.by = "CellType",
#'   plot.by = "feature",
#'   plot_type = "violin",
#'   stack = TRUE
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   features = c("AUCell_A", "GSVA_A"),
#'   xlab = "UMAP_1",
#'   ylab = "UMAP_2"
#' )
#'
#' GroupHeatmap(
#'  pancreas_sub,
#'  features = c("AUCell_A", "GSVA_A", "Sox9", "Anxa2", "Bicc1"),
#'  group.by = "CellType"
#' )
#'
#' data(panc8_sub)
#' panc8_sub <- integration_scop(
#'   panc8_sub,
#'   batch = "tech",
#'   integration_method = "Harmony"
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
#'   method = "AUCell",
#'   name = "GO",
#'   new_assay = TRUE
#' )
#'
#' panc8_sub <- integration_scop(
#'   panc8_sub,
#'   assay = "GO",
#'   batch = "tech",
#'   integration_method = "Harmony"
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
#'   method = "AUCell",
#'   name = "GO",
#'   new_assay = TRUE
#' )
#' pancreas_sub <- standard_scop(
#'   pancreas_sub,
#'   assay = "GO"
#' )
#'
#' pancreas_sub[["tech"]] <- "Mouse"
#' panc_merge <- integration_scop(
#'   srt_list = list(panc8_sub, pancreas_sub),
#'   assay = "GO",
#'   batch = "tech",
#'   integration_method = "Harmony"
#' )
#' CellDimPlot(
#'   srt = panc_merge,
#'   group.by = c("tech", "celltype", "SubCellType", "Phase")
#' )
#'
#' genenames <- make.unique(
#'   thisutils::capitalize(
#'     rownames(panc8_sub[["RNA"]]),
#'     force_tolower = TRUE
#'   )
#' )
#' names(genenames) <- rownames(panc8_sub)
#' panc8_sub <- RenameFeatures(
#'   panc8_sub,
#'   newnames = genenames,
#'   assay = "RNA"
#' )
#' panc_merge <- integration_scop(
#'   srt_list = list(panc8_sub, pancreas_sub),
#'   assay = "RNA",
#'   batch = "tech",
#'   integration_method = "Harmony"
#' )
#' CellDimPlot(
#'   srt = panc_merge,
#'   group.by = c("tech", "celltype", "SubCellType", "Phase")
#' )
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
  backend = c("cpp", "r"),
  cpp_strategy = c("sparse", "topk", "full"),
  classification = TRUE,
  name = "",
  new_assay = FALSE,
  store_metadata = NULL,
  seed = 11,
  cores = 1,
  verbose = TRUE,
  ...
) {
  log_message(
    "Start cell scoring",
    verbose = verbose
  )
  set.seed(seed)

  method <- unique(vapply(
    as.character(method),
    normalize_gene_set_scoring_method,
    FUN.VALUE = character(1),
    arg_name = "method"
  ))
  if (length(method) > 1L) {
    for (method_i in method) {
      name_i <- if (nzchar(name)) {
        paste(name, method_i, sep = "_")
      } else {
        method_i
      }
      srt <- CellScoring(
        srt = srt,
        features = features,
        layer = layer,
        assay = assay,
        split.by = split.by,
        IDtype = IDtype,
        species = species,
        db = db,
        termnames = termnames,
        db_update = db_update,
        db_version = db_version,
        convert_species = convert_species,
        Ensembl_version = Ensembl_version,
        mirror = mirror,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        method = method_i,
        backend = backend,
        cpp_strategy = cpp_strategy,
        classification = classification,
        name = name_i,
        new_assay = new_assay,
        store_metadata = store_metadata,
        seed = seed,
        cores = cores,
        verbose = verbose,
        ...
      )
    }
    return(srt)
  }
  method <- method[[1]]
  backend_missing <- missing(backend)
  backend <- match.arg(backend)
  cpp_strategy <- match.arg(cpp_strategy)
  cpp_supported_methods <- c("Seurat", "AUCell", "GSVA", "ssGSEA", "zscore", "PLAGE")
  if (!identical(backend, "r") && !method %in% cpp_supported_methods) {
    if (isTRUE(backend_missing)) {
      log_message(
        "{.arg method = {.val {method}}} does not have a C++ backend yet; using {.arg backend = 'r'} for this run.",
        message_type = "warning",
        verbose = verbose
      )
      backend <- "r"
    } else {
      log_message(
        "{.arg backend = 'cpp'} currently supports {.arg method} values {.val {cpp_supported_methods}} only",
        message_type = "error",
        verbose = verbose
      )
    }
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
        "Perform {.fn NormalizeData} with {.arg normalization.method = 'LogNormalize'} on {.arg srt}"
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
        "Perform {.fn NormalizeData} with {.arg normalization.method = 'LogNormalize'} on {.arg srt}"
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
  features_input_supplied <- !is.null(features)
  if (is.null(store_metadata)) {
    store_metadata <- isTRUE(features_input_supplied) || !isTRUE(new_assay)
  } else {
    if (!is.logical(store_metadata) || length(store_metadata) != 1L || is.na(store_metadata)) {
      log_message(
        "{.arg store_metadata} must be TRUE, FALSE, or NULL",
        message_type = "error",
        verbose = verbose
      )
    }
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
        mirror = mirror,
        ...
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
  dots <- list(...)
  kcdf <- dots[["kcdf"]] %||% if (identical(layer, "counts")) "Poisson" else "Gaussian"
  kcdf <- match.arg(kcdf, c("Gaussian", "Poisson"))
  cpp_chunk_size <- dots[["cpp_chunk_size"]] %||% NULL
  abs.ranking <- dots[["abs.ranking"]] %||% FALSE
  min.sz <- dots[["min.sz"]] %||% minGSSize
  max.sz <- dots[["max.sz"]] %||% maxGSSize
  mx.diff <- dots[["mx.diff"]] %||% TRUE
  tau <- dots[["tau"]] %||% 1
  ssgsea.norm <- dots[["ssgsea.norm"]] %||% TRUE
  expr_data <- GetAssayData5(
    srt,
    layer = layer,
    assay = assay
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
  names(features) <- make.unique(make.names(names(features)))
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
      if (identical(backend, "cpp")) {
        expr_sp <- GetAssayData5(
          srt_sp,
          layer = layer,
          assay = assay
        )
        module_scores <- run_seurat_module_scores(
          expr_data = expr_sp,
          features = features,
          pool = dots[["pool"]] %||% NULL,
          nbin = dots[["nbin"]] %||% 24,
          ctrl = dots[["ctrl"]] %||% 100,
          seed = seed
        )
        filtered <- names(features)[
          !names(features) %in% colnames(module_scores)
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
          as_matrix(module_scores)
        )[, names(features_keep), drop = FALSE]
      } else {
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
      }
    } else if (method == "UCell") {
      check_r("UCell", verbose = FALSE)
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
      expr_sp <- GetAssayData5(
        srt_sp,
        layer = layer,
        assay = assay
      )
      if (identical(backend, "cpp")) {
        auc_scores <- run_aucell_scores(
          expr_counts = expr_sp,
          gene_sets = features,
          strategy = cpp_strategy
        )
        filtered <- names(features)[
          !names(features) %in% colnames(auc_scores)
        ]
      } else {
        gene_set_scoring_require_namespace("AUCell")
        cell_rank <- AUCell::AUCell_buildRankings(
          as_matrix(expr_sp),
          plotStats = FALSE
        )
        cells_auc <- AUCell::AUCell_calcAUC(
          geneSets = features,
          rankings = cell_rank,
          ...
        )
        auc_scores <- Matrix::t(AUCell::getAUC(cells_auc))
        filtered <- names(features)[
          !names(features) %in% colnames(auc_scores)
        ]
      }
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
        as_matrix(auc_scores)
      )[, names(features_keep), drop = FALSE]
    } else if (method %in% c("GSVA", "ssGSEA", "zscore", "PLAGE")) {
      expr_sp <- GetAssayData5(
        srt_sp,
        layer = layer,
        assay = assay
      )
      min_size <- max(minGSSize, min.sz)
      max_size <- min(maxGSSize, max.sz)
      if (identical(backend, "cpp")) {
        if (method == "GSVA") {
          gs_scores <- run_gsva_scores(
            expr_counts = expr_sp,
            gene_sets = features,
            kcdf = kcdf,
            min_gs_size = min_size,
            max_gs_size = max_size,
            max_diff = mx.diff,
            abs_ranking = abs.ranking,
            tau = tau,
            chunk_size = cpp_chunk_size
          )
        } else if (method == "ssGSEA") {
          gs_scores <- run_ssgsea_scores(
            expr_counts = expr_sp,
            gene_sets = features,
            min_gs_size = min_size,
            max_gs_size = max_size,
            alpha = tau,
            normalize = ssgsea.norm
          )
        } else if (method == "zscore") {
          gs_scores <- run_zscore_scores(
            expr_counts = expr_sp,
            gene_sets = features,
            min_gs_size = min_size,
            max_gs_size = max_size
          )
        } else {
          gs_scores <- run_plage_scores(
            expr_counts = expr_sp,
            gene_sets = features,
            min_gs_size = min_size,
            max_gs_size = max_size
          )
        }
        filtered <- names(features)[
          !names(features) %in% colnames(gs_scores)
        ]
      } else {
        gene_set_scoring_require_namespace("GSVA")
        expr_mat <- as_matrix(expr_sp)
        expr_mat <- expr_mat[rowSums(expr_mat) > 0, , drop = FALSE]
        gene_sets_filtered <- lapply(features, function(gs) {
          intersect(gs, rownames(expr_mat))
        })
        keep <- lengths(gene_sets_filtered) >= min_size & lengths(gene_sets_filtered) <= max_size
        gene_sets_filtered <- gene_sets_filtered[keep]
        if (length(gene_sets_filtered) == 0L) {
          log_message(
            "No gene sets remain after intersecting with the expression matrix",
            message_type = "error",
            verbose = verbose
          )
        }
        if (method == "GSVA") {
          param <- GSVA::gsvaParam(
            exprData = expr_mat,
            geneSets = gene_sets_filtered,
            minSize = min_size,
            maxSize = max_size,
            kcdf = kcdf,
            tau = tau,
            maxDiff = mx.diff,
            absRanking = abs.ranking
          )
        } else if (method == "ssGSEA") {
          param <- GSVA::ssgseaParam(
            exprData = expr_mat,
            geneSets = gene_sets_filtered,
            minSize = min_size,
            maxSize = max_size,
            alpha = tau,
            normalize = ssgsea.norm,
            verbose = verbose
          )
        } else if (method == "zscore") {
          param <- GSVA::zscoreParam(
            exprData = expr_mat,
            geneSets = gene_sets_filtered,
            minSize = min_size,
            maxSize = max_size
          )
        } else {
          param <- GSVA::plageParam(
            exprData = expr_mat,
            geneSets = gene_sets_filtered,
            minSize = min_size,
            maxSize = max_size
          )
        }
        gs_scores <- GSVA::gsva(
          param = param,
          verbose = verbose
        )
        if (inherits(gs_scores, "SummarizedExperiment")) {
          gs_scores <- SummarizedExperiment::assay(gs_scores)
        }
        if (!is.matrix(gs_scores)) {
          gs_scores <- as.matrix(gs_scores)
        }
        if (identical(method, "PLAGE")) {
          gs_scores <- orient_plage_scores(
            scores = gs_scores,
            expr = expr_mat,
            gene_sets = gene_sets_filtered
          )
        }
        filtered <- names(features)[
          !names(features) %in% rownames(gs_scores)
        ]
        gs_scores <- Matrix::t(as_matrix(gs_scores))
      }
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
        as_matrix(gs_scores)
      )[, names(features_keep), drop = FALSE]
    } else if (method == "VISION") {
      expr_sp <- GetAssayData5(
        srt_sp,
        layer = layer,
        assay = assay
      )
      vision_scores <- run_vision_scores(
        expr_counts = expr_sp,
        gene_sets = features,
        scale_by_library = identical(layer, "counts")
      )
      filtered <- names(features)[
        !names(features) %in% colnames(vision_scores)
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
        as_matrix(vision_scores)
      )[, names(features_keep), drop = FALSE]
    }
    colnames(scores) <- gene_set_scoring_make_score_colnames(
      feature_names = names(features_nm),
      prefix = name
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
  }
  if (isTRUE(store_metadata)) {
    srt <- Seurat::AddMetaData(
      object = srt,
      metadata = scores_mat
    )
  }

  if (isTRUE(classification)) {
    feature_labels <- features_nm_used[colnames(scores_mat)]
    feature_labels[is.na(feature_labels)] <- colnames(scores_mat)[is.na(feature_labels)]
    assignments <- apply(
      scores_mat,
      MARGIN = 1,
      FUN = function(x) {
        x_valid <- x[!is.na(x)]
        if (length(x_valid) == 0L || all(x_valid < 0)) {
          return(NA)
        } else {
          hits <- which(!is.na(x) & x == max(x_valid))
          if (length(hits) > 1) {
            return(NA)
          } else {
            return(feature_labels[hits])
          }
        }
      }
    )
    class_col <- gene_set_scoring_make_classification_colname(
      prefix = name,
      fallback = method
    )
    srt[[class_col]] <- assignments[rownames(scores_mat)]
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
  seed = 11,
  search = FALSE,
  cores = 1,
  verbose = TRUE,
  ...
) {
  set.seed(seed = seed)
  assay_raw <- SeuratObject::DefaultAssay(object = object)
  assay <- assay %||% assay_raw
  SeuratObject::DefaultAssay(object = object) <- assay
  expr_data <- GetAssayData5(
    object = object,
    layer = layer
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
  ctrl_scores <- do.call(rbind, lapply(scores, `[[`, 1))
  features_scores <- do.call(rbind, lapply(scores, `[[`, 2))

  features_scores_use <- features_scores - ctrl_scores
  if (name == "") {
    rownames(features_scores_use) <- paste0("X", seq_len(cluster_length))
  } else {
    rownames(features_scores_use) <- paste0(name, seq_len(cluster_length))
  }
  features_scores_use <- as.data.frame(t(features_scores_use))
  rownames(features_scores_use) <- colnames(object)
  object[[colnames(features_scores_use)]] <- features_scores_use
  SeuratObject::CheckGC()
  SeuratObject::DefaultAssay(object = object) <- assay_raw

  return(object)
}

check_length <- function(
  values,
  cutoff = 0
) {
  vapply(
    X = values,
    FUN = function(x) {
      length(x) > cutoff
    },
    FUN.VALUE = logical(1)
  )
}
