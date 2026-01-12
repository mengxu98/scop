#' @title Run doublet-calling for single cell RNA-seq data.
#'
#' @md
#' @inheritParams standard_scop
#' @param assay The name of the assay to be used for doublet-calling.
#' Default is `"RNA"`.
#' @param db_method Method used for doublet-calling.
#' Can be one of `"scDblFinder"`, `"Scrublet"`, `"DoubletDetection"`,
#' `"scds_cxds"`, `"scds_bcds"`, `"scds_hybrid"`.
#' @param db_rate The expected doublet rate.
#' Default is calculated as `ncol(srt) / 1000 * 0.01`.
#' @param ... Additional arguments to be passed to the corresponding doublet-calling method.
#'
#' @return Returns a Seurat object with the doublet prediction results and prediction scores stored in the meta.data.
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunDoubletCalling(
#'   pancreas_sub,
#'   db_method = "scDblFinder"
#' )
#' CellDimPlot(
#'   pancreas_sub,
#'   reduction = "umap",
#'   group.by = "db.scDblFinder_class"
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   reduction = "umap",
#'   features = "db.scDblFinder_score"
#' )
RunDoubletCalling <- function(
    srt,
    assay = "RNA",
    db_rate = ncol(srt) / 1000 * 0.01,
    db_method = "scDblFinder",
    ...) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  status <- CheckDataType(srt, layer = "counts", assay = assay)
  if (status != "raw_counts") {
    log_message(
      "Data type is not raw counts",
      message_type = "error"
    )
  }
  db_methods <- c(
    "scDblFinder",
    "Scrublet",
    "DoubletDetection",
    "scds_cxds",
    "scds_bcds",
    "scds_hybrid"
  )
  if (db_method %in% db_methods) {
    methods <- unlist(strsplit(db_method, "_"))
    method1 <- methods[1]
    method2 <- methods[2]
    if (is.na(method2)) {
      args1 <- mget(names(formals()), sys.frame(sys.nframe()))
      args2 <- as.list(match.call())
    } else {
      args1 <- c(
        mget(names(formals()), sys.frame(sys.nframe())),
        method = method2
      )
      args2 <- c(as.list(match.call()), method = method2)
    }
    for (n in names(args2)) {
      args1[[n]] <- args2[[n]]
    }
    args1 <- args1[!names(args1) %in% c("db_method", "...")]
    tryCatch(
      expr = {
        srt <- do.call(
          what = paste0("db_", method1),
          args = args1
        )
      },
      error = function(e) {
        log_message(e, message_type = "error")
      }
    )
    return(srt)
  } else {
    log_message(
      "{.arg db_method} must be one of {.val {db_methods}}",
      message_type = "error"
    )
  }
}

#' @title Run doublet-calling with scDblFinder
#'
#' @md
#' @inheritParams RunDoubletCalling
#' @param ... Additional arguments to be passed to [scDblFinder::scDblFinder()].
#'
#' @export
db_scDblFinder <- function(
    srt,
    assay = "RNA",
    db_rate = ncol(srt) / 1000 * 0.01,
    ...) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  status <- CheckDataType(srt, layer = "counts", assay = assay)
  if (status != "raw_counts") {
    log_message(
      "Data type is not raw counts!",
      message_type = "error"
    )
  }
  check_r("scDblFinder", verbose = FALSE)
  sce <- Seurat::as.SingleCellExperiment(srt, assay = assay)
  sce <- scDblFinder::scDblFinder(sce, dbr = db_rate, verbose = FALSE, ...)
  srt[["db.scDblFinder_score"]] <- sce[["scDblFinder.score"]]
  srt[["db.scDblFinder_class"]] <- sce[["scDblFinder.class"]]
  return(srt)
}

#' @title Run doublet-calling with scds
#'
#' @md
#' @inheritParams RunDoubletCalling
#' @param method The method to be used for doublet-calling.
#' Options are `"hybrid"`, `"cxds"`, or `"bcds"`.
#' @param ... Additional arguments to be passed to [scds::cxds_bcds_hybrid()].
#'
#' @export
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- db_scds(pancreas_sub, method = "hybrid")
#' CellDimPlot(
#'   pancreas_sub,
#'   reduction = "umap",
#'   group.by = "db.scds_hybrid_class"
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   reduction = "umap",
#'   features = "db.scds_hybrid_score"
#' )
db_scds <- function(
    srt,
    assay = "RNA",
    db_rate = ncol(srt) / 1000 * 0.01,
    method = c("hybrid", "cxds", "bcds"),
    ...) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  status <- CheckDataType(srt, layer = "counts", assay = assay)
  if (status != "raw_counts") {
    log_message(
      "Data type is not raw counts!",
      message_type = "error"
    )
  }
  check_r("scds", verbose = FALSE)
  method <- match.arg(method)
  sce <- Seurat::as.SingleCellExperiment(srt, assay = assay)
  sce <- scds::cxds_bcds_hybrid(sce, ...)
  srt[["db.scds_cxds_score"]] <- sce[["cxds_score"]]
  srt[["db.scds_bcds_score"]] <- sce[["bcds_score"]]
  srt[["db.scds_hybrid_score"]] <- sce[["hybrid_score"]]
  ntop <- ceiling(db_rate * ncol(sce))
  db_qc <- names(sort(
    srt[[paste0("db.scds_", method, "_score"), drop = TRUE]],
    decreasing = TRUE
  )[1:ntop])
  srt[[paste0("db.scds_", method, "_class")]] <- "singlet"
  srt[[paste0("db.scds_", method, "_class")]][db_qc, ] <- "doublet"
  return(srt)
}

#' @title Run doublet-calling with Scrublet
#'
#' @md
#' @inheritParams RunDoubletCalling
#' @param ... Additional arguments to be passed to [scrublet.Scrublet](https://github.com/swolock/scrublet).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- db_Scrublet(pancreas_sub)
#' CellDimPlot(
#'   pancreas_sub,
#'   reduction = "umap",
#'   group.by = "db.Scrublet_class"
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   reduction = "umap",
#'   features = "db.Scrublet_score"
#' )
#' }
db_Scrublet <- function(
    srt,
    assay = "RNA",
    db_rate = ncol(srt) / 1000 * 0.01,
    ...) {
  PrepareEnv()
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  status <- CheckDataType(srt, layer = "counts", assay = assay)
  if (status != "raw_counts") {
    log_message(
      "Data type is not raw counts",
      message_type = "error"
    )
  }
  check_python("scrublet", verbose = FALSE)
  scr <- reticulate::import("scrublet")
  raw_counts <- Matrix::t(
    as_matrix(
      GetAssayData5(
        object = srt,
        assay = assay,
        layer = "counts"
      )
    )
  )
  scrub <- scr$Scrublet(raw_counts, expected_doublet_rate = db_rate, ...)
  res <- scrub$scrub_doublets()
  doublet_scores <- res[[1]]
  predicted_doublets <- res[[2]]

  srt[["db.Scrublet_score"]] <- doublet_scores
  srt[["db.Scrublet_class"]] <- sapply(
    predicted_doublets, function(i) {
      switch(as.character(i),
        "FALSE" = "singlet",
        "TRUE" = "doublet"
      )
    }
  )
  log_message(
    "{.pkg Scrublet} doublet calling completed",
    message_type = "success"
  )
  return(srt)
}

#' @title Run doublet-calling with DoubletDetection
#'
#' @md
#' @inheritParams RunDoubletCalling
#' @param ... Additional arguments to be passed to [doubletdetection.BoostClassifier](https://github.com/JonathanShor/DoubletDetection).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- db_DoubletDetection(pancreas_sub)
#' CellDimPlot(
#'   pancreas_sub,
#'   reduction = "umap",
#'   group.by = "db.DoubletDetection_class"
#' )
#'
#' FeatureDimPlot(
#'   pancreas_sub,
#'   reduction = "umap",
#'   features = "db.DoubletDetection_score"
#' )
#' }
db_DoubletDetection <- function(
    srt,
    assay = "RNA",
    db_rate = ncol(srt) / 1000 * 0.01,
    ...) {
  Sys.setenv(NUMBA_NUM_THREADS = "1")
  Sys.setenv(NUMBA_DISABLE_JIT = "0")
  PrepareEnv()

  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  status <- CheckDataType(srt, layer = "counts", assay = assay)
  if (status != "raw_counts") {
    log_message(
      "Data type is not raw counts",
      message_type = "error"
    )
  }
  check_python("doubletdetection")
  doubletdetection <- reticulate::import("doubletdetection")
  counts <- GetAssayData5(
    object = srt,
    assay = assay,
    layer = "counts"
  )
  clf <- doubletdetection$BoostClassifier(
    n_iters = as.integer(5),
    standard_scaling = TRUE,
    ...
  )
  labels <- clf$fit(Matrix::t(counts))$predict()
  scores <- clf$doublet_score()

  labels <- as.integer(reticulate::py_to_r(labels))
  scores <- as.numeric(reticulate::py_to_r(scores))

  n_cells <- ncol(srt)
  if (length(labels) != n_cells) {
    log_message(
      "Length of labels ({.val {length(labels)}}) does not match number of cells",
      message_type = "error"
    )
  }
  if (length(scores) != n_cells) {
    log_message(
      "Length of scores ({.val {length(scores)}}) does not match number of cells",
      message_type = "error"
    )
  }

  cell_names <- colnames(srt)

  if (!is.null(names(scores)) && !identical(names(scores), cell_names)) {
    scores <- scores[cell_names]
  }
  if (!is.null(names(labels)) && !identical(names(labels), cell_names)) {
    labels <- labels[cell_names]
  }

  scores <- unname(scores)
  labels <- unname(labels)

  class_labels <- sapply(
    labels, function(i) {
      switch(as.character(i),
        "0" = "singlet",
        "1" = "doublet"
      )
    }
  )

  # Add metadata directly to meta.data data frame
  srt@meta.data$db.DoubletDetection_score <- scores
  srt@meta.data$db.DoubletDetection_class <- class_labels

  log_message(
    "{.pkg DoubletDetection} doublet calling completed",
    message_type = "success"
  )
  return(srt)
}

#' @title Detect outliers using MAD (Median Absolute Deviation)
#'
#' @md
#' @param x Numeric vector.
#' @param nmads Number of median absolute deviations (MADs) from the median to define the boundaries for outliers.
#' Default is `2.5`.
#' @param constant Constant factor to convert the MAD to a standard deviation.
#' Default is `1.4826`, which is consistent with the MAD of a normal distribution.
#' @param type Type of outliers to detect.
#' Available options are `"both"`, `"lower"`, or `"higher"`.
#' If `type` is `"both"`, it detects both lower and higher outliers.
#' If `type` is `"lower"`, it detects only lower outliers.
#' If `type` is `"higher"`, it detects only higher outliers.
#'
#' @return Numeric vector of indices indicating the positions of outliers in `x`.
#'
#' @export
#'
#' @examples
#' x <- c(1, 2, 3, 4, 5, 100)
#' is_outlier(x) # returns 6
#'
#' x <- c(3, 4, 5, NA, 6, 7)
#' is_outlier(x, nmads = 1.5, type = "lower") # returns 4
#'
#' x <- c(10, 20, NA, 15, 35)
#' is_outlier(x, nmads = 2, type = "higher") # returns 3, 5
is_outlier <- function(
    x,
    nmads = 2.5,
    constant = 1.4826,
    type = c("both", "lower", "higher")) {
  type <- match.arg(type, c("both", "lower", "higher"))
  mad <- stats::mad(x, constant = constant, na.rm = TRUE)
  upper <- stats::median(x, na.rm = TRUE) + nmads * mad
  lower <- stats::median(x, na.rm = TRUE) - nmads * mad
  if (type == "both") {
    out <- which(x > upper | x < lower)
  }
  if (type == "lower") {
    out <- which(x < lower)
  }
  if (type == "higher") {
    out <- which(x > upper)
  }
  out <- c(which(is.na(x)), out)
  return(out)
}

#' @title Run cell-level quality control for single cell RNA-seq data.
#'
#' @md
#' @inheritParams CellDimPlot
#' @inheritParams RunDoubletCalling
#' @inheritParams standard_scop
#' @param return_filtered Logical indicating whether to return a cell-filtered Seurat object.
#' Default is `FALSE`.
#' @param qc_metrics A character vector specifying the quality control metrics to be applied.
#' Default is `c("doublets", "outlier", "umi", "gene", "mito", "ribo", "ribo_mito_ratio", "species")`.
#' @param outlier_threshold A character vector specifying the outlier threshold.
#' Default is `c("log10_nCount:lower:2.5", "log10_nCount:higher:5", "log10_nFeature:lower:2.5", "log10_nFeature:higher:5", "featurecount_dist:lower:2.5")`.
#' See [scuttle::isOutlier].
#' @param db_coefficient The coefficient used to calculate the doublet rate.
#' Default is `0.01`. Doublet rate is calculated as `ncol(srt) / 1000 * db_coefficient`.
#' @param outlier_n Minimum number of outlier metrics that meet the conditions for determining outlier cells.
#' Default is `1`.
#' @param UMI_threshold UMI number threshold. Cells that exceed this threshold will be considered as kept.
#' Default is `3000`.
#' @param gene_threshold Gene number threshold. Cells that exceed this threshold will be considered as kept.
#' Default is `1000`.
#' @param mito_threshold Percentage of UMI counts of mitochondrial genes. Cells that exceed this threshold will be considered as discarded.
#' Default is `20`.
#' @param mito_pattern Regex patterns to match the mitochondrial genes.
#' Default is `c("MT-", "Mt-", "mt-")`.
#' @param mito_gene A defined mitochondrial genes. If features provided, will ignore the `mito_pattern` matching.
#' Default is `NULL`.
#' @param ribo_threshold Percentage of UMI counts of ribosomal genes. Cells that exceed this threshold will be considered as discarded.
#' Default is `50`.
#' @param ribo_pattern Regex patterns to match the ribosomal genes.
#' Default is `c("RP[SL]\\d+\\w{0,1}\\d*$", "Rp[sl]\\d+\\w{0,1}\\d*$", "rp[sl]\\d+\\w{0,1}\\d*$")`.
#' @param ribo_gene A defined ribosomal genes. If features provided, will ignore the `ribo_pattern` matching.
#' Default is `NULL`.
#' @param ribo_mito_ratio_range A numeric vector specifying the range of ribosomal/mitochondrial gene expression ratios for ribo_mito_ratio outlier cells.
#' Default is `c(1, Inf)`.
#' @param species Species used as the suffix of the QC metrics. The first is the species of interest.
#' Default is `NULL`.
#' @param species_gene_prefix Species gene prefix used to calculate QC metrics for each species.
#' Default is `NULL`.
#' @param species_percent Percentage of UMI counts of the first species. Cells that exceed this threshold will be considered as kept.
#' Default is `95`.
#'
#' @return Returns Seurat object with the QC results stored in the meta.data layer.
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunCellQC(pancreas_sub)
#' CellStatPlot(
#'   pancreas_sub,
#'   stat.by = c(
#'     "db_qc", "outlier_qc",
#'     "umi_qc", "gene_qc",
#'     "mito_qc", "ribo_qc",
#'     "ribo_mito_ratio_qc", "species_qc"
#'   ),
#'   plot_type = "upset",
#'   stat_level = "Fail"
#' )
#' table(pancreas_sub$CellQC)
#'
#' data(ifnb_sub)
#' ifnb_sub <- RunCellQC(
#'   srt = ifnb_sub,
#'   split.by = "stim",
#'   UMI_threshold = 1000,
#'   gene_threshold = 550
#' )
#' CellStatPlot(
#'   srt = ifnb_sub,
#'   stat.by = c(
#'     "db_qc", "outlier_qc",
#'     "umi_qc", "gene_qc",
#'     "mito_qc", "ribo_qc",
#'     "ribo_mito_ratio_qc", "species_qc"
#'   ),
#'   plot_type = "upset",
#'   stat_level = "Fail"
#' )
#'
#' table(ifnb_sub$CellQC)
RunCellQC <- function(
    srt,
    assay = "RNA",
    split.by = NULL,
    return_filtered = FALSE,
    qc_metrics = c(
      "doublets",
      "outlier",
      "umi",
      "gene",
      "mito",
      "ribo",
      "ribo_mito_ratio",
      "species"
    ),
    db_method = "scDblFinder",
    db_rate = NULL,
    db_coefficient = 0.01,
    outlier_threshold = c(
      "log10_nCount:lower:2.5",
      "log10_nCount:higher:5",
      "log10_nFeature:lower:2.5",
      "log10_nFeature:higher:5",
      "featurecount_dist:lower:2.5"
    ),
    outlier_n = 1,
    UMI_threshold = 3000,
    gene_threshold = 1000,
    mito_threshold = 20,
    mito_pattern = c("MT-", "Mt-", "mt-"),
    mito_gene = NULL,
    ribo_threshold = 50,
    ribo_pattern = c(
      "RP[SL]\\d+\\w{0,1}\\d*$",
      "Rp[sl]\\d+\\w{0,1}\\d*$",
      "rp[sl]\\d+\\w{0,1}\\d*$"
    ),
    ribo_gene = NULL,
    ribo_mito_ratio_range = c(1, Inf),
    species = NULL,
    species_gene_prefix = NULL,
    species_percent = 95,
    seed = 11) {
  set.seed(seed)

  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  if (isFALSE(assay %in% SeuratObject::Assays(srt))) {
    log_message(
      "srt does not contain '", assay, "' assay.",
      message_type = "error"
    )
  }
  if (length(species) != length(species_gene_prefix)) {
    log_message(
      "'species_gene_prefix' must be the same length as 'species'.",
      message_type = "error"
    )
  }
  if (length(species) == 0) {
    species <- species_gene_prefix <- NULL
  }
  status <- CheckDataType(srt, layer = "counts", assay = assay)
  if (status != "raw_counts") {
    log_message(
      "Data type is not raw counts!",
      message_type = "warning"
    )
  }
  if (!paste0("nCount_", assay) %in% colnames(srt@meta.data)) {
    srt@meta.data[[paste0("nCount_", assay)]] <- Matrix::colSums(
      GetAssayData5(
        srt,
        assay = assay,
        layer = "counts"
      )
    )
  }
  if (!paste0("nFeature_", assay) %in% colnames(srt@meta.data)) {
    srt@meta.data[[paste0("nFeature_", assay)]] <- Matrix::colSums(
      GetAssayData5(
        srt,
        assay = assay,
        layer = "counts"
      ) > 0
    )
  }
  srt_raw <- srt
  if (!is.null(split.by)) {
    srt_list <- Seurat::SplitObject(srt, split.by = split.by)
  } else {
    srt_list <- list(srt)
  }

  for (i in seq_along(srt_list)) {
    srt <- srt_list[[i]]
    if (!is.null(split.by)) {
      log_message("Running QC for ", srt@meta.data[[split.by]][1])
    }
    ntotal <- ncol(srt)

    db_qc <- c()
    if ("doublets" %in% qc_metrics) {
      if (!is.null(db_method)) {
        if (is.null(db_rate)) {
          db_rate <- ncol(srt) / 1000 * db_coefficient
        }
        if (db_rate >= 1) {
          log_message(
            "The db_rate is equal to or greater than 1!",
            message_type = "error"
          )
        }
        for (dbm in db_method) {
          srt <- RunDoubletCalling(
            srt = srt,
            db_method = dbm,
            db_rate = db_rate
          )
          db_qc <- unique(c(
            db_qc,
            colnames(srt)[
              srt[[paste0("db.", dbm, "_class"), drop = TRUE]] == "doublet"
            ]
          ))
        }
      }
    }

    outlier_qc <- c()
    for (n in 1:length(species)) {
      if (n == 0) {
        break
      }
      sp <- species[n]
      prefix <- species_gene_prefix[n]
      sp_genes <- rownames(
        Seurat::GetAssay(
          srt,
          assay = assay
        )
      )[grep(
        pattern = paste0("^", prefix),
        x = rownames(
          Seurat::GetAssay(
            srt,
            assay = assay
          )
        )
      )]
      nCount <- srt[[paste0(
        c(paste0("nCount_", assay), sp),
        collapse = "."
      )]] <- Matrix::colSums(
        GetAssayData5(
          srt,
          assay = assay,
          layer = "counts"
        )[sp_genes, ]
      )
      nFeature <- srt[[paste0(
        c(paste0("nFeature_", assay), sp),
        collapse = "."
      )]] <- Matrix::colSums(
        GetAssayData5(
          srt,
          assay = assay,
          layer = "counts"
        )[sp_genes, ] > 0
      )
      percent.mito <- srt[[paste0(
        c("percent.mito", sp),
        collapse = "."
      )]] <- Seurat::PercentageFeatureSet(
        object = srt,
        assay = assay,
        pattern = paste0(
          "(",
          paste0("^", prefix, "-*", mito_pattern),
          ")",
          collapse = "|"
        ),
        features = mito_gene
      )
      percent.ribo <- srt[[paste0(
        c("percent.ribo", sp),
        collapse = "."
      )]] <- Seurat::PercentageFeatureSet(
        object = srt,
        assay = assay,
        pattern = paste0(
          "(",
          paste0("^", prefix, "-*", ribo_pattern),
          ")",
          collapse = "|"
        ),
        features = ribo_gene
      )
      percent.genome <- srt[[paste0(
        c("percent.genome", sp),
        collapse = "."
      )]] <- Seurat::PercentageFeatureSet(
        object = srt,
        assay = assay,
        pattern = paste0("^", prefix)
      )

      ribo.mito.ratio <- srt[[
        paste0(c("percent.ribo", sp), collapse = "."),
        drop = TRUE
      ]] /
        srt[[paste0(c("percent.mito", sp), collapse = "."), drop = TRUE]]
      ribo.mito.ratio[is.na(ribo.mito.ratio)] <- 1
      srt[[paste0(c("ribo.mito.ratio", sp), collapse = ".")]] <- ribo.mito.ratio

      if (n == 1) {
        if ("outlier" %in% qc_metrics) {
          log10_nFeature <- srt[[paste0(
            c(paste0("log10_nFeature_", assay), sp),
            collapse = "."
          )]] <- log10(nFeature)
          log10_nCount <- srt[[paste0(
            c(paste0("log10_nCount_", assay), sp),
            collapse = "."
          )]] <- log10(nCount)
          log10_nCount[is.infinite(log10_nCount)] <- NA
          log10_nFeature[is.infinite(log10_nFeature)] <- NA
          mod <- stats::loess(log10_nFeature ~ log10_nCount)
          pred <- stats::predict(
            mod,
            newdata = data.frame(log10_nCount = log10_nCount)
          )
          featurecount_dist <- srt[[paste0(
            c("featurecount_dist", sp),
            collapse = "."
          )]] <- log10_nFeature - pred

          var <- sapply(strsplit(outlier_threshold, ":"), function(x) x[[1]])
          var_valid <- var %in%
            colnames(srt@meta.data) |
            sapply(var, FUN = function(x) exists(x, where = environment()))
          if (any(!var_valid)) {
            log_message(
              "Variable ",
              paste0(names(var_valid)[!var_valid], collapse = ","),
              " is not found in the srt object.",
              message_type = "error"
            )
          }
          outlier <- lapply(
            strsplit(outlier_threshold, ":"), function(m) {
              colnames(srt)[is_outlier(
                get(m[1]),
                nmads = as.numeric(m[3]),
                type = m[2]
              )]
            }
          )
          names(outlier) <- outlier_threshold
          outlier_tb <- table(unlist(outlier))
          outlier_qc <- c(
            outlier_qc,
            names(outlier_tb)[outlier_tb >= outlier_n]
          )
          for (nm in names(outlier)) {
            srt[[make.names(nm)]] <- colnames(srt) %in% outlier[[nm]]
          }
        }
      }
    }

    umi_qc <- gene_qc <- mito_qc <- ribo_qc <- ribo_mito_ratio_qc <- species_qc <- c()
    if ("umi" %in% qc_metrics) {
      umi_qc <- colnames(srt)[which(
        srt[[
          paste0(c(paste0("nCount_", assay), species[1]), collapse = "."),
          drop = TRUE
        ]] <
          UMI_threshold
      )]
    }
    if ("gene" %in% qc_metrics) {
      gene_qc <- colnames(srt)[which(
        srt[[
          paste0(c(paste0("nFeature_", assay), species[1]), collapse = "."),
          drop = TRUE
        ]] <
          gene_threshold
      )]
    }
    if ("mito" %in% qc_metrics) {
      mito_qc <- colnames(srt)[which(
        srt[[
          paste0(c("percent.mito", species[1]), collapse = "."),
          drop = TRUE
        ]] >
          mito_threshold
      )]
    }
    if ("ribo" %in% qc_metrics) {
      ribo_qc <- colnames(srt)[which(
        srt[[
          paste0(c("percent.ribo", species[1]), collapse = "."),
          drop = TRUE
        ]] >
          ribo_threshold
      )]
    }
    if ("ribo_mito_ratio" %in% qc_metrics) {
      ribo_mito_ratio_qc <- colnames(srt)[which(
        srt[[
          paste0(c("ribo.mito.ratio", species[1]), collapse = "."),
          drop = TRUE
        ]] <
          ribo_mito_ratio_range[1] |
          srt[[
            paste0(c("ribo.mito.ratio", species[1]), collapse = "."),
            drop = TRUE
          ]] >
            ribo_mito_ratio_range[2]
      )]
    }
    if ("species" %in% qc_metrics) {
      species_qc <- colnames(srt)[which(
        srt[[
          paste0(c("percent.genome", species[1]), collapse = "."),
          drop = TRUE
        ]] <
          species_percent
      )]
    }

    CellQC <- unique(
      c(
        db_qc,
        outlier_qc,
        umi_qc,
        gene_qc,
        mito_qc,
        ribo_qc,
        ribo_mito_ratio_qc,
        species_qc
      )
    )
    log_message(">>> Total cells: ", ntotal)
    log_message(">>> Cells which are filtered out: ", length(CellQC))
    log_message(">>> ", length(db_qc), " potential doublets")
    log_message(">>> ", length(outlier_qc), " outlier cells")
    log_message(">>> ", length(umi_qc), "low-UMI cells")
    log_message(">>> ", length(gene_qc), "low-gene cells")
    log_message(">>> ", length(mito_qc), "high-mito cells")
    log_message(">>> ", length(ribo_qc), "high-ribo cells")
    log_message(">>> ", length(ribo_mito_ratio_qc), "ribo_mito_ratio outlier cells")
    log_message(">>> ", length(species_qc), "species-contaminated cells")
    log_message(">>> Remained cells after filtering: ", ntotal - length(CellQC))

    qc_nm <- c(
      "db_qc",
      "outlier_qc",
      "umi_qc",
      "gene_qc",
      "mito_qc",
      "ribo_qc",
      "ribo_mito_ratio_qc",
      "species_qc",
      "CellQC"
    )
    for (qc in qc_nm) {
      srt[[qc]] <- ifelse(colnames(srt) %in% get(qc), "Fail", "Pass")
      srt[[qc]] <- factor(srt[[qc, drop = TRUE]], levels = c("Pass", "Fail"))
    }

    if (return_filtered) {
      srt <- srt[, srt$CellQC == "Pass"]
      srt@meta.data[, intersect(qc_nm, colnames(srt@meta.data))] <- NULL
    }
    srt_list[[i]] <- srt
  }
  cells <- unlist(lapply(srt_list, colnames))
  srt_raw <- srt_raw[, cells]
  meta.data <- do.call(
    rbind.data.frame,
    unname(lapply(srt_list, function(x) x@meta.data))
  )
  srt_raw <- Seurat::AddMetaData(
    srt_raw,
    metadata = meta.data
  )
  return(srt_raw)
}
