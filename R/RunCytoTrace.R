#' @title Run CytoTRACE 2
#'
#' @description
#' Predicts cellular developmental potential from single-cell RNA-seq data
#' using the CytoTRACE 2 algorithm (Kang et al., 2025).
#' By default, this function uses the native `scop` R/C++ implementation. Set
#' `backend = "r"` to call the official `CytoTRACE2` R package.
#'
#' The algorithm consists of five stages:
#' \enumerate{
#'   \item \strong{Preprocessing}: Gene orthology mapping, feature selection,
#'     ranking, and log2-CPM transformation.
#'   \item \strong{GSBN Ensemble Prediction}: 19 pre-trained Gene Set Binary
#'     Network models predict a continuous developmental potency score (0-1)
#'     and a discrete potency category.
#'   \item \strong{Diffusion Smoothing}: A Markov random-walk-with-restart on
#'     a cell-cell similarity graph smooths the raw scores.
#'   \item \strong{Binning}: Within each potency category, cells are ranked
#'     and linearly scaled to corresponding segments of the unit interval.
#'   \item \strong{Adaptive kNN Smoothing}: PCA-based nearest-neighbor
#'     consensus refinement of the final scores.
#' }
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams standard_scop
#' @inheritParams GroupHeatmap
#' @param object An object.
#' This can be a Seurat object or a matrix-like object (genes as rows, cells as columns).
#' @param species The species of the input data.
#' Currently supported values are `"Homo_sapiens"` and `"Mus_musculus"`.
#' Default is `"Homo_sapiens"`.
#' @param batch_size The number of cells to process at once.
#' For datasets with more cells than this value, cells are randomly split
#' into batches and processed independently.
#' No batching if `NULL`.
#' Default is `10000`.
#' @param smooth_batch_size The number of cells per subsample for the
#' diffusion smoothing step.
#' No diffusion subsampling if `NULL`.
#' Default is `1000`.
#' @param compute_knn_smoothing Whether to run the final PCA-based adaptive
#' kNN smoothing step. Set to `FALSE` for a faster score using the pre-kNN
#' binned CytoTRACE2 output.
#' @param cores Number of cores for parallel processing.
#' Default is `1`.
#' @param backend Backend used to run CytoTRACE2. `"cpp"` uses the native
#' `scop` R/C++ backend and is the default. `"r"` calls the official
#' `CytoTRACE2::cytotrace2()` implementation.
#' @param seed Random seed for reproducibility. Default is `14`.
#' @param data_dir Path to the directory containing CytoTRACE2 model data files.
#' Used only by `backend = "cpp"`. If `NULL`, uses model data prepared by
#' `PrepareDB(db = "CytoTRACE2")`. Default is `NULL`.
#' @param ... Additional arguments passed to the official
#' `CytoTRACE2::cytotrace2()` call when `backend = "r"`.
#'
#' @rdname RunCytoTRACE
#' @export
#'
#' @return
#' When the input is a Seurat object,
#' the function returns a Seurat object with the following metadata columns added:
#' \itemize{
#'   \item \code{CytoTRACE2_Score}: The final predicted cellular potency score (0-1)
#'   \item \code{CytoTRACE2_Potency}: The final predicted cellular potency category
#'     (Differentiated, Unipotent, Oligopotent, Multipotent, Pluripotent, Totipotent)
#'   \item \code{CytoTRACE2_Relative}: The predicted relative order (normalized to 0-1)
#'   \item \code{preKNN_CytoTRACE2_Score}: The potency score before KNN smoothing
#'   \item \code{preKNN_CytoTRACE2_Potency}: The potency category before KNN smoothing
#' }
#'
#' When the input is a matrix or data.frame,
#' the function returns a data.frame with the same columns as above,
#' with cell IDs as row names.
#'
#' @references
#' Kang, M., Brown, E., Almagro Armenteros, J.J. et al.
#' "Improved reconstruction of single-cell developmental potential
#' with CytoTRACE 2." \emph{Nature Methods} (2025).
#' \doi{10.1038/s41592-025-02857-2}
#'
#' Model data: \url{https://github.com/mengxu98/datasets/tree/main/CytoTRACE2}
#'
#' @section License:
#' The CytoTRACE 2 model and associated data files are provided under the
#' Stanford Non-Commercial Software License Agreement.
#' Commercial entities wishing to use this software should contact
#' Stanford University's Office of Technology Licensing (docket S24-057).
#' See \url{https://github.com/mengxu98/datasets/blob/main/CytoTRACE2/LICENSE}
#' for complete terms.
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunCytoTRACE(
#'   pancreas_sub,
#'   species = "Mus_musculus"
#' )
#'
#' CytoTRACEPlot(
#'   pancreas_sub,
#'   xlab = "UMAP_1",
#'   ylab = "UMAP_2",
#'   ncol = 2
#' )
RunCytoTRACE <- function(object, ...) {
  UseMethod(generic = "RunCytoTRACE", object = object)
}

#' @rdname RunCytoTRACE
#' @method RunCytoTRACE Seurat
#' @export
RunCytoTRACE.Seurat <- function(
  object,
  assay = NULL,
  layer = c("counts", "data"),
  species = c("Homo_sapiens", "Mus_musculus"),
  batch_size = 10000,
  smooth_batch_size = 1000,
  compute_knn_smoothing = TRUE,
  cores = 1,
  backend = c("cpp", "r"),
  seed = 14,
  data_dir = NULL,
  verbose = TRUE,
  ...
) {
  log_message(
    "Running {.pkg CytoTRACE2}",
    message_type = "running",
    verbose = verbose
  )

  layer <- match.arg(layer)
  backend <- match.arg(backend)
  species <- match.arg(species)
  cores <- max(1L, as.integer(cores))

  assay <- assay %||% SeuratObject::DefaultAssay(object = object)

  log_message(
    "Extracting expression matrix from {.code assay = {assay}, layer = {layer}}",
    message_type = "info",
    verbose = verbose
  )

  if (identical(backend, "r")) {
    check_r("digitalcytometry/cytotrace2/cytotrace2_r", verbose = FALSE)
    if (!is.null(data_dir)) {
      log_message(
        "{.arg data_dir} is ignored by {.arg backend = 'r'} because the official {.pkg CytoTRACE2} package uses its installed model data.",
        message_type = "warning",
        verbose = verbose
      )
    }

    if (species == "Mus_musculus") {
      species <- "mouse"
    }
    if (species == "Homo_sapiens") {
      species <- "human"
    }

    cytotrace2_fun <- get_namespace_fun("CytoTRACE2", "cytotrace2")
    args <- utils::modifyList(
      list(
        input = object,
        species = species,
        is_seurat = identical(assay, "RNA"),
        slot_type = layer,
        batch_size = batch_size,
        smooth_batch_size = smooth_batch_size,
        ncores = cores,
        seed = seed
      ),
      list(...)
    )
    if (!identical(assay, "RNA")) {
      args$input <- SeuratObject::GetAssayData(
        object = object,
        assay = assay,
        layer = layer
      )
      args$is_seurat <- FALSE
    }
    args <- args[names(args) %in% names(formals(cytotrace2_fun))]
    result <- do.call(cytotrace2_fun, args)
    if (inherits(result, "Seurat")) {
      object <- result
    } else {
      result <- result[colnames(object), , drop = FALSE]
      object <- Seurat::AddMetaData(object = object, metadata = result)
    }
  } else {
    expression <- SeuratObject::GetAssayData(
      object = object,
      assay = assay,
      layer = layer
    )

    result_df <- RunCytoTRACE.default(
      object = expression,
      species = species,
      batch_size = batch_size,
      smooth_batch_size = smooth_batch_size,
      compute_knn_smoothing = compute_knn_smoothing,
      cores = cores,
      backend = "cpp",
      seed = seed,
      data_dir = data_dir,
      verbose = verbose
    )

    result_df <- result_df[colnames(object), , drop = FALSE]

    object <- Seurat::AddMetaData(object = object, metadata = result_df)
  }

  object <- Seurat::LogSeuratCommand(object = object)

  log_message(
    "{.pkg CytoTRACE2} computed successfully",
    message_type = "success",
    verbose = verbose
  )

  return(object)
}

#' @rdname RunCytoTRACE
#' @method RunCytoTRACE default
#' @export
RunCytoTRACE.default <- function(
  object,
  species = c("Homo_sapiens", "Mus_musculus"),
  batch_size = 10000,
  smooth_batch_size = 1000,
  compute_knn_smoothing = TRUE,
  cores = 1,
  backend = c("cpp", "r"),
  seed = 14,
  data_dir = NULL,
  verbose = TRUE,
  ...
) {
  species <- match.arg(species)
  backend <- match.arg(backend)
  cores <- max(1L, as.integer(cores))
  log_message(
    "Running {.pkg CytoTRACE2} with {.arg backend = {backend}}",
    message_type = "running",
    verbose = verbose
  )
  if (species == "Mus_musculus") {
    species <- "mouse"
  }
  if (species == "Homo_sapiens") {
    species <- "human"
  }

  if (!inherits(object, c("matrix", "Matrix", "data.frame"))) {
    log_message(
      "{.arg object} must be a {.cls Seurat}, matrix, or data.frame",
      message_type = "error"
    )
  }

  if (is.null(rownames(object))) {
    log_message(
      "{.arg object} must have row names (gene names)",
      message_type = "error"
    )
  }
  if (is.null(colnames(object))) {
    log_message(
      "{.arg object} must have column names (cell IDs)",
      message_type = "error"
    )
  }

  if (identical(backend, "r")) {
    check_r("digitalcytometry/cytotrace2/cytotrace2_r", verbose = FALSE)
    if (!is.null(data_dir)) {
      log_message(
        "{.arg data_dir} is ignored by {.arg backend = 'r'} because the official {.pkg CytoTRACE2} package uses its installed model data.",
        message_type = "warning",
        verbose = verbose
      )
    }
    cytotrace2_fun <- get_namespace_fun("CytoTRACE2", "cytotrace2")
    args <- utils::modifyList(
      list(
        input = object,
        species = species,
        is_seurat = FALSE,
        slot_type = "counts",
        batch_size = batch_size,
        smooth_batch_size = smooth_batch_size,
        ncores = cores,
        seed = seed
      ),
      list(...)
    )
    args <- args[names(args) %in% names(formals(cytotrace2_fun))]
    result_df <- do.call(cytotrace2_fun, args)
  } else {
    set.seed(seed)

    expression <- object
    model_data <- load_cytotrace2_data(data_dir, verbose)

    if (any(duplicated(rownames(expression)))) {
      log_message("Gene names must be unique", message_type = "error")
    }
    if (any(duplicated(colnames(expression)))) {
      log_message("Cell names must be unique", message_type = "error")
    }
    if (max(expression, na.rm = TRUE) < 20) {
      log_message(
        "Data may already be log-transformed (max < 20). Please provide raw counts or CPM/TPM normalized counts.",
        message_type = "warning",
        verbose = verbose
      )
    }

    log_message(
      "Dataset contains {.val {nrow(expression)}} genes and {.val {ncol(expression)}} cells.",
      message_type = "info",
      verbose = verbose
    )

    init_order <- colnames(expression)

    if (is.null(batch_size)) {
      batch_size <- ncol(expression)
    }
    if (batch_size > ncol(expression)) {
      batch_size <- ncol(expression)
    }

    chunk <- ceiling(ncol(expression) / batch_size)
    if (chunk > 1) {
      subsamples <- split(
        seq_len(ncol(expression)),
        sample(factor(seq_len(ncol(expression)) %% chunk))
      )
    } else {
      subsamples <- list(seq_len(ncol(expression)))
    }
    names(subsamples) <- paste0(
      "subsample ",
      seq_along(subsamples),
      "/",
      length(subsamples)
    )

    log_message(
      "Running on {.val {chunk}} subsample(s)",
      message_type = "info",
      verbose = verbose
    )

    outer_cores <- min(cores, length(subsamples))
    inner_cores <- max(1L, floor(cores / outer_cores))

    process_subsample <- function(subsample) {
      dt <- expression[, subsample, drop = FALSE]

      preprocessed <- cytotrace2_preprocess(
        data = dt,
        species = species,
        features = model_data$features,
        ortho_dict = model_data$ortho_dict,
        alias_dict = model_data$alias_dict,
        mouse_alias_dict = model_data$mouse_alias_dict,
        verbose = verbose
      )

      ranked_data <- preprocessed$ranked_data
      log2_data <- preprocessed$log2_data
      count_cells_few_genes <- preprocessed$count_cells_few_genes
      cell_names <- preprocessed$cell_names
      rm(dt, preprocessed)
      gc(verbose = FALSE)

      smooth_groups <- cytotrace2_smooth_groups(
        n_cells = nrow(log2_data),
        smooth_batch_size = smooth_batch_size,
        seed = seed
      )

      pca_coords <- matrix(0, nrow = nrow(log2_data), ncol = 0)
      if (
        isTRUE(compute_knn_smoothing) &&
          nrow(log2_data) > 100 &&
          stats::sd(as.vector(log2_data)) != 0
      ) {
        log2_scaled <- scale_rows_custom(as.matrix(log2_data))
        n_pcs <- min(30, nrow(log2_data) - 1)
        svd_res <- RSpectra::svds(
          log2_scaled,
          k = n_pcs,
          opts = list(center = TRUE, scale = FALSE)
        )
        pca_coords <- sweep(svd_res$u, 2, svd_res$d, "*")
        rownames(pca_coords) <- rownames(log2_data)
        rm(log2_scaled, svd_res)
        gc(verbose = FALSE)
      }

      result <- cytotrace2_main(
        rank_data = as.matrix(ranked_data),
        log2_data = as.matrix(log2_data),
        parameter_dict = model_data$parameter_dict,
        smooth_groups = smooth_groups,
        cores = inner_cores,
        seed = seed,
        pca_coords = pca_coords
      )

      result$count_cells_few_genes <- count_cells_few_genes
      result$cell_names <- cell_names
      result
    }

    results <- thisutils::parallelize_fun(
      x = subsamples,
      fun = process_subsample,
      cores = outer_cores,
      verbose = verbose
    )

    all_cell_names <- unlist(lapply(results, `[[`, "cell_names"))
    all_scores <- do.call(c, lapply(results, `[[`, "CytoTRACE2_Score"))
    all_potency <- do.call(c, lapply(results, `[[`, "CytoTRACE2_Potency"))
    all_preknn_score <- do.call(
      c,
      lapply(results, `[[`, "preKNN_CytoTRACE2_Score")
    )
    all_preknn_potency <- do.call(
      c,
      lapply(results, `[[`, "preKNN_CytoTRACE2_Potency")
    )

    all_count_cells_few_genes <- sum(unlist(lapply(
      results,
      `[[`,
      "count_cells_few_genes"
    )))

    frac_cells_few_genes <- all_count_cells_few_genes / ncol(expression)
    if (frac_cells_few_genes >= 0.2) {
      log_message(
        sprintf(
          "%.2f%% of input cells express fewer than %d genes. For best results, a minimum gene count of 500-1000 is recommended.",
          frac_cells_few_genes * 100,
          500
        ),
        message_type = "warning",
        verbose = verbose
      )
    }

    result_df <- data.frame(
      CytoTRACE2_Score = all_scores,
      CytoTRACE2_Potency = all_potency,
      CytoTRACE2_Relative = (all_scores - min(all_scores)) /
        (max(all_scores) - min(all_scores)),
      preKNN_CytoTRACE2_Score = all_preknn_score,
      preKNN_CytoTRACE2_Potency = all_preknn_potency,
      row.names = all_cell_names,
      stringsAsFactors = FALSE
    )

    result_df <- result_df[init_order, , drop = FALSE]
  }

  log_message(
    "{.pkg CytoTRACE2} computed successfully",
    message_type = "success",
    verbose = verbose
  )

  return(result_df)
}

.cytotrace2_data_cache <- new.env(parent = emptyenv())

load_cytotrace2_data <- function(data_dir, verbose) {
  if (!is.null(data_dir)) {
    if (!dir.exists(data_dir)) {
      log_message(
        "Directory {.path {data_dir}} does not exist",
        message_type = "error"
      )
    }
  } else {
    cyto_cache <- PrepareDB(db = "CytoTRACE2", verbose = verbose)[["CytoTRACE2"]]
    if (
      is.null(cyto_cache) ||
        is.null(cyto_cache$data_dir) ||
        !dir.exists(cyto_cache$data_dir)
    ) {
      log_message(
        "Prepared CytoTRACE2 model data were not found. Run {.code PrepareDB(db = 'CytoTRACE2')} first or provide {.arg data_dir}.",
        message_type = "error"
      )
    }
    data_dir <- cyto_cache$data_dir
  }
  cache_key <- normalizePath(data_dir, mustWork = FALSE)
  if (exists(cache_key, envir = .cytotrace2_data_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .cytotrace2_data_cache, inherits = FALSE))
  }

  log_message(
    "Loading model from {.path {data_dir}}",
    message_type = "info",
    verbose = verbose
  )

  params_candidates <- file.path(
    data_dir,
    c("parameter_dict_19.rds", "model_parameters.rds")
  )
  params_file <- params_candidates[file.exists(params_candidates)][1]
  if (is.na(params_file)) {
    log_message(
      "Model parameter file not found in {.path {data_dir}}",
      message_type = "error"
    )
  }
  parameter_dict <- readRDS(params_file)

  features_file <- file.path(data_dir, "features_model_training_17.csv")
  if (!file.exists(features_file)) {
    log_message(
      "Feature file not found: {.path {features_file}}",
      message_type = "error"
    )
  }
  features_tbl <- utils::read.csv(
    features_file,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (ncol(features_tbl) < 2) {
    log_message(
      "CytoTRACE2 feature file is malformed: {.path {features_file}}",
      message_type = "error"
    )
  }
  features <- as.character(features_tbl[[2]])
  features <- features[!is.na(features) & nzchar(features)]
  if (length(features) == 0) {
    log_message(
      "CytoTRACE2 feature file contains no feature names: {.path {features_file}}",
      message_type = "error"
    )
  }
  if (anyDuplicated(features) > 0) {
    log_message(
      "CytoTRACE2 feature file contains duplicated feature names. Re-run {.code PrepareDB(db = 'CytoTRACE2', db_update = TRUE)} to refresh the model data.",
      message_type = "error"
    )
  }

  ortho_file <- file.path(data_dir, "mt_dict_human_to_mouse.csv")
  ortho_dict <- NULL
  alias_dict <- NULL
  mouse_alias_dict <- NULL
  if (file.exists(ortho_file)) {
    ortho_dict <- suppressWarnings(utils::read.csv(
      ortho_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  alias_file <- file.path(data_dir, "mt_human_alias.csv")
  if (file.exists(alias_file)) {
    alias_dict <- suppressWarnings(read.csv(
      alias_file,
      header = TRUE,
      check.names = FALSE
    ))
  }

  mouse_alias_file <- file.path(data_dir, "mt_mouse_alias.csv")
  if (file.exists(mouse_alias_file)) {
    mouse_alias_dict <- suppressWarnings(read.csv(
      mouse_alias_file,
      header = TRUE,
      check.names = TRUE
    ))
  }

  model_data <- list(
    parameter_dict = parameter_dict,
    features = features,
    ortho_dict = ortho_dict,
    alias_dict = alias_dict,
    mouse_alias_dict = mouse_alias_dict
  )
  assign(cache_key, model_data, envir = .cytotrace2_data_cache)
  model_data
}

cytotrace2_smooth_groups <- function(n_cells, smooth_batch_size, seed) {
  if (is.null(smooth_batch_size)) {
    smooth_batch_size <- n_cells
  }
  if (smooth_batch_size > n_cells) {
    smooth_batch_size <- n_cells
  }

  chunk_number <- ceiling(n_cells / smooth_batch_size)
  if (chunk_number < 1) {
    chunk_number <- 1
  }

  set.seed(seed)
  as.integer(as.character(sample(factor(seq_len(n_cells) %% chunk_number))))
}

scale_rows_custom <- function(
  x,
  constant_threshold = 10 * .Machine$double.eps,
  mean_tolerance = 1e-10
) {
  row_means <- rowMeans(x, na.rm = TRUE)
  row_sds <- sqrt(rowMeans((x - row_means)^2, na.rm = TRUE))
  row_sds_handled <- ifelse(row_sds < constant_threshold, 1, row_sds)
  scaled_x <- (x - row_means) / row_sds_handled

  residual_means <- rowMeans(scaled_x, na.rm = TRUE)
  rows_to_adjust <- which(abs(residual_means) > mean_tolerance)
  if (length(rows_to_adjust) > 0) {
    scaled_x[rows_to_adjust, ] <- sweep(
      scaled_x[rows_to_adjust, , drop = FALSE],
      1,
      residual_means[rows_to_adjust],
      "-"
    )
  }
  return(scaled_x)
}

cytotrace2_preprocess <- function(
  data,
  species,
  features,
  ortho_dict,
  alias_dict,
  mouse_alias_dict,
  verbose
) {
  gene_names <- rownames(data)
  expression <- as.matrix(data)
  cell_names <- colnames(expression)

  if (species == "human") {
    if (is.null(ortho_dict)) {
      log_message(
        "Orthology dictionary required for human data",
        message_type = "error"
      )
    }

    idx <- match(gene_names, colnames(ortho_dict))
    mapping <- gene_names
    mapped_idx <- !is.na(idx)
    mapping[mapped_idx] <- as.character(ortho_dict[1, idx[mapped_idx]])

    if (!is.null(alias_dict)) {
      map_df <- data.frame(
        original_gene = gene_names,
        mapped_gene = mapping,
        mapped = ifelse(mapping %in% ortho_dict[1, ], "1", "0"),
        stringsAsFactors = FALSE
      )
      mapped_genes <- map_df$original_gene[map_df$mapped == "1"]
      unmapped_genes <- map_df$original_gene[map_df$mapped == "0"]
      is_alias <- intersect(unmapped_genes, alias_dict$alias)
      original_to_alias <- alias_dict$hsgene[
        alias_dict$alias %in% is_alias
      ]
      mapped_original_to_alias <- intersect(
        original_to_alias,
        mapped_genes
      )
      alias_to_unmapped <- alias_dict$alias[
        alias_dict$hsgene %in%
          setdiff(
            original_to_alias,
            mapped_original_to_alias
          )
      ]
      rownames(alias_dict) <- alias_dict$alias
      unmapped_with_alias <- map_df$mapped_gene[
        map_df$original_gene %in% alias_to_unmapped
      ]
      map_df$mapped_gene[
        map_df$original_gene %in% alias_to_unmapped
      ] <- alias_dict[unmapped_with_alias, "mmgene"]
      mapping <- map_df$mapped_gene
    }
  } else {
    mapping <- gene_names
    if (!is.null(mouse_alias_dict)) {
      unmapped_genes <- setdiff(mapping, features)
      unmapped_with_alias <- which(
        mapping %in% unmapped_genes[unmapped_genes %in% mouse_alias_dict$alias]
      )
      filtered <- mouse_alias_dict[
        (mouse_alias_dict$alias %in% mapping[unmapped_with_alias]) &
          !(mouse_alias_dict$mmgene %in% mapping),
      ]
      rownames(filtered) <- filtered$alias
      mapping[mapping %in% filtered$alias] <- filtered[
        mapping[mapping %in% filtered$alias],
        "mmgene"
      ]
    }
  }

  n_mapped <- length(base::intersect(mapping, features))
  log_message(
    "{.val {n_mapped}} input genes mapped to model genes.",
    message_type = "info",
    verbose = verbose
  )

  if (n_mapped < 9000) {
    log_message(
      "The number of input genes mapped to the model is low. Verify the input species is correct. Model performance might be compromised.",
      message_type = "warning",
      verbose = verbose
    )
  }

  rownames(expression) <- mapping

  common_genes <- intersect(mapping, features)
  expression_mapped <- expression[common_genes, , drop = FALSE]

  missing_features <- setdiff(features, common_genes)
  if (length(missing_features) > 0) {
    zero_mat <- matrix(
      0,
      nrow = length(missing_features),
      ncol = ncol(expression_mapped)
    )
    rownames(zero_mat) <- missing_features
    colnames(zero_mat) <- colnames(expression_mapped)
    expression_mapped <- rbind(expression_mapped, zero_mat)
  }

  expression_mapped <- expression_mapped[features, , drop = FALSE]

  preprocessed_numeric <- cytotrace2_preprocess_numeric(expression_mapped)
  ranked_data <- preprocessed_numeric$ranked_data
  log2_data <- preprocessed_numeric$log2_data
  colnames(ranked_data) <- features
  rownames(ranked_data) <- cell_names
  colnames(log2_data) <- features
  rownames(log2_data) <- cell_names

  list(
    ranked_data = ranked_data,
    log2_data = log2_data,
    cell_names = cell_names,
    count_cells_few_genes = preprocessed_numeric$count_cells_few_genes
  )
}
