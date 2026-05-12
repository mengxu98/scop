#' @title Run CytoTRACE 2
#'
#' @description
#' Predicts cellular developmental potential from single-cell RNA-seq data
#' using the CytoTRACE 2 algorithm (Kang et al., 2025).
#' This is a native scop implementation with C++ acceleration.
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
#' @param cores Number of cores for parallel processing.
#' Default is `1`.
#' @param seed Random seed for reproducibility. Default is `14`.
#' @param data_dir Path to the directory containing CytoTRACE2 model data files.
#' If `NULL`, uses data prepared by `PrepareDB()`, the user data cache, or
#' auto-downloads from the datasets repository. Default is `NULL`.
#' @param ... Additional arguments (reserved for future use).
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
  cores = 1,
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
  species <- match.arg(species)
  if (species == "Mus_musculus") {
    species <- "mouse"
  }
  if (species == "Homo_sapiens") {
    species <- "human"
  }

  assay <- assay %||% SeuratObject::DefaultAssay(object = object)

  log_message(
    "Extracting expression matrix from {.code assay = {assay}, layer = {layer}}",
    message_type = "info",
    verbose = verbose
  )

  expression <- SeuratObject::GetAssayData(
    object = object,
    assay = assay,
    layer = layer
  )

  result_df <- run_cytotrace2_impl(
    expression = expression,
    species = species,
    batch_size = batch_size,
    smooth_batch_size = smooth_batch_size,
    cores = cores,
    seed = seed,
    data_dir = data_dir,
    verbose = verbose
  )

  result_df <- result_df[colnames(object), , drop = FALSE]

  object <- Seurat::AddMetaData(object = object, metadata = result_df)

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
  cores = 1,
  seed = 14,
  data_dir = NULL,
  verbose = TRUE,
  ...
) {
  log_message(
    "Running {.pkg CytoTRACE2} (native scop implementation)",
    message_type = "running",
    verbose = verbose
  )

  species <- match.arg(species)
  if (species == "Mus_musculus") {
    species <- "mouse"
  }
  if (species == "Homo_sapiens") {
    species <- "human"
  }

  if (inherits(object, c("matrix", "Matrix"))) {
    object <- as.data.frame(object)
  }

  if (!is.data.frame(object)) {
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

  result_df <- run_cytotrace2_impl(
    expression = object,
    species = species,
    batch_size = batch_size,
    smooth_batch_size = smooth_batch_size,
    cores = cores,
    seed = seed,
    data_dir = data_dir,
    verbose = verbose
  )

  log_message(
    "{.pkg CytoTRACE2} computed successfully",
    message_type = "success",
    verbose = verbose
  )

  return(result_df)
}

run_cytotrace2_impl <- function(
  expression,
  species,
  batch_size,
  smooth_batch_size,
  cores,
  seed,
  data_dir,
  verbose
) {
  set.seed(seed)

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

  log_message(
    "Running on {.val {chunk}} subsample(s)",
    message_type = "info",
    verbose = verbose
  )

  results <- lapply(subsamples, function(subsample) {
    dt <- expression[, subsample, drop = FALSE]

    log_message(
      "Preprocessing subsample ({.val {ncol(dt)}} cells)",
      message_type = "info",
      verbose = verbose
    )
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
    cell_names <- preprocessed$cell_names
    count_cells_few_genes <- preprocessed$count_cells_few_genes

    log_message(
      "Running ensemble prediction and postprocessing",
      message_type = "info",
      verbose = verbose
    )

    smooth_groups <- cytotrace2_smooth_groups(
      n_cells = nrow(log2_data),
      smooth_batch_size = smooth_batch_size,
      seed = seed
    )

    pca_coords <- matrix(0, nrow = nrow(log2_data), ncol = 0)
    if (nrow(log2_data) > 100 && stats::sd(as.vector(log2_data)) != 0) {
      log_message(
        "Computing PCA for kNN smoothing",
        message_type = "info",
        verbose = verbose
      )
      log2_scaled <- scale_rows_custom(as.matrix(log2_data))
      n_pcs <- min(30, nrow(log2_data) - 1)
      svd_res <- RSpectra::svds(
        log2_scaled,
        k = n_pcs,
        opts = list(center = TRUE, scale = FALSE)
      )
      pca_coords <- svd_res$u %*% diag(svd_res$d)
      rownames(pca_coords) <- rownames(log2_data)
    }

    result <- cytotrace2_main(
      rank_data = as.matrix(ranked_data),
      log2_data = as.matrix(log2_data),
      parameter_dict = model_data$parameter_dict,
      smooth_groups = smooth_groups,
      cores = cores,
      seed = seed,
      pca_coords = pca_coords
    )

    result$count_cells_few_genes <- count_cells_few_genes
    result$cell_names <- cell_names
    return(result)
  })

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

  out <- data.frame(
    CytoTRACE2_Score = all_scores,
    CytoTRACE2_Potency = all_potency,
    CytoTRACE2_Relative = (all_scores - min(all_scores)) /
      (max(all_scores) - min(all_scores)),
    preKNN_CytoTRACE2_Score = all_preknn_score,
    preKNN_CytoTRACE2_Potency = all_preknn_potency,
    row.names = all_cell_names,
    stringsAsFactors = FALSE
  )

  out <- out[init_order, , drop = FALSE]

  return(out)
}

cytotrace2_files <- c(
  "model_parameters.rds",
  "features_model_training_17.csv",
  "mt_dict_human_to_mouse.csv",
  "mt_human_alias.csv",
  "mt_mouse_alias.csv"
)

load_cytotrace2_data <- function(data_dir, verbose) {
  data_dir <- resolve_cytotrace2_dir(data_dir, verbose)

  log_message(
    "Loading model from {.path {data_dir}}",
    message_type = "info",
    verbose = verbose
  )

  params_file <- file.path(data_dir, "model_parameters.rds")
  if (!file.exists(params_file)) {
    log_message(
      "Model parameter file not found: {.path {params_file}}",
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
  features <- read.csv(features_file, row.names = 1, check.names = FALSE)[[1]]

  ortho_file <- file.path(data_dir, "mt_dict_human_to_mouse.csv")
  ortho_dict <- NULL
  alias_dict <- NULL
  mouse_alias_dict <- NULL
  if (file.exists(ortho_file)) {
    ortho_dict <- utils::read.csv(
      ortho_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  alias_file <- file.path(data_dir, "mt_human_alias.csv")
  if (file.exists(alias_file)) {
    alias_dict <- read.csv(alias_file, header = TRUE, check.names = FALSE)
  }

  mouse_alias_file <- file.path(data_dir, "mt_mouse_alias.csv")
  if (file.exists(mouse_alias_file)) {
    mouse_alias_dict <- read.csv(
      mouse_alias_file,
      header = TRUE,
      check.names = TRUE
    )
  }

  list(
    parameter_dict = parameter_dict,
    features = features,
    ortho_dict = ortho_dict,
    alias_dict = alias_dict,
    mouse_alias_dict = mouse_alias_dict
  )
}

resolve_cytotrace2_dir <- function(data_dir, verbose) {
  if (!is.null(data_dir)) {
    if (!dir.exists(data_dir)) {
      log_message(
        "Directory {.path {data_dir}} does not exist",
        message_type = "error"
      )
    }
    return(data_dir)
  }

  cache_dir <- file.path(
    tools::R_user_dir("scop", "data"),
    "CytoTRACE2"
  )

  if (
    dir.exists(cache_dir) &&
      all(file.exists(file.path(cache_dir, cytotrace2_files)))
  ) {
    return(cache_dir)
  }

  cyto_cache_key <- list("1.1.0", "CytoTRACE2", "CytoTRACE2")
  if (requireNamespace("R.cache", quietly = TRUE)) {
    cyto_cached <- R.cache::loadCache(key = cyto_cache_key)
    cyto_cached_dir <- if (!is.null(cyto_cached$data_dir)) {
      normalizePath(cyto_cached$data_dir, mustWork = FALSE)
    } else {
      NULL
    }
    cache_dir_norm <- normalizePath(cache_dir, mustWork = FALSE)
    if (
      !is.null(cyto_cached) &&
        identical(cyto_cached_dir, cache_dir_norm) &&
        dir.exists(cyto_cached_dir) &&
        all(file.exists(file.path(
          cyto_cached_dir,
          cytotrace2_files
        )))
    ) {
      log_message(
        "Using CytoTRACE2 data from datasets cache: {.path {cyto_cached_dir}}",
        message_type = "info",
        verbose = verbose
      )
      return(cyto_cached_dir)
    } else if (
      !is.null(cyto_cached_dir) &&
        dir.exists(cyto_cached_dir)
    ) {
      log_message(
        "Ignoring legacy CytoTRACE2 cache outside the datasets cache: {.path {cyto_cached_dir}}",
        message_type = "info",
        verbose = verbose
      )
    }
  }

  log_message(
    "Downloading CytoTRACE2 model data from GitHub repository...",
    message_type = "info",
    verbose = verbose
  )

  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  old_timeout <- getOption("timeout")
  options(timeout = max(600, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  for (fname in cytotrace2_files) {
    url <- paste0(
      "https://raw.githubusercontent.com/mengxu98/datasets/main/CytoTRACE2",
      "/",
      fname
    )
    dest <- file.path(cache_dir, fname)
    if (!file.exists(dest)) {
      log_message(
        "  Downloading {.path {fname}} ...",
        message_type = "info",
        verbose = verbose
      )
      utils::download.file(
        url = url,
        destfile = dest,
        mode = "wb",
        quiet = !verbose
      )
    }
  }

  if (requireNamespace("R.cache", quietly = TRUE)) {
    cyto_cache <- list(
      data_dir = cache_dir,
      files = cytotrace2_files,
      version = "1.1.0"
    )
    R.cache::saveCache(
      cyto_cache,
      key = cyto_cache_key,
      comment = paste0(
        "1.1.0",
        " nterm:",
        length(cytotrace2_files),
        "|CytoTRACE2-CytoTRACE2"
      )
    )
  }

  return(cache_dir)
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
  pop_sd <- function(x, na.rm = TRUE) {
    sqrt(mean((x - mean(x, na.rm = na.rm))^2, na.rm = na.rm))
  }
  row_sds <- apply(x, 1, pop_sd)
  row_sds_handled <- ifelse(row_sds < constant_threshold, 1, row_sds)
  scaled_x <- sweep(x, 1, row_means, "-")
  scaled_x <- sweep(scaled_x, 1, row_sds_handled, "/")

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
  expression <- as.data.frame(data)
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

  expression_mapped <- expression
  rownames(expression_mapped) <- mapping

  common_genes <- intersect(mapping, features)
  expression_mapped <- expression_mapped[common_genes, , drop = FALSE]

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

  num_expressed_per_cell <- colSums(expression_mapped > 0)
  count_cells_few_genes <- sum(num_expressed_per_cell < 500)

  expr_mat <- as.matrix(expression_mapped)
  ranked_data <- t(Rfast::colRanks(
    expr_mat,
    descending = TRUE,
    method = "average"
  ))
  colnames(ranked_data) <- features
  rownames(ranked_data) <- cell_names

  col_sums <- colSums(expr_mat)
  log2_data <- log2((1e6 * t(expr_mat) / col_sums) + 1)
  colnames(log2_data) <- features
  rownames(log2_data) <- cell_names

  list(
    ranked_data = ranked_data,
    log2_data = log2_data,
    cell_names = cell_names,
    count_cells_few_genes = count_cells_few_genes
  )
}
