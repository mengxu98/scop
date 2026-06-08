#' @title Run Giotto cell proximity enrichment
#'
#' @description
#' Run Giotto spatial cell proximity enrichment as a temporary backend and
#' store the pairwise enrichment table in `srt@tools`.
#'
#' @md
#' @inheritParams RunGiottoCluster
#' @param group.by Seurat metadata column used as the cell or spot group.
#' @param network_method Spatial network method passed to
#' `Giotto::createSpatialNetwork()`.
#' @param number_of_simulations Number of random simulations.
#' @param adjust_method Multiple testing adjustment method.
#' @param enrichment_params Additional parameters passed to
#' `Giotto::cellProximityEnrichment()`.
#'
#' @return A `Seurat` object with Giotto proximity results in
#' `srt@tools[[tool_name]]`.
#'
#' @export
RunGiottoCellProximity <- function(
  srt,
  group.by,
  assay = NULL,
  layer = "data",
  image = NULL,
  coord.cols = c("x", "y"),
  network_method = "Delaunay",
  number_of_simulations = 1000,
  adjust_method = "fdr",
  tool_name = "GiottoCellProximity",
  store_giotto = FALSE,
  conversion_params = list(),
  network_params = list(),
  enrichment_params = list(),
  verbose = TRUE,
  seed = 11
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  giotto_validate_scalar_string(group.by, "group.by")
  giotto_validate_scalar_string(tool_name, "tool_name")
  giotto_validate_named_list(conversion_params, "conversion_params")
  giotto_validate_named_list(network_params, "network_params")
  giotto_validate_named_list(enrichment_params, "enrichment_params")
  if (!group.by %in% colnames(srt@meta.data)) {
    log_message("{.arg group.by} {.val {group.by}} is not present in {.cls Seurat} metadata", message_type = "error")
  }
  if (!is.numeric(number_of_simulations) || length(number_of_simulations) != 1L || is.na(number_of_simulations) || number_of_simulations < 1) {
    log_message("{.arg number_of_simulations} must be a positive number", message_type = "error")
  }
  network_method <- match.arg(network_method, c("Delaunay", "kNN"))
  adjust_method <- match.arg(adjust_method, c("none", "fdr", "bonferroni", "BH", "holm", "hochberg", "hommel", "BY"))

  giotto_require(verbose = verbose)
  giotto_old_options <- giotto_disable_r_only_options()
  on.exit(options(giotto_old_options), add = TRUE)

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  input <- giotto_prepare_input(
    srt = srt,
    assay = assay,
    layer = layer,
    features = NULL,
    image = image,
    coord.cols = coord.cols
  )
  gobject <- giotto_create_object(
    input = input,
    layer = layer,
    conversion_params = conversion_params,
    verbose = verbose
  )
  network_name <- giotto_spatial_network_name(network_method, network_params)
  gobject <- giotto_create_spatial_network(
    gobject = gobject,
    network_name = network_name,
    network_method = network_method,
    network_params = network_params,
    verbose = verbose
  )

  cellProximityEnrichment <- giotto_get_fun("cellProximityEnrichment")
  enrichment <- giotto_call(
    cellProximityEnrichment,
    giotto_merge_args(
      list(
        gobject = gobject,
        spatial_network_name = network_name,
        cluster_column = group.by,
        number_of_simulations = as.integer(number_of_simulations),
        adjust_method = adjust_method,
        set_seed = TRUE,
        seed_number = seed
      ),
      enrichment_params
    )
  )
  enrichment_table <- giotto_result_to_data_frame(enrichment)

  tool <- list(
    enrichment = enrichment_table,
    raw_result = enrichment,
    parameters = list(
      assay = assay,
      layer = layer,
      group.by = group.by,
      image = image,
      coord.cols = coord.cols,
      network_method = network_method,
      network_name = network_name,
      number_of_simulations = as.integer(number_of_simulations),
      adjust_method = adjust_method,
      conversion_params = conversion_params,
      network_params = network_params,
      enrichment_params = enrichment_params,
      seed = seed
    ),
    cells = input$cells,
    features = input$features
  )
  if (isTRUE(store_giotto)) {
    tool$giotto <- gobject
  }
  srt@tools[[tool_name]] <- tool
  log_message("{.pkg Giotto} cell proximity enrichment stored in {.code srt@tools[[tool_name]]}", message_type = "success", verbose = verbose)
  srt
}

#' @title Run Giotto spatial gene detection
#'
#' @description
#' Run Giotto `binSpect()` spatial gene detection and store the ranked result
#' table and top features in a `Seurat` object.
#'
#' @md
#' @inheritParams RunGiottoCellProximity
#' @param features Features to test. If `NULL`, variable features are used,
#' falling back to all assay features.
#' @param bin_method Binning method passed to `Giotto::binSpect()`.
#' @param set_variable_features Whether to update `VariableFeatures()`.
#' @param top_n Number of top spatial genes to store in `srt@misc`.
#' @param binSpect_params Additional parameters passed to `Giotto::binSpect()`.
#'
#' @return A `Seurat` object with spatial gene results in `srt@tools` and top
#' genes in `srt@misc[[tool_name]]`.
#'
#' @export
RunGiottoSpatialGenes <- function(
  srt,
  assay = NULL,
  layer = "data",
  features = NULL,
  image = NULL,
  coord.cols = c("x", "y"),
  network_method = "Delaunay",
  bin_method = c("kmeans", "rank"),
  set_variable_features = TRUE,
  top_n = 100,
  tool_name = "GiottoSpatialGenes",
  store_giotto = FALSE,
  conversion_params = list(),
  network_params = list(),
  binSpect_params = list(),
  verbose = TRUE,
  seed = 11
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  giotto_validate_scalar_string(tool_name, "tool_name")
  giotto_validate_named_list(conversion_params, "conversion_params")
  giotto_validate_named_list(network_params, "network_params")
  giotto_validate_named_list(binSpect_params, "binSpect_params")
  if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 1) {
    log_message("{.arg top_n} must be a positive number", message_type = "error")
  }
  network_method <- match.arg(network_method, c("Delaunay", "kNN"))
  bin_method <- match.arg(bin_method)

  giotto_require(verbose = verbose)
  giotto_old_options <- giotto_disable_r_only_options()
  on.exit(options(giotto_old_options), add = TRUE)

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  input <- giotto_prepare_input(
    srt = srt,
    assay = assay,
    layer = layer,
    features = features,
    image = image,
    coord.cols = coord.cols
  )
  gobject <- giotto_create_object(
    input = input,
    layer = layer,
    conversion_params = conversion_params,
    verbose = verbose
  )
  expression_values <- giotto_expression_values(layer)
  if (identical(expression_values, "raw")) {
    normalizeGiotto <- giotto_get_fun("normalizeGiotto")
    gobject <- giotto_call(normalizeGiotto, list(gobject = gobject, verbose = FALSE))
    expression_values <- "normalized"
  }
  network_name <- giotto_spatial_network_name(network_method, network_params)
  gobject <- giotto_create_spatial_network(
    gobject = gobject,
    network_name = network_name,
    network_method = network_method,
    network_params = network_params,
    verbose = verbose
  )

  binSpect <- giotto_get_fun("binSpect")
  result <- giotto_call(
    binSpect,
    giotto_merge_args(
      list(
        gobject = gobject,
        bin_method = bin_method,
        expression_values = expression_values,
        subset_feats = input$features,
        spatial_network_name = network_name,
        do_parallel = FALSE,
        cores = 1,
        seed = seed,
        verbose = verbose
      ),
      binSpect_params
    )
  )
  result_table <- giotto_result_to_data_frame(result)
  top_features <- giotto_top_features(result_table, n = as.integer(top_n))
  srt@misc[[tool_name]] <- top_features
  if (isTRUE(set_variable_features) && length(top_features) > 0L) {
    SeuratObject::VariableFeatures(srt, assay = assay) <- top_features
  }

  tool <- list(
    results = result_table,
    top_features = top_features,
    raw_result = result,
    parameters = list(
      assay = assay,
      layer = layer,
      image = image,
      coord.cols = coord.cols,
      network_method = network_method,
      network_name = network_name,
      bin_method = bin_method,
      set_variable_features = set_variable_features,
      top_n = as.integer(top_n),
      conversion_params = conversion_params,
      network_params = network_params,
      binSpect_params = binSpect_params,
      seed = seed
    ),
    cells = input$cells,
    features = input$features
  )
  if (isTRUE(store_giotto)) {
    tool$giotto <- gobject
  }
  srt@tools[[tool_name]] <- tool
  log_message("{.pkg Giotto} spatial genes stored in {.code srt@tools[[tool_name]]}", message_type = "success", verbose = verbose)
  srt
}

#' @title Run Giotto HMRF spatial domains
#'
#' @description
#' Run Giotto HMRF domain detection. This wrapper may require a configured
#' Giotto Python/HMRF environment.
#'
#' @md
#' @inheritParams RunGiottoCluster
#' @param spatial_genes Features used for HMRF.
#' @param k Number of HMRF states.
#' @param betas HMRF beta values.
#' @param dimensions_to_use PCA dimensions used by HMRF.
#' @param python_path Optional Python path for Giotto HMRF.
#' @param output_folder Output folder passed to `Giotto::doHMRF()`.
#' @param hmrf_params Additional parameters passed to `Giotto::doHMRF()`.
#'
#' @return A `Seurat` object with HMRF results in `srt@tools` and any
#' extractable domain assignments in metadata.
#'
#' @export
RunGiottoHMRF <- function(
  srt,
  assay = NULL,
  layer = "data",
  spatial_genes = NULL,
  image = NULL,
  coord.cols = c("x", "y"),
  k = 10,
  betas = c(0, 2, 50),
  dimensions_to_use = 1:10,
  cluster_colname = "Giotto_HMRF",
  tool_name = "GiottoHMRF",
  python_path = NULL,
  output_folder = tempdir(),
  store_giotto = FALSE,
  conversion_params = list(),
  network_params = list(),
  hmrf_params = list(),
  verbose = TRUE,
  seed = 11
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  giotto_validate_scalar_string(cluster_colname, "cluster_colname")
  giotto_validate_scalar_string(tool_name, "tool_name")
  giotto_validate_named_list(conversion_params, "conversion_params")
  giotto_validate_named_list(network_params, "network_params")
  giotto_validate_named_list(hmrf_params, "hmrf_params")
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 2) {
    log_message("{.arg k} must be a single number >= 2", message_type = "error")
  }
  if (!is.numeric(betas) || length(betas) == 0L || any(is.na(betas))) {
    log_message("{.arg betas} must be a numeric vector", message_type = "error")
  }
  dimensions_to_use <- unique(as.integer(dimensions_to_use))
  dimensions_to_use <- dimensions_to_use[is.finite(dimensions_to_use) & dimensions_to_use > 0L]
  if (length(dimensions_to_use) == 0L) {
    log_message("{.arg dimensions_to_use} must contain at least one positive dimension", message_type = "error")
  }

  giotto_require(verbose = verbose)
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (is.null(spatial_genes)) {
    spatial_genes <- SeuratObject::VariableFeatures(srt, assay = assay)
  }
  if (length(spatial_genes) == 0L) {
    log_message("{.arg spatial_genes} must be provided or available as variable features", message_type = "error")
  }
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  }
  input <- giotto_prepare_input(
    srt = srt,
    assay = assay,
    layer = layer,
    features = spatial_genes,
    image = image,
    coord.cols = coord.cols
  )
  max_dims <- min(nrow(input$expr), ncol(input$expr)) - 1L
  dimensions_to_use <- dimensions_to_use[dimensions_to_use <= max_dims]
  if (length(dimensions_to_use) == 0L) {
    log_message("No valid {.arg dimensions_to_use} remain for Giotto PCA", message_type = "error")
  }
  gobject <- giotto_create_object(
    input = input,
    layer = layer,
    conversion_params = conversion_params,
    verbose = verbose
  )
  expression_values <- giotto_expression_values(layer)
  if (identical(expression_values, "raw")) {
    normalizeGiotto <- giotto_get_fun("normalizeGiotto")
    gobject <- giotto_call(normalizeGiotto, list(gobject = gobject, verbose = FALSE))
    expression_values <- "normalized"
  }
  runPCA <- giotto_get_fun("runPCA")
  gobject <- giotto_call(
    runPCA,
    list(
      gobject = gobject,
      expression_values = expression_values,
      feats_to_use = input$features,
      name = "pca",
      ncp = max(dimensions_to_use),
      scale_unit = TRUE,
      center = TRUE
    )
  )
  network_name <- giotto_spatial_network_name("Delaunay", network_params)
  gobject <- giotto_create_spatial_network(
    gobject = gobject,
    network_name = network_name,
    network_method = "Delaunay",
    network_params = network_params,
    verbose = verbose
  )

  hmrf_name <- hmrf_params[["name"]] %||% cluster_colname
  doHMRF <- giotto_get_fun("doHMRF")
  hmrf_result <- giotto_call(
    doHMRF,
    giotto_merge_args(
      list(
        gobject = gobject,
        expression_values = expression_values,
        spatial_network_name = network_name,
        spatial_genes = input$features,
        dim_reduction_to_use = "pca",
        dim_reduction_name = "pca",
        dimensions_to_use = dimensions_to_use,
        seed = seed,
        name = hmrf_name,
        k = as.integer(k),
        betas = betas,
        python_path = python_path,
        output_folder = output_folder
      ),
      hmrf_params
    )
  )

  assignments <- giotto_extract_hmrf_assignments(hmrf_result, cells = input$cells)
  if (length(assignments) > 0L) {
    for (nm in names(assignments)) {
      colname <- if (length(assignments) == 1L) {
        cluster_colname
      } else {
        paste0(cluster_colname, "_", giotto_clean_colname(nm))
      }
      srt <- giotto_add_cluster_metadata(
        srt = srt,
        clusters = assignments[[nm]],
        cluster_colname = colname
      )
    }
  }

  tool <- list(
    assignments = assignments,
    tables = giotto_collect_data_frames(hmrf_result),
    raw_result = hmrf_result,
    parameters = list(
      assay = assay,
      layer = layer,
      image = image,
      coord.cols = coord.cols,
      spatial_genes = input$features,
      k = as.integer(k),
      betas = betas,
      dimensions_to_use = dimensions_to_use,
      cluster_colname = cluster_colname,
      hmrf_name = hmrf_name,
      python_path = python_path,
      output_folder = output_folder,
      conversion_params = conversion_params,
      network_params = network_params,
      hmrf_params = hmrf_params,
      seed = seed
    ),
    cells = input$cells,
    features = input$features
  )
  if (isTRUE(store_giotto)) {
    tool$giotto <- hmrf_result
  }
  srt@tools[[tool_name]] <- tool
  log_message("{.pkg Giotto} HMRF results stored in {.code srt@tools[[tool_name]]}", message_type = "success", verbose = verbose)
  srt
}

#' @title Run Giotto spatial co-expression modules
#'
#' @description
#' Run Giotto spatial correlation feature detection and clustering to identify
#' feature-level spatial co-expression modules.
#'
#' @md
#' @inheritParams RunGiottoSpatialGenes
#' @param cor_method Correlation method for `Giotto::detectSpatialCorFeats()`.
#' @param k Number of spatial correlation feature clusters.
#' @param detect_params Additional parameters passed to
#' `Giotto::detectSpatialCorFeats()`.
#' @param cluster_params Additional parameters passed to
#' `Giotto::clusterSpatialCorFeats()`.
#'
#' @return A `Seurat` object with module results in `srt@tools[[tool_name]]`.
#'
#' @export
RunGiottoSpatialModules <- function(
  srt,
  assay = NULL,
  layer = "data",
  features = NULL,
  image = NULL,
  coord.cols = c("x", "y"),
  network_method = "Delaunay",
  cor_method = c("pearson", "spearman", "kendall"),
  k = 10,
  tool_name = "GiottoSpatialModules",
  store_giotto = FALSE,
  conversion_params = list(),
  network_params = list(),
  detect_params = list(),
  cluster_params = list(),
  verbose = TRUE,
  seed = 11
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  giotto_validate_scalar_string(tool_name, "tool_name")
  giotto_validate_named_list(conversion_params, "conversion_params")
  giotto_validate_named_list(network_params, "network_params")
  giotto_validate_named_list(detect_params, "detect_params")
  giotto_validate_named_list(cluster_params, "cluster_params")
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 2) {
    log_message("{.arg k} must be a single number >= 2", message_type = "error")
  }
  network_method <- match.arg(network_method, c("Delaunay", "kNN"))
  cor_method <- match.arg(cor_method)

  giotto_require(verbose = verbose)
  giotto_old_options <- giotto_disable_r_only_options()
  on.exit(options(giotto_old_options), add = TRUE)

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  input <- giotto_prepare_input(
    srt = srt,
    assay = assay,
    layer = layer,
    features = features,
    image = image,
    coord.cols = coord.cols
  )
  gobject <- giotto_create_object(
    input = input,
    layer = layer,
    conversion_params = conversion_params,
    verbose = verbose
  )
  expression_values <- giotto_expression_values(layer)
  if (identical(expression_values, "raw")) {
    normalizeGiotto <- giotto_get_fun("normalizeGiotto")
    gobject <- giotto_call(normalizeGiotto, list(gobject = gobject, verbose = FALSE))
    expression_values <- "normalized"
  }
  network_name <- giotto_spatial_network_name(network_method, network_params)
  gobject <- giotto_create_spatial_network(
    gobject = gobject,
    network_name = network_name,
    network_method = network_method,
    network_params = network_params,
    verbose = verbose
  )

  detectSpatialCorFeats <- giotto_get_fun("detectSpatialCorFeats")
  spat_cor <- giotto_call(
    detectSpatialCorFeats,
    giotto_merge_args(
      list(
        gobject = gobject,
        method = "network",
        expression_values = expression_values,
        subset_feats = input$features,
        spatial_network_name = network_name,
        cor_method = cor_method
      ),
      detect_params
    )
  )
  clusterSpatialCorFeats <- giotto_get_fun("clusterSpatialCorFeats")
  modules <- giotto_call(
    clusterSpatialCorFeats,
    giotto_merge_args(
      list(
        spatCorObject = spat_cor,
        name = cluster_params[["name"]] %||% "spat_clus",
        k = as.integer(k),
        return_obj = TRUE
      ),
      cluster_params
    )
  )
  module_tables <- giotto_collect_data_frames(modules)

  tool <- list(
    spatial_cor = spat_cor,
    modules = modules,
    module_tables = module_tables,
    parameters = list(
      assay = assay,
      layer = layer,
      image = image,
      coord.cols = coord.cols,
      network_method = network_method,
      network_name = network_name,
      cor_method = cor_method,
      k = as.integer(k),
      conversion_params = conversion_params,
      network_params = network_params,
      detect_params = detect_params,
      cluster_params = cluster_params,
      seed = seed
    ),
    cells = input$cells,
    features = input$features
  )
  if (isTRUE(store_giotto)) {
    tool$giotto <- gobject
  }
  srt@tools[[tool_name]] <- tool
  log_message("{.pkg Giotto} spatial modules stored in {.code srt@tools[[tool_name]]}", message_type = "success", verbose = verbose)
  srt
}

giotto_disable_r_only_options <- function() {
  options(
    giotto.has_conda = FALSE,
    giotto.use_conda = FALSE,
    giotto.update_param = FALSE,
    giotto.no_python_warn = TRUE
  )
}

giotto_spatial_network_name <- function(network_method, network_params = list()) {
  network_params[["name"]] %||% paste0(network_method, "_network")
}

giotto_create_spatial_network <- function(
  gobject,
  network_name,
  network_method,
  network_params = list(),
  verbose = TRUE
) {
  createSpatialNetwork <- giotto_get_fun("createSpatialNetwork")
  giotto_call(
    createSpatialNetwork,
    giotto_merge_args(
      list(
        gobject = gobject,
        name = network_name,
        method = network_method,
        return_gobject = TRUE,
        verbose = verbose
      ),
      network_params
    )
  )
}

giotto_result_to_data_frame <- function(x) {
  if (is.null(x)) {
    return(data.frame())
  }
  if (inherits(x, "data.frame")) {
    return(as.data.frame(x, stringsAsFactors = FALSE))
  }
  if (inherits(x, "data.table")) {
    return(as.data.frame(x, stringsAsFactors = FALSE))
  }
  if (is.matrix(x)) {
    return(as.data.frame(x, stringsAsFactors = FALSE))
  }
  out <- tryCatch(as.data.frame(x, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(out)) {
    return(out)
  }
  data.frame(value = as.character(x), stringsAsFactors = FALSE)
}

giotto_collect_data_frames <- function(x, prefix = "result") {
  out <- list()
  if (inherits(x, "data.frame") || inherits(x, "data.table") || is.matrix(x)) {
    out[[prefix]] <- giotto_result_to_data_frame(x)
    return(out)
  }
  if (is.list(x)) {
    for (nm in names(x)) {
      key <- if (nzchar(nm)) paste0(prefix, ".", nm) else prefix
      out <- c(out, giotto_collect_data_frames(x[[nm]], prefix = key))
    }
  }
  out
}

giotto_top_features <- function(result_table, n = 100) {
  if (nrow(result_table) == 0L) {
    return(character())
  }
  feature_col <- giotto_pick_col(result_table, c("feat_ID", "feats", "feature", "gene", "gene_ID", "name"))
  if (is.null(feature_col)) {
    return(character())
  }
  p_col <- giotto_pick_col(result_table, c("adj.p.value", "p.adj", "p_adj", "pvalue_adj", "p.value", "pval", "p_val"))
  if (!is.null(p_col)) {
    result_table <- result_table[order(result_table[[p_col]], na.last = TRUE), , drop = FALSE]
  }
  unique(as.character(result_table[[feature_col]]))[seq_len(min(n, length(unique(as.character(result_table[[feature_col]])))))]
}

giotto_extract_hmrf_assignments <- function(x, cells) {
  frames <- giotto_collect_data_frames(x)
  out <- list()
  for (nm in names(frames)) {
    df <- frames[[nm]]
    ids <- tryCatch(giotto_metadata_cell_ids(df, cells = cells), error = function(e) NULL)
    if (is.null(ids)) {
      next
    }
    for (col in setdiff(colnames(df), c("cell_ID", "cell", "spat_ID", "spot", "barcode"))) {
      values <- df[[col]]
      if (length(values) != length(ids)) {
        next
      }
      named <- values
      names(named) <- ids
      if (all(cells %in% names(named)) && length(unique(stats::na.omit(named[cells]))) > 1L) {
        out[[paste0(nm, ".", col)]] <- named[cells]
      }
    }
  }
  out
}

giotto_pick_col <- function(x, candidates) {
  nm <- colnames(x)
  hit <- candidates[tolower(candidates) %in% tolower(nm)][1]
  if (is.na(hit)) {
    return(NULL)
  }
  nm[match(tolower(hit), tolower(nm))]
}

giotto_clean_colname <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}
