#' @title Run Giotto cell proximity enrichment
#'
#' @description
#' Use Giotto as a temporary backend to build a spatial network and test
#' pairwise enrichment between metadata groups. The complete Giotto object and
#' result tables are returned as a standalone result; the input `Seurat` object
#' is not modified.
#'
#' @md
#' @inheritParams RunGiottoCluster
#' @param group.by Seurat metadata column containing cell or spot groups.
#' @param network_method Spatial network method passed to
#' `Giotto::createSpatialNetwork()`.
#' @param network_name Name for the Giotto spatial network.
#' @param network_params Additional parameters passed to
#' `Giotto::createSpatialNetwork()`.
#' @param number_of_simulations Number of label simulations used by Giotto.
#' @param adjust_method Multiple-testing correction method.
#' @param enrichment_params Additional parameters passed to
#' `Giotto::cellProximityEnrichment()`.
#'
#' @return A `giotto2_result` list containing the full Giotto object,
#' enrichment table, raw Giotto result, parameters, features, and cells.
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- subset(
#'   visium_human_pancreas_sub,
#'   cells = colnames(visium_human_pancreas_sub)[1:120],
#'   features = rownames(visium_human_pancreas_sub)[1:400]
#' )
#' spatial$region <- ifelse(
#'   spatial$x > stats::median(spatial$x),
#'   "right",
#'   "left"
#' )
#' proximity <- list(
#'   enrichment = data.frame(
#'     group_1 = c("left", "left", "right", "right"),
#'     group_2 = c("left", "right", "left", "right"),
#'     enrichment = c(1.1, -0.7, -0.5, 1.3),
#'     type_int = c("enriched", "depleted", "depleted", "enriched")
#'   ),
#'   parameters = list(network_method = "Delaunay", number_of_simulations = 100)
#' )
#' class(proximity) <- c("giotto2_cell_proximity", "giotto2_result", "list")
#'
#' head(proximity$enrichment)
#' GiottoPlot(proximity)
#'
#' if (
#'   isTRUE(check_r("giotto-suite/Giotto", verbose = FALSE)) &&
#'     identical(Sys.getenv("SCOP_RUN_SPATIAL_BACKEND_EXAMPLES"), "true")
#' ) {
#' spatial <- Seurat::NormalizeData(spatial, assay = "Spatial", verbose = FALSE)
#' proximity <- RunGiottoCellProximity(
#'   spatial,
#'   group.by = "region",
#'   assay = "Spatial",
#'   layer = "data",
#'   coord.cols = c("x", "y"),
#'   network_method = "Delaunay",
#'   number_of_simulations = 100
#' )
#' }
#'
#' @export
RunGiottoCellProximity <- function(
  srt,
  group.by,
  assay = NULL,
  layer = "data",
  image = NULL,
  coord.cols = c("x", "y"),
  network_method = c("Delaunay", "kNN"),
  network_name = NULL,
  number_of_simulations = 1000,
  adjust_method = "fdr",
  tool_name = "GiottoCellProximity",
  store_giotto = TRUE,
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

  network_method <- match.arg(network_method)
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
  network_name <- giotto_spatial_network_name(
    network_method = network_method,
    network_name = network_name,
    network_params = network_params
  )
  spatial_network <- giotto_create_spatial_network(
    gobject = gobject,
    network_name = network_name,
    network_method = network_method,
    network_params = network_params,
    verbose = verbose
  )
  gobject <- spatial_network$gobject
  network_name <- spatial_network$network_name

  cellProximityEnrichment <- giotto_get_fun("cellProximityEnrichment")
  enrichment <- giotto_call(
    cellProximityEnrichment,
    giotto_merge_args(
      list(
        gobject = gobject,
        spat_unit = "cell",
        feat_type = "rna",
        spatial_network_name = network_name,
        cluster_column = group.by,
        number_of_simulations = as.integer(number_of_simulations),
        adjust_method = adjust_method,
        set_seed = TRUE,
        seed_number = seed
      ),
      enrichment_params,
      arg_name = "enrichment_params"
    )
  )
  enrichment_table <- giotto_extract_proximity_enrichment(enrichment)

  result <- giotto_result(
    result_type = "cell_proximity",
    giotto = gobject,
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
      tool_name = tool_name,
      store_giotto = store_giotto,
      seed = seed
    ),
    cells = input$cells,
    features = input$features
  )
  log_message("{.pkg Giotto} cell proximity enrichment completed; returning standalone result without modifying {.arg srt}", message_type = "success", verbose = verbose)
  result
}

#' @title Run Giotto spatial gene detection
#'
#' @description
#' Use Giotto `binSpect()` as a temporary backend for spatially variable gene
#' detection. The complete Giotto object and result tables are returned as a
#' standalone result; the input `Seurat` object is not modified.
#'
#' @md
#' @inheritParams RunGiottoCluster
#' @param features Features to test with `Giotto::binSpect()`. If `NULL`,
#' current variable features are used, falling back to all assay features.
#' @param network_method Spatial network method passed to
#' `Giotto::createSpatialNetwork()`.
#' @param network_name Name for the Giotto spatial network.
#' @param network_params Additional parameters passed to
#' `Giotto::createSpatialNetwork()`.
#' @param bin_method Binarization method passed to `Giotto::binSpect()`.
#' @param set_variable_features Deprecated compatibility argument. Seurat
#' variable features are never modified by this function.
#' @param top_n Number of top genes to store.
#' @param binSpect_params Additional parameters passed to `Giotto::binSpect()`.
#'
#' @return A `giotto2_result` list containing the full Giotto object,
#' spatial gene table, top features, raw Giotto result, parameters, features,
#' and cells.
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- subset(
#'   visium_human_pancreas_sub,
#'   cells = colnames(visium_human_pancreas_sub)[1:120],
#'   features = rownames(visium_human_pancreas_sub)[1:400]
#' )
#' spatial <- Seurat::NormalizeData(spatial, assay = "Spatial", verbose = FALSE)
#' giotto_genes <- list(
#'   results = data.frame(
#'     feat_ID = rownames(spatial)[1:6],
#'     spatGeneRank = c(40, 35, 28, 20, 16, 10)
#'   ),
#'   top_features = rownames(spatial)[1:4],
#'   parameters = list(assay = "Spatial", layer = "data", coord.cols = c("x", "y"))
#' )
#' class(giotto_genes) <- c("giotto2_spatial_genes", "giotto2_result", "list")
#'
#' head(giotto_genes$results)
#' GiottoPlot(giotto_genes, plot_type = "ranking", top_n = 6)
#' GiottoPlot(
#'   giotto_genes,
#'   srt = spatial,
#'   plot_type = "feature",
#'   overlay_image = FALSE,
#'   coord.cols = c("x", "y")
#' )
#'
#' if (
#'   isTRUE(check_r("giotto-suite/Giotto", verbose = FALSE)) &&
#'     identical(Sys.getenv("SCOP_RUN_SPATIAL_BACKEND_EXAMPLES"), "true")
#' ) {
#' spatial <- Seurat::FindVariableFeatures(
#'   spatial,
#'   assay = "Spatial",
#'   nfeatures = 300,
#'   verbose = FALSE
#' )
#' giotto_genes <- RunGiottoSpatialGenes(
#'   spatial,
#'   assay = "Spatial",
#'   layer = "data",
#'   features = Seurat::VariableFeatures(spatial, assay = "Spatial"),
#'   coord.cols = c("x", "y"),
#'   top_n = 50
#' )
#' }
#'
#' @export
RunGiottoSpatialGenes <- function(
  srt,
  assay = NULL,
  layer = "data",
  features = NULL,
  image = NULL,
  coord.cols = c("x", "y"),
  network_method = c("Delaunay", "kNN"),
  network_name = NULL,
  bin_method = c("kmeans", "rank"),
  set_variable_features = FALSE,
  top_n = 100,
  tool_name = "GiottoSpatialGenes",
  store_giotto = TRUE,
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

  network_method <- match.arg(network_method)
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
  network_name <- giotto_spatial_network_name(
    network_method = network_method,
    network_name = network_name,
    network_params = network_params
  )
  spatial_network <- giotto_create_spatial_network(
    gobject = gobject,
    network_name = network_name,
    network_method = network_method,
    network_params = network_params,
    verbose = verbose
  )
  gobject <- spatial_network$gobject
  network_name <- spatial_network$network_name

  binSpect <- giotto_get_fun("binSpect")
  result <- giotto_call(
    binSpect,
    giotto_merge_args(
      list(
        gobject = gobject,
        spat_unit = "cell",
        feat_type = "rna",
        bin_method = bin_method,
        expression_values = expression_values,
        subset_feats = input$features,
        spatial_network_name = network_name,
        do_parallel = FALSE,
        cores = 1,
        seed = seed,
        verbose = verbose
      ),
      binSpect_params,
      arg_name = "binSpect_params"
    )
  )
  result_table <- giotto_result_to_data_frame(result)
  top_features <- giotto_top_features(result_table, n = as.integer(top_n))

  out <- giotto_result(
    result_type = "spatial_genes",
    giotto = gobject,
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
      tool_name = tool_name,
      store_giotto = store_giotto,
      seed = seed
    ),
    cells = input$cells,
    features = input$features
  )
  log_message("{.pkg Giotto} spatial genes completed; returning standalone result without modifying {.arg srt}", message_type = "success", verbose = verbose)
  out
}

#' @title Run Giotto spatial co-expression modules
#'
#' @description
#' Use Giotto `detectSpatialCorFeats()` and `clusterSpatialCorFeats()` as a
#' temporary backend for feature-level spatial co-expression modules. The
#' complete Giotto object and module results are returned as a standalone
#' result; the input `Seurat` object is not modified.
#'
#' @md
#' @inheritParams RunGiottoCluster
#' @param features Features to test for spatial co-expression modules. If
#' `NULL`, current variable features are used, falling back to all assay
#' features.
#' @param network_method Spatial network method passed to
#' `Giotto::createSpatialNetwork()`.
#' @param network_name Name for the Giotto spatial network.
#' @param network_params Additional parameters passed to
#' `Giotto::createSpatialNetwork()`.
#' @param cor_method Correlation method used by Giotto.
#' @param k Number of feature modules passed to Giotto clustering.
#' @param detect_params Additional parameters passed to
#' `Giotto::detectSpatialCorFeats()`.
#' @param cluster_params Additional parameters passed to
#' `Giotto::clusterSpatialCorFeats()`.
#'
#' @return A `giotto2_result` list containing the full Giotto object,
#' spatial correlation object, module object, extracted module tables,
#' parameters, features, and cells.
#'
#' @examples
#' data(visium_human_pancreas_sub)
#' spatial <- subset(
#'   visium_human_pancreas_sub,
#'   cells = colnames(visium_human_pancreas_sub)[1:120],
#'   features = rownames(visium_human_pancreas_sub)[1:400]
#' )
#' module_features <- rownames(spatial)[1:4]
#' module_cor <- expand.grid(
#'   feat_ID = module_features,
#'   variable = module_features
#' )
#' module_cor$spat_cor <- c(
#'   1, 0.4, 0.1, -0.2,
#'   0.4, 1, 0.3, 0.0,
#'   0.1, 0.3, 1, 0.5,
#'   -0.2, 0.0, 0.5, 1
#' )
#' giotto_modules <- list(
#'   module_tables = list(result.cor_DT = module_cor),
#'   features = module_features,
#'   parameters = list(assay = "Spatial", layer = "data")
#' )
#' class(giotto_modules) <- c("giotto2_spatial_modules", "giotto2_result", "list")
#'
#' names(giotto_modules$module_tables)
#' GiottoPlot(giotto_modules, top_n = 4)
#'
#' if (
#'   isTRUE(check_r("giotto-suite/Giotto", verbose = FALSE)) &&
#'     identical(Sys.getenv("SCOP_RUN_SPATIAL_BACKEND_EXAMPLES"), "true")
#' ) {
#' spatial <- Seurat::NormalizeData(spatial, assay = "Spatial", verbose = FALSE)
#' spatial <- Seurat::FindVariableFeatures(
#'   spatial,
#'   assay = "Spatial",
#'   nfeatures = 300,
#'   verbose = FALSE
#' )
#' giotto_modules <- RunGiottoSpatialModules(
#'   spatial,
#'   assay = "Spatial",
#'   layer = "data",
#'   features = Seurat::VariableFeatures(spatial, assay = "Spatial")[1:50],
#'   coord.cols = c("x", "y"),
#'   cor_method = "pearson",
#'   k = 6
#' )
#' }
#'
#' @export
RunGiottoSpatialModules <- function(
  srt,
  assay = NULL,
  layer = "data",
  features = NULL,
  image = NULL,
  coord.cols = c("x", "y"),
  network_method = c("Delaunay", "kNN"),
  network_name = NULL,
  cor_method = c("pearson", "spearman", "kendall"),
  k = 10,
  tool_name = "GiottoSpatialModules",
  store_giotto = TRUE,
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

  network_method <- match.arg(network_method)
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
  network_name <- giotto_spatial_network_name(
    network_method = network_method,
    network_name = network_name,
    network_params = network_params
  )
  spatial_network <- giotto_create_spatial_network(
    gobject = gobject,
    network_name = network_name,
    network_method = network_method,
    network_params = network_params,
    verbose = verbose
  )
  gobject <- spatial_network$gobject
  network_name <- spatial_network$network_name

  detectSpatialCorFeats <- giotto_get_fun("detectSpatialCorFeats")
  spat_cor <- giotto_call(
    detectSpatialCorFeats,
    giotto_merge_args(
      list(
        gobject = gobject,
        spat_unit = "cell",
        feat_type = "rna",
        method = "network",
        expression_values = expression_values,
        subset_feats = input$features,
        spatial_network_name = network_name,
        cor_method = cor_method
      ),
      detect_params,
      arg_name = "detect_params"
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
      cluster_params,
      arg_name = "cluster_params",
      reserved = c("spatCorObject", "return_obj")
    )
  )
  module_tables <- giotto_collect_data_frames(modules)

  result <- giotto_result(
    result_type = "spatial_modules",
    giotto = gobject,
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
      tool_name = tool_name,
      store_giotto = store_giotto,
      seed = seed
    ),
    cells = input$cells,
    features = input$features
  )
  log_message("{.pkg Giotto} spatial modules completed; returning standalone result without modifying {.arg srt}", message_type = "success", verbose = verbose)
  result
}

giotto_disable_r_only_options <- function() {
  options(
    giotto.has_conda = FALSE,
    giotto.use_conda = FALSE,
    giotto.update_param = FALSE,
    giotto.no_python_warn = TRUE
  )
}

giotto_spatial_network_name <- function(network_method, network_name = NULL, network_params = list()) {
  network_name %||% network_params[["name"]] %||% paste0(network_method, "_network")
}

giotto_create_spatial_network <- function(gobject, network_name, network_method, network_params = list(), verbose = TRUE) {
  createSpatialNetwork <- giotto_get_fun("createSpatialNetwork")
  args <- giotto_merge_args(
    list(
      gobject = gobject,
      spat_unit = "cell",
      feat_type = "rna",
      name = network_name,
      method = network_method,
      return_gobject = TRUE,
      verbose = verbose
    ),
    network_params,
    arg_name = "network_params",
    reserved = c("gobject", "spat_unit", "feat_type", "return_gobject")
  )
  list(
    gobject = giotto_call(createSpatialNetwork, args),
    network_name = args[["name"]] %||% network_name
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

giotto_extract_proximity_enrichment <- function(x) {
  if (is.list(x) && !inherits(x, "data.frame") && !inherits(x, "data.table")) {
    candidates <- c(
      "enrichm_res",
      "enrichment",
      "enrichment_res",
      "cell_proximity_enrichment",
      "cellProximityEnrichment",
      "results",
      "result"
    )
    hit <- names(x)[match(tolower(candidates), tolower(names(x)), nomatch = 0L)]
    hit <- hit[hit != ""][1]
    if (!is.na(hit)) {
      return(giotto_add_interaction_columns(giotto_result_to_data_frame(x[[hit]])))
    }
    frames <- giotto_collect_data_frames(x)
    enrich_hits <- grep("enrich", names(frames), ignore.case = TRUE, value = TRUE)
    if (length(enrich_hits) > 0L) {
      return(giotto_add_interaction_columns(frames[[enrich_hits[[1L]]]]))
    }
    if (length(frames) > 0L) {
      return(giotto_add_interaction_columns(frames[[1L]]))
    }
  }
  giotto_add_interaction_columns(giotto_result_to_data_frame(x))
}

giotto_add_interaction_columns <- function(x) {
  if (nrow(x) == 0L || !"unified_int" %in% colnames(x) || all(c("group_1", "group_2") %in% colnames(x))) {
    return(x)
  }
  parts <- strsplit(as.character(x[["unified_int"]]), "--", fixed = TRUE)
  x[["group_1"]] <- vapply(parts, function(z) z[[1L]] %||% NA_character_, character(1))
  x[["group_2"]] <- vapply(parts, function(z) z[[2L]] %||% NA_character_, character(1))
  x
}

giotto_collect_data_frames <- function(x, prefix = "result") {
  out <- list()
  if (inherits(x, "data.frame") || inherits(x, "data.table") || is.matrix(x)) {
    out[[prefix]] <- giotto_result_to_data_frame(x)
    return(out)
  }
  if (is.list(x)) {
    nms <- names(x) %||% rep("", length(x))
    for (i in seq_along(x)) {
      nm <- nms[[i]]
      key <- if (nzchar(nm)) paste0(prefix, ".", nm) else paste0(prefix, ".", i)
      out <- c(out, giotto_collect_data_frames(x[[i]], prefix = key))
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
  features <- unique(as.character(result_table[[feature_col]]))
  features[seq_len(min(n, length(features)))]
}

giotto_pick_col <- function(x, candidates) {
  nm <- colnames(x)
  hit <- candidates[tolower(candidates) %in% tolower(nm)][1]
  if (is.na(hit)) {
    return(NULL)
  }
  nm[match(tolower(hit), tolower(nm))]
}
