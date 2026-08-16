# Official SpaTalk producer and stored-result plotting -------------------------

.spatalk_repository <- "ZJUFanLab/SpaTalk"
.spatalk_package <- "SpaTalk"
.spatalk_nnlm_repository <- "linxihui/NNLM"

spatalk_required_symbols <- function() {
  c("createSpaTalk", "dec_celltype", "find_lr_path", "dec_cci_all")
}

spatalk_check_r <- function() {
  package_status <- check_r(
    c("NNLM", .spatalk_package),
    install = FALSE,
    verbose = FALSE
  )
  package_status <- unname(unlist(package_status, use.names = FALSE))
  package_status <- suppressWarnings(as.logical(package_status))
  if (length(package_status) == 2L && !anyNA(package_status) && all(package_status)) {
    return(invisible(TRUE))
  }
  available <- check_r(c(.spatalk_nnlm_repository, .spatalk_repository), verbose = FALSE)
  available <- unname(unlist(available, use.names = FALSE))
  available <- suppressWarnings(as.logical(available))
  if (length(available) != 2L || anyNA(available) || !all(available)) {
    log_message(
      "Install the official {.pkg SpaTalk} backend from {.url https://github.com/ZJUFanLab/SpaTalk}",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

spatalk_get_fun <- function(fun) {
  if (!fun %in% spatalk_required_symbols()) {
    log_message("Unsupported {.pkg SpaTalk} function {.fn {fun}}", message_type = "error")
  }
  get_namespace_fun(.spatalk_package, fun)
}

spatalk_package_version <- function() {
  desc <- tryCatch(utils::packageDescription(.spatalk_package), error = function(e) NULL)
  if (is.null(desc)) NA_character_ else as.character(desc$Version %||% NA_character_)
}

spatalk_component <- function(object, name, default = NULL) {
  if (isS4(object) && name %in% methods::slotNames(object)) {
    return(methods::slot(object, name))
  }
  if (is.list(object) && name %in% names(object)) {
    return(object[[name]])
  }
  default
}

spatalk_load_data <- function(name, species) {
  env <- new.env(parent = emptyenv())
  utils::data(list = name, package = .spatalk_package, envir = env)
  if (!exists(name, envir = env, inherits = FALSE)) {
    log_message("The {.pkg SpaTalk} dataset {.val {name}} is unavailable", message_type = "error")
  }
  value <- get(name, envir = env, inherits = FALSE)
  value <- as.data.frame(value, stringsAsFactors = FALSE, check.names = FALSE)
  species_col <- intersect(c("species", "Species"), colnames(value))[1L]
  if (!is.na(species_col)) {
    keep <- tolower(as.character(value[[species_col]])) == tolower(species)
    value <- value[keep, , drop = FALSE]
  }
  if (nrow(value) == 0L) {
    log_message("No {.pkg SpaTalk} {.val {name}} rows are available for {.val {species}}", message_type = "error")
  }
  value
}

spatalk_call <- function(fun, args, extra = list()) {
  if (length(extra) > 0L && (is.null(names(extra)) || any(!nzchar(names(extra))))) {
    log_message("SpaTalk backend arguments in {.arg ...} must be named", message_type = "error")
  }
  formal_names <- names(formals(fun))
  if (!"..." %in% formal_names) {
    extra <- extra[intersect(names(extra), setdiff(formal_names, names(args)))]
  }
  do.call(fun, c(args, extra))
}

spatalk_matrix_values <- function(x) {
  if (methods::is(x, "sparseMatrix")) {
    return(methods::slot(x, "x"))
  }
  as.numeric(x)
}

spatalk_validate_matrix <- function(x, name) {
  if (!is.matrix(x) && !methods::is(x, "Matrix")) {
    log_message("{.arg {name}} must resolve to a matrix", message_type = "error")
  }
  if (is.null(rownames(x)) || is.null(colnames(x)) ||
      anyNA(rownames(x)) || anyNA(colnames(x)) ||
      any(!nzchar(rownames(x))) || any(!nzchar(colnames(x))) ||
      anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
    log_message("{.arg {name}} must have unique, non-missing row and column names", message_type = "error")
  }
  values <- spatalk_matrix_values(x)
  if (any(!is.finite(values)) || any(values < 0)) {
    log_message("{.arg {name}} must contain finite non-negative expression values", message_type = "error")
  }
  invisible(TRUE)
}

spatalk_input <- function(srt, group.by, assay, layer, image, coord.cols) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  if (!is.character(group.by) || length(group.by) != 1L || !group.by %in% colnames(srt[[]])) {
    log_message("{.arg group.by} must be one metadata column", message_type = "error")
  }
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message("Assay {.val {assay}} is absent from {.arg srt}", message_type = "error")
  }
  coords <- SpatialCoordinates(
    object = srt, image = image, coord.cols = coord.cols,
    space = "raw", image_policy = "strict"
  )
  cells <- as.character(coords$data$cell_id)
  if (length(cells) == 0L || anyDuplicated(cells) || any(!cells %in% colnames(srt))) {
    log_message("Spatial coordinates do not identify a valid unique subset of {.arg srt}", message_type = "error")
  }
  coords$data <- coords$data[match(cells, coords$data$cell_id), , drop = FALSE]
  if (any(!is.finite(coords$data$x)) || any(!is.finite(coords$data$y))) {
    log_message("SpaTalk requires finite raw x and y coordinates", message_type = "error")
  }
  expression <- GetAssayData5(srt, assay = assay, layer = layer)[, cells, drop = FALSE]
  spatalk_validate_matrix(expression, "srt expression")
  labels <- as.character(srt[[]][cells, group.by])
  if (anyNA(labels) || any(!nzchar(labels))) {
    log_message("{.arg group.by} contains missing or empty labels", message_type = "error")
  }
  list(
    expression = expression,
    coords = coords$data,
    source = coords$source,
    cells = cells,
    labels = labels,
    assay = assay
  )
}

spatalk_reference_input <- function(reference, reference.group.by, assay, layer) {
  if (!inherits(reference, "Seurat")) {
    log_message("Spot-mode {.arg reference} must be a {.cls Seurat} object", message_type = "error")
  }
  if (!is.character(reference.group.by) || length(reference.group.by) != 1L ||
      !reference.group.by %in% colnames(reference[[]])) {
    log_message("{.arg reference.group.by} must be one reference metadata column", message_type = "error")
  }
  assay <- assay %||% SeuratObject::DefaultAssay(reference)
  expression <- GetAssayData5(reference, assay = assay, layer = layer)
  spatalk_validate_matrix(expression, "reference expression")
  labels <- as.character(reference[[]][colnames(expression), reference.group.by])
  if (anyNA(labels) || any(!nzchar(labels))) {
    log_message("{.arg reference.group.by} contains missing or empty labels", message_type = "error")
  }
  list(expression = expression, labels = labels, assay = assay)
}

spatalk_normalize_ids <- function(x) {
  gsub("[ -]", "_", as.character(x))
}

spatalk_prepare_dec_result <- function(dec_result, spot_ids, reference_labels, source) {
  if (!is.matrix(dec_result) && !methods::is(dec_result, "Matrix")) {
    log_message(
      "SpaTalk {.arg dec_result} from {.val {source}} must be a numeric matrix",
      message_type = "error"
    )
  }
  dec_result <- as.matrix(dec_result)
  if (!is.numeric(dec_result) || nrow(dec_result) == 0L || ncol(dec_result) == 0L ||
      is.null(rownames(dec_result)) || is.null(colnames(dec_result)) ||
      anyNA(rownames(dec_result)) || anyNA(colnames(dec_result)) ||
      any(!nzchar(rownames(dec_result))) || any(!nzchar(colnames(dec_result))) ||
      anyDuplicated(rownames(dec_result)) || anyDuplicated(colnames(dec_result))) {
    log_message(
      "SpaTalk {.arg dec_result} from {.val {source}} must have unique spot and cell-type names",
      message_type = "error"
    )
  }
  if (any(!is.finite(dec_result)) || any(dec_result < 0)) {
    log_message(
      "SpaTalk {.arg dec_result} from {.val {source}} must contain finite non-negative weights",
      message_type = "error"
    )
  }
  missing_spots <- setdiff(spot_ids, rownames(dec_result))
  if (length(missing_spots) > 0L) {
    log_message(
      "SpaTalk {.arg dec_result} from {.val {source}} is missing {.val {length(missing_spots)}} selected spot{?s}",
      message_type = "error"
    )
  }
  dec_result <- dec_result[spot_ids, , drop = FALSE]
  if (any(rowSums(dec_result) <= 0)) {
    log_message(
      "SpaTalk {.arg dec_result} from {.val {source}} has selected spots with zero total weight",
      message_type = "error"
    )
  }
  normalized_types <- spatalk_normalize_ids(colnames(dec_result))
  normalized_reference <- spatalk_normalize_ids(unique(reference_labels))
  if (anyDuplicated(normalized_types) || anyDuplicated(normalized_reference)) {
    log_message(
      "SpaTalk cell-type names become duplicated after replacing spaces or hyphens with underscores",
      message_type = "error"
    )
  }
  unknown_types <- setdiff(normalized_types, normalized_reference)
  if (length(unknown_types) > 0L) {
    log_message(
      "SpaTalk {.arg dec_result} cell types are absent from {.arg reference.group.by}: {.val {paste(unknown_types, collapse = ', ')}}",
      message_type = "error"
    )
  }
  dec_result
}

spatalk_deconvolution_spec <- function(
  srt, input, reference_input, deconvolution, rctd.tool, dots
) {
  dec_result <- dots$dec_result %||% NULL
  dots$dec_result <- NULL
  if (identical(deconvolution, "NNLM")) {
    if (!is.null(dec_result)) {
      log_message(
        "Do not supply {.arg dec_result} with {.arg deconvolution = 'NNLM'}; use {.val RCTD} or {.val none}",
        message_type = "error"
      )
    }
    return(list(method = 1L, dec_result = NULL, source = "SpaTalk NNLM", dots = dots))
  }
  source <- "dec_result argument"
  if (identical(deconvolution, "RCTD") && is.null(dec_result)) {
    if (!is.character(rctd.tool) || length(rctd.tool) != 1L || is.na(rctd.tool) || !nzchar(rctd.tool)) {
      log_message("{.arg rctd.tool} must be one non-empty tool name", message_type = "error")
    }
    stored <- srt@tools[[rctd.tool]]
    dec_result <- stored$weights %||% NULL
    source <- paste0("srt@tools$", rctd.tool, "$weights")
    if (is.null(dec_result)) {
      log_message(
        "Stored RCTD weights are absent at {.code {source}}; run {.fn RunRCTD} first or provide {.arg dec_result}",
        message_type = "error"
      )
    }
  }
  if (is.null(dec_result)) {
    log_message(
      "{.arg deconvolution = 'none'} requires a named {.arg dec_result} matrix in {.arg ...}",
      message_type = "error"
    )
  }
  dec_result <- spatalk_prepare_dec_result(
    dec_result = dec_result,
    spot_ids = input$cells,
    reference_labels = reference_input$labels,
    source = source
  )
  list(method = 2L, dec_result = dec_result, source = source, dots = dots)
}

spatalk_pick_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0L) NULL else hit[[1L]]
}

spatalk_long_table <- function(lr_table) {
  lr_table <- as.data.frame(lr_table %||% data.frame(), stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(lr_table) == 0L) {
    log_message("The official {.pkg SpaTalk} backend returned no ligand-receptor interactions", message_type = "error")
  }
  required_map <- list(
    sender = c("celltype_sender", "sender", "source"),
    receiver = c("celltype_receiver", "receiver", "target"),
    ligand = "ligand", receptor = "receptor"
  )
  out <- lr_table
  for (target in names(required_map)) {
    source <- spatalk_pick_col(out, required_map[[target]])
    if (is.null(source)) {
      log_message("SpaTalk LR output is missing {.field {target}}", message_type = "error")
    }
    out[[target]] <- as.character(out[[source]])
  }
  score_col <- spatalk_pick_col(out, c("score", "lr_co_ratio", "lr_co_exp_num"))
  pvalue_col <- spatalk_pick_col(out, "lr_co_ratio_pvalue")
  pathway_col <- spatalk_pick_col(out, c("pathway_name", "pathway", "classification"))
  out$interaction_name <- paste(out$ligand, out$receptor, sep = "_")
  out$pathway_name <- if (is.null(pathway_col)) NA_character_ else as.character(out[[pathway_col]])
  out$score <- if (is.null(score_col)) NA_real_ else suppressWarnings(as.numeric(out[[score_col]]))
  out$pvalue <- if (is.null(pvalue_col)) NA_real_ else suppressWarnings(as.numeric(out[[pvalue_col]]))
  out$score_type <- if (identical(score_col, "score")) "spatalk_score" else (score_col %||% "not_reported")
  out$method <- "SpaTalk"
  standardize_long_df(out)
}

spatalk_pathway_table <- function(object, lr_table) {
  lr_path <- spatalk_component(object, "lr_path", list())
  pathways <- as.data.frame(lr_path$pathways %||% data.frame(), stringsAsFactors = FALSE)
  if (nrow(pathways) == 0L || nrow(lr_table) == 0L) {
    return(data.frame())
  }
  genes <- unique(c(as.character(lr_table$ligand), as.character(lr_table$receptor)))
  src <- spatalk_pick_col(pathways, c("src", "source", "from"))
  dest <- spatalk_pick_col(pathways, c("dest", "target", "to"))
  keep <- rep(FALSE, nrow(pathways))
  if (!is.null(src)) keep <- keep | as.character(pathways[[src]]) %in% genes
  if (!is.null(dest)) keep <- keep | as.character(pathways[[dest]]) %in% genes
  unique(pathways[keep, , drop = FALSE])
}

spatalk_validate_bundle <- function(bundle) {
  required <- c(
    "method", "results", "active_result", "long_table", "primary_table",
    "cells", "coordinates", "parameters", "summary", "provenance"
  )
  if (!is.list(bundle) || !all(required %in% names(bundle))) {
    log_message("SpaTalk result bundle is incomplete", message_type = "error")
  }
  if (!is.data.frame(bundle$long_table) || nrow(bundle$long_table) == 0L) {
    log_message("SpaTalk result bundle has no communication rows", message_type = "error")
  }
  invisible(TRUE)
}

#' @title Run official SpaTalk spatial communication analysis
#'
#' @description
#' Run the official SpaTalk R backend for single-cell spatial data or spot data
#' with a single-cell reference. Results are stored as a compact method bundle
#' and added to SCOP's unified CCC table.
#'
#' @param srt A `Seurat` spatial object.
#' @param group.by Metadata column containing cell labels for single-cell mode.
#' @param mode Analysis mode. `"auto"` selects spot mode when `reference` is
#'   supplied and single-cell mode otherwise.
#' @param reference Optional single-cell `Seurat` reference for spot mode.
#' @param reference.group.by Cell-type column in `reference`.
#' @param species SpaTalk species database.
#' @param assay,layer Spatial assay and raw-count layer.
#' @param reference.assay,reference.layer Reference assay and raw-count layer.
#' @param image Explicit spatial image. Required for multi-image objects.
#' @param coord.cols Metadata coordinate columns used when no image is present.
#' @param deconvolution Spot deconvolution mode. `"NNLM"` runs SpaTalk method 1;
#'   `"RCTD"` reuses `srt@tools[[rctd.tool]]$weights` from [RunRCTD()] (or an
#'   explicitly supplied `dec_result`); `"none"` requires a named external
#'   `dec_result` matrix in `...`. Imported matrices use SpaTalk method 2.
#' @param rctd.tool Tool key containing stored RCTD weights when
#'   `deconvolution = "RCTD"`.
#' @param result.name Stored result name.
#' @param store.object Whether to store the full upstream SpaTalk object.
#' @param overwrite Whether to replace an existing named result.
#' @param backend Backend used only for SCOP CCC result aggregation.
#' @param verbose Whether to print progress messages.
#' @param ... Named arguments forwarded to compatible official SpaTalk steps.
#'
#' @return The input `Seurat` object with `srt@tools$SpaTalk` and unified CCC
#'   results updated.
#' @references <https://github.com/ZJUFanLab/SpaTalk>
#' @concept spatial-producer
#' @export
#'
#' @examples
#' \dontrun{
#' available <- unname(unlist(thisutils::check_r(
#'   c("linxihui/NNLM", "ZJUFanLab/SpaTalk"), verbose = FALSE
#' )))
#' if (length(available) == 2L && all(available)) {
#'   spatial <- RunSpaTalk(
#'     spatial,
#'     group.by = "celltype",
#'     mode = "single_cell",
#'     image = "slice1"
#'   )
#' }
#' }
RunSpaTalk <- function(
  srt,
  group.by,
  mode = c("auto", "single_cell", "spot"),
  reference = NULL,
  reference.group.by = NULL,
  species = c("Human", "Mouse"),
  assay = NULL,
  layer = "counts",
  reference.assay = NULL,
  reference.layer = "counts",
  image = NULL,
  coord.cols = c("col", "row"),
  deconvolution = c("NNLM", "RCTD", "none"),
  rctd.tool = "RCTD",
  result.name = "default",
  store.object = c("minimal", "full"),
  overwrite = FALSE,
  backend = c("cpp", "r"),
  verbose = TRUE,
  ...
) {
  mode <- match.arg(mode)
  species <- match.arg(species)
  deconvolution <- match.arg(deconvolution)
  store.object <- match.arg(store.object)
  backend <- match.arg(backend)
  if (identical(mode, "auto")) mode <- if (is.null(reference)) "single_cell" else "spot"
  if (!is.character(result.name) || length(result.name) != 1L || is.na(result.name) || !nzchar(result.name)) {
    log_message("{.arg result.name} must be one non-empty string", message_type = "error")
  }
  validate_scalar_flag(overwrite, "overwrite")
  input <- spatalk_input(srt, group.by, assay, layer, image, coord.cols)
  existing <- srt@tools[["SpaTalk"]]
  if (!is.null(existing$results[[result.name]]) && !isTRUE(overwrite)) {
    log_message("SpaTalk result {.val {result.name}} already exists; set {.arg overwrite = TRUE}", message_type = "error")
  }
  reference_input <- NULL
  if (identical(mode, "spot")) {
    reference_input <- spatalk_reference_input(reference, reference.group.by, reference.assay, reference.layer)
    if (length(intersect(rownames(input$expression), rownames(reference_input$expression))) < 2L) {
      log_message("SpaTalk spot and reference inputs must share at least two genes", message_type = "error")
    }
  }
  dots <- list(...)
  spot_max_cell <- dots$spot_max_cell %||% if (identical(mode, "single_cell")) 1L else 30L
  deconvolution_spec <- NULL
  if (identical(mode, "spot")) {
    deconvolution_spec <- spatalk_deconvolution_spec(
      srt = srt, input = input, reference_input = reference_input,
      deconvolution = deconvolution, rctd.tool = rctd.tool, dots = dots
    )
    dots <- deconvolution_spec$dots
  }

  spatalk_check_r()
  create_fun <- spatalk_get_fun("createSpaTalk")
  dec_fun <- spatalk_get_fun("dec_celltype")
  find_fun <- spatalk_get_fun("find_lr_path")
  cci_fun <- spatalk_get_fun("dec_cci_all")
  lrpairs <- spatalk_load_data("lrpairs", species)
  pathways <- spatalk_load_data("pathways", species)

  meta <- data.frame(
    id = input$cells,
    x = as.numeric(input$coords$x),
    y = as.numeric(input$coords$y),
    stringsAsFactors = FALSE
  )
  colnames(meta)[1L] <- if (identical(mode, "single_cell")) "cell" else "spot"
  native <- spatalk_call(
    create_fun,
    list(
      st_data = input$expression,
      st_meta = meta,
      species = species,
      if_st_is_sc = identical(mode, "single_cell"),
      spot_max_cell = spot_max_cell,
      celltype = if (identical(mode, "single_cell")) input$labels else NULL
    ),
    dots
  )
  if (identical(mode, "spot")) {
    native <- spatalk_call(
      dec_fun,
      list(
        object = native,
        sc_data = reference_input$expression,
        sc_celltype = reference_input$labels,
        method = deconvolution_spec$method,
        dec_result = deconvolution_spec$dec_result
      ),
      dots
    )
  }
  native <- spatalk_call(find_fun, list(object = native, lrpairs = lrpairs, pathways = pathways), dots)
  native <- spatalk_call(cci_fun, list(object = native), dots)

  lr_table <- as.data.frame(spatalk_component(native, "lrpair", data.frame()), stringsAsFactors = FALSE)
  tf_table <- as.data.frame(spatalk_component(native, "tf", data.frame()), stringsAsFactors = FALSE)
  long_table <- spatalk_long_table(lr_table)
  pathway_table <- spatalk_pathway_table(native, lr_table)
  deconv <- spatalk_component(native, "coef", matrix())
  result <- list(
    lr_table = lr_table,
    pathway_table = pathway_table,
    tf_table = tf_table,
    deconvolution = deconv,
    long_table = long_table,
    native_object = if (identical(store.object, "full")) native else NULL
  )
  results <- existing$results %||% list()
  results[[result.name]] <- result
  parameters <- list(
    group.by = group.by, mode = mode, species = species, assay = input$assay,
    layer = layer, reference.assay = reference_input$assay %||% reference.assay,
    reference.layer = reference.layer, image = input$source$image %||% image,
    coord.cols = coord.cols, coordinate_space = "raw",
    deconvolution = deconvolution, result.name = result.name,
    deconvolution_source = deconvolution_spec$source %||% "single-cell labels",
    rctd.tool = if (identical(deconvolution, "RCTD")) rctd.tool else NULL,
    dec_result_dimensions = if (is.null(deconvolution_spec$dec_result)) NULL else dim(deconvolution_spec$dec_result),
    store.object = store.object, spot_max_cell = spot_max_cell,
    backend = backend, backend_scope = "unified CCC result aggregation",
    backend_args = dots
  )
  bundle <- list(
    method = "SpaTalk", results = results, active_result = result.name,
    long_table = long_table, primary_table = long_table,
    cells = input$cells, coordinates = input$coords,
    parameters = parameters,
    summary = list(
      result_name = result.name, mode = mode, n_interactions = nrow(long_table),
      n_pathway_edges = nrow(pathway_table), n_tf_rows = nrow(tf_table),
      store_object = store.object
    ),
    provenance = list(
      producer = "RunSpaTalk", backend_id = "spatalk",
      backend_version = spatalk_package_version(),
      repository = .spatalk_repository
    )
  )
  spatalk_validate_bundle(bundle)
  srt@tools[["SpaTalk"]] <- bundle
  srt <- ccc_update_unified_bundle(srt, method = "SpaTalk", bundle = bundle, backend = backend)
  log_message("SpaTalk analysis completed", message_type = "success", verbose = verbose)
  srt
}

spatalk_get_result <- function(object, result.name = NULL) {
  if (!inherits(object, "Seurat")) {
    log_message("{.arg object} must be a {.cls Seurat} object", message_type = "error")
  }
  bundle <- object@tools[["SpaTalk"]]
  if (is.null(bundle)) log_message("SpaTalk results are absent", message_type = "error")
  if (is.null(result.name) && length(bundle$results) > 1L) {
    log_message("Multiple SpaTalk results are stored; select {.arg result.name}", message_type = "error")
  }
  result.name <- result.name %||% bundle$active_result
  result <- bundle$results[[result.name]]
  if (is.null(result)) log_message("Unknown SpaTalk result {.val {result.name}}", message_type = "error")
  list(bundle = bundle, result = result, result.name = result.name)
}

spatalk_plot_object <- function(object, stored) {
  long_table <- stored$result$long_table %||% spatalk_long_table(stored$result$lr_table)
  bundle <- stored$bundle
  bundle$active_result <- stored$result.name
  bundle$long_table <- long_table
  bundle$primary_table <- long_table
  object@tools[["SpaTalk"]] <- bundle
  ccc_update_unified_bundle(object, method = "SpaTalk", bundle = bundle, backend = "r")
}

#' @title Plot stored SpaTalk results
#'
#' @description Plot a stored SpaTalk communication result without rerunning
#' the optional backend.
#'
#' @param object A `Seurat` object with SpaTalk results.
#' @param result.name Stored SpaTalk result name.
#' @param plot_type Network, ligand-receptor bubble, pathway bubble, or
#'   receptor-TF view.
#' @param ... Arguments passed to the corresponding SCOP CCC plot.
#'
#' @return A `ggplot` or compatible plot object.
#' @export
SpaTalkPlot <- function(
  object,
  result.name = NULL,
  plot_type = c("network", "bubble", "pathway", "tf"),
  ...
) {
  plot_type <- match.arg(plot_type)
  stored <- spatalk_get_result(object, result.name)
  if (identical(plot_type, "network")) {
    plot_object <- spatalk_plot_object(object, stored)
    return(do.call(CCCNetworkPlot, c(list(srt = plot_object, method = "SpaTalk", plot_type = "circle"), list(...))))
  }
  if (identical(plot_type, "bubble")) {
    plot_object <- spatalk_plot_object(object, stored)
    return(do.call(CCCHeatmap, c(list(srt = plot_object, method = "SpaTalk", plot_type = "bubble"), list(...))))
  }
  if (identical(plot_type, "pathway")) {
    plot_object <- spatalk_plot_object(object, stored)
    pathway <- plot_object@tools$SpaTalk$long_table$pathway_name
    if (all(is.na(pathway)) || all(!nzchar(pathway[!is.na(pathway)]))) {
      log_message("The stored SpaTalk result has no pathway labels", message_type = "error")
    }
    return(do.call(CCCHeatmap, c(list(srt = plot_object, method = "SpaTalk", plot_type = "pathway_bubble"), list(...))))
  }
  tf <- stored$result$tf_table
  if (!is.data.frame(tf) || nrow(tf) == 0L) {
    log_message("The stored SpaTalk result has no receptor-TF rows", message_type = "error")
  }
  receptor <- spatalk_pick_col(tf, "receptor")
  tf_col <- spatalk_pick_col(tf, c("tf", "target_tf", "transcription_factor"))
  receiver <- spatalk_pick_col(tf, c("receiver", "celltype_receiver"))
  score <- spatalk_pick_col(tf, c("score", "n_target", "n_hop"))
  if (is.null(receptor) || is.null(tf_col)) {
    log_message("SpaTalk receptor-TF rows are missing receptor or TF columns", message_type = "error")
  }
  plot_df <- data.frame(
    receptor = as.character(tf[[receptor]]),
    tf = as.character(tf[[tf_col]]),
    receiver = if (is.null(receiver)) "Receiver" else as.character(tf[[receiver]]),
    score = if (is.null(score)) 1 else suppressWarnings(as.numeric(tf[[score]])),
    stringsAsFactors = FALSE
  )
  plot_df$score[!is.finite(plot_df$score)] <- 1
  ggplot2::ggplot(plot_df, ggplot2::aes(
    x = .data[["receptor"]], y = .data[["tf"]],
    size = .data[["score"]], color = .data[["receiver"]]
  )) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::labs(x = "Receptor", y = "Transcription factor", color = "Receiver", size = "Score") +
    ggplot2::theme_bw()
}
