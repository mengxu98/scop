# Spatial CellChat producer and stored-result plotting -------------------------

.spatialcellchat_repository <- "jinworks/SpatialCellChat"
.spatialcellchat_package <- "SpatialCellChat"

spatialcellchat_required_symbols <- function(analysis.level = c("cell", "spot", "composition")) {
  analysis.level <- match.arg(analysis.level)
  common <- c(
    "createSpatialCellChat", "subsetDB", "subsetData", "preProcessing",
    "identifyOverExpressedGenes", "identifyOverExpressedInteractions",
    "computeCommunProb", "filterProbability", "filterCommunication",
    "computeCommunProbPathway", "aggregateNet", "subsetCommunication"
  )
  aggregate_symbol <- if (identical(analysis.level, "composition")) {
    "computeAvgCommunProb_Visium"
  } else {
    "computeAvgCommunProb"
  }
  c(common, aggregate_symbol)
}

spatialcellchat_check_r <- function(analysis.level, verbose = FALSE) {
  check_r(.spatialcellchat_repository, verbose = verbose)
  missing <- vapply(
    spatialcellchat_required_symbols(analysis.level),
    function(symbol) {
      tryCatch(
        {
          get_namespace_fun(.spatialcellchat_package, symbol)
          FALSE
        },
        error = function(e) TRUE
      )
    },
    logical(1)
  )
  if (any(missing)) {
    log_message(
      "The installed {.pkg SpatialCellChat} API is incompatible. Missing functions: {.val {names(missing)[missing]}}",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

spatialcellchat_get_fun <- function(symbol, analysis.level = NULL) {
  if (!is.null(analysis.level) && !symbol %in% spatialcellchat_required_symbols(analysis.level)) {
    log_message(
      "Function {.fn {symbol}} is not registered for SpatialCellChat analysis level {.val {analysis.level}}",
      message_type = "error"
    )
  }
  get_namespace_fun(.spatialcellchat_package, symbol)
}

spatialcellchat_remote_info <- function() {
  desc <- tryCatch(utils::packageDescription(.spatialcellchat_package), error = function(e) NULL)
  list(
    package_version = if (is.null(desc)) NA_character_ else as.character(desc$Version %||% NA_character_),
    remote_sha = if (is.null(desc)) NA_character_ else as.character(desc$RemoteSha %||% NA_character_),
    remote_repo = if (is.null(desc)) .spatialcellchat_repository else as.character(
      desc$RemoteRepo %||% .spatialcellchat_repository
    )
  )
}

spatialcellchat_validate_scalar <- function(x, name, positive = FALSE, allow_null = FALSE) {
  if (is.null(x) && isTRUE(allow_null)) return(invisible(TRUE))
  value <- suppressWarnings(as.numeric(x))
  ok <- length(value) == 1L && is.finite(value)
  if (isTRUE(positive)) ok <- ok && value > 0
  if (!ok) {
    qualifier <- if (isTRUE(positive)) "positive finite" else "finite"
    log_message("{.arg {name}} must be one {qualifier} number", message_type = "error")
  }
  invisible(TRUE)
}

spatialcellchat_validate_input <- function(
  srt,
  group.by,
  sample.by,
  assay,
  layer,
  result.name,
  composition
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  if (!is.character(group.by) || length(group.by) != 1L || !group.by %in% colnames(srt@meta.data)) {
    log_message("{.arg group.by} must be one metadata column", message_type = "error")
  }
  if (!is.null(sample.by) && (
    !is.character(sample.by) || length(sample.by) != 1L || !sample.by %in% colnames(srt@meta.data)
  )) {
    log_message("{.arg sample.by} must be one metadata column or NULL", message_type = "error")
  }
  if (!assay %in% names(srt@assays)) {
    log_message("Assay {.val {assay}} is absent from {.arg srt}", message_type = "error")
  }
  if (!is.character(layer) || length(layer) != 1L || is.na(layer) || !nzchar(layer)) {
    log_message("{.arg layer} must be one non-empty string", message_type = "error")
  }
  if (!is.character(result.name) || length(result.name) != 1L || is.na(result.name) || !nzchar(result.name)) {
    log_message("{.arg result.name} must be one non-empty string", message_type = "error")
  }
  labels <- as.character(srt@meta.data[[group.by]])
  if (anyNA(labels) || any(!nzchar(labels))) {
    log_message("{.arg group.by} contains missing or empty labels", message_type = "error")
  }
  if (!is.null(sample.by)) {
    samples <- as.character(srt@meta.data[[sample.by]])
    if (anyNA(samples) || any(!nzchar(samples))) {
      log_message("{.arg sample.by} contains missing or empty labels", message_type = "error")
    }
  }
  if (!is.null(composition) && !is.matrix(composition) && !is.data.frame(composition)) {
    log_message("{.arg composition} must be a matrix or data frame", message_type = "error")
  }
  invisible(TRUE)
}

spatialcellchat_validate_expression <- function(expression, expected_cells) {
  if (length(dim(expression)) != 2L || nrow(expression) == 0L || ncol(expression) == 0L) {
    log_message("SpatialCellChat expression input must be a non-empty feature-by-observation matrix", message_type = "error")
  }
  if (is.null(rownames(expression)) || anyNA(rownames(expression)) || anyDuplicated(rownames(expression))) {
    log_message("SpatialCellChat expression features must have unique non-missing names", message_type = "error")
  }
  if (is.null(colnames(expression)) || anyNA(colnames(expression)) || anyDuplicated(colnames(expression))) {
    log_message("SpatialCellChat expression observations must have unique non-missing names", message_type = "error")
  }
  if (!identical(as.character(colnames(expression)), as.character(expected_cells))) {
    log_message("SpatialCellChat expression columns do not align with the Seurat observations", message_type = "error")
  }
  values <- if (inherits(expression, "sparseMatrix")) expression@x else as.numeric(expression)
  if (any(!is.finite(values))) {
    log_message("SpatialCellChat expression input contains non-finite values", message_type = "error")
  }
  if (any(values < 0)) {
    log_message("SpatialCellChat requires non-negative normalized expression", message_type = "error")
  }
  if (length(values) == 0L || all(values == 0)) {
    log_message("SpatialCellChat expression input must contain at least one positive value", message_type = "error")
  }
  invisible(TRUE)
}

spatialcellchat_detect_technology <- function(srt, technology, image = NULL) {
  technology <- match.arg(
    technology,
    c("auto", "visium", "visium_hd", "xenium", "cosmx", "merfish", "generic")
  )
  if (!identical(technology, "auto")) return(technology)

  evidence <- character()
  if (!is.null(image) && image %in% tryCatch(SeuratObject::Images(srt), error = function(e) character())) {
    image_object <- srt[[image]]
    if (inherits(image_object, c("VisiumV1", "VisiumV2"))) evidence <- c(evidence, "visium")
  }
  values <- c(
    srt@misc$technology %||% character(),
    srt@misc$platform %||% character(),
    if ("technology" %in% colnames(srt@meta.data)) unique(srt@meta.data$technology) else character(),
    if ("platform" %in% colnames(srt@meta.data)) unique(srt@meta.data$platform) else character()
  )
  values <- tolower(as.character(values))
  mapping <- list(
    visium_hd = "visium[_ -]?hd",
    visium = "visium",
    xenium = "xenium",
    cosmx = "cosmx|cos mx",
    merfish = "merfish"
  )
  for (name in names(mapping)) {
    if (any(grepl(mapping[[name]], values))) evidence <- c(evidence, name)
  }
  evidence <- unique(evidence)
  if ("visium_hd" %in% evidence) evidence <- setdiff(evidence, "visium")
  if (length(evidence) != 1L) {
    log_message(
      "Cannot determine one spatial technology from the object. Supply {.arg technology} explicitly",
      message_type = "error"
    )
  }
  evidence[[1L]]
}

spatialcellchat_has_segmentation <- function(srt, image = NULL) {
  if (is.null(image) || !image %in% tryCatch(SeuratObject::Images(srt), error = function(e) character())) {
    return(FALSE)
  }
  boundaries <- tryCatch(SeuratObject::Boundaries(srt[[image]]), error = function(e) character())
  "segmentation" %in% boundaries || inherits(srt[[image]], "FOV")
}

spatialcellchat_detect_level <- function(srt, analysis.level, technology, composition, image = NULL) {
  analysis.level <- match.arg(analysis.level, c("auto", "cell", "spot", "composition"))
  if (!identical(analysis.level, "auto")) return(analysis.level)
  if (!is.null(composition)) return("composition")
  if (spatialcellchat_has_segmentation(srt, image = image) || technology %in% c("xenium", "cosmx", "merfish")) {
    return("cell")
  }
  if (technology %in% c("visium", "visium_hd")) return("spot")
  log_message(
    "Cannot determine {.arg analysis.level}; choose {.val cell}, {.val spot}, or {.val composition}",
    message_type = "error"
  )
}

spatialcellchat_image_map <- function(srt, sample_names, image = NULL) {
  images <- tryCatch(SeuratObject::Images(srt), error = function(e) character())
  if (length(images) == 0L) return(stats::setNames(rep(list(NULL), length(sample_names)), sample_names))
  if (is.null(image)) {
    if (length(images) > 1L) {
      log_message(
        "Multiple spatial images are available; provide a named {.arg image} vector mapping samples to images",
        message_type = "error"
      )
    }
    if (length(sample_names) > 1L) {
      log_message(
        "Multiple samples require an explicit named {.arg image} mapping",
        message_type = "error"
      )
    }
    return(stats::setNames(list(images[[1L]]), sample_names))
  }
  image <- as.character(image)
  if (length(sample_names) == 1L && length(image) == 1L) {
    if (!image %in% images) log_message("Unknown spatial image {.val {image}}", message_type = "error")
    return(stats::setNames(list(image), sample_names))
  }
  if (is.null(names(image)) || !all(sample_names %in% names(image))) {
    log_message("Multi-sample {.arg image} must be named for every sample", message_type = "error")
  }
  image <- image[sample_names]
  missing <- setdiff(unname(image), images)
  if (length(missing) > 0L) {
    log_message("Unknown spatial images: {.val {missing}}", message_type = "error")
  }
  stats::setNames(as.list(unname(image)), sample_names)
}

spatialcellchat_visium_spot_diameter <- function(srt, image) {
  candidates <- list(
    srt@misc$scalefactors_json[[image]]$spot_diameter_fullres,
    srt@misc$spatial[[image]]$scalefactors$spot_diameter_fullres,
    srt@misc$spatial[[image]]$spot_diameter_fullres
  )
  image_misc <- if (!is.null(image) && image %in% tryCatch(SeuratObject::Images(srt), error = function(e) character())) {
    tryCatch(srt[[image]]@misc, error = function(e) list())
  } else {
    list()
  }
  candidates <- c(
    candidates,
    list(image_misc$spot_diameter_fullres, image_misc$scalefactors$spot_diameter_fullres)
  )
  values <- suppressWarnings(as.numeric(unlist(candidates, use.names = FALSE)))
  values <- unique(values[is.finite(values) & values > 0])
  if (length(values) != 1L) return(NULL)
  values[[1L]]
}

spatialcellchat_metric_coordinates <- function(
  srt,
  cells,
  image,
  coord.cols,
  technology,
  coordinate.unit,
  ratio,
  tol
) {
  raw <- SpatialCoordinates(
    object = srt,
    image = image,
    coord.cols = coord.cols,
    space = "raw",
    image_policy = "strict"
  )
  display <- SpatialCoordinates(
    object = srt,
    image = image,
    coord.cols = coord.cols,
    space = "display",
    image_policy = "strict"
  )
  coords <- raw$data[raw$data$cell_id %in% cells, , drop = FALSE]
  coords <- coords[match(cells, coords$cell_id), , drop = FALSE]
  if (nrow(coords) != length(cells) || anyNA(coords$cell_id) || !identical(as.character(coords$cell_id), cells)) {
    log_message("Spatial coordinates do not align one-to-one with selected cells or spots", message_type = "error")
  }
  display_coords <- display$data[match(cells, display$data$cell_id), , drop = FALSE]
  if (
    nrow(display_coords) != length(cells) || anyNA(display_coords$cell_id) ||
      !identical(as.character(display_coords$cell_id), cells)
  ) {
    log_message("Display coordinates do not align one-to-one with selected cells or spots", message_type = "error")
  }
  coordinate.unit <- match.arg(coordinate.unit, c("auto", "pixel", "micron"))
  if (identical(coordinate.unit, "auto")) {
    coordinate.unit <- if (technology %in% c("visium", "visium_hd")) {
      "pixel"
    } else if (technology %in% c("xenium", "cosmx", "merfish")) {
      "micron"
    } else {
      log_message(
        "Generic spatial data require explicit {.arg coordinate.unit}",
        message_type = "error"
      )
    }
  }
  unit_source <- "user"
  if (identical(coordinate.unit, "micron")) {
    if (!is.null(ratio) && !isTRUE(all.equal(as.numeric(ratio), 1))) {
      log_message("Micron coordinates require {.arg ratio = 1}", message_type = "error")
    }
    ratio <- 1
    unit_source <- if (technology %in% c("xenium", "cosmx", "merfish")) "technology" else "user"
  } else if (is.null(ratio)) {
    if (!technology %in% c("visium", "visium_hd")) {
      log_message("Pixel coordinates require an explicit positive {.arg ratio}", message_type = "error")
    }
    diameter <- spatialcellchat_visium_spot_diameter(srt, image)
    if (is.null(diameter)) {
      log_message(
        "Visium pixel coordinates require {.arg ratio} or one trusted {.val spot_diameter_fullres} value",
        message_type = "error"
      )
    }
    ratio <- 65 / diameter
    unit_source <- "visium_scalefactors"
  }
  spatialcellchat_validate_scalar(ratio, "ratio", positive = TRUE)
  if (is.null(tol)) {
    if (technology %in% c("visium", "visium_hd")) {
      tol <- 32.5
    } else {
      log_message(
        "Cell-resolved or generic data require {.arg tol} in microns unless a trusted size is available",
        message_type = "error"
      )
    }
  }
  spatialcellchat_validate_scalar(tol, "tol", positive = TRUE)
  coords$x_raw <- coords$x
  coords$y_raw <- coords$y
  coords$x_display <- display_coords$x
  coords$y_display <- display_coords$y
  coords$x <- coords$x * as.numeric(ratio)
  coords$y <- coords$y * as.numeric(ratio)
  if (any(!is.finite(coords$x)) || any(!is.finite(coords$y))) {
    log_message("Coordinate conversion produced non-finite micron coordinates", message_type = "error")
  }
  source <- utils::modifyList(raw$source, list(
    coordinate_unit = coordinate.unit,
    unit_source = unit_source,
    scale_to_micron = as.numeric(ratio),
    analysis_unit = "micron",
    tol_um = as.numeric(tol)
  ))
  list(data = coords, source = source, spatial.factors = list(ratio = 1, tol = as.numeric(tol)))
}

spatialcellchat_validate_composition <- function(
  composition,
  cells,
  groups = NULL,
  normalize = TRUE
) {
  if (is.null(composition)) {
    log_message("{.arg composition} is required for composition-level analysis", message_type = "error")
  }
  composition <- as.matrix(composition)
  storage.mode(composition) <- "double"
  if (is.null(rownames(composition)) || anyNA(rownames(composition)) || anyDuplicated(rownames(composition))) {
    log_message("{.arg composition} must have unique spot identifiers as row names", message_type = "error")
  }
  if (is.null(colnames(composition)) || anyNA(colnames(composition)) || anyDuplicated(colnames(composition))) {
    log_message("{.arg composition} must have unique cell types as column names", message_type = "error")
  }
  missing <- setdiff(cells, rownames(composition))
  if (length(missing) > 0L) {
    log_message("Composition is missing selected spots: {.val {missing}}", message_type = "error")
  }
  composition <- composition[cells, , drop = FALSE]
  if (!is.null(groups)) {
    groups <- as.character(groups)
    missing_groups <- setdiff(groups, colnames(composition))
    extra_groups <- setdiff(colnames(composition), groups)
    if (length(missing_groups) > 0L || length(extra_groups) > 0L) {
      log_message(
        "Composition cell types must exactly match {.arg group.by}; missing: {.val {missing_groups}}, extra: {.val {extra_groups}}",
        message_type = "error"
      )
    }
    composition <- composition[, groups, drop = FALSE]
  }
  if (any(!is.finite(composition)) || any(composition < 0)) {
    log_message("{.arg composition} must contain finite non-negative values", message_type = "error")
  }
  totals <- rowSums(composition)
  if (any(!is.finite(totals)) || any(totals <= 0)) {
    log_message("Every composition row must have a positive finite sum", message_type = "error")
  }
  normalized <- FALSE
  if (!isTRUE(all.equal(totals, rep(1, length(totals)), tolerance = 1e-6))) {
    if (!isTRUE(normalize)) {
      log_message("Composition rows must sum to one when {.arg composition.normalize = FALSE}", message_type = "error")
    }
    composition <- composition / totals
    normalized <- TRUE
  }
  list(data = composition, normalized = normalized)
}

spatialcellchat_database <- function(species, database, custom.db = NULL) {
  if (identical(database, "custom")) {
    if (is.null(custom.db)) log_message("{.arg custom.db} is required for a custom database", message_type = "error")
    return(custom.db)
  }
  data_name <- switch(species,
    Homo_sapiens = "CellChatDB.human",
    Mus_musculus = "CellChatDB.mouse"
  )
  data_env <- new.env(parent = emptyenv())
  utils::data(list = data_name, package = .spatialcellchat_package, envir = data_env)
  db <- get0(data_name, envir = data_env, inherits = FALSE)
  if (is.null(db)) {
    log_message("Cannot load SpatialCellChat database {.val {data_name}}", message_type = "error")
  }
  if (identical(database, "protein")) {
    db <- spatialcellchat_get_fun("subsetDB")(
      db,
      search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"),
      non_protein = FALSE
    )
  }
  db
}

spatialcellchat_call <- function(symbol, args, analysis.level) {
  fun <- spatialcellchat_get_fun(symbol, analysis.level = analysis.level)
  formals_names <- names(formals(fun))
  unsupported <- setdiff(names(args), c(formals_names, if ("..." %in% formals_names) names(args) else character()))
  if (length(unsupported) > 0L && !"..." %in% formals_names) {
    log_message(
      "Installed SpatialCellChat function {.fn {symbol}} does not support required arguments: {.val {unsupported}}",
      message_type = "error"
    )
  }
  do.call(fun, args)
}

spatialcellchat_extract_table <- function(chat, sample, analysis.level, do.permutation) {
  table <- spatialcellchat_get_fun("subsetCommunication")(chat)
  table <- as.data.frame(table, stringsAsFactors = FALSE)
  if (nrow(table) == 0L) {
    log_message("SpatialCellChat returned no group-level communication records", message_type = "error")
  }
  table <- standardize_long_df(table)
  table$method <- "SpatialCellChat"
  table$sample <- sample
  table$modality <- "spatial"
  table$analysis_level <- analysis.level
  table$spatially_constrained <- TRUE
  table$significance_basis <- if (isTRUE(do.permutation)) "spatial_permutation" else "not_tested"
  table$passes_filter <- TRUE
  if (!isTRUE(do.permutation)) {
    table$pvalue <- NA_real_
    table$significant <- NA
    table$neglog10_pvalue <- NA_real_
  }
  table
}

spatialcellchat_extract_pathways <- function(table) {
  valid <- !is.na(table$pathway_name) & nzchar(as.character(table$pathway_name))
  table <- table[valid, , drop = FALSE]
  if (nrow(table) == 0L) return(data.frame())
  stats::aggregate(
    table$score,
    by = list(
      sender = table$sender,
      receiver = table$receiver,
      pathway_name = table$pathway_name,
      sample = table$sample
    ),
    FUN = sum,
    na.rm = TRUE
  ) |>
    stats::setNames(c("sender", "receiver", "pathway_name", "sample", "score"))
}

spatialcellchat_extract_network <- function(table) {
  if (nrow(table) == 0L) return(data.frame())
  score <- stats::aggregate(
    table$score,
    by = list(sender = table$sender, receiver = table$receiver, sample = table$sample),
    FUN = sum,
    na.rm = TRUE
  )
  count <- stats::aggregate(
    rep(1, nrow(table)),
    by = list(sender = table$sender, receiver = table$receiver, sample = table$sample),
    FUN = sum
  )
  names(score)[[4L]] <- "weight"
  names(count)[[4L]] <- "count"
  merge(score, count, by = c("sender", "receiver", "sample"), all = TRUE, sort = FALSE)
}

spatialcellchat_run_one <- function(
  expression,
  metadata,
  coordinates,
  spatial.factors,
  database.use,
  analysis.level,
  composition,
  interaction.range,
  contact.dependent,
  contact.range,
  scale.distance,
  min.cells,
  min.links,
  avg.type,
  do.permutation,
  nboot,
  seed.use
) {
  seed_was_present <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (seed_was_present) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (seed_was_present) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed.use))
  chat <- spatialcellchat_call("createSpatialCellChat", list(
    object = expression,
    meta = metadata,
    group.by = "labels",
    datatype = "spatial",
    coordinates = coordinates,
    spatial.factors = spatial.factors
  ), analysis.level)
  chat@DB <- database.use
  chat <- spatialcellchat_call("subsetData", list(object = chat), analysis.level)
  chat <- spatialcellchat_call("preProcessing", list(object = chat), analysis.level)
  chat <- spatialcellchat_call("identifyOverExpressedGenes", list(
    object = chat,
    selection.method = "meringue",
    do.grid = FALSE
  ), analysis.level)
  chat <- spatialcellchat_call("identifyOverExpressedInteractions", list(
    object = chat,
    variable.both = FALSE
  ), analysis.level)
  probability_args <- list(
    object = chat,
    raw.use = TRUE,
    distance.use = TRUE,
    scale.distance = scale.distance,
    interaction.range = interaction.range,
    contact.dependent = contact.dependent
  )
  if (isTRUE(contact.dependent)) probability_args$contact.range <- contact.range
  chat <- spatialcellchat_call("computeCommunProb", probability_args, analysis.level)
  chat <- spatialcellchat_call("filterProbability", list(
    object = chat,
    nboot = nboot,
    seed.use = seed.use,
    thresh = 0.05
  ), analysis.level)
  chat <- spatialcellchat_call("filterCommunication", list(
    object = chat,
    min.cells = NULL,
    min.links = min.links,
    min.cells.sr = min.cells
  ), analysis.level)
  if (identical(analysis.level, "composition")) {
    chat <- spatialcellchat_call("computeAvgCommunProb_Visium", list(
      object = chat,
      cell.type.decomposition = composition,
      avg.type = avg.type,
      do.permutation = do.permutation,
      nboot = nboot,
      seed.use = seed.use
    ), analysis.level)
  } else {
    chat <- spatialcellchat_call("computeAvgCommunProb", list(
      object = chat,
      group.by = "labels",
      avg.type = avg.type,
      min.cells.sr = min.cells,
      do.permutation = do.permutation,
      nboot = nboot,
      seed.use = seed.use
    ), analysis.level)
  }
  chat <- spatialcellchat_call("filterCommunication", list(
    object = chat,
    min.cells = min.cells,
    min.links = NULL,
    min.cells.sr = NULL
  ), analysis.level)
  chat <- spatialcellchat_call("computeCommunProbPathway", list(object = chat), analysis.level)
  spatialcellchat_call("aggregateNet", list(object = chat), analysis.level)
}

#' @title Run Spatial CellChat analysis
#'
#' @description
#' Run the SpatialCellChat v3 backend with explicit micron-scale coordinates,
#' sample isolation, truthful provenance, and schema-v1 storage. Cell-, spot-,
#' and composition-level results retain distinct interpretations.
#'
#' @param srt A `Seurat` object.
#' @param group.by Metadata column defining cell types or spot domains.
#' @param sample.by Optional metadata column defining independent spatial samples.
#' @param assay,layer Assay and normalized expression layer.
#' @param image Spatial image name, or a named character vector mapping samples
#'   to images. Multi-image objects require an explicit selection.
#' @param coord.cols Metadata coordinate columns when no image is present.
#' @param species Species database.
#' @param database Signaling database scope.
#' @param custom.db Custom SpatialCellChat database used with `database = "custom"`.
#' @param technology Spatial technology. Automatic detection is strict.
#' @param analysis.level Scientific interpretation of observations.
#' @param composition Spot-by-cell-type composition matrix for composition mode.
#' @param composition.normalize Normalize composition rows to one when needed.
#' @param coordinate.unit Unit of the raw coordinates.
#' @param ratio Raw-coordinate to micron multiplier.
#' @param tol Spatial tolerance in microns.
#' @param interaction.range Maximum signaling range in microns.
#' @param contact.dependent Whether to infer contact-dependent signaling.
#' @param contact.range Contact range in microns. Required when contact signaling
#'   is enabled for cell-level data.
#' @param scale.distance SpatialCellChat distance scaling parameter.
#' @param min.cells,min.links Minimum group sizes and individual links.
#' @param avg.type Group-level probability aggregation method.
#' @param do.permutation,nboot,seed.use Group permutation settings.
#' @param result.name Stored result name.
#' @param store.object Store only normalized results or also native backend objects.
#' @param overwrite Whether an existing named result may be replaced.
#' @param backend Backend for scop's unified CCC post-processing only.
#' @param verbose Whether to print progress messages.
#'
#' @return The input `Seurat` object with SpatialCellChat and unified CCC results.
#'
#' @references
#' [SpatialCellChat](https://github.com/jinworks/SpatialCellChat)
#'
#' @seealso [GetSpatialResult], [GetCCCObject], [SpatialCellChatPlot],
#' [CCCNetworkPlot], [CCCHeatmap], [CCCStatPlot]
#'
#' @export
RunSpatialCellChat <- function(
  srt,
  group.by,
  sample.by = NULL,
  assay = NULL,
  layer = "data",
  image = NULL,
  coord.cols = c("col", "row"),
  species = c("Homo_sapiens", "Mus_musculus"),
  database = c("protein", "all", "custom"),
  custom.db = NULL,
  technology = c("auto", "visium", "visium_hd", "xenium", "cosmx", "merfish", "generic"),
  analysis.level = c("auto", "cell", "spot", "composition"),
  composition = NULL,
  composition.normalize = TRUE,
  coordinate.unit = c("auto", "pixel", "micron"),
  ratio = NULL,
  tol = NULL,
  interaction.range = 250,
  contact.dependent = FALSE,
  contact.range = NULL,
  scale.distance = 0.2,
  min.cells = 10,
  min.links = 10,
  avg.type = c("avg", "sum"),
  do.permutation = TRUE,
  nboot = 100,
  seed.use = 1L,
  result.name = "default",
  store.object = c("minimal", "full"),
  overwrite = FALSE,
  backend = c("cpp", "r"),
  verbose = TRUE
) {
  assay <- assay %||% DefaultAssay(srt)
  species <- match.arg(species)
  database <- match.arg(database)
  technology <- match.arg(technology)
  analysis.level <- match.arg(analysis.level)
  coordinate.unit <- match.arg(coordinate.unit)
  avg.type <- match.arg(avg.type)
  store.object <- match.arg(store.object)
  backend <- match.arg(backend)
  spatialcellchat_validate_input(srt, group.by, sample.by, assay, layer, result.name, composition)
  spatialcellchat_validate_scalar(interaction.range, "interaction.range", positive = TRUE)
  spatialcellchat_validate_scalar(scale.distance, "scale.distance", positive = TRUE)
  spatialcellchat_validate_scalar(min.cells, "min.cells", positive = TRUE)
  spatialcellchat_validate_scalar(min.links, "min.links", positive = TRUE)
  spatialcellchat_validate_scalar(nboot, "nboot", positive = TRUE)
  if (!is.logical(contact.dependent) || length(contact.dependent) != 1L || is.na(contact.dependent)) {
    log_message("{.arg contact.dependent} must be TRUE or FALSE", message_type = "error")
  }
  if (!is.logical(do.permutation) || length(do.permutation) != 1L || is.na(do.permutation)) {
    log_message("{.arg do.permutation} must be TRUE or FALSE", message_type = "error")
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    log_message("{.arg overwrite} must be TRUE or FALSE", message_type = "error")
  }
  existing <- srt@tools[["SpatialCellChat"]]
  if (!is.null(existing$results[[result.name]]) && !isTRUE(overwrite)) {
    log_message(
      "SpatialCellChat result {.val {result.name}} already exists; set {.arg overwrite = TRUE} to replace it",
      message_type = "error"
    )
  }

  sample_values <- if (is.null(sample.by)) rep("ALL", ncol(srt)) else as.character(srt@meta.data[[sample.by]])
  names(sample_values) <- colnames(srt)
  sample_names <- unique(sample_values)
  image_map <- spatialcellchat_image_map(srt, sample_names, image = image)
  technology_requested <- technology
  if (identical(technology_requested, "auto")) {
    detected_technology <- unique(vapply(sample_names, function(sample_name) {
      spatialcellchat_detect_technology(srt, "auto", image = image_map[[sample_name]])
    }, character(1)))
    if (length(detected_technology) != 1L) {
      log_message(
        "Samples resolve to different spatial technologies; split the run or specify a consistent {.arg technology}",
        message_type = "error"
      )
    }
    technology <- detected_technology[[1L]]
  }
  analysis_level_requested <- analysis.level
  if (identical(analysis_level_requested, "auto")) {
    detected_level <- unique(vapply(sample_names, function(sample_name) {
      spatialcellchat_detect_level(
        srt,
        "auto",
        technology,
        composition,
        image = image_map[[sample_name]]
      )
    }, character(1)))
    if (length(detected_level) != 1L) {
      log_message(
        "Samples resolve to different analysis levels; split the run or specify {.arg analysis.level}",
        message_type = "error"
      )
    }
    analysis.level <- detected_level[[1L]]
  }
  if (identical(analysis.level, "composition") && is.null(composition)) {
    log_message("Composition-level analysis requires {.arg composition}", message_type = "error")
  }
  if (!identical(analysis.level, "composition") && !is.null(composition)) {
    log_message("{.arg composition} is only valid with composition-level analysis", message_type = "error")
  }
  if (identical(analysis.level, "spot") && isTRUE(contact.dependent)) {
    log_message(
      "Spot-level observations are not cells; contact-dependent signaling is disabled. Use cell-level segmentation for contact signaling",
      message_type = "error"
    )
  }
  if (isTRUE(contact.dependent)) {
    spatialcellchat_validate_scalar(contact.range, "contact.range", positive = TRUE)
  }
  spatialcellchat_check_r(analysis.level, verbose = FALSE)
  database.use <- spatialcellchat_database(species, database, custom.db = custom.db)
  expression_all <- tryCatch(
    GetAssayData5(srt, assay = assay, layer = layer),
    error = function(e) NULL
  )
  if (is.null(expression_all)) {
    log_message("Cannot read normalized expression from assay {.val {assay}} layer {.val {layer}}", message_type = "error")
  }
  spatialcellchat_validate_expression(expression_all, colnames(srt))

  log_message(
    "Start SpatialCellChat {.val {analysis.level}} analysis using {.val {technology}} coordinates",
    verbose = verbose
  )
  sample_results <- list()
  long_tables <- list()
  source_by_sample <- list()
  for (sample_name in sample_names) {
    cells <- names(sample_values)[sample_values == sample_name]
    metric <- spatialcellchat_metric_coordinates(
      srt = srt,
      cells = cells,
      image = image_map[[sample_name]],
      coord.cols = coord.cols,
      technology = technology,
      coordinate.unit = coordinate.unit,
      ratio = ratio,
      tol = tol
    )
    metric$source <- utils::modifyList(metric$source, list(
      image = image_map[[sample_name]],
      image_policy = "strict",
      technology = technology,
      analysis_level = analysis.level,
      interaction_range_um = as.numeric(interaction.range),
      contact_range_um = if (isTRUE(contact.dependent)) as.numeric(contact.range) else NA_real_,
      scale_distance = as.numeric(scale.distance)
    ))
    metadata <- data.frame(
      labels = as.character(srt@meta.data[cells, group.by]),
      samples = sample_name,
      row.names = cells,
      stringsAsFactors = FALSE
    )
    expression <- expression_all[, cells, drop = FALSE]
    if (!identical(colnames(expression), cells)) {
      log_message("Expression columns do not align with selected cells or spots", message_type = "error")
    }
    composition_i <- NULL
    composition_normalized <- FALSE
    if (identical(analysis.level, "composition")) {
      comp <- spatialcellchat_validate_composition(
        composition,
        cells = cells,
        groups = sort(unique(metadata$labels)),
        normalize = composition.normalize
      )
      composition_i <- comp$data
      composition_normalized <- comp$normalized
    }
    log_message(
      "Sample {.val {sample_name}}: coordinates converted from {.val {metric$source$coordinate_unit}} with {.val {metric$source$scale_to_micron}} micron per unit; interaction range {.val {interaction.range}} micron",
      verbose = verbose
    )
    chat <- spatialcellchat_run_one(
      expression = expression,
      metadata = metadata,
      coordinates = as.matrix(metric$data[, c("x", "y"), drop = FALSE]),
      spatial.factors = metric$spatial.factors,
      database.use = database.use,
      analysis.level = analysis.level,
      composition = composition_i,
      interaction.range = interaction.range,
      contact.dependent = contact.dependent,
      contact.range = contact.range,
      scale.distance = scale.distance,
      min.cells = as.integer(min.cells),
      min.links = as.integer(min.links),
      avg.type = avg.type,
      do.permutation = do.permutation,
      nboot = as.integer(nboot),
      seed.use = as.integer(seed.use)
    )
    table <- spatialcellchat_extract_table(chat, sample_name, analysis.level, do.permutation)
    coords_store <- metric$data[, c(
      "cell_id", "x", "y", "x_raw", "y_raw", "x_display", "y_display", "image"
    ), drop = FALSE]
    coords_store$label <- metadata[coords_store$cell_id, "labels"]
    native_object <- if (identical(store.object, "full")) chat else NULL
    if (!is.null(native_object)) {
      size <- as.numeric(utils::object.size(native_object))
      if (is.finite(size) && size > 500 * 1024^2) {
        log_message(
          "Native SpatialCellChat object for sample {.val {sample_name}} is {.val {round(size / 1024^2, 1)}} MiB and will enlarge serialized Seurat objects",
          message_type = "warning",
          verbose = verbose
        )
      }
    }
    sample_results[[sample_name]] <- list(
      interactions = table,
      pathways = spatialcellchat_extract_pathways(table),
      network = spatialcellchat_extract_network(table),
      coordinates = coords_store,
      native_object = native_object,
      diagnostics = list(
        interpretation = switch(analysis.level,
          cell = "cell-level communication",
          spot = "spot/domain-level communication",
          composition = "composition-aware spot communication"
        ),
        composition_normalized = composition_normalized,
        n_observations = length(cells),
        n_interactions = nrow(table)
      ),
      source = metric$source
    )
    long_tables[[sample_name]] <- table
    source_by_sample[[sample_name]] <- metric$source
  }
  long_table <- do.call(rbind, long_tables)
  rownames(long_table) <- NULL
  results <- existing$results %||% list()
  results[[result.name]] <- sample_results
  parameters <- list(
    group.by = group.by,
    sample.by = sample.by,
    assay = assay,
    layer = layer,
    image = image,
    coord.cols = coord.cols,
    species = species,
    database = database,
    technology = technology,
    analysis.level = analysis.level,
    coordinate.unit = coordinate.unit,
    ratio = ratio,
    tol = tol,
    interaction.range = interaction.range,
    contact.dependent = contact.dependent,
    contact.range = contact.range,
    scale.distance = scale.distance,
    min.cells = min.cells,
    min.links = min.links,
    avg.type = avg.type,
    do.permutation = do.permutation,
    nboot = nboot,
    seed.use = seed.use,
    result.name = result.name,
    store.object = store.object
  )
  remote <- spatialcellchat_remote_info()
  bundle <- spatial_result_build(
    bundle = list(
      results = results,
      active_result = result.name,
      long_table = long_table,
      cells = unique(unlist(lapply(sample_results, function(x) x$coordinates$cell_id), use.names = FALSE))
    ),
    method = "SpatialCellChat",
    result_type = "communication",
    source = list(
      image = unname(unlist(image_map, use.names = FALSE)),
      coordinate_space = "raw",
      analysis_unit = "micron",
      samples = source_by_sample,
      assay = assay,
      layer = layer
    ),
    provenance = list(
      producer = "RunSpatialCellChat",
      backend_id = "spatialcellchat",
      backend_versions = stats::setNames(remote$package_version, .spatialcellchat_package),
      remote_sha = remote$remote_sha,
      remote_repo = remote$remote_repo
    ),
    parameters = parameters,
    summary = list(
      result_name = result.name,
      samples = sample_names,
      analysis_level = analysis.level,
      n_interactions = nrow(long_table),
      store_object = store.object
    )
  )
  srt@tools[["SpatialCellChat"]] <- bundle
  srt <- ccc_update_unified_bundle(
    srt = srt,
    method = "SpatialCellChat",
    bundle = bundle,
    thresh = 0.05,
    backend = backend
  )
  log_message("SpatialCellChat analysis completed", message_type = "success", verbose = verbose)
  srt
}

#' @title Get a native CellChat-family object
#'
#' @param object A `Seurat` object.
#' @param method Either `"CellChat"` or `"SpatialCellChat"`.
#' @param result.name A CellChat condition or SpatialCellChat named result.
#' @param sample Spatial sample for SpatialCellChat results.
#'
#' @return A native backend object.
#' @export
GetCCCObject <- function(
  object,
  method = c("CellChat", "SpatialCellChat"),
  result.name = NULL,
  sample = NULL
) {
  if (!inherits(object, "Seurat")) {
    log_message("{.arg object} must be a {.cls Seurat} object", message_type = "error")
  }
  method <- match.arg(method)
  if (identical(method, "CellChat")) {
    if (!is.null(sample)) log_message("{.arg sample} is only used for SpatialCellChat", message_type = "error")
    return(get_single_cc_obj(object, condition = result.name))
  }
  bundle <- object@tools[["SpatialCellChat"]]
  if (is.null(bundle)) log_message("SpatialCellChat results are absent", message_type = "error")
  result.name <- result.name %||% bundle$active_result
  result <- bundle$results[[result.name]]
  if (is.null(result)) log_message("Unknown SpatialCellChat result {.val {result.name}}", message_type = "error")
  samples <- names(result)
  if (is.null(sample)) {
    if (length(samples) != 1L) {
      log_message("Multiple spatial samples are stored; select {.arg sample}: {.val {samples}}", message_type = "error")
    }
    sample <- samples[[1L]]
  }
  if (!sample %in% samples) log_message("Unknown SpatialCellChat sample {.val {sample}}", message_type = "error")
  native <- result[[sample]]$native_object
  if (is.null(native)) {
    log_message(
      "Native SpatialCellChat object was not stored. Rerun with {.arg store.object = 'full'}",
      message_type = "error"
    )
  }
  native
}

spatialcellchat_get_stored_sample <- function(object, result.name = NULL, sample = NULL) {
  bundle <- GetSpatialResult(object, method = "RunSpatialCellChat")
  result.name <- result.name %||% bundle$active_result
  result <- bundle$results[[result.name]]
  if (is.null(result)) log_message("Unknown SpatialCellChat result {.val {result.name}}", message_type = "error")
  samples <- names(result)
  if (is.null(sample)) {
    if (length(samples) != 1L) {
      log_message("Multiple samples are stored; select {.arg sample}: {.val {samples}}", message_type = "error")
    }
    sample <- samples[[1L]]
  }
  if (!sample %in% samples) log_message("Unknown sample {.val {sample}}", message_type = "error")
  list(result = result[[sample]], sample = sample, result.name = result.name)
}

#' @title Plot stored SpatialCellChat results
#'
#' @description
#' Plot normalized SpatialCellChat results without rerunning the backend. The
#' network view is a group-aggregated network over group centroids, not a
#' materialized individual cell-by-cell edge table.
#'
#' @param object A `Seurat` object with SpatialCellChat results.
#' @param result.name Stored result name.
#' @param sample Spatial sample. Required when multiple samples are stored.
#' @param plot_type Spatial result view.
#' @param signaling Optional pathway name.
#' @param pairLR.use Optional interaction name.
#' @param direction Score direction for spatial score plots.
#' @param top_n Maximum network edges.
#' @param point.size Coordinate point size.
#' @param palette,palcolor Color palette.
#' @param title Optional title.
#'
#' @return A `ggplot` object.
#' @export
SpatialCellChatPlot <- function(
  object,
  result.name = NULL,
  sample = NULL,
  plot_type = c("spatial_network", "lr_spatial", "pathway", "incoming", "outgoing"),
  signaling = NULL,
  pairLR.use = NULL,
  direction = c("outgoing", "incoming"),
  top_n = 30,
  point.size = 1.5,
  palette = "RdBu",
  palcolor = NULL,
  title = NULL
) {
  plot_type <- match.arg(plot_type)
  direction <- match.arg(direction)
  stored <- spatialcellchat_get_stored_sample(object, result.name, sample)
  result <- stored$result
  table <- result$interactions
  coords <- result$coordinates
  if (identical(plot_type, "pathway") && is.null(signaling)) {
    log_message("{.arg signaling} is required for a pathway spatial plot", message_type = "error")
  }
  if (identical(plot_type, "lr_spatial") && is.null(pairLR.use)) {
    log_message("{.arg pairLR.use} is required for an LR spatial plot", message_type = "error")
  }
  if (!is.null(signaling)) table <- table[table$pathway_name %in% signaling, , drop = FALSE]
  if (!is.null(pairLR.use)) table <- table[table$interaction_name %in% pairLR.use, , drop = FALSE]
  if (nrow(table) == 0L) log_message("No stored SpatialCellChat records match the requested filters", message_type = "error")
  colors <- thisplot::palette_colors(
    x = seq_len(100),
    n = 100,
    palette = palette,
    palcolor = palcolor,
    type = "continuous"
  )

  if (identical(plot_type, "spatial_network")) {
    network <- spatialcellchat_extract_network(table)
    network <- network[order(network$weight, decreasing = TRUE), , drop = FALSE]
    network <- utils::head(network, top_n)
    centers <- stats::aggregate(
      data.frame(x = coords$x_display, y = coords$y_display),
      by = list(label = coords$label),
      FUN = mean
    )
    from <- centers[match(network$sender, centers$label), , drop = FALSE]
    to <- centers[match(network$receiver, centers$label), , drop = FALSE]
    edges <- cbind(network, x = from$x, y = from$y, xend = to$x, yend = to$y)
    return(
      ggplot2::ggplot() +
        ggplot2::geom_segment(
          data = edges,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend, linewidth = weight),
          alpha = 0.55,
          arrow = grid::arrow(length = grid::unit(0.08, "inches"))
        ) +
        ggplot2::geom_point(
          data = centers,
          ggplot2::aes(x = x, y = y, color = label),
          size = point.size * 2
        ) +
        ggplot2::coord_fixed() +
        ggplot2::labs(
          title = title %||% "Aggregated spatial communication network",
          subtitle = paste(stored$sample, result$diagnostics$interpretation),
          color = "Group",
          linewidth = "Weight"
        ) +
        theme_scop()
    )
  }

  if (plot_type %in% c("incoming", "outgoing")) direction <- plot_type
  groups <- if (identical(direction, "outgoing")) table$sender else table$receiver
  values <- tapply(table$score, groups, sum, na.rm = TRUE)
  coords$communication_score <- as.numeric(values[as.character(coords$label)])
  if (all(!is.finite(coords$communication_score))) {
    log_message("No finite communication scores are available for plotting", message_type = "error")
  }
  ggplot2::ggplot(
    coords,
    ggplot2::aes(
      x = .data$x_display,
      y = .data$y_display,
      color = .data$communication_score
    )
  ) +
    ggplot2::geom_point(size = point.size) +
    ggplot2::scale_color_gradientn(colours = colors, na.value = "grey85") +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = title %||% paste("SpatialCellChat", direction, "score"),
      subtitle = paste(stored$sample, result$diagnostics$interpretation),
      color = "Score"
    ) +
    theme_scop()
}
