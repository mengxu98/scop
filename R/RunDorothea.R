#' @title Run DoRothEA transcription factor activity inference
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param layer Assay layer used as the expression matrix.
#' @param species Species used to select bundled DoRothEA regulons.
#' @param confidence DoRothEA confidence levels to keep.
#' @param regulons Optional regulon table with `tf`, `target`, `mor`, and
#' `confidence` columns. If `NULL`, bundled `dorothea_hs` or `dorothea_mm`
#' data are loaded from the `dorothea` package.
#' @param method Activity inference backend from `decoupleR`.
#' @param minsize Minimum regulon size passed to `decoupleR`.
#' @param options Additional named options passed to the selected `decoupleR`
#' function.
#' @param assay_name Name of the assay used to store TF activity scores.
#' @param new_assay Whether to store TF activity scores as a new assay.
#'
#' @return A `Seurat` object with DoRothEA results stored in
#' `srt@tools[["Dorothea"]]`.
#' @export
RunDorothea <- function(
  srt,
  assay = NULL,
  layer = "data",
  species = c("Homo_sapiens", "Mus_musculus"),
  confidence = c("A", "B", "C"),
  regulons = NULL,
  method = c("ulm", "viper", "wmean"),
  minsize = 5,
  options = list(),
  assay_name = "dorothea",
  new_assay = TRUE,
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  method <- match.arg(method)
  species <- match.arg(species)
  assay <- assay %||% SeuratObject::DefaultAssay(srt)

  check_r("decoupleR", verbose = FALSE)
  if (is.null(regulons)) {
    check_r("dorothea", verbose = FALSE)
    regulons <- dorothea_load_regulons(
      species = species,
      confidence = confidence
    )
  } else {
    regulons <- as.data.frame(regulons)
    regulons <- dorothea_filter_regulons(regulons, confidence = confidence)
  }
  dorothea_validate_regulons(regulons)

  expr <- GetAssayData5(srt, layer = layer, assay = assay)
  expr <- as.matrix(expr)
  expr <- expr[Matrix::rowSums(expr != 0) > 0, , drop = FALSE]
  if (nrow(expr) == 0 || ncol(expr) == 0) {
    log_message(
      "No expression values available for DoRothEA activity inference",
      message_type = "error"
    )
  }

  log_message(
    "Run {.pkg DoRothEA}/{.pkg decoupleR} with {.val {nrow(regulons)}} regulon edges",
    verbose = verbose
  )

  run_fun <- switch(
    method,
    ulm = getExportedValue("decoupleR", "run_ulm"),
    viper = getExportedValue("decoupleR", "run_viper"),
    wmean = getExportedValue("decoupleR", "run_wmean")
  )
  params <- c(
    list(
      mat = expr,
      network = regulons,
      .source = "tf",
      .target = "target",
      .mor = "mor",
      minsize = minsize
    ),
    options
  )
  res <- do.call(run_fun, params)
  scores <- dorothea_scores_to_matrix(res)
  missing_cells <- setdiff(colnames(srt), colnames(scores))
  if (length(missing_cells) > 0L) {
    log_message(
      "{.pkg decoupleR} did not return scores for all cells in {.arg srt}",
      message_type = "error"
    )
  }
  scores <- scores[, colnames(srt), drop = FALSE]

  if (isTRUE(new_assay)) {
    srt[[assay_name]] <- Seurat::CreateAssayObject(data = scores)
    srt[[assay_name]] <- Seurat::AddMetaData(
      object = srt[[assay_name]],
      metadata = data.frame(
        termnames = rownames(scores),
        row.names = rownames(scores),
        stringsAsFactors = FALSE
      )
    )
    log_message(
      "{.pkg DoRothEA} TF activity scores stored in assay {.val {assay_name}}",
      verbose = verbose
    )
  } else {
    meta_scores <- as.data.frame(t(scores), check.names = FALSE)
    colnames(meta_scores) <- make.names(
      paste(assay_name, colnames(meta_scores), sep = "_")
    )
    srt <- Seurat::AddMetaData(srt, metadata = meta_scores)
    log_message(
      "{.pkg DoRothEA} TF activity scores stored in {.cls Seurat} metadata",
      verbose = verbose
    )
  }

  srt@tools[["Dorothea"]] <- list(
    scores = scores,
    result = as.data.frame(res),
    regulon_summary = dorothea_regulon_summary(regulons),
    parameters = list(
      assay = assay,
      layer = layer,
      species = species,
      confidence = confidence,
      method = method,
      minsize = minsize,
      assay_name = assay_name,
      new_assay = new_assay,
      options = options
    )
  )
  srt
}

dorothea_load_regulons <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  confidence = c("A", "B", "C")
) {
  species <- match.arg(species)
  data_name <- switch(
    species,
    Homo_sapiens = "dorothea_hs",
    Mus_musculus = "dorothea_mm"
  )
  env <- new.env(parent = emptyenv())
  utils::data(list = data_name, package = "dorothea", envir = env)
  regulons <- get(data_name, envir = env)
  dorothea_filter_regulons(regulons, confidence = confidence)
}

dorothea_filter_regulons <- function(regulons, confidence = c("A", "B", "C")) {
  regulons <- as.data.frame(regulons)
  if (!is.null(confidence) && "confidence" %in% colnames(regulons)) {
    regulons <- regulons[regulons[["confidence"]] %in% confidence, , drop = FALSE]
  }
  regulons
}

dorothea_validate_regulons <- function(regulons) {
  required <- c("tf", "target", "mor")
  missing <- setdiff(required, colnames(regulons))
  if (length(missing) > 0) {
    log_message(
      "{.arg regulons} must contain columns: {.val {required}}",
      message_type = "error"
    )
  }
  if (nrow(regulons) == 0) {
    log_message(
      "No DoRothEA regulon edges remain after filtering",
      message_type = "error"
    )
  }
}

dorothea_scores_to_matrix <- function(res) {
  res <- as.data.frame(res)
  source_col <- intersect(c("source", "tf"), colnames(res))[1]
  condition_col <- intersect(c("condition", "sample", "cell"), colnames(res))[1]
  score_col <- intersect(c("score", "activity", "nes"), colnames(res))[1]
  if (any(is.na(c(source_col, condition_col, score_col)))) {
    log_message(
      "Unable to parse {.pkg decoupleR} result columns for DoRothEA scores",
      message_type = "error"
    )
  }
  sources <- unique(as.character(res[[source_col]]))
  conditions <- unique(as.character(res[[condition_col]]))
  scores <- matrix(
    NA_real_,
    nrow = length(sources),
    ncol = length(conditions),
    dimnames = list(sources, conditions)
  )
  idx <- cbind(
    match(as.character(res[[source_col]]), sources),
    match(as.character(res[[condition_col]]), conditions)
  )
  scores[idx] <- as.numeric(res[[score_col]])
  scores
}

dorothea_regulon_summary <- function(regulons) {
  data.frame(
    n_tfs = length(unique(regulons[["tf"]])),
    n_targets = length(unique(regulons[["target"]])),
    n_edges = nrow(regulons),
    confidence = if ("confidence" %in% colnames(regulons)) {
      paste(sort(unique(regulons[["confidence"]])), collapse = ",")
    } else {
      NA_character_
    },
    stringsAsFactors = FALSE
  )
}
