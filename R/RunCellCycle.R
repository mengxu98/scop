#' @title Run cell cycle scoring
#'
#' @description
#' Estimate cell cycle state with Seurat gene-set scoring, [scran::cyclone()],
#' or [tricycle](https://bioconductor.org/packages/tricycle).
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams CycGenePrefetch
#' @inheritParams thisutils::log_message
#' @param method Cell cycle estimation method. One of `"Seurat"`, `"cyclone"`,
#' or `"tricycle"`.
#' @param layer Data layer used by `cyclone` and `tricycle`. Default is `"counts"`.
#' @param name Prefix for metadata columns and tricycle reduction names.
#' Default is `"CellCycle"`.
#' @param phase_col Optional metadata column used to store the final phase call,
#' for example `"Phase"`. Default is `NULL`, which avoids writing a compatibility
#' phase column.
#' @param overwrite Whether to overwrite existing output columns. Default is `FALSE`.
#' @param ... Additional arguments passed to the selected method.
#'
#' @return A `Seurat` object with cell cycle metadata and, for `tricycle`,
#' a tricycle embedding reduction.
#'
#' @export
#'
#' @examples
#' data(pancreas_sub)
#' srt <- pancreas_sub[, 1:80]
#'
#' srt <- RunCellCycle(
#'   srt,
#'   method = "cyclone",
#'   species = "Mus_musculus",
#'   name = "Cyclone"
#' )
#'
#' srt <- RunCellCycle(
#'   srt,
#'   method = "tricycle",
#'   species = "Mus_musculus",
#'   name = "Tricycle"
#' )
#' if ("Cyclone_cyclone_Phase" %in% colnames(srt@meta.data)) {
#'     CellDimPlot(
#'       srt,
#'       reduction = "Tricycle_tricycleEmbedding",
#'       group.by = "Cyclone_cyclone_Phase"
#'     )
#'   }
#'   FeatureDimPlot(
#'     srt,
#'     reduction = "Tricycle_tricycleEmbedding",
#'     features = "Tricycle_tricyclePosition"
#'   )
RunCellCycle <- function(
  srt,
  method = c("Seurat", "cyclone", "tricycle"),
  assay = NULL,
  layer = "counts",
  species = "Homo_sapiens",
  name = "CellCycle",
  phase_col = NULL,
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  log_message(
    "Start cell cycle scoring",
    verbose = verbose
  )
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} is not a {.cls Seurat}",
      message_type = "error"
    )
  }
  method <- match.arg(method)
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message(
      "{.arg assay} must be one of {.val {SeuratObject::Assays(srt)}}",
      message_type = "error"
    )
  }
  if (!is.character(name) || length(name) != 1 || is.na(name) || name == "") {
    log_message(
      "{.arg name} must be a non-empty character scalar",
      message_type = "error"
    )
  }

  srt <- switch(method,
    Seurat = RunCellCycleSeurat(
      srt = srt,
      assay = assay,
      species = species,
      name = name,
      phase_col = phase_col,
      overwrite = overwrite,
      verbose = verbose,
      ...
    ),
    cyclone = RunCellCycleCyclone(
      srt = srt,
      assay = assay,
      layer = layer,
      species = species,
      name = name,
      phase_col = phase_col,
      overwrite = overwrite,
      verbose = verbose,
      ...
    ),
    tricycle = RunCellCycleTricycle(
      srt = srt,
      assay = assay,
      layer = layer,
      species = species,
      name = name,
      overwrite = overwrite,
      verbose = verbose,
      ...
    )
  )

  log_message(
    "Cell cycle scoring completed",
    message_type = "success",
    verbose = verbose
  )
  return(srt)
}

RunCellCycleSeurat <- function(
  srt,
  assay,
  species,
  name,
  phase_col = NULL,
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  output_cols <- c(
    paste0(name, "_S_score"),
    paste0(name, "_G2M_score"),
    paste0(name, "_Phase")
  )
  cellcycle_check_metadata(srt, c(output_cols, phase_col), overwrite = overwrite)

  ccgenes <- CycGenePrefetch(
    species = species,
    verbose = verbose
  )
  s_features <- intersect(ccgenes[["S"]], rownames(srt[[assay]]))
  g2m_features <- intersect(ccgenes[["G2M"]], rownames(srt[[assay]]))
  if (length(s_features) == 0 || length(g2m_features) == 0) {
    log_message(
      "No matched S or G2M cell cycle genes were found in {.arg srt}",
      message_type = "error"
    )
  }
  status <- suppressWarnings(
    CheckDataType(
      srt,
      layer = "data",
      assay = assay,
      verbose = FALSE
    )
  )
  if (status %in% c("raw_counts", "raw_normalized_counts", "unknown")) {
    log_message(
      "Perform {.fn NormalizeData} before cell cycle scoring",
      verbose = verbose
    )
    srt <- NormalizeData(
      object = srt,
      assay = assay,
      normalization.method = "LogNormalize",
      verbose = FALSE
    )
  }
  assay_raw <- SeuratObject::DefaultAssay(srt)
  SeuratObject::DefaultAssay(srt) <- assay
  srt_tmp <- Seurat::CellCycleScoring(
    object = srt,
    s.features = s_features,
    g2m.features = g2m_features,
    set.ident = FALSE,
    ...
  )
  SeuratObject::DefaultAssay(srt) <- assay_raw
  metadata <- data.frame(
    CellCycle_S_score = srt_tmp@meta.data[["S.Score"]],
    CellCycle_G2M_score = srt_tmp@meta.data[["G2M.Score"]],
    CellCycle_Phase = srt_tmp@meta.data[["Phase"]],
    row.names = colnames(srt)
  )
  colnames(metadata) <- output_cols
  srt <- cellcycle_add_metadata(
    srt = srt,
    metadata = metadata,
    phase_col = phase_col,
    phase = metadata[[paste0(name, "_Phase")]],
    overwrite = overwrite
  )
  return(srt)
}

RunCellCycleCyclone <- function(
  srt,
  assay,
  layer,
  species,
  name,
  phase_col = NULL,
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  check_r("scran", verbose = FALSE)
  output_cols <- c(
    paste0(name, "_cyclone_G1_score"),
    paste0(name, "_cyclone_S_score"),
    paste0(name, "_cyclone_G2M_score"),
    paste0(name, "_cyclone_Phase")
  )
  cellcycle_check_metadata(srt, c(output_cols, phase_col), overwrite = overwrite)

  pairs <- cellcycle_cyclone_pairs(
    species = species,
    verbose = verbose
  )
  expr <- GetAssayData5(srt, assay = assay, layer = layer)
  gene_names <- cellcycle_cyclone_gene_names(
    species = species,
    features = rownames(expr),
    pairs = pairs,
    verbose = verbose
  )
  assigned <- scran::cyclone(
    x = expr,
    pairs = pairs,
    gene.names = gene_names,
    ...
  )
  scores <- as.data.frame(assigned[["scores"]])
  if (is.null(rownames(scores))) {
    rownames(scores) <- colnames(srt)
  }
  scores <- scores[colnames(srt), , drop = FALSE]
  phases <- assigned[["phases"]]
  if (is.null(names(phases))) {
    names(phases) <- colnames(srt)
  }
  metadata <- data.frame(
    CellCycle_cyclone_G1_score = cellcycle_score_or_na(scores, "G1"),
    CellCycle_cyclone_S_score = cellcycle_score_or_na(scores, "S"),
    CellCycle_cyclone_G2M_score = cellcycle_score_or_na(scores, "G2M"),
    CellCycle_cyclone_Phase = phases[colnames(srt)],
    row.names = colnames(srt)
  )
  colnames(metadata) <- output_cols
  srt <- cellcycle_add_metadata(
    srt = srt,
    metadata = metadata,
    phase_col = phase_col,
    phase = metadata[[paste0(name, "_cyclone_Phase")]],
    overwrite = overwrite
  )
  return(srt)
}

RunCellCycleTricycle <- function(
  srt,
  assay,
  layer,
  species,
  name,
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  check_r(c("tricycle", "scuttle"), verbose = FALSE)
  output_cols <- paste0(name, "_tricyclePosition")
  cellcycle_check_metadata(srt, output_cols, overwrite = overwrite)
  reduction_name <- paste0(name, "_tricycleEmbedding")
  if (reduction_name %in% names(srt@reductions) && !isTRUE(overwrite)) {
    log_message(
      "Reduction {.val {reduction_name}} already exists. Set {.arg overwrite = TRUE} to replace it",
      message_type = "error"
    )
  }

  sce <- Seurat::as.SingleCellExperiment(srt, assay = assay)
  if (layer == "counts") {
    if (!"logcounts" %in% SummarizedExperiment::assayNames(sce)) {
      sce <- scuttle::logNormCounts(sce)
    }
  } else {
    SummarizedExperiment::assay(sce, "logcounts") <- GetAssayData5(
      srt,
      assay = assay,
      layer = layer
    )
  }
  project_args <- list(...)
  project_args$x <- sce
  project_args$species <- project_args$species %||% cellcycle_tricycle_species(species)
  project_args$gname.type <- project_args$gname.type %||% "SYMBOL"
  sce <- do.call(tricycle::project_cycle_space, project_args)
  sce <- tricycle::estimate_cycle_position(sce)

  position <- SummarizedExperiment::colData(sce)[["tricyclePosition"]]
  names(position) <- colnames(srt)
  metadata <- data.frame(
    CellCycle_tricyclePosition = position,
    row.names = colnames(srt)
  )
  colnames(metadata) <- output_cols
  srt <- cellcycle_add_metadata(
    srt = srt,
    metadata = metadata,
    phase_col = NULL,
    phase = NULL,
    overwrite = overwrite
  )

  if (!"tricycleEmbedding" %in% SingleCellExperiment::reducedDimNames(sce)) {
    log_message(
      "{.pkg tricycle} did not return {.val tricycleEmbedding}",
      message_type = "error"
    )
  }
  embeddings <- SingleCellExperiment::reducedDim(sce, "tricycleEmbedding")
  rownames(embeddings) <- colnames(srt)
  reduction_key <- cellcycle_reduction_key(name)
  colnames(embeddings) <- paste0(reduction_key, seq_len(ncol(embeddings)))
  srt[[reduction_name]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings,
    assay = assay,
    key = reduction_key,
    misc = list(
      method = "tricycle",
      species = species
    )
  )
  return(srt)
}

cellcycle_check_metadata <- function(srt, cols, overwrite = FALSE) {
  cols <- cols[!is.na(cols)]
  cols <- unique(cols[cols != ""])
  existing <- intersect(cols, colnames(srt@meta.data))
  if (length(existing) > 0 && !isTRUE(overwrite)) {
    log_message(
      "Metadata column{?s} already exist{?s}: {.val {existing}}. Set {.arg overwrite = TRUE} to replace existing result{?s}",
      message_type = "error"
    )
  }
}

cellcycle_add_metadata <- function(
  srt,
  metadata,
  phase_col = NULL,
  phase = NULL,
  overwrite = FALSE
) {
  if (isTRUE(overwrite)) {
    for (col in intersect(colnames(metadata), colnames(srt@meta.data))) {
      srt[[col]] <- NULL
    }
    if (!is.null(phase_col) && phase_col %in% colnames(srt@meta.data)) {
      srt[[phase_col]] <- NULL
    }
  }
  srt <- Seurat::AddMetaData(
    object = srt,
    metadata = metadata
  )
  if (!is.null(phase_col) && !phase_col %in% colnames(metadata)) {
    srt[[phase_col]] <- phase
  }
  return(srt)
}

cellcycle_cyclone_pairs <- function(
  species,
  verbose = TRUE
) {
  marker_file <- switch(species,
    Homo_sapiens = "human_cycle_markers.rds",
    Mus_musculus = "mouse_cycle_markers.rds",
    human = "human_cycle_markers.rds",
    mouse = "mouse_cycle_markers.rds",
    NULL
  )
  if (is.null(marker_file)) {
    log_message(
      "{.arg species} must be one of {.val Homo_sapiens}, {.val Mus_musculus}, {.val human}, or {.val mouse} for {.pkg scran::cyclone}",
      message_type = "error"
    )
  }
  marker_path <- system.file("exdata", marker_file, package = "scran")
  if (marker_path == "") {
    log_message(
      "Cannot find {.pkg scran} cyclone marker file {.file {marker_file}}",
      message_type = "error"
    )
  }
  readRDS(marker_path)
}

cellcycle_cyclone_gene_names <- function(
  species,
  features,
  pairs,
  verbose = TRUE
) {
  marker_genes <- unique(unlist(pairs, use.names = FALSE))
  if (length(intersect(marker_genes, features)) > 0) {
    return(features)
  }
  org_pkg <- switch(species,
    Homo_sapiens = "org.Hs.eg.db",
    human = "org.Hs.eg.db",
    Mus_musculus = "org.Mm.eg.db",
    mouse = "org.Mm.eg.db",
    NULL
  )
  if (is.null(org_pkg)) {
    log_message(
      "{.arg species} must be one of {.val Homo_sapiens}, {.val Mus_musculus}, {.val human}, or {.val mouse} for gene ID mapping in {.pkg scran::cyclone}",
      message_type = "error"
    )
  }
  check_r(c("AnnotationDbi", org_pkg), verbose = FALSE)
  orgdb <- get_namespace_fun(org_pkg, org_pkg)
  gene_map <- AnnotationDbi::select(
    x = orgdb,
    keys = features,
    keytype = "SYMBOL",
    columns = "ENSEMBL"
  )
  gene_map <- gene_map[!is.na(gene_map[["ENSEMBL"]]), , drop = FALSE]
  gene_map <- gene_map[!duplicated(gene_map[["SYMBOL"]]), , drop = FALSE]
  ensembl <- stats::setNames(gene_map[["ENSEMBL"]], gene_map[["SYMBOL"]])
  gene_names <- ensembl[features]
  gene_names[is.na(gene_names)] <- features[is.na(gene_names)]
  if (length(intersect(marker_genes, gene_names)) == 0) {
    log_message(
      "Unable to map input feature names to {.pkg scran::cyclone} marker IDs with {.pkg {org_pkg}}",
      message_type = "error"
    )
  }
  log_message(
    "Map input feature names to ENSEMBL IDs with {.pkg {org_pkg}} for {.pkg scran::cyclone}",
    verbose = verbose
  )
  return(gene_names)
}

cellcycle_tricycle_species <- function(species) {
  species_use <- switch(species,
    Homo_sapiens = "human",
    Mus_musculus = "mouse",
    human = "human",
    mouse = "mouse",
    NULL
  )
  if (is.null(species_use)) {
    log_message(
      "{.arg species} must be one of {.val Homo_sapiens}, {.val Mus_musculus}, {.val human}, or {.val mouse} for {.pkg tricycle}",
      message_type = "error"
    )
  }
  return(species_use)
}

cellcycle_score_or_na <- function(scores, phase) {
  if (phase %in% colnames(scores)) {
    return(scores[[phase]])
  }
  stats::setNames(rep(NA_real_, nrow(scores)), rownames(scores))
}

cellcycle_reduction_key <- function(name) {
  key <- gsub("[^A-Za-z0-9]", "", name)
  if (key == "") {
    key <- "CellCycle"
  }
  paste0(key, "TC_")
}
