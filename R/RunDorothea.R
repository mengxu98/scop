#' @title Run transcription factor activity inference
#'
#' @md
#' @inheritParams RunStandardWorkflow
#' @inheritParams thisutils::log_message
#' @param layer Assay layer used as the expression matrix.
#' @param species Species used to select the regulatory network. The bundled
#' DoRothEA and CollecTRI networks support human and mouse. For other input
#' species, set `input_species` and project expression values to this network
#' species through homologous gene conversion before activity inference.
#' @param input_species Species of the input expression features. If `NULL`,
#' the input is assumed to use the same gene namespace as `species`. When this
#' differs from `species`, expression features are converted with
#' [ConvertHomologs] before DoRothEA activity inference.
#' @param geneID_from_IDtype,geneID_to_IDtype Gene identifier types passed to
#' [ConvertHomologs] for cross-species projection. For bundled DoRothEA
#' regulons, `geneID_to_IDtype` should normally remain `"symbol"`.
#' @param homolog_params Additional named arguments passed to [ConvertHomologs]
#' when `input_species` differs from `species`, such as `Ensembl_version`,
#' `biomart`, `mirror`, `max_tries`, `multi_mapping`, and `collapse_fun`.
#' @param confidence DoRothEA confidence levels to keep when a `confidence`
#' column is present.
#' @param regulons Regulon source. `NULL` or `"dorothea"` loads the bundled
#' DoRothEA network, `"collectri"` loads CollecTRI, and a data-frame-like
#' object supplies a custom network. Custom networks may use `tf`, `TF`,
#' `source`, or `regulator` for the regulator column, `target` for targets,
#' and `mor`, `importance`, or `weight` for edge weights. `mor` is treated as
#' signed regulation; `importance` and `weight` are treated as non-negative
#' unsigned weights.
#' @param method Activity inference backend from `decoupleR`.
#' @param minsize Minimum regulon size passed to `decoupleR`.
#' @param options Additional named options passed to the selected `decoupleR`
#' function.
#' @param assay_name Name of the assay used to store TF activity scores.
#' @param new_assay Whether to store TF activity scores as a new assay.
#' @param add_meta Whether to also write TF activity scores to `srt@meta.data`
#' with the `assay_name` prefix for direct plotting with [FeatureDimPlot()].
#' @details
#' `RunGRN()` returns a `TF`, `target`, and `importance` table that can be
#' supplied directly through `regulons`. Such a network is unsigned: its
#' scores describe target-program activity and do not distinguish activation
#' from repression. For signed TF activity, supply a network with a `mor`
#' column instead.
#'
#' @return A `Seurat` object with TF activity results stored in
#' `srt@tools[["Dorothea"]]`, optionally TF activity scores stored in
#' `srt@meta.data`, and optionally a TF activity assay when
#' `new_assay = TRUE`. The normalized network is stored in
#' `srt@tools[["Dorothea"]]$regulons`, with the source and signedness recorded
#' in `network_info`. For cross-species runs, the homolog projection summary
#' is stored in `srt@tools[["Dorothea"]]$homolog_conversion`.
#' @export
#'
#' @references
#' Garcia-Alonso, L., Holland, C.H., Ibrahim, M.M., Turei, D.,
#' and Saez-Rodriguez, J. (2019). Benchmark and integration of
#' resources for the estimation of human transcription factor activities.
#' \emph{Genome Research}, 29, 1363-1375. \doi{10.1101/gr.240663.118}
#'
#' Badia-i-Mompel, P., Velez Santiago, J., Braunger, J., Geiss, C.,
#' Dimitrov, D., Muller-Dott, S., Taus, P., Dugourd, A., Holland, C.H.,
#' Ramirez Flores, R.O., and Saez-Rodriguez, J. (2022). decoupleR:
#' ensemble of computational methods to infer biological activities from
#' omics data. \emph{Bioinformatics Advances}, 2, vbac016.
#' \doi{10.1093/bioadv/vbac016}
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- RunStandardWorkflow(
#'   pancreas_sub,
#'   verbose = FALSE
#' )
#'
#' pancreas_sub <- RunDorothea(
#'   pancreas_sub,
#'   layer = "counts",
#'   species = "Mus_musculus",
#'   confidence = c("A", "B", "C"),
#'   method = "ulm",
#'   minsize = 5
#' )
#'
#' pancreas_sub@tools$Dorothea$regulon_summary
#' head(pancreas_sub@tools$Dorothea$result)
#'
#' tf_use <- intersect(
#'   c("Sox9", "Neurod1", "Pdx1", "Foxa2"),
#'   rownames(pancreas_sub@tools$Dorothea$scores)
#' )
#' FeatureDimPlot(
#'   pancreas_sub,
#'   assay = "dorothea",
#'   features = tf_use,
#'   ncol = 2
#' )
#' FeatureStatPlot(
#'   pancreas_sub,
#'   assay = "dorothea",
#'   stat.by = tf_use,
#'   group.by = "CellType",
#'   plot_type = "violin"
#' )
#'
#' DorotheaPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   features = "Sox9",
#'   plot_type = "dim"
#' )
#' DorotheaPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   group1 = "Endocrine",
#'   group2 = "Ductal",
#'   plot_type = "bar",
#'   top_n = 20
#' )
#' ht <- DorotheaPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   plot_type = "heatmap",
#'   top_n = 20
#' )
#' ht$plot
#' DorotheaPlot(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   group1 = "Endocrine",
#'   group2 = "Ductal",
#'   features = "Sox9",
#'   plot_type = "targets"
#' )
#'
#' # A RunGRN-compatible unsigned network can be passed directly.
#' grn_targets <- head(rownames(pancreas_sub), 5)
#' grn_edges <- data.frame(
#'   TF = rep("ExampleTF", length(grn_targets)),
#'   target = grn_targets,
#'   importance = seq_along(grn_targets) / length(grn_targets)
#' )
#' pancreas_sub <- RunDorothea(
#'   pancreas_sub,
#'   regulons = grn_edges,
#'   species = "Mus_musculus",
#'   minsize = 1,
#'   new_assay = FALSE,
#'   add_meta = FALSE,
#'   verbose = FALSE
#' )
#' pancreas_sub@tools$Dorothea$network_info
#' # A real RunGRN result can be connected in the same way:
#' # grn <- RunGRN(pancreas_sub, grn_method = "genie3")
#' # pancreas_sub <- RunDorothea(pancreas_sub, regulons = grn)
#' \dontrun{
#' pancreas_sub <- RunDorothea(
#'   pancreas_sub,
#'   regulons = "collectri",
#'   species = "Mus_musculus",
#'   minsize = 5
#' )
#' }
RunDorothea <- function(
  srt,
  assay = NULL,
  layer = "data",
  species = c("Homo_sapiens", "Mus_musculus"),
  input_species = NULL,
  geneID_from_IDtype = "symbol",
  geneID_to_IDtype = "symbol",
  homolog_params = list(),
  confidence = c("A", "B", "C"),
  regulons = NULL,
  method = c("ulm", "viper", "wmean"),
  minsize = 5,
  options = list(),
  assay_name = "dorothea",
  new_assay = TRUE,
  add_meta = TRUE,
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
  input_species <- input_species %||% species
  if (
    !is.character(input_species) ||
      length(input_species) != 1L ||
      is.na(input_species)
  ) {
    log_message(
      "{.arg input_species} must be a single species name",
      message_type = "error"
    )
  }
  assay <- assay %||% SeuratObject::DefaultAssay(srt)

  check_r("decoupleR", verbose = FALSE)
  regulons_input <- NULL
  if (is.null(regulons)) {
    network_name <- "dorothea"
  } else if (is.character(regulons) && length(regulons) == 1L) {
    network_name <- tolower(regulons)
    if (!network_name %in% c("dorothea", "collectri")) {
      log_message(
        "{.arg regulons} must be NULL, {.val 'dorothea'}, {.val 'collectri'}, or a network table",
        message_type = "error"
      )
    }
  } else {
    network_name <- "custom"
    regulons_input <- as.data.frame(regulons, stringsAsFactors = FALSE)
  }

  if (identical(network_name, "dorothea")) {
    check_r("dorothea", verbose = FALSE)
    data_name <- switch(species,
      Homo_sapiens = "dorothea_hs",
      Mus_musculus = "dorothea_mm"
    )
    env <- new.env(parent = emptyenv())
    utils::data(list = data_name, package = "dorothea", envir = env)
    regulons <- get(data_name, envir = env)
  } else if (identical(network_name, "collectri")) {
    organism <- switch(species,
      Homo_sapiens = "human",
      Mus_musculus = "mouse"
    )
    get_collectri <- get_namespace_fun("decoupleR", "get_collectri")
    if (!is.function(get_collectri)) {
      log_message(
        "{.pkg decoupleR} does not provide {.fun get_collectri}",
        message_type = "error"
      )
    }
    regulons <- get_collectri(
      organism = organism,
      split_complexes = FALSE
    )
  } else {
    regulons <- regulons_input
  }

  if (!is.null(confidence) && "confidence" %in% colnames(regulons)) {
    regulons <- regulons[
      regulons[["confidence"]] %in% confidence, ,
      drop = FALSE
    ]
  }
  network <- dorothea_normalize_network(
    regulons,
    source = network_name
  )
  regulons <- network$regulons
  network_info <- network$info
  network_info$species <- species
  network_info$input_species <- input_species
  if (identical(network_name, "collectri")) {
    network_info$split_complexes <- FALSE
  }

  expr <- GetAssayData5(srt, layer = layer, assay = assay)
  expr <- expr[Matrix::rowSums(expr != 0) > 0, , drop = FALSE]
  homolog_conversion <- NULL
  if (!identical(input_species, species)) {
    invalid_homolog_params <- length(homolog_params) > 0L &&
      (is.null(names(homolog_params)) || any(!nzchar(names(homolog_params))))
    if (!is.list(homolog_params) || invalid_homolog_params) {
      log_message(
        "{.arg homolog_params} must be a named list",
        message_type = "error"
      )
    }
    reserved <- c(
      "object",
      "species_from",
      "species_to",
      "geneID_from_IDtype",
      "geneID_to_IDtype",
      "keep_unmapped",
      "verbose"
    )
    duplicated_params <- intersect(names(homolog_params), reserved)
    if (length(duplicated_params) > 0L) {
      log_message(
        "{.arg homolog_params} cannot override {.val {duplicated_params}}",
        message_type = "error"
      )
    }
    log_message(
      "Project expression features from {.val {input_species}} to {.val {species}} homologs for {.val {network_info$label}}",
      verbose = verbose
    )
    expr <- do.call(
      ConvertHomologs.default,
      c(
        list(
          object = expr,
          species_from = input_species,
          species_to = species,
          geneID_from_IDtype = geneID_from_IDtype,
          geneID_to_IDtype = geneID_to_IDtype,
          keep_unmapped = FALSE,
          verbose = verbose
        ),
        homolog_params
      )
    )
    conversion <- attr(expr, "ConvertHomologs")
    if (!is.null(conversion)) {
      homolog_conversion <- data.frame(
        species_from = conversion$species_from,
        species_to = conversion$species_to,
        geneID_from_IDtype = paste(
          conversion$geneID_from_IDtype,
          collapse = ","
        ),
        geneID_to_IDtype = paste(
          conversion$geneID_to_IDtype,
          collapse = ","
        ),
        n_mapped_source_genes = length(unique(conversion$mapping$from_geneID)),
        n_target_homologs = length(unique(conversion$mapping$to_geneID)),
        n_unmapped_source_genes = length(conversion$unmapped),
        Ensembl_version = conversion$Ensembl_version %||% NA_character_,
        stringsAsFactors = FALSE
      )
    }
    expr <- expr[Matrix::rowSums(expr != 0) > 0, , drop = FALSE]
  }
  expr <- as.matrix(expr)
  if (nrow(expr) == 0 || ncol(expr) == 0) {
    log_message(
      "No expression values available for {.val {network_info$label}} activity inference",
      message_type = "error"
    )
  }

  log_message(
    "Run {.val {network_info$label}}/{.pkg decoupleR} with {.val {nrow(regulons)}} regulon edges",
    verbose = verbose
  )

  run_fun <- dorothea_get_run_fun(method)
  if (!is.function(run_fun)) {
    log_message(
      "The selected {.pkg decoupleR} method {.val {method}} is not available",
      message_type = "error"
    )
  }
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
  res_df <- as.data.frame(res)
  source_col <- intersect(c("source", "tf"), colnames(res_df))[1]
  condition_col <- intersect(
    c("condition", "sample", "cell"),
    colnames(res_df)
  )[1]
  score_col <- intersect(c("score", "activity", "nes"), colnames(res_df))[1]
  if (any(is.na(c(source_col, condition_col, score_col)))) {
    log_message(
      "Unable to parse {.pkg decoupleR} result columns for {.val {network_info$label}} scores",
      message_type = "error"
    )
  }
  sources <- unique(as.character(res_df[[source_col]]))
  conditions <- unique(as.character(res_df[[condition_col]]))
  scores <- matrix(
    NA_real_,
    nrow = length(sources),
    ncol = length(conditions),
    dimnames = list(sources, conditions)
  )
  idx <- cbind(
    match(as.character(res_df[[source_col]]), sources),
    match(as.character(res_df[[condition_col]]), conditions)
  )
  scores[idx] <- as.numeric(res_df[[score_col]])
  missing_cells <- setdiff(colnames(srt), colnames(scores))
  if (length(missing_cells) > 0L) {
    log_message(
      "{.pkg decoupleR} did not return scores for all cells in {.arg srt}",
      message_type = "error"
    )
  }
  scores <- scores[, colnames(srt), drop = FALSE]

  if (isTRUE(new_assay)) {
    srt <- dorothea_attach_assay(
      srt = srt,
      scores = scores,
      assay_name = assay_name
    )
    log_message(
      "{.val {network_info$label}} TF activity scores stored in assay {.val {assay_name}}",
      verbose = verbose
    )
  }
  if (isTRUE(add_meta)) {
    meta_scores <- as.data.frame(t(scores), check.names = FALSE)
    colnames(meta_scores) <- make.names(
      paste(assay_name, colnames(meta_scores), sep = "_")
    )
    srt <- Seurat::AddMetaData(srt, metadata = meta_scores)
    log_message(
      "{.val {network_info$label}} TF activity scores stored in {.cls Seurat} metadata",
      verbose = verbose
    )
  }

  srt@tools[["Dorothea"]] <- list(
    scores = scores,
    result = res_df,
    regulons = regulons,
    regulons_input = regulons_input,
    network_info = network_info,
    regulon_summary = data.frame(
      n_tfs = length(unique(regulons[["tf"]])),
      n_targets = length(unique(regulons[["target"]])),
      n_edges = nrow(regulons),
      confidence = if ("confidence" %in% colnames(regulons)) {
        paste(sort(unique(regulons[["confidence"]])), collapse = ",")
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    ),
    parameters = list(
      assay = assay,
      layer = layer,
      species = species,
      input_species = input_species,
      geneID_from_IDtype = geneID_from_IDtype,
      geneID_to_IDtype = geneID_to_IDtype,
      confidence = confidence,
      method = method,
      minsize = minsize,
      assay_name = assay_name,
      new_assay = new_assay,
      add_meta = add_meta,
      homolog_params = homolog_params,
      options = options
    ),
    homolog_conversion = homolog_conversion
  )
  srt
}

dorothea_get_run_fun <- function(method) {
  fun_name <- switch(method,
    ulm = "run_ulm",
    viper = "run_viper",
    wmean = "run_wmean"
  )
  get_namespace_fun("decoupleR", fun_name)
}

dorothea_normalize_network <- function(regulons, source = "custom") {
  regulons <- as.data.frame(regulons, stringsAsFactors = FALSE)
  if (nrow(regulons) == 0L) {
    log_message(
      "No {.val {source}} regulon edges remain after filtering",
      message_type = "error"
    )
  }
  if (is.null(colnames(regulons))) {
    log_message(
      "{.arg regulons} must have named regulator, target, and weight columns",
      message_type = "error"
    )
  }

  regulator_col <- dorothea_choose_network_column(
    regulons,
    c("tf", "TF", "source", "regulator"),
    "regulator"
  )
  target_col <- dorothea_choose_network_column(
    regulons,
    "target",
    "target"
  )
  weight_col <- dorothea_choose_network_column(
    regulons,
    c("mor", "importance", "weight"),
    "weight"
  )

  tf <- trimws(as.character(regulons[[regulator_col]]))
  target <- trimws(as.character(regulons[[target_col]]))
  weight <- suppressWarnings(as.numeric(regulons[[weight_col]]))
  if (anyNA(tf) || any(!nzchar(tf)) || anyNA(target) || any(!nzchar(target))) {
    log_message(
      "{.arg regulons} regulator and target columns must contain non-empty names",
      message_type = "error"
    )
  }
  if (anyNA(weight) || any(!is.finite(weight))) {
    log_message(
      "{.arg regulons} weight column must contain finite numeric values",
      message_type = "error"
    )
  }
  signed <- identical(weight_col, "mor")
  if (!signed && any(weight < 0)) {
    log_message(
      "Unsigned {.arg regulons} weights from {.val {weight_col}} must be non-negative",
      message_type = "error"
    )
  }
  if (any(weight == 0)) {
    log_message(
      "{.arg regulons} weights must be non-zero",
      message_type = "error"
    )
  }
  edge_id <- paste(tf, target, sep = "\r")
  if (anyDuplicated(edge_id)) {
    log_message(
      "{.arg regulons} must not contain duplicate regulator-target edges",
      message_type = "error"
    )
  }

  regulons$tf <- tf
  regulons$target <- target
  regulons$mor <- weight
  format <- if (identical(source, "dorothea")) {
    "dorothea"
  } else if (identical(source, "collectri")) {
    "collectri"
  } else if (identical(regulator_col, "TF") &&
      identical(weight_col, "importance")) {
    "scop_grn"
  } else if (identical(regulator_col, "source") &&
      identical(weight_col, "mor")) {
    "decoupler_network"
  } else {
    "custom_table"
  }
  label <- switch(source,
    dorothea = "DoRothEA",
    collectri = "CollecTRI",
    "custom TF network"
  )
  info <- list(
    source = source,
    label = label,
    format = format,
    signed = signed,
    regulator_column = regulator_col,
    target_column = target_col,
    weight_column = weight_col,
    weight_semantics = if (signed) {
      "signed_mor"
    } else {
      "unsigned_positive_weight"
    }
  )
  if (!signed) {
    log_message(
      "The supplied network is unsigned; {.val {weight_col}} is used as positive edge weight. Scores represent target-program activity, not signed TF activation or repression.",
      message_type = "info"
    )
  }
  list(regulons = regulons, info = info)
}

dorothea_choose_network_column <- function(regulons, candidates, role) {
  present <- candidates[candidates %in% colnames(regulons)]
  if (length(present) == 0L) {
    log_message(
      "Unable to identify {.val {role}} column in {.arg regulons}; supported columns are {.val {candidates}}",
      message_type = "error"
    )
  }
  chosen <- present[[1L]]
  if (length(present) > 1L) {
    reference <- as.character(regulons[[chosen]])
    same <- vapply(
      present[-1L],
      function(candidate) {
        value <- as.character(regulons[[candidate]])
        identical(reference, value) || all(
          ifelse(is.na(reference), "<NA>", reference) ==
            ifelse(is.na(value), "<NA>", value)
        )
      },
      logical(1)
    )
    if (any(!same)) {
      log_message(
        "Ambiguous {.val {role}} columns in {.arg regulons}: {.val {present}}",
        message_type = "error"
      )
    }
  }
  chosen
}

dorothea_attach_assay <- function(srt, scores, assay_name) {
  scores <- as.matrix(scores)
  missing_cells <- setdiff(colnames(srt), colnames(scores))
  if (length(missing_cells) > 0L) {
    pad <- matrix(
      NA_real_,
      nrow = nrow(scores),
      ncol = length(missing_cells),
      dimnames = list(rownames(scores), missing_cells)
    )
    scores <- cbind(scores, pad)
  }
  scores <- scores[, colnames(srt), drop = FALSE]
  assay_object <- suppressWarnings(Seurat::CreateAssayObject(data = scores))
  assay_object <- Seurat::AddMetaData(
    object = assay_object,
    metadata = data.frame(
      termnames = rownames(scores),
      row.names = rownames(scores),
      stringsAsFactors = FALSE
    )
  )
  suppressWarnings(srt[[assay_name]] <- assay_object)
  srt
}
