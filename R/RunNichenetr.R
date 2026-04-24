#' @title Run NicheNet analysis
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams standard_scop
#' @param group.by Metadata column defining cell types.
#' @param receiver Receiver cell type(s).
#' @param sender Sender cell type(s). Use `"all"` to use all non-receiver cell types.
#' @param condition.by Metadata column defining conditions.
#' @param condition_oi Condition of interest.
#' @param condition_reference Reference condition.
#' @param mode NicheNet wrapper mode. Supported values are `"aggregate"`,
#' `"aggregate_cluster_de"`, and `"custom"`.
#' @param assay Assay to use.
#' @param expression_pct Minimum fraction of cells expressing a gene.
#' @param geneset Optional target gene set for `mode = "custom"`.
#' @param background_expressed_genes Optional background genes for
#' `mode = "custom"`.
#' @param top_n_ligands Number of top ligands to keep in standardized summaries.
#' @param top_n_targets Number of top targets per ligand to keep in standardized summaries.
#' @param species Species for default NicheNet prior model loading.
#' @param lr_network Optional ligand-receptor prior model or path to an `.rds` file.
#' @param ligand_target_matrix Optional ligand-target prior model or path to an `.rds` file.
#' @param weighted_networks Optional weighted network list or path to an `.rds` file.
#' @param cutoff_visualization Cutoff passed to
#' `nichenetr::prepare_ligand_target_visualization()`.
#' @param lfc_cutoff Log fold change cutoff for DE analysis in `"aggregate"` and
#' `"aggregate_cluster_de"` modes. Default is `0.25`.
#' @param use_sender_agnostic_background Whether to use all sender cell types when
#' `sender = "all"`.
#' @return A Seurat object with standardized NicheNet results stored in
#' `srt@tools[["Nichenetr"]]`.
#' @export
RunNichenetr <- function(
  srt,
  group.by,
  receiver,
  sender = "all",
  condition.by = NULL,
  condition_oi = NULL,
  condition_reference = NULL,
  mode = c("aggregate", "aggregate_cluster_de", "custom"),
  assay = NULL,
  expression_pct = 0.10,
  geneset = NULL,
  background_expressed_genes = NULL,
  top_n_ligands = 30,
  top_n_targets = 200,
  species = c("Homo_sapiens", "Mus_musculus"),
  lr_network = NULL,
  ligand_target_matrix = NULL,
  weighted_networks = NULL,
  cutoff_visualization = 0.33,
  lfc_cutoff = 0.25,
  use_sender_agnostic_background = TRUE,
  verbose = TRUE
) {
  check_r("saeyslab/nichenetr", verbose = FALSE)

  mode <- match.arg(mode)
  species <- match.arg(species)
  assay <- assay %||% DefaultAssay(srt)

  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (!group.by %in% colnames(srt[[]])) {
    log_message(
      "{.arg group.by} must be a valid metadata column in {.cls Seurat}",
      message_type = "error"
    )
  }
  if (!is.null(condition.by)) {
    if (!condition.by %in% colnames(srt[[]])) {
      log_message(
        "{.arg condition.by} ({.val {condition.by}}) must be a valid metadata column in {.cls Seurat}",
        message_type = "error"
      )
    }
    condition_vals <- unique(as.character(srt[[condition.by]][, 1]))
    if (!is.null(condition_oi) && !condition_oi %in% condition_vals) {
      log_message(
        "{.arg condition_oi} ({.val {condition_oi}}) not found in column {.val {condition.by}}. Available values: {.val {condition_vals}}",
        message_type = "error"
      )
    }
    if (
      !is.null(condition_reference) && !condition_reference %in% condition_vals
    ) {
      log_message(
        "{.arg condition_reference} ({.val {condition_reference}}) not found in column {.val {condition.by}}. Available values: {.val {condition_vals}}",
        message_type = "error"
      )
    }
  }

  sender_use <- resolve_sender_groups(
    srt = srt,
    group.by = group.by,
    sender = sender,
    receiver = receiver,
    use_sender_agnostic_background = use_sender_agnostic_background
  )

  model <- load_nichenetr_models(
    species = species,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    weighted_networks = weighted_networks,
    verbose = verbose
  )

  old_idents <- Idents(srt)
  on.exit(
    try(Idents(srt) <- old_idents, silent = TRUE),
    add = TRUE
  )
  Idents(srt) <- srt[[group.by]][, 1]

  ns_fun <- function(name) getExportedValue("nichenetr", name)
  ligand_activities <- NULL
  ligand_target_df <- NULL
  ligand_receptor_df <- NULL
  raw_result <- NULL
  de_table <- NULL

  if (identical(mode, "aggregate")) {
    required <- c("condition.by", "condition_oi", "condition_reference")
    missing_req <- required[vapply(
      mget(required),
      is.null,
      logical(1)
    )]
    if (length(missing_req) > 0L) {
      log_message(
        "For {.val mode = 'aggregate'}, the following arguments are required: {.val {missing_req}}",
        message_type = "error"
      )
    }
    raw_result <- ns_fun("nichenet_seuratobj_aggregate")(
      seurat_obj = srt,
      receiver = receiver,
      condition_colname = condition.by,
      condition_oi = condition_oi,
      condition_reference = condition_reference,
      sender = sender_use,
      ligand_target_matrix = model$ligand_target_matrix,
      lr_network = model$lr_network,
      weighted_networks = model$weighted_networks,
      assay_oi = assay,
      expression_pct = expression_pct,
      lfc_cutoff = lfc_cutoff
    )
    ligand_activities <- raw_result$ligand_activities %||% NULL
    ligand_target_df <- raw_result$ligand_target_df %||% NULL
    ligand_receptor_df <- raw_result$ligand_receptor_df %||% NULL
    de_table <- raw_result$geneset_oi %||% NULL
  } else if (identical(mode, "aggregate_cluster_de")) {
    required <- c("condition.by", "condition_oi", "condition_reference")
    missing_req <- required[vapply(
      mget(required),
      is.null,
      logical(1)
    )]
    if (length(missing_req) > 0L) {
      log_message(
        "For {.val mode = 'aggregate_cluster_de'}, the following arguments are required: {.val {missing_req}}",
        message_type = "error"
      )
    }
    raw_result <- ns_fun("nichenet_seuratobj_aggregate_cluster_de")(
      seurat_obj = srt,
      receiver = receiver,
      condition_colname = condition.by,
      condition_oi = condition_oi,
      condition_reference = condition_reference,
      sender = sender_use,
      ligand_target_matrix = model$ligand_target_matrix,
      lr_network = model$lr_network,
      weighted_networks = model$weighted_networks,
      assay_oi = assay,
      expression_pct = expression_pct,
      lfc_cutoff = lfc_cutoff
    )
    ligand_activities <- raw_result$ligand_activities %||% NULL
    ligand_target_df <- raw_result$ligand_target_df %||% NULL
    ligand_receptor_df <- raw_result$ligand_receptor_df %||% NULL
    de_table <- raw_result$geneset_oi %||% NULL
  } else {
    if (is.null(geneset)) {
      log_message(
        "{.arg geneset} must be provided for {.val mode = 'custom'}",
        message_type = "error"
      )
    }
    get_expressed_genes <- ns_fun("get_expressed_genes")
    predict_ligand_activities <- ns_fun("predict_ligand_activities")
    get_weighted_ligand_target_links <- ns_fun(
      "get_weighted_ligand_target_links"
    )

    if (is.null(background_expressed_genes)) {
      background_expressed_genes <- get_expressed_genes(
        receiver,
        srt,
        pct = expression_pct,
        assay_oi = assay
      )
    }

    receiver_expr <- unique(background_expressed_genes)
    sender_expr_list <- stats::setNames(
      lapply(sender_use, function(x) {
        get_expressed_genes(x, srt, pct = expression_pct, assay_oi = assay)
      }),
      sender_use
    )
    sender_expr <- unique(unlist(sender_expr_list, use.names = FALSE))

    lr_df <- model$lr_network
    colnames(lr_df)[1:2] <- c("from", "to")
    candidate_ligands <- unique(lr_df$from[
      lr_df$from %in% sender_expr & lr_df$to %in% receiver_expr
    ])

    ligand_activities <- predict_ligand_activities(
      geneset = geneset,
      background_expressed_genes = background_expressed_genes,
      ligand_target_matrix = model$ligand_target_matrix,
      potential_ligands = candidate_ligands
    )
    ligand_activities <- ligand_activities[
      order(ligand_activities$pearson, decreasing = TRUE), ,
      drop = FALSE
    ]
    top_ligands <- utils::head(ligand_activities$test_ligand, top_n_ligands)

    ligand_target_df <- do.call(
      rbind,
      lapply(top_ligands, function(lig) {
        get_weighted_ligand_target_links(
          lig,
          geneset = geneset,
          ligand_target_matrix = model$ligand_target_matrix,
          n = top_n_targets
        )
      })
    )
    ligand_receptor_df <- lr_df[
      lr_df$from %in% top_ligands & lr_df$to %in% receiver_expr, ,
      drop = FALSE
    ]
    raw_result <- list(
      ligand_activities = ligand_activities,
      ligand_target_df = ligand_target_df,
      ligand_receptor_df = ligand_receptor_df,
      geneset_oi = geneset,
      background_expressed_genes = background_expressed_genes
    )
  }

  bundle <- standardize_nichenetr_result(
    raw_result = raw_result,
    ligand_activities = ligand_activities,
    ligand_target_df = ligand_target_df,
    ligand_receptor_df = ligand_receptor_df,
    sender_use = sender_use,
    receiver = receiver,
    cutoff_visualization = cutoff_visualization,
    top_n_ligands = top_n_ligands,
    top_n_targets = top_n_targets,
    mode = mode,
    group.by = group.by,
    condition.by = condition.by,
    condition_oi = condition_oi,
    condition_reference = condition_reference,
    assay = assay,
    species = species
  )

  srt@tools[["Nichenetr"]] <- bundle
  log_message(
    "{.pkg NicheNet} analysis completed",
    message_type = "success",
    verbose = verbose
  )
  srt
}

#' @title Run MultiNicheNet analysis
#'
#' @md
#' @inheritParams RunNichenetr
#' @param sample.by Metadata column defining biological samples.
#' @param sample_agnostic Whether to use the sample-agnostic MultiNicheNet wrapper.
#' @param contrast_tbl Optional contrast table passed to MultiNicheNet. If `NULL`,
#' a simple contrast table will be created from `condition_oi` and
#' `condition_reference`.
#' @param receiver_celltypes Receiver cell types of interest.
#' @param sender_celltypes Sender cell types of interest. Default is all available
#' cell types.
#' @param batches Optional metadata column(s) used as batches.
#' @param covariates Optional metadata column(s) used as covariates.
#' @param min_cells Minimum number of cells per sample-celltype combination.
#' @param fraction_cutoff Minimum expression fraction cutoff used by MultiNicheNet.
#' @param empirical_pval Whether to use empirical p-values.
#' @param top_n_interactions Number of top prioritized interactions kept in the
#' standardized long table.
#' @return A Seurat object with standardized MultiNicheNet results stored in
#' `srt@tools[["MultiNichenetr"]]`.
#' @export
RunMultiNichenetr <- function(
  srt,
  group.by,
  sample.by,
  condition.by,
  condition_oi,
  condition_reference,
  receiver_celltypes,
  sender_celltypes = NULL,
  assay = NULL,
  sample_agnostic = FALSE,
  contrast_tbl = NULL,
  batches = NULL,
  covariates = NULL,
  species = c("Homo_sapiens", "Mus_musculus"),
  lr_network = NULL,
  ligand_target_matrix = NULL,
  fraction_cutoff = 0.05,
  min_cells = 10,
  empirical_pval = TRUE,
  top_n_interactions = 250,
  verbose = TRUE
) {
  check_r(
    c("saeyslab/multinichenetr", "SingleCellExperiment", "muscat"),
    verbose = FALSE
  )

  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  needed_cols <- c(group.by, sample.by, condition.by)
  missing_cols <- setdiff(needed_cols, colnames(srt[[]]))
  if (length(missing_cols) > 0L) {
    log_message(
      "Missing required metadata columns: {.val {missing_cols}}",
      message_type = "error"
    )
  }

  assay <- assay %||% DefaultAssay(srt)
  species <- match.arg(species)
  if (is.null(sender_celltypes)) {
    sender_celltypes <- setdiff(
      unique(as.character(srt[[group.by]][, 1])),
      receiver_celltypes
    )
  }
  batches <- batches %||% NA
  covariates <- covariates %||% NA

  model <- load_nichenetr_models(
    species = species,
    lr_network = lr_network,
    ligand_target_matrix = ligand_target_matrix,
    weighted_networks = NULL,
    verbose = verbose
  )

  sce <- Seurat::as.SingleCellExperiment(srt, assay = assay)
  celltype_raw <- as.character(srt[[group.by]][, 1])
  celltype_safe <- make.names(celltype_raw)
  celltype_map <- stats::setNames(celltype_raw, celltype_safe)
  SummarizedExperiment::colData(sce)[[group.by]] <- celltype_safe
  SummarizedExperiment::colData(sce)[[
    sample.by
  ]] <- make.names(as.character(srt[[sample.by]][, 1]))
  SummarizedExperiment::colData(sce)[[
    condition.by
  ]] <- make.names(as.character(srt[[condition.by]][, 1]))
  sender_celltypes <- make.names(sender_celltypes)
  receiver_celltypes <- make.names(receiver_celltypes)
  condition_oi_safe <- make.names(condition_oi)
  condition_ref_safe <- make.names(condition_reference)
  contrast_tbl <- data.frame(
    contrast = c(
      paste0(condition_oi_safe, "-", condition_ref_safe),
      paste0(condition_ref_safe, "-", condition_oi_safe)
    ),
    group = c(condition_oi_safe, condition_ref_safe),
    stringsAsFactors = FALSE
  )
  contrasts_oi <- paste0("'", contrast_tbl$contrast, "'", collapse = ",")

  ns_fun <- function(name) getExportedValue("multinichenetr", name)
  raw_result <- tryCatch(
    {
      if (isTRUE(sample_agnostic)) {
        ns_fun("multi_nichenet_analysis_sampleAgnostic")(
          sce = sce,
          celltype_id = group.by,
          sample_id = sample.by,
          group_id = condition.by,
          batches = batches,
          covariates = covariates,
          lr_network = normalize_lr_network_for_multinichenetr(
            model$lr_network
          ),
          ligand_target_matrix = model$ligand_target_matrix,
          contrasts_oi = contrasts_oi,
          contrast_tbl = contrast_tbl,
          senders_oi = sender_celltypes,
          receivers_oi = receiver_celltypes,
          fraction_cutoff = fraction_cutoff,
          min_cells = min_cells
        )
      } else {
        ns_fun("multi_nichenet_analysis")(
          sce = sce,
          celltype_id = group.by,
          sample_id = sample.by,
          group_id = condition.by,
          batches = batches,
          covariates = covariates,
          lr_network = normalize_lr_network_for_multinichenetr(
            model$lr_network
          ),
          ligand_target_matrix = model$ligand_target_matrix,
          contrasts_oi = contrasts_oi,
          contrast_tbl = contrast_tbl,
          senders_oi = sender_celltypes,
          receivers_oi = receiver_celltypes,
          fraction_cutoff = fraction_cutoff,
          min_cells = min_cells,
          empirical_pval = empirical_pval
        )
      }
    },
    error = function(e) {
      msg <- cli::ansi_strip(conditionMessage(e))
      log_message(
        paste0(
          "MultiNicheNet analysis failed. ",
          "Common causes: (1) not enough samples per condition (>=2 required per group in {.arg sample.by}), ",
          "(2) some cell types only present in one condition, ",
          "(3) too few cells per sample-celltype combination (adjust {.arg min_cells}). ",
          "Original error: ",
          msg
        ),
        message_type = "error"
      )
    }
  )

  bundle <- standardize_multinichenetr_result(
    raw_result = raw_result,
    top_n_interactions = top_n_interactions,
    group.by = group.by,
    sample.by = sample.by,
    condition.by = condition.by,
    condition_oi = condition_oi,
    condition_reference = condition_reference,
    receiver_celltypes = receiver_celltypes,
    sender_celltypes = sender_celltypes,
    assay = assay,
    species = species,
    sample_agnostic = sample_agnostic
  )

  srt@tools[["MultiNichenetr"]] <- bundle
  log_message(
    "{.pkg MultiNicheNet} analysis completed",
    message_type = "success",
    verbose = verbose
  )
  srt
}

resolve_sender_groups <- function(
  srt,
  group.by,
  sender,
  receiver,
  use_sender_agnostic_background = TRUE
) {
  all_groups <- unique(as.character(srt[[group.by]][, 1]))
  receiver <- unique(as.character(receiver))
  if (length(sender) == 1L && identical(sender, "all")) {
    sender_use <- if (isTRUE(use_sender_agnostic_background)) {
      setdiff(all_groups, receiver)
    } else {
      all_groups
    }
  } else {
    sender_use <- unique(as.character(sender))
  }
  sender_use
}

normalize_lr_network_for_multinichenetr <- function(lr_network) {
  lr_network <- as.data.frame(lr_network)
  nm <- colnames(lr_network)
  if (all(c("from", "to") %in% nm)) {
    colnames(lr_network)[match(c("from", "to"), nm)] <- c("ligand", "receptor")
  }
  lr_network
}

load_nichenetr_models <- function(
  species = c("Homo_sapiens", "Mus_musculus"),
  lr_network = NULL,
  ligand_target_matrix = NULL,
  weighted_networks = NULL,
  verbose = TRUE
) {
  species <- match.arg(species)
  check_r(c("saeyslab/nichenetr", "R.cache"), verbose = FALSE)

  url_map <- list(
    Homo_sapiens = list(
      lr_network = "https://zenodo.org/records/3260758/files/lr_network.rds?download=1",
      ligand_target_matrix = "https://zenodo.org/records/3260758/files/ligand_target_matrix.rds?download=1",
      weighted_networks = "https://zenodo.org/records/3260758/files/weighted_networks.rds?download=1"
    ),
    Mus_musculus = list(
      lr_network = "https://zenodo.org/records/7074291/files/lr_network_mouse_21122021.rds?download=1",
      ligand_target_matrix = "https://zenodo.org/records/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds?download=1",
      weighted_networks = "https://zenodo.org/records/7074291/files/weighted_networks_nsga2r_final_mouse.rds?download=1"
    )
  )

  lr_network <- resolve_nichenetr_object(
    x = lr_network,
    package = "nichenetr",
    object_candidates = if (species == "Mus_musculus") {
      c("lr_network_mouse_21122021", "lr_network_mouse", "lr_network")
    } else {
      c("lr_network", "lr_network_human")
    },
    fallback_url = url_map[[species]][["lr_network"]],
    verbose = verbose
  )

  ligand_target_matrix <- resolve_nichenetr_object(
    x = ligand_target_matrix,
    package = "nichenetr",
    object_candidates = if (species == "Mus_musculus") {
      c(
        "ligand_target_matrix_nsga2r_final_mouse",
        "ligand_target_matrix_mouse",
        "ligand_target_matrix"
      )
    } else {
      c("ligand_target_matrix", "ligand_target_matrix_human")
    },
    fallback_url = url_map[[species]][["ligand_target_matrix"]],
    verbose = verbose
  )

  weighted_networks <- resolve_nichenetr_object(
    x = weighted_networks,
    package = "nichenetr",
    object_candidates = if (species == "Mus_musculus") {
      c(
        "weighted_networks_nsga2r_final_mouse",
        "weighted_networks_mouse",
        "weighted_networks"
      )
    } else {
      c("weighted_networks", "weighted_networks_human")
    },
    fallback_url = url_map[[species]][["weighted_networks"]],
    allow_null = TRUE,
    verbose = verbose
  )

  if (is.null(weighted_networks)) {
    sig_network <- resolve_nichenetr_object(
      x = NULL,
      package = "nichenetr",
      object_candidates = c("sig_network"),
      fallback_url = NULL,
      allow_null = TRUE,
      verbose = verbose
    )
    gr_network <- resolve_nichenetr_object(
      x = NULL,
      package = "nichenetr",
      object_candidates = c("gr_network"),
      fallback_url = NULL,
      allow_null = TRUE,
      verbose = verbose
    )
    source_weights_df <- resolve_nichenetr_object(
      x = NULL,
      package = "nichenetr",
      object_candidates = c("source_weights_df"),
      fallback_url = NULL,
      allow_null = TRUE,
      verbose = verbose
    )

    if (
      !is.null(sig_network) &&
        !is.null(gr_network) &&
        !is.null(source_weights_df)
    ) {
      construct_weighted_networks <- get0(
        "construct_weighted_networks",
        envir = asNamespace("nichenetr"),
        inherits = FALSE
      )
      if (!is.null(construct_weighted_networks)) {
        weighted_networks <- construct_weighted_networks(
          lr_network = lr_network,
          sig_network = sig_network,
          gr_network = gr_network,
          source_weights_df = source_weights_df
        )
      }
    }
  }

  list(
    lr_network = as.data.frame(lr_network),
    ligand_target_matrix = ligand_target_matrix,
    weighted_networks = weighted_networks
  )
}

resolve_nichenetr_object <- function(
  x,
  package,
  object_candidates,
  fallback_url = NULL,
  allow_null = FALSE,
  verbose = TRUE
) {
  if (!is.null(x)) {
    if (is.character(x) && length(x) == 1L && file.exists(x)) {
      return(readRDS(x))
    }
    return(x)
  }

  ns <- asNamespace(package)
  for (nm in object_candidates) {
    obj <- tryCatch(
      get0(nm, envir = ns, inherits = FALSE),
      error = function(e) NULL
    )
    if (!is.null(obj)) {
      return(obj)
    }
  }

  if (!is.null(fallback_url)) {
    cache_key <- list("nichenetr", object_candidates[1])
    cached <- tryCatch(
      R.cache::loadCache(key = cache_key),
      error = function(e) NULL
    )
    if (!is.null(cached)) {
      log_message(
        "Loading cached NicheNet prior model: {.val {object_candidates[1]}}",
        verbose = verbose
      )
      return(cached)
    }
    log_message(
      "Downloading NicheNet prior model from {.url {fallback_url}}",
      verbose = verbose
    )
    out <- tryCatch(readRDS(url(fallback_url)), error = function(e) NULL)
    if (!is.null(out)) {
      tryCatch(
        R.cache::saveCache(
          out,
          key = cache_key,
          comment = object_candidates[1]
        ),
        error = function(e) NULL
      )
      return(out)
    }
  }

  if (isTRUE(allow_null)) {
    return(NULL)
  }

  log_message(
    "Failed to resolve required NicheNet object. Checked candidates: {.val {object_candidates}}",
    message_type = "error"
  )
}

standardize_nichenetr_result <- function(
  raw_result,
  ligand_activities,
  ligand_target_df,
  ligand_receptor_df,
  sender_use,
  receiver,
  cutoff_visualization = 0.33,
  top_n_ligands = 30,
  top_n_targets = 200,
  mode,
  group.by,
  condition.by,
  condition_oi,
  condition_reference,
  assay,
  species
) {
  ligand_activities <- standardize_df(ligand_activities)
  ligand_target_df <- standardize_df(ligand_target_df)
  ligand_receptor_df <- standardize_df(ligand_receptor_df)

  ligand_col <- ccc_pick_col(
    ligand_receptor_df,
    c("from", "ligand", "test_ligand")
  )
  receptor_col <- ccc_pick_col(ligand_receptor_df, c("to", "receptor"))
  if (
    !is.null(ligand_receptor_df) &&
      !is.null(ligand_col) &&
      !is.null(receptor_col)
  ) {
    colnames(ligand_receptor_df)[match(
      c(ligand_col, receptor_col),
      colnames(ligand_receptor_df)
    )] <- c("ligand", "receptor")
  }

  act_ligand_col <- ccc_pick_col(ligand_activities, c("test_ligand", "ligand"))
  act_score_col <- ccc_pick_col(
    ligand_activities,
    c("pearson", "aupr_corrected", "aupr", "activity", "score")
  )
  lr_table <- data.frame()
  if (!is.null(ligand_receptor_df) && nrow(ligand_receptor_df) > 0L) {
    lr_table <- ligand_receptor_df
    if (!is.null(act_ligand_col) && !is.null(act_score_col)) {
      match_idx <- match(lr_table$ligand, ligand_activities[[act_ligand_col]])
      lr_table$score <- ligand_activities[[act_score_col]][match_idx]
    } else {
      lr_table$score <- NA_real_
    }
    lr_table$sender <- I(lapply(lr_table$ligand, function(x) sender_use))
    lr_table <- expand_sender_lr_table(lr_table)
    lr_table$receiver <- paste(receiver, collapse = ",")
    lr_table$interaction_name <- paste(
      lr_table$ligand,
      lr_table$receptor,
      sep = " - "
    )
    lr_table$method <- "Nichenetr"
  }

  ligand_target_table <- data.frame()
  if (!is.null(ligand_target_df) && nrow(ligand_target_df) > 0L) {
    ligand_target_table <- ligand_target_df
    lig_col <- ccc_pick_col(
      ligand_target_table,
      c("ligand", "from", "test_ligand")
    )
    tar_col <- ccc_pick_col(ligand_target_table, c("target", "to", "gene"))
    wt_col <- ccc_pick_col(
      ligand_target_table,
      c("weight", "regulatory_potential", "score")
    )
    if (!is.null(lig_col)) {
      colnames(ligand_target_table)[match(
        lig_col,
        colnames(ligand_target_table)
      )] <- "ligand"
    }
    if (!is.null(tar_col)) {
      colnames(ligand_target_table)[match(
        tar_col,
        colnames(ligand_target_table)
      )] <- "target"
    }
    if (!is.null(wt_col)) {
      colnames(ligand_target_table)[match(
        wt_col,
        colnames(ligand_target_table)
      )] <- "weight"
    }
    ligand_target_table <- ligand_target_table[
      order(ligand_target_table$weight, decreasing = TRUE), ,
      drop = FALSE
    ]
    ligand_target_table <- utils::head(
      ligand_target_table,
      top_n_targets *
        max(1L, min(top_n_ligands, length(unique(ligand_target_table$ligand))))
    )
  }

  pair_table <- aggregate_ccc_long(lr_table)

  list(
    method = "Nichenetr",
    results = raw_result,
    ligand_activities = ligand_activities,
    ligand_target_df = ligand_target_table,
    ligand_receptor_df = ligand_receptor_df,
    long_table = lr_table,
    pair_table = pair_table,
    parameters = list(
      sender_use = sender_use,
      receiver = receiver,
      mode = mode,
      group.by = group.by,
      condition.by = condition.by,
      condition_oi = condition_oi,
      condition_reference = condition_reference,
      assay = assay,
      species = species,
      cutoff_visualization = cutoff_visualization,
      top_n_ligands = top_n_ligands,
      top_n_targets = top_n_targets
    )
  )
}

expand_sender_lr_table <- function(df) {
  if (!"sender" %in% colnames(df)) {
    return(df)
  }
  out <- lapply(seq_len(nrow(df)), function(i) {
    senders_i <- df$sender[[i]] %||% NA_character_
    senders_i <- as.character(senders_i)
    if (length(senders_i) == 0L) {
      senders_i <- NA_character_
    }
    tmp <- df[rep(i, length(senders_i)), , drop = FALSE]
    tmp$sender <- senders_i
    tmp
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

standardize_multinichenetr_result <- function(
  raw_result,
  top_n_interactions = 250,
  group.by,
  sample.by,
  condition.by,
  condition_oi,
  condition_reference,
  receiver_celltypes,
  sender_celltypes,
  assay,
  species,
  sample_agnostic = FALSE
) {
  prior_tbl <- extract_multinichenetr_prior_table(raw_result)
  prior_tbl <- standardize_df(prior_tbl)
  lr_table <- standardize_multinichenetr_lr_table(prior_tbl)
  if (
    nrow(lr_table) > 0L &&
      top_n_interactions > 0L &&
      "score" %in% colnames(lr_table)
  ) {
    lr_table <- lr_table[
      order(lr_table$score, decreasing = TRUE), ,
      drop = FALSE
    ]
    lr_table <- utils::head(lr_table, top_n_interactions)
  }

  ligand_target_df <- extract_ligand_target_from_multinichenetr(raw_result)
  pair_table <- aggregate_ccc_long(lr_table)

  list(
    method = "MultiNichenetr",
    results = raw_result,
    prioritization_table = prior_tbl,
    ligand_target_df = ligand_target_df,
    long_table = lr_table,
    pair_table = pair_table,
    parameters = list(
      group.by = group.by,
      sample.by = sample.by,
      condition.by = condition.by,
      condition_oi = condition_oi,
      condition_reference = condition_reference,
      receiver_celltypes = receiver_celltypes,
      sender_celltypes = sender_celltypes,
      assay = assay,
      species = species,
      sample_agnostic = sample_agnostic
    )
  )
}

extract_multinichenetr_prior_table <- function(x) {
  if (is.data.frame(x)) {
    return(x)
  }
  if (is.list(x)) {
    if (!is.null(x$prioritization_tables)) {
      return(extract_multinichenetr_prior_table(x$prioritization_tables))
    }
    candidates <- c(
      "group_prioritization_tbl",
      "prioritization_tbl",
      "lr_network_prioritization_tbl",
      "prioritization_table"
    )
    for (nm in candidates) {
      if (!is.null(x[[nm]])) {
        return(extract_multinichenetr_prior_table(x[[nm]]))
      }
    }
    df_list <- lapply(x, extract_multinichenetr_prior_table)
    df_list <- Filter(function(el) is.data.frame(el) && nrow(el) > 0L, df_list)
    if (length(df_list) > 0L) {
      common <- Reduce(intersect, lapply(df_list, colnames))
      if (length(common) > 0L) {
        df_list <- lapply(df_list, function(el) el[, common, drop = FALSE])
      }
      return(do.call(rbind, df_list))
    }
  }
  data.frame()
}

extract_ligand_target_from_multinichenetr <- function(x) {
  if (is.data.frame(x)) {
    nm <- colnames(x)
    if (
      any(ccc_pick_col(x, c("ligand", "from", "test_ligand", "ligand_oi")) %in% nm) &&
        any(ccc_pick_col(x, c("target", "to", "gene", "target_gene", "target_genes", "gene_oi")) %in% nm)
    ) {
      out <- x
      lig_col <- ccc_pick_col(out, c("ligand", "from", "test_ligand", "ligand_oi"))
      tar_col <- ccc_pick_col(out, c("target", "to", "gene", "target_gene", "target_genes", "gene_oi"))
      snd_col <- ccc_pick_col(out, c("sender", "sender_celltype", "celltype_sender", "source_celltype"))
      rcv_col <- ccc_pick_col(out, c("receiver", "receiver_celltype", "celltype_receiver", "target_celltype"))
      wt_col <- ccc_pick_col(
        out,
        c(
          "weight",
          "score",
          "regulatory_potential",
          "pearson",
          "aupr_corrected",
          "aupr",
          "activity",
          "prioritization_score"
        )
      )
      if (!is.null(lig_col)) {
        colnames(out)[match(lig_col, colnames(out))] <- "ligand"
      }
      if (!is.null(tar_col)) {
        colnames(out)[match(tar_col, colnames(out))] <- "target"
      }
      if (!is.null(snd_col)) {
        colnames(out)[match(snd_col, colnames(out))] <- "sender"
      }
      if (!is.null(rcv_col)) {
        colnames(out)[match(rcv_col, colnames(out))] <- "receiver"
      }
      if (!is.null(wt_col)) {
        colnames(out)[match(wt_col, colnames(out))] <- "weight"
      } else {
        numeric_cols <- setdiff(
          colnames(out)[vapply(out, is.numeric, logical(1))],
          c("ligand", "target")
        )
        if (length(numeric_cols) > 0L) {
          colnames(out)[match(numeric_cols[1], colnames(out))] <- "weight"
        }
      }
      return(out)
    }
  }
  if (is.list(x)) {
    candidates <- c(
      "ligand_activities_targets_DEgenes",
      "ligand_target_df",
      "ligand_target_table"
    )
    for (nm in candidates) {
      if (!is.null(x[[nm]])) {
        out <- extract_ligand_target_from_multinichenetr(x[[nm]])
        if (is.data.frame(out) && nrow(out) > 0L) {
          return(out)
        }
      }
    }
    for (el in x) {
      out <- extract_ligand_target_from_multinichenetr(el)
      if (is.data.frame(out) && nrow(out) > 0L) {
        return(out)
      }
    }
  }
  data.frame()
}

standardize_multinichenetr_lr_table <- function(df) {
  df <- standardize_df(df)
  if (is.null(df) || nrow(df) == 0L) {
    return(data.frame())
  }
  out <- df
  out <- rename_by_candidates(
    out,
    "sender",
    c("sender", "from", "celltype_sender", "sender_celltype", "source")
  )
  out <- rename_by_candidates(
    out,
    "receiver",
    c("receiver", "to", "celltype_receiver", "receiver_celltype", "target")
  )
  out <- rename_by_candidates(
    out,
    "ligand",
    c("ligand", "from_ligand", "ligand_oi", "ligand_source")
  )
  out <- rename_by_candidates(
    out,
    "receptor",
    c("receptor", "to_receptor", "receptor_oi", "receptor_target")
  )
  out <- rename_by_candidates(
    out,
    "score",
    c(
      "prioritization_score",
      "scaled_ligand_receptor_score",
      "ligand_receptor_score",
      "score",
      "activity"
    )
  )
  out <- rename_by_candidates(
    out,
    "pvalue",
    c("pvalue", "p_val", "p_val_ligand_activity", "p_val_adj")
  )

  for (nm in c("sender", "receiver", "ligand", "receptor")) {
    if (!nm %in% colnames(out)) {
      out[[nm]] <- NA_character_
    }
  }
  if (!"score" %in% colnames(out)) {
    out$score <- NA_real_
  }
  if (!"pvalue" %in% colnames(out)) {
    out$pvalue <- NA_real_
  }
  out$interaction_name <- paste(out$ligand, out$receptor, sep = " - ")
  out$method <- "MultiNichenetr"
  out
}

rename_by_candidates <- function(df, new_name, candidates) {
  col <- ccc_pick_col(df, candidates)
  if (!is.null(col) && !new_name %in% colnames(df)) {
    colnames(df)[match(col, colnames(df))] <- new_name
  }
  df
}

standardize_df <- function(x) {
  if (is.null(x)) {
    return(data.frame())
  }
  if (inherits(x, "python.builtin.object")) {
    x <- py_to_r2(x)
  }
  if (is.data.frame(x)) {
    return(x)
  }
  if (is.matrix(x)) {
    return(as.data.frame(x))
  }
  data.frame()
}

ccc_pick_col <- function(df, candidates) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L && ncol(df) == 0L) {
    return(NULL)
  }
  normalize_name <- function(x) {
    tolower(gsub("[^[:alnum:]]+", "", x))
  }
  cn <- colnames(df)
  cn_norm <- normalize_name(cn)
  cand_norm <- normalize_name(candidates)
  hit <- match(cand_norm, cn_norm)
  hit <- hit[!is.na(hit)]
  if (length(hit) == 0L) {
    return(NULL)
  }
  cn[hit[1]]
}
