#' @title Run SCENICPlus-style eGRN analysis
#'
#' @description
#' Build a SCENICPlus-style eRegulon result from paired RNA and chromatin assays.
#' The native backend is a `scop` implementation and is not an exact re-run of
#' the official SCENIC+ Python workflow.
#'
#' @md
#' @inheritParams standard_scop
#' @param srt A Seurat object containing RNA and chromatin assays.
#' @param rna_assay RNA assay name.
#' @param atac_assay Chromatin assay name.
#' @param rna_layer RNA count layer used for TF-gene inference and AUC scoring.
#' @param atac_layer ATAC count layer used for peak-gene correlations.
#' @param backend Runtime backend. `"cpp"` uses the native package workflow.
#' `"python"` checks the official Python environment but does not silently
#' substitute the native result.
#' @param grn_method GRN method used for TF-gene inference. The native backend
#' can call `RunGRNBoost2()`, `RunGENIE3()`, or `RunGNIPLR()`;
#' RegDiffusion and GNIPLR are supported through the Python backend.
#' @param regulators Candidate transcription factors. If `NULL`, motif TF names
#' are used when available; otherwise all RNA genes are considered.
#' @param rna_expr Optional expression matrix used for native TF-gene and
#' region-gene inference. Rows should be genes and columns should be cells; a
#' cell-by-gene matrix is accepted and transposed when cell names make the
#' orientation unambiguous.
#' @param atac_expr Optional accessibility matrix used for native region-gene
#' inference. Rows should be regions and columns should be cells; a cell-by-region
#' matrix is accepted and transposed when cell names make the orientation
#' unambiguous.
#' @param tf_gene_prior Optional precomputed TF-gene table with columns `TF`,
#' `target`, and `importance`. When supplied with `backend = "cpp"`, native
#' GRN inference is skipped and the table is used as the TF-gene evidence layer.
#' @param region_gene_prior Optional precomputed region-gene table with columns
#' `region`, `gene`, and `score`. When supplied with `backend = "cpp"`, native
#' peak-to-gene correlation is skipped.
#' @param tf_region_prior Optional precomputed TF-region table with columns
#' `TF`, `region`, and `score`. When supplied with `backend = "cpp"`, existing
#' motif incidence is not required.
#' @param triplets_prior Optional precomputed TF-region-gene table with columns
#' `TF`, `region`, `gene`, and `score`. When supplied with `backend = "cpp"`,
#' triplets are read directly instead of assembled from the three edge layers.
#' @param eregulons_prior Optional precomputed eRegulon table with columns
#' `regulon` and `target`. When supplied with `backend = "cpp"`, eRegulon names
#' and gene sets are taken from this table.
#' @param auc_expr Optional expression matrix used only for eRegulon AUCell
#' scoring in the native backend. Rows should be genes and columns should be
#' cells; a cell-by-gene matrix is accepted and transposed when cell names make
#' the orientation unambiguous.
#' @param auc_rankings Optional precomputed 0-based ranking matrix used only for
#' eRegulon AUCell scoring in the native backend. Rows should be cells and
#' columns should be genes; a gene-by-cell matrix is accepted and transposed when
#' cell names make the orientation unambiguous. This is mainly useful for
#' official SCENIC+ parity tests because SCENIC+ uses seeded random tie-breaking.
#' @param region_gene_window Maximum peak-to-gene distance in base pairs.
#' @param region_gene_search_space Optional search-space table with columns
#' `region`/`Name` and `gene`/`Gene`. When supplied, native region-gene scoring
#' is limited to these candidate pairs.
#' @param region_gene_method Native region-gene scoring method. `"gbm"` follows
#' the official SCENIC+ region-to-gene strategy more closely by combining
#' boosting feature importance with Spearman correlation; `"correlation"` keeps
#' the older correlation-only approximation.
#' @param region_gene_n_rounds Number of boosting rounds for native
#' region-gene importance scoring.
#' @param max_region_gene Maximum region-gene links retained per gene.
#' @param max_tf_region Maximum motif-supported regions retained per TF.
#' @param min_eregulon_size Minimum genes per eRegulon retained for AUC scoring.
#' @param egrn_rho_threshold Absolute Spearman correlation threshold used to
#' split TF-gene and region-gene links into activating/repressing eGRN modules.
#' @param egrn_quantiles Region-gene importance quantiles used to build
#' candidate eModules before TF-gene leading-edge filtering.
#' @param egrn_top_n_region_gene Top-N region-gene links per gene used as
#' additional candidate eModules before TF-gene leading-edge filtering.
#' @param egrn_min_target_genes Minimum leading-edge target genes required for
#' an eRegulon module. The default follows the official SCENIC+ CLI workflow.
#' @param group.by Optional metadata column used to calculate eRegulon
#' specificity scores.
#' @param grn_top_targets Maximum TF-gene edges retained per target. The
#' default `Inf` keeps all positive links to match arboreto GRNBoost2 output.
#' @param max_grn_targets Maximum genes used as TF-gene GRN targets in the
#' native backend. With `grn_target_scope = "region_gene"`, genes are ranked by
#' strongest absolute peak-to-gene correlation. Set to `Inf` to keep the full
#' selected target universe.
#' @param grn_target_scope TF-gene target universe for the native backend.
#' `"all"` follows official SCENIC+ TF-to-gene inference by scoring every RNA
#' gene; `"region_gene"` restricts TF-gene inference to genes with peak links.
#' @param grn_n_rounds Number of native tree boosting rounds for TF-gene
#' inference. The default follows arboreto `SGBM_KWARGS`.
#' @param grn_learning_rate Native tree boosting learning rate.
#' @param grn_max_depth Maximum depth of each native TF-gene regression tree.
#' @param grn_max_features Fraction of candidate TFs sampled at each native
#' tree split.
#' @param grn_subsample Fraction of cells sampled for each native boosting
#' round.
#' @param grn_early_stop_window_length Out-of-bag improvement window used for
#' native GRNBoost2 early stopping. The RunSCENICPlus default is `0`, which
#' disables early stopping to match official SCENIC+ `SGBM_KWARGS`.
#' @param seed Random seed used by native stochastic boosting backends. The
#' default matches official SCENIC+ command-line wrappers.
#' @param cores Number of workers used by native TF-gene GRN inference.
#' @param assay_name Assay used to store eRegulon AUC scores.
#' @param tool_name Name of the `srt@tools` result entry.
#' @param python_result_dir Directory containing official SCENIC+ outputs already
#' exported as `tf_gene`, `region_gene`, `tf_region`, `triplets`, `eregulons`,
#' and `auc` tables. Used only for `backend = "python"`.
#' @param scplus_object Optional official SCENIC+ Python object path
#' (`.pkl`/`.pickle`) to extract into `python_result_dir` before
#' standardization. Used only for `backend = "python"`.
#' @param envname Python environment used when `backend = "python"`.
#' @param conda Conda-compatible executable used when `backend = "python"`.
#'
#' @return A Seurat object with SCENICPlus-style results.
#' @export
RunSCENICPlus <- function(
  srt,
  rna_assay = "RNA",
  atac_assay = "peaks",
  rna_layer = "counts",
  atac_layer = "counts",
  backend = c("cpp", "python"),
  grn_method = c("grnboost2", "regdiffusion", "genie3", "gniplr"),
  regulators = NULL,
  rna_expr = NULL,
  atac_expr = NULL,
  tf_gene_prior = NULL,
  region_gene_prior = NULL,
  tf_region_prior = NULL,
  triplets_prior = NULL,
  eregulons_prior = NULL,
  auc_expr = NULL,
  auc_rankings = NULL,
  region_gene_window = 250000,
  region_gene_search_space = NULL,
  region_gene_method = c("gbm", "correlation"),
  region_gene_n_rounds = 500,
  max_region_gene = 5,
  max_tf_region = 500,
  min_eregulon_size = 5,
  egrn_rho_threshold = 0.05,
  egrn_quantiles = c(0.85, 0.90, 0.95),
  egrn_top_n_region_gene = c(5L, 10L, 15L),
  egrn_min_target_genes = 10L,
  group.by = NULL,
  grn_top_targets = Inf,
  max_grn_targets = Inf,
  grn_target_scope = c("all", "region_gene"),
  grn_n_rounds = 5000,
  grn_learning_rate = 0.01,
  grn_max_depth = 3,
  grn_max_features = 0.1,
  grn_subsample = 0.9,
  grn_early_stop_window_length = 0,
  seed = 666,
  cores = 1,
  assay_name = "scenicplus",
  tool_name = "SCENICPlus",
  python_result_dir = NULL,
  scplus_object = NULL,
  envname = "scenicplus_env",
  conda = "auto",
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  backend <- match.arg(backend)
  grn_method <- match.arg(grn_method)
  region_gene_method <- match.arg(region_gene_method)
  grn_target_scope <- match.arg(grn_target_scope)
  if (identical(backend, "python")) {
    return(run_scenicplus_python(
      srt = srt,
      envname = envname,
      conda = conda,
      python_result_dir = python_result_dir,
      scplus_object = scplus_object,
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      assay_name = assay_name,
      tool_name = tool_name,
      group.by = group.by,
      verbose = verbose
    ))
  }
  if (!rna_assay %in% names(srt@assays)) {
    log_message(
      "{.arg rna_assay} is not present in {.arg srt}",
      message_type = "error"
    )
  }
  if (!atac_assay %in% names(srt@assays)) {
    log_message(
      "{.arg atac_assay} is not present in {.arg srt}",
      message_type = "error"
    )
  }
  if (!inherits(srt[[atac_assay]], "ChromatinAssay")) {
    log_message(
      "{.arg atac_assay} must refer to a {.cls ChromatinAssay}",
      message_type = "error"
    )
  }

  rna_counts <- GetAssayData5(srt, assay = rna_assay, layer = rna_layer)
  atac_counts <- GetAssayData5(srt, assay = atac_assay, layer = atac_layer)
  shared_cells <- intersect(colnames(rna_counts), colnames(atac_counts))
  if (length(shared_cells) < 3L) {
    log_message(
      "{.fn RunSCENICPlus} requires at least three shared cells between RNA and ATAC assays",
      message_type = "error"
    )
  }
  rna_counts <- rna_counts[, shared_cells, drop = FALSE]
  atac_counts <- atac_counts[, shared_cells, drop = FALSE]
  rna_counts <- scenicplus_prepare_rna_expr(
    rna_expr = rna_expr,
    fallback = rna_counts,
    cells = shared_cells
  )
  atac_counts <- scenicplus_prepare_atac_expr(
    atac_expr = atac_expr,
    fallback = atac_counts,
    cells = shared_cells
  )
  shared_cells <- intersect(colnames(rna_counts), colnames(atac_counts))
  rna_counts <- rna_counts[, shared_cells, drop = FALSE]
  atac_counts <- atac_counts[, shared_cells, drop = FALSE]
  auc_counts <- scenicplus_prepare_auc_expr(
    auc_expr = auc_expr,
    fallback = rna_counts,
    cells = shared_cells
  )
  auc_rankings <- scenicplus_prep_auc_rank(
    auc_rankings = auc_rankings,
    auc_counts = auc_counts
  )

  tf_region_prior <- scenicplus_standardize_prior(
    tf_region_prior,
    cols = c("TF", "region", "score"),
    label = "tf_region_prior"
  )
  region_gene_prior <- scenicplus_standardize_prior(
    region_gene_prior,
    cols = c("region", "gene", "score"),
    label = "region_gene_prior"
  )
  tf_gene_prior <- scenicplus_standardize_prior(
    tf_gene_prior,
    cols = c("TF", "target", "importance"),
    label = "tf_gene_prior"
  )
  triplets_prior <- scenicplus_standardize_prior(
    triplets_prior,
    cols = c("TF", "region", "gene", "score"),
    label = "triplets_prior"
  )
  eregulons_prior <- scenicplus_standardize_prior(
    eregulons_prior,
    cols = c("regulon", "target"),
    label = "eregulons_prior"
  )

  motif_links <- if (!is.null(tf_region_prior)) {
    scenicplus_filter_tfr_prior(
      tf_region_prior = tf_region_prior,
      peak_names = rownames(atac_counts),
      max_tf_region = max_tf_region
    )
  } else {
    scenicplus_tfr_motifs(
      srt = srt,
      atac_assay = atac_assay,
      peak_names = rownames(atac_counts),
      max_tf_region = max_tf_region
    )
  }
  motif_tfs <- unique(motif_links[["TF"]])
  regulators <- regulators %||% motif_tfs
  regulators <- intersect(regulators, rownames(rna_counts))
  if (length(regulators) == 0L) {
    log_message(
      "No candidate TFs are shared between RNA genes and motif annotations",
      message_type = "error"
    )
  }

  region_gene <- if (!is.null(region_gene_prior)) {
    scenicplus_filter_rg_prior(
      region_gene_prior = region_gene_prior,
      rna_counts = rna_counts,
      atac_counts = atac_counts,
      max_region_gene = max_region_gene
    )
  } else {
    scenicplus_region_gene_native(
      srt = srt,
      rna_counts = rna_counts,
      atac_counts = atac_counts,
      atac_assay = atac_assay,
      window = region_gene_window,
      max_region_gene = max_region_gene,
      search_space = region_gene_search_space,
      method = region_gene_method,
      n_rounds = region_gene_n_rounds,
      max_features = grn_max_features
    )
  }
  grn_targets <- scenicplus_select_grn_targets(
    region_gene = region_gene,
    rna_counts = rna_counts,
    max_grn_targets = max_grn_targets,
    target_scope = grn_target_scope
  )
  if (length(grn_targets) == 0L) {
    log_message(
      "No region-linked genes are available for TF-gene inference",
      message_type = "error"
    )
  }
  tf_gene <- if (!is.null(tf_gene_prior)) {
    scenicplus_filter_tfg_prior(
      tf_gene_prior = tf_gene_prior,
      regulators = regulators,
      targets = grn_targets,
      max_edges_per_target = grn_top_targets
    )
  } else {
    scenicplus_tf_gene_native(
      rna_counts = rna_counts,
      regulators = regulators,
      targets = grn_targets,
      grn_method = grn_method,
      max_edges_per_target = grn_top_targets,
      n_rounds = grn_n_rounds,
      learning_rate = grn_learning_rate,
      max_depth = grn_max_depth,
      max_features = grn_max_features,
      subsample = grn_subsample,
      early_stop_window_length = grn_early_stop_window_length,
      seed = seed,
      cores = cores,
      envname = envname,
      conda = conda,
      verbose = verbose
    )
  }
  tf_region <- motif_links[
    motif_links[["TF"]] %in% unique(tf_gene[["TF"]]), ,
    drop = FALSE
  ]
  triplets <- if (!is.null(triplets_prior)) {
    scenicplus_filter_tri_prior(
      triplets_prior = triplets_prior,
      tf_gene = tf_gene,
      region_gene = region_gene,
      tf_region = tf_region
    )
  } else {
    scenicplus_triplets(tf_gene, region_gene, tf_region)
  }
  if (!is.null(eregulons_prior)) {
    eregulon_info <- scenicplus_ereg_prior(
      eregulons_prior = eregulons_prior,
      triplets = triplets,
      min_size = min_eregulon_size
    )
    eregulons <- eregulon_info[["regulon_list"]]
    eregulon_table <- eregulon_info[["eregulons"]]
  } else {
    eregulon_info <- scenicplus_eregulon_build(
      triplets = triplets,
      tf_gene = tf_gene,
      region_gene = region_gene,
      rna_counts = rna_counts,
      min_size = min_eregulon_size,
      rho_threshold = egrn_rho_threshold,
      quantiles = egrn_quantiles,
      top_n_region_gene = egrn_top_n_region_gene,
      min_target_genes = egrn_min_target_genes
    )
    triplets <- eregulon_info[["triplets"]]
    eregulons <- eregulon_info[["regulon_list"]]
    eregulon_table <- data.frame(
      regulon = rep(names(eregulons), lengths(eregulons)),
      target = unlist(eregulons, use.names = FALSE),
      target_count = rep(as.integer(lengths(eregulons)), lengths(eregulons)),
      stringsAsFactors = FALSE
    )
  }
  if (length(eregulons) == 0L) {
    log_message(
      "No eRegulons remain after integrating TF-gene, region-gene, and TF-region links",
      message_type = "error"
    )
  }

  auc <- if (!is.null(auc_rankings)) {
    run_aucell_scores_from_rankings(
      rankings = auc_rankings,
      gene_sets = eregulons
    )
  } else {
    run_aucell_scores(
      expr_counts = auc_counts,
      gene_sets = eregulons,
      strategy = "full",
      algorithm = "ctxcore"
    )
  }
  auc <- as.data.frame(auc, check.names = FALSE)[
    colnames(auc_counts), ,
    drop = FALSE
  ]
  auc_mat <- as_matrix(auc)
  scores <- Matrix::t(methods::as(Matrix::Matrix(auc_mat, sparse = TRUE), "dgCMatrix"))
  rss <- scenicplus_calc_rss(
    auc = auc,
    srt = srt,
    cells = rownames(auc),
    group.by = group.by
  )

  result <- list(
    tf_gene = tf_gene,
    region_gene = region_gene,
    tf_region = tf_region,
    triplets = triplets,
    eregulons = eregulon_table,
    regulon_list = eregulons,
    auc = auc,
    rss = rss[["rss"]],
    rss_matrix = rss[["rss_matrix"]],
    scores = scores,
    parameters = list(
      backend = backend,
      method = "scenicplus_style_native_approx",
      grn_method = grn_method,
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      region_gene_window = region_gene_window,
      region_gene_method = region_gene_method,
      region_gene_n_rounds = region_gene_n_rounds,
      max_region_gene = max_region_gene,
      max_tf_region = max_tf_region,
      min_eregulon_size = min_eregulon_size,
      egrn_rho_threshold = egrn_rho_threshold,
      egrn_quantiles = egrn_quantiles,
      egrn_top_n_region_gene = egrn_top_n_region_gene,
      egrn_min_target_genes = egrn_min_target_genes,
      max_grn_targets = max_grn_targets,
      grn_target_scope = grn_target_scope,
      grn_target_count = length(grn_targets),
      seed = seed,
      prior_layers = c(
        if (!is.null(tf_gene_prior)) "tf_gene",
        if (!is.null(region_gene_prior)) "region_gene",
        if (!is.null(tf_region_prior)) "tf_region",
        if (!is.null(triplets_prior)) "triplets",
        if (!is.null(eregulons_prior)) "eregulons"
      ),
      auc_expr = if (is.null(auc_expr)) "rna_counts" else "custom",
      auc_rankings = !is.null(auc_rankings),
      group.by = group.by
    )
  )
  srt <- scenicplus_store_result(
    srt = srt,
    result = result,
    assay_name = assay_name,
    tool_name = tool_name
  )
  log_message(
    "{.pkg SCENICPlus}-style results stored in assay {.val {assay_name}} and tools slot {.val {tool_name}}",
    verbose = verbose
  )
  srt
}

run_scenicplus_python <- function(
  srt,
  envname,
  conda,
  python_result_dir = NULL,
  scplus_object = NULL,
  rna_assay = "RNA",
  atac_assay = "peaks",
  assay_name = "scenicplus",
  tool_name = "SCENICPlus",
  group.by = NULL,
  verbose = TRUE
) {
  needs_python_runtime <- !is.null(scplus_object)
  if (isTRUE(needs_python_runtime)) {
    PrepareEnv(
      envname = envname,
      conda = conda,
      modules = "scenicplus",
      verbose = verbose
    )
    check_python(
      packages = "scenicplus",
      envname = envname,
      conda = conda,
      verbose = verbose
    )
  }
  if (!is.null(scplus_object)) {
    if (is.null(python_result_dir)) {
      python_result_dir <- file.path(tempdir(), "scenicplus_export")
    }
    conda_resolved <- resolve_conda(conda)
    python_path <- conda_python(conda = conda_resolved, envname = envname)
    assert_python_runtime_switchable(python_path)
    configure_python_runtime(python_path)
    functions <- reticulate::import_from_path(
      "functions",
      path = system.file("python", package = "scop", mustWork = TRUE),
      convert = TRUE
    )
    functions$ExtractSCENICPlusTables(
      scplus_object = scplus_object,
      output_dir = python_result_dir,
      verbose = isTRUE(verbose)
    )
  }
  if (is.null(python_result_dir)) {
    log_message(
      "{.arg backend = 'python'} for {.fn RunSCENICPlus} standardizes official SCENIC+ outputs only. Provide {.arg python_result_dir} with exported official tables or {.arg scplus_object} to extract them; the native approximation is not used as a substitute.",
      message_type = "error"
    )
  }
  result <- scenicplus_read_python_result(
    result_dir = python_result_dir,
    srt = srt,
    cells = colnames(srt),
    group.by = group.by,
    parameters = list(
      backend = "python",
      method = "official_scenicplus_standardized",
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      python_result_dir = normalizePath(python_result_dir, mustWork = FALSE),
      scplus_object = scplus_object,
      group.by = group.by
    )
  )
  scenicplus_store_result(
    srt = srt,
    result = result,
    assay_name = assay_name,
    tool_name = tool_name
  )
}

scenicplus_store_result <- function(srt, result, assay_name, tool_name) {
  scores <- result[["scores"]]
  if (is.null(scores)) {
    auc <- result[["auc"]]
    auc_mat <- as_matrix(auc)
    scores <- Matrix::t(methods::as(Matrix::Matrix(auc_mat, sparse = TRUE), "dgCMatrix"))
    result[["scores"]] <- scores
  }
  regulon_list <- result[["regulon_list"]]
  target_count <- if (!is.null(regulon_list)) {
    lengths(regulon_list)[rownames(scores)]
  } else if (
    !is.null(result[["eregulons"]]) &&
      "target_count" %in% colnames(result[["eregulons"]])
  ) {
    stats::setNames(
      result[["eregulons"]][["target_count"]],
      result[["eregulons"]][["regulon"]]
    )[rownames(scores)]
  } else {
    rep(NA_integer_, nrow(scores))
  }
  score_termnames <- rownames(scores)
  srt[[assay_name]] <- Seurat::CreateAssayObject(data = scores)
  assay_features <- rownames(srt[[assay_name]])
  srt[[assay_name]] <- Seurat::AddMetaData(
    object = srt[[assay_name]],
    metadata = data.frame(
      termnames = score_termnames,
      target_count = as.integer(target_count),
      row.names = assay_features,
      stringsAsFactors = FALSE
    )
  )
  srt@tools[[tool_name]] <- result
  srt
}

scenicplus_read_python_result <- function(
  result_dir,
  srt,
  cells,
  group.by = NULL,
  parameters = list()
) {
  if (!dir.exists(result_dir)) {
    log_message(
      "Cannot find {.arg python_result_dir}: {.file {result_dir}}",
      message_type = "error"
    )
  }
  tf_gene <- scenicplus_read_required_table(
    result_dir,
    "tf_gene",
    c("TF", "target", "importance")
  )
  region_gene <- scenicplus_read_required_table(
    result_dir,
    "region_gene",
    c("region", "gene", "score")
  )
  tf_region <- scenicplus_read_required_table(
    result_dir,
    "tf_region",
    c("TF", "region", "score")
  )
  triplets <- scenicplus_read_required_table(
    result_dir,
    "triplets",
    c("TF", "region", "gene", "score")
  )
  eregulons_path <- scenicplus_table_path(result_dir, "eregulons")
  eregulons <- if (is.null(eregulons_path)) NULL else scenicplus_read_delim(eregulons_path)
  auc <- scenicplus_read_required_table(result_dir, "auc", character())
  auc <- scenicplus_normalize_auc_table(auc, cells = cells)
  regulon_list <- scenicplus_regulon_tbl(
    auc = auc,
    triplets = triplets,
    eregulons = eregulons
  )
  if (
    is.null(eregulons) ||
      !all(c("regulon", "target_count") %in% colnames(eregulons))
  ) {
    eregulons <- data.frame(
      regulon = names(regulon_list),
      target_count = lengths(regulon_list),
      stringsAsFactors = FALSE
    )
  }
  rss <- scenicplus_calc_rss(
    auc = auc,
    srt = srt,
    cells = rownames(auc),
    group.by = group.by
  )
  auc_mat <- as_matrix(auc)
  scores <- Matrix::t(methods::as(Matrix::Matrix(auc_mat, sparse = TRUE), "dgCMatrix"))
  list(
    tf_gene = tf_gene,
    region_gene = region_gene,
    tf_region = tf_region,
    triplets = triplets,
    eregulons = eregulons,
    regulon_list = regulon_list,
    auc = auc,
    rss = rss[["rss"]],
    rss_matrix = rss[["rss_matrix"]],
    scores = scores,
    parameters = parameters
  )
}

scenicplus_table_path <- function(result_dir, name) {
  candidates <- file.path(
    result_dir,
    paste0(name, c(".tsv", ".csv", ".txt"))
  )
  hits <- candidates[file.exists(candidates)]
  if (length(hits) > 0L) {
    return(hits[[1L]])
  }
  NULL
}

scenicplus_read_required_table <- function(result_dir, name, required_cols) {
  path <- scenicplus_table_path(result_dir, name)
  if (is.null(path)) {
    log_message(
      "Missing official SCENIC+ standardized table {.file {name}.tsv} in {.arg python_result_dir}",
      message_type = "error"
    )
  }
  out <- scenicplus_read_delim(path)
  missing_cols <- setdiff(required_cols, colnames(out))
  if (length(missing_cols) > 0L) {
    log_message(
      "SCENIC+ table {.file {basename(path)}} is missing column{?s} {.val {missing_cols}}",
      message_type = "error"
    )
  }
  out
}

scenicplus_read_delim <- function(path) {
  sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  utils::read.delim(
    path,
    sep = sep,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

scenicplus_prepare_rna_expr <- function(rna_expr, fallback, cells) {
  if (is.null(rna_expr)) {
    return(fallback)
  }
  if (!is.matrix(rna_expr) && !inherits(rna_expr, "Matrix")) {
    rna_expr <- as_matrix(rna_expr)
  }
  rna_expr <- methods::as(Matrix::Matrix(rna_expr, sparse = TRUE), "dgCMatrix")
  if (is.null(rownames(rna_expr)) || is.null(colnames(rna_expr))) {
    log_message(
      "{.arg rna_expr} must have gene and cell dimnames",
      message_type = "error"
    )
  }
  col_hits <- intersect(colnames(rna_expr), cells)
  row_hits <- intersect(rownames(rna_expr), cells)
  if (length(col_hits) == 0L && length(row_hits) > 0L) {
    rna_expr <- Matrix::t(rna_expr)
    col_hits <- intersect(colnames(rna_expr), cells)
  }
  if (length(col_hits) == 0L) {
    log_message(
      "{.arg rna_expr} has no cells matching the RNA/ATAC assays",
      message_type = "error"
    )
  }
  if (length(intersect(rownames(rna_expr), rownames(fallback))) == 0L) {
    log_message(
      "{.arg rna_expr} has no genes matching the RNA assay",
      message_type = "error"
    )
  }
  rna_expr[, col_hits, drop = FALSE]
}

scenicplus_prepare_atac_expr <- function(atac_expr, fallback, cells) {
  if (is.null(atac_expr)) {
    return(fallback)
  }
  if (!is.matrix(atac_expr) && !inherits(atac_expr, "Matrix")) {
    atac_expr <- as_matrix(atac_expr)
  }
  atac_expr <- methods::as(Matrix::Matrix(atac_expr, sparse = TRUE), "dgCMatrix")
  if (is.null(rownames(atac_expr)) || is.null(colnames(atac_expr))) {
    log_message(
      "{.arg atac_expr} must have region and cell dimnames",
      message_type = "error"
    )
  }
  rownames(atac_expr) <- sub(":", "-", rownames(atac_expr), fixed = TRUE)
  col_hits <- intersect(colnames(atac_expr), cells)
  row_hits <- intersect(rownames(atac_expr), cells)
  if (length(col_hits) == 0L && length(row_hits) > 0L) {
    atac_expr <- Matrix::t(atac_expr)
    rownames(atac_expr) <- sub(":", "-", rownames(atac_expr), fixed = TRUE)
    col_hits <- intersect(colnames(atac_expr), cells)
  }
  if (length(col_hits) == 0L) {
    log_message(
      "{.arg atac_expr} has no cells matching the RNA/ATAC assays",
      message_type = "error"
    )
  }
  if (length(intersect(rownames(atac_expr), rownames(fallback))) == 0L) {
    log_message(
      "{.arg atac_expr} has no regions matching the ATAC assay",
      message_type = "error"
    )
  }
  atac_expr[, col_hits, drop = FALSE]
}

scenicplus_standardize_prior <- function(x, cols, label) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) {
      log_message(
        "{.arg {label}} file does not exist: {.file {x}}",
        message_type = "error"
      )
    }
    x <- scenicplus_read_delim(x)
  }
  if (!is.data.frame(x)) {
    log_message(
      "{.arg {label}} must be a data frame or delimited file path",
      message_type = "error"
    )
  }
  missing_cols <- setdiff(cols, colnames(x))
  if (length(missing_cols) > 0L) {
    log_message(
      "{.arg {label}} is missing column{?s} {.val {missing_cols}}",
      message_type = "error"
    )
  }
  x <- x[, unique(c(cols, colnames(x))), drop = FALSE]
  for (col in intersect(cols, c("TF", "target", "region", "gene", "regulon"))) {
    x[[col]] <- as.character(x[[col]])
  }
  for (col in intersect(cols, c("importance", "score"))) {
    x[[col]] <- suppressWarnings(as.numeric(x[[col]]))
  }
  x <- x[stats::complete.cases(x[cols]), , drop = FALSE]
  rownames(x) <- NULL
  x
}

scenicplus_filter_rg_prior <- function(
  region_gene_prior,
  rna_counts,
  atac_counts,
  max_region_gene = 5
) {
  out <- region_gene_prior[
    region_gene_prior[["region"]] %in%
      rownames(atac_counts) &
      region_gene_prior[["gene"]] %in% rownames(rna_counts) &
      is.finite(region_gene_prior[["score"]]), ,
    drop = FALSE
  ]
  if (nrow(out) == 0L) {
    log_message(
      "{.arg region_gene_prior} has no links matching the RNA/ATAC assays",
      message_type = "error"
    )
  }
  max_region_gene <- suppressWarnings(as.numeric(max_region_gene))
  if (length(max_region_gene) == 1L && is.finite(max_region_gene)) {
    max_region_gene <- max(1L, as.integer(max_region_gene))
    out <- out[
      order(out[["gene"]], -abs(out[["score"]]), out[["region"]]), ,
      drop = FALSE
    ]
    out <- do.call(
      rbind,
      lapply(split(out, out[["gene"]]), utils::head, max_region_gene)
    )
  }
  first_cols <- c("region", "gene", "score")
  out <- out[,
    unique(c(first_cols, setdiff(colnames(out), first_cols))),
    drop = FALSE
  ]
  rownames(out) <- NULL
  out
}

scenicplus_filter_tfr_prior <- function(
  tf_region_prior,
  peak_names,
  max_tf_region = 500
) {
  out <- tf_region_prior[
    tf_region_prior[["region"]] %in%
      peak_names &
      is.finite(tf_region_prior[["score"]]), ,
    drop = FALSE
  ]
  if (nrow(out) == 0L) {
    log_message(
      "{.arg tf_region_prior} has no regions matching the ATAC assay",
      message_type = "error"
    )
  }
  max_tf_region <- suppressWarnings(as.numeric(max_tf_region))
  if (length(max_tf_region) == 1L && is.finite(max_tf_region)) {
    max_tf_region <- max(1L, as.integer(max_tf_region))
    out <- out[
      order(out[["TF"]], -out[["score"]], out[["region"]]), ,
      drop = FALSE
    ]
    out <- do.call(
      rbind,
      lapply(split(out, out[["TF"]]), utils::head, max_tf_region)
    )
  }
  first_cols <- c("TF", "region", "score")
  out <- out[,
    unique(c(first_cols, setdiff(colnames(out), first_cols))),
    drop = FALSE
  ]
  rownames(out) <- NULL
  out
}

scenicplus_filter_tfg_prior <- function(
  tf_gene_prior,
  regulators,
  targets,
  max_edges_per_target = Inf
) {
  out <- tf_gene_prior[
    tf_gene_prior[["TF"]] %in%
      regulators &
      tf_gene_prior[["target"]] %in% targets &
      is.finite(tf_gene_prior[["importance"]]), ,
    drop = FALSE
  ]
  if (nrow(out) == 0L) {
    log_message(
      "{.arg tf_gene_prior} has no links matching candidate TFs and targets",
      message_type = "error"
    )
  }
  max_edges_per_target <- suppressWarnings(as.numeric(max_edges_per_target))
  if (length(max_edges_per_target) == 1L && is.finite(max_edges_per_target)) {
    max_edges_per_target <- max(1L, as.integer(max_edges_per_target))
    out <- out[
      order(out[["target"]], -out[["importance"]], out[["TF"]]), ,
      drop = FALSE
    ]
    out <- do.call(
      rbind,
      lapply(split(out, out[["target"]]), utils::head, max_edges_per_target)
    )
  }
  first_cols <- c("TF", "target", "importance")
  out <- out[,
    unique(c(first_cols, setdiff(colnames(out), first_cols))),
    drop = FALSE
  ]
  rownames(out) <- NULL
  out
}

scenicplus_filter_tri_prior <- function(
  triplets_prior,
  tf_gene,
  region_gene,
  tf_region
) {
  tf_gene_keys <- paste(tf_gene[["TF"]], tf_gene[["target"]], sep = "\r")
  region_gene_keys <- paste(
    region_gene[["region"]],
    region_gene[["gene"]],
    sep = "\r"
  )
  tf_region_keys <- paste(tf_region[["TF"]], tf_region[["region"]], sep = "\r")
  out <- triplets_prior[
    paste(triplets_prior[["TF"]], triplets_prior[["gene"]], sep = "\r") %in%
      tf_gene_keys &
      paste(
        triplets_prior[["region"]],
        triplets_prior[["gene"]],
        sep = "\r"
      ) %in%
        region_gene_keys &
      paste(triplets_prior[["TF"]], triplets_prior[["region"]], sep = "\r") %in%
        tf_region_keys &
      is.finite(triplets_prior[["score"]]),
    c("TF", "region", "gene", "score"),
    drop = FALSE
  ]
  if (nrow(out) == 0L) {
    log_message(
      "{.arg triplets_prior} has no triplets matching the selected evidence layers",
      message_type = "error"
    )
  }
  out <- out[order(out[["TF"]], out[["gene"]], out[["region"]]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

scenicplus_prepare_auc_expr <- function(auc_expr, fallback, cells) {
  if (is.null(auc_expr)) {
    return(fallback)
  }
  if (!is.matrix(auc_expr) && !inherits(auc_expr, "Matrix")) {
    auc_expr <- as_matrix(auc_expr)
  }
  auc_expr <- methods::as(Matrix::Matrix(auc_expr, sparse = TRUE), "dgCMatrix")
  if (is.null(rownames(auc_expr)) || is.null(colnames(auc_expr))) {
    log_message(
      "{.arg auc_expr} must have gene and cell dimnames",
      message_type = "error"
    )
  }
  col_hits <- intersect(colnames(auc_expr), cells)
  row_hits <- intersect(rownames(auc_expr), cells)
  if (length(col_hits) == 0L && length(row_hits) > 0L) {
    auc_expr <- Matrix::t(auc_expr)
    col_hits <- intersect(colnames(auc_expr), cells)
  }
  if (length(col_hits) == 0L) {
    log_message(
      "{.arg auc_expr} has no cells matching the RNA/ATAC assays",
      message_type = "error"
    )
  }
  auc_expr <- auc_expr[, col_hits, drop = FALSE]
  genes <- intersect(rownames(auc_expr), rownames(fallback))
  if (length(genes) == 0L) {
    log_message(
      "{.arg auc_expr} has no genes matching the RNA assay",
      message_type = "error"
    )
  }
  auc_expr[genes, , drop = FALSE]
}

scenicplus_prep_auc_rank <- function(auc_rankings, auc_counts) {
  if (is.null(auc_rankings)) {
    return(NULL)
  }
  if (!is.matrix(auc_rankings)) {
    auc_rankings <- as_matrix(auc_rankings)
  }
  if (is.null(rownames(auc_rankings)) || is.null(colnames(auc_rankings))) {
    log_message(
      "{.arg auc_rankings} must have cell and gene dimnames",
      message_type = "error"
    )
  }
  cell_hits <- intersect(rownames(auc_rankings), colnames(auc_counts))
  transposed_cell_hits <- intersect(
    colnames(auc_rankings),
    colnames(auc_counts)
  )
  if (length(cell_hits) == 0L && length(transposed_cell_hits) > 0L) {
    auc_rankings <- t(auc_rankings)
    cell_hits <- intersect(rownames(auc_rankings), colnames(auc_counts))
  }
  if (length(cell_hits) == 0L) {
    log_message(
      "{.arg auc_rankings} has no cells matching the AUC expression matrix",
      message_type = "error"
    )
  }
  genes <- intersect(colnames(auc_rankings), rownames(auc_counts))
  if (length(genes) == 0L) {
    log_message(
      "{.arg auc_rankings} has no genes matching the AUC expression matrix",
      message_type = "error"
    )
  }
  auc_rankings[cell_hits, genes, drop = FALSE]
}

scenicplus_ereg_prior <- function(
  eregulons_prior,
  triplets,
  min_size = 5
) {
  out <- eregulons_prior[
    eregulons_prior[["target"]] %in% triplets[["gene"]],
    c("regulon", "target"),
    drop = FALSE
  ]
  if (nrow(out) == 0L) {
    log_message(
      "{.arg eregulons_prior} has no targets matching assembled triplets",
      message_type = "error"
    )
  }
  regulon_list <- lapply(split(out[["target"]], out[["regulon"]]), unique)
  regulon_list <- regulon_list[lengths(regulon_list) >= min_size]
  if (length(regulon_list) == 0L) {
    log_message(
      "No eRegulons remain after applying {.arg min_eregulon_size}",
      message_type = "error"
    )
  }
  out <- out[out[["regulon"]] %in% names(regulon_list), , drop = FALSE]
  counts <- stats::aggregate(
    out[["target"]],
    list(regulon = out[["regulon"]]),
    function(x) length(unique(x))
  )
  colnames(counts) <- c("regulon", "target_count")
  out <- merge(out, counts, by = "regulon", all.x = TRUE)
  out <- out[order(out[["regulon"]], out[["target"]]), , drop = FALSE]
  rownames(out) <- NULL
  list(regulon_list = regulon_list, eregulons = out)
}

scenicplus_normalize_auc_table <- function(auc, cells) {
  first_col <- colnames(auc)[[1L]]
  if (
    !is.null(first_col) &&
      !first_col %in% cells &&
      !is.numeric(auc[[first_col]])
  ) {
    rownames(auc) <- as.character(auc[[first_col]])
    auc[[first_col]] <- NULL
  }
  auc <- as.data.frame(auc, check.names = FALSE, stringsAsFactors = FALSE)
  auc[] <- lapply(auc, as.numeric)
  common_cells <- intersect(cells, rownames(auc))
  if (length(common_cells) == 0L && length(cells) == nrow(auc)) {
    rownames(auc) <- cells
    common_cells <- cells
  }
  if (length(common_cells) == 0L) {
    log_message(
      "SCENIC+ AUC table must have cell row names matching {.arg srt}",
      message_type = "error"
    )
  }
  auc[common_cells, , drop = FALSE]
}

scenicplus_regulon_tbl <- function(
  auc,
  triplets,
  eregulons = NULL
) {
  regulons <- colnames(auc)
  if (
    !is.null(eregulons) && all(c("regulon", "target") %in% colnames(eregulons))
  ) {
    out <- lapply(regulons, function(regulon) {
      unique(eregulons[eregulons[["regulon"]] == regulon, "target"])
    })
    names(out) <- regulons
    return(out)
  }
  if (!"TF" %in% colnames(triplets) || !"gene" %in% colnames(triplets)) {
    return(stats::setNames(vector("list", length(regulons)), regulons))
  }
  genes_by_tf <- lapply(split(triplets, triplets[["TF"]]), function(df) {
    unique(df[["gene"]])
  })
  out <- lapply(regulons, function(regulon) {
    tf <- sub("\\(.*$", "", regulon)
    genes_by_tf[[tf]] %||% character()
  })
  names(out) <- regulons
  out
}

scenicplus_calc_rss <- function(auc, srt, cells, group.by = NULL) {
  empty <- list(
    rss = data.frame(
      regulon = character(),
      group = character(),
      rss = numeric(),
      stringsAsFactors = FALSE
    ),
    rss_matrix = NULL
  )
  if (is.null(group.by)) {
    return(empty)
  }
  if (!group.by %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg group.by} column {.val {group.by}} is not present in {.arg srt}",
      message_type = "error"
    )
  }
  groups <- as.character(srt@meta.data[cells, group.by, drop = TRUE])
  auc_mat <- t(as_matrix(auc[cells, , drop = FALSE]))
  rss_matrix <- scenic_calc_rss_matrix(
    auc_mat = auc_mat,
    cell_annotation = groups
  )
  rss <- as.data.frame(as.table(rss_matrix), stringsAsFactors = FALSE)
  colnames(rss) <- c("regulon", "group", "rss")
  list(rss = rss, rss_matrix = rss_matrix)
}

scenicplus_tf_gene_native <- function(
  rna_counts,
  regulators,
  targets = NULL,
  grn_method = c("grnboost2", "regdiffusion", "genie3", "gniplr"),
  max_edges_per_target = Inf,
  n_rounds = 5000,
  learning_rate = 0.01,
  max_depth = 3,
  max_features = 0.1,
  subsample = 0.9,
  early_stop_window_length = 25,
  seed = 666,
  cores = 1,
  envname = "scenic_env",
  conda = "auto",
  verbose = TRUE
) {
  grn_method <- match.arg(grn_method)
  adjacency <- RunGRN(
    object = Matrix::t(rna_counts),
    regulators = regulators,
    targets = targets,
    genes_in = "columns",
    grn_method = grn_method,
    backend = "cpp",
    output_file = tempfile(pattern = "scenicplus_tf_gene_", fileext = ".tsv"),
    max_edges_per_target = max_edges_per_target,
    n_rounds = n_rounds,
    learning_rate = learning_rate,
    max_depth = max_depth,
    max_features = max_features,
    subsample = subsample,
    early_stop_window_length = early_stop_window_length,
    seed = seed,
    exclude_self = FALSE,
    cores = cores,
    envname = envname,
    conda = conda,
    force = TRUE,
    verbose = verbose
  )
  if (!identical(grn_method, "grnboost2")) {
    return(adjacency)
  }
  adjacency <- scenicplus_add_tfg_cor(
    adjacency,
    rna_counts = rna_counts
  )
  scenicplus_inject_tf_self(adjacency, rna_counts = rna_counts)
}

scenicplus_add_tfg_cor <- function(
  adjacency,
  rna_counts,
  rho_threshold = 0.03
) {
  if (nrow(adjacency) == 0L) {
    return(adjacency)
  }
  genes <- intersect(
    unique(c(adjacency[["TF"]], adjacency[["target"]])),
    rownames(rna_counts)
  )
  if (length(genes) == 0L) {
    adjacency[["regulation"]] <- NA_integer_
    adjacency[["rho"]] <- NA_real_
    return(adjacency)
  }
  expr <- as_matrix(rna_counts[genes, , drop = FALSE])
  corr <- suppressWarnings(stats::cor(
    t(expr),
    method = "pearson",
    use = "pairwise.complete.obs"
  ))
  rho <- vapply(
    seq_len(nrow(adjacency)),
    function(i) {
      tf <- adjacency[["TF"]][[i]]
      target <- adjacency[["target"]][[i]]
      if (!tf %in% rownames(corr) || !target %in% colnames(corr)) {
        return(NA_real_)
      }
      as.numeric(corr[tf, target])
    },
    numeric(1)
  )
  adjacency[["regulation"]] <- ifelse(
    is.finite(rho) & rho > rho_threshold,
    1L,
    ifelse(is.finite(rho) & rho < -rho_threshold, -1L, 0L)
  )
  adjacency[["rho"]] <- rho
  adjacency[["importance_x_rho"]] <- as.numeric(adjacency[["importance"]]) * rho
  adjacency[["importance_x_abs_rho"]] <- as.numeric(adjacency[["importance"]]) *
    abs(rho)
  adjacency
}

scenicplus_inject_tf_self <- function(
  adjacency,
  rna_counts,
  increase_importance_by = 0.00001
) {
  if (nrow(adjacency) == 0L) {
    return(adjacency)
  }
  tfs <- intersect(unique(adjacency[["TF"]]), rownames(rna_counts))
  if (length(tfs) == 0L) {
    return(adjacency)
  }
  max_importance <- tapply(
    as.numeric(adjacency[["importance"]]),
    adjacency[["TF"]],
    max,
    na.rm = TRUE
  )
  self <- data.frame(
    TF = tfs,
    target = tfs,
    importance = as.numeric(max_importance[tfs]) + increase_importance_by,
    stringsAsFactors = FALSE
  )
  self <- scenicplus_add_tfg_cor(self, rna_counts = rna_counts)
  out <- rbind(adjacency, self[, colnames(adjacency), drop = FALSE])
  rownames(out) <- NULL
  out
}

scenicplus_region_gene_native <- function(
  srt,
  rna_counts,
  atac_counts,
  atac_assay,
  window = 250000,
  max_region_gene = 5,
  search_space = NULL,
  method = c("gbm", "correlation"),
  n_rounds = 500,
  max_features = 0.1
) {
  method <- match.arg(method)
  hits <- if (!is.null(search_space)) {
    scenicplus_rg_search(
      search_space = search_space,
      rna_counts = rna_counts,
      atac_counts = atac_counts
    )
  } else {
    peaks <- scenicplus_peak_ranges(rownames(atac_counts))
    genes <- scenicplus_gene_ranges(srt, atac_assay, rownames(rna_counts))
    scenicplus_peak_gene_hits(peaks, genes, window = window)
  }
  if (nrow(hits) == 0L) {
    log_message(
      "No peak-gene pairs are available within {.arg region_gene_window}",
      message_type = "error"
    )
  }
  if (identical(method, "gbm")) {
    return(scenicplus_rg_gbm(
      hits = hits,
      rna_counts = rna_counts,
      atac_counts = atac_counts,
      max_region_gene = max_region_gene,
      n_rounds = n_rounds,
      max_features = max_features
    ))
  }
  hit_regions <- unique(hits[["region"]])
  hit_genes <- unique(hits[["gene"]])
  atac_mat <- as_matrix(atac_counts[hit_regions, , drop = FALSE])
  rna_mat <- as_matrix(rna_counts[hit_genes, , drop = FALSE])
  hits[["score"]] <- scenicplus_region_gene_cor(
    atac_log = atac_mat,
    rna_log = rna_mat,
    region_idx = as.integer(match(hits[["region"]], hit_regions)),
    gene_idx = as.integer(match(hits[["gene"]], hit_genes))
  )
  hits <- hits[is.finite(hits[["score"]]), , drop = FALSE]
  if (nrow(hits) == 0L) {
    log_message(
      "No finite peak-gene correlations are available for {.fn RunSCENICPlus}",
      message_type = "error"
    )
  }
  hits <- hits[
    order(hits[["gene"]], -abs(hits[["score"]]), hits[["region"]]), ,
    drop = FALSE
  ]
  hits <- do.call(
    rbind,
    lapply(split(hits, hits[["gene"]]), utils::head, max_region_gene)
  )
  rownames(hits) <- NULL
  hits[, c("region", "gene", "score"), drop = FALSE]
}

scenicplus_rg_search <- function(
  search_space,
  rna_counts,
  atac_counts
) {
  search_space <- scenicplus_std_search(search_space)
  out <- search_space[
    search_space[["region"]] %in%
      rownames(atac_counts) &
      search_space[["gene"]] %in% rownames(rna_counts), ,
    drop = FALSE
  ]
  if (nrow(out) == 0L) {
    log_message(
      "{.arg region_gene_search_space} has no candidate pairs matching RNA/ATAC assays",
      message_type = "error"
    )
  }
  out[, intersect(c("region", "gene", "distance"), colnames(out)), drop = FALSE]
}

scenicplus_std_search <- function(search_space) {
  if (is.character(search_space) && length(search_space) == 1L) {
    if (!file.exists(search_space)) {
      log_message(
        "{.arg region_gene_search_space} file does not exist: {.file {search_space}}",
        message_type = "error"
      )
    }
    search_space <- scenicplus_read_delim(search_space)
  }
  if (!is.data.frame(search_space)) {
    log_message(
      "{.arg region_gene_search_space} must be a data frame or delimited file path",
      message_type = "error"
    )
  }
  region_col <- intersect(
    c("region", "Region", "Name", "peak"),
    colnames(search_space)
  )
  gene_col <- intersect(
    c("gene", "Gene", "target", "Target"),
    colnames(search_space)
  )
  if (length(region_col) == 0L || length(gene_col) == 0L) {
    log_message(
      "{.arg region_gene_search_space} must contain region/Name and gene/Gene columns",
      message_type = "error"
    )
  }
  out <- data.frame(
    region = as.character(search_space[[region_col[[1L]]]]),
    gene = as.character(search_space[[gene_col[[1L]]]]),
    stringsAsFactors = FALSE
  )
  out[["region"]] <- sub(":", "-", out[["region"]], fixed = TRUE)
  distance_col <- intersect(c("distance", "Distance"), colnames(search_space))
  if (length(distance_col) > 0L) {
    out[["distance"]] <- search_space[[distance_col[[1L]]]]
  }
  out <- out[stats::complete.cases(out[c("region", "gene")]), , drop = FALSE]
  unique(out)
}

scenicplus_rg_gbm <- function(
  hits,
  rna_counts,
  atac_counts,
  max_region_gene = 5,
  n_rounds = 500,
  max_features = 0.1
) {
  hits <- hits[
    hits[["region"]] %in%
      rownames(atac_counts) &
      hits[["gene"]] %in% rownames(rna_counts), ,
    drop = FALSE
  ]
  if (nrow(hits) == 0L) {
    log_message(
      "No region-gene candidates remain for native GBM scoring",
      message_type = "error"
    )
  }
  atac_mat <- atac_counts[unique(hits[["region"]]), , drop = FALSE]
  rna_mat <- rna_counts[unique(hits[["gene"]]), , drop = FALSE]
  rows <- lapply(split(hits, hits[["gene"]]), function(df) {
    regions <- unique(df[["region"]])
    gene <- df[["gene"]][[1L]]
    regions <- intersect(regions, rownames(atac_mat))
    if (length(regions) == 0L || !gene %in% rownames(rna_mat)) {
      return(NULL)
    }
    expr <- cbind(
      as_matrix(Matrix::t(atac_mat[regions, , drop = FALSE])),
      target = as.numeric(rna_mat[gene, ])
    )
    edge_idx <- grnboost_tree(
      expr = expr,
      regulator_idx = seq_along(regions),
      target_idx = length(regions) + 1L,
      n_rounds = as.integer(max(1L, n_rounds)),
      learning_rate = 0.01,
      max_edges_per_target = 0L,
      max_depth = 3L,
      max_features = max_features,
      subsample = 1,
      early_stop_window_length = 25L,
      random_seed = 666L
    )
    if (nrow(edge_idx) == 0L) {
      importance <- stats::setNames(rep(0, length(regions)), regions)
    } else {
      importance <- stats::setNames(
        as.numeric(edge_idx[["importance"]]),
        regions[edge_idx[["regulator"]]]
      )
      importance <- importance[regions]
      importance[is.na(importance)] <- 0
      total <- sum(importance)
      if (is.finite(total) && total > 0) {
        importance <- importance / total
      }
    }
    rho <- vapply(
      regions,
      function(region) {
        suppressWarnings(stats::cor(
          as.numeric(atac_mat[region, ]),
          as.numeric(rna_mat[gene, ]),
          method = "spearman",
          use = "complete.obs"
        ))
      },
      numeric(1)
    )
    score <- as.numeric(importance) * abs(as.numeric(rho))
    out <- data.frame(
      region = regions,
      gene = gene,
      score = score,
      importance = as.numeric(importance),
      rho = as.numeric(rho),
      stringsAsFactors = FALSE
    )
    out[is.finite(out[["score"]]), , drop = FALSE]
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out) || nrow(out) == 0L) {
    log_message(
      "Native GBM region-gene scoring returned no finite links",
      message_type = "error"
    )
  }
  max_region_gene <- suppressWarnings(as.numeric(max_region_gene))
  if (length(max_region_gene) == 1L && is.finite(max_region_gene)) {
    max_region_gene <- max(1L, as.integer(max_region_gene))
    out <- out[
      order(out[["gene"]], -out[["score"]], out[["region"]]), ,
      drop = FALSE
    ]
    out <- do.call(
      rbind,
      lapply(split(out, out[["gene"]]), utils::head, max_region_gene)
    )
  }
  rownames(out) <- NULL
  out
}

scenicplus_select_grn_targets <- function(
  region_gene,
  rna_counts,
  max_grn_targets = Inf,
  target_scope = c("all", "region_gene")
) {
  target_scope <- match.arg(target_scope)
  genes <- if (identical(target_scope, "all")) {
    rownames(rna_counts)
  } else {
    intersect(unique(region_gene[["gene"]]), rownames(rna_counts))
  }
  if (length(genes) == 0L) {
    return(character())
  }
  max_grn_targets <- suppressWarnings(as.numeric(max_grn_targets))
  if (
    length(max_grn_targets) != 1L ||
      is.na(max_grn_targets) ||
      is.infinite(max_grn_targets)
  ) {
    return(genes)
  }
  max_grn_targets <- max(1L, as.integer(max_grn_targets))
  if (length(genes) <= max_grn_targets) {
    return(genes)
  }
  if (identical(target_scope, "all")) {
    return(genes[seq_len(min(max_grn_targets, length(genes)))])
  }
  score_by_gene <- tapply(
    abs(region_gene[["score"]]),
    region_gene[["gene"]],
    max,
    na.rm = TRUE
  )
  score_by_gene <- score_by_gene[intersect(names(score_by_gene), genes)]
  score_by_gene[!is.finite(score_by_gene)] <- 0
  names(sort(score_by_gene, decreasing = TRUE))[seq_len(min(
    max_grn_targets,
    length(score_by_gene)
  ))]
}

scenicplus_tfr_motifs <- function(
  srt,
  atac_assay,
  peak_names,
  max_tf_region = 500
) {
  motif_mat <- tryCatch(
    Signac::GetMotifData(object = srt[[atac_assay]], slot = "data"),
    error = function(...) NULL
  )
  if (is.null(motif_mat) || nrow(motif_mat) == 0L || ncol(motif_mat) == 0L) {
    log_message(
      "{.fn RunSCENICPlus} native backend requires existing motif annotations in the chromatin assay; de novo motif scanning is not performed.",
      message_type = "error"
    )
  }

  peak_rows <- intersect(rownames(motif_mat), peak_names)
  peak_cols <- intersect(colnames(motif_mat), peak_names)
  if (length(peak_rows) > 0L) {
    motif_mat <- motif_mat[peak_rows, , drop = FALSE]
    regions <- rownames(motif_mat)
    motifs <- colnames(motif_mat)
    motif_by_region <- Matrix::t(motif_mat)
  } else if (length(peak_cols) > 0L) {
    motif_mat <- motif_mat[, peak_cols, drop = FALSE]
    regions <- colnames(motif_mat)
    motifs <- rownames(motif_mat)
    motif_by_region <- motif_mat
  } else {
    log_message(
      "No motif peaks overlap the ATAC assay peaks",
      message_type = "error"
    )
  }

  motif_names <- tryCatch(
    Signac::GetMotifData(object = srt[[atac_assay]], slot = "motif.names"),
    error = function(...) NULL
  )
  if (is.null(motif_names) || length(motif_names) == 0L) {
    tf_names <- as.character(motifs)
  } else {
    motif_names <- unlist(motif_names, use.names = TRUE)
    tf_names <- motif_names[as.character(motifs)]
    tf_names[is.na(tf_names) | !nzchar(tf_names)] <- as.character(motifs)[
      is.na(tf_names) | !nzchar(tf_names)
    ]
    tf_names <- as.character(tf_names)
  }
  rownames(motif_by_region) <- tf_names
  colnames(motif_by_region) <- regions

  motif_by_region <- Matrix::Matrix(motif_by_region, sparse = TRUE)
  nz <- as.data.frame(Matrix::summary(motif_by_region))
  if (nrow(nz) == 0L) {
    log_message(
      "No TF-region motif links are available",
      message_type = "error"
    )
  }
  colnames(nz)[1:3] <- c("tf_idx", "region_idx", "score")
  nz <- nz[is.finite(nz[["score"]]) & nz[["score"]] > 0, , drop = FALSE]
  if (nrow(nz) == 0L) {
    log_message(
      "No positive TF-region motif links are available",
      message_type = "error"
    )
  }
  nz <- nz[
    order(nz[["tf_idx"]], -nz[["score"]], nz[["region_idx"]]), ,
    drop = FALSE
  ]
  nz <- do.call(
    rbind,
    lapply(split(nz, nz[["tf_idx"]]), utils::head, max_tf_region)
  )
  links <- data.frame(
    TF = tf_names[nz[["tf_idx"]]],
    region = colnames(motif_by_region)[nz[["region_idx"]]],
    score = nz[["score"]],
    stringsAsFactors = FALSE
  )
  if (is.null(links) || nrow(links) == 0L) {
    log_message(
      "No TF-region motif links are available",
      message_type = "error"
    )
  }
  rownames(links) <- NULL
  links
}

scenicplus_peak_ranges <- function(peak_names) {
  peaks <- tryCatch(
    Signac::StringToGRanges(peak_names, sep = c("-", "-")),
    error = function(err) {
      log_message(
        "Cannot parse ATAC peak names as genomic ranges: {conditionMessage(err)}",
        message_type = "error"
      )
    }
  )
  names(peaks) <- peak_names
  peaks
}

scenicplus_gene_ranges <- function(srt, atac_assay, genes_use) {
  check_r(c("GenomicRanges", "IRanges"), verbose = FALSE)
  annotation <- Signac::Annotation(srt[[atac_assay]])
  if (is.null(annotation) || length(annotation) == 0L) {
    log_message(
      "{.fn RunSCENICPlus} requires gene annotations on the chromatin assay",
      message_type = "error"
    )
  }
  gene_col <- intersect(
    c("gene_name", "gene", "symbol", "external_gene_name"),
    colnames(S4Vectors::mcols(annotation))
  )
  if (length(gene_col) == 0L) {
    log_message(
      "Chromatin assay annotation must contain a gene name column",
      message_type = "error"
    )
  }
  gene_names <- as.character(S4Vectors::mcols(annotation)[[gene_col[[1L]]]])
  keep <- gene_names %in% genes_use
  annotation <- annotation[keep]
  gene_names <- gene_names[keep]
  if (length(annotation) == 0L) {
    log_message(
      "No RNA genes overlap chromatin assay annotations",
      message_type = "error"
    )
  }
  starts <- ifelse(
    as.character(GenomicRanges::strand(annotation)) == "-",
    GenomicRanges::end(annotation),
    GenomicRanges::start(annotation)
  )
  genes <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(annotation),
    ranges = get_namespace_fun("IRanges", "IRanges")(start = starts, width = 1),
    strand = GenomicRanges::strand(annotation)
  )
  names(genes) <- gene_names
  genes <- genes[!duplicated(names(genes))]
  genes
}

scenicplus_peak_gene_hits <- function(peaks, genes, window = 250000) {
  genes_window <- GenomicRanges::resize(genes, width = 1)
  genes_window <- suppressWarnings(GenomicRanges::flank(
    genes_window,
    width = window,
    both = TRUE
  ))
  hits <- GenomicRanges::findOverlaps(peaks, genes_window, ignore.strand = TRUE)
  data.frame(
    region = names(peaks)[S4Vectors::queryHits(hits)],
    gene = names(genes)[S4Vectors::subjectHits(hits)],
    stringsAsFactors = FALSE
  )
}

scenicplus_triplets <- function(tf_gene, region_gene, tf_region) {
  out <- scenicplus_triplets_cpp(
    tf_gene_tf = as.character(tf_gene[["TF"]]),
    tf_gene_target = as.character(tf_gene[["target"]]),
    tf_gene_importance = as.numeric(tf_gene[["importance"]]),
    region_gene_region = as.character(region_gene[["region"]]),
    region_gene_gene = as.character(region_gene[["gene"]]),
    region_gene_score = as.numeric(region_gene[["score"]]),
    tf_region_tf = as.character(tf_region[["TF"]]),
    tf_region_region = as.character(tf_region[["region"]]),
    tf_region_score = as.numeric(tf_region[["score"]])
  )
  if (is.null(out) || nrow(out) == 0L) {
    log_message(
      "No SCENICPlus-style TF-region-gene triplets could be assembled",
      message_type = "error"
    )
  }
  rownames(out) <- NULL
  out
}

scenicplus_eregulon_build <- function(
  triplets,
  tf_gene = NULL,
  region_gene = NULL,
  rna_counts = NULL,
  min_size = 5,
  rho_threshold = 0.05,
  quantiles = c(0.85, 0.90, 0.95),
  top_n_region_gene = c(5L, 10L, 15L),
  min_target_genes = 10L
) {
  triplets <- scenicplus_add_eregulon_signs(
    triplets = triplets,
    tf_gene = tf_gene,
    region_gene = region_gene,
    rna_counts = rna_counts
  )
  filtered_triplets <- scenicplus_le_triplets(
    triplets = triplets,
    tf_gene = tf_gene,
    region_gene = region_gene,
    rho_threshold = rho_threshold,
    quantiles = quantiles,
    top_n_region_gene = top_n_region_gene,
    min_target_genes = min_target_genes
  )
  if (nrow(filtered_triplets) == 0L) {
    filtered_triplets <- triplets
  }
  split_key <- paste(
    filtered_triplets[["TF"]],
    filtered_triplets[["tf_sign"]],
    filtered_triplets[["r2g_sign"]],
    sep = "\r"
  )
  genes_by_key <- lapply(split(filtered_triplets, split_key), function(df) {
    unique(df[["gene"]])
  })
  genes_by_key <- genes_by_key[lengths(genes_by_key) >= min_size]
  if (length(genes_by_key) == 0L) {
    return(list(
      triplets = filtered_triplets[
        FALSE,
        c("TF", "region", "gene", "score"),
        drop = FALSE
      ],
      regulon_list = genes_by_key
    ))
  }
  retained_keys <- names(genes_by_key)
  names(genes_by_key) <- vapply(
    names(genes_by_key),
    function(key) {
      parts <- strsplit(key, "\r", fixed = TRUE)[[1L]]
      genes <- genes_by_key[[key]]
      sprintf(
        "%s_direct_%s/%s_(%dg)",
        parts[[1L]],
        parts[[2L]],
        parts[[3L]],
        length(genes)
      )
    },
    character(1)
  )
  keep_keys <- paste(
    filtered_triplets[["TF"]],
    filtered_triplets[["tf_sign"]],
    filtered_triplets[["r2g_sign"]],
    sep = "\r"
  ) %in%
    retained_keys
  filtered_triplets <- unique(filtered_triplets[
    keep_keys,
    c("TF", "region", "gene", "score"),
    drop = FALSE
  ])
  rownames(filtered_triplets) <- NULL
  list(triplets = filtered_triplets, regulon_list = genes_by_key)
}

scenicplus_le_triplets <- function(
  triplets,
  tf_gene,
  region_gene,
  rho_threshold = 0.05,
  quantiles = c(0.85, 0.90, 0.95),
  top_n_region_gene = c(5L, 10L, 15L),
  min_target_genes = 10L
) {
  if (
    is.null(tf_gene) ||
      is.null(region_gene) ||
      !"rho" %in% colnames(tf_gene) ||
      !"rho" %in% colnames(region_gene)
  ) {
    return(triplets)
  }
  region_gene <- region_gene[is.finite(region_gene[["rho"]]), , drop = FALSE]
  if (nrow(region_gene) == 0L) {
    return(triplets)
  }
  if (!"importance" %in% colnames(region_gene)) {
    region_gene[["importance"]] <- abs(as.numeric(region_gene[["score"]]))
  }
  tf_gene <- tf_gene[
    is.finite(tf_gene[["rho"]]) & is.finite(tf_gene[["importance"]]), ,
    drop = FALSE
  ]
  tf_regions <- split(triplets[["region"]], triplets[["TF"]])
  tf_pos <- split(
    tf_gene[tf_gene[["rho"]] > rho_threshold, , drop = FALSE],
    tf_gene[tf_gene[["rho"]] > rho_threshold, "TF"]
  )
  tf_neg <- split(
    tf_gene[tf_gene[["rho"]] < -rho_threshold, , drop = FALSE],
    tf_gene[tf_gene[["rho"]] < -rho_threshold, "TF"]
  )
  tf_pos <- tf_pos[vapply(tf_pos, nrow, integer(1)) >= min_target_genes]
  tf_neg <- tf_neg[vapply(tf_neg, nrow, integer(1)) >= min_target_genes]
  if (length(tf_pos) == 0L && length(tf_neg) == 0L) {
    return(triplets)
  }

  modules <- scenicplus_region_gene_modules(
    region_gene = region_gene,
    rho_threshold = rho_threshold,
    quantiles = quantiles,
    top_n_region_gene = top_n_region_gene
  )
  if (length(modules) == 0L) {
    return(triplets)
  }

  out <- vector("list", 0L)
  for (tf in intersect(names(tf_regions), unique(tf_gene[["TF"]]))) {
    regions <- unique(tf_regions[[tf]])
    for (module in modules) {
      emodule <- module[["data"]][
        module[["data"]][["region"]] %in% regions, ,
        drop = FALSE
      ]
      if (nrow(emodule) == 0L) {
        next
      }
      gene_set <- unique(emodule[["gene"]])
      if (tf %in% names(tf_pos)) {
        ranking <- stats::setNames(
          tf_pos[[tf]][["importance"]],
          tf_pos[[tf]][["target"]]
        )
        leading_edge <- scenicplus_gsea_leading_edge(ranking, gene_set)
        if (length(leading_edge) >= min_target_genes) {
          out[[length(out) + 1L]] <- scenicplus_fmt_le_tri(
            tf,
            "+",
            module[["sign"]],
            emodule,
            leading_edge
          )
        }
      }
      if (tf %in% names(tf_neg)) {
        ranking <- stats::setNames(
          tf_neg[[tf]][["importance"]],
          tf_neg[[tf]][["target"]]
        )
        leading_edge <- scenicplus_gsea_leading_edge(ranking, gene_set)
        if (length(leading_edge) >= min_target_genes) {
          out[[length(out) + 1L]] <- scenicplus_fmt_le_tri(
            tf,
            "-",
            module[["sign"]],
            emodule,
            leading_edge
          )
        }
      }
    }
  }
  if (length(out) == 0L) {
    return(triplets[FALSE, , drop = FALSE])
  }
  out <- unique(do.call(rbind, out))
  out <- out[
    order(
      out[["TF"]],
      out[["tf_sign"]],
      out[["r2g_sign"]],
      out[["gene"]],
      out[["region"]]
    ), ,
    drop = FALSE
  ]
  rownames(out) <- NULL
  out
}

scenicplus_region_gene_modules <- function(
  region_gene,
  rho_threshold = 0.05,
  quantiles = c(0.85, 0.90, 0.95),
  top_n_region_gene = c(5L, 10L, 15L)
) {
  modules <- vector("list", 0L)
  for (r2g_sign in c("-", "+")) {
    adj <- if (identical(r2g_sign, "+")) {
      region_gene[region_gene[["rho"]] > rho_threshold, , drop = FALSE]
    } else {
      region_gene[region_gene[["rho"]] < -rho_threshold, , drop = FALSE]
    }
    if (nrow(adj) == 0L) {
      next
    }
    for (threshold in quantiles) {
      passing <- as.logical(stats::ave(
        adj[["importance"]],
        adj[["gene"]],
        FUN = function(x) {
          x >
            as.numeric(stats::quantile(
              x,
              probs = threshold,
              names = FALSE,
              type = 7
            ))
        }
      ))
      if (any(passing)) {
        modules[[length(modules) + 1L]] <- list(
          sign = r2g_sign,
          data = adj[passing, , drop = FALSE]
        )
      }
    }
    for (n in top_n_region_gene) {
      n <- as.integer(n)
      if (!is.finite(n) || n < 1L) {
        next
      }
      passing <- as.logical(stats::ave(
        adj[["importance"]],
        adj[["gene"]],
        FUN = function(x) {
          x >= if (length(x) >= n) sort(x, decreasing = TRUE)[[n]] else min(x)
        }
      ))
      if (any(passing)) {
        modules[[length(modules) + 1L]] <- list(
          sign = r2g_sign,
          data = adj[passing, , drop = FALSE]
        )
      }
    }
  }
  modules
}

scenicplus_gsea_leading_edge <- function(ranking, gene_set) {
  ranking <- ranking[is.finite(ranking)]
  ranking <- sort(ranking, decreasing = TRUE)
  if (length(ranking) < 2L) {
    return(character())
  }
  hits <- names(ranking) %in% unique(gene_set)
  n_hits <- sum(hits)
  if (n_hits == 0L || n_hits == length(ranking)) {
    return(character())
  }
  hit_weights <- abs(as.numeric(ranking))
  hit_weights[!hits] <- 0
  total_hit_weight <- sum(hit_weights[hits])
  if (!is.finite(total_hit_weight) || total_hit_weight <= 0) {
    hit_weights[hits] <- 1 / n_hits
  } else {
    hit_weights[hits] <- hit_weights[hits] / total_hit_weight
  }
  miss_weights <- ifelse(hits, 0, 1 / (length(ranking) - n_hits))
  running_score <- cumsum(hit_weights - miss_weights)
  max_score <- max(running_score)
  min_score <- min(running_score)
  if (!is.finite(max_score) || max_score <= 0 || max_score < abs(min_score)) {
    return(character())
  }
  max_idx <- which.max(running_score)
  unique(names(ranking)[hits & seq_along(hits) <= max_idx])
}

scenicplus_fmt_le_tri <- function(
  tf,
  tf_sign,
  r2g_sign,
  emodule,
  leading_edge
) {
  out <- emodule[
    emodule[["gene"]] %in% leading_edge,
    c("region", "gene", "score"),
    drop = FALSE
  ]
  if (nrow(out) == 0L) {
    return(out)
  }
  data.frame(
    TF = tf,
    region = out[["region"]],
    gene = out[["gene"]],
    score = abs(as.numeric(out[["score"]])),
    tf_sign = tf_sign,
    r2g_sign = r2g_sign,
    stringsAsFactors = FALSE
  )
}

scenicplus_add_eregulon_signs <- function(
  triplets,
  tf_gene = NULL,
  region_gene = NULL,
  rna_counts = NULL
) {
  triplets[["tf_sign"]] <- "+"
  triplets[["r2g_sign"]] <- "+"
  if (
    !is.null(region_gene) &&
      all(c("region", "gene", "rho") %in% colnames(region_gene))
  ) {
    r2g_rho <- stats::setNames(
      as.numeric(region_gene[["rho"]]),
      paste(region_gene[["region"]], region_gene[["gene"]], sep = "\r")
    )
    rho <- r2g_rho[paste(triplets[["region"]], triplets[["gene"]], sep = "\r")]
    triplets[["r2g_sign"]] <- ifelse(is.finite(rho) & rho < 0, "-", "+")
  } else if (
    !is.null(region_gene) &&
      all(c("region", "gene", "score") %in% colnames(region_gene))
  ) {
    score <- stats::setNames(
      as.numeric(region_gene[["score"]]),
      paste(region_gene[["region"]], region_gene[["gene"]], sep = "\r")
    )
    r2g_score <- score[paste(
      triplets[["region"]],
      triplets[["gene"]],
      sep = "\r"
    )]
    triplets[["r2g_sign"]] <- ifelse(
      is.finite(r2g_score) & r2g_score < 0,
      "-",
      "+"
    )
  }
  tf_rho <- NULL
  if (
    !is.null(tf_gene) && all(c("TF", "target", "rho") %in% colnames(tf_gene))
  ) {
    tf_rho <- stats::setNames(
      as.numeric(tf_gene[["rho"]]),
      paste(tf_gene[["TF"]], tf_gene[["target"]], sep = "\r")
    )
  } else if (!is.null(rna_counts)) {
    tf_rho <- scenicplus_tf_gene_rho(
      pairs = unique(triplets[c("TF", "gene")]),
      rna_counts = rna_counts
    )
  }
  if (!is.null(tf_rho)) {
    rho <- tf_rho[paste(triplets[["TF"]], triplets[["gene"]], sep = "\r")]
    triplets[["tf_sign"]] <- ifelse(is.finite(rho) & rho < 0, "-", "+")
  }
  triplets
}

scenicplus_tf_gene_rho <- function(pairs, rna_counts) {
  pairs <- pairs[
    pairs[["TF"]] %in%
      rownames(rna_counts) &
      pairs[["gene"]] %in% rownames(rna_counts), ,
    drop = FALSE
  ]
  if (nrow(pairs) == 0L) {
    return(NULL)
  }
  expr <- log1p(as_matrix(rna_counts[
    unique(c(pairs[["TF"]], pairs[["gene"]])), ,
    drop = FALSE
  ]))
  out <- vapply(
    seq_len(nrow(pairs)),
    function(i) {
      suppressWarnings(stats::cor(
        as.numeric(expr[pairs[["TF"]][[i]], ]),
        as.numeric(expr[pairs[["gene"]][[i]], ]),
        method = "spearman",
        use = "complete.obs"
      ))
    },
    numeric(1)
  )
  stats::setNames(out, paste(pairs[["TF"]], pairs[["gene"]], sep = "\r"))
}
