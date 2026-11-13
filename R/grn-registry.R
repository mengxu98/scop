grn_method_spec <- function(method, backend = NULL) {
  key <- tolower(gsub("[^a-zA-Z0-9+]", "", method))
  aliases <- c("scenic+" = "scenicplus", "multicsn" = "multicsn", "ppcor" = "ppcor")
  if (key %in% names(aliases)) key <- unname(aliases[[key]])
  legacy_args <- c("max_edges_per_target", "n_rounds", "learning_rate", "max_depth", "max_features", "subsample", "early_stop_window_length", "exclude_self", "correlation_threshold", "lasso_degree", "lasso_alpha", "max_lag")
  registry <- list(
    grnboost2 = list(backends = c("cpp", "python", "native", "official"), default = "cpp", runner = "grn_run_rna", args = legacy_args, profile = "grnboost2"),
    regdiffusion = list(backends = c("python", "native"), default = "python", runner = "grn_run_rna", args = c("n_iter", "batch_size", "lr", "device", "max_edges_per_target"), profile = "regdiffusion"),
    genie3 = list(backends = c("cpp", "r", "native"), default = "cpp", runner = "grn_run_legacy", args = c("nTrees", "treeMethod", "K", "max_edges_per_target"), profile = "r"),
    gniplr = list(backends = c("cpp", "python", "native", "official"), default = "cpp", runner = "grn_run_rna", args = legacy_args, profile = "scientific"),
    ppcor = list(backends = c("r", "native"), default = "r", runner = "grn_run_rna", args = c("method"), profile = "r"),
    pidc = list(backends = c("r", "native"), default = "r", runner = "grn_run_rna", args = c("log_scale", "verbose"), profile = "r"),
    scorpion = list(backends = c("r", "native"), default = "r", runner = "grn_run_rna", args = c("n_pc", "gamma_value", "show_progress"), profile = "r"),
    scsgl = list(backends = c("python", "native"), default = "python", runner = "grn_run_rna", args = character(), profile = "scsgl"),
    leap = list(backends = c("r", "native"), default = "r", runner = "grn_run_time", args = character(), profile = "r", time = TRUE),
    scode = list(backends = c("r", "native"), default = "r", runner = "grn_run_time", args = c("latent_dim", "n_iter", "n_runs"), profile = "r", time = TRUE),
    sincerities = list(backends = c("r", "native"), default = "r", runner = "grn_run_time", args = c("sign", "n_bins", "min_cells", "smoothing", "lambda"), profile = "r", time = TRUE),
    multicsn = list(backends = "r", default = "r", runner = "grn_run_multicsn", prepare = "grn_prepare_multicsn", args = c("preset", "min_detected", "min_peak_cells", "upstream", "downstream", "tf_region_scope", "exclude_response_alias", "renormalize", "response_chunk", "target_chunk", "max_support_size", "min_improvement", "sort_regulators", "motif_backend", "exclude_exons", "align_annotation", "method", "lag_fraction", "lag_steps"), profile = "r", modes = c("rna", "multiome")),
    pando = list(backends = c("r", "native"), default = "r", runner = "grn_run_pando", args = character(), profile = "pando", modes = "multiome"),
    linger = list(backends = c("python", "native"), default = "python", runner = "grn_run_linger", args = character(), profile = "linger", modes = "multiome"),
    scenicplus = list(backends = c("python", "official", "native"), default = "python", runner = "grn_run_scenicplus", args = character(), profile = "scenicplus", modes = "multiome"),
    scenic = list(backends = c("python", "cpp"), default = "python", runner = "grn_run_scenic", args = character(), profile = "scenic")
  )
  if (length(key) != 1L || !key %in% names(registry)) stop("Unknown GRN method. Choose one of: ", paste(names(registry), collapse = ", "), call. = FALSE)
  spec <- registry[[key]]
  spec$method <- key
  spec$backend <- backend %||% spec$default
  if (length(spec$backend) != 1L || !spec$backend %in% spec$backends) stop("Unsupported backend for ", key, ". Choose: ", paste(spec$backends, collapse = ", "), call. = FALSE)
  spec$allowed_args <- if (identical(spec$runner, "grn_run_rna") && exists("grn_rna_method_args", mode = "function")) grn_rna_method_args(key) else spec$args
  if (key == "leap") spec$allowed_args <- c("max_lag_prop", "MAC_cutoff")
  if (key == "scode") spec$allowed_args <- c("trajectory", "latent_dim", "n_iter", "n_runs")
  if (key == "sincerities") spec$allowed_args <- c("trajectory", "max_time_points", "method", "noDIAG", "SIGN", "normalize")
  if (key == "linger") spec$allowed_args <- c("sample.by", "label.by", "min_genes", "min_cells", "genome", "activation", "patch_missing_keys", "resume_dir")
  if (key == "pando") spec$allowed_args <- c("exclude_exons", "nfeatures", "peak_to_gene_method", "upstream", "downstream", "extend", "only_tss", "tf_cor", "peak_cor", "family", "interaction_term", "adjust_method", "scale", "p_thresh", "rsq_thresh", "nvar_thresh", "min_genes_per_module")
  if (key == "scenicplus") spec$allowed_args <- c("gsea_n_perm", "min_target_genes", "data_preparation", "motif_enrichment", "inference", "resume_dir")
  spec$candidate_filter <- if (key == "scenicplus") "postfit_view" else "method_specific"
  if (key == "scenic") spec$allowed_args <- setdiff(names(formals(RunSCENIC)), c("object", "srt", "assay", "layer", "regulators", "targets", "work_dir", "backend", "cores", "seed", "return_seurat", "envname", "conda", "verbose", "ranking_dbs", "motif_annotations"))
  spec$allow_cell_filter <- identical(key, "linger")
  spec$input_modes <- spec$modes %||% "rna"
  if (key %in% c("grnboost2", "gniplr") && identical(spec$backend, "cpp")) {
    spec$runner <- "grn_run_legacy"
    spec$allowed_args <- legacy_args
    spec$profile <- "r"
  }
  if (key == "gniplr" && spec$runner != "grn_run_legacy") spec$profile <- "gniplr"
  if (key %in% c("genie3", "ppcor", "pidc", "scorpion", "leap", "scode", "sincerities", "multicsn")) spec$profile <- key
  spec
}
