#' Run mcRigor metacell partition assessment
#'
#' @description
#' `RunmcRigor()` wraps the upstream `JSB-UCLA/mcRigor` package to detect
#' dubious metacells for a supplied partition or to optimize across multiple
#' candidate metacell partitions. Results are stored in `srt@tools[[tool_name]]`
#' and the selected partition/status are written back to cell metadata.
#'
#' The upstream mcRigor package is installed at runtime when missing and is not
#' bundled with `scop`.
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param srt A `Seurat` object containing the original single-cell data.
#' @param cell_membership A data frame or matrix with cells in rows and one or
#' more metacell partitions in columns. Row names should be cell names. If row
#' names are missing and the row count equals `ncol(srt)`, cells are matched in
#' `colnames(srt)` order.
#' @param metacell.by Metadata column(s) in `srt` used as metacell partitions
#' when `cell_membership = NULL`.
#' @param mode mcRigor task. `"detect"` calls `mcRigor_DETECT()` for one
#' partition; `"optimize"` calls `mcRigor_OPTIMIZE()` across candidate
#' partitions.
#' @param tgamma Target partition/gamma for `"detect"`. Can be a membership
#' column name or the numeric gamma label used by mcRigor. If `NULL`, the first
#' membership column is used.
#' @param gamma_names Optional gamma labels for membership columns. mcRigor
#' requires numeric-like column labels; non-numeric labels are mapped internally
#' to `1:ncol(cell_membership)` and recorded in the stored result.
#' @param assay_type Assay type passed to mcRigor.
#' @param Gammas Candidate gamma labels for `"optimize"`. Can use original
#' membership column names or mapped mcRigor gamma labels.
#' @param aggregate_method Metacell aggregation method passed to mcRigor.
#' @param output_file Optional path where mcRigor writes the `TabMC` RDS file.
#' If `NULL`, a temporary file is used to avoid creating files in the working
#' directory.
#' @param Nrep Number of permutation repetitions used by mcRigor.
#' @param gene_filter,feature_use,cor_method,prePro,test_cutoff,thre_smooth,thre_bw Parameters forwarded to mcRigor.
#' @param D_bw,optim_method,weight,dub_rate Optimization parameters forwarded
#' to `mcRigor_OPTIMIZE()`.
#' @param draw,pur_metric,check_purity,fields,step_save Plotting, purity, and
#' intermediate-save parameters forwarded to mcRigor.
#' @param prefix Prefix for metadata columns written to `srt`.
#' @param tool_name Name of the `srt@tools` entry used to store results.
#'
#' @return A `Seurat` object with mcRigor metadata and a result list stored in
#' `srt@tools[[tool_name]]`.
#' @export
#'
#' @references
#' Liu, P. and Li, J.J. (2024). mcRigor: a statistical method to enhance the
#' rigor of metacell partitioning in single-cell data analysis.
#' \emph{bioRxiv}. \doi{10.1101/2024.10.30.621093}
#'
#' @examples
#' data(pancreas_sub)
#' set.seed(11)
#' pancreas_sub <- pancreas_sub[, seq_len(200)]
#' pancreas_sub <- standard_scop(
#'   pancreas_sub,
#'   nHVF = 500,
#'   linear_reduction_dims = 20,
#'   linear_reduction_dims_use = 1:20,
#'   nonlinear_reduction_dims = 2,
#'   verbose = FALSE
#' )
#' pancreas_sub <- RunMetaCell(
#'   pancreas_sub,
#'   method = "supercell",
#'   gamma = 25
#' )
#'
#' pancreas_sub <- RunmcRigor(
#'   pancreas_sub,
#'   metacell.by = "Metacell_id",
#'   Nrep = 1,
#'   feature_use = 100,
#'   draw = FALSE
#' )
#'
#' table(pancreas_sub$mcRigor_status)
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "Metacell_id"
#' )
#'
#' CellDimPlot(
#'   pancreas_sub,
#'   group.by = "mcRigor_status"
#' )
RunmcRigor <- function(
  srt,
  cell_membership = NULL,
  metacell.by = NULL,
  mode = c("detect", "optimize"),
  tgamma = NULL,
  gamma_names = NULL,
  assay_type = c("RNA", "ATAC"),
  Gammas = NULL,
  aggregate_method = c("mean", "sum", "geom"),
  output_file = NULL,
  Nrep = 1,
  gene_filter = 0.1,
  feature_use = 2000,
  cor_method = c("pearson", "spearman"),
  prePro = TRUE,
  test_cutoff = 0.01,
  thre_smooth = TRUE,
  thre_bw = 1 / 6,
  D_bw = 10,
  optim_method = c("tradeoff", "dub_rate_large", "dub_rate_small"),
  weight = 0.5,
  dub_rate = 0.1,
  draw = FALSE,
  pur_metric = NULL,
  check_purity = TRUE,
  fields = NULL,
  step_save = FALSE,
  prefix = "mcRigor",
  tool_name = "mcRigor",
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  mode <- match.arg(mode)
  assay_type <- match.arg(assay_type)
  aggregate_method <- match.arg(aggregate_method)
  cor_method <- match.arg(cor_method)
  optim_method <- match.arg(optim_method)

  for (name in c("prefix", "tool_name")) {
    value <- get(name)
    if (!is.character(value) || length(value) != 1L || !nzchar(value)) {
      log_message(
        "{.arg {name}} must be a single non-empty string",
        message_type = "error"
      )
    }
  }
  if (
    !is.numeric(Nrep) ||
      length(Nrep) != 1L ||
      is.na(Nrep) ||
      Nrep < 1 ||
      Nrep != as.integer(Nrep)
  ) {
    log_message(
      "{.arg Nrep} must be a single integer greater than or equal to 1",
      message_type = "error"
    )
  }
  if (
    !is.numeric(feature_use) ||
      length(feature_use) != 1L ||
      is.na(feature_use) ||
      feature_use < 1 ||
      feature_use != as.integer(feature_use)
  ) {
    log_message(
      "{.arg feature_use} must be a single integer greater than or equal to 1",
      message_type = "error"
    )
  }
  numeric_unit_params <- list(
    gene_filter = gene_filter,
    test_cutoff = test_cutoff,
    thre_bw = thre_bw,
    weight = weight,
    dub_rate = dub_rate
  )
  for (param_name in names(numeric_unit_params)) {
    value <- numeric_unit_params[[param_name]]
    if (
      !is.numeric(value) ||
        length(value) != 1L ||
        is.na(value) ||
        value < 0 ||
        value > 1
    ) {
      log_message(
        "{.arg {param_name}} must be a single number between 0 and 1",
        message_type = "error"
      )
    }
  }

  membership_raw <- if (is.null(cell_membership)) {
    if (is.null(metacell.by)) {
      log_message(
        "Provide {.arg cell_membership} or {.arg metacell.by}",
        message_type = "error"
      )
    }
    if (!all(metacell.by %in% colnames(srt[[]]))) {
      missing_columns <- setdiff(metacell.by, colnames(srt[[]]))
      log_message(
        "Missing metadata column{?s}: {.val {missing_columns}}",
        message_type = "error"
      )
    }
    as.data.frame(srt[[]][, metacell.by, drop = FALSE], check.names = FALSE)
  } else {
    as.data.frame(cell_membership, check.names = FALSE)
  }

  if (ncol(membership_raw) == 0L || nrow(membership_raw) == 0L) {
    log_message(
      "{.arg cell_membership} must contain at least one partition column",
      message_type = "error"
    )
  }
  default_rownames <- identical(
    rownames(membership_raw),
    as.character(seq_len(nrow(membership_raw)))
  )
  if (
    is.null(rownames(membership_raw)) ||
      any(!nzchar(rownames(membership_raw))) ||
      default_rownames
  ) {
    if (nrow(membership_raw) == ncol(srt)) {
      rownames(membership_raw) <- colnames(srt)
    } else {
      log_message(
        "{.arg cell_membership} must have cell names as row names",
        message_type = "error"
      )
    }
  }

  if (anyDuplicated(rownames(membership_raw)) > 0L) {
    log_message(
      "{.arg cell_membership} row names must be unique cell names",
      message_type = "error"
    )
  }
  common_cells <- intersect(colnames(srt), rownames(membership_raw))
  if (length(common_cells) == 0L) {
    log_message(
      "No cells are shared between {.arg srt} and {.arg cell_membership}",
      message_type = "error"
    )
  }
  if (length(common_cells) < ncol(srt)) {
    log_message(
      "Run {.pkg mcRigor} on {.val {length(common_cells)}} cells with supplied membership; unmatched cells receive {.val NA} metadata",
      message_type = "warning",
      verbose = verbose
    )
  }
  membership <- membership_raw[common_cells, , drop = FALSE]
  srt_use <- srt[, common_cells]

  original_names <- colnames(membership)
  if (is.null(original_names) || any(!nzchar(original_names))) {
    original_names <- paste0("partition_", seq_len(ncol(membership)))
  }
  if (!is.null(gamma_names)) {
    gamma_names <- as.character(gamma_names)
    if (
      !is.character(gamma_names) ||
        length(gamma_names) != ncol(membership) ||
        any(!nzchar(gamma_names))
    ) {
      log_message(
        "{.arg gamma_names} must be a character vector with one label per membership column",
        message_type = "error"
      )
    }
    original_names <- gamma_names
  }
  if (anyDuplicated(original_names) > 0L) {
    log_message(
      "{.arg gamma_names} and membership column names must be unique",
      message_type = "error"
    )
  }

  numeric_gamma <- suppressWarnings(as.numeric(original_names))
  use_original_gamma <- all(!is.na(numeric_gamma)) &&
    anyDuplicated(as.character(numeric_gamma)) == 0L
  mc_gamma <- if (use_original_gamma) {
    as.character(numeric_gamma)
  } else {
    as.character(seq_len(ncol(membership)))
  }
  colnames(membership) <- mc_gamma
  gamma_map <- data.frame(
    membership_column = original_names,
    mcRigor_gamma = mc_gamma,
    stringsAsFactors = FALSE
  )

  resolve_gamma <- function(x, arg_name) {
    if (is.null(x)) {
      return(NULL)
    }
    x <- as.character(x)
    out <- gamma_map$mcRigor_gamma[match(x, gamma_map$membership_column)]
    direct <- x %in% gamma_map$mcRigor_gamma
    out[direct] <- x[direct]
    if (any(is.na(out))) {
      log_message(
        "{.arg {arg_name}} must match membership columns or mcRigor gamma labels",
        message_type = "error"
      )
    }
    out
  }
  tgamma_use <- resolve_gamma(tgamma, "tgamma")
  Gammas_use <- resolve_gamma(Gammas, "Gammas")
  output_file_use <- output_file %||%
    tempfile(pattern = paste0(prefix, "_TabMC_"), fileext = ".rds")

  check_r("JSB-UCLA/mcRigor", verbose = FALSE)
  log_message(
    "Run {.pkg mcRigor} in {.val {mode}} mode with {.val {ncol(membership)}} partition{?s}",
    verbose = verbose
  )

  result <- if (identical(mode, "detect")) {
    get_namespace_fun("mcRigor", "mcRigor_DETECT")(
      obj_singlecell = srt_use,
      cell_membership = membership,
      tgamma = tgamma_use,
      assay_type = assay_type,
      aggregate_method = aggregate_method,
      output_file = output_file_use,
      Nrep = as.integer(Nrep),
      gene_filter = gene_filter,
      feature_use = as.integer(feature_use),
      cor_method = cor_method,
      prePro = isTRUE(prePro),
      test_cutoff = test_cutoff,
      thre_smooth = isTRUE(thre_smooth),
      thre_bw = thre_bw,
      draw = isTRUE(draw),
      pur_metric = pur_metric,
      check_purity = isTRUE(check_purity),
      fields = fields,
      step_save = isTRUE(step_save)
    )
  } else {
    get_namespace_fun("mcRigor", "mcRigor_OPTIMIZE")(
      obj_singlecell = srt_use,
      cell_membership = membership,
      assay_type = assay_type,
      Gammas = Gammas_use,
      aggregate_method = aggregate_method,
      output_file = output_file_use,
      Nrep = as.integer(Nrep),
      gene_filter = gene_filter,
      feature_use = as.integer(feature_use),
      cor_method = cor_method,
      prePro = isTRUE(prePro),
      test_cutoff = test_cutoff,
      thre_smooth = isTRUE(thre_smooth),
      thre_bw = thre_bw,
      D_bw = D_bw,
      optim_method = optim_method,
      weight = weight,
      dub_rate = dub_rate,
      draw = isTRUE(draw),
      pur_metric = pur_metric,
      check_purity = isTRUE(check_purity),
      fields = fields,
      step_save = isTRUE(step_save)
    )
  }

  selected_gamma <- if (identical(mode, "detect")) {
    tgamma_use %||% colnames(membership)[[1]]
  } else {
    as.character(result$best_granularity_level)
  }
  selected_column <- gamma_map$membership_column[
    match(selected_gamma, gamma_map$mcRigor_gamma)
  ]
  selected_label <- if (length(selected_column) == 1L && !is.na(selected_column)) {
    selected_column
  } else {
    selected_gamma
  }
  metacell_object <- result$obj_metacell %||% result$opt_metacell
  cell_status <- NULL
  if (!is.null(metacell_object@misc$cell_membership$mcRigor_sc)) {
    cell_status <- metacell_object@misc$cell_membership$mcRigor_sc
    names(cell_status) <- rownames(metacell_object@misc$cell_membership)
  }

  metadata <- data.frame(row.names = colnames(srt))
  metadata[[paste0(prefix, "_metacell")]] <- NA_character_
  metadata[[paste0(prefix, "_status")]] <- NA_character_
  metadata[[paste0(prefix, "_gamma")]] <- NA_character_
  selected_membership <- as.character(membership[[selected_gamma]])
  names(selected_membership) <- rownames(membership)
  metadata[names(selected_membership), paste0(prefix, "_metacell")] <-
    selected_membership
  metadata[names(selected_membership), paste0(prefix, "_gamma")] <-
    selected_label
  if (!is.null(cell_status)) {
    status_cells <- intersect(names(cell_status), rownames(metadata))
    metadata[status_cells, paste0(prefix, "_status")] <-
      as.character(cell_status[status_cells])
  }
  srt <- Seurat::AddMetaData(srt, metadata = metadata)

  srt@tools[[tool_name]] <- list(
    mode = mode,
    result = result,
    selected_gamma = selected_gamma,
    selected_column = selected_label,
    gamma_map = gamma_map,
    output_file = output_file_use,
    output_file_is_temporary = is.null(output_file),
    parameters = list(
      assay_type = assay_type,
      aggregate_method = aggregate_method,
      Nrep = as.integer(Nrep),
      gene_filter = gene_filter,
      feature_use = as.integer(feature_use),
      cor_method = cor_method,
      prePro = isTRUE(prePro),
      test_cutoff = test_cutoff,
      thre_smooth = isTRUE(thre_smooth),
      thre_bw = thre_bw,
      D_bw = D_bw,
      optim_method = optim_method,
      weight = weight,
      dub_rate = dub_rate,
      draw = isTRUE(draw),
      pur_metric = pur_metric,
      check_purity = isTRUE(check_purity),
      fields = fields,
      step_save = isTRUE(step_save),
      prefix = prefix
    )
  )

  log_message(
    "{.pkg mcRigor} results stored in {.code srt@tools[[{tool_name}]]}",
    message_type = "success",
    verbose = verbose
  )
  srt
}
