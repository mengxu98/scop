#' @title Run copy-number alteration inference
#'
#' @description
#' Run expression-based single-cell or spatial copy-number alteration backends
#' and store the results in a unified SCOP schema.
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param srt A `Seurat` object.
#' @param method CNA/CNV backend. Supported backends are `"copykat"`,
#' `"fastCNV"`, `"scevan"`, `"infercnv"`, `"numbat"`, and `"casper"`.
#' @param layer Assay layer used as the expression matrix.
#' @param group.by Optional metadata column forwarded to supported backends and
#' stored as cell annotation.
#' @param reference.by Metadata column identifying reference/normal cells.
#' Required for `"infercnv"` and `"fastCNV"`.
#' @param reference Reference labels in `reference.by`.
#' @param genome Reference genome label.
#' @param gene_order Gene coordinate table or a path to one. The table should
#' contain gene, chromosome, start, and end columns. If `NULL`, SCOP tries to
#' resolve these columns from assay feature metadata.
#' @param sample.by Optional sample metadata column.
#' @param allele_counts Allele count table for `"numbat"`. This is forwarded
#' to `numbat::run_numbat()` as `df_allele`.
#' @param reference_counts Reference expression profile for `"numbat"`. This
#' is forwarded to `numbat::run_numbat()` as `lambdas_ref`.
#' @param loh B-allele frequency/LOH signal for `"casper"`.
#' @param loh_name_mapping Optional CaSpER LOH-to-cell mapping table.
#' @param cytoband Cytoband table for `"casper"`.
#' @param output_dir Optional backend output directory.
#' @param prefix Prefix for metadata columns.
#' @param tool_name Name used for `srt@tools`.
#' @param store_matrix Whether to store the normalized CNV matrix in
#' `srt@tools[[tool_name]]`.
#' @param ... Additional parameters forwarded to the selected backend.
#'
#' @return A `Seurat` object with CNV metadata columns and a result bundle in
#' `srt@tools[[tool_name]]`.
#' @export
#'
#' @seealso [CNVPlot]
#'
#' @examples
#' \dontrun{
#' # copykat uses raw counts and can infer diploid/aneuploid cells directly.
#' srt <- RunCNV(
#'   srt,
#'   method = "copykat",
#'   genome = "hg38"
#' )
#'
#' # fastCNV and inferCNV require a normal/reference cell annotation.
#' srt <- RunCNV(
#'   srt,
#'   method = "fastCNV",
#'   reference.by = "celltype",
#'   reference = "Normal",
#'   genome = "hg38"
#' )
#'
#' gene_order <- data.frame(
#'   gene = rownames(srt),
#'   chr = "chr1",
#'   start = seq_len(nrow(srt)) * 1000,
#'   end = seq_len(nrow(srt)) * 1000 + 999
#' )
#' srt <- RunCNV(
#'   srt,
#'   method = "infercnv",
#'   reference.by = "celltype",
#'   reference = "Normal",
#'   gene_order = gene_order
#' )
#'
#' # Numbat and CaSpER can also be run when allele-aware preprocessing
#' # outputs are available.
#'
#' CNVPlot(srt, plot_type = "heatmap", group.by = "CNV_prediction")
#' CNVPlot(srt, plot_type = "dim", value = "CNV_prediction")
#' }
RunCNV <- function(
  srt,
  method = c("copykat", "fastCNV", "scevan", "infercnv", "numbat", "casper"),
  assay = NULL,
  layer = "counts",
  group.by = NULL,
  reference.by = NULL,
  reference = NULL,
  genome = c("hg38", "hg19", "mm10"),
  gene_order = NULL,
  sample.by = NULL,
  allele_counts = NULL,
  reference_counts = NULL,
  loh = NULL,
  loh_name_mapping = NULL,
  cytoband = NULL,
  output_dir = NULL,
  prefix = "CNV",
  tool_name = "CNV",
  store_matrix = TRUE,
  verbose = TRUE,
  ...
) {
  cnv_assert_seurat(srt)
  method <- cnv_match_method(method)
  genome <- match.arg(genome)
  cnv_assert_scalar_string(prefix, "prefix")
  cnv_assert_scalar_string(tool_name, "tool_name")

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  if (!assay %in% names(srt@assays)) {
    log_message(
      "{.arg assay} {.val {assay}} is not present in {.arg srt}",
      message_type = "error"
    )
  }
  cnv_validate_metadata(srt, group.by, "group.by")
  cnv_validate_metadata(srt, reference.by, "reference.by")
  cnv_validate_metadata(srt, sample.by, "sample.by")
  reference_cells <- cnv_reference_cells(
    srt = srt,
    method = method,
    reference.by = reference.by,
    reference = reference
  )
  gene_order_tbl <- cnv_resolve_gene_order(
    gene_order = gene_order,
    srt = srt,
    assay = assay,
    method = method,
    required = method %in% c("infercnv", "casper")
  )
  cnv_validate_backend_inputs(
    method = method,
    allele_counts = allele_counts,
    reference_counts = reference_counts,
    loh = loh,
    cytoband = cytoband,
    gene_order = gene_order_tbl
  )
  counts <- cnv_get_counts(
    srt = srt,
    assay = assay,
    layer = layer
  )

  if (!is.null(output_dir)) {
    output_dir <- normalizePath(
      path.expand(output_dir),
      mustWork = FALSE,
      winslash = "/"
    )
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
  }

  log_message(
    "Run {.pkg {method}} CNV backend with {.val {nrow(counts)}} features and {.val {ncol(counts)}} cells",
    verbose = verbose
  )
  backend_result <- cnv_run_backend(
    method = method,
    srt = srt,
    counts = counts,
    assay = assay,
    layer = layer,
    group.by = group.by,
    reference.by = reference.by,
    reference = reference,
    reference_cells = reference_cells,
    genome = genome,
    gene_order = gene_order_tbl,
    sample.by = sample.by,
    allele_counts = allele_counts,
    reference_counts = reference_counts,
    loh = loh,
    loh_name_mapping = loh_name_mapping,
    cytoband = cytoband,
    output_dir = output_dir,
    verbose = verbose,
    ...
  )
  method_bundle <- cnv_standardize_result(
    result = backend_result,
    method = method,
    cells = colnames(srt),
    features = rownames(counts),
    reference_cells = reference_cells,
    store_matrix = store_matrix,
    parameters = list(
      method = method,
      assay = assay,
      layer = layer,
      group.by = group.by,
      reference.by = reference.by,
      reference = reference,
      genome = genome,
      sample.by = sample.by,
      allele_counts = !is.null(allele_counts),
      reference_counts = !is.null(reference_counts),
      loh = !is.null(loh),
      cytoband = !is.null(cytoband),
      output_dir = output_dir,
      prefix = prefix,
      tool_name = tool_name
    )
  )

  srt <- cnv_add_metadata(
    srt = srt,
    bundle = method_bundle,
    method = method,
    prefix = prefix,
    tool_name = tool_name
  )

  old_store <- srt@tools[[tool_name]] %||% list()
  methods_store <- old_store$methods %||% list()
  methods_store[[method]] <- method_bundle
  srt@tools[[tool_name]] <- list(
    active_method = method,
    methods = methods_store,
    metadata = list(
      schema = "scop_cnv_v1",
      supported_methods = c("copykat", "fastCNV", "scevan", "infercnv", "numbat", "casper"),
      allele_aware_methods = c("numbat", "casper")
    )
  )
  srt <- Seurat::LogSeuratCommand(srt)

  log_message(
    "{.pkg {method}} CNV results stored in {.code srt@tools[[{tool_name}]]}",
    message_type = "success",
    verbose = verbose
  )
  srt
}

cnv_run_backend <- function(
  method,
  srt,
  counts,
  assay,
  layer,
  group.by,
  reference.by,
  reference,
  reference_cells,
  genome,
  gene_order,
  sample.by,
  allele_counts,
  reference_counts,
  loh,
  loh_name_mapping,
  cytoband,
  output_dir,
  verbose,
  ...
) {
  switch(method,
    copykat = cnv_run_copykat(
      counts = counts,
      group.by = group.by,
      reference_cells = reference_cells,
      genome = genome,
      sample.by = sample.by,
      output_dir = output_dir,
      verbose = verbose,
      ...
    ),
    fastCNV = cnv_run_fastcnv(
      srt = srt,
      assay = assay,
      layer = layer,
      group.by = group.by,
      reference.by = reference.by,
      reference = reference,
      genome = genome,
      sample.by = sample.by,
      output_dir = output_dir,
      verbose = verbose,
      ...
    ),
    scevan = cnv_run_scevan(
      counts = counts,
      genome = genome,
      output_dir = output_dir,
      verbose = verbose,
      ...
    ),
    infercnv = cnv_run_infercnv(
      counts = counts,
      reference.by = reference.by,
      reference = reference,
      reference_cells = reference_cells,
      gene_order = gene_order,
      output_dir = output_dir,
      verbose = verbose,
      ...
    ),
    numbat = cnv_run_numbat(
      counts = counts,
      allele_counts = allele_counts,
      reference_counts = reference_counts,
      genome = genome,
      output_dir = output_dir,
      verbose = verbose,
      ...
    ),
    casper = cnv_run_casper(
      counts = counts,
      reference_cells = reference_cells,
      gene_order = gene_order,
      loh = loh,
      loh_name_mapping = loh_name_mapping,
      cytoband = cytoband,
      genome = genome,
      output_dir = output_dir,
      verbose = verbose,
      ...
    )
  )
}

cnv_run_copykat <- function(
  counts,
  group.by = NULL,
  reference_cells = NULL,
  genome = "hg38",
  sample.by = NULL,
  output_dir = NULL,
  verbose = TRUE,
  ...
) {
  check_r("copykat", verbose = FALSE)
  copykat_fun <- get_namespace_fun("copykat", "copykat")
  genome_use <- switch(genome,
    hg38 = "hg20",
    hg19 = "hg19",
    mm10 = "mm10",
    genome
  )
  args <- utils::modifyList(
    list(
      rawmat = as.matrix(counts),
      id.type = "S",
      sam.name = sample.by %||% "scop_cnv",
      genome = genome_use,
      n.cores = 1,
      norm.cell.names = reference_cells %||% "",
      plot.genes = FALSE,
      output.seg = FALSE
    ),
    list(...)
  )
  oldwd <- NULL
  if (!is.null(output_dir)) {
    oldwd <- getwd()
    setwd(output_dir)
    on.exit(setwd(oldwd), add = TRUE)
  }
  result <- cnv_call_backend_fun(copykat_fun, args)
  cnv_extract_copykat(result)
}

cnv_run_fastcnv <- function(
  srt,
  assay,
  layer,
  group.by = NULL,
  reference.by = NULL,
  reference = NULL,
  genome = "hg38",
  sample.by = NULL,
  output_dir = NULL,
  verbose = TRUE,
  ...
) {
  check_r("fastCNV", verbose = FALSE)
  if (identical(genome, "mm10")) {
    log_message(
      "{.pkg fastCNV} currently supports human data only; use {.arg genome = 'hg38'} or {.arg genome = 'hg19'}",
      message_type = "error"
    )
  }
  fastcnv_fun <- get_namespace_fun("fastCNV", "fastCNV")
  sample_input <- cnv_fastcnv_sample_input(srt = srt, sample.by = sample.by)
  args <- utils::modifyList(
    list(
      seuratObj = sample_input$seuratObj,
      sampleName = sample_input$sampleName,
      referenceVar = reference.by,
      referenceLabel = reference,
      assay = assay,
      prepareCounts = FALSE,
      doPlot = FALSE,
      savePath = output_dir
    ),
    list(...)
  )
  result <- cnv_call_backend_fun(fastcnv_fun, args)
  cnv_extract_fastcnv(
    result,
    cells = colnames(srt),
    reference.by = reference.by,
    reference = reference
  )
}

cnv_run_scevan <- function(
  counts,
  genome = "hg38",
  output_dir = NULL,
  verbose = TRUE,
  ...
) {
  check_r("SCEVAN", verbose = FALSE)
  sample_name <- "scop_cnv"
  scevan_fun <- cnv_get_first_namespace_fun(
    "SCEVAN",
    c("pipelineCNA", "SCEVAN")
  )
  args <- utils::modifyList(
    list(
      count_mtx = as.matrix(counts),
      sample = sample_name,
      par_cores = 1,
      SUBCLONES = TRUE,
      plotTree = FALSE,
      organism = if (identical(genome, "mm10")) "mouse" else "human",
      out_dir = output_dir
    ),
    list(...)
  )
  run_dir <- output_dir %||% getwd()
  oldwd <- NULL
  if (!is.null(output_dir)) {
    oldwd <- getwd()
    setwd(output_dir)
    on.exit(setwd(oldwd), add = TRUE)
  }
  result <- cnv_call_backend_fun(scevan_fun, args)
  cnv_extract_scevan(
    result = result,
    cells = colnames(counts),
    sample = sample_name,
    output_dir = run_dir
  )
}

cnv_run_infercnv <- function(
  counts,
  reference.by,
  reference,
  reference_cells,
  gene_order,
  output_dir = NULL,
  verbose = TRUE,
  ...
) {
  check_r("infercnv", verbose = FALSE)
  log_message(
    "{.pkg infercnv} upstream README states the project is no longer supported; consider {.pkg copykat}, {.pkg fastCNV}, or {.pkg Numbat} for new analyses",
    message_type = "warning",
    verbose = verbose
  )
  if (is.null(gene_order) || nrow(gene_order) == 0L) {
    log_message(
      "{.arg gene_order} is required for {.arg method = 'infercnv'}",
      message_type = "error"
    )
  }
  out_dir <- output_dir %||% tempfile("scop_infercnv_")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  annotations <- data.frame(
    cell = colnames(counts),
    group = ifelse(colnames(counts) %in% reference_cells, "reference", "observation"),
    stringsAsFactors = FALSE
  )
  ann_file <- tempfile("infercnv_annotations_", fileext = ".txt")
  gene_file <- tempfile("infercnv_gene_order_", fileext = ".txt")
  utils::write.table(
    annotations,
    ann_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  utils::write.table(
    gene_order[, c("gene", "chr", "start", "end"), drop = FALSE],
    gene_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  create_obj <- get_namespace_fun("infercnv", "CreateInfercnvObject")
  run_fun <- get_namespace_fun("infercnv", "run")
  infer_obj <- create_obj(
    raw_counts_matrix = as.matrix(counts),
    annotations_file = ann_file,
    delim = "\t",
    gene_order_file = gene_file,
    ref_group_names = "reference"
  )
  args <- utils::modifyList(
    list(
      infercnv_obj = infer_obj,
      cutoff = 0.1,
      out_dir = out_dir,
      cluster_by_groups = FALSE,
      cluster_references = FALSE,
      denoise = FALSE,
      HMM = FALSE,
      analysis_mode = "samples",
      tumor_subcluster_partition_method = "qnorm",
      plot_steps = FALSE,
      resume_mode = FALSE,
      no_plot = TRUE,
      no_prelim_plot = TRUE,
      save_rds = FALSE,
      save_final_rds = FALSE,
      plot_probabilities = FALSE,
      write_expr_matrix = FALSE,
      write_phylo = FALSE,
      inspect_subclusters = FALSE,
      num_threads = cnv_default_threads()
    ),
    list(...)
  )
  result <- cnv_call_backend_fun(run_fun, args)
  cnv_extract_generic_result(result, method = "infercnv", cells = colnames(counts))
}

cnv_run_numbat <- function(
  counts,
  allele_counts,
  reference_counts,
  genome = "hg38",
  output_dir = NULL,
  verbose = TRUE,
  ...
) {
  check_r("numbat", verbose = FALSE)
  out_dir <- output_dir %||% tempfile("scop_numbat_")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  run_fun <- get_namespace_fun("numbat", "run_numbat")
  args <- utils::modifyList(
    list(
      count_mat = counts,
      lambdas_ref = reference_counts,
      df_allele = allele_counts,
      genome = genome,
      out_dir = out_dir,
      ncores = 1,
      plot = FALSE,
      verbose = verbose
    ),
    list(...)
  )
  result <- cnv_call_backend_fun(run_fun, args)
  cnv_extract_numbat(result = result, cells = colnames(counts), output_dir = out_dir)
}

cnv_run_casper <- function(
  counts,
  reference_cells,
  gene_order,
  loh,
  loh_name_mapping = NULL,
  cytoband,
  genome = "hg38",
  output_dir = NULL,
  verbose = TRUE,
  ...
) {
  check_r("CaSpER", verbose = FALSE)
  create_fun <- get_namespace_fun("CaSpER", "CreateCasperObject")
  run_fun <- get_namespace_fun("CaSpER", "runCaSpER")
  event_fun <- get_namespace_fun("CaSpER", "extractLargeScaleEvents")
  annotation <- cnv_gene_order_to_casper_annotation(gene_order)
  casper_args <- utils::modifyList(
    list(
      raw.data = as.matrix(counts),
      annotation = annotation,
      control.sample.ids = reference_cells,
      cytoband = cytoband,
      loh.name.mapping = loh_name_mapping,
      cnv.scale = 3,
      loh.scale = 3,
      method = "iterative",
      loh = loh,
      project = "scop_cnv",
      matrix.type = "raw",
      sequencing.type = "single-cell",
      log.transformed = FALSE,
      genomeVersion = if (identical(genome, "hg19")) "hg19" else "hg38"
    ),
    list(...)
  )
  object <- cnv_call_backend_fun(create_fun, casper_args)
  run_args <- utils::modifyList(
    list(
      object = object,
      removeCentromere = TRUE,
      cytoband = cytoband,
      method = casper_args$method %||% "iterative"
    ),
    list(...)
  )
  final_objects <- cnv_call_backend_fun(run_fun, run_args)
  event_args <- list(final.objects = final_objects)
  final_mat <- tryCatch(
    cnv_call_backend_fun(event_fun, event_args),
    error = function(e) NULL
  )
  cnv_extract_casper(
    result = list(final_objects = final_objects, large_scale_events = final_mat, object = object),
    cells = colnames(counts)
  )
}

cnv_extract_copykat <- function(result) {
  cnv_matrix <- result$CNAmat %||% result$CNA %||% result$cnv_matrix
  if (is.null(cnv_matrix)) {
    log_message(
      "{.pkg copykat} did not return {.field CNAmat}",
      message_type = "error"
    )
  }
  cnv_matrix <- as.data.frame(cnv_matrix, check.names = FALSE)
  coord_cols <- intersect(
    c("chrom", "chr", "chromosome", "Chromosome", "start", "end", "abspos", "chrompos", "cytoband", "gene", "genes"),
    colnames(cnv_matrix)
  )
  bin_info <- cnv_bin_info_from_matrix(cnv_matrix, coord_cols = coord_cols)
  value <- as.matrix(cnv_matrix[, setdiff(colnames(cnv_matrix), coord_cols), drop = FALSE])
  storage.mode(value) <- "double"

  pred <- result$prediction %||% result$pred
  cell_info <- NULL
  if (!is.null(pred)) {
    pred <- as.data.frame(pred, check.names = FALSE)
    cell_col <- cnv_first_col(pred, c("cell.names", "cell", "barcode", "Cell"))
    pred_col <- cnv_first_col(pred, c("copykat.pred", "prediction", "pred", "status"))
    if (!is.null(cell_col) && !is.null(pred_col)) {
      cell_info <- data.frame(
        cell = as.character(pred[[cell_col]]),
        prediction = as.character(pred[[pred_col]]),
        stringsAsFactors = FALSE
      )
      rownames(cell_info) <- cell_info$cell
    }
  }

  list(
    cnv_matrix = value,
    bin_info = bin_info,
    cell_info = cell_info,
    raw = list(
      class = class(result),
      names = names(result),
      prediction = pred
    )
  )
}

cnv_extract_generic_result <- function(result, method, cells) {
  if (inherits(result, "Seurat")) {
    return(cnv_extract_from_seurat_result(result, method = method, cells = cells))
  }
  candidates <- cnv_find_matrix_candidates(result, cells = cells)
  mat <- candidates$matrix
  if (is.null(mat)) {
    log_message(
      "{.pkg {method}} did not return a detectable CNV matrix",
      message_type = "error"
    )
  }
  bin_info <- cnv_bin_info_from_matrix(mat)
  mat <- cnv_drop_coordinate_columns(mat)
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  cell_info <- cnv_find_cell_info(result, cells = cells)
  list(
    cnv_matrix = mat,
    bin_info = bin_info,
    cell_info = cell_info,
    raw = cnv_light_raw(result)
  )
}

cnv_extract_numbat <- function(result, cells, output_dir = NULL) {
  loaded <- cnv_read_numbat_output(output_dir = output_dir)
  merged <- if (is.list(result)) {
    utils::modifyList(loaded, result)
  } else {
    loaded
  }
  if (isS4(result) || inherits(result, "R6")) {
    merged <- c(loaded, cnv_object_fields(result))
  }

  mat <- NULL
  bin_info <- NULL
  joint_post <- merged$joint_post %||% merged$joint.post
  if (!is.null(joint_post)) {
    from_post <- cnv_matrix_from_numbat_joint_post(joint_post, cells = cells)
    mat <- from_post$matrix
    bin_info <- from_post$bin_info
  }
  if (is.null(mat)) {
    candidates <- cnv_find_matrix_candidates(merged, cells = cells)
    mat <- candidates$matrix
  }
  if (is.null(mat)) {
    log_message(
      "{.pkg numbat} did not return a detectable CNV matrix. Check {.arg output_dir} for joint_post_*.tsv files.",
      message_type = "error"
    )
  }
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  cell_info <- cnv_cell_info_from_numbat_clone_post(merged$clone_post %||% merged$clone.post, cells = cells)
  if (is.null(cell_info)) {
    cell_info <- cnv_find_cell_info(merged, cells = cells)
  }
  list(
    cnv_matrix = mat,
    bin_info = bin_info %||% cnv_bin_info_from_matrix(mat),
    cell_info = cell_info,
    raw = cnv_light_raw(merged)
  )
}

cnv_extract_casper <- function(result, cells) {
  mat <- result$large_scale_events %||% result$finalChrMat %||% result$final_chr_mat
  if (is.null(mat)) {
    mat <- cnv_casper_matrix_from_objects(result$final_objects %||% result$objects %||% result)
  }
  if (is.null(mat)) {
    candidates <- cnv_find_matrix_candidates(result, cells = cells)
    mat <- candidates$matrix
  }
  if (is.null(mat)) {
    log_message(
      "{.pkg CaSpER} did not return a detectable CNV matrix",
      message_type = "error"
    )
  }
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  mat <- cnv_orient_matrix(mat, cells = cells)
  bin_info <- cnv_casper_bin_info(mat)
  score <- colMeans(abs(mat), na.rm = TRUE)
  prediction <- ifelse(score > 0, "aneuploid", "diploid")
  cell_info <- data.frame(
    cell = colnames(mat),
    score = as.numeric(score),
    prediction = prediction,
    stringsAsFactors = FALSE
  )
  rownames(cell_info) <- cell_info$cell
  list(
    cnv_matrix = mat,
    bin_info = bin_info,
    cell_info = cell_info,
    raw = cnv_light_raw(result)
  )
}

cnv_extract_fastcnv <- function(result, cells, reference.by = NULL, reference = NULL) {
  if (inherits(result, "Seurat")) {
    return(cnv_extract_from_seurat_result(
      result,
      method = "fastCNV",
      cells = cells,
      reference.by = reference.by,
      reference = reference
    ))
  }
  if (is.list(result) && length(result) > 0L && all(vapply(result, inherits, logical(1), what = "Seurat"))) {
    parts <- lapply(
      result,
      cnv_extract_from_seurat_result,
      method = "fastCNV",
      cells = cells,
      reference.by = reference.by,
      reference = reference
    )
    return(cnv_combine_extracted_results(parts, cells = cells, method = "fastCNV"))
  }
  cnv_extract_generic_result(result, method = "fastCNV", cells = cells)
}

cnv_extract_scevan <- function(result, cells, sample = "scop_cnv", output_dir = getwd()) {
  candidates <- cnv_find_matrix_candidates(result, cells = cells)
  mat <- candidates$matrix
  bin_info <- NULL
  raw <- cnv_light_raw(result)

  if (is.null(mat)) {
    rdata <- cnv_read_scevan_output(output_dir = output_dir, sample = sample, cells = cells)
    mat <- rdata$matrix
    bin_info <- rdata$bin_info
    raw$scevan_files <- rdata$files
  }

  if (is.null(mat)) {
    log_message(
      "{.pkg SCEVAN} did not return a detectable CNV matrix. Check the backend output directory for a *_CNAmtx.RData file.",
      message_type = "error"
    )
  }

  bin_info <- bin_info %||% cnv_bin_info_from_matrix(mat)
  mat <- cnv_drop_coordinate_columns(mat)
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  cell_info <- cnv_find_cell_info(result, cells = cells)
  list(
    cnv_matrix = mat,
    bin_info = bin_info,
    cell_info = cell_info,
    raw = raw
  )
}

cnv_extract_from_seurat_result <- function(
  result,
  method,
  cells,
  reference.by = NULL,
  reference = NULL
) {
  meta <- result@meta.data
  method_key <- gsub("[^A-Za-z0-9]+", "_", method)
  score_col <- cnv_first_col(meta, c(
    paste0("CNV_", method_key, "_score"),
    paste0("CNV_", method_key, "_scores"),
    paste0("CNV_", method_key, "_fraction"),
    paste0("CNV_", method_key, "_fractions"),
    paste0(method_key, "_score"),
    paste0(method_key, "_scores"),
    paste0(method_key, "_fraction"),
    paste0(method_key, "_fractions"),
    "cnv_fraction",
    "cnv_fractions",
    "CNV_fraction",
    "CNV_fractions",
    "cnv_score",
    "CNV_score"
  ))
  pred_col <- cnv_first_col(meta, c(
    paste0("CNV_", method_key, "_prediction"),
    paste0("CNV_", method_key, "_predictions"),
    paste0(method_key, "_prediction"),
    paste0(method_key, "_predictions"),
    "cnv_prediction",
    "CNV_prediction",
    "prediction",
    "pred",
    "malignant",
    "cell.assignment",
    "cell_assignment"
  ))
  cluster_col <- cnv_first_col(meta, c(
    paste0("CNV_", method_key, "_cluster"),
    paste0("CNV_", method_key, "_clusters"),
    paste0(method_key, "_cluster"),
    paste0(method_key, "_clusters"),
    "cnv_clusters",
    "cnv_cluster",
    "CNV_cluster",
    "subclone",
    "subclones",
    "Subclone",
    "Clone"
  ))
  cell_info <- data.frame(cell = rownames(meta), stringsAsFactors = FALSE)
  if (!is.null(score_col)) {
    cell_info$score <- suppressWarnings(as.numeric(meta[[score_col]]))
  }
  if (!is.null(pred_col)) {
    cell_info$prediction <- as.character(meta[[pred_col]])
  }
  if (
    identical(method, "fastCNV") &&
      "score" %in% colnames(cell_info) &&
      (!"prediction" %in% colnames(cell_info) || all(is.na(cell_info$prediction) | !nzchar(cell_info$prediction)))
  ) {
    cell_info$prediction <- cnv_fastcnv_prediction_from_score(
      score = cell_info$score,
      meta = meta,
      reference.by = reference.by,
      reference = reference
    )
  }
  if (!is.null(cluster_col)) {
    cell_info$cluster <- as.character(meta[[cluster_col]])
  }
  rownames(cell_info) <- cell_info$cell
  candidates <- cnv_find_matrix_candidates(result@tools, cells = cells)
  if (is.null(candidates$matrix)) {
    candidates <- cnv_find_seurat_assay_matrix(result, cells = cells)
  }
  list(
    cnv_matrix = candidates$matrix,
    bin_info = if (!is.null(candidates$matrix)) cnv_bin_info_from_matrix(candidates$matrix) else NULL,
    cell_info = cell_info,
    raw = list(class = class(result), meta_cols = colnames(meta), tool_names = names(result@tools))
  )
}

cnv_combine_extracted_results <- function(parts, cells, method) {
  matrices <- lapply(parts, `[[`, "cnv_matrix")
  matrices <- matrices[!vapply(matrices, is.null, logical(1))]
  mat <- NULL
  bin_info <- NULL
  if (length(matrices) > 0L) {
    bins <- unique(unlist(lapply(matrices, rownames), use.names = FALSE))
    mat <- matrix(NA_real_, nrow = length(bins), ncol = length(cells), dimnames = list(bins, cells))
    for (x in matrices) {
      x <- cnv_orient_matrix(x, cells = cells)
      mat[rownames(x), colnames(x)] <- x
    }
    bin_info_parts <- lapply(parts, `[[`, "bin_info")
    bin_info_parts <- bin_info_parts[!vapply(bin_info_parts, is.null, logical(1))]
    if (length(bin_info_parts) > 0L) {
      bin_info <- do.call(rbind, bin_info_parts)
      bin_info <- bin_info[!duplicated(bin_info$bin_id), , drop = FALSE]
      rownames(bin_info) <- bin_info$bin_id
      bin_info <- bin_info[intersect(rownames(mat), rownames(bin_info)), , drop = FALSE]
    }
  }
  cell_info <- do.call(rbind, lapply(parts, `[[`, "cell_info"))
  if (!is.null(cell_info) && nrow(cell_info) > 0L) {
    cell_info <- cell_info[!duplicated(cell_info$cell), , drop = FALSE]
    rownames(cell_info) <- cell_info$cell
  }
  list(
    cnv_matrix = mat,
    bin_info = bin_info,
    cell_info = cell_info,
    raw = list(class = "list", method = method, n_results = length(parts))
  )
}

cnv_read_numbat_output <- function(output_dir = NULL) {
  if (is.null(output_dir) || !dir.exists(output_dir)) {
    return(list())
  }
  joint_file <- cnv_latest_file(output_dir, "^joint_post_.*\\.tsv(\\.gz)?$")
  clone_file <- cnv_latest_file(output_dir, "^clone_post_.*\\.tsv(\\.gz)?$")
  seg_file <- cnv_latest_file(output_dir, "^segs_consensus_.*\\.tsv(\\.gz)?$")
  out <- list()
  if (!is.null(joint_file)) {
    out$joint_post <- cnv_read_table(joint_file)
  }
  if (!is.null(clone_file)) {
    out$clone_post <- cnv_read_table(clone_file)
  }
  if (!is.null(seg_file)) {
    out$segs_consensus <- cnv_read_table(seg_file)
  }
  files <- c(joint_post = joint_file, clone_post = clone_file, segs_consensus = seg_file)
  out$files <- files[!vapply(files, is.null, logical(1))]
  out
}

cnv_latest_file <- function(path, pattern) {
  files <- list.files(path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) {
    return(NULL)
  }
  files[order(file.info(files)$mtime, decreasing = TRUE)][[1L]]
}

cnv_read_table <- function(path) {
  utils::read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

cnv_object_fields <- function(x) {
  fields <- c(
    "joint_post", "exp_post", "allele_post", "clone_post",
    "segs_consensus", "gexp_roll_wide", "P"
  )
  out <- list()
  for (field in fields) {
    value <- tryCatch(x[[field]], error = function(e) NULL)
    if (!is.null(value)) {
      out[[field]] <- value
    }
  }
  out
}

cnv_matrix_from_numbat_joint_post <- function(joint_post, cells) {
  df <- as.data.frame(joint_post, check.names = FALSE)
  cell_col <- cnv_first_col(df, c("cell", "barcode", "cell_id", "Cell"))
  seg_col <- cnv_first_col(df, c("seg_label", "seg", "seg_cons", "component", "CNV"))
  state_col <- cnv_first_col(df, c("cnv_state_post", "cnv_state", "state_post", "state", "cnv_states"))
  if (is.null(cell_col) || is.null(seg_col)) {
    return(list(matrix = NULL, bin_info = NULL))
  }
  df <- df[as.character(df[[cell_col]]) %in% cells, , drop = FALSE]
  if (nrow(df) == 0L) {
    return(list(matrix = NULL, bin_info = NULL))
  }
  score <- if (!is.null(state_col)) {
    cnv_state_to_numeric(df[[state_col]])
  } else {
    rep(NA_real_, nrow(df))
  }
  prob_col <- cnv_first_col(df, c("p_cnv", "p_cnv_x", "p_cnv_y", "posterior", "prob", "p"))
  if (!is.null(prob_col)) {
    prob <- suppressWarnings(as.numeric(df[[prob_col]]))
    score <- score * ifelse(is.finite(prob), prob, 1)
  }
  if (all(!is.finite(score))) {
    score <- rep(0, nrow(df))
  }
  segs <- unique(as.character(df[[seg_col]]))
  mat <- matrix(NA_real_, nrow = length(segs), ncol = length(cells), dimnames = list(segs, cells))
  keys <- paste(as.character(df[[seg_col]]), as.character(df[[cell_col]]), sep = "\r")
  values <- stats::aggregate(score, by = list(key = keys), FUN = mean, na.rm = TRUE)
  split_key <- strsplit(values$key, "\r", fixed = TRUE)
  for (i in seq_len(nrow(values))) {
    mat[split_key[[i]][[1L]], split_key[[i]][[2L]]] <- values$x[[i]]
  }
  bin_info <- cnv_numbat_bin_info(df, seg_col = seg_col, segs = segs)
  list(matrix = mat, bin_info = bin_info)
}

cnv_state_to_numeric <- function(x) {
  state <- tolower(as.character(x))
  out <- rep(0, length(state))
  out[grepl("del|loss|bdel", state)] <- -1
  out[grepl("amp|gain|bamp", state)] <- 1
  out[grepl("loh|neu|normal|diploid", state)] <- 0
  out[is.na(state) | !nzchar(state)] <- NA_real_
  out
}

cnv_numbat_bin_info <- function(df, seg_col, segs) {
  chr_col <- cnv_first_col(df, c("CHROM", "chrom", "chr", "chromosome"))
  start_col <- cnv_first_col(df, c("seg_start", "start", "Start"))
  end_col <- cnv_first_col(df, c("seg_end", "end", "End"))
  out <- data.frame(
    bin_id = segs,
    chr = NA_character_,
    start = NA_real_,
    end = NA_real_,
    gene = segs,
    stringsAsFactors = FALSE
  )
  for (seg in segs) {
    rows <- which(as.character(df[[seg_col]]) == seg)
    if (length(rows) == 0L) {
      next
    }
    i <- match(seg, out$bin_id)
    if (!is.null(chr_col)) {
      out$chr[[i]] <- as.character(df[[chr_col]][[rows[[1L]]]])
    }
    if (!is.null(start_col)) {
      out$start[[i]] <- suppressWarnings(min(as.numeric(df[[start_col]][rows]), na.rm = TRUE))
    }
    if (!is.null(end_col)) {
      out$end[[i]] <- suppressWarnings(max(as.numeric(df[[end_col]][rows]), na.rm = TRUE))
    }
  }
  out$start[!is.finite(out$start)] <- NA_real_
  out$end[!is.finite(out$end)] <- NA_real_
  rownames(out) <- out$bin_id
  out
}

cnv_cell_info_from_numbat_clone_post <- function(clone_post, cells) {
  if (is.null(clone_post)) {
    return(NULL)
  }
  df <- as.data.frame(clone_post, check.names = FALSE)
  cell_col <- cnv_first_col(df, c("cell", "barcode", "cell_id", "Cell"))
  if (is.null(cell_col) || !any(as.character(df[[cell_col]]) %in% cells)) {
    return(NULL)
  }
  score_col <- cnv_first_col(df, c("p_cnv", "p_cnv_x", "p_cnv_y", "prob", "posterior"))
  pred_col <- cnv_first_col(df, c("compartment_opt", "compartment", "prediction", "classification"))
  cluster_col <- cnv_first_col(df, c("clone_opt", "clone", "GT_opt", "cluster", "subclone"))
  out <- data.frame(cell = as.character(df[[cell_col]]), stringsAsFactors = FALSE)
  if (!is.null(score_col)) {
    out$score <- suppressWarnings(as.numeric(df[[score_col]]))
  }
  if (!is.null(pred_col)) {
    out$prediction <- as.character(df[[pred_col]])
  }
  if (!is.null(cluster_col)) {
    out$cluster <- as.character(df[[cluster_col]])
  }
  rownames(out) <- out$cell
  out
}

cnv_casper_matrix_from_objects <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  objects <- if (is.list(x) && !isS4(x)) x else list(x)
  for (object in objects) {
    if (isS4(object) && "large.scale.cnv.events" %in% methods::slotNames(object)) {
      events <- methods::slot(object, "large.scale.cnv.events")
      mat <- cnv_casper_events_to_matrix(events)
      if (!is.null(mat)) {
        return(mat)
      }
    }
  }
  NULL
}

cnv_casper_events_to_matrix <- function(events) {
  if (is.null(events) || !is.data.frame(events) || nrow(events) == 0L) {
    return(NULL)
  }
  arms <- as.vector(rbind(paste0(seq_len(22), "p"), paste0(seq_len(22), "q")))
  mat <- matrix(0, nrow = nrow(events), ncol = length(arms), dimnames = list(rownames(events), arms))
  amp_col <- cnv_first_col(events, c("LargeScaleAmp", "largeScaleAmp", "amp"))
  del_col <- cnv_first_col(events, c("LargeScaleDel", "largeScaleDel", "del"))
  for (i in seq_len(nrow(events))) {
    if (!is.null(amp_col)) {
      amp <- unlist(strsplit(as.character(events[[amp_col]][[i]]), "\\s+"))
      mat[i, intersect(amp, arms)] <- 1
    }
    if (!is.null(del_col)) {
      del <- unlist(strsplit(as.character(events[[del_col]][[i]]), "\\s+"))
      mat[i, intersect(del, arms)] <- -1
    }
  }
  mat
}

cnv_casper_bin_info <- function(mat) {
  data.frame(
    bin_id = rownames(mat),
    chr = sub("([0-9XYM]+)[pq]$", "chr\\1", rownames(mat), ignore.case = TRUE),
    start = NA_real_,
    end = NA_real_,
    gene = rownames(mat),
    stringsAsFactors = FALSE,
    row.names = rownames(mat)
  )
}

cnv_standardize_result <- function(
  result,
  method,
  cells,
  features,
  reference_cells = NULL,
  store_matrix = TRUE,
  parameters = list()
) {
  mat <- result$cnv_matrix
  if (is.null(mat)) {
    mat <- matrix(
      nrow = 0L,
      ncol = length(cells),
      dimnames = list(character(), cells)
    )
  }
  mat <- cnv_orient_matrix(mat, cells = cells)
  cell_info <- cnv_standardize_cell_info(
    result$cell_info,
    mat = mat,
    cells = cells,
    method = method,
    reference_cells = reference_cells
  )
  bin_info <- cnv_standardize_bin_info(result$bin_info, mat = mat)
  if (!isTRUE(store_matrix)) {
    mat_store <- NULL
  } else {
    mat_store <- mat
  }
  list(
    method = method,
    cnv_matrix = mat_store,
    matrix_dim = dim(mat),
    bin_info = bin_info,
    cell_info = cell_info,
    parameters = parameters,
    raw = result$raw %||% list()
  )
}

cnv_add_metadata <- function(srt, bundle, method, prefix, tool_name) {
  cell_info <- bundle$cell_info
  rownames(cell_info) <- cell_info$cell
  cell_info <- cell_info[colnames(srt), , drop = FALSE]
  score <- cell_info$score
  prediction <- cell_info$prediction
  cluster <- cell_info$cluster
  meta <- data.frame(row.names = colnames(srt))
  meta[[paste0(prefix, "_score")]] <- score
  meta[[paste0(prefix, "_prediction")]] <- prediction
  meta[[paste0(prefix, "_cluster")]] <- cluster
  method_prefix <- paste0(prefix, "_", method)
  meta[[paste0(method_prefix, "_score")]] <- score
  meta[[paste0(method_prefix, "_prediction")]] <- prediction
  meta[[paste0(method_prefix, "_cluster")]] <- cluster
  Seurat::AddMetaData(srt, metadata = meta)
}

cnv_match_method <- function(method) {
  if (is.null(method) || length(method) == 0L) {
    method <- "copykat"
  }
  method <- method[[1L]]
  cnv_assert_scalar_string(method, "method")
  method_map <- c(
    copykat = "copykat",
    fastcnv = "fastCNV",
    scevan = "scevan",
    infercnv = "infercnv",
    numbat = "numbat",
    casper = "casper"
  )
  method_key <- tolower(method)
  if (!method_key %in% names(method_map)) {
    log_message(
      "{.arg method} must be one of {.val copykat}, {.val fastCNV}, {.val scevan}, {.val infercnv}, {.val numbat}, or {.val casper}",
      message_type = "error"
    )
  }
  unname(method_map[[method_key]])
}

cnv_assert_seurat <- function(srt) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

cnv_assert_scalar_string <- function(x, arg) {
  if (is.null(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    log_message(
      "{.arg {arg}} must be a single non-empty string",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

cnv_validate_metadata <- function(srt, column, arg) {
  if (!is.null(column) && !column %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg {arg}} {.val {column}} is not present in {.arg srt@meta.data}",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

cnv_get_counts <- function(srt, assay, layer) {
  mat <- tryCatch(
    GetAssayData5(srt, assay = assay, layer = layer),
    error = function(e) e
  )
  if (inherits(mat, "error")) {
    log_message(
      "Unable to read {.arg layer} {.val {layer}} from assay {.val {assay}}: {.val {conditionMessage(mat)}}",
      message_type = "error"
    )
  }
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    log_message(
      "CNV input matrix must contain feature and cell names",
      message_type = "error"
    )
  }
  if (!inherits(mat, "Matrix")) {
    mat <- Matrix::Matrix(
      if (is.data.frame(mat)) as.matrix(mat) else mat,
      sparse = TRUE
    )
  }
  mat <- methods::as(mat, "dgCMatrix")
  mat@x[!is.finite(mat@x) | mat@x < 0] <- 0
  Matrix::drop0(mat)
}

cnv_reference_cells <- function(srt, method, reference.by, reference) {
  if (method %in% c("infercnv", "fastCNV", "casper") && (is.null(reference.by) || is.null(reference))) {
    log_message(
      "{.arg reference.by} and {.arg reference} are required for {.arg method = {method}}",
      message_type = "error"
    )
  }
  if (is.null(reference.by) || is.null(reference)) {
    return(NULL)
  }
  ref_values <- as.character(srt@meta.data[[reference.by]])
  cells <- colnames(srt)[ref_values %in% as.character(reference)]
  if (length(cells) == 0L) {
    log_message(
      "No reference cells matched {.arg reference.by} and {.arg reference}",
      message_type = "error"
    )
  }
  attr(reference, "cells") <- cells
  cells
}

cnv_validate_backend_inputs <- function(
  method,
  allele_counts = NULL,
  reference_counts = NULL,
  loh = NULL,
  cytoband = NULL,
  gene_order = NULL
) {
  if (identical(method, "numbat")) {
    if (is.null(allele_counts)) {
      log_message(
        "{.arg allele_counts} is required for {.arg method = 'numbat'}",
        message_type = "error"
      )
    }
    if (is.null(reference_counts)) {
      log_message(
        "{.arg reference_counts} is required for {.arg method = 'numbat'}",
        message_type = "error"
      )
    }
  }
  if (identical(method, "casper")) {
    if (is.null(gene_order) || nrow(gene_order) == 0L) {
      log_message(
        "{.arg gene_order} is required for {.arg method = 'casper'}",
        message_type = "error"
      )
    }
    if (is.null(loh)) {
      log_message(
        "{.arg loh} is required for {.arg method = 'casper'}",
        message_type = "error"
      )
    }
    if (is.null(cytoband)) {
      log_message(
        "{.arg cytoband} is required for {.arg method = 'casper'}",
        message_type = "error"
      )
    }
  }
  invisible(TRUE)
}

cnv_fastcnv_sample_input <- function(srt, sample.by = NULL) {
  if (is.null(sample.by)) {
    return(list(seuratObj = srt, sampleName = "scop_cnv"))
  }
  samples <- as.character(srt@meta.data[[sample.by]])
  sample_levels <- unique(samples[!is.na(samples) & nzchar(samples)])
  if (length(sample_levels) <= 1L) {
    sample_name <- if (length(sample_levels) == 1L) sample_levels[[1L]] else "scop_cnv"
    return(list(seuratObj = srt, sampleName = sample_name))
  }
  split_obj <- Seurat::SplitObject(srt, split.by = sample.by)
  split_obj <- split_obj[names(split_obj) %in% sample_levels]
  list(
    seuratObj = split_obj,
    sampleName = names(split_obj)
  )
}

cnv_fastcnv_prediction_from_score <- function(score, meta, reference.by = NULL, reference = NULL) {
  reference_cells <- NULL
  if (!is.null(reference.by) && reference.by %in% colnames(meta) && !is.null(reference)) {
    reference_cells <- rownames(meta)[as.character(meta[[reference.by]]) %in% as.character(reference)]
  }
  cnv_prediction_from_score(
    score = score,
    cells = rownames(meta),
    reference_cells = reference_cells
  )
}

cnv_resolve_gene_order <- function(
  gene_order,
  srt,
  assay,
  method,
  required = FALSE
) {
  if (!is.null(gene_order)) {
    if (is.character(gene_order) && length(gene_order) == 1L) {
      path <- normalizePath(path.expand(gene_order), mustWork = TRUE, winslash = "/")
      gene_order <- if (grepl("\\.csv$", path, ignore.case = TRUE)) {
        utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
      } else {
        utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
      }
    }
    return(cnv_normalize_gene_order(gene_order, required = TRUE))
  }

  feature_meta <- tryCatch(GetFeaturesData(srt, assay = assay), error = function(e) NULL)
  if (!is.null(feature_meta) && ncol(feature_meta) > 0L) {
    feature_meta <- as.data.frame(feature_meta, check.names = FALSE)
    feature_meta$feature <- rownames(feature_meta)
    out <- tryCatch(
      cnv_normalize_gene_order(feature_meta, required = FALSE),
      error = function(e) NULL
    )
    if (!is.null(out) && nrow(out) > 0L) {
      return(out)
    }
  }

  if (isTRUE(required)) {
    log_message(
      paste(
        "{.arg gene_order} is required for {.arg method = {method}}.",
        "Provide a table/path with gene, chromosome, start, and end columns,",
        "or add those columns to assay feature metadata."
      ),
      message_type = "error"
    )
  }
  NULL
}

cnv_normalize_gene_order <- function(gene_order, required = TRUE) {
  if (!is.data.frame(gene_order)) {
    log_message(
      "{.arg gene_order} must be a data frame or readable file path",
      message_type = "error"
    )
  }
  nm <- colnames(gene_order)
  gene_col <- cnv_first_name(nm, c("gene", "Gene", "symbol", "features", "feature", "gene_name", "GeneSymbol"))
  chr_col <- cnv_first_name(nm, c("chr", "chrom", "chromosome", "Chromosome", "seqnames"))
  start_col <- cnv_first_name(nm, c("start", "Start", "txStart", "begin"))
  end_col <- cnv_first_name(nm, c("end", "End", "txEnd", "stop"))
  if (any(vapply(list(gene_col, chr_col, start_col, end_col), is.null, logical(1)))) {
    if (isTRUE(required)) {
      log_message(
        "{.arg gene_order} must contain gene, chromosome, start, and end columns",
        message_type = "error"
      )
    }
    return(NULL)
  }
  out <- data.frame(
    gene = as.character(gene_order[[gene_col]]),
    chr = as.character(gene_order[[chr_col]]),
    start = suppressWarnings(as.numeric(gene_order[[start_col]])),
    end = suppressWarnings(as.numeric(gene_order[[end_col]])),
    stringsAsFactors = FALSE
  )
  keep <- nzchar(out$gene) & nzchar(out$chr) & is.finite(out$start) & is.finite(out$end)
  out <- out[keep, , drop = FALSE]
  rownames(out) <- NULL
  if (nrow(out) == 0L && isTRUE(required)) {
    log_message(
      "{.arg gene_order} contains no valid rows after parsing",
      message_type = "error"
    )
  }
  out
}

cnv_gene_order_to_casper_annotation <- function(gene_order) {
  annotation <- data.frame(
    Gene = gene_order$gene,
    Chr = sub("^chr", "", gene_order$chr, ignore.case = TRUE),
    Start = gene_order$start,
    End = gene_order$end,
    Position = rowMeans(cbind(gene_order$start, gene_order$end), na.rm = TRUE),
    cytoband = NA_character_,
    stringsAsFactors = FALSE
  )
  rownames(annotation) <- annotation$Gene
  annotation
}

cnv_orient_matrix <- function(mat, cells) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  if (is.null(rownames(mat))) {
    rownames(mat) <- paste0("bin_", seq_len(nrow(mat)))
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- cells[seq_len(min(ncol(mat), length(cells)))]
  }
  row_match <- sum(rownames(mat) %in% cells)
  col_match <- sum(colnames(mat) %in% cells)
  if (row_match > col_match) {
    mat <- t(mat)
  }
  out <- matrix(
    NA_real_,
    nrow = nrow(mat),
    ncol = length(cells),
    dimnames = list(rownames(mat), cells)
  )
  common <- intersect(colnames(mat), cells)
  if (length(common) > 0L) {
    out[, common] <- mat[, common, drop = FALSE]
  } else if (ncol(mat) == length(cells)) {
    colnames(mat) <- cells
    out <- mat
  } else {
    log_message(
      "CNV matrix columns could not be aligned to cells",
      message_type = "error"
    )
  }
  out
}

cnv_standardize_cell_info <- function(
  cell_info,
  mat,
  cells,
  method = NULL,
  reference_cells = NULL
) {
  if (is.null(cell_info)) {
    cell_info <- data.frame(cell = cells, stringsAsFactors = FALSE)
  } else {
    cell_info <- as.data.frame(cell_info, check.names = FALSE)
    if (!"cell" %in% colnames(cell_info)) {
      cell_info$cell <- rownames(cell_info) %||% cells[seq_len(nrow(cell_info))]
    }
    cell_info$cell <- as.character(cell_info$cell)
  }
  rownames(cell_info) <- cell_info$cell
  out <- data.frame(cell = cells, stringsAsFactors = FALSE)
  rownames(out) <- out$cell
  common <- intersect(cells, rownames(cell_info))
  if (length(common) > 0L) {
    for (nm in setdiff(colnames(cell_info), "cell")) {
      out[common, nm] <- cell_info[common, nm]
    }
  }
  if (!"score" %in% colnames(out)) {
    out$score <- if (nrow(mat) == 0L) {
      rep(NA_real_, length(cells))
    } else {
      colMeans(abs(mat), na.rm = TRUE)
    }
  }
  out$score <- suppressWarnings(as.numeric(out$score))
  out$score[!is.finite(out$score)] <- NA_real_
  if (!"prediction" %in% colnames(out)) {
    out$prediction <- NA_character_
  }
  if (!"cluster" %in% colnames(out)) {
    out$cluster <- NA_character_
  }
  out$prediction <- as.character(out$prediction)
  missing_prediction <- is.na(out$prediction) | !nzchar(out$prediction)
  if (all(missing_prediction) && any(is.finite(out$score))) {
    out$prediction <- cnv_prediction_from_score(
      score = out$score,
      cells = out$cell,
      reference_cells = reference_cells
    )
  }
  out$cluster <- as.character(out$cluster)
  out[, c("cell", "score", "prediction", "cluster"), drop = FALSE]
}

cnv_prediction_from_score <- function(score, cells = NULL, reference_cells = NULL) {
  score <- suppressWarnings(as.numeric(score))
  out <- rep(NA_character_, length(score))
  finite <- is.finite(score)
  if (!any(finite)) {
    return(out)
  }

  reference_idx <- rep(FALSE, length(score))
  if (!is.null(cells) && !is.null(reference_cells) && length(reference_cells) > 0L) {
    reference_idx <- as.character(cells) %in% as.character(reference_cells)
  }

  cutoff <- 0
  reference_score <- score[reference_idx & finite]
  if (length(reference_score) > 0L) {
    cutoff <- stats::quantile(reference_score, probs = 0.99, na.rm = TRUE, names = FALSE)
    if (!is.finite(cutoff)) {
      cutoff <- 0
    }
    cutoff <- max(cutoff, 0)
  }

  out[finite] <- ifelse(score[finite] > cutoff, "aneuploid", "diploid")
  out[reference_idx & finite] <- "diploid"
  out
}

cnv_standardize_bin_info <- function(bin_info, mat) {
  if (is.null(bin_info)) {
    bin_info <- data.frame(bin_id = rownames(mat), stringsAsFactors = FALSE)
  } else {
    bin_info <- as.data.frame(bin_info, check.names = FALSE)
  }
  if (!"bin_id" %in% colnames(bin_info)) {
    bin_info$bin_id <- rownames(mat)[seq_len(nrow(bin_info))] %||% paste0("bin_", seq_len(nrow(bin_info)))
  }
  for (nm in c("chr", "start", "end", "gene")) {
    if (!nm %in% colnames(bin_info)) {
      bin_info[[nm]] <- if (nm %in% c("start", "end")) NA_real_ else NA_character_
    }
  }
  if (nrow(bin_info) != nrow(mat)) {
    bin_info <- data.frame(
      bin_id = rownames(mat),
      chr = NA_character_,
      start = NA_real_,
      end = NA_real_,
      gene = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  rownames(bin_info) <- rownames(mat)
  bin_info[, c("bin_id", "chr", "start", "end", "gene"), drop = FALSE]
}

cnv_bin_info_from_matrix <- function(mat, coord_cols = NULL) {
  df <- as.data.frame(mat, check.names = FALSE)
  coord_cols <- coord_cols %||% intersect(
    c("bin_id", "chrom", "chr", "chromosome", "Chromosome", "start", "Start", "end", "End", "abspos", "chrompos", "cytoband", "gene", "genes"),
    colnames(df)
  )
  if (length(coord_cols) == 0L) {
    return(data.frame(
      bin_id = rownames(df) %||% paste0("bin_", seq_len(nrow(df))),
      chr = NA_character_,
      start = NA_real_,
      end = NA_real_,
      gene = rownames(df) %||% NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  chr_col <- cnv_first_col(df, c("chr", "chrom", "chromosome", "Chromosome"))
  start_col <- cnv_first_col(df, c("start", "Start", "abspos", "chrompos"))
  end_col <- cnv_first_col(df, c("end", "End"))
  gene_col <- cnv_first_col(df, c("gene", "genes", "Gene"))
  bin_id <- rownames(df) %||% paste0("bin_", seq_len(nrow(df)))
  data.frame(
    bin_id = bin_id,
    chr = if (!is.null(chr_col)) as.character(df[[chr_col]]) else NA_character_,
    start = if (!is.null(start_col)) suppressWarnings(as.numeric(df[[start_col]])) else NA_real_,
    end = if (!is.null(end_col)) suppressWarnings(as.numeric(df[[end_col]])) else NA_real_,
    gene = if (!is.null(gene_col)) as.character(df[[gene_col]]) else bin_id,
    stringsAsFactors = FALSE
  )
}

cnv_drop_coordinate_columns <- function(mat) {
  df <- as.data.frame(mat, check.names = FALSE)
  coord_cols <- intersect(
    c("bin_id", "chrom", "chr", "chromosome", "Chromosome", "start", "Start", "end", "End", "abspos", "chrompos", "cytoband", "gene", "genes"),
    colnames(df)
  )
  df[, setdiff(colnames(df), coord_cols), drop = FALSE]
}

cnv_find_matrix_candidates <- function(x, cells) {
  if (is.data.frame(x)) {
    cell_info <- cnv_parse_cell_info_df(x, cells = cells)
    if (!is.null(cell_info) && ncol(cell_info) > 1L && !any(colnames(x) %in% cells) && ncol(x) != length(cells)) {
      return(list(matrix = NULL))
    }
  }
  if (inherits(x, c("matrix", "data.frame", "Matrix"))) {
    mat <- as.matrix(x)
    if (cnv_matrix_has_cells(mat, cells)) {
      return(list(matrix = x))
    }
  }
  if (isS4(x) && !inherits(x, "Seurat")) {
    preferred_slots <- intersect(
      c("expr.data", "count.data", "cnv_matrix", "CNAmat", "CNA", "CNA_profile"),
      methods::slotNames(x)
    )
    for (nm in preferred_slots) {
      candidate <- methods::slot(x, nm)
      if (inherits(candidate, c("matrix", "data.frame", "Matrix")) && cnv_matrix_has_cells(candidate, cells)) {
        return(list(matrix = candidate))
      }
    }
  }
  if (!is.list(x)) {
    return(list(matrix = NULL))
  }
  preferred <- c(
    "cnv_matrix", "CNAmat", "CNA", "CNA_profile", "cnv", "cnv_mat",
    "CNV", "expr.data", "cnv_profile", "cna_mat"
  )
  for (nm in intersect(preferred, names(x))) {
    candidate <- x[[nm]]
    if (inherits(candidate, c("matrix", "data.frame", "Matrix")) && cnv_matrix_has_cells(candidate, cells)) {
      return(list(matrix = candidate))
    }
  }
  for (nm in names(x)) {
    candidate <- x[[nm]]
    if (inherits(candidate, c("matrix", "data.frame", "Matrix")) && cnv_matrix_has_cells(candidate, cells)) {
      return(list(matrix = candidate))
    }
  }
  for (nm in names(x)) {
    if (is.list(x[[nm]])) {
      found <- cnv_find_matrix_candidates(x[[nm]], cells = cells)
      if (!is.null(found$matrix)) {
        return(found)
      }
    }
  }
  list(matrix = NULL)
}

cnv_matrix_has_cells <- function(mat, cells) {
  rn <- rownames(mat) %||% character(0)
  cn <- colnames(mat) %||% character(0)
  sum(rn %in% cells) > 0L || sum(cn %in% cells) > 0L || ncol(mat) == length(cells)
}

cnv_find_cell_info <- function(x, cells) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.data.frame(x)) {
    return(cnv_parse_cell_info_df(x, cells))
  }
  if (!is.list(x)) {
    return(NULL)
  }
  preferred <- c("cell_info", "prediction", "predictions", "classification", "classifications", "cell_annotation")
  for (nm in intersect(preferred, names(x))) {
    info <- cnv_find_cell_info(x[[nm]], cells = cells)
    if (!is.null(info)) {
      return(info)
    }
  }
  for (nm in names(x)) {
    info <- cnv_find_cell_info(x[[nm]], cells = cells)
    if (!is.null(info)) {
      return(info)
    }
  }
  NULL
}

cnv_parse_cell_info_df <- function(df, cells) {
  df <- as.data.frame(df, check.names = FALSE)
  cell_col <- cnv_first_col(df, c("cell", "cells", "cell.names", "barcode", "Cell", "cell_id"))
  if (is.null(cell_col) && !is.null(rownames(df)) && any(rownames(df) %in% cells)) {
    df$cell <- rownames(df)
    cell_col <- "cell"
  }
  if (is.null(cell_col) || !any(as.character(df[[cell_col]]) %in% cells)) {
    return(NULL)
  }
  out <- data.frame(cell = as.character(df[[cell_col]]), stringsAsFactors = FALSE)
  score_col <- cnv_first_col(df, c("score", "cnv_score", "CNV_score", "cnv_fraction", "CNV_fraction"))
  pred_col <- cnv_first_col(df, c("prediction", "pred", "copykat.pred", "classification", "class", "malignant", "status", "cell.assignment", "cell_assignment"))
  cluster_col <- cnv_first_col(df, c("cluster", "cnv_cluster", "CNV_cluster", "cnv_clusters", "subclone", "subclones", "Subclone", "Clone"))
  if (!is.null(score_col)) {
    out$score <- suppressWarnings(as.numeric(df[[score_col]]))
  }
  if (!is.null(pred_col)) {
    out$prediction <- as.character(df[[pred_col]])
  }
  if (!is.null(cluster_col)) {
    out$cluster <- as.character(df[[cluster_col]])
  }
  rownames(out) <- out$cell
  out
}

cnv_light_raw <- function(x) {
  if (inherits(x, "Seurat")) {
    return(list(class = class(x), meta_cols = colnames(x@meta.data), tool_names = names(x@tools)))
  }
  if (isS4(x)) {
    return(list(class = class(x), slots = methods::slotNames(x)))
  }
  if (is.list(x)) {
    return(list(class = class(x), names = names(x)))
  }
  list(class = class(x))
}

cnv_first_col <- function(df, candidates) {
  cnv_first_name(colnames(df), candidates)
}

cnv_first_name <- function(names, candidates) {
  hit <- candidates[candidates %in% names]
  if (length(hit) == 0L) {
    return(NULL)
  }
  hit[[1L]]
}

cnv_get_first_namespace_fun <- function(package, names) {
  for (nm in names) {
    fun <- tryCatch(get_namespace_fun(package, nm), error = function(e) NULL)
    if (is.function(fun)) {
      return(fun)
    }
  }
  log_message(
    "{.pkg {package}} does not expose one of {.val {names}}",
    message_type = "error"
  )
}

cnv_call_backend_fun <- function(fun, args) {
  formal_names <- names(formals(fun))
  if (is.null(formal_names)) {
    return(do.call(fun, args))
  }
  allowed <- setdiff(formal_names, "...")
  do.call(fun, args[names(args) %in% allowed])
}

cnv_default_threads <- function(max_threads = 4L) {
  cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (length(cores) == 0L || !is.finite(cores) || cores < 1L) {
    return(1L)
  }
  max(1L, min(as.integer(max_threads), as.integer(cores)))
}

cnv_find_seurat_assay_matrix <- function(result, cells) {
  assay_names <- intersect(
    c("genomicScores", "rawGenomicScores", "CNV", "CNA", "cnv", "cna"),
    names(result@assays)
  )
  for (assay in assay_names) {
    for (layer in c("data", "counts", "scale.data")) {
      mat <- tryCatch(
        GetAssayData5(result, assay = assay, layer = layer),
        error = function(e) NULL
      )
      if (!is.null(mat) && cnv_matrix_has_cells(mat, cells)) {
        return(list(matrix = mat))
      }
    }
  }
  list(matrix = NULL)
}

cnv_read_scevan_output <- function(output_dir, sample, cells) {
  dirs <- unique(normalizePath(
    c(output_dir, file.path(output_dir, "output")),
    mustWork = FALSE,
    winslash = "/"
  ))
  dirs <- dirs[dir.exists(dirs)]
  matrix_files <- unique(unlist(lapply(
    dirs,
    list.files,
    pattern = "CNAmtx.*\\.RData$",
    full.names = TRUE,
    recursive = TRUE
  ), use.names = FALSE))
  if (length(matrix_files) == 0L) {
    return(list(matrix = NULL, bin_info = NULL, files = character()))
  }
  preferred <- grep(paste0("^", sample, ".*CNAmtx.*\\.RData$"), basename(matrix_files), value = FALSE)
  if (length(preferred) > 0L) {
    matrix_files <- matrix_files[c(preferred, setdiff(seq_along(matrix_files), preferred))]
  }
  matrix_file <- matrix_files[[1L]]
  mat <- cnv_load_rdata_matrix(matrix_file, cells = cells)
  annot_file <- sub("CNAmtx", "count_mtx_annot", matrix_file, fixed = TRUE)
  bin_info <- NULL
  files <- matrix_file
  if (file.exists(annot_file)) {
    annot <- cnv_load_rdata_data_frame(annot_file)
    bin_info <- cnv_scevan_bin_info(annot, mat = mat)
    files <- c(files, annot_file)
  }
  list(matrix = mat, bin_info = bin_info, files = files)
}

cnv_load_rdata_matrix <- function(path, cells) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  objects <- mget(ls(env), envir = env)
  preferred <- intersect(c("CNA_mtx_relat", "CNA_mtx", "CNAmtx", "cnv_matrix"), names(objects))
  for (nm in c(preferred, setdiff(names(objects), preferred))) {
    x <- objects[[nm]]
    if (inherits(x, c("matrix", "data.frame", "Matrix")) && cnv_matrix_has_cells(as.matrix(x), cells)) {
      return(x)
    }
  }
  NULL
}

cnv_load_rdata_data_frame <- function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  objects <- mget(ls(env), envir = env)
  for (x in objects) {
    if (is.data.frame(x)) {
      return(x)
    }
  }
  NULL
}

cnv_scevan_bin_info <- function(annot, mat) {
  if (is.null(annot) || is.null(mat)) {
    return(NULL)
  }
  annot <- as.data.frame(annot, check.names = FALSE)
  chr_col <- cnv_first_col(annot, c("seqnames", "chr", "chrom", "chromosome", "Chromosome"))
  start_col <- cnv_first_col(annot, c("start", "Start", "txStart", "begin"))
  end_col <- cnv_first_col(annot, c("end", "End", "txEnd", "stop"))
  gene_col <- cnv_first_col(annot, c("gene", "genes", "Gene", "symbol", "gene_name"))
  data.frame(
    bin_id = rownames(mat) %||% paste0("bin_", seq_len(nrow(mat))),
    chr = if (!is.null(chr_col)) as.character(annot[[chr_col]]) else NA_character_,
    start = if (!is.null(start_col)) suppressWarnings(as.numeric(annot[[start_col]])) else NA_real_,
    end = if (!is.null(end_col)) suppressWarnings(as.numeric(annot[[end_col]])) else NA_real_,
    gene = if (!is.null(gene_col)) as.character(annot[[gene_col]]) else rownames(mat) %||% NA_character_,
    stringsAsFactors = FALSE
  )
}
