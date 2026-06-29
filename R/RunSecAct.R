#' @title Run SecAct secreted protein activity inference
#'
#' @description
#' Run the optional `SecAct` R package from `scop` to infer secreted protein
#' signaling activity from bulk, single-cell, or spatial transcriptomics
#' profiles. `SecAct` is not bundled with `scop`; it is checked and installed
#' from `data2intelligence/SecAct` when this function is called.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt Optional Seurat object. When `mode = "scRNAseq"`, this is passed
#' to `SecAct.activity.inference.scRNAseq`. When `mode = "matrix"`, expression
#' is extracted from `srt` if `inputProfile` is not supplied.
#' @param inputProfile Expression matrix, Seurat object, or SpaCET object.
#' Matrix input must be genes x samples, cells, or spatial spots.
#' @param inputProfile_control Optional control expression matrix or SpaCET
#' object passed to SecAct.
#' @param mode SecAct workflow. `auto` dispatches from the input class.
#' @param assay Assay used when extracting a matrix from `srt`; also sets the
#' default assay for Seurat input passed to SecAct.
#' @param layer Assay layer used when `mode = "matrix"` and `inputProfile` is
#' not supplied.
#' @param cellType_meta Metadata column containing cell type or state labels for
#' `mode = "scRNAseq"` when `is.singleCellLevel = FALSE`.
#' @param is.singleCellLevel Whether SecAct should return single-cell activity
#' rather than cell-state activity for Seurat input.
#' @param is.differential,is.paired,is.singleSampleLevel Parameters passed to
#' `SecAct.activity.inference` for matrix or bulk input.
#' @param scale.factor Spot-level scale factor passed to
#' `SecAct.activity.inference.ST`.
#' @param sigMatrix SecAct signature matrix name.
#' @param is.filter.sig,is.group.sig,is.group.cor,lambda,nrand,ncores,backend,rng_method
#' Parameters passed to SecAct activity inference.
#' @param batch_size,output_h5 Optional large-matrix controls passed only to
#' `SecAct.activity.inference`.
#' @param activity Activity matrix to store as a Seurat assay when possible.
#' @param assay_out Name of the Seurat assay used to store activity values.
#' @param store_assay Whether to store a SecAct activity matrix as a Seurat
#' assay when its columns match Seurat cells.
#' @param store_results Whether to store raw SecAct results in `srt@tools`.
#' @param tool_name Name used in `srt@tools`.
#' @param return_seurat Whether to return a Seurat object when `srt` is
#' supplied.
#'
#' @return A Seurat object, SpaCET object, or SecAct result list depending on
#' input and `return_seurat`.
#' @export
#'
#' @examples
#' \dontrun{
#' srt <- RunSecAct(
#'   srt,
#'   mode = "scRNAseq",
#'   cellType_meta = "celltype",
#'   is.singleCellLevel = TRUE
#' )
#'
#' res <- RunSecAct(
#'   inputProfile = expr_mat,
#'   mode = "matrix",
#'   is.differential = TRUE
#' )
#' }
RunSecAct <- function(
  srt = NULL,
  inputProfile = NULL,
  inputProfile_control = NULL,
  mode = c("auto", "matrix", "scRNAseq", "ST"),
  assay = NULL,
  layer = "counts",
  cellType_meta = NULL,
  is.singleCellLevel = FALSE,
  is.differential = FALSE,
  is.paired = FALSE,
  is.singleSampleLevel = FALSE,
  scale.factor = NULL,
  sigMatrix = "SecAct",
  is.filter.sig = FALSE,
  is.group.sig = TRUE,
  is.group.cor = 0.9,
  lambda = 5e5,
  nrand = 1000,
  ncores = 1L,
  backend = "auto",
  rng_method = "mt19937",
  batch_size = NULL,
  output_h5 = NULL,
  activity = c("zscore", "beta", "se", "pvalue"),
  assay_out = "SecAct",
  store_assay = !is.null(srt),
  store_results = TRUE,
  tool_name = "SecAct",
  return_seurat = !is.null(srt),
  verbose = TRUE
) {
  mode <- match.arg(mode)
  activity <- match.arg(activity)
  secact_check_r(verbose = verbose)

  if (identical(mode, "auto")) {
    mode <- secact_resolve_mode(srt = srt, inputProfile = inputProfile)
  }

  if (identical(mode, "scRNAseq")) {
    return(secact_run_scrnaseq(
      srt = srt,
      inputProfile = inputProfile,
      assay = assay,
      cellType_meta = cellType_meta,
      is.singleCellLevel = is.singleCellLevel,
      sigMatrix = sigMatrix,
      is.filter.sig = is.filter.sig,
      is.group.sig = is.group.sig,
      is.group.cor = is.group.cor,
      lambda = lambda,
      nrand = nrand,
      ncores = ncores,
      backend = backend,
      rng_method = rng_method,
      activity = activity,
      assay_out = assay_out,
      store_assay = store_assay,
      store_results = store_results,
      tool_name = tool_name,
      return_seurat = return_seurat,
      verbose = verbose
    ))
  }

  if (identical(mode, "ST")) {
    return(secact_run_st(
      inputProfile = inputProfile %||% srt,
      inputProfile_control = inputProfile_control,
      scale.factor = scale.factor %||% 1e5,
      sigMatrix = sigMatrix,
      is.filter.sig = is.filter.sig,
      is.group.sig = is.group.sig,
      is.group.cor = is.group.cor,
      lambda = lambda,
      nrand = nrand,
      ncores = ncores,
      backend = backend,
      rng_method = rng_method,
      verbose = verbose
    ))
  }

  secact_run_matrix(
    srt = srt,
    inputProfile = inputProfile,
    inputProfile_control = inputProfile_control,
    assay = assay,
    layer = layer,
    is.differential = is.differential,
    is.paired = is.paired,
    is.singleSampleLevel = is.singleSampleLevel,
    sigMatrix = sigMatrix,
    is.filter.sig = is.filter.sig,
    is.group.sig = is.group.sig,
    is.group.cor = is.group.cor,
    lambda = lambda,
    nrand = nrand,
    ncores = ncores,
    backend = backend,
    rng_method = rng_method,
    batch_size = batch_size,
    output_h5 = output_h5,
    activity = activity,
    assay_out = assay_out,
    store_assay = store_assay,
    store_results = store_results,
    tool_name = tool_name,
    return_seurat = return_seurat,
    verbose = verbose
  )
}

#' @title Run SecAct cell-cell communication analysis
#'
#' @description
#' Run SecAct cell-cell communication modules for Seurat scRNA-seq input or
#' SpaCET single-cell-resolution spatial transcriptomics input.
#'
#' @md
#' @inheritParams RunSecAct
#' @param mode `scRNAseq` or `scST`.
#' @param condition_meta Metadata column containing condition labels for
#' scRNA-seq CCC.
#' @param conditionCase,conditionControl Case and control labels for scRNA-seq
#' CCC.
#' @param act_diff_cutoff,exp_logFC_cutoff,exp_mean_all_cutoff,exp_fraction_case_cutoff,padj_cutoff
#' Cutoffs passed to `SecAct.CCC.scRNAseq`.
#' @param is.group.sig,is.group.cor,lambda,nrand Parameters passed to SecAct
#' signature grouping and randomization routines.
#' @param radius,ratio_cutoff,coreNo Parameters passed to `SecAct.CCC.scST`.
#'
#' @return A Seurat or SpaCET object with SecAct CCC results.
#' @export
RunSecActCCC <- function(
  srt = NULL,
  inputProfile = NULL,
  mode = c("scRNAseq", "scST"),
  cellType_meta,
  condition_meta = NULL,
  conditionCase = NULL,
  conditionControl = NULL,
  scale.factor = NULL,
  act_diff_cutoff = 2,
  exp_logFC_cutoff = 0.2,
  exp_mean_all_cutoff = 2,
  exp_fraction_case_cutoff = 0.1,
  padj_cutoff = 0.01,
  sigMatrix = "SecAct",
  is.group.sig = TRUE,
  is.group.cor = 0.9,
  lambda = 5e5,
  nrand = 1000,
  radius = 20,
  ratio_cutoff = 0.2,
  coreNo = 6,
  tool_name = "SecAct_CCC",
  store_results = TRUE,
  verbose = TRUE
) {
  mode <- match.arg(mode)
  secact_check_r(verbose = verbose)
  if (missing(cellType_meta) || is.null(cellType_meta) || !nzchar(cellType_meta)) {
    log_message("{.arg cellType_meta} is required", message_type = "error")
  }

  if (identical(mode, "scRNAseq")) {
    obj <- inputProfile %||% srt
    secact_assert_seurat(obj)
    secact_check_meta_columns(obj, c(cellType_meta, condition_meta))
    if (is.null(condition_meta)) {
      if (!is.null(conditionCase) || !is.null(conditionControl)) {
        warning(
          "conditionCase and conditionControl are ignored when condition_meta is NULL",
          call. = FALSE
        )
      }
      conditionCase <- NULL
      conditionControl <- NULL
    } else {
      secact_assert_scalar_string(condition_meta, "condition_meta")
      secact_assert_scalar_string(conditionCase, "conditionCase")
      secact_assert_scalar_string(conditionControl, "conditionControl")
    }
    fun <- get_namespace_fun("SecAct", "SecAct.CCC.scRNAseq")
    log_message("Running {.pkg SecAct} scRNA-seq CCC...", verbose = verbose)
    out <- fun(
      Seurat_obj = obj,
      cellType_meta = cellType_meta,
      condition_meta = condition_meta,
      conditionCase = conditionCase,
      conditionControl = conditionControl,
      scale.factor = scale.factor %||% 1e5,
      act_diff_cutoff = act_diff_cutoff,
      exp_logFC_cutoff = exp_logFC_cutoff,
      exp_mean_all_cutoff = exp_mean_all_cutoff,
      exp_fraction_case_cutoff = exp_fraction_case_cutoff,
      padj_cutoff = padj_cutoff,
      sigMatrix = sigMatrix,
      is.group.sig = is.group.sig,
      is.group.cor = is.group.cor,
      lambda = lambda,
      nrand = nrand
    )
    if (isTRUE(store_results)) {
      out@tools[[tool_name]] <- list(
        method = "SecAct.CCC.scRNAseq",
        parameters = as.list(match.call())[-1]
      )
    }
    return(out)
  }

  obj <- inputProfile %||% srt
  secact_assert_spacet(obj)
  fun <- get_namespace_fun("SecAct", "SecAct.CCC.scST")
  log_message("Running {.pkg SecAct} scST CCC...", verbose = verbose)
  fun(
    SpaCET_obj = obj,
    cellType_meta = cellType_meta,
    scale.factor = scale.factor %||% 1000,
    radius = radius,
    ratio_cutoff = ratio_cutoff,
    padj_cutoff = padj_cutoff,
    coreNo = coreNo
  )
}

#' @title Run SecAct spatial signaling pattern analysis
#'
#' @description
#' Run `SecAct.signaling.pattern` on a SpaCET object with existing SecAct
#' activity results.
#'
#' @md
#' @inheritParams RunSecAct
#' @param SpaCET_obj A SpaCET object.
#' @param radius Radius cutoff.
#' @param k Number of NMF patterns, or candidate pattern numbers.
#'
#' @return A SpaCET object with SecAct signaling pattern results.
#' @export
RunSecActSignalingPattern <- function(
  SpaCET_obj,
  scale.factor = 1e5,
  radius = 200,
  k,
  verbose = TRUE
) {
  secact_check_r(verbose = verbose)
  secact_assert_spacet(SpaCET_obj)
  if (missing(k) || is.null(k)) {
    log_message("{.arg k} is required", message_type = "error")
  }
  fun <- get_namespace_fun("SecAct", "SecAct.signaling.pattern")
  log_message("Running {.pkg SecAct} spatial signaling pattern analysis...", verbose = verbose)
  fun(
    SpaCET_obj = SpaCET_obj,
    scale.factor = scale.factor,
    radius = radius,
    k = k
  )
}

#' @title Extract SecAct pattern-associated secreted proteins
#'
#' @description
#' Run `SecAct.signaling.pattern.gene` on a SpaCET object after
#' [RunSecActSignalingPattern()].
#'
#' @md
#' @param SpaCET_obj A SpaCET object.
#' @param n Pattern index.
#' @inheritParams thisutils::log_message
#'
#' @return A matrix of pattern-associated secreted proteins.
#' @export
RunSecActPatternGenes <- function(
  SpaCET_obj,
  n,
  verbose = TRUE
) {
  secact_check_r(verbose = verbose)
  secact_assert_spacet(SpaCET_obj)
  if (missing(n) || length(n) != 1L || is.na(n)) {
    log_message("{.arg n} must be a single pattern index", message_type = "error")
  }
  fun <- get_namespace_fun("SecAct", "SecAct.signaling.pattern.gene")
  fun(SpaCET_obj = SpaCET_obj, n = n)
}

#' @title Run SecAct spatial signaling velocity
#'
#' @description
#' Run SecAct signaling velocity for spot-level ST or single-cell-resolution ST
#' SpaCET objects.
#'
#' @md
#' @inheritParams RunSecAct
#' @param SpaCET_obj A SpaCET object.
#' @param mode `spotST` calls `SecAct.signaling.velocity.spotST`; `scST` calls
#' `SecAct.signaling.velocity.scST`.
#' @param gene Secreted protein gene used by spot-level ST velocity.
#' @param signalMode `receiving` or `sending` for spot-level ST.
#' @param radius Spatial radius. Defaults to `200` for `spotST` and `20` for
#' `scST` when `NULL`.
#' @param contourMap,contourBins,animated Plot options for spot-level ST.
#' @param sender,secretedProtein,receiver,cellType_meta Parameters for scST
#' velocity.
#' @param CustomizedAreaCoordinates,show.coordinates,colors,pointSize,pointAlpha,legend.position,legend.size,arrow.color,arrow.width,arrow.size
#' Plot options for scST velocity.
#'
#' @return A ggplot object.
#' @export
RunSecActVelocity <- function(
  SpaCET_obj,
  mode = c("spotST", "scST"),
  scale.factor = 1e5,
  gene = NULL,
  signalMode = "receiving",
  radius = NULL,
  contourMap = FALSE,
  contourBins = 11,
  animated = FALSE,
  sender = NULL,
  secretedProtein = NULL,
  receiver = NULL,
  cellType_meta = NULL,
  CustomizedAreaCoordinates = NULL,
  show.coordinates = TRUE,
  colors = NULL,
  pointSize = 1,
  pointAlpha = 1,
  legend.position = "right",
  legend.size = 1,
  arrow.color = "#ff0099",
  arrow.width = 1,
  arrow.size = 0.3,
  verbose = TRUE
) {
  mode <- match.arg(mode)
  secact_check_r(verbose = verbose)
  secact_assert_spacet(SpaCET_obj)

  if (identical(mode, "spotST")) {
    signalMode <- match.arg(signalMode, c("receiving", "sending"))
    radius <- radius %||% 200
    secact_assert_scalar_string(gene, "gene")
    fun <- get_namespace_fun("SecAct", "SecAct.signaling.velocity.spotST")
    return(fun(
      SpaCET_obj = SpaCET_obj,
      scale.factor = scale.factor,
      gene = gene,
      signalMode = signalMode,
      radius = radius,
      contourMap = contourMap,
      contourBins = contourBins,
      animated = animated
    ))
  }

  secact_assert_scalar_string(sender, "sender")
  secact_assert_scalar_string(secretedProtein, "secretedProtein")
  secact_assert_scalar_string(receiver, "receiver")
  secact_assert_scalar_string(cellType_meta, "cellType_meta")
  radius <- radius %||% 20
  colors <- colors %||% c("#1f78b4", "#e31a1c", "#33a02c", "#ff7f00")
  fun <- get_namespace_fun("SecAct", "SecAct.signaling.velocity.scST")
  fun(
    SpaCET_obj = SpaCET_obj,
    sender = sender,
    secretedProtein = secretedProtein,
    receiver = receiver,
    cellType_meta = cellType_meta,
    scale.factor = scale.factor,
    CustomizedAreaCoordinates = CustomizedAreaCoordinates,
    radius = radius,
    show.coordinates = show.coordinates,
    colors = colors,
    pointSize = pointSize,
    pointAlpha = pointAlpha,
    legend.position = legend.position,
    legend.size = legend.size,
    arrow.color = arrow.color,
    arrow.width = arrow.width,
    arrow.size = arrow.size
  )
}

secact_check_r <- function(verbose = TRUE) {
  check_r("data2intelligence/SecAct", verbose = verbose)
  invisible(TRUE)
}

secact_resolve_mode <- function(srt = NULL, inputProfile = NULL) {
  obj <- inputProfile %||% srt
  if (inherits(obj, "Seurat")) {
    return("scRNAseq")
  }
  if (inherits(obj, "SpaCET")) {
    return("ST")
  }
  "matrix"
}

secact_run_scrnaseq <- function(
  srt,
  inputProfile,
  assay,
  cellType_meta,
  is.singleCellLevel,
  sigMatrix,
  is.filter.sig,
  is.group.sig,
  is.group.cor,
  lambda,
  nrand,
  ncores,
  backend,
  rng_method,
  activity,
  assay_out,
  store_assay,
  store_results,
  tool_name,
  return_seurat,
  verbose
) {
  obj <- inputProfile %||% srt
  secact_assert_seurat(obj)
  assay <- assay %||% SeuratObject::DefaultAssay(obj)
  secact_check_assay(obj, assay)
  secact_require_rna_assay(assay)
  DefaultAssay(obj) <- assay
  if (!isTRUE(is.singleCellLevel)) {
    secact_assert_scalar_string(cellType_meta, "cellType_meta")
    secact_check_meta_columns(obj, cellType_meta)
  }

  fun <- get_namespace_fun("SecAct", "SecAct.activity.inference.scRNAseq")
  log_message(
    "Running {.pkg SecAct} scRNA-seq activity inference...",
    verbose = verbose
  )
  out <- fun(
    inputProfile = obj,
    cellType_meta = cellType_meta,
    is.singleCellLevel = is.singleCellLevel,
    sigMatrix = sigMatrix,
    is.filter.sig = is.filter.sig,
    is.group.sig = is.group.sig,
    is.group.cor = is.group.cor,
    lambda = lambda,
    nrand = nrand,
    ncores = ncores,
    backend = backend,
    rng_method = rng_method
  )
  activity_res <- secact_extract_activity(out)
  if (!isTRUE(return_seurat)) {
    return(activity_res)
  }
  if (isTRUE(store_assay)) {
    out <- secact_store_activity_assay(
      srt = out,
      activity_res = activity_res,
      activity = activity,
      assay_out = assay_out,
      verbose = verbose
    )
  }
  if (isTRUE(store_results)) {
    out@tools[[tool_name]] <- list(
      method = "SecAct.activity.inference.scRNAseq",
      activity = activity_res,
      parameters = list(
        mode = "scRNAseq",
        assay = assay,
        cellType_meta = cellType_meta,
        is.singleCellLevel = is.singleCellLevel,
        sigMatrix = sigMatrix,
        is.filter.sig = is.filter.sig,
        is.group.sig = is.group.sig,
        is.group.cor = is.group.cor,
        lambda = lambda,
        nrand = nrand,
        ncores = ncores,
        backend = backend,
        rng_method = rng_method,
        assay_out = assay_out
      )
    )
  }
  out
}

secact_run_st <- function(
  inputProfile,
  inputProfile_control,
  scale.factor,
  sigMatrix,
  is.filter.sig,
  is.group.sig,
  is.group.cor,
  lambda,
  nrand,
  ncores,
  backend,
  rng_method,
  verbose
) {
  secact_assert_spacet(inputProfile)
  fun <- get_namespace_fun("SecAct", "SecAct.activity.inference.ST")
  log_message("Running {.pkg SecAct} spatial activity inference...", verbose = verbose)
  fun(
    inputProfile = inputProfile,
    inputProfile_control = inputProfile_control,
    scale.factor = scale.factor,
    sigMatrix = sigMatrix,
    is.filter.sig = is.filter.sig,
    is.group.sig = is.group.sig,
    is.group.cor = is.group.cor,
    lambda = lambda,
    nrand = nrand,
    ncores = ncores,
    backend = backend,
    rng_method = rng_method
  )
}

secact_run_matrix <- function(
  srt,
  inputProfile,
  inputProfile_control,
  assay,
  layer,
  is.differential,
  is.paired,
  is.singleSampleLevel,
  sigMatrix,
  is.filter.sig,
  is.group.sig,
  is.group.cor,
  lambda,
  nrand,
  ncores,
  backend,
  rng_method,
  batch_size,
  output_h5,
  activity,
  assay_out,
  store_assay,
  store_results,
  tool_name,
  return_seurat,
  verbose
) {
  expr <- inputProfile
  if (is.null(expr)) {
    secact_assert_seurat(srt)
    assay <- assay %||% SeuratObject::DefaultAssay(srt)
    secact_require_rna_assay(assay)
    expr <- GetAssayData5(srt, assay = assay, layer = layer)
  }
  expr <- secact_as_matrix(expr, "inputProfile")
  if (!is.null(inputProfile_control)) {
    inputProfile_control <- secact_as_matrix(inputProfile_control, "inputProfile_control")
  }

  fun <- get_namespace_fun("SecAct", "SecAct.activity.inference")
  log_message(
    "Running {.pkg SecAct} matrix activity inference on {.val {ncol(expr)}} profile{?s}...",
    verbose = verbose
  )
  res <- fun(
    inputProfile = expr,
    inputProfile_control = inputProfile_control,
    is.differential = is.differential,
    is.paired = is.paired,
    is.singleSampleLevel = is.singleSampleLevel,
    sigMatrix = sigMatrix,
    is.filter.sig = is.filter.sig,
    is.group.sig = is.group.sig,
    is.group.cor = is.group.cor,
    lambda = lambda,
    nrand = nrand,
    ncores = ncores,
    backend = backend,
    rng_method = rng_method,
    batch_size = batch_size,
    output_h5 = output_h5
  )
  if (!isTRUE(return_seurat)) {
    return(res)
  }
  secact_assert_seurat(srt)
  if (isTRUE(store_assay)) {
    srt <- secact_store_activity_assay(
      srt = srt,
      activity_res = res,
      activity = activity,
      assay_out = assay_out,
      verbose = verbose
    )
  }
  if (isTRUE(store_results)) {
    srt@tools[[tool_name]] <- list(
      method = "SecAct.activity.inference",
      activity = res,
      parameters = list(
        mode = "matrix",
        assay = assay,
        layer = layer,
        is.differential = is.differential,
        is.paired = is.paired,
        is.singleSampleLevel = is.singleSampleLevel,
        sigMatrix = sigMatrix,
        is.filter.sig = is.filter.sig,
        is.group.sig = is.group.sig,
        is.group.cor = is.group.cor,
        lambda = lambda,
        nrand = nrand,
        ncores = ncores,
        backend = backend,
        rng_method = rng_method,
        batch_size = batch_size,
        output_h5 = output_h5,
        assay_out = assay_out
      )
    )
  }
  srt
}

secact_extract_activity <- function(obj) {
  activity <- NULL
  if (inherits(obj, "Seurat")) {
    activity <- obj@misc$SecAct_output$SecretedProteinActivity
  } else if (inherits(obj, "SpaCET")) {
    activity <- obj@results$SecAct_output$SecretedProteinActivity
  } else if (is.list(obj) && all(c("beta", "se", "zscore", "pvalue") %in% names(obj))) {
    activity <- obj
  }
  if (!is.list(activity) || !"zscore" %in% names(activity)) {
    log_message(
      "{.pkg SecAct} did not return a valid SecretedProteinActivity result",
      message_type = "error"
    )
  }
  activity
}

secact_store_activity_assay <- function(
  srt,
  activity_res,
  activity,
  assay_out,
  verbose = TRUE
) {
  secact_assert_seurat(srt)
  if (!activity %in% names(activity_res)) {
    log_message(
      "Activity result {.val {activity}} was not returned by {.pkg SecAct}",
      message_type = "error"
    )
  }
  mat <- activity_res[[activity]]
  if (!inherits(mat, c("matrix", "Matrix"))) {
    log_message(
      "Cannot store {.pkg SecAct} activity as a Seurat assay because activity result {.val {activity}} is not matrix-like",
      message_type = "warning",
      verbose = verbose
    )
    return(srt)
  }
  mat <- as.matrix(mat)
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    log_message(
      "Cannot store {.pkg SecAct} activity as a Seurat assay because activity result {.val {activity}} lacks row or column names",
      message_type = "warning",
      verbose = verbose
    )
    return(srt)
  }
  cells <- colnames(srt)
  if (!all(cells %in% colnames(mat))) {
    log_message(
      "Cannot store {.pkg SecAct} activity as a Seurat assay because activity columns do not match all Seurat cells",
      message_type = "warning",
      verbose = verbose
    )
    return(srt)
  }
  mat <- mat[, cells, drop = FALSE]
  mat <- Matrix::Matrix(mat, sparse = TRUE)
  suppressWarnings({
    srt[[assay_out]] <- SeuratObject::CreateAssayObject(data = mat)
  })
  log_message(
    "{.pkg SecAct} {.val {activity}} stored in assay {.val {assay_out}}",
    message_type = "success",
    verbose = verbose
  )
  srt
}

secact_as_matrix <- function(x, arg) {
  if (inherits(x, "Seurat") || inherits(x, "SpaCET")) {
    log_message(
      "{.arg {arg}} is an object; use {.arg mode = 'scRNAseq'} or {.arg mode = 'ST'}",
      message_type = "error"
    )
  }
  if (inherits(x, "data.frame")) {
    x <- as.matrix(x)
  }
  if (!inherits(x, c("matrix", "Matrix"))) {
    log_message("{.arg {arg}} must be a matrix-like object", message_type = "error")
  }
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    log_message(
      "{.arg {arg}} must contain gene row names and profile column names",
      message_type = "error"
    )
  }
  x
}

secact_assert_seurat <- function(x) {
  if (!inherits(x, "Seurat")) {
    log_message("{.arg srt} must be a {.cls Seurat} object", message_type = "error")
  }
  invisible(TRUE)
}

secact_assert_spacet <- function(x) {
  if (!inherits(x, "SpaCET")) {
    log_message("{.arg SpaCET_obj} must be a {.cls SpaCET} object", message_type = "error")
  }
  invisible(TRUE)
}

secact_assert_scalar_string <- function(x, arg) {
  if (is.null(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    log_message(
      "{.arg {arg}} must be a single non-empty string",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

secact_check_assay <- function(srt, assay) {
  if (!assay %in% names(srt@assays)) {
    log_message("Assay {.val {assay}} was not found in {.arg srt}", message_type = "error")
  }
  invisible(TRUE)
}

secact_require_rna_assay <- function(assay) {
  if (!identical(assay, "RNA")) {
    log_message(
      "{.pkg SecAct} Seurat wrappers require {.arg assay = 'RNA'} to avoid reading a non-RNA assay silently",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

secact_check_meta_columns <- function(srt, cols) {
  cols <- cols[!is.null(cols)]
  missing_cols <- setdiff(cols, colnames(srt[[]]))
  if (length(missing_cols) > 0L) {
    log_message("Missing metadata column{?s}: {.val {missing_cols}}", message_type = "error")
  }
  invisible(TRUE)
}
