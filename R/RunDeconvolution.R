#' @title Run bulk or pseudobulk deconvolution
#'
#' @description
#' Estimate cell-type proportions from a bulk-like expression matrix stored in a
#' `SummarizedExperiment` object, using a `Seurat` reference.
#'
#' @md
#' @inheritParams standard_scop
#' @param object A `SummarizedExperiment` object containing bulk-like counts.
#' @param reference A `Seurat` reference object used to build cell-type
#' profiles. Not required for `"CIBERSORT"`.
#' @param method Deconvolution method. One of `"MuSiC"`, `"BisqueRNA"`,
#' `"BayesPrism"`, or `"CIBERSORT"`.
#' @param group.by Metadata column in `reference` defining reference cell
#' types.
#' @param sample.by Metadata column in `reference` defining biological
#' sample / donor IDs. Used by the `r` backends of `MuSiC` and `BisqueRNA`. If
#' `NULL`, SCOP will try to infer a suitable column automatically.
#' @param cellstate.by Metadata column in `reference` defining cell states for
#' the `r` backend of `BayesPrism`. If `NULL`, `group.by` is reused.
#' @param bulk_assay Assay name in `object` used as the bulk counts matrix.
#' @param ref_assay Assay name in `reference` used for the reference profiles.
#' @param ref_layer Layer name in `reference` used for reference counts.
#' @param backend Deconvolution engine backend. `"r"` uses the original method
#' package implementation. `"cpp"` is reserved for native SCOP implementations
#' when available.
#' @param ... Additional parameters forwarded to the internal deconvolution
#' backend.
#'
#' @return A `SummarizedExperiment` object with results stored in
#' `S4Vectors::metadata(object)[["Deconvolution"]]`.
#'
#' @seealso [DeconvolutionPlot]
#'
#' @export
#'
#' @examples
#' data(islet_bulk)
#' islet_bulk <- RunDeconvolution(
#'   islet_bulk,
#'   method = "CIBERSORT",
#'   backend = "cpp",
#'   perm = 0
#' )
#' DeconvolutionPlot(islet_bulk, plot_type = "bar")
#'
#' DeconvolutionPlot(
#'   islet_bulk,
#'   plot_type = "heatmap",
#'   sample_annotation = "condition",
#'   sample_split = "condition"
#' )
#'
#' DeconvolutionPlot(islet_bulk, plot_type = "box")
RunDeconvolution <- function(object, ...) {
  if (methods::is(object, "SummarizedExperiment")) {
    return(RunDeconvolution.SummarizedExperiment(object, ...))
  }
  UseMethod(generic = "RunDeconvolution", object = object)
}

#' @rdname RunDeconvolution
#' @export
RunDeconvolution.SummarizedExperiment <- function(
  object,
  reference = NULL,
  method = c("MuSiC", "BisqueRNA", "BayesPrism", "CIBERSORT"),
  group.by = NULL,
  sample.by = NULL,
  cellstate.by = NULL,
  bulk_assay = "counts",
  ref_assay = NULL,
  ref_layer = "counts",
  backend = c("cpp", "r"),
  verbose = TRUE,
  ...
) {
  method <- match.arg(method)
  if (missing(backend)) {
    backend <- if (identical(method, "BayesPrism")) "cpp" else "r"
  } else {
    backend <- match.arg(backend)
  }
  ctx <- build_context(
    mode = "pure_bulk",
    bulk_se = object,
    bulk_assay = bulk_assay
  )
  method_info <- canonical_method(method)
  if (!identical(method_info$module, "deconv")) {
    log_message(
      "{.arg method} must resolve to a supported deconvolution method.",
      message_type = "error"
    )
  }
  if (!identical(method_info$method, "deconv_CIBERSORT")) {
    ctx$reference <- build_reference_profiles(
      ref_srt = reference,
      group.by = group.by,
      sample.by = sample.by,
      cellstate.by = cellstate.by,
      assay = ref_assay,
      layer = ref_layer
    )
  }
  bundle <- run_deconv(
    ctx = ctx,
    method_name = method_info$method,
    deconv_args = utils::modifyList(
      list(backend = backend),
      list(...)
    ),
    verbose = verbose
  )
  bundle$results <- deconv_schema(bundle$results)

  store <- list(
    input = list(
      bulk_assay = bulk_assay,
      ref_assay = if (is.null(reference)) {
        ref_assay
      } else {
        ref_assay %||% SeuratObject::DefaultAssay(reference)
      },
      ref_layer = ref_layer,
      group.by = group.by,
      sample.by = sample.by,
      cellstate.by = cellstate.by,
      backend = backend
    ),
    active_method = method_info$method,
    methods = stats::setNames(list(bundle), method_info$method),
    results = bundle$results,
    parameters = bundle$parameters %||% list(),
    details = bundle$details %||% list(),
    status = list(
      method = method_info$method,
      status = bundle$status %||% "failed",
      reason = bundle$reason %||% NULL
    )
  )
  store_meta(object, "Deconvolution", store)
}

canonical_method <- function(method) {
  key <- tolower(as.character(method) %||% "")
  method_map <- list(
    music = list(method = "deconv_MuSiC", module = "deconv"),
    bisquerna = list(method = "deconv_BisqueRNA", module = "deconv"),
    bayesprism = list(method = "deconv_BayesPrism", module = "deconv"),
    cibersort = list(method = "deconv_CIBERSORT", module = "deconv")
  )
  out <- method_map[[key]]
  if (is.null(out)) {
    log_message(
      "{.arg method} must be one of {.val {c('MuSiC', 'BisqueRNA', 'BayesPrism', 'CIBERSORT')}}",
      message_type = "error"
    )
  }
  out
}

build_reference_profiles <- function(
  ref_srt,
  group.by,
  sample.by = NULL,
  cellstate.by = NULL,
  assay = NULL,
  layer = "counts"
) {
  if (!inherits(ref_srt, "Seurat")) {
    log_message(
      "{.arg reference} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (is.null(group.by) || !group.by %in% colnames(ref_srt@meta.data)) {
    log_message(
      "{.arg group.by} must be present in {.arg reference@meta.data}",
      message_type = "error"
    )
  }
  if (!is.null(sample.by) && !sample.by %in% colnames(ref_srt@meta.data)) {
    log_message(
      "{.arg sample.by} must be present in {.arg reference@meta.data}",
      message_type = "error"
    )
  }
  if (
    !is.null(cellstate.by) && !cellstate.by %in% colnames(ref_srt@meta.data)
  ) {
    log_message(
      "{.arg cellstate.by} must be present in {.arg reference@meta.data}",
      message_type = "error"
    )
  }
  list(
    object = ref_srt,
    group.by = group.by,
    sample.by = sample.by,
    cellstate.by = cellstate.by,
    assay = assay,
    layer = layer
  )
}
# Internal deconvolution and CSDE method implementations for
# RunDeconvolution() and related bulk-analysis helpers.
NULL

infer_ref_sample_col <- function(reference_meta, sample.by = NULL) {
  if (!is.null(sample.by)) {
    if (!sample.by %in% colnames(reference_meta)) {
      log_message(
        "{.arg sample.by} must be a valid metadata column in the reference object.",
        message_type = "error"
      )
    }
    return(sample.by)
  }

  candidates <- c(
    "sample",
    "subject",
    "subject_id",
    "donor",
    "patient",
    "orig.ident",
    "replicate",
    "dataset"
  )
  candidates <- intersect(candidates, colnames(reference_meta))
  for (nm in candidates) {
    values <- reference_meta[[nm]]
    values <- values[!is.na(values) & nzchar(as.character(values))]
    if (length(unique(as.character(values))) >= 2) {
      return(nm)
    }
  }

  log_message(
    paste(
      "Cannot infer a reference sample column automatically.",
      "Provide {.arg sample.by} explicitly."
    ),
    message_type = "error"
  )
}

prep_ref <- function(
  reference_srt,
  group.by,
  sample.by = NULL,
  cellstate.by = NULL,
  assay = NULL,
  layer = "counts",
  dense = TRUE
) {
  if (!inherits(reference_srt, "Seurat")) {
    log_message(
      "{.arg reference_srt} must be a {.cls Seurat} object.",
      message_type = "error"
    )
  }
  if (!group.by %in% colnames(reference_srt@meta.data)) {
    log_message(
      "{.arg group.by} must be present in {.arg reference_srt@meta.data}.",
      message_type = "error"
    )
  }

  assay_use <- assay %||% SeuratObject::DefaultAssay(reference_srt)
  ref_counts <- GetAssayData5(
    object = reference_srt,
    assay = assay_use,
    layer = layer
  )
  ref_meta <- reference_srt@meta.data
  sample_col_source <- infer_ref_sample_col(
    ref_meta,
    sample.by = sample.by
  )
  state_col <- cellstate.by %||% group.by
  ref_meta$celltype_scop <- as.character(ref_meta[[group.by]])
  ref_meta$sample_scop <- as.character(ref_meta[[sample_col_source]])
  if (!state_col %in% colnames(ref_meta)) {
    log_message(
      "{.arg cellstate.by} must be a valid metadata column in the reference object.",
      message_type = "error"
    )
  }
  ref_meta$cellstate_scop <- as.character(ref_meta[[state_col]])

  keep <- !is.na(ref_meta$celltype_scop) &
    nzchar(ref_meta$celltype_scop) &
    !is.na(ref_meta$sample_scop) &
    nzchar(ref_meta$sample_scop) &
    !is.na(ref_meta$cellstate_scop) &
    nzchar(ref_meta$cellstate_scop)
  if (sum(keep) == 0) {
    log_message(
      "No valid reference cells remain after filtering required metadata.",
      message_type = "error"
    )
  }

  ref_counts <- ref_counts[, keep, drop = FALSE]
  ref_meta <- ref_meta[keep, , drop = FALSE]
  ref_meta <- ref_meta[colnames(ref_counts), , drop = FALSE]
  if (isTRUE(dense)) {
    ref_counts <- as_matrix(ref_counts)
  }

  list(
    counts = ref_counts,
    meta = ref_meta,
    assay = assay_use,
    layer = layer,
    sample_col = "sample_scop",
    sample_col_source = sample_col_source,
    celltype_col = "celltype_scop",
    cellstate_col = "cellstate_scop",
    cellstate_col_source = state_col
  )
}

as_bulk_eset <- function(count_matrix) {
  count_matrix <- as_matrix(count_matrix)
  sample_meta <- data.frame(row.names = colnames(count_matrix))
  Biobase::ExpressionSet(
    assayData = count_matrix,
    phenoData = Biobase::AnnotatedDataFrame(sample_meta)
  )
}

as_ref_eset <- function(reference_info) {
  Biobase::ExpressionSet(
    assayData = as_matrix(reference_info$counts),
    phenoData = Biobase::AnnotatedDataFrame(reference_info$meta)
  )
}

match_backend <- function(
  backend,
  cpp_available = FALSE,
  method_label = "method"
) {
  backend <- match.arg(backend, c("cpp", "r"))
  if (identical(backend, "cpp") && !isTRUE(cpp_available)) {
    log_message(
      paste0(
        "{.val ",
        method_label,
        "} currently supports only {.arg backend = 'r'}."
      ),
      message_type = "error"
    )
  }
  backend
}

RunMuSiC <- function(
  count_matrix,
  reference_srt = NULL,
  group.by = NULL,
  sample.by = NULL,
  ref_assay = NULL,
  ref_layer = "counts",
  backend = c("cpp", "r"),
  verbose = TRUE
) {
  backend <- match_backend(
    backend = backend,
    cpp_available = FALSE,
    method_label = "MuSiC"
  )
  check_r(c("xuranw/MuSiC", "SingleCellExperiment"), verbose = FALSE)
  ref_info <- prep_ref(
    reference_srt = reference_srt,
    group.by = group.by,
    sample.by = sample.by,
    assay = ref_assay,
    layer = ref_layer,
    dense = FALSE
  )
  sc_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = ref_info$counts),
    colData = S4Vectors::DataFrame(ref_info$meta)
  )
  music_prop <- get_namespace_fun("MuSiC", "music_prop")
  prop_fit <- tryCatch(
    music_prop(
      bulk.mtx = as_matrix(count_matrix),
      sc.sce = sc_sce,
      clusters = ref_info$celltype_col,
      samples = ref_info$sample_col,
      verbose = verbose
    ),
    error = function(e) e
  )
  if (inherits(prop_fit, "error")) {
    return(list(
      status = "failed",
      reason = prop_fit$message,
      results = data.frame(),
      details = list(),
      parameters = list()
    ))
  }
  prop_matrix <- prop_fit$Est.prop.weighted %||% prop_fit$Est.prop.allgene

  list(
    status = "success",
    reason = NULL,
    results = prop_long(prop_matrix, method_name = "MuSiC"),
    details = list(
      engine = "MuSiC::music_prop",
      package_backend = TRUE,
      proportion_matrix = prop_matrix
    ),
    parameters = list(
      method = "MuSiC",
      backend = backend,
      sample.by = sample.by %||% ref_info$sample_col_source
    )
  )
}

RunBisqueRNA <- function(
  count_matrix,
  reference_srt = NULL,
  group.by = NULL,
  sample.by = NULL,
  ref_assay = NULL,
  ref_layer = "counts",
  backend = c("cpp", "r"),
  use.overlap = FALSE,
  old.cpm = TRUE,
  verbose = TRUE
) {
  backend <- match_backend(
    backend = backend,
    cpp_available = FALSE,
    method_label = "BisqueRNA"
  )
  check_r(c("BisqueRNA", "Biobase"), verbose = FALSE)

  ref_info <- prep_ref(
    reference_srt = reference_srt,
    group.by = group.by,
    sample.by = sample.by,
    assay = ref_assay,
    layer = ref_layer
  )
  bulk_eset <- as_bulk_eset(count_matrix)
  sc_eset <- as_ref_eset(ref_info)
  decompose <- get_namespace_fun(
    "BisqueRNA",
    "ReferenceBasedDecomposition"
  )
  prop_fit <- tryCatch(
    decompose(
      bulk.eset = bulk_eset,
      sc.eset = sc_eset,
      cell.types = ref_info$celltype_col,
      subject.names = ref_info$sample_col,
      use.overlap = use.overlap,
      old.cpm = old.cpm,
      verbose = verbose
    ),
    error = function(e) e
  )
  if (inherits(prop_fit, "error")) {
    return(list(
      status = "failed",
      reason = prop_fit$message,
      results = data.frame(),
      details = list(),
      parameters = list()
    ))
  }
  prop_matrix <- t(as_matrix(prop_fit$bulk.props))

  list(
    status = "success",
    reason = NULL,
    results = prop_long(prop_matrix, method_name = "BisqueRNA"),
    details = list(
      engine = "BisqueRNA::ReferenceBasedDecomposition",
      package_backend = TRUE,
      proportion_matrix = prop_matrix
    ),
    parameters = list(
      method = "BisqueRNA",
      backend = backend,
      sample.by = sample.by %||% ref_info$sample_col_source,
      use.overlap = use.overlap,
      old.cpm = old.cpm
    )
  )
}

RunBayesPrism <- function(
  count_matrix,
  reference_srt = NULL,
  group.by = NULL,
  sample.by = NULL,
  cellstate.by = NULL,
  ref_assay = NULL,
  ref_layer = "counts",
  backend = c("cpp", "r"),
  key = NULL,
  outlier.cut = 0.01,
  outlier.fraction = 0.1,
  pseudo.min = 1e-8,
  n.cores = 1,
  update.gibbs = TRUE,
  gibbs.control = list(),
  opt.control = list(),
  verbose = TRUE
) {
  n_cores_missing <- missing(n.cores)
  backend <- match_backend(
    backend = backend,
    cpp_available = TRUE,
    method_label = "BayesPrism"
  )
  check_r("BayesPrism", verbose = FALSE)
  suppressPackageStartupMessages(base::library(
    "BayesPrism",
    character.only = TRUE
  ))
  ref_info <- prep_ref(
    reference_srt = reference_srt,
    group.by = group.by,
    sample.by = sample.by,
    cellstate.by = cellstate.by,
    assay = ref_assay,
    layer = ref_layer
  )
  new_prism <- get_namespace_fun("BayesPrism", "new.prism")
  valid_gibbs_control <- get_namespace_fun("BayesPrism", "valid.gibbs.control")
  valid_opt_control <- get_namespace_fun("BayesPrism", "valid.opt.control")
  get_gibbs_idx <- get_namespace_fun("BayesPrism", "get.gibbs.idx")
  merge_k <- get_namespace_fun("BayesPrism", "mergeK")
  update_reference <- get_namespace_fun("BayesPrism", "updateReference")
  run_prism <- get_namespace_fun("BayesPrism", "run.prism")
  get_fraction <- get_namespace_fun("BayesPrism", "get.fraction")
  prism_obj <- tryCatch(
    new_prism(
      reference = t(as_matrix(ref_info$counts)),
      input.type = "count.matrix",
      cell.type.labels = ref_info$meta[[ref_info$celltype_col]],
      cell.state.labels = ref_info$meta[[ref_info$cellstate_col]],
      key = key,
      mixture = t(as_matrix(count_matrix)),
      outlier.cut = outlier.cut,
      outlier.fraction = outlier.fraction,
      pseudo.min = pseudo.min
    ),
    error = function(e) e
  )
  if (inherits(prism_obj, "error")) {
    return(list(
      status = "failed",
      reason = prism_obj$message,
      results = data.frame(),
      details = list(),
      parameters = list()
    ))
  }

  if (identical(backend, "r")) {
    bp_fit <- tryCatch(
      run_prism(
        prism = prism_obj,
        n.cores = n.cores,
        update.gibbs = update.gibbs,
        gibbs.control = gibbs.control,
        opt.control = opt.control
      ),
      error = function(e) e
    )
    if (inherits(bp_fit, "error")) {
      return(list(
        status = "failed",
        reason = bp_fit$message,
        results = data.frame(),
        details = list(),
        parameters = list()
      ))
    }
    prop_matrix <- get_fraction(
      bp_fit,
      which.theta = if (isTRUE(update.gibbs)) "final" else "first",
      state.or.type = "type"
    )
    details <- list(
      engine = "BayesPrism::run.prism",
      package_backend = TRUE,
      proportion_matrix = prop_matrix
    )
  } else {
    effective_n_cores <- if (isTRUE(n_cores_missing)) {
      max(
        1L,
        min(
          4L,
          as.integer(parallel::detectCores(logical = FALSE) %||% 1L)
        )
      )
    } else {
      as.integer(n.cores)
    }
    if (!"n.cores" %in% names(gibbs.control)) {
      gibbs.control$n.cores <- effective_n_cores
    }
    if (!"n.cores" %in% names(opt.control)) {
      opt.control$n.cores <- effective_n_cores
    }
    opt.control <- valid_opt_control(opt.control)
    gibbs.control <- valid_gibbs_control(gibbs.control)
    if (prism_obj@phi_cellState@pseudo.min == 0) {
      gibbs.control$alpha <- max(1, gibbs.control$alpha)
    }
    if (!all(is.na(prism_obj@key))) {
      return(list(
        status = "failed",
        reason = paste(
          "The BayesPrism cpp backend currently supports only references",
          "without a tumor key."
        ),
        results = data.frame(),
        details = list(),
        parameters = list()
      ))
    }

    gibbs_idx <- as.integer(get_gibbs_idx(gibbs.control))
    initial_cpp <- tryCatch(
      bayesprism_gibbs_initial_cpp(
        mixture = prism_obj@mixture,
        phi = prism_obj@phi_cellState@phi,
        gibbs_idx = gibbs_idx,
        alpha = gibbs.control$alpha,
        seed = as.integer(gibbs.control$seed),
        cores = as.integer(gibbs.control$n.cores)
      ),
      error = function(e) e
    )
    if (inherits(initial_cpp, "error")) {
      return(list(
        status = "failed",
        reason = initial_cpp$message,
        results = data.frame(),
        details = list(),
        parameters = list()
      ))
    }

    initial_z <- initial_cpp$Z
    dim(initial_z) <- c(
      nrow(prism_obj@mixture),
      ncol(prism_obj@mixture),
      nrow(prism_obj@phi_cellState@phi)
    )
    dimnames(initial_z) <- list(
      rownames(prism_obj@mixture),
      colnames(prism_obj@mixture),
      rownames(prism_obj@phi_cellState@phi)
    )
    initial_theta <- as_matrix(initial_cpp$theta)
    initial_theta_cv <- as_matrix(initial_cpp$theta_cv)
    rownames(initial_theta) <- rownames(prism_obj@mixture)
    colnames(initial_theta) <- rownames(prism_obj@phi_cellState@phi)
    rownames(initial_theta_cv) <- rownames(prism_obj@mixture)
    colnames(initial_theta_cv) <- rownames(prism_obj@phi_cellState@phi)
    joint_post_initial <- methods::new(
      "jointPost",
      Z = initial_z,
      theta = initial_theta,
      theta.cv = initial_theta_cv,
      constant = as.numeric(initial_cpp$constant %||% 0)
    )
    joint_post_type <- merge_k(
      jointPost.obj = joint_post_initial,
      map = prism_obj@map
    )

    if (isTRUE(update.gibbs)) {
      psi <- tryCatch(
        update_reference(
          Z = joint_post_type@Z,
          phi_prime = prism_obj@phi_cellType,
          map = prism_obj@map,
          key = prism_obj@key,
          opt.control = opt.control
        ),
        error = function(e) e
      )
      if (inherits(psi, "error")) {
        return(list(
          status = "failed",
          reason = psi$message,
          results = data.frame(),
          details = list(),
          parameters = list()
        ))
      }
      if (!methods::is(psi, "refPhi")) {
        return(list(
          status = "failed",
          reason = "The BayesPrism cpp backend currently supports only refPhi updates.",
          results = data.frame(),
          details = list(),
          parameters = list()
        ))
      }
      final_cpp <- tryCatch(
        bayesprism_gibbs_final_cpp(
          mixture = prism_obj@mixture,
          phi = psi@phi,
          gibbs_idx = gibbs_idx,
          alpha = gibbs.control$alpha,
          seed = as.integer(gibbs.control$seed),
          cores = as.integer(gibbs.control$n.cores)
        ),
        error = function(e) e
      )
      if (inherits(final_cpp, "error")) {
        return(list(
          status = "failed",
          reason = final_cpp$message,
          results = data.frame(),
          details = list(),
          parameters = list()
        ))
      }
      prop_matrix <- as_matrix(final_cpp$theta)
      rownames(prop_matrix) <- rownames(prism_obj@mixture)
      colnames(prop_matrix) <- rownames(psi@phi)
      details <- list(
        engine = "scop::bayesprism_gibbs_cpp",
        package_backend = FALSE,
        proportion_matrix = prop_matrix,
        gibbs_kept_iterations = gibbs_idx,
        updated_reference_n_cell_types = nrow(psi@phi),
        updated_reference_n_genes = ncol(psi@phi)
      )
    } else {
      prop_matrix <- joint_post_type@theta
      details <- list(
        engine = "scop::bayesprism_gibbs_cpp",
        package_backend = FALSE,
        proportion_matrix = prop_matrix,
        gibbs_kept_iterations = gibbs_idx
      )
    }
  }

  list(
    status = "success",
    reason = NULL,
    results = prop_long(prop_matrix, method_name = "BayesPrism"),
    details = details,
    parameters = list(
      method = "BayesPrism",
      backend = backend,
      sample.by = sample.by %||% ref_info$sample_col_source,
      cellstate.by = cellstate.by %||% ref_info$cellstate_col_source,
      key = key,
      outlier.cut = outlier.cut,
      outlier.fraction = outlier.fraction,
      pseudo.min = pseudo.min,
      n.cores = if (identical(backend, "cpp")) {
        as.integer(gibbs.control$n.cores)
      } else {
        n.cores
      },
      update.gibbs = update.gibbs
    )
  )
}

#' @title Run CIBERSORT deconvolution
#'
#' @description
#' Estimate immune cell proportions from a bulk expression matrix using the
#' external `CIBERSORT` package or the native `scop` C++ backend.
#' `sig_matrix = "LM22"` downloads the LM22 signature matrix from
#' `mengxu98/datasets` and caches it locally.
#'
#' @md
#' @inheritParams RunDeconvolution
#' @param object Optional `SummarizedExperiment` object or expression matrix.
#' When a `SummarizedExperiment` is provided, results are stored in
#' `metadata(object)[["Deconvolution"]]`.
#' @param count_matrix Optional expression matrix with genes in rows and samples
#' in columns. Used when `object` is not provided as a matrix.
#' @param sig_matrix Signature matrix, local file path, or `"LM22"`.
#' @param perm Number of CIBERSORT permutations.
#' @param QN Whether CIBERSORT should use quantile normalization.
#' @param absolute Passed to CIBERSORT when supported by the installed package.
#' The native C++ backend currently returns relative fractions.
#' @param backend CIBERSORT backend. `"r"` calls the external `CIBERSORT`
#' package and `"cpp"` uses the native `scop` LIBSVM implementation.
#' @param cores Number of CPU cores used by the C++ backend.
#' @param seed Random seed used by the C++ permutation backend.
#'
#' @return A deconvolution result bundle for matrix input, or the modified
#' `SummarizedExperiment` object for `SummarizedExperiment` input.
#'
#' @export
#'
#' @examples
#' data(islet_bulk)
#'
#' if (FALSE) {
#' # Run CIBERSORT
#' islet_bulk <- RunCIBERSORT(
#'   object = islet_bulk,
#'   sig_matrix = "LM22",
#'   bulk_assay = "counts",
#'   perm = 100,
#'   QN = TRUE
#' )
#'
#' # Immune abundance stacked bar plot
#' p1 <- ImmuneAbundancePlot(
#'   object = islet_bulk,
#'   plot_type = "bar",
#'   group.by = "condition"
#' )
#' p1
#'
#' # Immune cell correlation heatmap
#' p2 <- ImmuneAbundancePlot(
#'   object = islet_bulk,
#'   plot_type = "cor"
#' )
#' p2
#'
#' # Gene-immune correlation butterfly plot
#' p3 <- GeneImmuneCorPlot(
#'   object = islet_bulk,
#'   features = rownames(SummarizedExperiment::assay(islet_bulk, "counts"))[1:3]
#' )
#' p3
#' }
RunCIBERSORT <- function(
  object = NULL,
  count_matrix = NULL,
  sig_matrix = "LM22",
  bulk_assay = "counts",
  perm = 100,
  QN = TRUE,
  absolute = FALSE,
  backend = c("r", "cpp"),
  cores = 1L,
  seed = 123L,
  verbose = TRUE,
  ...
) {
  backend <- match.arg(backend)
  input_object <- object
  if (is.null(count_matrix)) {
    if (methods::is(object, "SummarizedExperiment")) {
      count_matrix <- SummarizedExperiment::assay(object, bulk_assay)
    } else if (inherits(object, c("matrix", "data.frame", "Matrix"))) {
      count_matrix <- object
      input_object <- NULL
    }
  }
  if (is.null(count_matrix)) {
    log_message(
      "{.arg object} must be a {.cls SummarizedExperiment} or expression matrix, or {.arg count_matrix} must be provided.",
      message_type = "error"
    )
  }

  bundle <- run_cibersort_bundle(
    count_matrix = count_matrix,
    sig_matrix = sig_matrix,
    perm = perm,
    QN = QN,
    absolute = absolute,
    backend = backend,
    cores = cores,
    seed = seed,
    verbose = verbose,
    ...
  )

  if (methods::is(input_object, "SummarizedExperiment")) {
    method_name <- "CIBERSORT"
    store <- list(
      input = list(
        bulk_assay = bulk_assay,
        backend = backend
      ),
      active_method = method_name,
      methods = stats::setNames(list(bundle), method_name),
      results = deconv_schema(bundle$results),
      parameters = bundle$parameters %||% list(),
      details = bundle$details %||% list(),
      status = list(
        method = method_name,
        status = bundle$status %||% "failed",
        reason = bundle$reason %||% NULL
      )
    )
    return(store_meta(input_object, "Deconvolution", store))
  }
  bundle
}

run_cibersort_bundle <- function(
  count_matrix,
  sig_matrix = "LM22",
  perm = 100,
  QN = TRUE,
  absolute = FALSE,
  backend = c("r", "cpp"),
  cores = 1L,
  seed = 123L,
  verbose = TRUE,
  ...
) {
  backend <- match.arg(backend)
  perm <- as.integer(perm)
  cores <- as.integer(cores)
  seed <- as.integer(seed)
  if (length(perm) != 1L || is.na(perm) || !is.finite(perm) || perm < 0L) {
    log_message(
      "{.arg perm} must be a non-negative integer.",
      message_type = "error"
    )
  }
  if (
    length(cores) != 1L ||
      is.na(cores) ||
      !is.finite(cores) ||
      cores < 1L
  ) {
    log_message(
      "{.arg cores} must be a positive integer.",
      message_type = "error"
    )
  }
  if (length(seed) != 1L || is.na(seed) || !is.finite(seed)) {
    log_message("{.arg seed} must be an integer.", message_type = "error")
  }
  count_matrix <- cibersort_check_matrix(count_matrix, "count_matrix")
  signature <- resolve_cibersort_signature(sig_matrix, verbose = verbose)
  sig_matrix_use <- cibersort_check_matrix(signature$matrix, "sig_matrix")
  common_genes <- sort(intersect(
    rownames(sig_matrix_use),
    rownames(count_matrix)
  ))
  if (length(common_genes) == 0L) {
    log_message(
      "No shared genes between {.arg sig_matrix} and {.arg count_matrix}.",
      message_type = "error"
    )
  }
  if (length(common_genes) < nrow(sig_matrix_use)) {
    log_message(
      "Use {.val {length(common_genes)}} shared genes for CIBERSORT",
      verbose = verbose
    )
  }
  sig_matrix_use <- sig_matrix_use[common_genes, , drop = FALSE]
  count_matrix <- count_matrix[common_genes, , drop = FALSE]

  if (identical(backend, "cpp")) {
    fit <- tryCatch(
      cibersort_cpp(
        signature = as_matrix(sig_matrix_use),
        mixture = as_matrix(count_matrix),
        perm = perm,
        QN = isTRUE(QN),
        absolute = isTRUE(absolute),
        cores = cores,
        seed = seed,
        verbose = verbose
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      return(list(
        status = "failed",
        reason = fit$message,
        results = data.frame(),
        details = list(
          engine = "scop::cibersort_cpp",
          package_backend = FALSE,
          signature_source = signature$source
        ),
        parameters = list(
          method = "CIBERSORT",
          backend = backend,
          perm = perm,
          QN = QN,
          absolute = absolute,
          cores = cores,
          seed = seed,
          sig_matrix = signature$label
        )
      ))
    }

    parsed <- parse_cibersort_result(cbind(
      as.data.frame(fit$proportion_matrix, check.names = FALSE),
      as.data.frame(fit$statistics, check.names = FALSE)
    ))
    return(list(
      status = "success",
      reason = NULL,
      results = prop_long(parsed$proportion_matrix, method_name = "CIBERSORT"),
      details = list(
        engine = "scop::cibersort_cpp",
        package_backend = FALSE,
        signature_source = signature$source,
        proportion_matrix = parsed$proportion_matrix,
        statistics = parsed$statistics,
        raw_results = parsed$raw_results
      ),
      parameters = list(
        method = "CIBERSORT",
        backend = backend,
        perm = perm,
        QN = QN,
        absolute = absolute,
        cores = cores,
        seed = seed,
        sig_matrix = signature$label
      )
    ))
  }

  if (!requireNamespace("CIBERSORT", quietly = TRUE)) {
    log_message(
      paste(
        "{.pkg CIBERSORT} is required for {.fn RunCIBERSORT} with {.arg backend = 'r'}.",
        "Install it with {.code devtools::install_github('Moonerss/CIBERSORT')}",
        "or use {.arg backend = 'cpp'}."
      ),
      message_type = "error"
    )
  }

  cibersort_fun <- get_namespace_fun("CIBERSORT", "cibersort")
  call_args <- utils::modifyList(
    list(
      sig_matrix = sig_matrix_use,
      mixture_file = as_matrix(count_matrix),
      perm = perm,
      QN = QN,
      absolute = absolute
    ),
    list(...)
  )
  call_args <- call_args[names(call_args) %in% names(formals(cibersort_fun))]
  fit <- tryCatch(
    do.call(cibersort_fun, call_args),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(
      status = "failed",
      reason = fit$message,
      results = data.frame(),
      details = list(
        engine = "CIBERSORT::cibersort",
        package_backend = TRUE,
        signature_source = signature$source
      ),
      parameters = list(
        method = "CIBERSORT",
        backend = backend,
        perm = perm,
        QN = QN,
        absolute = absolute,
        cores = cores,
        seed = seed,
        sig_matrix = signature$label
      )
    ))
  }

  parsed <- parse_cibersort_result(fit)
  list(
    status = "success",
    reason = NULL,
    results = prop_long(parsed$proportion_matrix, method_name = "CIBERSORT"),
    details = list(
      engine = "CIBERSORT::cibersort",
      package_backend = TRUE,
      signature_source = signature$source,
      proportion_matrix = parsed$proportion_matrix,
      statistics = parsed$statistics,
      raw_results = parsed$raw_results
    ),
    parameters = list(
      method = "CIBERSORT",
      backend = backend,
      perm = perm,
      QN = QN,
      absolute = absolute,
      cores = cores,
      seed = seed,
      sig_matrix = signature$label
    )
  )
}

cibersort_check_matrix <- function(x, arg_name) {
  if (is.null(x)) {
    log_message(
      "{.arg {arg_name}} must be provided.",
      message_type = "error"
    )
  }
  x <- as_matrix(x)
  dim_names <- dimnames(x)
  x <- suppressWarnings(matrix(
    as.numeric(x),
    nrow = nrow(x),
    ncol = ncol(x),
    dimnames = dim_names
  ))
  if (is.null(rownames(x)) || any(!nzchar(rownames(x)))) {
    log_message(
      "{.arg {arg_name}} must have gene names in rownames.",
      message_type = "error"
    )
  }
  if (is.null(colnames(x)) || any(!nzchar(colnames(x)))) {
    log_message(
      "{.arg {arg_name}} must have column names.",
      message_type = "error"
    )
  }
  x[!is.finite(x)] <- 0
  x
}

resolve_cibersort_signature <- function(sig_matrix = "LM22", verbose = TRUE) {
  if (inherits(sig_matrix, c("matrix", "data.frame", "Matrix"))) {
    return(list(
      matrix = sig_matrix,
      label = "custom",
      source = "user-provided matrix"
    ))
  }
  if (!is.character(sig_matrix) || length(sig_matrix) != 1L) {
    log_message(
      "{.arg sig_matrix} must be a matrix, data frame, local path, or {.val LM22}.",
      message_type = "error"
    )
  }
  if (file.exists(sig_matrix)) {
    mat <- if (grepl("\\.rds$", sig_matrix, ignore.case = TRUE)) {
      readRDS(sig_matrix)
    } else {
      utils::read.delim(sig_matrix, row.names = 1, check.names = FALSE)
    }
    return(list(
      matrix = mat,
      label = normalizePath(sig_matrix, mustWork = FALSE),
      source = normalizePath(sig_matrix, mustWork = FALSE)
    ))
  }
  if (!identical(toupper(sig_matrix), "LM22")) {
    log_message(
      "{.arg sig_matrix} is not a readable file and is not {.val LM22}.",
      message_type = "error"
    )
  }

  cache_dir <- file.path(tools::R_user_dir("scop", "data"), "CIBERSORT")
  cache_file <- file.path(cache_dir, "LM22.rds")
  url <- "https://raw.githubusercontent.com/mengxu98/datasets/main/CIBERSORT/LM22.rds"
  if (!file.exists(cache_file)) {
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    }
    tmp <- tempfile(fileext = ".rds")
    ok <- tryCatch(
      {
        utils::download.file(url, tmp, mode = "wb", quiet = !verbose)
        file.exists(tmp) && file.info(tmp)$size > 0
      },
      error = function(e) FALSE
    )
    if (!isTRUE(ok)) {
      log_message(
        paste(
          "Failed to download LM22 from {.url {url}}.",
          "Add {.file CIBERSORT/LM22.rds} to {.url https://github.com/mengxu98/datasets}",
          "or pass a local/custom {.arg sig_matrix}."
        ),
        message_type = "error"
      )
    }
    file.copy(tmp, cache_file, overwrite = TRUE)
  }
  list(
    matrix = readRDS(cache_file),
    label = "LM22",
    source = cache_file
  )
}

parse_cibersort_result <- function(fit) {
  df <- as.data.frame(fit, check.names = FALSE)
  if ("Mixture" %in% colnames(df)) {
    rownames(df) <- as.character(df$Mixture)
    df$Mixture <- NULL
  }
  if (
    is.null(rownames(df)) ||
      all(rownames(df) %in% as.character(seq_len(nrow(df))))
  ) {
    log_message(
      "{.pkg CIBERSORT} results must contain sample names in rownames or a {.field Mixture} column.",
      message_type = "error"
    )
  }
  stat_cols <- intersect(
    c("P-value", "P.value", "PValue", "Correlation", "RMSE"),
    colnames(df)
  )
  prop_cols <- setdiff(colnames(df), stat_cols)
  numeric_prop <- vapply(df[, prop_cols, drop = FALSE], is.numeric, logical(1))
  prop_cols <- prop_cols[numeric_prop]
  prop_matrix <- as_matrix(df[, prop_cols, drop = FALSE])
  prop_matrix[!is.finite(prop_matrix)] <- 0
  list(
    proportion_matrix = prop_matrix,
    statistics = df[, stat_cols, drop = FALSE],
    raw_results = df
  )
}

RunTOAST <- function(
  count_matrix,
  condition,
  proportions,
  condition1 = NULL,
  condition2 = NULL,
  backend = c("cpp", "r"),
  p.adjust.method = "BH",
  min_prop_variance = 1e-8,
  transform = c("log1p", "none"),
  var_shrinkage = TRUE,
  verbose = TRUE
) {
  backend <- match_backend(
    backend = backend,
    cpp_available = FALSE,
    method_label = "TOAST"
  )
  transform <- match.arg(transform)
  check_r("TOAST", verbose = FALSE)
  suppressPackageStartupMessages(base::library("TOAST", character.only = TRUE))

  pair <- tryCatch(
    resolve_condition_pair(
      condition = condition,
      condition1 = condition1,
      condition2 = condition2,
      strict_two_levels = TRUE
    ),
    error = function(e) e
  )
  if (inherits(pair, "error")) {
    return(list(
      status = "failed",
      reason = pair$message,
      results = data.frame(),
      details = list(),
      parameters = list()
    ))
  }

  cond <- as.character(condition)
  names(cond) <- names(condition) %||% colnames(count_matrix)
  sample_keep <- intersect(
    colnames(count_matrix)[
      cond[colnames(count_matrix)] %in% c(pair$condition1, pair$condition2)
    ],
    rownames(proportions)
  )
  if (length(sample_keep) < 4) {
    return(list(
      status = "failed",
      reason = "TOAST requires at least 4 matched samples between counts and proportions.",
      details = list(condition_pair = pair),
      results = data.frame(),
      parameters = list()
    ))
  }

  count_use <- as_matrix(count_matrix[, sample_keep, drop = FALSE])
  cond_use <- factor(
    cond[sample_keep],
    levels = c(pair$condition1, pair$condition2)
  )
  if (any(table(cond_use) < 2)) {
    return(list(
      status = "failed",
      reason = "TOAST requires at least 2 samples per condition.",
      details = list(condition_pair = pair),
      results = data.frame(),
      parameters = list()
    ))
  }
  prop_use <- as_matrix(proportions[sample_keep, , drop = FALSE])
  prop_use <- prop_use[,
    order(colMeans(prop_use, na.rm = TRUE), decreasing = TRUE),
    drop = FALSE
  ]
  celltype_original <- colnames(prop_use)
  celltype_safe <- make.names(celltype_original, unique = TRUE)
  colnames(prop_use) <- celltype_safe
  names(celltype_original) <- celltype_safe
  expr_use <- if (identical(transform, "log1p")) {
    log1p(count_use)
  } else {
    count_use
  }

  makeDesign <- get_namespace_fun("TOAST", "makeDesign")
  fitModel <- get_namespace_fun("TOAST", "fitModel")
  csTest <- get_namespace_fun("TOAST", "csTest")
  design_df <- data.frame(
    condition = factor(cond_use, levels = c(pair$condition1, pair$condition2))
  )
  prop_fit_use <- prop_use
  fit_error <- NULL
  dropped_cell_types <- character()
  repeat {
    fitted_model <- tryCatch(
      fitModel(
        Design_out = makeDesign(design = design_df, Prop = prop_fit_use),
        Y = expr_use
      ),
      error = function(e) e
    )
    if (!inherits(fitted_model, "error")) {
      break
    }
    fit_error <- fitted_model
    if (ncol(prop_fit_use) <= 1) {
      break
    }
    drop_ct <- colnames(prop_fit_use)[ncol(prop_fit_use)]
    dropped_cell_types <- c(
      dropped_cell_types,
      unname(celltype_original[[drop_ct]])
    )
    prop_fit_use <- prop_fit_use[,
      seq_len(ncol(prop_fit_use) - 1),
      drop = FALSE
    ]
  }
  if (inherits(fitted_model, "error")) {
    return(list(
      status = "failed",
      reason = fit_error$message,
      details = list(
        condition_pair = pair,
        dropped_cell_types = dropped_cell_types
      ),
      results = data.frame(),
      parameters = list()
    ))
  }
  csde_results <- list()
  for (ct in colnames(prop_fit_use)) {
    prop_vec <- as.numeric(prop_fit_use[, ct])
    if (stats::var(prop_vec, na.rm = TRUE) <= min_prop_variance) {
      next
    }
    tt <- tryCatch(
      suppressWarnings(
        csTest(
          fitted_model = fitted_model,
          coef = "condition",
          cell_type = ct,
          var_shrinkage = var_shrinkage,
          verbose = verbose,
          sort = FALSE
        )
      ),
      error = function(e) e
    )
    if (inherits(tt, "error")) {
      next
    }
    res_table <- if (is.data.frame(tt)) tt else tt$res_table
    if (is.null(res_table) || nrow(res_table) == 0) {
      next
    }
    out_df <- data.frame(
      gene = rownames(res_table),
      cell_type = unname(celltype_original[[ct]]),
      group1 = pair$condition2,
      group2 = pair$condition1,
      effect = res_table$effect_size %||% res_table$beta %||% NA_real_,
      p_val = res_table$p_value %||% NA_real_,
      p_val_adj = res_table$fdr %||%
        stats::p.adjust(res_table$p_value, method = p.adjust.method),
      method = "TOAST",
      stringsAsFactors = FALSE
    )
    csde_results[[ct]] <- out_df
  }

  if (length(csde_results) == 0) {
    return(list(
      status = "failed",
      reason = "No valid CSDE result was produced for method 'TOAST'.",
      details = list(condition_pair = pair),
      results = data.frame(
        gene = character(),
        cell_type = character(),
        group1 = character(),
        group2 = character(),
        effect = numeric(),
        p_val = numeric(),
        p_val_adj = numeric(),
        method = character(),
        stringsAsFactors = FALSE
      ),
      parameters = list()
    ))
  }

  res <- do.call(rbind, csde_results)
  rownames(res) <- NULL
  return(list(
    status = "success",
    reason = NULL,
    results = res,
    details = list(
      condition_pair = pair,
      engine = "TOAST::fitModel+csTest",
      package_backend = TRUE,
      retained_cell_types = unname(celltype_original[colnames(prop_fit_use)]),
      dropped_cell_types = dropped_cell_types
    ),
    parameters = list(
      backend = backend,
      p.adjust.method = p.adjust.method,
      min_prop_variance = min_prop_variance,
      transform = transform,
      var_shrinkage = var_shrinkage
    )
  ))
}

run_deconv <- function(
  ctx,
  method_name,
  deconv_args = list(),
  verbose = TRUE
) {
  method_map <- list(
    deconv_MuSiC = RunMuSiC,
    deconv_BisqueRNA = RunBisqueRNA,
    deconv_BayesPrism = RunBayesPrism,
    deconv_CIBERSORT = RunCIBERSORT
  )
  method_fun <- method_map[[method_name]]
  if (is.null(method_fun)) {
    log_message(
      paste0("Unsupported deconvolution method: ", method_name),
      message_type = "error"
    )
  }

  args <- utils::modifyList(
    list(
      count_matrix = ctx$counts_global,
      reference_srt = ctx$reference$object %||% NULL,
      group.by = ctx$reference$group.by %||% NULL,
      sample.by = ctx$reference$sample.by %||% NULL,
      cellstate.by = ctx$reference$cellstate.by %||% NULL,
      ref_assay = ctx$reference$assay %||% NULL,
      ref_layer = ctx$reference$layer %||% "counts",
      verbose = verbose
    ),
    deconv_args
  )
  bundle <- tryCatch(
    invoke_fun(
      method_fun,
      args[names(args) %in% names(formals(method_fun))]
    ),
    error = function(e) {
      list(
        status = "failed",
        reason = e$message,
        results = data.frame(),
        details = list(),
        parameters = list()
      )
    }
  )

  if (!is.list(bundle) || is.null(bundle$status)) {
    bundle <- list(
      status = "failed",
      reason = "Method did not return a valid bundle.",
      results = data.frame(),
      details = list(),
      parameters = list()
    )
  }
  bundle$parameters <- utils::modifyList(
    list(method = method_name),
    bundle$parameters %||% list()
  )
  bundle$results <- deconv_schema(bundle$results)
  if (nrow(bundle$results) > 0) {
    bundle$results$method <- method_name
  }
  bundle
}

run_csde <- function(
  ctx,
  method_name,
  deconv_bundle,
  condition1 = NULL,
  condition2 = NULL,
  csde_args = list(),
  verbose = TRUE
) {
  method_map <- list(
    csde_TOAST = RunTOAST
  )
  method_fun <- method_map[[method_name]]
  if (is.null(method_fun)) {
    log_message(
      paste0("Unsupported CSDE method: ", method_name),
      message_type = "error"
    )
  }

  prop_matrix <- deconv_bundle$details$proportion_matrix %||%
    deconv_mat(deconv_bundle$results)
  if (
    is.null(prop_matrix) || nrow(prop_matrix) == 0 || ncol(prop_matrix) == 0
  ) {
    return(list(
      status = "failed",
      reason = "CSDE requires successful deconvolution results.",
      parameters = list(method = method_name),
      results = data.frame(
        gene = character(),
        cell_type = character(),
        group1 = character(),
        group2 = character(),
        effect = numeric(),
        p_val = numeric(),
        p_val_adj = numeric(),
        method = character(),
        stringsAsFactors = FALSE
      ),
      details = list()
    ))
  }

  args <- utils::modifyList(
    list(
      count_matrix = ctx$counts_global,
      condition = ctx$condition_global,
      proportions = prop_matrix,
      condition1 = condition1,
      condition2 = condition2,
      verbose = verbose
    ),
    csde_args
  )
  bundle <- tryCatch(
    invoke_fun(
      method_fun,
      args[names(args) %in% names(formals(method_fun))]
    ),
    error = function(e) {
      list(
        status = "failed",
        reason = e$message,
        results = data.frame(),
        details = list(),
        parameters = list()
      )
    }
  )
  if (!is.list(bundle) || is.null(bundle$status)) {
    bundle <- list(
      status = "failed",
      reason = "Method did not return a valid bundle.",
      results = data.frame(),
      details = list(),
      parameters = list()
    )
  }
  bundle$parameters <- utils::modifyList(
    list(method = method_name),
    bundle$parameters %||% list()
  )
  bundle$results <- csde_schema(bundle$results)
  if (nrow(bundle$results) > 0) {
    bundle$results$method <- method_name
  }
  bundle
}

prop_long <- function(prop_matrix, method_name) {
  if (
    is.null(prop_matrix) || nrow(prop_matrix) == 0 || ncol(prop_matrix) == 0
  ) {
    return(data.frame(
      sample = character(),
      cell_type = character(),
      proportion = numeric(),
      method = character(),
      stringsAsFactors = FALSE
    ))
  }
  df <- as.data.frame(as.table(prop_matrix), stringsAsFactors = FALSE)
  colnames(df) <- c("sample", "cell_type", "proportion")
  df$method <- method_name
  deconv_schema(df)
}

deconv_mat <- function(deconv_results) {
  if (is.null(deconv_results) || nrow(deconv_results) == 0) {
    return(NULL)
  }
  df <- deconv_results[, c("sample", "cell_type", "proportion"), drop = FALSE]
  as_matrix(stats::xtabs(
    proportion ~ sample + cell_type,
    data = df
  ))
}

deconv_schema <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame(
      sample = character(),
      cell_type = character(),
      proportion = numeric(),
      method = character(),
      stringsAsFactors = FALSE
    ))
  }
  required <- c("sample", "cell_type", "proportion", "method")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    for (nm in missing) {
      if (identical(nm, "proportion")) {
        df[[nm]] <- NA_real_
      } else {
        df[[nm]] <- NA_character_
      }
    }
  }
  out <- df[, required, drop = FALSE]
  rownames(out) <- NULL
  out
}

csde_schema <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame(
      gene = character(),
      cell_type = character(),
      group1 = character(),
      group2 = character(),
      effect = numeric(),
      p_val = numeric(),
      p_val_adj = numeric(),
      method = character(),
      stringsAsFactors = FALSE
    ))
  }
  required <- c(
    "gene",
    "cell_type",
    "group1",
    "group2",
    "effect",
    "p_val",
    "p_val_adj",
    "method"
  )
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    for (nm in missing) {
      if (nm %in% c("effect", "p_val", "p_val_adj")) {
        df[[nm]] <- NA_real_
      } else {
        df[[nm]] <- NA_character_
      }
    }
  }
  out <- df[, required, drop = FALSE]
  rownames(out) <- NULL
  out
}
