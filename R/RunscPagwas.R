#' @title Run scPagwas
#'
#' @description
#' Run the optional `scPagwas` package from `scop` without bundling LD,
#' pathway, or block-annotation resources. The wrapper validates required GWAS
#' columns, normalizes output paths, and records provenance in Seurat tools or
#' a result attribute. It prefers the upstream `scPagwas_main2` runner and
#' applies Seurat 5 compatibility to local function copies without modifying
#' the installed backend namespace.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt Optional Seurat object used as single-cell input.
#' @param single_data Optional Seurat object or path to a Seurat `.rds` file.
#' @param gwas_data GWAS summary statistics as a data frame or delimited text
#' file. Required columns are `chrom`, `pos`, `rsid`, `se`, `beta`, and `maf`.
#' @param celltype_meta Optional Seurat metadata column used to set identities.
#' @param singlecell Whether to calculate single-cell results.
#' @param celltype Whether to calculate cell-type results.
#' @param assay Assay used by `scPagwas`. Defaults to the active assay for a
#' Seurat object and to `"RNA"` for other inputs.
#' @param block_annotation Genome build for bundled upstream annotations
#' (`"hg38"` or `"hg37"`) or a custom annotation path.
#' @param output.dirs Output directory passed to `scPagwas`.
#' @param cleanup_soar Deprecated compatibility argument. SOAR cleanup is
#' managed by the upstream `scPagwas` backend and is ignored by `scop`.
#' @param return_seurat Whether to return a Seurat object when one is available.
#' @param ... Additional arguments passed to the upstream `scPagwas` function
#' after filtering by its formal arguments.
#'
#' @return A Seurat object or upstream result list.
#' @seealso [PlotScPagwas()]
#' @export
RunscPagwas <- function(
  srt = NULL,
  single_data = NULL,
  gwas_data,
  celltype_meta = NULL,
  singlecell = TRUE,
  celltype = TRUE,
  assay = NULL,
  block_annotation = c("hg38", "hg37", "custom"),
  output.dirs = tempdir(),
  cleanup_soar = TRUE,
  return_seurat = !is.null(srt) || inherits(single_data, "Seurat"),
  verbose = TRUE,
  ...
) {
  block_annotation <- scpagwas_validate_block_annotation_selector(block_annotation)
  if (missing(gwas_data) || is.null(gwas_data)) {
    log_message("{.arg gwas_data} is required", message_type = "error")
  }
  gwas_data <- scpagwas_resolve_gwas_data(gwas_data)
  single_data <- scpagwas_resolve_single_data(srt = srt, single_data = single_data)
  if (is.null(assay)) {
    assay <- if (inherits(single_data, "Seurat")) {
      SeuratObject::DefaultAssay(single_data)
    } else {
      "RNA"
    }
  }
  output.dirs <- scpagwas_abs_path(output.dirs)
  dir.create(output.dirs, recursive = TRUE, showWarnings = FALSE)
  single_data <- scpagwas_prepare_seurat(single_data, celltype_meta, assay)
  scpagwas_check_r(verbose = verbose)
  block_annotation <- scpagwas_resolve_block_annotation(block_annotation)

  fun <- scpagwas_prepare_runner(scpagwas_find_runner())
  extra_args <- list(...)
  extra_args <- scpagwas_add_default_data_args(fun, extra_args)
  args <- list(
    Single_data = single_data,
    single_data = single_data,
    gwas_data = gwas_data,
    assay = assay,
    singlecell = singlecell,
    celltype = celltype,
    block_annotation = block_annotation,
    output.dirs = output.dirs,
    seurat_return = return_seurat
  )
  extra_args <- extra_args[setdiff(names(extra_args), names(args))]
  args <- c(args, extra_args)
  output_context <- scpagwas_output_context(output.dirs)
  args$output.dirs <- output_context$backend_dir
  on.exit(output_context$restore(), add = TRUE)
  old_r_local_cache <- Sys.getenv("R_LOCAL_CACHE", unset = NA_character_)
  on.exit(
    if (is.na(old_r_local_cache)) {
      Sys.unsetenv("R_LOCAL_CACHE")
    } else {
      Sys.setenv(R_LOCAL_CACHE = old_r_local_cache)
    },
    add = TRUE
  )
  res <- do.call(fun, scpagwas_filter_args(fun, args))

  meta <- list(
    method = "scPagwas",
    backend_runner = attr(fun, "scpagwas_runner", exact = TRUE),
    parameters = list(
      celltype_meta = celltype_meta,
      singlecell = singlecell,
      celltype = celltype,
      assay = assay,
      block_annotation = block_annotation,
      output.dirs = output.dirs,
      cleanup_soar = cleanup_soar
    )
  )
  if (isTRUE(return_seurat)) {
    if (!inherits(res, "Seurat")) {
      log_message(
        "The {.pkg scPagwas} backend did not return a Seurat object despite {.arg return_seurat = TRUE}",
        message_type = "error"
      )
    }
    res@tools$scPagwas <- meta
    return(res)
  }
  attr(res, "scPagwas") <- meta
  res
}

#' @title Deprecated scPagwas Alias
#'
#' @description
#' `RunscPaGWAS()` is retained for compatibility. Use [RunscPagwas()] for new
#' code.
#'
#' @md
#' @param ... Arguments passed to [RunscPagwas()].
#'
#' @return The value returned by [RunscPagwas()].
#' @name RunscPaGWAS-deprecated
#' @keywords internal
#' @export
RunscPaGWAS <- function(...) {
  .Deprecated("RunscPagwas")
  RunscPagwas(...)
}

#' @title Plot scPagwas Scores
#'
#' @description
#' Plot the score and adjusted p-value metadata produced by [RunscPagwas()] on
#' an existing Seurat reduction. The plots are returned and can optionally be
#' saved as PDF files.
#'
#' @md
#' @param srt A Seurat object returned by [RunscPagwas()].
#' @param reduction Reduction used for plotting, either `"umap"` or `"tsne"`.
#' @param features Numeric scPagwas metadata columns to plot. By default,
#' available gPAS, TRS, and down-TRS score columns are used.
#' @param p_threshold Adjusted p-value threshold used to identify significant
#' cells. Set to `NULL` to omit the significance plot.
#' @param output.dir Optional directory in which to save PDF files.
#' @param width,height PDF dimensions in inches.
#' @param point_size Point size passed to Seurat plotting functions.
#' @param low,high Colors for low and high score values.
#' @param do_plot Whether to print each plot.
#'
#' @return A named list of ggplot objects.
#' @export
PlotScPagwas <- function(
  srt,
  reduction = c("umap", "tsne"),
  features = NULL,
  p_threshold = 0.05,
  output.dir = NULL,
  width = 7,
  height = 7,
  point_size = NULL,
  low = "#000957",
  high = "#EBE645",
  do_plot = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message("{.arg srt} must be a Seurat object", message_type = "error")
  }
  reduction <- match.arg(reduction)
  if (!reduction %in% names(srt@reductions)) {
    log_message(
      "Reduction {.val {reduction}} is not present in {.arg srt}",
      message_type = "error"
    )
  }
  metadata <- srt[[]]
  if (is.null(features)) {
    features <- intersect(
      c(
        "scPagwas.gPAS.score",
        "scPagwas.TRS.Score1",
        "scPagwas.downTRS.Score2"
      ),
      colnames(metadata)
    )
  }
  missing_features <- setdiff(features, colnames(metadata))
  if (length(missing_features) > 0L) {
    log_message(
      "Missing scPagwas metadata column{?s}: {.val {missing_features}}",
      message_type = "error"
    )
  }
  if (length(features) < 1L) {
    log_message("No scPagwas score metadata was found", message_type = "error")
  }
  non_numeric <- features[!vapply(metadata[features], is.numeric, logical(1))]
  if (length(non_numeric) > 0L) {
    log_message(
      "scPagwas score column{?s} must be numeric: {.val {non_numeric}}",
      message_type = "error"
    )
  }

  plots <- lapply(features, function(feature) {
    Seurat::FeaturePlot(
      srt,
      features = feature,
      reduction = reduction,
      cols = c(low, high),
      pt.size = point_size,
      combine = FALSE
    )[[1]]
  })
  names(plots) <- features

  if (!is.null(p_threshold)) {
    if (length(p_threshold) != 1L || !is.finite(p_threshold) ||
      p_threshold <= 0 || p_threshold >= 1) {
      log_message(
        "{.arg p_threshold} must be a finite number between 0 and 1",
        message_type = "error"
      )
    }
    p_column <- "Random_Correct_BG_adjp"
    if (!p_column %in% colnames(metadata)) {
      log_message(
        "Missing scPagwas metadata column: {.val {p_column}}. Set {.arg p_threshold = NULL} to omit the significance plot.",
        message_type = "error"
      )
    }
    plot_srt <- srt
    plot_srt$.scPagwas_significant <- factor(
      ifelse(
        metadata[[p_column]] <= p_threshold,
        "Significant",
        "Not significant"
      ),
      levels = c("Not significant", "Significant")
    )
    plots$significant_cells <- Seurat::DimPlot(
      plot_srt,
      reduction = reduction,
      group.by = ".scPagwas_significant",
      cols = c("#E4DCCF", "#EA5455"),
      pt.size = point_size,
      combine = FALSE
    )[[1]] +
      ggplot2::ggtitle(
        paste0("Random_Correct_BG_adjp <= ", p_threshold)
      )
  }

  if (!is.null(output.dir)) {
    output.dir <- scpagwas_abs_path(output.dir)
    dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
    for (name in names(plots)) {
      ggplot2::ggsave(
        filename = file.path(
          output.dir,
          paste0(make.names(name), "_", reduction, ".pdf")
        ),
        plot = plots[[name]],
        width = width,
        height = height
      )
    }
  }
  if (isTRUE(do_plot)) {
    invisible(lapply(plots, print))
  }
  plots
}

scpagwas_check_r <- function(verbose = TRUE) {
  check_r("sulab-wmu/scPagwas", dependencies = NA, verbose = verbose)
  if (!is.function(scpagwas_find_runner(error = FALSE))) {
    log_message(
      "Failed to install or load {.pkg scPagwas}. Install it manually with {.code pak::pkg_install('sulab-wmu/scPagwas')}",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

scpagwas_get_fun <- function(fun, error = TRUE) {
  out <- tryCatch(
    suppressWarnings(get_namespace_fun("scPagwas", fun)),
    error = function(e) NULL
  )
  if (!is.function(out) && isTRUE(error)) {
    log_message(
      "Could not find {.pkg scPagwas} function {.val {fun}}",
      message_type = "error"
    )
  }
  out
}

scpagwas_find_runner <- function(error = TRUE) {
  candidates <- c("scPagwas_main2", "scPagwas_main", "scPagwas")
  for (fun in candidates) {
    runner <- scpagwas_get_fun(fun, error = FALSE)
    if (is.function(runner)) {
      attr(runner, "scpagwas_runner") <- fun
      return(runner)
    }
  }
  if (isTRUE(error)) {
    log_message("Could not find an upstream {.pkg scPagwas} runner", message_type = "error")
  }
  NULL
}

scpagwas_prepare_runner <- function(fun) {
  runner_name <- attr(fun, "scpagwas_runner", exact = TRUE)
  compat_env <- new.env(parent = environment(fun))
  for (name in c("Single_data_input", "Get_CorrectBg_p")) {
    helper <- scpagwas_get_fun(name, error = FALSE)
    if (is.function(helper)) {
      helper <- scpagwas_rewrite_seurat_calls(helper)
      environment(helper) <- compat_env
      assign(name, helper, envir = compat_env)
    }
  }
  fun <- scpagwas_rewrite_seurat_calls(fun)
  environment(fun) <- compat_env
  attr(fun, "scpagwas_runner") <- runner_name
  fun
}

scpagwas_rewrite_seurat_calls <- function(fun) {
  rewrite <- function(expr) {
    if (!is.call(expr)) {
      return(expr)
    }
    expr <- as.call(lapply(as.list(expr), rewrite))
    call_name <- scpagwas_call_name(expr)
    arg_names <- names(expr)
    if (call_name %in% c(
      "GetAssayData",
      "Seurat::GetAssayData",
      "SeuratObject::GetAssayData"
    )) {
      names(expr)[which(arg_names == "slot")] <- "layer"
    }
    if (identical(call_name, "Seurat::AddModuleScore")) {
      if (is.null(arg_names)) {
        arg_names <- rep("", length(expr))
      }
      feature_arg <- which(
        !nzchar(arg_names) &
          vapply(
            as.list(expr),
            identical,
            logical(1),
            y = as.name("scPagwas_topgenes")
          )
      )
      if (length(feature_arg) == 1L) {
        expr[[feature_arg]] <- call("list", expr[[feature_arg]])
        arg_names[[feature_arg]] <- "features"
        names(expr) <- arg_names
      }
    }
    expr
  }
  body(fun) <- rewrite(body(fun))
  fun
}

scpagwas_call_name <- function(expr) {
  head <- expr[[1]]
  if (is.symbol(head)) {
    return(as.character(head))
  }
  if (is.call(head) && identical(head[[1]], as.name("::"))) {
    return(paste0(as.character(head[[2]]), "::", as.character(head[[3]])))
  }
  ""
}

scpagwas_filter_args <- function(fun, args) {
  fun_formals <- names(formals(fun))
  if ("..." %in% fun_formals) {
    return(args)
  }
  args[intersect(names(args), fun_formals)]
}

scpagwas_add_default_data_args <- function(fun, args) {
  fun_formals <- names(formals(fun))
  if ("Pathway_list" %in% fun_formals && is.null(args$Pathway_list)) {
    args$Pathway_list <- scpagwas_package_data_raw("Genes_by_pathway_kegg")
  }
  if ("chrom_ld" %in% fun_formals && is.null(args$chrom_ld)) {
    args$chrom_ld <- scpagwas_package_data_raw("chrom_ld")
  }
  args
}

scpagwas_output_context <- function(output.dirs) {
  old_wd <- getwd()
  # scPagwas prefixes output.dirs with "./", so call it from the parent and
  # pass only the directory name to keep Windows drive paths valid.
  setwd(dirname(output.dirs))
  list(
    backend_dir = basename(output.dirs),
    restore = function() {
      setwd(old_wd)
      invisible(TRUE)
    }
  )
}

scpagwas_resolve_single_data <- function(srt = NULL, single_data = NULL) {
  if (!is.null(srt) && !is.null(single_data)) {
    log_message("Provide only one of {.arg srt} or {.arg single_data}", message_type = "error")
  }
  single_data <- single_data %||% srt
  if (is.null(single_data)) {
    log_message("Provide {.arg srt} or {.arg single_data}", message_type = "error")
  }
  if (is.character(single_data)) {
    if (length(single_data) != 1L || is.na(single_data) || !nzchar(single_data)) {
      log_message("{.arg single_data} path must be a single non-empty string", message_type = "error")
    }
    single_data <- scpagwas_abs_path(single_data)
    if (!grepl("\\.rds$", single_data, ignore.case = TRUE)) {
      log_message("{.arg single_data} path must point to an {.file .rds} file", message_type = "error")
    }
    if (!file.exists(single_data)) {
      log_message(
        "{.arg single_data} file does not exist: {.file {single_data}}",
        message_type = "error"
      )
    }
    single_data <- readRDS(single_data)
    if (!inherits(single_data, "Seurat")) {
      log_message(
        "{.arg single_data} must contain a Seurat object",
        message_type = "error"
      )
    }
  }
  single_data
}

scpagwas_prepare_seurat <- function(single_data, celltype_meta = NULL, assay = NULL) {
  if (!inherits(single_data, "Seurat")) {
    return(single_data)
  }
  assay <- assay %||% SeuratObject::DefaultAssay(single_data)
  if (!assay %in% names(single_data@assays)) {
    log_message("Missing assay: {.val {assay}}", message_type = "error")
  }
  data_layer <- suppressWarnings(
    tryCatch(
      GetAssayData5(single_data, assay = assay, layer = "data"),
      error = function(e) NULL
    )
  )
  if (is.null(data_layer) || any(dim(data_layer) == 0L)) {
    single_data <- Seurat::NormalizeData(single_data, assay = assay, verbose = FALSE)
  }
  if (is.null(celltype_meta)) {
    return(single_data)
  }
  if (!celltype_meta %in% colnames(single_data[[]])) {
    log_message("Missing metadata column: {.val {celltype_meta}}", message_type = "error")
  }
  SeuratObject::Idents(single_data) <- single_data[[celltype_meta]][, 1]
  single_data
}

scpagwas_resolve_gwas_data <- function(gwas_data) {
  if (inherits(gwas_data, "data.frame")) {
    scpagwas_check_gwas(gwas_data)
    return(gwas_data)
  }
  if (!is.character(gwas_data) || length(gwas_data) != 1L ||
    is.na(gwas_data) || !nzchar(gwas_data)) {
    log_message(
      "{.arg gwas_data} must be a data frame or a single file path",
      message_type = "error"
    )
  }
  gwas_data <- scpagwas_abs_path(gwas_data)
  if (!file.exists(gwas_data)) {
    log_message(
      "{.arg gwas_data} file does not exist: {.file {gwas_data}}",
      message_type = "error"
    )
  }
  columns <- scpagwas_gwas_file_columns(gwas_data)
  scpagwas_check_gwas_columns(columns)
  normalizePath(gwas_data, mustWork = TRUE, winslash = "/")
}

scpagwas_gwas_file_columns <- function(path) {
  plain_path <- sub("\\.gz$", "", path, ignore.case = TRUE)
  ext <- tolower(tools::file_ext(plain_path))
  header <- tryCatch(
    {
      if (identical(ext, "csv")) {
        utils::read.csv(
          path,
          nrows = 0L,
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      } else {
        utils::read.table(
          path,
          header = TRUE,
          nrows = 0L,
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
      }
    },
    error = function(e) {
      log_message(
        "Could not read the {.arg gwas_data} header from {.file {path}}: {conditionMessage(e)}",
        message_type = "error"
      )
    }
  )
  colnames(header)
}

scpagwas_check_gwas <- function(gwas_data) {
  if (!inherits(gwas_data, "data.frame")) {
    log_message("{.arg gwas_data} must be a data frame", message_type = "error")
  }
  scpagwas_check_gwas_columns(colnames(gwas_data))
  invisible(TRUE)
}

scpagwas_check_gwas_columns <- function(columns) {
  required <- c("chrom", "pos", "rsid", "se", "beta", "maf")
  missing_cols <- setdiff(required, columns)
  if (length(missing_cols) > 0L) {
    log_message(
      "{.arg gwas_data} is missing required column{?s}: {.val {missing_cols}}",
      message_type = "error"
    )
  }
  invisible(TRUE)
}

scpagwas_validate_block_annotation_selector <- function(block_annotation) {
  if (length(block_annotation) > 1L) {
    block_annotation <- block_annotation[[1]]
  }
  if (length(block_annotation) != 1L || is.na(block_annotation) || !nzchar(block_annotation)) {
    log_message("{.arg block_annotation} must be {.val hg38}, {.val hg37}, or a custom path", message_type = "error")
  }
  if (identical(block_annotation, "custom")) {
    log_message("{.arg block_annotation = 'custom'} requires a custom annotation path", message_type = "error")
  }
  block_annotation
}

scpagwas_resolve_block_annotation <- function(block_annotation) {
  if (identical(block_annotation, "hg38")) {
    return(scpagwas_package_data("block_annotation"))
  }
  if (identical(block_annotation, "hg37")) {
    return(scpagwas_package_data("block_annotation_hg37"))
  }
  scpagwas_read_block_annotation(block_annotation)
}

scpagwas_package_data <- function(name) {
  scpagwas_check_block_annotation(scpagwas_package_data_raw(name))
}

scpagwas_package_data_raw <- function(name) {
  env <- new.env(parent = emptyenv())
  utils::data(list = name, package = "scPagwas", envir = env)
  if (!exists(name, envir = env, inherits = FALSE)) {
    log_message("Could not load {.pkg scPagwas} data object {.val {name}}", message_type = "error")
  }
  get(name, envir = env, inherits = FALSE)
}

scpagwas_read_block_annotation <- function(path) {
  path <- scpagwas_abs_path(path)
  if (!file.exists(path)) {
    log_message("{.arg block_annotation} file does not exist: {.file {path}}", message_type = "error")
  }
  ext <- tolower(tools::file_ext(path))
  block_annotation <- switch(ext,
    rds = readRDS(path),
    rda = scpagwas_load_first_object(path),
    rdata = scpagwas_load_first_object(path),
    csv = utils::read.csv(path, stringsAsFactors = FALSE),
    tsv = utils::read.delim(path, stringsAsFactors = FALSE),
    txt = utils::read.delim(path, stringsAsFactors = FALSE),
    log_message(
      "{.arg block_annotation} custom files must be {.file .rds}, {.file .RData}, {.file .csv}, or tab-delimited text",
      message_type = "error"
    )
  )
  scpagwas_check_block_annotation(block_annotation)
}

scpagwas_load_first_object <- function(path) {
  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  if (length(loaded) < 1L) {
    log_message("{.arg block_annotation} file did not contain an R object", message_type = "error")
  }
  get(loaded[[1]], envir = env, inherits = FALSE)
}

scpagwas_check_block_annotation <- function(block_annotation) {
  if (!inherits(block_annotation, "data.frame")) {
    log_message("{.arg block_annotation} must resolve to a data frame", message_type = "error")
  }
  required <- c("chrom", "start", "end", "label")
  missing_cols <- setdiff(required, colnames(block_annotation))
  if (length(missing_cols) > 0L) {
    log_message(
      "{.arg block_annotation} is missing required column{?s}: {.val {missing_cols}}",
      message_type = "error"
    )
  }
  block_annotation
}

scpagwas_abs_path <- function(path) {
  path <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[\\/])", path)) {
    path <- file.path(getwd(), path)
  }
  normalizePath(path, mustWork = FALSE, winslash = "/")
}
