#' @title Run scPagwas/scPaGWAS
#'
#' @description
#' Run the optional `scPagwas` package from `scop` without bundling LD,
#' pathway, or block-annotation resources. The wrapper validates required GWAS
#' columns, normalizes output paths, and records provenance in Seurat tools or
#' a result attribute.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt Optional Seurat object used as single-cell input.
#' @param single_data Optional Seurat object or path to an `.rds` file used by
#' `scPagwas`.
#' @param gwas_data GWAS summary statistics as a data frame.
#' @param celltype_meta Optional Seurat metadata column used to set identities.
#' @param block_annotation Genome build for bundled upstream annotations
#' (`"hg38"` or `"hg37"`) or a custom annotation path.
#' @param output.dirs Output directory passed to `scPagwas`.
#' @param cleanup_soar Whether to remove SOAR objects after the run when the
#' optional SOAR package is available.
#' @param return_seurat Whether to return a Seurat object when one is available.
#' @param ... Additional arguments passed to the upstream `scPagwas` function
#' after filtering by its formal arguments.
#'
#' @return A Seurat object or upstream result list.
#' @export
RunscPagwas <- function(
  srt = NULL,
  single_data = NULL,
  gwas_data,
  celltype_meta = NULL,
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
  scpagwas_check_gwas(gwas_data)
  single_data <- scpagwas_resolve_single_data(srt = srt, single_data = single_data)
  output.dirs <- scpagwas_abs_path(output.dirs)
  single_data <- scpagwas_prepare_seurat(single_data, celltype_meta)
  scpagwas_check_r(verbose = verbose)
  block_annotation <- scpagwas_resolve_block_annotation(block_annotation)

  fun <- scpagwas_find_runner()
  args <- c(
    list(
      Single_data = single_data,
      single_data = single_data,
      gwas_data = gwas_data,
      block_annotation = block_annotation,
      output.dirs = output.dirs
    ),
    list(...)
  )
  res <- do.call(fun, scpagwas_filter_args(fun, args))
  if (isTRUE(cleanup_soar)) {
    scpagwas_cleanup_soar()
  }

  meta <- list(
    method = "scPagwas",
    parameters = list(
      celltype_meta = celltype_meta,
      block_annotation = block_annotation,
      output.dirs = output.dirs,
      cleanup_soar = cleanup_soar
    )
  )
  if (isTRUE(return_seurat)) {
    out <- if (inherits(res, "Seurat")) res else single_data
    if (!inherits(out, "Seurat")) {
      log_message("{.arg return_seurat = TRUE} requires Seurat input or Seurat output", message_type = "error")
    }
    out@tools$scPagwas <- c(meta, list(result = res))
    return(out)
  }
  attr(res, "scPagwas") <- meta
  res
}

#' @rdname RunscPagwas
#' @export
RunscPaGWAS <- RunscPagwas

scpagwas_check_r <- function(verbose = TRUE) {
  check_r("scPagwas", verbose = verbose)
  invisible(TRUE)
}

scpagwas_find_runner <- function() {
  candidates <- c("scPagwas_main", "scPagwas")
  for (fun in candidates) {
    runner <- tryCatch(get_namespace_fun("scPagwas", fun), error = function(e) NULL)
    if (is.function(runner)) {
      return(runner)
    }
  }
  log_message("Could not find an upstream {.pkg scPagwas} runner", message_type = "error")
}

scpagwas_filter_args <- function(fun, args) {
  fun_formals <- names(formals(fun))
  if ("..." %in% fun_formals) {
    return(args)
  }
  args[intersect(names(args), fun_formals)]
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
  }
  single_data
}

scpagwas_prepare_seurat <- function(single_data, celltype_meta = NULL) {
  if (is.null(celltype_meta) || !inherits(single_data, "Seurat")) {
    return(single_data)
  }
  if (!celltype_meta %in% colnames(single_data[[]])) {
    log_message("Missing metadata column: {.val {celltype_meta}}", message_type = "error")
  }
  SeuratObject::Idents(single_data) <- single_data[[celltype_meta]][, 1]
  single_data
}

scpagwas_check_gwas <- function(gwas_data) {
  if (!inherits(gwas_data, "data.frame")) {
    log_message("{.arg gwas_data} must be a data frame", message_type = "error")
  }
  required <- c("chrom", "pos", "rsid", "se", "beta", "maf")
  missing_cols <- setdiff(required, colnames(gwas_data))
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
  env <- new.env(parent = emptyenv())
  utils::data(list = name, package = "scPagwas", envir = env)
  if (!exists(name, envir = env, inherits = FALSE)) {
    log_message("Could not load {.pkg scPagwas} data object {.val {name}}", message_type = "error")
  }
  scpagwas_check_block_annotation(get(name, envir = env, inherits = FALSE))
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

scpagwas_cleanup_soar <- function() {
  if (!requireNamespace("SOAR", quietly = TRUE)) {
    return(invisible(FALSE))
  }
  ns <- asNamespace("SOAR")

  rm_fun <- get("RemoveAllObjects", envir = ns, inherits = FALSE)
  if (is.function(rm_fun)) {
    rm_fun()
    return(invisible(TRUE))
  }

  objects_fun <- get("Objects", envir = ns, inherits = FALSE)
  remove_fun <- get("Remove", envir = ns, inherits = FALSE)
  if (!is.function(objects_fun) || !is.function(remove_fun)) {
    return(invisible(TRUE))
  }

  object_ids <- tryCatch(objects_fun(), error = function(e) NULL)
  if (!is.character(object_ids) || length(object_ids) < 1L) {
    return(invisible(TRUE))
  }

  tryCatch(remove_fun(object_ids), error = function(e) NULL)
  invisible(TRUE)
}
