#' @title Run MDIC3 cell-cell communication inference
#'
#' @inheritParams thisutils::log_message
#' @param object A Seurat object or a gene-by-cell expression matrix.
#' @param group.by Metadata column used as cell labels when `object` is a Seurat
#' object.
#' @param grn Gene-by-gene GRN matrix aligned to expression genes, or a data
#' frame with columns `TF`, `target`, and `importance`.
#' @param grn_method GRN inference method used when `grn = NULL`.
#' `"grnboost2"` calls [RunGRNBoost2()]; `"gniplr"` calls [RunGNIPLR()];
#' `"correlation"` uses absolute expression correlation as a lightweight GRN
#' approximation.
#' @param grn_backend Backend passed to [RunGRNBoost2()] when
#' `grn_method = "grnboost2"`.
#' @param assay Assay used when `object` is a Seurat object.
#' @param layer Assay layer used when `object` is a Seurat object.
#' @param labels Cell labels used when `object` is a matrix.
#' @param genes_in Matrix orientation for matrix inputs. `"rows"` means genes x
#' cells; `"columns"` means cells x genes.
#' @param backend Runtime backend for the MDIC3 scoring step. Supports `"cpp"`
#' and `"python"`.
#' @param envname Python environment used when `backend = "python"`.
#' @param conda Conda-compatible executable used when `backend = "python"`.
#' @param ... Additional arguments passed to [RunGRNBoost2()] or [RunGNIPLR()]
#' when those methods are selected through `grn_method`.
#'
#' @return A Seurat object with results stored in `object@tools[["MDIC3"]]`, or
#' a list for matrix input.
#'
#' @examples
#' data(pancreas_sub)
#' expr <- GetAssayData5(
#'   pancreas_sub,
#'   assay = SeuratObject::DefaultAssay(pancreas_sub),
#'   layer = "counts"
#' )
#' expr <- as.matrix(expr[, seq_len(8)])
#' expr <- expr[
#'   names(sort(apply(expr, 1, stats::var), decreasing = TRUE))[seq_len(5)],
#' ]
#' labels <- as.character(pancreas_sub$SubCellType[seq_len(8)])
#'
#' mdic3 <- RunMDIC3(
#'   expr,
#'   labels = labels,
#'   grn_method = "gniplr",
#'   correlation_threshold = 0,
#'   lasso_degree = 1,
#'   max_lag = 1,
#'   max_edges_per_target = 2
#' )
#' mdic3$celltype_communication
#' @export
RunMDIC3 <- function(object, ...) {
  UseMethod("RunMDIC3", object)
}

#' @rdname RunMDIC3
#' @export
RunMDIC3.Seurat <- function(
  object,
  group.by,
  grn = NULL,
  grn_method = c("grnboost2", "gniplr", "correlation"),
  grn_backend = c("cpp", "python"),
  assay = NULL,
  layer = "data",
  backend = c("cpp", "python"),
  envname = "scop_env",
  conda = "auto",
  verbose = TRUE,
  ...
) {
  if (!inherits(object, "Seurat")) {
    log_message(
      "{.arg object} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (missing(group.by) || !group.by %in% colnames(object[[]])) {
    log_message(
      "{.arg group.by} must be a valid metadata column",
      message_type = "error"
    )
  }
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  expr <- GetAssayData5(object, assay = assay, layer = layer)
  res <- RunMDIC3.default(
    expr,
    labels = object[[group.by, drop = TRUE]],
    grn = grn,
    grn_method = grn_method,
    grn_backend = grn_backend,
    genes_in = "rows",
    backend = backend,
    envname = envname,
    conda = conda,
    verbose = verbose,
    ...
  )
  res$parameters$assay <- assay
  res$parameters$layer <- layer
  res$parameters$group.by <- group.by
  object@tools[["MDIC3"]] <- res
  ccc_update_unified_bundle(
    srt = object,
    method = "MDIC3",
    bundle = res,
    backend = match.arg(backend)
  )
}

#' @rdname RunMDIC3
#' @export
RunMDIC3.matrix <- function(object, ...) {
  RunMDIC3.default(object, ...)
}

#' @rdname RunMDIC3
#' @export
RunMDIC3.default <- function(
  object,
  labels,
  grn = NULL,
  grn_method = c("grnboost2", "gniplr", "correlation"),
  grn_backend = c("cpp", "python"),
  genes_in = c("rows", "columns"),
  backend = c("cpp", "python"),
  envname = "scop_env",
  conda = "auto",
  verbose = TRUE,
  ...
) {
  backend <- match.arg(backend)
  grn_method <- match.arg(grn_method)
  grn_backend <- match.arg(grn_backend)
  genes_in <- match.arg(genes_in)

  expr <- if (inherits(object, "Matrix")) {
    as.matrix(object)
  } else {
    as.matrix(object)
  }
  if (!is.numeric(expr)) {
    storage.mode(expr) <- "double"
  }
  if (identical(genes_in, "columns")) {
    expr <- t(expr)
  }
  if (is.null(rownames(expr)) || any(!nzchar(rownames(expr)))) {
    log_message(
      "{.arg object} must have gene names as row names",
      message_type = "error"
    )
  }
  if (is.null(colnames(expr)) || any(!nzchar(colnames(expr)))) {
    colnames(expr) <- paste0("cell", seq_len(ncol(expr)))
  }
  expr[!is.finite(expr)] <- 0

  labels <- as.character(labels)
  if (length(labels) != ncol(expr)) {
    log_message(
      "{.arg labels} length must match the number of cells",
      message_type = "error"
    )
  }
  if (anyNA(labels) || any(!nzchar(labels))) {
    log_message(
      "{.arg labels} must not contain missing or empty values",
      message_type = "error"
    )
  }
  group <- factor(labels, levels = unique(labels))
  group_ids <- as.integer(group)
  group_levels <- levels(group)
  genes <- rownames(expr)

  if (!is.null(grn)) {
    grn_method_used <- "external"
  } else if (identical(grn_method, "grnboost2")) {
    grn_method_used <- "grnboost2"
    grn <- RunGRNBoost2(
      expr,
      genes_in = "rows",
      backend = grn_backend,
      verbose = verbose,
      ...
    )
  } else if (identical(grn_method, "gniplr")) {
    grn_method_used <- "gniplr"
    grn_adjacency <- RunGNIPLR(
      expr,
      targets = genes,
      genes_in = "rows",
      envname = envname,
      conda = conda,
      verbose = verbose,
      ...
    )
    grn <- attr(grn_adjacency, "grn_matrix")
    if (is.null(grn)) {
      grn <- grn_adjacency
    }
  } else {
    grn_method_used <- "correlation"
    grn <- abs(suppressWarnings(stats::cor(
      t(expr),
      use = "pairwise.complete.obs"
    )))
    grn[!is.finite(grn)] <- 0
    diag(grn) <- 0
  }

  if (is.data.frame(grn)) {
    required <- c("TF", "target", "importance")
    if (!all(required %in% colnames(grn))) {
      log_message(
        "{.arg grn} data frame must contain columns {.field TF}, {.field target}, and {.field importance}",
        message_type = "error"
      )
    }
    grn_matrix <- matrix(
      0,
      nrow = length(genes),
      ncol = length(genes),
      dimnames = list(genes, genes)
    )
    from <- match(as.character(grn$TF), genes)
    to <- match(as.character(grn$target), genes)
    keep <- !is.na(from) & !is.na(to)
    grn_matrix[cbind(
      from[keep],
      to[keep]
    )] <- suppressWarnings(as.numeric(grn$importance[keep]))
  } else {
    grn_matrix <- as.matrix(grn)
    if (!is.numeric(grn_matrix)) {
      storage.mode(grn_matrix) <- "double"
    }
    if (nrow(grn_matrix) != ncol(grn_matrix)) {
      log_message(
        "{.arg grn} must be a square gene-by-gene matrix",
        message_type = "error"
      )
    }
    if (!is.null(rownames(grn_matrix)) && !is.null(colnames(grn_matrix))) {
      missing_genes <- setdiff(
        genes,
        intersect(rownames(grn_matrix), colnames(grn_matrix))
      )
      if (length(missing_genes) > 0L) {
        log_message(
          "{.arg grn} is missing {.val {length(missing_genes)}} expression genes",
          message_type = "error"
        )
      }
      grn_matrix <- grn_matrix[genes, genes, drop = FALSE]
    } else if (nrow(grn_matrix) != length(genes)) {
      log_message(
        "{.arg grn} without dimnames must have the same size as expression genes",
        message_type = "error"
      )
    } else {
      dimnames(grn_matrix) <- list(genes, genes)
    }
  }
  grn_matrix[!is.finite(grn_matrix)] <- 0

  log_message(
    "Running {.pkg MDIC3} with {.arg backend = {backend}} and {.arg grn_method = {grn_method_used}} on {.val {nrow(expr)}} genes and {.val {ncol(expr)}} cells",
    verbose = verbose
  )

  if (identical(backend, "cpp")) {
    raw <- mdic3_score_cpp(
      expression = as.matrix(expr),
      grn = grn_matrix,
      group = group_ids
    )
  } else {
    np <- reticulate::import("numpy", convert = TRUE)
    aa <- np$array(as.matrix(expr), dtype = "float64")
    grn_np <- np$array(grn_matrix, dtype = "float64")
    svd <- np$linalg$svd(aa, full_matrices = TRUE)
    singular_values <- svd[[1]]
    s_rect <- np$diag(singular_values)
    if (nrow(expr) < ncol(expr)) {
      s_rect <- np$hstack(list(
        s_rect,
        np$zeros(c(nrow(expr), ncol(expr) - nrow(expr)))
      ))
    } else if (nrow(expr) > ncol(expr)) {
      s_rect <- np$vstack(list(
        s_rect,
        np$zeros(c(nrow(expr) - ncol(expr), ncol(expr)))
      ))
    }
    m <- np$dot(np$linalg$pinv(np$dot(grn_np, s_rect)), aa)
    cellular <- matrix(0, nrow = nrow(m), ncol = ncol(m))
    pos <- m > 0 & is.finite(m)
    cellular[pos] <- m[pos]
    neg <- which(m < 0 & is.finite(m), arr.ind = TRUE)
    if (nrow(neg) > 0L) {
      cellular[cbind(neg[, 2], neg[, 1])] <- cellular[cbind(
        neg[, 2],
        neg[, 1]
      )] +
        abs(m[neg])
    }
    cellular_log <- matrix(0, nrow = nrow(cellular), ncol = ncol(cellular))
    idx <- is.finite(cellular) & cellular > 0
    logged <- log10(cellular[idx])
    cellular_log[idx] <- ifelse(logged > 0, logged, 0)
    type_raw <- matrix(
      0,
      nrow = length(group_levels),
      ncol = length(group_levels)
    )
    for (i in seq_len(nrow(cellular))) {
      gi <- group_ids[i]
      for (j in seq_len(ncol(cellular))) {
        value <- cellular[i, j]
        if (!is.finite(value) || value == 0) {
          next
        }
        gj <- group_ids[j]
        if (gi == gj) {
          type_raw[gi, gj] <- type_raw[gi, gj] + abs(value)
        } else if (value > 0) {
          type_raw[gi, gj] <- type_raw[gi, gj] + value
        } else {
          type_raw[gj, gi] <- type_raw[gj, gi] + abs(value)
        }
      }
    }
    type_log <- matrix(0, nrow = nrow(type_raw), ncol = ncol(type_raw))
    idx <- is.finite(type_raw) & type_raw > 0
    logged <- log10(type_raw[idx])
    type_log[idx] <- ifelse(logged > 0, logged, 0)
    raw <- list(
      mdic3_matrix = m,
      cellular_communication = cellular,
      cellular_communication_log = cellular_log,
      celltype_communication_raw = type_raw,
      celltype_communication = type_log
    )
  }

  raw <- lapply(raw, as.matrix)
  dimnames(raw$mdic3_matrix) <- list(colnames(expr), colnames(expr))
  dimnames(raw$cellular_communication) <- list(colnames(expr), colnames(expr))
  dimnames(raw$cellular_communication_log) <- list(
    colnames(expr),
    colnames(expr)
  )
  dimnames(raw$celltype_communication_raw) <- list(group_levels, group_levels)
  dimnames(raw$celltype_communication) <- list(group_levels, group_levels)

  long_idx <- which(is.finite(raw$celltype_communication), arr.ind = TRUE)
  long_table <- data.frame(
    sender = rownames(raw$celltype_communication)[long_idx[, 1]],
    receiver = colnames(raw$celltype_communication)[long_idx[, 2]],
    ligand = "",
    receptor = "",
    interaction_name = "MDIC3",
    interaction_label = "MDIC3",
    pair_lr = "MDIC3",
    pathway_name = "MDIC3",
    score = as.numeric(raw$celltype_communication[long_idx]),
    pvalue = 1,
    method = "MDIC3",
    stringsAsFactors = FALSE
  )
  c(
    raw,
    list(
      grn = grn_matrix,
      labels = data.frame(
        cell = colnames(expr),
        label = labels,
        stringsAsFactors = FALSE
      ),
      long_table = long_table,
      pair_table = aggregate_ccc_long(long_table, backend = "cpp"),
      parameters = list(
        method = "MDIC3",
        backend = backend,
        grn_method = grn_method_used,
        grn_backend = if (identical(grn_method_used, "grnboost2")) {
          grn_backend
        } else if (identical(grn_method_used, "gniplr")) {
          "python"
        } else {
          NA_character_
        }
      )
    )
  )
}
