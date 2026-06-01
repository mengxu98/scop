#' @title Infer gene regulatory networks with GENIE3
#'
#' @description
#' Run GENIE3 regulatory network inference and return a
#' standardized adjacency table with columns `TF`, `target`, and `importance`.
#'
#' @param object A Seurat object or expression matrix.
#' @param assay Assay used when `object` is a Seurat object.
#' @param layer Assay layer used when `object` is a Seurat object.
#' @param regulators Candidate transcription factor genes.
#' @param targets Optional target genes. If `NULL`, all genes are considered.
#' @param genes_in Matrix orientation for matrix inputs. `"rows"` means genes x
#' cells; `"columns"` means cells x genes.
#' @param max_edges_per_target Maximum incoming regulator edges retained per
#' target. The default `Inf` keeps all positive-importance links.
#' @param output_file Optional path where the adjacency table is written.
#' @param cores Number of workers used by GENIE3.
#' @param force Whether to rebuild existing `output_file`.
#' @param verbose Whether to print progress messages.
#' @param ... Additional backend-specific arguments.
#'
#' @return A data frame with columns `TF`, `target`, and `importance`.
#' @export
RunGENIE3 <- function(object, ...) {
  UseMethod("RunGENIE3", object)
}

#' @rdname RunGENIE3
#' @export
RunGENIE3.Seurat <- function(
  object,
  assay = NULL,
  layer = "counts",
  regulators = NULL,
  targets = NULL,
  max_edges_per_target = Inf,
  output_file = NULL,
  cores = 1,
  force = FALSE,
  verbose = TRUE,
  ...
) {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  grn_matrix <- GetAssayData5(object, assay = assay, layer = layer)
  RunGENIE3.default(
    grn_matrix,
    regulators = regulators,
    targets = targets,
    genes_in = "rows",
    max_edges_per_target = max_edges_per_target,
    output_file = output_file,
    cores = cores,
    force = force,
    verbose = verbose,
    ...
  )
}

#' @rdname RunGENIE3
#' @export
RunGENIE3.matrix <- function(object, ...) {
  RunGENIE3.default(object, ...)
}

#' @rdname RunGENIE3
#' @export
RunGENIE3.default <- function(
  object,
  regulators = NULL,
  targets = NULL,
  genes_in = c("rows", "columns"),
  max_edges_per_target = Inf,
  output_file = NULL,
  cores = 1,
  force = FALSE,
  verbose = TRUE,
  ...
) {
  genes_in <- match.arg(genes_in)
  if (identical(genes_in, "columns")) {
    object <- Matrix::t(object)
  }
  check_r("GENIE3", verbose = FALSE)
  grn_matrix <- object
  inputs <- scenic_normalize_grn_inputs(
    grn_matrix,
    regulators = regulators,
    targets = targets
  )
  expr <- Matrix::t(grn_matrix)
  weight_matrix <- GENIE3::GENIE3(
    exprMatrix = as.matrix(expr),
    regulators = inputs[["regulators"]],
    targets = inputs[["targets"]],
    nCores = as.integer(max(1L, cores)),
    ...
  )
  adjacency <- GENIE3::getLinkList(weightMatrix = weight_matrix)
  colnames(adjacency)[1:3] <- c("TF", "target", "importance")
  adjacency <- adjacency[, c("TF", "target", "importance"), drop = FALSE]
  adjacency <- scenic_cap_edges_per_target(
    adjacency,
    max_edges_per_target = max_edges_per_target
  )
  if (nrow(adjacency) == 0L) {
    log_message("GENIE3 returned no edges", message_type = "error")
  }
  scenic_write_grn_adjacency(
    adjacency,
    output_file = output_file,
    force = force
  )
}
