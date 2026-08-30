#' Aggregate expression by cell group
#'
#' Seurat already performs this operation as sparse matrix multiplication.
#' scop preserves its complete parameter and object contract.
#'
#' @inheritParams Seurat::AggregateExpression
#' @return A named list of aggregated matrices or a Seurat object.
#' @export
AggregateExpression <- function(
  object,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "ident",
  add.ident = NULL,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  margin = 1,
  verbose = TRUE,
  ...
) {
  Seurat::AggregateExpression(
    object = object,
    assays = assays,
    features = features,
    return.seurat = return.seurat,
    group.by = group.by,
    add.ident = add.ident,
    normalization.method = normalization.method,
    scale.factor = scale.factor,
    margin = margin,
    verbose = verbose,
    ...
  )
}

#' Average expression by cell group
#'
#' Seurat already performs this operation as sparse matrix multiplication.
#' scop preserves its complete parameter and object contract.
#'
#' @inheritParams Seurat::AverageExpression
#' @return A named list of averaged matrices or a Seurat object.
#' @export
AverageExpression <- function(
  object,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "ident",
  add.ident = NULL,
  layer = "data",
  slot = lifecycle::deprecated(),
  verbose = TRUE,
  ...
) {
  args <- list(
    object = object,
    assays = assays,
    features = features,
    return.seurat = return.seurat,
    group.by = group.by,
    add.ident = add.ident,
    layer = layer,
    verbose = verbose,
    ...
  )
  if (lifecycle::is_present(slot)) {
    args$slot <- slot
  }
  do.call(Seurat::AverageExpression, args)
}
