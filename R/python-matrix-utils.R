python_cells_by_features <- function(x) {
  if (inherits(x, "sparseMatrix")) {
    return(methods::as(Matrix::t(x), "CsparseMatrix"))
  }
  Matrix::t(as.matrix(x))
}

python_anndata_from_matrix <- function(
  x,
  embeddings = list(),
  anndata_module = NULL
) {
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    log_message(
      "Expression matrix must have feature and cell names",
      message_type = "error"
    )
  }
  anndata_module <- anndata_module %||% reticulate::import("anndata")
  adata <- anndata_module$AnnData(X = python_cells_by_features(x))
  adata$obs_names <- as.character(colnames(x))
  adata$var_names <- as.character(rownames(x))
  for (name in names(embeddings)) {
    embedding <- as.matrix(embeddings[[name]])
    if (!is.null(rownames(embedding))) {
      missing_cells <- setdiff(colnames(x), rownames(embedding))
      if (length(missing_cells) > 0) {
        log_message(
          "Embedding {.val {name}} is missing cells",
          message_type = "error"
        )
      }
      embedding <- embedding[colnames(x), , drop = FALSE]
    } else if (nrow(embedding) != ncol(x)) {
      log_message(
        "Embedding {.val {name}} has incompatible dimensions",
        message_type = "error"
      )
    }
    adata$obsm[[name]] <- embedding
  }
  adata
}
