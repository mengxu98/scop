mapQuery <- function(
    exp_query,
    metadata_query,
    ref_obj,
    vars = NULL,
    sigma = 0.1,
    verbose = TRUE) {
  log_message(
    "Scaling and synchronizing query gene expression",
    verbose = verbose
  )
  idx_shared_genes <- which(ref_obj$vargenes$symbol %in% rownames(exp_query))
  shared_genes <- ref_obj$vargenes$symbol[idx_shared_genes]
  log_message(
    "Found {.val {length(shared_genes)}} reference variable genes in query dataset",
    verbose = verbose
  )
  exp_query_scaled <- symphony::scaleDataWithStats(
    exp_query[shared_genes, ],
    ref_obj$vargenes$mean[idx_shared_genes],
    ref_obj$vargenes$stddev[idx_shared_genes],
    1
  )
  exp_query_scaled_sync <- matrix(
    0,
    nrow = length(ref_obj$vargenes$symbol),
    ncol = ncol(exp_query)
  )
  exp_query_scaled_sync[idx_shared_genes, ] <- exp_query_scaled
  rownames(exp_query_scaled_sync) <- ref_obj$vargenes$symbol
  colnames(exp_query_scaled_sync) <- colnames(exp_query)
  log_message(
    "Project query cells using reference gene loadings",
    verbose = verbose
  )
  Z_pca_query <- Matrix::t(ref_obj$loadings) %*% exp_query_scaled_sync
  log_message(
    "Clustering query cells to reference centroids",
    verbose = verbose
  )
  Z_pca_query_cos <- symphony:::cosine_normalize_cpp(
    V = Z_pca_query,
    dim = 2
  )
  R_query <- symphony:::soft_cluster(
    Y = ref_obj$centroids,
    Z = Z_pca_query_cos,
    sigma = sigma
  )
  log_message(
    "Correcting query batch effects",
    verbose = verbose
  )
  if (!is.null(vars)) {
    design <- droplevels(metadata_query)[, vars] %>% as.data.frame()
    onehot <- design %>%
      purrr::map(function(.x) {
        if (length(unique(.x)) == 1) {
          rep(1, length(.x))
        } else {
          stats::model.matrix(~ 0 + .x)
        }
      }) %>%
      purrr::reduce(cbind)
    Xq <- cbind(1, intercept = onehot) %>% Matrix::t()
  } else {
    Xq <- Matrix::Matrix(
      rbind(rep(1, ncol(Z_pca_query)), rep(1, ncol(Z_pca_query))),
      sparse = TRUE
    )
  }
  Zq_corr <- symphony:::moe_correct_ref(
    Zq = as_matrix(Z_pca_query),
    Xq = as_matrix(Xq),
    Rq = as_matrix(R_query),
    Nr = as_matrix(ref_obj$cache[[1]]),
    RrZtr = as_matrix(ref_obj$cache[[2]])
  )
  colnames(Z_pca_query) <- row.names(metadata_query)
  rownames(Z_pca_query) <- paste0("PC_", seq_len(nrow(Zq_corr)))
  colnames(Zq_corr) <- row.names(metadata_query)
  rownames(Zq_corr) <- paste0("harmony_", seq_len(nrow(Zq_corr)))

  return(list(Z_pca_query = Z_pca_query, Zq_corr = Zq_corr, R_query = R_query))
}
