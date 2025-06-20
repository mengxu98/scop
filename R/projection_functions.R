buildReferenceFromSeurat <- function(
    obj,
    assay = "RNA",
    pca = "pca",
    pca_dims = NULL,
    harmony = "harmony",
    umap = "umap") {
  if (!assay %in% c("RNA", "SCT")) {
    log_message(
      "Only supported assays are RNA or SCT.",
      message_type = "error"
    )
  }
  if (is.null(pca_dims)) {
    pca_dims <- seq_len(
      ncol(
        SeuratObject::Embeddings(obj, pca)
      )
    )
  }
  res <- list()
  ## TODO: check that these objects are all correctly initialized
  res$Z_corr <- Matrix::t(
    SeuratObject::Embeddings(obj, harmony)
  )
  res$Z_orig <- Matrix::t(
    SeuratObject::Embeddings(obj, pca)[, pca_dims, drop = FALSE]
  )
  log_message("Saved embeddings")

  res$R <- Matrix::t(obj[[harmony]]@misc$R)
  log_message("Saved soft cluster assignments")

  var_features <- SeuratObject::VariableFeatures(obj)

  if (assay == "RNA") {
    vargenes_means_sds <- data.frame(
      symbol = var_features,
      mean = Matrix::rowMeans(
        GetAssayData5(
          obj,
          assay = assay,
          layer = "data"
        )[
          var_features,
        ]
      )
    )

    vargenes_means_sds$stddev <- symphony::rowSDs(
      A = GetAssayData5(
        obj,
        assay = assay,
        layer = "data"
      )[var_features, ],
      row_means = vargenes_means_sds$mean
    )
  } else if (assay == "SCT") {
    vargenes_means_sds <- data.frame(
      symbol = var_features,
      mean = Matrix::rowMeans(
        GetAssayData5(
          obj,
          assay = assay,
          layer = "scale.data"
        )[
          var_features,
        ]
      )
    )
    asdgc <- Matrix::Matrix(
      GetAssayData5(
        obj,
        assay = assay,
        layer = "scale.data"
      )[var_features, ],
      sparse = TRUE
    )
    vargenes_means_sds$stddev <- symphony::rowSDs(
      asdgc,
      vargenes_means_sds$mean
    )
  }

  res$vargenes_means_sds <- vargenes_means_sds
  log_message(
    "Saved variable gene information for ",
    nrow(vargenes_means_sds),
    " genes."
  )

  res$loadings <- obj[[pca]]@feature.loadings[, pca_dims, drop = FALSE]
  log_message("Saved PCA loadings.")

  res$meta_data <- obj@meta.data
  log_message("Saved metadata.")

  if (is.null(obj[[umap]]@misc$model)) {
    log_message(
      "uwot model not initialiazed in Seurat object. Please do RunUMAP with umap.method='uwot', return.model=TRUE first.",
      message_type = "error"
    )
  }
  res$umap <- obj[[umap]]@misc$model

  ## Build Reference!
  log_message("Calculate final L2 normalized reference centroids (Y_cos)")
  res$centroids <- Matrix::t(
    symphony:::cosine_normalize_cpp(
      V = res$R %*% Matrix::t(res$Z_corr),
      dim = 1
    )
  )
  log_message("Calculate reference compression terms (Nr and C)")
  res$cache <- symphony:::compute_ref_cache(
    Rr = res$R,
    Zr = res$Z_corr
  )
  colnames(res$Z_orig) <- row.names(res$meta_data)
  rownames(res$Z_orig) <- paste0(
    SeuratObject::Key(
      obj[[pca]]
    ), seq_len(nrow(res$Z_corr))
  )
  colnames(res$Z_corr) <- row.names(res$meta_data)
  rownames(res$Z_corr) <- paste0(
    SeuratObject::Key(
      obj[[harmony]]
    ),
    seq_len(nrow(res$Z_corr))
  )
  log_message("Finished nicely.")
  return(res)
}

mapQuery <- function(
    exp_query,
    metadata_query,
    ref_obj,
    vars = NULL,
    sigma = 0.1,
    verbose = TRUE) {
  if (verbose) {
    log_message("Scaling and synchronizing query gene expression")
  }
  idx_shared_genes <- which(ref_obj$vargenes$symbol %in% rownames(exp_query))
  shared_genes <- ref_obj$vargenes$symbol[idx_shared_genes]
  if (verbose) {
    log_message(
      "Found ",
      length(shared_genes),
      " reference variable genes in query dataset"
    )
  }
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
  if (verbose) {
    log_message("Project query cells using reference gene loadings")
  }
  Z_pca_query <- Matrix::t(ref_obj$loadings) %*% exp_query_scaled_sync
  if (verbose) {
    log_message("Clustering query cells to reference centroids")
  }
  Z_pca_query_cos <- symphony:::cosine_normalize_cpp(
    V = Z_pca_query,
    dim = 2
  )
  R_query <- symphony:::soft_cluster(
    Y = ref_obj$centroids,
    Z = Z_pca_query_cos,
    sigma = sigma
  )
  if (verbose) {
    log_message("Correcting query batch effects")
  }
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
    Zq = as.matrix(Z_pca_query),
    Xq = as.matrix(Xq),
    Rq = as.matrix(R_query),
    Nr = as.matrix(ref_obj$cache[[1]]),
    RrZtr = as.matrix(ref_obj$cache[[2]])
  )
  colnames(Z_pca_query) <- row.names(metadata_query)
  rownames(Z_pca_query) <- paste0("PC_", seq_len(nrow(Zq_corr)))
  colnames(Zq_corr) <- row.names(metadata_query)
  rownames(Zq_corr) <- paste0("harmony_", seq_len(nrow(Zq_corr)))

  return(list(Z_pca_query = Z_pca_query, Zq_corr = Zq_corr, R_query = R_query))
}
