sct_scalar_pos <- function(x) {
  is.numeric(x) &&
    length(x) == 1L &&
    is.finite(x) &&
    x > 0
}

sct_clip_ok <- function(clip.range) {
  is.numeric(clip.range) &&
    length(clip.range) == 2L &&
    all(is.finite(clip.range)) &&
    clip.range[1] < clip.range[2]
}

sct_regression_supported <- function(vars.to.regress, latent.data, cell.attr, cells) {
  if (is.null(vars.to.regress) && is.null(latent.data)) {
    return(TRUE)
  }
  if (
    !is.character(cells) || anyNA(cells) || anyDuplicated(cells) ||
      (!is.null(vars.to.regress) && (
        !is.character(vars.to.regress) || !length(vars.to.regress) ||
          anyNA(vars.to.regress) || anyDuplicated(vars.to.regress) ||
          !is.data.frame(cell.attr) ||
          !all(vars.to.regress %in% colnames(cell.attr))
      ))
  ) {
    return(FALSE)
  }
  regressors <- if (is.null(vars.to.regress)) {
    NULL
  } else {
    cell.attr[, vars.to.regress, drop = FALSE]
  }
  if (!is.null(latent.data)) {
    latent <- tryCatch(as.data.frame(latent.data), error = function(e) NULL)
    if (is.null(latent) || !ncol(latent)) {
      return(FALSE)
    }
    latent <- latent[, !colnames(latent) %in% colnames(regressors), drop = FALSE]
    regressors <- if (is.null(regressors)) latent else cbind(regressors, latent)
  }
  if (is.null(regressors) || nrow(regressors) != length(cells)) {
    return(FALSE)
  }
  if (is.null(rownames(regressors))) {
    rownames(regressors) <- cells
  } else if (!setequal(rownames(regressors), cells)) {
    return(FALSE)
  }
  regressors <- regressors[cells, , drop = FALSE]
  if (any(!stats::complete.cases(regressors))) {
    return(FALSE)
  }
  names(regressors) <- paste0("v", seq_len(ncol(regressors)))
  design <- tryCatch(
    stats::model.matrix(
      stats::as.formula(paste("~", paste(names(regressors), collapse = " + "))),
      data = regressors
    ),
    error = function(e) NULL
  )
  !is.null(design) && all(is.finite(design))
}

sct_fast_path_supported <- function(
  reference.SCT.model,
  do.correct.umi,
  residual.features,
  conserve.memory,
  vst.flavor,
  do.scale,
  do.center,
  return.only.var.genes,
  extra_args,
  ncells,
  variable.features.n,
  variable.features.rv.th,
  clip.range,
  regression_ok = TRUE
) {
  is.null(reference.SCT.model) &&
    is.logical(do.correct.umi) && length(do.correct.umi) == 1L && !is.na(do.correct.umi) &&
    is.null(residual.features) &&
    isFALSE(conserve.memory) &&
    identical(vst.flavor, "v2") &&
    all(vapply(
      list(do.scale, do.center, return.only.var.genes),
      function(x) is.logical(x) && length(x) == 1L && !is.na(x),
      logical(1)
    )) &&
    length(extra_args) == 0L &&
    sct_scalar_pos(ncells) &&
    (is.null(variable.features.n) || sct_scalar_pos(variable.features.n)) &&
    sct_scalar_pos(variable.features.rv.th) &&
    sct_clip_ok(clip.range) &&
    isTRUE(regression_ok)
}

sct_regress_out <- function(mat, latent.df) {
  renamed <- latent.df
  names(renamed) <- paste0("v", seq_len(ncol(latent.df)))
  design <- stats::model.matrix(
    stats::as.formula(paste("~", paste(names(renamed), collapse = " + "))),
    data = renamed
  )
  # Normal equations are much faster for SCTransform's wide dense matrices
  # (genes x cells, cells >> genes) than lm.fit(t(mat)). They match the QR
  # residual path to numerical noise only for well-conditioned designs, so
  # keep QR for rank-deficient or ill-conditioned ones.
  xtx <- crossprod(design)
  xty <- t(mat %*% design)
  beta <- tryCatch(
    if (rcond(xtx) < 1e-8) NULL else solve(xtx, xty),
    error = function(e) NULL
  )
  if (is.null(beta) || any(!is.finite(beta))) {
    fit <- stats::lm.fit(design, t(mat))
    resid <- t(fit$residuals)
  } else {
    resid <- mat - t(beta) %*% t(design)
  }
  dimnames(resid) <- dimnames(mat)
  resid
}

sct_model_formula <- function(model_str) {
  stats::as.formula(sub("^\\s*y\\s*~", "~", model_str))
}


#' @export
SCTransform.default <- function(
  object,
  cell.attr,
  reference.SCT.model = NULL,
  do.correct.umi = TRUE,
  ncells = 5000,
  residual.features = NULL,
  variable.features.n = 3000,
  variable.features.rv.th = 1.3,
  vars.to.regress = NULL,
  latent.data = NULL,
  do.scale = FALSE,
  do.center = TRUE,
  clip.range = c(
    -sqrt(ncol(object) / 30),
    sqrt(ncol(object) / 30)
  ),
  vst.flavor = "v2",
  conserve.memory = FALSE,
  return.only.var.genes = TRUE,
  seed.use = 1448145,
  verbose = TRUE,
  cores = 1L,
  ...
) {
  extra_args <- list(...)
  sct_delegates <- !sct_fast_path_supported(
    reference.SCT.model = reference.SCT.model,
    do.correct.umi = do.correct.umi,
    residual.features = residual.features,
    conserve.memory = conserve.memory,
    vst.flavor = vst.flavor,
    do.scale = do.scale,
    do.center = do.center,
    return.only.var.genes = return.only.var.genes,
    extra_args = extra_args,
    ncells = ncells,
    variable.features.n = variable.features.n,
    variable.features.rv.th = variable.features.rv.th,
    clip.range = clip.range,
    regression_ok = sct_regression_supported(
      vars.to.regress = vars.to.regress,
      latent.data = latent.data,
      cell.attr = cell.attr,
      cells = colnames(object)
    )
  )
  if (sct_delegates) {
    log_message(
      "{.fn SCTransform} received arguments beyond the validated native path; delegating to Seurat.",
      message_type = "info"
    )
    seurat_sct <- utils::getFromNamespace("SCTransform.default", "Seurat")
    return(seurat_sct(
      object = object,
      cell.attr = cell.attr,
      reference.SCT.model = reference.SCT.model,
      do.correct.umi = do.correct.umi,
      ncells = ncells,
      residual.features = residual.features,
      variable.features.n = variable.features.n,
      variable.features.rv.th = variable.features.rv.th,
      vars.to.regress = vars.to.regress,
      latent.data = latent.data,
      do.scale = do.scale,
      do.center = do.center,
      clip.range = clip.range,
      vst.flavor = vst.flavor,
      conserve.memory = conserve.memory,
      return.only.var.genes = return.only.var.genes,
      seed.use = seed.use,
      verbose = verbose,
      ...
    ))
  }
  sct_latent_df <- NULL
  if (!is.null(vars.to.regress) && !is.null(cell.attr)) {
    sct_latent_df <- cell.attr[, vars.to.regress, drop = FALSE]
  }
  if (!is.null(latent.data)) {
    latent.df <- as.data.frame(latent.data)
    latent.df <- latent.df[, !colnames(latent.df) %in%
      (if (is.null(sct_latent_df)) character(0) else colnames(sct_latent_df)),
    drop = FALSE
    ]
    if (ncol(latent.df) > 0L) {
      if (is.null(sct_latent_df)) {
        sct_latent_df <- latent.df
      } else {
        if (is.null(rownames(latent.df))) {
          if (nrow(latent.df) != nrow(sct_latent_df)) {
            log_message(
              "latent.data has {nrow(latent.df)} rows but the cell attributes have {nrow(sct_latent_df)}.",
              message_type = "error"
            )
          }
          rownames(latent.df) <- rownames(sct_latent_df)
        } else if (!setequal(rownames(latent.df), rownames(sct_latent_df))) {
          log_message(
            "latent.data rownames do not match the cell attribute rownames.",
            message_type = "error"
          )
        }
        latent.df <- latent.df[rownames(sct_latent_df), , drop = FALSE]
        sct_latent_df <- cbind(sct_latent_df, latent.df)
      }
    }
  }
  if (!is.null(sct_latent_df)) {
    if (any(!stats::complete.cases(sct_latent_df))) {
      log_message(
        "Regression variables contain missing values; SCTransform requires complete cases.",
        message_type = "error"
      )
    }
    sct_cells <- colnames(SeuratObject::as.sparse(object))
    if (is.null(rownames(sct_latent_df))) {
      if (nrow(sct_latent_df) != length(sct_cells)) {
        log_message(
          "Regression variables have {nrow(sct_latent_df)} rows for {length(sct_cells)} cells.",
          message_type = "error"
        )
      }
      rownames(sct_latent_df) <- sct_cells
    } else if (!setequal(rownames(sct_latent_df), sct_cells)) {
      log_message(
        "Regression variable rownames do not match cell names.",
        message_type = "error"
      )
    }
    sct_latent_df <- sct_latent_df[sct_cells, , drop = FALSE]
  }
  check_r("sctransform", verbose = FALSE)
  check_r("glmGamPoi", verbose = FALSE)

  if (!is.null(seed.use)) {
    set.seed(seed.use)
  }
  vst.args <- list(...)
  object <- SeuratObject::as.sparse(object)
  umi <- object
  vst.args[["vst.flavor"]] <- vst.flavor
  vst.args[["umi"]] <- umi
  vst.args[["cell_attr"]] <- cell.attr
  vst.args[["verbosity"]] <- as.numeric(verbose) * 1
  vst.args[["return_cell_attr"]] <- TRUE
  vst.args[["return_gene_attr"]] <- FALSE
  vst.args[["return_corrected_umi"]] <- FALSE
  vst.args[["residual_type"]] <- "none"
  vst.args[["n_cells"]] <- min(ncells, ncol(umi))
  vst.out <- do.call(
    get_namespace_fun("sctransform", "vst"),
    args = vst.args
  )

  regressor_data_orig <- model.matrix(
    sct_model_formula(vst.out$model_str),
    vst.out$cell_attr
  )
  cell_attr_corr <- vst.out$cell_attr
  latent_var <- vst.out$arguments$latent_var
  cell_attr_corr[, latent_var] <- apply(
    cell_attr_corr[, latent_var, drop = FALSE],
    2,
    function(x) rep(median(x), length(x))
  )
  regressor_data_corr <- model.matrix(
    sct_model_formula(vst.out$model_str),
    cell_attr_corr
  )
  model_pars_fit <- vst.out$model_pars_fit
  sample_coefs <- model_pars_fit[1, -1, drop = FALSE]
  slope_col <- ncol(regressor_data_orig)
  corr_factor <- as.numeric(exp(
    (regressor_data_corr - regressor_data_orig) %*% t(sample_coefs)
  ))
  cell_mu_base <- as.numeric(exp(
    model_pars_fit[1, slope_col + 1] * regressor_data_orig[, slope_col]
  ))
  genes <- rownames(umi)[rownames(umi) %in% rownames(model_pars_fit)]
  if (identical(vst.out$arguments$min_variance, "umi_median")) {
    x_vals <- if (identical(genes, rownames(umi))) {
      umi@x
    } else {
      umi[genes, , drop = FALSE]@x
    }
    n <- length(x_vals)
    half <- n %/% 2L
    min_var <- if (n %% 2L == 1L) {
      (sort(x_vals, partial = half + 1L)[half + 1L] / 5)^2
    } else {
      sorted_partial <- sort(x_vals, partial = c(half, half + 1L))
      ((sorted_partial[half] + sorted_partial[half + 1L]) / 2 / 5)^2
    }
    rm(x_vals)
  } else {
    min_var <- as.numeric(vst.out$arguments$min_variance)
  }
  res.clip.range <- c(-sqrt(ncol(umi)), sqrt(ncol(umi)))
  col_names <- colnames(umi)
  all_gene_names <- rownames(umi)
  csr <- csc_to_csr(umi@i, umi@p, umi@x, nrow(umi), ncol(umi))
  bin_ind <- ceiling(seq_along(genes) / 7000L)
  max_bin <- max(bin_ind)
  corrected_list <- vector("list", max_bin)
  res_var <- numeric(length(genes))
  names(res_var) <- genes
  res_mean <- numeric(length(genes))
  names(res_mean) <- genes
  for (bin in seq_len(max_bin)) {
    genes_bin <- genes[bin_ind == bin]
    result <- sct_stats_correct_sparse(
      as.numeric(model_pars_fit[genes_bin, 2]),
      cell_mu_base,
      csr$row_ptr,
      csr$col_idx,
      csr$vals,
      match(genes_bin, all_gene_names) - 1L,
      model_pars_fit[genes_bin, 1],
      corr_factor,
      min_var,
      res.clip.range[1],
      res.clip.range[2],
      do.correct.umi
    )
    res_mean[genes_bin] <- result$res_mean
    res_var[genes_bin] <- result$res_var
    if (do.correct.umi) {
      corrected_list[[bin]] <- methods::new(
        "dgCMatrix",
        i = result$csc_i,
        p = result$csc_p,
        x = result$csc_x,
        Dim = c(length(genes_bin), ncol(umi)),
        Dimnames = list(genes_bin, col_names)
      )
    }
    rm(result)
  }
  vst.out$umi_corrected <- if (do.correct.umi) {
    do.call(rbind, corrected_list)
  } else {
    umi
  }

  gene_attr <- data.frame(
    residual_mean = res_mean,
    residual_variance = res_var,
    row.names = genes
  )
  vst.out$gene_attr <- gene_attr
  feature.variance <- sort(res_var, decreasing = TRUE)
  top.features <- if (is.null(variable.features.n)) {
    names(feature.variance)[feature.variance >= variable.features.rv.th]
  } else {
    names(feature.variance)[seq_len(min(variable.features.n, length(feature.variance)))]
  }
  output.features <- if (isFALSE(return.only.var.genes)) {
    genes
  } else {
    top.features
  }

  rm(umi)
  scale.data <- sct_fused_resid_center_sparse(
    as.numeric(model_pars_fit[output.features, 2]),
    cell_mu_base,
    csr$row_ptr,
    csr$col_idx,
    csr$vals,
    match(output.features, all_gene_names) - 1L,
    model_pars_fit[output.features, 1],
    min_var,
    res.clip.range[1],
    res.clip.range[2],
    clip.range[1],
    clip.range[2],
    isTRUE(do.center) && is.null(sct_latent_df),
    isTRUE(do.scale) && is.null(sct_latent_df)
  )
  dimnames(scale.data) <- list(output.features, col_names)
  rm(csr, corrected_list)
  if (!is.null(sct_latent_df)) {
    scale.data <- sct_regress_out(scale.data, sct_latent_df[col_names, , drop = FALSE])
    vst.out$arguments$sct.latent.vars <- colnames(sct_latent_df)
    scale.data <- sct_fastrowscale(scale.data, do.scale, do.center, Inf)
  }
  vst.out$y <- scale.data
  vst.out$variable_features <- top.features
  vst.out
}

#' @export
SCTransform.Seurat <- function(
  object,
  assay = "RNA",
  new.assay.name = "SCT",
  reference.SCT.model = NULL,
  do.correct.umi = TRUE,
  ncells = 5000,
  residual.features = NULL,
  variable.features.n = 3000,
  variable.features.rv.th = 1.3,
  vars.to.regress = NULL,
  do.scale = FALSE,
  do.center = TRUE,
  clip.range = c(
    -sqrt(ncol(object[[assay]]) / 30),
    sqrt(ncol(object[[assay]]) / 30)
  ),
  vst.flavor = "v2",
  conserve.memory = FALSE,
  return.only.var.genes = TRUE,
  seed.use = 1448145,
  verbose = TRUE,
  cores = 1L,
  ...
) {
  if (is.null(assay) || length(assay) != 1L || identical(assay, "SCT")) {
    log_message("SCTransform.Seurat requires a single non-SCT assay.", message_type = "error")
  }
  extra_args <- list(...)
  sct_delegates <- !sct_fast_path_supported(
    reference.SCT.model = reference.SCT.model,
    do.correct.umi = do.correct.umi,
    residual.features = residual.features,
    conserve.memory = conserve.memory,
    vst.flavor = vst.flavor,
    do.scale = do.scale,
    do.center = do.center,
    return.only.var.genes = return.only.var.genes,
    extra_args = extra_args,
    ncells = ncells,
    variable.features.n = variable.features.n,
    variable.features.rv.th = variable.features.rv.th,
    clip.range = clip.range,
    regression_ok = sct_regression_supported(
      vars.to.regress = vars.to.regress,
      latent.data = NULL,
      cell.attr = methods::slot(object, "meta.data"),
      cells = colnames(object[[assay]])
    )
  )
  if (sct_delegates) {
    log_message(
      "{.fn SCTransform} received arguments beyond the validated native path; delegating to Seurat.",
      message_type = "info"
    )
    seurat_sct <- utils::getFromNamespace("SCTransform.Seurat", "Seurat")
    return(seurat_sct(
      object = object,
      assay = assay,
      new.assay.name = new.assay.name,
      reference.SCT.model = reference.SCT.model,
      do.correct.umi = do.correct.umi,
      ncells = ncells,
      residual.features = residual.features,
      variable.features.n = variable.features.n,
      variable.features.rv.th = variable.features.rv.th,
      vars.to.regress = vars.to.regress,
      do.scale = do.scale,
      do.center = do.center,
      clip.range = clip.range,
      vst.flavor = vst.flavor,
      conserve.memory = conserve.memory,
      return.only.var.genes = return.only.var.genes,
      seed.use = seed.use,
      verbose = verbose,
      ...
    ))
  }
  if (!is.null(seed.use)) {
    set.seed(seed.use)
  }
  cell.attr <- methods::slot(object, "meta.data")[
    colnames(object[[assay]]),
  ]
  assay_object <- object[[assay]]
  umi <- tryCatch(
    SeuratObject::GetAssayData(assay_object, layer = "counts"),
    error = function(e) {
      if (!grepl("multiple layers", conditionMessage(e), fixed = TRUE)) {
        stop(e)
      }
      joined <- tryCatch(
        SeuratObject::JoinLayers(object = assay_object, layers = "counts"),
        error = function(e2) {
          log_message(
            "Failed to join Seurat v5 assay layers before reading {.val counts}: {conditionMessage(e2)}",
            message_type = "error"
          )
        }
      )
      SeuratObject::GetAssayData(joined, layer = "counts")
    }
  )
  vst.out <- SCTransform.default(
    object = umi,
    cell.attr = cell.attr,
    reference.SCT.model = reference.SCT.model,
    do.correct.umi = do.correct.umi,
    ncells = ncells,
    residual.features = residual.features,
    variable.features.n = variable.features.n,
    variable.features.rv.th = variable.features.rv.th,
    vars.to.regress = vars.to.regress,
    latent.data = NULL,
    do.scale = do.scale,
    do.center = do.center,
    clip.range = clip.range,
    vst.flavor = vst.flavor,
    conserve.memory = conserve.memory,
    return.only.var.genes = return.only.var.genes,
    seed.use = seed.use,
    verbose = verbose,
    cores = cores,
    ...
  )
  rm(umi)

  assay.out <- SeuratObject::CreateAssayObject(counts = vst.out$umi_corrected)
  data_layer <- vst.out$umi_corrected
  vst.out$umi_corrected <- NULL
  data_layer@x <- log1p(data_layer@x)
  SeuratObject::VariableFeatures(assay.out) <- vst.out$variable_features
  methods::slot(assay.out, "data") <- data_layer
  rm(data_layer)
  scale.data <- vst.out$y[
    rownames(assay.out)[rownames(assay.out) %in% rownames(vst.out$y)],
    ,
    drop = FALSE
  ]
  methods::slot(assay.out, "scale.data") <- scale.data
  rm(scale.data)
  vst.out$y <- NULL
  vst.out$arguments$sct.clip.range <- clip.range
  vst.out$arguments <- vst.out$arguments[
    !vapply(vst.out$arguments, is.null, logical(1))
  ]
  SeuratObject::Misc(assay.out, slot = "vst.out") <- vst.out
  rm(vst.out)
  old_validate <- getOption("Seurat.object.validate", default = TRUE)
  on.exit(options(Seurat.object.validate = old_validate), add = TRUE)
  options(Seurat.object.validate = FALSE)
  assay.out <- methods::as(assay.out, "SCTAssay")
  SCTAssay_fn <- utils::getFromNamespace("SCTAssay", "Seurat")
  assay.out <- SCTAssay_fn(assay.out, assay.orig = assay)
  methods::slot(
    methods::slot(assay.out, "SCTModel.list")[[1]],
    "umi.assay"
  ) <- assay
  SeuratObject::Key(assay.out) <- tolower(paste0(new.assay.name, "_"))
  assays_list <- methods::slot(object, "assays")
  assays_list[[new.assay.name]] <- assay.out
  methods::slot(object, "assays") <- assays_list
  rm(assays_list, assay.out)
  methods::slot(object, "active.assay") <- new.assay.name
  SeuratObject::LogSeuratCommand(object)
}

#' Apply SCTransform normalization
#'
#' @details
#' The validated sparse `vst.flavor = "v2"` path uses native corrected-count
#' and residual kernels. Reference models, specified residual features,
#' memory-conserving mode, unsupported regression designs, custom VST arguments,
#' and other flavors transparently delegate to `Seurat::SCTransform`.
#'
#' @param object Object containing count data.
#' @param ... Passed to methods.
#'
#' @return An object with SCTransform results.
#' @export
SCTransform <- function(object, ...) {
  UseMethod("SCTransform")
}
