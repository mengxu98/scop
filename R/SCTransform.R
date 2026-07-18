sct_is_default_clip <- function(clip.range, n_cells) {
  if (
    !is.numeric(clip.range) ||
      length(clip.range) != 2L ||
      any(!is.finite(clip.range))
  ) {
    return(FALSE)
  }
  isTRUE(all.equal(
    as.numeric(clip.range),
    c(-sqrt(n_cells / 30), sqrt(n_cells / 30)),
    tolerance = sqrt(.Machine$double.eps),
    check.attributes = FALSE
  ))
}

sct_model_formula <- function(model_str) {
  stats::as.formula(sub("^\\s*y\\s*~", "~", model_str))
}

sct_clip_matrix_values <- function(x, range) {
  if (!length(x)) {
    return(x)
  }
  x[x < range[1]] <- range[1]
  x[x > range[2]] <- range[2]
  x
}

sct_get_model_pars <- function(
  genes_step1,
  bin_size,
  umi,
  model_str,
  cells_step1,
  method,
  data_step1,
  theta_given,
  theta_estimation_fun,
  exclude_poisson = FALSE,
  fix_intercept = FALSE,
  fix_slope = FALSE,
  use_geometric_mean = TRUE,
  use_geometric_mean_offset = FALSE,
  verbosity = 0,
  cluster = NULL,
  cores = 1L
) {
  row_var <- get_namespace_fun("sctransform", "row_var")
  get_model_pars <- get_namespace_fun("sctransform", "get_model_pars")
  fit_glmGamPoi_offset <- get_namespace_fun(
    "sctransform",
    "fit_glmGamPoi_offset"
  )
  if (startsWith(method, "offset") || is.null(cluster)) {
    return(get_model_pars(
      genes_step1,
      bin_size,
      umi,
      model_str,
      cells_step1,
      method,
      data_step1,
      theta_given,
      theta_estimation_fun,
      exclude_poisson,
      fix_intercept,
      fix_slope,
      use_geometric_mean,
      use_geometric_mean_offset,
      verbosity
    ))
  }
  bin_ind <- ceiling(seq_along(genes_step1) / bin_size)
  max_bin <- max(bin_ind)
  model_pars <- list()
  cl <- cluster
  for (i in seq_len(max_bin)) {
    genes_bin_regress <- genes_step1[bin_ind == i]
    umi_bin <- as.matrix(umi[genes_bin_regress, cells_step1, drop = FALSE])
    n_genes_bin <- nrow(umi_bin)
    cores <- min(as.integer(cores), n_genes_bin)
    gene_indices <- split(
      seq_len(n_genes_bin),
      ceiling(
        seq_len(n_genes_bin) / (n_genes_bin / cores + 1e-10)
      )
    )
    chunk_list <- lapply(gene_indices, function(idx) {
      umi_bin[idx, , drop = FALSE]
    })
    task_env <- new.env(parent = globalenv())
    task_env$model_str <- model_str
    task_env$data_step1 <- data_step1
    task_env$exclude_poisson <- exclude_poisson
    task_env$fit_glmGamPoi_offset <- fit_glmGamPoi_offset
    task_fn <- function(chunk) {
      fit_glmGamPoi_offset(
        umi = chunk,
        model_str = model_str,
        data = data_step1,
        allow_inf_theta = exclude_poisson
      )
    }
    environment(task_fn) <- task_env
    par_results <- parallel::parLapply(cl, chunk_list, task_fn)
    model_pars[[i]] <- do.call(rbind, par_results)
  }
  model_pars <- do.call(rbind, model_pars)
  rownames(model_pars) <- genes_step1
  colnames(model_pars)[1] <- "theta"
  if (exclude_poisson) {
    umi_step1 <- umi[genes_step1, , drop = FALSE]
    genes_amean_step1 <- Matrix::rowMeans(umi_step1)
    genes_var_step1 <- row_var(umi_step1)
    predicted_theta <- genes_amean_step1^2 /
      (genes_var_step1 - genes_amean_step1)
    actual_theta <- model_pars[genes_step1, "theta"]
    diff_theta <- predicted_theta / actual_theta
    model_pars <- cbind(model_pars, diff_theta)
    diff_theta_index <- rownames(model_pars[
      model_pars[genes_step1, "diff_theta"] < 0.001,
    ])
    model_pars[diff_theta_index, 1] <- Inf
    model_pars <- model_pars[, -dim(model_pars)[2]]
  }
  model_pars
}

sct_vst <- function(
  umi,
  cell_attr = NULL,
  latent_var = c("log_umi"),
  batch_var = NULL,
  latent_var_nonreg = NULL,
  n_genes = 2000,
  n_cells = 5000,
  method = "poisson",
  do_regularize = TRUE,
  theta_given = NULL,
  theta_estimation_fun = "theta.ml",
  exclude_poisson = FALSE,
  use_geometric_mean = TRUE,
  use_geometric_mean_offset = FALSE,
  fix_intercept = FALSE,
  fix_slope = FALSE,
  scale_factor = NULL,
  vst.flavor = NULL,
  verbosity = 2,
  verbose = NULL,
  show_progress = NULL,
  residual_type = "pearson",
  return_cell_attr = FALSE,
  return_gene_attr = TRUE,
  return_corrected_umi = FALSE,
  min_cells = 5,
  gmean_eps = 1,
  theta_regularization = "od_factor",
  bin_size = 500,
  min_variance = -Inf,
  bw_adjust = 3,
  res_clip_range = c(-sqrt(ncol(umi)), sqrt(ncol(umi))),
  sct_cluster = NULL,
  sct_cores = 1L
) {
  make_cell_attr <- get_namespace_fun("sctransform", "make_cell_attr")
  reg_model_pars <- get_namespace_fun("sctransform", "reg_model_pars")
  row_gmean <- get_namespace_fun("sctransform", "row_gmean")
  row_var <- get_namespace_fun("sctransform", "row_var")
  get_model_pars <- function(
    genes_step1,
    bin_size,
    umi,
    model_str,
    cells_step1,
    method,
    data_step1,
    theta_given,
    theta_estimation_fun,
    exclude_poisson = FALSE,
    fix_intercept = FALSE,
    fix_slope = FALSE,
    use_geometric_mean = TRUE,
    use_geometric_mean_offset = FALSE,
    verbosity = 0
  ) {
    sct_get_model_pars(
      genes_step1 = genes_step1,
      bin_size = bin_size,
      umi = umi,
      model_str = model_str,
      cells_step1 = cells_step1,
      method = method,
      data_step1 = data_step1,
      theta_given = theta_given,
      theta_estimation_fun = theta_estimation_fun,
      exclude_poisson = exclude_poisson,
      fix_intercept = fix_intercept,
      fix_slope = fix_slope,
      use_geometric_mean = use_geometric_mean,
      use_geometric_mean_offset = use_geometric_mean_offset,
      verbosity = verbosity,
      cluster = sct_cluster,
      cores = sct_cores
    )
  }

  if (!is.null(vst.flavor)) {
    if (vst.flavor == "v2") {
      glmGamPoi_pkg <- paste0("glm", "GamPoi")
      check_r(glmGamPoi_pkg, verbose = FALSE)
      glmGamPoi_check <- TRUE
      method <- "glmGamPoi_offset"
      if (!glmGamPoi_check) {
        method <- "nb_offset"
      }
      exclude_poisson <- TRUE
      if (min_variance == -Inf) {
        min_variance <- "umi_median"
      }
      if (is.null(n_cells)) n_cells <- 2000
    }
  }
  arguments <- as.list(environment())
  arguments <- arguments[
    !names(arguments) %in%
      c(
        "umi",
        "cell_attr",
        "sct_cluster",
        "sct_cores"
      )
  ]
  if (startsWith(method, "offset")) {
    cell_attr <- NULL
    latent_var <- c("log_umi")
    batch_var <- NULL
    latent_var_nonreg <- NULL
    n_genes <- NULL
    n_cells <- NULL
    do_regularize <- FALSE
    if (is.null(theta_given)) {
      theta_given <- 100
    } else {
      theta_given <- theta_given[1]
    }
  }
  times <- list(start_time = Sys.time())
  cell_attr <- make_cell_attr(
    umi,
    cell_attr,
    latent_var,
    batch_var,
    latent_var_nonreg,
    verbosity
  )
  if (inherits(umi, "dgCMatrix")) {
    genes_cell_count <- tabulate(umi@i + 1L, nbins = nrow(umi))
    names(genes_cell_count) <- rownames(umi)
  } else {
    genes_cell_count <- Matrix::rowSums(umi >= 0.01)
  }
  genes <- rownames(umi)[genes_cell_count >= min_cells]
  umi <- umi[genes, ]
  if (use_geometric_mean) {
    genes_log_gmean <- log10(row_gmean(umi, eps = gmean_eps))
  } else {
    genes_log_gmean <- log10(Matrix::rowMeans(umi))
  }
  if (!do_regularize && !is.null(n_genes)) {
    n_genes <- NULL
  }
  if (!is.null(n_cells) && n_cells < ncol(umi)) {
    cells_step1 <- sample(colnames(umi), size = n_cells)
    if (inherits(umi, "dgCMatrix")) {
      umi_sub <- umi[, cells_step1]
      genes_cell_count_step1 <- tabulate(
        umi_sub@i + 1L,
        nbins = nrow(umi_sub)
      )
      names(genes_cell_count_step1) <- rownames(umi_sub)
      rm(umi_sub)
    } else {
      genes_cell_count_step1 <- Matrix::rowSums(umi[, cells_step1] > 0)
    }
    genes_step1 <- rownames(umi)[genes_cell_count_step1 >= min_cells]
    if (use_geometric_mean) {
      genes_log_gmean_step1 <- log10(row_gmean(
        umi[genes_step1, ],
        eps = gmean_eps
      ))
    } else {
      genes_log_gmean_step1 <- log10(Matrix::rowMeans(umi[genes_step1, ]))
    }
  } else {
    cells_step1 <- colnames(umi)
    genes_step1 <- genes
    genes_log_gmean_step1 <- genes_log_gmean
  }
  genes_amean <- NULL
  genes_var <- NULL
  if (do_regularize && exclude_poisson) {
    genes_amean <- Matrix::rowSums(umi) / ncol(umi)
    genes_var <- row_var(umi)
    overdispersion_factor <- genes_var - genes_amean
    overdispersion_factor_step1 <- overdispersion_factor[genes_step1]
    is_overdispersed <- overdispersion_factor_step1 > 0
    genes_step1 <- genes_step1[is_overdispersed]
    genes_log_gmean_step1 <- genes_log_gmean[genes_step1]
  }
  data_step1 <- cell_attr[cells_step1, , drop = FALSE]
  if (!is.null(n_genes) && n_genes < length(genes_step1)) {
    log_gmean_dens <- density(
      x = genes_log_gmean_step1,
      bw = "nrd",
      adjust = 1
    )
    sampling_prob <- 1 /
      (approx(
        x = log_gmean_dens$x,
        y = log_gmean_dens$y,
        xout = genes_log_gmean_step1
      )$y +
        .Machine$double.eps)
    genes_step1 <- sample(genes_step1, size = n_genes, prob = sampling_prob)
    if (use_geometric_mean) {
      genes_log_gmean_step1 <- log10(row_gmean(
        umi[genes_step1, ],
        eps = gmean_eps
      ))
    } else {
      genes_log_gmean_step1 <- log10(Matrix::rowMeans(umi[genes_step1, ]))
    }
  }
  model_str <- paste0("y ~ ", paste(latent_var, collapse = " + "))
  if (verbosity > 0) {
    message(
      "Variance stabilizing transformation of count matrix of size ",
      nrow(umi),
      " by ",
      ncol(umi)
    )
    message("Model formula is ", model_str)
  }
  times$get_model_pars <- Sys.time()
  model_pars <- get_model_pars(
    genes_step1,
    bin_size,
    umi,
    model_str,
    cells_step1,
    method,
    data_step1,
    theta_given,
    theta_estimation_fun,
    exclude_poisson,
    fix_intercept,
    fix_slope,
    use_geometric_mean,
    use_geometric_mean_offset,
    verbosity
  )
  min_theta <- 1e-7
  if (any(model_pars[, "theta"] < min_theta)) {
    model_pars[, "theta"] <- pmax(model_pars[, "theta"], min_theta)
  }
  times$reg_model_pars <- Sys.time()
  if (do_regularize) {
    model_pars_fit <- reg_model_pars(
      model_pars,
      genes_log_gmean_step1,
      genes_log_gmean,
      cell_attr,
      batch_var,
      cells_step1,
      genes_step1,
      umi,
      bw_adjust,
      gmean_eps,
      theta_regularization,
      genes_amean,
      genes_var,
      exclude_poisson,
      fix_intercept,
      fix_slope,
      use_geometric_mean,
      use_geometric_mean_offset,
      verbosity
    )
    model_pars_outliers <- attr(model_pars_fit, "outliers")
  } else {
    model_pars_fit <- model_pars
    model_pars_outliers <- rep(FALSE, nrow(model_pars))
  }
  regressor_data <- model.matrix(sct_model_formula(model_str), cell_attr)
  times$get_residuals <- Sys.time()
  res <- matrix(NA, nrow = 0, ncol = 0)
  rv <- list(
    y = res,
    model_str = model_str,
    model_pars = model_pars,
    model_pars_outliers = model_pars_outliers,
    model_pars_fit = model_pars_fit,
    model_str_nonreg = "",
    model_pars_nonreg = c(),
    arguments = arguments,
    genes_log_gmean_step1 = genes_log_gmean_step1,
    cells_step1 = cells_step1,
    cell_attr = cell_attr
  )
  rm(res)
  rv$y <- sct_clip_matrix_values(rv$y, res_clip_range)
  if (!return_cell_attr) {
    rv[["cell_attr"]] <- NULL
  }
  times$get_gene_attr <- Sys.time()
  times$done <- Sys.time()
  rv$times <- times
  rv
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
  if (
    !is.null(reference.SCT.model) ||
      !isTRUE(do.correct.umi) ||
      !is.numeric(ncells) ||
      length(ncells) != 1L ||
      !is.finite(ncells) ||
      ncells <= 0 ||
      !is.null(residual.features) ||
      is.null(variable.features.n) ||
      !is.numeric(variable.features.n) ||
      length(variable.features.n) != 1L ||
      !is.finite(variable.features.n) ||
      variable.features.n <= 0 ||
      !identical(variable.features.rv.th, 1.3) ||
      !is.null(vars.to.regress) ||
      !is.null(latent.data) ||
      !isFALSE(do.scale) ||
      !isTRUE(do.center) ||
      !sct_is_default_clip(clip.range, ncol(object)) ||
      !identical(vst.flavor, "v2") ||
      !isFALSE(conserve.memory) ||
      !isTRUE(return.only.var.genes) ||
      length(extra_args) != 0L
  ) {
    stop(
      "SCTransform.default received unsupported arguments for the scop implementation.",
      call. = FALSE
    )
  }
  check_r("sctransform", verbose = FALSE)
  check_r("glmGamPoi", verbose = FALSE)
  sct_cluster <- NULL
  sct_cores <- 1L
  requested_cores <- suppressWarnings(as.integer(cores))
  if (
    .Platform$OS.type == "windows" &&
      length(requested_cores) == 1L &&
      !is.na(requested_cores) &&
      requested_cores > 1L
  ) {
    available_cores <- as.integer(parallel::detectCores(logical = FALSE))
    if (is.na(available_cores) || available_cores < 1L) {
      available_cores <- 1L
    }
    if (available_cores > 1L) {
      sct_cores <- min(requested_cores, available_cores - 1L)
      sct_cluster <- tryCatch(
        {
          cl <- parallel::makeCluster(
            sct_cores,
            port = 20000L + sample.int(1000L, 1L)
          )
          message(sprintf(
            "[scop] SCTransform PSOCK ready (%d cores / %d available)",
            sct_cores,
            available_cores
          ))
          cl
        },
        error = function(e) {
          message(sprintf(
            "[scop] SCT PSOCK init failed (%s); GLM fitting will be sequential",
            conditionMessage(e)
          ))
          NULL
        }
      )
      if (is.null(sct_cluster)) {
        sct_cores <- 1L
      }
      if (!is.null(sct_cluster)) {
        on.exit(
          tryCatch(
            parallel::stopCluster(sct_cluster),
            error = function(e) NULL
          ),
          add = TRUE
        )
      }
    }
  }

  if (!is.null(seed.use)) {
    set.seed(seed.use)
  }
  vst.args <- list(...)
  object <- SeuratObject::as.sparse(object)
  umi <- object
  if (!is.null(vst.flavor) && vst.flavor == "v1") {
    vst.flavor <- NULL
  }
  vst.args[["vst.flavor"]] <- vst.flavor
  vst.args[["umi"]] <- umi
  vst.args[["cell_attr"]] <- cell.attr
  vst.args[["verbosity"]] <- as.numeric(verbose) * 1
  vst.args[["return_cell_attr"]] <- TRUE
  vst.args[["return_gene_attr"]] <- FALSE
  vst.args[["return_corrected_umi"]] <- FALSE
  vst.args[["residual_type"]] <- "none"
  vst.args[["n_cells"]] <- min(ncells, ncol(umi))
  vst.args[["bin_size"]] <- 15000
  vst.args[["sct_cluster"]] <- sct_cluster
  vst.args[["sct_cores"]] <- sct_cores

  orig_sct_vst <- get_namespace_fun("sctransform", "vst")
  orig_sct_rowgmean <- get_namespace_fun("sctransform", "row_gmean")
  row_gmean_cache <- new.env(parent = emptyenv())
  row_gmean_cache$full_result <- NULL
  row_gmean_cache$full_ncol <- NULL
  cached_row_gmean <- function(x, eps = 1) {
    if (!inherits(x, "dgCMatrix")) {
      return(orig_sct_rowgmean(x, eps = eps))
    }
    nc <- ncol(x)
    if (
      !is.null(row_gmean_cache$full_result) &&
        nc == row_gmean_cache$full_ncol
    ) {
      rn <- rownames(x)
      if (!is.null(rn) && all(rn %in% names(row_gmean_cache$full_result))) {
        return(row_gmean_cache$full_result[rn])
      }
    }
    value <- orig_sct_rowgmean(x, eps = eps)
    if (
      is.null(row_gmean_cache$full_result) ||
        length(value) > length(row_gmean_cache$full_result)
    ) {
      row_gmean_cache$full_result <- value
      row_gmean_cache$full_ncol <- nc
    }
    value
  }
  utils::assignInNamespace("vst", sct_vst, ns = "sctransform")
  utils::assignInNamespace("row_gmean", cached_row_gmean, ns = "sctransform")
  on.exit(
    {
      utils::assignInNamespace("vst", orig_sct_vst, ns = "sctransform")
      utils::assignInNamespace(
        "row_gmean",
        orig_sct_rowgmean,
        ns = "sctransform"
      )
    },
    add = TRUE
  )
  sct_vst_fun <- get_namespace_fun("sctransform", "vst")
  vst.out <- do.call(sct_vst_fun, args = vst.args)

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
    x_vals <- umi@x
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
  top.features <- names(feature.variance)[
    1:min(variable.features.n, length(feature.variance))
  ]
  feat_positions <- match(top.features, all_gene_names)
  top.features <- top.features[order(feat_positions)]

  rm(umi)
  scale.data <- sct_fused_resid_center_sparse(
    as.numeric(model_pars_fit[top.features, 2]),
    cell_mu_base,
    csr$row_ptr,
    csr$col_idx,
    csr$vals,
    match(top.features, all_gene_names) - 1L,
    model_pars_fit[top.features, 1],
    min_var,
    res.clip.range[1],
    res.clip.range[2],
    clip.range[1],
    clip.range[2]
  )
  dimnames(scale.data) <- list(top.features, col_names)
  rm(csr, corrected_list)
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
    stop("SCTransform.Seurat requires a single non-SCT assay.", call. = FALSE)
  }
  extra_args <- list(...)
  if (
    !is.null(reference.SCT.model) ||
      !isTRUE(do.correct.umi) ||
      !is.numeric(ncells) ||
      length(ncells) != 1L ||
      !is.finite(ncells) ||
      ncells <= 0 ||
      !is.null(residual.features) ||
      is.null(variable.features.n) ||
      !is.numeric(variable.features.n) ||
      length(variable.features.n) != 1L ||
      !is.finite(variable.features.n) ||
      variable.features.n <= 0 ||
      !identical(variable.features.rv.th, 1.3) ||
      !is.null(vars.to.regress) ||
      !isFALSE(do.scale) ||
      !isTRUE(do.center) ||
      !sct_is_default_clip(clip.range, ncol(object[[assay]])) ||
      !identical(vst.flavor, "v2") ||
      !isFALSE(conserve.memory) ||
      !isTRUE(return.only.var.genes) ||
      length(extra_args) != 0L
  ) {
    stop(
      "SCTransform.Seurat received unsupported arguments for the scop implementation.",
      call. = FALSE
    )
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
  methods::slot(assay.out, "scale.data") <- vst.out$y
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
#' @param object Object containing count data.
#' @param ... Passed to methods.
#'
#' @return An object with SCTransform results.
#' @export
SCTransform <- function(object, ...) {
  UseMethod("SCTransform")
}
