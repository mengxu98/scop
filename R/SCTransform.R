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
  fit <- stats::lm.fit(design, t(mat))
  resid <- t(fit$residuals)
  dimnames(resid) <- dimnames(mat)
  resid
}

sct_model_formula <- function(model_str) {
  stats::as.formula(sub("^\\s*y\\s*~", "~", model_str))
}

sct_single_offset_available <- function() {
  validated_glmgampoi <- c("1.18.0", "1.22.0")
  if (
    !requireNamespace("glmGamPoi", quietly = TRUE) ||
      !requireNamespace("sctransform", quietly = TRUE) ||
      !as.character(packageVersion("glmGamPoi")) %in% validated_glmgampoi ||
      packageVersion("sctransform") != "0.4.3"
  ) {
    return(FALSE)
  }
  required <- c(
    "get_groups_for_model_matrix",
    "estimate_dispersions_roughly",
    "estimate_betas_roughly_group_wise",
    "estimate_betas_group_wise",
    "estimate_betas_group_wise_optimize_helper",
    "calculate_mu",
    "overdispersion_mle",
    "overdispersion_shrinkage"
  )
  gp_ns <- asNamespace("glmGamPoi")
  all(vapply(required, exists, logical(1), envir = gp_ns, inherits = FALSE))
}

sct_fit_glmgampoi_offset_single <- function(umi, data) {
  gp_ns <- asNamespace("glmGamPoi")
  model_matrix <- matrix(1, nrow = ncol(umi), ncol = 1L)
  log_umi <- log(10^data$log_umi)
  offset_matrix <- matrix(
    rep(log_umi, each = nrow(umi)),
    nrow = nrow(umi),
    ncol = ncol(umi)
  )
  groups <- get("get_groups_for_model_matrix", gp_ns)(model_matrix)
  disp_init <- get("estimate_dispersions_roughly", gp_ns)(
    umi,
    model_matrix,
    offset_matrix
  )
  beta_group_init <- get("estimate_betas_roughly_group_wise", gp_ns)(
    umi,
    offset_matrix,
    groups
  )[, 1L]
  beta_res <- sct_fit_beta_intercept_offset(
    umi,
    log_offset = log_umi,
    dispersions = disp_init,
    beta_start = beta_group_init,
    return_mu = TRUE
  )
  if (any(!beta_res$converged)) {
    fallback <- which(!beta_res$converged)
    optimize_beta <- get(
      "estimate_betas_group_wise_optimize_helper",
      gp_ns
    )
    beta_res$beta[fallback] <- vapply(fallback, function(idx) {
      optimize_beta(umi[idx, ], log_umi, disp_init[idx])
    }, numeric(1))
    beta_res$mu[fallback, ] <- exp(
      beta_res$beta[fallback] + rep(log_umi, each = length(fallback))
    )
  }
  beta <- beta_res$beta
  mu <- beta_res$mu
  rm(offset_matrix, beta_group_init, groups)
  disp_est <- sct_overdispersion_mle_intercept(umi, mu)
  unstable_dispersion <- !is.finite(disp_est) |
    disp_est < 1e-3 |
    disp_est > 10
  if (any(unstable_dispersion)) {
    reference_dispersion <- get("overdispersion_mle", gp_ns)(
      umi[unstable_dispersion, , drop = FALSE],
      mu[unstable_dispersion, , drop = FALSE],
      model_matrix = model_matrix,
      do_cox_reid_adjustment = TRUE,
      subsample = FALSE
    )$estimate
    disp_est[unstable_dispersion] <- reference_dispersion
  }
  shrinkage <- get("overdispersion_shrinkage", gp_ns)(
    disp_est,
    gene_means = rowMeans(mu),
    df = ncol(umi) - ncol(model_matrix),
    ql_disp_trend = length(disp_est) >= 100L,
    npoints = max(0.1 * length(disp_est), 100),
    verbose = FALSE
  )
  beta_res <- sct_fit_beta_intercept_offset(
    umi,
    log_offset = log_umi,
    dispersions = shrinkage$dispersion_trend,
    beta_start = beta
  )
  if (any(!beta_res$converged)) {
    fallback <- which(!beta_res$converged)
    optimize_beta <- get(
      "estimate_betas_group_wise_optimize_helper",
      gp_ns
    )
    beta_res$beta[fallback] <- vapply(fallback, function(idx) {
      optimize_beta(umi[idx, ], log_umi, shrinkage$dispersion_trend[idx])
    }, numeric(1))
  }
  model_pars <- cbind(
    1 / disp_est,
    beta_res$beta,
    rep(log(10), nrow(umi))
  )
  dimnames(model_pars) <- list(
    rownames(umi),
    c("theta", "(Intercept)", "log_umi")
  )
  model_pars
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
  single_offset_fast <-
    is.null(cluster) &&
      identical(method, "glmGamPoi_offset") &&
      isTRUE(exclude_poisson) &&
      identical(model_str, "y ~ log_umi") &&
      !isTRUE(fix_intercept) &&
      !isTRUE(fix_slope) &&
      sct_single_offset_available()
  if (startsWith(method, "offset") || (is.null(cluster) && !single_offset_fast)) {
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
    if (single_offset_fast) {
      model_pars[[i]] <- sct_fit_glmgampoi_offset_single(
        umi = umi_bin,
        data = data_step1
      )
      next
    }
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
  if (length(res) > 0L) {
    res[res < res_clip_range[1L]] <- res_clip_range[1L]
    res[res > res_clip_range[2L]] <- res_clip_range[2L]
  }
  rv$y <- res
  rm(res)
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
  sct_cluster <- NULL
  sct_cores <- 1L

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
    x_vals <- umi[genes, , drop = FALSE]@x
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
