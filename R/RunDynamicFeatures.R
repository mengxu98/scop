#' @title Calculates dynamic features for lineages
#'
#' @md
#' @inheritParams thisutils::parallelize_fun
#' @inheritParams standard_scop
#' @inheritParams GroupHeatmap
#' @param lineages A character vector specifying the lineage names for which dynamic features should be calculated.
#' @param features A character vector of features to use.
#' If `NULL`, n_candidates must be provided.
#' @param suffix A character vector specifying the suffix to append to the output layer names for each lineage.
#' Default is the lineage names.
#' @param n_candidates A number of candidate features to select when features is `NULL`.
#' Default is `1000`.
#' @param minfreq An integer specifying the minimum frequency threshold for candidate features.
#' Features with a frequency less than minfreq will be excluded.
#' Default is `5`.
#' @param libsize A numeric or numeric vector specifying the library size correction factors for each cell.
#' If NULL, the library size correction factors will be calculated based on the expression matrix.
#' If length(libsize) is 1, the same value will be used for all cells.
#' Otherwise, libsize must have the same length as the number of cells in srt.
#' Default is `NULL`.
#' @param fit_method The method used for fitting features.
#' Either `"gam"` (generalized additive models) or `"pretsa"` (Pattern recognition in Temporal and Spatial Analyses).
#' Default is `"gam"`.
#' @param family A character or character vector specifying the family of distributions to use for the GAM.
#' If family is set to NULL, the appropriate family will be automatically determined based on the data.
#' If length(family) is 1, the same family will be used for all features.
#' Otherwise, family must have the same length as features.
#' @param knot For `fit_method = "pretsa"`: B-spline knots. `0` or `"auto"`.
#' Default is `0`.
#' @param max_knot_allowed For `fit_method = "pretsa"` when `knot = "auto"`: max knots.
#' Default is `10`.
#' @param padjust_method The method used for p-value adjustment.
#' Default is `"fdr"`.
#'
#' @return
#' Returns the modified Seurat object with the calculated dynamic features stored in the tools slot.
#'
#' @seealso
#' [DynamicHeatmap], [DynamicPlot], [RunDynamicEnrichment]
#'
#' @export
#'
#' @references
#' Zhuang, H., Ji, Z. PreTSA: computationally efficient modeling of temporal and spatial gene expression patterns.
#' Genome Biol (2026). https://doi.org/10.1186/s13059-026-03994-3
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunSlingshot(
#'   pancreas_sub,
#'   group.by = "SubCellType",
#'   reduction = "UMAP"
#' )
#'
#' pancreas_sub <- RunDynamicFeatures(
#'   pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   n_candidates = 200,
#'   fit_method = "gam"
#' )
#'
#' names(
#'   pancreas_sub@tools$DynamicFeatures_Lineage1
#' )
#' head(
#'   pancreas_sub@tools$DynamicFeatures_Lineage1$DynamicFeatures
#' )
#' ht <- DynamicHeatmap(
#'   pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   cell_annotation = "SubCellType",
#'   n_split = 3,
#'   reverse_ht = "Lineage1"
#' )
#' ht$plot
#'
#' DynamicPlot(
#'   pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   features = c("Arxes1", "Ncoa2"),
#'   group.by = "SubCellType",
#'   compare_lineages = TRUE,
#'   compare_features = FALSE
#' )
#'
#' pancreas_sub <- RunDynamicFeatures(
#'   pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   n_candidates = 200,
#'   fit_method = "pretsa"
#' )
#' head(
#'   pancreas_sub@tools$DynamicFeatures_Lineage1$DynamicFeatures
#' )
#' ht <- DynamicHeatmap(
#'   pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   cell_annotation = "SubCellType",
#'   n_split = 3,
#'   reverse_ht = "Lineage1"
#' )
#' ht$plot
#'
#' DynamicPlot(
#'   pancreas_sub,
#'   lineages = c("Lineage1", "Lineage2"),
#'   features = c("Arxes1", "Ncoa2"),
#'   group.by = "SubCellType",
#'   compare_lineages = TRUE,
#'   compare_features = FALSE
#' )
RunDynamicFeatures <- function(
    srt,
    lineages,
    features = NULL,
    suffix = lineages,
    n_candidates = 1000,
    minfreq = 5,
    family = NULL,
    layer = "counts",
    assay = NULL,
    libsize = NULL,
    fit_method = c("gam", "pretsa"),
    knot = 0,
    max_knot_allowed = 10,
    padjust_method = "fdr",
    cores = 1,
    verbose = TRUE,
    seed = 11) {
  set.seed(seed)
  assay <- assay %||% DefaultAssay(srt)
  fit_method <- match.arg(fit_method)

  log_message(
    "Start find dynamic features",
    verbose = verbose
  )

  if (fit_method == "gam") {
    check_r("mgcv", verbose = FALSE)
  }
  meta <- c()
  gene <- c()
  if (!is.null(features)) {
    gene <- features[features %in% rownames(srt[[assay]])]
    meta <- features[features %in% colnames(srt@meta.data)]
    isnum <- sapply(
      srt@meta.data[, meta, drop = FALSE], is.numeric
    )
    if (!all(isnum)) {
      log_message(
        "{.val {meta[!isnum]}} is not numeric and will be dropped",
        message_type = "warning",
        verbose = verbose
      )
      meta <- meta[isnum]
    }
    features <- c(gene, meta)
    if (length(features) == 0) {
      log_message(
        "No feature found in the srt object",
        message_type = "error"
      )
    }
  }

  y_mat <- GetAssayData5(
    srt,
    layer = layer,
    assay = assay
  )
  if (is.null(libsize)) {
    status <- CheckDataType(
      srt,
      assay = assay,
      layer = "counts",
      verbose = verbose
    )
    if (status != "raw_counts") {
      y_libsize <- stats::setNames(
        rep(1, ncol(srt)),
        colnames(srt)
      )
    } else {
      y_libsize <- Matrix::colSums(
        GetAssayData5(
          srt,
          assay = assay,
          layer = "counts"
        )
      )
    }
  } else {
    if (length(libsize) == 1) {
      y_libsize <- stats::setNames(
        rep(libsize, ncol(srt)),
        colnames(srt)
      )
    } else if (length(libsize) == ncol(srt)) {
      y_libsize <- stats::setNames(libsize, colnames(srt))
    } else {
      log_message(
        "{.arg libsize} must be length of 1 or the number of cells",
        message_type = "error"
      )
    }
  }

  if (length(meta) > 0) {
    y_mat <- rbind(
      y_mat, Matrix::t(srt@meta.data[, meta, drop = FALSE])
    )
  }

  features_list <- c()
  srt_sub_list <- list()
  for (l in lineages) {
    srt_sub <- subset(
      srt,
      cell = rownames(srt@meta.data)[is.finite(srt@meta.data[[l]])]
    )
    if (is.null(features)) {
      if (is.null(n_candidates)) {
        log_message(
          "{.arg features} or {.arg n_candidates} must be provided at least one",
          message_type = "error"
        )
      }
      HVF <- SeuratObject::VariableFeatures(
        Seurat::FindVariableFeatures(
          srt_sub,
          nfeatures = n_candidates,
          assay = assay,
          verbose = FALSE
        ),
        assay = assay
      )
      HVF_counts <- GetAssayData5(
        srt_sub,
        assay = assay,
        layer = "counts"
      )[HVF, , drop = FALSE]
      HVF <- HVF[
        apply(HVF_counts, 1, function(x) {
          length(unique(x))
        }) >=
          minfreq
      ]
      features_list[[l]] <- HVF
    } else {
      features_list[[l]] <- features
    }
    srt_sub_list[[l]] <- srt_sub
  }
  features <- unique(unlist(features_list))
  gene <- features[features %in% rownames(srt[[assay]])]
  meta <- features[features %in% colnames(srt@meta.data)]
  log_message(
    "Number of candidate features (union): {.val {length(features)}}",
    verbose = verbose
  )

  if (layer == "counts") {
    gene_status <- status
  }
  gene_status <- status <- CheckDataType(srt, assay = assay, layer = layer)
  meta_status <- sapply(meta, function(x) {
    CheckDataType(srt[[x]])
  })
  if (is.null(family)) {
    family <- rep("gaussian", length(features))
    names(family) <- features
    family[names(meta_status)[meta_status == "raw_counts"]] <- "nb"
    if (gene_status == "raw_counts") {
      family[gene] <- "nb"
    }
  } else {
    if (length(family) == 1) {
      family <- rep(family, length(features))
      names(family) <- features
    }
    if (length(family) != length(features)) {
      log_message(
        "{.arg family} must be one character or a vector of the same length as features",
        message_type = "error"
      )
    }
  }

  for (i in seq_along(lineages)) {
    l <- lineages[i]
    srt_sub <- srt_sub_list[[l]]
    t <- srt_sub[[l, drop = TRUE]]
    t <- t[is.finite(t)]
    t_ordered <- t[order(t)]
    y_ordered <- as_matrix(
      y_mat[features, names(t_ordered), drop = FALSE]
    )
    l_libsize <- y_libsize[names(t_ordered)]
    raw_matrix <- as_matrix(
      cbind(
        data.frame(pseudotime = t_ordered),
        Matrix::t(y_ordered)
      )
    )

    log_message(
      "Calculating dynamic features for {.val {l}}...",
      verbose = verbose
    )
    if (fit_method == "gam") {
      out <- dynamic_features_gam(
        y_ordered = y_ordered,
        t_ordered = t_ordered,
        features = features,
        gene = gene,
        meta = meta,
        family = family,
        layer = layer,
        y_libsize = y_libsize,
        padjust_method = padjust_method,
        cores = cores,
        verbose = verbose
      )
    } else {
      out <- dynamic_features_pretsa(
        y_ordered = y_ordered,
        t_ordered = t_ordered,
        features = features,
        gene = gene,
        meta = meta,
        gene_status = gene_status,
        layer = layer,
        y_libsize = y_libsize,
        family = family,
        knot = knot,
        max_knot_allowed = max_knot_allowed,
        padjust_method = padjust_method,
        verbose = verbose
      )
    }
    res <- list(
      DynamicFeatures = out$DynamicFeatures,
      raw_matrix = raw_matrix,
      fitted_matrix = out$fitted_matrix,
      upr_matrix = out$upr_matrix,
      lwr_matrix = out$lwr_matrix,
      libsize = l_libsize,
      lineages = l,
      family = family
    )
    srt@tools[[paste0("DynamicFeatures_", suffix[i])]] <- res
  }

  log_message(
    "Find dynamic features done",
    message_type = "success",
    verbose = verbose
  )

  return(srt)
}

dynamic_features_gam <- function(
    y_ordered,
    t_ordered,
    features,
    gene,
    meta,
    family,
    layer,
    y_libsize,
    padjust_method,
    cores,
    verbose) {
  l_libsize <- y_libsize[names(t_ordered)]
  gam_out <- parallelize_fun(
    rownames(y_ordered),
    function(n) {
      family_current <- family[n]
      if (
        min(y_ordered[n, ]) < 0 &&
          family_current %in% c("nb", "poisson", "binomial")
      ) {
        log_message(
          "Negative values detected. Replace family with {.pkg gaussian} for the feature: {.val {n}}",
          message_type = "warning",
          verbose = verbose
        )
        family_use <- "gaussian"
      } else {
        family_use <- family_current
      }
      if (layer == "counts" && family_use != "gaussian" && !n %in% meta) {
        l_use <- l_libsize
      } else {
        l_use <- rep(stats::median(y_libsize), ncol(y_ordered))
      }
      sizefactror <- stats::median(y_libsize) / l_use
      mod <- mgcv::gam(
        y ~ s(x, bs = "cs") + offset(log(l_use)),
        family = family_use,
        data = data.frame(
          y = y_ordered[n, ],
          x = t_ordered,
          l_use = l_use
        )
      )
      pre <- stats::predict(mod, type = "link", se.fit = TRUE)
      upr <- pre$fit + (2 * pre$se.fit)
      lwr <- pre$fit - (2 * pre$se.fit)
      upr <- mod$family$linkinv(upr)
      lwr <- mod$family$linkinv(lwr)
      res <- summary(mod)
      fitted <- fitted(mod)
      pvalue <- res$s.table[[4]]
      dev_expl <- res$dev.expl
      r_sq <- max(0, min(1, res$r.sq))
      fitted.values <- fitted * sizefactror
      upr.values <- upr * sizefactror
      lwr.values <- lwr * sizefactror
      exp_ncells <- sum(
        y_ordered[n, ] > min(y_ordered[n, ]),
        na.rm = TRUE
      )
      peaktime <- stats::median(
        t_ordered[
          fitted.values > stats::quantile(fitted.values, 0.99, na.rm = TRUE)
        ]
      )
      valleytime <- stats::median(
        t_ordered[
          fitted.values < stats::quantile(fitted.values, 0.01, na.rm = TRUE)
        ]
      )
      list(
        features = n,
        exp_ncells = exp_ncells,
        r.sq = r_sq,
        dev.expl = dev_expl,
        peaktime = peaktime,
        valleytime = valleytime,
        pvalue = pvalue,
        fitted.values = fitted.values,
        upr.values = upr.values,
        lwr.values = lwr.values
      )
    },
    cores = cores,
    verbose = verbose
  )
  names(gam_out) <- sapply(gam_out, `[[`, "features")
  gam_out <- gam_out[rownames(y_ordered)]
  fitted_matrix <- do.call(
    cbind,
    lapply(gam_out, function(x) x[["fitted.values"]])
  )
  colnames(fitted_matrix) <- rownames(y_ordered)
  fitted_matrix <- cbind(pseudotime = t_ordered, fitted_matrix)
  upr_matrix <- do.call(
    cbind,
    lapply(gam_out, function(x) x[["upr.values"]])
  )
  colnames(upr_matrix) <- rownames(y_ordered)
  upr_matrix <- cbind(pseudotime = t_ordered, upr_matrix)
  lwr_matrix <- do.call(
    cbind,
    lapply(gam_out, function(x) x[["lwr.values"]])
  )
  colnames(lwr_matrix) <- rownames(y_ordered)
  lwr_matrix <- cbind(pseudotime = t_ordered, lwr_matrix)
  DynamicFeatures <- as.data.frame(
    do.call(
      rbind.data.frame,
      lapply(
        gam_out,
        function(x) {
          x[!names(x) %in% c("fitted.values", "upr.values", "lwr.values")]
        }
      )
    )
  )
  char_var <- c("features")
  numb_var <- colnames(DynamicFeatures)[!colnames(DynamicFeatures) %in% char_var]
  DynamicFeatures[, char_var] <- lapply(
    DynamicFeatures[, char_var, drop = FALSE], as.character
  )
  DynamicFeatures[, numb_var] <- lapply(
    DynamicFeatures[, numb_var, drop = FALSE], as.numeric
  )
  rownames(DynamicFeatures) <- DynamicFeatures[["features"]]
  DynamicFeatures[, "padjust"] <- stats::p.adjust(
    DynamicFeatures[, "pvalue", drop = TRUE],
    method = padjust_method
  )
  list(
    DynamicFeatures = DynamicFeatures,
    fitted_matrix = fitted_matrix,
    upr_matrix = upr_matrix,
    lwr_matrix = lwr_matrix
  )
}

dynamic_features_pretsa <- function(
    y_ordered,
    t_ordered,
    features,
    gene,
    meta,
    gene_status,
    layer,
    y_libsize,
    family,
    knot,
    max_knot_allowed,
    padjust_method,
    verbose) {
  gene_sub <- features[features %in% gene]
  meta_sub <- features[features %in% meta]
  pseudotime_vec <- stats::setNames(t_ordered, names(t_ordered))
  if (length(gene_sub) == 0) {
    log_message(
      "PreTSA: fitting {.val {length(meta_sub)}} meta feature(s) only",
      verbose = verbose
    )
    expr_meta <- as_matrix(y_ordered[meta_sub, , drop = FALSE])
    out <- pretsa_one_block(
      expr_meta,
      meta_sub,
      t_ordered,
      pseudotime_vec,
      knot,
      max_knot_allowed,
      padjust_method
    )
    DynamicFeatures <- out$DynamicFeatures[features, ]
    fit_mat <- out$fit_mat
    fitted_matrix <- cbind(
      pseudotime = t_ordered,
      as_matrix(Matrix::t(fit_mat)[, meta_sub, drop = FALSE])
    )
    upr_matrix <- lwr_matrix <- fitted_matrix
    return(
      list(
        DynamicFeatures = DynamicFeatures,
        fitted_matrix = fitted_matrix,
        upr_matrix = upr_matrix,
        lwr_matrix = lwr_matrix
      )
    )
  }
  expr_norm <- as_matrix(
    y_ordered[gene_sub, , drop = FALSE]
  )
  if (gene_status == "raw_counts" || layer == "counts") {
    expr_norm <- log1p(expr_norm)
  }
  if (any(expr_norm < 0, na.rm = TRUE)) {
    log_message(
      "Negative values detected in expression. PreTSA fits B-spline on the given values",
      message_type = "warning",
      verbose = verbose
    )
  }
  if (length(meta_sub) > 0) {
    log_message(
      "PreTSA: fitting {.val {length(gene_sub)}} gene(s) and {.val {length(meta_sub)}} meta feature(s)",
      verbose = verbose
    )
  }
  out_gene <- pretsa_one_block(
    expr_norm,
    gene_sub,
    t_ordered,
    pseudotime_vec,
    knot,
    max_knot_allowed,
    padjust_method
  )
  DynamicFeatures_gene <- out_gene$DynamicFeatures
  fit_mat <- out_gene$fit_mat
  fitted_matrix <- cbind(
    pseudotime = t_ordered,
    as_matrix(Matrix::t(fit_mat)[, gene_sub, drop = FALSE])
  )
  upr_matrix <- lwr_matrix <- fitted_matrix
  if (length(meta_sub) > 0) {
    expr_meta <- as_matrix(
      y_ordered[meta_sub, , drop = FALSE]
    )
    out_meta <- pretsa_one_block(
      expr_meta,
      meta_sub,
      t_ordered,
      pseudotime_vec,
      knot,
      max_knot_allowed,
      padjust_method
    )
    DynamicFeatures <- rbind(
      DynamicFeatures_gene,
      out_meta$DynamicFeatures
    )[features, ]
    fitted_matrix <- cbind(
      fitted_matrix,
      as_matrix(
        Matrix::t(
          out_meta$fit_mat
        )[, meta_sub, drop = FALSE]
      )
    )

    upr_matrix <- cbind(
      upr_matrix,
      as_matrix(
        Matrix::t(out_meta$fit_mat)[, meta_sub, drop = FALSE]
      )
    )
    lwr_matrix <- cbind(
      lwr_matrix,
      as_matrix(
        Matrix::t(out_meta$fit_mat)[, meta_sub, drop = FALSE]
      )
    )
  } else {
    DynamicFeatures <- DynamicFeatures_gene[features, ]
  }
  list(
    DynamicFeatures = DynamicFeatures,
    fitted_matrix = fitted_matrix,
    upr_matrix = upr_matrix,
    lwr_matrix = lwr_matrix
  )
}

pretsa_one_block <- function(
    expr,
    feats,
    t_ordered,
    pseudotime_vec,
    knot,
    max_knot_allowed,
    padjust_method) {
  test_res <- pretsa_temporal(
    expr,
    pseudotime_vec,
    knot = knot,
    max_knot_allowed = max_knot_allowed,
    padjust_method = padjust_method
  )
  fit_mat <- pretsa_temporalFit(
    expr,
    pseudotime_vec,
    knot = knot,
    max_knot_allowed = max_knot_allowed
  )
  SSE <- Matrix::rowSums((expr - fit_mat)^2)
  SST <- Matrix::rowSums(sweep(expr, 1, Matrix::rowMeans(expr), "-")^2)
  r_sq_vec <- 1 - SSE / SST
  r_sq_vec[!is.finite(r_sq_vec) | r_sq_vec < 0] <- 0
  peaktime_vec <- apply(fit_mat, 1, function(v) {
    stats::median(t_ordered[v >= stats::quantile(v, 0.99, na.rm = TRUE)])
  })
  valleytime_vec <- apply(fit_mat, 1, function(v) {
    stats::median(t_ordered[v <= stats::quantile(v, 0.01, na.rm = TRUE)])
  })
  exp_ncells_vec <- apply(expr, 1, function(v) sum(v > min(v), na.rm = TRUE))
  DF <- data.frame(
    features = feats,
    exp_ncells = exp_ncells_vec,
    r.sq = r_sq_vec,
    dev.expl = r_sq_vec,
    peaktime = peaktime_vec,
    valleytime = valleytime_vec,
    pvalue = test_res[feats, "pval"],
    stringsAsFactors = FALSE
  )
  DF[, "padjust"] <- stats::p.adjust(
    DF[, "pvalue"],
    method = padjust_method
  )
  rownames(DF) <- DF[, "features"]
  list(
    DynamicFeatures = DF,
    fit_mat = fit_mat
  )
}

pretsa_Calbic <- function(
    numknot,
    Blist,
    expr) {
  B <- Blist[[as.character(numknot)]][["B"]]
  tBB <- Blist[[as.character(numknot)]][["tBB"]]
  beta <- as.matrix(tcrossprod(chol2inv(chol(tBB)), B) %*% expr)
  pred <- B %*% beta
  mse <- Matrix::colMeans((expr - pred)^2)
  bic <- nrow(B) * (1 + log(2 * pi) + log(mse)) + log(nrow(B)) * (ncol(B) + 1)
  return(bic)
}

pretsa_Calfstat <- function(
    expr,
    pseudotime,
    knot = 0,
    max_knot_allowed = 10) {
  expr <- expr[, names(pseudotime), drop = FALSE]
  if (knot != "auto") {
    knotnum <- rep(knot, nrow(expr))
    names(knotnum) <- rownames(expr)
    B <- splines::bs(pseudotime, intercept = FALSE, df = knot + 3)
    B <- B[, which(apply(B, 2, stats::sd) > 0), drop = FALSE]
    B <- cbind(1, B)
    rownames(B) <- colnames(expr)
    colnames(B) <- NULL
    tBB <- crossprod(B)
    rownames(tBB) <- colnames(tBB) <- NULL
    pred <- as.matrix(expr %*% B %*% tcrossprod(chol2inv(chol(tBB)), B))
    SSE <- Matrix::rowSums((expr - pred)^2)
    SST <- Matrix::rowSums(
      sweep(expr, 1, Matrix::rowMeans(expr), FUN = "-")^2
    )
    fstat <- ((SST - SSE) / (ncol(B) - 1)) / (SSE / (nrow(B) - ncol(B)))
    fstat[which(Matrix::rowSums(expr) == 0)] <- 0
    return(fstat)
  } else {
    knotnum0 <- 0:max_knot_allowed
    names(knotnum0) <- knotnum0
    Blist <- lapply(knotnum0, function(numknot) {
      B <- splines::bs(pseudotime, intercept = FALSE, df = numknot + 3)
      B <- B[, which(apply(B, 2, stats::sd) > 0), drop = FALSE]
      B <- cbind(1, B)
      rownames(B) <- colnames(expr)
      colnames(B) <- NULL
      tBB <- crossprod(B)
      rownames(tBB) <- colnames(tBB) <- NULL
      list(B = B, tBB = tBB)
    })
    names(Blist) <- as.character(knotnum0)
    testpos <- sapply(knotnum0, function(numknot) {
      tBB <- Blist[[as.character(numknot)]][["tBB"]]
      !"try-error" %in% class(try(chol(tBB), silent = TRUE))
    })
    if (mean(testpos) != 1) {
      maxknot <- which(testpos == FALSE)[1] - 2
      knotnum0 <- 0:maxknot
      names(knotnum0) <- knotnum0
    }
    expr <- Matrix::t(expr)
    bic <- sapply(knotnum0, pretsa_Calbic, Blist = Blist, expr = expr)
    knotnum <- knotnum0[apply(bic, 1, which.min)]
    names(knotnum) <- rownames(bic)
    fstat <- lapply(unique(knotnum), function(k) {
      B <- Blist[[as.character(k)]][["B"]]
      tBB <- Blist[[as.character(k)]][["tBB"]]
      beta <- as.matrix(
        tcrossprod(
          chol2inv(chol(tBB)), B
        ) %*% expr[, which(knotnum == k), drop = FALSE]
      )
      pred <- B %*% beta
      expr.sub <- expr[, colnames(pred), drop = FALSE]
      SSE <- Matrix::colSums((expr.sub - pred)^2)
      SST <- Matrix::colSums(
        sweep(expr.sub, 2, Matrix::colMeans(expr.sub), FUN = "-")^2
      )
      fstat <- ((SST - SSE) / (ncol(B) - 1)) / (SSE / (nrow(B) - ncol(B)))
      fstat[which(Matrix::colSums(expr.sub) == 0)] <- 0
      fstat
    })
    fstat <- unlist(fstat)
    fstat <- fstat[colnames(expr)]
    return(fstat)
  }
}

pretsa_temporal <- function(
    expr,
    pseudotime,
    knot,
    max_knot_allowed,
    padjust_method) {
  expr <- expr[, names(pseudotime), drop = FALSE]
  if (knot != "auto") {
    knotnum <- rep(knot, nrow(expr))
    names(knotnum) <- rownames(expr)
    B <- splines::bs(pseudotime, intercept = FALSE, df = knot + 3)
    B <- B[, which(apply(B, 2, stats::sd) > 0), drop = FALSE]
    B <- cbind(1, B)
    rownames(B) <- colnames(expr)
    colnames(B) <- NULL
    tBB <- crossprod(B)
    rownames(tBB) <- colnames(tBB) <- NULL
    pred <- as.matrix(expr %*% B %*% tcrossprod(chol2inv(chol(tBB)), B))
    SSE <- Matrix::rowSums((expr - pred)^2)
    SST <- Matrix::rowSums(sweep(expr, 1, Matrix::rowMeans(expr), FUN = "-")^2)
    fstat <- ((SST - SSE) / (ncol(B) - 1)) / (SSE / (nrow(B) - ncol(B)))
    fstat[which(Matrix::rowSums(expr) == 0)] <- 0
    pval <- stats::pf(
      q = fstat,
      df1 = ncol(B) - 1,
      df2 = nrow(B) - ncol(B),
      lower.tail = FALSE
    )
    logpval <- stats::pf(
      q = fstat,
      df1 = ncol(B) - 1,
      df2 = nrow(B) - ncol(B),
      lower.tail = FALSE,
      log.p = TRUE
    )
    res <- data.frame(
      fstat = fstat,
      pval = pval,
      logpval = logpval
    )
  } else {
    knotnum0 <- 0:max_knot_allowed
    names(knotnum0) <- knotnum0
    Blist <- lapply(knotnum0, function(numknot) {
      B <- splines::bs(pseudotime, intercept = FALSE, df = numknot + 3)
      B <- B[, which(apply(B, 2, stats::sd) > 0), drop = FALSE]
      B <- cbind(1, B)
      rownames(B) <- colnames(expr)
      colnames(B) <- NULL
      tBB <- crossprod(B)
      rownames(tBB) <- colnames(tBB) <- NULL
      list(B = B, tBB = tBB)
    })
    names(Blist) <- as.character(knotnum0)
    testpos <- sapply(knotnum0, function(numknot) {
      tBB <- Blist[[as.character(numknot)]][["tBB"]]
      !"try-error" %in% class(try(chol(tBB), silent = TRUE))
    })
    if (mean(testpos) != 1) {
      maxknot <- which(testpos == FALSE)[1] - 2
      knotnum0 <- 0:maxknot
      names(knotnum0) <- knotnum0
    }
    expr <- Matrix::t(expr)
    bic <- sapply(knotnum0, pretsa_Calbic, Blist = Blist, expr = expr)
    knotnum <- knotnum0[apply(bic, 1, which.min)]
    names(knotnum) <- rownames(bic)
    res <- lapply(unique(knotnum), function(k) {
      B <- Blist[[as.character(k)]][["B"]]
      tBB <- Blist[[as.character(k)]][["tBB"]]
      beta <- as.matrix(
        tcrossprod(
          chol2inv(chol(tBB)), B
        ) %*% expr[, which(knotnum == k), drop = FALSE]
      )
      pred <- B %*% beta
      expr.sub <- expr[, colnames(pred), drop = FALSE]
      SSE <- Matrix::colSums((expr.sub - pred)^2)
      SST <- Matrix::colSums(
        sweep(
          expr.sub, 2, Matrix::colMeans(expr.sub),
          FUN = "-"
        )^2
      )
      fstat <- ((SST - SSE) / (ncol(B) - 1)) / (SSE / (nrow(B) - ncol(B)))
      fstat[which(Matrix::colSums(expr.sub) == 0)] <- 0
      pval <- stats::pf(
        q = fstat,
        df1 = ncol(B) - 1,
        df2 = nrow(B) - ncol(B),
        lower.tail = FALSE
      )
      logpval <- stats::pf(
        q = fstat,
        df1 = ncol(B) - 1,
        df2 = nrow(B) - ncol(B),
        lower.tail = FALSE,
        log.p = TRUE
      )
      data.frame(
        fstat = fstat,
        pval = pval,
        logpval = logpval
      )
    })
    res <- do.call(rbind, res)
    res <- res[colnames(expr), ]
  }
  res$fdr <- stats::p.adjust(res$pval, method = padjust_method)
  res$knotnum <- knotnum
  res <- res[, c("fdr", "logpval", "pval", "fstat", "knotnum")]
  res <- res[order(res$pval, -res$fstat), ]
  return(res)
}

pretsa_temporalFit <- function(
    expr,
    pseudotime,
    knot = 0,
    max_knot_allowed = 10) {
  expr <- expr[, names(pseudotime), drop = FALSE]
  if (knot != "auto") {
    B <- splines::bs(pseudotime, intercept = FALSE, df = knot + 3)
    B <- B[, which(apply(B, 2, stats::sd) > 0), drop = FALSE]
    B <- cbind(1, B)
    rownames(B) <- colnames(expr)
    colnames(B) <- NULL
    tBB <- crossprod(B)
    rownames(tBB) <- colnames(tBB) <- NULL
    invB <- chol2inv(chol(tBB))
    pred <- as.matrix(expr %*% B %*% invB %*% t(B))
    return(pred)
  } else {
    knotnum0 <- 0:max_knot_allowed
    names(knotnum0) <- knotnum0
    Blist <- lapply(knotnum0, function(numknot) {
      B <- splines::bs(pseudotime, intercept = FALSE, df = numknot + 3)
      B <- B[, which(apply(B, 2, stats::sd) > 0), drop = FALSE]
      B <- cbind(1, B)
      rownames(B) <- colnames(expr)
      colnames(B) <- NULL
      tBB <- crossprod(B)
      rownames(tBB) <- colnames(tBB) <- NULL
      list(B = B, tBB = tBB)
    })
    names(Blist) <- as.character(knotnum0)
    testpos <- sapply(knotnum0, function(numknot) {
      tBB <- Blist[[as.character(numknot)]][["tBB"]]
      !"try-error" %in% class(try(chol(tBB), silent = TRUE))
    })
    if (mean(testpos) != 1) {
      maxknot <- which(testpos == FALSE)[1] - 2
      knotnum0 <- 0:maxknot
      names(knotnum0) <- knotnum0
    }
    expr_t <- Matrix::t(expr)
    bic <- sapply(knotnum0, pretsa_Calbic, Blist = Blist, expr = expr_t)
    knotnum <- knotnum0[apply(bic, 1, which.min)]
    names(knotnum) <- rownames(bic)
    pred_list <- list()
    for (k in unique(knotnum)) {
      B <- Blist[[as.character(k)]][["B"]]
      tBB <- Blist[[as.character(k)]][["tBB"]]
      invB <- chol2inv(chol(tBB))
      idx <- which(knotnum == k)
      beta <- as.matrix(
        tcrossprod(invB, B) %*% expr_t[, idx, drop = FALSE]
      )
      pred_k <- B %*% beta
      pred_list[[as.character(k)]] <- pred_k
    }
    pred <- do.call(cbind, pred_list)
    pred <- pred[, colnames(expr_t)]
    pred <- Matrix::t(pred)
    return(pred)
  }
}
