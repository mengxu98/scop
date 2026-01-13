#' @title RunDynamicFeatures
#'
#' @description
#' Calculates dynamic features for lineages in a single-cell RNA-seq dataset.
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
#' Features with a frequency less than minfreq will be excluded. Default is `5`.
#' @param family A character or character vector specifying the family of distributions to use for the generalized additive models (GAMs).
#' If family is set to NULL, the appropriate family will be automatically determined based on the data.
#' If length(family) is 1, the same family will be used for all features.
#' Otherwise, family must have the same length as features.
#' @param libsize A numeric or numeric vector specifying the library size correction factors for each cell.
#' If NULL, the library size correction factors will be calculated based on the expression matrix.
#' If length(libsize) is 1, the same value will be used for all cells.
#' Otherwise, libsize must have the same length as the number of cells in srt.
#' Default is `NULL`.
#'
#' @return Returns the modified Seurat object with the calculated dynamic features stored in the tools slot.
#'
#' @seealso
#' [DynamicHeatmap], [DynamicPlot], [RunDynamicEnrichment]
#'
#' @export
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
#'   n_candidates = 200
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
#'   n_split = 6,
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
    cores = 1,
    verbose = TRUE,
    seed = 11) {
  set.seed(seed)
  assay <- assay %||% DefaultAssay(srt)

  log_message(
    "Start find dynamic features",
    verbose = verbose
  )

  check_r("mgcv", verbose = FALSE)
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
    gam_out <- parallelize_fun(
      seq_len(nrow(y_ordered)),
      function(n) {
        feature_nm <- rownames(y_ordered)[n]
        family_current <- family[feature_nm]
        if (
          min(y_ordered[feature_nm, ]) < 0 &&
            family_current %in% c("nb", "poisson", "binomial")
        ) {
          log_message(
            "Negative values detected. Replace family with {.pkg gaussian} for the feature: {.val {feature_nm}}",
            message_type = "warning",
            verbose = verbose
          )
          family_use <- "gaussian"
        } else {
          family_use <- family_current
        }
        if (layer == "counts" && family_use != "gaussian" && !feature_nm %in% meta) {
          l_libsize <- l_libsize
        } else {
          l_libsize <- rep(stats::median(y_libsize), ncol(y_ordered))
        }
        sizefactror <- stats::median(y_libsize) / l_libsize
        mod <- mgcv::gam(
          y ~ s(x, bs = "cs") + offset(log(l_libsize)),
          family = family_use,
          data = data.frame(
            y = y_ordered[feature_nm, ],
            x = t_ordered,
            l_libsize = l_libsize
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
        r_sq <- res$r.sq
        fitted.values <- fitted * sizefactror
        upr.values <- upr * sizefactror
        lwr.values <- lwr * sizefactror
        exp_ncells <- sum(
          y_ordered[feature_nm, ] > min(y_ordered[feature_nm, ]),
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
          features = feature_nm,
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
    numb_var <- colnames(DynamicFeatures)[
      !colnames(DynamicFeatures) %in% char_var
    ]
    DynamicFeatures[, char_var] <- lapply(
      DynamicFeatures[, char_var, drop = FALSE],
      as.character
    )
    DynamicFeatures[, numb_var] <- lapply(
      DynamicFeatures[, numb_var, drop = FALSE],
      as.numeric
    )
    rownames(DynamicFeatures) <- DynamicFeatures[["features"]]
    DynamicFeatures[, "padjust"] <- stats::p.adjust(
      DynamicFeatures[, "pvalue", drop = TRUE]
    )

    res <- list(
      DynamicFeatures = DynamicFeatures,
      raw_matrix = raw_matrix,
      fitted_matrix = fitted_matrix,
      upr_matrix = upr_matrix,
      lwr_matrix = lwr_matrix,
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
