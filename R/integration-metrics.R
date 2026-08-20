scale_ilisi <- function(mean_lisi, n_labels) {
  if (!is.finite(mean_lisi) || !is.finite(n_labels) || n_labels <= 1) {
    return(NA_real_)
  }
  max(0, min(1, (mean_lisi - 1) / (n_labels - 1)))
}

scale_clisi <- function(mean_lisi, n_labels) {
  if (!is.finite(mean_lisi) || !is.finite(n_labels) || n_labels <= 1) {
    return(NA_real_)
  }
  max(0, min(1, 1 - (mean_lisi - 1) / (n_labels - 1)))
}

scale_asw_bio <- function(asw) {
  if (!is.finite(asw)) {
    return(NA_real_)
  }
  max(0, min(1, (asw + 1) / 2))
}

n_label_levels <- function(srt, col) {
  if (is.null(col) || !col %in% colnames(srt@meta.data)) {
    return(NA_integer_)
  }
  length(unique(stats::na.omit(as.character(srt[[col, drop = TRUE]]))))
}

resolve_integration_latent_reduction <- function(
  srt,
  method,
  linear_reduction = "pca"
) {
  reduc <- SeuratObject::Reductions(srt)
  candidates <- unique(c(
    method,
    paste0(method, linear_reduction),
    paste0(method, "pca"),
    paste0(method, "lsi"),
    paste0(method, "svd")
  ))
  hit <- candidates[candidates %in% reduc]
  if (length(hit) > 0L) {
    return(hit[[1L]])
  }
  umap <- paste0(method, "UMAP2D")
  if (umap %in% reduc) {
    return(umap)
  }
  NA_character_
}

resolve_integration_umap_reduction <- function(srt, method) {
  reduc <- SeuratObject::Reductions(srt)
  candidates <- c(
    paste0(method, "UMAP2D"),
    paste0(method, "umap"),
    paste0(method, "UMAP")
  )
  hit <- candidates[candidates %in% reduc]
  if (length(hit) > 0L) {
    return(hit[[1L]])
  }
  NA_character_
}

resolve_integration_cluster_col <- function(
  srt,
  method,
  linear_reduction = "pca"
) {
  candidates <- c(
    paste0(method, "clusters"),
    paste0(method, linear_reduction, "clusters"),
    paste0(method, "pcaclusters")
  )
  hit <- candidates[candidates %in% colnames(srt@meta.data)]
  if (length(hit) > 0L) {
    return(hit[[1L]])
  }
  NA_character_
}

integration_metric_meta <- function() {
  data.frame(
    metric = c(
      "iLISI",
      "cLISI",
      "batch_ASW_mixing",
      "celltype_ASW",
      "celltype_graph_connectivity",
      "celltype_NMI",
      "celltype_ARI",
      "celltype_purity"
    ),
    category = c(
      "batch",
      "bio",
      "batch",
      "bio",
      "bio",
      "bio",
      "bio",
      "bio"
    ),
    stringsAsFactors = FALSE
  )
}

format_integration_metrics <- function(
  raw_df,
  srt,
  batch_col,
  celltype_col = NULL,
  method = NULL
) {
  if (is.null(raw_df) || nrow(raw_df) == 0L) {
    return(raw_df[integer(), , drop = FALSE])
  }
  out <- raw_df
  out$metric <- as.character(out$metric)
  out$value <- as.numeric(out$value)
  n_batch <- n_label_levels(srt, batch_col)
  n_type <- n_label_levels(srt, celltype_col)
  extra <- list()
  if (!is.null(batch_col)) {
    batch_lisi <- paste0(batch_col, "_LISI_mean")
    if (batch_lisi %in% out$metric) {
      extra[[length(extra) + 1L]] <- data.frame(
        metric = "iLISI",
        value = out$value[match(batch_lisi, out$metric)],
        stringsAsFactors = FALSE
      )
    }
  }
  if (!is.null(celltype_col)) {
    type_lisi <- paste0(celltype_col, "_LISI_mean")
    if (type_lisi %in% out$metric) {
      extra[[length(extra) + 1L]] <- data.frame(
        metric = "cLISI",
        value = out$value[match(type_lisi, out$metric)],
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(extra) > 0L) {
    out <- rbind(out, do.call(rbind, extra))
  }
  meta <- integration_metric_meta()
  out$category <- meta$category[match(out$metric, meta$metric)]
  out$category[is.na(out$category)] <- "other"
  out$direction <- "higher"
  out$scaled <- NA_real_
  i_ilisi <- out$metric == "iLISI"
  out$scaled[i_ilisi] <- vapply(
    out$value[i_ilisi],
    scale_ilisi,
    numeric(1),
    n_labels = n_batch
  )
  i_clisi <- out$metric == "cLISI"
  out$scaled[i_clisi] <- vapply(
    out$value[i_clisi],
    scale_clisi,
    numeric(1),
    n_labels = n_type
  )
  i_asw <- out$metric == "celltype_ASW"
  out$scaled[i_asw] <- vapply(out$value[i_asw], scale_asw_bio, numeric(1))
  already_unit <- out$metric %in% c(
    "batch_ASW_mixing",
    "celltype_graph_connectivity",
    "celltype_NMI",
    "celltype_ARI",
    "celltype_purity"
  )
  out$scaled[already_unit] <- pmax(0, pmin(1, out$value[already_unit]))
  other <- is.na(out$scaled)
  out$scaled[other] <- out$value[other]
  if (!is.null(method)) {
    out$method <- method
  }
  out[, c(
    intersect(
      c("method", "metric", "category", "value", "scaled", "direction"),
      colnames(out)
    )
  ), drop = FALSE]
}

integration_overall_scores <- function(
  metrics_df,
  bio_weight = 0.6,
  batch_weight = 0.4
) {
  w_sum <- bio_weight + batch_weight
  if (!is.finite(w_sum) || w_sum <= 0) {
    bio_weight <- 0.6
    batch_weight <- 0.4
    w_sum <- 1
  }
  bio_weight <- bio_weight / w_sum
  batch_weight <- batch_weight / w_sum
  methods <- unique(as.character(metrics_df$method))
  rows <- lapply(methods, function(method) {
    sub <- metrics_df[metrics_df$method == method, , drop = FALSE]
    core <- sub[sub$metric %in% integration_metric_meta()$metric, , drop = FALSE]
    bio <- mean(core$scaled[core$category == "bio"], na.rm = TRUE)
    batch <- mean(core$scaled[core$category == "batch"], na.rm = TRUE)
    if (!is.finite(bio)) {
      bio <- NA_real_
    }
    if (!is.finite(batch)) {
      batch <- NA_real_
    }
    overall <- if (is.finite(bio) && is.finite(batch)) {
      bio_weight * bio + batch_weight * batch
    } else if (is.finite(bio)) {
      bio
    } else {
      batch
    }
    data.frame(
      method = method,
      bio = bio,
      batch = batch,
      overall = overall,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

integration_summary_wide <- function(metrics_df, overall_df, runs_df = NULL) {
  methods <- unique(as.character(metrics_df$method))
  pick <- function(method, metric) {
    hit <- metrics_df$scaled[
      metrics_df$method == method & metrics_df$metric == metric
    ]
    if (length(hit) == 0L) {
      return(NA_real_)
    }
    hit[[1L]]
  }
  out <- overall_df
  for (metric in c(
    "iLISI",
    "cLISI",
    "batch_ASW_mixing",
    "celltype_ASW",
    "celltype_graph_connectivity",
    "celltype_NMI",
    "celltype_ARI"
  )) {
    out[[metric]] <- vapply(out$method, pick, numeric(1), metric = metric)
  }
  if (!is.null(runs_df) && nrow(runs_df) > 0L) {
    out <- merge(out, runs_df, by = "method", all.x = TRUE, sort = FALSE)
  }
  out[match(methods, out$method), , drop = FALSE]
}
