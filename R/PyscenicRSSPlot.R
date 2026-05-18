#' @title Plot top regulon specificity scores from pySCENIC results
#'
#' @description
#' Calculate regulon specificity score (RSS) from pySCENIC regulon activity and
#' plot the top regulons for each group.
#'
#' @md
#' @param srt A Seurat object containing pySCENIC results from
#' [RunPyscenic()].
#' @param group.by Metadata column used as the cell group annotation.
#' @param tool_name Name of the `srt@tools` entry storing pySCENIC results.
#' @param assay Assay used as a fallback source of regulon activity.
#' @param layer Assay layer used as a fallback source of regulon activity.
#' @param top_n Number of top regulons labeled for each group.
#' @param combine Whether to combine group plots with [patchwork::wrap_plots()].
#' @param ncol Number of columns used when `combine = TRUE`.
#' @param return_data Whether to return RSS matrices and ranking tables together
#' with plots. If `FALSE`, only the plot object or plot list is returned.
#' @param title Optional title added to the combined plot.
#' @param point_color Color for all regulon rank points.
#' @param top_color Color for top regulon rank points.
#' @param point_size Point size.
#' @param point_alpha Alpha for all regulon rank points.
#' @param highlight_tf Optional TF or regulon names to highlight in every group
#' plot. Values can match either `TF` or `regulon`, for example `"Sox9"` or
#' `"Sox9(+)"`.
#' @param highlight_color Color for highlighted TF or regulon points and rank
#' lines.
#' @param highlight_point_size Point size for highlighted TFs or regulons.
#' @param highlight_linewidth Line width for highlighted TF or regulon rank
#' lines.
#' @param label_size Text size for top regulon labels.
#' @param verbose Whether to print messages.
#'
#' @return A list containing `rss_matrix`, `rank_table`, `top_table`, `plots`,
#' and `plot` when `return_data = TRUE`; otherwise a plot object or list of
#' plots.
#' @export
PyscenicRSSPlot <- function(
  srt,
  group.by,
  tool_name = "Pyscenic",
  assay = "pyscenic",
  layer = "data",
  top_n = 12,
  combine = TRUE,
  ncol = 3,
  return_data = TRUE,
  title = NULL,
  point_color = "#1F77B4",
  top_color = "#DC050C",
  point_size = 3,
  point_alpha = 0.4,
  highlight_tf = NULL,
  highlight_color = "#7A0177",
  highlight_point_size = 4,
  highlight_linewidth = 0.4,
  label_size = 3,
  verbose = TRUE
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (length(group.by) != 1 || !group.by %in% colnames(srt@meta.data)) {
    log_message(
      "{.arg group.by} must be one metadata column in {.arg srt}",
      message_type = "error"
    )
  }
  if (length(top_n) != 1 || is.na(top_n) || top_n < 1) {
    log_message(
      "{.arg top_n} must be a positive integer",
      message_type = "error"
    )
  }
  top_n <- as.integer(top_n)
  if (!is.null(highlight_tf)) {
    highlight_tf <- unique(as.character(highlight_tf))
    highlight_tf <- highlight_tf[!is.na(highlight_tf) & nzchar(highlight_tf)]
    if (length(highlight_tf) == 0) {
      highlight_tf <- NULL
    }
  }

  auc_mat <- pyscenic_get_rss_auc_matrix(
    srt = srt,
    tool_name = tool_name,
    assay = assay,
    layer = layer
  )
  group_annotation <- srt@meta.data[[group.by]]
  names(group_annotation) <- rownames(srt@meta.data)

  common_cells <- intersect(colnames(auc_mat), names(group_annotation))
  if (length(common_cells) == 0) {
    log_message(
      "No shared cells between pySCENIC regulon activity and {.arg srt} metadata",
      message_type = "error"
    )
  }
  auc_mat <- auc_mat[, common_cells, drop = FALSE]
  group_annotation <- group_annotation[common_cells]

  keep_cells <- !is.na(group_annotation)
  if (!any(keep_cells)) {
    log_message(
      "All cells have missing values in {.arg group.by}",
      message_type = "error"
    )
  }
  auc_mat <- auc_mat[, keep_cells, drop = FALSE]
  group_annotation <- group_annotation[keep_cells]
  group_names <- if (is.factor(group_annotation)) {
    levels(droplevels(group_annotation))
  } else {
    unique(as.character(group_annotation))
  }
  group_annotation <- as.character(group_annotation)

  log_message(
    "Calculating pySCENIC RSS for {.val {nrow(auc_mat)}} regulons across {.val {length(group_names)}} group{?s}",
    verbose = verbose
  )
  rss_matrix <- pyscenic_calc_rss_matrix(
    auc_mat = auc_mat,
    cell_annotation = group_annotation,
    cell_types = group_names
  )
  rss_matrix <- rss_matrix[stats::complete.cases(rss_matrix), , drop = FALSE]
  if (nrow(rss_matrix) == 0) {
    log_message(
      "No valid RSS values remain after removing missing scores",
      message_type = "error"
    )
  }

  rank_table <- do.call(
    rbind,
    lapply(colnames(rss_matrix), function(one_group) {
      specificity_score <- rss_matrix[, one_group]
      keep_regulons <- !is.na(specificity_score)
      regulons <- rownames(rss_matrix)[keep_regulons]
      specificity_score <- specificity_score[keep_regulons]
      regulon_order <- order(specificity_score, decreasing = TRUE)
      regulons <- regulons[regulon_order]
      specificity_score <- specificity_score[regulon_order]
      data.frame(
        group = one_group,
        regulon = regulons,
        TF = sub("\\(\\+\\)$", "", regulons),
        specificity_score = as.numeric(specificity_score),
        rank = seq_along(regulons),
        is_top = seq_along(regulons) <= min(top_n, length(regulons)),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(rank_table) <- NULL
  rank_table[["is_highlight"]] <- FALSE
  if (!is.null(highlight_tf)) {
    rank_table[["is_highlight"]] <- rank_table[["regulon"]] %in% highlight_tf |
      rank_table[["TF"]] %in% highlight_tf
  }
  top_table <- rank_table[rank_table[["is_top"]], , drop = FALSE]

  plots <- lapply(colnames(rss_matrix), function(one_group) {
    data_rank_plot <- rank_table[rank_table[["group"]] == one_group, , drop = FALSE]
    top_df <- top_table[top_table[["group"]] == one_group, , drop = FALSE]
    highlight_df <- data_rank_plot[data_rank_plot[["is_highlight"]], , drop = FALSE]
    label_df <- unique(rbind(top_df, highlight_df))

    plot_title <- paste0("Group ", one_group)
    if (!is.null(highlight_tf)) {
      highlight_title <- if (nrow(highlight_df) > 0) {
        paste0(highlight_df[["TF"]], " rank = ", highlight_df[["rank"]])
      } else {
        paste0(highlight_tf, " not found")
      }
      plot_title <- paste0(
        "Group ",
        one_group,
        "\n",
        paste(highlight_title, collapse = "; ")
      )
    }

    p <- ggplot2::ggplot(
      data_rank_plot,
      ggplot2::aes(x = .data[["rank"]], y = .data[["specificity_score"]])
    ) +
      ggplot2::geom_point(
        size = point_size,
        shape = 16,
        color = point_color,
        alpha = point_alpha
      ) +
      ggplot2::geom_point(
        data = top_df,
        size = point_size,
        color = top_color
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.title = ggplot2::element_text(colour = "black", size = 12),
        axis.text = ggplot2::element_text(colour = "black", size = 10),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      ) +
      ggplot2::labs(
        x = "Regulon rank",
        y = "Specificity Score",
        title = plot_title
      )

    if (nrow(highlight_df) > 0) {
      p <- p +
        ggplot2::geom_point(
          data = highlight_df,
          size = highlight_point_size,
          color = highlight_color
        ) +
        ggplot2::geom_vline(
          data = highlight_df,
          ggplot2::aes(xintercept = .data[["rank"]]),
          linetype = "dashed",
          color = highlight_color,
          linewidth = highlight_linewidth
        )
    }

    p +
      ggrepel::geom_text_repel(
        data = label_df,
        ggplot2::aes(label = .data[["TF"]]),
        color = "black",
        size = label_size,
        fontface = "italic",
        arrow = grid::arrow(ends = "first", length = grid::unit(0.01, "npc")),
        box.padding = 0.2,
        point.padding = 0.3,
        segment.color = "black",
        segment.size = 0.3,
        force = 1,
        max.iter = 3000
      )
  })
  names(plots) <- colnames(rss_matrix)

  plot <- plots
  if (isTRUE(combine)) {
    check_r("patchwork", verbose = FALSE)
    plot <- if (length(plots) == 1) {
      plots[[1]]
    } else {
      patchwork::wrap_plots(plotlist = plots, ncol = ncol)
    }
    if (!is.null(title)) {
      plot <- plot + patchwork::plot_annotation(title = title)
    }
  }

  if (isFALSE(return_data)) {
    return(plot)
  }

  list(
    rss_matrix = rss_matrix,
    rank_table = rank_table,
    top_table = top_table,
    plots = plots,
    plot = plot
  )
}

pyscenic_get_rss_auc_matrix <- function(
  srt,
  tool_name = "Pyscenic",
  assay = "pyscenic",
  layer = "data"
) {
  if (!is.null(srt@tools[[tool_name]][["scores_cells_by_regulon"]])) {
    scores_cells_by_regulon <- srt@tools[[tool_name]][["scores_cells_by_regulon"]]
    auc_mat <- as.matrix(scores_cells_by_regulon)
    if (is.null(rownames(auc_mat)) || is.null(colnames(auc_mat))) {
      log_message(
        "{.arg scores_cells_by_regulon} must have cell row names and regulon column names",
        message_type = "error"
      )
    }
    return(t(auc_mat))
  }

  if (!assay %in% SeuratObject::Assays(srt)) {
    log_message(
      "Cannot find pySCENIC results in tools slot {.val {tool_name}} or assay {.val {assay}}",
      message_type = "error"
    )
  }
  auc_mat <- GetAssayData5(srt, assay = assay, layer = layer)
  auc_mat <- as.matrix(auc_mat)
  if (is.null(rownames(auc_mat)) || is.null(colnames(auc_mat))) {
    log_message(
      "{.arg assay} regulon activity matrix must have regulon row names and cell column names",
      message_type = "error"
    )
  }
  auc_mat
}

pyscenic_calc_rss_matrix <- function(
  auc_mat,
  cell_annotation,
  cell_types = NULL
) {
  if (any(is.na(cell_annotation))) {
    log_message(
      "{.arg cell_annotation} contains missing values",
      message_type = "error"
    )
  }
  if (is.null(cell_types)) {
    cell_types <- unique(cell_annotation)
  }

  row_sums <- rowSums(auc_mat)
  norm_auc <- auc_mat / row_sums
  rss_list <- lapply(cell_types, function(this_type) {
    p_cell_type <- as.numeric(cell_annotation == this_type)
    p_cell_type <- p_cell_type / sum(p_cell_type)
    vapply(
      seq_len(nrow(norm_auc)),
      function(regulon_idx) {
        pyscenic_calc_one_rss(norm_auc[regulon_idx, ], p_cell_type)
      },
      numeric(1)
    )
  })
  rss <- do.call(cbind, rss_list)
  rownames(rss) <- rownames(auc_mat)
  colnames(rss) <- cell_types
  rss
}

pyscenic_calc_one_rss <- function(p_regulon, p_cell_type) {
  jsd <- pyscenic_calc_jsd(p_regulon, p_cell_type)
  1 - sqrt(jsd)
}

pyscenic_calc_jsd <- function(p_regulon, p_cell_type) {
  pyscenic_entropy((p_regulon + p_cell_type) / 2) -
    ((pyscenic_entropy(p_regulon) + pyscenic_entropy(p_cell_type)) / 2)
}

pyscenic_entropy <- function(p_vector) {
  p_vector <- p_vector[p_vector > 0]
  -sum(p_vector * log2(p_vector))
}
