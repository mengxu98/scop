#' @title Plot CellChat analysis results
#'
#' @description
#' [CellChatPlot] creates various visualizations for CellChat analysis results
#' stored in a Seurat object. It supports both single-condition CellChat objects
#' and merged comparison objects produced by [RunCellChat].
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param srt A Seurat object processed with [RunCellChat].
#' @param plot_type Type of plot to create.
#' @param condition Condition or comparison name to plot. For single-condition
#'   plotting this should match one entry in
#'   `CellChat results`. For comparison plotting this
#'   should match one entry in `CellChat comparisons`.
#' @param pathway Backward-compatible alias of `signaling`.
#' @param signaling Signaling pathway(s) to visualize.
#' @param pairLR.use Ligand-receptor pair(s) to visualize.
#' @param sources.use Source cell groups to include. Can be group names or indices.
#' @param targets.use Target cell groups to include. Can be group names or indices.
#' @param comparison Dataset indices or names for merged comparison plots.
#' @param dataset Dataset index or name to use when `condition` refers to a
#'   merged comparison but a single dataset-level plot is requested.
#' @param measure Communication summary to use for network-level plots.
#' @param pattern Pattern used in signaling role heatmaps.
#' @param idents.use Cell groups of interest for `plot_type = "role_change"`.
#' @param signaling.exclude Signaling pathways to exclude for
#'   `plot_type = "role_change"`.
#' @param slot.name Slot name passed to [CellChat::subsetCommunication].
#' @param thresh P-value threshold used when extracting or plotting
#'   communications.
#' @param return.data Whether to return the underlying communication table
#'   instead of a plot when supported.
#' @param remove.isolate Whether to remove isolated source-target pairs in bubble
#'   plots.
#' @param bubble_color.by Which statistic to map to bubble color in the robust
#'   table-based bubble plot.
#' @param bubble_size.range Size range of points in the robust bubble plot.
#' @param palette Color palette name. Available palettes can be found in
#'   `thisplot::show_palettes()`.
#' @param palcolor Custom colors used to create a color palette.
#' @param font.size Base text size passed to CellChat plotting functions when
#'   supported.
#' @param label.size Label size used in plots that explicitly expose label size.
#' @param title Plot title. If `NULL`, sensible defaults are used.
#' @param subtitle Plot subtitle.
#' @param title.name Backward-compatible alias of `title` for CellChat native
#'   plotting functions.
#' @param xlab X-axis label for ggplot-based outputs.
#' @param ylab Y-axis label for ggplot-based outputs.
#' @param angle.x Angle of x-axis text in custom bubble plots.
#' @param hjust.x Horizontal justification of x-axis text in custom bubble plots.
#' @param vjust.x Vertical justification of x-axis text in custom bubble plots.
#' @param aspect.ratio Aspect ratio for ggplot-based outputs.
#' @param legend.position Legend position.
#' @param legend.direction Legend direction.
#' @param legend.title Legend title.
#' @param theme_use Theme used for ggplot-based outputs. Can be a character
#'   string or a theme function.
#' @param theme_args Other arguments passed to `theme_use`.
#' @param combine Combine multiple ggplot outputs with `patchwork::wrap_plots`.
#' @param nrow Number of rows in combined ggplot output.
#' @param ncol Number of columns in combined ggplot output.
#' @param byrow Whether combined ggplot output should be filled by row.
#' @param layout Layout used for pathway or individual network plots.
#' @param top_n Number of top pathways to use when `signaling = NULL`.
#' @param reduce Whether to simplify the chord diagram by keeping only top links.
#' @param max.groups Maximum number of groups to retain in reduced chord plots.
#' @param small.gap Small gap between sectors in chord plots.
#' @param big.gap Large gap between sectors in chord plots.
#' @param lab.cex Label size in chord plots.
#' @param dirpath Directory to save plots.
#' @param output_format Output format. One of `"pdf"` or `"png"`.
#' @param base_height Height multiplier when auto-saving plots.
#' @param base_width Width multiplier when auto-saving plots.
#' @param res Resolution for PNG output.
#' @param ... Additional arguments passed to the underlying CellChat/Seurat
#'   plotting function for the chosen `plot_type`.
#'
#' @export
#'
#' @seealso
#' [RunCellChat]
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#' pancreas_sub <- RunCellChat(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   species = "Mus_musculus"
#' )
#'
#' CellChatPlot(pancreas_sub, plot_type = "aggregate")
#' CellChatPlot(pancreas_sub, plot_type = "bubble")
CellChatPlot <- function(
    srt,
    plot_type = c(
      "aggregate",
      "overview",
      "pathway",
      "comparison",
      "heatmap",
      "circle",
      "bubble",
      "gene",
      "diff",
      "ranknet",
      "chord",
      "lr_contribution",
      "individual_lr",
      "role_scatter",
      "role_heatmap",
      "role_change"
    ),
    condition = NULL,
    pathway = NULL,
    signaling = NULL,
    pairLR.use = NULL,
    sources.use = NULL,
    targets.use = NULL,
    comparison = c(1, 2),
    dataset = 1,
    measure = c("count", "weight"),
    pattern = c("outgoing", "incoming", "all"),
    idents.use = NULL,
    signaling.exclude = NULL,
    slot.name = "net",
    thresh = 0.05,
    return.data = FALSE,
    remove.isolate = TRUE,
    bubble_color.by = c("prob", "pval"),
    bubble_size.range = c(1.5, 8),
    palette = "Paired",
    palcolor = NULL,
    font.size = 10,
    label.size = 1,
    title = NULL,
    subtitle = NULL,
    title.name = NULL,
    xlab = NULL,
    ylab = NULL,
    angle.x = 45,
    hjust.x = 1,
    vjust.x = 1,
    aspect.ratio = NULL,
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = NULL,
    theme_use = "theme_scop",
    theme_args = list(),
    combine = TRUE,
    nrow = NULL,
    ncol = NULL,
    byrow = TRUE,
    layout = c("circle", "hierarchy", "chord"),
    top_n = 10,
    reduce = TRUE,
    max.groups = 8,
    small.gap = 1,
    big.gap = 8,
    lab.cex = 0.6,
    dirpath = NULL,
    output_format = c("pdf", "png"),
    base_height = 1,
    base_width = 1,
    res = 300,
    verbose = TRUE,
    ...
) {
  plot_type <- match.arg(plot_type)
  if (identical(plot_type, "overview")) {
    plot_type <- "aggregate"
  }
  measure <- match.arg(measure)
  pattern <- match.arg(pattern)
  bubble_color.by <- match.arg(bubble_color.by)
  layout <- match.arg(layout)
  output_format <- match.arg(output_format)

  store <- .cc_get_store(srt)
  if (is.null(store)) {
    log_message(
      "No CellChat results found. Please run {.fn RunCellChat} first",
      message_type = "error"
    )
  }

  signaling <- unique(stats::na.omit(c(signaling, pathway)))
  if (length(signaling) == 0L) {
    signaling <- NULL
  }

  if (!is.null(dirpath) && !dir.exists(dirpath)) {
    dir.create(dirpath, recursive = TRUE)
  }

  title.name <- title.name %||% title

  if (plot_type %in% c("aggregate", "pathway", "circle", "gene", "lr_contribution", "individual_lr", "chord", "role_scatter", "role_heatmap")) {
    obj_info <- .cc_get_dataset_object(
      srt = srt,
      condition = condition,
      dataset = dataset
    )
    cellchat_object <- obj_info$object
    plot_label <- obj_info$label
  }

  if (plot_type %in% c("comparison", "diff", "ranknet", "role_change")) {
    cmp <- .cc_get_cmp(srt = srt, condition = condition)
    merged_object <- cmp$comparison_object
    obj.list <- cmp$object.list
    comp_idx <- .cc_resolve_dataset_index(cmp, comparison = comparison)
  }

  if (plot_type == "aggregate") {
    p <- .cc_plot_aggregate(
      cellchat_object = cellchat_object,
      condition = plot_label,
      title = title.name,
      font.size = font.size,
      label.size = label.size,
      theme_use = theme_use,
      theme_args = theme_args,
      legend.position = legend.position,
      legend.direction = legend.direction,
      legend.title = legend.title,
      aspect.ratio = aspect.ratio,
      xlab = xlab,
      ylab = ylab,
      subtitle = subtitle,
      dirpath = dirpath,
      output_format = output_format,
      base_height = base_height,
      base_width = base_width,
      res = res,
      verbose = verbose,
      ...
    )
    log_message("Plot creation completed", message_type = "success", verbose = verbose)
    return(p)
  }

  if (plot_type == "pathway") {
    pathways_to_show <- signaling %||% cellchat_object@netP$pathways[seq_len(min(top_n, length(cellchat_object@netP$pathways)))]
    pathways_to_show <- pathways_to_show[!is.na(pathways_to_show)]
    if (length(pathways_to_show) == 0L) {
      log_message(
        "No signaling pathways available for plotting",
        message_type = "error"
      )
    }

    plots <- lapply(pathways_to_show, function(sig) {
      p <- .cc_capture_drawn_plot(CellChat::netVisual_aggregate(
        cellchat_object,
        signaling = sig,
        layout = layout,
        sources.use = .cc_resolve_group_index_single(cellchat_object, sources.use),
        targets.use = .cc_resolve_group_index_single(cellchat_object, targets.use),
        thresh = thresh,
        pt.title = font.size,
        vertex.label.cex = max(0.8, label.size),
        small.gap = small.gap,
        big.gap = big.gap,
        ...
      ))
      if (!is.null(dirpath)) {
        file <- file.path(dirpath, paste0(sig, "_", plot_label, "_pathway.", output_format))
        .cc_replay_recorded_plot(
          plot_obj = p,
          file = file,
          output_format = output_format,
          height = 5 * base_height,
          width = 5 * base_width,
          res = res
        )
      }
      p
    })

    plots <- .cc_simplify_plot_list(plots)
    log_message("Plot creation completed", message_type = "success", verbose = verbose)
    return(plots)
  }

  if (plot_type == "lr_contribution") {
    if (is.null(signaling)) {
      log_message(
        "{.arg signaling} must be provided for {.val plot_type = 'lr_contribution'}",
        message_type = "error"
      )
    }
    plots <- lapply(signaling, function(sig) {
      p <- CellChat::netAnalysis_contribution(
        cellchat_object,
        signaling = sig,
        title = title.name %||% sig,
        ...
      )
      p <- .cc_finalize_ggplot(
        p,
        title = title.name %||% sig,
        subtitle = subtitle,
        xlab = xlab,
        ylab = ylab,
        legend.position = legend.position,
        legend.direction = legend.direction,
        legend.title = legend.title,
        theme_use = theme_use,
        theme_args = theme_args,
        aspect.ratio = aspect.ratio,
        font.size = font.size
      )
      if (!is.null(dirpath)) {
        file <- file.path(dirpath, paste0(sig, "_", plot_label, "_lr_contribution.", output_format))
        .cc_safe_save_ggplot(p, file, width = 6 * base_width, height = 5 * base_height, res = res)
      }
      p
    })
    plots <- .cc_finalize_plot_list(plots, combine = combine, nrow = nrow, ncol = ncol, byrow = byrow)
    log_message("Plot creation completed", message_type = "success", verbose = verbose)
    return(plots)
  }

  if (plot_type == "individual_lr") {
    if (is.null(signaling)) {
      log_message(
        "{.arg signaling} must be provided for {.val plot_type = 'individual_lr'}",
        message_type = "error"
      )
    }
    if (is.null(pairLR.use)) {
      log_message(
        "{.arg pairLR.use} must be provided for {.val plot_type = 'individual_lr'}",
        message_type = "error"
      )
    }

    plots <- lapply(signaling, function(sig) {
      p <- .cc_capture_drawn_plot(CellChat::netVisual_individual(
        cellchat_object,
        signaling = sig,
        pairLR.use = pairLR.use,
        layout = layout,
        sources.use = .cc_resolve_group_index_single(cellchat_object, sources.use),
        targets.use = .cc_resolve_group_index_single(cellchat_object, targets.use),
        thresh = thresh,
        ...
      ))
      if (!is.null(dirpath)) {
        file <- file.path(dirpath, paste0(sig, "_", plot_label, "_individual_lr.", output_format))
        .cc_replay_recorded_plot(
          plot_obj = p,
          file = file,
          output_format = output_format,
          height = 5 * base_height,
          width = 5 * base_width,
          res = res
        )
      }
      p
    })
    plots <- .cc_simplify_plot_list(plots)
    log_message("Plot creation completed", message_type = "success", verbose = verbose)
    return(plots)
  }

  if (plot_type == "chord") {
    if (is.null(signaling)) {
      log_message(
        "{.arg signaling} must be provided for {.val plot_type = 'chord'}",
        message_type = "error"
      )
    }
    plots <- lapply(signaling, function(sig) {
      p <- .cc_chord_plot(
        cellchat_object = cellchat_object,
        signaling = sig,
        sources.use = sources.use,
        targets.use = targets.use,
        thresh = thresh,
        reduce = reduce,
        top_n = top_n,
        max.groups = max.groups,
        small.gap = small.gap,
        big.gap = big.gap,
        lab.cex = lab.cex,
        title.name = title.name %||% sig,
        dirpath = dirpath,
        output_format = output_format,
        fileprefix = paste0(sig, "_", plot_label, "_chord"),
        base_height = base_height,
        base_width = base_width,
        res = res,
        ...
      )
      p
    })
    plots <- .cc_simplify_plot_list(plots)
    log_message("Plot creation completed", message_type = "success", verbose = verbose)
    return(plots)
  }

  if (plot_type == "circle") {
    p <- .cc_plot_circle(
      cellchat_object = cellchat_object,
      condition = plot_label,
      measure = measure,
      label.size = label.size,
      dirpath = dirpath,
      output_format = output_format,
      base_height = base_height,
      base_width = base_width,
      res = res,
      ...
    )
    log_message("Plot creation completed", message_type = "success", verbose = verbose)
    return(p)
  }

  if (plot_type == "gene") {
    pathways_to_show <- signaling %||% cellchat_object@netP$pathways[seq_len(min(3, top_n, length(cellchat_object@netP$pathways)))]
    pathways_to_show <- pathways_to_show[!is.na(pathways_to_show)]
    if (length(pathways_to_show) == 0L) {
      log_message(
        "No signaling pathways available for gene plotting",
        message_type = "error"
      )
    }

    parameters <- .cc_get_store(srt)$parameters %||% list()
    plots <- lapply(pathways_to_show, function(sig) {
      p <- .cc_create_gene_plot(
        cellchat_object = cellchat_object,
        seurat_object = obj_info$seurat_object,
        signaling = sig,
        species = parameters$species,
        title = title.name %||% sig,
        theme_use = theme_use,
        theme_args = theme_args,
        legend.position = legend.position,
        legend.direction = legend.direction,
        legend.title = legend.title,
        subtitle = subtitle,
        xlab = xlab,
        ylab = ylab,
        aspect.ratio = aspect.ratio,
        ...
      )
      if (!is.null(dirpath)) {
        file <- file.path(dirpath, paste0(sig, "_", plot_label, "_gene.", output_format))
        .cc_safe_save_ggplot(p, file, width = 7 * base_width, height = 5 * base_height, res = res)
      }
      p
    })
    plots <- .cc_finalize_plot_list(plots, combine = combine, nrow = nrow, ncol = ncol, byrow = byrow)
    log_message("Plot creation completed", message_type = "success", verbose = verbose)
    return(plots)
  }

  if (plot_type == "bubble") {
    if (.cc_should_use_single_condition(srt, condition = condition)) {
      cond1 <- .cc_resolve_single_condition(srt, condition = condition)
      obj1 <- .cc_get_single_obj(srt, condition = cond1)

      df <- .cc_prepare_bubble_data(
        object = obj1,
        slot.name = slot.name,
        signaling = signaling,
        pairLR.use = pairLR.use,
        sources.use = .cc_resolve_group_index_single(obj1, sources.use),
        targets.use = .cc_resolve_group_index_single(obj1, targets.use),
        thresh = thresh,
        dataset = NULL,
        top_n = top_n,
        remove.isolate = remove.isolate
      )

      if (isTRUE(return.data)) {
        return(df)
      }

      p <- .cc_custom_bubble_plot(
        df = df,
        title = title.name %||% paste0(cond1, " bubble"),
        subtitle = subtitle,
        xlab = xlab,
        ylab = ylab,
        color.by = bubble_color.by,
        bubble_size.range = bubble_size.range,
        palette = palette,
        palcolor = palcolor,
        font.size = font.size,
        angle.x = angle.x,
        hjust.x = hjust.x,
        vjust.x = vjust.x,
        legend.position = legend.position,
        legend.direction = legend.direction,
        legend.title = legend.title,
        theme_use = theme_use,
        theme_args = theme_args,
        aspect.ratio = aspect.ratio,
        remove.isolate = remove.isolate
      )

      if (!is.null(dirpath)) {
        dims <- .cc_bubble_dims(df, base_height = base_height, base_width = base_width)
        file <- file.path(dirpath, paste0(cond1, "_bubble.", output_format))
        .cc_safe_save_ggplot(p, file, width = dims$width, height = dims$height, res = res)
      }
      log_message("Plot creation completed", message_type = "success", verbose = verbose)
      return(p)
    }

    cmp <- .cc_get_cmp(srt = srt, condition = condition)
    obj.list <- cmp$object.list
    comp_idx <- .cc_resolve_dataset_index(cmp, comparison = comparison)

    df_list <- lapply(names(obj.list)[comp_idx], function(ds) {
      obj <- obj.list[[ds]]
      .cc_prepare_bubble_data(
        object = obj,
        slot.name = slot.name,
        signaling = signaling,
        pairLR.use = pairLR.use,
        sources.use = .cc_resolve_group_index_single(obj, sources.use),
        targets.use = .cc_resolve_group_index_single(obj, targets.use),
        thresh = thresh,
        dataset = ds,
        top_n = top_n,
        remove.isolate = remove.isolate
      )
    })
    names(df_list) <- names(obj.list)[comp_idx]

    if (isTRUE(return.data)) {
      return(df_list)
    }

    df_list <- Filter(function(x) !is.null(x) && nrow(x) > 0L, df_list)
    if (length(df_list) == 0L) {
      log_message(
        "No communication records found under the current filtering conditions",
        message_type = "error"
      )
    }
    df <- do.call(rbind, df_list)

    p <- .cc_custom_bubble_plot(
      df = df,
      title = title.name %||% .cc_cmp_label(cmp, comp_idx),
      subtitle = subtitle,
      xlab = xlab,
      ylab = ylab,
      color.by = bubble_color.by,
      bubble_size.range = bubble_size.range,
      palette = palette,
      palcolor = palcolor,
      font.size = font.size,
      angle.x = angle.x,
      hjust.x = hjust.x,
      vjust.x = vjust.x,
      legend.position = legend.position,
      legend.direction = legend.direction,
      legend.title = legend.title,
      theme_use = theme_use,
      theme_args = theme_args,
      aspect.ratio = aspect.ratio,
      remove.isolate = remove.isolate
    )

    if (!is.null(dirpath)) {
      dims <- .cc_bubble_dims(df, base_height = base_height, base_width = base_width)
      file <- file.path(dirpath, paste0(.cc_cmp_label(cmp, comp_idx), "_bubble.", output_format))
      .cc_safe_save_ggplot(p, file, width = dims$width, height = dims$height, res = res)
    }
    log_message("Plot creation completed", message_type = "success", verbose = verbose)
    return(p)
  }

  if (plot_type == "comparison") {
    p <- CellChatPlot(
      srt = srt,
      plot_type = "bubble",
      condition = condition,
      pathway = pathway,
      signaling = signaling,
      pairLR.use = pairLR.use,
      sources.use = sources.use,
      targets.use = targets.use,
      comparison = comparison,
      dataset = dataset,
      measure = measure,
      pattern = pattern,
      idents.use = idents.use,
      signaling.exclude = signaling.exclude,
      slot.name = slot.name,
      thresh = thresh,
      return.data = return.data,
      remove.isolate = remove.isolate,
      bubble_color.by = bubble_color.by,
      bubble_size.range = bubble_size.range,
      palette = palette,
      palcolor = palcolor,
      font.size = font.size,
      label.size = label.size,
      title = title,
      subtitle = subtitle,
      title.name = title.name,
      xlab = xlab,
      ylab = ylab,
      angle.x = angle.x,
      hjust.x = hjust.x,
      vjust.x = vjust.x,
      aspect.ratio = aspect.ratio,
      legend.position = legend.position,
      legend.direction = legend.direction,
      legend.title = legend.title,
      theme_use = theme_use,
      theme_args = theme_args,
      combine = combine,
      nrow = nrow,
      ncol = ncol,
      byrow = byrow,
      layout = layout,
      top_n = top_n,
      reduce = reduce,
      max.groups = max.groups,
      small.gap = small.gap,
      big.gap = big.gap,
      lab.cex = lab.cex,
      dirpath = dirpath,
      output_format = output_format,
      base_height = base_height,
      base_width = base_width,
      res = res,
      verbose = verbose,
      ...
    )
    return(p)
  }

  if (plot_type == "diff") {
    p <- .cc_capture_drawn_plot(CellChat::netVisual_diffInteraction(
      merged_object,
      measure = measure,
      comparison = comp_idx,
      ...
    ))
    if (!is.null(dirpath)) {
      file <- file.path(dirpath, paste0(.cc_cmp_label(cmp, comp_idx), "_diff.", output_format))
      .cc_replay_recorded_plot(
        plot_obj = p,
        file = file,
        output_format = output_format,
        height = 5 * base_height,
        width = 6 * base_width,
        res = res
      )
    }
    log_message("Plot creation completed", message_type = "success", verbose = verbose)
    return(p)
  }

  if (plot_type == "heatmap") {
    if (.cc_should_use_single_condition(srt, condition = condition)) {
      cond1 <- .cc_resolve_single_condition(srt, condition = condition)
      obj1 <- .cc_get_single_obj(srt, condition = cond1)
      p <- CellChat::netVisual_heatmap(
        obj1,
        measure = measure,
        signaling = signaling,
        slot.name = slot.name,
        sources.use = .cc_resolve_group_index_single(obj1, sources.use),
        targets.use = .cc_resolve_group_index_single(obj1, targets.use),
        title.name = title.name %||% cond1,
        font.size = font.size,
        ...
      )
      if (!is.null(dirpath)) {
        file <- file.path(dirpath, paste0(cond1, "_heatmap.", output_format))
        .cc_draw_heatmap_to_file(
          heatmap_obj = p,
          file = file,
          output_format = output_format,
          height = 5.5 * base_height,
          width = 5.5 * base_width,
          res = res
        )
      }
      log_message("Plot creation completed", message_type = "success", verbose = verbose)
      return(p)
    }

    cmp <- .cc_get_cmp(srt = srt, condition = condition)
    merged_object <- cmp$comparison_object
    comp_idx <- .cc_resolve_dataset_index(cmp, comparison = comparison)

    p <- CellChat::netVisual_heatmap(
      merged_object,
      comparison = comp_idx,
      measure = measure,
      signaling = signaling,
      slot.name = slot.name,
      sources.use = .cc_resolve_group_index_merged(merged_object, sources.use),
      targets.use = .cc_resolve_group_index_merged(merged_object, targets.use),
      title.name = title.name %||% .cc_cmp_label(cmp, comp_idx),
      font.size = font.size,
      ...
    )
    if (!is.null(dirpath)) {
      file <- file.path(dirpath, paste0(.cc_cmp_label(cmp, comp_idx), "_heatmap.", output_format))
      .cc_draw_heatmap_to_file(
        heatmap_obj = p,
        file = file,
        output_format = output_format,
        height = 5.5 * base_height,
        width = 5.5 * base_width,
        res = res
      )
    }
    log_message("Plot creation completed", message_type = "success", verbose = verbose)
    return(p)
  }

  if (plot_type == "ranknet") {
    rank_colors <- .cc_get_palette_colors(
      palette = palette,
      palcolor = palcolor,
      n = max(2L, length(comp_idx))
    )
    p <- CellChat::rankNet(
      merged_object,
      mode = "comparison",
      comparison = comp_idx,
      color.use = rank_colors[seq_len(length(comp_idx))],
      stacked = TRUE,
      do.stat = TRUE,
      sources.use = .cc_resolve_group_index_merged(merged_object, sources.use),
      targets.use = .cc_resolve_group_index_merged(merged_object, targets.use),
      return.data = return.data,
      title = title.name %||% .cc_cmp_label(cmp, comp_idx),
      font.size = font.size,
      thresh = thresh,
      ...
    )
    if (isTRUE(return.data)) {
      return(p)
    }
    p <- .cc_finalize_patchwork(
      p,
      title = title.name %||% .cc_cmp_label(cmp, comp_idx),
      subtitle = subtitle,
      legend.position = legend.position,
      legend.direction = legend.direction,
      theme_use = theme_use,
      theme_args = theme_args,
      font.size = font.size
    )
    if (!is.null(dirpath)) {
      file <- file.path(dirpath, paste0(.cc_cmp_label(cmp, comp_idx), "_ranknet.", output_format))
      .cc_safe_save_ggplot(p, file, width = 8 * base_width, height = 6 * base_height, res = res)
    }
    log_message("Plot creation completed", message_type = "success", verbose = verbose)
    return(p)
  }

  if (plot_type == "role_scatter") {
    p <- CellChat::netAnalysis_signalingRole_scatter(
      cellchat_object,
      signaling = signaling,
      label.size = max(2.5, label.size),
      title = title.name %||% plot_label,
      font.size = font.size,
      ...
    )
    p <- .cc_finalize_ggplot(
      p,
      title = title.name %||% plot_label,
      subtitle = subtitle,
      xlab = xlab,
      ylab = ylab,
      legend.position = legend.position,
      legend.direction = legend.direction,
      legend.title = legend.title,
      theme_use = theme_use,
      theme_args = theme_args,
      aspect.ratio = aspect.ratio,
      font.size = font.size
    )
    if (!is.null(dirpath)) {
      file <- file.path(dirpath, paste0(plot_label, "_role_scatter.", output_format))
      .cc_safe_save_ggplot(p, file, width = 6 * base_width, height = 5 * base_height, res = res)
    }
    log_message("Plot creation completed", message_type = "success", verbose = verbose)
    return(p)
  }

  if (plot_type == "role_heatmap") {
    p <- CellChat::netAnalysis_signalingRole_heatmap(
      cellchat_object,
      signaling = signaling,
      pattern = pattern,
      title = title.name %||% plot_label,
      font.size = font.size,
      ...
    )
    if (!is.null(dirpath)) {
      file <- file.path(dirpath, paste0(plot_label, "_role_heatmap_", pattern, ".", output_format))
      .cc_draw_heatmap_to_file(
        heatmap_obj = p,
        file = file,
        output_format = output_format,
        height = 6 * base_height,
        width = 6 * base_width,
        res = res
      )
    }
    log_message("Plot creation completed", message_type = "success", verbose = verbose)
    return(p)
  }

  if (plot_type == "role_change") {
    if (is.null(idents.use)) {
      log_message(
        "{.arg idents.use} must be provided for {.val plot_type = 'role_change'}",
        message_type = "error"
      )
    }
    p <- CellChat::netAnalysis_signalingChanges_scatter(
      merged_object,
      idents.use = idents.use,
      comparison = comp_idx,
      signaling = signaling,
      signaling.exclude = signaling.exclude,
      label.size = max(2.5, label.size),
      title = title.name %||% .cc_cmp_label(cmp, comp_idx),
      font.size = font.size,
      ...
    )
    p <- .cc_finalize_ggplot(
      p,
      title = title.name %||% .cc_cmp_label(cmp, comp_idx),
      subtitle = subtitle,
      xlab = xlab,
      ylab = ylab,
      legend.position = legend.position,
      legend.direction = legend.direction,
      legend.title = legend.title,
      theme_use = theme_use,
      theme_args = theme_args,
      aspect.ratio = aspect.ratio,
      font.size = font.size
    )
    if (!is.null(dirpath)) {
      file <- file.path(dirpath, paste0(.cc_cmp_label(cmp, comp_idx), "_role_change.", output_format))
      .cc_safe_save_ggplot(p, file, width = 6 * base_width, height = 5 * base_height, res = res)
    }
    log_message("Plot creation completed", message_type = "success", verbose = verbose)
    return(p)
  }

  log_message(
    "Unknown {.arg plot_type}: {.val {plot_type}}",
    message_type = "error"
  )
}

.cc_get_store <- function(x) {
  if (!inherits(x, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  x@tools[["CellChat"]]
}

.cc_result_names <- function(srt) {
  store <- .cc_get_store(srt)
  if (is.null(store) || is.null(store$results)) {
    return(character(0))
  }
  names(store$results)
}

.cc_comparison_names <- function(srt) {
  store <- .cc_get_store(srt)
  if (is.null(store) || is.null(store$comparisons)) {
    return(character(0))
  }
  names(store$comparisons)
}

.cc_resolve_single_condition <- function(srt, condition = NULL) {
  result_names <- .cc_result_names(srt)
  if (!is.null(condition)) {
    if (!condition %in% result_names) {
      log_message(
        "Condition {.val {condition}} not found in CellChat results",
        message_type = "error"
      )
    }
    return(condition)
  }
  if (length(result_names) == 1L) {
    return(result_names[1])
  }
  if ("ALL" %in% result_names) {
    return("ALL")
  }
  log_message(
    "Multiple CellChat results found. Please specify {.arg condition}",
    message_type = "error"
  )
}

.cc_should_use_single_condition <- function(srt, condition = NULL) {
  result_names <- .cc_result_names(srt)
  cmp_names <- .cc_comparison_names(srt)
  if (!is.null(condition)) {
    if (condition %in% result_names) {
      return(TRUE)
    }
    if (condition %in% cmp_names) {
      return(FALSE)
    }
    log_message(
      "{.arg condition} must be one of CellChat result names or comparison names",
      message_type = "error"
    )
  }
  if (length(result_names) == 1L && length(cmp_names) == 0L) {
    return(TRUE)
  }
  if ("ALL" %in% result_names && length(cmp_names) == 0L) {
    return(TRUE)
  }
  if (length(cmp_names) == 1L && length(result_names) == 0L) {
    return(FALSE)
  }
  log_message(
    "The requested plot is ambiguous because both single-condition results and/or comparison results are available. Please specify {.arg condition}",
    message_type = "error"
  )
}

.cc_get_single_obj <- function(srt, condition = NULL) {
  store <- .cc_get_store(srt)
  condition <- .cc_resolve_single_condition(srt, condition = condition)
  store$results[[condition]]$cellchat_object
}

.cc_get_cmp <- function(srt, condition = NULL) {
  store <- .cc_get_store(srt)
  cmp_names <- .cc_comparison_names(srt)
  if (is.null(condition)) {
    if (length(cmp_names) == 1L) {
      condition <- cmp_names[1]
    } else if (length(cmp_names) == 0L) {
      log_message(
        "No CellChat comparisons found",
        message_type = "error"
      )
    } else {
      log_message(
        "Multiple CellChat comparisons found. Please specify {.arg condition}",
        message_type = "error"
      )
    }
  }
  if (!condition %in% cmp_names) {
    log_message(
      "Comparison {.val {condition}} not found in CellChat results",
      message_type = "error"
    )
  }
  store$comparisons[[condition]]
}

.cc_get_dataset_object <- function(srt, condition = NULL, dataset = 1) {
  store <- .cc_get_store(srt)
  result_names <- .cc_result_names(srt)
  cmp_names <- .cc_comparison_names(srt)

  if (!is.null(condition) && condition %in% result_names) {
    return(list(
      object = store$results[[condition]]$cellchat_object,
      seurat_object = store$results[[condition]]$seurat_object,
      label = condition,
      source = "single"
    ))
  }

  if (!is.null(condition) && condition %in% cmp_names) {
    cmp <- .cc_get_cmp(srt, condition = condition)
    ds_name <- .cc_pick_dataset_name(cmp, dataset)
    seu_object <- NULL
    if (!is.null(store$results[[ds_name]]) && !is.null(store$results[[ds_name]]$seurat_object)) {
      seu_object <- store$results[[ds_name]]$seurat_object
    }
    return(list(
      object = cmp$object.list[[ds_name]],
      seurat_object = seu_object,
      label = ds_name,
      source = "comparison"
    ))
  }

  if (is.null(condition) && length(result_names) == 1L) {
    return(list(
      object = store$results[[result_names[1]]]$cellchat_object,
      seurat_object = store$results[[result_names[1]]]$seurat_object,
      label = result_names[1],
      source = "single"
    ))
  }

  if (is.null(condition) && "ALL" %in% result_names) {
    return(list(
      object = store$results[["ALL"]]$cellchat_object,
      seurat_object = store$results[["ALL"]]$seurat_object,
      label = "ALL",
      source = "single"
    ))
  }

  if (is.null(condition) && length(cmp_names) == 1L) {
    cmp <- .cc_get_cmp(srt, condition = cmp_names[1])
    ds_name <- .cc_pick_dataset_name(cmp, dataset)
    seu_object <- NULL
    if (!is.null(store$results[[ds_name]]) && !is.null(store$results[[ds_name]]$seurat_object)) {
      seu_object <- store$results[[ds_name]]$seurat_object
    }
    return(list(
      object = cmp$object.list[[ds_name]],
      seurat_object = seu_object,
      label = ds_name,
      source = "comparison"
    ))
  }

  log_message(
    "Unable to determine which CellChat object to plot. Please specify {.arg condition}",
    message_type = "error"
  )
}

.cc_resolve_dataset_index <- function(cmp, comparison = c(1, 2)) {
  ds_names <- names(cmp$object.list)
  if (is.numeric(comparison)) {
    idx <- as.integer(comparison)
    if (any(is.na(idx)) || any(idx < 1L) || any(idx > length(ds_names))) {
      log_message(
        "comparison indices out of range. Available datasets: {.val {seq_along(ds_names)}}",
        message_type = "error"
      )
    }
    return(idx)
  }
  idx <- match(as.character(comparison), ds_names)
  if (any(is.na(idx))) {
    log_message(
      "comparison names not found: {.val {as.character(comparison)[is.na(idx)]}}. Available datasets: {.val {ds_names}}",
      message_type = "error"
    )
  }
  idx
}

.cc_pick_dataset_name <- function(cmp, dataset = 1) {
  nm <- names(cmp$object.list)
  if (is.character(dataset)) {
    if (!dataset %in% nm) {
      log_message(
        "dataset {.val {dataset}} not found. Available: {.val {nm}}",
        message_type = "error"
      )
    }
    return(dataset)
  }
  if (length(dataset) != 1L || dataset < 1L || dataset > length(nm)) {
    log_message(
      "dataset index out of range. Available: {.val {seq_along(nm)}}",
      message_type = "error"
    )
  }
  nm[dataset]
}

.cc_resolve_group_index_single <- function(object, groups = NULL) {
  if (is.null(groups)) {
    return(NULL)
  }
  lev <- levels(object@idents)
  if (is.numeric(groups)) {
    return(as.integer(groups))
  }
  idx <- match(as.character(groups), lev)
  if (any(is.na(idx))) {
    log_message(
      "These cell groups were not found: {.val {as.character(groups)[is.na(idx)]}}. Available groups: {.val {lev}}",
      message_type = "error"
    )
  }
  idx
}

.cc_resolve_group_index_merged <- function(object, groups = NULL) {
  if (is.null(groups)) {
    return(NULL)
  }
  lev <- if (is.list(object@idents) && !is.null(object@idents$joint)) {
    levels(object@idents$joint)
  } else {
    levels(object@idents)
  }
  if (is.numeric(groups)) {
    return(as.integer(groups))
  }
  idx <- match(as.character(groups), lev)
  if (any(is.na(idx))) {
    log_message(
      "These cell groups were not found in merged object: {.val {as.character(groups)[is.na(idx)]}}. Available groups: {.val {lev}}",
      message_type = "error"
    )
  }
  idx
}

.cc_cmp_label <- function(cmp, comp_idx = c(1, 2)) {
  nm <- names(cmp$object.list)[comp_idx]
  if (length(nm) <= 1L) {
    return(nm[1])
  }
  paste0(nm[1], "_vs_", nm[2])
}

.cc_get_palette_colors <- function(palette = "Paired", palcolor = NULL, n = 10) {
  if (!is.null(palcolor)) {
    return(grDevices::colorRampPalette(palcolor)(n))
  }

  cols <- tryCatch(
    get_colors(palette),
    error = function(e) NULL
  )

  if (!is.null(cols)) {
    if (is.data.frame(cols) && "hex" %in% colnames(cols)) {
      cols <- cols$hex
    }
    cols <- as.character(cols)
    cols <- cols[!is.na(cols)]
    if (length(cols) > 0L) {
      return(grDevices::colorRampPalette(cols)(n))
    }
  }

  pal_fallback <- if (palette %in% grDevices::hcl.pals()) palette else "viridis"
  grDevices::hcl.colors(n, palette = pal_fallback)
}

.cc_theme_object <- function(theme_use = "theme_scop", theme_args = list()) {
  if (is.null(theme_use)) {
    return(NULL)
  }
  if (inherits(theme_use, "theme")) {
    return(theme_use)
  }
  theme_fun <- if (is.character(theme_use)) {
    tryCatch(get(theme_use, mode = "function", inherits = TRUE), error = function(e) NULL)
  } else if (is.function(theme_use)) {
    theme_use
  } else {
    NULL
  }

  if (is.null(theme_fun)) {
    theme_fun <- ggplot2::theme_bw
  }

  do.call(theme_fun, theme_args)
}

.cc_finalize_ggplot <- function(
    p,
    title = NULL,
    subtitle = NULL,
    xlab = NULL,
    ylab = NULL,
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = NULL,
    theme_use = "theme_scop",
    theme_args = list(),
    aspect.ratio = NULL,
    font.size = NULL
) {
  if (!inherits(p, c("gg", "ggplot"))) {
    return(p)
  }
  theme_obj <- .cc_theme_object(theme_use = theme_use, theme_args = theme_args)
  p <- p +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab,
      color = legend.title,
      fill = legend.title,
      size = legend.title
    ) +
    ggplot2::theme(
      legend.position = legend.position,
      legend.direction = legend.direction
    )
  if (!is.null(theme_obj)) {
    p <- p + theme_obj
  }
  if (!is.null(aspect.ratio)) {
    p <- p + ggplot2::theme(aspect.ratio = aspect.ratio)
  }
  if (!is.null(font.size)) {
    p <- p + ggplot2::theme(
      axis.text = ggplot2::element_text(size = font.size),
      axis.title = ggplot2::element_text(size = font.size + 1),
      plot.title = ggplot2::element_text(size = font.size + 2, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = font.size),
      legend.text = ggplot2::element_text(size = font.size),
      legend.title = ggplot2::element_text(size = font.size + 1)
    )
  }
  p
}

.cc_finalize_patchwork <- function(
    p,
    title = NULL,
    subtitle = NULL,
    legend.position = "right",
    legend.direction = "vertical",
    theme_use = "theme_scop",
    theme_args = list(),
    font.size = NULL
) {
  if (!inherits(p, "patchwork")) {
    return(.cc_finalize_ggplot(
      p,
      title = title,
      subtitle = subtitle,
      legend.position = legend.position,
      legend.direction = legend.direction,
      theme_use = theme_use,
      theme_args = theme_args,
      font.size = font.size
    ))
  }
  theme_obj <- .cc_theme_object(theme_use = theme_use, theme_args = theme_args)
  if (!is.null(title) || !is.null(subtitle)) {
    p <- p + patchwork::plot_annotation(title = title, subtitle = subtitle)
  }
  if (!is.null(theme_obj)) {
    p <- p & theme_obj & ggplot2::theme(
      legend.position = legend.position,
      legend.direction = legend.direction
    )
  }
  p
}

.cc_finalize_plot_list <- function(plot_list, combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE) {
  plot_list <- Filter(Negate(is.null), plot_list)
  if (length(plot_list) == 0L) {
    return(NULL)
  }
  if (length(plot_list) == 1L) {
    return(plot_list[[1]])
  }
  if (!isTRUE(combine)) {
    return(plot_list)
  }
  if (all(vapply(plot_list, function(x) inherits(x, c("gg", "ggplot", "patchwork")), logical(1)))) {
    return(patchwork::wrap_plots(plot_list, nrow = nrow, ncol = ncol, byrow = byrow))
  }
  plot_list
}

.cc_simplify_plot_list <- function(plot_list) {
  plot_list <- Filter(Negate(is.null), plot_list)
  if (length(plot_list) == 1L) {
    return(plot_list[[1]])
  }
  plot_list
}

.cc_capture_drawn_plot <- function(x) {
  out <- x
  if (inherits(out, c("gg", "ggplot", "patchwork", "recordedplot", "Heatmap", "HeatmapList"))) {
    return(out)
  }
  tryCatch(grDevices::recordPlot(), error = function(e) out)
}

create_output_device <- function(filename, output_format = "pdf", height = 6, width = 6, res = 300, units = "in") {
  if (output_format == "png") {
    grDevices::png(filename, height = height, width = width, res = res, units = units)
  } else if (output_format == "pdf") {
    grDevices::pdf(filename, height = height, width = width)
  }
}

save_ggplot <- function(plot, filename, output_format = "pdf", height = 6, width = 6, res = 300) {
  if (output_format == "png") {
    ggsave(filename, plot = plot, height = height, width = width, dpi = res)
  } else if (output_format == "pdf") {
    ggsave(filename, plot = plot, height = height, width = width)
  }
}

.cc_safe_save_ggplot <- function(plot, filename, width = 6, height = 6, res = 300) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ext <- tolower(tools::file_ext(filename))
  save_ggplot(
    plot = plot,
    filename = filename,
    output_format = if (ext %in% c("png", "pdf")) ext else "pdf",
    height = height,
    width = width,
    res = res
  )
}

.cc_replay_recorded_plot <- function(plot_obj, file, output_format = "pdf", height = 6, width = 6, res = 300) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  create_output_device(
    filename = file,
    output_format = output_format,
    height = height,
    width = width,
    res = res
  )
  grDevices::replayPlot(plot_obj)
  grDevices::dev.off()
  invisible(file)
}

.cc_draw_heatmap_to_file <- function(heatmap_obj, file, output_format = "pdf", height = 6, width = 6, res = 300) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  create_output_device(
    filename = file,
    output_format = output_format,
    height = height,
    width = width,
    res = res
  )
  ComplexHeatmap::draw(heatmap_obj)
  grDevices::dev.off()
  invisible(file)
}

.cc_plot_aggregate <- function(
    cellchat_object,
    condition,
    title = NULL,
    font.size = 10,
    label.size = 1,
    theme_use = "theme_scop",
    theme_args = list(),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = NULL,
    aspect.ratio = NULL,
    xlab = NULL,
    ylab = NULL,
    subtitle = NULL,
    dirpath = NULL,
    output_format = "pdf",
    base_height = 1,
    base_width = 1,
    res = 300,
    verbose = TRUE,
    ...
) {
  p_count <- .cc_plot_circle(
    cellchat_object = cellchat_object,
    condition = condition,
    measure = "count",
    title.name = "Number of interactions",
    label.size = label.size,
    dirpath = dirpath,
    output_format = output_format,
    base_height = base_height,
    base_width = base_width,
    res = res,
    ...
  )

  p_weight <- .cc_plot_circle(
    cellchat_object = cellchat_object,
    condition = condition,
    measure = "weight",
    title.name = "Interaction weights/strength",
    label.size = label.size,
    dirpath = dirpath,
    output_format = output_format,
    base_height = base_height,
    base_width = base_width,
    res = res,
    ...
  )

  p_scatter <- CellChat::netAnalysis_signalingRole_scatter(
    cellchat_object,
    title = title %||% condition,
    label.size = max(2.5, label.size),
    font.size = font.size
  )
  p_scatter <- .cc_finalize_ggplot(
    p_scatter,
    title = title %||% condition,
    subtitle = subtitle,
    xlab = xlab,
    ylab = ylab,
    legend.position = legend.position,
    legend.direction = legend.direction,
    legend.title = legend.title,
    theme_use = theme_use,
    theme_args = theme_args,
    aspect.ratio = aspect.ratio,
    font.size = font.size
  )
  if (!is.null(dirpath)) {
    file <- file.path(dirpath, paste0(condition, "_aggregate_role_scatter.", output_format))
    .cc_safe_save_ggplot(p_scatter, file, width = 5.5 * base_width, height = 4.5 * base_height, res = res)
  }

  list(
    circle_count = p_count,
    circle_weight = p_weight,
    role_scatter = p_scatter
  )
}

.cc_plot_circle <- function(
    cellchat_object,
    condition,
    measure = c("count", "weight"),
    title.name = NULL,
    label.size = 1,
    dirpath = NULL,
    output_format = "pdf",
    base_height = 1,
    base_width = 1,
    res = 300,
    ...
) {
  measure <- match.arg(measure)
  group_size <- as.numeric(table(cellchat_object@idents))
  mat <- if (identical(measure, "count")) cellchat_object@net$count else cellchat_object@net$weight
  ng <- length(group_size)

  height_in <- max(4.5, 3.2 + 0.15 * ng) * base_height
  width_in <- max(5.5, 3.8 + 0.2 * ng) * base_width

  if (!is.null(dirpath)) {
    file <- file.path(dirpath, paste0(condition, "_circle_", measure, ".", output_format))
    create_output_device(file, output_format = output_format, height = height_in, width = width_in, res = res)
  }
  CellChat::netVisual_circle(
    mat,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    vertex.label.cex = max(0.8, label.size),
    title.name = title.name %||% paste(condition, measure),
    ...
  )
  p <- tryCatch(grDevices::recordPlot(), error = function(e) NULL)
  if (!is.null(dirpath)) {
    grDevices::dev.off()
  }
  p
}

.cc_subset_table <- function(
    object,
    slot.name = "net",
    signaling = NULL,
    pairLR.use = NULL,
    sources.use = NULL,
    targets.use = NULL,
    thresh = 0.05,
    dataset = NULL
) {
  df <- CellChat::subsetCommunication(
    object = object,
    slot.name = slot.name,
    sources.use = sources.use,
    targets.use = targets.use,
    signaling = signaling,
    pairLR.use = pairLR.use,
    thresh = thresh
  )
  df <- as.data.frame(df)
  if (!is.null(dataset) && nrow(df) > 0L) {
    df$dataset <- dataset
  }
  df
}

.cc_prepare_bubble_data <- function(
    object,
    slot.name = "net",
    signaling = NULL,
    pairLR.use = NULL,
    sources.use = NULL,
    targets.use = NULL,
    thresh = 0.05,
    dataset = NULL,
    top_n = 10,
    remove.isolate = TRUE
) {
  df <- .cc_subset_table(
    object = object,
    slot.name = slot.name,
    signaling = signaling,
    pairLR.use = pairLR.use,
    sources.use = sources.use,
    targets.use = targets.use,
    thresh = thresh,
    dataset = dataset
  )
  if (is.null(df) || nrow(df) == 0L) {
    return(df)
  }

  if (isTRUE(remove.isolate)) {
    df <- df[df$prob > 0, , drop = FALSE]
  }
  if (nrow(df) == 0L) {
    return(df)
  }

  if (!is.null(top_n) && is.numeric(top_n) && length(top_n) == 1L && top_n > 0L && is.null(pairLR.use)) {
    if (length(unique(df$interaction_name)) > top_n) {
      score_interaction <- sort(
        tapply(df$prob, as.character(df$interaction_name), sum, na.rm = TRUE),
        decreasing = TRUE
      )
      keep_interaction <- names(score_interaction)[seq_len(min(top_n, length(score_interaction)))]
      df <- df[df$interaction_name %in% keep_interaction, , drop = FALSE]
    }

    max_pairs <- max(15L, as.integer(top_n) * 2L)
    pair_vec <- paste(df$source, "->", df$target)
    if (length(unique(pair_vec)) > max_pairs && is.null(sources.use) && is.null(targets.use)) {
      score_pair <- sort(
        tapply(df$prob, pair_vec, sum, na.rm = TRUE),
        decreasing = TRUE
      )
      keep_pair <- names(score_pair)[seq_len(min(max_pairs, length(score_pair)))]
      df <- df[paste(df$source, "->", df$target) %in% keep_pair, , drop = FALSE]
    }
  }

  df
}

.cc_bubble_dims <- function(df, base_height = 1, base_width = 1) {
  n_pair <- length(unique(paste(df$source, "->", df$target, if ("dataset" %in% colnames(df)) paste0("[", df$dataset, "]") else "")))
  n_lr <- length(unique(df$interaction_name))
  list(
    width = min(14, max(6, 2.5 + 0.28 * n_pair)) * base_width,
    height = min(12, max(4.5, 2.5 + 0.16 * n_lr)) * base_height
  )
}

.cc_custom_bubble_plot <- function(
    df,
    title = NULL,
    subtitle = NULL,
    xlab = NULL,
    ylab = NULL,
    color.by = c("prob", "pval"),
    bubble_size.range = c(1.5, 8),
    palette = "Paired",
    palcolor = NULL,
    font.size = 10,
    angle.x = 45,
    hjust.x = 1,
    vjust.x = 1,
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = NULL,
    theme_use = "theme_scop",
    theme_args = list(),
    aspect.ratio = NULL,
    remove.isolate = TRUE
) {
  color.by <- match.arg(color.by)
  if (is.null(df) || nrow(df) == 0L) {
    log_message(
      "No communication records available for bubble plot",
      message_type = "error"
    )
  }

  req_cols <- c("source", "target", "interaction_name", "prob", "pval")
  miss <- setdiff(req_cols, colnames(df))
  if (length(miss) > 0L) {
    log_message(
      "Missing required columns for bubble plot: {.val {miss}}",
      message_type = "error"
    )
  }

  if (isTRUE(remove.isolate)) {
    df <- df[df$prob > 0, , drop = FALSE]
  }
  if (nrow(df) == 0L) {
    log_message(
      "No non-zero communication records remain after filtering",
      message_type = "error"
    )
  }

  df$pair <- paste(df$source, "->", df$target)
  if ("dataset" %in% colnames(df)) {
    df$pair <- paste0(df$pair, " [", df$dataset, "]")
  }
  df$interaction_plot <- as.character(df$interaction_name)

  pair_order <- sort(tapply(df$prob, df$pair, sum, na.rm = TRUE), decreasing = TRUE)
  interaction_order <- sort(tapply(df$prob, df$interaction_plot, sum, na.rm = TRUE), decreasing = TRUE)
  df$pair <- factor(df$pair, levels = names(pair_order))
  df$interaction_plot <- factor(df$interaction_plot, levels = rev(names(interaction_order)))

  cols <- .cc_get_palette_colors(palette = palette, palcolor = palcolor, n = 9)

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = pair,
      y = interaction_plot,
      size = prob,
      color = if (identical(color.by, "pval")) -log10(pval + 1e-300) else prob
    )
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_size_continuous(range = bubble_size.range) +
    ggplot2::scale_color_gradientn(colours = cols) +
    ggplot2::labs(
      x = xlab,
      y = ylab,
      title = title,
      subtitle = subtitle,
      size = "prob",
      color = legend.title %||% if (identical(color.by, "pval")) "-log10(pval)" else "prob"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = angle.x, hjust = hjust.x, vjust = vjust.x),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey90"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      color = ggplot2::guide_colorbar(order = 1),
      size = ggplot2::guide_legend(order = 2)
    )

  p <- .cc_finalize_ggplot(
    p,
    title = title,
    subtitle = subtitle,
    xlab = xlab,
    ylab = ylab,
    legend.position = legend.position,
    legend.direction = legend.direction,
    legend.title = legend.title,
    theme_use = theme_use,
    theme_args = theme_args,
    aspect.ratio = aspect.ratio,
    font.size = font.size
  )
  p
}

.cc_chord_plot <- function(
    cellchat_object,
    signaling,
    sources.use = NULL,
    targets.use = NULL,
    thresh = 0.05,
    reduce = TRUE,
    top_n = 30,
    max.groups = 8,
    small.gap = 1,
    big.gap = 8,
    lab.cex = 0.6,
    title.name = NULL,
    dirpath = NULL,
    output_format = "pdf",
    fileprefix = "cellchat_chord",
    base_height = 1,
    base_width = 1,
    res = 300,
    ...
) {
  src_idx <- .cc_resolve_group_index_single(cellchat_object, sources.use)
  tgt_idx <- .cc_resolve_group_index_single(cellchat_object, targets.use)
  all_groups <- levels(cellchat_object@idents)

  df <- tryCatch(
    CellChat::subsetCommunication(
      object = cellchat_object,
      signaling = signaling,
      sources.use = src_idx,
      targets.use = tgt_idx,
      thresh = thresh
    ),
    error = function(e) NULL
  )

  if (is.null(df) || nrow(df) == 0L) {
    log_message(
      "There is no significant communication for {.val {signaling}}",
      message_type = "error"
    )
  }

  df <- as.data.frame(df)
  df <- df[!is.na(df$prob), , drop = FALSE]
  if (nrow(df) == 0L) {
    log_message(
      "No valid communication records remain for the chord plot",
      message_type = "error"
    )
  }

  if (isTRUE(reduce)) {
    df <- df[order(df$prob, decreasing = TRUE), , drop = FALSE]
    if (nrow(df) > top_n) {
      df <- df[seq_len(top_n), , drop = FALSE]
    }

    group_score <- tapply(df$prob, as.character(df$source), sum, na.rm = TRUE)
    target_score <- tapply(df$prob, as.character(df$target), sum, na.rm = TRUE)
    score <- tapply(
      c(group_score, target_score),
      c(names(group_score), names(target_score)),
      sum,
      na.rm = TRUE
    )
    score <- sort(score, decreasing = TRUE)
    keep_groups <- names(score)[seq_len(min(max.groups, length(score)))]
    df <- df[df$source %in% keep_groups & df$target %in% keep_groups, , drop = FALSE]
  }

  if (nrow(df) == 0L) {
    log_message(
      "No communication records remain after reduction for the chord plot",
      message_type = "error"
    )
  }

  src_use2 <- match(unique(as.character(df$source)), all_groups)
  tgt_use2 <- match(unique(as.character(df$target)), all_groups)
  src_use2 <- src_use2[!is.na(src_use2)]
  tgt_use2 <- tgt_use2[!is.na(tgt_use2)]

  if (!is.null(dirpath)) {
    file <- file.path(dirpath, paste0(fileprefix, ".", output_format))
    create_output_device(file, output_format = output_format, height = 7 * base_height, width = 7 * base_width, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  circlize::circos.clear()
  on.exit(try(circlize::circos.clear(), silent = TRUE), add = TRUE)

  p <- CellChat::netVisual_chord_cell(
    object = cellchat_object,
    signaling = signaling,
    sources.use = src_use2,
    targets.use = tgt_use2,
    lab.cex = lab.cex,
    small.gap = small.gap,
    big.gap = big.gap,
    title.name = title.name %||% paste(signaling, collapse = ", "),
    ...
  )
  p
}

.cc_create_gene_plot <- function(
    cellchat_object,
    seurat_object,
    signaling,
    species,
    title = NULL,
    theme_use = "theme_scop",
    theme_args = list(),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = NULL,
    subtitle = NULL,
    xlab = NULL,
    ylab = NULL,
    aspect.ratio = NULL,
    ...
) {
  pairLR <- CellChat::extractEnrichedLR(
    cellchat_object,
    signaling = signaling,
    geneLR.return = FALSE
  )
  lr_genes <- unique(unlist(strsplit(split = "_", x = pairLR$interaction_name)))
  genes_upper <- toupper(rownames(seurat_object))
  genes1 <- lr_genes[lr_genes %in% genes_upper]
  genes2 <- lr_genes[!(lr_genes %in% genes_upper)]
  genes22 <- sub(".?.?", "", genes2)
  lr_genes <- c(genes1, genes22[genes22 %in% genes_upper])

  if (identical(species, "mouse")) {
    genes <- rownames(cellchat_object@data)
    indices <- match(lr_genes, toupper(genes))
    lr_genes <- genes[indices]
  }
  lr_genes <- unique(stats::na.omit(lr_genes))
  if (length(lr_genes) == 0L) {
    log_message(
      "No ligand/receptor genes found in the Seurat object for signaling {.val {signaling}}",
      message_type = "error"
    )
  }

  p <- VlnPlot(
    object = seurat_object,
    features = lr_genes,
    pt.size = -1,
    stack = TRUE,
    ...
  )
  p <- .cc_finalize_patchwork(
    p,
    title = title %||% signaling,
    subtitle = subtitle,
    legend.position = legend.position,
    legend.direction = legend.direction,
    theme_use = theme_use,
    theme_args = theme_args,
    font.size = 10
  )
  p
}
