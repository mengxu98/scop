ExpressionStatPlot <- function(
    exp.data,
    meta.data,
    stat.by,
    group.by = NULL,
    split.by = NULL,
    bg.by = NULL,
    plot.by = c("group", "feature"),
    fill.by = c("group", "feature", "expression"),
    cells = NULL,
    keep_empty = FALSE,
    individual = FALSE,
    plot_type = c("violin", "box", "bar", "dot", "col"),
    palette = "Paired",
    palcolor = NULL,
    alpha = 1,
    bg_palette = "Paired",
    bg_palcolor = NULL,
    bg_alpha = 0.2,
    add_box = FALSE,
    box_color = "black",
    box_width = 0.1,
    box_ptsize = 2,
    add_point = FALSE,
    pt.color = "grey30",
    pt.size = NULL,
    pt.alpha = 1,
    jitter.width = 0.4,
    jitter.height = 0.1,
    add_trend = FALSE,
    trend_color = "black",
    trend_linewidth = 1,
    trend_ptsize = 2,
    add_stat = c("none", "mean", "median"),
    stat_color = "black",
    stat_size = 1,
    stat_stroke = 1,
    stat_shape = 25,
    add_line = NULL,
    line_color = "red",
    line_size = 1,
    line_type = 1,
    cells.highlight = NULL,
    cols.highlight = "red",
    sizes.highlight = 1,
    alpha.highlight = 1,
    calculate_coexp = FALSE,
    same.y.lims = FALSE,
    y.min = NULL,
    y.max = NULL,
    y.trans = "identity",
    y.nbreaks = 5,
    sort = FALSE,
    stack = FALSE,
    flip = FALSE,
    comparisons = NULL,
    ref_group = NULL,
    pairwise_method = "wilcox.test",
    multiplegroup_comparisons = FALSE,
    multiple_method = "kruskal.test",
    sig_label = c("p.signif", "p.format"),
    sig_labelsize = 3.5,
    aspect.ratio = NULL,
    title = NULL,
    subtitle = NULL,
    xlab = NULL,
    ylab = "Expression level",
    legend.position = "right",
    legend.direction = "vertical",
    theme_use = "theme_scop",
    theme_args = list(),
    force = FALSE,
    seed = 11) {
  set.seed(seed)

  plot.by <- match.arg(plot.by)
  plot_type <- match.arg(plot_type)
  fill.by <- match.arg(fill.by)
  sig_label <- match.arg(sig_label)
  add_stat <- match.arg(add_stat)
  if (!is.null(add_line)) {
    stopifnot(is.numeric(add_line))
  }

  if (missing(exp.data)) {
    exp.data <- matrix(
      0,
      nrow = 1,
      ncol = nrow(meta.data),
      dimnames = list("", rownames(meta.data))
    )
  }

  allfeatures <- rownames(exp.data)
  allcells <- rownames(meta.data)

  if (plot_type == "col") {
    if (
      isTRUE(add_box) ||
        isTRUE(add_point) ||
        isTRUE(add_trend) ||
        isTRUE(add_stat != "none")
    ) {
      log_message(
        "Cannot add other layers when plot_type is 'col'",
        message_type = "warning"
      )
      add_box <- add_point <- add_trend <- FALSE
    }
  }
  if (
    (isTRUE(multiplegroup_comparisons) || length(comparisons) > 0) &&
      plot_type %in% c("col")
  ) {
    log_message(
      "Cannot add comparison when plot_type is 'col'",
      message_type = "warning"
    )
    multiplegroup_comparisons <- FALSE
    comparisons <- NULL
  }
  if (isTRUE(comparisons) && is.null(split.by)) {
    log_message(
      "'split.by' must provided when comparisons=TRUE",
      message_type = "error"
    )
  }

  if (nrow(meta.data) == 0) {
    log_message(
      "meta.data is empty.",
      message_type = "error"
    )
  }
  if (is.null(group.by)) {
    group.by <- "All.groups"
    xlab <- ""
    meta.data[[group.by]] <- factor("")
  }
  if (is.null(split.by)) {
    split.by <- "All.groups"
    meta.data[[split.by]] <- factor("")
  }
  if (group.by == split.by && group.by == "All.groups") {
    legend.position <- "none"
  }
  for (i in unique(c(group.by, split.by, bg.by))) {
    if (!i %in% colnames(meta.data)) {
      log_message(
        paste0(i, " is not in the meta.data."),
        message_type = "error"
      )
    }
    if (!is.factor(meta.data[[i]])) {
      meta.data[[i]] <- factor(meta.data[[i]], levels = unique(meta.data[[i]]))
    }
  }
  bg_map <- NULL
  if (!is.null(bg.by)) {
    for (g in group.by) {
      df_table <- table(meta.data[[g]], meta.data[[bg.by]])
      if (max(Matrix::rowSums(df_table > 0), na.rm = TRUE) > 1) {
        log_message(
          "'group.by' must be a part of 'bg.by'",
          message_type = "error"
        )
      } else {
        bg_map[[g]] <- stats::setNames(
          colnames(df_table)[apply(df_table, 1, function(x) which(x > 0))],
          rownames(df_table)
        )
      }
    }
  } else {
    for (g in group.by) {
      bg_map[[g]] <- stats::setNames(levels(meta.data[[g]]), levels(meta.data[[g]]))
    }
  }
  if (!is.null(cells.highlight) && isFALSE(cells.highlight)) {
    if (!any(cells.highlight %in% allcells)) {
      log_message(
        "No cells in 'cells.highlight' found.",
        message_type = "error"
      )
    }
    if (!all(cells.highlight %in% allcells)) {
      log_message(
        "Some cells in 'cells.highlight' not found.",
        message_type = "warning"
      )
    }
    cells.highlight <- intersect(cells.highlight, allcells)
  }
  if (isTRUE(cells.highlight)) {
    cells.highlight <- allcells
  }
  if (!is.null(cells.highlight) && isFALSE(add_point)) {
    log_message(
      "'cells.highlight' is valid only when add_point=TRUE.",
      message_type = "warning"
    )
  }
  if (isTRUE(stack) & isTRUE(sort)) {
    log_message(
      "Set sort to FALSE when stack is TRUE",
      message_type = "warning"
    )
    sort <- FALSE
  }
  if (isTRUE(multiplegroup_comparisons) || length(comparisons) > 0) {
    check_r("ggpubr")
    ncomp <- sapply(comparisons, length)
    if (any(ncomp > 2)) {
      log_message(
        "'comparisons' must be a list in which all elements must be vectors of length 2",
        message_type = "error"
      )
    }
  }

  stat.by <- unique(stat.by)
  features_drop <- stat.by[
    !stat.by %in% c(rownames(exp.data), colnames(meta.data))
  ]
  if (length(features_drop) > 0) {
    log_message(
      paste0(features_drop, collapse = ","),
      " are not found.",
      message_type = "warning"
    )
    stat.by <- stat.by[!stat.by %in% features_drop]
  }

  features_gene <- stat.by[stat.by %in% rownames(exp.data)]
  features_meta <- stat.by[stat.by %in% colnames(meta.data)]
  if (length(intersect(features_gene, features_meta)) > 0) {
    log_message(
      "Features appear in both gene names and metadata names: ",
      paste0(intersect(features_gene, features_meta), collapse = ","),
      message_type = "warning"
    )
  }

  if (isTRUE(calculate_coexp) && length(features_gene) > 1) {
    if (length(features_meta) > 0) {
      log_message(
        paste(features_meta, collapse = ","),
        "is not used when calculating co-expression",
        message_type = "warning"
      )
    }
    status <- CheckDataType(exp.data)
    log_message(
      "Data type: ", status
    )
    if (status %in% c("raw_counts", "raw_normalized_counts")) {
      meta.data[["CoExp"]] <- apply(
        exp.data[features_gene, , drop = FALSE],
        2,
        function(x) exp(mean(log(x)))
      )
    } else if (status == "log_normalized_counts") {
      meta.data[["CoExp"]] <- apply(
        expm1(exp.data[features_gene, , drop = FALSE]),
        2,
        function(x) log1p(exp(mean(log(x))))
      )
    } else {
      log_message(
        "Can not determine the data type.",
        message_type = "error"
      )
    }
    stat.by <- c(stat.by, "CoExp")
    features_meta <- c(features_meta, "CoExp")
  }
  if (length(features_gene) > 0) {
    if (all(allfeatures %in% features_gene)) {
      dat_gene <- Matrix::t(exp.data)
    } else {
      dat_gene <- Matrix::t(exp.data[features_gene, , drop = FALSE])
    }
  } else {
    dat_gene <- matrix(nrow = length(allcells), ncol = 0)
  }
  if (length(features_meta) > 0) {
    dat_meta <- as_matrix(meta.data[, features_meta, drop = FALSE])
  } else {
    dat_meta <- matrix(nrow = length(allcells), ncol = 0)
  }
  dat_exp <- cbind(dat_gene, dat_meta)
  stat.by <- unique(stat.by[stat.by %in% c(features_gene, features_meta)])

  if (!is.numeric(dat_exp) && !inherits(dat_exp, "Matrix")) {
    log_message(
      "'stat.by' must be type of numeric variable.",
      message_type = "error"
    )
  }
  dat_group <- meta.data[,
    unique(c("cells", group.by, bg.by, split.by)),
    drop = FALSE
  ]
  dat_use <- cbind(dat_group, dat_exp[row.names(dat_group), , drop = FALSE])
  if (!is.null(cells)) {
    dat_group <- dat_group[
      intersect(rownames(dat_group), cells), ,
      drop = FALSE
    ]
    dat_use <- dat_use[intersect(rownames(dat_use), cells), , drop = FALSE]
  }
  if (nrow(dat_group) == 0) {
    log_message(
      "No specified cells found.",
      message_type = "error"
    )
  }

  if (is.null(pt.size)) {
    pt.size <- min(3000 / nrow(dat_group), 0.5)
  }

  nlev <- sapply(dat_group, nlevels)
  nlev <- nlev[nlev > 100]
  if (length(nlev) > 0 && isFALSE(force)) {
    log_message(
      paste(names(nlev), sep = ","),
      " have more than 100 levels.",
      message_type = "warning"
    )
    answer <- utils::askYesNo("Are you sure to continue?", default = FALSE)
    if (isFALSE(answer)) {
      return(invisible(NULL))
    }
  }

  if (isTRUE(same.y.lims)) {
    valus <- as_matrix(dat_use[,
      stat.by,
      drop = FALSE
    ])[is.finite(as_matrix(dat_use[, stat.by, drop = FALSE]))]
    if (is.null(y.max)) {
      y.max <- max(valus, na.rm = TRUE)
    } else if (is.character(y.max)) {
      q.max <- as.numeric(sub("(^q)(\\d+)", "\\2", y.max)) / 100
      y.max <- stats::quantile(values, q.max, na.rm = TRUE)
    }
    if (is.null(y.min)) {
      y.min <- min(valus, na.rm = TRUE)
    } else if (is.character(y.min)) {
      q.min <- as.numeric(sub("(^q)(\\d+)", "\\2", y.min)) / 100
      y.min <- stats::quantile(values, q.min, na.rm = TRUE)
    }
  }

  plist <- list()

  comb_list <- list()
  comb <- expand.grid(
    group_name = group.by,
    stat_name = stat.by,
    stringsAsFactors = FALSE
  )
  if (isTRUE(individual)) {
    for (g in group.by) {
      comb_list[[g]] <- merge(
        comb,
        expand.grid(
          group_name = g,
          group_element = levels(dat_use[[g]]),
          split_name = levels(dat_use[[split.by]]),
          stringsAsFactors = FALSE
        ),
        by = "group_name",
        all = FALSE
      )
    }
  } else {
    for (g in group.by) {
      comb_list[[g]] <- merge(
        comb,
        expand.grid(
          group_name = g,
          group_element = list(levels(dat_use[[g]])),
          split_name = list(levels(dat_use[[split.by]])),
          stringsAsFactors = FALSE
        ),
        by = "group_name",
        all = FALSE
      )
    }
  }
  comb <- do.call(rbind, comb_list)
  rownames(comb) <- paste0(
    comb[["stat_name"]],
    ":",
    comb[["group_name"]],
    ":",
    sapply(comb[["group_element"]], function(x) paste0(x, collapse = ",")),
    ":",
    sapply(comb[["split_name"]], function(x) paste0(x, collapse = ","))
  )

  plist <- lapply(stats::setNames(rownames(comb), rownames(comb)), function(i) {
    g <- comb[i, "group_name"]
    f <- comb[i, "stat_name"]
    single_group <- comb[[i, "group_element"]]
    sp <- comb[[i, "split_name"]]
    xlab <- xlab %||% g
    ylab <- ylab %||% "Expression level"
    if (identical(theme_use, "theme_blank")) {
      theme_args[["xlab"]] <- xlab
      theme_args[["ylab"]] <- ylab
    }
    if (fill.by == "feature") {
      colors <- palette_colors(stat.by, palette = palette, palcolor = palcolor)
    }
    if (fill.by == "group") {
      if (split.by != "All.groups") {
        colors <- palette_colors(
          levels(dat_use[[split.by]]),
          palette = palette,
          palcolor = palcolor
        )
      } else {
        colors <- palette_colors(
          levels(dat_use[[g]]),
          palette = palette,
          palcolor = palcolor
        )
      }
    }
    if (fill.by == "expression") {
      median_values <- stats::aggregate(
        dat_use[, stat.by, drop = FALSE],
        by = list(dat_use[[g]], dat_use[[split.by]]),
        FUN = stats::median
      )
      rownames(median_values) <- paste0(
        median_values[, 1],
        "-",
        median_values[, 2]
      )
      colors <- palette_colors(
        unlist(median_values[, stat.by]),
        type = "continuous",
        palette = palette,
        palcolor = palcolor
      )
      colors_limits <- range(median_values[, stat.by])
    }

    dat <- dat_use[
      dat_use[[g]] %in% single_group & dat_use[[split.by]] %in% sp,
      c(colnames(dat_group), f)
    ]
    dat[[g]] <- factor(
      dat[[g]],
      levels = levels(dat[[g]])[levels(dat[[g]]) %in% dat[[g]]]
    )
    if (!is.null(bg.by)) {
      bg <- bg.by
      bg_color <- palette_colors(
        levels(dat[[bg]]),
        palette = bg_palette,
        palcolor = bg_palcolor
      )
    } else {
      bg <- g
      bg_color <- palette_colors(
        levels(dat[[bg]]),
        palcolor = bg_palcolor %||%
          rep(c("transparent", "grey85"), nlevels(dat[[bg]]))
      )
    }
    dat[["bg.by"]] <- dat[[bg]]
    dat[["value"]] <- dat[[f]]
    dat[["group.by"]] <- dat[[g]]
    dat[["split.by"]] <- dat[[split.by]]
    if (split.by == g) {
      dat[["split.by"]] <- dat[["group.by"]]
    }

    dat[, "features"] <- rep(f, nrow(dat))
    if (
      nrow(dat) > 0 && ((is.character(x = sort) && nchar(x = sort) > 0) || sort)
    ) {
      df_sort <- stats::aggregate(
        dat[, "value", drop = FALSE],
        by = list(dat[["group.by"]]),
        FUN = stats::median
      )
      if (is.character(sort) && sort == "increasing") {
        decreasing <- FALSE
      } else {
        decreasing <- TRUE
      }
      sortlevel <- as.character(df_sort[
        order(df_sort[["value"]], decreasing = decreasing),
        1
      ])
      dat[, "group.by"] <- factor(dat[, "group.by"], levels = sortlevel)
    }
    if (fill.by == "feature") {
      dat[, "fill.by"] <- rep(f, nrow(dat))
      keynm <- "Features"
    }
    if (fill.by == "group") {
      dat[, "fill.by"] <- if (split.by == "All.groups") {
        dat[, "group.by"]
      } else {
        dat[, "split.by"]
      }
      keynm <- ifelse(split.by == "All.groups", g, split.by)
    }
    if (fill.by == "expression") {
      dat[, "fill.by"] <- median_values[
        paste0(dat[["group.by"]], "-", dat[["split.by"]]),
        f
      ]
      keynm <- "Median expression"
    }
    if (split.by != "All.groups") {
      levels_order <- levels(dat[["split.by"]])
    } else {
      levels_order <- levels(dat[["group.by"]])
    }
    if (fill.by == "feature") {
      levels_order <- unique(stat.by)
    }

    group_comb <- expand.grid(
      x = levels(dat[["split.by"]]),
      y = levels(dat[["group.by"]])
    )
    dat[["group.unique"]] <- utils::head(
      factor(
        paste("sp", dat[["split.by"]], "gp", dat[["group.by"]], sep = "-"),
        levels = paste("sp", group_comb[[1]], "gp", group_comb[[2]], sep = "-")
      ),
      nrow(dat)
    )
    dat <- dat[order(dat[["group.unique"]]), , drop = FALSE]

    values <- dat[, "value"][is.finite(x = dat[, "value"])]
    if (is.null(y.max)) {
      y_max_use <- max(values, na.rm = TRUE)
    } else if (is.character(y.max)) {
      q.max <- as.numeric(sub("(^q)(\\d+)", "\\2", y.max)) / 100
      y_max_use <- stats::quantile(values, q.max, na.rm = TRUE)
    } else {
      y_max_use <- y.max
    }
    if (is.null(y.min)) {
      y_min_use <- min(values, na.rm = TRUE)
    } else if (is.character(y.min)) {
      q.min <- as.numeric(sub("(^q)(\\d+)", "\\2", y.min)) / 100
      y_min_use <- stats::quantile(values, q.min, na.rm = TRUE)
    } else {
      y_min_use <- y.min
    }

    if (isTRUE(flip)) {
      dat[["group.by"]] <- factor(
        dat[["group.by"]],
        levels = rev(levels(dat[["group.by"]]))
      )
      aspect.ratio <- 1 / aspect.ratio
      if (length(aspect.ratio) == 0 || is.na(aspect.ratio)) {
        aspect.ratio <- NULL
      }
    }

    if (plot_type == "col") {
      if (isTRUE(flip)) {
        dat[["cell"]] <- rev(seq_len(nrow(dat)))
      } else {
        dat[["cell"]] <- seq_len(nrow(dat))
      }
      p <- ggplot(
        dat,
        aes(
          x = .data[["cell"]],
          y = .data[["value"]],
          fill = .data[["fill.by"]]
        )
      )
    } else {
      p <- ggplot(
        dat,
        aes(
          x = .data[["group.by"]],
          y = .data[["value"]],
          fill = .data[["fill.by"]]
        )
      )
    }

    if (isFALSE(individual)) {
      if (plot_type == "col") {
        x_index <- split(dat[["cell"]], dat[["group.by"]])
        bg_data <- as.data.frame(Matrix::t(sapply(x_index, range)))
        colnames(bg_data) <- c("xmin", "xmax")
        bg_data[["group.by"]] <- names(x_index)
        bg_data[["xmin"]] <- ifelse(
          bg_data[["xmin"]] == min(bg_data[["xmax"]]),
          -Inf,
          bg_data[["xmin"]] - 0.5
        )
        bg_data[["xmax"]] <- ifelse(
          bg_data[["xmax"]] == max(bg_data[["xmax"]]),
          Inf,
          bg_data[["xmax"]] + 0.5
        )
        bg_data[["ymin"]] <- -Inf
        bg_data[["ymax"]] <- Inf
        bg_data[["fill"]] <- bg_color[bg_map[[g]][as.character(bg_data[[
          "group.by"
        ]])]]
      } else {
        bg_data <- unique(dat[, "group.by", drop = FALSE])
        bg_data[["x"]] <- as.numeric(bg_data[["group.by"]])
        bg_data[["xmin"]] <- ifelse(
          bg_data[["x"]] == min(bg_data[["x"]]),
          -Inf,
          bg_data[["x"]] - 0.5
        )
        bg_data[["xmax"]] <- ifelse(
          bg_data[["x"]] == max(bg_data[["x"]]),
          Inf,
          bg_data[["x"]] + 0.5
        )
        bg_data[["ymin"]] <- -Inf
        bg_data[["ymax"]] <- Inf
        bg_data[["fill"]] <- bg_color[bg_map[[g]][as.character(bg_data[[
          "group.by"
        ]])]]
      }
      bg_layer <- geom_rect(
        data = bg_data,
        xmin = bg_data[["xmin"]],
        xmax = bg_data[["xmax"]],
        ymin = bg_data[["ymin"]],
        ymax = bg_data[["ymax"]],
        fill = bg_data[["fill"]],
        alpha = bg_alpha,
        inherit.aes = FALSE
      )
      p <- p + bg_layer
    }

    if (plot_type %in% c("bar", "col")) {
      p <- p + geom_hline(yintercept = 0, linetype = 2)
    }
    if (plot_type == "violin") {
      p <- p +
        geom_violin(
          scale = "width",
          trim = TRUE,
          alpha = alpha,
          position = position_dodge()
        )
    }
    if (plot_type == "box") {
      add_box <- FALSE
      p <- p +
        geom_boxplot(
          mapping = aes(group = .data[["group.unique"]]),
          position = position_dodge(width = 0.9),
          color = "black",
          width = 0.8,
          outlier.shape = NA
        ) +
        stat_summary(
          fun = stats::median,
          geom = "point",
          mapping = aes(group = .data[["split.by"]]),
          position = position_dodge(width = 0.9),
          color = "black",
          fill = "white",
          size = 1.5,
          shape = 21,
        )
    }
    if (plot_type == "bar") {
      p <- p +
        stat_summary(
          fun = mean,
          geom = "col",
          mapping = aes(group = .data[["split.by"]]),
          position = position_dodge(width = 0.9),
          width = 0.8,
          color = "black"
        ) +
        stat_summary(
          fun.data = mean_sdl,
          fun.args = list(mult = 1),
          geom = "errorbar",
          mapping = aes(group = .data[["split.by"]]),
          position = position_dodge(width = 0.9),
          width = 0.2,
          color = "black"
        )
      y_min_use <- layer_scales(p)$y$range$range[1]
    }
    if (plot_type == "dot") {
      bins <- cut(
        dat$value,
        breaks = seq(min(dat$value), max(dat$value), length.out = 15),
        include.lowest = TRUE
      )
      bins_median <- sapply(
        strsplit(levels(bins), ","),
        function(x) {
          stats::median(
            as.numeric(gsub("\\(|\\)|\\[|\\]", "", x)),
            na.rm = TRUE
          )
        }
      )
      names(bins_median) <- levels(bins)
      dat[["bins"]] <- bins_median[bins]
      p <- p +
        geom_count(
          data = dat,
          aes(y = bins),
          shape = 21,
          alpha = alpha,
          position = position_dodge(width = 0.9)
        ) +
        scale_size_area(name = "Count", max_size = 6, n.breaks = 4) +
        guides(
          size = guide_legend(
            override.aes = list(fill = "grey30", shape = 21),
            order = 2
          )
        )
    }
    if (plot_type == "col") {
      p <- p + geom_col()
      if (flip) {
        p <- p +
          theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      } else {
        p <- p +
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      }
      if (isFALSE(individual) && isTRUE(nlevels(dat[["group.by"]]) > 1)) {
        x_index <- split(dat[["cell"]], dat[["group.by"]])
        border_data <- as.data.frame(sapply(x_index, min) - 0.5)
        colnames(border_data) <- "xintercept"
        border_data <- border_data[2:nrow(border_data), , drop = FALSE]
        border_layer <- geom_vline(
          xintercept = border_data[["xintercept"]],
          linetype = 2,
          alpha = 0.5
        )
        p <- p + border_layer
      }
    }

    if (length(comparisons) > 0) {
      if (isTRUE(comparisons)) {
        group_use <- names(which(
          Matrix::rowSums(table(dat[["group.by"]], dat[["split.by"]]) >= 2) >= 2
        ))
        if (
          any(Matrix::rowSums(table(dat[["group.by"]], dat[["split.by"]]) >= 2) >= 3)
        ) {
          log_message(
            "Detected more than 2 groups. Use multiple_method for comparison"
          )
          method <- multiple_method
        } else {
          method <- pairwise_method
        }
        p <- p +
          ggpubr::stat_compare_means(
            data = dat[dat[["group.by"]] %in% group_use, , drop = FALSE],
            mapping = aes(
              x = .data[["group.by"]],
              y = .data[["value"]],
              group = .data[["group.unique"]]
            ),
            label = sig_label,
            label.y = y_max_use,
            size = sig_labelsize,
            step.increase = 0.1,
            tip.length = 0.03,
            vjust = 1,
            method = method
          )

        y_max_use <- layer_scales(p)$y$range$range[2]
      } else {
        p <- p +
          ggpubr::stat_compare_means(
            mapping = aes(
              x = .data[["group.by"]],
              y = .data[["value"]],
              group = .data[["group.unique"]]
            ),
            label = sig_label,
            label.y = y_max_use,
            size = sig_labelsize,
            step.increase = 0.1,
            tip.length = 0.03,
            vjust = 0,
            comparisons = comparisons,
            ref.group = ref_group,
            method = pairwise_method
          )
        y_max_use <- layer_scales(p)$y$range$range[1] +
          (layer_scales(p)$y$range$range[2] -
            layer_scales(p)$y$range$range[1]) *
            1.15
      }
    }
    if (isTRUE(multiplegroup_comparisons)) {
      p <- p +
        ggpubr::stat_compare_means(
          aes(
            x = .data[["group.by"]],
            y = .data[["value"]],
            group = .data[["group.unique"]]
          ),
          method = multiple_method,
          label = sig_label,
          label.y = y_max_use,
          size = sig_labelsize,
          vjust = 1.2,
          hjust = 0
        )
      y_max_use <- layer_scales(p)$y$range$range[1] +
        (layer_scales(p)$y$range$range[2] - layer_scales(p)$y$range$range[1]) *
          1.15
    }

    if (isTRUE(add_point)) {
      suppressWarnings(
        p <- p +
          geom_point(
            aes(
              x = .data[["group.by"]],
              y = .data[["value"]],
              linetype = rep(f, nrow(dat)),
              group = .data[["group.unique"]]
            ),
            inherit.aes = FALSE,
            color = pt.color,
            size = pt.size,
            alpha = pt.alpha,
            position = position_jitterdodge(
              jitter.width = jitter.width,
              jitter.height = jitter.height,
              dodge.width = 0.9,
              seed = 11
            ),
            show.legend = FALSE
          )
      )
      if (!is.null(cells.highlight)) {
        cell_df <- subset(p$data, rownames(p$data) %in% cells.highlight)
        if (nrow(cell_df) > 0) {
          p <- p +
            geom_point(
              data = cell_df,
              aes(
                x = .data[["group.by"]],
                y = .data[["value"]],
                linetype = rep(f, nrow(cell_df)),
                group = .data[["group.unique"]]
              ),
              inherit.aes = FALSE,
              color = cols.highlight,
              size = sizes.highlight,
              alpha = alpha.highlight,
              position = position_jitterdodge(
                jitter.width = jitter.width,
                jitter.height = jitter.height,
                dodge.width = 0.9,
                seed = 11
              ),
              show.legend = FALSE
            )
        }
      }
    }
    if (isTRUE(add_box)) {
      p <- p +
        geom_boxplot(
          aes(group = .data[["group.unique"]]),
          position = position_dodge(width = 0.9),
          color = box_color,
          fill = box_color,
          width = box_width,
          show.legend = FALSE,
          outlier.shape = NA
        ) +
        stat_summary(
          fun = stats::median,
          geom = "point",
          mapping = aes(group = .data[["split.by"]]),
          position = position_dodge(width = 0.9),
          color = "black",
          fill = "white",
          size = box_ptsize,
          shape = 21,
        )
    }
    if (isTRUE(add_trend)) {
      if (plot_type %in% c("violin", "box")) {
        if (nlevels(dat[["split.by"]]) > 1) {
          point_layer <- stat_summary(
            fun = stats::median,
            geom = "point",
            mapping = aes(
              group = .data[["split.by"]],
              color = .data[["group.by"]]
            ),
            position = position_dodge(width = 0.9),
            fill = "white",
            size = trend_ptsize,
            shape = 21
          )
          p_data <- p + point_layer
          p <- p +
            geom_line(
              data = layer_data(p_data, length(p_data$layers)),
              aes(x = x, y = y, group = colour),
              color = trend_color,
              linewidth = trend_linewidth,
              inherit.aes = FALSE
            ) +
            stat_summary(
              fun = stats::median,
              geom = "point",
              mapping = aes(group = .data[["split.by"]]),
              position = position_dodge(width = 0.9),
              color = "black",
              fill = "white",
              size = trend_ptsize,
              shape = 21
            )
        } else {
          p <- p +
            stat_summary(
              fun = stats::median,
              geom = "line",
              mapping = aes(group = .data[["split.by"]]),
              position = position_dodge(width = 0.9),
              color = trend_color,
              linewidth = trend_linewidth
            ) +
            stat_summary(
              fun = stats::median,
              geom = "point",
              mapping = aes(group = .data[["split.by"]]),
              position = position_dodge(width = 0.9),
              color = "black",
              fill = "white",
              size = trend_ptsize,
              shape = 21
            )
        }
      }
      if (plot_type %in% c("bar")) {
        if (nlevels(dat[["split.by"]]) > 1) {
          point_layer <- stat_summary(
            fun = mean,
            geom = "point",
            mapping = aes(
              group = .data[["split.by"]],
              color = .data[["group.by"]]
            ),
            position = position_dodge(width = 0.9),
            fill = "white",
            size = trend_ptsize,
            shape = 21
          )
          p_data <- p + point_layer
          p <- p +
            geom_line(
              data = layer_data(p_data, length(p_data$layers)),
              aes(x = x, y = y, group = colour),
              color = trend_color,
              linewidth = trend_linewidth,
              inherit.aes = FALSE
            ) +
            stat_summary(
              fun = mean,
              geom = "point",
              mapping = aes(group = .data[["split.by"]]),
              position = position_dodge(width = 0.9),
              color = "black",
              fill = "white",
              size = trend_ptsize,
              shape = 21
            )
        } else {
          p <- p +
            stat_summary(
              fun = mean,
              geom = "line",
              mapping = aes(group = .data[["split.by"]]),
              position = position_dodge(width = 0.9),
              color = trend_color,
              linewidth = trend_linewidth,
            ) +
            stat_summary(
              fun = mean,
              geom = "point",
              mapping = aes(group = .data[["split.by"]]),
              position = position_dodge(width = 0.9),
              color = "black",
              fill = "white",
              size = trend_ptsize,
              shape = 21
            )
        }
      }
    }
    if (add_stat != "none") {
      p <- p +
        stat_summary(
          fun = add_stat,
          geom = "point",
          mapping = aes(group = .data[["split.by"]], shape = stat_shape),
          position = position_dodge(width = 0.9),
          color = stat_color,
          fill = stat_color,
          size = stat_size,
          stroke = stat_stroke,
        ) +
        scale_shape_identity()
    }
    if (!is.null(add_line)) {
      p <- p +
        geom_hline(
          yintercept = add_line,
          color = line_color,
          linetype = line_type,
          linewidth = line_size
        )
    }

    if (nrow(dat) == 0) {
      p <- p + facet_null()
    } else {
      if (isTRUE(stack) && isFALSE(flip)) {
        p <- p +
          facet_grid(features ~ .) +
          theme(strip.text.y = element_text(angle = 0))
      } else {
        p <- p + facet_grid(. ~ features)
      }
    }
    p <- p + labs(title = title, subtitle = subtitle, x = xlab, y = ylab)
    if (nrow(dat) != 0) {
      p <- p + scale_x_discrete(drop = !keep_empty)
    }

    if (isTRUE(flip)) {
      if (isTRUE(stack)) {
        p <- p +
          do.call(theme_use, theme_args) +
          theme(
            aspect.ratio = aspect.ratio,
            axis.text.x = element_text(angle = 90, hjust = 1),
            strip.text.x = element_text(angle = 90),
            panel.grid.major.x = element_line(color = "grey", linetype = 2),
            legend.position = legend.position,
            legend.direction = legend.direction
          ) +
          coord_flip(ylim = c(y_min_use, y_max_use))
      } else {
        p <- p +
          do.call(theme_use, theme_args) +
          theme(
            aspect.ratio = aspect.ratio,
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            strip.text.y = element_text(angle = 0),
            panel.grid.major.x = element_line(color = "grey", linetype = 2),
            legend.position = legend.position,
            legend.direction = legend.direction
          ) +
          coord_flip(ylim = c(y_min_use, y_max_use))
      }
    } else {
      p <- p +
        do.call(theme_use, theme_args) +
        theme(
          aspect.ratio = aspect.ratio,
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          strip.text.y = element_text(angle = 0),
          panel.grid.major.y = element_line(color = "grey", linetype = 2),
          legend.position = legend.position,
          legend.direction = legend.direction
        ) +
        coord_cartesian(ylim = c(y_min_use, y_max_use))
    }

    if (isTRUE(stack)) {
      p <- p +
        scale_y_continuous(
          trans = y.trans,
          breaks = c(y_min_use, y_max_use),
          labels = c(round(y_min_use, 1), round(y_max_use, 1))
        )
    } else {
      p <- p + scale_y_continuous(trans = y.trans, n.breaks = y.nbreaks)
    }

    if (fill.by != "expression") {
      if (isTRUE(stack)) {
        p <- p +
          scale_fill_manual(
            name = paste0(keynm, ":"),
            values = colors,
            breaks = levels_order,
            limits = levels_order,
            drop = FALSE
          ) +
          scale_color_manual(
            name = paste0(keynm, ":"),
            values = colors,
            breaks = levels_order,
            limits = levels_order,
            drop = FALSE
          )
      } else {
        p <- p +
          scale_fill_manual(
            name = paste0(keynm, ":"),
            values = colors,
            breaks = levels_order,
            drop = FALSE
          ) +
          scale_color_manual(
            name = paste0(keynm, ":"),
            values = colors,
            breaks = levels_order,
            drop = FALSE
          )
      }
      p <- p +
        guides(
          fill = guide_legend(
            title.hjust = 0,
            order = 1,
            override.aes = list(size = 4, color = "black", alpha = 1)
          )
        )
    } else {
      p <- p +
        scale_fill_gradientn(
          name = paste0(keynm, ":"),
          colours = colors,
          limits = colors_limits
        ) +
        guides(
          fill = guide_colorbar(
            frame.colour = "black",
            ticks.colour = "black",
            title.hjust = 0,
            order = 1
          )
        )
    }
  })

  return(plist)
}
