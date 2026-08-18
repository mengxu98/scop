scenic_plot_coverage <- function(
  srt,
  tool_name,
  group.by,
  group_annotation,
  group_names,
  features = NULL,
  top_table = NULL,
  atac_assay = NULL,
  extend.upstream = 50000,
  extend.downstream = 50000,
  palette = "RdYlBu",
  palcolor = NULL,
  group_palette = "Chinese",
  group_palcolor = NULL,
  coverage_args = list(),
  title = NULL,
  verbose = TRUE
) {
  links <- scenic_get_region_gene_links(srt, tool_name)
  genes <- scenic_resolve_coverage_genes(
    features = features,
    links = links,
    top_table = top_table
  )
  atac_assay <- scenic_chromatin_assay(srt, atac_assay)
  extend.upstream <- max(0, as.integer(extend.upstream))
  extend.downstream <- max(0, as.integer(extend.downstream))
  fill_colors <- palette_colors(type = "continuous", palette = palette, palcolor = palcolor)
  group_cols <- palette_colors(
    group_names,
    palette = group_palette,
    palcolor = group_palcolor
  )

  plots <- lapply(genes, function(gene) {
    scenic_plot_one_coverage_locus(
      srt = srt,
      gene = gene,
      links = links,
      group.by = group.by,
      group_annotation = group_annotation,
      group_names = group_names,
      atac_assay = atac_assay,
      extend.upstream = extend.upstream,
      extend.downstream = extend.downstream,
      fill_colors = fill_colors,
      group_cols = group_cols,
      group_palette = group_palette,
      group_palcolor = group_palcolor,
      coverage_args = coverage_args,
      title = title,
      verbose = verbose
    )
  })
  names(plots) <- genes
  plot <- if (length(plots) == 1L) {
    plots[[1L]]
  } else {
    patchwork::wrap_plots(plots, ncol = 1)
  }
  list(
    plot = plot,
    plots = plots,
    data = links[links[["gene"]] %in% genes, , drop = FALSE]
  )
}

scenic_get_region_gene_links <- function(srt, tool_name) {
  result <- srt@tools[[tool_name]]
  if (is.null(result)) {
    log_message(
      "Cannot find SCENIC results in tools slot {.val {tool_name}}",
      message_type = "error"
    )
  }
  links <- result[["region_gene"]]
  if (is.null(links) || nrow(links) == 0L) {
    triplets <- result[["triplets"]]
    if (!is.null(triplets) && nrow(triplets) > 0L) {
      links <- as.data.frame(triplets, stringsAsFactors = FALSE)
    }
  } else {
    links <- as.data.frame(links, stringsAsFactors = FALSE)
  }
  if (is.null(links) || nrow(links) == 0L) {
    log_message(
      "{.val coverage} plots need region-gene links in {.code srt@tools[['{tool_name}']]$region_gene} or {.code $triplets}",
      message_type = "error"
    )
  }
  if (!"region" %in% colnames(links)) {
    log_message(
      "Region-gene table must contain a {.val region} column",
      message_type = "error"
    )
  }
  gene_col <- if ("gene" %in% colnames(links)) "gene" else NULL
  if (is.null(gene_col)) {
    log_message(
      "Region-gene table must contain a {.val gene} column",
      message_type = "error"
    )
  }
  score <- if ("score" %in% colnames(links)) {
    as.numeric(links[["score"]])
  } else {
    rep(1, nrow(links))
  }
  tf <- if ("TF" %in% colnames(links)) {
    as.character(links[["TF"]])
  } else {
    rep(NA_character_, nrow(links))
  }
  out <- data.frame(
    region = as.character(links[["region"]]),
    gene = as.character(links[[gene_col]]),
    score = score,
    TF = tf,
    stringsAsFactors = FALSE
  )
  parsed <- scenic_parse_genomic_region(out[["region"]])
  out[["seqnames"]] <- parsed[["seqnames"]]
  out[["start"]] <- parsed[["start"]]
  out[["end"]] <- parsed[["end"]]
  out[["region_mid"]] <- (out[["start"]] + out[["end"]]) / 2
  out <- out[is.finite(out[["start"]]) & is.finite(out[["end"]]), , drop = FALSE]
  if (nrow(out) == 0L) {
    log_message(
      "No parseable genomic regions were found for {.val coverage} plots",
      message_type = "error"
    )
  }
  out
}

scenic_resolve_coverage_genes <- function(features, links, top_table = NULL) {
  genes <- unique(as.character(links[["gene"]]))
  tfs <- unique(as.character(links[["TF"]]))
  tfs <- tfs[!is.na(tfs) & nzchar(tfs)]
  if (is.null(features) || length(features) == 0L) {
    if (!is.null(top_table) && nrow(top_table) > 0L) {
      tf <- scenic_tf_from_regulon(top_table[["regulon"]][[1L]])
      tf_links <- links[links[["TF"]] == tf | links[["gene"]] == tf, , drop = FALSE]
      if (nrow(tf_links) > 0L) {
        tf_links <- tf_links[order(abs(tf_links[["score"]]), decreasing = TRUE), , drop = FALSE]
        return(tf_links[["gene"]][[1L]])
      }
    }
    ranked <- links[order(abs(links[["score"]]), decreasing = TRUE), , drop = FALSE]
    return(ranked[["gene"]][[1L]])
  }
  resolved <- unique(unlist(lapply(as.character(features), function(feature) {
    feature_tf <- scenic_tf_from_regulon(feature)
    hit_genes <- genes[tolower(genes) %in% tolower(c(feature, feature_tf))]
    if (length(hit_genes) > 0L) {
      return(hit_genes)
    }
    hit_tf <- tfs[tolower(tfs) %in% tolower(c(feature, feature_tf))]
    if (length(hit_tf) > 0L) {
      tf_links <- links[links[["TF"]] %in% hit_tf, , drop = FALSE]
      tf_links <- tf_links[order(abs(tf_links[["score"]]), decreasing = TRUE), , drop = FALSE]
      return(unique(tf_links[["gene"]]))
    }
    character(0)
  }), use.names = FALSE))
  if (length(resolved) == 0L) {
    log_message(
      "None of {.arg features} matched region-gene TFs or target genes",
      message_type = "error"
    )
  }
  resolved
}

scenic_chromatin_assay <- function(srt, atac_assay = NULL) {
  assays <- SeuratObject::Assays(srt)
  if (!is.null(atac_assay)) {
    if (!atac_assay %in% assays) {
      log_message(
        "{.arg atac_assay} {.val {atac_assay}} is not present in {.arg srt}",
        message_type = "error"
      )
    }
    return(atac_assay)
  }
  for (nm in c("peaks", assays)) {
    if (nm %in% assays && inherits(srt[[nm]], "ChromatinAssay")) {
      return(nm)
    }
  }
  NULL
}

scenic_parse_genomic_region <- function(region) {
  region <- as.character(region)
  region_clean <- gsub(",", "", region, fixed = TRUE)
  matches <- regexec(
    "^([^:[:space:]-]+)[:_-]([0-9.]+(?:[eE][+-]?[0-9]+)?)[:_-]([0-9.]+(?:[eE][+-]?[0-9]+)?)$",
    region_clean
  )
  parts <- regmatches(region_clean, matches)
  out <- data.frame(
    region = region,
    seqnames = NA_character_,
    start = NA_real_,
    end = NA_real_,
    stringsAsFactors = FALSE
  )
  ok <- lengths(parts) == 4L
  if (any(ok)) {
    parsed <- do.call(rbind, parts[ok])
    out[["seqnames"]][ok] <- parsed[, 2]
    out[["start"]][ok] <- as.numeric(parsed[, 3])
    out[["end"]][ok] <- as.numeric(parsed[, 4])
    swap <- !is.na(out[["start"]]) & !is.na(out[["end"]]) & out[["start"]] > out[["end"]]
    if (any(swap)) {
      tmp <- out[["start"]][swap]
      out[["start"]][swap] <- out[["end"]][swap]
      out[["end"]][swap] <- tmp
    }
  }
  out
}

scenic_annotation_gene_coords <- function(srt, atac_assay, genes) {
  if (is.null(atac_assay) || !atac_assay %in% SeuratObject::Assays(srt)) {
    return(NULL)
  }
  if (!inherits(srt[[atac_assay]], "ChromatinAssay")) {
    return(NULL)
  }
  ann <- tryCatch(Signac::Annotation(srt[[atac_assay]]), error = function(...) NULL)
  if (is.null(ann) || length(ann) == 0L) {
    return(NULL)
  }
  ann_df <- as.data.frame(ann, stringsAsFactors = FALSE)
  name_col <- scenic_first(intersect(c("gene_name", "gene_id"), colnames(ann_df)))
  if (is.null(name_col)) {
    return(NULL)
  }
  keep <- tolower(as.character(ann_df[[name_col]])) %in% tolower(genes)
  if (!any(keep)) {
    return(NULL)
  }
  ann_df <- ann_df[keep, , drop = FALSE]
  data.frame(
    gene = as.character(ann_df[[name_col]]),
    seqnames = as.character(ann_df[["seqnames"]]),
    start = as.numeric(ann_df[["start"]]),
    end = as.numeric(ann_df[["end"]]),
    strand = as.character(ann_df[["strand"]]),
    stringsAsFactors = FALSE
  )
}

scenic_infer_gene_coords <- function(gene_links) {
  do.call(rbind, lapply(split(gene_links, gene_links[["gene"]]), function(df) {
    data.frame(
      gene = df[["gene"]][[1L]],
      seqnames = df[["seqnames"]][[1L]],
      start = min(df[["start"]], na.rm = TRUE),
      end = max(df[["end"]], na.rm = TRUE),
      strand = "*",
      stringsAsFactors = FALSE
    )
  }))
}

scenic_has_fragments <- function(srt, atac_assay) {
  if (is.null(atac_assay)) {
    return(FALSE)
  }
  tryCatch(
    length(Signac::Fragments(srt[[atac_assay]])) > 0L,
    error = function(...) FALSE
  )
}

scenic_peak_coverage_data <- function(
  srt,
  atac_assay,
  window,
  group_annotation,
  group_names
) {
  peaks <- rownames(srt[[atac_assay]])
  peak_gr <- scenic_parse_genomic_region(peaks)
  keep <- !is.na(peak_gr[["seqnames"]]) &
    peak_gr[["seqnames"]] == window[["seqnames"]] &
    peak_gr[["end"]] >= window[["start"]] &
    peak_gr[["start"]] <= window[["end"]]
  if (!any(keep)) {
    return(NULL)
  }
  peak_gr <- peak_gr[keep, , drop = FALSE]
  mat <- scenic_get_expression_matrix(srt, assay = atac_assay, layer = "counts")
  if (nrow(mat) == 0L) {
    mat <- scenic_get_expression_matrix(srt, assay = atac_assay, layer = "data")
  }
  peak_ids <- intersect(peak_gr[["region"]], rownames(mat))
  if (length(peak_ids) == 0L) {
    peak_ids <- intersect(
      gsub(":", "-", peak_gr[["region"]], fixed = TRUE),
      rownames(mat)
    )
  }
  if (length(peak_ids) == 0L) {
    return(NULL)
  }
  cells <- intersect(colnames(mat), names(group_annotation))
  if (length(cells) == 0L) {
    return(NULL)
  }
  mat <- mat[peak_ids, cells, drop = FALSE]
  avg <- scenic_group_average_matrix(
    auc_mat = mat,
    group_annotation = group_annotation[cells],
    group_names = group_names
  )
  long <- scenic_matrix_to_long(avg, "region", "group", "accessibility")
  coords <- peak_gr
  coords[["match"]] <- coords[["region"]]
  alt <- gsub(":", "-", coords[["region"]], fixed = TRUE)
  long <- merge(
    long,
    coords[, c("match", "seqnames", "start", "end")],
    by.x = "region",
    by.y = "match",
    all.x = TRUE
  )
  if (any(is.na(long[["start"]]))) {
    long_alt <- merge(
      scenic_matrix_to_long(avg, "region", "group", "accessibility"),
      data.frame(
        match = alt,
        seqnames = coords[["seqnames"]],
        start = coords[["start"]],
        end = coords[["end"]],
        stringsAsFactors = FALSE
      ),
      by.x = "region",
      by.y = "match",
      all.x = TRUE
    )
    fill <- is.na(long[["start"]])
    long[fill, c("seqnames", "start", "end")] <- long_alt[
      match(long[["region"]][fill], long_alt[["region"]]),
      c("seqnames", "start", "end")
    ]
  }
  long[["group"]] <- factor(long[["group"]], levels = group_names)
  long[!is.na(long[["start"]]), , drop = FALSE]
}

scenic_try_coverage_track_plot <- function(
  srt,
  atac_assay,
  window,
  group.by,
  group_palette,
  group_palcolor,
  coverage_args,
  verbose
) {
  if (!scenic_has_fragments(srt, atac_assay)) {
    return(NULL)
  }
  locus <- paste0(window[["seqnames"]], "-", window[["start"]], "-", window[["end"]])
  args <- list(
    srt = srt,
    region = locus,
    assay = atac_assay,
    group.by = group.by,
    palette = group_palette,
    palcolor = group_palcolor,
    extend.upstream = 0,
    extend.downstream = 0,
    annotation = TRUE,
    peaks = TRUE,
    links = TRUE,
    verbose = verbose
  )
  extra <- coverage_args[setdiff(names(coverage_args), names(args))]
  args[names(extra)] <- extra
  tryCatch(
    do.call(CoverageTrackPlot, args),
    error = function(...) NULL
  )
}

scenic_set_signac_links <- function(srt, atac_assay, locus_links, gene_coords) {
  if (is.null(atac_assay) || !inherits(srt[[atac_assay]], "ChromatinAssay")) {
    return(invisible(NULL))
  }
  if (!requireNamespace("GenomicRanges", quietly = TRUE) ||
    !requireNamespace("IRanges", quietly = TRUE)) {
    return(invisible(NULL))
  }
  gene_pos <- stats::setNames(
    (gene_coords[["start"]] + gene_coords[["end"]]) / 2,
    gene_coords[["gene"]]
  )
  keep <- locus_links[["gene"]] %in% names(gene_pos)
  if (!any(keep)) {
    return(invisible(NULL))
  }
  df <- locus_links[keep, , drop = FALSE]
  df[["gene_pos"]] <- unname(gene_pos[df[["gene"]]])
  gr <- GenomicRanges::GRanges(
    seqnames = df[["seqnames"]],
    ranges = IRanges::IRanges(
      start = pmin(df[["region_mid"]], df[["gene_pos"]]),
      end = pmax(df[["region_mid"]], df[["gene_pos"]])
    ),
    score = df[["score"]]
  )
  tryCatch(
    Signac::Links(srt[[atac_assay]]) <- gr,
    error = function(...) NULL
  )
  invisible(NULL)
}

scenic_coverage_theme <- function() {
  theme_scop() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

scenic_plot_one_coverage_locus <- function(
  srt,
  gene,
  links,
  group.by,
  group_annotation,
  group_names,
  atac_assay,
  extend.upstream,
  extend.downstream,
  fill_colors,
  group_cols,
  group_palette = "Chinese",
  group_palcolor = NULL,
  coverage_args,
  title,
  verbose
) {
  locus_links <- links[links[["gene"]] == gene, , drop = FALSE]
  if (nrow(locus_links) == 0L) {
    log_message(
      "No region-gene links remain for {.val {gene}}",
      message_type = "error"
    )
  }
  gene_coords <- scenic_annotation_gene_coords(srt, atac_assay, gene)
  if (is.null(gene_coords) || nrow(gene_coords) == 0L) {
    gene_coords <- scenic_infer_gene_coords(locus_links)
    span <- max(500L, as.integer(diff(range(c(gene_coords[["start"]], gene_coords[["end"]]))) * 0.05))
    mid <- (gene_coords[["start"]] + gene_coords[["end"]]) / 2
    gene_coords[["start"]] <- mid - span
    gene_coords[["end"]] <- mid + span
  } else {
    gene_coords <- data.frame(
      gene = gene,
      seqnames = gene_coords[["seqnames"]][[1L]],
      start = min(gene_coords[["start"]], na.rm = TRUE),
      end = max(gene_coords[["end"]], na.rm = TRUE),
      strand = gene_coords[["strand"]][[1L]],
      stringsAsFactors = FALSE
    )
  }
  chr <- unique(c(locus_links[["seqnames"]], gene_coords[["seqnames"]]))[[1L]]
  locus_links <- locus_links[locus_links[["seqnames"]] == chr, , drop = FALSE]
  window <- list(
    seqnames = chr,
    start = min(c(locus_links[["start"]], gene_coords[["start"]]), na.rm = TRUE) - extend.upstream,
    end = max(c(locus_links[["end"]], gene_coords[["end"]]), na.rm = TRUE) + extend.downstream
  )
  window[["start"]] <- max(1, as.integer(window[["start"]]))
  window[["end"]] <- as.integer(window[["end"]])
  x_limits <- c(window[["start"]], window[["end"]])
  gene_coords[["tss"]] <- ifelse(
    gene_coords[["strand"]] == "-",
    gene_coords[["end"]],
    gene_coords[["start"]]
  )
  locus_links[["gene_pos"]] <- gene_coords[["tss"]][[1L]]

  log_message(
    "Drawing coverage / region-gene track for {.val {gene}} on {.val {chr}}:{window$start}-{window$end}",
    verbose = verbose
  )

  scenic_set_signac_links(srt, atac_assay, locus_links, gene_coords)
  signac_plot <- scenic_try_coverage_track_plot(
    srt = srt,
    atac_assay = atac_assay,
    window = window,
    group.by = group.by,
    group_palette = group_palette,
    group_palcolor = group_palcolor,
    coverage_args = coverage_args,
    verbose = verbose
  )

  tracks <- list()
  cov_data <- NULL
  if (!is.null(atac_assay) && is.null(signac_plot)) {
    cov_data <- scenic_peak_coverage_data(
      srt = srt,
      atac_assay = atac_assay,
      window = window,
      group_annotation = group_annotation,
      group_names = group_names
    )
  }
  if (!is.null(cov_data) && nrow(cov_data) > 0L) {
    tracks[["coverage"]] <- ggplot2::ggplot(cov_data) +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = .data[["start"]],
          xmax = .data[["end"]],
          ymin = 0,
          ymax = .data[["accessibility"]],
          fill = .data[["group"]]
        ),
        color = NA
      ) +
      ggplot2::scale_fill_manual(values = group_cols, guide = "none") +
      ggplot2::facet_grid(rows = ggplot2::vars(.data[["group"]]), scales = "free_y") +
      ggplot2::scale_x_continuous(limits = x_limits, expand = c(0.01, 0)) +
      scenic_coverage_theme() +
      ggplot2::labs(y = "ATAC", title = title %||% paste(gene, "locus"))
  } else if (is.null(signac_plot)) {
    tracks[["coverage"]] <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = mean(x_limits),
        y = 0.5,
        label = "No chromatin assay / peak coverage in this window"
      ) +
      ggplot2::scale_x_continuous(limits = x_limits, expand = c(0.01, 0)) +
      ggplot2::scale_y_continuous(limits = c(0, 1), breaks = NULL) +
      scenic_coverage_theme() +
      ggplot2::labs(y = "ATAC", title = title %||% paste(gene, "locus"))
  }

  peak_df <- unique(locus_links[, c("region", "start", "end", "score", "TF"), drop = FALSE])
  tracks[["peaks"]] <- ggplot2::ggplot(peak_df) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = .data[["start"]],
        xmax = .data[["end"]],
        ymin = 0.25,
        ymax = 0.75,
        fill = .data[["score"]]
      ),
      color = "black",
      linewidth = 0.2
    ) +
    scenic_fill_scale(fill_colors, "Score") +
    ggplot2::scale_x_continuous(limits = x_limits, expand = c(0.01, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = NULL) +
    scenic_coverage_theme() +
    ggplot2::labs(y = "Peaks")

  tracks[["links"]] <- ggplot2::ggplot(locus_links) +
    ggplot2::geom_curve(
      ggplot2::aes(
        x = .data[["region_mid"]],
        xend = .data[["gene_pos"]],
        y = 0,
        yend = 0,
        color = .data[["score"]],
        linewidth = abs(.data[["score"]])
      ),
      curvature = -0.45,
      ncp = 20
    ) +
    scenic_color_scale(fill_colors, "Score") +
    ggplot2::scale_linewidth(range = c(0.3, 1.4), guide = "none") +
    ggplot2::scale_x_continuous(limits = x_limits, expand = c(0.01, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = NULL) +
    scenic_coverage_theme() +
    ggplot2::labs(y = "Links")

  tracks[["genes"]] <- ggplot2::ggplot(gene_coords) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data[["start"]], xend = .data[["end"]], y = 0.5, yend = 0.5),
      linewidth = 2,
      color = "grey20"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data[["tss"]], y = 0.5),
      size = 2.5,
      color = "grey20"
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(x = (.data[["start"]] + .data[["end"]]) / 2, y = 0.75, label = .data[["gene"]]),
      size = 3,
      fontface = "italic",
      max.overlaps = Inf
    ) +
    ggplot2::scale_x_continuous(
      limits = x_limits,
      expand = c(0.01, 0),
      labels = function(x) format(x, big.mark = ",", scientific = FALSE)
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = NULL) +
    theme_scop() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(colour = "black")
    ) +
    ggplot2::labs(x = window[["seqnames"]], y = "Genes")

  ggplot_tracks <- Filter(Negate(is.null), tracks)
  heights <- c(
    coverage = 2.4,
    peaks = 0.55,
    links = 0.9,
    genes = 0.7
  )[names(ggplot_tracks)]
  ggplot_plot <- patchwork::wrap_plots(
    ggplot_tracks,
    ncol = 1,
    heights = unname(heights)
  )
  if (is.null(signac_plot)) {
    return(ggplot_plot)
  }
  patchwork::wrap_plots(signac_plot, ggplot_plot, ncol = 1, heights = c(2.5, 2))
}
