#' Shorten and offset the segment
#'
#' This function takes a data frame representing segments in a plot and shortens
#' and offsets them based on the provided arguments.
#'
#' @param data A data frame containing the segments.
#' It should have columns 'x', 'y', 'xend', and 'yend' representing the start and end points of each segment.
#' @param shorten_start The amount to shorten the start of each segment by.
#' @param shorten_end The amount to shorten the end of each segment by.
#' @param offset The amount to offset each segment by.
#'
#' @return The modified data frame with the shortened and offset segments.
#'
#' @export
#'
#' @examples
#' library(ggplot2)
#' temp_nodes <- data.frame(
#'   "x" = c(10, 40),
#'   "y" = c(10, 30)
#' )
#' data <- data.frame(
#'   "x" = c(10, 40),
#'   "y" = c(10, 30),
#'   "xend" = c(40, 10),
#'   "yend" = c(30, 10)
#' )
#'
#' ggplot(temp_nodes, aes(x = x, y = y)) +
#'   geom_point(size = 12) +
#'   xlim(0, 50) +
#'   ylim(0, 50) +
#'   geom_segment(
#'     data = data,
#'     aes(x = x, xend = xend, y = y, yend = yend)
#'   )
#'
#' ggplot(temp_nodes, aes(x = x, y = y)) +
#'   geom_point(size = 12) +
#'   xlim(0, 50) +
#'   ylim(0, 50) +
#'   geom_segment(
#'     data = segements_df(
#'       data,
#'       shorten_start = 2,
#'       shorten_end = 3,
#'       offset = 1
#'     ),
#'     aes(x = x, xend = xend, y = y, yend = yend)
#'   )
segements_df <- function(
    data,
    shorten_start,
    shorten_end, offset) {
  data$dx <- data$xend - data$x
  data$dy <- data$yend - data$y
  data$dist <- sqrt(data$dx^2 + data$dy^2)
  data$px <- data$dx / data$dist
  data$py <- data$dy / data$dist

  data$x <- data$x + data$px * shorten_start
  data$y <- data$y + data$py * shorten_start
  data$xend <- data$xend - data$px * shorten_end
  data$yend <- data$yend - data$py * shorten_end
  data$x <- data$x - data$py * offset
  data$xend <- data$xend - data$py * offset
  data$y <- data$y + data$px * offset
  data$yend <- data$yend + data$px * offset

  return(data)
}

fc_matrix <- function(matrix) {
  matrix / rowMeans(matrix)
}

zscore_matrix <- function(matrix, ...) {
  Matrix::t(scale(Matrix::t(matrix), ...))
}

log2fc_matrix <- function(matrix) {
  log2(matrix / rowMeans(matrix))
}

log1p_matrix <- function(matrix) {
  log1p(matrix)
}

matrix_process <- function(
    matrix,
    method = c("raw", "zscore", "fc", "log2fc", "log1p"),
    ...) {
  if (is.function(method)) {
    matrix_processed <- method(matrix, ...)
  } else if (method == "raw") {
    matrix_processed <- matrix
  } else if (method == "fc") {
    matrix_processed <- fc_matrix(matrix)
  } else if (method == "zscore") {
    matrix_processed <- zscore_matrix(matrix, ...)
  } else if (method == "log2fc") {
    matrix_processed <- log2fc_matrix(matrix)
  } else if (method == "log1p") {
    matrix_processed <- log1p_matrix(matrix)
  }
  if (!identical(dim(matrix_processed), dim(matrix))) {
    stop("The dimensions of the matrix are changed after processing")
  }
  return(matrix_processed)
}

extractgrobs <- function(vlnplots, x_nm, y_nm, x, y) {
  grobs <- vlnplots[paste0(x_nm[x], ":", y_nm[y])]
  if (length(grobs) == 1) {
    grobs <- grobs[[1]]
  }
  return(grobs)
}

grid_draw <- function(groblist, x, y, width, height) {
  if (grid::is.grob(groblist)) {
    groblist <- list(groblist)
  }
  for (i in seq_along(groblist)) {
    groblist[[i]]$vp <- grid::viewport(
      x = x[i],
      y = y[i],
      width = width[i],
      height = height[i]
    )
    grid::grid.draw(groblist[[i]])
  }
}

cluster_within_group2 <- function(mat, factor) {
  check_r("dendextend")
  if (!is.factor(factor)) {
    factor <- factor(factor, levels = unique(factor))
  }
  dend_list <- list()
  order_list <- list()
  for (le in unique(levels(factor))) {
    m <- mat[, factor == le, drop = FALSE]
    if (ncol(m) == 1) {
      order_list[[le]] <- which(factor == le)
      dend_list[[le]] <- structure(
        which(factor == le),
        class = "dendrogram",
        leaf = TRUE, # height = 0,
        label = 1,
        members = 1
      )
    } else if (ncol(m) > 1) {
      hc1 <- stats::hclust(
        stats::as.dist(
          proxyC::dist(Matrix::t(m))
        )
      )
      dend_list[[le]] <- stats::as.dendrogram(hc1)
      order_list[[le]] <- which(factor == le)[stats::order.dendrogram(dend_list[[le]])]
      stats::order.dendrogram(dend_list[[le]]) <- order_list[[le]]
    }
    attr(dend_list[[le]], ".class_label") <- le
  }
  parent <- stats::as.dendrogram(
    stats::hclust(
      stats::as.dist(
        proxyC::dist(
          Matrix::t(sapply(
            order_list,
            function(x) rowMeans(mat[, x, drop = FALSE])
          ))
        )
      )
    )
  )
  dend_list <- lapply(dend_list, function(dend) {
    stats::dendrapply(
      dend,
      function(node) {
        if (is.null(attr(node, "height"))) {
          attr(node, "height") <- 0
        }
        node
      }
    )
  })

  dend <- ComplexHeatmap::merge_dendrogram(parent, dend_list)
  stats::order.dendrogram(dend) <- unlist(
    order_list[stats::order.dendrogram(parent)])
  return(dend)
}

heatmap_enrichment <- function(
    geneID,
    geneID_groups,
    feature_split_palette = "simspec",
    feature_split_palcolor = NULL,
    ha_right = NULL,
    flip = FALSE,
    anno_terms = FALSE,
    anno_keys = FALSE,
    anno_features = FALSE,
    terms_width = grid::unit(4, "in"),
    terms_fontsize = 8,
    keys_width = grid::unit(2, "in"),
    keys_fontsize = c(6, 10),
    features_width = grid::unit(2, "in"),
    features_fontsize = c(6, 10),
    IDtype = "symbol",
    species = "Homo_sapiens",
    db_update = FALSE,
    db_combine = FALSE,
    db_version = "latest",
    convert_species = FALSE,
    Ensembl_version = 103,
    mirror = NULL,
    db = "GO_BP",
    TERM2GENE = NULL,
    TERM2NAME = NULL,
    minGSSize = 10,
    maxGSSize = 500,
    GO_simplify = FALSE,
    GO_simplify_cutoff = "p.adjust < 0.05",
    simplify_method = "Wang",
    simplify_similarityCutoff = 0.7,
    pvalueCutoff = NULL,
    padjustCutoff = 0.05,
    topTerm = 5,
    show_termid = FALSE,
    topWord = 20,
    words_excluded = NULL) {
  res <- NULL
  words_excluded <- words_excluded %||% scop::words_excluded

  if (isTRUE(anno_keys) || isTRUE(anno_features) || isTRUE(anno_terms)) {
    if (isTRUE(flip)) {
      stop(
        "anno_keys, anno_features and anno_terms can only be used when flip is FALSE."
      )
    }
    if (all(is.na(geneID_groups))) {
      geneID_groups <- rep(1, length(geneID))
    }
    if (!is.factor(geneID_groups)) {
      geneID_groups <- factor(geneID_groups, levels = unique(geneID_groups))
    }
    fill_split <- palette_scop(
      levels(geneID_groups),
      type = "discrete",
      palette = feature_split_palette,
      palcolor = feature_split_palcolor
    )[levels(geneID_groups) %in% geneID_groups]
    res <- RunEnrichment(
      geneID = geneID,
      geneID_groups = geneID_groups,
      IDtype = IDtype,
      species = species,
      db_update = db_update,
      db_version = db_version,
      db_combine = db_combine,
      convert_species = convert_species,
      Ensembl_version = Ensembl_version,
      mirror = mirror,
      db = db,
      TERM2GENE = TERM2GENE,
      TERM2NAME = TERM2NAME,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      GO_simplify = GO_simplify,
      GO_simplify_cutoff = GO_simplify_cutoff,
      simplify_method = simplify_method,
      simplify_similarityCutoff = simplify_similarityCutoff
    )
    if (isTRUE(db_combine)) {
      db <- "Combined"
    }
    if (isTRUE(GO_simplify) && any(db %in% c("GO_BP", "GO_CC", "GO_MF"))) {
      db[db %in% c("GO_BP", "GO_CC", "GO_MF")] <- paste0(
        db[db %in% c("GO_BP", "GO_CC", "GO_MF")],
        "_sim"
      )
    }
    if (nrow(res$enrichment) == 0) {
      warning("No enrichment result found.", immediate. = TRUE)
    } else {
      metric <- ifelse(is.null(padjustCutoff), "pvalue", "p.adjust")
      metric_value <- ifelse(
        is.null(padjustCutoff),
        pvalueCutoff,
        padjustCutoff
      )
      pvalueCutoff <- ifelse(is.null(pvalueCutoff), 1, pvalueCutoff)
      padjustCutoff <- ifelse(is.null(padjustCutoff), 1, padjustCutoff)

      df <- res$enrichment
      df <- df[df[["Database"]] %in% db, , drop = FALSE]
      df <- df[df[[metric]] < metric_value, , drop = FALSE]
      df <- df[order(df[[metric]]), , drop = FALSE]
      if (nrow(df) == 0) {
        warning(
          "No term enriched using the threshold: ",
          paste0("pvalueCutoff = ", pvalueCutoff),
          "; ",
          paste0("padjustCutoff = ", padjustCutoff),
          immediate. = TRUE
        )
      } else {
        df_list <- split.data.frame(df, ~ Database + Groups)
        df_list <- df_list[lapply(df_list, nrow) > 0]

        for (enrich in db) {
          nm <- strsplit(names(df_list), "\\.")
          subdf_list <- df_list[
            unlist(lapply(nm, function(x) x[[1]])) %in% enrich
          ]
          if (length(subdf_list) == 0) {
            warning(
              "No ",
              enrich,
              " term enriched using the threshold: ",
              paste0("pvalueCutoff = ", pvalueCutoff),
              "; ",
              paste0("padjustCutoff = ", padjustCutoff),
              immediate. = TRUE
            )
            next
          }
          nm <- strsplit(names(subdf_list), "\\.")

          ha_terms <- NULL
          if (isTRUE(anno_terms)) {
            terms_list <- lapply(subdf_list, function(df) {
              if (isTRUE(show_termid)) {
                terms <- paste(
                  utils::head(df$ID, topTerm),
                  utils::head(df$Description, topTerm)
                )
              } else {
                terms <- utils::head(df$Description, topTerm)
                terms <- capitalize(terms)
              }
              df_out <- data.frame(keyword = terms)
              df_out[["col"]] <- palette_scop(
                -log10(utils::head(df[, metric], topTerm)),
                type = "continuous",
                palette = "Spectral",
                matched = TRUE
              )
              df_out[["col"]] <- sapply(
                df_out[["col"]],
                function(x) blendcolors(c(x, "black"))
              )
              df_out[["fontsize"]] <- rep(terms_fontsize, nrow(df_out))
              return(df_out)
            })
            names(terms_list) <- unlist(lapply(nm, function(x) x[[2]]))
            if (length(intersect(geneID_groups, names(terms_list))) > 0) {
              ha_terms <- ComplexHeatmap::HeatmapAnnotation(
                "terms_empty" = ComplexHeatmap::anno_empty(
                  width = grid::unit(0.05, "in"),
                  border = FALSE,
                  which = "row"
                ),
                "terms_split" = ComplexHeatmap::anno_block(
                  gp = grid::gpar(fill = fill_split),
                  width = grid::unit(0.1, "in"),
                  which = "row"
                ),
                "terms" = ComplexHeatmap::anno_textbox(
                  align_to = geneID_groups,
                  text = terms_list,
                  max_width = terms_width,
                  word_wrap = TRUE,
                  add_new_line = TRUE,
                  background_gp = grid::gpar(fill = "grey98", col = "black"),
                  round_corners = TRUE,
                  which = "row"
                ),
                which = "row",
                gap = grid::unit(0, "points")
              )
              names(ha_terms) <- paste0(names(ha_terms), "_", enrich)
            }
          }

          ha_keys <- NULL
          if (isTRUE(anno_keys)) {
            check_r("simplifyEnrichment")
            keys_list <- lapply(subdf_list, function(df) {
              if (all(df$Database %in% c("GO", "GO_BP", "GO_CC", "GO_MF"))) {
                df0 <- simplifyEnrichment::keyword_enrichment_from_GO(df[[
                  "ID"
                ]])
                if (nrow(df0) > 0) {
                  df <- df0 %>%
                    dplyr::reframe(
                      keyword = .data[["keyword"]],
                      score = -(log10(.data[["padj"]])),
                      count = .data[["n_term"]],
                      Database = df[["Database"]][1],
                      Groups = df[["Groups"]][1]
                    ) %>%
                    dplyr::filter(
                      !grepl(pattern = "\\[.*\\]", x = .data[["keyword"]])
                    ) %>%
                    dplyr::filter(nchar(.data[["keyword"]]) >= 1) %>%
                    dplyr::filter(
                      !tolower(.data[["keyword"]]) %in% tolower(words_excluded)
                    ) %>%
                    dplyr::distinct() %>%
                    dplyr::mutate(
                      angle = 90 *
                        sample(c(0, 1), dplyr::n(), replace = TRUE, prob = c(60, 40))
                    ) %>%
                    as.data.frame()
                  df <- df[
                    utils::head(order(df[["score"]], decreasing = TRUE), topWord), ,
                    drop = FALSE
                  ]
                } else {
                  df <- NULL
                }
              } else {
                df <- df %>%
                  dplyr::mutate(
                    keyword = strsplit(
                      tolower(as.character(.data[["Description"]])),
                      " "
                    )
                  ) %>%
                  unnest(cols = "keyword") %>%
                  dplyr::group_by(.data[["keyword"]], Database, Groups) %>%
                  dplyr::reframe(
                    keyword = .data[["keyword"]],
                    score = sum(-(log10(.data[[metric]]))),
                    count = dplyr::n(),
                    Database = .data[["Database"]],
                    Groups = .data[["Groups"]],
                    .groups = "keep"
                  ) %>%
                  dplyr::filter(
                    !grepl(pattern = "\\[.*\\]", x = .data[["keyword"]])
                  ) %>%
                  dplyr::filter(nchar(.data[["keyword"]]) >= 1) %>%
                  dplyr::filter(
                    !tolower(.data[["keyword"]]) %in% tolower(words_excluded)
                  ) %>%
                  dplyr::distinct() %>%
                  dplyr::mutate(
                    angle = 90 *
                      sample(c(0, 1), dplyr::n(), replace = TRUE, prob = c(60, 40))
                  ) %>%
                  as.data.frame()
                df <- df[
                  utils::head(order(df[["score"]], decreasing = TRUE), topWord), ,
                  drop = FALSE
                ]
              }
              if (isTRUE(nrow(df) > 0)) {
                df[["col"]] <- palette_scop(
                  df[, "score"],
                  type = "continuous",
                  palette = "Spectral",
                  matched = TRUE
                )
                df[["col"]] <- sapply(
                  df[["col"]],
                  function(x) blendcolors(c(x, "black"))
                )
                df[["fontsize"]] <- rescale(df[, "count"], to = keys_fontsize)
                return(df)
              } else {
                return(NULL)
              }
            })
            names(keys_list) <- unlist(lapply(nm, function(x) x[[2]]))
            keys_list <- keys_list[lapply(keys_list, length) > 0]
            if (length(intersect(geneID_groups, names(keys_list))) > 0) {
              ha_keys <- ComplexHeatmap::HeatmapAnnotation(
                "keys_empty" = ComplexHeatmap::anno_empty(
                  width = grid::unit(0.05, "in"),
                  border = FALSE,
                  which = "row"
                ),
                "keys_split" = ComplexHeatmap::anno_block(
                  gp = grid::gpar(fill = fill_split),
                  width = grid::unit(0.1, "in"),
                  which = "row"
                ),
                "keys" = ComplexHeatmap::anno_textbox(
                  align_to = geneID_groups,
                  text = keys_list,
                  max_width = keys_width,
                  background_gp = grid::gpar(fill = "grey98", col = "black"),
                  round_corners = TRUE,
                  which = "row"
                ),
                which = "row",
                gap = grid::unit(0, "points")
              )
              names(ha_keys) <- paste0(names(ha_keys), "_", enrich)
            }
          }

          ha_features <- NULL
          if (isTRUE(anno_features)) {
            features_list <- lapply(subdf_list, function(df) {
              df <- df %>%
                dplyr::mutate(
                  keyword = strsplit(as.character(.data[["geneID"]]), "/")
                ) %>%
                unnest(cols = "keyword") %>%
                dplyr::group_by(.data[["keyword"]], Database, Groups) %>%
                dplyr::reframe(
                  keyword = .data[["keyword"]],
                  score = sum(-(log10(.data[[metric]]))),
                  count = dplyr::n(),
                  Database = .data[["Database"]],
                  Groups = .data[["Groups"]],
                  .groups = "keep"
                ) %>%
                dplyr::distinct() %>%
                dplyr::mutate(
                  angle = 90 *
                    sample(c(0, 1), dplyr::n(), replace = TRUE, prob = c(60, 40))
                ) %>%
                as.data.frame()
              df <- df[
                utils::head(order(df[["score"]], decreasing = TRUE), topWord), ,
                drop = FALSE
              ]
              df[["col"]] <- palette_scop(
                df[, "score"],
                type = "continuous",
                palette = "Spectral",
                matched = TRUE
              )
              df[["col"]] <- sapply(
                df[["col"]],
                function(x) blendcolors(c(x, "black"))
              )
              df[["fontsize"]] <- rescale(df[, "count"], to = features_fontsize)
              return(df)
            })
            names(features_list) <- unlist(lapply(nm, function(x) x[[2]]))
            if (length(intersect(geneID_groups, names(features_list))) > 0) {
              ha_features <- ComplexHeatmap::HeatmapAnnotation(
                "features_empty" = ComplexHeatmap::anno_empty(
                  width = grid::unit(0.05, "in"),
                  border = FALSE,
                  which = "row"
                ),
                "features_split" = ComplexHeatmap::anno_block(
                  gp = grid::gpar(fill = fill_split),
                  width = grid::unit(0.1, "in"),
                  which = "row"
                ),
                "features" = ComplexHeatmap::anno_textbox(
                  align_to = geneID_groups,
                  text = features_list,
                  max_width = features_width,
                  background_gp = grid::gpar(fill = "grey98", col = "black"),
                  round_corners = TRUE,
                  which = "row"
                ),
                which = "row",
                gap = grid::unit(0, "points")
              )
              names(ha_features) <- paste0(names(ha_features), "_", enrich)
            }
          }

          ha_enrichment <- list(ha_terms, ha_keys, ha_features)
          ha_enrichment <- ha_enrichment[sapply(ha_enrichment, length) > 0]
          ha_enrichment <- do.call(c, ha_enrichment)

          if (is.null(ha_right)) {
            ha_right <- ha_enrichment
          } else {
            ha_right <- c(ha_right, ha_enrichment)
          }
        }
      }
    }
  }
  return(list(ha_right = ha_right, res = res))
}

heatmap_rendersize <- function(
    width,
    height,
    units,
    ha_top_list,
    ha_left,
    ha_right,
    ht_list,
    legend_list,
    flip) {
  width_annotation <- height_annotation <- 0
  if (isTRUE(flip)) {
    width_sum <- width[1] %||%
      grid::convertWidth(
        grid::unit(1, "in"),
        units,
        valueOnly = TRUE
      )
    height_sum <- sum(
      height %||% grid::convertHeight(
        grid::unit(1, "in"),
        units,
        valueOnly = TRUE
      )
    )
    if (length(ha_top_list) > 0) {
      width_annotation <- grid::convertWidth(
        grid::unit(width_annotation, units) +
          ComplexHeatmap::width.HeatmapAnnotation(ha_top_list[[1]]),
        units,
        valueOnly = TRUE
      )
    }
    if (!is.null(ha_left)) {
      height_annotation <- grid::convertHeight(
        grid::unit(height_annotation, units) + ComplexHeatmap::height.HeatmapAnnotation(ha_left),
        units,
        valueOnly = TRUE
      )
    }
    if (!is.null(ha_right)) {
      height_annotation <- grid::convertHeight(
        grid::unit(height_annotation, units) + ComplexHeatmap::height.HeatmapAnnotation(ha_right),
        units,
        valueOnly = TRUE
      )
    }
  } else {
    width_sum <- sum(
      width %||% grid::convertWidth(
        grid::unit(1, "in"),
        units,
        valueOnly = TRUE
      )
    )
    height_sum <- height[1] %||%
      grid::convertHeight(grid::unit(1, "in"), units, valueOnly = TRUE)
    if (length(ha_top_list) > 0) {
      height_annotation <- grid::convertHeight(
        grid::unit(height_annotation, units) +
          ComplexHeatmap::height.HeatmapAnnotation(ha_top_list[[1]]),
        units,
        valueOnly = TRUE
      )
    }
    if (!is.null(ha_left)) {
      width_annotation <- grid::convertWidth(
        grid::unit(width_annotation, units) + ComplexHeatmap::width.HeatmapAnnotation(ha_left),
        units,
        valueOnly = TRUE
      )
    }
    if (!is.null(ha_right)) {
      width_annotation <- grid::convertWidth(
        grid::unit(width_annotation, units) + ComplexHeatmap::width.HeatmapAnnotation(ha_right),
        units,
        valueOnly = TRUE
      )
    }
  }
  dend_width <- name_width <- NULL
  dend_height <- name_height <- NULL
  if (inherits(ht_list, "HeatmapList")) {
    for (nm in names(ht_list@ht_list)) {
      ht <- ht_list@ht_list[[nm]]
      dend_width <- max(ht@row_dend_param$width, dend_width)
      dend_height <- max(ht@column_dend_param$height, dend_height)
      name_width <- max(ht@row_names_param$max_width, name_width)
      name_height <- max(ht@column_names_param$max_height, name_height)
    }
  } else if (inherits(ht_list, "Heatmap")) {
    ht <- ht_list
    dend_width <- max(ht@row_dend_param$width, dend_width)
    dend_height <- max(ht@column_dend_param$height, dend_height)
    name_width <- max(ht@row_names_param$max_width, name_width)
    name_height <- max(ht@column_names_param$max_height, name_height)
  } else {
    stop("ht_list is not a class of HeatmapList or Heatmap.")
  }

  lgd_width <- grid::convertWidth(
    grid::unit(
      unlist(lapply(legend_list, ComplexHeatmap::width.Legends)),
      grid::unitType(
        ComplexHeatmap::width.Legends(legend_list[[1]])
      )
    ),
    unitTo = units,
    valueOnly = TRUE
  )
  width_sum <- grid::convertWidth(
    grid::unit(width_sum, units) +
      grid::unit(width_annotation, units) +
      dend_width +
      name_width,
    units,
    valueOnly = TRUE
  ) +
    sum(lgd_width)
  height_sum <- max(
    grid::convertHeight(
      grid::unit(height_sum, units) +
        grid::unit(height_annotation, units) +
        dend_height +
        name_height,
      units,
      valueOnly = TRUE
    ),
    grid::convertHeight(
      grid::unit(0.95, "npc"),
      units,
      valueOnly = TRUE
    )
  )
  return(
    list(
      width_sum = width_sum,
      height_sum = height_sum
    )
  )
}

standardise <- function(data) {
  data[] <- t(apply(data, 1, scale))
  return(data)
}

mestimate <- function(data) {
  N <- nrow(data)
  D <- ncol(data)
  m.sj <- 1 +
    (1418 / N + 22.05) * D^(-2) +
    (12.33 / N + 0.243) *
      D^(-0.0406 * log(N) - 0.1134)
  return(m.sj)
}

adjustlayout <- function(
    graph,
    layout,
    width,
    height = 2,
    scale = 100,
    iter = 100) {
  w <- width / 2
  layout[, 1] <- layout[, 1] / diff(range(layout[, 1])) * scale
  layout[, 2] <- layout[, 2] / diff(range(layout[, 2])) * scale

  adjusted <- c()
  for (v in order(igraph::degree(graph), decreasing = TRUE)) {
    adjusted <- c(adjusted, v)
    neighbors <- as.numeric(
      igraph::neighbors(
        graph, igraph::V(graph)[v]
      )
    )
    neighbors <- setdiff(neighbors, adjusted)
    x <- layout[v, 1]
    y <- layout[v, 2]
    r <- w[v]
    for (neighbor in neighbors) {
      nx <- layout[neighbor, 1]
      ny <- layout[neighbor, 2]
      ndist <- sqrt((nx - x)^2 + (ny - y)^2)
      nr <- w[neighbor]
      expect <- r + nr
      if (ndist < expect) {
        dx <- (x - nx) * (expect - ndist) / ndist
        dy <- (y - ny) * (expect - ndist) / ndist
        layout[neighbor, 1] <- nx - dx
        layout[neighbor, 2] <- ny - dy
        adjusted <- c(adjusted, neighbor)
      }
    }
  }

  for (i in seq_len(iter)) {
    dist_matrix <- Matrix::as.matrix(
      proxyC::dist(layout)
    )
    nearest_neighbors <- apply(
      dist_matrix,
      2,
      function(x) which(x == min(x[x > 0])),
      simplify = FALSE
    )

    for (v in sample(seq_len(nrow(layout)))) {
      neighbors <- unique(nearest_neighbors[[v]])
      x <- layout[v, 1]
      y <- layout[v, 2]
      r <- w[v]
      for (neighbor in neighbors) {
        nx <- layout[neighbor, 1]
        ny <- layout[neighbor, 2]
        nr <- w[neighbor]
        if (abs(nx - x) < (r + nr) && abs(ny - y) < height) {
          dx <- r + nr - (nx - x)
          dy <- height - (ny - y)
          if (sample(c(1, 0), 1) == 1) {
            dx <- 0
          } else {
            dy <- 0
          }
          layout[neighbor, 1] <- nx - dx
          layout[neighbor, 2] <- ny - dy
        }
      }
    }
  }
  return(layout)
}
