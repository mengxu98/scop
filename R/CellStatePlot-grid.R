cellstate_encode <- function(data, min.pct) {
  markers <- data$markers
  valid <- is.finite(markers$mean_z) & is.finite(markers$proportion)
  markers$show <- valid & markers$proportion >= min.pct
  domain <- if (any(markers$show)) range(markers$proportion[markers$show]) else c(0, 1)
  if (!any(markers$show) || diff(domain) == 0) domain <- c(min.pct, 1)
  if (diff(domain) == 0) domain <- c(0, 1)
  marker.radius <- function(p) {
    u <- pmax(0, pmin(1, (p - domain[1]) / diff(domain)))
    sqrt(0.0354 * (20 + 130 * u^1.8))
  }
  pathways <- data$pathways
  upper <- max(c(0, pathways$neg_log_p), na.rm = TRUE)
  if (upper == 0) upper <- 1
  pathway.radius <- function(p) {
    3.259 * sqrt(pmax(0, pmin(1, p / upper)))
  }
  markers$radius <- marker.radius(markers$proportion)
  pathways$radius <- pathway.radius(pathways$neg_log_p)
  list(
    markers = markers, pathways = pathways, marker.valid = valid,
    marker.radius = marker.radius, pathway.radius = pathway.radius,
    size.ticks = c(.2, .6, 1) * upper
  )
}

cellstate_palettes <- function(data, palcolor) {
  clusters <- data$clusters
  classes <- unique(clusters$class[order(clusters$cluster)])
  subclasses <- unique(clusters$subclass)
  base <- c("#9DC75D", "#BFBCBC", "#9F522C", "#35ADE2", "#BA8D8D", "#79A38F", "#CF8C7D")
  head <- c("#EDEEC4", "#E6E4E4", "#D8B794", "#9EDDF5", "#DAC6C6", "#C7D9D1", "#ECC1BA")
  body <- c("#EEF6E2", "#E6E4E4", "#F3EADF", "#E2F5FC", "#F6EEEE", "#EDF4F0", "#FAECE8")
  mix <- function(color, amount) {
    rgb <- grDevices::col2rgb(color) / 255
    grDevices::rgb(t(rgb * (1 - amount) + amount))
  }
  if (length(classes) > length(base)) {
    base <- grDevices::hcl.colors(length(classes), "Dark 3")
    head <- vapply(base, mix, character(1), amount = .5)
    body <- vapply(base, mix, character(1), amount = .87)
  }
  base <- stats::setNames(base[seq_along(classes)], classes)
  if (!is.null(palcolor)) {
    if (!is.list(palcolor) ||
        any(!names(palcolor) %in% c("class", "subclass", "program", "sample"))) {
      stop("palcolor must be a list with class, subclass, program or sample colors.")
    }
    for (kind in names(palcolor)) {
      colors <- palcolor[[kind]]
      if (is.null(names(colors)) || anyNA(colors) || anyDuplicated(names(colors))) {
        stop("palcolor entries must be uniquely named color vectors.")
      }
      grDevices::col2rgb(colors)
    }
  }
  override <- function(colors, custom) {
    if (!is.null(custom)) {
      hit <- intersect(names(custom), names(colors))
      colors[hit] <- custom[hit]
    }
    colors
  }
  custom <- if (is.null(palcolor)) list() else palcolor
  base <- override(base, custom$class)
  head <- stats::setNames(head[seq_along(classes)], classes)
  body <- stats::setNames(body[seq_along(classes)], classes)
  if (!is.null(custom$class)) {
    use <- intersect(names(custom$class), classes)
    head[use] <- vapply(base[use], mix, character(1), amount = .5)
    body[use] <- vapply(base[use], mix, character(1), amount = .87)
  }
  subcolors <- stats::setNames(rep(NA_character_, length(subclasses)), subclasses)
  for (class in classes) {
    ids <- unique(clusters$subclass[clusters$class == class])
    shades <- seq(0, .18, length.out = length(ids))
    for (j in seq_along(ids)) {
      if (is.na(subcolors[ids[j]])) subcolors[ids[j]] <- mix(base[class], shades[j])
    }
  }
  categories <- function(kind, initial, fallback) {
    ids <- unique(as.character(data[[paste0(kind, "s")]]$category))
    colors <- if (length(ids) <= length(initial)) initial[seq_along(ids)] else
      grDevices::hcl.colors(length(ids), fallback)
    override(stats::setNames(colors, ids), custom[[kind]])
  }
  list(
    header = head, body = body,
    subclass = override(subcolors, custom$subclass),
    program = categories("program", c("#2077B4", "#F47E20", "#2D9F48", "#D52928", "#9168AB", "#8B574C"), "Dark 3"),
    sample = categories("sample", c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5"), "Set 3")
  )
}

cellstate_grob <- function(data, tree, features, feature.labels, pathways,
                           pathway.labels, encoded, palcolor, family, title,
                           panel_label) {
  if (length(family) != 1L || is.na(family)) stop("family must be one font family.")
  if (length(title) > 1L || length(panel_label) > 1L ||
      anyNA(title) || anyNA(panel_label)) stop("title and panel_label must be single strings or NULL.")
  measure <- function(s, size = 6, face = "plain") {
    if (!length(s)) return(numeric())
    if (requireNamespace("systemfonts", quietly = TRUE)) {
      return(vapply(s, function(value) {
        systemfonts::string_width(
          value, family = family, size = size, res = 7200,
          weight = if (face == "bold") "bold" else "normal", italic = face == "italic"
        ) / 100
      }, numeric(1)))
    }
    nchar(s, type = "width") * size * .55
  }
  wrap <- function(s, width) {
    if (grepl("\n", s, fixed = TRUE)) return(s)
    words <- strsplit(s, " ", fixed = TRUE)[[1]]
    lines <- words[1]
    if (length(words) > 1L) {
      for (word in words[-1]) {
        j <- length(lines)
        if (measure(paste(lines[j], word)) <= width) {
          lines[j] <- paste(lines[j], word)
        } else {
          lines <- c(lines, word)
        }
      }
    }
    paste(lines, collapse = "\n")
  }
  palette <- cellstate_palettes(data, palcolor)
  n <- nrow(data$clusters)
  xs <- 28.2732 + (seq_len(n) - 1) * 7.07006
  left <- xs[1] - 5.7712
  right <- tail(xs, 1) + 6.3442
  label.x <- right + 6.0072
  runs <- rle(data$clusters$class)
  ends <- cumsum(runs$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)
  edges <- c(left, (xs[head(ends, -1L)] + xs[starts[-1L]]) / 2 - .35, right)
  headers <- mapply(wrap, gsub("_", " ", runs$values), diff(edges) - 2, USE.NAMES = FALSE)
  rows <- lengths(strsplit(headers, "\n", fixed = TRUE))
  header.starts <- integer(length(headers))
  occupied <- numeric()
  for (i in seq_along(headers)) {
    width <- max(measure(strsplit(headers[i], "\n", fixed = TRUE)[[1]]))
    center <- mean(edges[i + 0:1])
    row <- 0L
    repeat {
      use <- row + seq_len(rows[i])
      free <- vapply(use, function(j) {
        j > length(occupied) || occupied[j] + 2 <= center - width / 2
      }, logical(1))
      if (all(free)) break
      row <- row + 1L
    }
    header.starts[i] <- row
    length(occupied) <- max(length(occupied), max(use))
    occupied[is.na(occupied)] <- -Inf
    occupied[use] <- center + width / 2
  }
  title.lines <- if (is.null(title) || !nzchar(title)) character() else
    strsplit(wrap(title, max(110, right - left)), "\n", fixed = TRUE)[[1]]
  extra <- max(0, length(occupied) - 2L) * 7.002 +
    max(0, length(title.lines) - 1L) * 7.002
  display <- function(s) sub("^t cell$", "T cell", gsub("_", " ", s, fixed = TRUE))
  cluster.labels <- display(data$clusters$cluster)
  program.y <- 119.13574 + extra + max(0, max(measure(cluster.labels, 5)) - 42.6)
  sub.top <- 53.3329 + extra
  capacity <- max(1L, floor((program.y - 14 - sub.top) / 6.474) + 1L)
  sub.rows <- min(length(palette$subclass), capacity)
  sub.columns <- ceiling(length(palette$subclass) / sub.rows)
  sub.widths <- vapply(seq_len(sub.columns), function(j) {
    ids <- which((seq_along(palette$subclass) - 1L) %/% sub.rows == j - 1L)
    max(measure(display(names(palette$subclass)[ids]))) + 11.1
  }, numeric(1))
  sub.offsets <- c(0, head(cumsum(sub.widths + 5), -1L))
  program.y <- max(program.y, sub.top + (sub.rows - 1L) * 6.474 + 14)
  program.key.y <- program.y + 3.686
  sample.y <- max(program.y + 36.109,
                  program.key.y + (ceiling(length(palette$program) / 2) - 1L) * 6.4082 + 14)
  sample.key.y <- sample.y - .8684
  sample.end <- sample.key.y + (ceiling(length(palette$sample) / 2) - 1L) * 6.58332
  missing <- any(!encoded$marker.valid) || any(!is.finite(data$pathways$neg_log_p))
  gene.start <- max(sample.y + 40.69437, sample.end + if (missing) 15.2 else 8)
  gene.y <- gene.start + (seq_along(features) - 1L) * 5.418856
  path.lines <- lengths(strsplit(pathway.labels, "\n", fixed = TRUE))
  path.step <- max(9.483429, (max(path.lines) - 1L) * 5.25 + 4)
  path.y <- tail(gene.y, 1) + 10.67364 + (seq_along(pathways) - 1L) * path.step
  bottom <- tail(path.y, 1) + max(5.244, max(path.lines) * 2.625)
  legend.y <- bottom + 6.8023
  category.gap <- function(colors, minimum) {
    max(minimum, max(measure(head(names(colors), ceiling(length(colors) / 2)))) + 15)
  }
  program.gap <- category.gap(palette$program, 32.25)
  sample.gap <- category.gap(palette$sample, 25.885)
  category.width <- function(colors, gap) {
    nr <- ceiling(length(colors) / 2)
    max(measure(names(colors)) + 11.1 + ((seq_along(colors) - 1L) %/% nr) * gap)
  }
  key.end <- function(label, values, x) x + measure(label) + 6 + 44 + 4 + measure(tail(values, 1))
  size.labels <- format(round(encoded$size.ticks, 1), trim = TRUE)
  key.left.end <- max(key.end("Expression proportion", c(".40", ".60", ".80"), 6.8),
                      key.end("z-normalized expression", c("0", "1", "2"), 1))
  key.right.x <- max(146.85, key.left.end + 8)
  W <- max(244, label.x + 5 +
    max(76.484, max(measure(unlist(strsplit(c(feature.labels, pathway.labels), "\n", fixed = TRUE)), 5)) + 3,
        tail(sub.offsets + sub.widths, 1),
        category.width(palette$program, program.gap),
        category.width(palette$sample, sample.gap)),
    key.end("-ln(P)", size.labels, key.right.x) + 4,
    key.end("GSEA NES", c("-2.5", "0", "2.5"), key.right.x) + 4)
  H <- legend.y + 10.9507 + 5.177
  color.file <- system.file("extdata", "cell-state-colors.csv", package = "scop")
  colors <- utils::read.csv(color.file, stringsAsFactors = FALSE)
  layout <- list(
    W = W, H = H, xs = xs, left = left, right = right, label.x = label.x,
    runs = runs, edges = edges, headers = headers, header.starts = header.starts,
    occupied = occupied, extra = extra, cluster.labels = cluster.labels,
    title.lines = title.lines, panel_label = panel_label, family = family,
    program.y = program.y, program.key.y = program.key.y,
    sample.y = sample.y, sample.key.y = sample.key.y, sample.end = sample.end,
    sub.top = sub.top, sub.rows = sub.rows, sub.offsets = sub.offsets,
    program.gap = program.gap, sample.gap = sample.gap, missing = missing,
    gene.y = gene.y, path.y = path.y, bottom = bottom, legend.y = legend.y,
    key.right.x = key.right.x, measure = measure, display = display
  )
  plot <- grid::gTree(
    data = data, tree = tree, features = features, feature.labels = feature.labels,
    pathways = pathways, pathway.labels = pathway.labels, encoded = encoded,
    palette = palette, colors = colors, layout = layout,
    name = "cell-state-plot", cl = "scop_cell_state_plot"
  )
  attr(plot, "width_mm") <- W / 72 * 25.4
  attr(plot, "height_mm") <- H / 72 * 25.4
  plot
}

#' @export
print.scop_cell_state_plot <- function(x, ...) {
  grid::grid.newpage()
  grid::grid.draw(x)
  invisible(x)
}

#' @export
#' @method makeContent scop_cell_state_plot
#' @importFrom grid makeContent
makeContent.scop_cell_state_plot <- function(x) {
  a <- x$layout
  data <- x$data
  palette <- x$palette
  scale <- min(
    grid::convertWidth(grid::unit(1, "npc"), "inches", valueOnly = TRUE) * 72 / a$W,
    grid::convertHeight(grid::unit(1, "npc"), "inches", valueOnly = TRUE) * 72 / a$H
  )
  children <- list()
  add <- function(grob) children[[length(children) + 1L]] <<- grob
  native <- function(value) grid::unit(value, "native")
  rect <- function(xx, yy, w, h, fill, col = NA, lw = .25) {
    add(grid::rectGrob(native(xx), native(yy),
                       grid::unit(w / a$W, "npc"), grid::unit(h / a$H, "npc"),
                       just = c("left", "top"),
                       gp = grid::gpar(fill = fill, col = col, lwd = lw * scale / .75)))
  }
  line <- function(xx, yy, col = "#D3D3D3", lw = .259) {
    add(grid::linesGrob(native(xx), native(yy),
      gp = grid::gpar(col = col, lwd = lw * scale / .75, lineend = "square")))
  }
  text <- function(s, xx, yy, size = 6, face = "plain", just = "left",
                   rot = 0, col = "#231F20") {
    if (!length(s)) return(invisible(NULL))
    add(grid::textGrob(s, native(xx), native(yy), just = c(just, "bottom"), rot = rot,
      gp = grid::gpar(fontfamily = a$family, fontsize = size * scale,
                      fontface = face, col = col)))
  }
  dots <- function(xx, yy, radius, fill, lw = .25) {
    if (!length(xx)) return(invisible(NULL))
    yy <- rep(yy, length.out = length(xx))
    fill <- rep(fill, length.out = length(xx))
    keep <- is.finite(radius) & radius > 0
    if (!any(keep)) return(invisible(NULL))
    add(grid::circleGrob(native(xx[keep]), native(yy[keep]), native(radius[keep]),
      gp = grid::gpar(fill = fill[keep], col = "black", lwd = lw * scale / .75)))
  }
  cross <- function(xx, yy) {
    for (i in seq_along(xx)) {
      line(xx[i] + c(-.7, .7), yy[i] + c(-.7, .7), "#777777", .35)
      line(xx[i] + c(-.7, .7), yy[i] + c(.7, -.7), "#777777", .35)
    }
  }
  continuous <- function(values, kind) {
    if (!length(values)) return(character())
    p <- x$colors[x$colors$scale == kind, ]
    p <- p[order(p$value), ]
    rgb <- grDevices::col2rgb(p$color)
    channels <- matrix(vapply(seq_len(3), function(k) {
      stats::approx(p$value, rgb[k, ], xout = values, rule = 2)$y
    }, numeric(length(values))), nrow = length(values), ncol = 3)
    grDevices::rgb(channels[, 1], channels[, 2], channels[, 3], maxColorValue = 255)
  }
  ids <- data$clusters$cluster
  n <- length(ids)
  for (i in seq_along(a$runs$values)) {
    class <- a$runs$values[i]
    rect(a$edges[i], 64.229 + a$extra, diff(a$edges)[i],
         a$bottom - 64.229 - a$extra, palette$body[class])
    rect(a$edges[i], 27.918 + a$extra, diff(a$edges)[i], 35.872, palette$header[class])
    words <- strsplit(a$headers[i], "\n", fixed = TRUE)[[1]]
    title.extra <- max(0, length(a$title.lines) - 1L) * 7.002
    yy <- if (length(a$occupied) == 1L) 21.4217 + title.extra else
      17.9207 + title.extra + (a$header.starts[i] + seq_along(words) - 1L) * 7.002
    text(words, mean(a$edges[i + 0:1]), yy, just = "centre", col = "black")
  }
  text(a$panel_label, 1.7457, 8.1712, 8, "bold")
  text(a$title.lines, mean(c(a$left, a$right)), 8.1668 +
         (seq_along(a$title.lines) - 1L) * 7.002, just = "centre")
  leaf.x <- stats::setNames(a$xs, ids)
  tree.x <- tree.y <- numeric(n - 1L)
  height <- max(x$tree$height)
  for (i in seq_len(n - 1L)) {
    xy <- lapply(x$tree$merge[i, ], function(k) {
      if (k < 0) c(leaf.x[x$tree$labels[-k]], 62.271545 + a$extra) else
        c(tree.x[k], tree.y[k])
    })
    yy <- 62.271545 + a$extra -
      if (height > 0) 32.612 * x$tree$height[i] / height else 0
    line(c(xy[[1]][1], xy[[1]][1], xy[[2]][1], xy[[2]][1]),
         c(xy[[1]][2], yy, yy, xy[[2]][2]), "#555555", 1.242)
    tree.x[i] <- mean(c(xy[[1]][1], xy[[2]][1]))
    tree.y[i] <- yy
  }
  for (j in seq_len(n)) {
    rect(a$xs[j] - 2.2436, 65.7556 + a$extra, 3.536, 6.4534,
         palette$subclass[data$clusters$subclass[j]])
    text(a$cluster.labels[j], a$xs[j] + 1.063, 74.7436 + a$extra +
           a$measure(a$cluster.labels[j], 5), 5, rot = 90, col = "#58595B")
  }
  text("Classes", 16.7994, 56.4034 + a$extra, rot = 90)
  text("Subclasses", 6.7194, 84.5314 + a$extra, rot = 90)
  text("Clusters", 16.7994, 103.3474 + a$extra, rot = 90)
  line(c(12.55, 20.55), rep(69.1 + a$extra, 2), "#231F20", .5)
  stack <- function(table, colors, y) {
    for (j in seq_len(n)) {
      cursor <- y + 32.23
      for (category in names(colors)) {
        value <- table$fraction[table$cluster == ids[j] & table$category == category]
        if (!length(value) || value == 0) next
        h <- 32.23 * value
        cursor <- cursor - h
        rect(a$xs[j] - 1.8405, cursor, 3.532, h, colors[category], "black")
      }
    }
  }
  stack(data$programs, palette$program, a$program.y)
  stack(data$samples, palette$sample, a$sample.y)
  for (xx in a$xs) {
    line(c(xx, xx), c(a$gene.y[1] - 5.6911, tail(a$gene.y, 1) + 5.6899))
    line(c(xx, xx), c(tail(a$gene.y, 1) + 5.9113, a$bottom - 2.8724))
  }
  for (yy in c(a$gene.y, a$path.y)) line(c(a$left - 3.7733, a$right + 3.2515), rep(yy, 2))
  markers <- data$markers
  mx <- a$xs[match(markers$cluster, ids)]
  my <- a$gene.y[match(markers$feature, x$features)]
  keep <- markers$show
  dots(mx[keep], my[keep], markers$radius[keep], continuous(markers$mean_z[keep], "expression"))
  cross(mx[!x$encoded$marker.valid], my[!x$encoded$marker.valid])
  gsea <- data$pathways
  gx <- a$xs[match(gsea$cluster, ids)]
  gy <- a$path.y[match(gsea$feature, x$pathways)]
  keep <- is.finite(gsea$neg_log_p) & is.finite(gsea$nes)
  dots(gx[keep], gy[keep], gsea$radius[keep], continuous(gsea$nes[keep], "enrichment"))
  cross(gx[!keep], gy[!keep])
  text(x$feature.labels, a$label.x, a$gene.y + 1.5384, 5, "italic", col = "black")
  for (i in seq_along(x$pathway.labels)) {
    words <- strsplit(x$pathway.labels[i], "\n", fixed = TRUE)[[1]]
    yy <- a$path.y[i] + 1.54 - (length(words) - 1L) * 2.625 +
      (seq_along(words) - 1L) * 5.25
    text(words, a$label.x, yy, 5, col = "black")
  }
  text(c("Marker gene", "expression"), c(9.8267, 16.8287),
       mean(range(a$gene.y)) + c(7.735, 5.509), rot = 90)
  text(c("Pathway", "enrichment"), c(9.8267, 16.8287),
       mean(range(a$path.y)) + c(7.045, 10.501), rot = 90)
  text(c("Cell population", "annotations"), a$label.x, 38.0432 + a$extra + c(0, 7.002))
  for (i in seq_along(palette$subclass)) {
    xx <- a$label.x + a$sub.offsets[(i - 1L) %/% a$sub.rows + 1L]
    yy <- a$sub.top + ((i - 1L) %% a$sub.rows) * 6.474
    rect(xx, yy - 3.3619, 8.651, 3.028, palette$subclass[i])
    text(a$display(names(palette$subclass)[i]), xx + 11.0013, yy)
  }
  category.key <- function(colors, y, gap) {
    nr <- ceiling(length(colors) / 2)
    for (i in seq_along(colors)) {
      xx <- a$label.x + ((i - 1L) %/% nr) * gap
      yy <- y + ((i - 1L) %% nr) * 6.5
      rect(xx, yy - 3.3437, 8.65, 3.03, colors[i])
      text(names(colors)[i], xx + 11.1, yy)
    }
  }
  text("Program assignment", a$label.x, a$program.y - 4.2046)
  category.key(palette$program, a$program.key.y, a$program.gap)
  text("Donors / samples", a$label.x, a$sample.y - 8.8364)
  category.key(palette$sample, a$sample.key.y, a$sample.gap)
  if (a$missing) {
    cross(a$label.x + 1, a$sample.end + 5.6)
    text("Not evaluable", a$label.x + 4, a$sample.end + 7.2, 5, col = "#666666")
  }
  dot.key <- function(label, labels, radii, fill, xx, yy) {
    text(label, xx, yy + 1.84, col = "black")
    dx <- xx + a$measure(label) + 6 + c(0, 22, 44)
    dots(dx, rep(yy, 3), radii, fill)
    text(labels, dx + 4, yy + 1.84, col = "black")
  }
  dot.key("Expression proportion", c("0.40", "0.60", "0.80"),
          x$encoded$marker.radius(c(.4, .6, .8)), "#D3D3D3", 6.8, a$legend.y)
  dot.key("z-normalized expression", c("0", "1", "2"), rep(1.505, 3),
          continuous(c(0, 1, 2), "expression"), 1, a$legend.y + 10.9507)
  dot.key("-ln(P)", format(round(x$encoded$size.ticks, 1), trim = TRUE),
          x$encoded$pathway.radius(x$encoded$size.ticks), "#333333",
          a$key.right.x, a$legend.y)
  dot.key("GSEA NES", c("-2.5", "0", "2.5"), rep(1.129, 3),
          continuous(c(-2.5, 0, 2.5), "enrichment"),
          a$key.right.x, a$legend.y + 10.9507)
  viewport <- grid::viewport(
    width = grid::unit(a$W * scale, "bigpts"), height = grid::unit(a$H * scale, "bigpts"),
    xscale = c(0, a$W), yscale = c(a$H, 0)
  )
  child <- grid::gTree(children = do.call(grid::gList, children), vp = viewport)
  grid::setChildren(x, grid::gList(child))
}
