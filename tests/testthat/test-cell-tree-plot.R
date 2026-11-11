test_that("annotation paths determine leaf counts and percentages", {
  skip_if_not_installed("ggforce")
  hierarchy <- c("major", "minor")
  annotations <- data.frame(
    major = c("B", "A", "A", "B", "A", "A/B"),
    minor = c("Shared", "Shared", "Shared", "Only", "B/C", "C")
  )
  result <- scop::CellTreePlot(annotations, hierarchy, return_data = TRUE)

  expect_s3_class(result$plot, "ggplot")
  expect_null(result$tracks)
  expect_equal(result$leaves$major, c("A", "A", "A/B", "B", "B"))
  expect_equal(result$leaves$minor, c("B/C", "Shared", "C", "Only", "Shared"))
  expect_equal(result$leaves$.count, c(1, 2, 1, 1, 1))
  expect_equal(result$leaves$.proportion, c(1, 2, 1, 1, 1) / 6)
  expect_no_warning(ggplot2::ggplotGrob(result$plot))

  summary <- data.frame(
    major = c("B", "A"), minor = c("Shared", "Shared"),
    cells = c(0, 30)
  )
  weighted <- CellTreePlot(summary, hierarchy,
    count = "cells",
    return_data = TRUE
  )
  expect_equal(weighted$leaves$.count, c(30, 0))
  expect_equal(weighted$leaves$.proportion, c(1, 0))
})

test_that("rings use complete annotation paths and preserve missing values", {
  skip_if_not_installed("ggforce")
  hierarchy <- c("major", "minor")
  annotations <- data.frame(
    major = c("B", "A", "A"),
    minor = c("Shared", "Other", "Shared")
  )
  values <- data.frame(
    major = c("B", "A"), minor = c("Shared", "Shared"),
    RNA = c(0, 0.3), ATAC = c(-0.5, NA_real_)
  )
  result <- CellTreePlot(annotations, hierarchy,
    tracks = values,
    return_data = TRUE
  )
  expect_equal(result$tracks$major, c("A", "A", "B"))
  expect_equal(result$tracks$minor, c("Other", "Shared", "Shared"))
  expect_equal(result$tracks$RNA, c(NA_real_, 0.3, 0))
  expect_equal(result$tracks$ATAC, c(NA_real_, NA_real_, -0.5))
  fill_scale <- ggplot2::ggplot_build(result$plot)$plot$scales$get_scales("fill")
  expect_length(unique(fill_scale$map(c(NA_real_, 0))), 2L)
  breaks <- fill_scale$get_breaks()
  expect_length(breaks, 3L)
  expect_lt(breaks[[1L]], 0)
  expect_equal(breaks, c(-breaks[[3L]], 0, breaks[[3L]]))
  expect_no_warning(ggplot2::ggplotGrob(result$plot))

  shuffled <- CellTreePlot(annotations[3:1, ], hierarchy,
    tracks = values[2:1, ], return_data = TRUE
  )
  expect_equal(shuffled$leaves, result$leaves)
  expect_equal(shuffled$tracks, result$tracks)
  expect_error(CellTreePlot(annotations, hierarchy,
    tracks = values[c(1, 1), ]
  ))
  expect_error(CellTreePlot(annotations, hierarchy,
    tracks = transform(values, minor = "Unknown")
  ))
  expect_error(CellTreePlot(annotations, hierarchy,
    tracks = transform(values, RNA = Inf)
  ))
})

test_that("named custom colors match annotation groups", {
  skip_if_not_installed("ggforce")
  annotations <- data.frame(group = c("B", "A"))
  colors <- c(B = "#2040AA", A = "#AA4020")
  plot <- scop::CellTreePlot(annotations, "group", palcolor = colors)
  scale <- ggplot2::ggplot_build(plot)$plot$scales$get_scales("colour")
  expect_identical(scale$map(c("A", "B")), unname(colors[c("A", "B")]))
})

test_that("ring colors remain distinct across extreme finite values", {
  skip_if_not_installed("ggforce")
  annotations <- data.frame(group = c("A", "B", "C", "D"))
  values <- c(-1e308, 0, 1e308, NA_real_)
  plot <- scop::CellTreePlot(
    annotations, "group",
    tracks = transform(annotations, score = values)
  )
  built <- ggplot2::ggplot_build(plot)
  ring <- Filter(function(layer) "fill" %in% names(layer), built$data)[[1L]]
  tiles <- unique(ring[c("group", "fill")])
  expect_identical(
    tiles$fill[order(tiles$group)], c("#2C6B9E", "#F8F7F4", "#BD4E42", "#A9ADB4")
  )
  key <- ggplot2::get_guide_data(plot, "fill")
  expect_true(all(is.finite(key$.value)) && all(diff(key$.value) > 0))
  expect_equal(sign(as.numeric(key$.label)), c(-1, 0, 1))
  expect_no_warning(ggplot2::ggplotGrob(plot))
})

test_that("each leaf has a broad-group tile inside the numeric rings", {
  skip_if_not_installed("ggforce")
  annotations <- data.frame(major = c("A", "A", "B"), minor = c("x", "y", "x"))
  radius <- function(layer) sqrt(layer$x^2 + layer$y^2)

  for (tracks in list(NULL, transform(annotations, score = 1:3))) {
    plot <- CellTreePlot(annotations, c("major", "minor"), tracks = tracks)
    built <- ggplot2::ggplot_build(plot)
    arcs <- which(vapply(plot$layers, function(layer) {
      inherits(layer$geom, "GeomArcBar")
    }, logical(1)))
    band <- built$data[[arcs[[2L]]]]
    fill <- grep("^fill", names(band), value = TRUE)
    tiles <- unique(band[c("group", fill)])
    expect_equal(nrow(tiles), nrow(annotations))
    expect_identical(duplicated(tiles[[fill]]), c(FALSE, TRUE, FALSE))
    if (!is.null(tracks)) {
      expect_lt(max(radius(band)), min(radius(built$data[[arcs[[3L]]]])))
    }
  }
})

test_that("one-leaf and deeper hierarchies render", {
  skip_if_not_installed("ggforce")
  one <- data.frame(label = "NA")
  one_ring <- data.frame(label = "NA", score = 0)
  p <- CellTreePlot(one, "label", tracks = one_ring, return_data = TRUE)
  expect_equal(p$leaves$.proportion, 1)
  expect_no_warning(ggplot2::ggplotGrob(p$plot))

  nested <- data.frame(
    broad = c("A", "A", "B"),
    fine = c("Family", "Family", "Family"),
    leaf = c("x", "y", "x")
  )
  result <- CellTreePlot(nested, c("broad", "fine", "leaf"),
    return_data = TRUE
  )
  expect_equal(result$leaves$.proportion, rep(1 / 3, 3))
  expect_null(result$plot$labels$title)
  expect_no_warning(ggplot2::ggplotGrob(result$plot))
  text_labels <- function(plot) {
    unlist(lapply(ggplot2::ggplot_build(plot)$data, function(layer) layer$label))
  }
  expect_false("Family" %in% text_labels(result$plot))
  labeled <- CellTreePlot(
    nested, c("broad", "fine", "leaf"),
    title = "My tree", show_internal_labels = TRUE
  )
  expect_identical(labeled$labels$title, "My tree")
  expect_true("Family" %in% text_labels(labeled))
  expect_error(CellTreePlot(
    transform(nested, leaf = NA_character_),
    c("broad", "fine", "leaf")
  ))
})
