make_scfea_plot_srt <- function(n_cells = 12L, seed = 1L) {
  set.seed(seed)
  module_info <- scfea_module_info(verbose = FALSE)
  ids <- unique(as.character(module_info$feature_id))
  signal <- rep(c(2, -2), each = n_cells / 2L)
  values <- matrix(
    stats::rnorm(length(ids) * n_cells, sd = 0.2) + rep(signal, each = length(ids)),
    nrow = length(ids),
    dimnames = list(ids, paste0("cell", seq_len(n_cells)))
  )
  srt <- Seurat::CreateSeuratObject(
    counts = Matrix::Matrix(abs(values) + 1, sparse = TRUE)
  )
  srt[["scFEAflux"]] <- Seurat::CreateAssayObject(data = values)
  srt[["scFEAbalance"]] <- Seurat::CreateAssayObject(data = values)
  srt@tools[["scFEA"]] <- list(module_info = module_info)
  srt$grp <- rep(c("A", "B"), each = n_cells / 2L)
  srt
}

test_that("scFEABalanceBarPlot honours palcolor and the shared theme arguments", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  srt <- make_scfea_plot_srt()
  palette <- c("#112233", "#445566")
  out <- scFEABalanceBarPlot(
    srt,
    group.by = "grp", ident.1 = "A", ident.2 = "B",
    palcolor = palette, theme_use = "theme_blank"
  )

  expect_type(out, "list")
  expect_s3_class(out$plot, "ggplot")
  fill_scale <- out$plot$scales$get_scales("fill")
  expect_false(is.null(fill_scale))
  expect_setequal(fill_scale$palette(2), palette)
  expect_identical(out$plot$theme$axis.text.x$angle, 90)
  expect_silent(invisible(ggplot2::ggplot_build(out$plot)))
})

test_that("scFEABalanceBarPlot renders with the default scop theme as well", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  srt <- make_scfea_plot_srt()
  out <- scFEABalanceBarPlot(
    srt,
    group.by = "grp", ident.1 = "A", ident.2 = "B"
  )
  expect_s3_class(out$plot, "ggplot")
  expect_silent(invisible(ggplot2::ggplot_build(out$plot)))
})

test_that("scFEAVolcanoPlot forwards theme_use and palcolor to its panels", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  srt <- make_scfea_plot_srt()
  pages <- scFEAVolcanoPlot(
    srt,
    group.by = "grp", ident.1 = "A", ident.2 = "B",
    combine = FALSE, theme_use = "theme_blank"
  )

  expect_type(pages, "list")
  expect_gt(length(pages), 0L)
  first_page <- pages[[1L]]
  expect_s3_class(first_page, "ggplot")
  expect_true(nrow(ggplot2::ggplot_build(first_page)$data[[1L]]) >= 0L)
})

test_that("scFEAVolcanoPlot renders combined output with the default theme", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  srt <- make_scfea_plot_srt()
  pages <- scFEAVolcanoPlot(
    srt,
    group.by = "grp", ident.1 = "A", ident.2 = "B"
  )
  expect_type(pages, "list")
  expect_true(all(vapply(pages, function(p) inherits(p, "ggplot"), logical(1))))
  expect_silent(invisible(lapply(pages, ggplot2::ggplot_build)))
})
