spatial_pie_dependency_object <- function() {
  counts <- Matrix::Matrix(
    matrix(1, nrow = 2, ncol = 3,
      dimnames = list(c("Gene1", "Gene2"), paste0("Spot", 1:3))),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts)
  object$x <- c(1, 2, 3)
  object$y <- c(1, 2, 1)
  object
}

test_that("empty pies render without looking up the optional backend", {
  object <- spatial_pie_dependency_object()
  dependency_calls <- character()
  testthat::local_mocked_bindings(
    check_r = function(packages, ...) {
      dependency_calls <<- c(dependency_calls, packages)
      stop("Unexpected pie dependency lookup")
    },
    .package = "scop"
  )
  values <- list(
    zero = matrix(0, 3, 2),
    missing = matrix(NA_real_, 3, 2),
    incomplete = cbind(c(1, NA, 0), c(NA, 1, 0))
  )
  for (value in values) {
    dimnames(value) <- list(colnames(object), c("A", "B"))
    plot <- scop::SpatialSpotPlot(
      object, values = value, plot_type = "pie", overlay_image = FALSE,
      legend.title = "Empty proportions", theme_use = NULL
    )
    expect_s3_class(plot, "ggplot")
    expect_identical(plot$data$label, "No spots with positive pie values")
    expect_identical(plot$labels$title, "Empty proportions")
    expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
    expect_identical(dependency_calls, character())
  }
})

test_that("pie filtering respects selected cells before requiring a backend", {
  object <- spatial_pie_dependency_object()
  dependency_calls <- character()
  testthat::local_mocked_bindings(
    check_r = function(packages, ...) {
      dependency_calls <<- c(dependency_calls, packages)
      stop("Unexpected pie dependency lookup")
    },
    .package = "scop"
  )
  values <- matrix(c(1, 0, 0, 0, 0, 0), 3, 2,
    dimnames = list(colnames(object), c("A", "B")))
  plot <- scop::SpatialSpotPlot(
    object, values = values, cells = c("Spot2", "Spot3"),
    plot_type = "pie", overlay_image = FALSE, theme_use = NULL
  )
  expect_identical(plot$data$label, "No spots with positive pie values")
  expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
  expect_identical(dependency_calls, character())
})

test_that("invalid pie values fail validation before a dependency lookup", {
  object <- spatial_pie_dependency_object()
  dependency_calls <- character()
  testthat::local_mocked_bindings(
    check_r = function(packages, ...) {
      dependency_calls <<- c(dependency_calls, packages)
      stop("Unexpected pie dependency lookup")
    },
    .package = "scop"
  )
  for (invalid in c(-1, Inf, -Inf)) {
    values <- matrix(NA_real_, 3, 2,
      dimnames = list(colnames(object), c("A", "B")))
    values[1, 1] <- invalid
    expect_error(
      scop::SpatialSpotPlot(
        object, values = values, plot_type = "pie", overlay_image = FALSE
      ),
      "Pie values must be non-negative and cannot be infinite"
    )
    expect_identical(dependency_calls, character())
  }
})

test_that("nonempty pies still require the scatterpie backend", {
  object <- spatial_pie_dependency_object()
  dependency_calls <- list()
  testthat::local_mocked_bindings(
    check_r = function(packages, verbose, ...) {
      dependency_calls[[length(dependency_calls) + 1L]] <<-
        list(packages = packages, verbose = verbose)
      stop("Pie backend is unavailable")
    },
    .package = "scop"
  )
  values <- matrix(0, 3, 2,
    dimnames = list(colnames(object), c("A", "B")))
  values[2, ] <- c(2, 3)
  expect_error(
    scop::SpatialSpotPlot(
      object, values = values, plot_type = "pie", overlay_image = FALSE
    ),
    "Pie backend is unavailable"
  )
  expect_identical(dependency_calls, list(list(packages = "scatterpie", verbose = FALSE)))
})
