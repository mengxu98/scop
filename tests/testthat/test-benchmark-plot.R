make_scop_benchmark_plot_result <- function(status = c("success", "success", "failed")) {
  summary <- data.frame(
    method = c("BayesSpace", "BANKSY", "A very long unavailable method name"),
    ARI = c(0.82, 0.63, NA),
    NMI = c(0.76, 0.59, NA),
    purity = c(0.88, 0.71, NA),
    runtime_s = c(1.2, 24, NA),
    baseline_memory_mb = c(90, 120, NA),
    peak_memory_mb = c(180, 2400, NA),
    memory_delta_mb = c(90, 2280, NA),
    n_evaluated = c(100L, 96L, 0L),
    n_clusters = c(4L, 5L, NA_integer_),
    status = status,
    error = c("", "", "optional backend is unavailable"),
    stringsAsFactors = FALSE
  )
  metric_names <- c("ARI", "NMI", "purity", "runtime_s", "peak_memory_mb")
  metrics <- do.call(rbind, lapply(seq_len(nrow(summary)), function(i) {
    data.frame(
      method = summary$method[[i]],
      metric = metric_names,
      value = unlist(summary[i, metric_names], use.names = FALSE),
      workflow = "SpatialDomain",
      direction = c("higher", "higher", "higher", "lower", "lower"),
      status = summary$status[[i]],
      stringsAsFactors = FALSE
    )
  }))
  structure(
    list(
      summary = summary,
      metrics = metrics,
      predictions = data.frame(),
      runs = data.frame(),
      objects = list(),
      parameters = list(metrics = c("ARI", "NMI"))
    ),
    class = c("scop_benchmark", "list")
  )
}

test_that("BenchmarkPlot selects publication benchmark views", {
  result <- make_scop_benchmark_plot_result()
  expect_s3_class(BenchmarkPlot(data = result, plot_type = "quality"), "ggplot")
  expect_s3_class(BenchmarkPlot(data = result, plot_type = "efficiency"), "ggplot")
  expect_s3_class(BenchmarkPlot(data = result, plot_type = "heatmap"), "ggplot")
  expect_s3_class(BenchmarkPlot(data = result), "patchwork")
  expect_s3_class(plot(result), "patchwork")
})

test_that("quality ordering and automatic log scales are explicit", {
  result <- make_scop_benchmark_plot_result()
  expect_identical(
    benchmark_plot_method_levels(result, c("ARI", "NMI"), "quality")[[1]],
    "BayesSpace"
  )
  efficiency <- BenchmarkPlot(data = result, plot_type = "efficiency")
  expect_true(inherits(efficiency$scales$get_scales("x")$trans, "transform"))
  expect_identical(efficiency$labels$x, "Runtime (s, log10)")
  expect_match(efficiency$labels$y, "log10")
  labels <- unlist(lapply(ggplot2::ggplot_build(efficiency)$data, `[[`, "label"))
  expect_true(any(grepl("s \\| .*(MB|GB)", labels)))
})

test_that("quality view uses score tracks and explicit direction", {
  quality <- BenchmarkPlot(
    data = make_scop_benchmark_plot_result(),
    plot_type = "quality"
  )
  expect_identical(quality$labels$x, "Agreement score (higher is better)")
  expect_true(any(vapply(quality$layers, function(layer) {
    inherits(layer$geom, "GeomSegment")
  }, logical(1))))
})

test_that("legacy benchmark tiers are visible without changing method identity", {
  result <- make_scop_benchmark_plot_result(rep("success", 3))
  result$summary$tier <- c("stable", "stable", "legacy")
  labels <- benchmark_plot_method_labels(result, result$summary$method)
  expect_identical(labels[[1]], "BayesSpace")
  expect_match(labels[[3]], "\\[legacy\\]")
  expect_no_error(ggplot2::ggplot_build(BenchmarkPlot(data = result, plot_type = "quality")))
})

test_that("heatmap normalization respects metric direction", {
  expect_equal(benchmark_normalize_metric(c(1, 3), "higher"), c(0, 1))
  expect_equal(benchmark_normalize_metric(c(1, 3), "lower"), c(1, 0))
  expect_equal(benchmark_normalize_metric(c(2, 2), "higher"), c(0.5, 0.5))

  heatmap <- BenchmarkPlot(
    data = make_scop_benchmark_plot_result(),
    plot_type = "heatmap"
  )
  expect_gte(length(heatmap$layers), 2L)
  expect_match(heatmap$labels$subtitle, "labels show raw values")
  expect_identical(levels(heatmap$data$metric)[4:5], c("Runtime", "Peak memory"))
})

test_that("status strip appears only for incomplete runs", {
  incomplete <- make_scop_benchmark_plot_result()
  complete <- make_scop_benchmark_plot_result(rep("success", 3))
  incomplete_plot <- BenchmarkPlot(data = incomplete)
  complete_plot <- BenchmarkPlot(data = complete)
  expect_identical(incomplete_plot[[2]]$labels$title, "Incomplete runs")
  expect_identical(complete_plot[[2]]$labels$title, "Computational efficiency")
})

test_that("legacy benchmark bar mode remains available", {
  data <- data.frame(
    method = c("A", "B"), metric = c("ARI", "ARI"), value = c(0.8, 0.6)
  )
  expect_s3_class(BenchmarkPlot(data = data, plot_type = "bar"), "ggplot")
})

test_that("all-unavailable benchmark results render truthful empty panels", {
  result <- make_scop_benchmark_plot_result(rep("unavailable", 3))
  result$summary[, c(
    "ARI", "NMI", "purity", "runtime_s", "baseline_memory_mb",
    "peak_memory_mb", "memory_delta_mb"
  )] <- NA_real_
  result$summary$n_evaluated <- 0L
  result$summary$n_clusters <- NA_integer_
  result$summary$error <- "optional backend is unavailable"

  quality <- BenchmarkPlot(data = result, plot_type = "quality")
  efficiency <- BenchmarkPlot(data = result, plot_type = "efficiency")
  heatmap <- BenchmarkPlot(data = result, plot_type = "heatmap")
  overview <- BenchmarkPlot(data = result)

  expect_s3_class(quality, "ggplot")
  expect_s3_class(efficiency, "ggplot")
  expect_s3_class(heatmap, "ggplot")
  expect_s3_class(overview, "patchwork")
  expect_no_error(ggplot2::ggplot_build(quality))
  expect_no_error(ggplot2::ggplot_build(efficiency))
  expect_no_error(ggplot2::ggplot_build(heatmap))
  expect_no_error(patchwork::patchworkGrob(overview))
})

test_that("zero or missing resource measurements do not fabricate efficiency points", {
  result <- make_scop_benchmark_plot_result(rep("success", 3))
  result$summary$runtime_s <- c(0, NA, Inf)
  result$summary$peak_memory_mb <- c(100, NA, 200)
  plot <- BenchmarkPlot(data = result, plot_type = "efficiency")
  built <- ggplot2::ggplot_build(plot)
  expect_s3_class(plot, "ggplot")
  expect_match(built$data[[1]]$label[[1]], "No successful run")
  expect_match(built$data[[1]]$label[[1]], "both\\nruntime")
})
