grn_plot_fixture <- function(method = "fixture", directed = TRUE, signed = TRUE, membership = FALSE) {
  tf_gene <- data.frame(regulator = c("TF1", "TF2", "TF1"), target = c("G1", "G2", "G2"), weight = c(0.9, -0.8, 0.4), beta = c(0.2, -0.3, 0.1))
  tf_region <- data.frame(regulator = c("TF1", "TF2"), region = c("chr1-100-200", "chr1-300-400"), weight = c(0.8, -0.7), beta = c(0.2, -0.3))
  region_gene <- data.frame(region = c("chr1-100-200", "chr1-300-400"), target = c("G1", "G2"), weight = c(-0.6, 0.5), beta = c(-0.4, 0.3))
  triplets <- data.frame(regulator = c("TF1", "TF2"), region = c("chr1-100-200", "chr1-300-400"), target = c("G1", "G2"), weight = c(-0.9, -0.6))
  spec <- function(path, source, target, source_type, target_type) list(
    status = "computed", table_path = path, source = source, target = target,
    score = "weight", source_type = source_type, target_type = target_type,
    score_type = if (membership) "membership" else "coefficient",
    directed = directed, sign_column = if (signed && !membership) "beta" else NULL,
    sign_basis = if (signed && !membership) "coefficient" else "unknown", rank_by = "absolute_score"
  )
  structure(list(
    schema_version = "1.0", method = method, implementation = "fixture",
    package_versions = character(), parameters = list(), input_summary = list(),
    resources = list(), execution = list(status = "completed"),
    groups = list(A = list(
      cells = c("C1", "C2"), native = list(tf_gene = tf_gene, tf_region = tf_region, region_gene = region_gene, triplets = triplets),
      layer_specs = list(tf_gene = spec("tf_gene", "regulator", "target", if (directed) "TF" else "gene", "gene"), tf_region = spec("tf_region", "regulator", "region", "TF", "region"), region_gene = spec("region_gene", "region", "target", "region", "gene")),
      chains = list(status = "computed", table_path = "triplets"), gene_sets = list(status = "not_supported"),
      parameters = list(), implementation = "fixture", package_versions = character(), execution = list(status = "completed")
    ))
  ), class = c("GRNResult", "list"))
}

test_that("GRNPlot keeps signed native evidence and limits only its view", {
  result <- grn_plot_fixture()
  original <- serialize(result, NULL)
  plot <- GRNPlot(result, group = "A", top_n = 2)
  view <- attr(plot, "grn_data")
  expect_s3_class(plot, "ggplot")
  expect_equal(nrow(view$edges), 2L)
  expect_setequal(view$edges$sign, c(1L, -1L))
  expect_setequal(as.character(view$edges$evidence), c("Positive", "Negative"))
  expect_match(plot$labels$subtitle, "2 of 3")
  expect_identical(serialize(result, NULL), original)
  expect_silent(ggplot2::ggplot_build(plot))
})

test_that("layered views use each segment sign and actual native chains", {
  plot <- GRNPlot(grn_plot_fixture(), group = "A", type = "layered", top_n = 1)
  view <- attr(plot, "grn_data")
  expect_equal(nrow(view$chains), 1L)
  expect_equal(nrow(view$edges), 2L)
  expect_equal(view$edges$sign[view$edges$layer == "tf_region"], 1L)
  expect_equal(view$edges$sign[view$edges$layer == "region_gene"], -1L)
  expect_setequal(view$nodes$node_type, c("TF", "region", "gene"))
  expect_equal(length(unique(view$nodes$x)), 3L)
  expect_silent(ggplot2::ggplot_build(plot))
  invalid <- grn_plot_fixture()
  invalid$groups$A$native$triplets$regulator[[1L]] <- "Absent"
  expect_error(GRNPlot(invalid, group = "A", type = "layered"), "absent")
})

test_that("unknown, membership and undirected views retain their semantics", {
  unknown <- GRNPlot(grn_plot_fixture(signed = FALSE), group = "A")
  expect_true(all(attr(unknown, "grn_data")$edges$evidence == "Unknown"))
  membership <- GRNPlot(grn_plot_fixture(membership = TRUE), group = "A", type = "layered")
  expect_true(all(attr(membership, "grn_data")$edges$evidence == "Membership"))
  expect_true(all(is.na(attr(membership, "grn_data")$edges$sign)))
  undirected <- GRNPlot(grn_plot_fixture(directed = FALSE), group = "A")
  expect_true(all(!attr(undirected, "grn_data")$edges$directed))
  expect_null(undirected$layers[[1L]]$geom_params$arrow)
  expect_equal(nrow(undirected$layers[[2L]]$data), 0L)
  expect_silent(ggplot2::ggplot_build(undirected))
})

test_that("overlap uses a common entity scope and never adds method scores", {
  first <- grn_plot_fixture(method = "one")
  second <- grn_plot_fixture(method = "two")
  second$groups$A$native$tf_gene$weight <- c(10, -9, 3)
  plot <- GRNPlot(first, group = "A", type = "overlap", compare = second, top_n = 2)
  view <- attr(plot, "grn_data")
  expect_true(all(view$overlap$shared_edges == 2L))
  expect_true(all(view$counts$shown == 2L))
  expect_true(all(view$counts$total == 3L))
  second <- grn_plot_fixture(method = "two", directed = FALSE)
  separate <- attr(GRNPlot(first, group = "A", type = "overlap", compare = second), "grn_data")$overlap
  expect_true(all(separate$shared_edges[as.character(separate$method_1) != as.character(separate$method_2)] == 0L))
  first <- grn_plot_fixture(method = "one", directed = FALSE)
  second <- grn_plot_fixture(method = "two", directed = FALSE)
  second$groups$A$native$tf_gene[c("regulator", "target")] <- second$groups$A$native$tf_gene[c("target", "regulator")]
  reversed <- attr(GRNPlot(first, group = "A", type = "overlap", compare = second), "grn_data")$overlap
  expect_true(all(reversed$shared_edges == 3L))
  expect_silent(ggplot2::ggplot_build(plot))
})

test_that("undirected duplicates do not use up the plotting limit", {
  result <- grn_plot_fixture(directed = FALSE)
  result$groups$A$native$tf_gene <- data.frame(
    regulator = c("A", "D", "B"), target = c("D", "A", "A"),
    weight = c(-0.9, -0.9, 0.8), beta = c(-0.9, -0.9, 0.8)
  )
  edges <- attr(GRNPlot(result, top_n = 2L), "grn_data")$edges
  expect_equal(edges$raw_score, c(-0.9, 0.8))
  expect_equal(nrow(GetGRN(result)), 3L)
  view <- attr(GRNPlot(result, type = "overlap", compare = result, top_n = 2L), "grn_data")
  expect_true(all(view$counts$total == 2L))
  expect_true(all(view$counts$shown == 2L))
  expect_true(all(view$overlap$shared_edges == 2L))
})

test_that("overlap matches endpoints even when identifiers contain separators", {
  first <- grn_plot_fixture(method = "one")
  second <- grn_plot_fixture(method = "two")
  first$groups$A$native$tf_gene[c("regulator", "target")] <- list(
    c("a", "a\rgene\rb", "p"), c("b\rgene\rc", "q", "c")
  )
  second$groups$A$native$tf_gene[c("regulator", "target")] <- list(
    c("a\rgene\rb", "a", "p"), c("c", "q", "b\rgene\rc")
  )
  overlap <- attr(GRNPlot(first, type = "overlap", compare = second), "grn_data")$overlap
  expect_equal(overlap$shared_edges[as.character(overlap$method_1) != as.character(overlap$method_2)], c(0L, 0L))
})

test_that("empty plots and saved results remain drawable", {
  result <- grn_plot_fixture()
  result$groups$A$native$tf_gene <- result$groups$A$native$tf_gene[FALSE, ]
  plot <- GRNPlot(result, group = "A")
  expect_equal(nrow(attr(plot, "grn_data")$edges), 0L)
  expect_silent(ggplot2::ggplot_build(plot))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(grn_plot_fixture(), path)
  original <- GRNPlot(grn_plot_fixture(), group = "A")
  restored <- GRNPlot(readRDS(path), group = "A")
  expect_identical(attr(original, "grn_data"), attr(restored, "grn_data"))
  set.seed(91)
  seed <- .Random.seed
  GRNPlot(grn_plot_fixture(), group = "A", layout = "fr")
  expect_identical(.Random.seed, seed)
  expect_error(GRNPlot(grn_plot_fixture(), top_n = 0), "top_n")
  expect_error(GRNPlot(grn_plot_fixture(), type = "overlap"), "compare")
})
