test_that("cross kNN matches the raw cosine and euclidean rankings", {
  reference <- rbind(
    c(1.2, 0.1, 0.4),
    c(0.2, 1.3, 0.6),
    c(0.3, 0.5, 1.4),
    c(1.1, 1.2, 0.2),
    c(0.7, 0.4, 0.9)
  )
  query <- rbind(
    c(1.0, 0.2, 0.3),
    c(0.4, 1.1, 0.5),
    c(0.2, 0.6, 1.2)
  )

  for (metric in c("cosine", "euclidean")) {
    raw_distance <- if (identical(metric, "cosine")) {
      1 - proxyC::simil(reference, query, method = metric, use_nan = TRUE)
    } else {
      proxyC::dist(reference, query, method = metric, use_nan = TRUE)
    }
    expected <- thisutils::run_dense_topk_by_column(
      as.matrix(raw_distance),
      k = 3L,
      decreasing = FALSE
    )
    actual <- cross_knn_f32(
      reference = reference,
      query = query,
      k = 3L,
      metric = metric,
      cores = 2L
    )

    expect_equal(actual[["idx"]], expected[["idx"]])
    expect_equal(actual[["distance"]], expected[["value"]], tolerance = 1e-5)
  }
})

test_that("native cross kNN supports exact Pearson and Spearman rankings", {
  reference <- rbind(
    c(1, 0, 3),
    c(0, 3, 1),
    c(2, 1, 3),
    c(4, 1, 0)
  )
  query <- rbind(c(1, 0, 2), c(0, 2, 1))

  expect_null(knn_cross_topk_native(
    Matrix::Matrix(reference, sparse = TRUE), query, 2L, "cosine"
  ))
  expect_null(knn_cross_topk_native(
    reference, rbind(c(0, 0, 0), query), 2L, "cosine"
  ))

  for (metric in c("pearson", "spearman")) {
    expected_input <- if (identical(metric, "spearman")) {
      list(
        reference = t(apply(reference, 1L, rank)),
        query = t(apply(query, 1L, rank))
      )
    } else {
      list(reference = reference, query = query)
    }
    expected <- thisutils::run_dense_topk_by_column(
      1 - proxyC::simil(
        expected_input$reference,
        expected_input$query,
        method = "correlation",
        use_nan = TRUE
      ),
      k = 2L,
      decreasing = FALSE
    )
    actual <- knn_cross_topk_native(reference, query, 2L, metric)
    expect_equal(actual[["idx"]], expected[["idx"]])
    expect_equal(actual[["distance"]], expected[["value"]], tolerance = 1e-5)
  }

  expect_null(knn_cross_topk_native(
    rbind(reference, c(1, 1, 1)), query, 2L, "pearson"
  ))
})

test_that("fast row ranks install matrixStats and preserve Spearman preprocessing", {
  x <- rbind(
    c(3, 1, 3, 2),
    c(0, 0, 2, 1),
    c(-1, 4, 4, 4)
  )
  calls <- list()
  testthat::local_mocked_bindings(
    check_r = function(packages, verbose = TRUE, ...) {
      calls[[length(calls) + 1L]] <<- list(
        packages = packages,
        verbose = verbose
      )
      TRUE
    }
  )

  expect_identical(knn_rank_rows(x), t(apply(x, 1L, rank)))
  sparse_x <- methods::as(Matrix::Matrix(x, sparse = TRUE), "dgCMatrix")
  expect_identical(knn_rank_rows(sparse_x), t(apply(sparse_x, 1L, rank)))
  expect_length(calls, 2L)
  expect_true(all(vapply(calls, `[[`, character(1), "packages") == "matrixStats"))
  expect_true(all(!vapply(calls, `[[`, logical(1), "verbose")))
})

test_that("cross kNN preserves raw tie ordering", {
  reference <- rbind(c(1, 0), c(0, 1), c(1, 1), c(1, 0))
  query <- rbind(c(1, 0), c(0, 1))

  for (metric in c("cosine", "euclidean")) {
    raw_distance <- if (identical(metric, "cosine")) {
      1 - proxyC::simil(reference, query, method = metric, use_nan = TRUE)
    } else {
      proxyC::dist(reference, query, method = metric, use_nan = TRUE)
    }
    expected <- thisutils::run_dense_topk_by_column(
      as.matrix(raw_distance),
      k = 4L,
      decreasing = FALSE
    )
    actual <- cross_knn_f32(reference, query, 4L, metric, cores = 2L)

    expect_identical(actual[["idx"]], expected[["idx"]])
  }
})
