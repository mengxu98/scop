time_methods <- list(
  leap = list(), scode = list(n_iter = 10L, n_runs = 2L, latent_dim = 2L),
  sincerities = list(max_time_points = 8L, SIGN = 0L)
)

for (method in names(time_methods)) {
  test_that(paste(method, "preserves native score and sign contracts"), {
    skip_if(Sys.getenv("SCOP_TEST_GRN_TIME") != "true", "Set SCOP_TEST_GRN_TIME=true for real R time backends")
    source_dir <- Sys.getenv("SCOP_GRN_SOURCE_DIR")
    if (method != "leap") {
      skip_if(!nzchar(source_dir), "Set SCOP_GRN_SOURCE_DIR to the external method sources")
    }
    for (package in switch(method,
      leap = "LEAP",
      scode = "MASS",
      sincerities = c("glmnet", "ppcor")
    )) {
      skip_if_not_installed(package)
    }
    set.seed(41)
    expression <- matrix(stats::rpois(8L * 80L, lambda = 4),
      nrow = 8L,
      dimnames = list(paste0("g", 1:8), paste0("c", 1:80))
    )
    time <- matrix(seq(0, 1, length.out = 80L),
      ncol = 1L,
      dimnames = list(colnames(expression), "trajectory")
    )
    input <- list(
      expression = expression, pseudotime = time,
      cells = colnames(expression), group = "all", regulators = rownames(expression)[1:4],
      targets = rownames(expression)[5:8]
    )
    resources <- if (method == "leap") list() else list(source_dir = source_dir)
    result <- RunGRN(expression, methods = method, pseudotime = time,
      regulators = input$regulators, targets = input$targets, resources = resources,
      method_args = time_methods[[method]], return_type = "result", seed = 42L, verbose = FALSE
    )
    shuffled_time <- time[rev(seq_len(nrow(time))), , drop = FALSE]
    shuffled <- RunGRN(expression, methods = method, pseudotime = shuffled_time,
      regulators = input$regulators, targets = input$targets, resources = resources,
      method_args = time_methods[[method]], return_type = "result", seed = 42L, verbose = FALSE
    )
    expect_identical(GetGRN(shuffled), GetGRN(result))
    fit <- result$groups$all
    expect_identical(fit$cells, input$cells)
    expect_identical(fit$execution$status, "completed")
    expect_null(fit$execution$process)
    expect_true(all(fit$native$edges$regulator %in% input$regulators))
    expect_true(all(fit$native$edges$target %in% input$targets))
    expect_true(all(c("weight", "signed_weight") %in% names(fit$native$edges)))
    expect_identical(fit$layer_specs$tf_gene$sign_column, if (method == "sincerities") NULL else "signed_weight")
    expect_true(is.data.frame(fit$native$edges))
  })
}
