time_methods <- list(leap = list(), scode = list(n_iter = 10L, n_runs = 2L, latent_dim = 2L),
                     sincerities = list(max_time_points = 8L, SIGN = 0L))

for (method in names(time_methods)) {
  test_that(paste(method, "preserves native score and sign contracts"), {
    skip_if(Sys.getenv("SCOP_TEST_GRN_TIME") != "true", "Set SCOP_TEST_GRN_TIME=true for real R time backends")
    source_dir <- Sys.getenv("SCOP_GRN_SOURCE_DIR")
    if (method != "leap") {
      skip_if(!nzchar(source_dir), "Set SCOP_GRN_SOURCE_DIR to the external method sources")
    }
    for (package in grn_runtime_profiles()[[method]]$packages) skip_if_not_installed(package)
    set.seed(41)
    expression <- matrix(stats::rpois(8L * 80L, lambda = 4), nrow = 8L,
                         dimnames = list(paste0("g", 1:8), paste0("c", 1:80)))
    time <- matrix(seq(0, 1, length.out = 80L), ncol = 1L,
                   dimnames = list(colnames(expression), "trajectory"))
    input <- list(expression = expression, pseudotime = time,
                  cells = colnames(expression), group = "all", regulators = rownames(expression)[1:4],
                  targets = rownames(expression)[5:8])
    resources <- if (method == "leap") list() else list(source_dir = source_dir)
    config <- list(method = method, backend = "r", method_args = time_methods[[method]], resources = resources,
                   runtime = list(rscript = file.path(R.home("bin"), "Rscript")),
                   work_dir = tempfile(paste0("grn-", method, "-")), cores = 1L,
                   seed = 42L, force = FALSE, verbose = FALSE)
    fit <- grn_run_time(input, config)
    expect_identical(fit$cells, input$cells)
    expect_identical(fit$execution$status, "completed")
    expect_true(all(fit$native$edges$regulator %in% input$regulators))
    expect_true(all(fit$native$edges$target %in% input$targets))
    expect_true(all(c("weight", "signed_weight") %in% names(fit$native$edges)))
    expect_identical(fit$layer_specs$tf_gene$sign_column, if (method == "sincerities") NULL else "signed_weight")
    expect_true(is.data.frame(fit$native$edges))
  })
}
