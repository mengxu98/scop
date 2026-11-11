test_that("motif logos use explicit IDs from the selected assay", {
  x <- make_atac_test_object()
  incidence <- Matrix::Matrix(matrix(1, nrow(x), 1, dimnames = list(rownames(x), "M1")), sparse = TRUE)
  pwm <- list(M1 = matrix(seq_len(20), 4, 5, dimnames = list(c("A", "C", "G", "T"), NULL)))
  Signac::Motifs(x[["peaks"]]) <- Signac::CreateMotifObject(data = incidence, pwm = pwm)
  expect_error(ATACMotifPlot(x, "unknown"), "unique motif IDs")
  expect_error(ATACMotifPlot(x, c("M1", "M1")), "unique motif IDs")
  skip_if_not_installed("ggseqlogo")
  if (utils::packageVersion("ggseqlogo") == "0.2.2") {
    withr::local_options(lifecycle_verbosity = "warning")
    withCallingHandlers(
      {
        actual <- ATACMotifPlot(x, "M1", use.names = FALSE)
        expected <- Signac::MotifPlot(x, "M1", use.names = FALSE)
      },
      warning = function(w) {
        message <- conditionMessage(w)
        if (inherits(w, "lifecycle_warning_deprecated") &&
          grepl("argument of .guides.*cannot be .FALSE.|aes_string.*deprecated", message)) {
          invokeRestart("muffleWarning")
        }
      }
    )
  } else {
    actual <- ATACMotifPlot(x, "M1", use.names = FALSE)
    expected <- Signac::MotifPlot(x, "M1", use.names = FALSE)
  }
  expect_equal(actual$data, expected$data)
})
