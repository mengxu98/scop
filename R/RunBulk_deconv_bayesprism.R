RunBulk_deconv_bayesprism <- function(
  count_matrix,
  reference_matrix,
  backend = c("internal", "native"),
  transform = "none",
  min_overlap = 10,
  verbose = TRUE
) {
  backend <- match.arg(backend)
  if (identical(backend, "native")) {
    return(.bulk_method_failed(
      reason = paste(
        "The native BayesPrism backend is not implemented yet.",
        "Use backend = 'internal' for SCOP internal profile fitting."
      ),
      parameters = list(backend = backend)
    ))
  }
  log_message(
    "{.val deconv_BayesPrism} currently uses SCOP internal profile fitting.",
    verbose = verbose
  )

  prop_matrix <- tryCatch(
    .bulk_fit_proportions(
      count_matrix = count_matrix,
      reference_matrix = reference_matrix,
      transform = transform,
      min_overlap = min_overlap
    ),
    error = function(e) e
  )
  if (inherits(prop_matrix, "error")) {
    return(.bulk_method_failed(reason = prop_matrix$message))
  }

  .bulk_method_success(
    results = .bulk_prop_matrix_to_long(prop_matrix, method_name = "BayesPrism"),
    details = list(
      engine = "internal_profile_fit",
      package_backend = FALSE,
      note = paste(
        "This runner currently uses SCOP internal profile fitting rather",
        "than BayesPrism's native model."
      ),
      proportion_matrix = prop_matrix
    ),
    parameters = list(
      method = "BayesPrism",
      backend = backend,
      transform = transform,
      min_overlap = min_overlap
    )
  )
}
