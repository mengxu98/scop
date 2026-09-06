# Verified migration of two pinned metadata-only assets, not a general v2
# upgrader. Their numeric raw x/y input is identical under contracts 2 and 3;
# no image axes, factors, or display-unit conversion can affect these results.
library(scop)
stopifnot(get(".spatial_coordinate_contract_version", asNamespace("scop")) == 3L)
pins <- c(visium_human_pancreas_results_sub = "8049c6851567618e5222ad86cd6ce3ab",
          visium_human_pancreas_pair_sub = "2260392ad1765d084bc232f76323e888")
stamp <- function(x, from, to) {
  if (is.list(x)) for (i in seq_along(x)) {
    key <- names(x)[i]
    if (identical(key, "coordinate_contract_version")) {
      stopifnot(all(x[[i]] == from))
      x[i] <- list(rep.int(as.integer(to), length(x[[i]])))
    } else {
      if (identical(key, "coordinate_space")) stopifnot(identical(x[[i]], "raw"))
      x[i] <- list(stamp(x[[i]], from, to))
    }
  }
  for (key in c("coordinate_contract_version", "spatial_source", "spatial_transform")) {
    value <- attr(x, key, exact = TRUE)
    if (!is.null(value)) {
      if (key == "coordinate_contract_version") {
        stopifnot(identical(value, as.integer(from)))
        attr(x, key) <- as.integer(to)
      } else attr(x, key) <- stamp(value, from, to)
    }
  }
  x
}
for (name in names(pins)) {
  path <- file.path("data", paste0(name, ".rda"))
  stopifnot(unname(tools::md5sum(path)) == pins[[name]])
  env <- new.env(); load(path, env)
  before <- env[[name]]
  stopifnot(length(before@images) == 0L, is.numeric(before$x), is.numeric(before$y),
            all(is.finite(before$x)), all(is.finite(before$y)))
  current <- SpatialCoordinates(before, coord.cols = c("x", "y"))$data
  stopifnot(identical(current$cell_id, colnames(before)),
            identical(current$x, unname(before$x)), identical(current$y, unname(before$y)))
  after <- before
  after@tools <- stamp(before@tools, 2L, 3L)
  # Reverse the metadata-only edits to verify that no scientific payload,
  # cell identity, assay, coordinate, or original backend provenance changed.
  restored <- after; restored@tools <- stamp(after@tools, 3L, 2L)
  stopifnot(identical(restored, before))
  after@tools$SCOPExamples$coordinate_contract_migration <- list(
    from = 2L, to = 3L, source_md5 = pins[[name]],
    validation = "Pinned metadata-only numeric raw x/y; v2 and v3 input coordinates identical; all scientific payloads unchanged",
    backend_rerun = FALSE)
  env[[name]] <- after
  save(list = name, file = path, envir = env, compress = "xz")
}
