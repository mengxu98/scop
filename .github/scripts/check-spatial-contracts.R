source(file.path("R", "SpatialRegistry.R"), local = environment())

registry <- spatial_method_registry()
implementation_paths <- file.path("R", registry$implementation_files)
missing_files <- unique(implementation_paths[!file.exists(implementation_paths)])
if (length(missing_files) > 0L) {
  stop(
    "Spatial registry implementation files are missing: ",
    paste(missing_files, collapse = ", "),
    call. = FALSE
  )
}

namespace <- readLines("NAMESPACE", warn = FALSE)
exports <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", namespace, value = TRUE))
missing_exports <- setdiff(registry$method, exports)
if (length(missing_exports) > 0L) {
  stop(
    "Spatial registry methods are not exported: ",
    paste(missing_exports, collapse = ", "),
    call. = FALSE
  )
}

dist_calls <- function(file) {
  parsed <- parse(file, keep.source = TRUE)
  tokens <- utils::getParseData(parsed)
  calls <- tokens[
    tokens$token == "SYMBOL_FUNCTION_CALL" & tokens$text == "dist",
    c("line1", "col1"),
    drop = FALSE
  ]
  if (nrow(calls) == 0L) return(character())
  paste0(file, ":", calls$line1, ":", calls$col1)
}

sparse_files <- unique(
  implementation_paths[registry$scalability == "sparse_required"]
)
dense_calls <- unlist(lapply(sparse_files, dist_calls), use.names = FALSE)
if (length(dense_calls) > 0L) {
  stop(
    "Dense dist() calls remain in sparse-required spatial paths: ",
    paste(dense_calls, collapse = ", "),
    call. = FALSE
  )
}

contract_files <- file.path(
  "R",
  c(
    "SpatialRegistry.R", "RunSpatialNetwork.R", "SpatialCellPlot.R",
    "SpatialFrameworkConvert.R", "SpatialCore.R"
  )
)
contract_text <- unlist(lapply(contract_files, readLines, warn = FALSE), use.names = FALSE)
forbidden <- grep(
  "\\b(library|require|requireNamespace)\\s*\\(",
  contract_text,
  perl = TRUE,
  value = TRUE
)
if (length(forbidden) > 0L) {
  stop(
    "Forbidden package loading call in spatial contract files: ",
    paste(forbidden, collapse = " | "),
    call. = FALSE
  )
}

message("Spatial source contracts passed for ", nrow(registry), " methods.")
