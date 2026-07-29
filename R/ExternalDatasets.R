#' @title List SCOP external datasets
#'
#' @description
#' Read a dataset manifest from the external `mengxu98/datasets` repository or
#' a local mirror. This keeps example data assets out of the SCOP package while
#' still making them discoverable and reproducible.
#'
#' @md
#' @param collection Dataset collection directory, for example `"Xenium"`.
#' @param datasets_base_url Base URL or local directory containing SCOP dataset
#' collections.
#'
#' @return A data frame parsed from `manifest.tsv`.
#' @export
#'
#' @examples
#' \dontrun{
#' ListExternalDatasets("Xenium")
#' }
ListExternalDatasets <- function(
  collection = "Xenium",
  datasets_base_url = "https://raw.githubusercontent.com/mengxu98/datasets/main"
) {
  validate_scalar_string(collection, "collection")
  validate_scalar_string(datasets_base_url, "datasets_base_url")
  dataset_read_manifest(
    collection = collection,
    datasets_base_url = datasets_base_url
  )
}

#' @title Load a SCOP external dataset
#'
#' @description
#' Download a dataset listed in `mengxu98/datasets`, validate its size and
#' sha256 checksum when the manifest provides them, cache it under
#' `tools::R_user_dir("scop", "data")`, and return the R object.
#'
#' @md
#' @param dataset Dataset id from the collection manifest.
#' @param collection Dataset collection directory, for example `"Xenium"`.
#' @param cache_dir Directory used to cache downloaded files. If `NULL`, uses
#' `tools::R_user_dir("scop", "data")/datasets/<collection>`.
#' @param datasets_base_url Base URL or local directory containing SCOP dataset
#' collections.
#' @param update Whether to redownload the file even when a valid cached copy is
#' available.
#' @param return_path Whether to return the cached file path instead of reading
#' the R object.
#' @param verbose Whether to print progress messages.
#'
#' @return The loaded R object, or a file path when `return_path = TRUE`.
#' @export
#'
#' @examples
#' \dontrun{
#' xenium <- LoadExternalDataset("xenium_human_pancreas_sub", collection = "Xenium")
#' SpatialSpotPlot(xenium, group.by = "nCount_Xenium")
#' }
LoadExternalDataset <- function(
  dataset,
  collection = "Xenium",
  cache_dir = NULL,
  datasets_base_url = "https://raw.githubusercontent.com/mengxu98/datasets/main",
  update = FALSE,
  return_path = FALSE,
  verbose = TRUE
) {
  validate_scalar_string(dataset, "dataset")
  validate_scalar_string(collection, "collection")
  validate_scalar_string(datasets_base_url, "datasets_base_url")
  cache_dir <- cache_dir %||%
    file.path(tools::R_user_dir("scop", "data"), "datasets", collection)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  manifest <- dataset_read_manifest(
    collection = collection,
    datasets_base_url = datasets_base_url
  )
  row <- dataset_select_row(manifest, dataset = dataset)
  path <- dataset_download_file(
    row = row,
    collection = collection,
    cache_dir = cache_dir,
    datasets_base_url = datasets_base_url,
    update = update,
    verbose = verbose
  )
  if (isTRUE(return_path)) {
    return(path)
  }
  readRDS(path)
}

dataset_read_manifest <- function(collection, datasets_base_url) {
  manifest_ref <- resource_ref(
    datasets_base_url,
    file.path(collection, "manifest.tsv")
  )
  manifest <- tryCatch(
    utils::read.delim(manifest_ref, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) {
      log_message(
        "Failed to read SCOP dataset manifest from {.val {manifest_ref}}: {.val {e$message}}",
        message_type = "error"
      )
    }
  )
  required <- c("dataset", "file", "object", "source", "sha256", "size_bytes")
  missing <- setdiff(required, colnames(manifest))
  if (length(missing) > 0L) {
    log_message(
      "SCOP dataset manifest is missing required column{?s}: {.val {missing}}",
      message_type = "error"
    )
  }
  if (nrow(manifest) == 0L) {
    log_message(
      "SCOP dataset manifest {.val {manifest_ref}} is empty",
      message_type = "error"
    )
  }
  manifest
}

dataset_select_row <- function(manifest, dataset) {
  idx <- which(manifest$dataset == dataset)
  if (length(idx) != 1L) {
    candidates <- manifest[
      ,
      intersect(c("dataset", "file", "object", "description"), colnames(manifest)),
      drop = FALSE
    ]
    log_message(
      paste0(
        "Could not find a unique SCOP dataset {.val {dataset}}. ",
        "Available datasets: {.val {utils::capture.output(print(candidates))}}"
      ),
      message_type = "error"
    )
  }
  manifest[idx, , drop = FALSE]
}

dataset_download_file <- function(
  row,
  collection,
  cache_dir,
  datasets_base_url,
  update = FALSE,
  verbose = TRUE
) {
  key <- row$file
  dest <- file.path(cache_dir, key)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (
    file.exists(dest) &&
      !isTRUE(update) &&
      dataset_validate_file(dest, sha256 = row$sha256, size = row$size_bytes)
  ) {
    log_message("Use cached SCOP dataset: {.file {dest}}", verbose = verbose)
    return(dest)
  }

  source <- resource_ref(
    datasets_base_url,
    file.path(collection, row$file)
  )
  log_message(
    "Download SCOP dataset {.file {key}} to {.file {cache_dir}}",
    verbose = verbose
  )
  tmp <- paste0(dest, ".tmp")
  if (file.exists(tmp)) {
    unlink(tmp)
  }
  old_timeout <- getOption("timeout", 60)
  options(timeout = max(600, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  status <- if (file.exists(source)) {
    file.copy(source, tmp, overwrite = TRUE)
  } else {
    identical(
      utils::download.file(
        url = source,
        destfile = tmp,
        mode = "wb",
        quiet = !isTRUE(verbose)
      ),
      0L
    )
  }
  if (!isTRUE(status) || !file.exists(tmp)) {
    unlink(tmp)
    log_message(
      "Failed to download SCOP dataset from {.val {source}}",
      message_type = "error"
    )
  }
  if (!dataset_validate_file(tmp, sha256 = row$sha256, size = row$size_bytes)) {
    unlink(tmp)
    log_message(
      "Downloaded SCOP dataset failed size or sha256 validation: {.file {key}}",
      message_type = "error"
    )
  }
  if (file.exists(dest)) {
    unlink(dest)
  }
  if (!file.rename(tmp, dest)) {
    unlink(tmp)
    log_message(
      "Failed to cache SCOP dataset at {.file {dest}}",
      message_type = "error"
    )
  }
  dest
}

dataset_validate_file <- function(path, sha256 = NA_character_, size = NA_real_) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  size <- suppressWarnings(as.numeric(size))
  if (!is.na(size) && file.info(path)$size != size) {
    return(FALSE)
  }
  if (!is.na(sha256) && nzchar(sha256)) {
    observed <- unname(tools::sha256sum(path))
    if (!identical(tolower(observed), tolower(sha256))) {
      return(FALSE)
    }
  }
  TRUE
}
