#' @title Run tAge transcriptomic aging-clock prediction
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param object A `Seurat` object, `ExpressionSet`, or expression matrix-like
#' object with genes in rows and pseudobulk/bulk samples in columns.
#' @param model_paths Named list or character vector of local tAge model files
#' (`.rds` for `backend = "cpp"`, `.pkl` for `backend = "python"`).
#' Names must match one or more of `"scaled"`, `"scaled_diff"`, `"yugene"`,
#' and `"yugene_diff"`. If `NULL`, model files are downloaded from
#' `mengxu98/datasets` for the C++ backend or Zenodo for the Python backend.
#' @param backend Model prediction backend. `"cpp"` uses converted Elastic Net
#' models and C++ prediction. `"python"` uses `tAge::predict_tAge()`.
#' @param clock tAge clock family used when `model_paths = NULL`.
#' @param model_species Model species scope used when `model_paths = NULL`.
#' `"auto"` chooses `"Mouse"` for `species = "mouse"`, `"Rodents"` for
#' `species = "rat"`, and `"Multispecies"` otherwise.
#' @param model_tissue Model tissue scope used when `model_paths = NULL`.
#' @param model_preprocessing Preprocessing-specific model file(s) used when
#' `model_paths = NULL`. Values map to tAge model filename suffixes.
#' @param model_cache_dir Directory used to cache downloaded tAge models. If
#' `NULL`, uses `tools::R_user_dir("scop", "data")`.
#' @param datasets_base_url Base URL or local directory containing converted
#' tAge EN models from `mengxu98/datasets`.
#' @param zenodo_record Zenodo record ID for tAge model files.
#' @param max_model_size Maximum model file size, in bytes, allowed for
#' automatic download. Increase this or use `Inf` for large Bayesian Ridge
#' models.
#' @param metadata Sample metadata for matrix input. Row names must match
#' columns of `object`. If `NULL`, minimal `sample_id` metadata is created.
#' @param group.by Metadata columns used to form pseudobulk groups for `Seurat`
#' input. If `NULL`, cells are aggregated sequentially across the whole object.
#' `split.by` and `control_group_column` are automatically included in the
#' aggregation metadata when provided.
#' @param split.by Optional pseudobulk metadata column used to run tAge
#' separately by group, such as tissue.
#' @param assay,layer Assay and layer used as raw counts for `Seurat` input.
#' @param species Species passed to tAge. Supported values are `"mouse"`,
#' `"human"`, `"rat"`, and `"monkey"`.
#' @param mode tAge model mode: `"EN"` for Elastic Net or `"BR"` for Bayesian
#' Ridge.
#' @param gene_mapping_type Gene identifier type passed to tAge preprocessing:
#' `"Gene.Symbol"` or `"Ensembl"`.
#' @param coverage_threshold Minimum cumulative read count per pseudobulk
#' sample for `Seurat` input.
#' @param count_threshold,percent_threshold Gene filtering thresholds passed to
#' `tAge::tAge_preprocessing()`.
#' @param control_group_column,control_group_label Optional control group used
#' by tAge control subtraction.
#' @param remove_outliers Whether to run `tAge::remove_outliers()` on the
#' pseudobulk `ExpressionSet`.
#' @param outlier.by Optional metadata column used to split samples before
#' outlier detection. Defaults to `split.by`.
#' @param outlier_n_components,outlier_threshold_quantile,outlier_min_samples
#' Parameters passed to `tAge::remove_outliers()`.
#' @param min_samples Minimum pseudobulk samples per `split.by` group for
#' `tAge::tAge_by_group()`.
#' @param shuffle,seed Whether to shuffle cells during pseudobulk aggregation
#' and the random seed used for that shuffle.
#' @param check_python Whether to verify that the active `reticulate` Python can
#' import `joblib`, `pandas`, and `sklearn` before tAge prediction.
#' @param tool_name Name of the Seurat tool entry used to store results.
#' @param store_eset,store_processed Whether to store the pseudobulk
#' `ExpressionSet` and processed tAge `ExpressionSet` list in the returned
#' Seurat tool entry or list result.
#'
#' @return A `Seurat` object with tAge results stored in `object@tools`, or a
#' list with `predictions`, metadata, and parameters for non-Seurat input.
#' @export
#'
#' @references
#' Tyshkovskiy, A., Glubokov, D., Moliere, A., et al. (2026). Universal
#' transcriptomic hallmarks of mammalian ageing and mortality.
#' \emph{Nature}. \doi{10.1038/s41586-026-10542-3}
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- RuntAge(
#'   pancreas_sub,
#'   group.by = "CellType",
#'   species = "mouse",
#'   mode = "EN"
#' )
#' head(pancreas_sub@tools$tAge$predictions)
#' }
RuntAge <- function(
  object,
  model_paths = NULL,
  backend = c("python", "cpp"),
  clock = c("Chronoage", "NormalizedAge", "Mortality"),
  model_species = c("auto", "Multispecies", "Mouse", "Rodents"),
  model_tissue = "Multitissue",
  model_preprocessing = "scaled_diff",
  model_cache_dir = NULL,
  datasets_base_url = "https://raw.githubusercontent.com/mengxu98/datasets/main/tAge/EN",
  zenodo_record = "18763485",
  max_model_size = 1024^3,
  metadata = NULL,
  group.by = NULL,
  split.by = NULL,
  assay = NULL,
  layer = "counts",
  species = c("mouse", "human", "rat", "monkey"),
  mode = c("EN", "BR"),
  gene_mapping_type = c("Gene.Symbol", "Ensembl"),
  coverage_threshold = 1e6,
  count_threshold = 10,
  percent_threshold = 20,
  control_group_column = NULL,
  control_group_label = NULL,
  remove_outliers = FALSE,
  outlier.by = split.by,
  outlier_n_components = 10,
  outlier_threshold_quantile = 0.99,
  outlier_min_samples = 10,
  min_samples = 5,
  shuffle = FALSE,
  seed = NULL,
  check_python = TRUE,
  tool_name = "tAge",
  store_eset = FALSE,
  store_processed = FALSE,
  verbose = TRUE
) {
  species <- match.arg(species)
  mode <- match.arg(mode)
  backend <- match.arg(backend)
  clock <- match.arg(clock)
  model_species <- match.arg(model_species)
  gene_mapping_type <- match.arg(gene_mapping_type)
  if (identical(backend, "cpp") && !identical(mode, "EN")) {
    log_message(
      "The C++ tAge backend currently supports only {.val EN} models; use {.code backend = 'python'} for {.val BR}.",
      message_type = "error"
    )
  }

  if (!is.numeric(coverage_threshold) || length(coverage_threshold) != 1L || is.na(coverage_threshold) || coverage_threshold <= 0) {
    log_message(
      "{.arg coverage_threshold} must be a positive number",
      message_type = "error"
    )
  }
  if (!is.character(tool_name) || length(tool_name) != 1L || !nzchar(tool_name)) {
    log_message(
      "{.arg tool_name} must be a non-empty string",
      message_type = "error"
    )
  }

  check_r(c("Gladyshev-Lab/tAge", "Biobase", "edgeR"), verbose = FALSE)
  if (identical(backend, "python")) {
    PrepareEnv(modules = "tage")
    configure_python_thread_env()
    if (isTRUE(check_python)) {
      check_tage_python(verbose = verbose)
    }
  }

  is_seurat <- inherits(object, "Seurat")
  if (is.null(model_paths)) {
    model_paths <- if (identical(backend, "cpp")) {
      fetch_tage_r_model_paths(
        mode = mode,
        clock = clock,
        model_species = model_species,
        species = species,
        model_tissue = model_tissue,
        model_preprocessing = model_preprocessing,
        model_cache_dir = model_cache_dir,
        datasets_base_url = datasets_base_url,
        verbose = verbose
      )
    } else {
      fetch_tage_model_paths(
        mode = mode,
        clock = clock,
        model_species = model_species,
        species = species,
        model_tissue = model_tissue,
        model_preprocessing = model_preprocessing,
        model_cache_dir = model_cache_dir,
        zenodo_record = zenodo_record,
        max_model_size = max_model_size,
        verbose = verbose
      )
    }
  } else {
    model_paths <- normalize_tage_model_paths(model_paths)
  }

  log_message(
    "Run {.pkg tAge} with {.val {length(model_paths)}} model{?s} using {.val {backend}} backend",
    message_type = "running",
    verbose = verbose
  )

  eset <- make_tage_eset(
    object = object,
    metadata = metadata,
    group.by = group.by,
    split.by = split.by,
    assay = assay,
    layer = layer,
    coverage_threshold = coverage_threshold,
    control_group_column = control_group_column,
    shuffle = shuffle,
    seed = seed,
    verbose = verbose
  )

  if (isTRUE(remove_outliers)) {
    if (!is.null(outlier.by) && !outlier.by %in% colnames(Biobase::pData(eset))) {
      log_message(
        "{.arg outlier.by} must be a column in tAge pseudobulk metadata",
        message_type = "error"
      )
    }
    eset <- get_namespace_fun("tAge", "remove_outliers")(
      eset,
      n_components = outlier_n_components,
      threshold_quantile = outlier_threshold_quantile,
      split_by = outlier.by,
      min_samples = outlier_min_samples,
      verbose = verbose
    )
  }

  processed <- NULL
  if (!is.null(split.by)) {
    if (!split.by %in% colnames(Biobase::pData(eset))) {
      log_message(
        "{.arg split.by} must be a column in tAge pseudobulk metadata",
        message_type = "error"
      )
    }
    predictions <- run_tage_by_group_scop(
      eset = eset,
      split_by = split.by,
      model_paths = model_paths,
      species = species,
      mode = mode,
      backend = backend,
      gene_mapping_type = gene_mapping_type,
      control_group_column = control_group_column,
      control_group_label = control_group_label,
      count_threshold = count_threshold,
      percent_threshold = percent_threshold,
      min_samples = min_samples,
      verbose = verbose
    )
  } else {
    processed <- get_namespace_fun("tAge", "tAge_preprocessing")(
      eset = eset,
      species = species,
      gene_mapping_type = gene_mapping_type,
      verbose = verbose,
      control_group_column = control_group_column,
      control_group_label = control_group_label,
      count_threshold = count_threshold,
      percent_threshold = percent_threshold
    )
    predictions <- predict_tage_models(
      processed = processed,
      model_paths = model_paths,
      species = species,
      mode = mode,
      backend = backend
    )
  }

  predictions <- as.data.frame(predictions, check.names = FALSE)
  parameters <- list(
    group.by = group.by,
    split.by = split.by,
    assay = assay %||% if (is_seurat) SeuratObject::DefaultAssay(object) else NULL,
    layer = layer,
    species = species,
    mode = mode,
    backend = backend,
    clock = clock,
    model_species = model_species,
    model_tissue = model_tissue,
    model_preprocessing = model_preprocessing,
    datasets_base_url = datasets_base_url,
    gene_mapping_type = gene_mapping_type,
    model_paths = model_paths,
    coverage_threshold = coverage_threshold,
    count_threshold = count_threshold,
    percent_threshold = percent_threshold,
    control_group_column = control_group_column,
    control_group_label = control_group_label,
    remove_outliers = remove_outliers,
    outlier.by = outlier.by,
    min_samples = min_samples,
    shuffle = shuffle,
    seed = seed
  )

  result <- list(
    predictions = predictions,
    pseudobulk_metadata = Biobase::pData(eset),
    parameters = parameters
  )
  if (isTRUE(store_eset)) {
    result$pseudobulk_eset <- eset
  }
  if (isTRUE(store_processed) && !is.null(processed)) {
    result$processed <- processed
  }

  if (is_seurat) {
    object@tools[[tool_name]] <- result
    object <- Seurat::LogSeuratCommand(object = object)
    log_message(
      "{.pkg tAge} predictions stored in {.code object@tools[[{tool_name}]]}",
      message_type = "success",
      verbose = verbose
    )
    return(object)
  }

  result
}

normalize_tage_model_paths <- function(model_paths) {
  if (is.null(model_paths) || length(model_paths) == 0L) {
    log_message(
      "{.arg model_paths} must provide at least one local tAge model file",
      message_type = "error"
    )
  }
  model_paths <- as.list(model_paths)
  model_paths <- lapply(model_paths, function(path) {
    if (length(path) != 1L || is.na(path) || !nzchar(path)) {
      log_message(
        "Each {.arg model_paths} entry must be one non-empty file path",
        message_type = "error"
      )
    }
    normalizePath(path, mustWork = FALSE)
  })
  missing_files <- vapply(model_paths, function(path) !file.exists(path), logical(1))
  if (any(missing_files)) {
    log_message(
      "Missing tAge model file{?s}: {.file {unlist(model_paths[missing_files])}}",
      message_type = "error"
    )
  }

  valid_names <- c("scaled", "scaled_diff", "yugene", "yugene_diff")
  model_names <- names(model_paths)
  if (is.null(model_names) || any(!nzchar(model_names))) {
    model_names <- vapply(model_paths, infer_tage_model_name, character(1))
  }
  invalid_names <- setdiff(model_names, valid_names)
  if (length(invalid_names) > 0L || any(is.na(model_names))) {
    log_message(
      "{.arg model_paths} names must be one or more of {.val {valid_names}}",
      message_type = "error"
    )
  }
  names(model_paths) <- model_names
  model_paths
}

fetch_tage_r_model_paths <- function(
  mode,
  clock,
  model_species = "auto",
  species = "mouse",
  model_tissue = "Multitissue",
  model_preprocessing = "scaled_diff",
  model_cache_dir = NULL,
  datasets_base_url = "https://raw.githubusercontent.com/mengxu98/datasets/main/tAge/EN",
  verbose = TRUE
) {
  if (!identical(mode, "EN")) {
    log_message(
      "Converted R-native tAge models are available only for {.val EN}",
      message_type = "error"
    )
  }
  model_species <- resolve_tage_model_species(
    model_species = model_species,
    species = species
  )
  model_preprocessing <- normalize_tage_model_preprocessing(model_preprocessing)
  model_cache_dir <- model_cache_dir %||%
    file.path(tools::R_user_dir("scop", "data"), "tAge", "datasets-EN")
  dir.create(model_cache_dir, recursive = TRUE, showWarnings = FALSE)

  manifest <- read_tage_datasets_manifest(datasets_base_url = datasets_base_url)
  model_paths <- lapply(model_preprocessing, function(preprocessing) {
    row <- select_tage_r_model_row(
      manifest = manifest,
      mode = mode,
      clock = clock,
      model_species = model_species,
      model_tissue = model_tissue,
      preprocessing = preprocessing
    )
    download_tage_r_model_file(
      row = row,
      model_cache_dir = model_cache_dir,
      datasets_base_url = datasets_base_url,
      verbose = verbose
    )
  })
  names(model_paths) <- model_preprocessing
  normalize_tage_model_paths(model_paths)
}

read_tage_datasets_manifest <- function(datasets_base_url) {
  manifest_ref <- tage_resource_ref(datasets_base_url, "manifest.tsv")
  tryCatch(
    utils::read.delim(manifest_ref, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) {
      log_message(
        "Failed to read tAge EN manifest from {.val {manifest_ref}}: {.val {e$message}}",
        message_type = "error"
      )
    }
  )
}

select_tage_r_model_row <- function(
  manifest,
  mode,
  clock,
  model_species,
  model_tissue,
  preprocessing
) {
  idx <- which(
    manifest$mode == mode &
      manifest$clock == clock &
      manifest$model_species == model_species &
      manifest$tissue == model_tissue &
      manifest$preprocessing == preprocessing
  )
  if (length(idx) != 1L) {
    candidates <- manifest[
      manifest$mode == mode &
        manifest$clock == clock &
        manifest$preprocessing == preprocessing,
      c("rds_file", "model_species", "tissue", "preprocessing"),
      drop = FALSE
    ]
    log_message(
      paste0(
        "Could not find a unique converted tAge EN model for ",
        "{.val {clock}}, {.val {model_species}}, {.val {model_tissue}}, ",
        "{.val {preprocessing}}. Matching candidates: {.val {utils::capture.output(print(candidates))}}"
      ),
      message_type = "error"
    )
  }
  manifest[idx, , drop = FALSE]
}

download_tage_r_model_file <- function(
  row,
  model_cache_dir,
  datasets_base_url,
  verbose = TRUE
) {
  key <- basename(row$rds_file)
  dest <- file.path(model_cache_dir, key)
  expected_size <- row$rds_size
  expected_md5 <- row$rds_md5
  if (file.exists(dest) && validate_tage_cached_file(dest, checksum = paste0("md5:", expected_md5), size = expected_size)) {
    log_message(
      "Use cached converted tAge model: {.file {dest}}",
      verbose = verbose
    )
    return(dest)
  }
  source <- tage_resource_ref(datasets_base_url, row$rds_file)
  log_message(
    "Download converted tAge model {.file {key}} to {.file {model_cache_dir}}",
    verbose = verbose
  )
  tmp <- paste0(dest, ".tmp")
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
      "Failed to download converted tAge model from {.val {source}}",
      message_type = "error"
    )
  }
  if (!validate_tage_cached_file(tmp, checksum = paste0("md5:", expected_md5), size = expected_size)) {
    unlink(tmp)
    log_message(
      "Downloaded converted tAge model failed size or checksum validation: {.file {key}}",
      message_type = "error"
    )
  }
  if (file.exists(dest)) {
    unlink(dest)
  }
  file.rename(tmp, dest)
  dest
}

tage_resource_ref <- function(base, path) {
  path <- gsub("^/+", "", path)
  if (dir.exists(base)) {
    return(file.path(base, path))
  }
  paste0(sub("/+$", "", base), "/", path)
}

fetch_tage_model_paths <- function(
  mode,
  clock,
  model_species = "auto",
  species = "mouse",
  model_tissue = "Multitissue",
  model_preprocessing = "scaled_diff",
  model_cache_dir = NULL,
  zenodo_record = "18763485",
  max_model_size = 1024^3,
  verbose = TRUE
) {
  model_species <- resolve_tage_model_species(
    model_species = model_species,
    species = species
  )
  model_preprocessing <- normalize_tage_model_preprocessing(model_preprocessing)
  model_cache_dir <- model_cache_dir %||%
    file.path(tools::R_user_dir("scop", "data"), "tAge", paste0("zenodo-", zenodo_record))
  dir.create(model_cache_dir, recursive = TRUE, showWarnings = FALSE)

  metadata <- fetch_tage_zenodo_metadata(zenodo_record = zenodo_record)
  files <- metadata[["files"]]
  if (length(files) == 0L) {
    log_message(
      "No files were found in Zenodo record {.val {zenodo_record}}",
      message_type = "error"
    )
  }

  model_paths <- lapply(model_preprocessing, function(preprocessing) {
    file_info <- select_tage_model_file(
      files = files,
      mode = mode,
      clock = clock,
      model_species = model_species,
      model_tissue = model_tissue,
      preprocessing = preprocessing
    )
    download_tage_model_file(
      file_info = file_info,
      model_cache_dir = model_cache_dir,
      max_model_size = max_model_size,
      verbose = verbose
    )
  })
  names(model_paths) <- model_preprocessing
  normalize_tage_model_paths(model_paths)
}

resolve_tage_model_species <- function(model_species, species) {
  if (!identical(model_species, "auto")) {
    return(model_species)
  }
  switch(species,
    mouse = "Mouse",
    rat = "Rodents",
    human = "Multispecies",
    monkey = "Multispecies",
    "Multispecies"
  )
}

normalize_tage_model_preprocessing <- function(model_preprocessing) {
  model_preprocessing <- unique(as.character(model_preprocessing))
  valid_names <- c("scaled_diff", "yugene_diff")
  invalid <- setdiff(model_preprocessing, valid_names)
  if (length(invalid) > 0L) {
    log_message(
      "{.arg model_preprocessing} must be one or more of {.val {valid_names}}",
      message_type = "error"
    )
  }
  model_preprocessing
}

fetch_tage_zenodo_metadata <- function(zenodo_record = "18763485") {
  check_r("jsonlite", verbose = FALSE)
  url <- paste0("https://zenodo.org/api/records/", zenodo_record)
  metadata <- tryCatch(
    {
      con <- url(url, open = "rb")
      on.exit(close(con), add = TRUE)
      get_namespace_fun("jsonlite", "fromJSON")(
        paste(readLines(con, warn = FALSE), collapse = "\n"),
        simplifyVector = FALSE
      )
    },
    error = function(e) {
      log_message(
        "Failed to retrieve tAge Zenodo metadata from {.url {url}}: {.val {e$message}}",
        message_type = "error"
      )
    }
  )
  metadata
}

select_tage_model_file <- function(
  files,
  mode,
  clock,
  model_species,
  model_tissue,
  preprocessing
) {
  expected <- normalize_tage_model_key(
    paste0(mode, "_", clock, "_", model_species, "_", model_tissue, "_", preprocessing, ".pkl")
  )
  keys <- vapply(files, function(file) file[["key"]] %||% "", character(1))
  normalized_keys <- normalize_tage_model_key(keys)
  idx <- which(normalized_keys == expected)
  if (length(idx) != 1L) {
    candidates <- keys[
      grepl(paste0("^", mode, "_", clock), keys, ignore.case = TRUE) &
        grepl(preprocessing_filename_token(preprocessing), keys, ignore.case = TRUE)
    ]
    log_message(
      paste0(
        "Could not find a unique tAge model for ",
        "{.val {mode}}, {.val {clock}}, {.val {model_species}}, ",
        "{.val {model_tissue}}, {.val {preprocessing}}. ",
        "Matching candidates: {.val {head(candidates, 12)}}"
      ),
      message_type = "error"
    )
  }
  files[[idx]]
}

normalize_tage_model_key <- function(x) {
  tolower(gsub("[^A-Za-z0-9]", "", x))
}

preprocessing_filename_token <- function(preprocessing) {
  switch(preprocessing,
    scaled_diff = "scaleddiff",
    yugene_diff = "yugenediff",
    preprocessing
  )
}

download_tage_model_file <- function(
  file_info,
  model_cache_dir,
  max_model_size = 1024^3,
  verbose = TRUE
) {
  key <- file_info[["key"]]
  size <- file_info[["size"]] %||% NA_real_
  checksum <- file_info[["checksum"]] %||% NA_character_
  url <- file_info[["links"]][["self"]]
  dest <- file.path(model_cache_dir, key)

  if (file.exists(dest) && validate_tage_cached_file(dest, checksum = checksum, size = size)) {
    log_message(
      "Use cached tAge model: {.file {dest}}",
      verbose = verbose
    )
    return(dest)
  }

  if (is.finite(max_model_size) && is.numeric(size) && !is.na(size) && size > max_model_size) {
    log_message(
      paste0(
        "tAge model {.file {key}} is {.val {round(size / 1024^3, 2)}} GB, ",
        "above {.arg max_model_size}. Set {.code max_model_size = Inf} ",
        "or pass a local {.arg model_paths} file to download it explicitly."
      ),
      message_type = "error"
    )
  }

  log_message(
    "Download tAge model {.file {key}} to {.file {model_cache_dir}}",
    verbose = verbose
  )
  tmp <- paste0(dest, ".tmp")
  status <- utils::download.file(
    url = url,
    destfile = tmp,
    mode = "wb",
    quiet = !isTRUE(verbose)
  )
  if (!identical(status, 0L) || !file.exists(tmp)) {
    unlink(tmp)
    log_message(
      "Failed to download tAge model from {.url {url}}",
      message_type = "error"
    )
  }
  if (!validate_tage_cached_file(tmp, checksum = checksum, size = size)) {
    unlink(tmp)
    log_message(
      "Downloaded tAge model failed size or checksum validation: {.file {key}}",
      message_type = "error"
    )
  }
  if (file.exists(dest)) {
    unlink(dest)
  }
  file.rename(tmp, dest)
  dest
}

validate_tage_cached_file <- function(path, checksum = NA_character_, size = NA_real_) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  if (!is.na(size) && file.info(path)$size != size) {
    return(FALSE)
  }
  if (!is.na(checksum) && grepl("^md5:", checksum)) {
    expected <- sub("^md5:", "", checksum)
    observed <- unname(tools::md5sum(path))
    if (!identical(tolower(observed), tolower(expected))) {
      return(FALSE)
    }
  }
  TRUE
}

infer_tage_model_name <- function(path) {
  name <- tolower(basename(path))
  if (grepl("scaled[_\\.-]*diff|scaleddiff", name)) {
    return("scaled_diff")
  }
  if (grepl("yugene[_\\.-]*diff|yugenediff", name)) {
    return("yugene_diff")
  }
  if (grepl("yugene", name)) {
    return("yugene")
  }
  if (grepl("scaled", name)) {
    return("scaled")
  }
  NA_character_
}

check_tage_python <- function(verbose = TRUE) {
  missing_modules <- vapply(
    c("joblib", "pandas", "sklearn"),
    function(pkg) !isTRUE(reticulate::py_module_available(pkg)),
    logical(1)
  )
  if (any(missing_modules)) {
    log_message(
      paste0(
        "Missing Python module{?s} for {.pkg tAge}: {.pkg {names(missing_modules)[missing_modules]}}. ",
        "Run {.code PrepareEnv(modules = 'tage')} or set {.envvar RETICULATE_PYTHON} ",
        "to an environment with {.pkg joblib}, {.pkg pandas}, and {.pkg scikit-learn}."
      ),
      message_type = "error",
      verbose = verbose
    )
  }
}

predict_tage_models <- function(
  processed,
  model_paths,
  species,
  mode,
  backend = c("python", "cpp")
) {
  backend <- match.arg(backend)
  if (identical(backend, "python")) {
    return(get_namespace_fun("tAge", "predict_tAge")(
      tAge_eset = processed,
      model_paths = model_paths,
      species = species,
      mode = mode
    ))
  }

  valid_names <- intersect(names(processed), names(model_paths))
  valid_names <- valid_names[valid_names %in% c("scaled", "scaled_diff", "yugene", "yugene_diff")]
  if (length(valid_names) == 0L) {
    log_message(
      "No overlapping tAge processed data and model names were found",
      message_type = "error"
    )
  }

  result <- NULL
  for (name in valid_names) {
    model <- readRDS(model_paths[[name]])
    pred <- predict_tage_cpp_one(
      eset = processed[[name]],
      model = model,
      species = species,
      prefix = paste0(name, "_", mode, "_")
    )
    if (is.null(result)) {
      result <- pred
    } else {
      new_cols <- setdiff(colnames(pred), colnames(result))
      result <- cbind(result, pred[, new_cols, drop = FALSE])
    }
  }
  result
}

predict_tage_cpp_one <- function(eset, model, species, prefix = "EN_") {
  if (!identical(model$model_type, "ElasticNet")) {
    log_message(
      "Unsupported converted tAge model type: {.val {model$model_type}}",
      message_type = "error"
    )
  }
  expr <- Biobase::exprs(eset)
  features <- as.character(model$feature_names)
  selected <- model$selected
  if (is.logical(selected)) {
    selected_idx <- which(selected)
  } else {
    selected_idx <- as.integer(selected)
  }
  features_use <- features[selected_idx]
  feature_match <- match(features_use, rownames(expr))
  prediction <- tage_elastic_net_predict_cpp(
    expr = expr,
    feature_match = feature_match,
    imputer = model$imputer_statistics[selected_idx],
    center = model$center[selected_idx],
    scale = model$scale %||% numeric(),
    coef = model$coef,
    intercept = model$intercept
  )

  adjustment <- model$species_adjustment[[species]]
  if (!is.null(adjustment) && !is.na(adjustment)) {
    prediction <- prediction * adjustment
  }

  annotation <- Biobase::pData(eset)
  annotation[[paste0(prefix, "tAge")]] <- prediction
  annotation
}

make_tage_eset <- function(
  object,
  metadata = NULL,
  group.by = NULL,
  split.by = NULL,
  assay = NULL,
  layer = "counts",
  coverage_threshold = 1e6,
  control_group_column = NULL,
  shuffle = FALSE,
  seed = NULL,
  verbose = TRUE
) {
  if (inherits(object, "ExpressionSet")) {
    return(object)
  }

  if (inherits(object, "Seurat")) {
    assay <- assay %||% SeuratObject::DefaultAssay(object)
    aggregate_by <- unique(c(group.by, split.by, control_group_column))
    aggregate_by <- aggregate_by[!is.na(aggregate_by) & nzchar(aggregate_by)]
    if (length(aggregate_by) == 0L) {
      return(get_namespace_fun("tAge", "aggregate_pseudobulk")(
        seurat_obj = object,
        coverage_threshold = coverage_threshold,
        assay = assay,
        layer = layer,
        shuffle = shuffle,
        seed = seed,
        verbose = verbose
      ))
    }
    missing_cols <- setdiff(aggregate_by, colnames(object@meta.data))
    if (length(missing_cols) > 0L) {
      log_message(
        "Missing Seurat metadata column{?s}: {.field {missing_cols}}",
        message_type = "error"
      )
    }
    return(get_namespace_fun("tAge", "aggregate_on_obs_columns")(
      seurat_obj = object,
      obs_column_names = aggregate_by,
      coverage_threshold = coverage_threshold,
      assay = assay,
      layer = layer,
      shuffle = shuffle,
      seed = seed,
      verbose = verbose
    ))
  }

  expr <- object
  if (!inherits(expr, c("matrix", "Matrix", "data.frame"))) {
    log_message(
      "{.arg object} must be a {.cls Seurat}, {.cls ExpressionSet}, matrix, or data.frame",
      message_type = "error"
    )
  }
  expr <- as.matrix(expr)
  if (is.null(rownames(expr)) || is.null(colnames(expr))) {
    log_message(
      "Matrix input must have gene row names and sample column names",
      message_type = "error"
    )
  }
  if (is.null(metadata)) {
    metadata <- data.frame(
      sample_id = colnames(expr),
      row.names = colnames(expr),
      stringsAsFactors = FALSE
    )
  } else {
    metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
    if (is.null(rownames(metadata))) {
      log_message(
        "{.arg metadata} must have row names matching matrix columns",
        message_type = "error"
      )
    }
    missing_samples <- setdiff(colnames(expr), rownames(metadata))
    if (length(missing_samples) > 0L) {
      log_message(
        "{.arg metadata} is missing sample{?s}: {.val {missing_samples}}",
        message_type = "error"
      )
    }
    metadata <- metadata[colnames(expr), , drop = FALSE]
  }
  get_namespace_fun("tAge", "make_ExpressionSet")(
    exprs_data = expr,
    phenodata = metadata,
    verbose = verbose
  )
}

run_tage_by_group_scop <- function(
  eset,
  split_by,
  model_paths,
  species = "mouse",
  mode = "EN",
  backend = c("python", "cpp"),
  gene_mapping_type = "Gene.Symbol",
  control_group_column = NULL,
  control_group_label = NULL,
  count_threshold = 10,
  percent_threshold = 20,
  min_samples = 5,
  verbose = TRUE
) {
  backend <- match.arg(backend)
  groups <- unique(Biobase::pData(eset)[[split_by]])
  groups <- groups[!is.na(groups)]
  results_list <- list()

  for (group in groups) {
    mask <- Biobase::pData(eset)[[split_by]] == group
    eset_sub <- eset[, mask]
    if (ncol(eset_sub) < min_samples) {
      log_message(
        "Skip {.val {group}}: {.val {ncol(eset_sub)}} samples < {.arg min_samples}",
        verbose = verbose
      )
      next
    }

    log_message(
      "Run {.pkg tAge} for {.field {split_by}} = {.val {group}}",
      verbose = verbose
    )
    result <- tryCatch(
      {
        processed <- get_namespace_fun("tAge", "tAge_preprocessing")(
          eset = eset_sub,
          species = species,
          gene_mapping_type = gene_mapping_type,
          verbose = verbose,
          control_group_column = control_group_column,
          control_group_label = control_group_label,
          count_threshold = count_threshold,
          percent_threshold = percent_threshold
        )
        prediction <- predict_tage_models(
          processed = processed,
          model_paths = model_paths,
          species = species,
          mode = mode,
          backend = backend
        )
        prediction[[split_by]] <- group
        prediction
      },
      error = function(e) {
        log_message(
          "Skip {.val {group}} after tAge error: {.val {e$message}}",
          message_type = "warning",
          verbose = verbose
        )
        NULL
      }
    )
    if (!is.null(result)) {
      results_list[[as.character(group)]] <- result
    }
  }

  if (length(results_list) == 0L) {
    log_message(
      "No {.arg split.by} group produced valid tAge predictions",
      message_type = "error"
    )
  }

  all_cols <- unique(unlist(lapply(results_list, colnames)))
  results_aligned <- lapply(results_list, function(result) {
    missing_cols <- setdiff(all_cols, colnames(result))
    for (col in missing_cols) {
      result[[col]] <- NA
    }
    result[, all_cols, drop = FALSE]
  })
  do.call(rbind, results_aligned)
}
