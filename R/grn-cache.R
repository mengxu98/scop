grn_cached_run <- function(input, config, compute) {
  if (!is.function(compute)) stop("GRN cache compute must be a function.", call. = FALSE)
  implementation_functions <- c(config$runner, config$prepare)
  implementations <- lapply(implementation_functions, function(name) {
    fun <- get(name, mode = "function", inherits = TRUE)
    list(name = name, formals = formals(fun), body = body(fun))
  })
  identity <- list(schema_version = "1.0", method = config$method, backend = config$backend,
                   implementation = config$implementation, implementation_functions = implementations,
                   package_versions = config$package_versions,
                   package_code = grn_package_code_identity(),
                   parameters = config[c("method_args", "assay", "layer", "atac_assay", "atac_layer", "group.by", "cores", "seed", "runtime_isolation")],
                   input = grn_input_identity(input, config), runtime = config$runtime,
                   resources = config$resource_identity %||% grn_resource_identity(config$resources))
  key <- grn_hash(identity)
  directory <- file.path(config$work_dir %||% tempdir(), "cache", config$method, key)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  directory <- normalizePath(directory, winslash = "/", mustWork = TRUE)
  manifest_file <- file.path(directory, "completed.rds")
  read_completed <- function() {
    if (!file.exists(manifest_file)) return(NULL)
    manifest <- tryCatch(readRDS(manifest_file), error = function(e) NULL)
    if (!is.list(manifest) || !identical(manifest$status, "completed") || !identical(manifest$key, key) ||
        !is.character(manifest$result_file) || length(manifest$result_file) != 1L ||
        !grepl("^attempt-[A-Za-z0-9]+/result[.]rds$", manifest$result_file)) return(NULL)
    result_file <- file.path(directory, manifest$result_file)
    if (!file.exists(result_file) ||
        !identical(manifest$result_sha256, grn_hash(file = result_file))) return(NULL)
    value <- tryCatch(readRDS(result_file), error = function(e) NULL)
    if (!inherits(value, "GRNResult") || !identical(value$schema_version, "1.0")) return(NULL)
    artifacts <- manifest$artifacts
    if (length(artifacts) && !all(vapply(artifacts, function(artifact) {
      is.list(artifact) && !is.null(artifact$path) && length(artifact$path) == 1L && file.exists(artifact$path) && !dir.exists(artifact$path) &&
        isTRUE(tryCatch(identical(artifact$sha256, grn_hash(file = artifact$path)), error = function(e) FALSE))
    }, logical(1)))) return(NULL)
    value$execution$cache <- list(key = key, hit = TRUE, directory = directory)
    value
  }
  if (!isTRUE(config$force)) {
    value <- read_completed()
    if (!is.null(value)) return(value)
  }
  lock <- file.path(directory, "running.lock")
  if (!dir.create(lock, showWarnings = FALSE)) {
    stop("Another GRN run owns this cache key. Lock: ", lock, call. = FALSE)
  }
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  if (!isTRUE(config$force)) {
    value <- read_completed()
    if (!is.null(value)) return(value)
  }
  stage <- tempfile("attempt-", tmpdir = directory)
  dir.create(stage)
  status_file <- file.path(stage, "status.rds")
  status <- list(status = "running", key = key, pid = Sys.getpid(), started_at = Sys.time(), identity = identity)
  saveRDS(status, status_file)
  committed <- FALSE
  on.exit({
    if (!committed && !status$status %in% c("failed", "cancelled")) {
      status$status <- "failed"
      status$finished_at <- Sys.time()
      saveRDS(status, status_file)
    }
  }, add = TRUE)
  value <- tryCatch(compute(), interrupt = function(e) {
    status$status <<- "cancelled"
    status$finished_at <<- Sys.time()
    saveRDS(status, status_file)
    stop(e)
  }, error = function(e) {
    status$status <<- "failed"
    status$error <<- conditionMessage(e)
    status$finished_at <<- Sys.time()
    saveRDS(status, status_file)
    stop(e)
  })
  if (!inherits(value, "GRNResult") || !identical(value$schema_version, "1.0") ||
      !identical(value$execution$status, "completed")) {
    status$status <- "failed"
    status$error <- "The method did not return a completed GRNResult."
    saveRDS(status, status_file)
    stop(status$error, call. = FALSE)
  }
  value$execution$cache <- list(key = key, hit = FALSE, directory = directory)
  staged_result <- file.path(stage, "result.rds")
  saveRDS(value, staged_result)
  status$status <- "completed"
  status$finished_at <- Sys.time()
  status$result_file <- paste0(basename(stage), "/result.rds")
  status$result_sha256 <- grn_hash(file = staged_result)
  artifacts <- value$execution$artifacts %||% list()
  status$artifacts <- unlist(lapply(artifacts, function(artifact) {
    if (is.character(artifact) && length(artifact) == 1L) artifact <- list(path = artifact)
    if (is.null(artifact$path) || length(artifact$path) != 1L || !file.exists(artifact$path)) stop("GRN artifacts must identify existing files or directories.", call. = FALSE)
    path <- normalizePath(artifact$path, winslash = "/", mustWork = TRUE)
    files <- if (dir.exists(path)) list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE) else path
    files <- files[!dir.exists(files)]
    files <- normalizePath(files, winslash = "/", mustWork = TRUE)
    relative <- if (dir.exists(path)) substring(files, nchar(path) + 2L) else basename(files)
    files <- files[!startsWith(files, paste0(directory, "/")) &
      !grepl("(^|/)(cache|[.]tmp)(/|$)|[.]tmp($|/)", relative)]
    if (!length(files)) stop("GRN artifact directory contains no completed files: ", path, call. = FALSE)
    lapply(sort(files), function(file) list(path = file, sha256 = grn_hash(file = file)))
  }), recursive = FALSE)
  if (length(status$artifacts)) status$artifacts <- status$artifacts[!duplicated(vapply(status$artifacts, `[[`, "", "path"))]
  staged_manifest <- file.path(stage, "completed.rds")
  saveRDS(status, staged_manifest)
  if (.Platform$OS.type == "windows" && file.exists(manifest_file)) unlink(manifest_file)
  if (!file.rename(staged_manifest, manifest_file)) {
    stop("Unable to commit the completed GRN cache. Attempt: ", stage, call. = FALSE)
  }
  committed <- TRUE
  saveRDS(status, status_file)
  value
}
