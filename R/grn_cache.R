grn_cached_run <- function(input, config, compute) {
  if (!is.function(compute)) stop("GRN cache compute must be a function.", call. = FALSE)
  identity <- list(
    schema_version = "1.0", method = config$method, backend = config$backend,
    implementation = config$implementation,
    package_versions = config$package_versions,
    package_code = grn_package_code_identity(),
    parameters = config[c("method_args", "assay", "layer", "atac_assay", "atac_layer", "group.by", "cores", "seed", "runtime_isolation")],
    input = grn_input_identity(input, config), runtime = config$runtime,
    resources = config$resource_identity %||% grn_resource_identity(config$resources)
  )
  key <- grn_hash(identity)
  directory <- file.path(config$work_dir %||% tempdir(), "cache", config$method, key)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  directory <- normalizePath(directory, winslash = "/", mustWork = TRUE)
  manifest_file <- file.path(directory, "completed.rds")
  read_completed <- function() {
    if (!file.exists(manifest_file)) {
      return(NULL)
    }
    manifest <- tryCatch(readRDS(manifest_file), error = function(e) NULL)
    if (!is.list(manifest) || !identical(manifest$status, "completed") || !identical(manifest$key, key) ||
      !is.character(manifest$result_file) || length(manifest$result_file) != 1L ||
      !grepl("^attempt-[A-Za-z0-9]+/result[.]rds$", manifest$result_file)) {
      return(NULL)
    }
    result_file <- file.path(directory, manifest$result_file)
    if (!file.exists(result_file) ||
      !identical(manifest$result_sha256, grn_hash(file = result_file))) {
      return(NULL)
    }
    value <- tryCatch(readRDS(result_file), error = function(e) NULL)
    if (!inherits(value, "GRNResult") || !identical(value$schema_version, "1.0")) {
      return(NULL)
    }
    artifacts <- manifest$artifacts
    if (length(artifacts) && !all(vapply(artifacts, function(artifact) {
      is.list(artifact) && !is.null(artifact$path) && length(artifact$path) == 1L && file.exists(artifact$path) && !dir.exists(artifact$path) &&
        isTRUE(tryCatch(identical(artifact$sha256, grn_hash(file = artifact$path)), error = function(e) FALSE))
    }, logical(1)))) {
      return(NULL)
    }
    value$execution$cache <- list(key = key, hit = TRUE, directory = directory)
    value
  }
  if (!isTRUE(config$force)) {
    value <- read_completed()
    if (!is.null(value)) {
      return(value)
    }
  }
  lock <- file.path(directory, "running.lock")
  if (!dir.create(lock, showWarnings = FALSE)) {
    stop("Another GRN run owns this cache key. Lock: ", lock, call. = FALSE)
  }
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  if (!isTRUE(config$force)) {
    value <- read_completed()
    if (!is.null(value)) {
      return(value)
    }
  }
  stage <- tempfile("attempt-", tmpdir = directory)
  dir.create(stage)
  status_file <- file.path(stage, "status.rds")
  status <- list(status = "running", key = key, pid = Sys.getpid(), started_at = Sys.time(), identity = identity)
  saveRDS(status, status_file)
  committed <- FALSE
  on.exit(
    {
      if (!committed && !status$status %in% c("failed", "cancelled")) {
        status$status <- "failed"
        status$finished_at <- Sys.time()
        saveRDS(status, status_file)
      }
    },
    add = TRUE
  )
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

grn_hash <- function(value = NULL, file = NULL) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("Install digest before using GRN cache identities.", call. = FALSE)
  if (!is.null(file)) {
    return(digest::digest(file = file, algo = "sha256"))
  }
  digest::digest(serialize(value, NULL, xdr = FALSE, version = 3L), algo = "sha256", serialize = FALSE)
}

grn_resource_identity <- function(resource) {
  if (is.null(resource)) {
    return(NULL)
  }
  type <- class(resource)
  is_genome <- isS4(resource) && !is.null(attr(resource, "pkgname")) && !is.null(attr(resource, "seqinfo"))
  exact_genome <- identical(type, structure("BSgenome", package = "BSgenome"))
  if (is_genome && (exact_genome || inherits(resource, "BSgenome"))) {
    package <- attr(resource, "pkgname")
    return(list(
      class = type, package = package,
      version = as.character(utils::packageVersion(package)),
      seqinfo = if (exact_genome) attr(resource, "seqinfo") else GenomeInfoDb::seqinfo(resource)
    ))
  }
  if (is.character(resource) && length(resource) && all(file.exists(resource))) {
    paths <- normalizePath(resource, winslash = "/", mustWork = TRUE)
    return(lapply(paths, function(path) {
      files <- if (dir.exists(path)) sort(list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)) else path
      files <- files[!dir.exists(files)]
      list(path = path, files = stats::setNames(vapply(files, function(file) grn_hash(file = file), ""), if (dir.exists(path)) substring(files, nchar(path) + 2L) else basename(files)))
    }))
  }
  if (is.list(resource)) {
    return(lapply(resource, grn_resource_identity))
  }
  list(class = type, sha256 = grn_hash(resource))
}

grn_input_identity <- function(input, config) {
  object <- input$object
  chromatin <- NULL
  rna <- NULL
  metadata <- NULL
  if (!is.null(object) && identical(input$input_mode %||% config$input_mode, "multiome")) {
    assay <- config$atac_assay
    assay_layers <- function(assay, requested) {
      available <- SeuratObject::Layers(object[[assay]], search = NA)
      layers <- unique(c(intersect(requested, available), grep("^(counts|data)($|[.])", available, value = TRUE)))
      stats::setNames(lapply(layers, function(layer) SeuratObject::LayerData(object, assay = assay, layer = layer)[, input$cells, drop = FALSE]), layers)
    }
    chromatin <- list(layers = assay_layers(assay, config$atac_layer %||% "data"))
    rna <- assay_layers(config$assay, config$layer %||% "data")
    metadata <- object[[]][input$cells, , drop = FALSE]
    state <- object@misc[["CSNObject"]]
    chromatin$candidates <- if (is.list(state)) state[intersect(names(state), c("params", "attributes", "regions", "tfs", "celltypes"))] else state
    chromatin$input_annotation <- object@misc$multiCSN_annotations[[assay]] %||% object@misc$multiCSN_input$annotation
    if (inherits(object[[assay]], "ChromatinAssay")) {
      chromatin$annotation <- Signac::Annotation(object[[assay]])
      chromatin$ranges <- Signac::granges(object[[assay]])
      chromatin$motifs <- methods::slot(object[[assay]], "motifs")
    }
  }
  list(
    expression = grn_hash(input$expression), rna = grn_hash(rna), chromatin = grn_hash(chromatin), metadata = grn_hash(metadata),
    cells = input$cells, groups = input$groups %||% config$groups,
    regulators = input$regulators, targets = input$targets,
    pseudotime = input$pseudotime, input_mode = input$input_mode %||% config$input_mode
  )
}

grn_package_code_identity <- function(package = find.package("scop")) {
  directories <- file.path(package, c("R", "libs", "grn", "python"))
  files <- sort(unique(c(file.path(package, "DESCRIPTION"), unlist(lapply(directories[dir.exists(directories)], list.files,
    recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE
  ), use.names = FALSE))))
  files <- files[file.exists(files) & !dir.exists(files)]
  relative <- substring(files, nchar(package) + 2L)
  files <- files[!grepl("(^|/)(__pycache__|[.]DS_Store)(/|$)|[.](pyc|pyo)$", relative)]
  stats::setNames(vapply(files, function(path) grn_hash(file = path), ""), substring(files, nchar(package) + 2L))
}
