grn_hash <- function(value = NULL, file = NULL) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("Install digest before using GRN cache identities.", call. = FALSE)
  if (!is.null(file)) return(digest::digest(file = file, algo = "sha256"))
  digest::digest(serialize(value, NULL, xdr = FALSE, version = 3L), algo = "sha256", serialize = FALSE)
}

grn_resource_identity <- function(resource) {
  if (is.null(resource)) return(NULL)
  type <- class(resource)
  is_genome <- isS4(resource) && !is.null(attr(resource, "pkgname")) && !is.null(attr(resource, "seqinfo"))
  exact_genome <- identical(type, structure("BSgenome", package = "BSgenome"))
  if (is_genome && (exact_genome || inherits(resource, "BSgenome"))) {
    package <- attr(resource, "pkgname")
    return(list(class = type, package = package,
                version = as.character(utils::packageVersion(package)),
                seqinfo = if (exact_genome) attr(resource, "seqinfo") else GenomeInfoDb::seqinfo(resource)))
  }
  if (is.character(resource) && length(resource) && all(file.exists(resource))) {
    paths <- normalizePath(resource, winslash = "/", mustWork = TRUE)
    return(lapply(paths, function(path) {
      files <- if (dir.exists(path)) sort(list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)) else path
      files <- files[!dir.exists(files)]
      list(path = path, files = stats::setNames(vapply(files, function(file) grn_hash(file = file), ""), if (dir.exists(path)) substring(files, nchar(path) + 2L) else basename(files)))
    }))
  }
  if (is.list(resource)) return(lapply(resource, grn_resource_identity))
  list(class = type, sha256 = grn_hash(resource))
}

grn_hash_reference_files <- function(paths, config, algorithm) {
  script <- system.file("grn", "hash_references.py", package = "scop")
  if (!nzchar(script)) stop("The installed GRN resource checker is missing; reinstall scop.", call. = FALSE)
  process <- grn_run_process(config$runtime$python, c(script, algorithm, unname(paths)), config)
  hashes <- jsonlite::read_json(process$stdout, simplifyVector = TRUE)
  if (!is.character(hashes) || length(hashes) != length(paths) || anyNA(hashes) ||
      any(!grepl("^[0-9a-f]{64}$", hashes))) {
    stop("GRN resource verification returned invalid fingerprints.", call. = FALSE)
  }
  hashes
}

grn_scenicplus_resource_identity <- function(config) {
  resources <- config$resources
  paths <- vapply(resources[c("ranking_db", "scores_db")], normalizePath, "", winslash = "/", mustWork = TRUE)
  large <- file.info(paths)$size >= 64 * 1024^2
  if (!any(large)) return(grn_resource_identity(resources))
  paths <- paths[large]
  algorithm <- "sha256-tree-64m-v1"
  hashes <- grn_hash_reference_files(paths, config, algorithm)
  result <- lapply(resources[setdiff(names(resources), names(paths))], grn_resource_identity)
  for (i in seq_along(paths)) {
    result[[names(paths)[[i]]]] <- list(list(path = paths[[i]], algorithm = algorithm,
      files = stats::setNames(hashes[[i]], basename(paths[[i]]))))
  }
  result[names(resources)]
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
  list(expression = grn_hash(input$expression), rna = grn_hash(rna), chromatin = grn_hash(chromatin), metadata = grn_hash(metadata),
       cells = input$cells, groups = input$groups %||% config$groups,
       regulators = input$regulators, targets = input$targets,
       pseudotime = input$pseudotime, input_mode = input$input_mode %||% config$input_mode)
}

grn_package_code_identity <- function(package = find.package("scop")) {
  directories <- file.path(package, c("R", "libs", "grn", "python"))
  files <- sort(unique(c(file.path(package, "DESCRIPTION"), unlist(lapply(directories[dir.exists(directories)], list.files,
    recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE), use.names = FALSE))))
  files <- files[file.exists(files) & !dir.exists(files)]
  relative <- substring(files, nchar(package) + 2L)
  files <- files[!grepl("(^|/)(__pycache__|[.]DS_Store)(/|$)|[.](pyc|pyo)$", relative)]
  stats::setNames(vapply(files, function(path) grn_hash(file = path), ""), substring(files, nchar(package) + 2L))
}

grn_external_sources <- function(method, resources) {
  files <- switch(method,
    gniplr = c(entry = "code/methods/GNIPLR.py"),
    pidc = c(entry = "code/vendor/PIDC/PIDC.R"),
    scode = c(entry = "code/vendor/SCODE/SCODE.R"),
    sincerities = c(entry = "code/methods/SINCERITIES.R"),
    scsgl = c(
      entry = "code/vendor/SCSGL/run_scSGL.py",
      paste0("code/vendor/SCSGL/pysrc/", c(
        "__init__.py", "associations/__init__.py", "associations/correlation.py",
        "graphlearning/__init__.py", "graphlearning/signed.py", "graphlearning/unsigned.py",
        "utils/__init__.py", "utils/rowsum_matrix.py"
      ))
    ),
    NULL
  )
  if (is.null(files)) return(NULL)
  root <- resources$source_dir
  if (!is.character(root) || length(root) != 1L || is.na(root) || !dir.exists(root)) {
    stop(method, " requires resources$source_dir pointing to a multiCSN_manuscript checkout. See ?RunGRN.", call. = FALSE)
  }
  paths <- file.path(root, files)
  missing <- !file.exists(paths) | dir.exists(paths)
  if (any(missing)) {
    stop("Missing ", method, " upstream source files: ", paste(files[missing], collapse = ", "), call. = FALSE)
  }
  stats::setNames(normalizePath(paths, winslash = "/", mustWork = TRUE), names(files))
}
