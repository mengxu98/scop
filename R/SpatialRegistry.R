# Spatial method discovery and result inspection --------------------------------

.spatial_backend_cache <- new.env(parent = emptyenv())

spatial_registry_entry <- function(
  method,
  kind,
  task,
  implementation_files,
  tool_key = NA_character_,
  backend_id = "core",
  status = "stable",
  coordinate_space_current = "none",
  coordinate_space_target = coordinate_space_current,
  coordinate_requirement = "none",
  scalability = "standard",
  plot_function = NA_character_,
  documentation = method,
  backend_requirement = "all"
) {
  data.frame(
    method = method,
    kind = kind,
    task = task,
    implementation_files = implementation_files,
    tool_key = tool_key,
    backend_id = backend_id,
    status = status,
    coordinate_space_current = coordinate_space_current,
    coordinate_space_target = coordinate_space_target,
    coordinate_requirement = coordinate_requirement,
    scalability = scalability,
    plot_function = plot_function,
    documentation = documentation,
    backend_requirement = backend_requirement,
    stringsAsFactors = FALSE
  )
}

spatial_method_registry <- function() {
  entry <- spatial_registry_entry
  rows <- list(
    entry("RunSpotQC", "analysis", "quality_control", "RunSpotQC.R", backend_id = "core", coordinate_space_current = "none", coordinate_space_target = "none", coordinate_requirement = "identity_only", plot_function = "SpatialSpotPlot"),
    entry("RunSpaNorm", "analysis", "normalization", "RunSpaNorm.R", "SpaNorm", "spanorm", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "SpatialSpotPlot"),
    entry("RunSpatialQM", "analysis", "quality_control", "RunSpatialQM.R", "SpatialQM", "spatialqm", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "SpatialSpotPlot"),
    entry("RunSpotSweeper", "analysis", "quality_control", "RunSpotSweeper.R", "SpotSweeper", "spotsweeper", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "SpatialSpotPlot"),

    entry("srt_to_giotto", "bridge", "framework_bridge", "SpatialFrameworkConvert.R", backend_id = "giotto;giotto_class", coordinate_space_current = "raw", coordinate_requirement = "backend_managed"),
    entry("giotto_to_srt", "bridge", "framework_bridge", "SpatialFrameworkConvert.R", backend_id = "giotto_class", coordinate_space_current = "raw", coordinate_requirement = "backend_managed"),
    entry("srt_to_spata2", "bridge", "framework_bridge", "SpatialFrameworkConvert.R", backend_id = "spata2", coordinate_space_current = "raw", coordinate_requirement = "backend_managed"),
    entry("spata2_to_srt", "bridge", "framework_bridge", "SpatialFrameworkConvert.R", backend_id = "spata2", coordinate_space_current = "raw", coordinate_requirement = "backend_managed"),
    entry("SeuratToScopGiotto", "bridge", "framework_bridge", "GiottoObject.R", "Giotto", "giotto;giotto_class", "legacy", "legacy_display", "legacy_display", "backend_managed"),
    entry("RunGiottoWorkflow", "workflow", "framework_workflow", "GiottoObject.R", "Giotto", "giotto;giotto_class", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("GiottoPreprocess", "analysis", "framework_workflow", "GiottoObject.R", "Giotto", "giotto", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("GiottoReduce", "analysis", "framework_workflow", "GiottoObject.R", "Giotto", "giotto", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("GiottoCluster", "analysis", "framework_workflow", "GiottoObject.R", "Giotto", "giotto", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("GiottoSpatialNetwork", "analysis", "framework_workflow", "GiottoObject.R", "Giotto", "giotto", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("GiottoSpatialGenes", "analysis", "framework_workflow", "GiottoObject.R", "Giotto", "giotto", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("GiottoSpatialModules", "analysis", "framework_workflow", "GiottoObject.R", "Giotto", "giotto", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("GiottoCellProximity", "analysis", "framework_workflow", "GiottoObject.R", "Giotto", "giotto", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("GiottoHMRF", "analysis", "framework_workflow", "GiottoObject.R", "Giotto", "giotto", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("AddGiottoToSeurat", "bridge", "framework_bridge", "GiottoObject.R", "Giotto", "giotto;giotto_class", "legacy", "legacy_display", "legacy_display", "backend_managed"),
    entry("RunGiottoCluster", "analysis", "framework_workflow", "RunGiottoCluster.R", "GiottoCluster", "giotto;giotto_class", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("RunGiottoCellProximity", "analysis", "framework_workflow", "RunGiottoSpatialMethods.R", "GiottoCellProximity", "giotto;giotto_class", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("RunGiottoSpatialGenes", "analysis", "framework_workflow", "RunGiottoSpatialMethods.R", "GiottoSpatialGenes", "giotto;giotto_class", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),
    entry("RunGiottoSpatialModules", "analysis", "framework_workflow", "RunGiottoSpatialMethods.R", "GiottoSpatialModules", "giotto;giotto_class", "legacy", "legacy_display", "legacy_display", "backend_managed", plot_function = "GiottoPlot"),

    entry("RunDeconvolution", "analysis", "deconvolution", "RunDeconvolution.R", backend_id = "music;bisquerna;bayesprism;cibersort", coordinate_space_current = "none", coordinate_space_target = "none", coordinate_requirement = "identity_only", plot_function = "DeconvolutionPlot", backend_requirement = "any"),
    entry("RunRCTD", "analysis", "deconvolution", "RunRCTD.R", "RCTD", "spacexr", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "DeconvolutionPlot"),
    entry("RunCSIDE", "analysis", "deconvolution", "RunCSIDE.R", "CSIDE", "spacexr", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "DeconvolutionPlot"),
    entry("RunCARD", "analysis", "deconvolution", "RunCARD.R", "CARD", "card", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "DeconvolutionPlot"),
    entry("RunSTdeconvolve", "analysis", "deconvolution", "RunSTdeconvolve.R", "STdeconvolve", "stdeconvolve", coordinate_space_current = "none", coordinate_space_target = "none", coordinate_requirement = "identity_only", plot_function = "STdeconvolvePlot"),
    entry("RunSPOTlight", "analysis", "deconvolution", "RunSPOTlight.R", "SPOTlight", "spotlight", coordinate_space_current = "none", coordinate_space_target = "none", coordinate_requirement = "identity_only", plot_function = "DeconvolutionPlot"),
    entry("RunSpatialDWLS", "analysis", "deconvolution", "RunSpatialDWLS.R", "SpatialDWLS", "core", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "DeconvolutionPlot"),
    entry("RunSpatialEcoTyper", "analysis", "ecotype", "RunSpatialEcoTyper.R", "SpatialEcoTyper", "spatialecotyper", coordinate_space_current = "none", coordinate_space_target = "none", coordinate_requirement = "identity_only", plot_function = "SpatialEcoTyperSpatialPlot"),

    entry("RunBayesSpace", "analysis", "domain", "RunBayesSpace.R", "BayesSpace", "bayesspace", coordinate_space_current = "raw", coordinate_requirement = "distance_sensitive", plot_function = "SpatialSpotPlot"),
    entry("RunBANKSY", "analysis", "domain", "RunBANKSY.R", "BANKSY", "banksy", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "SpatialSpotPlot"),
    entry("RunCytoSPACE", "analysis", "mapping", "RunCytoSPACE.R", "CytoSPACE", "core", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "SpatialSpotPlot"),
    entry("RunSmoothClust", "analysis", "domain", "RunSmoothClust.R", "SmoothClust", "smoothclust", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", scalability = "sparse_required", plot_function = "SpatialSpotPlot"),
    entry("RunMERINGUE", "analysis", "feature_pattern", "RunMERINGUE.R", "MERINGUE", "meringue", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "SpatialVariableFeaturePlot"),
    entry("RunSpatialVariableFeatures", "analysis", "feature_pattern", "RunSpatialVariableFeatures.R", "SpatialVariableFeatures", "core;sparkx;nnsvg", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", scalability = "sparse_required", plot_function = "SpatialVariableFeaturePlot", backend_requirement = "any"),
    entry("RunSpatialGradientFeatures", "analysis", "feature_pattern", "RunSpatialGradientFeatures.R", "SpatialGradientFeatures", "core;spata2", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "SpatialGradientPlot", backend_requirement = "any"),

    entry("RunSpatialNetwork", "analysis", "network", "RunSpatialNetwork.R", "SpatialNetwork", "biocneighbors", coordinate_space_current = "raw", coordinate_requirement = "distance_sensitive", scalability = "sparse", plot_function = "SpatialNetworkPlot"),
    entry("RunSpatialNeighborhood", "analysis", "neighborhood", "RunSpatialNeighborhood.R", "SpatialNeighborhood", "spicyr", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", scalability = "sparse_required", plot_function = "SpatialNeighborhoodPlot"),
    entry("RunStatialKontextual", "analysis", "neighborhood", "RunStatialKontextual.R", "StatialKontextual", "statial", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive"),
    entry("RunSpatialIntegration", "analysis", "integration", "RunSpatialIntegration.R", "SpatialIntegration", "precast;bass;spatialmnn", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive", plot_function = "SpatialIntegrationPlot", backend_requirement = "any"),
    entry("RunMistyR", "analysis", "neighborhood", "RunMistyR.R", "MistyR", "mistyr", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive"),
    entry("RunSemlaSpatialNetwork", "analysis", "neighborhood", "RunSemla.R", "SemlaSpatialNetwork", "semla", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive"),
    entry("RunSemlaLocalG", "analysis", "neighborhood", "RunSemla.R", "SemlaLocalG", "semla", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive"),
    entry("RunSemlaRadialDistance", "analysis", "neighborhood", "RunSemla.R", "SemlaRadialDistance", "semla", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive"),
    entry("RunSemlaRegionNeighbors", "analysis", "neighborhood", "RunSemla.R", "SemlaRegionNeighbors", "semla", coordinate_space_current = "legacy_display", coordinate_space_target = "raw", coordinate_requirement = "distance_sensitive"),

    entry("GiottoPlot", "plot", "visualization", "GiottoPlot.R", backend_id = "giotto;giotto_class", status = "legacy", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("DeconvolutionPlot", "plot", "visualization", "DeconvolutionPlot.R", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("SpatialEcoTyperCompositionPlot", "plot", "visualization", "RunSpatialEcoTyper.R", coordinate_space_current = "none"),
    entry("SpatialEcoTyperSpatialPlot", "plot", "visualization", "RunSpatialEcoTyper.R", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("SpatialGradientPlot", "plot", "visualization", "RunSpatialGradientFeatures.R", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("SpatialIntegrationPlot", "plot", "visualization", "RunSpatialIntegration.R", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("SpatialNeighborhoodPlot", "plot", "visualization", "RunSpatialNeighborhood.R", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("SpatialNetworkPlot", "plot", "visualization", "RunSpatialNetwork.R", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("SpatialCellPlot", "plot", "visualization", "SpatialCellPlot.R", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("SpatialSpotPlot", "plot", "visualization", "SpatialSpotPlot.R", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("SpatialVariableFeaturePlot", "plot", "visualization", "RunSpatialVariableFeatures.R", coordinate_space_current = "display", coordinate_requirement = "display_only"),
    entry("STdeconvolvePlot", "plot", "visualization", "RunSTdeconvolve.R", backend_id = "stdeconvolve", coordinate_space_current = "none"),
    entry("standard_scop", "workflow", "recommended_workflow", "standard_scop.R", backend_id = "core", coordinate_space_current = "mixed", coordinate_space_target = "mixed", coordinate_requirement = "backend_managed")
  )
  registry <- do.call(rbind, rows)
  rownames(registry) <- NULL
  registry
}

spatial_giotto_symbol_registry <- function() {
  class_symbols <- c(
    "createExprObj", "createGiottoObject", "createNearestNetwork",
    "createSpatialNetwork", "fDataDT", "getDimReduction",
    "getSpatialNetwork", "pDataDT", "setExpression"
  )
  giotto_symbols <- c(
    "normalizeGiotto", "runPCA", "runUMAP", "doLeidenCluster",
    "doLouvainCluster", "binSpect", "detectSpatialCorFeats",
    "clusterSpatialCorFeats", "cellProximityEnrichment", "doHMRF",
    "addHMRF", "addStatistics", "calculateHVF"
  )
  data.frame(
    symbol = c(class_symbols, giotto_symbols),
    package = c(
      rep("GiottoClass", length(class_symbols)),
      rep("Giotto", length(giotto_symbols))
    ),
    required = !c(
      rep(FALSE, length(class_symbols)),
      giotto_symbols %in% c("addStatistics", "calculateHVF")
    ),
    stringsAsFactors = FALSE
  )
}

spatial_giotto_converter_name <- function() {
  seurat_major <- suppressWarnings(as.integer(
    strsplit(as.character(utils::packageVersion("Seurat")), "\\.")[[1L]][1L]
  ))
  if (!is.na(seurat_major) && seurat_major >= 5L) {
    "giottoToSeuratV5"
  } else {
    "giottoToSeuratV4"
  }
}

spatial_backend_registry <- function() {
  backend <- function(id, package, repository = package, symbols = character()) {
    list(id = id, package = package, repository = repository, symbols = symbols)
  }
  giotto_symbols <- spatial_giotto_symbol_registry()
  list(
    core = backend("core", "scop"),
    biocneighbors = backend("biocneighbors", "BiocNeighbors", symbols = c("findKNN", "findNeighbors")),
    spanorm = backend("spanorm", "SpaNorm", symbols = "SpaNorm"),
    spatialqm = backend("spatialqm", "SpatialQM"),
    spotsweeper = backend("spotsweeper", "SpotSweeper"),
    giotto = backend(
      "giotto", "Giotto", "drieslab/Giotto",
      giotto_symbols$symbol[giotto_symbols$package == "Giotto" & giotto_symbols$required]
    ),
    giotto_class = backend(
      "giotto_class", "GiottoClass", "drieslab/GiottoClass",
      giotto_symbols$symbol[giotto_symbols$package == "GiottoClass" & giotto_symbols$required]
    ),
    spata2 = backend("spata2", "SPATA2", "theMILOlab/SPATA2", c("asSPATA2", "asSeurat")),
    music = backend("music", "MuSiC", "xuranw/MuSiC", "music_prop"),
    bisquerna = backend("bisquerna", "BisqueRNA", symbols = "ReferenceBasedDecomposition"),
    bayesprism = backend("bayesprism", "BayesPrism", symbols = c("new.prism", "run.prism")),
    cibersort = backend("cibersort", "CIBERSORT", "Moonerss/CIBERSORT", "cibersort"),
    spacexr = backend("spacexr", "spacexr", "dmcable/spacexr", c("createRctd", "runRctd")),
    card = backend("card", "CARD", "YingMa0107/CARD"),
    stdeconvolve = backend("stdeconvolve", "STdeconvolve", "JEFworks-Lab/STdeconvolve", c("cleanCounts", "fitLDA")),
    spotlight = backend("spotlight", "SPOTlight", symbols = "SPOTlight"),
    spatialecotyper = backend("spatialecotyper", "SpatialEcoTyper", "digitalcytometry/SpatialEcoTyper"),
    bayesspace = backend("bayesspace", "BayesSpace"),
    banksy = backend("banksy", "Banksy", symbols = "computeBanksy"),
    smoothclust = backend("smoothclust", "smoothclust", "lmweber/smoothclust"),
    meringue = backend("meringue", "MERINGUE"),
    sparkx = backend("sparkx", "SPARK"),
    nnsvg = backend("nnsvg", "nnSVG"),
    spicyr = backend("spicyr", "spicyR", symbols = "spicy"),
    statial = backend("statial", "Statial", symbols = "Kontextual"),
    precast = backend("precast", "PRECAST", "feiyoung/PRECAST", "CreatePRECASTObject"),
    bass = backend("bass", "BASS", symbols = "createBASSObject"),
    spatialmnn = backend("spatialmnn", "spatialMNN", "Pixel-Dream/spatialMNN", "spatialMNN"),
    mistyr = backend("mistyr", "mistyR", symbols = c("create_initial_view", "run_misty")),
    semla = backend("semla", "semla")
  )
}

spatial_registry_split_backends <- function(x) {
  unique(trimws(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)))
}

#' @title List spatial methods
#'
#' @description
#' List the public spatial analysis, visualization, framework bridge, and
#' workflow entry points registered by scop. This function does not inspect
#' optional packages unless `available` is `TRUE` or `FALSE`.
#'
#' @param task,kind,backend,status Optional exact filters.
#' @param available Optional logical availability filter. `NULL` avoids backend
#'   inspection.
#' @param pattern Optional case-insensitive pattern matched against method,
#'   task, backend, and documentation fields.
#'
#' @return A data frame with one row per registered public method.
#'
#' @export
ListSpatialMethods <- function(
  task = NULL,
  kind = NULL,
  backend = NULL,
  status = NULL,
  available = NULL,
  pattern = NULL
) {
  registry <- spatial_method_registry()
  exact_filter <- function(data, column, value) {
    if (is.null(value)) return(data)
    data[data[[column]] %in% value, , drop = FALSE]
  }
  registry <- exact_filter(registry, "task", task)
  registry <- exact_filter(registry, "kind", kind)
  registry <- exact_filter(registry, "status", status)
  if (!is.null(backend)) {
    keep <- vapply(
      registry$backend_id,
      function(x) any(spatial_registry_split_backends(x) %in% backend),
      logical(1)
    )
    registry <- registry[keep, , drop = FALSE]
  }
  if (!is.null(pattern)) {
    if (!is.character(pattern) || length(pattern) != 1L || is.na(pattern) || !nzchar(pattern)) {
      log_message("{.arg pattern} must be one non-empty string", message_type = "error")
    }
    haystack <- apply(
      registry[, c("method", "task", "backend_id", "documentation"), drop = FALSE],
      1,
      paste,
      collapse = " "
    )
    registry <- registry[grepl(pattern, haystack, ignore.case = TRUE), , drop = FALSE]
  }
  if (!is.null(available)) {
    if (!is.logical(available) || length(available) != 1L || is.na(available)) {
      log_message("{.arg available} must be TRUE, FALSE, or NULL", message_type = "error")
    }
    status_table <- SpatialBackendStatus(api_check = TRUE)
    ready <- stats::setNames(status_table$availability == "available", status_table$backend_id)
    method_ready <- vapply(seq_len(nrow(registry)), function(i) {
      ids <- spatial_registry_split_backends(registry$backend_id[[i]])
      values <- ready[ids]
      values[is.na(values)] <- FALSE
      if (identical(registry$backend_requirement[[i]], "any")) {
        any(values)
      } else {
        all(values)
      }
    }, logical(1))
    registry <- registry[method_ready == available, , drop = FALSE]
  }
  rownames(registry) <- NULL
  registry
}

spatial_backend_cache_key <- function(api_check, backend_ids, method) {
  paste(
    "v1",
    paste(normalizePath(.libPaths(), winslash = "/", mustWork = FALSE), collapse = "|"),
    isTRUE(api_check),
    paste(sort(backend_ids), collapse = "|"),
    paste(sort(method %||% character()), collapse = "|"),
    sep = "::"
  )
}

spatial_backend_required_symbols <- function(spec, method = NULL) {
  symbols <- spec$symbols
  if (
    identical(spec$id, "giotto_class") &&
      !is.null(method) &&
      "giotto_to_srt" %in% method
  ) {
    symbols <- c(symbols, spatial_giotto_converter_name())
  }
  unique(symbols)
}

#' @title Inspect spatial backend availability
#'
#' @description
#' Diagnose optional spatial backend installation and exported API status. The
#' function is read-only: it never installs or updates a package. Actual method
#' calls remain the authoritative runtime check for system libraries, S4
#' registration, and data-specific requirements.
#'
#' @param backend Optional backend identifier filter.
#' @param method Optional registered spatial method filter.
#' @param api_check Whether to inspect required namespace exports.
#' @param refresh Whether to refresh the session-local diagnostic cache.
#'
#' @return A data frame with one row per backend.
#'
#' @export
SpatialBackendStatus <- function(
  backend = NULL,
  method = NULL,
  api_check = TRUE,
  refresh = FALSE
) {
  if (!is.logical(api_check) || length(api_check) != 1L || is.na(api_check)) {
    log_message("{.arg api_check} must be TRUE or FALSE", message_type = "error")
  }
  if (!is.logical(refresh) || length(refresh) != 1L || is.na(refresh)) {
    log_message("{.arg refresh} must be TRUE or FALSE", message_type = "error")
  }
  backends <- spatial_backend_registry()
  backend_ids <- names(backends)
  if (!is.null(method)) {
    registry <- spatial_method_registry()
    registry <- registry[registry$method %in% method, , drop = FALSE]
    backend_ids <- intersect(
      backend_ids,
      spatial_registry_split_backends(registry$backend_id)
    )
  }
  if (!is.null(backend)) {
    backend_ids <- intersect(backend_ids, backend)
  }
  backends <- backends[backend_ids]
  if (length(backends) == 0L) {
    return(data.frame(
      backend_id = character(), package = character(), repository = character(),
      installed = logical(), api_checked = logical(), required_symbols = character(),
      missing_symbols = character(), availability = character(), namespace_error = character(),
      stringsAsFactors = FALSE
    ))
  }
  cache_key <- spatial_backend_cache_key(api_check, backend_ids, method)
  if (isTRUE(refresh) && exists(cache_key, envir = .spatial_backend_cache, inherits = FALSE)) {
    rm(list = cache_key, envir = .spatial_backend_cache)
  }
  if (!exists(cache_key, envir = .spatial_backend_cache, inherits = FALSE)) {
    installed <- rownames(utils::installed.packages())
    rows <- lapply(backends, function(spec) {
      package <- spec$package
      symbols <- spatial_backend_required_symbols(spec, method = method)
      is_core <- identical(package, "scop")
      is_installed <- is_core || package %in% installed
      namespace_error <- NA_character_
      missing_symbols <- character()
      if (!is_installed) {
        availability <- "missing"
      } else if (isTRUE(api_check) && length(symbols) > 0L) {
        exports <- tryCatch(
          getNamespaceExports(package),
          error = function(e) {
            namespace_error <<- conditionMessage(e)
            NULL
          }
        )
        if (is.null(exports)) {
          availability <- "namespace_error"
        } else {
          missing_symbols <- setdiff(symbols, exports)
          availability <- if (length(missing_symbols) == 0L) "available" else "api_incompatible"
        }
      } else {
        availability <- "available"
      }
      data.frame(
        backend_id = spec$id,
        package = package,
        repository = spec$repository,
        installed = is_installed,
        api_checked = isTRUE(api_check) && length(symbols) > 0L,
        required_symbols = paste(symbols, collapse = ","),
        missing_symbols = paste(missing_symbols, collapse = ","),
        availability = availability,
        namespace_error = namespace_error,
        stringsAsFactors = FALSE
      )
    })
    assign(cache_key, do.call(rbind, rows), envir = .spatial_backend_cache)
  }
  out <- get(cache_key, envir = .spatial_backend_cache, inherits = FALSE)
  rownames(out) <- NULL
  out
}

spatial_result_payload_size <- function(x) {
  if (is.null(x)) return(0L)
  if (is.data.frame(x) || is.matrix(x)) return(nrow(x))
  if (is.atomic(x)) return(length(x))
  if (is.list(x)) return(length(x))
  1L
}

spatial_result_state <- function(bundle, object_cells = NULL) {
  metadata_names <- c("method", "schema_version", "source", "parameters", "summary", "active_graph")
  if (!is.list(bundle)) {
    return(list(state = if (length(bundle) > 0L) "ready" else "empty", n_items = length(bundle), reason = NA_character_))
  }
  if (!is.null(bundle$graphs) && length(bundle$graphs) > 0L) {
    node_ids <- unique(unlist(lapply(bundle$graphs, function(graph) {
      nodes <- graph$nodes
      if (is.data.frame(nodes) && "cell_id" %in% colnames(nodes)) as.character(nodes$cell_id) else character()
    }), use.names = FALSE))
    stale <- length(node_ids) > 0L && !is.null(object_cells) && any(!node_ids %in% object_cells)
    return(list(
      state = if (stale) "stale" else "ready",
      n_items = length(bundle$graphs),
      reason = if (stale) "stored graph nodes are absent from object" else NA_character_
    ))
  }
  payload_names <- setdiff(names(bundle), metadata_names)
  payload_sizes <- vapply(bundle[payload_names], spatial_result_payload_size, integer(1))
  n_items <- sum(payload_sizes)
  if (n_items == 0L) {
    return(list(state = "empty", n_items = 0L, reason = "no logical result payload"))
  }
  source <- bundle$source
  partial <- is.null(bundle$method) && is.null(source) && is.null(bundle$parameters)
  list(
    state = if (partial) "partial" else "ready",
    n_items = n_items,
    reason = if (partial) "result provenance is incomplete" else NA_character_
  )
}

#' @title Inspect stored spatial results
#'
#' @description
#' Read spatial result bundles from a Seurat object's `tools` slot without
#' migrating, renaming, or modifying legacy results.
#'
#' @param object A `Seurat` object.
#' @param method Optional registered method filter.
#' @param tool_name Optional stored tool key filter.
#' @param include_empty Whether to include recognized keys without a logical
#'   result payload.
#'
#' @return A data frame describing recognized stored spatial results.
#'
#' @export
SpatialResultInfo <- function(
  object,
  method = NULL,
  tool_name = NULL,
  include_empty = FALSE
) {
  if (!inherits(object, "Seurat")) {
    log_message("{.arg object} must be a {.cls Seurat} object", message_type = "error")
  }
  if (!is.logical(include_empty) || length(include_empty) != 1L || is.na(include_empty)) {
    log_message("{.arg include_empty} must be TRUE or FALSE", message_type = "error")
  }
  registry <- spatial_method_registry()
  registry <- registry[!is.na(registry$tool_key) & nzchar(registry$tool_key), , drop = FALSE]
  if (!is.null(method)) registry <- registry[registry$method %in% method, , drop = FALSE]
  if (!is.null(tool_name)) registry <- registry[registry$tool_key %in% tool_name, , drop = FALSE]
  registry <- registry[order(is.na(registry$plot_function)), , drop = FALSE]
  registry <- registry[!duplicated(registry$tool_key), , drop = FALSE]
  cells <- colnames(object)
  rows <- lapply(seq_len(nrow(registry)), function(i) {
    key <- registry$tool_key[[i]]
    if (is.null(object@tools[[key]])) return(NULL)
    bundle <- object@tools[[key]]
    state <- spatial_result_state(bundle, object_cells = cells)
    if (!include_empty && identical(state$state, "empty")) return(NULL)
    source <- if (is.list(bundle)) bundle$source else NULL
    data.frame(
      method = if (is.list(bundle) && !is.null(bundle$method)) as.character(bundle$method[[1L]]) else registry$method[[i]],
      tool_name = key,
      schema_version = if (is.list(bundle) && !is.null(bundle$schema_version)) as.integer(bundle$schema_version[[1L]]) else NA_integer_,
      result_state = state$state,
      n_items = state$n_items,
      coordinate_space = if (!is.null(source$coordinate_space)) as.character(source$coordinate_space[[1L]]) else registry$coordinate_space_current[[i]],
      image = if (!is.null(source$image)) as.character(source$image[[1L]]) else NA_character_,
      parameters_available = is.list(bundle) && !is.null(bundle$parameters),
      summary_available = is.list(bundle) && !is.null(bundle$summary),
      plot_function = registry$plot_function[[i]],
      empty_reason = state$reason,
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(data.frame(
      method = character(), tool_name = character(), schema_version = integer(),
      result_state = character(), n_items = integer(), coordinate_space = character(),
      image = character(), parameters_available = logical(), summary_available = logical(),
      plot_function = character(), empty_reason = character(), stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
