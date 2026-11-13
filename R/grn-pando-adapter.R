grn_pando_defaults <- function() {
  list(
    exclude_exons = TRUE, nfeatures = 2000L,
    peak_to_gene_method = "GREAT", upstream = 100000, downstream = 0,
    extend = 1000000, only_tss = FALSE, tf_cor = 0.1, peak_cor = 0,
    family = "gaussian", interaction_term = ":", adjust_method = "fdr",
    scale = FALSE, p_thresh = 0.05, rsq_thresh = 0.1,
    nvar_thresh = 10L, min_genes_per_module = 5L
  )
}

grn_pando_preserve_session <- function() {
  previous_search <- search()
  previous_options <- options()
  had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_rng) previous_rng <- get(".Random.seed", envir = .GlobalEnv)
  pando_loaded <- "Pando" %in% loadedNamespaces()
  s3_tables <- lapply(c(stats = "stats", base = "base"), function(package) {
    get(".__S3MethodsTable__.", asNamespace(package), inherits = FALSE)
  })
  s3_keys <- list(
    stats = paste0("coef.", c("GRNData", "Network", "RegulatoryNetwork")),
    base = paste0("print.", c("Modules", "Network", "Regions", "RegulatoryNetwork"))
  )
  s3_methods <- Map(function(table, keys) {
    stats::setNames(lapply(keys, get0, envir = table, inherits = FALSE), keys)
  }, s3_tables, s3_keys)
  registrations <- lapply(c(foreach = ".foreachGlobals", doParallel = ".options"), function(name) NULL)
  for (package in names(registrations)) {
    if (package %in% loadedNamespaces()) {
      name <- if (package == "foreach") ".foreachGlobals" else ".options"
      registrations[[package]] <- as.list(get(name, asNamespace(package)), all.names = TRUE)
    }
  }
  function() {
    for (name in setdiff(search(), previous_search)) detach(name, character.only = TRUE, unload = FALSE)
    if (!pando_loaded && "Pando" %in% loadedNamespaces()) {
      unloadNamespace("Pando")
      for (package in names(s3_methods)) {
        table <- s3_tables[[package]]
        for (key in names(s3_methods[[package]])) {
          method <- s3_methods[[package]][[key]]
          if (is.null(method)) {
            if (exists(key, envir = table, inherits = FALSE)) rm(list = key, envir = table)
          } else {
            assign(key, method, envir = table)
          }
        }
      }
    }
    extra_options <- setdiff(names(options()), names(previous_options))
    if (length(extra_options)) options(stats::setNames(rep(list(NULL), length(extra_options)), extra_options))
    options(previous_options)
    if (had_rng) assign(".Random.seed", previous_rng, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    for (package in names(registrations)) {
      if (!package %in% loadedNamespaces()) next
      name <- if (package == "foreach") ".foreachGlobals" else ".options"
      registration <- get(name, asNamespace(package))
      previous <- registrations[[package]]
      extra <- setdiff(ls(registration, all.names = TRUE), names(previous))
      if (length(extra)) rm(list = extra, envir = registration)
      for (key in names(previous)) {
        if (!exists(key, envir = registration, inherits = FALSE) || !identical(get(key, envir = registration, inherits = FALSE), previous[[key]])) {
          assign(key, previous[[key]], envir = registration)
        }
      }
    }
    invisible(NULL)
  }
}

grn_pando_current_runtime <- function(config) {
  spec <- grn_runtime_profiles()[["pando"]]
  packages <- unique(c(spec$packages, grn_pando_runtime_packages(config)))
  genome <- config$resources$genome
  if (isS4(genome) && identical(attr(class(genome), "package"), "BSgenome")) genome <- attr(genome, "pkgname")
  if (is.character(genome) && length(genome) == 1L && startsWith(genome, "BSgenome.")) packages <- unique(c(packages, genome))
  spec$packages <- packages
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Pando requires package(s) in the current R session: ", paste(missing, collapse = ", "), ". Install them or select method_args$runtime$rscript and r_library.", call. = FALSE)
  executable <- file.path(R.home("bin"), "Rscript")
  resolved <- list(kind = "r", executable = executable, python = NULL, rscript = executable,
    r_library = .libPaths(), timeout = Inf, version = as.character(getRversion()),
    package_versions = stats::setNames(vapply(packages, function(package) as.character(getNamespaceVersion(package)), character(1)), packages),
    environment = NULL)
  structure(resolved, class = c("GRNRuntime", "list"), profile = grn_hash(spec))
}

grn_run_pando <- function(input, config) {
  if (!inherits(input$object, "Seurat")) {
    stop("Pando requires a Seurat object with RNA and chromatin assays.", call. = FALSE)
  }
  if (is.null(config$resources$genome)) {
    stop("Pando requires resources$genome matching the chromatin coordinates.", call. = FALSE)
  }
  supplied <- config$method_args
  if (length(supplied) && (is.null(names(supplied)) || any(!nzchar(names(supplied))) || anyDuplicated(names(supplied)))) {
    stop("Pando method_args must have unique names.", call. = FALSE)
  }
  defaults <- grn_pando_defaults()
  unknown <- setdiff(names(supplied), names(defaults))
  if (length(unknown)) stop("Unknown Pando argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  config$method_args <- utils::modifyList(defaults, supplied)
  if (!identical(config$method_args$interaction_term, ":")) {
    stop("Pando GRN layers require interaction_term=':' to identify TF-region terms.", call. = FALSE)
  }
  if (!identical(as.numeric(config$method_args$extend), 1000000)) {
    stop("Pando 1.1.1 does not forward a custom extend value; supply resources$gene_domains for explicit gene domains.", call. = FALSE)
  }
  isolated <- if (!is.null(config$runtime_isolation)) isTRUE(config$runtime_isolation) else !inherits(config$runtime, "GRNRuntime") && any(c("rscript", "r_library", "timeout") %in% names(config$runtime))
  script <- system.file("grn", "pando", if (isolated) "run.R" else "native.R", package = "scop")
  if (!nzchar(script)) {
    stop("The installed Pando runner is missing; reinstall scop.", call. = FALSE)
  }
  dir.create(config$work_dir, recursive = TRUE, showWarnings = FALSE)
  object <- subset(input$object, cells = input$cells)
  request <- list(
    object = object,
    expression = input$expression,
    cells = input$cells,
    regulators = input$regulators,
    targets = input$targets,
    config = config
  )
  if (!isolated) {
    if (is.null(config$runtime_isolation)) {
      cleanup <- grn_pando_preserve_session()
      on.exit(cleanup(), add = TRUE)
    }
    runner <- new.env(parent = environment(grn_run_pando))
    sys.source(script, envir = runner)
    native <- runner$grn_pando_fit(request)
    return(grn_pando_group(native, input, list(mode = "current_session", r_version = as.character(getRversion()), r_library = .libPaths())))
  }
  runtime <- grn_resolve_runtime(config, profile = "pando")
  config$runtime <- runtime
  request$config <- config
  request_path <- tempfile("pando-request-", tmpdir = config$work_dir, fileext = ".rds")
  result_path <- tempfile("pando-result-", tmpdir = config$work_dir, fileext = ".rds")
  saveRDS(request, request_path)
  process <- grn_run_process(
    runtime$rscript,
    c("--vanilla", script, request_path, result_path),
    config,
    env = if (length(runtime$r_library)) c(R_LIBS_USER = paste(runtime$r_library, collapse = .Platform$path.sep)) else character()
  )
  if (!file.exists(result_path)) {
    stop("Pando completed without writing its result.", call. = FALSE)
  }
  native <- readRDS(result_path)
  process$mode <- "subprocess"
  process$r_version <- runtime$version
  process$r_library <- runtime$r_library
  grn_pando_group(native, input, process)
}

grn_pando_group <- function(native, input, process = list()) {
  required <- c("coefficients", "goodness_of_fit", "modules", "triplets", "tf_region", "region_gene", "gene_sets", "parameters", "package_versions")
  if (!is.list(native) || !all(required %in% names(native))) {
    stop("Pando returned an incomplete native result.", call. = FALSE)
  }
  tables <- list(
    tf_gene = c("modules", "tf", "target", "estimate", "TF", "gene", "aggregated_interaction_coefficient"),
    tf_region = c("tf_region", "tf", "region", "membership", "TF", "region", "membership"),
    region_gene = c("region_gene", "region", "target", "membership", "region", "gene", "membership")
  )
  specs <- lapply(names(tables), function(layer) {
    columns <- tables[[layer]]
    table <- native[[columns[[1L]]]]
    if (!is.data.frame(table) || !all(columns[2:4] %in% names(table))) {
      stop("Pando returned an invalid ", layer, " table.", call. = FALSE)
    }
    membership <- identical(columns[[7L]], "membership")
    list(
      status = "computed", table_path = columns[[1L]],
      source = columns[[2L]], target = columns[[3L]], score = columns[[4L]],
      source_type = columns[[5L]], target_type = columns[[6L]],
      score_type = columns[[7L]], directed = TRUE,
      sign_column = if (membership) NULL else "estimate",
      sign_basis = if (membership) "unknown" else "aggregated_interaction_coefficient",
      rank_by = if (membership) "score" else "absolute_score"
    )
  })
  names(specs) <- names(tables)
  list(
    cells = input$cells,
    native = native,
    layer_specs = specs,
    chains = list(status = "computed", table_path = "triplets"),
    gene_sets = list(status = "computed", table_path = "gene_sets"),
    parameters = native$parameters,
    implementation = "Pando::infer_grn/glm",
    package_versions = native$package_versions,
    execution = utils::modifyList(process, list(status = "completed", artifacts = list()))
  )
}
