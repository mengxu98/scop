grn_scenicplus_check_resources <- function(resources) {
  required <- c("cistopic", "gex", "region_sets", "ranking_db", "scores_db", "motif_annotation", "genome_annotation", "chromsizes")
  absent <- required[!vapply(required, function(name) {
    path <- resources[[name]]
    is.character(path) && length(path) == 1L && !is.na(path) && file.exists(path)
  }, logical(1))]
  if (length(absent)) {
    stop("Missing SCENIC+ resources: ", paste(absent, collapse = ", "), call. = FALSE)
  }
  files <- setdiff(required, "region_sets")
  invalid <- files[vapply(files, function(name) {
    path <- resources[[name]]
    dir.exists(path) || is.na(file.info(path)$size) || file.info(path)$size < 1
  }, logical(1))]
  if (length(invalid)) stop("SCENIC+ resources must be nonempty files: ", paste(invalid, collapse = ", "), call. = FALSE)
  if (!dir.exists(resources$region_sets)) stop("SCENIC+ region_sets must be a directory.", call. = FALSE)
  if (length(Sys.glob(file.path(resources$region_sets, "*", "*.bed"))) < 2L) {
    stop("Official SCENIC+ DEM requires at least two region sets to define foreground and background.", call. = FALSE)
  }
  if (!is.null(resources$workflow)) {
    workflow <- resources$workflow
    if (!is.character(workflow) || length(workflow) != 1L || is.na(workflow) ||
        !all(file.exists(file.path(workflow, c("Snakefile", "config.yaml"))))) {
      stop("SCENIC+ workflow must contain the official Snakefile and config.yaml.", call. = FALSE)
    }
  }
  invisible(required)
}

grn_run_scenicplus <- function(input, config) {
  required <- grn_scenicplus_check_resources(config$resources)
  runtime <- if (!is.null(config$runtime$executable)) config$runtime else grn_resolve_runtime(config, "scenicplus")
  config$runtime <- runtime
  script <- system.file("grn", "scenicplus.py", package = "scop")
  if (!nzchar(script)) stop("The installed official SCENIC+ runner is missing; reinstall scop.", call. = FALSE)
  dir.create(config$work_dir, recursive = TRUE, showWarnings = FALSE)
  directory <- normalizePath(config$work_dir, winslash = "/", mustWork = TRUE)
  expression_path <- file.path(directory, "scenicplus-expression.mtx")
  feature_path <- file.path(directory, "scenicplus-features.tsv")
  if (!identical(colnames(input$expression), input$cells)) {
    stop("SCENIC+ expression columns must match the requested cells.", call. = FALSE)
  }
  if (any(grepl("[\r\n]", rownames(input$expression)))) {
    stop("SCENIC+ feature identifiers cannot contain line breaks.", call. = FALSE)
  }
  Matrix::writeMM(Matrix::Matrix(input$expression, sparse = TRUE), expression_path)
  writeLines(rownames(input$expression), feature_path, useBytes = TRUE)
  request_path <- file.path(directory, "scenicplus-request.json")
  output <- config$method_args$resume_dir %||% file.path(directory, "official")
  if (!is.character(output) || length(output) != 1L || is.na(output) || !nzchar(output)) {
    stop("SCENIC+ resume_dir must be one directory path.", call. = FALSE)
  }
  resources <- config$resources[intersect(c(required, "workflow"), names(config$resources))]
  resources <- lapply(resources, normalizePath, winslash = "/", mustWork = TRUE)
  jsonlite::write_json(list(resources = resources, cells = unname(input$cells), group = input$group,
    cores = config$cores, seed = config$seed, method_args = config$method_args,
    resource_identity = config$resource_identity,
    expression = expression_path, features = feature_path, output_dir = output),
    request_path, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
  launcher <- "import runpy,sys; sys.path.pop(0) if sys.path[:1] == [''] else None; sys.argv = sys.argv[1:]; runpy.run_path(sys.argv[0], run_name='__main__')"
  process <- grn_run_process(runtime$python, c("-c", launcher, script, request_path), config,
    env = c(PYTHONHASHSEED = as.character(config$seed), OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1"))
  grn_scenicplus_group(output, input, config, process)
}

grn_scenicplus_group <- function(directory, input, config, process = list()) {
  manifest <- file.path(directory, "execution.json")
  if (!file.exists(manifest)) stop("Official SCENIC+ did not produce its execution receipt.", call. = FALSE)
  receipt <- jsonlite::read_json(manifest, simplifyVector = TRUE)
  if (!identical(receipt$status, "completed") || !setequal(receipt$cells, input$cells)) {
    stop("Official SCENIC+ execution is incomplete or has different cells.", call. = FALSE)
  }
  filenames <- c(tf_gene = "tf_gene.tsv", tf_region = "tf_region.tsv", region_gene = "region_gene.tsv",
    direct = "direct.tsv", extended = "extended.tsv", chains = "chains.tsv")
  native <- lapply(filenames, function(filename) {
    path <- file.path(directory, filename)
    if (!file.exists(path)) stop("Missing official SCENIC+ result: ", filename, call. = FALSE)
    header <- names(utils::read.delim(path, nrows = 0L, check.names = FALSE))
    identifiers <- intersect(header, c("TF", "target", "region", "Region", "Gene", "eRegulon_name", "Gene_signature_name", "Region_signature_name", "annotation_type"))
    classes <- rep(NA_character_, length(header))
    classes[header %in% identifiers] <- "character"
    utils::read.delim(path, check.names = FALSE, colClasses = classes, na.strings = NULL)
  })
  native$gene_sets <- jsonlite::read_json(file.path(directory, "gene_sets.json"), simplifyVector = TRUE)
  native$auc <- lapply(receipt$auc_files, function(filename) {
    path <- file.path(directory, filename)
    header <- names(utils::read.delim(path, nrows = 0L, check.names = FALSE))
    classes <- c("character", rep(NA_character_, length(header) - 1L))
    table <- utils::read.delim(path, colClasses = classes, check.names = FALSE, na.strings = NULL)
    rownames(table) <- table[[1L]]
    table <- table[-1L]
    if (!setequal(rownames(table), input$cells)) stop("Official SCENIC+ AUC cells differ from the requested cells.", call. = FALSE)
    table[input$cells, , drop = FALSE]
  })
  names(native$auc) <- tools::file_path_sans_ext(receipt$auc_files)
  regulators <- input$regulators %||% unique(native$tf_gene$TF)
  targets <- input$targets %||% unique(native$tf_gene$target)
  selected_rows <- list(
    tf_gene = which(native$tf_gene$TF %in% regulators & native$tf_gene$target %in% targets),
    tf_region = which(native$tf_region$TF %in% regulators),
    region_gene = which(native$region_gene$target %in% targets)
  )
  chain_rows <- which(native$chains$TF %in% regulators & native$chains$Gene %in% targets)
  native$gene_sets_view <- if ("Gene_signature_name" %in% names(native$chains)) {
    lapply(split(native$chains$Gene[chain_rows], native$chains$Gene_signature_name[chain_rows]), unique)
  } else {
    lapply(native$gene_sets, intersect, targets)
  }
  native$candidate_scope <- list(selection = "postfit_view", requested_regulators = regulators, requested_targets = targets,
    fitted_regulators = unique(native$tf_gene$TF), fitted_targets = unique(native$tf_gene$target))
  native$parameters <- receipt$parameters
  native$execution <- receipt
  columns <- list(tf_gene = c("TF", "target", "importance", "TF", "gene", "rho", "feature_importance"),
    tf_region = c("TF", "Region", "membership", "TF", "region", "", "membership"),
    region_gene = c("region", "target", "importance", "region", "gene", "rho", "feature_importance"))
  specs <- lapply(names(columns), function(layer) {
    fields <- columns[[layer]]
    if (!all(fields[1:3] %in% names(native[[layer]]))) {
      stop("Invalid official SCENIC+ ", layer, " schema.", call. = FALSE)
    }
    membership <- identical(layer, "tf_region")
    list(status = "computed", table_path = layer, rows = selected_rows[[layer]], source = fields[[1L]], target = fields[[2L]], score = fields[[3L]],
      source_type = fields[[4L]], target_type = fields[[5L]], score_type = fields[[7L]], directed = TRUE,
      sign_column = if (membership) NULL else fields[[6L]], sign_basis = if (membership) "unknown" else "spearman_correlation",
      rank_by = "score")
  })
  names(specs) <- names(columns)
  list(cells = input$cells, native = native, layer_specs = specs,
    chains = list(status = "computed", table_path = "chains", rows = chain_rows, source = "TF", region = "Region", target = "Gene", score = "triplet_rank", rank_by = "ascending_score"),
    gene_sets = list(status = "computed", table_path = "gene_sets_view"),
    parameters = receipt$parameters, implementation = "SCENIC+::official_snakemake",
    package_versions = unlist(receipt$package_versions, use.names = TRUE),
    execution = c(list(status = "completed", artifacts = list(object = file.path(directory, "scplusmdata.h5mu"),
      manifest = manifest, stdout = file.path(directory, "snakemake.stdout.log"), stderr = file.path(directory, "snakemake.stderr.log"))),
      process[setdiff(names(process), "status")]))
}
