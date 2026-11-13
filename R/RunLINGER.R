grn_run_linger <- function(input, config) {
  arguments <- config$method_args
  resume_dir <- arguments$resume_dir
  if (!is.null(resume_dir)) {
    if (!is.character(resume_dir) || length(resume_dir) != 1L || is.na(resume_dir) ||
      !file.exists(file.path(resume_dir, "prepare.completed.json")) ||
      dir.exists(file.path(resume_dir, "prepare.completed.json"))) {
      stop("LINGER resume_dir must contain a completed prepare checkpoint.", call. = FALSE)
    }
    resume_dir <- normalizePath(resume_dir, winslash = "/", mustWork = TRUE)
  }
  activation <- arguments$activation %||% "ReLU"
  if (!is.character(activation) || length(activation) != 1L || is.na(activation) || !activation %in% c("ReLU", "sigmoid", "tanh")) {
    stop("Unsupported LINGER activation.", call. = FALSE)
  }
  patch_missing_keys <- arguments$patch_missing_keys %||% TRUE
  if (!is.logical(patch_missing_keys) || length(patch_missing_keys) != 1L || is.na(patch_missing_keys)) {
    stop("LINGER patch_missing_keys must be TRUE or FALSE.", call. = FALSE)
  }
  min_genes <- arguments$min_genes %||% 200L
  min_cells <- arguments$min_cells %||% 3L
  for (value in list(min_genes, min_cells)) {
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value) || value < 0 || value != floor(value)) {
      stop("LINGER min_genes and min_cells must be nonnegative integers.", call. = FALSE)
    }
  }
  reference <- grn_linger_check_resources(input, config)
  selection <- reference$selection
  rna <- grn_linger_counts(input$object, config$assay, config$layer, input$cells)
  atac <- grn_linger_counts(input$object, config$atac_assay, config$atac_layer, input$cells)
  metadata <- grn_linger_metadata(input, config)
  if (ncol(rna) < 50L) {
    stop("LINGER requires at least 50 paired cells and 50 native pseudobulks for training.", call. = FALSE)
  }
  input_dir <- file.path(config$work_dir, "linger-input")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  Matrix::writeMM(rna, file.path(input_dir, "rna.mtx"))
  Matrix::writeMM(atac, file.path(input_dir, "atac.mtx"))
  writeLines(rownames(rna), file.path(input_dir, "genes.txt"))
  writeLines(rownames(atac), file.path(input_dir, "regions.txt"))
  utils::write.table(metadata, file.path(input_dir, "cells.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
  output_dir <- file.path(config$work_dir, "linger-output")
  request <- list(
    input_dir = normalizePath(input_dir, winslash = "/", mustWork = TRUE),
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
    bulk_reference = reference$reference,
    min_genes = as.integer(min_genes), min_cells = as.integer(min_cells),
    genome = selection$genome, activation = activation,
    patch_missing_keys = patch_missing_keys,
    seed = as.integer(config$seed), cores = as.integer(config$cores),
    regulators = as.list(input$regulators), targets = as.list(input$targets),
    expected_reference_files = as.list(selection$files),
    reference_selection = list(
      schema = selection$schema, genome = selection$genome,
      chromosomes = as.list(selection$chromosomes), manifest_sha256 = selection$manifest_sha256
    ),
    resume_dir = resume_dir,
    checkpoint_prepare = !startsWith(normalizePath(config$work_dir, winslash = "/", mustWork = TRUE), paste0(normalizePath(tempdir(), winslash = "/", mustWork = TRUE), "/")) || !is.null(resume_dir)
  )
  request_path <- file.path(input_dir, "request.json")
  runner_write_json(request, request_path)
  script <- system.file("grn", "linger.py", package = "scop", mustWork = TRUE)
  python <- config$runtime$python %||% conda_python(
    envname = config$envname %||% env_module_environment("linger"),
    conda = config$conda %||% "auto"
  )
  process <- runner_run_process(python, c(normalizePath(script), request_path),
    work_dir = config$work_dir, timeout = config$runtime$timeout %||% Inf,
    env = c(
      PYTHONHASHSEED = as.character(config$seed), OMP_NUM_THREADS = as.character(config$cores),
      OPENBLAS_NUM_THREADS = as.character(config$cores), CUDA_VISIBLE_DEVICES = ""
    )
  )
  receipt <- runner_read_json(file.path(output_dir, "receipt.json"))
  if (!identical(receipt$status, "completed")) {
    stop("LINGER did not produce a completed execution receipt.", call. = FALSE)
  }
  read_edges <- function(layer) {
    table <- utils::read.delim(gzfile(file.path(output_dir, "networks", paste0(layer, ".tsv.gz"))),
      colClasses = c("character", "character", "numeric"), check.names = FALSE, stringsAsFactors = FALSE
    )
    if (any(!is.finite(table$score))) {
      stop("LINGER returned non-finite network scores.", call. = FALSE)
    }
    table
  }
  native <- list(
    tf_gene = read_edges("tf_gene"), tf_region = read_edges("tf_region"),
    region_gene = read_edges("region_gene"), receipt = receipt,
    training_coverage = utils::read.delim(file.path(output_dir, "training-coverage.tsv"), stringsAsFactors = FALSE)
  )
  specs <- list()
  for (layer in c("tf_gene", "tf_region", "region_gene")) {
    specs[[layer]] <- list(
      status = "computed", table_path = layer,
      source = if (identical(layer, "region_gene")) "region" else "regulator",
      target = if (identical(layer, "tf_region")) "region" else "target",
      score = "score", source_type = if (identical(layer, "region_gene")) "region" else "TF",
      target_type = if (identical(layer, "tf_region")) "region" else "gene",
      score_type = "native_linger_score", directed = TRUE,
      sign_column = NULL, sign_basis = "unknown", rank_by = "score"
    )
  }
  retained <- utils::read.delim(file.path(output_dir, "retained_cells.tsv"), colClasses = "character", stringsAsFactors = FALSE)
  list(
    cells = as.character(retained$cell_id), native = native, layer_specs = specs,
    chains = list(status = "not_supported"), gene_sets = list(status = "not_supported"),
    parameters = c(request, list(
      training_seed = 42L, candidate_filter = "preprocess_gene_union_and_postfit",
      input_cells = input$cells, excluded_cells = setdiff(input$cells, retained$cell_id),
      cell_filter = "native min_genes/min_cells in both modalities"
    )),
    implementation = receipt$implementation,
    package_versions = unlist(receipt$versions, use.names = TRUE),
    execution = list(
      status = "completed", artifacts = list(output_dir = output_dir),
      stdout = process$stdout, stderr = process$stderr, runtime = list(python = python)
    )
  )
}

.linger_reference_spec <- list(
  schema = "linger-1.110-reference-v1",
  chromosomes = paste0("chr", c(seq_len(22L), "X")),
  common = c(
    "TFName.txt", "Match2.txt", "bulk_gene_all.txt", "all_hg19.txt",
    "motifWeight.txt", "RE_gene_corr_hg19.bed"
  ),
  genome_pair = c(hg19 = "hg19_hg38_pair.bed", hg38 = "hg38_hg19_pair.bed"),
  per_chromosome = c(
    "MotifTarget_Matrix_{chromosome}.txt",
    "MotifTarget_matrix_{chromosome}.bed", "{genome}_Peaks_{chromosome}.bed",
    "{chromosome}_gene.txt", "{chromosome}_index.txt", "{chromosome}_index_all.txt",
    "fisher_{chromosome}.pt", "all_models_{chromosome}.pt"
  )
)

grn_linger_reference_files <- function(input, config) {
  spec <- .linger_reference_spec
  genome <- config$method_args$genome %||% "hg38"
  if (!is.character(genome) || length(genome) != 1L || is.na(genome) || !genome %in% names(spec$genome_pair)) {
    stop("LINGER supports genome='hg19' or genome='hg38'.", call. = FALSE)
  }
  features <- rownames(input$object[[config$atac_assay]])
  if (!length(features) || anyNA(features)) {
    stop("LINGER accessibility peaks require complete coordinates.", call. = FALSE)
  }
  pattern <- "^([^:-]+)[:-]([0-9]+)-([0-9]+)$"
  matches <- regmatches(features, regexec(pattern, features, perl = TRUE))
  if (any(lengths(matches) != 4L)) {
    stop("LINGER requires valid primary-chromosome peak coordinates.", call. = FALSE)
  }
  chromosomes <- vapply(matches, `[[`, "", 2L)
  starts <- as.numeric(vapply(matches, `[[`, "", 3L))
  ends <- as.numeric(vapply(matches, `[[`, "", 4L))
  if (any(!chromosomes %in% spec$chromosomes) || any(!is.finite(starts) | !is.finite(ends) | ends <= starts)) {
    stop("LINGER requires valid primary-chromosome peak coordinates.", call. = FALSE)
  }
  chromosomes <- sort(unique(chromosomes), method = "radix")
  files <- c(spec$common, spec$genome_pair[[genome]], unlist(lapply(chromosomes, function(chromosome) {
    gsub("{genome}", genome, gsub("{chromosome}", chromosome, spec$per_chromosome, fixed = TRUE), fixed = TRUE)
  }), use.names = FALSE))
  list(
    files = files, schema = spec$schema, genome = genome, chromosomes = chromosomes,
    manifest_sha256 = digest::digest(serialize(spec, NULL, xdr = FALSE, version = 3L), algo = "sha256", serialize = FALSE)
  )
}

grn_linger_check_resources <- function(input, config) {
  reference <- config$resources$bulk_reference
  if (!is.character(reference) || length(reference) != 1L || !dir.exists(reference)) {
    stop("LINGER requires resources$bulk_reference pointing to the data_bulk directory.", call. = FALSE)
  }
  reference <- normalizePath(reference, winslash = "/", mustWork = TRUE)
  selection <- grn_linger_reference_files(input, config)
  paths <- file.path(reference, selection$files)
  missing <- selection$files[!file.exists(paths) | dir.exists(paths)]
  if (length(missing)) {
    stop("Missing LINGER bulk reference files: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  list(reference = reference, selection = selection)
}

grn_linger_metadata <- function(input, config) {
  arguments <- config$method_args
  metadata <- input$object[[]]
  cells <- as.character(input$cells)
  if (!all(cells %in% rownames(metadata))) {
    stop("LINGER cell identifiers are missing from Seurat metadata.", call. = FALSE)
  }
  metadata <- metadata[cells, , drop = FALSE]
  read_column <- function(column, fallback) {
    if (is.null(column)) {
      return(rep(fallback, length(cells)))
    }
    if (length(column) != 1L || !is.character(column) || !column %in% colnames(metadata)) {
      stop("LINGER sample.by and label.by must name existing metadata columns.", call. = FALSE)
    }
    values <- as.character(metadata[[column]])
    if (anyNA(values) || any(!nzchar(values))) {
      stop("LINGER sample and label metadata must be complete.", call. = FALSE)
    }
    values
  }
  data.frame(
    cell_id = cells,
    sample = read_column(arguments$sample.by, "sample1"),
    label = read_column(arguments$label.by, input$group %||% "all"),
    stringsAsFactors = FALSE
  )
}

grn_linger_counts <- function(object, assay, layer, cells) {
  if (!identical(layer, "counts")) {
    stop("LINGER requires RNA and ATAC count layers; set layer and atac_layer to 'counts'.", call. = FALSE)
  }
  if (is.null(assay) || !assay %in% SeuratObject::Assays(object)) {
    stop("LINGER requires the selected RNA and ATAC assays.", call. = FALSE)
  }
  available <- SeuratObject::Layers(object[[assay]], search = NA)
  layers <- available[grepl("^counts($|\\.)", available)]
  if (!identical(layers, "counts")) {
    stop("LINGER requires one counts layer per assay. Join split layers with SeuratObject::JoinLayers() first.", call. = FALSE)
  }
  matrix <- SeuratObject::LayerData(object, assay = assay, layer = layer)
  if (is.null(rownames(matrix)) || anyDuplicated(rownames(matrix)) || anyNA(rownames(matrix)) || any(!nzchar(rownames(matrix)))) {
    stop("LINGER count features must have unique complete identifiers.", call. = FALSE)
  }
  if (!all(cells %in% colnames(matrix))) {
    stop("LINGER RNA and ATAC cells must match the selected group.", call. = FALSE)
  }
  matrix <- methods::as(matrix[, cells, drop = FALSE], "dgCMatrix")
  if (any(!is.finite(matrix@x)) || any(matrix@x < 0) || any(matrix@x != floor(matrix@x))) {
    stop("LINGER requires finite nonnegative integer RNA and ATAC counts.", call. = FALSE)
  }
  matrix
}
