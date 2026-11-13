grn_multicsn_arguments <- function() {
  list(
    prepare = c("motif_backend", "exclude_exons", "align_annotation"),
    multiome = c(
      "preset", "min_detected", "min_peak_cells", "upstream", "downstream",
      "tf_region_scope", "exclude_response_alias", "renormalize", "response_chunk",
      "target_chunk", "max_support_size", "min_improvement", "sort_regulators"
    ),
    rna = c("method", "max_support_size", "lag_fraction", "lag_steps")
  )
}

grn_prepare_multicsn <- function(input, config) {
  if (!requireNamespace("multiCSN", quietly = TRUE) ||
      utils::packageVersion("multiCSN") < "0.0.2") {
    stop("multiCSN >= 0.0.2 is required. Install the current independent multiCSN package.", call. = FALSE)
  }
  mode <- input$input_mode
  if (is.null(mode)) mode <- config$input_mode
  if (!mode %in% c("rna", "multiome")) {
    stop("multiCSN input_mode must be 'rna' or 'multiome'.", call. = FALSE)
  }
  argument_sets <- grn_multicsn_arguments()
  allowed <- if (identical(mode, "rna")) {
    argument_sets$rna
  } else {
    c(argument_sets$prepare, argument_sets$multiome)
  }
  unknown <- setdiff(names(config$method_args), allowed)
  if (length(unknown)) {
    stop("Unsupported multiCSN ", mode, " argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  input$input_mode <- mode
  if (identical(mode, "rna")) return(input)
  if (!inherits(input$object, "Seurat")) {
    stop("multiCSN multiome input requires a Seurat object.", call. = FALSE)
  }
  if (!identical(config$layer, "data") || !identical(config$atac_layer, "data")) {
    stop("multiCSN multiome inference uses the RNA and ATAC 'data' layers; select layer = 'data' and atac_layer = 'data'.", call. = FALSE)
  }
  if (is.null(config$resources$genome)) {
    stop("multiCSN multiome inference requires resources$genome.", call. = FALSE)
  }
  object <- input$object[, input$cells]
  groups <- input$groups
  if (is.null(groups)) groups <- config$groups
  grouped_cells <- unlist(groups, use.names = FALSE)
  if (is.null(names(groups)) || anyDuplicated(grouped_cells) ||
      !setequal(grouped_cells, input$cells)) {
    stop("multiCSN groups must assign every selected cell exactly once.", call. = FALSE)
  }
  labels <- stats::setNames(rep(NA_character_, length(input$cells)), input$cells)
  for (group in names(groups)) labels[groups[[group]]] <- group
  if (anyNA(labels) || !identical(names(labels), input$cells)) {
    stop("multiCSN groups must assign every selected cell exactly once.", call. = FALSE)
  }
  group_column <- ".scop_grn_group"
  while (group_column %in% colnames(object[[]])) group_column <- paste0(group_column, "_")
  object[[group_column]] <- unname(labels[colnames(object)])
  if (!is.null(config$resources$annotation)) {
    Signac::Annotation(object[[config$atac_assay]]) <- config$resources$annotation
  }
  prepare_args <- config$method_args[intersect(names(config$method_args), argument_sets$prepare)]
  if (!is.null(prepare_args$motif_backend)) {
    prepare_args$backend <- prepare_args$motif_backend
    prepare_args$motif_backend <- NULL
  }
  object <- do.call(multiCSN::prepare_candidates, c(list(
    object = object, group.by = group_column, genome = config$resources$genome,
    rna_assay = config$assay, peak_assay = config$atac_assay,
    targets = input$targets, pfm = config$resources$motif,
    motif_tfs = config$resources$motif_tfs, verbose = config$verbose
  ), prepare_args))
  input$object <- object
  candidate_regulators <- multiCSN::get_attribute(object, attribute = "tfs")
  if (!is.null(input$regulators)) {
    input$regulators <- intersect(input$regulators, candidate_regulators)
    if (!length(input$regulators)) {
      stop("No requested regulators have motif-supported TF candidates in the multiome input.", call. = FALSE)
    }
  }
  input$multicsn_preparation <- list(
    rebuilt = TRUE, group_column = group_column, cells = input$cells,
    rna_assay = config$assay, peak_assay = config$atac_assay,
    motif_backend = if (is.null(prepare_args$backend)) "motifmatchr" else prepare_args$backend,
    targets = input$targets,
    candidate_regulators = candidate_regulators
  )
  input
}

grn_multicsn_layer <- function(table_path, source, target,
                               sign_column = "standardized_beta",
                               sign_basis = "regression_coefficient") {
  list(
    status = "computed", table_path = table_path,
    source = source, target = target, score = "weight",
    source_type = if (identical(source, "region")) "region" else "TF",
    target_type = if (identical(target, "region")) "region" else "gene",
    score_type = "signed_ordinal_weight", directed = TRUE,
    sign_column = sign_column, sign_basis = sign_basis,
    rank_by = "absolute_score"
  )
}

grn_run_multicsn <- function(input, config) {
  arguments <- grn_multicsn_arguments()
  mode <- input$input_mode
  if (is.null(mode)) mode <- config$input_mode
  versions <- vapply(c("multiCSN", "inferCSN"), function(package) {
    as.character(utils::packageVersion(package))
  }, character(1))
  if (identical(mode, "rna")) {
    method_args <- config$method_args[intersect(names(config$method_args), arguments$rna)]
    expression <- if (inherits(input$expression, "Matrix")) {
      Matrix::t(input$expression)
    } else {
      t(as.matrix(input$expression))
    }
    native <- do.call(multiCSN::inferCSN, c(list(
      object = expression, pseudotime = input$pseudotime,
      regulators = input$regulators, targets = input$targets,
      cores = config$cores, verbose = config$verbose
    ), method_args))
    parameters <- utils::modifyList(list(
      method = "greedy_l0", max_support_size = NULL, lag_fraction = 0.05, lag_steps = NULL
    ), method_args, keep.null = TRUE)
    return(list(
      cells = input$cells, native = native,
      layer_specs = list(tf_gene = grn_multicsn_layer(
        character(), "regulator", "target", "weight", "signed_ordinal_weight"
      )),
      chains = list(status = "not_supported"),
      gene_sets = list(status = "not_computed"),
      parameters = c(list(input_mode = "rna", pseudotime = !is.null(input$pseudotime)), parameters),
      implementation = "inferCSN::inferCSN",
      package_versions = versions,
      execution = list(status = "completed", artifacts = list())
    ))
  }
  if (is.null(input$multicsn_preparation)) {
    stop("Prepare multiCSN candidates before running a multiome group.", call. = FALSE)
  }
  method_args <- config$method_args[intersect(names(config$method_args), arguments$multiome)]
  native <- do.call(multiCSN::fit_layered_network, c(list(
    object = input$object, cell_group = input$group, cells = input$cells,
    regulators = input$regulators, targets = input$targets,
    cores = config$cores, store = FALSE, verbose = config$verbose
  ), method_args))
  list(
    cells = native$cells, native = native,
    layer_specs = list(
      tf_gene = grn_multicsn_layer("tf_gene", "regulator", "target"),
      tf_region = grn_multicsn_layer("tf_region", "regulator", "region"),
      region_gene = grn_multicsn_layer("region_gene", "region", "target"),
      mediated_tf_gene = grn_multicsn_layer(
        "mediated_tf_gene", "regulator", "target", "weight", "strongest_supported_path"
      )
    ),
    chains = list(status = "computed", table_path = "triplets"),
    gene_sets = list(status = "not_computed"),
    parameters = c(native$settings, list(
      input_mode = "multiome", preparation = input$multicsn_preparation
    )),
    implementation = "multiCSN::fit_layered_network",
    package_versions = versions,
    execution = list(status = "completed", artifacts = list())
  )
}
