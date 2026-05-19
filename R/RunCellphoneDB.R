#' @title Run CellphoneDB analysis
#'
#' @md
#' @inheritParams thisutils::log_message
#' @inheritParams standard_scop
#' @inheritParams CellDimPlot
#' @param species Species of the input data. CellphoneDB is human-centric; when a
#' non-human species is supplied and `convert_to_human = TRUE`, gene symbols will
#' be converted to human ortholog symbols before running the analysis.
#' @param assay Assay to use.
#' @param layer Layer to use as CellphoneDB input. Default is `"data"` because
#' CellphoneDB generally expects normalized expression.
#' @param counts_data Gene identifier type used by CellphoneDB. One of
#' `"hgnc_symbol"`, `"ensembl"`, or `"gene_name"`.
#' @param method CellphoneDB analysis mode. One of `"statistical_analysis"`,
#' `"analysis"`, or `"degs_analysis"`.
#' @param cpdb_file_path Path to `cellphonedb.zip`. If `NULL`, the wrapper will try
#' common locations and download it automatically when not found.
#' @param convert_to_human Whether to convert non-human gene symbols to human
#' ortholog symbols before running CellphoneDB.
#' @param gene_id_from_IDtype Input gene identifier type passed to [GeneConvert]
#' when `convert_to_human = TRUE`.
#' @param microenvs Optional CellphoneDB microenvironment specification. Can be a
#' `data.frame` or a file path.
#' @param active_tfs Optional CellphoneDB active TF specification. Can be a
#' `data.frame` or a file path.
#' @param degs Optional DEG table for `method = "degs_analysis"`. Can be a
#' `data.frame` or a file path.
#' @param score_interactions Whether to compute CellphoneDB interaction scores.
#' @param threshold Minimum fraction of cells expressing a gene for the gene to be
#' considered.
#' @param pvalue P-value threshold used by CellphoneDB.
#' @param iterations Number of permutations for statistical analysis.
#' @param cores Number of cores used by CellphoneDB.
#' @param separator Separator used by CellphoneDB for sender/receiver pair columns.
#' @param debug Whether to ask CellphoneDB to save intermediate debug tables.
#' @param debug_seed Random seed used by CellphoneDB.
#' @param result_precision Result precision used by CellphoneDB.
#' @param subsampling Whether to enable CellphoneDB subsampling.
#' @param subsampling_log Whether to log-transform during CellphoneDB subsampling.
#' @param subsampling_num_pc Number of PCs used during CellphoneDB subsampling.
#' @param subsampling_num_cells Number of cells retained during CellphoneDB subsampling.
#' @param output_path Optional directory to keep CellphoneDB output files.
#' @param output_suffix Optional output suffix passed to CellphoneDB.
#' @param keep_output Whether to keep temporary output files when `output_path` is
#' not supplied.
#'
#' @return A Seurat object with results stored in `srt@tools[["CellphoneDB"]]`.
#' @export
RunCellphoneDB <- function(
  srt,
  group.by,
  species = c("Homo_sapiens", "Mus_musculus"),
  assay = NULL,
  layer = "data",
  counts_data = c("hgnc_symbol", "ensembl", "gene_name"),
  method = c("statistical_analysis", "analysis", "degs_analysis"),
  cpdb_file_path = NULL,
  convert_to_human = TRUE,
  gene_id_from_IDtype = "symbol",
  microenvs = NULL,
  active_tfs = NULL,
  degs = NULL,
  score_interactions = TRUE,
  threshold = 0.1,
  pvalue = 0.05,
  iterations = 1000,
  cores = 1,
  separator = "|",
  debug = FALSE,
  debug_seed = 42,
  result_precision = 3,
  subsampling = FALSE,
  subsampling_log = FALSE,
  subsampling_num_pc = 100,
  subsampling_num_cells = 1000,
  output_path = NULL,
  output_suffix = NULL,
  keep_output = FALSE,
  verbose = TRUE
) {
  PrepareEnv(modules = "cellphonedb")
  check_python(c("cellphonedb==5.0.1", "scanpy"), verbose = verbose)

  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  if (!group.by %in% colnames(srt[[]])) {
    log_message(
      "{.arg group.by} must be a valid metadata column in {.cls Seurat}",
      message_type = "error"
    )
  }

  species <- match.arg(species)
  counts_data <- match.arg(counts_data)
  method <- match.arg(method)
  assay <- assay %||% DefaultAssay(srt)
  if (isTRUE(convert_to_human) && species != "Homo_sapiens") {
    counts_data <- "hgnc_symbol"
  }

  log_message(
    "Running {.pkg CellphoneDB} ({.val {method}})...",
    verbose = verbose
  )

  adata <- build_cpdb_adata(
    srt = srt,
    assay = assay,
    layer = layer,
    group.by = group.by,
    species = species,
    convert_to_human = convert_to_human,
    gene_id_from_IDtype = gene_id_from_IDtype,
    verbose = verbose
  )

  functions <- reticulate::import_from_path(
    "functions",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )

  args <- mget(names(formals()))
  args <- lapply(args, function(x) {
    if (is.numeric(x)) {
      ifelse(grepl("\\.", as.character(x)), as.double(x), as.integer(x))
    } else {
      x
    }
  })
  call_envir <- parent.frame(1)
  args <- lapply(args, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call_envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call_envir)
    } else {
      arg
    }
  })
  args[["adata"]] <- adata
  args[["celltype_key"]] <- group.by
  args[["threads"]] <- args[["cores"]]
  args <- args[c(
    "adata", "celltype_key", "method", "cpdb_file_path", "counts_data",
    "microenvs", "active_tfs", "degs", "score_interactions", "iterations",
    "threshold", "pvalue", "threads", "separator", "debug", "debug_seed",
    "result_precision", "subsampling", "subsampling_log", "subsampling_num_pc",
    "subsampling_num_cells", "output_path", "output_suffix", "keep_output",
    "verbose"
  )]

  res <- do.call(functions$CellphoneDB, args)
  res <- tryCatch(check_python_element(res), error = function(e) res)

  bundle <- standardize_cellphonedb_result(
    res = res,
    group.by = group.by,
    species = species,
    assay = assay,
    layer = layer,
    method = method,
    separator = separator,
    converted_to_human = isTRUE(convert_to_human) && species != "Homo_sapiens"
  )

  if (is.null(bundle$long_table) || nrow(bundle$long_table) == 0L) {
    log_message(
      "{.pkg CellphoneDB} returned no detectable interactions for the current input and parameters. Empty result tables were stored in {.code srt@tools[['CellphoneDB']]}",
      message_type = "warning",
      verbose = verbose
    )
  }

  srt@tools[["CellphoneDB"]] <- bundle
  srt <- ccc_update_unified_bundle(
    srt = srt,
    method = "CellphoneDB",
    bundle = bundle
  )

  log_message(
    "{.pkg CellphoneDB} analysis completed",
    message_type = "success",
    verbose = verbose
  )
  srt
}

build_cpdb_adata <- function(
  srt,
  assay,
  layer,
  group.by,
  species,
  convert_to_human = FALSE,
  gene_id_from_IDtype = "symbol",
  verbose = TRUE
) {
  mat <- GetAssayData5(srt, assay = assay, layer = layer)
  obs <- srt[[]]
  obs[[group.by]] <- as.character(obs[[group.by]])
  obs <- as.data.frame(obs)

  converted <- FALSE
  if (isTRUE(convert_to_human) && !identical(species, "Homo_sapiens")) {
    log_message(
      "Converting {.val {species}} genes to human ortholog symbols for {.pkg CellphoneDB}",
      verbose = verbose
    )
    mat <- ConvertHomologs(
      object = mat,
      geneID_from_IDtype = gene_id_from_IDtype,
      geneID_to_IDtype = "symbol",
      species_from = species,
      species_to = "Homo_sapiens",
      verbose = verbose
    )
    converted <- TRUE
  }

  adata <- matrix_to_adata(mat = mat, obs = obs)
  adata$uns[["scop_cellphonedb"]] <- list(
    group_by = group.by,
    converted_to_human = converted,
    source_species = species
  )
  adata
}

matrix_to_adata <- function(mat, obs) {
  sc <- reticulate::import("scanpy", convert = FALSE)

  if (!inherits(mat, "Matrix")) {
    mat <- methods::as(mat, "dgCMatrix")
  }

  obs <- as.data.frame(obs)
  for (i in seq_len(ncol(obs))) {
    if (is.logical(obs[, i])) {
      obs[, i] <- factor(as.character(obs[, i]), levels = c("TRUE", "FALSE"))
    }
  }

  adata <- sc$AnnData(
    X = reticulate::r_to_py(Matrix::t(mat)),
    obs = obs,
    var = data.frame(features = rownames(mat), row.names = rownames(mat))
  )
  adata$var_names <- rownames(mat)
  adata
}

standardize_cellphonedb_result <- function(
  res,
  group.by,
  species,
  assay,
  layer,
  method,
  separator = "|",
  converted_to_human = FALSE
) {
  results <- res[["results"]] %||% res
  tables <- lapply(results, function(x) {
    if (is.data.frame(x)) {
      x
    } else if (inherits(x, "python.builtin.object")) {
      py_to_r2(x)
    } else {
      x
    }
  })

  long_table <- cpdb_tables_to_long(tables = tables, separator = separator)
  pair_table <- aggregate_ccc_long(long_table)

  list(
    method = "CellphoneDB",
    results = tables,
    long_table = long_table,
    pair_table = pair_table,
    parameters = list(
      group.by = group.by,
      species = species,
      assay = assay,
      layer = layer,
      method = method,
      separator = separator,
      converted_to_human = converted_to_human
    ),
    output_path = res[["output_path"]] %||% NULL,
    cpdb_file_path = res[["cpdb_file_path"]] %||% NULL
  )
}

cpdb_tables_to_long <- function(tables, separator = "|") {
  means_df <- tables[["means"]] %||% tables[["means_result"]]
  if (is.null(means_df) || !is.data.frame(means_df)) {
    return(data.frame())
  }

  pvalues_df <- tables[["pvalues"]] %||% tables[["pvalues_result"]]
  scores_df <- tables[["interaction_scores"]] %||%
    tables[["interaction_scores_result"]]
  sig_df <- tables[["significant_means"]] %||%
    tables[["significant_means_result"]]

  pair_cols <- names(means_df)[vapply(
    names(means_df),
    function(x) grepl(separator, x, fixed = TRUE),
    logical(1)
  )]
  if (length(pair_cols) == 0L) {
    return(data.frame())
  }

  info_cols <- setdiff(names(means_df), pair_cols)
  n_interactions <- nrow(means_df)
  n_pairs <- length(pair_cols)
  pair_parts <- strsplit(pair_cols, split = separator, fixed = TRUE)
  sender <- vapply(pair_parts, function(x) x[1] %||% NA_character_, character(1))
  receiver <- vapply(pair_parts, function(x) x[2] %||% NA_character_, character(1))
  sender[is.na(sender)] <- pair_cols[is.na(sender)]
  receiver[is.na(receiver)] <- pair_cols[is.na(receiver)]

  table_values <- function(df, default = NA_real_) {
    out <- matrix(default, nrow = n_interactions, ncol = n_pairs)
    if (!is.null(df)) {
      available <- intersect(pair_cols, colnames(df))
      if (length(available) > 0L) {
        out[, match(available, pair_cols)] <- as.matrix(df[, available, drop = FALSE])
      }
    }
    as.vector(out)
  }

  out <- means_df[
    rep(seq_len(n_interactions), times = n_pairs),
    info_cols,
    drop = FALSE
  ]
  out[["sender"]] <- rep(sender, each = n_interactions)
  out[["receiver"]] <- rep(receiver, each = n_interactions)
  out[["pair"]] <- rep(pair_cols, each = n_interactions)
  out[["means"]] <- table_values(means_df)
  out[["score"]] <- out[["means"]]
  out[["interaction_score"]] <- table_values(scores_df)
  out[["pvalue"]] <- table_values(pvalues_df)
  out[["significant_mean"]] <- table_values(sig_df)
  rownames(out) <- NULL

  ligand_col <- intersect(
    c("gene_a", "ligand", "partner_a", "multidata_1_id", "interactor_1"),
    colnames(out)
  )[1]
  receptor_col <- intersect(
    c("gene_b", "receptor", "partner_b", "multidata_2_id", "interactor_2"),
    colnames(out)
  )[1]
  interaction_col <- intersect(
    c("interacting_pair", "interaction_name", "id_cp_interaction"),
    colnames(out)
  )[1]

  out[["ligand"]] <- if (!is.na(ligand_col)) {
    out[[ligand_col]]
  } else {
    NA_character_
  }
  out[["receptor"]] <- if (!is.na(receptor_col)) {
    out[[receptor_col]]
  } else {
    NA_character_
  }
  out[["interaction_name"]] <- if (!is.na(interaction_col)) {
    out[[interaction_col]]
  } else {
    paste(out[["ligand"]], out[["receptor"]], sep = " - ")
  }
  out[["method"]] <- "CellphoneDB"
  out
}

aggregate_ccc_long <- function(df) {
  if (is.null(df) || nrow(df) == 0L) {
    return(data.frame())
  }
  score_col <- if ("score" %in% colnames(df)) {
    "score"
  } else if ("means" %in% colnames(df)) {
    "means"
  } else if ("prob" %in% colnames(df)) {
    "prob"
  } else {
    NULL
  }
  if (is.null(score_col) || !all(c("sender", "receiver") %in% colnames(df))) {
    return(data.frame())
  }

  pair_key <- paste(
    df$sender,
    df$receiver,
    sep = "
"
  )
  score <- df[[score_col]]
  score[is.na(score)] <- 0
  significant <- if ("significant" %in% colnames(df)) {
    suppressWarnings(as.numeric(df$significant))
  } else if ("pvalue" %in% colnames(df)) {
    pvalue <- suppressWarnings(as.numeric(df$pvalue))
    ifelse(is.finite(pvalue), as.numeric(pvalue < 0.05), as.numeric(is.finite(score) & score > 0))
  } else {
    as.numeric(is.finite(score) & score > 0)
  }
  significant[!is.finite(significant)] <- 0

  sum_vec <- tapply(score, pair_key, sum, na.rm = TRUE)
  mean_vec <- tapply(score, pair_key, mean, na.rm = TRUE)
  max_vec <- tapply(score, pair_key, max, na.rm = TRUE)
  count_vec <- tapply(significant, pair_key, sum, na.rm = TRUE)

  keys <- names(sum_vec)
  parts <- strsplit(
    keys,
    split = "
",
    fixed = TRUE
  )
  out <- data.frame(
    sender = vapply(parts, `[`, character(1), 1),
    receiver = vapply(parts, `[`, character(1), 2),
    sum = as.numeric(sum_vec),
    mean = as.numeric(mean_vec),
    max = as.numeric(max_vec),
    count = as.numeric(count_vec),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}
