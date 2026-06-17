#' @title Run cisTarget motif enrichment on a GRN adjacency table
#'
#' @description
#' cisTarget performs motif enrichment analysis on gene regulatory network
#' adjacency tables. For each transcription factor, it identifies target
#' genes whose regulatory regions are enriched for the TF's binding motifs,
#' producing regulons (TF + enriched target gene sets).
#'
#' @md
#' @inheritParams RunSCENIC
#' @param adj A data frame with columns `TF`, `target`, and optionally
#' `importance`, or a path to a TSV adjacency file. This is typically the
#' output of [RunGRNBoost2()], [RunGENIE3()], or [RunSCENIC()] step 1.
#' @param backend cisTarget runtime backend.
#' \describe{
#'   \item{`"python"`}{Uses the pySCENIC/ctxcore Python pipeline. This is the
#'     most tested backend and produces results identical to official SCENIC.}
#'   \item{`"r"`}{Uses the `RcisTarget` Bioconductor package. Requires the
#'     cisTarget ranking databases and motif annotations to be available in
#'     the same format as the Python backend (feather files and motif2tf
#'     table).}
#' }
#' @param ranking_dbs Character vector of cisTarget ranking feather files.
#' @param motif_annotations Motif annotation table path (motif2tf).
#' @param expression_mtx Optional expression matrix path (CSV). When
#' `NULL` and `backend = "python"`, the expression matrix is reconstructed
#' from the unique genes in `adj` as a minimal stub.
#' @param ctx_output Optional output file path for the cisTarget result.
#' @param min_regulon_size Minimum number of target genes per regulon.
#' @param gmt_output Optional output path for the regulon GMT file.
#' @param txt_output Optional output path for the regulon TXT file.
#' @param work_dir Working directory used by Python backend.
#' @param prefix Prefix for output files.
#' @param envname Python environment name (Python backend only).
#' @param conda Conda-compatible executable (Python backend only).
#' @param prepare_env Whether to prepare the Python environment.
#' @param cores Number of workers.
#' @param force Whether to rebuild existing outputs.
#' @param verbose Whether to print progress messages.
#' @param ... Additional arguments passed to the backend.
#'
#' @return A named list with components `regulons` (list of gene vectors),
#' `ctx_file` (path to raw cisTarget output), `gmt_file`, and `txt_file`.
#' @export
#'
#' @examples
#' \dontrun{
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(pancreas_sub)
#'
#' # First run GRN inference
#' grn <- RunGRNBoost2(
#'   pancreas_sub,
#'   regulators = c("Neurod1", "Arx", "Pax6"),
#'   backend = "cpp"
#' )
#'
#' # Then run cisTarget (Python backend)
#' regulons <- RunCisTarget(
#'   grn,
#'   species = "Mus_musculus",
#'   backend = "python"
#' )
#' }
RunCisTarget <- function(
  adj,
  species = c("Homo_sapiens", "Mus_musculus", "Drosophila_melanogaster"),
  backend = c("python", "r"),
  ranking_dbs = NULL,
  motif_annotations = NULL,
  expression_mtx = NULL,
  ctx_output = NULL,
  min_regulon_size = 10,
  gmt_output = NULL,
  txt_output = NULL,
  work_dir = NULL,
  prefix = "cisTarget",
  data_dir = NULL,
  envname = NULL,
  conda = "auto",
  prepare_env = TRUE,
  cores = 1,
  force = FALSE,
  verbose = TRUE,
  ...
) {
  backend <- match.arg(backend)
  species <- match.arg(species)

  work_dir <- work_dir %||% "."
  if (!dir.exists(work_dir)) {
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Resolve ranking databases and motif annotations from species
  reference_data <- scenic_reference(
    species = species,
    data_dir = data_dir,
    ranking_dbs = ranking_dbs,
    motif_annotations = motif_annotations,
    verbose = verbose
  )
  ranking_dbs <- reference_data[["ranking_dbs"]]
  motif_annotations <- reference_data[["motif_annotations"]]

  # Normalize adjacency input
  if (is.character(adj) && length(adj) == 1 && file.exists(adj)) {
    adj_file <- adj
    adj <- utils::read.table(
      adj_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE
    )
  } else if (is.data.frame(adj)) {
    adj_file <- file.path(work_dir, paste0(prefix, "_adj.tsv"))
    utils::write.table(
      adj, adj_file, sep = "\t", row.names = FALSE, quote = FALSE
    )
  } else {
    log_message(
      "{.arg adj} must be a data.frame or a path to a TSV file",
      message_type = "error"
    )
  }

  required_cols <- c("TF", "target")
  missing_cols <- setdiff(required_cols, colnames(adj))
  if (length(missing_cols) > 0) {
    log_message(
      "{.arg adj} missing required columns: {.val {missing_cols}}",
      message_type = "error"
    )
  }

  # Output files
  ctx_output <- ctx_output %||% file.path(
    work_dir, paste0(prefix, "_ctx.csv")
  )
  gmt_output <- gmt_output %||% file.path(
    work_dir, paste0(prefix, "_regulons.gmt")
  )
  txt_output <- txt_output %||% file.path(
    work_dir, paste0(prefix, "_regulons.txt")
  )

  switch(
    backend,
    python = cisTarget_python(
      adj = adj,
      adj_file = adj_file,
      ranking_dbs = ranking_dbs,
      motif_annotations = motif_annotations,
      expression_mtx = expression_mtx,
      ctx_output = ctx_output,
      gmt_output = gmt_output,
      txt_output = txt_output,
      min_regulon_size = min_regulon_size,
      work_dir = work_dir,
      prefix = prefix,
      envname = envname,
      conda = conda,
      prepare_env = prepare_env,
      cores = cores,
      force = force,
      verbose = verbose,
      ...
    ),
    r = cisTarget_r(
      adj = adj,
      ranking_dbs = ranking_dbs,
      motif_annotations = motif_annotations,
      ctx_output = ctx_output,
      gmt_output = gmt_output,
      txt_output = txt_output,
      min_regulon_size = min_regulon_size,
      force = force,
      verbose = verbose,
      ...
    )
  )
}

# ── Python backend ─────────────────────────────────────────────────────────

cisTarget_python <- function(
  adj,
  adj_file,
  ranking_dbs,
  motif_annotations,
  expression_mtx,
  ctx_output,
  gmt_output,
  txt_output,
  min_regulon_size,
  work_dir,
  prefix,
  envname,
  conda,
  prepare_env,
  cores,
  force,
  verbose,
  ...
) {
  check_r("reticulate", verbose = FALSE)

  if (isTRUE(prepare_env)) {
    PrepareEnv(envname = envname, modules = "scenic")
  }

  functions <- reticulate::import_from_path(
    "functions",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )

  # Build expression matrix CSV if not provided
  if (is.null(expression_mtx)) {
    expr_csv <- file.path(work_dir, paste0(prefix, "_expr.csv"))
    # Collect all genes referenced anywhere in the adjacency
    all_genes <- unique(unlist(lapply(adj, function(col) {
      if (is.character(col) || is.factor(col)) as.character(col) else NULL
    })))
    all_genes <- all_genes[nzchar(all_genes) & !is.na(all_genes)]
    # Build a minimal 3-cell stub so pySCENIC correlation works
    n_cells <- min(3L, length(all_genes))
    stub <- as.data.frame(matrix(
      1.0,
      nrow = length(all_genes),
      ncol = n_cells,
      dimnames = list(all_genes, paste0("Cell", seq_len(n_cells)))
    ))
    utils::write.csv(stub, expr_csv)
    expression_mtx <- expr_csv
  }

  # Convert motif_annotations to string (handles data.frame input)
  if (is.data.frame(motif_annotations)) {
    motif_file <- file.path(work_dir, paste0(prefix, "_motif2tf.tbl"))
    utils::write.table(
      motif_annotations, motif_file,
      sep = "\t", row.names = FALSE, quote = FALSE
    )
    motif_annotations <- motif_file
  }

  log_message(
    "Running cisTarget (Python ctxcore) on {.val {nrow(adj)}} edges...",
    verbose = verbose
  )

  functions$RunSCENICCtx(
    expression_mtx = expression_mtx,
    ranking_dbs = as.list(ranking_dbs),
    motif_annotations = motif_annotations,
    adj_output = adj_file,
    ctx_output = ctx_output,
    cores = as.integer(cores),
    force = isTRUE(force),
    verbose = isTRUE(verbose)
  )

  # Convert ctx output to regulon files
  functions$SCENICRegulonsToFiles(
    ctx_output,
    gmt_output,
    txt_output,
    min_regulon_size = as.integer(min_regulon_size)
  )

  regulons <- read_regulons_from_gmt(gmt_output, min_regulon_size)

  log_message(
    "{.fn RunCisTarget} (python) produced {.val {length(regulons)}} regulons",
    message_type = "success",
    verbose = verbose
  )

  list(
    regulons = regulons,
    ctx_file = ctx_output,
    gmt_file = gmt_output,
    txt_file = txt_output
  )
}

# ── R backend (RcisTarget) ─────────────────────────────────────────────────

cisTarget_r <- function(
  adj,
  ranking_dbs,
  motif_annotations,
  ctx_output,
  gmt_output,
  txt_output,
  min_regulon_size,
  force,
  verbose,
  ...
) {
  check_r("RcisTarget", verbose = FALSE)
  importRankings <- get_namespace_fun("RcisTarget", "importRankings")
  cisTarget <- get_namespace_fun("RcisTarget", "cisTarget")

  # RcisTarget requires the ranking databases in feather format
  # and the motif annotations as a data.table
  log_message(
    "Running cisTarget (RcisTarget) on {.val {nrow(adj)}} edges...",
    verbose = verbose
  )

  motif_annotations_dt <- data.table::fread(motif_annotations)

  # Import ranking databases for RcisTarget
  db_paths <- ranking_dbs
  if (length(db_paths) > 1) {
    log_message(
      "RcisTarget does not natively support multiple ranking databases. Using first: {.file {db_paths[1]}}",
      message_type = "warning",
      verbose = verbose
    )
  }

  # Build gene list per TF from adjacency
  tf_targets <- split(adj[["target"]], adj[["TF"]])
  tf_targets <- lapply(tf_targets, unique)

  # Run RcisTarget
  motif_rankings <- tryCatch(
    importRankings(db_paths[1]),
    error = function(e) {
      log_message(
        "Failed to import rankings: {.val {conditionMessage(e)}}",
        message_type = "error"
      )
    }
  )

  motif_enrichment <- cisTarget(
    tf_targets,
    motif_rankings,
    motifAnnot = motif_annotations_dt
  )

  # Build regulons from RcisTarget output
  regulons <- list()
  if (is.data.frame(motif_enrichment)) {
    for (tf_name in unique(motif_enrichment[["geneSet"]])) {
      enriched <- motif_enrichment[motif_enrichment[["geneSet"]] == tf_name, , drop = FALSE]
      if (!nrow(enriched)) {
        next
      }
      top_row <- enriched[which.max(enriched[["NES"]]), , drop = FALSE]
      target_genes <- unique(unlist(strsplit(top_row[["enrichedGenes"]], ";")))
      target_genes <- target_genes[nzchar(target_genes)]
      if (length(target_genes) >= min_regulon_size) {
        regulons[[tf_name]] <- target_genes
      }
    }
  } else {
    for (tf_name in names(motif_enrichment)) {
      enriched <- motif_enrichment[[tf_name]]
      if (is.null(enriched) || nrow(enriched) == 0) {
        next
      }
      # Take top motif per TF
      top_row <- enriched[which.max(enriched[["NES"]]), ]
      target_genes <- unique(unlist(strsplit(top_row[["enrichedGenes"]], ";")))
      target_genes <- target_genes[nzchar(target_genes)]
      if (length(target_genes) >= min_regulon_size) {
        regulons[[tf_name]] <- target_genes
      }
    }
  }

  # Write GMT
  write_regulons_to_gmt(regulons, gmt_output)
  write_regulons_to_txt(regulons, txt_output)

  log_message(
    "{.fn RunCisTarget} (r) produced {.val {length(regulons)}} regulons",
    message_type = "success",
    verbose = verbose
  )

  list(
    regulons = regulons,
    ctx_file = ctx_output,
    gmt_file = gmt_output,
    txt_file = txt_output
  )
}

# ── Helpers ────────────────────────────────────────────────────────────────

read_regulons_from_gmt <- function(gmt_file, min_regulon_size = 10) {
  if (!file.exists(gmt_file)) {
    return(list())
  }
  lines <- readLines(gmt_file, warn = FALSE)
  regulons <- list()
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) < 3) next
    tf_name <- gsub("\\(\\d+g\\)$", "", parts[1])
    genes <- parts[-(1:2)]
    genes <- genes[nzchar(genes)]
    if (length(genes) >= min_regulon_size) {
      regulons[[tf_name]] <- genes
    }
  }
  regulons
}

write_regulons_to_gmt <- function(regulons, gmt_file) {
  lines <- character(length(regulons))
  for (i in seq_along(regulons)) {
    tf <- names(regulons)[i]
    genes <- regulons[[i]]
    lines[i] <- paste0(
      tf, "(", length(genes), "g)", "\tna\t",
      paste(genes, collapse = "\t")
    )
  }
  writeLines(lines, gmt_file)
}

write_regulons_to_txt <- function(regulons, txt_file) {
  lines <- character(length(regulons))
  for (i in seq_along(regulons)) {
    tf <- names(regulons)[i]
    genes <- regulons[[i]]
    lines[i] <- paste0(
      tf, "(", length(genes), "g)", "\tna\t",
      paste(genes, collapse = ",")
    )
  }
  writeLines(lines, txt_file)
}
