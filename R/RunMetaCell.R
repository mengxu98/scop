#' @title Run metacell partitioning for single-cell data
#'
#' @md
#' @inheritParams standard_scop
#' @inheritParams thisutils::log_message
#' @param method Metacell construction method. One of `"supercell"`,
#' `"seacells"`, or `"metacell"`.
#' @param assay Assay to use for metacell construction.
#' @param layer Assay layer used to extract the count matrix.
#' @param gamma Metacell granularity parameter. For SuperCell, larger values
#' produce fewer metacells (typical range 10--50). For MetaCell, this is the
#' K parameter controlling the number of metacells. For SEACells, the
#' comparable parameter is passed via `...` (e.g. `n_metacells`).
#' @param group.by Optional metadata column used to build metacells within
#' each group independently (e.g. by sample or cell type), preventing
#' metacells from crossing group boundaries.
#' @param envname Python environment name (SEACells only). Passed to
#' `reticulate::use_condaenv()` when `method = "seacells"`.
#' @param conda Conda executable path (SEACells only).
#' @param prefix Prefix for metadata columns written to `srt`.
#' @param tool_name Name of the `srt@tools` entry.
#' @param ... Additional arguments passed to the underlying metacell method.
#'
#' @return A metacell-level `Seurat` object. The original single-cell Seurat
#' is stored in `@misc[["original_srt"]]` and the cell-to-metacell membership
#' vector in `@misc[["cell_membership"]]`. The returned object can be passed
#' directly to any scop function (`standard_scop()`, `CellDimPlot()`, etc.).
#' @export
#'
#' @references
#' Baran, Y. et al. (2019). MetaCell: analysis of single-cell RNA-seq data
#' using K-nn graph partitions. \emph{Genome Biology}.
#'
#' Bilous, M. et al. (2022). SuperCell: a versatile tool for single-cell
#' data analysis. \emph{Genome Biology}.
#'
#' Persad, S. et al. (2023). SEACells infers transcriptional and epigenomic
#' cellular states from single-cell genomics data. \emph{Nature Biotechnology}.
#'
#' @examples
#' data(pancreas_sub)
#' pancreas_sub <- standard_scop(
#'   pancreas_sub,
#'   nHVF = 500,
#'   linear_reduction_dims = 20,
#'   linear_reduction_dims_use = 1:20,
#'   nonlinear_reduction_dims = 2,
#'   verbose = FALSE
#' )
#'
#' mc1 <- RunMetaCell(
#'   pancreas_sub,
#'   method = "supercell",
#'   gamma = 20
#' )
#'
#' MetaCellPlot(mc1, group.by = "CellType")
#'
#' mc2 <- RunMetaCell(
#'   pancreas_sub,
#'   method = "metacell",
#'   gamma = 20
#' )
#'
#' MetaCellPlot(mc2, group.by = "CellType")
RunMetaCell <- function(
  srt,
  method = c("supercell", "seacells", "metacell"),
  assay = NULL,
  layer = "counts",
  gamma = 20,
  group.by = NULL,
  envname = NULL,
  conda = "auto",
  prefix = "Metacell",
  tool_name = "Metacell",
  verbose = TRUE,
  ...
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }
  method <- match.arg(method)

  for (name in c("prefix", "tool_name")) {
    value <- get(name)
    if (!is.character(value) || length(value) != 1L || !nzchar(value)) {
      log_message(
        "{.arg {name}} must be a single non-empty string",
        message_type = "error"
      )
    }
  }
  if (!is.numeric(gamma) || length(gamma) != 1L || is.na(gamma) || gamma <= 0) {
    log_message(
      "{.arg gamma} must be a single positive number",
      message_type = "error"
    )
  }

  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  counts <- GetAssayData5(srt, assay = assay, layer = layer)
  if (nrow(counts) == 0 || ncol(counts) == 0) {
    log_message(
      "No expression values available for metacell construction",
      message_type = "error"
    )
  }

  group_df <- NULL
  if (!is.null(group.by)) {
    missing_cols <- setdiff(group.by, colnames(srt@meta.data))
    if (length(missing_cols) > 0) {
      log_message(
        "{.arg group.by} columns not found: {.val {missing_cols}}",
        message_type = "error"
      )
    }
    group_df <- srt@meta.data[, group.by, drop = FALSE]
    if (any(!stats::complete.cases(group_df))) {
      log_message(
        "{.arg group.by} contains NA values; fill or remove missing annotations",
        message_type = "error"
      )
    }
  }

  mc_result <- switch(
    method,
    supercell = metacell_supercell(
      counts = counts,
      gamma = gamma,
      group_df = group_df,
      verbose = verbose,
      ...
    ),
    seacells = metacell_seacells(
      srt = srt,
      counts = counts,
      gamma = gamma,
      group_df = group_df,
      assay = assay,
      layer = layer,
      envname = envname,
      conda = conda,
      verbose = verbose,
      ...
    ),
    metacell = metacell_metacell(
      counts = counts,
      gamma = gamma,
      group_df = group_df,
      verbose = verbose,
      ...
    ),
    log_message(
      "Method {.val {method}} is not yet implemented",
      message_type = "error"
    )
  )

  raw_membership <- mc_result[["membership"]]
  mc_sc <- mc_result[["SC"]]
  mc_counts <- mc_result[["counts"]]

  raw_membership <- as.character(raw_membership)
  names(raw_membership) <- colnames(srt)

  unique_raw <- unique(raw_membership)
  mc_labels <- paste0("MC", seq_along(unique_raw))
  names(mc_labels) <- unique_raw

  membership <- mc_labels[raw_membership]
  names(membership) <- colnames(srt)
  metacell_sizes <- table(membership)

  log_message(
    "{.fn RunMetaCell} ({.val {method}}) built {.val {length(metacell_sizes)}} metacells from {.val {ncol(srt)}} cells",
    verbose = verbose
  )
  log_message(
    "Metacell size summary: min {.val {min(metacell_sizes)}}, median {.val {stats::median(metacell_sizes)}}, mean {.val {round(mean(metacell_sizes), 2)}}, max {.val {max(metacell_sizes)}} cells",
    verbose = verbose
  )

  raw_to_mc <- stats::setNames(mc_labels, unique_raw)
  mc_count_ids <- colnames(mc_counts)
  if (
    is.null(mc_count_ids) ||
      any(is.na(mc_count_ids)) ||
      any(!nzchar(mc_count_ids))
  ) {
    if (ncol(mc_counts) != length(unique_raw)) {
      log_message(
        "Metacell count matrix has {.val {ncol(mc_counts)}} columns but {.val {length(unique_raw)}} metacell memberships",
        message_type = "error"
      )
    }
    mc_count_ids <- unique_raw
  } else {
    mc_count_ids <- as.character(mc_count_ids)
  }
  if (!all(mc_count_ids %in% names(raw_to_mc))) {
    if (ncol(mc_counts) == length(unique_raw)) {
      mc_count_ids <- unique_raw
    } else {
      log_message(
        "Metacell count column names do not match membership IDs",
        message_type = "error"
      )
    }
  }
  colnames(mc_counts) <- unname(raw_to_mc[mc_count_ids])

  mc_seurat <- Seurat::CreateSeuratObject(
    counts = mc_counts,
    assay = assay
  )

  mc_seurat[["metacell_size"]] <- as.integer(metacell_sizes[colnames(
    mc_seurat
  )])
  mc_seurat[["metacell_method"]] <- method

  cell_meta <- srt@meta.data
  transfer_cols <- setdiff(
    names(cell_meta),
    c("nCount_RNA", "nFeature_RNA")
  )
  for (col_name in transfer_cols) {
    vals <- cell_meta[[col_name]]
    names(vals) <- colnames(srt)
    split_vals <- split(vals, membership)
    if (is.numeric(vals)) {
      mapped <- vapply(
        split_vals,
        function(x) mean(x, na.rm = TRUE),
        numeric(1)
      )
      mapped[is.nan(mapped)] <- NA_real_
      mc_seurat@meta.data[[col_name]] <- mapped[colnames(mc_seurat)]
    } else {
      mapped <- vapply(
        split_vals,
        function(x) {
          x <- x[!is.na(x)]
          if (length(x) == 0) {
            return(NA_character_)
          }
          names(sort(table(x), decreasing = TRUE))[1]
        },
        character(1)
      )
      mapped <- mapped[colnames(mc_seurat)]
      if (is.factor(vals)) {
        mc_seurat@meta.data[[col_name]] <- factor(mapped, levels = levels(vals))
      } else if (is.logical(vals)) {
        mc_seurat@meta.data[[col_name]] <- mapped == "TRUE"
      } else {
        mapped[is.na(mapped)] <- "Unknown"
        mc_seurat@meta.data[[col_name]] <- mapped
      }
    }
  }

  mc_seurat@misc[["original_srt"]] <- srt
  mc_seurat@misc[["cell_membership"]] <- membership

  mc_seurat@tools[[tool_name]] <- list(
    method = method,
    membership = membership,
    metacell_sizes = as.integer(metacell_sizes[colnames(mc_seurat)]),
    SC = mc_sc,
    counts = mc_counts,
    parameters = list(
      assay = assay,
      layer = layer,
      gamma = gamma,
      group.by = group.by,
      prefix = prefix,
      tool_name = tool_name
    )
  )

  log_message(
    "{.fn RunMetaCell} returned metacell Seurat with {.val {ncol(mc_seurat)}} metacells. Original cells in {.code @misc[[\"original_srt\"]]}",
    message_type = "success",
    verbose = verbose
  )

  mc_seurat
}

metacell_supercell <- function(
  counts,
  gamma,
  k.knn = 5,
  n.var.genes = 1000,
  group_df = NULL,
  verbose = TRUE,
  ...
) {
  check_r("GfellerLab/SuperCell", verbose = FALSE)

  counts <- methods::as(counts, "CsparseMatrix")
  ge <- lognorm_counts(counts)

  if (is.null(group_df)) {
    log_message(
      "Running {.pkg SuperCell} with gamma = {.val {gamma}}, k.knn = {.val {k.knn}} on {.val {ncol(counts)}} cells",
      verbose = verbose
    )

    SC <- get_namespace_fun("SuperCell", "SCimplify")(
      ge,
      k.knn = k.knn,
      gamma = gamma,
      n.var.genes = n.var.genes,
      ...
    )

    SC.GE <- get_namespace_fun("SuperCell", "supercell_GE")(
      ge,
      SC$membership
    )

    return(list(
      membership = as.character(SC$membership),
      SC = SC,
      counts = SC.GE
    ))
  }

  group_interaction <- do.call(
    interaction,
    c(group_df, list(drop = TRUE, sep = "."))
  )
  group_levels <- levels(group_interaction)

  log_message(
    "Running {.pkg SuperCell} within {.val {length(group_levels)}} group{?s}",
    verbose = verbose
  )

  membership <- character(ncol(counts))

  for (grp in group_levels) {
    grp_cells <- which(group_interaction == grp)
    if (length(grp_cells) == 0) {
      next
    }

    grp_ge <- ge[, grp_cells, drop = FALSE]
    log_message(
      "  Group {.val {grp}}: {.val {length(grp_cells)}} cells",
      verbose = verbose
    )

    SC_grp <- get_namespace_fun("SuperCell", "SCimplify")(
      grp_ge,
      k.knn = k.knn,
      gamma = gamma,
      n.var.genes = n.var.genes,
      ...
    )

    grp_membership <- paste0(grp, ".", SC_grp[["membership"]])
    membership[grp_cells] <- grp_membership
  }

  SC.GE <- get_namespace_fun("SuperCell", "supercell_GE")(
    ge,
    membership
  )

  list(
    membership = membership,
    SC = NULL,
    counts = SC.GE
  )
}

metacell_seacells <- function(
  srt,
  counts,
  gamma,
  group_df = NULL,
  assay = NULL,
  layer = "counts",
  envname = NULL,
  conda = "auto",
  verbose = TRUE,
  ...
) {
  check_r(
    c("reticulate", "Matrix"),
    verbose = FALSE
  )
  assay <- assay %||% SeuratObject::DefaultAssay(srt)

  PrepareEnv(modules = c("scanpy", "seacells"))

  functions <- reticulate::import_from_path(
    "functions",
    path = system.file("python", package = "scop", mustWork = TRUE),
    convert = TRUE
  )

  pca_reduction <- tryCatch(
    DefaultReduction(srt, pattern = "pca"),
    error = function(e) NULL
  )
  if (is.null(pca_reduction)) {
    log_message(
      "{.pkg SEACells} requires PCA. Run {.fn standard_scop} or {.fn RunPCA} first to compute PCA on the original Seurat.",
      message_type = "error"
    )
  }

  adata_full <- srt_to_adata(
    srt = srt,
    assay_x = assay,
    layer_x = layer,
    reductions = pca_reduction,
    verbose = verbose
  )

  pca_ok <- tryCatch(
    {
      adata_full$obsm[pca_reduction]
      TRUE
    },
    error = function(e) FALSE
  )
  if (!isTRUE(pca_ok)) {
    log_message(
      "{.arg {pca_reduction}} not found in AnnData obsm after conversion",
      message_type = "error"
    )
  }

  extra_args <- list(...)

  run_seacells <- function(adata, n_meta, cell_names, log_prefix = "") {
    n_cells <- reticulate::py_to_r(adata$n_obs)
    if (n_cells < 3) {
      mem <- stats::setNames(rep("1", n_cells), cell_names)
      return(list(membership = mem))
    }
    n_meta_use <- max(2L, min(n_meta, as.integer(n_cells / 2)))

    py_args <- list(
      adata = adata,
      n_SEACells = as.integer(n_meta_use),
      build_kernel_on = pca_reduction,
      n_waypoint_eigs = as.integer(min(10L, n_meta_use - 1L)),
      convergence_epsilon = 1e-5,
      min_iter = 10L,
      max_iter = 100L,
      verbose = isTRUE(verbose)
    )
    for (nm in names(extra_args)) {
      py_args[[nm]] <- extra_args[[nm]]
    }

    log_message(
      "{log_prefix}SEACells: {.val {n_meta_use}} archetypes on {.val {n_cells}} cells",
      verbose = verbose
    )

    py_result <- do.call(functions$RunSEACells, py_args)

    membership_py <- py_result[["membership"]]
    n_meta_py <- py_result[["n_metacells"]]
    if (is.null(membership_py) || length(membership_py) == 0) {
      log_message(
        "{log_prefix}SEACells returned no assignments, falling back to random",
        message_type = "warning",
        verbose = verbose
      )
      mem <- stats::setNames(
        as.character(sample(seq_len(n_meta_use), n_cells, replace = TRUE)),
        cell_names
      )
      return(list(membership = mem))
    }

    membership <- as.character(membership_py)
    names(membership) <- cell_names
    list(membership = membership)
  }

  aggregate_counts <- function(counts, membership) {
    unique_mc <- unique(membership)
    mc_counts <- Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      dims = c(nrow(counts), length(unique_mc))
    )
    rownames(mc_counts) <- rownames(counts)
    colnames(mc_counts) <- unique_mc
    for (mc_name in unique_mc) {
      mc_cells <- which(membership == mc_name)
      mc_counts[, mc_name] <- Matrix::rowSums(
        counts[, mc_cells, drop = FALSE]
      )
    }
    mc_counts
  }

  if (!is.null(group_df)) {
    group_interaction <- do.call(
      interaction,
      c(group_df, list(drop = TRUE, sep = "."))
    )
    group_levels <- levels(group_interaction)
    n_metacells <- as.integer(gamma)

    log_message(
      "Running {.pkg SEACells} within {.val {length(group_levels)}} group{?s}",
      verbose = verbose
    )

    membership <- character(ncol(counts))
    names(membership) <- colnames(counts)

    for (grp in group_levels) {
      grp_cells <- which(group_interaction == grp)
      if (length(grp_cells) < 3) {
        membership[grp_cells] <- paste0(grp, ".1")
        next
      }

      n_meta_grp <- max(
        2L,
        as.integer(ceiling(
          n_metacells * length(grp_cells) / ncol(counts)
        ))
      )

      grp_indices <- reticulate::r_to_py(as.integer(grp_cells - 1L))
      grp_adata <- adata_full[grp_indices, ]

      result <- run_seacells(
        adata = grp_adata,
        n_meta = n_meta_grp,
        cell_names = colnames(counts)[grp_cells],
        log_prefix = paste0("  Group ", grp, ": ")
      )

      grp_membership <- paste0(grp, ".", result[["membership"]])
      membership[grp_cells] <- grp_membership
    }

    mc_counts <- aggregate_counts(counts, membership)

    list(
      membership = membership,
      SC = NULL,
      counts = mc_counts
    )
  } else {
    n_cells <- ncol(counts)
    n_metacells <- as.integer(gamma)
    if (n_metacells < 2) {
      n_metacells <- max(2L, as.integer(ceiling(n_cells / 50)))
    }
    if (n_metacells >= n_cells) {
      n_metacells <- max(2L, as.integer(n_cells / 2))
    }

    result <- run_seacells(
      adata = adata_full,
      n_meta = n_metacells,
      cell_names = colnames(counts),
      log_prefix = ""
    )

    membership <- result[["membership"]]
    mc_counts <- aggregate_counts(counts, membership)

    list(
      membership = membership,
      SC = NULL,
      counts = mc_counts
    )
  }
}

metacell_metacell <- function(
  counts,
  gamma,
  group_df = NULL,
  verbose = TRUE,
  ...
) {
  check_r(
    c("Matrix", "RANN"),
    verbose = FALSE
  )
  counts <- methods::as(counts, "CsparseMatrix")
  ge <- lognorm_counts(counts)

  extra_args <- list(...)
  k_knn <- extra_args[["k_knn"]] %||% max(5L, as.integer(gamma))
  n_pcs <- extra_args[["n_pcs"]] %||% min(50L, ncol(counts) - 1L)
  n_var_genes <- extra_args[["n_var_genes"]] %||% min(2000L, nrow(counts))

  gene_vars <- apply(ge, 1, stats::var)
  gene_vars <- sort(gene_vars, decreasing = TRUE)
  top_genes <- names(gene_vars)[seq_len(min(n_var_genes, length(gene_vars)))]
  ge_sub <- ge[top_genes, , drop = FALSE]

  ge_scaled <- Matrix::t(ge_sub)
  if (ncol(ge_scaled) > 1000) {
    pca_res <- tryCatch(
      irlba::irlba(ge_scaled, nv = min(n_pcs, ncol(ge_scaled))),
      error = function(e) {
        log_message(
          "irlba PCA failed: {.val {conditionMessage(e)}}. Falling back to prcomp.",
          message_type = "warning",
          verbose = verbose
        )
        NULL
      }
    )
    if (!is.null(pca_res)) {
      pca_emb <- pca_res$u
    } else {
      pca_emb <- stats::prcomp(
        as.matrix(ge_scaled),
        rank. = min(n_pcs, ncol(ge_scaled)),
        scale. = TRUE
      )[["x"]]
    }
  } else {
    pca_emb <- stats::prcomp(
      as.matrix(ge_scaled),
      rank. = min(n_pcs, ncol(ge_scaled)),
      scale. = TRUE
    )[["x"]]
  }
  colnames(pca_emb) <- paste0("PC", seq_len(ncol(pca_emb)))

  run_metacell_knn <- function(emb, k, grp_label = NULL) {
    n <- nrow(emb)
    k_use <- min(k, n - 1L)
    if (k_use < 2) {
      return(stats::setNames(as.character(seq_len(n)), rownames(emb)))
    }

    knn <- RANN::nn2(emb, k = k_use + 1L)
    knn_idx <- knn[["nn.idx"]][, -1L, drop = FALSE]

    adj <- Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = c(n, n))
    for (i in seq_len(n)) {
      adj[i, knn_idx[i, ]] <- 1
    }
    adj <- adj + Matrix::t(adj)
    adj@x <- rep(1, length(adj@x))

    if (!requireNamespace("igraph", quietly = TRUE)) {
      clusters <- as.integer(stats::cutree(
        stats::hclust(stats::dist(emb), method = "ward.D2"),
        k = max(2L, as.integer(n / k_use))
      ))
    } else {
      g <- igraph::graph_from_adjacency_matrix(
        adj,
        mode = "undirected",
        weighted = NULL
      )
      clusters <- igraph::membership(
        igraph::cluster_louvain(g)
      )
      clusters <- as.integer(clusters)
    }

    membership <- as.character(clusters)
    if (!is.null(grp_label)) {
      membership <- paste0(grp_label, ".", membership)
    }
    names(membership) <- rownames(emb)
    membership
  }

  if (is.null(group_df)) {
    log_message(
      "Running {.pkg MetaCell}-style KNN partitioning with k = {.val {k_knn}} on {.val {ncol(counts)}} cells",
      verbose = verbose
    )

    membership <- run_metacell_knn(
      emb = pca_emb,
      k = k_knn
    )
  } else {
    group_interaction <- do.call(
      interaction,
      c(group_df, list(drop = TRUE, sep = "."))
    )
    group_levels <- levels(group_interaction)

    log_message(
      "Running {.pkg MetaCell}-style partitioning within {.val {length(group_levels)}} group{?s}",
      verbose = verbose
    )

    membership <- character(ncol(counts))
    names(membership) <- colnames(counts)

    for (grp in group_levels) {
      grp_cells <- which(group_interaction == grp)
      if (length(grp_cells) < 3) {
        membership[grp_cells] <- paste0(grp, ".1")
        next
      }

      grp_emb <- pca_emb[grp_cells, , drop = FALSE]
      log_message(
        "  Group {.val {grp}}: {.val {length(grp_cells)}} cells",
        verbose = verbose
      )

      membership[grp_cells] <- run_metacell_knn(
        emb = grp_emb,
        k = k_knn,
        grp_label = grp
      )
    }
  }

  unique_mc <- unique(membership)
  mc_counts <- Matrix::sparseMatrix(
    i = integer(0),
    j = integer(0),
    dims = c(nrow(counts), length(unique_mc))
  )
  rownames(mc_counts) <- rownames(counts)
  colnames(mc_counts) <- unique_mc
  for (mc_name in unique_mc) {
    mc_cells <- which(membership == mc_name)
    mc_counts[, mc_name] <- Matrix::rowSums(
      counts[, mc_cells, drop = FALSE]
    )
  }

  list(
    membership = membership,
    SC = NULL,
    counts = mc_counts
  )
}

lognorm_counts <- function(counts, scale_factor = 10000) {
  lib_sizes <- Matrix::colSums(counts)
  lib_sizes[lib_sizes == 0] <- 1
  log1p(sweep(counts, 2, lib_sizes, "/") * scale_factor)
}

as.MetaCellSeurat <- function(
  srt,
  fields = NULL,
  assay = "RNA",
  ...
) {
  if (!inherits(srt, "Seurat")) {
    log_message(
      "{.arg srt} must be a {.cls Seurat} object",
      message_type = "error"
    )
  }

  mc_tool <- srt@tools[["Metacell"]]
  if (is.null(mc_tool)) {
    log_message(
      "No metacell data found. Run {.fn RunMetaCell} first.",
      message_type = "error"
    )
  }

  mc_sc <- mc_tool[["SC"]]
  mc_counts <- mc_tool[["counts"]]

  if (is.null(mc_counts)) {
    log_message(
      "Metacell counts not available. Re-run {.fn RunMetaCell} or ensure method supports count output.",
      message_type = "error"
    )
  }

  check_r("GfellerLab/SuperCell", verbose = FALSE)

  sc_fields <- list()
  if (!is.null(fields)) {
    for (f in fields) {
      if (f %in% colnames(srt@meta.data)) {
        sc_fields[[f]] <- srt@meta.data[[f]]
        names(sc_fields[[f]]) <- colnames(srt)
      }
    }
  }

  if (!is.null(mc_sc)) {
    mc_sc[names(sc_fields)] <- lapply(sc_fields, function(vals) {
      get_namespace_fun("SuperCell", "supercell_assign")(
        vals,
        supercell_membership = mc_sc$membership,
        method = "jaccard"
      )
    })

    mc_seurat <- get_namespace_fun("SuperCell", "supercell_2_Seurat")(
      SC.GE = as.matrix(mc_counts),
      SC = mc_sc,
      fields = names(sc_fields),
      ...
    )
  } else {
    mc_seurat <- Seurat::CreateSeuratObject(
      counts = mc_counts,
      assay = assay
    )
    membership <- mc_tool[["membership"]]
    for (f in names(sc_fields)) {
      vals <- sc_fields[[f]]
      mc_seurat[[f]] <- vapply(
        split(vals, membership),
        function(x) names(sort(table(x), decreasing = TRUE))[1],
        character(1)
      )
    }
    mc_seurat[["metacell_size"]] <- as.integer(
      table(membership)[colnames(mc_seurat)]
    )
  }

  log_message(
    "{.fn as.MetaCellSeurat} created Seurat object with {.val {ncol(mc_seurat)}} metacells",
    message_type = "success"
  )

  mc_seurat
}
