#' @title Simulate Single-Cell Data for Proportion Testing
#'
#' @description
#' Generate a synthetic `Seurat` object with sample-level replicates, condition
#' labels, and controllable cell-type composition shifts.
#' The output is designed to be directly compatible with common `scop` workflows,
#' especially [RunProportionTest].
#' Default settings provide richer cell-type diversity (including a rare group),
#' multiple batches, and stronger sample-level heterogeneity for DA regression.
#'
#' @md
#' @inheritParams thisutils::log_message
#' @param n_genes Number of genes to simulate.
#' @param samples_per_condition Number of biological samples per condition.
#' @param cells_per_sample Number of cells per sample.
#' @param conditions Condition labels.
#' Default is `c("Control", "Treatment")`.
#' @param cell_types Cell type labels.
#' @param base_props Baseline cell-type proportions for the first condition.
#' Must have the same length as `cell_types`.
#' @param condition_prop_shift Additive shift applied to `base_props` for the
#' second condition when `condition_prop_matrix = NULL`.
#' @param condition_prop_matrix Optional matrix of condition-specific proportions.
#' Rows correspond to `conditions`, columns correspond to `cell_types`.
#' When provided, `base_props` and `condition_prop_shift` are ignored.
#' @param prop_concentration Dirichlet concentration for sample-level proportion
#' variability. Higher values produce less variability among samples.
#' @param marker_genes_per_type Number of marker genes for each cell type.
#' @param marker_fc Fold-change multiplier for marker genes.
#' @param condition_de_genes Number of genes with global condition-level effect
#' (for non-reference conditions).
#' @param condition_de_fc Fold-change multiplier for global condition-level genes.
#' @param dispersion Negative binomial size parameter.
#' @param sample_prefix Prefix for simulated sample IDs.
#' @param batch_levels Batch labels assigned to samples in a round-robin manner.
#' @param assay Assay name used in [Seurat::CreateSeuratObject].
#' @param seed Random seed.
#'
#' @return
#' A `Seurat` object with:
#' - metadata columns: `Sample`, `Condition`, `Batch`, `CellType`, `SubCellType`, `Phase`
#' - simulation metadata stored in `srt@tools[["SimulateProportionData"]]`
#'
#' @export
#'
#' @examples
#' \dontrun{
#' srt_sim <- SimulateProportionData(
#'   n_genes = 1500,
#'   samples_per_condition = 4,
#'   cells_per_sample = 250,
#'   seed = 42
#' )
#'
#' # Compatible with proportion testing
#' srt_sim <- RunProportionTest(
#'   srt_sim,
#'   group.by = "CellType",
#'   split.by = "Condition",
#'   sample.by = "Sample",
#'   proportion_method = "propeller",
#'   comparison = list(c("Control", "Treatment"))
#' )
#' }
SimulateProportionData <- function(
    n_genes = 2000L,
    samples_per_condition = 5L,
    cells_per_sample = 350L,
    conditions = c("Control", "Treatment"),
    cell_types = c("Tcell", "Bcell", "Myeloid", "NK", "DC", "Mono", "Epithelial", "Fibro", "Mast"),
    base_props = c(0.20, 0.14, 0.16, 0.10, 0.08, 0.10, 0.12, 0.07, 0.03),
    condition_prop_shift = c(0.07, -0.03, -0.03, -0.01, -0.01, -0.01, 0.01, 0.00, 0.01),
    condition_prop_matrix = NULL,
    prop_concentration = 35,
    marker_genes_per_type = 30L,
    marker_fc = 4,
    condition_de_genes = 100L,
    condition_de_fc = 1.5,
    dispersion = 2,
    sample_prefix = "S",
    batch_levels = c("Batch1", "Batch2", "Batch3"),
    assay = "RNA",
    seed = 11,
    verbose = TRUE) {
  set.seed(seed)

  n_genes <- as.integer(n_genes)
  samples_per_condition <- as.integer(samples_per_condition)
  cells_per_sample <- as.integer(cells_per_sample)
  marker_genes_per_type <- as.integer(marker_genes_per_type)
  condition_de_genes <- as.integer(condition_de_genes)

  if (n_genes < 100) {
    log_message(
      "{.arg n_genes} must be >= 100",
      message_type = "error"
    )
  }
  if (samples_per_condition < 2) {
    log_message(
      "{.arg samples_per_condition} must be >= 2 for DA analysis",
      message_type = "error"
    )
  }
  if (cells_per_sample < 50) {
    log_message(
      "{.arg cells_per_sample} should be >= 50",
      message_type = "error"
    )
  }
  if (length(conditions) < 2) {
    log_message(
      "{.arg conditions} must contain at least two groups",
      message_type = "error"
    )
  }
  if (length(cell_types) < 2) {
    log_message(
      "{.arg cell_types} must contain at least two cell types",
      message_type = "error"
    )
  }

  n_types <- length(cell_types)
  n_conditions <- length(conditions)

  if (is.null(condition_prop_matrix)) {
    if (length(base_props) != n_types) {
      log_message(
        "Length of {.arg base_props} must equal length of {.arg cell_types}",
        message_type = "error"
      )
    }
    if (length(condition_prop_shift) != n_types) {
      log_message(
        "Length of {.arg condition_prop_shift} must equal length of {.arg cell_types}",
        message_type = "error"
      )
    }

    condition_prop_matrix <- matrix(
      rep(base_props, times = n_conditions),
      nrow = n_conditions,
      byrow = TRUE
    )
    if (n_conditions >= 2) {
      condition_prop_matrix[2, ] <- base_props + condition_prop_shift
    }
    if (n_conditions > 2) {
      for (i in 3:n_conditions) {
        condition_prop_matrix[i, ] <- base_props + stats::rnorm(n_types, mean = 0, sd = 0.03)
      }
    }
  } else {
    condition_prop_matrix <- as.matrix(condition_prop_matrix)
    if (!all(dim(condition_prop_matrix) == c(n_conditions, n_types))) {
      log_message(
        "{.arg condition_prop_matrix} must have dimensions length(conditions) x length(cell_types)",
        message_type = "error"
      )
    }
  }

  condition_prop_matrix <- pmax(condition_prop_matrix, 1e-4)
  condition_prop_matrix <- condition_prop_matrix / rowSums(condition_prop_matrix)
  rownames(condition_prop_matrix) <- conditions
  colnames(condition_prop_matrix) <- cell_types

  n_samples <- n_conditions * samples_per_condition
  sample_ids <- paste0(sample_prefix, seq_len(n_samples))
  sample_condition <- rep(conditions, each = samples_per_condition)
  sample_batch <- rep(batch_levels, length.out = n_samples)

  log_message(
    "Simulating {.val {n_samples}} samples and {.val {n_samples * cells_per_sample}} cells",
    verbose = verbose
  )

  sample_prop_list <- vector("list", length = n_samples)
  names(sample_prop_list) <- sample_ids

  cell_meta <- vector("list", length = n_samples)
  for (i in seq_len(n_samples)) {
    cond_i <- sample_condition[i]
    alpha <- condition_prop_matrix[cond_i, ] * prop_concentration
    p_sample <- .rdirichlet1(alpha)
    p_sample <- p_sample / sum(p_sample)

    ct_counts <- as.integer(stats::rmultinom(1, size = cells_per_sample, prob = p_sample))
    names(ct_counts) <- cell_types
    sample_prop_list[[i]] <- p_sample

    ct_labels <- rep(cell_types, times = ct_counts)
    n_cells_i <- length(ct_labels)
    cell_ids <- paste0(sample_ids[i], "_Cell", seq_len(n_cells_i))
    phase_levels <- c("G1", "S", "G2M")

    cell_meta[[i]] <- data.frame(
      cell_id = cell_ids,
      Sample = sample_ids[i],
      Condition = cond_i,
      Batch = sample_batch[i],
      CellType = ct_labels,
      SubCellType = paste0(ct_labels, "_sub"),
      Phase = sample(phase_levels, size = n_cells_i, replace = TRUE, prob = c(0.55, 0.25, 0.20)),
      stringsAsFactors = FALSE
    )
  }

  meta_df <- do.call(rbind, cell_meta)
  rownames(meta_df) <- meta_df$cell_id
  meta_df$cell_id <- NULL
  meta_df$Condition <- factor(meta_df$Condition, levels = conditions)
  meta_df$CellType <- factor(meta_df$CellType, levels = cell_types)
  meta_df$SubCellType <- factor(meta_df$SubCellType)
  meta_df$Phase <- factor(meta_df$Phase, levels = c("G1", "S", "G2M"))
  meta_df$Batch <- factor(meta_df$Batch, levels = unique(sample_batch))

  cell_ids_all <- rownames(meta_df)
  n_cells <- length(cell_ids_all)
  gene_ids <- sprintf("Gene%04d", seq_len(n_genes))

  baseline_expr <- stats::rgamma(n_genes, shape = 2, rate = 1)
  baseline_expr <- pmax(baseline_expr, 0.01)

  marker_genes_per_type <- min(marker_genes_per_type, floor(n_genes / n_types))
  marker_idx <- split(seq_len(n_types * marker_genes_per_type), rep(cell_types, each = marker_genes_per_type))
  marker_idx <- lapply(marker_idx, function(x) x[x <= n_genes])

  de_idx <- seq_len(min(condition_de_genes, n_genes))
  cond_fc <- matrix(1, nrow = n_conditions, ncol = n_genes)
  rownames(cond_fc) <- conditions
  if (n_conditions >= 2 && length(de_idx) > 0) {
    cond_fc[2, de_idx] <- condition_de_fc
  }
  if (n_conditions > 2 && length(de_idx) > 0) {
    for (i in 3:n_conditions) {
      cond_fc[i, de_idx] <- 1 + (condition_de_fc - 1) * (i - 1) / max(1, n_conditions - 1)
    }
  }

  type_fc <- matrix(1, nrow = n_types, ncol = n_genes)
  rownames(type_fc) <- cell_types
  for (ct in cell_types) {
    idx <- marker_idx[[ct]]
    if (length(idx) > 0) {
      type_fc[ct, idx] <- marker_fc
    }
  }

  counts <- matrix(0L, nrow = n_genes, ncol = n_cells, dimnames = list(gene_ids, cell_ids_all))

  type_levels <- as.character(meta_df$CellType)
  cond_levels <- as.character(meta_df$Condition)
  sample_levels <- as.character(meta_df$Sample)
  unique_groups <- unique(paste(type_levels, cond_levels, sample_levels, sep = "||"))

  for (grp in unique_groups) {
    idx_cells <- which(paste(type_levels, cond_levels, sample_levels, sep = "||") == grp)
    if (length(idx_cells) == 0) {
      next
    }
    parts <- strsplit(grp, "\\|\\|")[[1]]
    ct <- parts[1]
    cond <- parts[2]

    mu_gene <- baseline_expr * type_fc[ct, ] * cond_fc[cond, ]
    lib_factor <- stats::rlnorm(length(idx_cells), meanlog = 0, sdlog = 0.35)
    mu_mat <- outer(mu_gene, lib_factor, "*")

    counts_grp <- stats::rnbinom(
      n = length(mu_mat),
      size = dispersion,
      mu = as.vector(mu_mat)
    )
    counts[, idx_cells] <- matrix(
      as.integer(counts_grp),
      nrow = n_genes,
      ncol = length(idx_cells),
      byrow = FALSE
    )
  }

  counts_sparse <- Matrix::Matrix(counts, sparse = TRUE)
  srt <- Seurat::CreateSeuratObject(
    counts = counts_sparse,
    meta.data = meta_df,
    assay = assay
  )

  marker_gene_names <- lapply(marker_idx, function(x) gene_ids[x])
  sample_props_mat <- do.call(
    rbind,
    lapply(sample_prop_list, function(x) {
      stats::setNames(as.numeric(x), cell_types)
    })
  )
  rownames(sample_props_mat) <- sample_ids

  srt@tools[["SimulateProportionData"]] <- list(
    seed = seed,
    condition_props = condition_prop_matrix,
    sample_props = sample_props_mat,
    marker_genes = marker_gene_names,
    parameters = list(
      n_genes = n_genes,
      samples_per_condition = samples_per_condition,
      cells_per_sample = cells_per_sample,
      conditions = conditions,
      cell_types = cell_types,
      marker_genes_per_type = marker_genes_per_type,
      marker_fc = marker_fc,
      condition_de_genes = condition_de_genes,
      condition_de_fc = condition_de_fc,
      dispersion = dispersion,
      prop_concentration = prop_concentration
    )
  )

  log_message(
    "Synthetic dataset generated: {.val {ncol(srt)}} cells x {.val {nrow(srt)}} genes",
    message_type = "success",
    verbose = verbose
  )

  srt
}

.rdirichlet1 <- function(alpha) {
  x <- stats::rgamma(length(alpha), shape = alpha, rate = 1)
  x / sum(x)
}
