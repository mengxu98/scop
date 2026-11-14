infer_scorpion <- function(expression, arguments, cores) {
  if (!requireNamespace("SCORPION", quietly = TRUE)) {
    stop("Install the SCORPION package before running scorpion.", call. = FALSE)
  }
  original_genes <- rownames(expression)
  pad_small <- arguments$pad_small_input %||% TRUE
  if (nrow(expression) < 10L && isTRUE(pad_small)) {
    scale <- stats::sd(as.numeric(expression))
    if (!is.finite(scale) || scale <= 0) {
      scale <- 1
    }
    signal <- colMeans(expression)
    padding <- t(vapply(
      seq_len(10L - nrow(expression)),
      function(i) signal + i * seq_along(signal) * scale * 1e-9,
      numeric(ncol(expression))
    ))
    rownames(padding) <- paste0(".SCORPION_PAD_", seq_len(nrow(padding)))
    colnames(padding) <- colnames(expression)
    expression <- rbind(expression, padding)
  }
  arguments <- utils::modifyList(
    list(
      n_pc = min(25L, nrow(expression) - 1L, ncol(expression) - 1L),
      gamma_value = 10,
      show_progress = FALSE,
      tf_motifs = NULL,
      ppi_network = NULL,
      pad_small_input = pad_small
    ),
    arguments
  )
  n_pc <- arguments$n_pc
  gamma_value <- arguments$gamma_value
  if (!is.numeric(n_pc) || length(n_pc) != 1L || !is.finite(n_pc) ||
    n_pc < 1L || n_pc != as.integer(n_pc) || n_pc >= min(dim(expression))) {
    stop("SCORPION n_pc must be a positive integer smaller than both expression dimensions.", call. = FALSE)
  }
  if (!is.numeric(gamma_value) || length(gamma_value) != 1L ||
    !is.finite(gamma_value) || gamma_value <= 0) {
    stop("SCORPION gamma_value must be positive and finite.", call. = FALSE)
  }
  fit <- SCORPION::scorpion(
    tfMotifs = arguments$tf_motifs,
    gexMatrix = expression,
    ppiNet = arguments$ppi_network,
    nCores = as.integer(cores),
    nPC = as.integer(n_pc),
    gammaValue = gamma_value,
    showProgress = arguments$show_progress
  )
  scores <- if (is.list(fit)) fit$regNet %||% fit else fit
  extra <- if (is.list(fit)) fit[setdiff(names(fit), "regNet")] else list()
  extra$input_preprocessing <- list(
    original_genes = original_genes,
    padded_genes = setdiff(rownames(expression), original_genes)
  )
  list(
    edges = grn_rna_matrix_edges(as.matrix(scores)),
    extra = extra,
    parameters = arguments,
    implementation = "SCORPION::scorpion",
    package_versions = c(SCORPION = as.character(utils::packageVersion("SCORPION")))
  )
}
