# Benchmark RunCIBERSORT backend = "r" vs backend = "cpp".
# Run from the package root after installing CIBERSORT and scop.

make_mock_cibersort_data <- function(
  n_genes = 10,
  n_cell_types = 3,
  n_samples = 5,
  seed = 1
) {
  set.seed(seed)
  signature <- matrix(
    stats::runif(n_genes * n_cell_types, 50, 500),
    nrow = n_genes,
    ncol = n_cell_types,
    dimnames = list(
      sprintf("Gene%04d", seq_len(n_genes)),
      sprintf("Cell%02d", seq_len(n_cell_types))
    )
  )
  fractions <- matrix(
    stats::runif(n_cell_types * n_samples),
    nrow = n_cell_types,
    ncol = n_samples,
    dimnames = list(colnames(signature), sprintf("Sample%03d", seq_len(n_samples)))
  )
  fractions <- sweep(fractions, 2, colSums(fractions), "/")
  mixture <- signature %*% fractions
  rownames(mixture) <- rownames(signature)
  colnames(mixture) <- colnames(fractions)
  list(signature = signature, mixture = mixture, truth = t(fractions))
}

make_lm22_simulation <- function(n_samples = 20, seed = 1) {
  signature <- scop:::resolve_cibersort_signature("LM22", verbose = FALSE)$matrix
  signature <- as.matrix(signature)
  signature <- signature[rowSums(signature) > 0, , drop = FALSE]
  set.seed(seed)
  fractions <- matrix(
    stats::rgamma(ncol(signature) * n_samples, shape = 1),
    nrow = ncol(signature),
    ncol = n_samples,
    dimnames = list(colnames(signature), sprintf("LM22_%03d", seq_len(n_samples)))
  )
  fractions <- sweep(fractions, 2, colSums(fractions), "/")
  mixture <- signature %*% fractions
  rownames(mixture) <- rownames(signature)
  colnames(mixture) <- colnames(fractions)
  list(signature = signature, mixture = mixture, truth = t(fractions))
}

extract_fraction <- function(result) {
  x <- result$details$proportion_matrix
  rn <- rownames(x)
  cn <- colnames(x)
  x <- as.matrix(x)
  matrix(
    as.numeric(x),
    nrow = nrow(x),
    ncol = ncol(x),
    dimnames = list(rn, cn)
  )
}

fraction_metrics <- function(old, new) {
  new <- new[rownames(old), colnames(old), drop = FALSE]
  rank_spearman <- vapply(seq_len(nrow(old)), function(i) {
    rho <- suppressWarnings(stats::cor(old[i, ], new[i, ], method = "spearman"))
    if (is.na(rho)) {
      old_order <- order(old[i, ], decreasing = TRUE)
      new_order <- order(new[i, ], decreasing = TRUE)
      if (identical(old_order, new_order)) {
        return(1)
      }
    }
    rho
  }, numeric(1))
  top_k <- min(3L, ncol(old))
  top_k_jaccard <- vapply(seq_len(nrow(old)), function(i) {
    old_top <- order(old[i, ], decreasing = TRUE)[seq_len(top_k)]
    new_top <- order(new[i, ], decreasing = TRUE)[seq_len(top_k)]
    length(intersect(old_top, new_top)) / top_k
  }, numeric(1))
  data.frame(
    fraction_pearson = stats::cor(as.vector(old), as.vector(new)),
    max_abs_diff = max(abs(old - new)),
    mae = mean(abs(old - new)),
    rmse = sqrt(mean((old - new)^2)),
    mean_rank_spearman = mean(rank_spearman, na.rm = TRUE),
    top3_jaccard = mean(top_k_jaccard, na.rm = TRUE),
    major_cell_consistency = mean(max.col(old) == max.col(new))
  )
}

pvalue_metrics <- function(old_stats, new_stats) {
  common <- intersect(rownames(old_stats), rownames(new_stats))
  old_p <- as.numeric(old_stats[common, "P-value"])
  new_p <- as.numeric(new_stats[common, "P-value"])
  bins <- c(-Inf, 0.01, 0.05, 0.1, 1, Inf)
  pvalue_bin_consistency <- if (length(common) == 0L) {
    NA_real_
  } else {
    mean(cut(old_p, bins) == cut(new_p, bins), na.rm = TRUE)
  }
  pvalue_spearman <- if (length(common) < 3 ||
      length(unique(old_p)) < 2 ||
      length(unique(new_p)) < 2) {
    NA_real_
  } else {
    stats::cor(old_p, new_p, method = "spearman")
  }
  data.frame(
    pvalue_spearman = pvalue_spearman,
    pvalue_bin_consistency = pvalue_bin_consistency
  )
}

time_run <- function(fn, repeats = 1L) {
  repeats <- as.integer(repeats)
  if (length(repeats) != 1L || is.na(repeats) || repeats < 1L) {
    repeats <- 1L
  }
  elapsed <- numeric(repeats)
  value <- NULL
  for (i in seq_len(repeats)) {
    gc()
    start <- proc.time()
    value <- fn()
    elapsed[i] <- as.numeric((proc.time() - start)["elapsed"])
  }
  elapsed <- pmax(elapsed, 1e-3)
  list(
    value = value,
    elapsed = stats::median(elapsed),
    elapsed_min = min(elapsed),
    elapsed_sd = if (repeats > 1L) stats::sd(elapsed) else 0,
    elapsed_all = paste(signif(elapsed, 4), collapse = ";")
  )
}

benchmark_one <- function(
  label,
  dat,
  perm,
  QN,
  n_threads = 1L,
  seed = 123L,
  time_repeats = 1L
) {
  r_run <- time_run(function() {
    set.seed(seed)
    scop::RunCIBERSORT(
      count_matrix = dat$mixture,
      sig_matrix = dat$signature,
      backend = "r",
      perm = perm,
      QN = QN,
      verbose = FALSE
    )
  }, repeats = time_repeats)
  cpp_run <- time_run(function() {
    scop::RunCIBERSORT(
      count_matrix = dat$mixture,
      sig_matrix = dat$signature,
      backend = "cpp",
      perm = perm,
      QN = QN,
      n_threads = n_threads,
      seed = seed,
      verbose = FALSE
    )
  }, repeats = time_repeats)

  old_fraction <- extract_fraction(r_run$value)
  new_fraction <- extract_fraction(cpp_run$value)
  new_fraction <- new_fraction[rownames(old_fraction), colnames(old_fraction), drop = FALSE]
  metrics <- cbind(
    data.frame(
      dataset = label,
      samples = ncol(dat$mixture),
      perm = perm,
      QN = QN,
      n_threads = n_threads,
      time_r = r_run$elapsed,
      time_cpp = cpp_run$elapsed,
      speedup = r_run$elapsed / cpp_run$elapsed,
      time_repeats = as.integer(time_repeats),
      time_r_min = r_run$elapsed_min,
      time_cpp_min = cpp_run$elapsed_min,
      time_r_sd = r_run$elapsed_sd,
      time_cpp_sd = cpp_run$elapsed_sd,
      time_r_all = r_run$elapsed_all,
      time_cpp_all = cpp_run$elapsed_all
    ),
    fraction_metrics(old_fraction, new_fraction),
    pvalue_metrics(r_run$value$details$statistics, cpp_run$value$details$statistics)
  )

  fraction_table <- data.frame(
    dataset = label,
    perm = perm,
    QN = QN,
    sample = rep(rownames(old_fraction), times = ncol(old_fraction)),
    cell_type = rep(colnames(old_fraction), each = nrow(old_fraction)),
    r_fraction = as.vector(old_fraction),
    cpp_fraction = as.vector(new_fraction)
  )
  old_stats <- r_run$value$details$statistics
  new_stats <- cpp_run$value$details$statistics[rownames(old_stats), , drop = FALSE]
  pvalue_table <- data.frame(
    dataset = label,
    perm = perm,
    QN = QN,
    sample = rownames(old_stats),
    r_pvalue = as.numeric(old_stats[["P-value"]]),
    cpp_pvalue = as.numeric(new_stats[["P-value"]])
  )
  ranking_table <- data.frame(
    dataset = label,
    perm = perm,
    QN = QN,
    sample = rownames(old_fraction),
    r_major = colnames(old_fraction)[max.col(old_fraction)],
    cpp_major = colnames(new_fraction)[max.col(new_fraction)]
  )
  ranking_table$consistent <- ranking_table$r_major == ranking_table$cpp_major

  list(
    metrics = metrics,
    fractions = fraction_table,
    pvalues = pvalue_table,
    ranking = ranking_table
  )
}

save_cibersort_benchmark <- function(result, output_dir) {
  if (is.null(output_dir)) {
    return(invisible(result))
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  utils::write.csv(result$metrics, file.path(output_dir, "metrics.csv"), row.names = FALSE)
  utils::write.csv(result$fractions, file.path(output_dir, "fractions.csv"), row.names = FALSE)
  utils::write.csv(result$pvalues, file.path(output_dir, "pvalues.csv"), row.names = FALSE)
  utils::write.csv(result$ranking, file.path(output_dir, "ranking.csv"), row.names = FALSE)

  if (!is.null(result$plots) && requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::ggsave(
      file.path(output_dir, "runtime.png"),
      result$plots$runtime,
      width = 8,
      height = 5,
      dpi = 200
    )
    ggplot2::ggsave(
      file.path(output_dir, "speedup.png"),
      result$plots$speedup,
      width = 7,
      height = 4,
      dpi = 200
    )
    ggplot2::ggsave(
      file.path(output_dir, "fraction_parity.png"),
      result$plots$fraction_parity,
      width = 8,
      height = 5,
      dpi = 200
    )
    ggplot2::ggsave(
      file.path(output_dir, "pvalue_ranking.png"),
      result$plots$pvalue_ranking,
      width = 8,
      height = 5,
      dpi = 200
    )
  }
  invisible(result)
}

run_cibersort_benchmark <- function(
  perm_values = c(0, 100),
  QN_values = c(FALSE, TRUE),
  n_threads = 4L,
  include_islet = TRUE,
  small_perm_values = 0,
  time_repeats = 1L,
  seed = 123L,
  output_dir = NULL
) {
  if (!requireNamespace("CIBERSORT", quietly = TRUE)) {
    stop("Install CIBERSORT before running the R-vs-C++ benchmark.")
  }
  datasets <- list(
    small = make_mock_cibersort_data(),
    lm22_20 = make_lm22_simulation(n_samples = 20)
  )
  islet_loaded <- tryCatch({
    utils::data("islet_bulk", package = "scop", envir = environment())
    exists("islet_bulk", inherits = FALSE)
  }, error = function(e) FALSE)
  if (isTRUE(include_islet) && isTRUE(islet_loaded) &&
      requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    counts <- SummarizedExperiment::assay(islet_bulk, "counts")
    datasets$islet_bulk <- make_lm22_simulation(n_samples = min(20, ncol(counts)))
    datasets$islet_bulk$mixture <- counts[
      intersect(rownames(datasets$islet_bulk$signature), rownames(counts)),
      seq_len(min(20, ncol(counts))),
      drop = FALSE
    ]
    datasets$islet_bulk$signature <- datasets$islet_bulk$signature[
      rownames(datasets$islet_bulk$mixture),
      ,
      drop = FALSE
    ]
  }

  out <- list()
  i <- 1L
  for (label in names(datasets)) {
    dataset_perm_values <- if (identical(label, "small")) {
      small_perm_values
    } else {
      perm_values
    }
    for (perm in dataset_perm_values) {
      for (QN in QN_values) {
        out[[i]] <- benchmark_one(
          label,
          datasets[[label]],
          perm,
          QN,
          n_threads = n_threads,
          seed = seed,
          time_repeats = time_repeats
        )
        i <- i + 1L
      }
    }
  }
  result <- list(
    metrics = do.call(rbind, lapply(out, `[[`, "metrics")),
    fractions = do.call(rbind, lapply(out, `[[`, "fractions")),
    pvalues = do.call(rbind, lapply(out, `[[`, "pvalues")),
    ranking = do.call(rbind, lapply(out, `[[`, "ranking"))
  )

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    runtime_plot <- ggplot2::ggplot(result$metrics, ggplot2::aes(samples, time_r, color = "r")) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::geom_point(ggplot2::aes(y = time_cpp, color = "cpp")) +
      ggplot2::geom_line(ggplot2::aes(y = time_cpp, color = "cpp")) +
      ggplot2::facet_grid(QN ~ perm, labeller = ggplot2::label_both) +
      ggplot2::scale_y_log10() +
      ggplot2::labs(x = "Samples", y = "Seconds", color = "Backend")

    speedup_plot <- ggplot2::ggplot(result$metrics, ggplot2::aes(factor(perm), dataset, fill = speedup)) +
      ggplot2::geom_tile() +
      ggplot2::facet_wrap(~QN, labeller = ggplot2::label_both) +
      ggplot2::labs(x = "Permutation", y = "Dataset", fill = "Speedup")

    parity_plot <- ggplot2::ggplot(
      result$fractions,
      ggplot2::aes(r_fraction, cpp_fraction, color = dataset)
    ) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.3) +
      ggplot2::geom_point(alpha = 0.55, size = 0.9) +
      ggplot2::facet_grid(QN ~ perm, labeller = ggplot2::label_both) +
      ggplot2::labs(x = "R backend fraction", y = "C++ backend fraction")

    pvalue_plot <- ggplot2::ggplot(
      merge(
        result$pvalues,
        result$ranking[, c("dataset", "perm", "QN", "sample", "consistent")],
        by = c("dataset", "perm", "QN", "sample")
      ),
      ggplot2::aes(r_pvalue, cpp_pvalue, color = consistent)
    ) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.3) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::facet_grid(QN ~ perm, labeller = ggplot2::label_both) +
      ggplot2::labs(
        x = "R backend P value",
        y = "C++ backend P value",
        color = "Major cell match"
      )

    print(runtime_plot)
    print(speedup_plot)
    print(parity_plot)
    print(pvalue_plot)
    result$plots <- list(
      runtime = runtime_plot,
      speedup = speedup_plot,
      fraction_parity = parity_plot,
      pvalue_ranking = pvalue_plot
    )
  }

  save_cibersort_benchmark(result, output_dir)
  result
}

if (sys.nframe() == 0L) {
  print(run_cibersort_benchmark())
}
