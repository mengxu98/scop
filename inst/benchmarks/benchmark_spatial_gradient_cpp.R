# Benchmark RunSpatialGradientFeatures backend = "spata2" vs backend = "cpp".
# Run from the package root after installing SPATA2 and scop.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

default_thresholds <- list(
  top_k_jaccard = 0.70,
  rank_spearman = 0.80,
  fdr_spearman = 0.90,
  p_value_spearman = 0.90,
  tot_var_spearman = 0.80,
  estimate_curve_pearson = 0.90,
  value_curve_pearson = 0.80
)

parse_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list(
    output_dir = file.path("inst", "benchmarks", "spatial_gradient_cpp"),
    reference = "trajectory",
    n_spots = 80L,
    n_variables = 6L,
    n_random = 20L,
    n_bins = 25L,
    resolution = 25,
    seed = 1L,
    run_spata2 = TRUE,
    strict = FALSE,
    verbose = FALSE
  )
  for (arg in args) {
    if (identical(arg, "--skip-spata2") || identical(arg, "--cpp-only")) {
      out$run_spata2 <- FALSE
    } else if (identical(arg, "--strict")) {
      out$strict <- TRUE
    } else if (identical(arg, "--verbose")) {
      out$verbose <- TRUE
    } else if (grepl("^--[^=]+=", arg)) {
      key <- sub("^--([^=]+)=.*$", "\\1", arg)
      value <- sub("^--[^=]+=", "", arg)
      key <- chartr("-", "_", key)
      out[[key]] <- value
    }
  }
  int_keys <- c("n_spots", "n_variables", "n_random", "n_bins", "seed")
  for (key in int_keys) {
    out[[key]] <- as.integer(out[[key]])
  }
  out$resolution <- as.numeric(out$resolution)
  out$reference <- match.arg(out$reference, c("trajectory", "annotation"))
  out
}

load_scop_for_benchmark <- function(package_root = ".") {
  loaded <- FALSE
  if (file.exists(file.path(package_root, "DESCRIPTION")) &&
      requireNamespace("pkgload", quietly = TRUE)) {
    loaded <- tryCatch({
      pkgload::load_all(package_root, quiet = TRUE)
      TRUE
    }, error = function(e) FALSE)
  }
  if (!isTRUE(loaded)) {
    suppressPackageStartupMessages(library(scop))
  }
  invisible(TRUE)
}

format_seconds <- function(x) {
  if (!is.finite(x)) {
    return("NA")
  }
  if (x < 0.01) {
    sprintf("%.3fs", x)
  } else if (x < 10) {
    sprintf("%.2fs", x)
  } else {
    sprintf("%.1fs", x)
  }
}

format_number <- function(x, digits = 3) {
  if (!is.finite(x)) {
    return("NA")
  }
  format(signif(x, digits), trim = TRUE, scientific = FALSE)
}

safe_cor <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x, y, method = method))
}

safe_mae <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (!any(keep)) {
    return(NA_real_)
  }
  mean(abs(x[keep] - y[keep]))
}

time_run <- function(fn) {
  warnings <- character()
  start <- proc.time()
  value <- tryCatch(
    withCallingHandlers(
      fn(),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  elapsed <- as.numeric((proc.time() - start)[["elapsed"]])
  if (inherits(value, "error")) {
    list(ok = FALSE, value = NULL, elapsed = elapsed, error = conditionMessage(value), warnings = warnings)
  } else {
    list(ok = TRUE, value = value, elapsed = elapsed, error = NA_character_, warnings = warnings)
  }
}

make_spatial_gradient_fixture <- function(n_spots = 80L, n_variables = 6L, seed = 1L) {
  data(visium_human_pancreas_sub, package = "scop", envir = environment())
  srt <- visium_human_pancreas_sub
  SeuratObject::DefaultAssay(srt) <- "Spatial"
  if (is.finite(n_spots) && n_spots > 0L && n_spots < ncol(srt)) {
    meta <- srt@meta.data
    set.seed(seed)
    panin_cells <- rownames(meta)[is.finite(meta$coda_panin) & meta$coda_panin > 0]
    keep_panin <- head(panin_cells[order(meta[panin_cells, "coda_panin"], decreasing = TRUE)], min(length(panin_cells), max(1L, n_spots %/% 5L)))
    other_cells <- setdiff(colnames(srt), keep_panin)
    keep_other <- sample(other_cells, min(length(other_cells), max(0L, n_spots - length(keep_panin))))
    cells <- unique(c(keep_panin, keep_other))
    srt <- subset(srt, cells = cells)
  }

  candidate_variables <- c(
    "TMSB4X", "COL1A1", "COL3A1", "FOS", "B2M", "IGFBP7",
    "SPP1", "KRT19", "MALAT1", "ACTB"
  )
  variables <- intersect(candidate_variables, rownames(srt))
  if (length(variables) < n_variables) {
    fallback <- setdiff(rownames(srt), variables)
    variables <- c(variables, head(fallback, n_variables - length(variables)))
  }
  variables <- head(unique(variables), n_variables)

  if (all(c("x", "y") %in% colnames(srt@meta.data))) {
    coords <- data.frame(
      x = srt@meta.data[["x"]],
      y = srt@meta.data[["y"]],
      row.names = rownames(srt@meta.data),
      stringsAsFactors = FALSE
    )
  } else {
    coords <- spatial_dim_coords(
      srt = srt,
      image = NULL,
      coord.cols = c("x", "y"),
      overlay_image = FALSE
    )$data
  }
  coords <- coords[colnames(srt), c("x", "y"), drop = FALSE]
  meta_use <- srt@meta.data[rownames(coords), , drop = FALSE]
  if ("coda_panin" %in% colnames(meta_use) && any(meta_use$coda_panin > 0, na.rm = TRUE)) {
    start_idx <- which.max(meta_use$coda_panin)
  } else {
    start_idx <- which.min(coords$x)
  }
  start <- as.numeric(coords[start_idx, c("x", "y")])
  sq_dist <- (coords$x - start[1])^2 + (coords$y - start[2])^2
  end <- as.numeric(coords[which.max(sq_dist), c("x", "y")])

  list(
    srt = srt,
    variables = variables,
    start = start,
    end = end,
    annotation.variable = "coda_panin",
    annotation.threshold = 0
  )
}

extract_gradient_result <- function(srt, result_name) {
  x <- srt@tools[["SpatialGradientFeatures"]][[result_name]]
  if (is.null(x)) {
    stop("Stored spatial gradient result was not found: ", result_name, call. = FALSE)
  }
  x
}

run_gradient_backend <- function(
  fixture,
  backend,
  reference,
  result_name,
  n_random,
  n_bins,
  resolution,
  seed,
  verbose
) {
  if (identical(backend, "spata2") && !requireNamespace("SPATA2", quietly = TRUE)) {
    return(list(
      backend = backend,
      status = "skipped",
      elapsed = NA_real_,
      error = "SPATA2 is not installed",
      warnings = character(),
      result = NULL
    ))
  }
  if (identical(backend, "spata2")) {
    suppressPackageStartupMessages(library(SPATA2))
  }

  base_args <- list(
    srt = fixture$srt,
    reference = reference,
    backend = backend,
    result_name = result_name,
    assay = "Spatial",
    layer = "counts",
    platform = "VisiumSmall",
    variables = fixture$variables,
    n_random = n_random,
    seed = seed,
    sign_threshold = 1,
    nfeatures = length(fixture$variables),
    verbose = verbose
  )
  if (identical(reference, "trajectory")) {
    base_args$start <- fixture$start
    base_args$end <- fixture$end
  } else {
    base_args$annotation.variable <- fixture$annotation.variable
    base_args$annotation.threshold <- fixture$annotation.threshold
  }
  if (identical(backend, "cpp")) {
    base_args$n_bins <- n_bins
    base_args$min_spots <- 1L
  } else {
    base_args$resolution <- resolution
    base_args$force_comp <- TRUE
    base_args$rm_zero_infl <- FALSE
    base_args$estimate_R2 <- FALSE
  }

  timed <- time_run(function() {
    out <- do.call(scop::RunSpatialGradientFeatures, base_args)
    extract_gradient_result(out, result_name)
  })
  list(
    backend = backend,
    status = if (isTRUE(timed$ok)) "ok" else "error",
    elapsed = timed$elapsed,
    error = timed$error,
    warnings = timed$warnings,
    result = timed$value
  )
}

table_nrow <- function(result, name) {
  if (is.null(result) || !is.data.frame(result[[name]])) {
    return(NA_integer_)
  }
  nrow(result[[name]])
}

schema_ok <- function(result) {
  expected <- c("screening", "significance", "model_fits", "top_variables", "parameters")
  if (is.null(result) || !identical(names(result), expected)) {
    return(FALSE)
  }
  all(vapply(result, is.data.frame, logical(1)))
}

contains_spata2_object <- function(x) {
  if (is.null(x)) {
    return(FALSE)
  }
  if (methods::is(x, "SPATA2")) {
    return(TRUE)
  }
  if (is.list(x)) {
    return(any(vapply(x, contains_spata2_object, logical(1))))
  }
  FALSE
}

result_summary <- function(run, case) {
  data.frame(
    case = case,
    backend = run$backend,
    status = run$status,
    elapsed = run$elapsed,
    error = run$error,
    warnings = paste(unique(run$warnings), collapse = " | "),
    screening_rows = table_nrow(run$result, "screening"),
    significance_rows = table_nrow(run$result, "significance"),
    model_fits_rows = table_nrow(run$result, "model_fits"),
    top_variables_rows = table_nrow(run$result, "top_variables"),
    schema_ok = schema_ok(run$result),
    no_spata2_object = !contains_spata2_object(run$result),
    stringsAsFactors = FALSE
  )
}

merge_by_variable <- function(ref, cmp, table_name, cols) {
  if (is.null(ref) || is.null(cmp) ||
      !is.data.frame(ref[[table_name]]) || !is.data.frame(cmp[[table_name]])) {
    return(data.frame())
  }
  a <- ref[[table_name]]
  b <- cmp[[table_name]]
  if (!"variable" %in% colnames(a) || !"variable" %in% colnames(b)) {
    return(data.frame())
  }
  keep_a <- unique(c("variable", intersect(cols, colnames(a))))
  keep_b <- unique(c("variable", intersect(cols, colnames(b))))
  a <- a[, keep_a, drop = FALSE]
  b <- b[, keep_b, drop = FALSE]
  merge(a, b, by = "variable", suffixes = c("_spata2", "_cpp"))
}

compare_significance <- function(ref, cmp) {
  merged <- merge_by_variable(ref, cmp, "significance", c("tot_var", "p_value", "fdr"))
  if (nrow(merged) == 0L) {
    return(data.frame(
      common_variables = 0,
      fdr_spearman = NA_real_,
      p_value_spearman = NA_real_,
      tot_var_spearman = NA_real_,
      fdr_mae = NA_real_,
      p_value_mae = NA_real_,
      significance_call_agreement = NA_real_
    ))
  }
  fdr_ref <- as.numeric(merged$fdr_spata2)
  fdr_cpp <- as.numeric(merged$fdr_cpp)
  p_ref <- as.numeric(merged$p_value_spata2)
  p_cpp <- as.numeric(merged$p_value_cpp)
  tot_ref <- as.numeric(merged$tot_var_spata2)
  tot_cpp <- as.numeric(merged$tot_var_cpp)
  data.frame(
    common_variables = nrow(merged),
    fdr_spearman = safe_cor(fdr_ref, fdr_cpp, method = "spearman"),
    p_value_spearman = safe_cor(p_ref, p_cpp, method = "spearman"),
    tot_var_spearman = safe_cor(tot_ref, tot_cpp, method = "spearman"),
    fdr_mae = safe_mae(fdr_ref, fdr_cpp),
    p_value_mae = safe_mae(p_ref, p_cpp),
    significance_call_agreement = mean((fdr_ref <= 0.05) == (fdr_cpp <= 0.05), na.rm = TRUE)
  )
}

compare_top_variables <- function(ref, cmp, k = 10L) {
  ref_top <- ref$top_variables$variable %||% character()
  cmp_top <- cmp$top_variables$variable %||% character()
  ref_top <- ref_top[!is.na(ref_top) & nzchar(ref_top)]
  cmp_top <- cmp_top[!is.na(cmp_top) & nzchar(cmp_top)]
  k <- min(k, length(ref_top), length(cmp_top))
  if (k < 1L) {
    return(data.frame(top_k = 0, top_k_jaccard = NA_real_, rank_spearman = NA_real_, top1_match = NA_real_))
  }
  ref_k <- head(ref_top, k)
  cmp_k <- head(cmp_top, k)
  union_k <- union(ref_k, cmp_k)
  common <- intersect(ref_top, cmp_top)
  ref_rank <- match(common, ref_top)
  cmp_rank <- match(common, cmp_top)
  data.frame(
    top_k = k,
    top_k_jaccard = length(intersect(ref_k, cmp_k)) / length(union_k),
    rank_spearman = safe_cor(ref_rank, cmp_rank, method = "spearman"),
    top1_match = identical(ref_top[1], cmp_top[1])
  )
}

curve_vector <- function(df, variable, column, grid) {
  if (!all(c("variable", "distance", column) %in% colnames(df))) {
    return(NULL)
  }
  x <- df[df$variable == variable, c("distance", column), drop = FALSE]
  distance <- suppressWarnings(as.numeric(x$distance))
  value <- suppressWarnings(as.numeric(x[[column]]))
  keep <- is.finite(distance) & is.finite(value)
  distance <- distance[keep]
  value <- value[keep]
  if (length(distance) < 2L || diff(range(distance)) <= 0) {
    return(NULL)
  }
  distance <- (distance - min(distance)) / diff(range(distance))
  agg <- stats::aggregate(value, list(distance = distance), mean)
  if (nrow(agg) < 2L) {
    return(NULL)
  }
  stats::approx(agg$distance, agg$x, xout = grid, rule = 2, ties = mean)$y
}

compare_screening_curves <- function(ref, cmp, grid_n = 50L) {
  if (!is.data.frame(ref$screening) || !is.data.frame(cmp$screening)) {
    return(data.frame(
      curve_variables = 0,
      estimate_curve_pearson = NA_real_,
      estimate_curve_spearman = NA_real_,
      value_curve_pearson = NA_real_,
      value_curve_spearman = NA_real_
    ))
  }
  variables <- intersect(ref$screening$variable, cmp$screening$variable)
  grid <- seq(0, 1, length.out = grid_n)
  out <- list()
  for (column in c("estimate", "value")) {
    ref_values <- numeric()
    cmp_values <- numeric()
    used <- character()
    for (variable in variables) {
      a <- curve_vector(ref$screening, variable, column, grid)
      b <- curve_vector(cmp$screening, variable, column, grid)
      if (!is.null(a) && !is.null(b)) {
        ref_values <- c(ref_values, a)
        cmp_values <- c(cmp_values, b)
        used <- c(used, variable)
      }
    }
    out[[paste0(column, "_variables")]] <- length(unique(used))
    out[[paste0(column, "_curve_pearson")]] <- safe_cor(ref_values, cmp_values, method = "pearson")
    out[[paste0(column, "_curve_spearman")]] <- safe_cor(ref_values, cmp_values, method = "spearman")
  }
  data.frame(
    curve_variables = max(out$estimate_variables, out$value_variables),
    estimate_curve_pearson = out$estimate_curve_pearson,
    estimate_curve_spearman = out$estimate_curve_spearman,
    value_curve_pearson = out$value_curve_pearson,
    value_curve_spearman = out$value_curve_spearman
  )
}

compare_backends <- function(spata2_result, cpp_result) {
  if (is.null(spata2_result) || is.null(cpp_result)) {
    return(data.frame(
      common_variables = 0,
      fdr_spearman = NA_real_,
      p_value_spearman = NA_real_,
      tot_var_spearman = NA_real_,
      fdr_mae = NA_real_,
      p_value_mae = NA_real_,
      significance_call_agreement = NA_real_,
      top_k = 0,
      top_k_jaccard = NA_real_,
      rank_spearman = NA_real_,
      top1_match = NA,
      curve_variables = 0,
      estimate_curve_pearson = NA_real_,
      estimate_curve_spearman = NA_real_,
      value_curve_pearson = NA_real_,
      value_curve_spearman = NA_real_
    ))
  }
  cbind(
    compare_significance(spata2_result, cpp_result),
    compare_top_variables(spata2_result, cpp_result),
    compare_screening_curves(spata2_result, cpp_result)
  )
}

metric_status <- function(value, threshold, higher_is_better = TRUE) {
  if (!is.finite(value) || !is.finite(threshold)) {
    return("na")
  }
  pass <- if (isTRUE(higher_is_better)) value >= threshold else value <= threshold
  if (isTRUE(pass)) "pass" else "fail"
}

metric_row <- function(case, method, metric, group, value, score, label, status, threshold = NA_real_, higher_is_better = TRUE) {
  data.frame(
    case = case,
    method = method,
    metric = metric,
    group = group,
    value = value,
    score = pmin(pmax(score, 0), 1),
    label = label,
    status = status,
    threshold = threshold,
    higher_is_better = higher_is_better,
    stringsAsFactors = FALSE
  )
}

make_summary_metrics <- function(case, runs, pairwise) {
  run_df <- do.call(rbind, lapply(runs, result_summary, case = case))
  ok_times <- run_df$elapsed[run_df$status == "ok" & is.finite(run_df$elapsed)]
  min_time <- if (length(ok_times) > 0L) min(ok_times) else NA_real_
  spata_time <- run_df$elapsed[match("spata2", run_df$backend)]
  rows <- list()
  i <- 1L

  for (run in runs) {
    method <- if (identical(run$backend, "spata2")) "SPATA2" else "SCOP C++"
    row <- run_df[run_df$backend == run$backend, , drop = FALSE]
    ok <- identical(row$status, "ok")
    runtime_score <- if (ok && is.finite(min_time) && is.finite(row$elapsed)) min_time / row$elapsed else 0
    speedup <- if (identical(run$backend, "cpp") && ok && is.finite(spata_time) && is.finite(row$elapsed)) {
      spata_time / row$elapsed
    } else if (identical(run$backend, "spata2") && ok) {
      1
    } else {
      NA_real_
    }
    speedup_score <- if (is.finite(speedup)) min(log1p(speedup) / log1p(max(speedup, 1)), 1) else 0

    rows[[i]] <- metric_row(case, method, "Runtime", "Runtime", row$elapsed, runtime_score, format_seconds(row$elapsed), row$status)
    i <- i + 1L
    rows[[i]] <- metric_row(case, method, "Speedup", "Runtime", speedup, speedup_score, if (is.finite(speedup)) paste0(format_number(speedup, 3), "x") else "NA", row$status)
    i <- i + 1L
    rows[[i]] <- metric_row(case, method, "Schema", "Integrity", as.numeric(row$schema_ok), as.numeric(row$schema_ok), if (row$schema_ok) "pass" else "fail", if (row$schema_ok) "pass" else "fail", threshold = 1)
    i <- i + 1L
    rows[[i]] <- metric_row(case, method, "No S4 object", "Integrity", as.numeric(row$no_spata2_object), as.numeric(row$no_spata2_object), if (row$no_spata2_object) "pass" else "fail", if (row$no_spata2_object) "pass" else "fail", threshold = 1)
    i <- i + 1L

    agreement_values <- if (identical(run$backend, "spata2") && ok) {
      list(
        "TopK Jaccard" = 1,
        "Rank rho" = 1,
        "FDR rho" = 1,
        "P-value rho" = 1,
        "Curve rho" = 1
      )
    } else if (identical(run$backend, "cpp") && ok) {
      list(
        "TopK Jaccard" = pairwise$top_k_jaccard,
        "Rank rho" = pairwise$rank_spearman,
        "FDR rho" = pairwise$fdr_spearman,
        "P-value rho" = pairwise$p_value_spearman,
        "Curve rho" = pairwise$estimate_curve_pearson
      )
    } else {
      list(
        "TopK Jaccard" = NA_real_,
        "Rank rho" = NA_real_,
        "FDR rho" = NA_real_,
        "P-value rho" = NA_real_,
        "Curve rho" = NA_real_
      )
    }
    thresholds <- c(
      "TopK Jaccard" = default_thresholds$top_k_jaccard,
      "Rank rho" = default_thresholds$rank_spearman,
      "FDR rho" = default_thresholds$fdr_spearman,
      "P-value rho" = default_thresholds$p_value_spearman,
      "Curve rho" = default_thresholds$estimate_curve_pearson
    )
    for (metric in names(agreement_values)) {
      value <- as.numeric(agreement_values[[metric]])
      status <- if (identical(run$backend, "spata2") && ok) {
        "reference"
      } else {
        metric_status(value, thresholds[[metric]], higher_is_better = TRUE)
      }
      rows[[i]] <- metric_row(
        case = case,
        method = method,
        metric = metric,
        group = "Agreement",
        value = value,
        score = if (is.finite(value)) value else 0,
        label = format_number(value, 3),
        status = status,
        threshold = thresholds[[metric]]
      )
      i <- i + 1L
    }
  }
  do.call(rbind, rows)
}

write_backend_tables <- function(run, output_dir) {
  if (is.null(run$result)) {
    return(invisible(FALSE))
  }
  for (name in names(run$result)) {
    utils::write.csv(
      run$result[[name]],
      file.path(output_dir, paste0(run$backend, "_", name, ".csv")),
      row.names = FALSE
    )
  }
  invisible(TRUE)
}

save_benchmark <- function(runs, pairwise, summary_metrics, output_dir, case) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  run_df <- do.call(rbind, lapply(runs, result_summary, case = case))
  utils::write.csv(run_df, file.path(output_dir, "backend_runs.csv"), row.names = FALSE)
  utils::write.csv(pairwise, file.path(output_dir, "pairwise_metrics.csv"), row.names = FALSE)
  utils::write.csv(summary_metrics, file.path(output_dir, "summary_metrics.csv"), row.names = FALSE)
  for (run in runs) {
    write_backend_tables(run, output_dir)
  }
  invisible(TRUE)
}

strict_check <- function(runs, pairwise) {
  status <- vapply(runs, `[[`, character(1), "status")
  if (!all(status == "ok")) {
    stop("Strict benchmark failed because at least one backend did not finish.", call. = FALSE)
  }
  failures <- character()
  for (name in names(default_thresholds)) {
    value <- as.numeric(pairwise[[name]])
    if (!is.finite(value) || value < default_thresholds[[name]]) {
      failures <- c(failures, sprintf("%s=%s", name, format_number(value, 4)))
    }
  }
  if (length(failures) > 0L) {
    stop("Strict benchmark failed: ", paste(failures, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

run_spatial_gradient_benchmark <- function(
  output_dir = file.path("inst", "benchmarks", "spatial_gradient_cpp"),
  reference = c("trajectory", "annotation"),
  n_spots = 80L,
  n_variables = 6L,
  n_random = 20L,
  n_bins = 25L,
  resolution = 25,
  seed = 1L,
  run_spata2 = TRUE,
  strict = FALSE,
  verbose = FALSE
) {
  reference <- match.arg(reference)
  load_scop_for_benchmark(".")
  fixture <- make_spatial_gradient_fixture(n_spots = n_spots, n_variables = n_variables, seed = seed)
  case <- paste0(reference, "_", length(fixture$variables), "genes_", ncol(fixture$srt), "spots")

  cpp_run <- run_gradient_backend(
    fixture = fixture,
    backend = "cpp",
    reference = reference,
    result_name = paste0("benchmark_", reference, "_cpp"),
    n_random = n_random,
    n_bins = n_bins,
    resolution = resolution,
    seed = seed,
    verbose = verbose
  )
  spata2_run <- if (isTRUE(run_spata2)) {
    run_gradient_backend(
      fixture = fixture,
      backend = "spata2",
      reference = reference,
      result_name = paste0("benchmark_", reference, "_spata2"),
      n_random = n_random,
      n_bins = n_bins,
      resolution = resolution,
      seed = seed,
      verbose = verbose
    )
  } else {
    list(
      backend = "spata2",
      status = "skipped",
      elapsed = NA_real_,
      error = "Skipped by --skip-spata2",
      warnings = character(),
      result = NULL
    )
  }
  runs <- list(spata2 = spata2_run, cpp = cpp_run)
  pairwise <- compare_backends(spata2_run$result, cpp_run$result)
  pairwise$case <- case
  pairwise$reference <- reference
  pairwise <- pairwise[, c("case", "reference", setdiff(colnames(pairwise), c("case", "reference"))), drop = FALSE]
  summary_metrics <- make_summary_metrics(case = case, runs = runs, pairwise = pairwise)
  save_benchmark(runs, pairwise, summary_metrics, output_dir, case)
  if (isTRUE(strict)) {
    strict_check(runs, pairwise)
  }
  list(runs = runs, pairwise = pairwise, summary_metrics = summary_metrics, output_dir = output_dir)
}

if (sys.nframe() == 0L) {
  args <- parse_cli_args()
  result <- do.call(run_spatial_gradient_benchmark, args)
  cat("output_dir=", normalizePath(result$output_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("backend_runs=", normalizePath(file.path(result$output_dir, "backend_runs.csv"), winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("pairwise_metrics=", normalizePath(file.path(result$output_dir, "pairwise_metrics.csv"), winslash = "/", mustWork = FALSE), "\n", sep = "")
  cat("summary_metrics=", normalizePath(file.path(result$output_dir, "summary_metrics.csv"), winslash = "/", mustWork = FALSE), "\n", sep = "")
}
