#!/usr/bin/env Rscript
# Benchmark RunMonocle2 R vs C++ backends: speed, memory, consistency

suppressPackageStartupMessages({
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools required")
  }
})

repo_root <- normalizePath(Sys.getenv("SCOP_ROOT", getwd()), wins = FALSE, mustWork = TRUE)
setwd(repo_root)

cat("Loading scop from source...\n")
devtools::load_all(repo_root, quiet = TRUE, export_all = TRUE)

required <- c("Seurat", "monocle", "DDRTree", "SeuratObject", "Matrix")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing packages: ", paste(missing, collapse = ", "))
}

benchmark_run <- function(expr) {
  gc(verbose = FALSE, reset = TRUE)
  t <- system.time(result <- expr)
  mem_peak_mb <- gc(verbose = FALSE)[2, 7]
  list(
    elapsed = unname(t["elapsed"]),
    user = unname(t["user"]),
    mem_peak_mb = mem_peak_mb,
    result = result
  )
}

compare_results <- function(out_r, out_cpp, label) {
  pt_r <- as.numeric(out_r$Monocle2_Pseudotime)
  pt_cpp <- as.numeric(out_cpp$Monocle2_Pseudotime)
  st_r <- as.character(out_r$Monocle2_State)
  st_cpp <- as.character(out_cpp$Monocle2_State)

  pt_max_diff <- if (all(is.finite(pt_r)) && all(is.finite(pt_cpp))) {
    max(abs(pt_r - pt_cpp))
  } else {
    NA_real_
  }
  pt_cor <- if (sum(is.finite(pt_r) & is.finite(pt_cpp)) >= 2L) {
    suppressWarnings(stats::cor(pt_r, pt_cpp, use = "complete.obs"))
  } else {
    NA_real_
  }
  state_match <- mean(st_r == st_cpp, na.rm = TRUE)

  cat(sprintf(
    "\n[%s] Consistency (root_state=1):\n  pseudotime max|diff| = %.2e, cor = %.6f\n  state match rate = %.1f%%\n",
    label, pt_max_diff, pt_cor, 100 * state_match
  ))

  list(
    pt_max_diff = pt_max_diff,
    pt_cor = pt_cor,
    state_match = state_match,
    identical_pt = isTRUE(all.equal(pt_r, pt_cpp, tolerance = 1e-8)),
    identical_state = identical(st_r, st_cpp)
  )
}

make_synthetic <- function(n_genes, n_cells, seed = 1L) {
  set.seed(seed)
  counts <- matrix(
    stats::rnbinom(n_genes * n_cells, mu = 5, size = 1),
    nrow = n_genes,
    dimnames = list(paste0("g", seq_len(n_genes)), paste0("c", seq_len(n_cells)))
  )
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  srt <- Seurat::CreateSeuratObject(counts = counts)
  SeuratObject::VariableFeatures(srt) <- rownames(srt)
  srt
}

clone_srt <- function(srt) {
  counts <- SeuratObject::GetAssayData(srt, layer = "counts")
  out <- Seurat::CreateSeuratObject(
    counts = counts,
    meta.data = srt@meta.data
  )
  vf <- SeuratObject::VariableFeatures(srt)
  if (length(vf)) {
    SeuratObject::VariableFeatures(out) <- vf
  }
  out
}

run_scenario <- function(name, srt, n_features, root_state = 1L, ddrtree_maxIter = NULL) {
  features <- rownames(srt)[seq_len(min(n_features, nrow(srt)))]
  n_cells <- ncol(srt)
  n_genes <- nrow(srt)

  cat(sprintf(
    "\n========== %s (%d genes x %d cells, %d features) ==========\n",
    name, n_genes, n_cells, length(features)
  ))

  run_args <- list(
    features = features,
    root_state = root_state,
    verbose = FALSE,
    ddrtree_maxIter = ddrtree_maxIter
  )

  bench_r <- benchmark_run({
    do.call(RunMonocle2, c(list(srt = clone_srt(srt)), run_args, list(backend = "r")))
  })
  bench_cpp <- benchmark_run({
    do.call(RunMonocle2, c(list(srt = clone_srt(srt)), run_args, list(backend = "cpp")))
  })
  out_r <- bench_r$result
  out_cpp <- bench_cpp$result

  speedup <- bench_r$elapsed / bench_cpp$elapsed
  cat(sprintf(
    "  R backend:   %.2f s (user %.2f s, peak mem %.0f MB)\n",
    bench_r$elapsed, bench_r$user, bench_r$mem_peak_mb
  ))
  cat(sprintf(
    "  C++ backend: %.2f s (user %.2f s, peak mem %.0f MB)\n",
    bench_cpp$elapsed, bench_cpp$user, bench_cpp$mem_peak_mb
  ))
  cat(sprintf(
    "  Speedup (R/C++): %.2fx  |  Peak mem delta (C++ - R): %+.0f MB\n",
    speedup, bench_cpp$mem_peak_mb - bench_r$mem_peak_mb
  ))

  consistency <- compare_results(out_r, out_cpp, name)

  data.frame(
    scenario = name,
    n_genes = n_genes,
    n_cells = n_cells,
    n_features = length(features),
    time_r = bench_r$elapsed,
    time_cpp = bench_cpp$elapsed,
    speedup = speedup,
    mem_r_mb = bench_r$mem_peak_mb,
    mem_cpp_mb = bench_cpp$mem_peak_mb,
    mem_delta_mb = bench_cpp$mem_peak_mb - bench_r$mem_peak_mb,
    pt_max_diff = consistency$pt_max_diff,
    pt_cor = consistency$pt_cor,
    state_match_pct = 100 * consistency$state_match,
    stringsAsFactors = FALSE
  )
}

results <- list()

# Small: matches unit test size
results[[1L]] <- run_scenario(
  "small (test size)",
  make_synthetic(50L, 120L),
  n_features = 30L
)

# Medium: synthetic
results[[2L]] <- run_scenario(
  "medium (synthetic)",
  make_synthetic(200L, 500L),
  n_features = 100L
)

# Real data if available
panc <- tryCatch({
  data(pancreas_sub, package = "scop", envir = environment())
  get("pancreas_sub", envir = environment())
}, error = function(e) NULL)
if (!is.null(panc)) {
  if (!length(SeuratObject::VariableFeatures(panc))) {
    SeuratObject::VariableFeatures(panc) <- rownames(panc)
  }
  results[[length(results) + 1L]] <- run_scenario(
    "pancreas_sub (1000 cells)",
    panc,
    n_features = min(500L, nrow(panc))
  )
} else {
  cat("\n(pancreas_sub not found in package data, skipping)\n")
}

# Faster DDRTree for large exploratory timing
results[[length(results) + 1L]] <- run_scenario(
  "medium + ddrtree_maxIter=5",
  make_synthetic(200L, 500L),
  n_features = 100L,
  ddrtree_maxIter = 5L
)

summary_df <- do.call(rbind, results)
cat("\n\n========== SUMMARY ==========\n")
print(summary_df, row.names = FALSE)

cat("\nNote: Most runtime is in Monocle2 reduceDimension (DDRTree), shared by both backends.\n")
cat("C++ speedup applies mainly to orderCells / MST ordering step.\n")

# Ordering-only micro-benchmark (after shared reduceDimension)
cat("\n\n========== ORDERING-ONLY (shared reduceDimension, timed separately) ==========\n")
set.seed(1)
srt_order <- make_synthetic(200L, 500L)
features_order <- rownames(srt_order)[1:100]

prep_cds <- function(backend) {
  out <- RunMonocle2(
    clone_srt(srt_order),
    features = features_order,
    backend = backend,
    root_state = 1L,
    verbose = FALSE
  )
  out@tools$Monocle2$cds
}

cds_r <- prep_cds("r")
cds_cpp <- prep_cds("cpp")

order_r <- function() {
  get_namespace_fun("monocle", "orderCells")(cds_r, root_state = 1L)
}
order_cpp <- function() {
  scop:::monocle2_cpp_order_cells(cds_cpp, root_state = 1L, initial_root_state = 1L)
}

bench_order_r <- benchmark_run(order_r())
bench_order_cpp <- benchmark_run(order_cpp())
cat(sprintf(
  "  monocle orderCells (R):     %.3f s, peak mem %.0f MB\n",
  bench_order_r$elapsed, bench_order_r$mem_peak_mb
))
cat(sprintf(
  "  monocle2_cpp_order_cells:   %.3f s, peak mem %.0f MB\n",
  bench_order_cpp$elapsed, bench_order_cpp$mem_peak_mb
))
cat(sprintf(
  "  Ordering speedup (R/C++):   %.1fx\n",
  bench_order_r$elapsed / bench_order_cpp$elapsed
))
