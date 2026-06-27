#!/usr/bin/env Rscript
# Expanded RunMonocle2 C++ backend validation:
# multi-state parity, parameter coverage, plot/tool compatibility, and edge behavior.

suppressPackageStartupMessages({
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools required", call. = FALSE)
  }
})

repo_root <- normalizePath(Sys.getenv("SCOP_ROOT", getwd()), wins = FALSE, mustWork = TRUE)
setwd(repo_root)
devtools::load_all(repo_root, quiet = TRUE, export_all = TRUE)

required <- c("Seurat", "SeuratObject", "Matrix", "monocle", "DDRTree")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

make_branching_srt <- function(n_genes = 90L, n_per_branch = 50L, seed = 1L) {
  set.seed(seed)
  n_cells <- n_per_branch * 3L
  branch <- rep(c("root", "branch_a", "branch_b"), each = n_per_branch)
  branch_time <- rep(seq(0, 1, length.out = n_per_branch), 3L)
  mu <- matrix(1, nrow = n_genes, ncol = n_cells)
  mu[1:20, branch == "root"] <- 5 + branch_time[branch == "root"] * 2
  mu[21:45, branch == "branch_a"] <- 3 + branch_time[branch == "branch_a"] * 8
  mu[46:70, branch == "branch_b"] <- 3 + branch_time[branch == "branch_b"] * 8
  mu[71:n_genes, ] <- 2
  counts <- matrix(
    rnbinom(n_genes * n_cells, mu = as.vector(mu), size = 2),
    nrow = n_genes,
    dimnames = list(paste0("g", seq_len(n_genes)), paste0("c", seq_len(n_cells)))
  )
  srt <- Seurat::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))
  srt$branch <- branch
  SeuratObject::VariableFeatures(srt) <- rownames(srt)
  srt
}

clone_counts_srt <- function(srt) {
  assay <- SeuratObject::DefaultAssay(srt)
  counts <- SeuratObject::GetAssayData(srt, assay = assay, layer = "counts")
  out <- Seurat::CreateSeuratObject(counts = counts, meta.data = srt@meta.data, assay = assay)
  vf <- SeuratObject::VariableFeatures(srt, assay = assay)
  if (length(vf)) {
    SeuratObject::VariableFeatures(out, assay = assay) <- vf
  }
  out
}

timed_run <- function(expr) {
  gc(verbose = FALSE, reset = TRUE)
  time <- system.time(value <- force(expr))
  list(
    value = value,
    elapsed = unname(time[["elapsed"]]),
    user = unname(time[["user.self"]]),
    sys = unname(time[["sys.self"]]),
    gc_vcell_peak_mb = gc(verbose = FALSE)[2, 7]
  )
}

compare_outputs <- function(out_r, out_cpp) {
  pt_r <- as.numeric(out_r$Monocle2_Pseudotime)
  pt_cpp <- as.numeric(out_cpp$Monocle2_Pseudotime)
  st_r <- as.character(out_r$Monocle2_State)
  st_cpp <- as.character(out_cpp$Monocle2_State)
  finite <- is.finite(pt_r) & is.finite(pt_cpp)
  data.frame(
    cells = length(pt_r),
    states_r = length(unique(st_r)),
    states_cpp = length(unique(st_cpp)),
    pt_max_abs_diff = max(abs(pt_r[finite] - pt_cpp[finite])),
    pt_mean_abs_diff = mean(abs(pt_r[finite] - pt_cpp[finite])),
    pt_cor = suppressWarnings(stats::cor(pt_r, pt_cpp, use = "complete.obs")),
    pt_spearman = suppressWarnings(stats::cor(pt_r, pt_cpp, method = "spearman", use = "complete.obs")),
    state_match_pct = 100 * mean(st_r == st_cpp, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

run_pair <- function(name, srt, args) {
  cat("\n== ", name, " ==\n", sep = "")
  r <- tryCatch(
    timed_run(do.call(RunMonocle2, c(list(srt = clone_counts_srt(srt), backend = "r"), args))),
    error = function(e) e
  )
  cpp <- tryCatch(
    timed_run(do.call(RunMonocle2, c(list(srt = clone_counts_srt(srt), backend = "cpp"), args))),
    error = function(e) e
  )
  if (inherits(r, "error") || inherits(cpp, "error")) {
    out <- data.frame(
      scenario = name,
      r_status = if (inherits(r, "error")) paste("error:", conditionMessage(r)) else "ok",
      cpp_status = if (inherits(cpp, "error")) paste("error:", conditionMessage(cpp)) else "ok",
      stringsAsFactors = FALSE
    )
    print(out, row.names = FALSE)
    return(out)
  }
  cmp <- compare_outputs(r$value, cpp$value)
  perf <- data.frame(
    scenario = name,
    time_r = r$elapsed,
    time_cpp = cpp$elapsed,
    speedup = r$elapsed / cpp$elapsed,
    mem_r_mb = r$gc_vcell_peak_mb,
    mem_cpp_mb = cpp$gc_vcell_peak_mb,
    mem_delta_mb = cpp$gc_vcell_peak_mb - r$gc_vcell_peak_mb,
    stringsAsFactors = FALSE
  )
  out <- cbind(perf, cmp)
  print(out, row.names = FALSE)
  out
}

branching <- make_branching_srt()
branching_args <- list(
  features = rownames(branching)[1:80],
  root_state = 1,
  expressionFamily = "uninormal",
  norm_method = "none",
  ddrtree_maxIter = 5,
  verbose = FALSE
)

results <- list(
  run_pair("multi_state_explicit_features", branching, branching_args),
  run_pair("multi_state_hvf_features", branching, modifyList(branching_args, list(features = NULL, feature_type = "HVF"))),
  run_pair("multi_state_group_root", branching, modifyList(branching_args, list(group.by = "branch", root_state = "branch_a")))
)

cat("\n== cpp plot/tool compatibility ==\n")
plot_out <- RunMonocle2(
  clone_counts_srt(branching),
  features = rownames(branching)[1:70],
  backend = "cpp",
  group.by = "branch",
  root_state = 1,
  expressionFamily = "uninormal",
  norm_method = "none",
  ddrtree_maxIter = 5,
  show_plot = TRUE,
  verbose = FALSE
)
plot_checks <- data.frame(
  has_reduction = "DDRTree" %in% names(plot_out@reductions),
  has_state = "Monocle2_State" %in% colnames(plot_out@meta.data),
  has_pseudotime = "Monocle2_Pseudotime" %in% colnames(plot_out@meta.data),
  graph_rows = nrow(plot_out@tools$Monocle2$graph),
  engine = plot_out@tools$Monocle2$parameters$engine,
  stringsAsFactors = FALSE
)
print(plot_checks, row.names = FALSE)

cat("\n== cpp missing root_state behavior ==\n")
missing_root <- tryCatch(
  RunMonocle2(
    clone_counts_srt(branching),
    features = rownames(branching)[1:70],
    backend = "cpp",
    root_state = "missing_state",
    expressionFamily = "uninormal",
    norm_method = "none",
    ddrtree_maxIter = 5,
    verbose = FALSE
  ),
  error = function(e) e
)
cat(if (inherits(missing_root, "error")) conditionMessage(missing_root) else "unexpected success", "\n")

real_rds <- Sys.getenv("MONOCLE2_CPP_REAL_RDS", "")
if (nzchar(real_rds) && file.exists(real_rds)) {
  cat("\n== real RDS ==\n")
  real_srt <- readRDS(real_rds)
  real_features <- SeuratObject::VariableFeatures(real_srt)
  if (!length(real_features)) {
    real_features <- rownames(real_srt)
  }
  real_args <- list(
    features = real_features[seq_len(min(500L, length(real_features)))],
    root_state = 1,
    verbose = FALSE,
    show_plot = FALSE
  )
  results[[length(results) + 1L]] <- run_pair(basename(real_rds), real_srt, real_args)
} else {
  cat("\n(real RDS skipped; set MONOCLE2_CPP_REAL_RDS=/path/to/object.rds)\n")
}

invisible(results)
