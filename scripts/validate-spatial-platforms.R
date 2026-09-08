# Run from the package root, after compiling with pkgload::load_all().
# Rscript scripts/validate-spatial-platforms.R DATA_ROOT OUTPUT_DIR [CASE ...]
# File provenance and downloads are documented in spatial-platform-data.md.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Supply DATA_ROOT and OUTPUT_DIR")
root <- normalizePath(".", winslash = "/", mustWork = TRUE)
data_root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[2], winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cases <- list(
  visium = function(data_root) {
    data("visium_human_pancreas_sub", package = "scop")
    raw <- visium_human_pancreas_sub
    cells <- colnames(raw)[unique(round(seq(1, ncol(raw), length.out = 120)))]
    features <- names(sort(Matrix::rowSums(GetAssayData5(raw, assay = "Spatial")), decreasing = TRUE))[1:200]
    x <- subset(raw, cells = cells, features = features)
    out <- RunStandardWorkflow(x, workflow = "spatial", assay = "Spatial", image = "slice1",
      spot_qc_params = list(qc_metrics = c("umi", "gene")), nHVF = 80,
      spatial_variable_features_params = list(nfeatures = 20, min_spots = 3),
      linear_reduction_dims = 5, linear_reduction_dims_use = 1:5,
      nonlinear_reduction = character(), verbose = FALSE)
    stopifnot(identical(colnames(out), colnames(x)), identical(out@tools$run_standard_spatial_workflow$status, "completed"))
    list(object = out, plot = SpatialSpotPlot(out, group.by = "SpotQC", image = "slice1"),
      detail = sprintf("%s real Visium spots; QC + preprocessing + native SVG", ncol(out)))
  },
  banksy = function(data_root) {
    data("visium_human_pancreas_sub", package = "scop")
    x <- visium_human_pancreas_sub[, 1:120]
    out <- RunStandardWorkflow(x, workflow = "spatial", assay = "Spatial", image = "slice1",
      do_spot_qc = FALSE, do_spatial_variable_features = FALSE,
      do_spatial_cluster = TRUE, spatial_cluster_method = "BANKSY",
      spatial_cluster_params = list(features = rownames(x)[1:80], npcs = 5,
        k_geom = 6, k_neighbors = 10, algo = "louvain"),
      nHVF = 80, linear_reduction_dims = 5, linear_reduction_dims_use = 1:5,
      nonlinear_reduction = character(), verbose = FALSE)
    stopifnot(!anyNA(out$BANKSY_cluster))
    list(object = out, plot = SpatialSpotPlot(out, group.by = "BANKSY_cluster", image = "slice1"),
      detail = "120 real Visium spots; installed BANKSY producer through standard workflow")
  },
  smoothclust = function(data_root) {
    data("visium_human_pancreas_sub", package = "scop")
    x <- visium_human_pancreas_sub[, 1:120]
    out <- RunStandardWorkflow(x, workflow = "spatial", assay = "Spatial", image = "slice1",
      do_spot_qc = FALSE, do_spatial_variable_features = FALSE, do_spatial_cluster = TRUE,
      spatial_cluster_method = "SmoothClust", spatial_q = 3,
      spatial_cluster_params = list(nfeatures = 80, n_pcs = 5, smooth_method = "knn", k = 6),
      nHVF = 80, linear_reduction_dims = 5, linear_reduction_dims_use = 1:5,
      nonlinear_reduction = character(), verbose = FALSE)
    stopifnot(!anyNA(out$SmoothClust_cluster), length(unique(out$SmoothClust_cluster)) == 3)
    list(object = out, plot = SpatialSpotPlot(out, group.by = "SmoothClust_cluster", image = "slice1"),
      detail = "120 real Visium spots; installed smoothclust producer through standard workflow")
  },
  spanorm = function(data_root) {
    data("visium_human_pancreas_sub", package = "scop")
    raw <- visium_human_pancreas_sub
    features <- names(sort(Matrix::rowSums(GetAssayData5(raw, assay = "Spatial")), decreasing = TRUE))[1:50]
    x <- subset(raw, cells = colnames(raw)[unique(round(seq(1, ncol(raw), length.out = 240)))], features = features)
    counts <- GetAssayData5(x, assay = "Spatial")
    out <- RunStandardWorkflow(x, workflow = "spatial", assay = "Spatial", image = "slice1",
      do_spot_qc = FALSE, normalization_method = "SpaNorm", spanorm_params = list(sample.p = 1),
      nHVF = 30, linear_reduction_dims = 5, linear_reduction_dims_use = 1:5,
      spatial_variable_features_params = list(nfeatures = 10, min_spots = 3),
      nonlinear_reduction = character(), verbose = FALSE)
    stopifnot(identical(counts, GetAssayData5(out, assay = "Spatial")),
      identical(out@tools$run_standard_spatial_workflow$parameters$analysis_assay, "SpaNorm"))
    list(object = out, plot = SpatialSpotPlot(out, features = rownames(out[["SpaNorm"]])[1],
      assay = "SpaNorm", layer = "data", image = "slice1"),
      detail = "240 real Visium spots; live SpaNorm + preprocessing + SVG; original counts identical")
  },
  spotsweeper = function(data_root) {
    data("visium_human_pancreas_sub", package = "scop")
    x <- visium_human_pancreas_sub[, 1:120]
    out <- RunStandardWorkflow(x, workflow = "spatial", assay = "Spatial", image = "slice1",
      do_spot_qc = FALSE, do_spatial_qc = TRUE,
      spatial_qc_params = list(run_artifact = FALSE, n_neighbors = 6),
      do_spatial_variable_features = FALSE, nHVF = 80, linear_reduction_dims = 5,
      linear_reduction_dims_use = 1:5, nonlinear_reduction = character(), verbose = FALSE)
    stopifnot(identical(colnames(out), colnames(x)), identical(out@tools$SpotSweeper$status, "completed"))
    list(object = out, plot = SpatialSpotPlot(out, group.by = "SpotSweeper_QC", image = "slice1"),
      detail = "120 real Visium spots; live local SpotSweeper QC; artifact detection explicitly off")
  },
  hd = function(data_root) {
    x <- ReadSpatialData(file.path(data_root, "hd"), "visium_hd", sample_id = "tiny_mouse")
    counts <- GetAssayData5(x, assay = "Spatial.008um")
    out <- RunSpatialSketch(x, assay = "Spatial.008um", image = "slice1.008um",
      ncells = 300, nfeatures = 150, npcs = 8, method = "Uniform", verbose = FALSE)
    stopifnot(identical(counts, GetAssayData5(out, assay = "Spatial.008um")),
      length(out@tools$SpatialSketch$projected_cells) == ncol(counts),
      identical(out@misc$scop_spatial_input[["slice1.016um"]]$resolution_um, 16))
    list(object = out, plot = SpatialSpotPlot(out, group.by = "SpatialSketch_projected", image = "slice1.008um"),
      detail = sprintf("10x developer HD fixture: 300 sampled, %s full bins; both 8/16 um assays retained", ncol(counts)))
  },
  xenium = function(data_root) {
    x <- ReadSpatialData(file.path(data_root, "xenium"), "xenium", sample_id = "tiny_ileum")
    qc <- SpatialSegmentationQC(x, assay = "Xenium", image = "fov", min_counts = 1)
    stopifnot(identical(qc$cell_id, colnames(x)), all(is.finite(qc$area)), all(qc$area > 0),
      identical(SpatialDataInfo(x, assay = "Xenium", image = "fov")$coordinate_units, "micron"))
    x$segmentation_qc <- stats::setNames(qc$status, qc$cell_id)
    x@misc$segmentation_qc <- qc
    list(object = x, plot = SpatialCellPlot(x, image = "fov", group.by = "segmentation_qc"),
      detail = sprintf("%s real Xenium cells; vendor loader, true polygons and count/area QC", ncol(x)))
  },
  xenium_analysis = function(data_root) {
    x <- readRDS(file.path(data_root, "xenium_human_pancreas_sub.rds"))
    counts <- GetAssayData5(x, assay = "Xenium")
    out <- RunSpatialSketch(x, assay = "Xenium", ncells = 300, nfeatures = 100,
      npcs = 8, method = "LeverageScore", resolution = 1, verbose = FALSE)
    stopifnot(identical(counts, GetAssayData5(out, assay = "Xenium")), !anyNA(out$SpatialSketch_projected))
    list(object = out, plot = SpatialSpotPlot(out, group.by = "SpatialSketch_projected", overlay_image = FALSE),
      detail = sprintf("%s real Xenium pancreas cells; LeverageScore sketch/full annotation; centroid-only fixture", ncol(x)))
  },
  subjects = function(data_root) {
    env <- new.env()
    load(file.path(data_root, "diabetesData.rda"), envir = env)
    meta <- as.data.frame(SummarizedExperiment::colData(env$diabetesData))
    rownames(meta) <- meta$imageCellID
    summary <- SpatialSampleSummary(meta, "cellType", "imageID", "stage", "case")
    out <- SpatialSampleComparison(summary, c("Non-diabetic", "Onset"))
    stopifnot(all(out$comparisons$n_reference == 4), all(out$comparisons$n_treatment == 4),
      all(out$comparisons$status %in% c("tested", "not_tested")))
    list(object = out, plot = SpatialSamplePlot(out),
      detail = sprintf("Real Damond IMC design; %s source subjects; contrast uses 4+4 patients, not cells or images as replicates", length(unique(meta$case))))
  }
)

backends <- c(banksy = "Banksy", smoothclust = "smoothclust", spanorm = "SpaNorm", spotsweeper = "SpotSweeper")
selected <- if (length(args) > 2L) args[-c(1, 2)] else names(cases)
if (!all(selected %in% names(cases))) stop("Unknown case")
records <- list()
for (name in selected) {
  backend <- backends[name]
  if (!is.na(backend) && !nzchar(system.file(package = backend))) {
    records[[name]] <- data.frame(case = name, status = "unavailable", seconds = 0,
      sampled_peak_mb = NA_real_, detail = paste("Missing optional backend", backend))
    next
  }
  cat("RUN", name, "\n")
  started <- proc.time()[["elapsed"]]
  process <- callr::r_bg(function(fun, root, data_root, output_dir, name) {
    suppressPackageStartupMessages(pkgload::load_all(root, quiet = TRUE, compile = FALSE))
    value <- withCallingHandlers(fun(data_root), error = function(e) print(sys.calls()))
    file <- file.path(output_dir, paste0(name, ".rds"))
    saveRDS(value$object, file)
    reloaded <- readRDS(file)
    stopifnot(identical(value$object, reloaded))
    ggplot2::ggsave(file.path(output_dir, paste0(name, ".png")), value$plot,
      width = 9, height = 6, dpi = 110)
    writeLines(capture.output(sessionInfo()), file.path(output_dir, paste0(name, "-session.txt")))
    value$detail
  }, args = list(fun = cases[[name]], root = root, data_root = data_root,
    output_dir = output_dir, name = name),
    stdout = file.path(output_dir, paste0(name, "-stdout.log")),
    stderr = file.path(output_dir, paste0(name, "-stderr.log")))
  peak <- 0
  while (process$is_alive()) {
    rss <- tryCatch({
      h <- ps::ps_handle(process$get_pid())
      handles <- c(list(h), ps::ps_children(h, recursive = TRUE))
      sum(vapply(handles, function(x) tryCatch(unname(ps::ps_memory_info(x)["rss"]), error = function(e) 0), numeric(1)))
    }, error = function(e) 0)
    peak <- max(peak, rss)
    if (proc.time()[["elapsed"]] - started > 600) {
      process$kill_tree(); break
    }
    Sys.sleep(0.25)
  }
  result <- tryCatch(process$get_result(), error = identity)
  records[[name]] <- data.frame(case = name,
    status = if (inherits(result, "error")) "failed" else "passed",
    seconds = proc.time()[["elapsed"]] - started, sampled_peak_mb = peak / 1024^2,
    detail = if (inherits(result, "error")) conditionMessage(result) else result)
  utils::write.csv(do.call(rbind, records), file.path(output_dir, "acceptance.csv"), row.names = FALSE)
  cat(records[[name]]$status, name, "\n")
}
if (any(vapply(records, function(x) x$status != "passed", logical(1)))) quit(status = 1L)
