# Rebuild plots from saved, reloaded results; no analysis backend is rerun.
# Rscript scripts/verify-spatial-artifacts.R DATA_ROOT
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply DATA_ROOT")
root <- args[1]
pkgload::load_all(".", quiet = TRUE, compile = FALSE)
dirs <- list.dirs(root, recursive = FALSE, full.names = FALSE)
dirs <- dirs[file.exists(file.path(root, dirs, "acceptance.csv"))]
dirs <- dirs[order(file.info(file.path(root, dirs, "acceptance.csv"))$mtime)]
if (!length(dirs)) stop("No acceptance results found")
tables <- lapply(dirs, function(d) {
  x <- utils::read.csv(file.path(root, d, "acceptance.csv"), stringsAsFactors = FALSE)
  x$artifact_dir <- d
  x
})
tab <- do.call(rbind, tables)
tab <- tab[!duplicated(tab$case, fromLast = TRUE), ]
if (any(tab$status != "passed")) stop("The latest run of at least one case did not pass")
outdir <- file.path(root, "verified")
dir.create(outdir, showWarnings = FALSE)
for (i in seq_len(nrow(tab))) {
  name <- tab$case[i]
  file <- file.path(root, tab$artifact_dir[i], paste0(name, ".rds"))
  object <- readRDS(file)
  plot <- switch(name,
    visium = SpatialSpotPlot(object, group.by = "SpotQC", image = "slice1"),
    banksy = SpatialSpotPlot(object, group.by = "BANKSY_cluster", image = "slice1"),
    smoothclust = SpatialSpotPlot(object, group.by = "SmoothClust_cluster", image = "slice1"),
    spanorm = SpatialSpotPlot(object, features = rownames(object[["SpaNorm"]])[1], assay = "SpaNorm", layer = "data", image = "slice1"),
    spotsweeper = SpatialSpotPlot(object, group.by = "SpotSweeper_QC", image = "slice1"),
    hd = SpatialSpotPlot(object, group.by = "SpatialSketch_projected", image = "slice1.008um"),
    xenium = SpatialCellPlot(object, image = "fov", group.by = "segmentation_qc"),
    xenium_analysis = SpatialSpotPlot(object, group.by = "SpatialSketch_projected", overlay_image = FALSE),
    subjects = SpatialSamplePlot(object))
  stopifnot(inherits(plot, c("ggplot", "patchwork")))
  ggplot2::ggsave(file.path(outdir, paste0(name, ".png")), plot, width = 9, height = 6, dpi = 110)
  file.copy(file, file.path(outdir, basename(file)), overwrite = TRUE)
}
stopifnot(nrow(tab) == 9L)
tab$reload_plot <- "passed"
utils::write.csv(tab, file.path(root, "acceptance-final.csv"), row.names = FALSE)
sources <- c("visium-hd-tiny.zip", "xenium-tiny.zip", "diabetesData.rda", "xenium_human_pancreas_sub.rds")
manifest <- data.frame(file = sources, md5 = unname(tools::md5sum(file.path(root, sources))))
utils::write.csv(manifest, file.path(root, "source-checksums.csv"), row.names = FALSE)
