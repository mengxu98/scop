suppressPackageStartupMessages(library(Seurat))

if (!requireNamespace("SeuratData", quietly = TRUE)) {
  stop("Install SeuratData before rebuilding visium_mouse_brain_slices_sub.", call. = FALSE)
}

if (!requireNamespace("stxBrain.SeuratData", quietly = TRUE)) {
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  # stxBrain.SeuratData is about 136 MB, so the default 60 second timeout is
  # often too short for a one-time data rebuild.
  options(timeout = max(600, old_timeout))
  install.packages(
    "stxBrain.SeuratData",
    repos = "https://seurat.nygenome.org",
    type = "source"
  )
}

set.seed(98)

load_stx_slice <- function(type, nspots = 1000) {
  obj <- SeuratData::LoadData("stxBrain", type = type)
  obj$sample <- type

  coords <- SeuratObject::GetTissueCoordinates(obj)
  rownames(coords) <- coords$cell
  obj$x <- coords[colnames(obj), "x"]
  obj$y <- coords[colnames(obj), "y"]

  cells <- sample(colnames(obj), size = min(nspots, ncol(obj)))
  obj[, cells]
}

brain1 <- load_stx_slice("anterior1")
brain2 <- load_stx_slice("anterior2")

visium_mouse_brain_slices_sub <- merge(
  brain1,
  y = brain2,
  add.cell.ids = c("anterior1", "anterior2")
)

joined <- SeuratObject::JoinLayers(visium_mouse_brain_slices_sub[["Spatial"]])
counts <- SeuratObject::GetAssayData(joined, layer = "counts")
feature_totals <- Matrix::rowSums(counts)
features <- names(sort(feature_totals, decreasing = TRUE))[
  seq_len(min(4000L, length(feature_totals)))
]

visium_mouse_brain_slices_sub <- visium_mouse_brain_slices_sub[features, ]
visium_mouse_brain_slices_sub[["Spatial"]] <- SeuratObject::JoinLayers(
  visium_mouse_brain_slices_sub[["Spatial"]]
)
visium_mouse_brain_slices_sub$sample <- factor(
  visium_mouse_brain_slices_sub$sample,
  levels = c("anterior1", "anterior2")
)

save(
  visium_mouse_brain_slices_sub,
  file = file.path("data", "visium_mouse_brain_slices_sub.rda"),
  compress = "xz"
)
