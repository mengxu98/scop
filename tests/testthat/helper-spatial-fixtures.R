make_spatial_test_grid <- function() {
  counts <- matrix(
    c(3, 1, 0, 2, 0, 4, 1, 0, 2, 1, 3, 0),
    nrow = 3,
    dimnames = list(paste0("gene", 1:3), paste0("spot", 1:4))
  )
  srt <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  srt$col <- c(0, 1, 0, 1)
  srt$row <- c(0, 0, 1, 1)
  srt
}

make_spatial_test_object <- function() {
  counts <- matrix(
    c(
      5, 4, 0, 0, 0,
      0, 0, 4, 5, 4,
      1, 1, 1, 1, 1,
      3, 0, 3, 0, 3
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("Gene", 1:4), paste0("Spot", 1:5))
  )
  counts <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
  srt <- Seurat::CreateSeuratObject(counts = counts)
  srt$x <- c(1, 2, 3, 4, 5)
  srt$y <- c(1, 1, 2, 2, 3)
  srt
}

make_visium_test_subset <- function(n, seed = 42) {
  data(visium_human_pancreas_sub)
  srt <- visium_human_pancreas_sub
  srt <- Seurat::NormalizeData(srt, verbose = FALSE)
  srt@images <- list()
  set.seed(seed)
  cells <- sample(colnames(srt), n)
  srt <- srt[, cells]
  srt
}

make_panc8_test_subset <- function(n, seed = 7) {
  data(panc8_sub)
  set.seed(seed)
  panc8_sub[, sample(colnames(panc8_sub), n)]
}
