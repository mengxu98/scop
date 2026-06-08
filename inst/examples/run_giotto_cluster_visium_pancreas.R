# Minimal smoke test for scop::RunGiottoCluster() on official scop data.
#
# Giotto has several optional/suggested dependencies that are not needed for
# nearest-network clustering. The default RunGiottoCluster() path uses
# Giotto's R/igraph Leiden implementation, so no Giotto Python environment is
# required. Keep dependencies restricted for a lightweight install.
if (!requireNamespace("Giotto", quietly = TRUE)) {
  if (!requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak", repos = "https://cloud.r-project.org")
  }
  pak::pkg_install(
    "giotto-suite/Giotto",
    ask = FALSE,
    upgrade = FALSE,
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}

# Keep BLAS thread allocation modest on Windows smoke-test machines.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)

library(scop)

set.seed(1)
data(visium_human_pancreas_sub, package = "scop")

# Full official spatial example: 5000 genes x 1986 spots.
spatial <- visium_human_pancreas_sub

# For a faster smoke test, uncomment:
# spatial <- subset(spatial, cells = colnames(spatial)[1:400])

spatial <- Seurat::NormalizeData(
  spatial,
  assay = "Spatial",
  verbose = FALSE
)

spatial <- RunGiottoCluster(
  spatial,
  assay = "Spatial",
  features = rownames(spatial)[1:200],
  dims = 1:10,
  k = 10,
  resolution = 0.4,
  verbose = TRUE
)

print(table(spatial[["Giotto_cluster"]][, 1], useNA = "ifany"))
print(head(spatial@tools[["GiottoCluster"]][["clusters"]]))

p <- SpatialSpotPlot(
  spatial,
  group.by = "Giotto_cluster",
  overlay_image = FALSE
)
print(p)

saveRDS(spatial, file = "run_giotto_cluster_visium_pancreas_result.rds")
