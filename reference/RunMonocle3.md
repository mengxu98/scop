# Run Monocle3 analysis

Run Monocle3 analysis

## Usage

``` r
RunMonocle3(
  srt,
  group.by = NULL,
  assay = NULL,
  layer = "counts",
  reduction = NULL,
  clusters = NULL,
  graph = NULL,
  partition_qval = 0.05,
  k = 50,
  cluster_method = "louvain",
  num_iter = 2,
  resolution = NULL,
  use_partition = NULL,
  close_loop = TRUE,
  root_pr_nodes = NULL,
  root_cells = NULL,
  seed = 11,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- group.by:

  Name of one or more meta.data columns to group (color) cells by.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used.

- layer:

  Which layer to use. Default is `data`.

- reduction:

  Which dimensionality reduction to use. If not specified, will use the
  reduction returned by
  [DefaultReduction](https://mengxu98.github.io/scop/reference/DefaultReduction.md).

- clusters:

  The cluster variable in the Seurat object to use for analysis.
  Defaults to NULL, in which case use Monocle clusters is used.

- graph:

  The name of the graph slot in the Seurat object to use for analysis.
  Defaults to NULL, in which case Monocle graph is used.

- partition_qval:

  The q-value threshold for partitioning cells. Defaults to 0.05.

- k:

  The number of nearest neighbors to consider for clustering. Defaults
  to 50.

- cluster_method:

  The clustering method to use. Defaults to "louvain".

- num_iter:

  The number of iterations for clustering. Defaults to 2.

- resolution:

  The resolution parameter for clustering. Defaults to NULL.

- use_partition:

  Whether to use partitions to learn disjoint graph in each partition.
  If not specified, user will be prompted for input. Defaults to NULL.

- close_loop:

  Whether to close loops in the graph. Defaults to TRUE.

- root_pr_nodes:

  The root nodes to order cells. If not specified, user will be prompted
  for input. Defaults to NULL.

- root_cells:

  The root cells to order cells. If not specified, user will be prompted
  for input. Defaults to NULL.

- seed:

  Random seed for reproducibility. Default is `11`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Examples

``` r
if (interactive()) {
  data(pancreas_sub)
  pancreas_sub <- standard_scop(pancreas_sub)
  pancreas_sub <- RunMonocle3(
    pancreas_sub,
    reduction = "UMAP"
  )

  pancreas_sub <- RunMonocle3(
    pancreas_sub,
    reduction = "UMAP",
    group.by = "CellType"
  )
  names(pancreas_sub@tools$Monocle3)
  trajectory <- pancreas_sub@tools$Monocle3$trajectory
  milestones <- pancreas_sub@tools$Monocle3$milestones

  CellDimPlot(
    pancreas_sub,
    group.by = "Monocle3_partitions",
    reduction = "UMAP",
    label = TRUE,
    theme_use = "theme_blank"
  ) +
    trajectory +
    milestones
  CellDimPlot(
    pancreas_sub,
    group.by = "Monocle3_clusters",
    reduction = "UMAP",
    label = TRUE,
    theme_use = "theme_blank"
  ) +
    trajectory
  FeatureDimPlot(
    pancreas_sub,
    features = "Monocle3_Pseudotime",
    reduction = "UMAP",
    theme_use = "theme_blank"
  ) +
    trajectory

  if (FALSE) {
    # Select the lineage using monocle3::choose_graph_segments
    cds <- pancreas_sub@tools$Monocle3$cds
    cds_sub <- thisutils::get_namespace_fun(
      "monocle3", "choose_graph_segments"
    )(
      cds,
      starting_pr_node = NULL,
      ending_pr_nodes = NULL
    )
    pancreas_sub$Lineages_1 <- NA
    pancreas_sub$Lineages_1[colnames(
      cds_sub
    )] <- pancreas_sub$Monocle3_Pseudotime[colnames(cds_sub)]
    CellDimPlot(
      pancreas_sub,
      group.by = "SubCellType",
      lineages = "Lineages_1",
      lineages_span = 0.1,
      theme_use = "theme_blank"
    )
  }

  # Use Seurat clusters to infer the trajectories
  pancreas_sub <- standard_scop(pancreas_sub)
  CellDimPlot(
    pancreas_sub,
    group.by = c("Standardclusters", "CellType"),
    label = TRUE,
    theme_use = "theme_blank"
  )

  pancreas_sub <- RunMonocle3(
    pancreas_sub,
    clusters = "Standardclusters"
  )

  trajectory <- pancreas_sub@tools$Monocle3$trajectory
  CellDimPlot(
    pancreas_sub,
    group.by = "Monocle3_partitions",
    reduction = "StandardUMAP2D",
    label = TRUE, theme_use = "theme_blank"
  ) + trajectory

  CellDimPlot(
    pancreas_sub,
    group.by = "Monocle3_clusters",
    reduction = "StandardUMAP2D",
    label = TRUE, theme_use = "theme_blank"
  ) + trajectory
  FeatureDimPlot(
    pancreas_sub,
    features = "Monocle3_Pseudotime",
    reduction = "StandardUMAP2D",
    theme_use = "theme_blank"
  ) + trajectory

  # Use custom graphs and cell clusters to infer
  # the partitions and trajectories, respectively
  pancreas_sub <- standard_scop(
    pancreas_sub,
    cluster_resolution = 5
  )
  CellDimPlot(
    pancreas_sub,
    group.by = c("Standardclusters", "CellType"),
    label = TRUE
  )
  pancreas_sub <- RunMonocle3(
    pancreas_sub,
    clusters = "Standardclusters",
    graph = "Standardpca_SNN"
  )
  trajectory <- pancreas_sub@tools$Monocle3$trajectory
  CellDimPlot(
    pancreas_sub,
    group.by = "Monocle3_partitions",
    reduction = "StandardUMAP2D",
    label = TRUE, theme_use = "theme_blank"
  ) + trajectory
  CellDimPlot(
    pancreas_sub,
    group.by = "Monocle3_clusters",
    reduction = "StandardUMAP2D",
    label = TRUE, theme_use = "theme_blank"
  ) + trajectory
  FeatureDimPlot(
    pancreas_sub,
    features = "Monocle3_Pseudotime",
    reduction = "StandardUMAP2D",
    theme_use = "theme_blank"
  ) + trajectory
}
```
