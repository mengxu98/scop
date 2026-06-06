# Co-embed reference and query cells

Impute reference expression into query cells with `TransferData`, merge
reference and query cells, and compute PCA/UMAP on the shared expression
space. The current implementation mirrors the Seurat ATAC integration
vignette for ATAC query and RNA reference.

## Usage

``` r
RunCoEmbedding(
  srt,
  reference,
  assay = NULL,
  reference_assay = NULL,
  reference_reduction = "pca",
  reference_dims = 1:30,
  gene_activity_assay = "ACTIVITY",
  weight_reduction = NULL,
  dims = 2:30,
  genes_use = NULL,
  imputed_assay = "RNA",
  modality_col = "modality",
  coembed_prefix = "CoEmbed",
  npcs = 30,
  umap_dims = 1:30,
  k.weight = 100,
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- srt:

  A Seurat object.

- reference:

  RNA reference `Seurat` object used for label transfer.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- reference_assay:

  Assay used in the reference object.

- reference_reduction:

  Reduction used in the reference object.

- reference_dims:

  Dimensions used from the reference reduction.

- gene_activity_assay:

  Name of the gene activity assay used for mapping.

- weight_reduction:

  Reduction in `srt` used to weight transferred labels. If `NULL`, an
  ATAC linear reduction is resolved automatically from
  `ATAC_default_linear_reduction`, `{prefix}lsi`, `{prefix}svd`, or the
  current default reduction.

- dims:

  Query reduction dimensions used by `TransferData`.

- genes_use:

  Genes/features used for expression imputation and co-embedding. If
  `NULL`, reference variable features are used.

- imputed_assay:

  Assay name for imputed RNA expression in ATAC cells.

- modality_col:

  Metadata column storing RNA/ATAC origin after merge.

- coembed_prefix:

  Prefix for co-embedding reductions.

- npcs:

  Number of PCs to compute.

- umap_dims:

  PC dimensions used for UMAP.

- k.weight:

  Number of neighbors used when weighting transfer anchors.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed used for UMAP.

## Value

A merged `Seurat` object containing RNA reference and ATAC query cells.

## Examples

``` r
data(pbmcmultiome_sub)
pbmcmultiome_sub <- standard_scop(
  pbmcmultiome_sub,
  assay = c("RNA", "peaks"),
  linear_reduction_dims = 20
)
coembed <- RunCoEmbedding(
  srt = pbmcmultiome_sub,
  reference = pbmcmultiome_sub,
  assay = "peaks",
  reference_assay = "RNA",
  gene_activity_assay = "RNA",
  reference_reduction = "RNApca",
  reference_dims = 1:10,
  dims = 2:10,
  umap_dims = 1:10
)

CellDimPlot(
  coembed,
  group.by = c("modality", "CellType"),
  xlab = "CoEmbedUMAP_1",
  ylab = "CoEmbedUMAP_2"
)
```
