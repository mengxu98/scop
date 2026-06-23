# Cell scoring

This function performs cell scoring on a Seurat object. It calculates
scores for a given set of features and stores them in `meta.data` and/or
a score assay, depending on `new_assay` and `store_metadata`.

## Usage

``` r
CellScoring(
  srt,
  features = NULL,
  layer = "data",
  assay = NULL,
  split.by = NULL,
  IDtype = "symbol",
  species = "Homo_sapiens",
  db = "GO_BP",
  termnames = NULL,
  db_update = FALSE,
  db_version = "latest",
  convert_species = TRUE,
  Ensembl_version = NULL,
  mirror = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  method = "Seurat",
  backend = c("cpp", "r"),
  cpp_strategy = c("sparse", "aucell", "topk", "full"),
  classification = TRUE,
  name = "",
  new_assay = FALSE,
  store_metadata = NULL,
  seed = 11,
  cores = 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- srt:

  A Seurat object.

- features:

  A named list of feature lists for scoring. If `NULL`, `db` will be
  used to create features sets.

- layer:

  Which layer to use. Default is `data`.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- split.by:

  Name of a column in meta.data column to split plot by. Default is
  `NULL`.

- IDtype:

  A character vector specifying the type of gene IDs in the `srt` object
  or `geneID` argument. This argument is used to convert the gene IDs to
  a different type if `IDtype` is different from `result_IDtype`.

- species:

  A character vector specifying the species for which the gene
  annotation databases should be prepared. Can be `"Homo_sapiens"` or
  `"Mus_musculus"`.

- db:

  A character vector specifying the annotation sources to be included in
  the gene annotation databases. Can be one or more of
  `"GO", "GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway", "Reactome", "CORUM", "MP", "DO", "HPO", "PFAM", "CSPA", "Surfaceome", "SPRomeDB", "VerSeDa", "TFLink", "hTFtarget", "TRRUST", "JASPAR", "ENCODE", "MSigDB", "CellTalk", "CellChat", "Chromosome", "GeneType", "Enzyme", "TF", "CytoTRACE2"`.
  MSigDB subcollections can be requested as `"MSigDB_<collection>"`,
  such as `"MSigDB_H"` for human Hallmark and `"MSigDB_MH"` for mouse
  Hallmark. Note: `"CytoTRACE2"` is species-independent and downloads
  pre-trained model data required by
  [RunCytoTRACE](https://mengxu98.github.io/scop/reference/RunCytoTRACE.md).

- termnames:

  A vector of term names to be used from the database. Default is
  `NULL`, in which case all features from the database are used.

- db_update:

  Whether the gene annotation databases should be forcefully updated. If
  set to FALSE, the function will attempt to load the cached databases
  instead. Default is `FALSE`.

- db_version:

  A character vector specifying the version of the gene annotation
  databases to be retrieved. Default is `"latest"`.

- convert_species:

  Whether to use a species-converted database when the annotation is
  missing for the specified species. Default is `TRUE`.

- Ensembl_version:

  An integer specifying the Ensembl version. Default is `NULL`. If
  `NULL`, the latest version will be used.

- mirror:

  Specify an Ensembl mirror to connect to. The valid options here are
  `"www"`, `"uswest"`, `"useast"`, `"asia"`.

- minGSSize:

  The minimum size of a gene set to be considered in the enrichment
  analysis.

- maxGSSize:

  The maximum size of a gene set to be considered in the enrichment
  analysis.

- method:

  The method to use for scoring. Can be `"Seurat"`, `"AUCell"`,
  `"UCell"`, `"GSVA"`, `"ssGSEA"`, `"zscore"`, `"PLAGE"`, or `"VISION"`.
  Multiple methods can be supplied at once; in that case each method is
  run separately and stored with a method suffix such as `"GO_AUCell"`
  or `"GO_GSVA"`. Default is `"Seurat"`.

- backend:

  Scoring backend. `"cpp"` is the default for supported methods. `"r"`
  uses the original package implementation. `"cpp"` currently supports
  `method = "Seurat"`, `method = "AUCell"`, `method = "GSVA"`,
  `method = "ssGSEA"`, `method = "zscore"`, and `method = "PLAGE"`.
  `method = "UCell"` and `method = "VISION"` fall back to `"r"` when
  `backend` is not explicitly set.

- cpp_strategy:

  AUCell scoring strategy used when `backend = "cpp"`. `"sparse"` ranks
  non-zero genes and approximates zero ties without densifying the
  expression matrix. `"aucell"` calls the official
  [`AUCell::AUCell_buildRankings()`](https://rdrr.io/pkg/AUCell/man/AUCell_buildRankings.html)
  and
  [`AUCell::AUCell_calcAUC()`](https://rdrr.io/pkg/AUCell/man/AUCell_calcAUC.html)
  path for exact consistency with the R backend, `"topk"` ranks only
  genes that can contribute to AUCell AUC, and `"full"` ranks all genes.

- classification:

  Whether to perform classification based on the scores. Default is
  `TRUE`.

- name:

  The name of the assay to store the scores in. Only used if new_assay
  is TRUE. Default is `""`.

- new_assay:

  Whether to create a new assay for storing the scores. Default is
  `FALSE`.

- store_metadata:

  Whether to also store score columns in `meta.data`. When `NULL`,
  manual `features = list(...)` input is stored in `meta.data` by
  default, while database-derived results stay assay-only when
  `new_assay = TRUE`.

- seed:

  Random seed for reproducibility. Default is `11`.

- cores:

  The number of cores to use for parallelization with
  [foreach::foreach](https://rdrr.io/pkg/foreach/man/foreach.html).
  Default is `1`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- ...:

  Additional arguments to be passed to the scoring methods.

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
features_all <- rownames(pancreas_sub)
pancreas_sub <- CellScoring(
  pancreas_sub,
  features = list(
    A = features_all[1:100],
    B = features_all[101:200]
  ),
  method = "AUCell",
  name = "test"
)
CellDimPlot(pancreas_sub, "test_classification")

FeatureDimPlot(
  pancreas_sub,
  features = "test_A"
)

pancreas_sub <- CellScoring(
  pancreas_sub,
  features = list(A = features_all[1:100]),
  method = c("AUCell", "GSVA")
)
FeatureStatPlot(
  pancreas_sub,
  stat.by = c("AUCell_A", "GSVA_A"),
  group.by = "CellType",
  plot.by = "feature",
  plot_type = "violin",
  stack = TRUE
)

FeatureDimPlot(
  pancreas_sub,
  features = c("AUCell_A", "GSVA_A"),
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

GroupHeatmap(
  pancreas_sub,
  features = c("AUCell_A", "GSVA_A", "Sox9", "Anxa2", "Bicc1"),
  group.by = "CellType"
)

data(panc8_sub)
panc8_sub <- integration_scop(
  panc8_sub,
  batch = "tech",
  integration_method = "Harmony"
)

panc8_sub <- CellScoring(
  panc8_sub,
  layer = "data",
  assay = "RNA",
  db = "GO_BP",
  species = "Homo_sapiens",
  minGSSize = 10,
  maxGSSize = 100,
  method = "AUCell",
  name = "GO",
  new_assay = TRUE
)

panc8_sub <- integration_scop(
  panc8_sub,
  assay = "GO",
  batch = "tech",
  integration_method = "Harmony"
)
CellDimPlot(
  panc8_sub,
  group.by = c("tech", "celltype")
)

pancreas_sub <- CellScoring(
  pancreas_sub,
  layer = "data",
  assay = "RNA",
  db = "GO_BP",
  species = "Mus_musculus",
  termnames = panc8_sub[["GO"]]@meta.features[, "termnames"],
  method = "AUCell",
  name = "GO",
  new_assay = TRUE
)
pancreas_sub <- standard_scop(
  pancreas_sub,
  assay = "GO"
)

pancreas_sub[["tech"]] <- "Mouse"
panc_merge <- integration_scop(
  srt_list = list(panc8_sub, pancreas_sub),
  assay = "GO",
  batch = "tech",
  integration_method = "Harmony"
)
CellDimPlot(
  srt = panc_merge,
  group.by = c("tech", "celltype", "SubCellType", "Phase")
)

genenames <- make.unique(
  thisutils::capitalize(
    rownames(panc8_sub[["RNA"]]),
    force_tolower = TRUE
  )
)
names(genenames) <- rownames(panc8_sub)
panc8_sub <- RenameFeatures(
  panc8_sub,
  newnames = genenames,
  assay = "RNA"
)
panc_merge <- integration_scop(
  srt_list = list(panc8_sub, pancreas_sub),
  assay = "RNA",
  batch = "tech",
  integration_method = "Harmony"
)
CellDimPlot(
  srt = panc_merge,
  group.by = c("tech", "celltype", "SubCellType", "Phase")
)
```
