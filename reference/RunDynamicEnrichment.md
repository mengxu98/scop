# RunDynamicEnrichment

This function calculates gene-set scores from the specified database
(`db`) for each lineage using the specified scoring method
(`score_method`). It then treats these scores as expression values and
uses them as input to the RunDynamicFeatures function to identify
dynamically enriched terms along the lineage.

## Usage

``` r
RunDynamicEnrichment(
  srt,
  lineages,
  score_method = "AUCell",
  layer = "data",
  assay = NULL,
  min_expcells = 20,
  r.sq = 0.2,
  dev.expl = 0.2,
  padjust = 0.05,
  IDtype = "symbol",
  species = "Homo_sapiens",
  db = "GO_BP",
  db_update = FALSE,
  db_version = "latest",
  convert_species = TRUE,
  Ensembl_version = NULL,
  mirror = NULL,
  TERM2GENE = NULL,
  TERM2NAME = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  cores = 1,
  verbose = TRUE,
  seed = 11
)
```

## Arguments

- srt:

  A Seurat object containing the results of differential expression
  analysis (RunDEtest). If specified, the genes and groups will be
  extracted from the Seurat object automatically. If not specified, the
  `geneID` and `geneID_groups` arguments must be provided.

- lineages:

  A character vector specifying the lineages to plot.

- score_method:

  The method to use for scoring. Can be `"Seurat"`, `"AUCell"`, or
  `"UCell"`. Default is `"Seurat"`.

- layer:

  Which layer to use. Default is `"counts"`.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used.

- min_expcells:

  The minimum number of expected cells. Default is `20`.

- r.sq:

  The R-squared threshold. Default is `0.2`.

- dev.expl:

  The deviance explained threshold. Default is `0.2`.

- padjust:

  The p-value adjustment threshold. Default is `0.05`.

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
  `"GO", "GO_BP", "GO_CC", "GO_MF", "KEGG", "WikiPathway", "Reactome", "CORUM", "MP", "DO", "HPO", "PFAM", "CSPA", "Surfaceome", "SPRomeDB", "VerSeDa", "TFLink", "hTFtarget", "TRRUST", "JASPAR", "ENCODE", "MSigDB", "CellTalk", "CellChat", "Chromosome", "GeneType", "Enzyme", "TF"`.

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

- TERM2GENE:

  A data frame specifying the gene-term mapping for a custom database.
  The first column should contain the term IDs, and the second column
  should contain the gene IDs.

- TERM2NAME:

  A data frame specifying the term-name mapping for a custom database.
  The first column should contain the term IDs, and the second column
  should contain the corresponding term names.

- minGSSize:

  The minimum size of a gene set to be considered in the enrichment
  analysis.

- maxGSSize:

  The maximum size of a gene set to be considered in the enrichment
  analysis.

- cores:

  The number of cores to use for parallelization with
  [foreach::foreach](https://rdrr.io/pkg/foreach/man/foreach.html).
  Default is `1`.

- verbose:

  Whether to print the message. Default is `TRUE`.

- seed:

  Random seed for reproducibility. Default is `11`.

## See also

[RunDynamicFeatures](https://mengxu98.github.io/scop/reference/RunDynamicFeatures.md),
[DynamicHeatmap](https://mengxu98.github.io/scop/reference/DynamicHeatmap.md)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(pancreas_sub)
pancreas_sub <- RunSlingshot(
  pancreas_sub,
  group.by = "SubCellType",
  reduction = "UMAP"
)

pancreas_sub <- RunDynamicFeatures(
  pancreas_sub,
  lineages = "Lineage1",
  n_candidates = 200
)
#> ⠙ [2026-01-20 07:52:27] Running for Gcg [1/199] ■                              …
#> ⠹ [2026-01-20 07:52:27] Running for Hmgb2 [25/199] ■■■■■                       …
#> ⠸ [2026-01-20 07:52:27] Running for Ptn [110/199] ■■■■■■■■■■■■■■■■■■           …
#> ⠼ [2026-01-20 07:52:27] Running for Notch2 [193/199] ■■■■■■■■■■■■■■■■■■■■■■■■■■…
#> ✔ [2026-01-20 07:52:27] Completed 199 tasks in 7s
#> 
ht1 <- DynamicHeatmap(
  pancreas_sub,
  lineages = "Lineage1",
  cell_annotation = "SubCellType",
  n_split = 4
)

ht1$plot


pancreas_sub <- RunDynamicEnrichment(
  pancreas_sub,
  lineages = "Lineage1",
  score_method = "UCell",
  db = "GO_BP",
  species = "Mus_musculus"
)
#> ⠙ [2026-01-20 07:55:47] Running for GO-BP-2..deoxyribonucleotide.biosynthetic.p…
#> ⠹ [2026-01-20 07:55:47] Running for GO-BP-amyloid.precursor.protein.metabolic.p…
#> ⠸ [2026-01-20 07:55:47] Running for GO-BP-coenzyme.A.metabolic.process [442/280…
#> ⠼ [2026-01-20 07:55:47] Running for GO-BP-hydrogen.peroxide.catabolic.process […
#> ⠴ [2026-01-20 07:55:47] Running for GO-BP-myelin.maintenance [1068/2800] ■■■■■■…
#> ⠦ [2026-01-20 07:55:47] Running for GO-BP-neurotransmitter.receptor.localizatio…
#> ⠧ [2026-01-20 07:55:47] Running for GO-BP-positive.regulation.of.microtubule.po…
#> ⠇ [2026-01-20 07:55:47] Running for GO-BP-regulation.of.amyloid.beta.clearance …
#> ⠏ [2026-01-20 07:55:47] Running for GO-BP-regulation.of.heart.growth [2154/2800…
#> ⠋ [2026-01-20 07:55:47] Running for GO-BP-renal.system.process [2462/2800] ■■■■…
#> ⠙ [2026-01-20 07:55:47] Running for GO-BP-ventricular.cardiac.muscle.cell.diffe…
#> ✔ [2026-01-20 07:55:47] Completed 2800 tasks in 28.6s
#> 
ht2 <- DynamicHeatmap(
  pancreas_sub,
  assay = "GO_BP",
  lineages = "Lineage1_GO_BP",
  cell_annotation = "SubCellType",
  n_split = 4,
  split_method = "kmeans-peaktime"
)

ht2$plot
```
