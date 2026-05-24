# Run DoRothEA transcription factor activity inference

Run DoRothEA transcription factor activity inference

## Usage

``` r
RunDorothea(
  srt,
  assay = NULL,
  layer = "data",
  species = c("Homo_sapiens", "Mus_musculus"),
  input_species = NULL,
  geneID_from_IDtype = "symbol",
  geneID_to_IDtype = "symbol",
  homolog_params = list(),
  confidence = c("A", "B", "C"),
  regulons = NULL,
  method = c("ulm", "viper", "wmean"),
  minsize = 5,
  options = list(),
  assay_name = "dorothea",
  new_assay = TRUE,
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- assay:

  Which assay to use. If `NULL`, the default assay of the Seurat object
  will be used. When the object also contains `ChromatinAssay`, the
  default assay and additional `ChromatinAssay` will be preprocessed
  sequentially.

- layer:

  Assay layer used as the expression matrix.

- species:

  Species used to select bundled DoRothEA regulons. DoRothEA only
  provides human and mouse regulons. For other input species, set
  `input_species` and project expression values to this regulon species
  through homologous gene conversion before activity inference.

- input_species:

  Species of the input expression features. If `NULL`, the input is
  assumed to use the same gene namespace as `species`. When this differs
  from `species`, expression features are converted with
  [ConvertHomologs](https://mengxu98.github.io/scop/reference/ConvertHomologs.md)
  before DoRothEA activity inference.

- geneID_from_IDtype, geneID_to_IDtype:

  Gene identifier types passed to
  [ConvertHomologs](https://mengxu98.github.io/scop/reference/ConvertHomologs.md)
  for cross-species projection. For bundled DoRothEA regulons,
  `geneID_to_IDtype` should normally remain `"symbol"`.

- homolog_params:

  Additional named arguments passed to
  [ConvertHomologs](https://mengxu98.github.io/scop/reference/ConvertHomologs.md)
  when `input_species` differs from `species`, such as
  `Ensembl_version`, `biomart`, `mirror`, `max_tries`, `multi_mapping`,
  and `collapse_fun`.

- confidence:

  DoRothEA confidence levels to keep.

- regulons:

  Optional regulon table with `tf`, `target`, `mor`, and `confidence`
  columns. If `NULL`, bundled `dorothea_hs` or `dorothea_mm` data are
  loaded from the `dorothea` package.

- method:

  Activity inference backend from `decoupleR`.

- minsize:

  Minimum regulon size passed to `decoupleR`.

- options:

  Additional named options passed to the selected `decoupleR` function.

- assay_name:

  Name of the assay used to store TF activity scores.

- new_assay:

  Whether to store TF activity scores as a new assay.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A `Seurat` object with DoRothEA results stored in
`srt@tools[["Dorothea"]]`. For cross-species runs, the homolog
projection summary is stored in
`srt@tools[["Dorothea"]]$homolog_conversion`.

## References

Garcia-Alonso, L., Holland, C.H., Ibrahim, M.M., Turei, D., and
Saez-Rodriguez, J. (2019). Benchmark and integration of resources for
the estimation of human transcription factor activities. *Genome
Research*, 29, 1363-1375.
[doi:10.1101/gr.240663.118](https://doi.org/10.1101/gr.240663.118)

Badia-i-Mompel, P., Velez Santiago, J., Braunger, J., Geiss, C.,
Dimitrov, D., Muller-Dott, S., Taus, P., Dugourd, A., Holland, C.H.,
Ramirez Flores, R.O., and Saez-Rodriguez, J. (2022). decoupleR: ensemble
of computational methods to infer biological activities from omics data.
*Bioinformatics Advances*, 2, vbac016.
[doi:10.1093/bioadv/vbac016](https://doi.org/10.1093/bioadv/vbac016)

## Examples

``` r
data(pancreas_sub)
pancreas_sub <- standard_scop(
  pancreas_sub,
  verbose = FALSE
)
#> ℹ [2026-05-24 16:03:22] Skip `log1p()` because `layer = data` is not "counts"

pancreas_sub <- RunDorothea(
  pancreas_sub,
  layer = "counts",
  species = "Mus_musculus",
  confidence = c("A", "B", "C"),
  method = "ulm",
  minsize = 5,
  new_assay = FALSE
)
#> ℹ [2026-05-24 16:03:32] Run DoRothEA/decoupleR with 12895 regulon edges
#> ℹ [2026-05-24 16:03:42] DoRothEA TF activity scores stored in <Seurat> metadata

pancreas_sub@tools$Dorothea$regulon_summary
#>   n_tfs n_targets n_edges confidence
#> 1   273      5130   12895      A,B,C
head(pancreas_sub@tools$Dorothea$result)
#>   statistic source        condition     score      p_value
#> 1       ulm     Ar AAACCTGAGCCTTGAT 1.9230738 0.0544885134
#> 2       ulm     Ar AAACCTGGTAAGTGGC 3.8375556 0.0001247429
#> 3       ulm     Ar AAACGGGAGATATGGT 2.3175899 0.0204841737
#> 4       ulm     Ar AAACGGGCAAAGAATC 0.9565534 0.3388071679
#> 5       ulm     Ar AAACGGGGTACAGTTC 0.6864447 0.4924426684
#> 6       ulm     Ar AAACGGGTCAGCTCTC 1.0246312 0.3055527415

activity_cols <- head(
  grep("^dorothea_", colnames(pancreas_sub@meta.data), value = TRUE),
  2
)
head(pancreas_sub@meta.data[, activity_cols, drop = FALSE])
#>                  dorothea_Ar dorothea_Arid2
#> AAACCTGAGCCTTGAT   1.9230738     -0.5482112
#> AAACCTGGTAAGTGGC   3.8375556     -0.2029179
#> AAACGGGAGATATGGT   2.3175899     -0.1763966
#> AAACGGGCAAAGAATC   0.9565534     -0.5784380
#> AAACGGGGTACAGTTC   0.6864447     -0.2284842
#> AAACGGGTCAGCTCTC   1.0246312     -0.6278209

FeatureDimPlot(
  pancreas_sub,
  features = activity_cols,
  reduction = "StandardUMAP2D",
  ncol = 2
)
```
