# NMF similarity heatmap

NMF similarity heatmap

## Usage

``` r
NMFHeatmap(
  srt,
  plot_type = c("cells", "features"),
  reduction = "nmf",
  dims = NULL,
  cells = NULL,
  features = NULL,
  similarity_metric = "cosine",
  cell_annotation = NULL,
  feature_annotation = NULL,
  assay = NULL,
  border = TRUE,
  heatmap_border = NULL,
  cell_annotation_border = NULL,
  feature_annotation_border = NULL,
  heatmap_border_palcolor = "black",
  cell_annotation_border_palcolor = "black",
  feature_annotation_border_palcolor = "black",
  heatmap_border_size = 1,
  cell_annotation_border_size = 1,
  feature_annotation_border_size = 1,
  show_row_names = FALSE,
  row_names_wrap = NULL,
  show_column_names = FALSE,
  row_names_side = "left",
  column_names_side = "top",
  row_names_rot = 0,
  column_names_rot = 90,
  row_title = NULL,
  column_title = NULL,
  anno_terms = FALSE,
  anno_keys = FALSE,
  anno_features = FALSE,
  terms_width = grid::unit(4, "in"),
  terms_stat_width = grid::unit(1.35, "in"),
  terms_fontsize = 8,
  terms_stat = "none",
  terms_stat_digits = 2,
  terms_stat_label = "value",
  terms_stat_axis = FALSE,
  terms_stat_background_palcolor = NULL,
  terms_stat_border = NULL,
  terms_stat_border_palcolor = NULL,
  terms_stat_border_size = NULL,
  terms_stat_label_palcolor = NULL,
  terms_group_background = FALSE,
  terms_background_palcolor = "grey98",
  terms_background_alpha = 1,
  terms_border = TRUE,
  terms_border_palcolor = "black",
  terms_border_size = 0.8,
  terms_text_palcolor = NULL,
  terms_bar_palcolor = NULL,
  keys_width = grid::unit(2, "in"),
  keys_fontsize = c(6, 10),
  features_width = grid::unit(2, "in"),
  features_fontsize = c(6, 10),
  IDtype = "symbol",
  species = "Homo_sapiens",
  db_update = FALSE,
  db_version = "latest",
  db_combine = FALSE,
  convert_species = FALSE,
  Ensembl_version = NULL,
  mirror = NULL,
  db = "GO_BP",
  TERM2GENE = NULL,
  TERM2NAME = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  GO_simplify = FALSE,
  GO_simplify_cutoff = "p.adjust < 0.05",
  simplify_method = "Wang",
  simplify_similarityCutoff = 0.7,
  pvalueCutoff = NULL,
  padjustCutoff = 0.05,
  topTerm = 5,
  show_termid = FALSE,
  topWord = 20,
  words_excluded = NULL,
  heatmap_palette = "simspec",
  heatmap_palcolor = c("#ffffe5", "#d9f0d3", "#74add1", "#2166ac"),
  heatmap_limits = NULL,
  cluster_palette = "simspec",
  cluster_palcolor = NULL,
  cell_annotation_palette = "Chinese",
  cell_annotation_palcolor = NULL,
  feature_annotation_palette = "Dark2",
  feature_annotation_palcolor = NULL,
  use_raster = NULL,
  raster_device = "png",
  raster_by_magick = FALSE,
  height = NULL,
  width = NULL,
  units = "inch",
  cores = 1,
  seed = 11,
  legend.position = "right",
  ht_params = list(),
  verbose = TRUE
)
```

## Arguments

- srt:

  A Seurat object containing an NMF dimensional reduction.

- plot_type:

  Plot type. `"cells"` plots cell/spot similarity from NMF embeddings.
  `"features"` plots feature similarity from NMF loadings.

- reduction:

  Name of the NMF reduction.

- dims:

  Dimensions/components from the NMF reduction to use. If `NULL`, all
  available dimensions are used.

- cells:

  Cells/spots to include when `plot_type = "cells"`.

- features:

  Features to include when `plot_type = "features"`. If `NULL`, variable
  features shared with the loading matrix are used; if none are found,
  all features in the loading matrix are used.

- similarity_metric:

  Similarity metric.

- cell_annotation:

  Metadata columns to show as column annotations in cell mode.

- feature_annotation:

  Feature metadata columns to show as column annotations in feature
  mode.

- assay:

  Assay to use. `NULL` uses the default assay.

- border:

  Draw borders. Kept for compatibility; more specific `*_border`
  arguments inherit this when `NULL`.

- heatmap_border, cell_annotation_border, feature_annotation_border:

  Borders for the heatmap body, annotations, and their matching legends.
  `NULL` inherits `border`.

- heatmap_border_palcolor, cell_annotation_border_palcolor,
  feature_annotation_border_palcolor:

  Border colors for the matching body, annotation, and legend.

- heatmap_border_size, cell_annotation_border_size,
  feature_annotation_border_size:

  Border line widths for the matching body, annotation, and legend.

- show_row_names:

  Whether to draw row/column names for the heatmap body.

- row_names_wrap:

  Maximum number of characters per displayed row-name line. When set to
  a positive number, underscores are displayed as spaces and labels are
  wrapped without changing the underlying item identifiers. `NULL`
  disables wrapping.

- show_column_names:

  Whether to draw row/column names for the heatmap body.

- row_names_side, column_names_side, row_names_rot, column_names_rot:

  Name placement.

- row_title:

  The title for the row names in the heatmap. If not provided, the
  default is to use the query grouping variable.

- column_title:

  The title for the column names in the heatmap. Default is to use the
  reference grouping variable.

- anno_terms, anno_keys, anno_features:

  Enrichment annotations.

- terms_width, terms_stat_width, terms_fontsize:

  Term annotation size.

- terms_stat:

  Enrichment statistic for term bars: `"none"`, `"score"` (`-log10` of
  the active p-value), or an enrichment column such as `"p.adjust"`.

- terms_stat_digits, terms_stat_label, terms_stat_axis:

  Statistic labels (`"none"`, `"value"`, `"significance"`, `"both"`) and
  shared axis.

- terms_stat_background_palcolor, terms_stat_border,
  terms_stat_border_palcolor, terms_stat_border_size,
  terms_stat_label_palcolor:

  Statistic-panel appearance. `NULL` inherits the matching `terms_*`
  setting.

- terms_group_background, terms_background_palcolor,
  terms_background_alpha, terms_border, terms_border_palcolor,
  terms_border_size, terms_text_palcolor, terms_bar_palcolor:

  Term-block appearance. `terms_text_palcolor = NULL` maps text to
  enrichment significance; `terms_bar_palcolor = NULL` matches bar color
  to term text.

- keys_width, keys_fontsize, features_width, features_fontsize:

  Key and feature annotations.

- IDtype, species, db_combine, mirror, db, TERM2GENE, TERM2NAME,
  minGSSize, maxGSSize:

  Gene-set database (see
  [PrepareDB](https://mengxu98.github.io/scop/reference/PrepareDB.md)).

- db_update:

  Force a refresh. `FALSE` loads the cache when available.

- db_version:

  Database version to retrieve.

- convert_species:

  Use a species-converted database when the annotation is missing for
  `species`.

- Ensembl_version:

  Ensembl version. `NULL` uses the latest.

- GO_simplify, GO_simplify_cutoff, simplify_method,
  simplify_similarityCutoff:

  GO simplification.

- pvalueCutoff, padjustCutoff, topTerm, show_termid, topWord,
  words_excluded:

  Enrichment filters.

- heatmap_palette:

  Palette used for CNV heatmap values.

- heatmap_palcolor:

  Palette used for CNV heatmap values.

- heatmap_limits:

  Numeric breaks for the heatmap color scale. If `NULL`, defaults to
  `c(0, 0.35, 0.75, 1)`.

- cluster_palette:

  Palette used for NMF cluster/program annotations.

- cluster_palcolor:

  Optional custom colors for NMF cluster/program annotations.

- cell_annotation_palette:

  Color palette for cell-type annotations.

- cell_annotation_palcolor:

  Custom colors for cell-type annotations.

- feature_annotation_palette:

  Color palette for feature annotations.

- feature_annotation_palcolor:

  Custom colors for feature annotations.

- use_raster, raster_device, raster_by_magick:

  Raster device (`NULL` chooses automatically).

- width, height, units:

  Heatmap size. `NULL` sizes from matrix dimensions.

- cores:

  The number of worker processes to use for parallelization. Default is
  `1`.

- seed:

  Optional integer seed. When supplied, every input receives a
  deterministic independent L'Ecuyer-CMRG random-number stream, making
  results reproducible across worker counts and scheduling order. The
  caller's random number state is restored when the call finishes.

- legend.position:

  Legend side (`"right"`, `"left"`, `"top"`, `"bottom"`). Gap to the
  heatmap grows automatically when long row names are on the right.

- ht_params:

  Extra arguments passed to
  [ComplexHeatmap::Heatmap](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html),
  overriding defaults.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A list with the following elements:

- `plot`: The heatmap plot as a patchwork/ggplot object.

- `similarity_matrix`: The ordered similarity matrix used for plotting.

- `nmf_cluster`: The ordered NMF cluster/program assignment.

- `order`: The ordered row/column names.

- `metadata`: Ordered cell or feature metadata used for annotations.

- `enrichment`: Enrichment results for feature mode when requested,
  otherwise `NULL`.

## See also

[RunNMF](https://mengxu98.github.io/scop/reference/RunNMF.md)

## Examples

``` r
library(Matrix)
data(pancreas_sub)
pancreas_sub <- NormalizeData(pancreas_sub)
pancreas_sub <- FindVariableFeatures(
  pancreas_sub,
  nfeatures = 1000
)
pancreas_sub <- RunNMF(
  pancreas_sub,
  features = SeuratObject::VariableFeatures(pancreas_sub),
  nbes = 5,
  maxit = 50
)
#> ℹ [2026-09-06 21:39:22] Running NMF...
#> ℹ BE_ 1 
#> ℹ Positive:  Spp1, Clu, Ttr, Krt18, Ptma, Rpl12, Sparc, Dbi, Gapdh, Mt1 
#> ℹ      Cd24a, Mgst1, H19, Pebp1, Myl12a, Cldn3, Clps, Atp1b1, Sox4, Gnas 
#> ℹ      Vim, Ambp, Cdkn1c, Jun, Mdk, Serpinh1, Eno1, Anxa2, Acot1, Tmsb4x 
#> ℹ Negative:  Fgf8, Mapk12, Lrrc6, Spock1, Lrrc9, Fam71b, Il1r2, Serpini1, Gng4, Cdca2 
#> ℹ      Sulf2, Pgf, Dusp26, Ucn3, Entpd3, Gm13373, Megf11, Acvr1c, Krtap16-1, Nrp2 
#> ℹ      Mmel1, Pabpn1l, Sept3, Hepacam2, Rnf138rt1, Scn9a, Tex36, Syt13, Bace2, Igsf21 
#> ℹ BE_ 2 
#> ℹ Positive:  Gnas, Pyy, Rbp4, Chgb, Chga, Cpe, Slc25a5, Hmgn3, Ttr, Pcsk1n 
#> ℹ      Bex2, Isl1, Aplp1, Fam183b, Rap1b, Glud1, Lrpprc, Fev, Slc38a5, Mid1ip1 
#> ℹ      Akr1c19, Cck, Gch1, Tm4sf4, Ptma, Map1b, Sec61b, Clps, Tuba1a, 1700086L19Rik 
#> ℹ Negative:  Fam124a, 1810034E14Rik, 5730507C01Rik, Tmem100, Fam71b, Sycp3, Fscn1, Cdca2, Traip, Gm8113 
#> ℹ      C2cd4c, Plpp2, Sulf2, Ucn3, Col1a1, Megf11, Bcl2, Gm28875, Ugt2b35, Ugt2b36 
#> ℹ      A730098A19Rik, Serpinb6b, Eya2, AA986860, Palmd, Vps8, Crybb1, Pabpn1l, Il18, Gjb1 
#> ℹ BE_ 3 
#> ℹ Positive:  Tmsb4x, Neurog3, Mdk, Cck, Sox4, Btg2, Btbd17, Gadd45a, Ptma, Selm 
#> ℹ      Krt7, Gnas, Hn1, Cdkn1a, Hes6, Cd24a, Smarcd2, Camk2n1, Rpl12, Cotl1 
#> ℹ      Cldn6, Map1b, Clps, Aplp1, Tubb3, Pax4, Slc25a5, Gpx2, Igfbpl1, Nkx6-1 
#> ℹ Negative:  Mapk12, Lrrc6, Ccl28, Spock1, Hoxb2, Tmem100, Il1r2, Wdr86, Angptl4, Cdca2 
#> ℹ      Traip, Slc16a10, Dusp26, Ucn3, Entpd3, Gmfg, Acvr1c, Col27a1, Kctd8, Mdm1 
#> ℹ      Adora2b, C530044C16Rik, Gm28875, Ugt2b35, Ugt2b36, Fosb, A730098A19Rik, Serpinb6b, Pax6os1, Nrsn1 
#> ℹ BE_ 4 
#> ℹ Positive:  Iapp, Pyy, Nnat, Ins2, Ins1, Rbp4, Gnas, Ttr, Dlk1, Sec61b 
#> ℹ      Pcsk2, Calr, Ppp1r1a, Hspa5, Pdia6, Sdf2l1, Hsp90b1, Gng12, Tuba1a, Pcsk1n 
#> ℹ      Hadh, Cpe, Clps, Ptma, Mafb, Chgb, Scg2, Gapdh, Fkbp2, Chga 
#> ℹ Negative:  Fgf8, Mapk12, Ccl28, Lrrc9, Hoxb2, Tmem100, Serpini1, Lrrn1, Angptl4, Traip 
#> ℹ      Cmtm3, Pgf, Gm13373, Gmfg, Nrp2, Mmel1, Ugt2b35, Ugt2b36, A730098A19Rik, Serpinb6b 
#> ℹ      Tyrobp, AA986860, Palmd, Lmo4, Gjb1, Pdlim1, Scara3, Tex36, P2ry14, Cdc42ep1 
#> ℹ BE_ 5 
#> ℹ Positive:  Tuba1b, Hmgb2, Tubb5, 2810417H13Rik, Ptma, H2afz, Ran, H2afx, Ranbp1, Tubb4b 
#> ℹ      Birc5, Cks1b, Mif, Slc25a5, H1f0, Spc24, Hn1, Gapdh, Cks2, Mdk 
#> ℹ      Rpl12, Spp1, Cdk1, Dut, Hmgb1, Snrpd1, Anp32b, Ldha, Hspe1, Cdca3 
#> ℹ Negative:  Lrrc6, Ccl28, 1810034E14Rik, Tmem100, Fam71b, Sycp3, Kiss1r, Gng4, Prodh2, Lingo1 
#> ℹ      Mpzl1, Angptl4, Slc16a10, Gm8113, C2cd4c, Dpysl3, Entpd3, Col1a1, Rem2, Acvr1c 
#> ℹ      Bcl2, Kctd8, Adora2b, C530044C16Rik, Gm28875, C1qa, Serpinb6b, Camk2n1, Pax6os1, Nrsn1 
#> ✔ [2026-09-06 21:44:13] NMF compute completed
ht_cells <- NMFHeatmap(
  pancreas_sub,
  plot_type = "cells",
  cell_annotation = "CellType"
)
#> ℹ [2026-09-06 21:44:13] `NMFHeatmap()` input: 1000 cells x 5 NMF dimensions. Computing a 1000 x 1000 similarity matrix (~0.01 GiB dense numeric matrix).
#> ℹ [2026-09-06 21:44:13] Ordering `NMFHeatmap()` rows and columns ...
#> ℹ [2026-09-06 21:44:13] Building ComplexHeatmap object for `NMFHeatmap()` ...
#> ℹ [2026-09-06 21:44:13] Calculating `NMFHeatmap()` render size ...
#> ℹ [2026-09-06 21:44:13] Drawing `NMFHeatmap()`; this can take time for large similarity matrices ...
#> ℹ [2026-09-06 21:44:14] Assembling `NMFHeatmap()` plot object ...
ht_cells$plot


ht_features <- NMFHeatmap(
  pancreas_sub,
  plot_type = "features"
)
#> ℹ [2026-09-06 21:44:15] `NMFHeatmap()` input: 1000 features x 5 NMF dimensions. Computing a 1000 x 1000 similarity matrix (~0.01 GiB dense numeric matrix).
#> ℹ [2026-09-06 21:44:15] Ordering `NMFHeatmap()` rows and columns ...
#> ℹ [2026-09-06 21:44:15] Building ComplexHeatmap object for `NMFHeatmap()` ...
#> ℹ [2026-09-06 21:44:15] Calculating `NMFHeatmap()` render size ...
#> ℹ [2026-09-06 21:44:15] Drawing `NMFHeatmap()`; this can take time for large similarity matrices ...
#> ℹ [2026-09-06 21:44:16] Assembling `NMFHeatmap()` plot object ...
ht_features$plot
```
