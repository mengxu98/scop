# Run SecAct cell-cell communication analysis

Run SecAct cell-cell communication modules for Seurat scRNA-seq input or
SpaCET single-cell-resolution spatial transcriptomics input.

## Usage

``` r
RunSecActCCC(
  srt = NULL,
  inputProfile = NULL,
  mode = c("scRNAseq", "scST"),
  cellType_meta,
  condition_meta = NULL,
  conditionCase = NULL,
  conditionControl = NULL,
  scale.factor = NULL,
  act_diff_cutoff = 2,
  exp_logFC_cutoff = 0.2,
  exp_mean_all_cutoff = 2,
  exp_fraction_case_cutoff = 0.1,
  padj_cutoff = 0.01,
  sigMatrix = "SecAct",
  is.group.sig = TRUE,
  is.group.cor = 0.9,
  lambda = 5e+05,
  nrand = 1000,
  radius = 20,
  ratio_cutoff = 0.2,
  coreNo = 6,
  tool_name = "SecAct_CCC",
  store_results = TRUE,
  verbose = TRUE
)
```

## Arguments

- srt:

  Optional Seurat object. When `mode = "scRNAseq"`, this is passed to
  `SecAct.activity.inference.scRNAseq`. When `mode = "matrix"`,
  expression is extracted from `srt` if `inputProfile` is not supplied.

- inputProfile:

  Expression matrix, Seurat object, or SpaCET object. Matrix input must
  be genes x samples, cells, or spatial spots.

- mode:

  `scRNAseq` or `scST`.

- cellType_meta:

  Metadata column containing cell type or state labels for
  `mode = "scRNAseq"` when `is.singleCellLevel = FALSE`.

- condition_meta:

  Metadata column containing condition labels for scRNA-seq CCC.

- conditionCase, conditionControl:

  Case and control labels for scRNA-seq CCC.

- scale.factor:

  Spot-level scale factor passed to `SecAct.activity.inference.ST`.

- act_diff_cutoff, exp_logFC_cutoff, exp_mean_all_cutoff,
  exp_fraction_case_cutoff, padj_cutoff:

  Cutoffs passed to `SecAct.CCC.scRNAseq`.

- sigMatrix:

  SecAct signature matrix name.

- is.group.sig, is.group.cor, lambda, nrand:

  Parameters passed to SecAct signature grouping and randomization
  routines.

- radius, ratio_cutoff, coreNo:

  Parameters passed to `SecAct.CCC.scST`.

- tool_name:

  Name used in `srt@tools`.

- store_results:

  Whether to store raw SecAct results in `srt@tools`.

- verbose:

  Whether to print the message. Default is `TRUE`.

## Value

A Seurat or SpaCET object with SecAct CCC results.
