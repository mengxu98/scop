# Convert a Seurat object to an `.h5ad` file

Convert a Seurat object to an `.h5ad` file

## Usage

``` r
srt_to_h5ad(
  object,
  path,
  features = NULL,
  assay_x = "RNA",
  layer_x = "counts",
  assay_y = c("spliced", "unspliced"),
  layer_y = "counts",
  reductions = NULL,
  graphs = NULL,
  neighbors = NULL,
  convert_tools = FALSE,
  convert_misc = FALSE,
  overwrite = FALSE,
  verbose = TRUE,
  srt = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- path:

  Path to the output `.h5ad` file.

- features:

  Optional vector of features to include in the anndata object. Default
  is all features in `assay_x`.

- assay_x:

  Assay to convert as the main data matrix in the anndata object.

- layer_x:

  Layer name for assay_x in the Seurat object.

- assay_y:

  Assays, or extra layers of `assay_x`, to convert as layers in the
  anndata object. Each name is matched against the assays of the object
  first, then against the layers of `assay_x`, which covers velocity
  matrices such as `spliced` and `unspliced` stored as layers of a
  Seurat v5 `RNA` assay.

- layer_y:

  Layer names for the assays in `assay_y`. It is ignored for names that
  are matched as layers of `assay_x`, where the layer name itself is
  used.

- reductions:

  Character vector specifying which Seurat reductions to convert into
  `obsm`. Default is `NULL`, which converts all available reductions.

- graphs:

  Character vector specifying which Seurat graphs to convert into
  `obsp`. Default is `NULL`, which converts all available graphs.

- neighbors:

  Character vector specifying which Seurat neighbor objects to convert
  into `obsp`. Default is `NULL`, which converts all available neighbor
  objects.

- convert_tools:

  Whether to convert the tool-specific data.

- convert_misc:

  Whether to convert the miscellaneous data.

- overwrite:

  Whether to overwrite an existing file.

- verbose:

  Whether to print the message. Default is `TRUE`.

- srt:

  Deprecated alias for `object`; supply exactly one of the two. It will
  be removed in scop 1.0.0.

## Value

Invisibly returns the normalized `path` of the written file.

## See also

[srt_to_adata](https://mengxu98.github.io/scop/reference/srt_to_adata.md),
[h5ad_to_srt](https://mengxu98.github.io/scop/reference/h5ad_to_srt.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(pancreas_sub)
srt_to_h5ad(pancreas_sub, "pancreas_sub.h5ad")
} # }
```
