# Add Giotto results back to Seurat

Add Giotto results back to Seurat

## Usage

``` r
AddGiottoToSeurat(
  srt,
  x,
  result = c("cluster", "hmrf"),
  name = NULL,
  tool_name = "Giotto",
  store_result = TRUE
)
```

## Arguments

- srt:

  A Seurat object.

- x:

  A \`giotto2\` workflow object.

- result:

  Giotto result to copy back.

- name:

  Metadata column name to write. If \`NULL\`, a default name is used.

- tool_name:

  Name used to store the Giotto workflow object in \`srt@tools\`.

- store_result:

  Whether to store the Giotto workflow object in
  \`srt@tools\[\[tool_name\]\]\`.

## Value

A Seurat object.
